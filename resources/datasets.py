from flask_restful import reqparse, Resource
from flask import send_from_directory
import os, pandas

from resources import auth
from models import datasets, genes
from resources.errors import DatasetIdNotFoundError, DatasetIsPrivateError, DatasetGeneIdNotInExpressionError, UserNotAuthenticatedError, KeyNotFoundError

def protectedDataset(datasetId):
    """For many classes here where we check if a dataset is private or not before proceding, this convenience function peforms
    the taks and returns the datasets.Dataset instance.
    """
    try:
        ds = datasets.Dataset(datasetId)
        if ds.isPrivate() and not auth.AuthUser().username():
            raise DatasetIsPrivateError
        return ds
    except datasets.DatasetIdNotFoundError:
        raise DatasetIdNotFoundError

# ----------------------------------------------------------
# Working on a single specified dataset
# ----------------------------------------------------------

class DatasetMetadata(Resource):
    def get(self, datasetId):
        """Return dataset metadata for a dataset with datasetId.
        """
        ds = protectedDataset(datasetId)
        return ds.metadata()

class DatasetSamples(Resource):
    def get(self, datasetId):
        """Return sample table for a dataset with datasetId.
        Parameters:
            orient: same options used in pandas.DataFrame.to_dict(orient=...). Default "records".
            na: string to replace na or null values with. Default "". 
        """
        parser = reqparse.RequestParser()
        parser.add_argument('orient', type=str, required=False, default="records")
        parser.add_argument('na', type=str, required=False, default="")
        parser.add_argument('as_file', type=str, required=False, default="false")  
            # - if set to boolean, any string is parsed as true
        args = parser.parse_args()

        ds = protectedDataset(datasetId)
        
        # These fields are used internally by the model and not useful for API
        hideKeys = ["dataset_id"]
        df = ds.samples()
        if len(df)>0:
            df = df.drop(hideKeys, axis=1).fillna(args.get("na"))

            if args.get('as_file').lower().startswith('t'):
                from tempfile import NamedTemporaryFile
                from flask import send_file
                with NamedTemporaryFile() as temp_file:
                    df.to_csv(temp_file.name, sep='\t')
                    return send_file(temp_file.name, as_attachment=True, attachment_filename="stemformatics_dataset_%s.samples.tsv" % datasetId)
            else:
                if args.get('orient')=='records':  # include index
                    df = df.reset_index()
                return df.to_dict(orient=args.get("orient")) 
        else:
            return {}

class DatasetExpression(Resource):
    def get(self, datasetId):
        """Return expression table for a dataset with datasetId and gene id(s).
        """
        parser = reqparse.RequestParser()
        parser.add_argument('gene_id', type=str, required=False)  # Ensembl ids. May be a comma separated list
        parser.add_argument('key', type=str, required=False, default="raw") # ['raw','genes','cpm']
        parser.add_argument('log2', type=str, required=False, default="False") # will apply numpy.log2(exp+1) if true (only to RNASeq data)
        parser.add_argument('orient', type=str, required=False, default="records")  # from pandas
        parser.add_argument('as_file', type=str, required=False, default="False")  
            # - if set to boolean, any string is parsed as true
        args = parser.parse_args()

        ds = protectedDataset(datasetId)

        if args.get('as_file').lower().startswith('t'):  # file download for entire expression matrix - ignore gene_id
            filename = "stemformatics_dataset_%s.%s.tsv" % (datasetId, args.get('key'))  # user will see this name for download
            if ds.metadata()['platform_type']=='RNASeq' and args.get('key')=='cpm': # we don't have cpm saved on file - calculate it and return it
                from tempfile import NamedTemporaryFile
                from flask import send_file
                with NamedTemporaryFile() as temp_file:
                    df = ds.expressionMatrix(key=args.get('key'))
                    df.to_csv(temp_file.name, sep='\t')
                    return send_file(temp_file.name, as_attachment=True, attachment_filename=filename)
            else:
                filepath = ds.expressionFilePath(args.get('key'))
                return send_from_directory(os.path.dirname(filepath), os.path.basename(filepath), as_attachment=True, attachment_filename=filename)
        else:
            geneIds = args.get('gene_id').split(',') if args.get('gene_id') is not None else []
            df = ds.expressionMatrix(key=args.get('key'), applyLog2=args.get('log2').lower().startswith('t'))
            geneIds = df.index.intersection(geneIds)
            df = df.loc[geneIds]
            
            if args.get('orient')=='records':
                df = df.reset_index()
            if len(df)>0:
                return df.to_dict(orient=args.get('orient'))
            else:
                raise DatasetGeneIdNotInExpressionError

class DatasetPca(Resource):
    def get(self, datasetId):
        """Return pca data for a dataset with datasetId.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('orient', type=str, required=False, default="records")
        parser.add_argument('dims', type=int, required=False, default=20)
        args = parser.parse_args()

        ds = protectedDataset(datasetId)
        coords = ds.pcaCoordinates()
        attributes = ds.pcaAttributes()

        # some datasets (eg 6646) may have really large number of PCA dimensions, but we don't need so many
        coords = coords.iloc[:,:args.get('dims')]
        attributes = attributes.iloc[:args.get('dims'),:]

        return {'coordinates': coords.to_dict(orient=args.get('orient')), 
                'attributes': attributes.to_dict(orient=args.get('orient'))}

# ----------------------------------------------------------
# Working on groups of datasets - searching, fetching values
# ----------------------------------------------------------

class DatasetSearch(Resource):
    def get(self):
        """Return matching dataset + sample info based on query.
        Note that in the current implementation, if none of the parameters have been specified or other parameters
        not recognised here have been specified, this will fetch data for all datasets.

        Note that multiple parameters work with "AND" operator, so that platform_type=Microarray&projects=dc_atlas will
        only fetch datasets with both of these conditions satisfied - same applies for specific dataset ids supplied
        through the dataset_id parameter. The only exception to this is when query_string and include_samples_query are
        both specified. In this case, query_string is run against both datasets and samples collections and the union
        of the results is returned rather than an intersection.

        """
        parser = reqparse.RequestParser()
        # commonly used query keys
        parser.add_argument('dataset_id', type=str, required=False) # comma separated list of dataset ids
        parser.add_argument('query_string', type=str, required=False) # arbitrary query string - will perform text search

        # specific fields found in datasets and samples
        parser.add_argument('platform_type', type=str, required=False)
        parser.add_argument('projects', type=str, required=False)
        parser.add_argument('name', type=str, required=False)  # fetch dataset by name
        parser.add_argument('organism', type=str, required=False, default="homo sapiens")  # use 'all' to fetch all organisms

        # output control
        parser.add_argument('format', type=str, required=False)  # ['sunburst1','sunburst2']
        parser.add_argument('limit', type=int, required=False)  # limit the total number of search results

        # filtering
        parser.add_argument('include_samples_query', type=str, required=False, default="False")
        parser.add_argument('filter_Project', type=str, required=False)
        parser.add_argument('filter_platform_type', type=str, required=False)
        parser.add_argument('filter_cell_type', type=str, required=False)
        parser.add_argument('filter_tissue_of_origin', type=str, required=False)

        # sorting
        parser.add_argument('sort_field', type=str, required=False, default="name")
        parser.add_argument('sort_ascending', type=str, required=False, default="True")

        # pagination
        parser.add_argument('pagination_limit', type=int, required=False)  # leave as None if not requiring pagination
        parser.add_argument('pagination_start', type=int, required=False, default=0)  # start page
        args = parser.parse_args()

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        df = datasets.datasetMetadataFromQuery(dataset_id=args.get("dataset_id").split(',') if args.get('dataset_id') is not None else [],
                                               name=args.get("name"),
                                               query_string=args.get("query_string"),
                                               platform_type=args.get("platform_type"),
                                               projects=args.get("projects"),
                                               organism=args.get("organism"),
                                               limit=args.get("limit"),
                                               public_only=publicOnly,
                                               include_samples_query=args.get("include_samples_query").lower().startswith('t'))
        samples = datasets.samplesFromDatasetIds(df.index.tolist())

        if args.get('format')=='sunburst1':  # Returns sunburst plot data, rather than dataset + samples data
            df = datasets.sunburstData(samples)
            return df.reset_index().to_dict(orient='list')
        elif args.get('format')=='sunburst2': 
            df = datasets.sunburstData(samples, parentKey='tissue_of_origin', childKey='cell_type')
            return df.reset_index().to_dict(orient='list')
        
        if len(df)==0: # no matching result found
            return []
            
        # Add sample related columns
        samples = samples.fillna('[unassigned]')
        df['samples'] = [len(samples[samples['dataset_id']==index]) for index in df.index]  # number of samples in each dataset
        df['cell_type'] = [','.join(samples[samples['dataset_id']==index]['cell_type'].unique().tolist()) for index in df.index]
        df['tissue_of_origin'] = [','.join(samples[samples['dataset_id']==index]['tissue_of_origin'].unique().tolist()) for index in df.index]

        # Add some derived columns for convenience
        displayNames, pubmedIds, years = [], [], []
        for name in df["name"]:
            items = name.split("_")
            pubmedIds.append(items[2])
            try:
                years.append(int(items[1]))
                displayNames.append("{} ({})".format(items[0],items[1]))
            except ValueError:  # assume wrong year string to convert to int
                years.append(0)
                displayNames.append("{}".format(items[0]))
        df["display_name"] = displayNames
        df["pubmed_id"] = pubmedIds
        df["year"] = years
        
        # Sort
        df = df.sort_values(args.get('sort_field'), ascending=args.get('sort_ascending').lower().startswith('t'))
        df = df.fillna("[unassigned]")

        # Pagination
        if args.get('pagination_limit') is not None: # output is in different format
            # also need to provide numbers in each category
            counts = {}
            df['Project'] = [','.join(sorted(item)) if len(item)>0 else '[unassigned]' for item in df['projects']]
            counts['Project'] = pandas.Series(','.join(df['Project'].tolist()).split(',')).value_counts().to_dict()
            counts['Platform Type'] = df['platform_type'].value_counts().to_dict()
            counts['cell_type'] = pandas.Series(','.join(df['cell_type'].tolist()).split(',')).value_counts().to_dict()
            counts['tissue_of_origin'] = pandas.Series(','.join(df['tissue_of_origin'].tolist()).split(',')).value_counts().to_dict()
            total = len(df)

            # filter_xxx params are applied here, so that counts can still count values before these filters are applied
            if args.get('filter_Project'):
                df = df[[any(value in item for value in args.get('filter_Project').split(',')) for item in df['Project']]]
            if args.get('filter_platform_type'):
                df = df[df['platform_type']==args.get('filter_platform_type')]
            if args.get('filter_cell_type'):
                df = df[[any(value in item for value in args.get('filter_cell_type').split(',')) for item in df['cell_type']]]
            if args.get('filter_tissue_of_origin'):
                df = df[[any(value in item for value in args.get('filter_tissue_of_origin').split(',')) for item in df['tissue_of_origin']]]

            start = args.get('pagination_start')
            limit = args.get('pagination_limit')
            if len(df)>start*limit:
                subset = df.iloc[start*limit:start*limit+limit]
                return {'total':total, 'filtered_total':len(df), 'counts':counts, 'results':subset.fillna("").reset_index().to_dict(orient="records")}
            else:  # wrong combination of start and limit specified, so just return everything
                return {'total':total, 'filtered_total':len(df), 'counts':counts, 'results':df.fillna("").reset_index().to_dict(orient="records")}
        else:
            return df.reset_index().to_dict(orient="records")

class SampleSearch(Resource):
    def get(self):
        """Return matching sample info based on query.
        Note that in the current implementation, if none of the parameters have been specified or other parameters
        not recognised here have been specified, this will fetch data for all samples but a limit of 50 is imposed.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('dataset_id', type=str, required=False) # comma separated list of dataset ids
        parser.add_argument('query_string', type=str, required=False)
        parser.add_argument('field', type=str, required=False, default='')  # comma separated list of fields to include
        parser.add_argument('limit', type=int, required=False, default=50)
        parser.add_argument('orient', type=str, required=False, default='records')
        parser.add_argument('organism', type=str, required=False, default="homo sapiens")  # use 'all' to fetch all organisms
        args = parser.parse_args()

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        df = datasets.datasetMetadataFromQuery(dataset_id=args.get("dataset_id").split(',') if args.get('dataset_id') is not None else [],
                                               query_string=args.get("query_string"),
                                               limit=args.get("limit"),
                                               organism=args.get("organism"),                                               
                                               public_only=publicOnly)
        samples = datasets.samplesFromDatasetIds(df.index.tolist())

        # subset columns of samples if specified
        commonCols = set(args.get('field').split(',')).intersection(set(samples.columns))
        if commonCols:
            samples = samples[[item for item in samples.columns if item in commonCols]]

        if args.get('orient')=='records':  # include index
            samples = samples.reset_index()

        return samples.fillna("").to_dict(orient=args.get('orient'))

class Values(Resource):
    def get(self, collection, key):
        """Return all values for a key (=field) in the datasets collection.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('include_count', type=str, required=False, default="false")
        parser.add_argument('organism', type=str, required=False, default="homo sapiens")  # use 'all' to fetch all organisms
        args = parser.parse_args()

        if collection not in ['datasets','samples']:
            raise KeyNotFoundError

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        if args.get('include_count').lower().startswith('t'):
            values = datasets.allValues(collection, key, includeCount=True, public_only=publicOnly, 
                                        excludeDatasets=_exclude_some_datasets, organism=args.get('organism'))
            if values is None:
                raise KeyNotFoundError
            else:
                return values.to_dict() 
        else:
            values = datasets.allValues(collection, key, public_only=publicOnly, 
                                        excludeDatasets=_exclude_some_datasets,organism=args.get('organism'))
            if values is None:
                raise KeyNotFoundError
            else:
                return sorted(values)
        
class Download(Resource):
    def get(self):
        """Return all dataset files as one zip file.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('dataset_id', type=str, required=True) # comma separated list of dataset ids
        parser.add_argument('exclude_some_datasets', type=str, required=False, default="True")
        parser.add_argument('sep', type=str, required=False, default="\t")  # separator to use for each file
        args = parser.parse_args()

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        datasetIds = map(int, args.get('dataset_id').split(','))
        if not args.get('exclude_some_datasets').lower().startswith('f'):
            datasetIds = [id for id in datasetIds if id not in datasets._exclude_list]
        filepath = datasets.dataAsZipfile(datasetIds, publicOnly=publicOnly)
        filename = "Stemformatics_downloaded_datasets.zip"
        return send_from_directory(os.path.dirname(filepath), os.path.basename(filepath), as_attachment=True, attachment_filename=filename)
