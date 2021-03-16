from flask_restful import reqparse, Resource

from models import datasets
from resources.errors import DatasetIdNotFoundError, DatasetIsPrivateError

class DatasetMetadata(Resource):
    def get(self, datasetId):
        """Return dataset metadata for a dataset with datasetId.
        """
        try:
            ds = datasets.Dataset(datasetId)
            if ds and ds.isPrivate():
                raise DatasetIsPrivateError
            return ds.metadata()
        except datasets.DatasetIdNotFoundError:
            raise DatasetIdNotFoundError

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
        args = parser.parse_args()

        ds = datasets.Dataset(datasetId)
        if ds and ds.isPrivate():
            return DatasetIsPrivateError
        
        # These fields are used internally by the model and not useful for API
        hideKeys = ["dataset_id"]
        df = ds.samples()
        if len(df)>0:
            df = df.drop(hideKeys, axis=1).fillna(args.get("na"))
            return df.reset_index().to_dict(orient=args.get("orient"))
        else:
            return {}

class DatasetExpression(Resource):
    def get(self, datasetId):
        """Return expression table for a dataset with datasetId and gene id(s).
        """
        parser = reqparse.RequestParser()
        parser.add_argument('gene_id', type=str, required=False, default="", action="append")
        parser.add_argument('key', type=str, required=False, default="raw")
        parser.add_argument('orient', type=str, required=False, default="records")
        args = parser.parse_args()

        ds = datasets.Dataset(datasetId)
        if ds and ds.isPrivate():
            return DatasetIsPrivateError
        return ds.expressionMatrix().loc[args.get('gene_id')].to_dict(orient=args.get('orient'))

class DatasetGovernance(Resource):
    def get(self, datasetId):
        """Return governmence data for a dataset with datasetId.
        """
        return {}

class DatasetPca(Resource):
    def get(self, datasetId):
        """Return pca data for a dataset with datasetId.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('orient', type=str, required=False, default="records")
        args = parser.parse_args()

        ds = datasets.Dataset(datasetId)
        if ds and ds.isPrivate():
            return DatasetIsPrivateError
        return {'coordinates': ds.pcaCoordinates().to_dict(orient=args.get('orient')), 
                'attributes': ds.pcaAttributes().to_dict(orient=args.get('orient'))}

class DatasetSearch(Resource):
    def get(self):
        """Return matching dataset + sample info based on query.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('dataset_id', type=str, required=False, action='append') # list of dataset ids
        parser.add_argument('query_string', type=str, required=False)
        parser.add_argument('platform_type', type=str, required=False)
        parser.add_argument('projects', type=str, required=False)
        parser.add_argument('format', type=str, required=False)
        parser.add_argument('limit', type=int, required=False)
        args = parser.parse_args()

        df = datasets.datasetMetadataFromQuery(dataset_id=args.get("dataset_id"),
                                               query_string=args.get("query_string"),
                                               platform_type=args.get("platform_type"),
                                               projects=args.get("projects"),
                                               limit=args.get("limit"))
        samples = datasets.samplesFromDatasetIds(df.index.tolist())

        if args.get('format')=='sunburst1':  # Returns sunburst plot data, rather than dataset + samples data
            df = datasets.sunburstData(samples)
            return df.reset_index().to_dict(orient='list')
        elif args.get('format')=='sunburst2':  # Returns sunburst plot data, rather than dataset + samples data
            df = datasets.sunburstData(samples, parentKey='tissue_of_origin', childKey='cell_type')
            return df.reset_index().to_dict(orient='list')
        
        # Add sample related columns
        samples = samples.fillna('')
        df['samples'] = [len(samples[samples['dataset_id']==index]) for index in df.index]
        df['cell_type'] = [','.join(samples[samples['dataset_id']==index]['cell_type'].unique().tolist()) for index in df.index]

        # Add some derived columns for convenience
        displayNames, pubmedIds = [], []
        for name in df["name"]:
            items = name.split("_")
            displayNames.append("{} ({})".format(items[0],items[1]))
            pubmedIds.append(items[2])
        df["display_name"] = displayNames
        df["pubmed_id"] = pubmedIds
        
        return df.fillna("").reset_index().to_dict(orient="records")


class ValuesDatasets(Resource):
    def get(self, key):
        """Return all values for a key (=field) in the datasets collection.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('include_count', type=bool, required=False, default=False)
        args = parser.parse_args()

        if args.get('include_count'):
            return datasets.allValues("datasets", key, includeCount=True).to_dict()
        else:
            return sorted(datasets.allValues("datasets", key))

class ValuesSamples(Resource):
    def get(self, key):
        """Return all values for a key (=field) in the samples collection.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('include_count', type=bool, required=False, default=False)
        args = parser.parse_args()

        if args.get('include_count'):
            return datasets.allValues("samples", key, includeCount=True).to_dict()
        else:
            return sorted(datasets.allValues("samples", key))
