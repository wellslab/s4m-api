from flask_restful import reqparse, Resource

from resources import auth
from models import datasets
from resources.errors import DatasetIdNotFoundError, DatasetIsPrivateError, UserNotAuthenticatedError

class DatasetMetadata(Resource):
    def get(self, datasetId):
        """Return dataset metadata for a dataset with datasetId.
        """
        try:
            ds = datasets.Dataset(datasetId)
            if ds.isPrivate() and not auth.AuthUser().username():
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

        try:
            ds = datasets.Dataset(datasetId)
            if ds.isPrivate() and not auth.AuthUser().username():
                raise DatasetIsPrivateError
        except datasets.DatasetIdNotFoundError:
            raise DatasetIdNotFoundError
        
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

        try:
            ds = datasets.Dataset(datasetId)
            if ds.isPrivate() and not auth.AuthUser().username():
                raise DatasetIsPrivateError
        except datasets.DatasetIdNotFoundError:
            raise DatasetIdNotFoundError

        return ds.expressionMatrix().loc[args.get('gene_id')].to_dict(orient=args.get('orient'))

class DatasetPca(Resource):
    def get(self, datasetId):
        """Return pca data for a dataset with datasetId.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('orient', type=str, required=False, default="records")
        args = parser.parse_args()

        try:
            ds = datasets.Dataset(datasetId)
            if ds.isPrivate() and not auth.AuthUser().username():
                raise DatasetIsPrivateError
        except datasets.DatasetIdNotFoundError:
            raise DatasetIdNotFoundError

        return {'coordinates': ds.pcaCoordinates().to_dict(orient=args.get('orient')), 
                'attributes': ds.pcaAttributes().to_dict(orient=args.get('orient'))}

class DatasetSearch(Resource):
    def get(self):
        """Return matching dataset + sample info based on query.
        Note that in the current implementation, if none of the parameters have been specified or other parameters
        not recognised here have been specified, this will fetch data for all datasets.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('dataset_id', type=str, required=False, action='append') # list of dataset ids
        parser.add_argument('query_string', type=str, required=False)
        parser.add_argument('platform_type', type=str, required=False)
        parser.add_argument('projects', type=str, required=False)
        parser.add_argument('name', type=str, required=False)
        parser.add_argument('format', type=str, required=False)  # ['sunburst1','sunburst2']
        parser.add_argument('limit', type=int, required=False)
        args = parser.parse_args()

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        df = datasets.datasetMetadataFromQuery(dataset_id=args.get("dataset_id"),
                                               name=args.get("name"),
                                               query_string=args.get("query_string"),
                                               platform_type=args.get("platform_type"),
                                               projects=args.get("projects"),
                                               limit=args.get("limit"),
                                               public_only=publicOnly)
        samples = datasets.samplesFromDatasetIds(df.index.tolist())

        if args.get('format')=='sunburst1':  # Returns sunburst plot data, rather than dataset + samples data
            df = datasets.sunburstData(samples)
            return df.reset_index().to_dict(orient='list')
        elif args.get('format')=='sunburst2': 
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

class SampleSearch(Resource):
    def get(self):
        """Return matching sample info based on query.
        Note that in the current implementation, if none of the parameters have been specified or other parameters
        not recognised here have been specified, this will fetch data for all samples but a limit of 50 is imposed.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('dataset_id', type=str, required=False, action='append') # list of dataset ids
        parser.add_argument('query_string', type=str, required=False)
        parser.add_argument('field', type=str, required=False, action='append')  # list of fields to include
        parser.add_argument('limit', type=int, required=False, default=50)
        parser.add_argument('orient', type=str, required=False, default='records')
        args = parser.parse_args()

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        df = datasets.datasetMetadataFromQuery(dataset_id=args.get("dataset_id"),
                                               query_string=args.get("query_string"),
                                               limit=args.get("limit"),
                                               public_only=publicOnly)
        samples = datasets.samplesFromDatasetIds(df.index.tolist())

        # subset columns of samples if specified
        fieldsToReturn = [item for item in args.get('field') if item in samples.columns]
        if len(fieldsToReturn)>0:
            samples = samples[fieldsToReturn]

        if args.get('orient')=='records':  # include index
            samples = sample.reset_index()

        return samples.fillna("").to_dict(orient=args.get('orient'))


class ValuesDatasets(Resource):
    def get(self, key):
        """Return all values for a key (=field) in the datasets collection.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('include_count', type=bool, required=False, default=False)
        args = parser.parse_args()

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        if args.get('include_count'):
            return datasets.allValues("datasets", key, includeCount=True, public_only=publicOnly).to_dict()
        else:
            return sorted(datasets.allValues("datasets", key, public_only=publicOnly))

class ValuesSamples(Resource):
    def get(self, key):
        """Return all values for a key (=field) in the samples collection.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('include_count', type=bool, required=False, default=False)
        args = parser.parse_args()

        publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
        if args.get('include_count'):
            return datasets.allValues("samples", key, includeCount=True, public_only=publicOnly).to_dict()
        else:
            return sorted(datasets.allValues("samples", key, public_only=publicOnly))

