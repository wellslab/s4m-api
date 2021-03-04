from flask_restful import reqparse, Resource

from models import datasets
from resources.errors import DatasetIdNotFoundError

class DatasetMetadata(Resource):
    def get(self, datasetId):
        """Return dataset metadata for a dataset with datasetId.
        """
        try:
            ds = datasets.Dataset(datasetId)
            if ds and ds.isPrivate():
                return None
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
            return None
        
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
            return None
        return {'coordinates': ds.pcaCoordinates().to_dict(orient=args.get('orient')), 
                'attributes': ds.pcaAttributes().to_dict(orient=args.get('orient'))}

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
