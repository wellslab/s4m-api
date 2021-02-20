from flask_restful import reqparse, Resource

from models import datasets
#from bson.json_util import loads, dumps

class DatasetMetadata(Resource):
    def get(self, datasetId):
        """Return dataset metadata for a dataset with datasetId.
        """        
        ds = datasets.datasetFromDatasetId(datasetId)
        if ds and ds.isPrivate():
            return None
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
        args = parser.parse_args()

        ds = datasets.datasetFromDatasetId(datasetId)
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
        """Return expression table for a dataset with datasetId.
        """
        return {}

class DatasetGovernance(Resource):
    def get(self, datasetId):
        """Return governmence data for a dataset with datasetId.
        """
        return {}

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
