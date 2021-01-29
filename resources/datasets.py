from flask_restful import reqparse, Resource

from models import datasets
from bson.json_util import loads, dumps

class DatasetMetadata(Resource):
    def get(self, datasetId):
        """Return dataset metadata for a dataset with datasetId.
        """        
        ds = datasets.datasetFromDatasetId(datasetId)
        if ds and ds.isPrivate():
            return None
        return ds.metadata()

class Samples(Resource):
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

class Metadata(Resource):
    def get(self, item):
        """Return metadata such as platform_types.
        """
        if item=="platform_types":
            return datasets.Dataset.platform_types
        return {}
