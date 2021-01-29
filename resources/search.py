from flask_restful import reqparse, Resource

import pandas
from models import datasets

class Search(Resource):
    def get(self):
        
        parser = reqparse.RequestParser()
        parser.add_argument('platform_type', type=str, required=False)
        args = parser.parse_args()

        datasetIds = datasets.datasetIdsFromQuery(platform_type=args.get("platform_type"))
        df = pandas.DataFrame.from_records([datasets.datasetFromDatasetId(dsId).metadata() for dsId in datasetIds])
        
        # Add some derived columns for convenience
        displayNames, pubmedIds = [], []
        for name in df["name"]:
            items = name.split("_")
            displayNames.append("{} ({})".format(items[0],items[1]) if len(items)>1 else "")
            pubmedIds.append(items[2] if len(items)>2 else "")
        df["display_name"] = displayNames
        df["pubmed_id"] = pubmedIds
        
        return df.to_dict(orient="records")