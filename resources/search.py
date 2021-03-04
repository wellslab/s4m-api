from flask_restful import reqparse, Resource

import pandas
from models import datasets

class SearchDatasets(Resource):
    def get(self):
        
        parser = reqparse.RequestParser()
        parser.add_argument('platform_type', type=str, required=False)
        parser.add_argument('limit', type=int, required=False)
        args = parser.parse_args()

        datasetIds = datasets.datasetIdsFromQuery(platform_type=args.get("platform_type"), limit=args.get("limit"))
        df = pandas.DataFrame.from_records([datasets.Dataset(dsId).metadata() for dsId in datasetIds])
        
        # Add some derived columns for convenience
        displayNames, pubmedIds = [], []
        for name in df["name"]:
            items = name.split("_")
            displayNames.append("{} ({})".format(items[0],items[1]) if len(items)>1 else "")
            pubmedIds.append(items[2] if len(items)>2 else "")
        df["display_name"] = displayNames
        df["pubmed_id"] = pubmedIds
        
        res = df.fillna("").to_dict(orient="records")
        return res

class SearchSamples(Resource):
    def get(self):
        
        parser = reqparse.RequestParser()
        parser.add_argument('cell_type', type=str, required=False)
        parser.add_argument('format', type=str, required=False)
        args = parser.parse_args()

        if args.get('format')=='sunburst':
            import pandas, os
            df = pandas.read_csv(os.path.join(os.environ["ATLAS_FILEPATH"],'sunburst1.tsv'), sep='\t').fillna("")
            return df.to_dict(orient='list')

        df = datasets.samplesFromQuery(cell_type=args.get("cell_type"))
        df = df.fillna("")        
        
        return df.to_dict(orient="records")