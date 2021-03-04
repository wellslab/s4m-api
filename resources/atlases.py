from flask_restful import reqparse, Resource
from flask import send_from_directory

from models import atlases
import os

class Atlas(Resource):
    def get(self, atlasType, item):
        """Returns different objects from atlasType depending on item.
        """
        # Note that these arguments will not be case sensitive, coming from URL
        parser = reqparse.RequestParser()
        parser.add_argument('orient', type=str, required=False, default="records")
        parser.add_argument('filtered', type=str, required=False, default=True)
        parser.add_argument('query_string', type=str, required=False, default="")
        parser.add_argument('gene_id', type=str, required=False, default="", action="append")  # list of strings
        args = parser.parse_args()

        atlas = atlases.Atlas(atlasType)
        if item=="coordinates":
            return atlas.pcaCoordinates().to_dict(orient=args.get("orient"))

        elif item=="samples":
            return atlas.sampleMatrix().fillna('').to_dict(orient=args.get("orient"))

        elif item=="expression-file":  # this is served as a file download
            filepath = atlas.expressionFilePath(args.get('filtered'))
            return send_from_directory(os.path.dirname(filepath), os.path.basename(filepath), as_attachment=True, 
                attachment_filename="stemformatics_atlas_%s.%s.%s" % (atlasType, atlas.version, os.path.basename(filepath)))
        
        elif item=="expression-values":
            df = atlas.expressionValues(args.get("gene_id"))
            return df.to_dict(orient=args.get("orient"))
        
        elif item=="colours-and-ordering":
            return atlas.coloursAndOrdering()
        
        elif item=="possible-genes":
            df = atlas.geneInfo().fillna("")
            df = df[df["symbol"].str.lower().str.startswith(args.get("query_string").lower())]
            return df.sort_values(["inclusion","symbol"], ascending=[False,True]).reset_index().to_dict(orient="records")
        
        else:
            return []
