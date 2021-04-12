from flask_restful import reqparse, Resource

from models import genes

class Geneset(Resource):
    def get(self):
        """Return attibutes of genes from query parameters. The returned dictionary looks like
            [
                {
                    "gene_id": "ENSG00000133055",
                    "gene_name": "MYBPH",
                    "gene_description": "myosin binding protein H [Source:HGNC Symbol;Acc:HGNC:7552]"
                },
                ...
            ]
        """
        parser = reqparse.RequestParser()
        parser.add_argument('gene_id', type=str, required=False, default="", action="append")  # list of strings
        parser.add_argument('search_string', type=str, required=False, default="")
        parser.add_argument('limit', type=int, required=False, default=50, action="append")
        parser.add_argument('orient', type=str, required=False, default="records")
        args = parser.parse_args()
        
        df = genes.geneset(geneIds=args.get("gene_id"), searchString=args.get("search_string"), limit=args.get("limit"))
        if args.get('orient')=='records':
            df = df.reset_index()
        return df.fillna("").to_dict(orient=args.get("orient"))

