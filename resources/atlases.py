from flask_restful import reqparse, Resource
from flask import send_from_directory

from models import atlases, datasets
from resources import auth, errors

import os, werkzeug, pandas, json

class Atlas(Resource):
    def get(self, atlasType, item):
        """Returns different objects from atlasType depending on item.
        """
        # Note that these arguments will not be case sensitive, coming from URL
        parser = reqparse.RequestParser()
        parser.add_argument('version', type=str, required=False)
        parser.add_argument('orient', type=str, required=False, default="records")
        parser.add_argument('filtered', type=str, required=False, default='false')
        parser.add_argument('query_string', type=str, required=False, default="")
        parser.add_argument('gene_id', type=str, required=False)  # comma separated list of strings
        parser.add_argument('as_file', type=bool, required=False, default=False)
        args = parser.parse_args()

        atlas = atlases.Atlas(atlasType, version=args.get('version'))
        filtered = args.get('filtered').lower().startswith('t')
        if item=="coordinates":
            df = atlas.pcaCoordinates()
            filepath = os.path.join(atlas.atlasFilePath, "coordinates.tsv")

        elif item=="samples":
            df = atlas.sampleMatrix().fillna('')
            filepath = os.path.join(atlas.atlasFilePath, "samples.tsv")

        elif item=="expression-values":  # subset expression matrix on gene ids specified
            geneIds = args.get('gene_id').split(',') if args.get('gene_id') is not None else []
            df = atlas.expressionMatrix(filtered=filtered).loc[geneIds]
        
        elif item=="expression-file":  # this is served as a file download regardless of as_file flag
            filepath = atlas.expressionFilePath(filtered=filtered)
            args['as_file'] = True

        elif item=="genes":  # also served as a file download regardless of as_file flag
            filepath = os.path.join(atlas.atlasFilePath, "genes.tsv")
            args['as_file'] = True

        elif item=="colours-and-ordering":  # note this is actually a dictionary, not data frame, so just return it here if not serving a file
            colours = atlas.coloursAndOrdering()
            if not args.get('as_file'):
                return colours
            filepath = os.path.join(atlas.atlasFilePath, "colours.json")
        
        elif item=="possible-genes":  # no option to return a file here - just return a dictionary
            df = atlas.geneInfo().fillna("")
            df = df[df["symbol"].str.lower().str.startswith(args.get("query_string").lower())]
            return df.sort_values(["inclusion","symbol"], ascending=[False,True]).reset_index().to_dict(orient="records")
        
        else:
            return []

        if args.get('as_file'):
            return send_from_directory(os.path.dirname(filepath), os.path.basename(filepath), as_attachment=True, 
                attachment_filename="stemformatics_atlas_%s.%s.%s" % (atlasType, atlas.version, os.path.basename(filepath)))
        else:
            if args.get('orient')=='records':
                df = df.reset_index()
            return df.to_dict(orient=args.get("orient"))

class AtlasTypes(Resource):
    def get(self):
        """Return a dictionary of available atlas types and versions.
        """
        return atlases.atlasTypes()

class AtlasProjection(Resource):
    def post(self, atlasType, dataSource):
        """Project data onto the atlas of atlasType. dataSource is one of  ['Stemformatics','User'].
        """
        try:
            parser = reqparse.RequestParser()
            parser.add_argument('name', type=str, required=False)  # selected Stemformatics dataset name (eg. Helft_2017_28723558)
            parser.add_argument('test_name', type=str, required=False, default="test-data")  # user projected dataset name
            parser.add_argument('test_sample_column', type=str, required=False)  # user projected sample column to use for mapping
            # following returns a FileStorage object in this key
            parser.add_argument('test_expression', type=werkzeug.datastructures.FileStorage, location='files')
            parser.add_argument('test_samples', type=werkzeug.datastructures.FileStorage, location='files')
            args = parser.parse_args()

            if dataSource=="Stemformatics":
                name = args.get('name').split("_")[0] # only take the author name for this
                column = 'cell_type'
                # Find dataset with same name from Stemformatics
                publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
                dsIds = datasets.datasetMetadataFromQuery(name=args.get('name'), public_only=publicOnly, ids_only=True)
                ds = datasets.Dataset(dsIds[0])  # should be one matching dataset
                df = ds.expressionMatrix(key="genes" if ds.metadata()['platform_type']=='Microarray' else "raw")
                samples = ds.samples()

            else:
                name = args.get('test_name')
                df = pandas.read_csv(args.get('test_expression'), sep='\t', index_col=0)
                samples = pandas.read_csv(args.get('test_samples'), sep='\t', index_col=0)
                # Some validation on user supplied data
                if len(df)==0:
                    return {'error': 'The expression matrix came back as zero length. Check its format.'}
                elif len(samples)==0:
                    return {'error': 'The sample table came back as zero length. Check its format and ensure its row index match columns of expression matrix.'}
                samples = samples.loc[df.columns]
                column = args.get('test_sample_column')
                if column not in samples.columns: column = samples.columns[0]
                
            # Create atlas data instance
            atlas = atlases.Atlas(atlasType)

            # Perform projection
            result = atlas.projection(name, df, includeCombinedCoords=False)
            if result["error"] !="": # Returning empty data frame may cause exception when trying to parse as json, so just return error string
                return {"error": result["error"]}

            # Prepare the dictionary to return - each object must be JSON serializable (so don't return data frame).
            result["coords"] = result["coords"].to_dict(orient="records")
            result["samples"] = samples.reset_index().fillna('').to_dict(orient="records")
            result["sampleIds"] = ["%s_%s" % (name, item) for item in samples.index]
            result["column"] = column
            if "combinedCoords" in result:
                result["combinedCoords"] = result["combinedCoords"].to_dict(orient="split")
                
            return result
            
        except:
            raise errors.DatasetProjectionFailedError
