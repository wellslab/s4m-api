from flask_restful import reqparse, Resource
from flask import send_from_directory

from models import atlases, datasets
from resources import auth, errors

import os, werkzeug, pandas, json #, uuid, pathlib, csv

class AtlasTypes(Resource):
    def get(self):
        """Return a dictionary of available atlas types and versions.
        """
        return atlases.atlasTypes()

class Atlas(Resource):
    def get(self, atlasType, item):
        """Returns different objects from atlasType depending on item.
        """
        # Note that these arguments will not be case sensitive, coming from URL
        parser = reqparse.RequestParser()
        parser.add_argument('version', type=str, required=False)    # specify to get info from a specific version of the atlas
        parser.add_argument('orient', type=str, required=False, default="records")  # when returning data frames
        parser.add_argument('filtered', type=str, required=False, default='false')  # for expression values (use filtered genes)
        parser.add_argument('query_string', type=str, required=False, default="")   # for possible-genes (all genes which belong to the atlas)
        parser.add_argument('gene_id', type=str, required=False)  # for expression-values and heatmap-data: comma separated list of strings
        parser.add_argument('as_file', type=bool, required=False, default=False)  # download some data frames as file

        # More specific to heatmap-data:
        parser.add_argument('cluster_columns', type=str, required=False, default='true')
        parser.add_argument('groupby', type=str, required=False)  # eg. "treatment"
        parser.add_argument('subsetby', type=str, required=False)  # eg. "time"
        parser.add_argument('subsetby_item', type=str, required=False)  # eg. "2hr"
        parser.add_argument('relative_value', type=str, required=False, default="zscore")  # calculate row values relative to this

        args = parser.parse_args()
        filtered = args.get('filtered').lower().startswith('t')
        clusterColumns = args.get('cluster_columns').lower().startswith('t')

        atlas = atlases.Atlas(atlasType, version=args.get('version'))
        
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
        
        elif item=="heatmap-data":
            geneIds = args.get('gene_id').split(',') if args.get('gene_id') is not None else []
            result = atlas.heatmapData(geneIds=geneIds, clusterColumns=clusterColumns, groupby=args.get('groupby'), 
                                       subsetby=args.get('subsetby'), subsetbyItem=args.get('subsetby_item'), relativeValue=args.get('relative_value'))
            result['dataframe'] = result['dataframe'].to_dict(orient='split')
            return result
            
        else:
            return []

        if args.get('as_file'):
            return send_from_directory(os.path.dirname(filepath), os.path.basename(filepath), as_attachment=True, 
                attachment_filename="stemformatics_atlas_%s.%s.%s" % (atlasType, atlas.version, os.path.basename(filepath)))
        else:
            if args.get('orient')=='records':
                df = df.reset_index()
            return df.to_dict(orient=args.get("orient"))

class AtlasProjection(Resource):
    # Check the bulk data file format
    def check_format(self, expression, samples, sample_column):
        # Check number of columns in expression matrix and number of rows in sample matrix is the same
        if len(expression.columns) != len(samples):
            return {'error': 'The number of columns in the expression matrix and rows in the sample matrix need to match. Check the file formats.'}
        
        # Check the expression matrix columns are present in the row IDs in the sample matrix
        for column in expression.columns:
            if column not in samples.index:
                return {'error': 'The column headings in the expression matrix need to be present in the row IDs in the sample matrix. Check the file formats.'}
            
        # Check the expression matrix uses Ensembl IDs
        for index in expression.index:
            if index[0:4] != 'ENSG':
                return {'error': 'The expression matrix must use Ensembl IDs as row IDs. Check the file format.'}
        
        # Check the sample column given exists in the sample matrix
        if sample_column not in samples.columns and sample_column != "":
            return {'error': 'The sample column must be a column that exists in the sample matrix. Check the sample column.'}
        
        return {'error': ""}


    def post(self, atlasType, dataSource):
        """Project data onto the atlas of atlasType. dataSource is one of  ['stemformatics','user'].
        The word 'query' would be better to use rather than 'test' for the data being projected,
        but it got stuck here for historical reasons.
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

            if dataSource.lower()=="stemformatics":
                name = args.get('name').split("_")[0] # only take the author name for this
                # Find dataset with same name from Stemformatics
                publicOnly = auth.AuthUser().username()==None  # public datasets only if authenticated username returns None
                dsId = datasets.datasetIdFromName(args.get('name'), publicOnly=publicOnly)
                ds = datasets.Dataset(dsId)
                df = ds.expressionMatrix(key="genes" if ds.metadata()['platform_type']=='Microarray' else "raw")
                samples = ds.samples()
                # Select column to use - best column may not be cell_type though
                column = 'cell_type'
                for col in ['cell_type','sample_type','final_cell_type']:
                    if len(samples[col].unique())>=2:  # use this
                        column = col
                        break
                samples = samples.fillna('[not assigned]')

            else:
                name = args.get('test_name')
                df = pandas.read_csv(args.get('test_expression'), sep='\t', index_col=0)
                samples = pandas.read_csv(args.get('test_samples'), sep='\t', index_col=0)
                
                # Check validation on the uploaded data files format
                check = self.check_format(df, samples, args.get('test_sample_column'))
                if check["error"] != "":
                    return {"error": check["error"]}
                
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
            if "capybara" in result:
                for col in result["capybara"]:
                    result["capybara"][col] = result["capybara"][col].to_dict(orient="split")

            return result
        
        except:
            raise errors.DatasetProjectionFailedError

# Get - encaapsulates all parameters in URL, Post - from form submitted to server
# Return JSON object to UI server
class AtlasProjectionResults(Resource):
    
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument("id", type=str, required=False)
        args = parser.parse_args()
        id = args.get("id")

        # Return JSON object from userID folder
        path = "/../mnt/stemformatics-data/user_projection_data/{}/output.json".format(id)
        with open(path, 'r') as f:
            result = json.load(f)

        return result

