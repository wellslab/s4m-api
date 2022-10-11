from flask_restful import reqparse, Resource
from flask import send_from_directory

from models import atlases, datasets
from resources import auth, errors
from datetime import datetime
import pandas as pd

import os, werkzeug, pandas, json, uuid, pathlib, csv
# import pickle

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

class AtlasProjection(Resource):
    def post(self, atlasType, dataSource):
        """Project data onto the atlas of atlasType. dataSource is one of  ['stemformatics','user','user-single'].
        """

        try:
            # Project single-cell-data - user upload
            if dataSource.lower()=="user-single":
                parser = reqparse.RequestParser()
                parser.add_argument('email', type=str, required=False)
                parser.add_argument('data', type=werkzeug.datastructures.FileStorage, location='files')
                args = parser.parse_args()

                now = datetime.now()

                email = args.get('email')
                df = pandas.read_csv(args.get('data'), sep='\t', index_col=0)
                if len(df)==0:
                    print('error')
                    return {'error': 'The data matrix came back as zero length. Check its format.'}
                # print(df)

                # Read samples.tsv saved in folder

                # Create unique ID
                id = uuid.uuid4()

                # Check for unique ID in ids.tsv 
                try:
                    path = "/mnt/stemformatics-data/user_projection_data/ids.csv"

                    ids_df = pandas.read_csv(path, index_col=0)

                    unique = False
                    # Loop until unique ID
                    while not unique:
                        if id in ids_df.index:
                            id = uuid.uuid4()
                        else:
                            unique = True
                    
                    # Create file called ID
                    new_path = "/mnt/stemformatics-data/user_projection_data/{}".format(id)
                    if not os.path.exists(new_path):
                        os.makedirs(new_path)
                        
                    # Save df to input
                    df.to_csv(new_path + '/input.tsv', sep='\t')

                    # Append details to ids.csv
                    now = datetime.now()
                    try:
                        data = [id, email, atlasType, now.strftime('%d/%m/%Y %H:%M:%S'),'','']
                        with open(path, 'a', newline='') as fd:
                            writer_object = csv.writer(fd)
                            writer_object.writerow(data)
                            fd.close()
                    except FileNotFoundError:
                        print('File not found.')
                    except IOError:
                        print('File IO error.')
                    except Exception as e:
                        print(e)
                    
                    return
                    

                except IOError:
                    print("Error: file does not exist.")
                    # print(pathlib.Path.cwd())
                    return

            else:
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
                # print(atlas.expressionMatrix(filtered=True))
                # print(atlas.geneInfo())

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

                # print(type(result))
                # output_path = "/mnt/stemformatics-data/user_projection_data/test.json"
                # with open(output_path, 'w') as f:
                #     json.dump(result, f, ensure_ascii=False, indent=4)

                # print(result)
                # print(type(result["coords"]))
                return result
            
        except:
            raise errors.DatasetProjectionFailedError

# Get - encaapsulates all parameters in URL, Post - from form submitted to server
# Return JSON object to UI server
class AtlasProjectionResults(Resource):

    # def get(userID):
    #     parser = reqparse.RequestParser()
    #     parser.add_argument('userID', type=str, required=False)
    #     # parser.add_argument('atlasType', type=str, required=False)
    #     args = parser.parse_args()
    #     id = args.get("userID")
    #     # atlasTypee = args.get("atlasType")
    #     # print(atlasTypee)
    #     atlasType = "myeloid"

    #     # Return JSON object from userID folder
    #     input_path = "/mnt/stemformatics-data/user_projection_data/{}/input.tsv".format(id)
    #     df = pandas.read_csv(input_path, sep='\t', index_col=0)

    #     # Create atlas data instance
    #     atlas = atlases.Atlas(atlasType)
    #     # print(atlasType)
    #     # print(atlas.expressionMatrix(filtered=True))
    #     # print(atlas.geneInfo())

    #     # Perform projection
    #     result = atlas.projection(id, df, includeCombinedCoords=False)
    #     if result["error"] !="": # Returning empty data frame may cause exception when trying to parse as json, so just return error string
    #         return {"error": result["error"]}

    #     # Read samples.tsv from folder
    #     samples_path = "/mnt/stemformatics-data/user_projection_data/ae962dbd-aff5-484a-9cdc-94885aeafe28/input.tsv"
    #     samples = pandas.read_csv(samples_path, sep='\t', index_col=0)
    #     # print(samples)

    #     if len(df)==0:
    #         return {'error': 'The expression matrix came back as zero length. Check its format.'}
    #     elif len(samples)==0:
    #         return {'error': 'The sample table came back as zero length. Check its format and ensure its row index match columns of expression matrix.'}

    #     samples = samples.loc[df.columns]
    #     # column = args.get('test_sample_column')
    #     column = 'Sample Type'
    #     if column not in samples.columns: column = samples.columns[0]

    #     # Prepare the dictionary to return - each object must be JSON serializable (so don't return data frame).
    #     name = "notta"
    #     result["coords"] = result["coords"].to_dict(orient="records")
    #     result["samples"] = samples.reset_index().fillna('').to_dict(orient="records")
    #     result["sampleIds"] = ["%s_%s" % (name, item) for item in samples.index]
    #     result["column"] = column
    #     if "combinedCoords" in result:
    #         result["combinedCoords"] = result["combinedCoords"].to_dict(orient="split")
    #     if "capybara" in result:
    #         for col in result["capybara"]:
    #             result["capybara"][col] = result["capybara"][col].to_dict(orient="split")
        
    #     return result

    def jsonKeys2int(x):
        if isinstance(x, dict):
            return {int(k):v for k,v in x.items()}
        return x
    
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument("id", type=str, required=False)
        args = parser.parse_args()
        id = args.get("id")
        # print(id)

        # Return JSON object from userID folder
        path = "/../mnt/stemformatics-data/user_projection_data/{}/output.json".format(id)
        with open(path, 'r') as f:
            # result = json.load(f, object_hook=jsonKeys2int)
            result = json.load(f)


        # print("Type after loading...")
        # print(type(result))

        return result

    # def get(userID):
    #     # print(userID)
    #     parser = reqparse.RequestParser()
    #     parser.add_argument('userID', type=str, required=False)
    #     args = parser.parse_args()
    #     id = args.get("userID")
    #     print(id)

    #     # Return JSON object from userID folder
    #     path = "/../mnt/stemformatics-data/user_projection_data/{}/output.json".format(id)
    #     with open(path, 'r') as f:
    #         # result = json.load(f, object_hook=jsonKeys2int)
    #         result = json.load(f)


    #     print("Type after loading...")
    #     print(type(result))
    #     # print(type(result["coords"]))
    #     # return userID

    #     # input_path = "/mnt/stemformatics-data/user_projection_data/{}/input.tsv".format(id)
    #     # input_df = pandas.read_csv(input_path, sep='\t', index_col=0)

    #     # samples_path = "/mnt/stemformatics-data/user_projection_data/ae962dbd-aff5-484a-9cdc-94885aeafe28/input.tsv"
    #     # samples = pandas.read_csv(samples_path, sep='\t', index_col=0)
    #     # samples = samples.loc[input_df.columns]
    #     # # column = args.get('test_sample_column')
    #     # column = 'Sample Type'
    #     # if column not in samples.columns: column = samples.columns[0]


    #     # # Prepare the dictionary to return - each object must be JSON serializable (so don't return data frame).
    #     # name = id
    #     # column = 'Sample Type'
    #     # result["coords"] = result["coords"].to_dict(orient="records")
    #     # result["samples"] = samples.reset_index().fillna('').to_dict(orient="records")
    #     # result["sampleIds"] = ["%s_%s" % (name, item) for item in samples.index]
    #     # result["column"] = column
    #     # if "combinedCoords" in result:
    #     #     result["combinedCoords"] = result["combinedCoords"].to_dict(orient="split")
    #     # if "capybara" in result:
    #     #     for col in result["capybara"]:
    #     #         result["capybara"][col] = result["capybara"][col].to_dict(orient="split")

    #     # print(result)
    #     return result
