"""
Examples of how to run this script (ensure you're in the application directory):
(Requires environment variables EXPRESSION_FILEPATH and ATLAS_FILEPATH, which point to where the expression and atlas files are)
(s4m-api) [ec2-user@api-dev s4m-api]$ python -m scripts.check_data_consistency -a

"""

import os, sys, pandas, re, argparse

sys.path.append(os.path.join(sys.path[0]))
from models import utilities, datasets, atlases

def checkMissingData(publicDatasetsOnly=False, atlasDatasetsOnly=False, checkSampleIds=False):
    """Function to check for inconsistencies, where dataset metadata may have datasets with no corresponding
    sample table, for example.
    """
    if "EXPRESSION_FILEPATH" not in os.environ:
        print("No EXPRESSION_FILEPATH in environment")
        return
    if "ATLAS_FILEPATH" not in os.environ:
        print("No ATLAS_FILEPATH in environment")
        return

    database = utilities.mongoClient()["dataportal"]
    option = {'private':False} if publicDatasetsOnly else {}
    if atlasDatasetsOnly:
        option["projects"] = {"$in":["%s_atlas" % atlasType for atlasType in atlases.Atlas.all_atlas_types]}

    datasetIdsFromMetadata = set([item["dataset_id"] for item in database["datasets"].find(option)])
    datasetIdsFromSamples = set([item["dataset_id"] for item in database["samples"].find({})])
    datasetIdsFromExpression = set()
    for item in os.listdir(os.environ.get("EXPRESSION_FILEPATH")):
        match = re.findall("^\d{4}", item)
        if len(match)>0: datasetIdsFromExpression.add(int(match[0]))
    datasetIdsFromAtlasFiles = set()
    for atlasType in atlases.Atlas.all_atlas_types:
        datasetIdsFromAtlasFiles =  datasetIdsFromAtlasFiles.union(set(atlases.Atlas(atlasType).datasetIds()))

    print("Number of dataset ids from metadata", len(datasetIdsFromMetadata))
    print("Number of dataset ids from samples", len(datasetIdsFromSamples))
    print("Number of dataset ids from expression", len(datasetIdsFromExpression))
    print("Number of dataset ids from atlas files", len(datasetIdsFromAtlasFiles))

    diff1 = datasetIdsFromMetadata.difference(datasetIdsFromSamples)
    diff2 = datasetIdsFromMetadata.difference(datasetIdsFromExpression)
    diff3 = datasetIdsFromAtlasFiles.difference(datasetIdsFromMetadata)
    print("\nDataset ids from metadata missing from samples", len(diff1), sorted(diff1))
    print("Dataset ids from metadata missing from expression", len(diff2), sorted(diff2))
    print("Dataset ids from atlas files missing from metadata", len(diff3), sorted(diff3))

    if checkSampleIds:    # Check for consistency between sample ids - this takes longer since we need to open each expression file
        print("\n")
        for datasetId in sorted(datasetIdsFromMetadata):
            ds = datasets.Dataset(datasetId)
            sampleIdsFromMetadata = set(ds.samples().index)
            sampleIdsFromExpression = set(ds.expressionMatrix().columns)
            if sampleIdsFromMetadata!=sampleIdsFromExpression:
                print("Non-matching sample ids for dataset %s" % (datasetId))


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="only look at public datasets", action="store_true")
    parser.add_argument("-a", help="only look at atlas datasets", action="store_true")
    parser.add_argument("-s", help="check sample ids too", action="store_true")
    args = parser.parse_args()
    
    checkMissingData(publicDatasetsOnly=args.p, atlasDatasetsOnly=args.a, checkSampleIds=args.s)
