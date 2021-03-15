"""
Examples of how to run this script (ensure you're in the application directory):
(Requires environment variable EXPRESSION_FILEPATH, which points to where the expression files are)
(s4m-api) [ec2-user@api-dev s4m-api]$ python -m scripts.check_data_consistency checkMissingData

"""

import os, sys, pandas, re, argparse

sys.path.append(os.path.join(sys.path[0]))
from models import utilities, datasets, atlases

def checkMissingData(publicDatasetsOnly=False, atlasDatasetsOnly=False):
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

    return
    diff1 = datasetIdsFromSamples.difference(datasetIdsFromMetadata)
    diff2 = datasetIdsFromSamples.difference(datasetIdsFromExpression)
    print("\nDataset ids from samples missing from metadata", len(diff1), list(diff1)[:3])
    print("Dataset ids from samples missing from expression", len(diff2), list(diff2)[:3])

    diff1 = datasetIdsFromExpression.difference(datasetIdsFromMetadata)
    diff2 = datasetIdsFromExpression.difference(datasetIdsFromSamples)
    print("\nDataset ids from expression missing from metadata", len(diff1), list(diff1)[:3])
    print("Dataset ids from expression missing from samples", len(diff2), list(diff2)[:3])

def checkMismatchedSampleIds():
    """Some datasets have different sample ids in the expression file compared to samples table.
    """
    return

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="only look at public datasets", action="store_true")
    parser.add_argument("-a", help="only look at atlas datasets", action="store_true")
    args = parser.parse_args()
    
    checkMissingData(publicDatasetsOnly=args.p, atlasDatasetsOnly=args.a)
