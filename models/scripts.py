"""
This is a collection functions mainly used to manipulate the underlying data.

Notes: 2020-12-31.
MDAP produced a cleaned version of a subset of the sample table from stemformatics biosamples_metadata table.
This file is called biosamples_atlas_v3.1_clean.txt and was used as input to samples collection in mongodb using
_setSamplesCollectionFromFile function below.
Since this is only a subset, other sample data also need to be added to the samples collection in mongodb.
Catherine has done some cleaning of this other data, and the file called biosamples_metadata-edited.tsv represents
this data (shared by her on google drive).
Note that cleaned data may also be not complete per dataset - ie. some samples from a dataset may be in the cleaned
set, while others from the same dataset may be in the uncleaned set. This is the breakdown of datasets where this
is observed:
    dataset_id  samples_from_biosamples_metadata-edited samples_from_biosamples_atlas_v3.1_clean
    1000    50  9
    6610    24  12
    6612    19  2
    6646    1448    124
Hence for these datasets, we need to clean up the rest of the samples before they become consistent with their clean
siblings. Eg. For dataset 1000, we have 9 cleaned up samples and 41 unclean samples.

"""
import pandas, sys, os
from models import utilities, datasets

database = utilities.mongoClient()["dataportal"]

def convertBiosamplesMetadata(infile, outfile):
    """Given a dump of the old stemformatics biosamples_metadata table, in the form of a tab separated file (infile),
    create a tab separated file which is ready to be inserted into the new mongodb as samples collection.
    
    The biosamples_metadata table dump looks like this (it may also contain extra sample_id column which looks like "6003_GSM396481",
    as in the case of the cleaned up file from MDAP):

    chip_type     chip_id                   md_name                              md_value  ds_id
    205  GSM2064216        replicate_group_id  Clec4e-/- microglia I/R, replicate 1   6731
    205  GSM2064216               sample_type               Clec4e-/- microglia I/R   6731
    205  GSM2064216          sample_type_long               Clec4e-/- microglia I/R   6731

    The output file will look like this:

    sample_id   cell_line   parental_cell_type  ...    treatment    external_source_id
    6731_GSM2064216   normal    microglia   ...     I/R                            
    6731_GSM2064205   normal    microglia   ...     sham treatment

    So the workflow is to clean up biosamples_metadata file using OpenRefine, and then use the cleaned up version of the file
    as input to this function to create the outfile, which is then used to insert into the mongodb.

    """
    # columns we will end up in the outfile matrix (sample_id will be the index)
    columns = datasets.Dataset.sample_fields

    # Read infile, which has the long format (md_name, md_value)
    df = pandas.read_csv(infile, sep="\t")

    # Create sample_id column if it doesn't exist already
    if "sample_id" not in df.columns:
        df["sample_id"] = ["%s_%s" % (row['ds_id'], row['chip_id']) for index,row in df.iterrows()]

    # Define new data frame to be filled in
    newdf = pandas.DataFrame(index=df["sample_id"].unique(), columns=columns[1:])
    newdf.index.name = "sample_id"

    # Loop through each row of df and set value in newdf - map some old columns to new
    for index,row in df.iterrows():
        targetColumn = None
        if row["md_name"] in columns:
            targetColumn = row["md_name"]
        elif "tissue" in row["md_name"]:
            targetColumn = "tissue_of_origin"
        elif "media" in row["md_name"]:
            targetColumn = "media"
        elif "age" in row["md_name"] or "time" in row["md_name"]:
            targetColumn = "age_time"
        if targetColumn is None:    # can't find matching column, so skip this record
            continue
        newdf.at[row['sample_id'], targetColumn] = row['md_value']
        newdf.at[row["sample_id"], "dataset_id"] = row["ds_id"]
    
    newdf.to_csv(outfile, sep="\t")

def checkMissingData(publicDatasetsOnly=False):
    """Function to check for inconsistencies, where dataset metadata may have datasets with no corresponding
    sample table, for example.
    """
    import re
    option = {'private':False} if publicDatasetsOnly else {}
    datasetIdsFromMetadata = set([item["dataset_id"] for item in database["datasets"].find(option)])
    datasetIdsFromSamples = set([item["dataset_id"] for item in database["samples"].find({})])
    datasetIdsFromExpression = set()
    for item in os.listdir(os.environ.get("EXPRESSION_FILEPATH")):
        match = re.findall("^\d{4}", item)
        if len(match)>0: datasetIdsFromExpression.add(int(match[0]))

    print("Number of dataset ids from metadata", len(datasetIdsFromMetadata))
    print("Number of dataset ids from samples", len(datasetIdsFromSamples))
    print("Number of dataset ids from expression", len(datasetIdsFromExpression))

    diff1 = datasetIdsFromMetadata.difference(datasetIdsFromSamples)
    diff2 = datasetIdsFromMetadata.difference(datasetIdsFromExpression)
    print("\nDataset ids from metadata missing from samples", len(diff1), list(diff1)[:3])
    print("Dataset ids from metadata missing from expression", len(diff2), list(diff2)[:3])

    diff1 = datasetIdsFromSamples.difference(datasetIdsFromMetadata)
    diff2 = datasetIdsFromSamples.difference(datasetIdsFromExpression)
    print("\nDataset ids from samples missing from metadata", len(diff1), list(diff1)[:3])
    print("Dataset ids from samples missing from expression", len(diff2), list(diff2)[:3])

    diff1 = datasetIdsFromExpression.difference(datasetIdsFromMetadata)
    diff2 = datasetIdsFromExpression.difference(datasetIdsFromSamples)
    print("\nDataset ids from expression missing from metadata", len(diff1), list(diff1)[:3])
    print("Dataset ids from expression missing from samples", len(diff2), list(diff2)[:3])


# ----------------------------------------------------------
# Collection of mostly one-off functions, kept here for history
# ----------------------------------------------------------
def combineSamplesMatricesAndInsertIntoMongo_20210103(fileCath, fileMdap):
    """Read the two tsv files produced by the convertBiosamplesMetadata function above, where one file comes from
    MDAP's cleaned up version of biosamples_metadata, and the other comes from Catherine, who worked on other parts
    not worked on by MDAP. Combine these two files and insert the records into mongodb.
    """
    dfCath = pandas.read_csv(fileCath, sep="\t")
    dfMdap = pandas.read_csv(fileMdap, sep="\t")

    # MDAP samples are better annotated, so drop these from dfCath
    dfCath = dfCath[~dfCath['sample_id'].isin(dfMdap['sample_id'])]

    # combine
    df = pandas.concat([dfMdap, dfCath], axis=0)
    
    # insert
    collection = mongoClient()["dataportal"]["samples"]
    collection.drop()
    collection.insert_many(df.to_dict("records"))

def fixDatasetCollection_20201231():
    """dataset collection was built initially by Jack, and there are some minor fixes required, such as converting
    boolean to proper boolean types instead of using strings, and removing some unwanted keys.
    """
    collection = database["datasets"]
    collection.update_many({'private': "False"},{'$set': {'private':False}}, upsert=False)
    collection.update_many({'private': "True"},{'$set': {'private':True}}, upsert=False)
    collection.update_many({},{'$set': {'version':'1.0'}}, upsert=False)
    collection.update_many({},{'$unset': {'annotator':1}}, upsert=False)
    collection.update_many({},{'$unset': {'can_annotate':1}}, upsert=False)
    collection.update_many({},{'$unset': {'number_of_samples':1}}, upsert=False)
    ds = datasetFromDatasetId(2000)
    print(ds.metadata())

def checkDataFromCatherine_20210101(filepath):
    """Check the data in filepath (given to us by Catherine, Dec 2020), for datasets where only some of the samples
    got cleaned up by MDAP. File content looks like:

    chip_type     chip_id                   md_name                              md_value  ds_id
    205  GSM2064216        replicate_group_id  Clec4e-/- microglia I/R, replicate 1   6731
    205  GSM2064216               sample_type               Clec4e-/- microglia I/R   6731
    205  GSM2064216          sample_type_long               Clec4e-/- microglia I/R   6731

    See notes at the top of this file for the result of running this function.
    """
    # Read file
    df = pandas.read_csv(filepath, sep="\t")

    # Read current samples collection and get sample ids
    currentSampleIds = [item["sample_id"] for item in database["samples"].find({})]
    
    # Subset df on datasets found in currentSampleIds
    currentDatasetIds = set([int(item.split("_")[0]) for item in currentSampleIds])
    subset = df[df["ds_id"].isin(currentDatasetIds)]

    # Group sample ids by dataset id for samples currently in mongo
    currentSampleIdsFromDatasetId = {}
    for item in currentSampleIds:
        dsId,sampleId = item.split("_")[:2]  # remember dsId is string here
        if int(dsId) not in currentSampleIdsFromDatasetId:
            currentSampleIdsFromDatasetId[int(dsId)] = []
        currentSampleIdsFromDatasetId[int(dsId)].append(sampleId)

    # Show datasets where not all samples in the dataset got cleaned up
    for dsId,sampleIds in subset.groupby("ds_id")["chip_id"]:
        if len(set(sampleIds).difference(set(currentSampleIdsFromDatasetId[dsId])))>0:
            print(dsId, len(set(sampleIds)), len(set(currentSampleIdsFromDatasetId[dsId])))


def compareExpressionMatrixLoadingTimes_20210102():
    """Compare the loading time between loading from tsv file vs mongo record.
    Could not complete this function. Kept saying "Killed" for the process, perhaps due to a lack of memory
    during the insert_many process (changing this to a loop of insert_one didn't make any difference).
    """
    # Choose dataset 6646, which is large and has 1448 samples.
    import time
    import datasets
    t0 = time.time()
    df = datasets.Dataset(6646).expressionMatrix("raw")
    print(time.time()-t0)

    df.index.name = "index"
    collection = database["expression"]
    collection.drop()
    collection.insert_many(df.reset_index().to_dict("records"))

    t1 = time.time()
    df = pandas.DataFrame(database["expression"].find({},{"_id":0})).set_index("index")
    print(time.time()-t1)
    print(df.shape)


def updateDataset_20210118():
    """Noticed a dataset which should be public is private. Fixing this.
    """
    result = database["datasets"].find_one({"dataset_id": 7179}, {"_id":0})
    print(result)
    result = database["datasets"].update_one({"dataset_id": 7179}, 
        {"$set": {"private":False, "pubmed_id":"32415101"}}, upsert=False)


if __name__=="__main__":
    if sys.argv[1]=="convertBiosamplesMetadata":
        convertBiosamplesMetadata(sys.argv[2], sys.argv[3])
    elif sys.argv[1]=="checkMissingData":
        checkMissingData(publicDatasetsOnly=len(sys.argv)>2 and sys.argv[2]=='publicOnly')

    elif sys.argv[1]=="combineSamplesMatricesAndInsertIntoMongo_20210103":
        combineSamplesMatricesAndInsertIntoMongo_20210103(sys.argv[2], sys.argv[3])
    elif sys.argv[1]=="fixDatasetCollection_20201231":
        fixDatasetCollection_20201231()
    elif sys.argv[1]=="checkDataFromCatherine_20210101":
        checkDataFromCatherine_20210101(sys.argv[2])
    elif sys.argv[1]=="compareExpressionMatrixLoadingTimes_20210102":
        compareExpressionMatrixLoadingTimes_20210102()
    elif sys.argv[1]=="updateDataset_20210118":
        updateDataset_20210118()
    elif sys.argv[1]=="getMissingExpressionDataFiles_20210120":
        getMissingExpressionDataFiles_20210120()
    else:
        print("Unknown function.")
