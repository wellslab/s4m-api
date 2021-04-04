"""
Examples of how to run this script (ensure you're in the application directory):
(s4m-api) [ec2-user@api-dev s4m-api]$ python -m scripts.backup_and_restore backupCollectionToCSV dataportal datasets /mnt/stemformatics-data/backups/datasets_20210129.tsv
(s4m-api) [ec2-user@api-dev s4m-api]$ python -m scripts.backup_and_restore createCollectionFromCSV dataportal datasets /mnt/stemformatics-data/backups/datasets_20210129.tsv

In the first example, a tsv file is created from dataportal.datasets collection, while the reverse happens in the second example.
Since this will delete the existing collection first, you have to confirm this if the collection exits. 
"""

import os, sys, pandas

sys.path.append(os.path.join(sys.path[0]))
from models import utilities, datasets

def backupCollectionToCSV(database, collection, filepath):
    """Make a backup of a collection in database to a tsv file at filepath. Example:
    backupCollectionToSCV("dataportal", "samples", "samples_20210102.tsv")
    """
    coll = utilities.mongoClient()[database][collection]
    df = pandas.DataFrame(coll.find({},{"_id":0}))
    df.to_csv(filepath, sep="\t", index=False)

def createCollectionFromCSV(database, collection, filepath):
    """Create a collection in database, given tsv file in filepath. Example:
    createCollectionFromCSV("dataportal", "samples", "samples_20210102.tsv").
    If collection already exists, you'll have to confirm that you want to delete all records in it first.
    """
    coll = utilities.mongoClient()[database][collection]
    df = pandas.read_csv(filepath, sep="\t")

    # For any column which is a list, parse it, otherwise it will go into mongo as a string
    import ast
    listColumns = ['projects']  # all columns where list is expected
    for column in listColumns:
        if column in df.columns:
            df[column] = [ast.literal_eval(item) for item in df[column]]

    if coll.count_documents({})>0:
        confirm = input("This collection contains documents. Are you sure you want to delete all before inserting new documents? (y/[n])\n")
        if confirm!="y": return

    coll.drop()
    coll.insert_many(df.to_dict("records"))
    createTextIndex(database, collection)

def createTextIndex(database, collection):
    """Create text index on a collection.
    """
    if collection=='datasets':
        index_cols = ['title', 'authors', 'description', 'platform']
    elif collection=='samples':
        index_cols = [item for item in datasets.Dataset.sample_fields if item not in ['sample_id']]

    coll = utilities.mongoClient()[database][collection]
    coll.create_index([(item,'text') for item in index_cols])


if __name__=="__main__":
    if sys.argv[1]=="backupCollectionToCSV":
        backupCollectionToCSV(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1]=="createCollectionFromCSV":
        createCollectionFromCSV(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1]=="createTextIndex":
        createTextIndex(sys.argv[2], sys.argv[3])
