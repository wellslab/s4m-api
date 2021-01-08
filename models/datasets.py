"""
Main interface to stemformatics data. Dataset class defined here handles most interactions with underlying
data, including sample tables and expression matrix files. Atlas data are handled separately in the atlas.py.
"""
import pymongo, os, pandas
from models.utilities import mongoClient

database = mongoClient()["dataportal"]

# ----------------------------------------------------------
# Functions
# ----------------------------------------------------------

def datasetFromDatasetId(datasetId):
    """Return a Dataset instance given datasetId. Returns None if no matching Dataset was found
    in the database.
    """
    # Unset _id field before returning the dict, since this field's value is ObjectId class which can't be
    # serialised by json, which can create issues downstream, and this field is not needed outside mongo anyway.
    result = database["datasets"].find_one({"dataset_id": datasetId}, {"_id":0})
    return Dataset(result["dataset_id"], metadata=result) if result else None


# ----------------------------------------------------------
# Dataset class
# ----------------------------------------------------------

class DatasetIdNotFoundError(Exception):
    """Use this when Dataset object has been created with an arbitrary id but there's actually no
    corresponding records in the database.
    """
    pass

class Dataset(object):
    """Main class to encapsulate a dataset in Stemformatics. 
    Each dataset has a unique id, and we can associate 4 types of entities to a dataset:
    1. dataset metadata: information about the dataset, such as description, pubmed_id, etc.
    2. sample metadata: come from sample annotation, such as cell type, tissue, organism, etc.
    3. expression data: may be raw counts, cpm, etc.
    4. processing details: information about how the dataset was processed.
    """

    def __init__(self, datasetId, **kwargs):
        """Initialise a dataset with Id. Note that id is an integer, and will be coherced into one.
        """
        self.datasetId = int(datasetId)
        self._metadata = kwargs.get('metadata')
        self._samples = None
        self._expressionMatrix = {}

        # Make a query to database for dataset metadata now, so if we try to create an instance with 
        # no matching dataset id, we can throw an exception
        if self._metadata is None:  # make a query to db and fetch this, throw exception if no matching datasetId found
            result = database["datasets"].find_one({"dataset_id": self.datasetId}, {"_id":0})
            if not result:
                raise DatasetIdNotFoundError("No matching dataset id found in database:<%s>" % self.datasetId)
            self._metadata = result

    def __repr__(self):
        return "<Dataset id={0.datasetId}>".format(self)

    def __eq__(self, other):
        return self.datasetId==other.datasetId

    def isPrivate(self):
        return self.metadata()["private"]

    # dataset metadata -------------------------------------
    def metadata(self):
        """
        Return dataset metadata, such as description and pubmed id as a dictionary.
        """
        return self._metadata

    # sample metadata -------------------------------------
    def samples(self):
        """
        Return samples in the dataset as a pandas DataFrame object.
        """
        if self._samples is None:   # make a query, construct the DataFrame and cache it
            cursor = database["samples"].find({"dataset_id": self.datasetId}, {"_id":0})
            self._samples = pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()!=0 else pandas.DataFrame()

        return self._samples

    # expression matrix -------------------------------------
    def expressionMatrix(self, key):
        """Return expression matrix for this dataset as a pandas DataFrame.
        """
        if key not in self._expressionMatrix:
            self._expressionMatrix[key] = pandas.read_csv(self.expressionFilePath(key), sep="\t", index_col=0)
        return self._expressionMatrix[key]

    def expressionFilePath(self, key):
        """Return the full path to the expression file.
        """
        return os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s.%s.tsv" % (self.datasetId, key))


# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------
# May have to use export MONGO_URI='xxxx' before running these tests, in order to set environment variables. 
# See .env file for a full list of variables to set.

def test_datasetMetadata():
    assert datasetFromDatasetId(0) is None
    assert datasetFromDatasetId(2000).metadata()["name"] == "Matigian_2010_20699480"

def test_samples():
    df = datasetFromDatasetId(6003).samples()
    print(df.iloc[:3,:4])
