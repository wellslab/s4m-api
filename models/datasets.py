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

def datasetIdsFromQuery(**kwargs):
    """Return a list of dataset ids which match a query.
    """
    limit = kwargs.get("limit")
    platform_type = kwargs.get("platform_type")
    params = {"private": False} # return only public datasets currently
    if platform_type:
        params['platform_type'] = platform_type

    if limit:
        cursor = database["datasets"].find(params, {"dataset_id":1, "_id":0}).limit(limit)
    else:
        cursor = database["datasets"].find(params, {"dataset_id":1, "_id":0})

    return [item["dataset_id"] for item in cursor]

def samplesFromQuery(**kwargs):
    """Return DataFrame of samples which match a query.
    """
    limit = kwargs.get("limit")
    cell_type = kwargs.get("cell_type")

    # Only returning samples which belong to public datasets currently
    dsIds = datasetIdsFromQuery()
    params = {"dataset_id": {"$in": dsIds}}

    if cell_type:
        params['cell_type'] = cell_type

    if limit:
        cursor = database["samples"].find(params, {"_id":0}).limit(limit)
    else:
        cursor = database["samples"].find(params, {"_id":0})

    df = pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()!=0 else pandas.DataFrame()
    return df

def allValues(collection, key, includeCount=False):
    """Return a set of all the values for a key in a collection.
    Eg: allValues("samples", "sex") returns {'', 'female', 'male'}
    If includeCount is True, returns a pandas Series that includes count of each value
    """
    # Only returning datasets/samples which belong to public datasets currently
    if collection=="datasets":
        params = {"private": False}
    elif collection=="samples":
        params = {"dataset_id": {"$in": datasetIdsFromQuery()}}
    else: # unknown collection
        return set.Set()

    cursor = database[collection].find(params, {key:1, "_id":0})
    values = ["" if pandas.isnull(item[key]) else item[key] for item in cursor]
    if includeCount:
        return pandas.Series(values).value_counts()
    else:
        return set(values)

# ----------------------------------------------------------
# Dataset class
# ----------------------------------------------------------

class DatasetIdNotFoundError(Exception):
    """Use this when Dataset object has been created with an arbitrary id but there's actually no
    corresponding records in the database.
    """
    pass

class ExpressionFilePathNotFoundError(Exception):
    """Use this when EXPRESSION_FILEPATH is not found in os.environ when trying to access expression data.
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
    
    # This is the full list of fields associated with dataset metadata
    dataset_fields = ["dataset_id", "name", "title", "authors", "description", "platform_type", "platform", 
                      "private", "pubmed_id", "accession", "version"]

    # This is the full list of fields associated with sample metadata
    sample_fields = ["sample_id", "dataset_id", "cell_type", "parental_cell_type", "final_cell_type", "disease_state", 
               "organism", "sample_type", "tissue_of_origin", "sample_name_long", "media", "cell_line", "facs_profile", 
               "sample_description", "age_time", "sex", "reprogramming_method", "genetic_modification",
               "labelling", "developmental_stage", "treatment", "external_source_id"]

    # All available platform_type
    platform_types = ["Microarray", "RNASeq", "scRNASeq", "other"]

    def __init__(self, datasetId):
        """Initialise a dataset with Id. Note that id is an integer, and will be coherced into one.
        """
        self.datasetId = int(datasetId)
        self._samples = None
        self._expressionMatrix = {}

        # Make a query to database for dataset metadata now, so if we try to create an instance with 
        # no matching dataset id, we can throw an exception
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
    def expressionMatrix(self, key="raw"):
        """Return expression matrix for this dataset as a pandas DataFrame.
        key is one of ["raw","gene"]
        """
        if key not in self._expressionMatrix:
            self._expressionMatrix[key] = pandas.read_csv(self.expressionFilePath(key=key), sep="\t", index_col=0)
            # Until we fix columns of expression matrix to match sample_id from database, we need to prefix dataset id
            self._expressionMatrix[key].columns = ["%s_%s" % (self.datasetId, col) for col in self._expressionMatrix[key].columns]
        return self._expressionMatrix[key]

    def expressionFilePath(self, key="raw"):
        """Return the full path to the expression file.
        """
        if "EXPRESSION_FILEPATH" not in os.environ:
            raise ExpressionFilePathNotFoundError("EXPRESSION_FILEPATH not found in os.environ.")
        return os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s/%s.%s.tsv" % (self.datasetId, self.datasetId, key))

    # pca data -------------------------------------
    def pcaCoordinates(self):
        """Return PCA coordinates as a pandas DataFrame object.
        """
        filepath = os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s/%s.pca.tsv" % (self.datasetId, self.datasetId))
        df = pandas.read_csv(filepath, sep="\t", index_col=0) if os.path.exists(filepath) else pandas.DataFrame()

        # Until we fix columns of expression matrix to match sample_id from database, we need to prefix dataset id
        df.index = ["%s_%s" % (self.datasetId, index) for index in df.index]

        return df

    def pcaAttributes(self):
        """Return PCA attributes, such as amount of variance explained by each component, as a pandas DataFrame object.
        """
        filepath = os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s/%s.pca_attributes.tsv" % (self.datasetId, self.datasetId))
        return pandas.read_csv(filepath, sep="\t", index_col=0) if os.path.exists(filepath) else pandas.DataFrame()


# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------
# May have to use export MONGO_URI='xxxx' before running these tests, in order to set environment variables. 
# See .env file for a full list of variables to set.

def test_metadata():
    assert datasetFromDatasetId(0) is None
    assert datasetFromDatasetId(2000).metadata()["name"] == "Matigian_2010_20699480"

def test_samples():
    df = datasetFromDatasetId(6003).samples()
    assert df.shape==(9,21)
    assert df.at["6003_GSM396481", "cell_type"] == "hESC-derived monocyte"

def test_datasetIdsFromQuery():
    print(datasetIdsFromQuery(limit=2))

def test_samplesFromQuery():
    print(samplesFromQuery(limit=5))

def test_allValues():
    print(allValues("samples", "sex"))