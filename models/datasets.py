"""
Main interface to stemformatics data. Dataset class defined here handles most interactions with underlying
data, including sample tables and expression matrix files. Atlas data are handled separately in the atlas.py.

Note that model here does not enforce rules on private and public datasets, hence it's up to the controller
to ensure that correct parameter is used to call the functions here, if trying to ensure only public datasets
are returned, for example.
"""
import pymongo, os, pandas
from models.utilities import mongoClient
from models.atlases import Atlas

database = mongoClient()["dataportal"]

# ----------------------------------------------------------
# Functions
# ----------------------------------------------------------

def datasetMetadataFromQuery(**kwargs):
    """Return DataFrame of dataset metadata which match a query. Rows will have dataset ids,
    while columns will be attributes of dataset metadata. Use this instead of Dataset instance
    for fetching large numbers of datasets.
    If ids_only=True, only a list of dataset ids will be returned, instead of a DataFrame.
    Note that currently query_string is searching only in datasets collection, rather than including samples.
    This should be changed in future once samples is better annotated.
    """
    limit = kwargs.get("limit")
    ids_only = kwargs.get('ids_only', False)
    public_only = kwargs.get("public_only", True)

    dataset_id = kwargs.get("dataset_id")
    name = kwargs.get("name")

    query_string = kwargs.get("query_string")
    platform_type = kwargs.get("platform_type")
    projects = kwargs.get("projects")
    status = kwargs.get("status")

    params = {}
    attributes = {"dataset_id":1, "_id":0} if ids_only==True else {"_id":0}

    if public_only:
        params['private'] = False
    if dataset_id:
        params['dataset_id'] = {"$in": [int(item) for item in dataset_id]}
    if platform_type:
        params['platform_type'] = platform_type
    if projects:
        if projects=='atlas':  # any atlas project
            params["projects"] = {"$in":["%s_atlas" % atlasType for atlasType in Atlas.all_atlas_types]}
        else:
            params["projects"] = {"$in":[projects]}
    if status:
        params["status"] = status
    if name:
        params["name"] = name
    if query_string:
        params['$text'] = {"$search": query_string}

    if limit:
        cursor = database["datasets"].find(params, attributes).limit(limit)
    else:
        cursor = database["datasets"].find(params, attributes)
    
    if ids_only:
        return [item["dataset_id"] for item in cursor]
    else:
        return pandas.DataFrame(cursor).set_index("dataset_id") if cursor.count()!=0 else pandas.DataFrame()

def samplesFromDatasetIds(datasetIds):
    """Return DataFrame of samples which belong to datasets with datasetIds.
    """
    params = {"dataset_id": {"$in": datasetIds}}
    cursor = database["samples"].find(params, {"_id":0})
    return pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()!=0 else pandas.DataFrame()

def allValues(collection, key, includeCount=False, public_only=True):
    """Return a set of all the values for a key in a collection.
    Eg: allValues("samples", "sex") returns {'', 'female', 'male'}
    If includeCount is True, returns a pandas Series that includes count of each value
    """
    params = {}
    if public_only:
        params = {"dataset_id": {"$in": datasetMetadataFromQuery(ids_only=True, public_only=True)}}

    cursor = database[collection].find(params, {key:1, "_id":0})
    values = ["" if pandas.isnull(item[key]) else item[key] for item in cursor]
    if includeCount:
        return pandas.Series(values).value_counts()
    else:
        return set(values)

def sunburstData(dataFrame, childKey='final_cell_type', parentKey='parental_cell_type',
                   parentCutoff=12, childCutoff=16, sep='_', includeOther=False):
    """Return a pandas DataFrame that can be used as input to sunburst plot.
    """
    topParents = dataFrame[parentKey].value_counts()[:parentCutoff].index
    
    # Find children with highest numbers from topParent
    df = dataFrame[dataFrame[parentKey].isin(topParents)]
    children = df[childKey].value_counts()[:childCutoff].index
                
    subset = df if includeOther else df[df[childKey].isin(children)]
    # If a parent has only one child, remove this key
    s = subset.groupby(parentKey)[childKey].unique()  
        # - creates a series with parentKey value as index, and list of uniquely matching childKey values
    parentsToKeep = [key for key,val in s.items() if len(val)>1]
    subset = subset[subset[parentKey].isin(parentsToKeep)]

    # Create a dictionary that counts number of samples for each child, which is actually
    # a combination of parent_child, because child may be duplicated (same child occurs under a different parent)
    # and each child must have a unique id.
    combo = {}
    for index,row in subset.iterrows():
        child = row[childKey] if not includeOther or row[childKey] in children else 'other'
        key = "%s%s%s" % (row[parentKey], sep, child)
        if key not in combo:
            combo[key] = []
        combo[key].append(index)

    # combo.keys() gives us the id to use for sunburst labels. 
    # Create a data frame we can use for sunburst plot.
    df = pandas.DataFrame(index=list(combo.keys()), columns=['labels','parents','values'])
    df.index.name = 'ids'
    for key,val in combo.items():
        df.at[key,'labels'] = key.split(sep)[1]
        df.at[key,'parents'] = key.split(sep)[0]
        df.at[key,'values'] = val
    
    # if a parent is not found existing childKey though, we need to add them also
    parentsWithoutLabels = [item for item in df['parents'].unique() if item not in df.index]
    for item in parentsWithoutLabels:
        df.at[item] = [item, '', sum(df[df['parents']==item]['values'].tolist(),[])]

    return df

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
        cursor = database["samples"].find({"dataset_id": self.datasetId}, {"_id":0})
        return pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()!=0 else pandas.DataFrame()

    # expression matrix -------------------------------------
    def expressionMatrix(self, key="raw"):
        """Return expression matrix for this dataset as a pandas DataFrame.
        key is one of ["raw","gene"] for microarray and ["raw","cpm"] for rna-seq.
        """
        if self.metadata()['platform_type']=='RNASeq' and key=='cpm': # get raw and calculate cpm
            df = pandas.read_csv(self.expressionFilePath(key='raw'), sep="\t", index_col=0)
            df = (df * 1000000).div(df.sum(axis = 0), axis = 1)    # counts per million
        else:
            df = pandas.read_csv(self.expressionFilePath(key=key), sep="\t", index_col=0)

        # Until we fix columns of expression matrix to match sample_id from database, we need to prefix dataset id
        # and remove .CEL suffixes
        df.columns = ["%s_%s" % (self.datasetId, col.replace(".CEL","")) for col in df.columns]

        return df

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

def test_datasetMetadataFromQuery():
    """Compare times for bulk query in mongo vs constructing a data frame after individual queries
    (457, 12) 0.012085914611816406
    (457, 13) 0.2130753993988037
    """
    import time
    t0 = time.time()
    df = datasetMetadataFromQuery(limit=500)
    print(df.shape, time.time()-t0)
    t1 = time.time()
    print(pandas.DataFrame.from_records([Dataset(datasetId).metadata() for datasetId in df.index]).shape, time.time()-t1)

def test_samplesFromQuery():
    print(samplesFromQuery(limit=5))

def test_allValues():
    print(allValues("samples", "sex"))