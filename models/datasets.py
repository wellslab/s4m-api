"""
Main interface to stemformatics data. Dataset class defined here handles most interactions with underlying
data, including sample tables and expression matrix files. Atlas data are handled separately in the atlas.py.

Note that model here does not enforce rules on private and public datasets, hence it's up to the controller
to ensure that correct parameter is used to call the functions here, if trying to ensure only public datasets
are returned, for example.
"""
import pymongo, os, pandas, numpy
from models.utilities import mongoClient
from models.atlases import Atlas

database = mongoClient()["dataportal"]

# ----------------------------------------------------------
# Functions
# ----------------------------------------------------------
def cpm(df):
    """Return the counts per million version of data frame.
    """
    return (df * 1000000).div(df.sum(axis = 0), axis = 1)

def datasetMetadataFromQuery(**kwargs):
    """Return DataFrame of dataset metadata which match a query. Rows will have dataset ids,
    while columns will be attributes of dataset metadata. Use this instead of Dataset instance
    for fetching large numbers of datasets.
    If ids_only=True, only a list of dataset ids will be returned, instead of a DataFrame.
    Note that query_string will search samples collection as well if include_samples_query=true.
    """
    limit = kwargs.get("limit")
    ids_only = kwargs.get('ids_only', False)
    public_only = kwargs.get("public_only", True)
    include_samples_query = kwargs.get("include_samples_query", False)

    dataset_id = kwargs.get("dataset_id",[]) # list of dataset ids specified in the query
    if dataset_id is None: dataset_id = []
    name = kwargs.get("name")
    query_string = kwargs.get("query_string")
    platform_type = kwargs.get("platform_type")
    projects = kwargs.get("projects")
    organism = kwargs.get("organism")
    status = kwargs.get("status")
    
    params = {}
    attributes = {"dataset_id":1, "_id":0} if ids_only==True else {"_id":0}

    if public_only:
        params['private'] = False

    datasetIds = []  # this is additional dataset ids to search, based on sample search
    if include_samples_query and query_string:  
        # perform text search in both datasets and samples and use union
        sampleSearch = database["samples"].find({'$text': {'$search':query_string}}, {'dataset_id':1})
        datasetIds = [item['dataset_id'] for item in sampleSearch]
        datasetsSearch = database["datasets"].find({'$text': {'$search':query_string}}, {'dataset_id':1})
        datasetIds = list(set(datasetIds).union(set([item['dataset_id'] for item in datasetsSearch])))

    if organism and organism!='all':  # restrict datasets to samples with this organism
        sampleSearch = database["samples"].find({'organism': organism}, {'dataset_id':1})
        datasetIds = list(set(datasetIds).union(set([item['dataset_id'] for item in sampleSearch])))

    if len(dataset_id)>0 and len(datasetIds)>0:  # find common dataset ids
        datasetIds = list(set(datasetIds).intersection(set([int(item) for item in dataset_id])))
    elif len(dataset_id)>0 and len(datasetIds)==0: # just specified by parameter
        datasetIds = [int(item) for item in dataset_id]
    
    if len(datasetIds)>0:
        params['dataset_id'] = {"$in": datasetIds}

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
    if query_string and not include_samples_query:  # otherwise it's been done already above
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

def allValues(collection, key, includeCount=False, public_only=True, excludeDatasets=[], organism='homo sapiens'):
    """Return a set of all the values for a key in a collection.
    Eg: allValues("samples", "sex") returns {'', 'female', 'male'}
    If includeCount is True, returns a pandas Series that includes count of each value.
    excludeDatasets can be a list of dataset ids to exclude explicitly from getting these values.
    """
    params = {}
    if public_only:
        params = {"dataset_id": {"$in": datasetMetadataFromQuery(ids_only=True, public_only=True, organism=organism)}}

    cursor = database[collection].find(params, {key:1, "dataset_id":1, "_id":0})

    # Deal with excludeDatasets
    values = [item.get(key) for item in cursor if len(excludeDatasets)==0 or item['dataset_id'] not in excludeDatasets]

    # Deal with arrays
    values = [','.join(item) if isinstance(item, list) else item for item in values]

    # Deal with nulls
    values = ["" if pandas.isnull(item) else item for item in values]
    if len(set(values))==1 and values[0]=="": # assume no matching key
        return None

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

def dataAsZipfile(datasetIds, publicOnly=True):
    """Return a zip file that contains all relevant files for a list of datasetIds.
    """
    import string, zipfile

    # Create a name for the zip file using string and numpy.random
    randString = ''.join(numpy.random.choice(list(string.ascii_lowercase), size=5))

    # Dataset ids to focus on
    datasetIds = datasetMetadataFromQuery(dataset_id=datasetIds, ids_only=True, public_only=True)

    # Write to file
    filepath = '/tmp/s4m_zipfile_%s.zip' % randString
    with zipfile.ZipFile(filepath, 'w') as zf:
        for datasetId in datasetIds:
            ds = Dataset(datasetId)
            metadata = pandas.DataFrame.from_dict(ds.metadata(), orient='index', columns=['value'])
            metadata.index.name = 'key'
            zf.writestr("%s_samples.tsv" % datasetId, ds.samples().to_csv(sep="\t"))
            zf.writestr("%s_metadata.tsv" % datasetId, metadata.to_csv(sep="\t"))
            if ds.platformType()=='Microarray':
                zf.write(ds.expressionFilePath(key='raw'), "%s_expression_probes.tsv" % datasetId)
                zf.write(ds.expressionFilePath(key='genes'), "%s_expression_genes.tsv" % datasetId)
            else:
                zf.write(ds.expressionFilePath(key='raw'), "%s_expression_counts.tsv" % datasetId)
    return filepath


def sampleSummaryTable():
    """Return a pandas DataFrame which summaries sample data we have in the system.
    Work in progress...
    """
    # Focus on these datasets
    datasetIds = datasetMetadataFromQuery(ids_only=True, public_only=True)

    # Part 1: Count most commonly used cell type values. Since we don't have a very complete sample annotation
    # use the following fields as synonyms to collect together
    keysToUse = ['cell_type','parental_cell_type','final_cell_type','sample_type']
    samples = samplesFromDatasetIds(datasetIds)
    samples = samples[keysToUse]
    
    # For each unique value in this matrix, create a set of sample ids
    valuesToIgnore = ['E7','Day_7','E6','Day_6','E5','Day_5','LPS 24hr','IFNg 24hr','LPS 2hr']
    valuesToCollapse = {'induced pluripotent stem cell':'iPSC', 'mesenchymal stromal cell':'MSC', 'embryonic stem cell':'ESC', 'hiPSC':'iPSC'}
    values = {}
    for column in samples.columns:
        s = samples[column]
        for item in s.unique():
            if item in valuesToIgnore: continue
            sampleIds = set(s[s==item].index.tolist())
            itemToUSe = valuesToCollapse.get(item, item)
            if itemToUSe in values:
                values[itemToUSe] = values[itemToUSe].union(sampleIds)
            else:
                values[itemToUSe] = sampleIds

    # Sort values by number of sample ids and show the highest
    values = [(len(val),key) for key,val in values.items()]
    for item in sorted(values, reverse=True)[:15]:
        print("{sample:'%s', number_of_samples:%s}" % (item[1],item[0]))

    return 'done'

    keysToIgnore = ['sample_id','dataset_id','organism','media','sample_type_long','sample_description',
                    'external_source_id','treatment','sex','labelling','reprogramming_method','sample_name_long',
                    'facs_profile']
    for key in Dataset.sample_fields:
        if key in keysToIgnore: continue
        values = allValues('samples', key, includeCount=True)
        values = values.loc[[index for index in values.index if index!='']].sort_values(ascending=False)
        print(key, values[:3].index.tolist(), values[:3].tolist())

    focusParentKey = 'tissue_of_origin'
    focusParentValues = ['blood','bone marrow']
    focusChildKey = 'cell_type'
    result = {}

    params = {"dataset_id": {"$in": datasetIds}}
    for parentValue in focusParentValues:
        params[focusParentKey] = parentValue
        cursor = database['samples'].find(params, {focusChildKey:1, "_id":0})
        values = [item[focusChildKey] for item in cursor if pandas.notnull(item[focusChildKey])]
        values = pandas.Series(values).value_counts().sort_values(ascending=False)
        result[parentValue] = values[:6]

    return result


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

    # All available platform_type values
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

    def platformType(self):
        return self.metadata()['platform_type']

    # sample metadata -------------------------------------
    def samples(self):
        """
        Return samples in the dataset as a pandas DataFrame object.
        """
        cursor = database["samples"].find({"dataset_id": self.datasetId}, {"_id":0})
        return pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()!=0 else pandas.DataFrame()

    # expression matrix -------------------------------------
    def expressionMatrix(self, key="raw", applyLog2=False):
        """Return expression matrix for this dataset as a pandas DataFrame.
        key is one of ["raw","genes"] for microarray and ["raw","cpm"] for RNASeq.
        Using 'genes' for RNASeq will still work, and fetch 'raw' in that case.
        Using 'cpm' for Microarray will still work, and fetch 'genes' in that case.
        applyLog2 will apply log2(df+1) if platform_type is RNASeq and max value is greater than 100.
        """
        if self.platformType()=='RNASeq' and key=='cpm': # get raw and calculate cpm
            df = pandas.read_csv(self.expressionFilePath(key='raw'), sep="\t", index_col=0)
            df = cpm(df)
        elif self.platformType()=='Microarray' and key=='cpm':
            df = pandas.read_csv(self.expressionFilePath(key='genes'), sep="\t", index_col=0)
        else:
            df = pandas.read_csv(self.expressionFilePath(key=key), sep="\t", index_col=0)

        if applyLog2 and self.platformType()=='RNASeq' and df.max().max()>100: # 
            df = numpy.log2(df+1)

        return df

    def expressionFilePath(self, key="raw"):
        """Return the full path to the expression file.
        key is one of ["raw","genes"] for microarray and ["raw"] for RNASeq.
        Using 'genes' for RNASeq will still work, and fetch 'raw' in that case.
        Using 'cpm' for Microarray will still work, and fetch 'genes' in that case.
        """
        if "EXPRESSION_FILEPATH" not in os.environ:
            raise ExpressionFilePathNotFoundError("EXPRESSION_FILEPATH not found in os.environ.")
        if self.platformType()=='RNASeq' and key=='genes': # make it same as raw
            key = 'raw'
        elif self.platformType()=='Microarray' and key=='cpm': # make it same as genes
            key = 'genes'
        return os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s/%s.%s.tsv" % (self.datasetId, self.datasetId, key))

    # pca data -------------------------------------
    def pcaCoordinates(self):
        """Return PCA coordinates as a pandas DataFrame object.
        """
        filepath = os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s/%s.pca.tsv" % (self.datasetId, self.datasetId))
        df = pandas.read_csv(filepath, sep="\t", index_col=0) if os.path.exists(filepath) else pandas.DataFrame()

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
    assert Dataset(0) is None
    assert Dataset(2000).metadata()["name"] == "Matigian_2010_20699480"
    print()

def test_samples():
    df = datasetFromDatasetId(6003).samples()
    assert df.shape==(9,21)
    assert df.at["6003_GSM396481", "cell_type"] == "hESC-derived monocyte"

def test_datasetMetadataFromQuery():
    df = datasetMetadataFromQuery()
    print(df.shape)
    df['projects'] = ['atlas' if len(item)>0 else '' for item in df['projects']]
    #df = df[df['projects']!='']
    #print(df.shape)
    s = df.groupby('platform_type').size()
    for key,val in s.items():
        print('{platform_type:%s, number_of_datasets:%s}' % (key,val))
    print(df.groupby(['projects','platform_type']).size().to_dict())

def test_datasetMetadataVsDatasetLoadingTime():
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

