"""
Main interface to stemformatics data. Dataset class defined here handles most interactions with underlying
data, including sample tables and expression matrices. Atlas data are handled separately in the atlas.py.

Note that model here does not enforce rules on private and public datasets, hence it's up to the controller
to ensure that correct parameter is used to call the functions here, if trying to ensure only public datasets
are returned, for example.

The following datasets had expression files but no metadata, which means they're likely to be irrelevant
for the new system (it implies that we copied the expression files across from various sources but there was
no entry for them in the old sql tables). But we keep a record of them here to examine later. When we converted
the expression from text files to anndata objects, these were not included, so their expression values are
only held in text files.
[5010, 5037, 6031, 6054, 6060, 6097, 6105, 6121, 6159, 6165, 6191, 6192, 6204, 6209, 6211, 6212, 6214, 6215, 
 6217, 6219, 6220, 6234, 6235, 6240, 6241, 6243, 6246, 6256, 6257, 6259, 6262, 6294, 6299, 6301, 6304, 6305, 
 6312, 6314, 6330, 6333, 6341, 6352, 6356, 6361, 6365, 6372, 6374, 6384, 6387, 6397, 6400, 6401, 6415, 6418, 
 6419, 6421, 6422, 6423, 6424, 6425, 6426, 6427, 6428, 6437, 6441, 6442, 6460, 6464, 6469, 6472, 6488, 6493, 
 6494, 6504, 6508, 6510, 6511, 6519, 6520, 6524, 6525, 6535, 6537, 6538, 6539, 6547, 6548, 6549, 6550, 6551, 
 6553, 6554, 6563, 6564, 6568, 6581, 6586, 6587, 6595, 6616, 6626, 6630, 6636, 6640, 6645, 6647, 6648, 6650, 
 6651, 6667, 6668, 6681, 6689, 6690, 6691, 6692, 6697, 6707, 6709, 6720, 6722, 6726, 6728, 6738, 6739, 6743, 
 6744, 6759, 6800, 6803, 6833, 6840, 6841, 6856, 6888, 6962, 6964, 6985, 7012, 7069, 7070, 7097, 7110, 7156, 
 7157, 7190, 7255, 7257, 7285, 7286, 7292, 7327, 7346, 7377]

"""
import pymongo, os, pandas, numpy, anndata
from models.utilities import mongoClient
from models.atlases import Atlas

database = mongoClient()["dataportal"]

# These datasets have entries in the metadata table but are not ready to be exposed to the public - such as 6131, 
# which is not in the normal format. Hence by default, datasetMetadataFromDatasetIds() function here 
# will exclude this hard coded list of datasets. You can still access these by using Dataset() initialiser.
_exclude_list = [5002, 6056, 6127, 6130, 6131, 6149, 6150, 6151, 6155, 6187, 6197, 6198, 6368, 6655, 6701, 6754, 6776, 6948, 7012, 7115, 7209, 7217, 7218, 7250, 7311, 7401]

# ----------------------------------------------------------
# Functions
# ----------------------------------------------------------
def cpm(df):
    """Return the counts per million version of data frame.
    """
    return (df * 1000000).div(df.sum(axis = 0), axis = 1)

def datasetMetadataFromDatasetIds(datasetIds, publicOnly=True):
    """Given a list of dataset ids, return a DataFrame where rows will be dataset ids,
    while columns will be attributes of dataset metadata. Use this instead of Dataset instance
    for fetching large numbers of datasets. Also the most "primitive" of all datasetMetadata
    functions here, in the sense that all the others will call this after working out the
    dataset ids first.
    """
    params = {'dataset_id': {"$nin": _exclude_list, "$in":datasetIds}}
    if publicOnly:
        params['private'] = False
    cursor = database["datasets"].find(params, {"_id":0})
    return pandas.DataFrame(cursor).set_index("dataset_id") if cursor.count()>0 else pandas.DataFrame()

def datasetIdsFromQuery(query_string, include_samples_query=False):
    """Return dataset ids which match a query.
    Note that samples collection will be searched as well if include_samples_query=true.
    Use datasetIdsFromQuery('*') to fetch all dataset ids (including private dataset ids)
    """
    # params for find function (ie. fetch all records matching params) and attributes for what to return
    params = {} if query_string=='*' else {'$text': {'$search':query_string}}
    datasetIds = []

    if include_samples_query:
        # perform text search in both datasets and samples and use union
        sampleSearch = database["samples"].find(params, {'dataset_id':1})
        datasetIds = [item['dataset_id'] for item in sampleSearch]

    datasetsSearch = database["datasets"].find(params, {'dataset_id':1})
    return list(set(datasetIds).union(set([item['dataset_id'] for item in datasetsSearch])))

def datasetIdsFromFields(platform_type=[], projects=[], organism=['homo sapiens'], status=[], publicOnly=True):
    """Return dataset ids which match values specified. The query is an 'and' query for all fields.
    Call the function with default values to get all public human datasets.
    """
    # Get dataset ids matching organism first
    datasetIds = None
    if len(organism)>0 and 'all' not in organism:  # restrict datasets to samples with this organism
        sampleSearch = database["samples"].find({'organism': {'$in':organism}}, {'dataset_id':1})
        datasetIds = set([item['dataset_id'] for item in sampleSearch])
    
    if datasetIds is None or len(datasetIds)>0: # Get dataset ids matching other fields
        params = {}
        if len(platform_type)>0:
            params['platform_type'] = {"$in": platform_type}
        if len(projects)>0:
            if projects==['atlas']:  # any atlas project
                params["projects"] = {"$in":["%s_atlas" % atlasType for atlasType in Atlas.all_atlas_types]}
            else:
                params["projects"] = {"$in":projects}
        if publicOnly:
            params['private'] = False

        if params:
            cursor = database["datasets"].find(params, {"dataset_id":1, "_id":0})
            ids = set([item["dataset_id"] for item in cursor])
            datasetIds = ids if datasetIds is None else datasetIds.intersection(ids)

    return [] if datasetIds is None else list([dsId for dsId in datasetIds if dsId not in _exclude_list])

def datasetIdFromName(name, publicOnly=True):
    params = {'name':name}
    if publicOnly:
        params['private'] = False
    cursor = database["datasets"].find(params, {"dataset_id":1, "_id":0})
    return [item['dataset_id'] for item in cursor][0] if cursor.count()>0 else None

def samplesFromQuery(datasetIds=[], queryString="", organism=["homo sapiens"], limit=100, publicOnly=True):
    """Return a pandas DataFrame containing sample table which match the query.
    """    
    params = {}
    if queryString and queryString!='*':
        params['$text'] = {'$search':queryString}
    if datasetIds:
        params['dataset_id'] = {'$in':datasetIds}
    if organism and 'all' not in organism:
        params['organism'] = {'$in':organism}

    if limit:
        cursor = database["samples"].find(params, {'_id':0}).limit(limit)
    else:
        cursor = database["samples"].find(params, {'_id':0})

    df = pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()>0 else pandas.DataFrame()
    
    if publicOnly and len(df)>0:  # subset df based on dataset property
        cursor = database["datasets"].find({'private':False}, {"dataset_id":1, "_id":0})
        ids = [item['dataset_id'] for item in cursor]
        df = df[df['dataset_id'].isin(ids)]

    return df

def samplesFromDatasetIds(datasetIds):
    """Return DataFrame of samples which belong to datasets with datasetIds.
    """
    params = {"dataset_id": {"$in": datasetIds}}
    cursor = database["samples"].find(params, {"_id":0})
    return pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()!=0 else pandas.DataFrame()

def allValues(collection, key, includeCount=False, public_only=True, organism='homo sapiens', excludeDatasets=[]):
    """Return a set of all the values for a key in a collection.
    Eg: allValues("samples", "sex") returns {'', 'female', 'male'}
    If includeCount is True, returns a pandas Series that includes count of each value.
    excludeDatasets can be a list of dataset ids to exclude explicitly from getting these values.
    """
    # datasets with matching public_only and organism
    datasetIds = datasetIdsFromFields(publicOnly=public_only, organism=organism.split(','))
    if len(datasetIds)==0: # must be due to organism not matching
        return pandas.Series([]) if includeCount else set()
    params = {'dataset_id':{'$in':datasetIds}} if organism else {}

    # Make query for key
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

def sunburstData(samples, childKey='final_cell_type', parentKey='parental_cell_type',
                   parentCutoff=12, childCutoff=8, sep='_', includeOther=False):
    """Return a pandas DataFrame that can be used as input to sunburst plot.
    samples is a samples table (in same format as Dataset.samples()).
    """
    topParents = samples[parentKey].value_counts().sort_values(ascending=False)[:parentCutoff].index
    
    # Find children with highest numbers from topParent
    df = samples[samples[parentKey].isin(topParents)]
    children = df[childKey].value_counts().sort_values(ascending=False)[:childCutoff].index
                
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

    # Filter out any
    if publicOnly:
        datasetIds = set(datasetIds).intersection(set(datasetIdsFromFields(organism=['all'], publicOnly=publicOnly)))

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
    datasetIds = datasetIdsFromFields()

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
               "sample_description", "experiment_time", "sex", "reprogramming_method", "genetic_modification",
               "sample_source", "developmental_stage", "treatment", "external_source_id"]

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
        key may be one of ['raw','genes','cpm'].

        For Microarray data, 'cpm' is treated as 'genes', as there is no cpm calculation done to them.
        For RNASeq data, 'raw' and 'genes' are the same, while 'cpm' calculates cpm values.

        applyLog2 will apply log2(df+1) if platform_type is RNASeq and max value is greater than 100.
        """
        # First get filepath to the expression matrix - always fetch h5 file if we can for speed
        isMicroarray = self.platformType()=='Microarray'
        filepath = ''
        if isMicroarray: # fetch genes.h5 unless key='raw'
            if key=='cpm': key = 'genes'
            filepath = self.expressionFilePath(key=key, hdf5=key=='genes')
        else:
            if key=='raw': key = 'genes'
            filepath = self.expressionFilePath(hdf5=True)

        if filepath.endswith('.h5'):
            df = pandas.read_hdf(filepath, key='genes')
        else:
            df = pandas.read_csv(filepath, sep="\t", index_col=0)

        if not isMicroarray and key=='cpm':
            df = cpm(df)

        if applyLog2 and not isMicroarray and df.max().max()>100: # 
            df = numpy.log2(df+1)

        return df

        # ada = anndata.read_h5ad(f'/mnt/stemformatics-data/expression-files/{self.datasetId}_1.0.h5ad')
        # df = ada.to_df().transpose()
        # if 'RNASeq' in self.platformType() and key=='cpm': # get raw and calculate cpm
        #     df = cpm(df)
        # if applyLog2 and 'RNASeq' in self.platformType() and df.max().max()>100: # 
        #     df = numpy.log2(df+1)
        # return df
        
    def expressionFilePath(self, key="genes", hdf5=False):
        """Return the full path to the expression file. key may be one of ["raw","genes"].
        For RNASeq data, key is ignored, so key=genes is same as raw.
        If hdf5=True, will return the filepath to the .h5 file instead of .tsv (only for key=genes)
        """
        if "EXPRESSION_FILEPATH" not in os.environ:
            raise ExpressionFilePathNotFoundError("EXPRESSION_FILEPATH not found in os.environ.")
        
        filesuffix = ""
        if self.platformType()=='Microarray':
            if key=='raw' and hdf5:
                raise Exception("There is no h5 version of raw file for microarray data.")
            elif key=='genes' and hdf5:
                filesuffix = 'genes.h5'
            else:
                filesuffix = f'{key}.tsv'
        else:
            filesuffix = 'genes.h5' if hdf5 else 'raw.tsv'

        return os.path.join(os.environ["EXPRESSION_FILEPATH"], f"{self.datasetId}/{self.datasetId}.{filesuffix}")

    # pca data -------------------------------------
    def pcaCoordinates(self):
        """Return PCA coordinates as a pandas DataFrame.
        """
        filepath = os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s/%s.pca.tsv" % (self.datasetId, self.datasetId))
        df = pandas.read_csv(filepath, sep="\t", index_col=0) if os.path.exists(filepath) else pandas.DataFrame()

        return df

    def pcaAttributes(self):
        """Return PCA attributes, such as amount of variance explained by each component, as a pandas DataFrame object.
        """
        filepath = os.path.join(os.environ["EXPRESSION_FILEPATH"], "%s/%s.pca_attributes.tsv" % (self.datasetId, self.datasetId))
        return pandas.read_csv(filepath, sep="\t", index_col=0) if os.path.exists(filepath) else pandas.DataFrame()

    # Analysis functions -------------------------------------
    def correlatedGenes(self, geneId, cutoff=30):
        """Return correlation as a pandas Series.
        """
        df = self.expressionMatrix(key='cpm')
        if geneId not in df.index:
            return None
        values = df.loc[geneId]
        
        # Speed up the query by filtering out genes with low variance?
        # What if geneId is one of the genes with low variance?
        #var = df.var(axis=1)
        #df = df.loc[var[var>1].index]

        corr = df.corrwith(values, axis=1).sort_values(ascending=False)
        return corr[:cutoff]

    def ttest(self, geneId, sampleGroup, sampleGroupItems):
        """Return the result of running T-test between elements of sampleGroupItems.
        The result looks like: {statistic:-1.0418722798394733, pvalue:0.3220043225434629}
        """
        from scipy.stats import ttest_ind
        df = self.expressionMatrix(key='cpm', applyLog2=True)
        if geneId not in df.index:
            return None
        samples = self.samples()
        values = df.loc[geneId, samples.index]
        sampleIds1 = samples[samples[sampleGroup]==sampleGroupItems[0]].index
        sampleIds2 = samples[samples[sampleGroup]==sampleGroupItems[1]].index
        result = ttest_ind(values[sampleIds1], values[sampleIds2])

        return {'statistic': result.statistic, 'pvalue': result.pvalue}

# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------
# May have to use export MONGO_URI='xxxx' before running these tests, in order to set environment variables. 
# See .env file for a full list of variables to set.

def test_metadata():
    assert Dataset(2000).metadata()["name"] == "Matigian_2010_20699480"
    assert Dataset(6864).metadata()['platform'] == "Illumina MouseRef-8 V2 (GPL6885 and A-MEXP-1174)"

def test_samples():
    df = Dataset(6003).samples()
    assert set(df.reset_index().columns)==set(Dataset.sample_fields)
    assert df.shape==(9,21)
    assert df.at["6003_GSM396481", "cell_type"]=='monocyte'

def test_expressionFilepath():
    ds = Dataset(6003)
    assert os.path.basename(ds.expressionFilePath(key='raw'))=='6003.raw.tsv'
    assert os.path.basename(ds.expressionFilePath(key='genes'))=='6003.genes.tsv'
    assert os.path.basename(ds.expressionFilePath(key='genes', hdf5=True))=='6003.genes.h5'
    ds = Dataset(7419)
    assert os.path.basename(ds.expressionFilePath(key='raw'))=='7419.raw.tsv'
    assert os.path.basename(ds.expressionFilePath(key='genes'))=='7419.raw.tsv'
    assert os.path.basename(ds.expressionFilePath(key='raw', hdf5=True))=='7419.genes.h5'
    assert os.path.basename(ds.expressionFilePath(hdf5=True))=='7419.genes.h5'

def test_datasetMetadataFromDatasetIds():
    df = datasetMetadataFromDatasetIds([2000, _exclude_list[0]])
    assert df.shape==(1,12)
    df = datasetMetadataFromDatasetIds([])
    assert df.shape==(0,0)

def test_datasetIdsFromQuery():
    datasetIds = datasetIdsFromQuery(query_string='xdfdfdx', include_samples_query=True)
    assert len(datasetIds)==0
    datasetIds = datasetIdsFromQuery('')
    assert len(datasetIds)==0
    datasetIds = datasetIdsFromQuery('*')
    assert len(datasetIds)==661
    datasetIds = datasetIdsFromQuery(query_string='abud', include_samples_query=True)
    assert datasetIds==[7268]

def test_datasetIdsFromFields():
    datasetIds = datasetIdsFromFields()
    assert len(datasetIds)==338
    datasetIds = datasetIdsFromFields(organism=['mus musculus'])
    assert len(datasetIds)==129

def test_datasetIdFromName():
    assert datasetIdFromName('Abud_2017_28426964')==7268

def test_expression():
    df = Dataset(2000).expressionMatrix(key='raw')
    assert df.index[0].startswith('ILMN')
    df = Dataset(2000).expressionMatrix(key='cpm')
    assert df.index[0].startswith('ENSG')
    df = Dataset(7419).expressionMatrix(key='raw')
    assert df.index[0].startswith('ENSG') and df.iloc[0,0]==61
    df = Dataset(7419).expressionMatrix(key='cpm', applyLog2=True)
    assert df.max().max()<15

def test_ttest():
    ds = Dataset(8144) 
    assert ds.ttest('ENSG00000197576', 'cell_type', ['conventional dendritic cell', 'macrophage'])['pvalue']<0.03
    assert ds.ttest('ENSG00000150048', 'sample_type', ['CD14+ cell', 'cDC1'])['statistic']<-16
    assert ds.ttest('ENSG00000172954', 'sample_type', ['CD14+ cell', 'cDC1'])['pvalue']>0.3

def test_datasetMetadataVsDatasetLoadingTime():
    """Compare times for bulk query in mongo vs constructing a data frame after individual queries
    (457, 12) 0.012085914611816406
    (457, 13) 0.2130753993988037
    """
    return
    import time
    t0 = time.time()
    df = datasetMetadataFromDatasetIds(datasetIdsFromQuery())
    print(df.shape, time.time()-t0)
    t1 = time.time()
    print(pandas.DataFrame.from_records([Dataset(datasetId).metadata() for datasetId in df.index]).shape, time.time()-t1)

