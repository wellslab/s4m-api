"""
Main interface to gene expression data. Even though we can consider gene expression as a part of dataset,
performing gene expression analyses across multiple datasets happens often, so it gets its own package here.

# How to use mygene python package
    mg = mygene.MyGeneInfo()
    res = mg.querymany(geneSymbols, ensemblonly=True, scopes='symbol', fields=['symbol','ensembl.gene'], species='human')
    # Note that some ensembl entries are multiple matches
    geneIds = {}
    for item in res:
        if isinstance(item['ensembl'], list):
            geneIds[item['symbol']] = [val['gene'] for val in item['ensembl']]
        else:
            geneIds[item['symbol']] = [item['ensembl']['gene']]
    return geneIds

"""
import requests, os, pandas, numpy, json, anndata
from models import datasets
from models.utilities import mongoClient

def sampleGroupToGenes(sampleGroup, sampleGroupItem, sampleGroupItem2=None, cutoff=20):
    """Given a sampleGroup and sampleGroupItem, (eg cell_type=monocyte), loop through each dataset which
    contains this sample and calculate genes for high expression.
    """
    if sampleGroupItem2=='': sampleGroupItem2 = None

    # Find records in samples collection matching sampleGroupItem - just get dataset ids for now
    cursor = mongoClient()["dataportal"]["samples"].find({sampleGroup: sampleGroupItem}, {"_id":0, "dataset_id":1})
    datasetIds = [item['dataset_id'] for item in cursor]

    # Restrict to public human datasets, Microarray and RNASeq only
    datasetIds = list(set(datasets.publicDatasetIds()).intersection(set(datasetIds)))

    # Get all samples for these datasetIds - quicker to make one mongo query than to loop through Dataset object
    cursor = mongoClient()["dataportal"]["samples"].find({'dataset_id': {'$in':datasetIds}}, {"_id":0})
    allSamples = pandas.DataFrame(cursor).set_index("sample_id") if cursor.count()!=0 else pandas.DataFrame()

    # This will hold genes as index and rank score for each gene in each dataset
    rankScore = pandas.Series(dtype=float)
    datasetIds = {}
    uniqueDatasetIds = set()  # keep track of all dataset ids used for scoring

    for datasetId in allSamples['dataset_id'].unique():
        samples = allSamples[allSamples['dataset_id']==datasetId]
        
        # ignore dataset if there's only one sampleGroupItem as we can't make comparisons
        if len(samples[sampleGroup].unique())==1: continue
        
        if not sampleGroupItem2 is None: # further ignore if this isn't found in the dataset
            # we could do this in the mongo search above but that first search is fast enough
            if not sampleGroupItem2 in samples[sampleGroup].tolist(): continue

        exp = datasets.Dataset(datasetId).expressionMatrix(key='cpm', applyLog2=True)
        exp = exp[samples.index]
        
        # Calculate the difference between mean of sampleGroupItem samples vs max of other in sampleGroup
        df1 = exp[samples[samples[sampleGroup]==sampleGroupItem].index].mean(axis=1)
        if sampleGroupItem2 is None: 
            df2 = exp[samples[samples[sampleGroup]!=sampleGroupItem].index].max(axis=1)
        else:
            df2 = exp[samples[samples[sampleGroup]==sampleGroupItem2].index].max(axis=1)
        diff = df1 - df2

        # Only keep +ve scores
        diff = diff[diff>0]

        # Remember dataset id for all the genes in diff
        for geneId in diff.index:
            if geneId not in datasetIds: datasetIds[geneId] = []
            datasetIds[geneId].append(datasetId)
            uniqueDatasetIds.add(datasetId)

        # Row concatenate this Series for all datasets after converting diff into normalised ranks 
        # (so 0 is lowest diff value, 1 is highest)
        rankScore = pandas.concat([rankScore, diff.rank()/len(diff)])

    # Create data frame of average rank values
    df = pandas.DataFrame({'meanRank': rankScore.groupby(rankScore.index).mean()})

    if len(df)==0:
        return {'rankScore':pandas.DataFrame(), 'totalDatasets':len(uniqueDatasetIds)}

    # Add datasetIds and count of them
    df['datasetIds'] = [','.join(map(str,datasetIds[geneId])) for geneId in df.index]
    df['count'] = [len(datasetIds[geneId]) for geneId in df.index]
    df.index.name = "geneId"

    # Apply cutoff - this is first applied to each combination of geneId-count
    if cutoff:
        df = pandas.concat([df[df['count']==count].sort_values('meanRank', ascending=False).iloc[:cutoff,:] for count in df['count'].unique()])

    df = df.sort_values(['count','meanRank'], ascending=False)

    # apply cutoff again, just to make it less confusing for the output
    if cutoff:
        df = df.iloc[:cutoff,:]
     
    return {'rankScore':df, 'totalDatasets':len(uniqueDatasetIds)}

def geneToSampleGroups(geneId, sampleGroup='cell_type'):
    # filepath = '/mnt/stemformatics-data/backups/gene_to_sample_groups.h5ad'
    # ada = anndata.AnnData(obsm=pandas.DataFrame(columns=['cell_type', 'datasetId']))
    # ada.write(filepath)
    # return
    # ada = anndata.read_h5ad(filepath)

    # Restrict to public human datasets, Microarray and RNASeq only
    allDatasetIds = datasets.publicDatasetIds()

    # Also no need to look at cell types with only a few samples assigned to them
    #sampleCount = datasets.allValues('samples', sampleGroup, includeCount=True, excludeDatasets=datasets._exclude_list)
    #sampleGroupItems =  sampleCount[sampleCount>10].index.tolist()

    # DataFrame to hold the result
    result = pandas.DataFrame()

    for datasetId in allDatasetIds:
        ds = datasets.Dataset(datasetId)
        samples = ds.samples()
        samples = samples[samples[sampleGroup].notnull()] # focus on only non-null cell types
        #samples = samples[samples[sampleGroup].isin(sampleGroupItems)]  # other sampleGroupItems are too infrequent
        if len(samples[sampleGroup].unique())<2:  # we need at least 2 different cell types which aren't null
            continue
        df = ds.expressionMatrix(key='cpm')
        if geneId not in df.index or len(df.columns.intersection(samples.index))==0:
            continue

        # get values of the gene aligned to samples
        values = df.loc[geneId][samples.index]

        mean = values.groupby(samples[sampleGroup]).mean().round()
        var = mean.var()
        if var<1: continue

        # only keep those above median
        diff = mean - mean.median()
        diff = diff[diff>0]

        # rank = mean.rank()/len(mean)  # {'fibroblast': 0.5, 'peripheral blood mononuclear cell': 1.0}
        # rank = rank.sort_values(ascending=False)

        # Create a data frame with all the info and row concatenate this for all datasets
        df = pandas.DataFrame({'score':diff, 'datasetIds':[datasetId for item in diff.index]})
        result = pandas.concat([result,df])

    # Transform result by grouping sampleGroupItem values
    result = result.groupby(result.index).agg({'score':list, 'datasetIds':list})

    # Add other properties
    result['count'] = [len(item) for item in result['datasetIds']]

    return result

def createGeneset(msigName=None, name=None, geneSymbols=[], sampleGroupItem='cell_type'):
    """Create a geneset from some source. If msigName, use name of MsigDB from Broad Institute.
    """
    if "GENESET_FILEPATH" not in os.environ:
        print("GENESET_FILEPATH not found in os.environ.")
        return

    if msigName is not None:
        # First use the api from msigdb to get all gene symbols from the name
        name = msigName
        url = 'https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=%s&fileType=txt' % name
        r = requests.get(url)
        geneSymbols = r.text.split('\n')[2:]
    elif len(geneSymbols)==0:
        print("No MsigDB geneset name or a list of gene symbols provided")
        return

    # Now use mygene.info to get ensembl ids
    params = {'q':','.join(geneSymbols),'scopes':'symbol', 'species':'human', 'fields':'symbol,ensembl.gene', 'ensemblonly':True}
    res = requests.post('http://mygene.info/v3/query', params)

    # Create dict of gene ids keyed on symbols. Note that some ensembl entries are multiple matches
    geneIdsFromSymbol = {}
    for item in res.json():
        if 'ensembl' not in item:
            continue
        if isinstance(item['ensembl'], list):
            geneIdsFromSymbol[item['symbol']] = [val['gene'] for val in item['ensembl']]
        else:
            geneIdsFromSymbol[item['symbol']] = [item['ensembl']['gene']]

    # Loop through each dataset in the system
    scores = pandas.DataFrame(columns=['var','high','low','diff'])
    scores.index.name = 'datasetId'
    allDatasets = datasets.datasetMetadataFromQuery(ids_only=True, organism='homo sapiens')
    for datasetId in allDatasets:
        ds = datasets.Dataset(datasetId)
        if not (ds.platformType()=='RNASeq' or ds.platformType()=='Microarray'): continue
        samples = ds.samples()
        if len(samples[sampleGroupItem].unique())==1:  # some datasets may have one sampleGroupItem specified only
            continue
        # Remove any samples where there is only one sample for sampleGroupItem
        vc = samples[sampleGroupItem].value_counts()
        goodSamples = vc[vc>1].index
        samples = samples[samples[sampleGroupItem].isin(goodSamples)]
        
        exp = ds.expressionMatrix(key='cpm', applyLog2=True)
        # Take mean of geneset and mean across sampleGroupItem, and evaluate variance of these means
        df = exp.loc[exp.index.intersection(sum(list(geneIdsFromSymbol.values()),[]))]
        if len(samples.index.intersection(df.columns))==0: continue
        if len(df)>0:
            group = df[samples.index].groupby(samples[sampleGroupItem], axis=1)
            mean = group.mean().mean()
            variance = mean.var()
            if pandas.notnull(variance):
                scores.loc[datasetId,'var'] = variance
                meanSorted = mean.sort_values()
                scores.loc[datasetId,'high'] = meanSorted.index[-1]
                scores.loc[datasetId,'low'] = meanSorted.index[0]
                scores.loc[datasetId,'diff'] = meanSorted[-1] - meanSorted[0]
    
    if len(scores)==0:
        print("No dataset found with sufficient scores")
        return
    else:
        print("Datasets with scores:", len(scores))

    # Save results to file
    result = {'geneIdsFromSymbol':geneIdsFromSymbol, 'scores':scores.sort_values('var', ascending=False).reset_index().to_dict(orient='records')}
    with open(os.path.join(os.environ['GENESET_FILEPATH'], '%s.json' % name), 'w') as f:
        f.write(json.dumps(result))

def allGenesets():
    return [item.split('.')[0] for item in os.listdir(os.environ['GENESET_FILEPATH']) if item.endswith('.json')]

def genesetFromName(name):
    filepath = os.path.join(os.environ['GENESET_FILEPATH'],"%s.json" % name)
    if os.path.exists(filepath):
        with open(filepath) as f:
            return json.loads(f.read())

def createGenesets():
    createGeneset(msigName='NABA_COLLAGENS')
    createGeneset(msigName='REACTOME_COLLAGEN_FORMATION')
    return
    df = pandas.read_csv(os.path.join(os.environ['EXPRESSION_FILEPATH'],'../received/zahra/VillaniGenesets.tsv'), sep='\t')
    for column in df.columns:
        name = df.columns[2]
        createGeneset(name=name, geneSymbols=df[name].dropna().tolist())

# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------
def test_sampleGroupToGenes():
    print(sampleGroupToGenes('cell_type','fibroblast'))
    #['rankScore'].head())

def test_geneToSampleGroups():
    print(geneToSampleGroups('ENSG00000102145'))

def test_geneset():
    gs = genesetFromName('NABA_COLLAGENS')
    df = pandas.DataFrame.from_records(gs['scores']).set_index('index')
    print(sorted(df['diff']))
    return
    print(df['high'].value_counts())
    print(df['low'].value_counts())

# def geneset(geneIds=[], searchString="", limit=100):
#     """Return a pandas DataFrame which contains attributes of genes matching the parameters.
#     """
#     geneIds = list(set([item for item in geneIds if item.startswith("ENS")]))

#     # Don't include _id field, since this its value is ObjectId class which can't be
#     # serialised by json, which can create issues downstream, and this field is not needed outside mongo anyway.
#     if len(geneIds)>0:
#         params = {"gene_id": {"$in": geneIds}}
#     else:
#         params = {'gene_name': {'$regex': searchString, '$options': 'i'}}

#     if limit:
#         cursor = database["genes"].find(params, {"_id":0}).limit(limit)
#     else:
#         cursor = database["genes"].find(params, {"_id":0})

#     return pandas.DataFrame(cursor).set_index("gene_id") if cursor.count()!=0 else pandas.DataFrame()

