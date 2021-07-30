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
import requests, os, pandas, json, anndata
from models import datasets

def sampleGroupToGenes(sampleGroup, sampleGroupItem, public_only=True, cutoff=10):
    import time
    t0 = time.time()
    allDatasetIds = datasets.datasetMetadataFromQuery(ids_only=True, public_only=public_only)

    # This will hold genes as index and rank score for each gene in each dataset
    rankScore = pandas.Series(dtype=float)
    datasetIds = {}
    uniqueDatasetIds = set()  # keep track of all dataset ids used for scoring

    for datasetId in allDatasetIds:
        ds = datasets.Dataset(datasetId)
        samples = ds.samples()

        # ignore dataset if sampleGroupItem isn't found or there's only one sampleGroupItem
        if sampleGroupItem not in samples[sampleGroup] or len(samples[sampleGroup].unique())==1:
            continue
                
        exp = ds.expressionMatrix(key='cpm', applyLog2=True)
        exp = exp[samples.index]
        
        # Calculate the difference between mean of sampleGroupItem samples vs max of other in sampleGroup
        df1 = exp[samples[samples[sampleGroup]==sampleGroupItem].index].mean(axis=1)
        df2 = exp[samples[samples[sampleGroup]!=sampleGroupItem].index].max(axis=1)
        diff = df1 - df2

        # Only keep +ve scores
        diff = diff[diff>0]

        # Remember dataset id for all the genes in diff
        for geneId in diff.index:
            if geneId not in datasetIds: datasetIds[geneId] = []
            datasetIds[geneId].append(ds.datasetId)
            uniqueDatasetIds.add(ds.datasetId)

        # Row concatenate this Series for all datasets after converting diff into normalised ranks 
        # (so 0 is lowest diff value, 1 is highest)
        rankScore = pandas.concat([rankScore, diff.rank()/len(diff)])

    # Count the number of times each gene occurs in this series
    df = rankScore.index.value_counts().to_frame(name='count')
    if len(df)==0:
        return {'rankScore':pandas.DataFrame(), 'totalDatasets':len(uniqueDatasetIds)}

    # Add average rank values
    df['meanRank'] = [rankScore.loc[geneId].mean() for geneId in df.index]
    print(df.head())
    # Apply cutoff - this is applied to each combination of geneId-count
    if not 'cutoff' is None:
        df = pandas.concat([df[df['count']==count].sort_values('meanRank', ascending=False).iloc[:cutoff,:] for count in df['count'].unique()])

    df = df.sort_values(['count','meanRank'], ascending=False)
    
    # Add datasetIds
    df['datasetIds'] = [','.join(map(str,datasetIds[geneId])) for geneId in df.index]
    df.index.name = "geneId"
    
    print(time.time()-t0)
    return {'rankScore':df, 'totalDatasets':len(uniqueDatasetIds)}

def createGeneToSampleGroupsData():
    filepath = '/mnt/stemformatics-data/backups/gene_to_sample_groups.h5ad'
    ada = anndata.AnnData(obsm=pandas.DataFrame(columns=['cell_type', 'datasetId']))
    ada.write(filepath)
    return
    ada = anndata.read_h5ad(filepath)

    geneId = 'ENSG00000102145'
    sampleGroup = 'cell_type'  # perhaps turn this into a parameter later

    expressionFilepath = "/mnt/stemformatics-data/expression-files"
    os.chdir(expressionFilepath)
    for path in os.listdir(expressionFilepath):
        datasetId = int(path.split('_')[0])
        ds = datasets.Dataset(datasetId)
        ad = anndata.read_h5ad(path)

        samples = ds.samples()
        samples = samples.loc[samples.index.intersection(ad.obs_names)]
        if len(samples)<4: # we need at least 2 replicates in each cell type
            continue
        elif len(samples[sampleGroup].unique())<2:  # we need at least 2 different cell types
            continue
        samples = samples.fillna("[unspecified]")

        print(datasetId)
        if geneId in ad.var_names:
            series = ad[samples.index,geneId].to_df()[geneId]
            mean = series.groupby(samples[sampleGroup]).mean().round()
            df = pandas.DataFrame({'var':mean.var(), 'rank':mean.rank()})
            
            print(df)
        return

        if len(result)>30: break


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
    print(sampleGroupToGenes('cell_type','macrophage'))

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

