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
import requests, os, pandas, json
from models import datasets

def createGenesetFromMsigDB(name, sampleGroupItem='cell_type'):
    """Create a geneset from name of MsigDB from Broad Institute.
    """
    if "GENESET_FILEPATH" not in os.environ:
        print("GENESET_FILEPATH not found in os.environ.")
        return

    # First use the api from msigdb to get all gene symbols from the name
    url = 'https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=%s&fileType=txt' % name
    r = requests.get(url)

    # Now use mygene.info to get ensembl ids
    params = {'q':','.join(r.text.split('\n')[2:]),'scopes':'symbol', 'species':'human', 'fields':'symbol,ensembl.gene', 'ensemblonly':True}
    res = requests.post('http://mygene.info/v3/query', params)
    
    # Create dict of gene symbols keyed on ids. Note that some ensembl entries are multiple matches
    #geneSymbolsFromId = {}
    geneIdsFromSymbol = {}
    for item in res.json():
        if isinstance(item['ensembl'], list):
            geneIdsFromSymbol[item['symbol']] = [val['gene'] for val in item['ensembl']]
            # for val in item['ensembl']:
            #     if val['gene'] not in geneSymbolsFromId: geneSymbolsFromId[val['gene']] = []
            #     geneSymbolsFromId[val['gene']].append(item['symbol'])
        else:
            geneIdsFromSymbol[item['symbol']] = [item['ensembl']['gene']]
            # if item['ensembl']['gene'] not in geneSymbolsFromId: geneSymbolsFromId[item['ensembl']['gene']] = []
            # geneSymbolsFromId[item['ensembl']['gene']].append(item['symbol'])

    # Loop through each dataset in the system
    vars = {}
    allDatasets = datasets.datasetMetadataFromQuery(ids_only=True, organism='homo sapiens')
    for datasetId in allDatasets:
        ds = datasets.Dataset(datasetId)
        if not (ds.platformType()=='RNASeq' or ds.platformType()=='Microarray'): continue
        samples = ds.samples()
        if len(samples[sampleGroupItem].unique())==1:  # some datasets may have one sampleGroupItem specified only
            continue
        
        exp = ds.expressionMatrix(key='cpm', applyLog2=True)
        # Take mean of geneset and mean across sampleGroupItem, and evaluate variance of these means
        #df = exp.loc[exp.index.intersection(list(geneSymbolsFromId.keys()))]
        df = exp.loc[exp.index.intersection(sum(list(geneIdsFromSymbol.values()),[]))]
        if len(samples.index.intersection(df.columns))==0: continue  # shouldn't happen but will for some datasets
        if len(df)>0:
            variance = df[samples.index].groupby(samples[sampleGroupItem], axis=1).mean().mean().var()
            if pandas.notnull(variance):
                vars[datasetId] = variance
                
    # Save results to file
    # result = {'geneSymbolsFromId':geneSymbolsFromId, 'vars':pandas.Series(vars).sort_values(ascending=False).to_dict()}
    result = {'geneIdsFromSymbol':geneIdsFromSymbol, 'vars':pandas.Series(vars).sort_values(ascending=False).to_dict()}
    with open(os.path.join(os.environ['GENESET_FILEPATH'], '%s.json' % name), 'w') as f:
        f.write(json.dumps(result))

def allGenesets():
    return [item.split('.')[0] for item in os.listdir(os.environ['GENESET_FILEPATH']) if item.endswith('.json')]

def allGenesetsOld():
    genesets = []
    for item in os.listdir(os.environ['GENESET_FILEPATH']):
        if item.endswith('.json'):
            name = item.split('.')[0]
            vars = pandas.Series(genesetFromName(name)['vars'])
            vars = vars[vars>0.1].sort_values(ascending=False)
            genesets.append({'name':name, 'datasetIds':vars.index.tolist(), 'scores':vars.tolist()})
    return genesets

def genesetFromName(name):
    filepath = os.path.join(os.environ['GENESET_FILEPATH'],"%s.json" % name)
    if os.path.exists(filepath):
        with open(filepath) as f:
            return json.loads(f.read())

# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------

def test_createGenesetFromMsigDB():
    createGenesetFromMsigDB('REACTOME_COLLAGEN_FORMATION')

def test_genesetFromName():
    print(genesetFromName('NABA_COLLAGENS'))

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

