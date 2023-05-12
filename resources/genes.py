"""
Gene expression analysis related resources here.
"""
from flask_restful import reqparse, Resource
import os, pandas
from models import genes, atlases
from resources.datasets import protectedDataset, DatasetSearch

# def geneSymbolsFromGeneIds(geneIds):
#     headers = {'content-type': 'application/x-www-form-urlencoded'}
#     params = 'ids=%s&fields=symbol,ensembl.gene' % ','.join(geneIds)
#     r = requests.post('http://mygene.info/v3/gene', data=params, headers=headers)
#     result = {}
#     for item in r.json():
#         result[item['query']] = item['symbol'] if 'symbol' in item else item['query']
#     return result

class SampleGroupToGenes(Resource):
    def get(self):
        """Returns highly expressed genes in sampleGroup = sampleGroupItem. Returned object is a dictionary
        containing 'rankScore' key, which is a pandas DataFrame that contains genes and dataset ids and their scores.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('sample_group', type=str, required=True)
        parser.add_argument('sample_group_item', type=str, required=True)
        parser.add_argument('sample_group_item2', type=str, required=False)
        parser.add_argument('cutoff', type=int, required=False, default=20)  # use 0 or None to get all results
        parser.add_argument('scoring_method', type=str, required=False, default='max')
        args = parser.parse_args()
        
        result = genes.sampleGroupToGenes(args.get('sample_group'), args.get('sample_group_item'), 
            sampleGroupItem2=args.get('sample_group_item2'), cutoff=args.get('cutoff'), scoringMethod=args.get('scoring_method'))
        return {'sampleGroup':args.get('sample_group'), 'sampleGroupItem':args.get('sample_group_item'), 'sampleGroupItem2':args.get('sample_group_item2'),
                'rankScore':result['rankScore'].reset_index().to_dict(orient='records'), 'totalDatasets':result['totalDatasets']}

class GeneToSampleGroups(Resource):
    def get(self):
        """Returns highly expressed sample group items under sample group, given gene_id.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('gene_id', type=str, required=True)
        parser.add_argument('sample_group', type=str, required=False, default='cell_type')
        args = parser.parse_args()
                
        result = genes.geneToSampleGroups(args.get('gene_id'), args.get('sample_group'))
        return result.to_dict(orient='index')

"""Return index of df after hierarchical clustering the rows"""
def hclusteredRows(df):
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist, squareform
    clust = linkage(squareform(pdist(df.values, metric='euclidean')), method='complete')
    dendro = dendrogram(clust, no_plot=True)
    return df.index[dendro['leaves']]

class GenesetTable(Resource):
    def get(self):
        #from scipy.stats import zscore
        parser = reqparse.RequestParser()
        parser.add_argument('geneset_group', type=str, required=False, default='DE genes') # ['DE genes','Hallmark','WGCNA']
        parser.add_argument('geneset', type=str, required=False)
        parser.add_argument('timepoint', type=str, required=False, default='2h')
        parser.add_argument('cluster_columns', type=str, required=False, default='false')
        args = parser.parse_args()

        gs = genes.genesetTable(genesetGroup=args.get('geneset_group')).fillna('')
        genesetName = args.get('geneset') if args.get('geneset') else gs.index[0]
        timepoint = args.get('timepoint')
        genelist = gs.at[genesetName, 'genes'].split(',')
        clusterColumns = args.get('cluster_columns').lower().startswith('t')

        # Read sample table and subset on timepoint
        atl = atlases.Atlas('ma')
        samples = atl.sampleMatrix()
        #samples = samples[(samples['treatment']==treatment) | (samples['treatment']=='NS')]
        samples = samples[samples['time']==timepoint]

        # Read expression matrix and subset on genes and samples
        df = atl.expressionMatrix().loc[genelist, samples.index]

        # group by treatment (at this timepoint) and work out the mean
        df = df[samples.index].groupby(samples['treatment'], axis=1).mean()#.apply(zscore, axis=1)
        df = df.sub(df['NS'], axis=0)
        df = df.loc[hclusteredRows(df)]

        # Substitute gene symbols
        geneSymbolFromId = atl.geneInfo().loc[genelist, 'symbol'].fillna('').to_dict()
        geneSymbols = [geneSymbolFromId.get(geneId, geneId) for geneId in df.index]

        # sort columns, either by treatment type or by hcluster. NS is always first
        orderedColumns = ['NS']
        if clusterColumns:
            for item in hclusteredRows(df.drop(columns=['NS']).transpose()):
                orderedColumns.append(item)
        else:
            for treatmentType in samples['type'].unique():
                for item in samples[samples['type']==treatmentType]['treatment'].unique():
                    if item!='NS': orderedColumns.append(item)
        df = df[orderedColumns]

        #print(samples.head())
        #print(df.shape)
        return {'genesets':gs[['description']].reset_index().to_dict(orient='records'), 
                'expression':df.to_dict(orient='split'),
                'geneSymbols':geneSymbols}

        # Change the format of the data
        df = pandas.DataFrame(columns=['2h','6h','16h'])
        df.index.name = 'name'

        treatments = list(set(['_'.join(item.split('_')[:-1]) for item in gs.index]))
        for treatment in treatments:
            for col in df.columns:
                genelist = gs.at[f"{treatment}_{col}",'genes']
                df.at[treatment, col] = genelist.split(',') if len(genelist)>0 else []

        return df.reset_index().to_dict(orient='records')
    
# Not used for now
# class GenesetCollection(Resource):
#     def get(self):
#         parser = reqparse.RequestParser()
#         parser.add_argument('name', type=str, required=False)
#         parser.add_argument('dataset_id', type=int, required=False)
#         parser.add_argument('sample_group', type=int, required=False, default='cell_type')
#         args = parser.parse_args()

#         if args.get('name') is None: # return a list of current geneset names
#             return genes.allGenesets()
#         elif args.get('dataset_id') is None: # return dataset ids for named geneset
#             return genes.genesetFromName(args.get('name'))
#         else:   # return details about particular geneset and dataset combo
#             geneset = genes.genesetFromName(args.get('name'))
#             #geneIds = list(geneset['geneSymbolsFromId'].keys())
#             ds = protectedDataset(args.get('dataset_id'))
#             exp = ds.expressionMatrix(key='cpm', applyLog2=True)

#             # Now we subset exp on geneset. There may be multiple gene symbols per gene id, 
#             # but we just choose the one with most variance
#             selectedGeneIds, selectedGeneSymbols = [], []
#             for geneSymbol,geneIds in geneset['geneIdsFromSymbol'].items():
#                 df = exp.loc[exp.index.intersection(geneIds)]
#                 if len(df)>0:
#                     selectedGeneIds.append(df.var(axis=1).sort_values().index[-1])
#                     selectedGeneSymbols.append(geneSymbol)

#             df = exp.loc[selectedGeneIds]
#             samples = ds.samples().fillna('')

#             from scipy.stats import zscore
#             zscores = pandas.DataFrame([zscore(df.loc[rowId]) for rowId in df.index], index=df.index, columns=df.columns).fillna(0)

#             # cluster rows and columns based on zscore
#             from scipy.spatial.distance import pdist, squareform
#             import scipy.cluster.hierarchy as hc

#             rowDist = squareform(pdist(zscores.to_numpy()))
#             rowOrdering = hc.leaves_list(hc.linkage(rowDist, method='centroid'))

#             colDist = squareform(pdist(zscores.transpose().to_numpy()))
#             colOrdering = hc.leaves_list(hc.linkage(colDist, method='centroid'))

#             zscores = zscores.iloc[rowOrdering, colOrdering]

#             # Relabel columns to 'xxx_1','xxx_2', etc, where xxx is the sample group item
#             sampleGroupItems = samples.loc[zscores.columns][args.get('sample_group')]
#             uniqueItems = dict([(item,[]) for item in sampleGroupItems.unique()])
#             columns = []
#             for item in sampleGroupItems:
#                 columns.append("%s_%s" % (item, len(uniqueItems[item])))
#                 uniqueItems[item].append(item)
#             zscores.columns = columns
#             zscores.index = selectedGeneSymbols
#             return zscores.to_dict(orient='split')


# class Geneset(Resource):
#     def get(self):
#         """Return attibutes of genes from query parameters. The returned dictionary looks like
#             [
#                 {
#                     "gene_id": "ENSG00000133055",
#                     "gene_name": "MYBPH",
#                     "gene_description": "myosin binding protein H [Source:HGNC Symbol;Acc:HGNC:7552]"
#                 },
#                 ...
#             ]
#         """
#         parser = reqparse.RequestParser()
#         parser.add_argument('gene_id', type=str, required=False, default="", action="append")  # list of strings
#         parser.add_argument('search_string', type=str, required=False, default="")
#         parser.add_argument('limit', type=int, required=False, default=50, action="append")
#         parser.add_argument('orient', type=str, required=False, default="records")
#         args = parser.parse_args()
        
#         df = genes.geneset(geneIds=args.get("gene_id"), searchString=args.get("search_string"), limit=args.get("limit"))
#         if args.get('orient')=='records':
#             df = df.reset_index()
#         return df.fillna("").to_dict(orient=args.get("orient"))

