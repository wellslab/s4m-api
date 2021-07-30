"""
Gene expression analysis related resources here.
"""
from flask_restful import reqparse, Resource
import os, pandas
from models import genes, datasets
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
        """Returns highly expressed genes in sampleGroup = sampleGroupItem. Returned object is a pandas DataFrame
        which contains scores and dataset ids.
        Bulk of this code should be in the model, but it's easier to filter out dataset we don't want to look at
        by using DatasetSearch class from resources for now. Once these filters don't need to be applied, we can move
        this code into models.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('sample_group', type=str, required=True)
        parser.add_argument('sample_group_item', type=str, required=True)
        parser.add_argument('cutoff', type=int, required=False, default=10)  # use 0 or None to get all results
        parser.add_argument('orient', type=str, required=False, default="records")
        args = parser.parse_args()
        
        allDatasets = DatasetSearch().get()
        sampleGroup = args.get("sample_group")
        sampleGroupItem = args.get('sample_group_item')

        # This will hold genes as index and rank score for each gene in each dataset
        rankScore = pandas.Series(dtype=float)
        datasetIds = {}
        uniqueDatasetIds = set()  # keep track of all dataset ids used for scoring

        for datasetDict in allDatasets:
            if sampleGroupItem not in datasetDict[sampleGroup] or len(datasetDict[sampleGroup].split(","))==1:
                continue
            
            ds = datasets.Dataset(datasetDict['dataset_id'])
            samples = ds.samples()
            groupItems = [item if item==sampleGroupItem else 'other' for item in samples[sampleGroup]]
            if len(set(groupItems))==1:  # some datasets may have one sampleGroupItem specified only
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
            return {'sampleGroup':sampleGroup, 'sampleGroupItem':sampleGroupItem, 'rankScore':pandas.DataFrame().to_dict(orient=args.get('orient')),
                    'totalDatasets':len(uniqueDatasetIds)}

        # Add average rank values
        df['meanRank'] = [rankScore.loc[geneId].mean() for geneId in df.index]
        
        # Apply cutoff - this is applied to each combination of geneId-count
        if args.get('cutoff'):
            df = pandas.concat([df[df['count']==count].sort_values('meanRank', ascending=False).iloc[:args.get('cutoff'),:] for count in df['count'].unique()])

        df = df.sort_values(['count','meanRank'], ascending=False)
        
        # Add datasetIds
        df['datasetIds'] = [','.join(map(str,datasetIds[geneId])) for geneId in df.index]
        df.index.name = "geneId"
        
        # # Add gene symbols
        # geneSymbols = geneSymbolsFromGeneIds(df.index.tolist())
        # df['geneSymbol'] = [geneSymbols.get(geneId,geneId) for geneId in df.index]

        if args.get('orient')=='records':
            df = df.reset_index()
        return {'sampleGroup':sampleGroup, 'sampleGroupItem':sampleGroupItem, 'rankScore':df.to_dict(orient=args.get('orient')), 'totalDatasets':len(uniqueDatasetIds)}

class GeneToSampleGroups(Resource):
    def get(self):
        """Returns highly expressed genes in sampleGroup = sampleGroupItem. Returned object is a pandas DataFrame
        which contains scores and dataset ids.
        Bulk of this code should be in the model, but it's easier to filter out dataset we don't want to look at
        by using DatasetSearch class from resources for now. Once these filters don't need to be applied, we can move
        this code into models.
        """
        parser = reqparse.RequestParser()
        parser.add_argument('gene_id', type=str, required=True)
        args = parser.parse_args()
        
        allDatasets = DatasetSearch().get()
        geneId = args.get('gene_id')
        sampleGroup = 'cell_type'  # perhaps turn this into a parameter later
        result = pandas.DataFrame(columns=[sampleGroup, 'datasetIds', 'var', 'ranks']).set_index(sampleGroup)  # result to fill in
        for datasetDict in allDatasets:
            ds = datasets.Dataset(datasetDict['dataset_id'])
            exp = ds.expressionMatrix(key='cpm', applyLog2=True)
            if geneId not in exp.index:
                continue

            samples = ds.samples()
            samples = samples.loc[samples.index.intersection(exp.columns)]
            if len(samples)<4: # we need at least 2 replicates in each cell type
                continue
            elif len(samples[sampleGroup].unique())<2:  # we need at least 2 different cell types
                continue

            # Find mean of gene expression across cell types
            values = exp.loc[geneId, samples.index]
            mean = values.groupby(samples[sampleGroup]).mean().round()
            rank = mean.rank()
            for celltype,value in rank.items():
                if celltype not in result:
                    result.loc[celltype] = [[], 0, []]
                result.loc[celltype,'datasetIds'].append(ds.datasetId)
                result.loc[celltype,'ranks'].append(value)
                result.loc[celltype,'var'] = mean.var()

            if len(result)>30: break
        result.to_csv('%s.tsv' % geneId, sep='\t')
        return result.reset_index().to_dict(orient='records')          

class GenesetCollection(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('name', type=str, required=False)
        parser.add_argument('dataset_id', type=int, required=False)
        parser.add_argument('sample_group', type=int, required=False, default='cell_type')
        args = parser.parse_args()

        if args.get('name') is None: # return a list of current geneset names
            return genes.allGenesets()
        elif args.get('dataset_id') is None: # return dataset ids for named geneset
            return genes.genesetFromName(args.get('name'))
        else:   # return details about particular geneset and dataset combo
            geneset = genes.genesetFromName(args.get('name'))
            #geneIds = list(geneset['geneSymbolsFromId'].keys())
            ds = protectedDataset(args.get('dataset_id'))
            exp = ds.expressionMatrix(key='cpm', applyLog2=True)

            # Now we subset exp on geneset. There may be multiple gene symbols per gene id, 
            # but we just choose the one with most variance
            selectedGeneIds, selectedGeneSymbols = [], []
            for geneSymbol,geneIds in geneset['geneIdsFromSymbol'].items():
                df = exp.loc[exp.index.intersection(geneIds)]
                if len(df)>0:
                    selectedGeneIds.append(df.var(axis=1).sort_values().index[-1])
                    selectedGeneSymbols.append(geneSymbol)

            df = exp.loc[selectedGeneIds]
            samples = ds.samples().fillna('')

            from scipy.stats import zscore
            zscores = pandas.DataFrame([zscore(df.loc[rowId]) for rowId in df.index], index=df.index, columns=df.columns).fillna(0)

            # cluster rows and columns based on zscore
            from scipy.spatial.distance import pdist, squareform
            import scipy.cluster.hierarchy as hc

            rowDist = squareform(pdist(zscores.to_numpy()))
            rowOrdering = hc.leaves_list(hc.linkage(rowDist, method='centroid'))

            colDist = squareform(pdist(zscores.transpose().to_numpy()))
            colOrdering = hc.leaves_list(hc.linkage(colDist, method='centroid'))

            zscores = zscores.iloc[rowOrdering, colOrdering]

            # Relabel columns to 'xxx_1','xxx_2', etc, where xxx is the sample group item
            sampleGroupItems = samples.loc[zscores.columns][args.get('sample_group')]
            uniqueItems = dict([(item,[]) for item in sampleGroupItems.unique()])
            columns = []
            for item in sampleGroupItems:
                columns.append("%s_%s" % (item, len(uniqueItems[item])))
                uniqueItems[item].append(item)
            zscores.columns = columns
            zscores.index = selectedGeneSymbols
            return zscores.to_dict(orient='split')


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

