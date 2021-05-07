"""
2021-04-13. Jarny.

To run this script (ensure you're in the application directory):
(s4m-api) [ec2-user@api-dev s4m-api]$ python -m scripts.create_gene_ids_for_microarray -d 1000

where 1000 is the dataset id example.
This script creates .genes version of the expression matrix for a microarray dataset,
by using maximum value for all probes that map to the same gene.
It uses probe id to gene id mapping file which comes from psql dump of the stemformatics dbs,
more specifically stemformatics.feature_mappings table.
"""

import os, sys, pandas, argparse

# append current directoy to path so modules can be imported from models
sys.path.append(os.path.join(sys.path[0]))
from models import datasets, probes

def createFile(datasetId, report_only=False):
    ds = datasets.Dataset(datasetId)
    if ds.metadata()['platform_type']!='Microarray':
        print("Not microarray dataset")
        return

    # Read raw expression matrix - some probe ids look like integers, but we'll convert to string
    df = ds.expressionMatrix()
    df.index = df.index.astype(str)

    # Read gene id probe id mapping
    mapping = probes.probeMappingMatrix(datasetId)
    if mapping is None:
        print("Can't find probe mapping file for this dataset")
        return
    
    if report_only:
        print(mapping.head())
        print(df.shape, len(set(df.index)), mapping.shape)
        print(len(set(df.index).intersection(set(mapping['probeId']))))
        mapping = mapping[mapping['probeId'].isin(df.index)]
        print(mapping.shape, len(set(mapping['geneId'])))
        return
        #     probeId           geneId
        # 0  17112607  ENSG00000000003
        # 1  17105321  ENSG00000000005
        # 2  16920287  ENSG00000000419
        # 3  16696276  ENSG00000000457
        # 4  16673557  ENSG00000000460
        # (53617, 10) 53617 (60366, 2)
        # 39697
        # (60366, 2) 40820

    # Subset mapping using probes only found in df
    # then create a dictionary of probe ids keyed on gene id
    probeIdsFromGeneIds = {}
    mapping = mapping[mapping['probeId'].isin(df.index)]
    for index,row in mapping.iterrows():
        probeId = row['probeId']
        geneId = row['geneId']
        if geneId not in probeIdsFromGeneIds:
            probeIdsFromGeneIds[geneId] = set()
        probeIdsFromGeneIds[geneId].add(probeId)
    
    # Each key of the dictionary will form a row of new dataframe
    # with all matching probe rows aggregated using max.
    genes = []
    data = []
    for geneId,probeIds in probeIdsFromGeneIds.items():
        data.append(df.loc[list(probeIds)].max().tolist())
        genes.append(geneId)

    if len(data)<0.5*len(df):
        print("Less than 50% of probes recovered through mapping. Stopping.")
        return

    # write to file
    print(len(df), len(data))
    pandas.DataFrame(data, index=genes, columns=df.columns).to_csv(ds.expressionFilePath(key='genes'), sep='\t')

def createFileOld(datasetId, report_only=False):
    ds = datasets.Dataset(datasetId)
    if ds.metadata()['platform_type']!='Microarray':
        print("Not microarray dataset")
        return

    # Read raw expression matrix
    df = ds.expressionMatrix()

    # Read gene id probe id mapping
    mapping = pandas.read_csv('/mnt/stemformatics-data/received/www1/feature_mappings.tsv', sep='\t', header=None)
    mapping.columns = ['db_id','mapping_id','from_type','from_id','to_type','to_id']
    
    if report_only:
        print(mapping.head())
        print(df.shape, len(set(df.index)), mapping.shape)
        print(len(set(df.index).intersection(set(mapping['from_id']))))
        print(len(set(df.index).intersection(set(mapping['to_id']))))
        mapping = mapping[mapping['to_id'].isin(df.index)]
        print(mapping.shape, len(set(mapping['from_id'])))
        print(mapping.groupby('db_id').size())
        return
    #
    #     db_id  mapping_id from_type             from_id to_type                      to_id
    # 0     56         121      Gene    ENSG000001154151   Probe            ENSG00000115415
    # 1     46         113      Gene  ENSMUSG00000056912   Probe  chr10:100056156-100056207
    # 2     46         113      Gene  ENSMUSG00000044921   Probe  chr10:101974535-101974639
    # 3     46         113      Gene  ENSMUSG00000019892   Probe  chr10:102698896-102699016
    # 4     46         113      Gene  ENSMUSG00000019894   Probe  chr10:102830118-102830211
    # (46693, 50) 46693 (4136182, 6)
    # 0
    # 28579
    # (194758, 6) 28108
    # db_id
    # 59    194758

    # Subset mapping using probes only found in df
    # then create a dictionary of probe ids keyed on gene id
    probeIdsFromGeneIds = {}
    mapping = mapping[mapping['to_id'].isin(df.index)]
    for index,row in mapping.iterrows():
        probeId = row['to_id']
        geneId = row['from_id']
        if geneId not in probeIdsFromGeneIds:
            probeIdsFromGeneIds[geneId] = set()
        probeIdsFromGeneIds[geneId].add(probeId)
    
    # Each key of the dictionary will form a row of new dataframe
    # with all matching probe rows aggregated using max.
    genes = []
    data = []
    for geneId,probeIds in probeIdsFromGeneIds.items():
        data.append(df.loc[list(probeIds)].max().tolist())
        genes.append(geneId)

    # write to file
    pandas.DataFrame(data, index=genes, columns=df.columns).to_csv(ds.expressionFilePath(key='genes'), sep='\t')

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", help="dataset id")
    parser.add_argument("-r", help="show report only", default=False, action="store_true")
    args = parser.parse_args()
    createFile(args.d, report_only=args.r)
