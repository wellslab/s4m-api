"""Main interface for microarray probe annotation and mapping data.

The files used here come from the current/previous Stemformatics system, and you can see the list
of files online at https://www.stemformatics.org/contents/download_mappings, where the URLs
point to /var/www/html/mappings/, on the server (www1 for example).

Each file has the format of mapping_xx.txt, where xx is the chip_type, as found in the
assay_platforms table in the psql db:

[jarny@w1-s4m ~]$ psql -U portaladmin portal_prod

portal_prod=# select * from assay_platforms where platform like 'HG-U133%';
 chip_type | species | manufacturer |    platform     | version | platform_type | min_y_axis |  y_axis_label   | y_axis_
label_description | mapping_id | log_2 | probe_name |  default_graph_title
-----------+---------+--------------+-----------------+---------+---------------+------------+-----------------+--------
------------------+------------+-------+------------+-----------------------
        16 | Human   | Affymetrix   | HG-U133_2       |         | microarray    |          0 | Log2 Expression |
                  |         16 | t     | Probe      | Gene Expression Graph
        51 | Human   | Affymetrix   | HG-U133A        |         | microarray    |          0 | Log2 Expression |
                  |         51 | t     | Probe      | Gene Expression Graph
        56 | Human   | Affymetrix   | HG-U133A        | V2      | microarray    |          0 | Log2 Expression |
                  |         56 | t     | Probe      | Gene Expression Graph
       136 | Human   | Affymetrix   | HG-U133B        |         | microarray    |          0 | Log2 Expression |
                  |        136 | t     | Probe      | Gene Expression Graph
       142 | Human   | Affymetrix   | HG-U133 A+B Set |         | microarray    |          0 | Log2 Expression |
                  |        142 | t     | Probe      | Gene Expression Graph

Hence the combination of this table and all the .txt files will enable us to map
any probe id to gene id in the system. So after dumping the contents of the sql table into a tsv file:
copy (select * from assay_platforms) to '/tmp/assay_platforms.tsv' with (format csv, header, delimiter E'\t');

These were copied to PROBE_MAPPING_FILEPATH in the new system.
"""
import os, pandas

def createMappingFile():
    """The sql dump file assay_platforms.tsv has extra columns we don't need, so trim these and save
    a different version that will be used in this model.
    (Used nosetests -s models/probes.py:createMappingFile to run this one time)
    """
    os.chdir(os.environ['PROBE_MAPPING_FILEPATH'])
    df = pandas.read_csv('assay_platforms.tsv', sep='\t', index_col=0)

    # Check the match between chip_type of this table and .txt file existing
    ids = [int(item.split("_")[1].replace('.txt','')) for item in os.listdir() if item.startswith('mapping')]
    print(set(ids).issubset(set(df.index)), set(df.index).issubset(set(ids)))
    # prints True, False, so finding chip_type in the table does not guarantee that file exists, though inverse is true

    # Only keep these columns and write to file
    df = df[['species','manufacturer','platform','version','platform_type']]
    print(df.head())
    df.to_csv('assay_platforms_cleaned.tsv', sep='\t')

def probeMappingFilepath(**kwargs):
    """Return the full path to the probe mapping file corresponding to the query.
    """
    datasetId = kwargs.get('datasetId')
    topDir = os.environ['PROBE_MAPPING_FILEPATH']

    os.chdir(topDir)
    df = pandas.read_csv('assay_platforms_cleaned.tsv', sep='\t', index_col=0)
    subset = pandas.DataFrame()  # will be assigned as the matching subset of df

    def _filepathFromSubset(subset):
        if len(subset)==1:  # exact match
            filepath = 'mapping_%s.txt' % subset.index[0]
            if os.path.exists(filepath):
                return os.path.join(topDir, filepath)
        return None

    if datasetId:  # try to work out platform based on what's in the dataset metadata
        from models import datasets
        platform = datasets.Dataset(datasetId).metadata()['platform']
        # First we just try some pre-defined matches
        platformRep = {'Affymetrix Human Genome U133 Plus 2.0 Array [HG-U133_Plus_2]':{'manufacturer':'Affymetrix', 'platform':'HG-U133_2'},
                       'Illumina HumanHT-12 v4.0 Expression BeadChip':{'manufacturer':'Illumina', 'platform':'HumanHT-12', 'version':'V4'},
                       'Affymetrix Human Exon 1.0 ST Array [transcript (gene) version] [HuEx-1_0-st]':{'manufacturer':'Affymetrix', 'platform':'HuEx-1_0-ST', 'version':'V2'},
                       'Illumina MouseRef-8 v2.0 Expression BeadChip':{'manufacturer':'Illumina', 'platform':'MouseRef-8', 'version':'V2'},
                       'Agilent-014850 Whole Human Genome Microarray 4x44K G4112F':{'manufacturer':'Agilent', 'platform':'4x44 014850 G4112F'},
                       'Illumina Human-6 v2.0 Expression BeadChip':{'manufacturer':'Illumina', 'platform':'HumanWG-6', 'version':'V2'},
                       'Illumina HumanRef-8 v2.0 Expression BeadChip':{'manufacturer':'Illumina', 'platform':'HumanRef-8', 'version':'V2'},
                       'Affymetrix HT Human Genome U133A Array [HT_HG-U133A]':{'manufacturer':'Affymetrix', 'platform':'HT-HG-U133A',},
                      }
        #print(platform)
        if platform in platformRep:
            subset = df[(df['manufacturer']==platformRep[platform]['manufacturer']) & (df['platform']==platformRep[platform]['platform'])]
            if 'version' in platformRep[platform]:
                subset = subset[subset['version']==platformRep[platform]['version']]
            return _filepathFromSubset(subset)

        platform = platform.lower()  # lowercase version for easier pattern matching
        # This is usually a combination of several fields
        values = platform.split(' ')
        manufacturer = values[0]
        platform = ' '.join(values[1:])

    # 
    for col in df.columns:  # use all lowercase for better matching
        df[col] = df[col].str.lower()
    subset = df[(df['manufacturer']==manufacturer) & (df['platform']==platform)]
    if len(subset)==0:  # try to match just platform?
        subset = df[df['platform']==platform]
    if len(subset)==0:  # some platforms have square brackets eg. "human gene 2.0 st array [hugene-2_0-st]"
        import re
        search = re.search(r"\[(.*)\]", platform)  # match stuff in the square brackets
        if search:
            subset = df[df['platform']==search.group(1)]
    if len(subset)==0:  # try just first word in platform
        subset = df[(df['manufacturer']==manufacturer) & (df['platform']==platform.split(' ')[0])]

    return _filepathFromSubset(subset)

def probeMappingMatrix(datasetId):
    filepath = probeMappingFilepath(datasetId=datasetId)
    if not filepath: return None

    df = pandas.read_csv(filepath, sep='\t', header=None, names=['probeId','geneId'])
    # some probe ids look like integers
    df['probeId'] = df['probeId'].astype(str)
    return df

def test_probeMapping():
    df = probeMappingMatrix(6051)
    print(df.head() if df is not None else "mapping file not found")