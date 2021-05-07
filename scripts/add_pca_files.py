"""
2021-02-18. Jarny.

To run this script (ensure you're in the application directory):
(s4m-api) [ec2-user@api-dev s4m-api]$ python -m scripts.add_pca_files

Create pca coordinates for all the expression files.
"""

import os, sys, pandas, numpy, shutil
from sklearn.decomposition import PCA

# append current directoy to path so modules can be imported from models
sys.path.append(os.path.join(sys.path[0]))
from models import datasets
from models import utilities

expressionFilepath = os.environ['EXPRESSION_FILEPATH']

# These datasets gave Segmentation fault (core dumped) error on api-dev server, apart from
# 7139,7077,6128,6645,6992,7204,7194,6309,6572,6599 which give ValueError: Input contains NaN, infinity or a value too large for dtype('float64')
datasetsToSkip = ['6121', '6128', '6190', '6256', '6271', '6278', '6329', '6353', '6400', 
                  '6401', '6418', '6425', '6459', '6461', '6468', '6491', '6518', '6580', 
                  '6586', '6602', '6639', '6645', '6646', '6667', '6668', '6728', '6730', 
                  '6735', '6739', '6741', '6748', '6991', '6992', '7032', '7077', '7128', 
                  '7129', '7135', '7139', '7142', '7168', '7169', '7170', '7192', '7194', 
                  '7200', '7204', '7239', '7243', '7253', '7268', '7327', '7346', '7378', 
                  '7379', '6214', '6231', '6286', '6309', '6385', '6460', '6572', '6599',
                  '6932', '6936', '7064', '7171', '7254', '7274', '7357', '7387']

def createPCAFiles(datasetId):
    """Given datasetId, create the pca files and place them where they should go.
    """
    os.chdir(expressionFilepath)
    df = pandas.read_csv(os.path.join(datasetId, "%s.raw.tsv" % datasetId), sep="\t", index_col=0)
    print(datasetId, df.shape, df.max().max())
    if df.max().max()>100: # log this first
        df = numpy.log2(df+1)

    # Auto select components using this. From docs (https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html):
    # If 0 < n_components < 1 and svd_solver == 'full', select the number of components such that the
    # amount of variance that needs to be explained is greater than the percentage specified by n_components.
    pca = PCA(n_components=0.99, svd_solver='full')
    coords = pca.fit_transform(df.values.T)
    pandas.DataFrame(coords, index=df.columns).to_csv(os.path.join(datasetId, "%s.pca.tsv" % datasetId), sep="\t")
    pandas.DataFrame({"explained_variance_":pca.explained_variance_,
                      "explained_variance_ratio_":pca.explained_variance_ratio_,
                      "singular_values_":pca.singular_values_}).to_csv(os.path.join(datasetId, "%s.pca_attributes.tsv" % datasetId), sep="\t")

def createPcaCoordinates():
    os.chdir(expressionFilepath)
    for dirname in sorted(os.listdir(expressionFilepath)):
        if os.path.isfile(dirname) or "_" in dirname or \
            os.path.exists(os.path.join(dirname, "%s.pca.tsv" % dirname)) or dirname in datasetsToSkip: 
        #if os.path.isfile(dirname) or "_" in dirname: 
            continue  # ignore files or actual directory (use symlink) or if pca file already exists
        createPCAFiles(dirname)

def copyDatasetsToSkip():
    """Copy the files for datasets to skip to gadi so that they can be processed there.
    We'll just copy the files to a temporary directory so that 
    """
    destination = "pca_files_todo"  # subdirectory of expressionFilepath
    os.chdir(expressionFilepath)
    for dirname in sorted(os.listdir(expressionFilepath)):
        if os.path.isfile(dirname) or "_" in dirname or dirname not in datasetsToSkip:
            continue  # ignore files or actual directory (use symlink) or if dataset not in datasetsToSkip
        shutil.copyfile(os.path.join(dirname, "%s.raw.tsv" % dirname), os.path.join(destination, "%s.raw.tsv" % dirname))

def movePCAFiles():
    """2021-03-02. Move the pca files which were done in gadi to where they're supposed to go.
    """
    os.chdir("/mnt/stemformatics-data/received/pca_files_todo")
    for filename in sorted(os.listdir()):
        datasetId = filename.split(".")[0]
        shutil.move(filename, os.path.join(expressionFilepath, "%s/" % datasetId))

def addMorePCAFiles():
    """2021-03-11. Create more pca files for some datasets we found.
    """
    datasetIds = ['6003', '6530']
    os.chdir(expressionFilepath)
    for datasetId in datasetIds:
        fileToUse = os.path.join(datasetId, "%s.raw.tsv" % datasetId) # apply PCA to this file
        df = pandas.read_csv(fileToUse, sep="\t", index_col=0)
        print(datasetId, df.shape, df.max().max())
        if df.max().max()>100: # log this first
            df = numpy.log2(df+1)

        pca = PCA(n_components=0.99, svd_solver='full')
        coords = pca.fit_transform(df.values.T)
        pandas.DataFrame(coords, index=df.columns).to_csv(os.path.join(datasetId, "%s.pca.tsv" % datasetId), sep="\t")
        pandas.DataFrame({"explained_variance_":pca.explained_variance_,
                          "explained_variance_ratio_":pca.explained_variance_ratio_,
                          "singular_values_":pca.singular_values_}).to_csv(os.path.join(datasetId, "%s.pca_attributes.tsv" % datasetId), sep="\t")

def fixPCAFiles():
    """2021-04-04. We first created PCA files using raw data for RNASeq, but we should have used cpm values. So re-do these.
    """
    os.chdir(expressionFilepath)
    # These are dc atlas datasets where cpm has already been used to generate the pca, so can ignore them.
    datasetsToIgnore = ['4135_1.0', '3082_1.0', '8144_1.0', '3120_1.0', '9747_1.0', '8507_1.0', '2494_1.0', '1611_1.0', '9002_1.0', '3559_1.0', '2865_1.0',
                        '6309_1.0']
    for dirname in sorted(os.listdir()):
        if os.path.isfile(dirname) or "_" not in dirname or dirname in datasetsToIgnore: continue  # ignore files or symlinks
        datasetId = dirname.split("_")[0]
        df = None
        try:
            ds = datasets.Dataset(datasetId)
            if ds.metadata()['platform_type']=='RNASeq':
                df = ds.expressionMatrix(key='cpm') # apply PCA to this file
        except datasets.DatasetIdNotFoundError:  # dataset not in metadata - we can still perform pca though
            df = pandas.read_csv(os.path.join(dirname, "%s.raw.tsv" % datasetId), sep="\t", index_col=0)
            if df.min().min()>0: df = None  # assume microarray if there are no zeros

        if df is not None:
            print("Working on", dirname, df.shape, df.max().max())
            if df.min().min()>=0:
                df = numpy.log2(df+1)

            # Auto select components using this. From docs (https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html):
            # If 0 < n_components < 1 and svd_solver == 'full', select the number of components such that the
            # amount of variance that needs to be explained is greater than the percentage specified by n_components.
            pca = PCA(n_components=0.99, svd_solver='full')
            coords = pca.fit_transform(df.values.T)
            pandas.DataFrame(coords, index=df.columns).to_csv(os.path.join(dirname, "%s.pca.tsv" % datasetId), sep="\t")
            pandas.DataFrame({"explained_variance_":pca.explained_variance_,
                              "explained_variance_ratio_":pca.explained_variance_ratio_,
                              "singular_values_":pca.singular_values_}).to_csv(os.path.join(dirname, "%s.pca_attributes.tsv" % datasetId), sep="\t")
            print("Done", dirname)

if __name__=="__main__":
    #createPcaCoordinates()
    #copyDatasetsToSkip()
    #movePCAFiles()
    #addMorePCAFiles()
    #fixPCAFiles()
    createPCAFiles(sys.argv[1])
    
