"""Model to handle atlas data in Stemformatics. Atlas data refers to an object of integrated data, 
and comprises of expression and sample matrices, as well as supporting data such as colours and genes.

colours.json
coordinates.tsv
expression.filtered.tsv
expression.tsv
genes.tsv
Readme.txt
samples.tsv

Examle usage:
--------------

# Note that we need to provide the path to the files which contains the atlas data.
> export ATLAS_FILEPATH=/path/to/atlas-files
atlas = Atlas('myeloid')
print(atlas.pcaCoordinates().head())
"""
import os, pandas, json, numpy, random

# ----------------------------------------------------------
# Functions
# ----------------------------------------------------------
def rankTransform(df):
    """Return a rank transformed version of data frame, where smallest value is 0 and largest is 1.
    For tied values, they are given same values, rather than one given higher value arbitrarily.

    For gene expression matrices, it's likely to be better to perform rank transform on the raw counts
    rather than log normalised values, because log values may differ very slightly when running the
    program and two genes may end up with different rank normalised values even though they may have
    started with same integer counts. Of course this doesn't apply to microarray data.

    The implementation here doesn't actually produce 0 as smallest value, but 1/df.shape[0].
    This is for historical reasons in our lab. When working with gene expression matrices, this lowest value
    will be small (large number of genes), and there will also be many tied values for RNA-seq data.
    Use this form to get guaranteed range of [0,1]:
        (df.rank(axis=0, ascending=True, na_option='bottom')-1)/(df.shape[0]-1) 
    """
    return (df.shape[0] - df.rank(axis=0, ascending=False, na_option='bottom')+1)/df.shape[0]

def hclusteredRows(df):
    """Return index of df after hierarchical clustering the rows.
    So use df = df.loc[hclusteredRows(df)] to change index after hclust.
    Seem to get "The symmetric non-negative hollow observation matrix looks suspiciously like an uncondensed distance matrix" warning
    for linkage function, even though squereform is called beforehand. So just ignoring this warning here.
    """
    from scipy.cluster.hierarchy import linkage, dendrogram, ClusterWarning
    from scipy.spatial.distance import pdist, squareform
    from warnings import simplefilter
    simplefilter("ignore", ClusterWarning)

    #clust = linkage(squareform(pdist(df.values, metric='euclidean')), method='complete')
    sf = squareform(pdist(df.values, metric='euclidean'))
    clust = linkage(sf, method='complete')
    dendro = dendrogram(clust, no_plot=True)
    return df.index[dendro['leaves']]

def atlasTypes():
    """Return a dictionary of available atlas types and versions.
        {'dc': {'versions': ['1.3', '1.2', '1.1'], 'current_version': '1.3', 'release_notes': 'Updated xxx...'}, ...}
    Note that 'versions' will be reverse sorted.
    """
    os.chdir(os.environ['ATLAS_FILEPATH'])
    dictToReturn = {}
    filelist = os.listdir()
    for atype in Atlas.all_atlas_types:
        dirlist = sorted([item for item in filelist if item.startswith('%s_' % atype)], reverse=True) # so we can ignore symlink
        dictToReturn[atype] = {'current_version': os.readlink(atype).split('_')[1],
                               'versions': [item.split('_')[1] for item in dirlist],
                               'release_notes': [open(os.path.join(item, 'Readme.txt')).read() for item in dirlist]}
    return dictToReturn

def projectSingleCellData(scData):
    """
    """
    return pandas.DataFrame()

# ----------------------------------------------------------
# Atlas class
# ----------------------------------------------------------
class Atlas(object):
    """Main class representing an atlas in Stemformatics.
    The atlas type specified at initialisation will be used to access the text files
    which reside under the correponding os.environ['ATLAS_FILEPATH'] directory.
    """

    # Full list of current atlas types
    all_atlas_types = ['myeloid','blood','dc','activation','ma']

    def __init__(self, atlasType, version=None):
        # each atlas type is under its own directory under ATLAS_FILEPATH
        directory = atlasType if not version else '%s_%s' % (atlasType, version)
        self.atlasFilePath = os.path.join(os.environ["ATLAS_FILEPATH"], directory)

        # Work out version based on link
        self.version = os.readlink(self.atlasFilePath).split("_")[1] if version is None else version

        self.atlasType = atlasType

    def pcaCoordinates(self):
        """Return a pandas DataFrame object, specifying the PCA coordinates. The data frame will have sample ids as index.
        Columns will be named as ['0','1',...] (as strings).
        """
        df = pandas.read_csv(os.path.join(self.atlasFilePath, "coordinates.tsv"), sep="\t", index_col=0)
        df.columns = [str(i) for i in range(len(df.columns))]
        return df

    def expressionFilePath(self, filtered=False):
        """Return the full file path to the expression matrix file on disk.
        """
        filename = "expression.filtered.tsv" if filtered else "expression.tsv"
        return os.path.join(self.atlasFilePath, filename)

    def expressionMatrix(self, filtered=False):
        """Return a pandas DataFrame of expression matrix after reading from file.
        If filtered=False, this is the "full" expression matrix that includes all genes. 
        The expression values are from rank normalised values.
        """
        return pandas.read_csv(self.expressionFilePath(filtered=filtered), sep="\t", index_col=0)

    def datasetIds(self):
        """Return all dataset ids in this atlas as a list. Note that each element will be integer type.
         Relies on sample ids in the form of "datasetId_sampleId", matching the format of sampleId in sample data.
        """
        return [int(item.split("_")[0]) for item in self.sampleMatrix().index]

    def sampleMatrix(self):
        """Return sample annotation matrix for the atlas as a pandas DataFrame object. The shape of the data frame
        will be number_of_samples x number_of_columns, with sample ids as index.
        NA values shouldn't happen in atlas samples, since they are more carefully annotated, but if they do exist,
        downstream analyses may be affected and it's probably best to assign them something.
        There are already 'unknown' values, but it may be a case that a particular column doesn't apply (eg. cell line)
        """
        return pandas.read_csv(os.path.join(self.atlasFilePath, "samples.tsv"), sep="\t", index_col=0).fillna('[NA]')

    def geneInfo(self):
        """Return a pandas DataFrame of information about all genes in the atlas, after reading the genes.tsv
        file in the atlas file directory. Ensembl ids form the index.
        """
        df = pandas.read_csv(os.path.join(self.atlasFilePath, "genes.tsv"), sep="\t", index_col=0)
        df.index.name = 'ensembl'
        return df

    def coloursAndOrdering(self):
        """Return dictionaries of colours and ordering of sample type items based on "colours.json" file inside the
        atlas file directory. Return empty dictionaries if such a file doesn't exist.
        Note that it's possible for colours and ordering to not exist for certain sample columns.
        """
        filepath = os.path.join(self.atlasFilePath, "colours.json")
        return json.loads(open(filepath).read()) if os.path.exists(filepath) else {}

    def projection(self, name, queryData, includeCombinedCoords=True, includeCapybara=True, includeRandomSamples=True):
        """Perform projection of queryData onto this atlas and return a dictionary of objects.
        Params:
            name (str): Name of the queryData, used as prefix for projected points if includeCombinedCoords is True.
            queryData (DataFrame): Expression matrix of query data to be projected onto this atlas.
            includeCombinedCoords (bool): Should atlas coords and projected coords be returned together, rather than just projections?
            includeRandomSamples (bool): If true, 3 random columns from queryData will have their row ids randomised, to provide "null" background. 

        Returns a dictionary with following keys and values:
            coords (DataFrame): projected coordinates of queryData.
            combinedCoords (DataFrame): data frame of atlas coordinates + projected coords if includeCombinedCoords is True.
                The projected points will have index in the format of "{name}_sample1", etc, so that these points
                can be distinguished from atlas points.
            error (str): Error message. Empty string if there was no error.
            name (str): Same as the value used as input.
        """
        result = {"error":"", "coords":pandas.DataFrame(), "name":name, "combinedCoords":pandas.DataFrame()}

        # Some validation before projecting
        if len(queryData)==0:
            result["error"] = "Data to project appears to have 0 rows. Check format of the file."

        # Read expression matrix - we only need filtered version. 
        df = self.expressionMatrix(filtered=True)
        genes = self.geneInfo()

        commonGenes = queryData.index.intersection(df.index)  # common index between test and atlas

        if len(queryData)!=len(set(queryData.index)):
            result["error"] = "This expression data contain duplicate index. Please remove these first."

        if len(commonGenes)==0:
            result["error"] = "No genes common between query data and atlas, likely due to row ids not being Ensembl ids."
            
        elif len(commonGenes)/len(genes[genes["inclusion"]])<0.5:
            result["error"] = f"Less than 50% of genes in query data are common with atlas ({len(commonGenes)} common). This may lead to an unreliable result."

        if result["error"]!="":
            return result

        # We reindex queryData on df.index, not on commonGenes, since pca is done on df and not on commonGenes. 
        # This means any genes in queryData not found in df will be dropped, and any genes in df not found in queryData will be assigned Nan
        # - we convert these to zero and live with this, as long as there aren't so many!
        dfQuery = rankTransform(queryData.reindex(df.index).fillna(0))

        if includeRandomSamples:
            n = 3 if len(dfQuery.columns)>2 else len(dfQuery.columns)  # max 3 random columns chosen from dfQuery
            for i,col in enumerate(random.sample(dfQuery.columns.tolist(), n)):  # add new column based on values of these columns
                dfQuery[f"random_{i}"] = random.sample(dfQuery[col].tolist(), len(dfQuery))

        # perform pca on atlas
        from sklearn.decomposition import PCA
        pca = PCA(n_components=10, svd_solver='full')
        coords = pandas.DataFrame(pca.fit_transform(df.values.T), index=df.columns)  # can also just run fit
        
        # make projection
        result["coords"] = pandas.DataFrame(pca.transform(dfQuery.values.T)[:,:3], index=dfQuery.columns)

        if includeCombinedCoords:   # also return row concatenated data frame of atlas+projection.
            projectedCoords = result['coords']
            projectedCoords.index = ["%s_%s" % (name, item) for item in projectedCoords.index]
            coords = coords.iloc[:,:3]
            coords.columns = projectedCoords.columns
            result['combinedCoords'] = pandas.concat([coords, projectedCoords])

        if includeCapybara:
            cutoff = 0.01   # drop any columns with max less than this, as these aren't useful.
            
            # Not sure if this works that well. Leave out for now.
            # if includeRandomSamples:
            #     for i,col in enumerate(random.sample(queryData.columns.tolist(), n)):  # add new column based on values of these columns
            #         queryData[f"random_{i}"] = random.sample(queryData[col].tolist(), len(queryData))

            capy = self.capybara(queryData)
            # drop columns with very low values
            for column,df in capy.items():
                df = df[[col for col in df.columns if df[col].max()>cutoff]]
            result['capybara'] = capy

        return result

    def capybara(self, query, rankNormalise=True):
        """Calculates capybara score for each sample in query against each column of atlas samples table
        and return the result as a dictionary, keyed on sample column and valued as a pandas DataFrame. 
        Implementation of https://doi.org/10.1101/2020.02.17.947390
        
        > self.capybara(external_expression_dataframe).
        The returned dictionary has same keys as columns of self.sampleMatrix().
        The value of dictionary is a DataFrame with columns of query as rows, and unique values of the atlas sample group.
        So each value is the score of that query column against each unique value of atlas sample group.

        Parameters
        ----------         
        query
            Gene expression dataframe of query data that will be classified.
        """
        
        from cvxopt import matrix, solvers
        solvers.options['show_progress'] = False
        #solvers.options['abstol'] = 1.e-10
        #solvers.options['reltol'] = 1.e-10
        #solvers.options['feastol'] = 1.e-10

        # Get reference and query matrices and subset both on columns of
        df = self.expressionMatrix(filtered=True)
        samples = self.sampleMatrix()
        df = df[samples.index]  # subset columns on index of samples (shouldn't have to be done for atlas, but just in case)
        query = query.loc[df.index.intersection(query.index)]
        df = df.loc[query.index]
    
        if rankNormalise:
            query = rankTransform(query)
    
        result = {}
        for column in samples.columns:
            anno = samples[column]
            n = anno.unique().shape[0]
            reference = pandas.DataFrame(index=anno.unique(), columns=df.index)
            for i_celltype in anno.unique():
                reference.loc[i_celltype,:] = df.loc[:,anno==i_celltype].mean(axis=1)
            P = matrix(reference.dot(reference.transpose()).values.astype(numpy.double))
            G = matrix(0.0, (n,n))
            G[::n+1] = -1.0
            h = matrix(0.0, (n,1))
            A = matrix(1.0, (1,n))
            b = matrix(1.0)
            ##########
            # New code for Jarny here, PA
            G = matrix(numpy.vstack((G,A)).astype(numpy.double))
            h = matrix(numpy.vstack((h,b)).astype(numpy.double)) 
            ##########
            cell_score = pandas.DataFrame(columns=anno.unique(), index=query.columns)
            for i in range(query.shape[1]):
                y = -1.0*query.iloc[:,i].values.astype(numpy.double)
                q = matrix(numpy.matmul(reference.values, y).astype(numpy.double))
                ##########
                # New code for Jarny here, PA
                # I included a numpy.round because the tiny negative residuals were annoying me
                solution = solvers.qp(P, q, G, h, show_progress=False)
                cell_score.iloc[i,:] = numpy.round(numpy.array(solution['x']).reshape(-1),6)
                ##########
                # Old code for comparison
                #solution = solvers.qp(P, q, G, h, A, b, show_progress=False)
                #cell_score.iloc[i,:] = numpy.array(solution['x']).reshape(-1)
            result[column] = cell_score
    
        return result

    def heatmapData(self, geneIds=[], clusterColumns=False, groupby=None, subsetby=None, subsetbyItem=None, relativeValue='zscore'):
        """Return a data frame which is suitable for rendering a heatmap (though it could be used in other cases) by specifying geneIds (Ensembl ids). 
        Columns are sorted according to specified ordering in the atlas or hierarchically clustered if clusterColumns is True. 
        """
        geneInfo = self.geneInfo()

        if len(geneIds)==0: return {}
        if not geneIds[0].startswith("ENS"): # assume that gene ids are symbols - we can actually work out ids from geneInfo
            geneIds = geneInfo[geneInfo['symbol'].isin(geneIds)].index.tolist()
            
        # Split geneIds into genes found/not found in the atlas, then filtered/unfiltered
        filtered = set(geneIds).intersection(set(geneInfo[geneInfo['inclusion']].index))
        unfiltered = set(geneIds).intersection(set(geneInfo[~geneInfo['inclusion']].index))
        notFound = set(geneIds).difference(filtered.union(unfiltered))

        # Work out subset of expression matrix to use
        df = self.expressionMatrix(filtered=True).loc[list(filtered)]
        samples = self.sampleMatrix()

        # Apply any subsetby
        if subsetby is not None:
            samples = samples[samples[subsetby]==subsetbyItem]
            df = df[samples.index]

        # Apply any groupby (eg. if groupby='treatment', mean of samples for each treatment will be calculated)
        if groupby is not None:
            df = df[samples.index].groupby(samples[groupby], axis=1).mean()

        if len(df.columns)<=1:  # application of subset + groupby reduced the matrix too much
            return {'dataframe':pandas.DataFrame(), 'error':'Not enough samples to render a heatmap.'}
        
        # Transform values according to relativeValue
        if relativeValue=='rowAverage':
            df = df.sub(df.mean(axis=1), axis=0)
        elif relativeValue=='zscore':  # use zscore on each row
            from scipy.stats import zscore
            df = df.apply(zscore, axis=1)
        elif relativeValue in df.columns:  # subtract this for each value
            df = df.sub(df[relativeValue], axis=0)
        
        # Ordering of the columns of df may be specified or clustered
        ordered = self.coloursAndOrdering()['ordering']
        if clusterColumns:
            df = df[hclusteredRows(df.transpose())]
        elif groupby in ordered: # use ordering specified
            df = df[ordered[groupby]]

        if relativeValue in df.columns: # Ensure that relativeValue is the first column
            columns = df.columns.tolist()
            columns.remove(relativeValue)
            df = df[sum([[relativeValue], columns], [])]

        # Get gene symbols
        geneSymbolFromId = geneInfo.loc[df.index, 'symbol'].fillna('').to_dict()
        geneSymbols = [geneSymbolFromId.get(geneId, geneId) for geneId in df.index]

        return {'filtered':list(filtered), 'unfiltered':list(unfiltered), 'notFound':list(notFound), 'geneSymbols':geneSymbols, 'dataframe':df}

# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------
def test_atlas():
    """Run tests for the current files in the atlas.
    Since we can't pass an argument to this function through nosetests, use environ variables to test for a particular atlas type.
    > export ATLAS_TYPE=dc
    (don't forget to also export ATLAS_FILEPATH=/path/to/atlas/files)
    """
    if 'ATLAS_TYPE' not in os.environ:
        print("test_atlas works on a single atlas specified by environ variable, for example:")
        print("> export ATLAS_TYPE=dc")
        print("Try again after setting this variable")
        return
    atlas = Atlas(os.environ['ATLAS_TYPE'])
    
    # Check that all files exist
    filelist = ['colours.json','coordinates.tsv','expression.tsv','expression.filtered.tsv','genes.tsv','samples.tsv','Readme.txt']
    assert set(filelist).issubset(set(os.listdir(atlas.atlasFilePath)))

    # Read all the files
    samples = atlas.sampleMatrix()
    exp = atlas.expressionMatrix()
    expFilt = atlas.expressionMatrix(filtered=True)
    pca = atlas.pcaCoordinates()
    colours = atlas.coloursAndOrdering()['colours']
    ordering = atlas.coloursAndOrdering()['ordering']
    genes = atlas.geneInfo()

    # Check that sample ids match in all files
    assert set(samples.index)==set(exp.columns)
    assert set(samples.index)==set(expFilt.columns)
    assert set(samples.index)==set(pca.index)

    # Check that gene ids match between exp and genes
    assert set(genes.index)==set(exp.index)
    assert set(genes[genes['inclusion']].index)==set(expFilt.index)

    # Check that colours and ordering have values which match those in samples
    assert set(colours.keys()).issubset(set(samples.columns))
    for key,val in colours.items():
        assert set(val.keys()).issubset(set(samples[key]))

    assert set(ordering.keys()).issubset(set(samples.columns))
    for key,val in ordering.items():
        assert set(val).issubset(set(samples[key]))

def test_projection():
    from models import datasets
    atlas = Atlas('dc')
    ds = datasets.Dataset(2000)
    df = ds.expressionMatrix(key='genes')
    res = atlas.projection('test', df, includeCombinedCoords=False)
    assert res['coords'].shape == (127,3)  # dataset has 124 samples but we added 3 random samples
    assert round(res['capybara']['Cell Type'].at['2000_1699538155_A','cDC1']*1000)==351

def test_capybara():
    from models import datasets
    atlas = Atlas('dc')
    ds = datasets.Dataset(2000)
    df = ds.expressionMatrix(key='genes')
    import time
    t0 = time.time()
    output = atlas.capybara(df)
    print(len(output), time.time()-t0)
    assert output['Cell Type'].at['2000_1699538155_A','monocyte']<0.001
    assert round(output['Cell Type'].at['2000_1699538155_C','cDC1']*100)==48

def test_heatmapData():
    atlas = Atlas('ma')
    atlas.heatmapData(genesetCollection='Hallmark', genesetName='HALLMARK_ADIPOGENESIS', groupby='treatment', relativeValue='NS')

def test_rankTransform():
    #df = pandas.DataFrame([[0, 0, 0, 3, 5, 7, 11]]).transpose()
    df = pandas.DataFrame([[2,2,5,10]]).transpose()
    print((df.rank(axis=0, ascending=True, na_option='bottom')-1)/(df.shape[0]-1))
    print(df.shape[0] - df.rank(axis=0, ascending=False, na_option='bottom')+1)
    print(rankTransform(df))
