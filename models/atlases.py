"""Model to handle atlas data in Stemformatics. Atlas data refers to an object of integrated data, 
and comprises of expression and sample matrices, as well as supporting data such as colours and genes. 

Examle usage:
--------------

# Note that we need to provide the path to the files which contains the atlas data.
> export ATLAS_FILEPATH=/path/to/atlas-files
atlas = Atlas('myeloid')
print(atlas.pcaCoordinates().head())
"""
import os, pandas, json

class Atlas(object):

    # Full list of current atlas types
    all_atlas_types = ['myeloid','blood','dc']

    def __init__(self, atlasType):
        # each atlas type is under its own directory under ATLAS_FILEPATH
        self.atlasFilePath = os.path.join(os.environ["ATLAS_FILEPATH"], atlasType)

        # Work out version based on link
        self.version = os.readlink(self.atlasFilePath).split("_")[1]

        self.atlasType = atlasType

    def pcaCoordinates(self):
        """Return a pandas DataFrame object, specifying the PCA coordinates. The data frame will have sample ids as index.
        Columns will be named as ['0','1',...] (as strings).
        """
        df = pandas.read_csv(os.path.join(self.atlasFilePath, "coordinates.tsv"), sep="\t", index_col=0)
        df.columns = [str(i) for i in range(len(df.columns))]
        return df

    def expressionFilePath(self, filtered=False):
        filename = "expression.filtered.tsv" if filtered else "expression.tsv"
        return os.path.join(self.atlasFilePath, filename)

    def expressionValues(self, geneIds, filtered=False):
        """Return expression values for geneIds as a pandas DataFrame. 
        If filtered=False, this is the "full" expression matrix that includes all genes. 
        The expression values are from rank normalised values.
        """
        df = pandas.read_csv(self.expressionFilePath(filtered=filtered), sep="\t", index_col=0)
        return df.loc[geneIds]

    def datasetIds(self):
        """Return all dataset ids in this atlas as a list. Note that each element will be integer type.
         Relies on sample ids in the form of "datasetId_sampleId", matching the format of sampleId in sample data.
        """
        return [int(item.split("_")[0]) for item in self.sampleMatrix().index]

    def sampleMatrix(self):
        """Return sample annotation matrix for the atlas as a pandas DataFrame object. The shape of the data frame
        will be number_of_samples x number_of_columns, with sample ids as index.
        """
        return pandas.read_csv(os.path.join(self.atlasFilePath, "samples.tsv"), sep="\t", index_col=0)   

    def geneInfo(self):
        """Return a pandas DataFrame of information about all genes in the atlas, after reading the genes.tsv
        file in the atlas file directory. Ensembl ids form the index.
        """
        df = pandas.read_csv(os.path.join(self.atlasFilePath, "genes.tsv"), sep="\t", index_col=0)
        df.index.name = 'ensembl'
        return df

    def coloursAndOrdering(self):
        """Return dictionaries of colours and ordering of sample type items based on "colours.tsv" file inside the
        atlas file directory. Return empty dictionaries if such a file doesn't exist.
        """
        filepath = os.path.join(self.atlasFilePath, "colours.tsv")
        return json.loads(open(filepath).read()) if os.path.exists(filepath) else {}

    def projection(self, name, testData, includeCombinedCoords=True):
        """Perform projection of testData onto this atlas and return a dictionary of objects.
        Params:
            name (str): Name of the testData, used as prefix for projected points if includeCombinedCoords is True.
            testData (DataFrame): Expression matrix of test data to be projected onto this atlas.
            includeCombinedCoords (bool): Should atlas coords and projected coords be returned together, rather than just projections?

        Returns a dictionary with following keys and values:
            coords (DataFrame): projected coordinates of testData.
            combinedCoords (DataFrame): data frame of atlas coordinates + projected coords if includeCombinedCoords is True.
                The projected points will have index in the format of "{name}_sample1", etc, so that these points
                can be distinguished from atlas points.
            error (str): Error message. Empty string if there was no error.
            name (str): Same as the value used as input.
        """
        result = {"error":"", "coords":pandas.DataFrame(), "name":name, "combinedCoords":pandas.DataFrame()}

        # Some validation before projecting
        if len(testData)==0:
            result["error"] = "Data to project appears to have 0 rows. Check format of the file."

        # Read expression matrix - we only need filtered version. 
        df = self.expressionTable(filtered=True)
        genes = self.geneInfo()

        commonGenes = testData.index.intersection(df.index)  # common index between test and atlas
        if len(commonGenes)==0:
            result["error"] = "No genes common between test data and atlas, likely due to row ids not in Ensembl ids."
            
        elif len(commonGenes)/len(genes[genes["inclusion"]])<0.5:
            result["error"] = "Less than 50% of genes in test data are common with atlas ({} common)".format(len(commonGenes))

        if result["error"]!="":
            return result

        # We reindex testData on df.index, not on commonGenes, since pca is done on df. This means any genes in df not found in 
        # testData will gene None assigned - we will live with this, as long as there aren't so many.
        dfTest = rankTransform(testData.reindex(df.index))
        expression = pandas.concat([df, dfTest], axis=1)

        # perform pca on atlas
        from sklearn.decomposition import PCA
        pca = PCA(n_components=10, svd_solver='full')
        coords = pandas.DataFrame(pca.fit_transform(df.values.T), index=df.columns)  # can also just run fit
        
        # make projection
        result["coords"] = pandas.DataFrame(pca.transform(dfTest.values.T)[:,:3], index=dfTest.columns, columns=['x','y','z'])

        if includeCombinedCoords:   # also return row concatenated data frame of atlas+projection.
            projectedCoords = result['coords']
            projectedCoords.index = ["%s_%s" % (name, item) for item in projectedCoords.index]
            coords = coords.iloc[:,:3]
            coords.columns = projectedCoords.columns
            result['combinedCoords'] = pandas.concat([coords, projectedCoords])

        return result

# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------
def test_samples():
    atlas = Atlas("myeloid")
    print(atlas.pcaCoordinates().head())