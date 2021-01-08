Project template described [here](https://flask-restful.readthedocs.io/en/latest/intermediate-usage.html) has been used to create the structure in this project.

We'll also need /expression/[datasetId]?type=raw&geneId=ENSG00000023242 to return one row of values from the expression matrix, where rowId is one of the row index elements of the matrix (usually a gene id like ENSG0000034353). Now we have to see how quick this query is - currently with each expression matrix stored as tsv file, this API call will need to read this file using pandas.read_csv(), then return the values using .loc, so the time it takes for this query will be the same as /expression/raw/[datasetId] even though we're only returning one row. If too slow, we may have to look at ways of optimising this later.


conda install waitress
conda install -c conda-forge flask-restful
conda install pymongo
conda install nose
conda install -c conda-forge python-dotenv
conda install pandas