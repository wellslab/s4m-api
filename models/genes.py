"""
Main interface to gene expression data. Even though we can consider gene expression as a part of dataset,
performing gene expression analyses across multiple datasets happens often, so it gets its own package here.
"""
import pymongo, os, pandas

def geneset(geneIds=[], searchString="", limit=100):
    """Return a pandas DataFrame which contains attributes of genes matching the parameters.
    """
    geneIds = list(set([item for item in geneIds if item.startswith("ENS")]))

    # Don't include _id field, since this its value is ObjectId class which can't be
    # serialised by json, which can create issues downstream, and this field is not needed outside mongo anyway.
    if len(geneIds)>0:
        params = {"gene_id": {"$in": geneIds}}
    else:
        params = {'gene_name': {'$regex': searchString, '$options': 'i'}}

    if limit:
        cursor = database["genes"].find(params, {"_id":0}).limit(limit)
    else:
        cursor = database["genes"].find(params, {"_id":0})

    return pandas.DataFrame(cursor).set_index("gene_id") if cursor.count()!=0 else pandas.DataFrame()

