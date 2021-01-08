"""
Collection of mostly small functions for convenience used by all models.
"""
import pymongo, os

def mongoClient():
    return pymongo.MongoClient(os.environ.get("MONGO_URI"))
