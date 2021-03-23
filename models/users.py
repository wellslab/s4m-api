"""
Main interface to user data.
"""
import pymongo, os, pandas
from flask_bcrypt import Bcrypt

from models.utilities import mongoClient

database = mongoClient()["dataportal"]

def isValid(username, password):
    result = database['users'].find_one({'username':username}, {"_id":0})
    if result:
        return Bcrypt().check_password_hash(result.get('password'), password)
    return False