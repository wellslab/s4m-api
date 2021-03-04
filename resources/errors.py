from werkzeug.exceptions import HTTPException

class DatasetIdNotFoundError(HTTPException):
    pass

errors = {
    "DatasetIdNotFoundError": {
        "message": "Dataset id not found in the system",
        "status": 400
    }
}