from werkzeug.exceptions import HTTPException

class DatasetIdNotFoundError(HTTPException):
    pass

class DatasetIsPrivateError(HTTPException):
    pass

errors = {
    "DatasetIdNotFoundError": {
        "message": "Dataset id not found in the system.",
        "status": 400
    },
    "DatasetIsPrivateError": {
        "message": "Dataset is private and cannot be accessed without authentication.",
        "status": 400
    }
}