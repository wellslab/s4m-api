from werkzeug.exceptions import HTTPException

class DatasetIdNotFoundError(HTTPException):
    pass

class DatasetIsPrivateError(HTTPException):
    pass

class UserNotAuthenticatedError(HTTPException):
    pass

class DatasetQCFilesMissingError(HTTPException):
    pass

class DatasetProjectionFailedError(HTTPException):
    pass

errors = {
    "DatasetIdNotFoundError": {
        "message": "Dataset id not found in the system.",
        "status": 400
    },
    "DatasetIsPrivateError": {
        "message": "Dataset is private and cannot be accessed without authentication.",
        "status": 400
    },
    "UserNotAuthenticatedError": {
        "message": "User authentication failed.",
        "status": 400
    },
    "DatasetQCFilesMissingError": {
        "message": "QC files are missing for this dataset.",
        "status": 400
    },
    "DatasetProjectionFailedError": {
        "message": "Projection of dataset failed. Check file formats.",
        "stauts": 400
    }
}