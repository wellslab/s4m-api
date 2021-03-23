from flask_restful import reqparse, Resource

from resources.errors import UserNotAuthenticatedError

# ----------------------------------------------------------
# Dataset access after authentication
# ----------------------------------------------------------
class DatasetSummary(Resource):
    def get(self):
        from flask import send_from_directory
        return send_from_directory("/mnt/stemformatics-data/received/qc","multiqc_report.html")

class DatasetReport(Resource):
    def get(self, datasetId):
        parser = reqparse.RequestParser()
        parser.add_argument('Authorization', location='headers', default='')  # eg {'Authorization':'Bearer yJ0eXAiOiJKV1QiLCJhbGciOiJIU'}
        args = parser.parse_args()
        try:
            authorizationHeader = args.get('Authorization').split(' ')
            if decode_auth_token(authorizationHeader[1]):
                from models import datasets
                ds = datasets.Dataset(datasetId)
                return ds.metadata()
        except:
            raise UserNotAuthenticatedError

def test_tokens():
    token = encode_auth_token('jarny')
    print(token)
    print(decode_auth_token(token))