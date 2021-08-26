import os, re, pandas
from flask_restful import reqparse, Resource

from resources import auth
from resources.errors import UserNotAuthenticatedError, DatasetQCFilesMissingError
from models import datasets

# ----------------------------------------------------------
# Governance related dataset access after authentication
# ----------------------------------------------------------
"""All classes here inherit from Governance, which checks for authentication.
There are other ways of implementing this, but this works fine, and also enables
definition of other common methods.
"""
class Governance(Resource):
    def __init__(self):
        if not auth.AuthUser().username():
            raise UserNotAuthenticatedError

    def directoryFromDataset(self):
        """Return a dictionary of path to QC files, keyed on dataset id (integer)
        eg. {3331:'/path/to/qc/files/',...}.
        """
        # Fetch qc results. These are retrieved from data-source if there's not a local copy already.
        # (to-do: how do we update the local copy if the version at data-source changes?)
        dfd = {}
        for item in os.listdir(os.getenv('QC_FILEPATH')):
            match = re.findall("^\d{4}", item)
            if len(match)>0: dfd[int(match[0])] = os.path.join(os.getenv("QC_FILEPATH"), item)
        return dfd

class DatasetSummary(Governance):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str, required=False, default='pending')
        args = parser.parse_args()

        status = args.get('status').split(',') if args.get('status') else []
        df = datasets.datasetMetadataFromDatasetIds(datasets.datasetIdsFromFields(status=status, publicOnly=False), publicOnly=False)
        return df.sort_index().reset_index().fillna('').to_dict(orient='records')

class DatasetReport(Governance):
    def get(self, datasetId):
        # Fetch governance report on this dataset
        # to-do

        return {'history':['history of dataset processing']}

class DatasetQCHtml(Governance):
    def get(self, datasetId):

        parser = reqparse.RequestParser()
        parser.add_argument('type', type=str, required=False, default='multiqc')
        args = parser.parse_args()

        if args.get('type')=='multiqc':
            filename = 'multiqc_report.html'

        dfd = self.directoryFromDataset()        
        if dfd.get(datasetId):
            from flask import send_from_directory
            return send_from_directory(dfd[datasetId], filename, as_attachment=True)
        else:
            raise DatasetQCFilesMissingError

class DatasetQCPCA(Governance):
    def get(self, datasetId):
        # Read pca coordinates file and return it along with samples table
        dfd = self.directoryFromDataset()        
        if dfd.get(datasetId):
            pca = pandas.read_csv(os.path.join(dfd[datasetId],'pcas/%s.pca.tsv' % datasetId), sep='\t', index_col=0)
            pca.index.name = 'sampleId'
            return {'pca':pca.to_dict(orient='dict')}
        else:
            raise DatasetQCFilesMissingError


def test_tokens():
    token = encode_auth_token('jarny')
    print(token)
    print(decode_auth_token(token))