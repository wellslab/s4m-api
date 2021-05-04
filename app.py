from flask import Flask
from flask_restful import Api
from flask_cors import CORS
import os, logging

# Load environment vars in .env file. Even though load_dotenv function call is not even necessary
# when the app is called directly by python app.py, it is necessary for nohup.
from dotenv import load_dotenv
load_dotenv()

from resources import datasets, atlases, genes, auth, governance
from resources.errors import errors

app = Flask(__name__)
api = Api(app, errors=errors)
#cors = CORS(app)

if os.getenv('FLASK_ENV')!='development':
    logging.basicConfig(filename='app.log', level=logging.ERROR, format=f'%(asctime)s %(levelname)s %(name)s : %(message)s')

# Get tables for a dataset with id
api.add_resource(datasets.DatasetMetadata, '/datasets/<int:datasetId>/metadata')
api.add_resource(datasets.DatasetSamples, '/datasets/<int:datasetId>/samples')
api.add_resource(datasets.DatasetExpression, '/datasets/<int:datasetId>/expression')
api.add_resource(datasets.DatasetPca, '/datasets/<int:datasetId>/pca')

# Dataset and sample search
api.add_resource(datasets.DatasetSearch, '/search/datasets')
api.add_resource(datasets.SampleSearch, '/search/samples')

# Get available values
api.add_resource(datasets.ValuesDatasets, '/values/datasets/<key>')
api.add_resource(datasets.ValuesSamples, '/values/samples/<key>')

# Atlas data
api.add_resource(atlases.Atlas, '/atlases/<atlasType>/<item>')
api.add_resource(atlases.AtlasTypes, '/atlas-types')
api.add_resource(atlases.AtlasProjection, '/atlas-projection/<atlasType>/<dataSource>')

# Gene annotation data - use mygene.info instead of using my own
# api.add_resource(genes.Geneset, '/genes')

# API for authentication
api.add_resource(auth.AuthLogin, '/auth/login')
api.add_resource(auth.AuthLogout, '/auth/logout')
api.add_resource(auth.AuthUser, '/auth/user')

# Dataset governance pages (require authentication)
api.add_resource(governance.DatasetSummary, '/governance/summary')
api.add_resource(governance.DatasetReport, '/governance/<int:datasetId>/report')
api.add_resource(governance.DatasetQCHtml, '/governance/<int:datasetId>/html')
api.add_resource(governance.DatasetQCPCA, '/governance/<int:datasetId>/pca')

if __name__ == '__main__':
    app.run(debug=True, port=5000)
