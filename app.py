from flask import Flask
from flask_restful import Api
import logging

# Load environment vars in .env file. Even though load_dotenv function call is not even necessary
# when the app is called directly by python app.py, it is necessary for nohup.
from dotenv import load_dotenv
load_dotenv()

from resources.datasets import DatasetMetadata, DatasetSamples, DatasetExpression, DatasetGovernance, DatasetPca, DatasetSearch,\
                               ValuesDatasets, ValuesSamples
from resources.atlases import Atlas
from resources.genes import Geneset
from resources.errors import errors

app = Flask(__name__)
api = Api(app, errors=errors)

logging.basicConfig(filename='app.log', level=logging.ERROR, format=f'%(asctime)s %(levelname)s %(name)s : %(message)s')

# Get tables for a dataset with id
api.add_resource(DatasetMetadata, '/datasets/<int:datasetId>/metadata')
api.add_resource(DatasetSamples, '/datasets/<int:datasetId>/samples')
api.add_resource(DatasetExpression, '/datasets/<int:datasetId>/expression')
api.add_resource(DatasetGovernance, '/datasets/<int:datasetId>/governance')
api.add_resource(DatasetPca, '/datasets/<int:datasetId>/pca')

# Perform search
api.add_resource(DatasetSearch, '/search')

# Get available values
api.add_resource(ValuesDatasets, '/values/datasets/<key>')
api.add_resource(ValuesSamples, '/values/samples/<key>')

# Atlas data
api.add_resource(Atlas, '/atlases/<atlasType>/<item>')

# Gene annotation data
api.add_resource(Geneset, '/genes')

if __name__ == '__main__':
    app.run(debug=True, port=5000)
