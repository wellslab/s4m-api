from flask import Flask
from flask_restful import Api
import logging

from resources.datasets import DatasetMetadata, DatasetSamples, DatasetExpression, DatasetGovernance, ValuesDatasets, ValuesSamples
from resources.search import SearchDatasets, SearchSamples
from resources.atlases import Atlas

app = Flask(__name__)
api = Api(app)

logging.basicConfig(filename='app.log', level=logging.ERROR, format=f'%(asctime)s %(levelname)s %(name)s %(threadName)s : %(message)s')

# Get tables for a dataset with id
api.add_resource(DatasetMetadata, '/datasets/<int:datasetId>/metadata')
api.add_resource(DatasetSamples, '/datasets/<int:datasetId>/samples')
api.add_resource(DatasetExpression, '/datasets/<int:datasetId>/expression')
api.add_resource(DatasetGovernance, '/datasets/<int:datasetId>/governance')

# Perform search
api.add_resource(SearchDatasets, '/search/datasets')
api.add_resource(SearchSamples, '/search/samples')

# Get available values
api.add_resource(ValuesDatasets, '/values/datasets/<key>')
api.add_resource(ValuesSamples, '/values/samples/<key>')

# Atlas data
api.add_resource(Atlas, '/atlases/<atlasType>/<item>')

if __name__ == '__main__':
    app.run(debug=True, port=5000)
