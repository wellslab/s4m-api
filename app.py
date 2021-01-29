from flask import Flask
from flask_restful import Api

from resources.datasets import DatasetMetadata, Samples, Metadata
from resources.search import Search

app = Flask(__name__)
api = Api(app)

# Link resources with routes
api.add_resource(DatasetMetadata, '/dataset/<int:datasetId>/metadata', '/dataset/<int:datasetId>/metadata/', 
                                  '/dataset/<int:datasetId>', '/dataset/<int:datasetId>/')
api.add_resource(Samples, '/dataset/<int:datasetId>/samples', '/dataset/<int:datasetId>/samples/',
                          '/samples/<int:datasetId>', '/samples/<int:datasetId>/')
#api.add_resource(Expression, '/dataset/<int:datasetId>/expression', '/expression/<int:datasetId>')
api.add_resource(Search, '/search', '/search/')
api.add_resource(Metadata, '/metadata/<item>','/metadata/<item>/')
#api.add_resource(Atlas, '/atlas/myeloid')

if __name__ == '__main__':
    app.run(debug=True, port=5000)
