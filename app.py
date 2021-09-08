from flask import Flask, render_template, request

from flask_restful import Api
#from flask_cors import CORS
import os, logging, datetime

# Load environment vars in .env file. Even though load_dotenv function call is not even necessary
# when the app is called directly by python app.py, it is necessary for nohup.
from dotenv import load_dotenv
load_dotenv()

from resources import datasets, atlases, genes, auth, governance
from resources.errors import errors

app = Flask(__name__)
api = Api(app, errors=errors)
#cors = CORS(app)

# Added this to enable larger file uploads (this makes the limit 100Mb).
# According to this https://stackoverflow.com/questions/31873989/rejecting-files-greater-than-a-certain-amount-with-flask-uploads
# not setting this should allow any size upload, but the server returns 413 error for larger file uploads.
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024
app.config['JSONIFY_PRETTYPRINT_REGULAR'] = True  # seems to ignore this unless debug mode is on

#logging.basicConfig(filename='app.log', level=logging.ERROR, format=f'%(asctime)s %(levelname)s : %(message)s', datefmt='%Y-%m-%d_%H:%M')
logger = logging.getLogger('werkzeug') # grabs underlying WSGI logger
fh = logging.FileHandler('app.log')
fh.setFormatter(logging.Formatter('%(asctime)s %(levelname)s : %(message)s', '%Y-%m-%d_%H:%M'))
logger.setLevel(logging.WARNING)
logger.addHandler(fh)

# Get tables for a dataset with id
api.add_resource(datasets.DatasetMetadata, '/datasets/<int:datasetId>/metadata')
api.add_resource(datasets.DatasetSamples, '/datasets/<int:datasetId>/samples')
api.add_resource(datasets.DatasetExpression, '/datasets/<int:datasetId>/expression')
api.add_resource(datasets.DatasetPca, '/datasets/<int:datasetId>/pca')
api.add_resource(datasets.DatasetCorrelatedGenes, '/datasets/<int:datasetId>/correlated-genes')
api.add_resource(datasets.DatasetTtest, '/datasets/<int:datasetId>/ttest')

# Dataset and sample search
api.add_resource(datasets.DatasetSearch, '/search/datasets')
api.add_resource(datasets.SampleSearch, '/search/samples')

# Get available values eg: /values/samples/cell_type
api.add_resource(datasets.Values, '/values/<collection>/<key>')

# Download multiple datasets
api.add_resource(datasets.Download, '/download')

# Gene expression analyses
api.add_resource(genes.SampleGroupToGenes, '/genes/sample-group-to-genes')
api.add_resource(genes.GeneToSampleGroups, '/genes/gene-to-sample-groups')
#api.add_resource(genes.GenesetCollection, '/genes/geneset-collection')

# Atlas data
api.add_resource(atlases.AtlasTypes, '/atlas-types')
api.add_resource(atlases.Atlas, '/atlases/<atlasType>/<item>')
api.add_resource(atlases.AtlasProjection, '/atlas-projection/<atlasType>/<dataSource>')

# API for authentication
api.add_resource(auth.AuthLogin, '/auth/login')
api.add_resource(auth.AuthLogout, '/auth/logout')
api.add_resource(auth.AuthUser, '/auth/user')

# Dataset governance pages (require authentication)
api.add_resource(governance.DatasetSummary, '/governance/summary')
api.add_resource(governance.DatasetReport, '/governance/<int:datasetId>/report')
api.add_resource(governance.DatasetQCHtml, '/governance/<int:datasetId>/html')
api.add_resource(governance.DatasetQCPCA, '/governance/<int:datasetId>/pca')

# Static html at the home page level
@app.route("/")
def homepage():
    return render_template('index.html')

# Record URLs before each request
@app.before_request
def request_tracer():
    pathsToIgnore = ['','/auth/user','favicon.ico']
    if request.path not in pathsToIgnore:
        with open("app.log","a") as filehandle:
            filehandle.write(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")} {request.remote_addr} {request.full_path}\n')

if __name__ == '__main__':
    app.run(debug=True, port=5000)
