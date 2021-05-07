"""
Script that imports newly created expression files from data-source.
Assumes that dataset and sample metadata are already present in the current system.

Examples of how to run this script (ensure you're in the application directory):
(Requires environment variable EXPRESSION_FILEPATH, which points to where the expression files are)
(s4m-api) [ec2-user@api-dev s4m-api]$ python -m scripts.import_expression_data 6122
"""

import os, sys, pandas, scp, argparse
from paramiko import SSHClient
from scp import SCPClient

sys.path.append(os.path.join(sys.path[0]))
from models import datasets

def importExpressionData(datasetId):
    """
    """
    if "EXPRESSION_FILEPATH" not in os.environ:
        print("No EXPRESSION_FILEPATH in environment")
        return

    ds = datasets.Dataset(datasetId)
    print("This will import expression files for dataset %s, %s [%s] [version %s] [status %s]" % \
        (datasetId, ds.metadata()['name'], ds.metadata()['platform_type'], ds.metadata()['version'], ds.metadata()['status']))
    if ds.metadata()['status']!='passed':
        print("This dataset does not have passed status.")
        return
        
    samples = ds.samples()
    print("Sample metadata shows %s samples" % len(samples))
    
    # Check if expression file exists - modify later to cope with version change rather than a completley new dataset
    fileExists = os.path.exists(ds.expressionFilePath())
    if fileExists:
        print("Expression file already exists at", ds.expressionFilePath())
        return
    
    # Make directory locally
    os.chdir(os.getenv('EXPRESSION_FILEPATH'))
    dirname = "%s_%s" % (datasetId, ds.metadata()['version'])
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    if not os.path.exists(str(datasetId)):
        os.symlink(dirname, str(datasetId), target_is_directory=True)

    # Copy relevant files from data-source
    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.connect('data-source.stemformatics.org')
    scp = SCPClient(ssh.get_transport())
    scp.get('/mnt/data/pending/%s/qc/pcas/%s.pca.tsv' % (dirname, datasetId), local_path=dirname)
    scp.get('/mnt/data/pending/%s/qc/pcas/%s.pca_attributes.tsv' % (dirname, datasetId), local_path=dirname)
    if ds.metadata()['platform_type']=='RNASeq':
        scp.get('/mnt/data/pending/%s/feature_counts/%s.raw.tsv' % (dirname, datasetId), local_path=dirname)
    else:
        scp.get('/mnt/data/pending/%s/tables/%s.raw.tsv' % (dirname, datasetId), local_path=dirname)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", help="dataset id")
    #parser.add_argument("-a", help="only look at atlas datasets", action="store_true")
    args = parser.parse_args()
    
    importExpressionData(datasetId=args.d)
