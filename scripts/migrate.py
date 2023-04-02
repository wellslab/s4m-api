"""
Script to migrate code and data from one server to another. Use it like this on the destination server
(so it should be run under the s4m-api environment)
    (s4m-api) [ec2-user@api-test s4m-api]$ python -m scripts.migrate

It will ask for the source server, which will be used to copy data from, if that is required.
For updating code, it will use git pull, which will ignore the source server.
"""

import os, subprocess, scp
from paramiko import SSHClient
from scp import SCPClient
from datetime import date

source = None
_ssh = None

def _getSSH():
    """Make a ssh connection to the source server and return the SSHClient object
    """
    global _ssh
    if _ssh is None:
        _ssh = SSHClient()
        _ssh.load_system_host_keys()
        _ssh.connect(f'{source}.stemformatics.org')
    return _ssh

def rsyncFiles():
    """Expression files are copied from source server using rsync. 
    Example rsync usage on command line (note trailing / on the source directory!):
    rsync -avz dev.stemformatics.org:/mnt/stemformatics-data/atlas/ /mnt/stemformatics-data/atlas
    """
    answer = input("rsync expression files? [N]/y ")
    if (answer=='y'):
        filepath = os.environ['EXPRESSION_FILEPATH']
        command = ["rsync","-avz",f"{source}.stemformatics.org:{filepath}/", filepath]
        print(subprocess.list2cmdline(command))
        process = subprocess.run(command)

    answer = input("rsync atlas files? [N]/y ")
    if (answer=='y'):
        filepath = os.environ['ATLAS_FILEPATH']
        command = ["rsync","-avz",f"{source}.stemformatics.org:{filepath}/", filepath]
        print(subprocess.list2cmdline(command)) 
        process = subprocess.run(command)

def copyMongoData():
    """Metadata are copied from source server using dump to text then read from text.
    """
    for key in ['datasets','samples']:
        answer = input(f"Copy {key} metadata files? [N]/y ")
        if (answer=='y'):
            # Make a backup of metadata at source and copy here
            filename = f"{key}_{date.today().strftime('%Y%m%d')}"    # eg. datasets_20210902

            # First check if backup file already exists
            # Runs a bash command like this on source: test -f /mnt/stemformatics-data/expression_files/../backups/datasets_20210902 && echo 'true'
            # which will return 'true' if the file exists
            filepath = f"{os.environ['EXPRESSION_FILEPATH'].replace('expression_files','backups')}/{filename}.tsv"
            command = f"test -f {filepath} && echo 'true'" 
            print(command)
            ssh = _getSSH()
            stdin, stdout, stderr = ssh.exec_command(command)
            output = stdout.readlines()
            createBackupFile = True
            if output==['true\n']:
                if input(f"File at source ({filepath}) already exists. Use this? [N]/y ")=='y':
                    createBackupFile = False
                else:   # just quit the program and handle the existing file first manually
                    return

            if createBackupFile:    # create backup
                command = f"cd; cd s4m-api; conda activate s4m-api; python -m scripts.backup_and_restore backupCollectionToCSV dataportal {key} {filepath}"
                print(command)
                stdin, stdout, stderr = ssh.exec_command(command)

            # scp here
            localfile = f"{os.environ['EXPRESSION_FILEPATH'].replace('expression_files','received')}/{source}/{filename}.tsv"
            scp = SCPClient(ssh.get_transport())
            scp.get(filepath, local_path=localfile)

            # Now run backup_and_restore locally
            from scripts.backup_and_restore import createCollectionFromCSV
            createCollectionFromCSV('dataportal',key,localfile)

def gitPull():
    """Code update is done through git pull (easier for public repositories)
    """
    for key in ["s4m-api", "s4m-ui"]:
        answer = input(f"git pull {key}? [N]/y ")
        if (answer=='y'):
            process = subprocess.run(f"cd; cd {key}; git pull", shell=True)

def restartServers():
    """Restart server after updates.
    """
    answer = input(f"restart s4m-api server? [N]/y ")
    if (answer=='y'):
        # find pid and kill it
        ps = subprocess.run(['ps', '-u'], stdout=subprocess.PIPE, universal_newlines=True)
        cols = [line for line in ps.stdout.splitlines() if 'waitress-serve' in line]
        if len(cols)>0:
            subprocess.run(['kill', cols[0].split()[1]])
        # restart
        subprocess.run("nohup waitress-serve --port=5000 app:app > waitress.log 2>&1 &", shell=True)

    answer = input(f"restart s4m-ui server? [N]/y ")
    if (answer=='y'):
        subprocess.run("cd; cd s4m-ui; source /mnt/miniconda3/bin/activate s4m-ui; npm run build; pm2 stop ecosystem.config.js; pm2 start", shell=True)

def main():
    rsyncFiles()
    print("\n")
    copyMongoData()
    print("\n")
    gitPull()
    print("\n")
    restartServers()

if __name__=="__main__":
    source = input("Set source of migration [dev,test,prod1,prod2]: ")
    if source in ['test','dev','prod1','prod2']:
        print(f"Source set to {source}.\n")
        main()
    else:
        print("Source should be either test or dev")
