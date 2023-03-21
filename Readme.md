# Stemformatics API Server

This is a flask-restful based API server, designed to serve [Stemformatics](http://stemformatics.org) data to the UI server.

## Building the environment
conda is used for the environment and the necessary packages are:

```bash
conda install -c conda-forge flask-restful
conda install flask-cors
conda install flask-bcrypt
conda install pymongo
conda install nose
conda install -c conda-forge python-dotenv
conda install pandas
conda install pyjwt
conda install scp
conda install requests
conda install waitress
conda install scikit-learn
conda install -c conda-forge cvxopt
```

## Running the server
For development, you can run from command line, which will print any errors to the command line:
```bash
python app.py
```

For test/production, use waitress to run it more permanently
```bash
nohup waitress-serve --port=5000 app:app > waitress.log 2>&1 &

# manually kill the server by finding it id:
ps -fu ec2-user | grep waitress
```

## Notes

Project template described [here](https://flask-restful.readthedocs.io/en/latest/intermediate-usage.html) has been used to create the structure in this project.

API via python requests
https://realpython.com/python-requests/#ssl-certificate-verification

Rotating logs
https://stackoverflow.com/questions/42797276/flask-how-to-write-werkzeug-logs-to-log-file-using-rotatingfilehandler

Authentication including token generation
https://scotch.io/tutorials/build-a-restful-api-with-flask-the-tdd-way-part-2
https://realpython.com/token-based-authentication-with-flask/#jwt-setup

To access a private resource using auth token, copy the token from /datasets/governance and add auth_token key into the URL (value starts after "Bearer"):
https://api-dev.stemformatics.org/datasets/9762/samples?as_file=true&auth_token=eyJ0eXAiOiJKV1QiL...


If the api suddenly starts returning 400 bad request, it could be that werkzeug got upgraded. See [here](https://stackoverflow.com/questions/72157708/flask-restx-request-parser-returns-400-bad-request/72174259#72174259). This can be fixed by:
```bash
conda install werkzeug=2.0.*
```


Return a string as file in response
```
from flask import make_response

response = make_response(df.to_json(orient=args.get('orient')))
response.headers['content-type'] = 'application/octet-stream'
return response
```