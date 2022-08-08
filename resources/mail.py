from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from dotenv import load_dotenv
import smtplib, ssl, logging, os

# ----------------------------------------------------------
# Steps to allow mail.py to work
# ----------------------------------------------------------
# -- (Organisation Gmail Account) Allow less secure apps in Gmail settings --
# 1. Go to Google account settings.
# 2. Go to 'Security'.
# 3. Find the 'Less secure app access' tab.
# 4. Turn access 'On'.

# -- (Personal Gmail Account) Obtain app password from Gmail settings --
# 1. Go to Google account settings.
# 2. Go to 'Security'.
# 3. Click on 'App passwords' under 'Signing in to Google'.
# 4. Select an app ('Mail'), and select a device '(Windows Computer' or 'Other').
# 5. Copy and paste the generated password below under the password variable, or save to .env.

# -- Sending limits --
# There is a sending limit of either 500 emails per day, or 500 recipients in a single email.
# Per day is defined as every rolling 24 hours.

# -- Todos for mail.py --
# 1. Check the security issues that arise when allowing less secure app access.

# Load .env variables
load_dotenv()

# @param list recipients The recipients email addresses as a list of strings.
# A bare string will be interpreted as a list with 1 address.
# Must be either a list of strings, or a single address as a bare string, any other format will not work properly.
# @param str subject The subject of the email.
# @param str body The text body of the email.
def send_mail(recipients, subject, body):

  # Set up logger
  LOG_FILENAME = 'mail.log'
  logger = logging.getLogger(LOG_FILENAME)
  logger.setLevel(logging.INFO)
  file_handler = logging.FileHandler(LOG_FILENAME)
  formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%d/%m/%Y %H:%M:%S')
  file_handler.setFormatter(formatter)
  logger.addHandler(file_handler)

  # Use Gmail SMTP server, with port 587 for TLS
  smtp_server = 'smtp.gmail.com'
  port = 587 

  # Sender email and password
  sender_email = os.environ["EMAIL_USER"]
  password = os.environ["EMAIL_PASS"]

  # Check if input parameter is a bare string, then put in a list
  if(isinstance(recipients,str)):
    recipients = [recipients]

  # Create message
  msg = MIMEMultipart()
  msg["Subject"] = subject
  msg["From"] = 'Stemformatics' + f' <sender_email>' # Set custom email 'From' ID 
  msg['To'] = ", ".join(recipients)

  # Attach body to message
  body_text = MIMEText(body, 'plain')
  msg.attach(body_text)

  context = ssl.create_default_context()

  # Try to log in to server and send email
  try:
    server = smtplib.SMTP(smtp_server, port)
    server.ehlo()  # check connection
    server.starttls(context=context)  # Secure the connection
    server.ehlo()  # check connection
    server.login(sender_email, password)

    # Send email
    response = server.sendmail(sender_email, recipients, msg.as_string())

  except Exception as e:
    # Log any error messages to mail.log
    logger.info(e)
  finally:
    # Log successful send mail to mail.log
    logger.info(recipients)
    server.quit()
    return response

# ----------------------------------------------------------
# tests: eg. $nosetests -s <filename>:ClassName.func_name
# ----------------------------------------------------------
# May have to use export MONGO_URI='xxxx' before running these tests, in order to set environment variables. 
# See .env file for a full list of variables to set.

# send_mail() should return a dictionary with one entry for each recipient that was refused - empty if all successful
def test_mail():

  # response = send_mail(['jbarry1@student.unimelb.edu.au', 'jake.barry95@gmail.com'], 'Testing Email', 'This is a testing email. Please do not reply.')
  response = send_mail('jbarry1@student.unimelb.edu.au', 'Testing Email', 'This is a testing email. Please do not reply.')
  assert response == {}
