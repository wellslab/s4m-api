from flask_restful import reqparse, Resource
import jwt, datetime, os

from resources.errors import UserNotAuthenticatedError
from models import users

# ----------------------------------------------------------
# Token encoding and decoding
# ----------------------------------------------------------
def encode_auth_token(username):
    """
    Generates the Authentication Token.
    :return: string
    """
    payload = {
        'exp': datetime.datetime.utcnow() + datetime.timedelta(hours=2),
        'iat': datetime.datetime.utcnow(),
        'sub': username
    }
    return jwt.encode(payload, os.getenv('SECRET_KEY'), algorithm='HS256')

def decode_auth_token(auth_token):
    """
    Decodes the auth token
    :param auth_token:
    :return: integer|string
    """
    try:
        payload = jwt.decode(auth_token, os.getenv('SECRET_KEY'), algorithms=['HS256'])
        return payload['sub']
    except jwt.ExpiredSignatureError:
        raise Exception('Signature expired. Please log in again.')
    except jwt.InvalidTokenError:
        raise Exception('Invalid token. Please log in again.')

# ----------------------------------------------------------
# Authentication interactions
# ----------------------------------------------------------
class AuthLogin(Resource):
    def post(self):
        """Login with username and password
        """
        parser = reqparse.RequestParser()
        parser.add_argument('username', type=str)
        parser.add_argument('password', type=str)
        args = parser.parse_args()
        
        try:
            if users.isValid(args.get('username'), args.get('password')):
                return {'token':encode_auth_token(args.get('username'))}
            else:
                raise UserNotAuthenticatedError
        except:
            raise UserNotAuthenticatedError

class AuthLogout(Resource):
    def delete(self):
        """Logout
        """
        return {}

class AuthUser(Resource):
    """
    Use username method to check if user is authenticated before allowing access to a restricted resource by checking
    the token. Token may be present either as a parameter with key 'auth_token' or in the header as key 'Authorization'.
    """
    def username(self):
        parser = reqparse.RequestParser()
        parser.add_argument('auth_token', type=str, required=False)
        parser.add_argument('Authorization', type=str, location='headers', default='')  # eg {'Authorization':'Bearer yJ0eXAiOiJKV1QiLCJhbGciOiJIU'}
        args = parser.parse_args()

        token = None
        try:
            if args.get('auth_token'): # token present, so try to decode it
                token = args.get('auth_token').split(' ')[1] if ' ' in args.get('auth_token') else args.get('auth_token')
            else: # look at header
                token = args.get('Authorization').split(' ')[1]
            return decode_auth_token(token)
        except:
            return None

    def get(self):
        username = self.username()
        if username:
            return {'username': username}
        else:
            raise UserNotAuthenticatedError


def test_tokens():
    token = encode_auth_token('jarny')
    print(token)
    print(decode_auth_token(token))