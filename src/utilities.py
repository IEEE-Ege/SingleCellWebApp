# utilities.py
import jwt
from datetime import datetime, timedelta


# Function to create a JWT token
def create_jwt_token(username: str, SECRET_KEY) -> str:
    # Payload for the JWT token
    payload = {
        "user": username,
        "exp": datetime.datetime.utcnow() + datetime.timedelta(hours=1), # Token expires in 1 hour
    }
    # Encode the payload with the secret key and HS256 algorithm
    return jwt.encode(payload, SECRET_KEY, algorithm="HS256")