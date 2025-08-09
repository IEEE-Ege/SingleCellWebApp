# utilities.py
import jwt
from datetime import datetime, timedelta
import os
from dotenv import load_dotenv

SECRET_KEY = os.getenv("SECRET_KEY")

# Function to create a JWT token
def create_jwt_token(username: str) -> str:
    # Payload for the JWT token
    payload = {
        "user": username,
        "exp": datetime.datetime.utcnow() + datetime.timedelta(hours=1), # Token expires in 1 hour
    }
    # Encode the payload with the secret key and HS256 algorithm
    return jwt.encode(payload, SECRET_KEY, algorithm="HS256")