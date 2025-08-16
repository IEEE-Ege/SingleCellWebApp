# utilities.py
import jwt
import os
from datetime import datetime, timedelta
from dotenv import load_dotenv

# load .env file to get the SECRET_KEY
load_dotenv()
SECRET_KEY = os.getenv("SECRET_KEY")

# Function to create a JWT token
def create_jwt_token(username: str) -> str:
    """
    Creates a JWT token.
    Uses the module-level SECRET_KEY.
    """
    if not SECRET_KEY:
        raise ValueError("SECRET_KEY environment variable not set.")
    
    payload = {
        "username": username,  # changed "user" to "username" for clarity
        "exp": datetime.utcnow() + timedelta(hours=1),
    }
    return jwt.encode(payload, SECRET_KEY, algorithm="HS256")

def is_token_expired(token: str, secret_key: str) -> bool:
    """Checks if a JWT token has expired."""
    try:
        payload = jwt.decode(token, secret_key, algorithms=["HS256"])
        exp = payload.get('exp')
        if exp is None:
            return True # Expiration time not found, treat as expired.
        return datetime.utcnow() > datetime.fromtimestamp(exp)
    except jwt.PyJWTError: # Catch specific JWT errors
        return True