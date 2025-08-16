# utilities.py
import jwt
import os
from datetime import datetime, timedelta
from dotenv import load_dotenv

# .env dosyasını yükle

# Function to create a JWT token
def create_jwt_token(username: str, SECRET_KEY) -> str:
    if not SECRET_KEY:
        raise ValueError("SECRET_KEY environment variable not set.")
    
    payload = {
        "user": username,
        "exp": datetime.utcnow() + timedelta(hours=1),  # Token expires in 1 hour
    }
    return jwt.encode(payload, SECRET_KEY, algorithm="HS256")

def is_token_expired(token: str, secret_key: str) -> bool:
    try:
        payload = jwt.decode(token, secret_key, algorithms=["HS256"])
        exp = payload.get('exp')
        if exp is None:
            return True
        return datetime.utcnow() > datetime.fromtimestamp(exp)
    except:
        return True