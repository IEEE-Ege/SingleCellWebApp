# dependencies.py
from db import SessionLocal
from sqlalchemy.orm import Session
import jwt
from db_crud import get_user_by_username
from utilities import SECRET_KEY # import secret key from utilities

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

def get_current_user(token: str, db: Session):
    """
    Decodes the JWT token and fetches the user from the database.
    This acts as a dependency to get the currently authenticated user.
    """
    if not token:
        return None
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=["HS256"])
        username = payload.get("username")
        if username is None:
            return None
    except jwt.PyJWTError:
        return None
    
    user = get_user_by_username(db, username=username)
    return user