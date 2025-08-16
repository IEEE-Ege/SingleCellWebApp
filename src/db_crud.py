from sqlalchemy import select, or_
from sqlalchemy.orm import Session
from db_models import User
import bcrypt

# CREATE
def create_user(db: Session, username: str, email: str, password: str):
    hashed = bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()
    user = User(username=username, email=email, password_hash=hashed)
    db.add(user)
    db.commit()
    db.refresh(user)
    return user

# READ
def get_user_by_username(db: Session, username: str):
    stmt = select(User).where(User.username == username)
    return db.execute(stmt).scalar_one_or_none()

def get_user_by_username_or_email(db: Session, username: str, email: str):
    """
    Fetches a user if EITHER the username or email exists. 
    Used for registration check.
    """
    stmt = select(User).where(or_(User.username == username, User.email == email))
    return db.execute(stmt).scalar_one_or_none()

# --- NEW SECURE FUNCTION ---
def get_user_by_username_and_email(db: Session, username: str, email: str):
    """
    Fetches a user only if BOTH the username and email match a single record.
    This is used for secure operations like password reset verification.
    """
    stmt = select(User).where(User.username == username, User.email == email)
    return db.execute(stmt).scalar_one_or_none()
# ---------------------------

# UPDATE
def update_user_password(db: Session, username: str, email: str, new_password: str):
    """
    Updates a user's password after secure verification and checks against the old password.
    """
    # Use the new secure function for verification before updating
    user = get_user_by_username_and_email(db, username, email)
    if user:
        # Check if new password is the same as the old one
        if bcrypt.checkpw(new_password.encode(), user.password_hash.encode()):
            return None, "Your new password cannot be the same as your old password. Please choose a different password."
        
        hashed = bcrypt.hashpw(new_password.encode(), bcrypt.gensalt()).decode()
        user.password_hash = hashed
        db.commit()
        return user, None
    return None, "User not found."

