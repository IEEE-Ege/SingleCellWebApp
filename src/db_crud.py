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
    stmt = select(User).where(or_(User.username == username, User.email == email))
    return db.execute(stmt).scalar_one_or_none()

# UPDATE
def update_user_password(db: Session, username: str, email: str, new_password: str):
    stmt = select(User).where(User.username == username, User.email == email)
    user = db.execute(stmt).scalar_one_or_none()
    if user:
        hashed = bcrypt.hashpw(new_password.encode(), bcrypt.gensalt()).decode()
        user.password_hash = hashed
        db.commit()
    return user

# DELETE