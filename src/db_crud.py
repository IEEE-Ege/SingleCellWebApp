from sqlalchemy import select, or_
from sqlalchemy.ext.asyncio import AsyncSession
from db_models import User
import bcrypt

# CREATE
async def create_user(db: AsyncSession, username: str, email: str, password: str):
    hashed = bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()
    user = User(username=username, email=email, password_hash=hashed)
    db.add(user)
    await db.commit()
    await db.refresh(user)
    return user

# READ
async def get_user_by_username(db: AsyncSession, username: str):
    stmt = select(User).where(User.username == username)
    result = await db.execute(stmt)
    return result.scalar_one_or_none()

async def get_user_by_username_or_email(db: AsyncSession, username: str, email: str):
    stmt = select(User).where(or_(User.username == username, User.email == email))
    result = await db.execute(stmt)
    return result.scalar_one_or_none()

async def get_user_by_username_and_email(db: AsyncSession, username: str, email: str):
    stmt = select(User).where(User.username == username, User.email == email)
    result = await db.execute(stmt)
    return result.scalar_one_or_none()

# UPDATE
async def update_user_password(db: AsyncSession, username: str, email: str, new_password: str):
    user = await get_user_by_username_and_email(db, username, email)
    if user:
        if bcrypt.checkpw(new_password.encode(), user.password_hash.encode()):
            return None, "Your new password cannot be the same as your old password. Please choose a different password."
        
        hashed = bcrypt.hashpw(new_password.encode(), bcrypt.gensalt()).decode()
        user.password_hash = hashed
        await db.commit()
        return user, None
    return None, "User not found."
