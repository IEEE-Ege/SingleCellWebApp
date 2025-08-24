from sqlalchemy.ext.asyncio import create_async_engine, async_sessionmaker
import os
from dotenv import load_dotenv
from db_models import Base

# Database Configuration
# Create database tables if they don't exist.

load_dotenv()

DATABASE_URL = os.getenv("DATABASE_URL")
engine = create_async_engine(DATABASE_URL, echo=False)

# Async session factory
SessionLocal = async_sessionmaker(bind=engine, expire_on_commit=False)

# Async context manager for session
from contextlib import asynccontextmanager

@asynccontextmanager
async def get_session():
    session = SessionLocal()
    try:
        yield session
    finally:
        await session.close()

async def init_db():
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)