from sqlalchemy.ext.asyncio import create_async_engine, AsyncSession
from sqlalchemy.orm import sessionmaker
import os
from dotenv import load_dotenv
from db_models import Base

# .env dosyasını yükle
load_dotenv()

# Örn: postgresql+asyncpg://username:password@localhost:5432/mydb
DATABASE_URL = os.getenv("DATABASE_URL")

# Asenkron engine
engine = create_async_engine(DATABASE_URL, echo=False, future=True)

# Asenkron session factory
SessionLocal = sessionmaker(
    bind=engine,
    expire_on_commit=False,
    class_=AsyncSession
)

# Database init (asenkron)
async def init_db():
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)
