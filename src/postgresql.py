from shiny import App, ui, reactive, render
from htmltools import head_content
from sqlalchemy import create_engine, Integer, String, select, or_
from sqlalchemy.orm import Mapped, mapped_column, Session, DeclarativeBase
import bcrypt
import jwt
import datetime

# SQLAlchemy 2.0 ORM Base
class Base(DeclarativeBase):
    pass

# Define the User model for SQLAlchemy
class User(Base):
    __tablename__ = "users"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, index=True)
    username: Mapped[str] = mapped_column(String, unique=True, index=True, nullable=False)
    email: Mapped[str] = mapped_column(String, unique=True, nullable=False)
    password_hash: Mapped[str] = mapped_column(String, nullable=False)

# Database Configuration
# IMPORTANT: Ensure your PostgreSQL server is running and the 'postgres' user
# with password '1234' has access to the 'AZRA_SCA_DEMO' database.
# Create database tables if they don't exist.
#DATABASE_URL = "postgresql://postgres:1234@localhost:5432/AZRA_SCA_DEMO"
DATABASE_URL = "postgresql+psycopg2://postgres:oyku@localhost/sca"
engine = create_engine(DATABASE_URL, echo=False)
Base.metadata.create_all(engine)

# JWT Secret Key
# IMPORTANT: In a real application, use a strong, randomly generated
# key stored securely (e.g., in environment variables), not hardcoded.
SECRET_KEY = "supersecretkey123"

# Function to create a JWT token
def create_jwt_token(username: str) -> str:
    # Payload for the JWT token
    payload = {
        "user": username,
        "exp": datetime.datetime.utcnow() + datetime.timedelta(hours=1), # Token expires in 1 hour
    }
    # Encode the payload with the secret key and HS256 algorithm
    return jwt.encode(payload, SECRET_KEY, algorithm="HS256")

# Define the User Interface (UI) for the Shiny app
app_ui = ui.page_fluid(
    # Custom CSS styling for the app's appearance
    ui.tags.head(
        ui.tags.link(rel="stylesheet", href="frontend/styles.css"),
        ui.tags.script(src="frontend/app.js"),
    ),


    # Main application div container
    ui.div(
        ui.h2("SCA-Web"), # Application title
        ui.output_ui("message_text"), # Displays messages to the user
        ui.output_ui("main_ui"), # UI for login/register/welcome screen
        ui.hr(), # Horizontal rule
        ui.output_ui("protected_content"), # Content visible only when logged in
        class_="container" # Apply the CSS container class
    ),
)

# Define the Server logic for the Shiny app
def server(input, output, session):
    # Reactive values to manage application state
    logged_in = reactive.Value(False)
    current_user = reactive.Value(None)
    jwt_token = reactive.Value(None)
    message = reactive.Value("")
    message_type = reactive.Value("error")
    # Default page state is "register"
    page_state = reactive.Value("register") # Can be "register", "login", "forgot_password_initiate", "forgot_password_reset"
    show_token = reactive.Value(False)
    
    # Store temporary user info for password reset flow
    reset_username = reactive.Value(None)
    reset_email = reactive.Value(None)

    # Render reactive message text to the UI with dynamic styling
    @output
    @render.ui
    def message_text():
        if message_type() == "success":
            return ui.tags.div(message(), class_="message-success")
        return ui.tags.div(message())

    # Dynamically render the main UI based on login state and page state
    @output
    @render.ui
    def main_ui():
        if logged_in():
            return ui.div(
                ui.h3(f"Welcome, {current_user()}!"),
                ui.p("You are successfully logged in and can now access protected content."),
                ui.input_action_button("btn_toggle_token", "Toggle JWT Token Display", class_="btn"),
                ui.output_ui("token_text"),
                ui.input_action_button("btn_logout", "Log out", class_="btn"),
            )
        elif page_state() == "register":
            return ui.div(
                ui.h3("Register New Account"),
                ui.input_text("reg_username", "Username", placeholder="Enter your desired username"),
                ui.input_text("reg_email", "Email", placeholder="Enter your email address"),
                ui.input_password("reg_password", "Password", placeholder="Create a strong password"),
                ui.input_action_button("btn_register", "Register", class_="btn"),
                ui.a(
                    "Already have an account? Log in here.",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_login', Math.random())",
                    class_="form-link"
                ),
            )
        elif page_state() == "login":
            return ui.div(
                ui.h3("Log In to Your Account"),
                ui.input_text("login_username", "Username", placeholder="Enter your username"),
                ui.input_password("login_password", "Password", placeholder="Enter your password"),
                ui.input_action_button("btn_login", "Log In", class_="btn"),
                ui.a(
                    "Don't have an account? Register here.",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_register', Math.random())",
                    class_="form-link"
                ),
                ui.a(
                    "Forgot Password?",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_forgot_password_initiate', Math.random())",
                    class_="form-link"
                )
            )
        elif page_state() == "forgot_password_initiate":
            return ui.div(
                ui.h3("Reset Your Password"),
                ui.p("Please enter your username and email to proceed."),
                ui.input_text("reset_username_input", "Username", placeholder="Enter your username"),
                ui.input_text("reset_email_input", "Email", placeholder="Enter your email address"),
                ui.input_action_button("btn_reset_password_initiate", "Continue", class_="btn"),
                ui.a(
                    "Back to Login",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_login', Math.random())",
                    class_="form-link"
                )
            )
        elif page_state() == "forgot_password_reset":
            return ui.div(
                ui.h3("Set New Password"),
                ui.p(f"Setting new password for {reset_username()}."),
                ui.input_password("new_password", "New Password", placeholder="Enter your new strong password"),
                ui.input_password("confirm_new_password", "Confirm New Password", placeholder="Confirm your new password"),
                ui.input_action_button("btn_reset_password_final", "Reset Password", class_="btn"),
                ui.a(
                    "Back to Login",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_login', Math.random())",
                    class_="form-link"
                )
            )

    # Render protected content based on JWT token presence
    @output
    @render.ui
    def protected_content():
        if jwt_token():
            return ui.div(
                ui.h4("🔒 Protected Application Dashboard"),
                ui.p("Welcome to your secure dashboard! This area is only accessible after successful authentication."),
                ui.p("You can integrate your secure analytics, personalized reports, or confidential data visualizations here."),
                ui.tags.ul(
                    ui.tags.li("View real-time data analytics."),
                    ui.tags.li("Manage user settings."),
                    ui.tags.li("Access exclusive features."),
                ),
            )
        else:
            return ui.div(
                ui.h4("🚫 Access Restricted"),
                ui.p("Please log in to view the protected content and features."),
                ui.p("Register if you don't have an account yet.")
            )

    # Render the JWT token text with styling for readability
    @output
    @render.ui
    def token_text():
        if show_token() and jwt_token():
           return ui.div(
               jwt_token(),
               class_="token-display"
           )
        return None

    # Reactive effect to toggle the visibility of the JWT token
    @reactive.Effect
    @reactive.event(input.btn_toggle_token)
    def toggle_token():
        show_token.set(not show_token())

    # Reactive effect to switch to the login page
    @reactive.Effect
    @reactive.event(input.go_to_login)
    def go_to_login_event():
        message.set("")
        message_type.set("error")
        page_state.set("login")

    # Reactive effect to switch to the register page
    @reactive.Effect
    @reactive.event(input.go_to_register)
    def go_to_register_event():
        message.set("")
        message_type.set("error")
        page_state.set("register")

    # Reactive effect to switch to the forgot password initiation page
    @reactive.Effect
    @reactive.event(input.go_to_forgot_password_initiate)
    def go_to_forgot_password_initiate_event():
        message.set("")
        message_type.set("error")
        page_state.set("forgot_password_initiate")

    # Reactive effect to handle user registration
    @reactive.Effect
    @reactive.event(input.btn_register)
    def register():
        username = input.reg_username()
        email = input.reg_email()
        password = input.reg_password()

        if not username or not email or not password:
            message_type.set("error")
            message.set("Please fill in all registration fields.")
            return
        
        with Session(engine) as db:
            stmt = select(User).where(or_(User.username == username, User.email == email))
            existing = db.execute(stmt).scalar_one_or_none()

            if existing:
                message_type.set("error")
                message.set("Username or Email already registered. Please log in or use different credentials.")
                page_state.set("login")
                return

            hashed = bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()
            user = User(username=username, email=email, password_hash=hashed)
            db.add(user)
            db.commit()

            message_type.set("success")
            message.set("Registration successful! You can now log in.")
            page_state.set("login")          

    # Reactive effect to handle user login
    @reactive.Effect
    @reactive.event(input.btn_login)
    def login():
        username = input.login_username()
        password = input.login_password()

        if not username or not password:
            message_type.set("error")
            message.set("Please enter both username and password.")
            return

        with Session(engine) as db:
            stmt = select(User).where(User.username == username)
            user = db.execute(stmt).scalar_one_or_none()

        if not user:
            message_type.set("error")
            message.set("User not found or incorrect username.")
            return

        if bcrypt.checkpw(password.encode(), user.password_hash.encode()):
            logged_in.set(True)
            current_user.set(username)
            token = create_jwt_token(username)
            jwt_token.set(token)
            message_type.set("success")
            message.set("Login successful! Welcome.")
        else:
            message_type.set("error")
            message.set("Incorrect password.")

    # Reactive effect to handle initiation of password reset
    @reactive.Effect
    @reactive.event(input.btn_reset_password_initiate)
    def reset_password_initiate():
        username = input.reset_username_input()
        email = input.reset_email_input()

        if not username or not email:
            message_type.set("error")
            message.set("Please provide both username and email.")
            return

        with Session(engine) as db:
            stmt = select(User).where(User.username == username, User.email == email)
            user = db.execute(stmt).scalar_one_or_none()

        if user:
            reset_username.set(username)
            reset_email.set(email)
            message_type.set("success")
            message.set("Username and email verified. Please set your new password.")
            page_state.set("forgot_password_reset")
        else:
            message_type.set("error")
            message.set("Username or email is incorrect. Please try again.")

    # Reactive effect to handle final password reset
    @reactive.Effect
    @reactive.event(input.btn_reset_password_final)
    def reset_password_final():
        new_password = input.new_password()
        confirm_new_password = input.confirm_new_password()

        if not new_password or not confirm_new_password:
            message_type.set("error")
            message.set("Please enter and confirm your new password.")
            return
        
        if new_password != confirm_new_password:
            message_type.set("error")
            message.set("Passwords do not match. Please try again.")
            return

        if not reset_username() or not reset_email():
            message_type.set("error")
            message.set("Error: User information not found for password reset. Please start again.")
            page_state.set("forgot_password_initiate")
            return

        with Session(engine) as db:
            stmt = select(User).where(
                User.username == reset_username(),
                User.email == reset_email()
            )
            user = db.execute(stmt).scalar_one_or_none()

            if user:
                hashed = bcrypt.hashpw(new_password.encode(), bcrypt.gensalt()).decode()
                user.password_hash = hashed
                db.commit()
                message_type.set("success")
                message.set("Your password has been successfully reset. You can now log in with your new password.")
                page_state.set("login")
                # Clear reset state
                reset_username.set(None) 
                reset_email.set(None)            
            if not user:
                message_type.set("error")
                message.set("User not found during password reset. Please try again.")
                page_state.set("forgot_password_initiate")
                return

            user.password_hash = bcrypt.hashpw(new_password.encode(), bcrypt.gensalt()).decode()
            db.add(user)
            db.commit()

        message_type.set("success")
        message.set("Password reset successfully.")
        reset_username.set(None)
        reset_email.set(None)
        page_state.set("login")

    # Reactive effect to handle user logout
    @reactive.Effect
    @reactive.event(input.btn_logout)
    def logout():
        logged_in.set(False)
        current_user.set(None)
        jwt_token.set(None)
        show_token.set(False)
        page_state.set("login")
        message_type.set("success")
        message.set("You have been successfully logged out.")

# Create the Shiny App instance
app = App(app_ui, server)

# This line starts the Shiny web server. 
if __name__ == "__main__":
    app.run()