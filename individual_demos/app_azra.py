from shiny import App, ui, reactive, render
from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy.orm import sessionmaker, declarative_base
import bcrypt
import jwt
import datetime

# Database Configuration
# IMPORTANT: Replace 'AZRA_SCA_DEMO' with your actual database name if it's different.
# Ensure the PostgreSQL server is running and the 'postgres' user with password '1234' exists
# and has access to this database.
DATABASE_URL = "postgresql://postgres:1234@localhost:5432/AZRA_SCA_DEMO"
engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(bind=engine)
Base = declarative_base()

# Define the User model for SQLAlchemy
class User(Base):
    __tablename__ = "users"
    id = Column(Integer, primary_key=True, index=True)
    username = Column(String, unique=True, index=True, nullable=False)
    email = Column(String, unique=True, nullable=False)
    password_hash = Column(String, nullable=False)

# Create database tables if they don't exist
Base.metadata.create_all(bind=engine)

# JWT Secret Key - IMPORTANT: In a real application, use a strong, randomly generated
# key stored securely (e.g., in environment variables), not hardcoded.
SECRET_KEY = "supersecretkey123"

# Function to create a JWT token
def create_jwt_token(username):
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
    ui.tags.style("""
        /* General body styling for a softer background and font */
        body {
            background-color: #eef2f7; /* Light blue-grey background */
            font-family: 'Inter', 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
            color: #333;
            line-height: 1.6;
            margin: 0;
            padding: 0;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh; /* Full viewport height */
        }
        /* Container styling for a modern card-like look */
        .container {
            max-width: 450px; /* Slightly wider for better spacing */
            width: 90%; /* Responsive width */
            padding: 30px 35px;
            border-radius: 12px; /* More rounded corners */
            background-color: #ffffff; /* Pure white background */
            box-shadow: 0 8px 25px rgba(0,0,0,0.12); /* Stronger, softer shadow */
            box-sizing: border-box; /* Include padding in element's total width and height */
        }
        /* Title styling */
        h2 {
            text-align: center;
            color: #2c3e50; /* Darker title color */
            margin-bottom: 30px;
            font-size: 2.2rem; /* Larger title */
            font-weight: 700; /* Bolder font */
        }
        /* Form headings styling */
        h3 {
            color: #34495e;
            margin-bottom: 25px; /* More space below heading */
            text-align: center;
            font-size: 1.6rem; /* Larger form headings */
            font-weight: 600;
        }
        /* Input field styling */
        input[type="text"], input[type="password"], input[type="email"] {
            width: 100%;
            padding: 12px 15px; /* More padding */
            margin-bottom: 20px; /* More space between inputs */
            border: 1px solid #dcdcdc; /* Lighter border */
            border-radius: 8px; /* More rounded input fields */
            font-size: 1.05rem; /* Slightly larger text in inputs */
            box-sizing: border-box;
            transition: border-color 0.3s ease, box-shadow 0.3s ease;
        }
        /* Input focus effect */
        input[type="text"]:focus, input[type="password"]:focus, input[type="email"]:focus {
            border-color: #007bff; /* Brighter blue on focus */
            outline: none;
            box-shadow: 0 0 0 3px rgba(0, 123, 255, 0.25); /* Glow effect */
        }
        /* Button styling */
        button, .btn {
            background-color: #007bff; /* Primary blue color */
            color: white;
            padding: 12px 0; /* More vertical padding */
            width: 100%;
            font-size: 1.1rem; /* Larger font for buttons */
            font-weight: 600; /* Bolder button text */
            border: none;
            border-radius: 8px; /* More rounded buttons */
            cursor: pointer;
            transition: background-color 0.3s ease, transform 0.2s ease;
            letter-spacing: 0.5px; /* Slightly spaced letters */
        }
        /* Button hover effect */
        button:hover, .btn:hover {
            background-color: #0056b3; /* Darker blue on hover */
            transform: translateY(-2px); /* Slight lift effect */
        }
        button:active, .btn:active {
            transform: translateY(0); /* Press down effect */
        }
        /* Link styling for switch between forms */
        a.form-link {
            display: block;
            text-align: center;
            margin-top: 20px;
            color: #007bff;
            cursor: pointer;
            text-decoration: none;
            font-size: 0.95rem;
            transition: color 0.3s ease;
        }
        /* Link hover effect */
        a.form-link:hover {
            text-decoration: underline;
            color: #0056b3;
        }
        /* Message area styling (for error/success messages) */
        #message_text {
            text-align: center;
            font-weight: 600;
            margin-top: 15px; /* Adjusted margin */
            margin-bottom: 25px; /* Adjusted margin */
            color: #dc3545; /* Red for error messages */
            font-size: 0.95rem;
        }
        /* Styling for success messages */
        .message-success {
            color: #28a745 !important; /* Green for success */
        }
        /* Horizontal rule styling */
        hr {
            border: 0;
            height: 1px;
            background: #eee; /* Lighter divider */
            margin: 30px 0;
        }
        /* Protected content styling */
        #protected_content {
            background-color: #f8f9fa;
            border: 1px solid #e9ecef;
            border-radius: 8px;
            padding: 20px;
            margin-top: 25px;
            text-align: center;
            color: #495057;
        }
        #protected_content h4 {
            color: #007bff;
            margin-bottom: 15px;
            font-size: 1.3rem;
        }
        /* JWT Token display styling */
        .token-display {
            white-space: pre-wrap;
            word-wrap: break-word;
            max-height: 150px;
            overflow-y: auto;
            background-color: #e9ecef; /* Slightly darker background for token */
            padding: 15px;
            border: 1px solid #ced4da;
            border-radius: 8px;
            font-family: 'Fira Code', 'Cascadia Code', monospace; /* Modern monospace font */
            font-size: 0.85rem;
            margin-top: 15px;
            color: #495057;
        }
    """),
    # Main application div container
    ui.div(
        ui.h2("SCA-Web"), # Application title
        ui.output_ui("message_text"), # Displays messages to the user (CHANGED TO output_ui)
        ui.output_ui("main_ui"), # UI for login/register/welcome screen
        ui.hr(), # Horizontal rule
        ui.output_ui("protected_content"), # Content visible only when logged in
        class_="container" # Apply the CSS container class
    ),
)

# Define the Server logic for the Shiny app
def server(input, output, session):
    # Reactive values to manage application state
    logged_in = reactive.Value(False) # True if a user is logged in
    current_user = reactive.Value(None) # Stores the username of the logged-in user
    jwt_token = reactive.Value(None) # Stores the JWT token
    message = reactive.Value("") # Stores messages to display to the user
    message_type = reactive.Value("error") # "error" or "success" for message styling
    page_state = reactive.Value("register") # Controls whether to show register or login form
    show_token = reactive.Value(False) # Controls visibility of the JWT token

    # Render reactive message text to the UI
    @output
    @render.ui # CHANGED FROM render.text TO render.ui
    def message_text():
        # Apply success class if message_type is "success"
        if message_type() == "success":
            return ui.tags.div(message(), class_="message-success")
        return ui.tags.div(message()) # Returns a div for error messages too

    # Dynamically render the main UI based on login state and page state
    @output
    @render.ui
    def main_ui():
        if logged_in():
           # UI for logged-in users
           return ui.div(
                ui.h3(f"Welcome, {current_user()}!"),
                ui.p("You are successfully logged in and can now access protected content."),
                ui.input_action_button("btn_toggle_token", "Toggle JWT Token Display", class_="btn"),
                ui.output_ui("token_text"), # Displays the JWT token if toggled
                ui.input_action_button("btn_logout", "Log out", class_="btn"),
            )
        elif page_state() == "register":
            # UI for the registration form
            return ui.div(
                ui.h3("Register New Account"),
                ui.input_text("reg_username", "Username", placeholder="Enter your desired username"),
                ui.input_text("reg_email", "Email", placeholder="Enter your email address"),
                ui.input_password("reg_password", "Password", placeholder="Create a strong password"),
                ui.input_action_button("btn_register", "Register", class_="btn"),
                ui.a(
                    "Already have an account? Log in here.",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_login', Math.random())", # Use JS to trigger input event
                    class_="form-link"
                ),
            )
        elif page_state() == "login":
            # UI for the login form
            return ui.div(
                ui.h3("Log In to Your Account"),
                ui.input_text("login_username", "Username", placeholder="Enter your username"),
                ui.input_password("login_password", "Password", placeholder="Enter your password"),
                ui.input_action_button("btn_login", "Log In", class_="btn"),
                ui.a(
                    "Don't have an account? Register here.",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_register', Math.random())", # New JS trigger for going back to register
                    class_="form-link"
                ),
            )

    # Render protected content based on JWT token presence
    @output
    @render.ui
    def protected_content():
        if jwt_token():
            # Content visible when logged in
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
            # Content visible when not logged in
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
                jwt_token(), # The actual JWT token string
                class_="token-display" # Apply the new CSS class for token display
            )
        return None # Hide the token if not toggled or no token exists

    # Reactive effect to toggle the visibility of the JWT token
    @reactive.Effect
    @reactive.event(input.btn_toggle_token)
    def toggle_token():
        show_token.set(not show_token()) # Invert the current state of show_token

    # Reactive effect to switch to the login page
    @reactive.Effect
    @reactive.event(input.go_to_login)
    def go_to_login_event(): # Renamed to avoid confusion with go_to_login() function
        message.set("") # Clear any previous messages
        page_state.set("login") # Change the page state to "login"

    # NEW: Reactive effect to switch to the register page
    @reactive.Effect
    @reactive.event(input.go_to_register)
    def go_to_register_event():
        message.set("") # Clear any previous messages
        page_state.set("register") # Change the page state to "register"


    # Reactive effect to handle user registration
    @reactive.Effect
    @reactive.event(input.btn_register)
    def register():
        username = input.reg_username()
        email = input.reg_email()
        password = input.reg_password()

        # Input validation
        if not username or not email or not password:
            message_type.set("error")
            message.set("Please fill in all registration fields.")
            return

        db = SessionLocal() # Get a new database session
        try:
            # Check if username or email already exists
            if db.query(User).filter((User.username == username) | (User.email == email)).first():
                message_type.set("error")
                message.set("Username or Email already registered. Please log in or use different credentials.")
                page_state.set("login")
                return

            # Hash the password using bcrypt
            hashed = bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()
            # Create a new User object
            new_user = User(username=username, email=email, password_hash=hashed)
            db.add(new_user) # Add the new user to the session
            db.commit() # Commit the changes to the database
            message_type.set("success")
            message.set("Registration successful! You can now log in.")
            page_state.set("login") # Switch to login page after successful registration
        except Exception as e:
            db.rollback() # Rollback in case of error
            message_type.set("error")
            message.set(f"An error occurred during registration: {e}")
        finally:
            db.close() # Close the database session

    # Reactive effect to handle user login
    @reactive.Effect
    @reactive.event(input.btn_login)
    def login():
        username = input.login_username()
        password = input.login_password()

        # Input validation
        if not username or not password:
            message_type.set("error")
            message.set("Please enter both username and password.")
            return

        db = SessionLocal() # Get a new database session
        user = db.query(User).filter(User.username == username).first() # Find user by username
        db.close() # Close the database session

        if not user:
            message_type.set("error")
            message.set("User not found or incorrect username.")
            return

        # Verify password using bcrypt
        if bcrypt.checkpw(password.encode(), user.password_hash.encode()):
            logged_in.set(True) # Set logged_in state to True
            current_user.set(username) # Store current username
            token = create_jwt_token(username) # Create JWT token
            jwt_token.set(token) # Store JWT token
            message_type.set("success")
            message.set("Login successful! Welcome.")
        else:
            message_type.set("error")
            message.set("Incorrect password.")

    # Reactive effect to handle user logout
    @reactive.Effect
    @reactive.event(input.btn_logout)
    def logout():
        # Reset all login-related reactive values
        logged_in.set(False)
        current_user.set(None)
        jwt_token.set(None)
        show_token.set(False)
        page_state.set("login") # Return to login page
        message_type.set("success")
        message.set("You have been successfully logged out.")

# Create the Shiny App instance
app = App(app_ui, server)

# This line starts the Shiny web server. It should be at the end of your script.
if __name__ == "__main__":
    app.run()