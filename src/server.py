from shiny import App, ui, reactive, render
from htmltools import head_content
from sqlalchemy import create_engine, Integer, String, select, or_
from sqlalchemy.orm import Mapped, mapped_column, Session, DeclarativeBase
from utilities import create_jwt_token
import bcrypt
import jwt
import datetime
import os
from dotenv import load_dotenv #pip install python-dotenv
from db import SessionLocal, init_db
from db_crud import create_user, get_user_by_username, get_user_by_username_or_email, update_user_password

# Database Configuration
init_db()

load_dotenv()
SECRET_KEY = os.getenv("SECRET_KEY")

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
                ui.h3(f"Welcome, {current_user().username}!"),
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
        
        db = SessionLocal()
        existing = get_user_by_username_or_email(db, username, email)
        if existing:
            message_type.set("error")
            message.set("Username or Email already registered. Please log in or use different credentials.")
            page_state.set("login")
            return

        create_user(db, username, email, password)
        message_type.set("success")
        message.set("Registration successful! You can now log in.")
        page_state.set("login") 
        
    # Reactive effect to handle user login
    @reactive.Effect #get userdan login statusunu alıp eger loginse print bişe bişe
    @reactive.event(input.btn_login)
    def login():
        username = input.login_username()
        password = input.login_password()

        if not username or not password:
            message_type.set("error")
            message.set("Please enter both username and password.")
            return
        
        db = SessionLocal()
        user = get_user_by_username(db, username)

        if not user or not bcrypt.checkpw(password.encode(), user.password_hash.encode()):
            message_type.set("error")
            message.set("User not found or incorrect username.")
            return

        token = create_jwt_token({"username": user.username}, SECRET_KEY)
        jwt_token.set(token)
        current_user.set(user)
        logged_in.set(True)
        message_type.set("success")
        message.set("Login successful! Welcome!")

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

        db = SessionLocal()
        user = get_user_by_username_or_email(db, username, email)

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

        db = SessionLocal()
        user = update_user_password(db, reset_username(), reset_email(), new_password)
         
        if not user:
                message_type.set("error")
                message.set("User not found during password reset. Please try again.")
                page_state.set("forgot_password_initiate")
                return

        message_type.set("success")
        message.set("Your password has been successfully reset. You can now log in with your new password.")
        page_state.set("login")
        reset_username.set(None)
        reset_email.set(None)

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

    # Protected function: Accessible only to logged-in users
    @reactive.Effect
    def protected_action():
        if logged_in() and jwt_token():
            print("A message only visible to logged-in users.")
        else:
            print("Unauthorized access attempt.")