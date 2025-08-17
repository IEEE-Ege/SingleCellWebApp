from shiny import App, ui, reactive, render
from htmltools import head_content
import bcrypt
import datetime

# Imports from our own modules
import asyncio
from db import init_db
from db_crud import (
    create_user,
    get_user_by_username,
    get_user_by_username_or_email,
    update_user_password,
    get_user_by_username_and_email
)
from utilities import create_jwt_token, is_token_expired, SECRET_KEY
from dependencies import get_db, get_current_user

# Initialize the database
init_db()

# Define the Server logic for the Shiny app
def server(input, output, session):
    # Reactive values to manage application state
    logged_in = reactive.Value(False)
    current_user = reactive.Value(None)
    jwt_token = reactive.Value(None)
    message = reactive.Value("")
    message_type = reactive.Value("error")
    page_state = reactive.Value("register")
    show_token = reactive.Value(False)
    
    # Store temporary user info for password reset flow
    reset_username = reactive.Value(None)
    reset_email = reactive.Value(None)
    last_activity = reactive.Value(None)

    # --- UI Rendering Functions ---
    @output
    @render.ui
    def message_text():
        if message_type() == "success":
            return ui.tags.div(message(), class_="message-success")
        return ui.tags.div(message(), class_="message-error")

    @output
    @render.ui
    def main_ui():
        if logged_in():
            username = current_user().username if current_user() else "User"
            return ui.div(
                ui.h3(f"Welcome, {username}!"),
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
                ui.a("Already have an account? Log in here.", href="#", onclick="Shiny.setInputValue('go_to_login', Math.random())", class_="form-link"),
            )
        elif page_state() == "login":
            return ui.div(
                ui.h3("Log In to Your Account"),
                ui.input_text("login_username", "Username", placeholder="Enter your username"),
                ui.input_password("login_password", "Password", placeholder="Enter your password"),
                ui.input_action_button("btn_login", "Log In", class_="btn"),
                ui.a("Don't have an account? Register here.", href="#", onclick="Shiny.setInputValue('go_to_register', Math.random())", class_="form-link"),
                ui.a("Forgot Password?", href="#", onclick="Shiny.setInputValue('go_to_forgot_password_initiate', Math.random())", class_="form-link")
            )
        elif page_state() == "forgot_password_initiate":
            return ui.div(
                ui.h3("Reset Your Password"),
                ui.p("Please enter your username and email to proceed."),
                ui.input_text("reset_username_input", "Username", placeholder="Enter your username"),
                ui.input_text("reset_email_input", "Email", placeholder="Enter your email address"),
                ui.input_action_button("btn_reset_password_initiate", "Continue", class_="btn"),
                ui.a("Back to Login", href="#", onclick="Shiny.setInputValue('go_to_login', Math.random())", class_="form-link")
            )
        elif page_state() == "forgot_password_reset":
            return ui.div(
                ui.h3("Set New Password"),
                ui.p(f"Setting new password for {reset_username()}."),
                ui.input_password("new_password", "New Password", placeholder="Enter your new strong password"),
                ui.input_password("confirm_new_password", "Confirm New Password", placeholder="Confirm your new password"),
                ui.input_action_button("btn_reset_password_final", "Reset Password", class_="btn"),
                ui.a("Back to Login", href="#", onclick="Shiny.setInputValue('go_to_login', Math.random())", class_="form-link")
            )

    @output
    @render.ui
    def protected_content():
        if jwt_token():
            return ui.div(
                ui.h4("🔒 Protected Application Dashboard"),
                ui.p("Welcome to your secure dashboard! This area is only accessible after successful authentication."),
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
            )

    @output
    @render.ui
    def token_text():
        if show_token() and jwt_token():
           return ui.div(jwt_token(), class_="token-display")
        return None
    
    # --- Page Transitions ---
    @reactive.Effect
    @reactive.event(input.btn_toggle_token)
    def toggle_token():
        show_token.set(not show_token())

    @reactive.Effect
    @reactive.event(input.go_to_login)
    def go_to_login_event():
        message.set("")
        page_state.set("login")

    @reactive.Effect
    @reactive.event(input.go_to_register)
    def go_to_register_event():
        message.set("")
        page_state.set("register")

    @reactive.Effect
    @reactive.event(input.go_to_forgot_password_initiate)
    def go_to_forgot_password_initiate_event():
        message.set("")
        page_state.set("forgot_password_initiate")

    # --- Main Logic Functions ---
    @reactive.Effect
    @reactive.event(input.btn_register)
    async def register():
        username = input.reg_username()
        email = input.reg_email()
        password = input.reg_password()

        if not all([username, email, password]):
            message.set("Please fill in all registration fields.")
            message_type.set("error")
            return
        
        db_generator = get_db()
        db = await anext(db_generator)
        try:
            existing = await get_user_by_username_or_email(db, username, email)
            if existing:
                message.set("Username or Email already registered.")
                message_type.set("error")
                page_state.set("login")
                return

            await create_user(db, username, email, password)
            message.set("Registration successful! You can now log in.")
            message_type.set("success")
            page_state.set("login")
        finally:
            await db.close()
        
    @reactive.Effect
    @reactive.event(input.btn_login)
    async def login():
        username = input.login_username()
        password = input.login_password()

        if not all([username, password]):
            message.set("Please enter both username and password.")
            message_type.set("error")
            return
        
        db_generator = get_db()
        db = await anext(db_generator)
        try:
            user = await get_user_by_username(db, username)
            if not user or not bcrypt.checkpw(password.encode(), user.password_hash.encode()):
                message.set("Invalid username or password.")
                message_type.set("error")
                return

            token = create_jwt_token(username=user.username)
            jwt_token.set(token)
            current_user.set(user)
            logged_in.set(True)
            last_activity.set(datetime.datetime.now())
            message.set("Login successful! Welcome!")
            message_type.set("success")
        finally:
            await db.close()

    # --- MODIFIED FOR SECURITY ---
    @reactive.Effect
    @reactive.event(input.btn_reset_password_initiate)
    async def reset_password_initiate():
        username = input.reset_username_input()
        email = input.reset_email_input()

        if not all([username, email]):
            message.set("Please provide both username and email.")
            message_type.set("error")
            return

        db_generator = get_db()
        db = await anext(db_generator)
        try:
            user = await get_user_by_username_and_email(db, username, email)
            if user:
                reset_username.set(username)
                reset_email.set(email)
                message.set("User verified. Please set your new password.")
                message_type.set("success")
                page_state.set("forgot_password_reset")
            else:
                message.set("The username and email combination is incorrect.")
                message_type.set("error")
        finally:
            await db.close()
    # ---------------------------

    @reactive.Effect
    @reactive.event(input.btn_reset_password_final)
    async def reset_password_final():
        new_password = input.new_password()
        confirm_new_password = input.confirm_new_password()

        if not new_password or not confirm_new_password:
            message.set("Please enter and confirm your new password.")
            message_type.set("error")
            return
        
        if new_password != confirm_new_password:
            message.set("Passwords do not match.")
            message_type.set("error")
            return

        db_generator = get_db()
        db = await anext(db_generator)
        try:
            user, error = await update_user_password(db, reset_username(), reset_email(), new_password)
            if error:
                message.set(error)
                message_type.set("error")
                return
            if not user:
                message.set("User not found during password reset.")
                message_type.set("error")
                page_state.set("forgot_password_initiate")
                return

            message.set("Password has been successfully reset.")
            message_type.set("success")
            page_state.set("login")
            reset_username.set(None)
            reset_email.set(None)
        finally:
            await db.close()

    @reactive.Effect
    @reactive.event(input.btn_logout)
    def logout():
        logged_in.set(False)
        current_user.set(None)
        jwt_token.set(None)
        show_token.set(False)
        page_state.set("login")
        message.set("You have been successfully logged out.")
        message_type.set("success")

    # --- Background Checks ---
    @reactive.Effect
    def check_session_timeout():
        reactive.invalidate_later(60)
        
        if not logged_in() or not last_activity():
            return
            
        is_expired = (
            (datetime.datetime.now() - last_activity()).total_seconds() > 3600 or
            (jwt_token() and is_token_expired(jwt_token(), SECRET_KEY))
        )
        
        if is_expired:
            message.set("Your session has expired. Please log in again.")
            message_type.set("error")
            logout()

    @reactive.Effect
    async def protected_action():
        if logged_in() and jwt_token():
            db_generator = get_db()
            db = await anext(db_generator)
            try:
                user = await get_current_user(token=jwt_token(), db=db)
                if user:
                    last_activity.set(datetime.datetime.now())
                    print(f"Activity refreshed for user: {user.username}")
                else:
                    print("User from token not found in DB. Logging out.")
                    message.set("Your account could not be verified. Please log in again.")
                    message_type.set("error")
                    logout()
            finally:
                await db.close()