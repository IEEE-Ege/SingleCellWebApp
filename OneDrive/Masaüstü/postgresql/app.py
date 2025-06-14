from shiny import App, ui, reactive, render
import psycopg
import bcrypt
import jwt
import datetime

DB_PARAMS = {
    "host": "localhost",
    "dbname": "sca",
    "user": "postgres",
    "password": "oyku",
}

SECRET_KEY = "supersecretkey123"

def get_db_connection():
    return psycopg.connect(**DB_PARAMS)

def create_jwt_token(username):
    payload = {
        "user": username,
        "exp": datetime.datetime.utcnow() + datetime.timedelta(hours=1),
    }
    return jwt.encode(payload, SECRET_KEY, algorithm="HS256")

app_ui = ui.page_fluid(
    ui.tags.style("""
        /* General container center and max width */
        .container {
            max-width: 400px;
            margin: 40px auto;
            padding: 20px 25px;
            border-radius: 10px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            background-color: #fafafa;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }
        /* Title */
        h2 {
            text-align: center;
            color: #2c3e50;
            margin-bottom: 25px;
        }
        /* Form headings */
        h3 {
            color: #34495e;
            margin-bottom: 15px;
            text-align: center;
        }
        /* Inputs */
        input[type="text"], input[type="password"], input[type="email"] {
            width: 100%;
            padding: 10px 12px;
            margin-bottom: 15px;
            border: 1px solid #ccc;
            border-radius: 6px;
            font-size: 1rem;
            box-sizing: border-box;
            transition: border-color 0.3s ease;
        }
        input[type="text"]:focus, input[type="password"]:focus, input[type="email"]:focus {
            border-color: #3498db;
            outline: none;
        }
        /* Buttons */
        button, .btn {
            background-color: #3498db;
            color: white;
            padding: 10px 0;
            width: 100%;
            font-size: 1rem;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            transition: background-color 0.3s ease;
        }
        button:hover, .btn:hover {
            background-color: #2980b9;
        }
        /* Link */
        a {
            display: block;
            text-align: center;
            margin-top: 12px;
            color: #3498db;
            cursor: pointer;
            text-decoration: none;
            font-size: 0.9rem;
        }
        a:hover {
            text-decoration: underline;
            color: #2980b9;
        }
        /* Message area */
        #message_text {
            text-align: center;
            font-weight: 600;
            margin-bottom: 20px;
            color: #e74c3c; /* red error color */
        }
    """),
    ui.div(
        ui.h2("SCA-Web"),
        ui.output_ui("main_ui"),
        ui.output_text("message_text"),
        ui.hr(),
        ui.output_ui("protected_content"),
        class_="container"
    ),
)

def server(input, output, session):
    logged_in = reactive.Value(False)
    current_user = reactive.Value(None)
    jwt_token = reactive.Value(None)
    message = reactive.Value("")
    page_state = reactive.Value("register")

    @output
    @render.text
    def message_text():
        return message()

    @output
    @render.ui
    def main_ui():
        if logged_in():
           return ui.div(
        ui.h3(f"Welcome, {current_user()}!"),
        ui.p("You are successfully logged in."),
        ui.p("Welcome to the application! We are glad to have you here."),
        ui.input_action_button("btn_toggle_token", "Show token", class_="btn"),
        ui.output_ui("token_text"),
        ui.input_action_button("btn_logout", "Log out", class_="btn"),
    )


        elif page_state() == "register":
            return ui.div(
                ui.h3("Register"),
                ui.input_text("reg_username", "Username"),
                ui.input_text("reg_email", "Email"),
                ui.input_password("reg_password", "Password"),
                ui.input_action_button("btn_register", "Register", class_="btn"),
                ui.hr(),
                ui.a(
                    "Already registered? Click here to log in.",
                    href="#",
                    onclick="Shiny.setInputValue('go_to_login', Math.random())",
                ),
            )
        elif page_state() == "login":
            return ui.div(
                ui.h3("Log In"),
                ui.input_text("login_username", "Username"),
                ui.input_password("login_password", "Password"),
                ui.input_action_button("btn_login", "Log In", class_="btn"),
            )

    @output
    @render.ui
    def protected_content():
        if jwt_token():
            return ui.div(
                ui.h4("Protected Content"),
                ui.p("This page is visible only to logged-in users."),
                ui.p("You can add graphs, dashboards, etc. here."),
            )
        else:
            return ui.div(
                ui.h4("Access Denied"),
                ui.p("Please log in to see this content.")
            )
    
    show_token = reactive.Value(False)

    @output
    @render.ui
    def token_text():
        if show_token() and jwt_token():
           return ui.div(
                jwt_token(),
                style=(
                   "white-space: pre-wrap; "
                   "word-wrap: break-word; "
                   "max-height: 150px; "
                   "overflow-y: auto; "
                   "background-color: #f8f9fa; "
                   "padding: 10px; "
                   "border: 1px solid #ccc; "
                   "border-radius: 6px; "
                   "font-family: monospace; "
                   "margin-top: 10px;"
                )
            )
        return None

    
    @reactive.Effect
    @reactive.event(input.btn_toggle_token)
    def toggle_token():
        show_token.set(not show_token())


    @reactive.Effect
    @reactive.event(input.go_to_login)
    def go_to_login():
        page_state.set("login")

    @reactive.Effect
    @reactive.event(input.btn_register)
    def register_effect():
        username = input.reg_username()
        email = input.reg_email()
        password = input.reg_password()

        if not username or not email or not password:
            message.set("Please fill in all registration fields.")
            return

        hashed = bcrypt.hashpw(password.encode("utf-8"), bcrypt.gensalt())

        try:
            with get_db_connection() as conn:
                with conn.cursor() as cur:
                    cur.execute("SELECT * FROM users WHERE username = %s OR email = %s", (username, email))
                    existing = cur.fetchone()
                    if existing:
                        message.set("You are already registered. Please log in.")
                        page_state.set("login")
                        return

                    cur.execute(
                        "INSERT INTO users (username, email, password_hash) VALUES (%s, %s, %s)",
                        (username, email, hashed.decode()),
                    )
                    conn.commit()

            message.set("Registration successful! Please log in.")
            page_state.set("login")
        except Exception as e:
            message.set(f"An error occurred during registration: {e}")

    @reactive.Effect
    @reactive.event(input.btn_login)
    def login_effect():
        username = input.login_username()
        password = input.login_password()

        if not username or not password:
            message.set("Please enter both username and password.")
            return

        try:
            with get_db_connection() as conn:
                with conn.cursor() as cur:
                    cur.execute("SELECT password_hash FROM users WHERE username = %s", (username,))
                    result = cur.fetchone()
                    if not result:
                        message.set("User not found.")
                        return

                    stored_hash = result[0]
                    if bcrypt.checkpw(password.encode("utf-8"), stored_hash.encode("utf-8")):
                        logged_in.set(True)
                        current_user.set(username)
                        token = create_jwt_token(username)
                        jwt_token.set(token)
                        message.set("Login successful.")
                    else:
                        message.set("Incorrect password.")
        except Exception as e:
            message.set(f"An error occurred during login: {e}")

    @reactive.Effect
    @reactive.event(input.btn_logout)
    def logout_effect():
        logged_in.set(False)
        current_user.set(None)
        jwt_token.set(None)
        page_state.set("login")
        message.set("You have been logged out.")

app = App(app_ui, server)
