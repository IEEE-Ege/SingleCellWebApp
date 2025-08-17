from shiny import App, ui, reactive, render
from htmltools import head_content
from server import server



# Define the User Interface (UI) for the Shiny app
app_ui = ui.page_fluid(
    # Custom CSS styling for the app's appearance
    ui.head_content(ui.include_css("www/frontend/styles.css"),
                    ui.tags.script(src="www/frontend/app.js")),


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

# Create the Shiny App instance
app = App(app_ui, server)

# This line starts the Shiny web server. 
if __name__ == "__main__":
    app.run()