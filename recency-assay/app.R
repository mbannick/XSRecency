#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Recency Assay"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            textInput("express", "Phi Formula", value="1 - pgamma(t, 1, 2)"),
            sliderInput("bigT", "T*", min=0, max=10, value=2)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("phi.plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$phi.plot <- renderPlot({
        # generate bins based on input$bins from ui.R
        t <- seq(0, 12, by=0.01)
        phi <- eval(parse(text=input$express))

        plot(x=t, y=phi, type='l')
        abline(v=input$bigT, lty='dashed')
    })
}

# Run the application
shinyApp(ui = ui, server = server)
