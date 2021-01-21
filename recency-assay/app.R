#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(magrittr)
library(data.table)
library(tidyr)
source("external-study-data.R")
source("external-study-sim.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Recency Assay"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            HTML("First, choose a functional form for phi. When you run the simulation, it uses
                 whatever is currently showing up in this formula box. <br/><br/>"),
            textInput("express", "Phi Formula", value="1 - pgamma(t, 1, 2.7)"),
            HTML("The following choices are some standard characteristic functions.
                 When you select each one, it will allow you to choose a mean window period
                 and a shadow period (and a baseline false recency rate if there is FRR). <br/><br/>"),
            HTML("In order for the phi function to be updated with your choices of the selections,
                 click Update Phi Function. This will calculate the gamma parameters based on your
                 selections. <br/><br/>"),
            selectInput("standard", "Standard Choices",
                        choices=c("No FRR", "Constant FRR", "Non-Constant FRR"),
                        selected="No FRR"),
            conditionalPanel(condition="input.standard == 'Constant FRR' || input.standard == 'Non-Constant FRR'",
                             numericInput("frr", "FRR", value=0.015)),
            conditionalPanel(condition="!is.null(input.standard)",
                             numericInput("window", "Mean Window Period (Days)", value=142)),
            conditionalPanel(condition="!is.null(input.standard)",
                             numericInput("shadow", "Shadow Period (Days)", value=150)),
            actionButton("update_phi", "Update Phi Function"),
            HTML("<br/><br/>For the simulations, you may choose which T* to integrate to for MDRI
                 and which tau to integrate to for mu."),

            sliderInput("bigT", "T*", min=0, max=20, value=2),
            sliderInput("maxT", "tau", min=0, max=20, value=12),
            HTML("Choose the number of simulations to run. Then it will show each line
                 in the top plot for the estimated phi, plus estimates and confidence intervals
                 for each of the recency assay quantities including mu, omega, and beta.
                 The <strong>black line</strong> is the true values based on integrating the true phi curve
                 and the <strong>dashed</strong> line is the mean of the simulations.<br/><br/>"),
            sliderInput("n_sims", "Number of Simulations", min=1, max=25, value=3),
            actionButton("simulate", "Run Simulation"),
            width=6
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("phi.plot"),
           plotOutput("est.plot"),
           width=6
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    observeEvent(input$update_phi, {

        new_text <- ""

        window <- input$window/(365.25)
        shadow <- input$shadow/(365.25)

        alpha <- window ** 2 / (2 * shadow - window ** 2)
        beta <- window / (2 * shadow - window ** 2)

        if(input$standard == "No FRR"){
            new_text <- paste0("1 - pgamma(t, ", round(alpha, 4), ", ", round(beta, 3), ")")
        } else if(input$standard == "Constant FRR"){
            new_text <- paste0("(1 - pgamma(t, ", round(alpha, 4), ", ", round(beta, 4), ")) * (1 - ", input$frr, ") + ", frr)
        } else {
            new_text <- paste0("(1 - pgamma(t, ", round(alpha, 4), ", ", round(beta, 4), ")) * (1 - ", input$frr, ") + ", frr, " + dnorm(t-7, mean=0, sd=1) / 8")
        }

        updateTextInput(session, "express", value=new_text)

    })

    toListen <- reactive({
        list(input$bigT, input$maxT, input$express)
    })

    phi.function <- reactive({
        v$estimators <- NULL
        v$phi <- NULL
        func <- eval(parse(text=paste0("function(t) ", input$express)))
        return(func)
    })

    phi.data <- reactive({
        t <- seq(0, 12, by=0.01)
        phi <- phi.function()(t=t)
        df <- data.frame(t=t, phi=phi)
        return(df)
    })

    v <- reactiveValues(estimators=NULL,
                        phi=NULL)

    observeEvent(input$simulate, {
        sims <- assay.properties.nsim(n_sims=input$n_sims, phi.func=phi.function(),
                                      max_T=input$maxT, big_T=input$bigT)
        v$estimators <- sims$estimators
        v$phi <- sims$phi
    })

    observeEvent(toListen(), {
        v$estimators <- NULL
        v$phi <- NULL
    })

    # simulants <- eventReactive(input$simulate, {
    #
    #     return(sims)
    # })

    output$phi.plot <- renderPlot({
        p <- ggplot() + geom_line(data=phi.data(), mapping=aes(x=t, y=phi)) +
            geom_vline(xintercept=input$bigT, linetype=2) + theme_bw() +
            ggtitle("Phi Function")

        if(!is.null(v$estimators)){
            p <- p + geom_line(data=v$phi, mapping=aes(x=time, y=value, color=sim))
        }
        p
    })

    output$est.plot <- renderPlot({
        if(is.null(v$estimators)){
            p <- ggplot()
        } else {
            phi.function()
            small <- v$estimators[, lapply(.SD, mean), .SDcols=c("estimate", "truth"), by=c("quantity")]

            p <- ggplot(data=v$estimators, mapping=aes(x=sim, y=estimate)) +
                geom_pointrange(aes(ymin=lower, ymax=upper)) +
                geom_hline(data=small, aes(yintercept=truth)) +
                geom_hline(data=small, aes(yintercept=estimate), linetype='dashed',
                           lwd=0.5) +
                theme_minimal() +
                theme(legend.position="right") +
                theme(panel.grid.minor = element_line(size = 0.25),
                      panel.grid.major = element_line(size = 0.25)) +
                # scale_color_brewer(palette="Set1") +
                # labs("Coverage") + xlab("Simulation") + ylab("Estimate") +
                facet_wrap(~ quantity, nrow=3, strip.position="right", scales="free")
        }
        p
    })

}

# Run the application
shinyApp(ui = ui, server = server)
