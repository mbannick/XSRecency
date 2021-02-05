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
source("data-generator.R")
source("estimators.R")

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
            HTML("<br/><br/>"),
            sliderInput("bigT", "T*", min=0, max=20, value=2),
            sliderInput("maxT", "tau", min=0, max=20, value=12),
            HTML("Choose the number of simulations to run. Then it will show each line
                 in the top plot for the estimated phi, plus estimates and confidence intervals
                 for each of the recency assay quantities including mu, omega, and beta.
                 The <strong>black line</strong> is the true values based on integrating the true phi curve
                 and the <strong>dashed</strong> line is the mean of the simulations.<br/><br/>"),
            sliderInput("n_sims", "Number of Simulations", min=1, max=50, value=25),
            actionButton("simulate", "Run Asssay Simulation"),
            HTML("<br/><br/>Now you can change the following trial and epidemiological settings
                 without re-running the recency assay simulation."),
            HTML("<br/><br/>"),
            numericInput("inc", "Baseline Incidence", min=0, max=1, value=0.032, step=0.01),
            numericInput("prev", "Prevalence", min=0, max=1, value=0.29, step=0.01),
            radioButtons("type", "Incidence Function", choices=c("constant",
                                                                 "linear",
                                                                 "exponential"),
                         selected="constant", inline=TRUE),
            conditionalPanel(condition="input.type == 'exponential' || input.type == 'linear'",
                             numericInput("rho", "Incidence Decay", value=0.0028)),
            numericInput("n", "Number Screened", value=5000, min=1, max=10000),
            width=6
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("phi.plot"),
           plotOutput("est.plot"),
           plotOutput("inf.plot"),
           plotOutput("res.plot"),
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
            new_text <- paste0("(1 - pgamma(t, ", round(alpha, 4), ", ", round(beta, 4), ")) * (1 - ", input$frr, ") + ", input$frr)
        } else {
            new_text <- paste0("(1 - pgamma(t, ", round(alpha, 4), ", ", round(beta, 4), ")) * (1 - ", input$frr, ") + ", input$frr, " + dnorm(t-7, mean=0, sd=1) / 8")
        }

        updateTextInput(session, "express", value=new_text)

    })

    observeEvent(input$type, {
        if(input$type == "linear"){
            updateNumericInput(session, "rho", value=0.0028)
        }
        if(input$type == "exponential"){
            updateNumericInput(session, "rho", value=0.07)
        }
    })

    # This is what will update the phi simulation
    toListen <- reactive({
        list(input$bigT, input$maxT, input$express, input$n_sims)
    })

    inc.function <- reactive({
        if(input$type == "constant"){
            return(c.incidence)
        } else if(input$type == "linear"){
            return(l.incidence)
        } else if(input$type == "exponential"){
            return(e.incidence)
        }
    })

    inf.function <- reactive({
        if(input$type == "constant"){
            return(c.infections)
        } else if(input$type == "linear"){
            return(l.infections)
        } else if(input$type == "exponential"){
            return(e.infections)
        }
    })

    phi.function <- reactive({
        v$estimators <- NULL
        v$phi <- NULL
        func <- eval(parse(text=paste0("function(t) ", input$express)))
        return(func)
    })

    phi.data <- reactive({
        t <- seq(0, input$maxT, by=0.01)
        phi <- phi.function()(t=t)
        df <- data.frame(t=t, phi=phi)
        return(df)
    })

    inf.data <- reactive({
        e <- runif(1000)
        print(input$rho)
        infect.times <- inf.function()(e, t=0, p=input$prev, lambda_0=input$inc, rho=input$rho)
        t.data <- data.table(times=infect.times)
        return(t.data)
    })

    inc.data <- reactive({
        t <- -1*seq(0, input$maxT, by=1)
        print(inc.function())
        print(input$rho)
        print(t)
        inc <- inc.function()(t=t, lambda_0=input$inc, rho=input$rho)
        print(inc)
        return(data.table(t=t, inc=inc))
    })

    data.data <- reactive({
        df <- generate.data(n=input$n, n_sims=input$n_sims, phi.func=phi.function(),
                            infection.function=inf.function(),
                            baseline_incidence=input$inc, prevalence=input$prev,
                            rho=input$rho)
        print(df)
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

    estimator.data <- reactive({
        trial.data <- data.data()
        rec.data <- v$estimators

        mu <- copy(rec.data[quantity == 'mu'])
        mu.est <- mu$estimate %>% c
        mu.var <- mu$variance %>% c

        beta <- copy(rec.data[quantity == 'beta'])
        beta.est <- beta$estimate %>% c
        beta.var <- beta$variance %>% c

        omega <- copy(rec.data[quantity == 'omega'])
        omega.est <- omega$estimate %>% c
        omega.var <- omega$variance %>% c

        snap <- get.snapshot(
            n_r=trial.data$n_r, n_n=trial.data$n_n, n_p=trial.data$n_p,
            n=trial.data$n, mu=mu.est, mu_var=mu.var
        )
        adj <- get.adjusted(
            n_r=trial.data$n_r, n_n=trial.data$n_n,
            n_p=trial.data$n_p, n=trial.data$n,
            omega=omega.est, omega_var=omega.var,
            beta=beta.est, beta_var=beta.var,
            big_T=input$bigT
        )
        df <- data.table(
            truth=rep(input$inc, input$n_sims),
            snap_est=snap$est,
            snap_var=snap$var,
            adj_est=adj$est,
            adj_var=adj$var
        )
        df[, num := .I]
        print(df)
        melt.est <- data.table::melt(df, id.vars=c("truth", "num"),
                                     measure.vars=patterns("est$", cols=names(df)),
                                     variable.name="estimator", value.name="estimate")
        print(melt.est)
        melt.est[, estimator := lapply(.SD, function(x) gsub("_est$", "", x)), .SDcols="estimator"]
        melt.var <- data.table::melt(df, id.vars=c("truth", "num"),
                                     measure.vars=patterns("var$", cols=names(df)),
                                     variable.name="estimator", value.name="variance")
        melt.var[, estimator := lapply(.SD, function(x) gsub("_var$", "", x)), .SDcols="estimator"]
        melted <- merge(melt.est, melt.var, by=c("truth", "num", "estimator"))
        print(melted)
        melted[, lower := estimate - qnorm(0.975) * variance ** 0.5]
        melted[, upper := estimate + qnorm(0.975) * variance ** 0.5]
        melted[, cover := (lower < truth) & (upper > truth)]
        print(melted)
        return(melted)
    })

    estimator.data.summ <- reactive({
        melted.b <- estimator.data()[, lapply(.SD, mean),
                                     .SDcols="estimate",
                                     by=c("truth", "estimator")]
        return(melted.b)
    })

    output$res.plot <- renderPlot({

        if(is.null(v$estimators)){
            pp <- ggplot()
        } else {
            truth <- unique(estimator.data.summ()$truth)
            total <- nrow(estimator.data())
            all.true <- estimator.data()[cover == TRUE] %>% nrow
            print(all.true)
            if(all.true == total){
                colors <- c("#3669FB")
            } else if(all.true == 0){
                colors <- c("#DF340F")
            } else {
                colors <- c("#DF340F", "#3669FB")
            }

            pp <- ggplot(data=estimator.data(), mapping=aes(x=num, y=estimate, color=cover)) +
                geom_pointrange(aes(ymin=lower, ymax=upper, color=cover), size=0.3, fatten=0.5) +
                geom_hline(yintercept=truth) +
                geom_hline(data=estimator.data.summ(), aes(yintercept=estimate), linetype='dashed',
                           lwd=0.5) +
                theme_minimal() +
                theme(legend.position="none") +
                theme(panel.grid.minor = element_line(size = 0.25),
                      panel.grid.major = element_line(size = 0.25)) +
                scale_color_manual(values=colors) +
                labs("Coverage") + xlab("Simulation") + ylab("Estimate") +
                facet_wrap(~ factor(estimator,
                                    levels=c("snap", "adj")),
                           nrow=1, ncol=2, strip.position="right")
        }
        pp
    })

    output$phi.plot <- renderPlot({
        p <- ggplot() + geom_line(data=phi.data(), mapping=aes(x=t, y=phi)) +
            geom_vline(xintercept=input$bigT, linetype=2) + theme_bw() +
            ggtitle("Phi Function")

        if(!is.null(v$estimators)){
            p <- p + geom_line(data=v$phi, mapping=aes(x=time, y=value, color=sim))
        }
        p
    })

    output$inf.plot <- renderPlot({
        par(mfrow=c(1, 2))
        plot(inc.data()$t, inc.data()$inc, type='l', xlab="-t", ylab="incidence",
             main="Incidence Function")
        hist(inf.data()$times, main="Infection Time Distribution", xlab="-t",
             freq=FALSE, col='white')
        dens <- density(inf.data()$times)
        lines(dens)
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
                facet_wrap(~ quantity, nrow=3, strip.position="right", scales="free")
        }
        p
    })

}

# Run the application
shinyApp(ui = ui, server = server)
