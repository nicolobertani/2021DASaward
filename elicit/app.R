suppressMessages({
  library(shiny)
  library(purrr)
  library(linprog)
  library(quadprog)
  library(MCMCpack)
  setwd("~/OneDrive - ucp.pt/Shared Working Directory - ED, NB/ELICIT EXPERIMENT/R code and deliverables")
  source("../Helper functions/M-spline and I-spline.R")
  source("../Helper functions/Helper functions for data analysis.R")
})

cols <- c(
  rgb(.8,0,0),
  rgb(26/255,48/255,132/255),
  rgb(0,.6,0))

helper <- seq(0, 1, .01)
chosen.xi <- c(.1, .9)
order <- 3
m <- length(M_knots_sequence(order, chosen.xi)) - order

sim.answers <- data.frame(matrix(NA, nrow = 1, ncol = 4))
colnames(sim.answers) <- c('p.x', 'z', 'w.p', 's')
x <- 120
y <- 10
epsilon <- Inf
epsilon_threshold <- .1
iteration <- 1
# keep.questioning <- T
A1 <- matrix(rep(rep(0, m), m), ncol = m)
diag(A1) <- 1
lower_bound <- I_spline(helper, order, interior_knots = chosen.xi, individual = T)[, m]
upper_bound <- I_spline(helper, order, interior_knots = chosen.xi, individual = T)[, 1]

sim.answers[iteration, 1] <- .9
sim.answers[iteration, 2] <- 65
sim.answers[iteration, 3] <- ((sim.answers[iteration, 2] - y) / (x - y))

choice.text <- c(
  paste0('winning €', sim.answers[iteration, 2], ' for sure.'),
  paste0('A ', (1 - sim.answers[iteration, 1]) * 100, '% chance of winning €', y, ' and a ', sim.answers[iteration, 1] * 100, '% chance of winning €', x, ".")
)
choice.values <- c(1, 0)
shuffle.order <- sample(1:2)
choice.text <- choice.text[shuffle.order]
choice.values <- choice.values[shuffle.order]


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  fluidRow(
    
    column(5, 
           titlePanel("Experimental input:"),
           
           uiOutput(outputId = 'ui.buttons'),
           
           uiOutput(outputId = 'ui.confirm'),
           
           uiOutput(outputId = 'ui.done')
    ),
    
    column(7, 
           titlePanel("Workings of the procedure:"),
           
           plotOutput(outputId = 'bounds')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  choice <- eventReactive(input$confirm, input$choice)
  rv <- reactiveValues(keep.questioning = T)
  
  output$ui.buttons <- renderUI({
    if (!rv$keep.questioning) return()
    inputPanel(
      radioButtons(inputId = 'choice',
                   choiceValues = choice.values,
                   choiceNames = choice.text,
                   selected = "",
                   label = 'What do you prefer?',
                   width = 'auto')
    )
  })

  output$ui.confirm <- renderUI({
    if (!rv$keep.questioning) return()
    actionButton(inputId = 'confirm',
                 'Confirm my choice')
  })
  
  output$bounds <- renderPlot({
    
    if (input$confirm == 0) {
      separator <- 8.3/10
      par(fig = c(0, 1, separator, 1), new = F,
          mar = c(0,0,0,0), oma = c(0,0,0,0), xpd = T)
      plot(0:1, 0:1, type = 'n', axes = F, xlab = NA, ylab = NA)
      legend(x = .2, y = 1,
             legend = c("Questions:", "Chose sure", "Chose lottery", 'Next question'),
             pch = c(NA, 19, 19, 19),
             col = c(NA, cols[2:1], gray(.5)),
             text.font = c(1, 3, 3, 3),
             cex = 1, box.lty = "blank", y.intersp = 1
      )
      legend(x = .5, y = 1,
             legend = c("Lines - bounds:", "Upper bound", "Lower bound"),
             lty = c(1, 1),
             lwd = c(2, 2),
             col = c(NA, cols[2:1]),
             text.font = c(1, 3, 3),
             cex = 1, box.lty = "blank", y.intersp = 1
      )
      
      par(fig = c(0, 1, 0, separator), new = T, 
          oma = c(0, 0, 0, 0), mar = c(2, 4, 1, 1) + .1, xaxs = "i", yaxs = "i", pty = "s")
      plot(0:1, 0:1, type = 'n', axes = F, xlab = NA, ylab = NA)
      par(fig = c(0, 1, 0, separator), new = T, 
          oma = c(0, 0, 0, 0), mar = c(3, 3, 0.6, 0) + .1, xaxs = "i", yaxs = "i", pty = "s")
      s()
      # add next question
      points(sim.answers[iteration, ]$p.x, sim.answers[iteration, ]$w.p, col = gray(.5), pch = 19)
      ## add bounds
      polygon(x = c(1, helper),
              y = c(0, lower_bound), 
              border = NA, col = ggplot2::alpha(cols[1], .2))
      polygon(x = c(0, helper),
              y = c(1, upper_bound), 
              border = NA, col = ggplot2::alpha(cols[2], .2))
      
      points(helper, lower_bound, type = 'l', 
             col = ggplot2::alpha(cols[1], .8), lwd = 2, lty = 1)
      points(helper, upper_bound, type = 'l', 
             col = ggplot2::alpha(cols[2], .8), lwd = 2, lty = 1)
      # axis labels 
      mtext(expression(p), 1, 2, cex = 1.2)
      mtext(expression(w(p)), 2, 2, cex = 1.2)
      
      
    } else {
      
      sim.answers[iteration, 4] <<- as.numeric(choice())
      iteration <<- iteration + 1
      
      ## update upper and lower bounds  
      input_lower <- sim.answers[sim.answers$s == 0, ]
      input_upper <- sim.answers[sim.answers$s == 1, ]
      
      ## both positive number of constraints
      if(nrow(input_lower) != 0 && nrow(input_upper) != 0) {
        A2 <- t(rbind(I_spline(input_lower$p.x, order, interior_knots = chosen.xi, individual = T),
                      I_spline(input_upper$p.x, order, interior_knots = chosen.xi, individual = T)))
        
        b <- c(1,
               rep(0, m),
               input_lower$w.p,
               input_upper$w.p
        )
        
        constraint_signs <- c("==", 
                              rep(">=", m),
                              rep(">=", length(input_lower$w.p)),
                              rep("<=", length(input_upper$w.p)))
      }
      
      ## if no lower bound constraints
      if(nrow(input_lower) == 0 && nrow(input_upper) != 0) {
        A2 <- t(I_spline(input_upper$p.x, order, interior_knots = chosen.xi, individual = T))
        
        b <- c(1,
               rep(0, m),
               input_upper$w.p
        )
        
        constraint_signs <- c("==", 
                              rep(">=", m),
                              rep("<=", length(input_upper$w.p)))
        
        lower_bound <- I_spline(seq(.01, .99, by = .01), order, interior_knots = chosen.xi, individual = T)[, m]
      }
      
      ## if no upper bound constraints
      if(nrow(input_lower) != 0 && nrow(input_upper) == 0) {
        A2 <- t(I_spline(input_lower$p.x, order, interior_knots = chosen.xi, individual = T))
        
        b <- c(1,
               rep(0, m),
               input_lower$w.p
        )
        
        constraint_signs <- c("==", 
                              rep(">=", m),
                              rep(">=", length(input_lower$w.p)))
        
        upper_bound <- I_spline(seq(.01, .99, by = .01), order, interior_knots = chosen.xi, individual = T)[, 1]
      }
      
      if(1 %in% dim(A2)) A2 <- t(A2)
      A <- t(cbind(
        rep(1, m),
        A1,
        A2
      ))
      
      ## calculate lower bound
      lower_bound <- sapply(helper, function(local.x) {
        c <- I_spline(local.x, 3, interior_knots = chosen.xi, individual = T)
        sol <- solveLP(c, b, A,
                       maximum = F,
                       const.dir = constraint_signs,
                       lpSolve = T)
        sol$solution %*% c
      })
      
      ## calculate upper bound
      upper_bound <- sapply(helper, function(local.x) {
        c <- I_spline(local.x, order, interior_knots = chosen.xi, individual = T)
        sol <- solveLP(c, b, A, 
                       maximum = T,
                       const.dir = constraint_signs,
                       lpSolve = T)
        sol$solution %*% c
      })
      
      ## calculate max difference based on updated bounds
      D <- upper_bound - lower_bound
      if (max(D) < .1) {
        rv$keep.questioning <- F
      }
      new.w.p <- ((upper_bound + lower_bound)[D == max(D)] / 2)[1]
      if (length(D[D == max(D)]) == 1) {
        sim.answers[iteration, 1] <<- helper[D == max(D)]
      } else {
        warning('multiple optimal bisection points')      
        sim.answers[iteration, 1] <<- helper[D == max(D)][which.max(abs(helper[D == max(D)] - .5))]
      }
      sim.answers[iteration, 2] <<- (new.w.p) * (x - y) + y
      sim.answers[iteration, 3] <<- new.w.p
      
      choice.text <- c(
        paste0('winning €', round(sim.answers[iteration, 2]), ' for sure.'),
        paste0('A ', round((1 - sim.answers[iteration, 1]) * 100), '% chance of winning €', y, ' and a ', round(sim.answers[iteration, 1] * 100), '% chance of winning €', x, ".")
      )
      choice.values <- c(1, 0)
      shuffle.order <- sample(1:2)
      choice.text <- choice.text[shuffle.order]
      choice.values <- choice.values[shuffle.order]
      updateRadioButtons(
        session,
        'choice',
        choiceValues = choice.values,
        choiceNames = choice.text,
        selected = ""
      )
      
      # plot ------
      separator <- 8.3/10
      par(fig = c(0, 1, separator, 1), new = F,
          mar = c(0,0,0,0), oma = c(0,0,0,0), xpd = T)
      plot(0:1, 0:1, type = 'n', axes = F, xlab = NA, ylab = NA)
      legend(x = .2, y = 1,
             legend = c("Questions:", "Chose sure", "Chose lottery", 'Next question'),
             pch = c(NA, 19, 19, 19),
             col = c(NA, cols[2:1], gray(.5)),
             text.font = c(1, 3, 3, 3),
             cex = 1, box.lty = "blank", y.intersp = 1
      )
      legend(x = .5, y = 1,
             legend = c("Lines - bounds:", "Upper bound", "Lower bound"),
             lty = c(1, 1),
             lwd = c(2, 2),
             col = c(NA, cols[2:1]),
             text.font = c(1, 3, 3),
             cex = 1, box.lty = "blank", y.intersp = 1
      )
      
      par(fig = c(0, 1, 0, separator), new = T, 
          oma = c(0, 0, 0, 0), mar = c(2, 4, 1, 1) + .1, xaxs = "i", yaxs = "i", pty = "s")
      plot(0:1, 0:1, type = 'n', axes = F, xlab = NA, ylab = NA)
      par(fig = c(0, 1, 0, separator), new = T, 
          oma = c(0, 0, 0, 0), mar = c(3, 3, 0.6, 0) + .1, xaxs = "i", yaxs = "i", pty = "s")
      s()
      points(sim.answers$p.x, sim.answers$w.p,
             xlim = c(0,1), ylim = c(0,1), pch = 19,
             col = cols[1 + as.numeric(sim.answers$s)])
      # add next question
      if (rv$keep.questioning) {
        points(sim.answers[iteration, ]$p.x, sim.answers[iteration, ]$w.p, col = gray(.5), pch = 19)
      }
      ## add bounds
      polygon(x = c(1, helper),
              y = c(0, lower_bound), 
              border = NA, col = ggplot2::alpha(cols[1], .2))
      polygon(x = c(0, helper),
              y = c(1, upper_bound), 
              border = NA, col = ggplot2::alpha(cols[2], .2))
      
      points(helper, lower_bound, type = 'l', 
             col = ggplot2::alpha(cols[1], .8), lwd = 2, lty = 1)
      points(helper, upper_bound, type = 'l', 
             col = ggplot2::alpha(cols[2], .8), lwd = 2, lty = 1)
      # axis labels 
      mtext(expression(p), 1, 2, cex = 1.2)
      mtext(expression(w(p)), 2, 2, cex = 1.2)
      
    }
    
  })
  
  output$ui.done <- renderUI({
    if (rv$keep.questioning) return()
    actionButton("done", "Done! Select one shape.")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
