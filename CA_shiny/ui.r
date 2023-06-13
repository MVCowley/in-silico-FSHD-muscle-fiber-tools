library(shiny)
library("tidyverse")
library('deSolve')
library(lattice)
library('ggplot2')
load("files_for_CA_MT_size_shiny_not_refined.rd")
ui <- fluidPage(
  
  # Application title
  titlePanel("Cellular Automaton Model of DUX4 in a Myotube"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      actionButton("go", "Run Automaton"),
      sliderInput("vD",
                  "DUX4 transcription rate:",
                  min = 0.0001,
                  max = 0.1,
                  value = vD),
      sliderInput("d0",
                  "DUX4 mRNA degredation rate:",
                  min = 0.01,
                  max = 1,
                  value = d0),
      sliderInput("vT",
                  "DUX4 target transcription rate:",
                  min = 0.1,
                  max = 10,
                  value = vT),
      sliderInput("TD",
                  "Translation rate from DUX4 mRNA to active protien:",
                  min = 0.001,
                  max = 0.5,
                  value = TD),
      sliderInput("death_rate",
                  "Death rate of DUX4 Target+ cells:",
                  min = 0.001,
                  max = 0.1,
                  value = death_rate),
      sliderInput("diffusion",
                  "Diffusion rate of DUX4 to neighboring cells:",
                  min = 0.001,
                  max = 0.1,
                  value = 0.04023596),
      sliderInput("mt_c",
                  "Myotube Circumference:",
                  min = 5,
                  max = 100,
                  value = 9,
                  round = T),
      sliderInput("mt_l",
                  "Myotube Length:",
                  min = 5,
                  max = 100,
                  value = 44,
                  round = T),
      sliderInput("n_iter",
                  "Number of hours to run:",
                  min = 3,
                  max = 100,
                  value = n_iter,
                  round = T),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plot1"),plotOutput("plot2")
    )
  )
)