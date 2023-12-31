# This file is part of a set of shiny applications for the analysis of in silico
# FSHD muscle fibres.
# Copyright © 2023 C. R. S. Banerji
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

load("files_for_KM_plot_shiny_not_refined.rd")
library(shiny)
library("tidyverse")
library('deSolve')
library('survival')
library("survminer")
library('ggplot2')
ui <- fluidPage(
  
  # Application title
  titlePanel("Comparing Empirical and Test Parameters"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      actionButton("go", "Go"),
      numericInput("vD",
                   "DUX4 transcription rate:",
                   min = 0.0001,
                   max = 0.1,step=0.001,
                   value = vD),
      numericInput("d0",
                   "DUX4 mRNA degredation rate:",
                   min = 0.01,
                   max = 1,step=0.01,
                   value = d0),
      numericInput("vT",
                   "DUX4 target transcription rate:",
                   min = 1,
                   max = 10,step=1,
                   value = vT),
      numericInput("TD",
                   "Translation rate from DUX4 mRNA to active protien:",
                   min = 0.001,
                   max = 0.5,step=0.01,
                   value = TD),
      numericInput("death_rate",
                   "Death rate of DUX4 Target+ cells:",
                   min = 0.001,
                   max = 0.1,step=0.01,
                   value = death_rate),
      numericInput("diffusion",
                   "Diffusion rate of DUX4 to neighboring cells:",
                   min = 0.001,
                   max = 0.1,step=0.01,
                   value = 0.0402),
      sliderInput("mt_c",
                  "Myotube Circumference:",
                  min = 5,
                  max = 100,
                  value = 9,
                  round=T
      ),
      sliderInput("mt_l",
                  "Myotube Length:",
                  min = 5,
                  max = 100,
                  value = 44,
                  round=T
      ),
      sliderInput("n_iter",
                  "Number of hours to run CA:",
                  min = 3,
                  max = 100,
                  value = 48,
                  round=T),
    ),
    
    # Show a plot of the generated distribution
    # mainPanel(fluidRow(
    #   column(4,
    #          plotOutput("plot1")),
    #   column(4,
    #          plotOutput("plot2")),
    #   column(4,
    #          plotOutput("plot3"))
    # ),
    # fluidRow(
    #   column(4,
    #          plotOutput("plot4")),
    #   column(4,
    #          plotOutput("plot5")),
    #   column(4,
    #          plotOutput("plot6"))
    # )
    # )
    mainPanel(fluidRow(plotOutput('plot1'),plotOutput('plot4')))
  )
)