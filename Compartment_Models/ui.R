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

library(shiny)
library("tidyverse")
library('deSolve')
library('ggplot2')
load("files_for_SEIR_shiny_not_refined.rd")
ui <- fluidPage(
  
  # Application title
  titlePanel("SEIR Models of DUX4 in myocytes and myotubes"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("vD",
                  "DUX4 transcription rate:",
                  min = 0.0001,
                  max = 0.1,
                  value = vD,animate = T),
      sliderInput("d0",
                  "DUX4 mRNA degredation rate:",
                  min = 0.01,
                  max = 1,
                  value = d0,animate = T),
      sliderInput("vT",
                  "DUX4 target transcription rate:",
                  min = 0.1,
                  max = 10,
                  value = vT,animate = T),
      sliderInput("TD",
                  "Translation rate from DUX4 mRNA to active protien:",
                  min = 0.001,
                  max = 0.5,
                  value = TD,animate = T),
      sliderInput("death_rate",
                  "Death rate of DUX4 Target+ cells:",
                  min = 0.001,
                  max = 0.1,
                  value = death_rate,animate = T),
      sliderInput("diffusion",
                  "Diffusion rate of DUX4 - equal mixing of cells:",
                  min = 0.0005,
                  max = 0.001,
                  value = SEIR_diffusion_GA,animate = T),
      sliderInput("n_iter",
                  "Number of hours to run:",
                  min = 3,
                  max = 24*100,
                  value =24*100),
      sliderInput("n_sc",
                  "Number of Starting Single Myocytes:",
                  min = 100,
                  max = 10000,
                  value =5133*(1+per_D_at_D3_sc_adj_trans)),
      sliderInput("n_sn",
                  "Number of Starting Single Myonuclei:",
                  min = 50,
                  max = 5000,
                  value =ceiling(n))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(
        column(6,
               plotOutput("plot1")),
        column(6,
               plotOutput("plot2"))
      ),
      fluidRow(
        column(6,
               plotOutput("plot3")),
        column(6,
               plotOutput("plot4"))
      )
    )
  )
)
