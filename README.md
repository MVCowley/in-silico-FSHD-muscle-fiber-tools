# in-silico-FSHD-muscle-fiber-tools

Shiny tools implementing the model from the eLife article: An in silico FSHD muscle fibre for modelling DUX4 dynamics and predicting the impact of therapy (https://elifesciences.org/articles/88345).
These tools allow the pre-screening of theraputic stratergies via the modification of model parameters.

Web applications of the tools are hosted at:

1. Compartment Models: https://crsbanerji.shinyapps.io/compartment_models/  
2. Cellular Automaton: https://crsbanerji.shinyapps.io/ca_shiny/ 
3. Survival Analysis: https://crsbanerji.shinyapps.io/survival_sim/ 

Alternatively, clone this repository and open an R terminal in the directory of the shiny app of interest.
Run `ui.r` and `server.r` to create the objects `ui` and `server`. Then run `shinyApp(ui, server)` to open a local server.
Use `install.packages()` to install any required packages.
Tested on R version 4.2.3.
