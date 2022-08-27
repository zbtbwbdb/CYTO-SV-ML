# CYTO-SV-ML
#### CYTO-SV-ML PIPELINE WORKFLOW:
![CYTO-SV-ML PIPELINE WORKFLOW](Workflow.png)


# Run Shiny Application Locally

## Install R Packages

EFS Calculator Shiny Application requires the following R packages to be installed.
```
RUN install2.r devtools DT ggplot2 ggforce gridExtra remotes DALEX hash cluster stringr shinydashboard plotly shinyWidgets lineupjs shinyjs
RUN R -e 'remotes::install_github("jcheng5/d3scatter")'
```

# Run App

### Change directory to `cyto-sv-ml`
```
cd cyto-sv-ml
```

### Start Shiny Application from R Console

```
shiny::runApp()
```
