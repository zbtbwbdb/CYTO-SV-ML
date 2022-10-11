# CYTO-SV-ML
#### CYTO-SV-ML PIPELINE WORKFLOW:
![CYTO-SV-ML PIPELINE WORKFLOW](Workflow.png)


# Run Shiny Application Locally

## Install R Packages

CYTO-SV-ML Shiny Application requires the following R packages to be installed.
```
RUN install2.r remotes
RUN R -e 'install.packages("lineupjs",dependencies=TRUE, repos="https://cran.r-project.org/")'
#RUN R -e 'remotes::install_github("tzhang-nmdp/lineup_htmlwidget")'
RUN install2.r ggforce plotly DT ggplot2 stringr shinydashboard shinyjs shinyWidgets dplyr Hmisc
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

### Online web portal
http://cyto-sv-ml.b12x.org/
