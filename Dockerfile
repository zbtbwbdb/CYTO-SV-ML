FROM dockerhub.nmdp.org:8443/nmdp-shiny-base:3.6.1


MAINTAINER "Tao Zhang tzhang@nmdp.org"

RUN apt-get update && apt-get install -y libxtst6
RUN apt-get update && apt-get install -y openjdk-8-jdk libssl-dev && apt-get autoremove

RUN R CMD javareconf

RUN install2.r DT ggplot2 ggforce gridExtra remotes DALEX hash cluster stringr shinydashboard plotly shinyjs shinyWidgets lineupjs dplyr Hmisc
RUN R -e 'remotes::install_github("jcheng5/d3scatter")'

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY /cyto-sv-ml /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server/

EXPOSE 80

COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]

