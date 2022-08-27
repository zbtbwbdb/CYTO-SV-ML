FROM dockerhub.nmdp.org:8443/nmdp-shiny-base:3.6.1


MAINTAINER "Tao Zhang tzhang@nmdp.org"

RUN apt-get update && apt-get install -y libxtst6
RUN apt-get update && apt-get install -y openjdk-8-jdk libssl-dev && apt-get autoremove

RUN R CMD javareconf
RUN install2.r remotes
#RUN R -e 'remotes::install_github("jcheng5/d3scatter")'
RUN R -e 'remotes::install_github("lineupjs/lineup_htmlwidget")'
RUN install2.r DT ggplot2 remotes stringr shinydashboard plotly shinyjs shinyWidgets dplyr Hmisc

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY /cyto-sv-ml /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server/

EXPOSE 80

COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]

