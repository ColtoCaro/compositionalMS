FROM jrnold/rstan

WORKDIR /usr/src

#RUN apt-get update

COPY . /usr/src/compositionalMS

# Add Michael Rutter's c2d4u3.5 PPA (and rrutter3.5 for CRAN builds too)
#RUN sudo add-apt-repository -y "ppa:marutter/rrutter3.5"
#RUN sudo add-apt-repository -y "ppa:marutter/c2d4u3.5"
#RUN sudo apt update
#RUN apt install r-cran-rstan
#debian testing user
#RUN apt-get install r-cran-rstan

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-6.tar.gz', repos = NULL, dependencies = TRUE)"
RUN R -e "install.packages('car', repos='http://cran.rstudio.com/', dependencies = TRUE)"


#RUN R CMD INSTALL compositionalMS --no-staged-install
RUN R CMD INSTALL compositionalMS

# finish
RUN echo "compMS installed!"


