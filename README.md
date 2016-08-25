# SRA2R
SRA2R, a package to import SRA data directly into R

# Developing on the package

Right now, Docker is needed to work on the package.  This is because the NCBI NGS SDK is challenging to
install.  We have a docker image based on Bioconductor/sequencing_release.

In a terminal:

```sh
docker run -d -p 8787:8787 seandavi/sra2r       ## run docker in background and mirror site running RStudio
docker-machine ip                               ## obtain IP address of Docker machine (Mac/Windows only)
docker exec -ti DOCKER_PROCESS_ID /bin/bash     ## run terminal in Docker machine
                                                ##   only needed if you need command-line access to R
```

RStudio should be running on the docker image, so using the address should connect to it 
(details depend on docker environment).vIn a browser, navigate to http://DOCKER_IP:8787/. 
DOCKER_IP should be replaced with:
- the output of the `docker-machine ip` on windows or Mac
- `localhost` or `127.0.0.1` on linux

Checkout the package using git (or Rstudio) and change the working directory to the 
SRA2R directory (with the DESCRIPTION file in it).  Then, install the package for local development
using devtools.

```R
install.packages('devtools')
devtools::document()
devtools::load_all()
```

# SRA ToolKit examples

## sra-stat example 

http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-stat

```
(known good but large)
sra-stat --quick --xml SRR390728
(smaller file)
sra-stat --quick --xml SRR2971307
(small and no alignment)
sra-stat --quick --xml ERR1162649
```

## sra-pileup

http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-pileup

Command-line

```
sra-pileup -r chr20:1530960-1540960 SRR2971307
```

In R

```
x = getPileUp('SRR390728','chr20',1530960,1540960)
```


# ncbi ngs SDK details

- /usr/include/ngs (interfaces for C++ ngs)
- /usr/include/ncbi-vdb (NGS.hpp)
- /usr/local/share/doc/ngs (javadoc)
- LD_LIBRARY_PATH = /usr/local/ngs/ngs-sdk/lib64:/usr/local/ncbi/ncbi-vdb/lib64:

# R and Rcpp documentation of interest

- http://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
- http://adv-r.had.co.nz/Rcpp.html
- http://r-pkgs.had.co.nz/
- http://statr.me/rcpp-note/index.html
- https://cran.r-project.org/web/packages/Rcpp/vignettes/

