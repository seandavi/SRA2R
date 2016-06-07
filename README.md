# SRA2R
SRA2R, a package to import SRA data directly into R

## The Vignette

http://watson.nci.nih.gov/~sdavis/sra.html

### git setup for github

Before using git with github, you'll want to do this on the AWS instance.

```
ssh-keygen
cat ~/.ssh/id_rsa.pub
```

copy the contents of that string to the github settings/ssh keys location

```
git config --global user.name 'My Name'
git config --global user.email 'myemail@email.com'
```

### Developing on the package

Checkout the package using git (or Rstudio) and change the working directory to the 
SRA2R directory (with the DESCRIPTION file in it).

```
install.packages('devtools')
devtools::document()
devtools::load_all()
```

## SRA ToolKit examples

### sra-stat example 

http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-stat

```
(known good but large)
sra-stat --quick --xml SRR390728
(smaller file)
sra-stat --quick --xml SRR2971307
(small and no alignment)
sra-stat --quick --xml ERR1162649
```

### sra-pileup

http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=sra-pileup

Command-line

```
sra-pileup -r chr20:1530960-1540960 SRR2971307
```

In R

```
x = read.delim(text = system('~/sratoolkit.2.5.7-ubuntu64/bin/sra-pileup -r chr20:1530960-1540960 SRR2971307',intern=TRUE))
```


## ncbi ngs SDK details

- /usr/include/ngs (interfaces for C++ ngs)
- /usr/include/ncbi-vdb (NGS.hpp)
- /usr/local/share/doc/ngs (javadoc)
- LD_LIBRARY_PATH = /usr/local/ngs/ngs-sdk/lib64:/usr/local/ncbi/ncbi-vdb/lib64:


## R and Rcpp documentation of interest

- http://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
- http://adv-r.had.co.nz/Rcpp.html
- http://r-pkgs.had.co.nz/
- http://statr.me/rcpp-note/index.html
- https://cran.r-project.org/web/packages/Rcpp/vignettes/

## Link to presentation to a small example of using read fastq from SRA

- http://ec2-52-90-90-120.compute-1.amazonaws.com/presentation/view#/

## Setting up Docker on local machine
In terminal:
```
docker run -d -p 8787:8787 seandavi/sra2r       ## run docker in background and mirror site running RStudio

docker-machine ip                               ## obtain IP address of Docker machine

docker exec -ti e264 /bin/bash                  ## run terminal in Docker machine
```
In browser:
```
192.168.99.100:8787                             ## navigate to RStudio running on Docker
```