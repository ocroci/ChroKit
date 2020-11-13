### find and, if necessary, install, required packages for Chrokit
x=rownames(installed.packages())
#bioCversion will define the version for installing BiocManager packages
#change this if you have an R version that is not supported. For example, 
#bioCversion="3.8" for R 3.5, bioCversion="3.9" for R 3.6 and so on...
bioCversion="3.9"
cranRepo='http://cran.us.r-project.org'
## IMPORTANT: on linux (tested on Ubuntu), additional libraries are rquired:

## libcurl4-openssl-dev is required for RCurl package; install with:
# sudo apt install libcurl4-openssl-dev
## libxml2-dev is required for XML package; install with:
# sudo apt install libxml2-dev
## libssl-dev is required for openssl package; install with:
# sudo apt install libssl-dev
## libmysqlclient-dev is required for RMySQL package; install with:
# sudo apt install libmysqlclient-dev (or default-libmysqlclient-dev)
## for zlib.h (package data.table) install also the following library:
# sudo apt-get install libz-dev

## If you find problems in installing the above dependencies, try also with:
# sudo apt install libgcrypt11-dev (libgcrypt20-dev or other versions) libgnutls-dev (libgnutls28-dev or newer) librtmp-dev



## to correct undefined symbol: __atomic_fetch_add_8 error (tested on raspbian buster, raspberry pi 4)

## add this line to .bashrc file in home directory
# export LD_PRELOAD=/usr/lib/arm-linux-gnueabihf/libatomic.so.1.2.0
## ...and source the file with:
# source .bashrc


#shiny
if(! ("shiny" %in% x)){
	print("Installing shiny package...")
	install.packages('shiny', repos=cranRepo)
}else{
	print("shiny package already installed...")
}

#shinyFiles
if(! ("shinyFiles" %in% x)){
	print("Installing shinyFiles package...")
	install.packages('shinyFiles', repos=cranRepo)
}else{
	print("shinyFiles package already installed...")
}

#shinydashboard
if(! ("shinydashboard" %in% x)){
	print("Installing shinydashboard package...")
	install.packages('shinydashboard', repos=cranRepo)
}else{
	print("shinydashboard package already installed...")
}

#shinyWidgets
if(! ("shinyWidgets" %in% x)){
	print("Installing shinyWidgets package...")
	install.packages('shinyWidgets', repos=cranRepo)
}else{
	print("shinyWidgets package already installed...")
}

#fastcluster
if(! ("fastcluster" %in% x)){
	print("Installing fastcluster package...")
	install.packages('fastcluster', repos=cranRepo)
}else{
	print("fastcluster package already installed...")
}


#VennDiagram
if(! ("VennDiagram" %in% x)){
	print("Installing VennDiagram package...")
	install.packages('VennDiagram', repos=cranRepo)
}else{
	print("VennDiagram package already installed...")
}

#data.table
if(! ("data.table" %in% x)){
	print("Installing data.table package...")
	install.packages('data.table', repos=cranRepo)
}else{
	print("data.table package already installed...")
}


#RColorBrewer
if(! ("RColorBrewer" %in% x)){
	print("Installing RColorBrewer package...")
	install.packages('RColorBrewer', repos=cranRepo)
}else{
	print("RColorBrewer package already installed...")
}


#ppcor
if(! ("ppcor" %in% x)){
	print("Installing ppcor package...")
	install.packages('ppcor', repos=cranRepo)
}else{
	print("ppcor package already installed...")
}

#inline
if(! ("inline" %in% x)){
	print("Installing inline package...")
	install.packages('inline', repos=cranRepo)
}else{
	print("inline package already installed...")
}


###determine the version of R. if R>=3.5, install BiocManager package
###and install packages from bioconductor from this package
Rversion=strsplit(version$version.string,split=" ")[[1]][3]
Rversion=strsplit(Rversion,split="\\.")[[1]]
Rversion_main=as.numeric(Rversion[1])
Rversion_submain=as.numeric(Rversion[2])
if(Rversion_main==3){
	if(Rversion_submain>=5){
		R35=TRUE
	}else{
		R35=FALSE
	}
}else if(Rversion_main>3){
	R35=TRUE
}else{
	R35=FALSE
}

#BiocManager if R version> 3.5
if(! ("BiocManager" %in% x) & R35){
	print("Installing BiocManager package for R > 3.5.0 ...")
	install.packages("BiocManager")
	BiocManager::install(version=bioCversion)
}else{
	print("BiocManager package already installed...")
}



#GenomicRanges
if(! ("GenomicRanges" %in% x)){
	print("Installing GenomicRanges package...")
	install.packages('RCurl', repos=cranRepo)
	if(R35){
		BiocManager::install("GenomicRanges", version = bioCversion,ask=FALSE)
	}else{
		source("http://bioconductor.org/biocLite.R"); biocLite("GenomicRanges",ask=FALSE)
	}
}else{
	print("GenomicRanges package already installed...")
}


#Rsamtools
if(! ("Rsamtools" %in% x)){
	print("Installing Rsamtools package...")
	if(R35){
		BiocManager::install("Rsamtools", version = bioCversion,ask=FALSE)
	}else{
		source("http://bioconductor.org/biocLite.R"); biocLite("Rsamtools",ask=FALSE)
	}
	
}else{
	print("Rsamtools package already installed...")
}



#rtracklayer
if(! ("rtracklayer" %in% x)){
	print("Installing rtracklayer package...")
	if(R35){
		BiocManager::install("rtracklayer", version = bioCversion,ask=FALSE)
	}else{
		source("http://bioconductor.org/biocLite.R"); biocLite("rtracklayer",ask=FALSE)
	}
	
}else{
	print("rtracklayer package already installed...")
}




#################################################################
# with CONDA
#################################################################

### create "environment" in which install all the dependencies in conda
### maybe R version will be downgraded

# conda create --name environment 
# source activate environment
# conda install -c r r-base
# conda install -c conda-forge r-inline r-shinyfiles r-shinyDashboard r-fastcluster r-shinywidgets r-ppcor r-lattice
# conda install -c bioconda bioconductor-genomicranges bioconductor-rsamtools bioconductor-rtracklayer bioconductor-graph bioconductor-rbgl r-venndiagram
# conda install r-rcolorbrewer r-gtools r-reshape r-data.table r-RCurl r-xml

#rtracklayer may not work (shared objects not found). To solve this, use BiocManager if R>3.5 from within R





