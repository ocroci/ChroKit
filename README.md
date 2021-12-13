# ChroKit
ChroKit (The **Chro**matin Tool**Kit**) is a Shiny-based framework to analyze and visualize interactively genomic data.

This bioinformatics tool can help the researchers to process the data from next generation sequencing (NGS) experiments, such as Chip-Seq, ATAC-Seq and any other NGS experiments aimed at analyzing the enrichment of specific genomic features in particular regions of interest. 

As input, ChroKit takes aligned reads (BAM or WIG files), genomic ranges (BED or GTF/GFF formats) and list of genes and performs a several operations on them, such as the refinement of ranges boundaries, extraction of sequence patterns, gene ontologies on the annotated genes and the calculation of reads enrichment. The user can then perform logical operation on the genomic regions and unsupervised clustering to create further subsets of regions.
The wide variety of interactive plots offered by ChroKit can be modified simply by playing with mouse cursor and downloaded as pdf files. Working sessions can be exported as RDS files, allowing data sharing and reproducibility of the analyses.


<img src="https://github.com/ocroci/ChroKit/blob/master/logo2.png" height="50%" width="50%">

# Run the program

## From Docker image (recommended)
Docker image is available with pre-installed libraries of human and mouse genome assemblies (https://hub.docker.com/r/ocroci/chrokit). Be sure to have Docker installed on your system.

To pull the image, use the following command from the terminal:\
```sudo docker pull ocroci/chrokit:HumanMouseGenomes```

To run the program, type this command from the terminal:\
```sudo docker run -v /home/:/mnt/ -p <port>:6060 -it ocroci/chrokit:HumanMouseGenomes```

In this case the /mnt folder inside the container will mount the /home folder of the host system; change these folders according to your needs.
\<port\> is an arbitrary port on the host system to use for accessing the docker image; this port must be free.

To use the application, open a web browser and
  - if you are using a personal computer, go to:\
    ```127.0.0.1:<port>```
  - if you are using a remote machine, go to:\
    ```<IP>:<port> ```\
    where \<IP\> is the IP address of the remote machine in which the Docker container is running and the \<port\> is the port selected when running the image.

For further instructions, go to https://hub.docker.com/r/ocroci/chrokit


## From source (using R interpreter)

### Install dependencies 
- Download and install the R interpreter (suggested version 3.5 or higher) on your computer or on a remote machine
- Download the ChroKit source code in this gitHub page. Unzip the folder if necessary 
- Make sure the OS-specific requirements are satisfied:

  #### Linux users
  For Linux users, make sure to install the required system packages; this can be done from a terminal with the following command:\
  ```sudo apt install libcurl4-openssl-dev libxml2-dev libssl-dev libz-dev```\
  In case of further errors, try to follow the suggestions at the beginning of the *installChrokitDependencies.R* script.


  #### MacOS users
  To run ChroKit from within R in a MacOS machine, make sure the Xcode command line tools are properly installed and updated.
  This could be done by simply typing the following command on a terminal:
  ```xcode-select --install```\


  #### Windows users
  Make sure to install the R interpreter in a directory path without spaces: when prompted the path for installation choose C:\R\ as the path.
  After installing the R interpreter, install Rtools; then, modify the PATH variable to include also all the binaries of Rtools.
  To install Rtools and modify the PATH variable, follow the instructions at the link: https://cran.r-project.org/bin/windows/Rtools/

  Note for windows users: only BAM file association is allowed (WIG files not supported); moreover, only one core is allowed, due to the use of “parallel”  library, which works only on UNIX operating systems.


- Open the R interpreter and go to the main source directory of the program; type:\
 ``` > setwd("/path/to/the/ChroKit/folder")```\
  where "/path/to/the/ChroKit/folder" is the path on the system in which installChrokitDependencies.R and shinyapp.R scripts are located. For example, if the    ChroKit source code has been downloaded in "/Users/ocroci/Downloads" directory, just type:\
  ``` > setwd("/Users/ocroci/Downloads/ChroKit-master/source/")```\
  in the R console.
- Run the script and wait for all dependencies to be downloaded from internet, by typing the following command in the R console:\
 ``` > source ("installChrokitDependencies.R")```\
 
 
 Alternatively, make sure the following R libraries are installed:
  - shiny
  - shinyFiles
  - shinydashboard
  - shinyWidgets
  - fastcluster
  - VennDiagram
  - GenomicRanges
  - rtracklayer
  - data.table
  - RColorBrewer
  - Rsamtools
  - ppcor
  - inline

Note: if the error *"Bioconductor version X.X requires R version Y.Y"* is displayed in the R console, it means that the bioconductor version is not compatible with the version of the installed R interpreter. To solve this, edit the *installChrokitDependencies.R* text file and change the variable **bioCversion** to insert the version number compatible with your R interpreter. 
For example, R 4.0 is compatible with bioconductor version 3.11. In that case, change the default value:\
```bioCversion="3.8"```\
in\
```bioCversion="3.11"```\
and source the *installChrokitDependencies.R* script again:\
``` > source ("installChrokitDependencies.R")```\

### Basic setup
Some parameters could be set in the **shinyapp.r** script, such as the listening port or the number of cores, as well as the colors available for the heatmaps.
- The variable **Port** specify the listening port of the program. Default: 6060.
- The variable **nc** specify the number of cores that will be used for computation. The higher the number, the faster the program will be, but it will require more RAM. Windows users will always use 1 single core for operations because of technical issues.
- The variable **ColsArray** specify all colors available in the palettes for heatmaps (gradient from white)
- The variable **bioCversion** specify the appropriate version of Bioconductor for your R interpreter for the download of databases

### launch the program
- Launch the application using\
  ``` > source("shinyapp.r")```
- If the message "Listening on http://0.0.0.0:6060" is displayed in the R console, the application is running properly. Open the application using your web browser.
  - if you are using a personal computer, go to:\
    ```127.0.0.1:6060```
  - if you are using a remote machine, go to:\
    ```<IP>:6060 ```\
    where \<IP\> is the IP address of the remote machine, and 6060 is the listening port used by ChroKit.


**IMPORTANT**: to carry out gene ontology analyses, you must put gene signatures under the appContent/signatures directory. Those signatures must be in gmt format, and their file name must end with **\_symbols.gmt**. Signatures from MSigDB (Molecular Signature Database) (https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb) are already preloaded. 


## Tutorial
To learn the basics on how to analyse NGS data with ChroKit, please follow the tutorials at this page: https://ocroci.github.io/ChroKit/ using sample data.


## Credits
- The function for drawing color bars was adapted from John Colby (stackoverflow) http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
- The function to generate a number of most distinctive colors in R was taken from http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r (Megatron)
- "distanceFromTSS3", "countOverlapsInBins", "GRbaseCoverage2" and "summitFromBaseCoverage" functions were taken (or adapted) from the compEpiTools R package:\
**Kishore K, de Pretis S, Lister R, Morelli MJ, Bianchi V, Amati B, Ecker JR, Pelizzola M (2015). “methylPipe and compEpiTools: a suite of R packages for the integrative analysis of epigenomics data.” BMC Bioinformatics. doi: 10.1186/s12859-015-0742-6.**
- Pre-loaded genesets for gene ontology analyses were downloaded from MSigDB v.7.4 (copyright (c) 2004-2020 Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California) (http://www.gsea-msigdb.org/gsea/downloads.jsp) and are under Creative Commons Attribution 4.0 International License (https://creativecommons.org/licenses/by/4.0/legalcode)


## Citation
If you use this framework for your project, please acknowledge Ottavio Croci, PhD at Center for Genomic Science of IIT@CGS, Fondazione Istituto Italiano di Tecnologia. 


## License
Copyright (c) 2020 Ottavio Croci\
[The GNU GENERAL PUBLIC LICENSE v3.0](LICENSE)

