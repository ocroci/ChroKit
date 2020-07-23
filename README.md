# ChroKit
ChroKit (The **Chro**matin Tool**kit**) is a Shiny-based framework to analyze and visualize interactively genomic data.

<img src="https://github.com/ocroci/ChroKit/blob/master/logo2.png" height="50%" width="50%">

## Install dependencies and launch the program
- Download and install the R interpreter on your computer or on a remote machine (suggested version 3.5 or higher)
- Open the R interpreter and go into the main source directory of the program
- Run the script and wait for all dependencies to be downloaded from internet:\
 ``` > source ("installChrokitDependencies.R")```\
 Alternatively, make sure to the following R libraries installed:
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
  
If the error *"Bioconductor version X.X requires R version Y.Y"* appears, it means that the bioconductor version is not compatible with the version of R interpreter installed in the computer. To solve this, edit the *installChrokitDependencies.R* file and change the variable **bioCversion** to a different number, representing a bioconductor version compatible with your current R interpreter. For example, R 4.0 is compatible with bioconductor version 3.11: in that case, change:\
```bioCversion="3.8"```\
in\
```bioCversion="3.11"```\
and source the *installChrokitDependencies.R* script again
- Launch the application using\
  ``` > source("shinyapp.r")```
- Go with your browser and open the application.
  - if you are using your computer, go to:\
    ```127.0.0.1:6060```
  - if you are using a remote machine, go to:\
    ```<IP>:6060 ```\
    where \<IP\> is the IP address of the remote machine

**IMPORTANT**: to carry out gene ontology analyses, you must put gene signatures under appContent/signatures directory. Those signatures must be in gmt format, and their file name must end with **\_symbols.gmt**. Suggested: download MSigDB (Molecular Signature Database) signatures from https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb. 

### Windows users
Make sure to install the R interpreter in a directory path without spaces: when prompted the path for installation choose C:\R\ as the path.
After installing the R interpreter, install Rtools; then, modify the PATH variable to include also all the binaries of Rtools.
To install Rtools and modify the PATH variable, follow the instructions at the link: https://cran.r-project.org/bin/windows/Rtools/

Note for windows users: WIG file association is not supported and only BAM file association is allowed; moreover, only one core is allowed, due to the use of “parallel” library, which works only on UNIX operating systems.


## Basic setup
Some parameters could be set in the **shinyapp.r** script, such as the listening port or the number of cores, as well as the colors available for the heatmaps.
- The variable **Port** specify the listening port of the program
- The variable **nc** specify the number of cores that will be used for computation. The higher the number, the faster the program will be, but it will require more RAM. Windows users will always use 1 single core for operations.
- The variable **ColsArray** specify all colors available for palettes in heatmaps (gradient from white)


## Credits
- Function for drawing color bars was adapted from John Colby (stackoverflow) http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
- Function to generate a number of most distinctive colors in R was taken from http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r (Megatron)
- "distanceFromTSS3", "countOverlapsInBins", "GRbaseCoverage2" and "summitFromBaseCoverage" functions were taken (or adapted) from compEpiTools R package:\
**Kishore K, de Pretis S, Lister R, Morelli MJ, Bianchi V, Amati B, Ecker JR, Pelizzola M (2015). “methylPipe and compEpiTools: a suite of R packages for the integrative analysis of epigenomics data.” BMC Bioinformatics. doi: 10.1186/s12859-015-0742-6.**


## Citation
If you use this framework for your project, please acknowledge Ottavio Croci, PhD at Center for Genomic Science of IIT@SEMM


## License
Copyright (c) 2020 Ottavio Croci\
[The GNU GENERAL PUBLIC LICENSE v3.0](LICENSE)
