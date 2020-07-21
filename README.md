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
- Launch the application using\
  ``` > source("shinyapp.r")```
- Go with your browser and open the application.
  - if you are using your computer, go to:\
    ```127.0.0.1:6060```
  - if you are using a remote machine, go to:\
    ```<IP>:6060 ```\
    where \<IP\> is the IP address of the remote machine

## Basic setup
Some parameters could be set in the **shinyapp.r** script, such as the listening port or the number of cores, as well as the colors available for the heatmaps.

## Windows users
Make sure to install the R interpreter in a directory path without spaces: when prompted the path for installation choose C:\R\ as the path.
After installing the R interpreter, install Rtools; then, find the path in which the gcc.exe executable had been installed (usually C:\Rtools\mingw_32\bin) and the directory in which all the Rtools binaries had been installed (usually under C:\Rtools\bin). Finally, add these two directories to you Windows PATH variable.
To temporarily set these paths only for the current R session, from the R interpreter, type the following lines of code:
```
rtools<-"C:\\Rtools\\bin"
gcc<-"C:\\Rtools\\mingw_32\\bin"
path <- strsplit(Sys.getenv("PATH"),";")[[1]]
new_path<-c(rtools,gcc,path)
new_path<-new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH=paste(new_path,collapse=";"))
```
Where rtools is the path in which all Rtools binaries had been installed and gcc is the directory in which gcc.exe executable had been installed. 
Note for windows users: WIG file association is not supported and only BAM file association is allowed; moreover, only one core is allowed, due to the use of “parallel” library, which works only on UNIX operating systems.

## Credits
- Function for drawing color bars was adapted from John Colby (stackoverflow) http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
- Function to generate a number of most distinctive colors in R was taken from http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r (Megatron)
- "distanceFromTSS3", "countOverlapsInBins", "GRbaseCoverage2" and "summitFromBaseCoverage" functions were taken (or adapted) from compEpiTools R package:\
**Kishore K, de Pretis S, Lister R, Morelli MJ, Bianchi V, Amati B, Ecker JR, Pelizzola M (2015). “methylPipe and compEpiTools: a suite of R packages for the integrative analysis of epigenomics data.” BMC Bioinformatics. doi: 10.1186/s12859-015-0742-6.**
- The instructions to use the program under Windows operating system were taken from https://stackoverflow.com/questions/23141982/inline-function-code-doesnt-compile (Richie Cotton)

## Citation
If you use this framework for your project, please acknowledge Ottavio Croci, PhD at Center for Genomic Science of IIT@SEMM

## License
Copyright (c) 2020 Ottavio Croci\
[The GNU GENERAL PUBLIC LICENSE v3.0](LICENSE)
