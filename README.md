# ChroKit
ChroKit (The **Chro**matin Tool**kit**) is a Shiny-based framework to analyze and visualize interactively genomic data.

<img src="https://github.com/ocroci/ChroKit/blob/master/logo2.png" height="50%" width="50%">

## Installation
- Download and install the R interpreter on your computer or on a remote machine (suggested version 3.5 or higher)
- Open the R interpreter and go into the main source directory of the program
- Run the script and wait for all dependencies to be downloaded from internet:\
 ``` > source ("installChrokitDependencies.R")```
- Launch the application using\
  ``` > source("shinyapp.r")```
- Go with your browser and open the application.
  - if you are using your computer, go to:\
    ```127.0.0.1:6060```
  - if you are using a remote machine, go to:\
    ```<IP>:6060 ```\
    where \<IP\> is the IP address of the remote machine
    
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
