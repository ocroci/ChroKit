options(install.packages.compile.from.source = "always")
############################################
#General parameters
############################################
#select the listening port
Port=6061
### select the number of cores for execution
### on windows systems, it will put nc=1
nc=4
#define colors for the palettes (from white)
ColsArray=c("red","#e202b2","hotpink","deepskyblue","darkorange","darkorchid",
			"#00CC59","#00ba96","#f46036","#031d44","#138a36","#3891a6","#748386","#fff05a",
			"#824c71","#7a306c","#04395e","#ffaa5a",
			"#69306d","#8B0000","#595959","#4169e1","#f03a47","blue","#c8ad4d","#85cb33")
#define the bioCversion version for installing packages with BiocManager
#R 3.5 wants bioCversion="3.8", R 3.6 wants bioCversion="3.9" and so on...
numeric_version=as.character(tools:::.BioC_version_associated_with_R_version())
bioCversion=numeric_version




############################################
# Parameters (port, user, bioCversion and number of cores) can be 
# introduced also as command line arguments
############################################
### command line arguments (port, USER)
### if multiple users are using the program inside the same machine, USER
### variable can be used to keep trace of the USER
### if port is not set, default is 6060
args = commandArgs(trailingOnly=TRUE)
print (args)
if (!is.na(args[1])){
	port=as.numeric(args[1])
}else{
	port=Port
}

if(!is.na(args[2])){
	USER=args[2]
}else{
	USER=NULL
}

if(!is.na(as.character(args[3]))){
	bioCversion=args[3]
}

if(!is.na(args[4])){
	nc=args[4]
}









#set the root directory for the shinyFiles buttons according to the OS type in use
#if windows, nc will be forced to 1

if(.Platform$OS.type=="unix"){
	separation=rootdir=rootsavedir="/"
}else{
	rootdir="C:\\Users\\"
	rootsavedir="C:/Users"
	separation="\\"
	nc=1
}

############################################
# Loading required libraries
############################################

#shiny, shinyFiles, shinydashboard and shinyWidgets for the GUI
library(shiny)
library(shinyFiles)
library(shinydashboard)
library(shinyWidgets)
#fastcluster for efficient clusterings
library(fastcluster)
#VennDiagram for the pairwise overlaps venn
library(VennDiagram)
#rtracklayer for the Wig file support
library(rtracklayer)
#GenomicRanges for all ROI operations
library(GenomicRanges)
#data.table for efficient tables management
library(data.table)
#RColorBrewer for colors and color palettes
library(RColorBrewer)
#Rsamtools for BAM files support
library(Rsamtools)
#ppcor for calculating partial correlations
library(ppcor)
#inline, for compiling Rcpp functions
library(inline)
#bamsignals, for efficient pileups/coverage computation from BAM files
library(bamsignals)




#DEFINE all possible loadable or downloadable databases: GLOBAL variables
#define downloadable and available databases from current session
x=rownames(installed.packages())
##define TxDb transcripts packages of bioconductor for each genome assembly (KnownGene, refGene or ensGene)
##and name each element with the corresponding annotation libraries of the same organism
all_avail_assemblies=c(
"TxDb.Btaurus.UCSC.bosTau8.refGene",  
"TxDb.Celegans.UCSC.ce11.refGene",  
"TxDb.Celegans.UCSC.ce6.ensGene",  #
"TxDb.Cfamiliaris.UCSC.canFam3.refGene",
"TxDb.Dmelanogaster.UCSC.dm3.ensGene",  #
"TxDb.Dmelanogaster.UCSC.dm6.ensGene",  #
"TxDb.Drerio.UCSC.danRer10.refGene",
"TxDb.Ggallus.UCSC.galGal4.refGene",  
"TxDb.Ggallus.UCSC.galGal5.refGene",      
"TxDb.Hsapiens.UCSC.hg18.knownGene",  
"TxDb.Hsapiens.UCSC.hg19.knownGene",  
"TxDb.Hsapiens.UCSC.hg38.knownGene",
"TxDb.Mmulatta.UCSC.rheMac3.refGene",  
"TxDb.Mmulatta.UCSC.rheMac8.refGene",     
"TxDb.Mmusculus.UCSC.mm10.knownGene",  
"TxDb.Mmusculus.UCSC.mm9.knownGene"  ,
"TxDb.Ptroglodytes.UCSC.panTro4.refGene",  
"TxDb.Ptroglodytes.UCSC.panTro5.refGene",    
"TxDb.Rnorvegicus.UCSC.rn4.ensGene",  #
"TxDb.Rnorvegicus.UCSC.rn5.refGene",  
"TxDb.Rnorvegicus.UCSC.rn6.refGene"
)    

names(all_avail_assemblies)=c(
"org.Bt.eg.db",   
"org.Ce.eg.db",
"org.Ce.eg.db",
"org.Cf.eg.db",
"org.Dm.eg.db",
"org.Dm.eg.db",
"org.Dr.eg.db",
"org.Gg.eg.db", 
"org.Gg.eg.db",
"org.Hs.eg.db",
"org.Hs.eg.db",
"org.Hs.eg.db",
"org.Mmu.eg.db",
"org.Mmu.eg.db",
"org.Mm.eg.db",
"org.Mm.eg.db",
"org.Pt.eg.db",
"org.Pt.eg.db",
"org.Rn.eg.db",
"org.Rn.eg.db",
"org.Rn.eg.db"
)

#here define packages already installed, and create avail_assemblies vector.
#both org and txdb must be present to be added to avail_assemblies (variable x for installed packages).
#then, immediately put the variable in a reactive variable that can change

#after that, we know which can be served in the interface as "loadable"

#define functions
source("appContent/_functions.R")

#open help buttons (?) text messages:
source("appContent/_help_messages.R")

#define the assemblies available at th startup
availASSEMBLIES=getExistingDB(all_avail_assemblies)$assemblies_we_have
missingASSEMBLIES=getExistingDB(all_avail_assemblies)$assemblies_we_donthave


#define genesets (MSigDB) in GMT format as global variable to show to the user for GO analyses:
totMsigDB=dir(paste(getwd(),"/appContent/signatures/",sep=""),full.names=TRUE)
GenesetsGMT=grep("*_symbols.gmt$",totMsigDB,value=TRUE)
#attrib names to GenesetsGMT (clean base names)
bN=basename(GenesetsGMT)
GenesetsGMT=as.list(GenesetsGMT)
sN=strsplit(bN,split="_symbol")
bN=sapply(sN,"[[",1)
names(GenesetsGMT)=bN
temporary_GMTstorage=as.list(rep(NA,length(GenesetsGMT)))
names(temporary_GMTstorage)=bN




#run the app
runApp("appContent",host="0.0.0.0",port=port,launch.browser=FALSE)


