#biocversion
#define the bioCversion version for installing packages with BiocManager
#R 3.5 wants bioCversion="3.8", R 3.6 wants bioCversion="3.9" and so on...
numeric_version=as.character(tools:::.BioC_version_associated_with_R_version())
bioCversion=numeric_version


#source _functions
source("appContent/_functions.R")
library(GenomicRanges)

#define all available assemblies. Update if necessary
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

miss=getExistingDB(all_avail_assemblies)
miss=miss$assemblies_we_donthave

#download missing dbs
if (length(miss)>0){
	for (i in 1:length(miss)){downloadDB(miss[i],all_avail_assemblies)}
}


#pre-calculate tables of each DB
dbs=getExistingDB(all_avail_assemblies)
dbs=dbs$assemblies_we_have

if (!file.exists("appContent/assemblies")){
	dir.create("appContent/assemblies")
}

for (i in 1:length(dbs)){
	assembly=dbs[i]
	tab__=extractFromDB(assembly,avail_assemblies=all_avail_assemblies)
	df=as.data.frame(tab__$transcripts)
	write.table(df,file=paste0("appContent/assemblies/",assembly),col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE )
}




