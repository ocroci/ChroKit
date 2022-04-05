#####FUNCTIONS#####

##functions for get and set the enrichment blocks in the reactiveVariables ##
# those values are not part of the ROI object anymore, but contained in
# Enrichlist<-reactiveValues(rawcoverage=NULL,normfactlist=NULL)
#these functions cannot read or modify reactiveValues (outside the scope) but give useful info
#on how to do that

getEnrichList<-function(ROIname) {
	nome=ROIname
	pos=match(nome,names(Enrichlist$rawcoverage))
	rawvals=Enrichlist$rawcoverage[[pos]]
	normvals=Enrichlist$normfactlist[[pos]]
	return(list(rawvals,normvals))

}
updateEnrichList<-function(ROIname,enrichlist,normfact) {
	nome=ROIname
	pos=match(nome,names(Enrichlist$rawcoverage))
	#are we sure absence of mem leaks? Return value?
	Enrichlist$rawcoverage[[pos]]=enrichlist
	Enrichlist$normfactlist[[pos]]=normfact
}

#update enrich list before removing an entire ROI
removeelementEnrichList<-function(ROIname) {
	nome=ROIname
	pos=match(nome,names(Enrichlist$rawcoverage))
	Enrichlist$rawcoverage[pos]<-NULL
	Enrichlist$normfactlist[pos]<-NULL
}
renameROIEnrichlist<-function(ROIname,newname) {
	nome=ROIname
	pos=match(nome,names(Enrichlist$rawcoverage))
	names(Enrichlist$rawcoverage)[pos]<-newname
	names(Enrichlist$normfactlist)[pos]<-newname
}


#function to display the help buttons, here choose color, symbol, ID and title of the box
#thanks to Roman Hillje, useful code taken from cerebroApp (https://github.com/romanhaa/cerebroApp/blob/master/inst/shiny/overview/UI.R)
boxHelp<-function(ID,title,col="#FF8080",symbol="?") {
  #col: background color for the help button
  #symbol: symbol to show in the button (default "?")
  #ID: ID to use for triggering the server part
  #title: title of the box
  tagList(
    p(title,style = "padding-right: 5px; display: inline"),
    actionButton(
      inputId = ID,
      label = symbol,
      icon = NULL,
      class = "btn-xs",
      title = "Show help",
      style = paste("margin-right: 5px; color: black; background-color: ",col,sep="")
    )
  )

}


#function to trigger the help button (server side).
#parameter is the object in help_messages.R script which is a list of title,text
#thanks to Roman Hillje, useful code taken from cerebroApp (https://github.com/romanhaa/cerebroApp/blob/master/inst/shiny/overview/UI.R)
boxHelpServer<-function(object) {
  #object: defined in help_messages.R, list containing a title and text (HTML format)
  showModal(
    modalDialog(
      object$text,
      title = object$title,
      easyClose = TRUE,
      footer = NULL
    )
  )
}





#checks whether an input is valid or set (different from "", NA, NULL...)
#valid only for single values
isvalid<-function(inputfield) {
  innocentUntilProvenGuilty=TRUE
  #if length of input is 0, refuse
  if(length(inputfield)==0){
    innocentUntilProvenGuilty=FALSE
  }
  #if length of input is 1, checks...

  if(length(inputfield)==1){
    if(!is.null(inputfield)){
      #if input is ""
      if(nchar(inputfield)==0 | is.na(nchar(inputfield))) {
        innocentUntilProvenGuilty=FALSE
      }
      #if input is NA
      if(is.na(inputfield)){
        innocentUntilProvenGuilty=FALSE
      }    
    }else{
      innocentUntilProvenGuilty=FALSE
    }
  }

  return(innocentUntilProvenGuilty)
}


##reads BED or GTF/GFF and automatically recognise the file format
readBEDGTFF<-function(bedpath,Header=TRUE,Skip=0){
  #bedpath= path to the file to be opened
  #Header= whether to include the header
  #Skip= number defining how many lines of the file to skip 
  if (!file.exists(bedpath)){
      stop("file does not exist")
    }
    rt=suppressWarnings(read.table(bedpath,header=Header,sep="\t",skip=Skip))
    ncl=ncol(rt)
    if (ncl<3){
      stop("too few columns...")
    }
    #check col 2 and 3.
    numeric2=any(is.na(suppressWarnings(as.numeric(as.character(rt[,2])))))
    numeric3=any(is.na(suppressWarnings(as.numeric(as.character(rt[,3])))))   
    # #here check if column1 (seqnames / chr names) are in the correct UCSC format.
    # #must start with "chr...", otherwise format is not valid
    # seqnamecheck=as.character(rt[,1])
    # startseqnames=substr(seqnamecheck,start=1,stop=3)
    # if (all(startseqnames=="chr")){
    #   print("All seqnames starts correctly with chr")
    # }else if (!all(startseqnames=="chr")&"chr"%in% startseqnames){
    #   print("Warning: some of the seqnames do not start with chr")
    # }else if (!"chr"%in% startseqnames){
    #   print("None of the seqnames start with chr. Check the format.")
    #   stop("None of the seqnames start with chr. Check the format.")
    # }

    if (ncl <8){
      if (!numeric2&!numeric3){
        if (ncl>=4){
          #check if 4th column is a strand
          sign=as.character(rt[,4])
          logic=  sign=="+" | sign=="-" | sign=="*" # | sign=="." ??
          logicneg=!logic
          if(any(logicneg)){
            return(rt[,1:3])
          }else{
            return(rt[,1:4])
          }

        }else{
          return(rt[,1:3])
        }         
      }else{
        stop("seems a BED, but column 2 or 3 are not numeric...")
      }
      #if numeric2 is FALSE and numeric3 is false (all numeric) or 
      #check if there is a 4th column:

    }else{
      #can be BED or GTF. check strand at 7th column.
      numeric4=any(is.na(suppressWarnings(as.numeric(as.character(rt[,4])))))
      numeric5=any(is.na(suppressWarnings(as.numeric(as.character(rt[,5])))))  
      sign=as.character(rt[,7])
      logic=  sign=="+" | sign=="-" | sign=="*" | sign =="."
      logicneg=!logic
      if (!numeric4&!numeric5 & !any(logicneg)){
        #if both numeric, and 7th is strand, => GTF
        df=data.frame(rt[,1],rt[,4],rt[,5],rt[,7])
        colnames(df)=colnames(rt)[c(1,4,5,7)]
        return(df)
      }else{
        #in this case can be a BED with multiple columns. check if cols 2 and 3 are numeric
        if(!numeric2&!numeric3){
        sign=as.character(rt[,4])
        logic=  sign=="+" | sign=="-" | sign=="*" # | sign=="." ??
        logicneg=!logic
        if(any(logicneg)){
          return(rt[,1:3])
        }else{
          return(rt[,1:4])
        }
        }else{
          stop("it seems not a BED nor a GTF... some of columns 2/3/4/5 not numeric...")
        }
      }
        
    }

}




####################################################################
####################################################################
#DATABASES
####################################################################
####################################################################
checkBiocConnection<-function(url="bioconductor.org/packages/") {
  library(RCurl)
  connectivity=FALSE
  tryCatch({ 
    x=getURL(url)
    connectivity=TRUE
  },
  warning = function( w ){
    connectivity=FALSE
  },
  error = function( err ){
    connectivity=FALSE
  })
  return(connectivity)
}
downloadDB<-function(assembly,avail_assemblies) {
  #assembly: character with the name of the assembly for which txdb and org DBs 
  #   have to be installed
  #avail assemblies= named charactr vector with all existing Txdb and names are
  #   the relative org db
  av_packages=rownames(installed.packages())
  if(class(assembly)!="character"){
    stop("'assembly' must be a named character")
  }  
  if(class(avail_assemblies)!="character"){
    stop("'avail_assemblies' must be a named character")
  }  
  #search, from the assembly, the couple txdb/org inside avail_assemblies vector:
  pos=grep(assembly,avail_assemblies)
  if(length(pos)!=1){
    stop("'assembly' not found in 'avail_assemblies' or more than 1 match...")
  }

  todownload_txdb=avail_assemblies[pos]
  todownload_org=names(todownload_txdb)
  todownload_txdb=unname(todownload_txdb)
  #check if we have both, one of them or none
  is_txdb=todownload_txdb%in%av_packages
  is_org=todownload_org%in%av_packages

  #if both, exit form the function: we cannot do anything useful...
  #but check, because this function should be called only when at least one DB is not available
  if(is_txdb&is_org){
    print(paste("both 'org' and 'txdb' for the",assembly,"assembly were available!"))
    return()
  }
  #for R versions <=3.4. For R >= 3.5 you can download the package BiocManager
  #and use BiocManager::install("package")  

  #here, check if R>3.5 and BiocManager is installed
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
  if(! ("BiocManager" %in% av_packages) & R35){
    print("Installing BiocManager package for R > 3.5.0 ...")
    install.packages("BiocManager")
  }else{
    #print("BiocManager package already installed...")
  }

  if(R35){
    #upgrading bioconductor packages to the correct version bioCversion detected with .BioC_version_associated_with_R_version() function
    print(paste("updating bioconductor packages to the version",bioCversion))
    BiocManager::install(version = bioCversion,ask=FALSE,force=TRUE)
    if(is_txdb&!is_org){
      print(paste("Downloading",todownload_org,"..."))
      BiocManager::install(todownload_org, version = bioCversion,ask=FALSE,force=TRUE)
    }
    if(!is_txdb&is_org){
      print(paste("Downloading",todownload_txdb,"..."))
      BiocManager::install(todownload_txdb, version = bioCversion,ask=FALSE,force=TRUE)
    }
    if(!is_txdb&!is_org){
      print(paste("Downloading",todownload_org,"..."))
      BiocManager::install(todownload_org, version = bioCversion,ask=FALSE,force=TRUE)
      print(paste("Downloading",todownload_txdb,"..."))
      BiocManager::install(todownload_txdb, version = bioCversion,ask=FALSE,force=TRUE)
    }  
  }else{
    source("https://bioconductor.org/biocLite.R")
    if(is_txdb&!is_org){
      print(paste("Downloading",todownload_org,"..."))
      biocLite(todownload_org,suppressUpdates=TRUE,ask=FALSE)
    }
    if(!is_txdb&is_org){
      print(paste("Downloading",todownload_txdb,"..."))
      biocLite(todownload_txdb,suppressUpdates=TRUE,ask=FALSE)
    }
    if(!is_txdb&!is_org){
      print(paste("Downloading",todownload_txdb,"..."))
      biocLite(todownload_txdb,suppressUpdates=TRUE,ask=FALSE)
      print(paste("Downloading",todownload_org,"..."))
      biocLite(todownload_org,suppressUpdates=TRUE,ask=FALSE)
    }      
  }



}


#from all possible assemblies, return those for which we have libraries
getExistingDB<-function(avail_assemblies) {
  #avail_assemblies is a named array; elements are character strings of TxDb bioC annotation packages,
  #  while the names are the org.XX.eg.db of the species associated to that genome assembly TxDb.
  #  (defined in shinyapp.r initial script)
  if(class(avail_assemblies)!="character"){
    stop("'avail_assemblies' must be a named character")
  }
  av_packages=rownames(installed.packages())
  query=paste("TxDb\\..+\\.UCSC\\..+\\..+Gene$",sep="")
  #find the position of installed packages in which we have txdb and org databases
  pos=grep(query,av_packages)
  current_txdbs=av_packages[pos]
  #now detect if we have the txdb
  pos=match(current_txdbs,avail_assemblies) 
  #retrieve txdbs in avail_assemblies found in avail packages
  current_txdbs=avail_assemblies[pos]
  current_txdbs=current_txdbs[!is.na(current_txdbs)]
  #now check which of them have also the relative org DB
  orgs=names(current_txdbs)
  pos=orgs%in%av_packages
  final_availables=current_txdbs[pos]

  #now extract the available assembly codes:
  assemblies_we_have=sapply(strsplit(final_availables,split="\\."),"[[",4)
  notpresent=avail_assemblies[!avail_assemblies%in%final_availables]
  assemblies_we_donthave=sapply(strsplit(notpresent,split="\\."),"[[",4)
  dblist=list(unname(assemblies_we_have),unname(assemblies_we_donthave))
  names(dblist)=c("assemblies_we_have","assemblies_we_donthave")
  return(dblist)
}


#takes the name of the assembly (eg. "hg19") and the array
#of all possible libraries of all assemblies; return a list
#of transcripts with all possible IDs/symbols and their fraction
#out of the total amount of transcripts
extractFromDB<-function(assembly,avail_assemblies) {
  #assembly is a character string representing the genome assembly of the species 
  # for example, mm10, hg19, rheMac3, rn6, which is part of the name in TxDb libraries
  #avail_assemblies is a named array; elements are character strings of TxDb bioC annotation packages,
  # while the names are the org.XX.eg.db of the species associated to that genome assembly TxDb
  if(class(assembly)!="character"){
    stop("'assembly' must be a character of the genome assembly (example: mm9,hg19)")
  }
  if(length(assembly)!=1){
    stop("'assembly' length must be 1")
  }
  if(class(avail_assemblies)!="character"){
    stop("'avail_assemblies' must be a named character")
  }
  if(!any(grepl(assembly,avail_assemblies))) {
    stop("'assembly' must be in the avail_assemblies")
  }

  #extract the correct libraries (txdb and org) and import them with library()
  query=paste("TxDb\\..+\\.UCSC\\.",assembly,"\\..+Gene$",sep="")
  pos=grep(query,avail_assemblies)
  if(length(pos)==0){
    stop("'assembly must be in 'avail_assemblies'")
  }
  txdb_lib=avail_assemblies[pos]
  org_lib=names(txdb_lib)
  txdb_lib=unname(txdb_lib)
  
  #verify if both txdb and org libs are available as packages, otherwise stop function:
  #here, it will tell you which packages are missing
  av_packages=rownames(installed.packages())
  txdb_confirm=TRUE
  org_confirm=TRUE
  if(! (org_lib %in% av_packages)){
    org_confirm=FALSE
  }
  if(! (txdb_lib %in% av_packages)){
    txdb_confirm=FALSE
  }
  if(!org_confirm &!txdb_confirm ){
    stop(paste(org_lib,"and",txdb_lib,"packages not found, maybe not installed?"))
  }
  if(!org_confirm&txdb_confirm){
    stop(paste(org_lib,"package not found, maybe not installed?"))
  }
  if(org_confirm&!txdb_confirm){
    stop(paste(txdb_lib,"package not found, maybe not installed?"))
  }

  ###if we are here, both txdb and org packages are found###
  #import txdb package
  library(txdb_lib,character.only=TRUE)
  #import org.XXX package
  library(org_lib,character.only=TRUE)
  transc=transcripts(eval(parse(text = txdb_lib)),columns=c("gene_id","tx_id"))
  transc=unique(transc)
  #check if txdb ends with "knowngene","refGene" or "ensGene" and do the 
  #various translations appropriately
  tocheck=strsplit(txdb_lib,split="\\.")[[1]][5]  #5 is the last part of name txdb.XX.UCSC.XX.??Gene
  baseNameOrg=strsplit(org_lib,split="\\.")[[1]]
  baseNameOrg=paste(baseNameOrg[1:(length(baseNameOrg)-1)],collapse=".")
  if(tocheck=="knownGene" | tocheck=="refGene"){
    #we have already entrez IDs. Load SYMBOL, then ENSEMBL, then REFSEQ (and GENENAME).
    #then, append everything as metadata in the transcripts
    dictSymbol=paste(baseNameOrg,"SYMBOL",sep="")
    dictEnsembl=paste(baseNameOrg,"ENSEMBL",sep="")
    dictRefseq=paste(baseNameOrg,"REFSEQ",sep="")
    #dictGenename=paste(baseNameOrg,"GENENAME",sep="")

    dictSymbol=as.list(eval(parse(text=dictSymbol))[mappedkeys(eval(parse(text=dictSymbol)))])
    dictEnsembl=as.list(eval(parse(text=dictEnsembl))[mappedkeys(eval(parse(text=dictEnsembl)))])
    dictRefseq=as.list(eval(parse(text=dictRefseq))[mappedkeys(eval(parse(text=dictRefseq)))])
    #dictGenename=as.list(eval(parse(text=dictGenename))[mappedkeys(eval(parse(text=dictGenename)))])
    #use all these dictionaries to populate transcripts table granges

    totalSymbols=dictSymbol[as.character(transc$gene_id)]
    totalEnsembl=dictEnsembl[as.character(transc$gene_id)]
    totalRefseq=dictRefseq[as.character(transc$gene_id)]
    #totalGenename=dictSymbol[as.character(transc$gene_id)]

    #if more than 2 match per single ENTREZ found, pick the first one:
    totalSymbols=lapply(totalSymbols,"[[",1)
    totalEnsembl=lapply(totalEnsembl,"[[",1)
    totalRefseq=lapply(totalRefseq,"[[",1)
    #totalGenename=lapply(totalGenename,"[[",1)

    #transform NULL in NA
    totalSymbols[sapply(totalSymbols, is.null)]=NA
    totalEnsembl[sapply(totalEnsembl, is.null)]=NA
    totalRefseq[sapply(totalRefseq, is.null)]=NA
    #totalGenename[sapply(totalGenename, is.null)]=NA

    totalSymbols=unlist(totalSymbols)
    totalEnsembl=unlist(totalEnsembl)
    totalRefseq=unlist(totalRefseq)
    #totalGenename=unlist(totalGenename)

    #finally append to transcript table
    df=data.frame(gene_id=as.character(transc$gene_id),
            symbol=totalSymbols,
            ensembl_id=totalEnsembl,
            refSeq_id=totalRefseq)

    elementMetadata(transc)=df

  }else if (tocheck=="ensGene"){
    #further subset if C. elegans (ce6)(ENSEMBL transcript name) or the others (real ENSEMBL IDs)
    if(assembly=="ce6"){
      #we have ensembl transcripts names for ce6. We don't know why!
      #we need ENSEMBLTRANS to retrieve entrez ID. then all the others (SYMBOL,REFSEQ,ENSEMBL)
      dictEnsembltrans=paste(baseNameOrg,"ENSEMBLTRANS2EG",sep="")
      
      dictSymbol=paste(baseNameOrg,"SYMBOL",sep="")
      dictEnsembl=paste(baseNameOrg,"ENSEMBL",sep="")
      dictRefseq=paste(baseNameOrg,"REFSEQ",sep="")

      dictEnsembltrans=as.list(eval(parse(text=dictEnsembltrans))[mappedkeys(eval(parse(text=dictEnsembltrans)))])
      dictSymbol=as.list(eval(parse(text=dictSymbol))[mappedkeys(eval(parse(text=dictSymbol)))])
      dictEnsembl=as.list(eval(parse(text=dictEnsembl))[mappedkeys(eval(parse(text=dictEnsembl)))])
      dictRefseq=as.list(eval(parse(text=dictRefseq))[mappedkeys(eval(parse(text=dictRefseq)))])
      
      #ensembl transcript name in transc are real names plus .*: we have to 
      #remove the dot and the following part
      ensembltranscname=strsplit(as.character(transc$gene_id),split="\\.")
      ensembltranscname=lapply(ensembltranscname,"[[",1)
      ensembltranscname[sapply(ensembltranscname, is.null)]=NA
      ensembltranscname=unlist(ensembltranscname)

      #first: convert the ensembl trancript names in ENTREZ:
      totalEntrez=dictEnsembltrans[ensembltranscname]
      totalEntrez=lapply(totalEntrez,"[[",1)
      totalEntrez[sapply(totalEntrez, is.null)]=NA
      totalEntrez=unlist(totalEntrez)

      #now, using the ENTREZ, find symbols, ensembl and refseq:
      totalSymbols=dictSymbol[as.character(totalEntrez)]
      totalEnsembl=dictEnsembl[as.character(totalEntrez)]
      totalRefseq=dictRefseq[as.character(totalEntrez)]
      #if more than 2 match per single ENTREZ found, pick the first one:
      totalSymbols=lapply(totalSymbols,"[[",1)
      totalEnsembl=lapply(totalEnsembl,"[[",1)
      totalRefseq=lapply(totalRefseq,"[[",1)
      #transform NULL in NA
      totalSymbols[sapply(totalSymbols, is.null)]=NA
      totalEnsembl[sapply(totalEnsembl, is.null)]=NA
      totalRefseq[sapply(totalRefseq, is.null)]=NA
      #finally,unlist
      totalSymbols=unlist(totalSymbols)
      totalEnsembl=unlist(totalEnsembl)
      totalRefseq=unlist(totalRefseq)
      #finally append to transcript table
      df=data.frame(gene_id=totalEntrez,
            symbol=totalSymbols,
            ensembl_id=totalEnsembl,
            refSeq_id=totalRefseq)
      elementMetadata(transc)=df


    }else{
      #We have ensembls. Import ENSEMBL2EG, then use SYMBOL and REFSEQ on ENTREZ
      dictSymbol=paste(baseNameOrg,"SYMBOL",sep="")
      dictEnsembl2eg=paste(baseNameOrg,"ENSEMBL2EG",sep="")
      dictRefseq=paste(baseNameOrg,"REFSEQ",sep="")

      dictSymbol=as.list(eval(parse(text=dictSymbol))[mappedkeys(eval(parse(text=dictSymbol)))])
      dictEnsembl2eg=as.list(eval(parse(text=dictEnsembl2eg))[mappedkeys(eval(parse(text=dictEnsembl2eg)))])
      dictRefseq=as.list(eval(parse(text=dictRefseq))[mappedkeys(eval(parse(text=dictRefseq)))])
      
      #first: convert the ensembl in ENTREZ:
      totalEntrez=dictEnsembl2eg[as.character(transc$gene_id)]
      totalEntrez=lapply(totalEntrez,"[[",1)
      totalEntrez[sapply(totalEntrez, is.null)]=NA
      totalEntrez=unlist(totalEntrez)

      #now, using the ENTREZ, find symbols and refseq:
      totalSymbols=dictSymbol[as.character(totalEntrez)]
      totalRefseq=dictRefseq[as.character(totalEntrez)]
          
      totalSymbols=lapply(totalSymbols,"[[",1)    
      totalRefseq=lapply(totalRefseq,"[[",1)

      totalSymbols[sapply(totalSymbols, is.null)]=NA
      totalRefseq[sapply(totalRefseq, is.null)]=NA
          
      totalSymbols=unlist(totalSymbols) 
      totalRefseq=unlist(totalRefseq)
      #finally append to transcript table
      df=data.frame(gene_id=totalEntrez,
            symbol=totalSymbols,
            ensembl_id=as.character(transc$gene_id),
            refSeq_id=totalRefseq)
      elementMetadata(transc)=df


    }
    

  }
  recovery=sapply(1:ncol(df),function(x){tab=table(!is.na(df[,x]));   
          if(!is.na(tab["TRUE"]) & !is.na(tab["FALSE"])){
            return(round((tab["TRUE"]/(tab["FALSE"]+tab["TRUE"]))*100,2))
          }else if(!is.na(tab["TRUE"]) &is.na(tab["FALSE"]) ){
            return(100)
          }else if (is.na(tab["TRUE"]) &!is.na(tab["FALSE"]) ){
            return(0)
          }
        })
  names(recovery)=c("gene_id","symbol","ensembl_id","refSeq_id")
  print(recovery)
  res_=list(transc,recovery)
  names(res_)=c("transcripts","translation_recovery")
  return(res_)
}









#from characters of symbols or gene IDs and the range of transcripts,
#(and the maximum width threshold), returns the positions of transcripts 
#that were found and pass the threshold
findPositionFromGene<-function(genelist,annotatedrange,kindofID="entrez",thresh=250000) {
  #genelist: is a vector containing symbols or IDs
  # can be symbols, ENTREZ ids, refseq ids, ensembl ids
  #annotatedrange: are the genomic ranges of the transcripts for the current assembly,
  # complete with metatata containing all IDs and symbols (the result from extractFromDB function)
  # or any annotated ROI with the same columns
  #kindofID: whether genelist given to the function is an "entrez", "ensembl", "symbol" or "refseq"
  #thresh: threshold width for the transcripts or annotatedranges
  if(class("genelist")!="character"){
    stop("'genelist' must be a character vector of genes...")
  }
  if(length(genelist)==0){
    stop("'genelist' must not be 0 length")
  }
  if(class(annotatedrange)!="GRanges" ){
    stop("'annotatedrange' must be genomic ranges")
  }
  if(length(thresh)!=1 | ( class(thresh)!="numeric" & class(thresh)!="integer") ){
    stop("'thresh' must be a number")
  }
  if(thresh<=0){
    stop("'thresh' must be positive")
  }
  if(!(kindofID  %in% c("entrez","ensembl","symbol","refseq"))){
    stop("'kindofID' must be 'entrez','ensembl','symbol' or 'refseq'")
  }

  ids=genelist
  ids=toupper(ids)
  #extract element metadata from annotatedrange range: 
  #columns are gene_id symbol ensembl_id refSeq_id
  df=elementMetadata(annotatedrange)

  if(kindofID=="entrez"){
    coltouse=toupper(as.character(df$gene_id))
  }else if (kindofID=="ensembl"){
    coltouse=toupper(as.character(df$ensembl_id))
  }else if (kindofID=="symbol"){
    coltouse=toupper(as.character(df$symbol))
  }else if (kindofID=="refseq"){
    coltouse=toupper(as.character(df$refSeq_id))
  }

  
  totake=coltouse %in% ids
  lostinmatch=table(ids %in% as.character(coltouse))["FALSE"]
  if(is.na(lostinmatch)){
    lostinmatch=0
  }

  #find those genes not found at all in the annotatedranges
  geneslostinmatch=ids[is.na(match(ids,coltouse))]

  #filter annotated range based on a specific threshold width
  position_width=width(annotatedrange)<thresh

  totake2=which(totake &position_width)


  #if some match found (totake2), find which genes/isoforms were excluded
  if (length(totake2)>0 & !all(is.na(totake2))){
    lostpos=totake & (!position_width)
    lost_genes_ID=as.character(coltouse[which(lostpos)])
    lost_genes_ID=unique(lost_genes_ID)

    #catch only those whose alternative transcripts are lost:
    truly_lost=lost_genes_ID[!lost_genes_ID %in% coltouse[totake2]]
    truly_lost=unique(truly_lost)

    #genes with lost isoform but still alive in the filtering 
    isoform_lost=lost_genes_ID[!(lost_genes_ID %in% truly_lost)]

    #take lost_genes_ID (symbols) and truly_lost (symbols)
    #the first are isoforms lost, while the second are the complete lost
    pos_truly_lost=match(coltouse,truly_lost)
    pos_isoform_lost=match(coltouse,isoform_lost)

    #find symbols associated with those positions
    symbols_truly_lost=unique(as.character(df[!is.na(pos_truly_lost),]$symbol))
    symbols_isoform_lost=unique(as.character(df[!is.na(pos_isoform_lost),]$symbol))

    lres=list(totake2,geneslostinmatch,symbols_truly_lost,symbols_isoform_lost)
    names(lres)=c("totake","notfound","losttoolarge","lostisoformstoolarge")
    return(lres)
  }else{
    return(NULL)
  }

}


#given a range and TSSs (width=1), find distance and annotation. This function was taken from
#compEpiTools package
#modified from the previous distanceFromTSS2, because we have standardized the way in which TSS are build
#now, each TSS will have both ENTREZ ID, symbol, REFSEQ ID, ENSEMBL ID
#shuld be run 1) every time a change in database is made, for all ROIs 2) the first time of the database
#3) new ROI created, but only for : summit, union/intersections, resize; otherwise, inherit annotation from previus ROI
#tested on mac 4 core i5, 8Gb RAM; tested 50 reps of same ROI with 80k ranges each
#with mclapply 4 cores takes about 20 seconds. Feasible. 
distanceFromTSS3<-function(Object, Tss,criterion="midpoint") {
    #Object: GRanges from which calculate the distance with Tss and to be annotated
    #Tss: GRanges of length =1 (calculated with ENTREZ, symbol, refseq, ensembl, from function extractFromDB)
    #criterion: criterion to which callculate the distance; if not set, 
    #     the distance will be calculated from any of the points inside the Object,
    #     as below:
    #     TSS:    |             |                |
    #     Object:           ---------------
    #     => the distance will be 0
    #     TSS:    |             |                |
    #     Object:                    ---------  
    #     The distance will be:               ****   <- this one!          
    #     if set to "midpoint", the distance will be calculated
    #     from the midpoint of Object, as below:
    #     TSS:   |               |            |
    #     Object:                    ------|------ 
    #     The distance will be:             ***   <- this one! 
    if(class(Object)!="GRanges"){
      stop("Object must be a GRange")
    }
    if(class(Tss)!="GRanges"){
      stop("Tss must be a GRange")
    }
    if(!all(width(Tss))==1){
      stop("Tss must be transcription start site of width == 1")
    }

    if(length(Tss$gene_id)==0 | length(Tss$symbol)==0 | length(Tss$ensembl_id)==0 | length(Tss$refSeq_id)==0){
      stop("Tss must have gene_id, symbol, ensembl_id, refSeq_id column included...")
    }

    #resize Object to 1 (midpoint)
    #do not resize if you want the nearest from the boundaries of ranges and not the midpoint
    ######################################################
    if(criterion=="midpoint"){
      Object=resize(Object,width=1,fix="center")
    }

    TSSpos <- Tss
    #extract only TSS ranges with
    #gene_id
    #ensembl_id
    #symbol
    #refSeq_id
    
    ##could be avoided: in theory, they are all >0, we save computing time...
    # pos=sapply(TSSpos$gene_id,length)>0
    # TSSpos=TSSpos[pos]

    #match the ranges (slow)
    nearestInd <- nearest(Object, TSSpos)
    nonNAinds <- which(!is.na(nearestInd))
    nearestENTREZId= nearestENSEMBLId=nearestsymbol=nearestRefSeqId= rep(NA, length(Object))
    nearestDist <- rep(NA, length(Object))

    nearestENTREZId[nonNAinds] <- as.character(TSSpos[nearestInd[nonNAinds]]$gene_id)
    nearestENSEMBLId[nonNAinds] <- as.character(TSSpos[nearestInd[nonNAinds]]$ensembl_id)
    nearestsymbol[nonNAinds] <- as.character(TSSpos[nearestInd[nonNAinds]]$symbol)
    nearestRefSeqId[nonNAinds] <- as.character(TSSpos[nearestInd[nonNAinds]]$refSeq_id)

    suppressWarnings(nearestDist[nonNAinds] <- distance(Object[nonNAinds],TSSpos[nearestInd[nonNAinds]]))
    
    df=data.frame(gene_id=nearestENTREZId,symbol=nearestsymbol,ensembl_id=nearestENSEMBLId,refSeq_id=nearestRefSeqId,distance_from_TSS=nearestDist)
    return(df)  

}







#it returns the fraction of 1 and/or 2 that overlaps if a threshold on enrichment was set
test_ov<-function(ov,bam1_1,bam2_2,quantile) {
  #ov=object of overlap resulting from ov=findOverlaps(range1,range2)
  #bam1_1= coverage of range1 with bam1
  #bam2_2= coverage of range2 with bam2
  #quantile= quantile to threshold for the overlap
  if (is.null(bam1_1) & is.null(bam2_2)){
    stop("At least bam1_1 or bam2_2 must be set")
  }

  ov2=as.matrix(ov)
  #depending on if we have both bam1 and 2 or only one of them,
  #returns the fraction of overlap of 1 and 2
  #if both bams are present, the fractions are plotted when BOTH
  #ranges are filtered on their quantiles and not each of them separately


  if(!is.null(bam1_1) & is.null(bam2_2)){
    pos1=which(bam1_1>quantile(bam1_1,quantile))
    regions_total_1=length(pos1)
    regions_overlapping_1=length(  unique(ov2[,1][ov2[,1] %in% pos1]))
    perc1only=regions_overlapping_1/regions_total_1
    perc2only=NULL
    perc1combination=NULL
    perc2combination=NULL

  }else if(is.null(bam1_1) & !is.null(bam2_2)){
    pos2=which(bam2_2>quantile(bam2_2,quantile))
    regions_total_2=length(pos2)
    regions_overlapping_2=length(  unique(ov2[,2][ov2[,2] %in% pos2]))
    perc1only=NULL
    perc2only=regions_overlapping_2/regions_total_2
    perc1combination=NULL
    perc2combination=NULL

  }else if(!is.null(bam1_1) & !is.null(bam2_2)){
    pos1=which(bam1_1>quantile(bam1_1,quantile))
    regions_total_1=length(pos1)
    pos2=which(bam2_2>quantile(bam2_2,quantile))
    regions_total_2=length(pos2)  

    idx1=ov2[,1] %in% pos1
    idx2=ov2[,2] %in% pos2
    regions_overlapping_1_only=length(  unique(ov2[,1][idx1]))
    regions_overlapping_2_only=length(  unique(ov2[,2][idx2]))
    #find overlapping AREAS that pass the threshold. 2 are subjectHits, 1 are queryHits

    regions_overlapping_TOT= ov2[ idx1 &  idx2 ,,drop=FALSE]
    regions_overlapping_2=length(  unique(regions_overlapping_TOT[,2]))
    regions_overlapping_1=length(  unique(regions_overlapping_TOT[,1])) 
    perc1only=regions_overlapping_1_only/regions_total_1
    perc2only=regions_overlapping_2_only/regions_total_2
    perc1combination=regions_overlapping_1/regions_total_1
    perc2combination=regions_overlapping_2/regions_total_2
    
  }
  return(list(perc1only,perc2only,perc1combination,perc2combination))
}

shinyFileName<-function(obj) {
  #shinyFileChoose object
  return(unlist(obj[[1]])[length(unlist(obj[[1]]))])
}

shinyFilePath<-function(obj) {
  #shinyFileChoose object
  #returns complete path of the file selected
  partial=paste0(unlist(obj[[1]]),collapse=separation)
  if(.Platform$OS.type=="windows"){
    partial=paste("C:\\Users",partial,sep="")
  }
  return(partial)
}


#draw an empty plot
emptyplot<-function() {
  par(mar = rep(0, 4))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')  
}


# Function to plot color bar.
# adapted from John Colby (stackoverflow) http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=2), title='') {
  scale = (length(lut)-1)/(max-min)
  par(mar=c(2.1,1,2.1,1))
  plot(c(min,max), c(0,2), type='n', bty='n', xaxt='n', yaxt='n', ylab='', main=title,xlab="n",cex.lab=0.6)
  axis(3, ticks, las=1,cex=0.3)
  title(xlab="Quantile norm. rpm/bp", line=0.4, cex=0.4)
  for (i in 1:(length(lut)-1)) {
    x = (i-1)/scale + min
    rect(x,0,x+1/scale,2, col=lut[i], border=NA)
  }
}


color.bar2 <- function(lut, min, max=-min, margins=c(2,2,2,2),nticks=11, ticks=seq(min, max, len=2), title='') {
  scale = (length(lut)-1)/(max-min)
  par(mar=margins)
  plot(c(min,max), c(0,2), type='n', bty='n', xaxt='n', yaxt='n', ylab='', main=title,xlab="n",cex.lab=0.6)
  axis(3, ticks, las=1,cex=0.3)
  title(xlab="-log10 padj", line=0.4, cex=0.4)
  for (i in 1:(length(lut)-1)) {
    x = (i-1)/scale + min
    rect(x,0,x+1/scale,2, col=lut[i], border=NA)
  }
}


#this function checks the maximum number of bins in which a ROI can be divided into
checkMaxBins<-function(roiobject){
  if(class(roiobject)!="RegionOfInterest"){
    stop("roiobject must be a 'RegionOfInterest' object...")
  }
  rangeSelected=getRange(roiobject)
  rangeSelected=granges(rangeSelected)
  maxtoshow=min(width(rangeSelected))
  return(maxtoshow)
}


#this function calculate the pileup for each base pair for each range of the ROI
#is a little improvement of the GRbaseCoverage function from compEpiTools (lapply intead of for)
#a better strategy to save RAM could be to keep the norm. factor associated to each BAM, and use it when needed.
#in that case it is possible to store everything as integer, more than 2 billions values 
#https://stackoverflow.com/questions/23660094/whats-the-difference-between-integer-class-and-numeric-class-in-r
GRbaseCoverage2<-function(Object, signalfile,signalfileNorm=NULL,signalControl=NULL,signalControlSpike=NULL,multiplFactor=1e+06)
{
    #Object: genomic range in which calculate the base coverage
    # signalfile= path to the signal file file (the associated signal file index must be present in the same directory if bam file)
    #      or the bw/bigwig file 
    #signalfileNorm= path to the signal file (BAM) for which normalize the signalfile (can be the same file for library normalization,
    #             or the spike in bam, for spike in normalization)
    #signalControl=path to the signal file (BAM) of the input for spike-in normalization (ChIP-Seq). If this option
    #         is provided, signalControlSpike must be also provided
    #signalControlSpike=path to the signal file (BAM) of the spike-in in the input (ChIP-Seq). If this option is
    #         provided, also signalControl must be provided
    #multiplFactor= constant factor to multiply the final normalization coefficient, usually 1 million
    if (!is.character(signalfile)) {
      stop("signalfile has to be a file path of class character...")
    }
    if(!file.exists(signalfile)){
      stop("'signalfile' doesn't exist...")
    }

    if(!is.null(signalfileNorm)){
      if(!is.character(signalfileNorm)){
        stop("'signalfileNorm' has to be a file path of class character...")
      }
      if(!file.exists(signalfileNorm)){
        stop("'signalfileNorm' doesn't exist...")
      }
    }

    #check the existence of control and spikein
    if(!is.null(signalControl)){
      if(!is.character(signalControl)){
        stop("'signalControl' has to be a file path of class character...")
      }
      if(!file.exists(signalControl)){
        stop("'signalControl' doesn't exist...")
      }
    }

    if(!is.null(signalControlSpike)){
      if(!is.character(signalControlSpike)){
        stop("'signalControlSpike' has to be a file path of class character...")
      }
      if(!file.exists(signalControlSpike)){
        stop("'signalControlSpike' doesn't exist...")
      }
    }   
    #if provide signalControl you must provide also signalControlSpike
    if (xor(is.null(signalControl),is.null(signalControlSpike))){
      stop("If 'signalControl' is provided, also 'signalControlSpike' must be provided, and vice versa...")
    }

    #if you provide signalControlSpike and signalControl you must also have provided signalfileNorm
    if(!is.null(signalControl) & is.null(signalfileNorm)){
      stop("If 'signalControl' and 'signalControlSpike' are provided, also 'signalfileNorm' must be provided")
    }


    #select if bam:
    if(substring(signalfile,nchar(signalfile)-3,nchar(signalfile))==".bam"){
      BAMseqs <- names(scanBamHeader(signalfile)[[1]]$targets)
      matchingSeqs <- which(as.character(seqnames(Object)) %in% BAMseqs)
      if (length(matchingSeqs) == 0) 
          return(sapply(width(Object), function(x) rep(0, x)))
      #use ss parameter for splitting strands
      #already integer!
      covList=as.list(bamCoverage(bampath=signalfile, gr=Object[matchingSeqs], verbose=FALSE))
      # param <- ApplyPileupsParam(which = Object[matchingSeqs],what = "seq")
      # pileupFiles=PileupFiles(signalfile)
      # coverage <- applyPileups(pileupFiles, FUN = function(x) x,param = param)
      # rm(pileupFiles)
      # widths <- width(Object[matchingSeqs])
      # covList <- list()
      # starts <- start(Object[matchingSeqs])
      # covList=lapply(1:length(Object[matchingSeqs]), function(i) {
      #     covx <- coverage[[i]]
      #     cvec <- rep(0, widths[i])
      #     inds <- covx$pos - starts[i] + 1
      #     #sometims max(inds) is > than length of cvec. => NAs in some positions. Why?
      #     cvec[inds] <- colSums(covx$seq)
      #     return(cvec)
      # })
      # covList<-lapply(covList,function(i)as.integer(i))



      ####################################################################
      # normalization
      ####################################################################

      if(!is.null(signalfileNorm)){
          param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
          nreads <- countBam(signalfileNorm, param = param)$records
          norm_factor=multiplFactor/nreads

          if(!is.null(signalControl)) {
            nreadCtr=countBam(signalControl, param = param)$records
            nreadCtrSpike=countBam(signalControlSpike, param = param)$records
            norm_factor=norm_factor* nreadCtrSpike/  nreadCtr         
          }
          print (paste("Norm. factor=",norm_factor))  
          #covList <- lapply(covList, function(x) (x*norm_factor) )
      }else{
      	norm_factor=1
      }

      

      if (length(matchingSeqs) < length(Object)) {
          coverageTot <- sapply(width(Object), function(x) list(rep(0, 
              x)))
          coverageTot[matchingSeqs] <- covList
      }
      else coverageTot <- covList

      rm(covList)
      return(list(coverageTot,norm_factor))


    #select if wig:
    }else if (substring(signalfile,nchar(signalfile)-2,nchar(signalfile))==".bw" | tolower(substring(signalfile,nchar(signalfile)-6,nchar(signalfile)))==".bigwig"){
    
      #open wig file only for names and lengths of chromosomes, to find
      #those matching with "Object"
      print ("    Coverage of WIG file...")
      chromosomes=seqlevelsInUse(Object)
      #just open the wig file, without opening all chromosome content. This serves just to know 
      #how long is each chromosome
      grfake=GRanges(Rle(chromosomes[1]),IRanges(1,1))
      wig=import(signalfile,which=grfake,as = 'Rle')
      chrswig=lengths(wig)
      rm(wig)
      common_chromosomes=intersect(chromosomes,names(chrswig))
      #initialize 0s for all ranges (if wig do not have a chr, the remaining are 0s)
      widths=width(Object)
      coverageTot=lapply(1:length(Object),function(i){integer(widths[i])})
      for (i in 1:length(common_chromosomes)){
        #length(as.character(seqnames(Object))) should be equal to length(Object) =>
        #find positions of Object that have that chromosome
        #opening chromosome by chromosome
        pos=as.character(seqnames(Object))==common_chromosomes[i]
        wigtempchr=import(signalfile,which=Object[pos],as = 'Rle')
        #count coverage for that chromosome in wig in the Object ranges in the same chromosome
        counts=Views(unlist(wigtempchr[[common_chromosomes[i]]]),ranges(Object[pos]))
        rm(wigtempchr)
        #extract the base coverage for all those ranges for ith chromosome.
        #used "as.integer" to save space: all numbers are pure integers
        tmp=viewApply(counts, as.integer)
        
        if(length(counts)==1){
          tmp=list(as.vector(tmp))
        }
        if(any(class(tmp)!="list") ){
          ##check whether the order is correct
          tmp=lapply(seq_len(ncol(tmp)), function(i) as.integer(tmp[,i]))
        }
        rm(counts)
        ##IMPORTANT##
        #############
        #change here if numbers can not be integer, but real
        #############
        #tmp<-lapply(tmp,function(i)as.integer(i))
        ########################################################
        coverageTot[pos] <- tmp


      }
      #return coverage and norm factor (in case of WIG, simply 1) to keep the same result structure 
      #for the function

      return(list(coverageTot,1))      
    }else{
      stop("Error in retrieving file/extension...")
    }

}















#from a range and one of its associated BAM (tipically the BAM from which this range
#derives), the range containing the summit (width=1) is returned by this function
#this function was taken and adapted from "GRcoverageSummit" function of compEpiTools package
summitFromBaseCoverage<-function(Object,baseCoverageOutput) {

  #Object: GRanges from which derive the summits
  #baseCoverageOutput: list of output of GRbaseCoverage2 function

  if(class(Object)!="GRanges"){
    stop("'Object' must be a GRange...")
  }
  if(class(baseCoverageOutput)!="list"){
    stop("'baseCoverageOutput' must be a list")
  }

  blcov=baseCoverageOutput

  #if many maxPos in very different indexes, maybe the peak is not sharp
  #and doesn't have a good shape...
  maxPos <- lapply(blcov, function(x) which(x == max(x)))
  starts <- start(Object)
  maxPos2=lapply(1:length(maxPos),function(i) {maxPos[[i]]+starts[i] - 1})
  maxPos=maxPos2

  maxPosL <- sapply(maxPos, length)
  if (max(maxPosL) > 1) {
      maxPos[maxPosL > 1] <- sapply(maxPos[maxPosL > 1], function(x) sample(x, 
          1))
      maxPos <- unlist(maxPos)
  }
  else maxPos <- unlist(maxPos)
  start(Object) <- maxPos
  end(Object) <- maxPos
  return(Object)
}


#given a range and TSSs (width=1), find distance. This function was taken from
#compEpiTools package
distanceFromTSS2<-function (Object, Tss,criterion="midpoint") 
{
    #Object: GRanges from which calculate the distance with Tss
    #Tss: GRanges of length =1 
    #criterion: criterion to which callculate the distance; if not set, 
    #     the distance will be calculated from the boundaries of Object
    #     otherwise, if set to "midpoint", the distance will be calculated
    #     from the midpoint of Object
    if(class(Object)!="GRanges"){
      stop("Object must be a GRange")
    }
    if(class(Tss)!="GRanges"){
      stop("Tss must be a GRange")
    }
    if(!all(width(Tss))==1){
      stop("Tss must be transcription start site of width == 1")
    }

    if(length(Tss$gene_id)==0){
      stop("Tss must have gene_id column included...")
    }

    TSSpos <- Tss
    #extract only TSS ranges with gene_id
    pos=sapply(TSSpos$gene_id,length)>0
    TSSpos=TSSpos[pos]
    #resize Object to 1 (midpoint)
    #do not resize if you want the nearest from the boundaries of ranges and not the midpoint
    ######################################################
    if(criterion=="midpoint"){
      Object=resize(Object,width=1,fix="center")
    }
    ######################################################
    nearestInd <- nearest(Object, TSSpos)
    nonNAinds <- which(!is.na(nearestInd))
    nearestId <- rep(NA, length(Object))
    nearestId[nonNAinds] <- unlist(TSSpos[nearestInd[nonNAinds]]$gene_id)
    nearestDist <- rep(NA, length(Object))
    suppressWarnings(nearestDist[nonNAinds] <- distance(Object[nonNAinds],TSSpos[nearestInd[nonNAinds]]))
    df=data.frame(gene_id=nearestId,distance_from_TSS=nearestDist)
    return(df)
}



#functions that takes GRbaseCoverage and output GRcoverageInbins output
#taken from GRcoverageInbins function from compEpiTools
# and re-implemented in Rcpp
makeMatrixFrombaseCoverageCPP <- cxxfunction(signature(GRbaseCoverageOutput='List',Nbins='integer',Snorm="integer"), plugin='Rcpp', body = '  
     Rcpp::List xlist(GRbaseCoverageOutput); 
     int nbins = Rcpp::as<int>(Nbins);
     int snormbool=Rcpp::as<int>(Snorm);
     int n = xlist.size(); 
     //std::vector<double> res(n);  
     double res;
     //define final matrix output: ncol=nbin, nrow=n
     //double mat[n][nbins];
     Rcpp::NumericMatrix mat( n , nbins );

     std::vector<int> lengths(n);

     for(int i=0; i<n; i++) {     
         SEXP ll = xlist[i]; 
         Rcpp::NumericVector y(ll);  
         int m=y.size();   
         //divide m (size) by the number of bins: how many elements to be summed for each bin?
         int goodpart=m/nbins;
         //printf("size of a bin: %d\\n",goodpart);

         // for each bin, sum elements in each bin
         for(int k=0; k<nbins; k++){
          res=0;
          for(int j=k*goodpart; j<(k+1)*goodpart; j++){     
              res=res+y[j]; 
          }   
          mat(i,k)=res;
         }

         //populate the array of lengths if Snorm=TRUE
         lengths[i]=m;
     }
     //if Snorm=TRUE (Snorm>0), divide each matrix value for m
     if(snormbool>0){
      for(int i=0; i<n; i++){
        for(int k=0; k<nbins; k++){
          mat(i,k) = mat(i,k)/lengths[i];
        }
      }
     }
       
   return mat;  
') 

#wrapper for makeMatrixFrombaseCoverageCPP function
makeMatrixFrombaseCoverage <-function(GRbaseCoverageOutput,Nbins,Snorm=FALSE,norm_factor=1) {
    # GRbaseCoverageOutput: output from GRbaseCoverage2 function
    # Nbins: the number of bins to ivide the coverage for each range into
    # Snorm: whether to normalize the coverage for each bin for the length of the range (TRUE/FALSE)
    if(Snorm==TRUE){
      norm=1
    }else{
      norm=0
    }
    #here introduce normalization step. We could have done easyly in the CPP function,
    #but it's fast anyway
    result=makeMatrixFrombaseCoverageCPP(GRbaseCoverageOutput,Nbins,Snorm=norm)
    if(norm_factor!=1){
      #here the result should still be a matrix
    	result <- result*norm_factor
    }
    return(result)
}

#cut the transcripts ranges in specific indexes and sum the coverage inside the internal part (remaining range),
#the input is always a list of kind "GRbaseCoverageOutput", from GRbaseCoverage2 function
#StartingPositions and EndingPositions are the array of positions for each element of the GRbaseCoverageOutput list
#it returns an integer vector, that is the sums of the elements of the list within the cut part (be careful to the 0-index of C!)
cutAndSumTranscriptsCPP <-cxxfunction(signature(GRbaseCoverageOutput='List',StartingPositions="vector",EndingPositions="vector"), plugin='Rcpp', body = '  
     Rcpp::List xlist(GRbaseCoverageOutput); 
     Rcpp::IntegerVector starts(StartingPositions);
     Rcpp::IntegerVector ends(EndingPositions);

     int n = xlist.size(); 
     double res;
     //define final result
     //std::vector<double> results(n);
     Rcpp::NumericVector results(n);
     for(int i=0; i<n; i++) {     
         SEXP ll = xlist[i]; 
         Rcpp::NumericVector y(ll);  
         //int m=y.size();   
         res=0;
         // for each bin, sum elements in each bin
         for(int k=starts[i]; k<= ends[i]; k++){
          //printf("cumulative: %f\\n",res);
          //sum the cut element
          res=res+y[k-1]; 
         }

         //populate the final array with the cut sum
         results[i]=res;
     }
     //return the numeric array as result
     return (results); 
') 
#wrapper for cutAndSumTranscriptsCPP function
cutAndSumTranscripts <-function(GRbaseCoverageOutput,StartingPositions,EndingPositions,norm_factor=1) {
    # GRbaseCoverageOutput: output from GRbaseCoverage2 function
    #here introduce normalization step. We could have done easyly in the CPP function,
    #but it's fast anyway
    result=cutAndSumTranscriptsCPP(GRbaseCoverageOutput,StartingPositions,EndingPositions)
    if(norm_factor!=1){
      #here the result should still be a matrix
      result <- result*norm_factor
    }
    return(result)
}


#OLD function (matrix frm base coverage). New functions implemented in Rcpp (much faster)
# #functions that takes GRbaseCoverage and output GRcoverageInins output
# #to parallelize externally among multiple bams
# #taken from GRcoverageInbins
# makeMatrixFrombaseCoverage <-function(GRbaseCoverageOutput,Nbins,Snorm=FALSE) {
#   #GRbaseCoverageOutput: output from GRbasecoverage function (list ov coverages of each single bp)
#   #Nbins:number of bins to split
#   #Snorm: normalize for length of grange

#   if (class(Nbins)!="integer" & class(Nbins)!="numeric"){
#     stop("Nbins must be an integer or numeric...")
#   }
#   if (class(GRbaseCoverageOutput)!="list"){
#     stop("GRbaseCoverageOutput must be a list...")
#   }
#   if (class(GRbaseCoverageOutput[[1]])!="numeric" & class(GRbaseCoverageOutput[[1]])!="integer"){
#       stop("GRbaseCoverageOutput elements must be numeric...")
#   }    

#   widths=unlist(lapply(GRbaseCoverageOutput,length))
#   if (min(widths)<Nbins){
#     print("Warning: Nbins is greater than the length of some ranges, they will be NA")
#   }
#   #split
#   coverageMat <- matrix(NA, length(GRbaseCoverageOutput), Nbins)
#   binsize <- floor(widths / Nbins)
#   remains = widths%%Nbins
#   #find end position to keep (excluding remains)
#   endpos_list=widths-remains

#   #check if ROI is "fixed size". In this case, try with operation on matrix, much faster
#   logicfixed=length(table(widths))==1
#   #create factor specific for each bin (will be used in tapply later)
  

#   #too much concnetrated on the last bin. If you want to distribute the
#   #remaining of the division you should do:
#   # binsize_integer <- rep( floor(widths / Nbins),Nbins)
#   # binsize_remain<- widths %% Nbins
#   #and sum, for each binsize_integer elements, +1 to the last binsize_remain elements
#   #and below substitute binsize with binsize[[i]]

#   #integer part of the division for Nbins-1 times, the last
#   #one is current position -> end

#   #if elements of the coverage list are of different lengths (variable peak size for example)
#   #different elements have different binSize
#   if (logicfixed==FALSE){
    
#     print ("not fixed calc")
#     #Nbins threshold: 150: above this, the tapply method is more convenient
#     #with an "average" ROI
#     if (Nbins>150){
#       factor_binning=lapply(1:length(widths),function (x) {rep(1:Nbins,each=binsize[x])})
#       provv_calc=lapply(1:length(GRbaseCoverageOutput),function(gr) {
#               element=GRbaseCoverageOutput[[gr]][1:endpos_list[gr]]
#               return(tapply(element,INDEX=factor_binning[[gr]],sum))
#             }) 
#       coverageMat <- matrix(unlist(provv_calc), ncol = Nbins, byrow = TRUE)
#     }else{
#       for(bin in 1:Nbins) {
#         startPos <- binsize * {bin - 1} + 1
#         #if(bin == Nbins) endPos <- widths
#         endPos <- startPos + binsize - 1

#         # for(gr in 1:length(GRbaseCoverageOutput)) coverageMat[gr,bin] <-
#         #     sum(GRbaseCoverageOutput[[gr]][startPos[gr]:endPos[gr]])
#         row_sum=sapply(1:length(GRbaseCoverageOutput), function(gr) {
#             return(sum(  GRbaseCoverageOutput[[gr]][startPos[gr]:endPos[gr]]   ))
#         })
#         coverageMat[,bin]=row_sum
#       }  
#     }
    
#   #if, however, ROI is fixed size, GRbaseCoverageOutput can be treated as matrix:
#   #operations are muc more efficient (~8 times faster!)  
#   }else{
#     print("fixed calc")
#     #transform list in matrix 

#     matToUse <- matrix(unlist(GRbaseCoverageOutput), ncol = widths[1], byrow = TRUE)
#     #calculate remove the rest of last bin from the matrix columns (thrown). 
#     remain=widths[1]%%Nbins
#     #ncolumns=ncol(matToUse)
#     #startCol=ncolumns-remain
#     #matToUse=matToUse[,-c( startCol  ,endCol )]
#     #for each submatrix, do the colsum and add to the final coverageMat
#     binsize=binsize[1]
#     for(bin in 1:Nbins) {
#       startPos <- binsize * {bin - 1} + 1
#       #if(bin == Nbins) endPos <- widths
#       endPos <- startPos + binsize - 1
#       submatrix=matToUse[,startPos:endPos]
#       #colSums. As integer approximates (wrong!)
#       coverageMat[,bin]=as.numeric(rowSums(submatrix))
#     }
#   }

#   if(Snorm){
#     coverageMat=apply(coverageMat,2,function(i){i/widths})
#   }
#   coverageMat[Nbins>widths]<- NA
#   return(coverageMat)
# } 









#countOverlapsInBin function from compEpiTools package
countOverlapsInBins<-function (query, subject, nbins,strandspecific=FALSE) 
{
    #query: GRanges in which calculate the overlaps
    #subject: GRanges to use for overlap computation on query
    #nbins: the number of bins to divide the overlap into
    #strandspecific: if consider the strand specificity (+ strand ranges will overlap only with 
            #+ strand ranges) or not (all strands will be *, therefore + and - strands can overlap)
    if (!is(subject, "GRanges")) 
        stop("subject has to be of class GRanges ..")
    if (!is.numeric(nbins)) 
        stop("nbins has to be of class numeric ..")
    countMat <- matrix(0, length(query), nbins)
    binsize <- floor(width(query)/nbins)
    starts <- start(query)
    if(strandspecific==TRUE){
      strands <- strand(query)
    }else{
      strands <- rep("*",length(query))
    }
    for (bin in 1:nbins) {
        startPos <- starts + binsize * (bin - 1)
        if (bin == nbins) 
            endPos <- end(query)
        else endPos <- startPos + binsize - 1
        queryBin <- GRanges(seqnames = seqnames(query), ranges = IRanges(start = startPos, 
            end = endPos),strand=strands)
        countMat[, bin] <- countOverlaps(queryBin, subject, maxgap = 0L, 
             type = "any")
    }
    countMat[countMat > 1] <- 1
    return(countMat)
}



#cluster a matlist with hierarchical method according to parameters (clust. methods, indexes of list...). 
#return matrix ordered and the order
#hclust function must be overwritten with the hclust function of fastcluster package
#CANBERRA distance method gives problems on apparently good matrixes
clusterMatrix<-function(matlist,distmethod,clustmethod,clustinds){
  #matlist: list of matrixes for which the rows must be clustered
  #distmethod: method to calculate the distance
  #clustmethod: method for cluster calculation
  #clustinds: indexes rpresenting the position of the matrixes inside matlist 
  #         that will drive the clustering. All the other matrixes will be reordered accordingly
  if (!is.null(clustinds)){
    if(class(clustinds)!="integer"){
      stop("'clustinds' must be an integer or NULL...")
    }
  }
  if(class(distmethod)!="character" | class(clustmethod) != "character"){
    stop("'distmethod' and 'clustmethod' must be characters...")
  }
  if(class(matlist)!="list"){
    stop("'matlist' must be a list...")
  }

  for(i in 1:length(matlist)){
    if(class(matlist[[i]])!="matrix"){
      stop("each element of 'matlist' must be a matrix...")
    }
  }

  mat <- NULL
  for (i in 1:length(matlist)) {
    mat <- cbind(mat, matlist[[i]])
  }
  if(!is.null(clustinds)){
    clmat = NULL      
    for (i in clustinds){
      clmat <- cbind(clmat, matlist[[i]])
    } 
    NAcounts <- apply(mat, 1, function(x) length(which(is.na(x))))
    NAinds <- which(NAcounts > ncol(mat) * 0.3)
    if (length(NAinds) == nrow(mat)) 
        return(NULL)
    if (length(NAinds) > 0) {
        mat <- mat[-NAinds,,drop=FALSE]
        clmat <- clmat[-NAinds,,drop=FALSE]
    }
    clustobject <- hclust(dist(clmat,method=distmethod), method = clustmethod)
    ord= clustobject$order
    #order the matrix
    mat =mat[ord,,drop=FALSE]
  }else{
    ord=1:nrow(mat)
    clustobject=NULL
  }
  resul=list(mat,ord,clustobject)
  names(resul)=c("mat","ord","clustobject")
  return(resul)
}


#cluster a matlist using k-means clustering. Number of clusters decided a priori
clusterMatrixKmeans<-function(matlist,clustinds,numberclusters,startingpoints,iter){
  #matlist: list of matrixes for which the rows must be clustered
  #clustinds: indexes rpresenting the position of the matrixes inside matlist 
  #         that will drive the clustering. All the other matrixes will be reordered accordingly
  #numberclusters: number of clusters for K-means
  #startingpoints: number of starting points for the kmeans clustering
  #iter: number of iteration for k-mean clustering

  #k-means is not deterministic! make it reproducible across different runs
  set.seed(123)
  if (!is.null(clustinds)){
    if(class(clustinds)!="integer"){
      stop("'clustinds' must be an integer or NULL...")
    }
  }
  if(class(matlist)!="list"){
    stop("'matlist' must be a list...")
  }
  for(i in 1:length(matlist)){
    if(class(matlist[[i]])!="matrix"){
      stop("each element of 'matlist' must be a matrix...")
    }
  }
  if (class(numberclusters)!="numeric" & class(numberclusters)!="integer"){
    stop("'numberclusters' must be a number...")
  }
  if(class(startingpoints)!="numeric" & class(startingpoints)!="integer"){
    stop("'startingpoints' must be a number...")
  }else{
    if(startingpoints<=0){
      stop("'startingpoints' must be > 0...")
    }
  }
  if(class(iter)!="numeric" & class(iter)!="integer"){
    stop("'iter' must be a number...")
  }else{
    if(iter<=0){
      stop("'iter' must be > 0...")
    }
  }

  #initialize the big matrx from all the matrixes
  mat <- NULL
  for (i in 1:length(matlist)) {
    mat <- cbind(mat, matlist[[i]])
  }

  #if indexes of clustering are not null,
  if(!is.null(clustinds) & numberclusters>0 & numberclusters<434){
    ###SHOULD implemented: if nclustering is < than single data point
    ###should be very rare... in this case, try to do unique of rowsums of the matrix
    #and check whether the length of the unique values is < than number of clustering.
    #if so, no clustering....
    #prepare the matrix composed only by the matrixes that will drive the clustering
    clmat = NULL      
    for (i in clustinds){
      clmat <- cbind(clmat, matlist[[i]])
    } 
    #treat NAs (if all NAs, return NULL)
    NAcounts <- apply(mat, 1, function(x) length(which(is.na(x))))
    NAinds <- which(NAcounts > ncol(mat) * 0.3)
    if (length(NAinds) == nrow(mat)) 
        return(NULL)
    if (length(NAinds) > 0) {
        mat <- mat[-NAinds, ,drop=FALSE]
        clmat <- clmat[-NAinds, ,drop=FALSE]
    }

    #check clmat: no more clusters than 2^n, where n is all the possible combinations 
    #of the rows of the matrix
    maxclusterallowed=table(!duplicated(clmat))["TRUE"]
    if(numberclusters > maxclusterallowed ){
      numberclusters=maxclusterallowed
      print(paste("warning: max clustering allowed is",maxclusterallowed))
    }
    #HERE modify or add kmeans parameters, or alternatively set them by the user from UI
    #if very fast and you want more precise, increase iter.max (10 is default)
    #nstart is the initial random configuration (default of original function=1)
    clustobject <- kmeans(clmat, centers = numberclusters,nstart=startingpoints,iter.max=iter)
    groups=clustobject$centers
    clusters= clustobject$cluster
    ord=order(clusters)
    #you can improve ordering/separation of clusters between them according to groups matrix
    #keep more sparate those more different
    #order the matrix
    mat =mat[ord,,drop=FALSE]
  
  #otherwise cluster object is null and the order is not given (no clustering, no ranking)
  }else{
    ord=1:nrow(mat)
    clustobject=NULL
  }
  resul=list(mat,ord,clustobject)
  names(resul)=c("mat","ord","clustobject")
  return(resul)
}




# function for generating new ROIs starting from existing ROIs and 
# using exclusion/inclusion criteria with other ROIs
generateROI<-function(selectedlist,selectedfix=NULL,overlaplist,notoverlaplist,method,criterion1,criterion2,bamlist,minbp=1,strandSpecific=FALSE){
  #selectedlist: list of GR to start from (primary ROI or custom ROI)
  #selectedfix: the fix for the ROI selected. If selectedlist is composed by only one ROI, preserve the 
  # fix of the old ROI (maybe fix not simmetric to the center and strand specific!!), otherwise, if
  #  multiple ROIs are going to be joined in a union or intersection, the fix will be the simple midpoint
  #  in fact, you cannot predict where the original fix was if multiple ROIs are fused together
  #overlaplist: list of GR. selectedlist must overlap with these GRs.
  #notoverlaplist: list of GR. selectedlist must NOT overlap with these GRs.
  #method: combine method for starting selectedlist: union or intersection
  #strandSpecific: TRUE or FALSE. If TRUE, - will overlap with other - or * (but not +)
  #                 if FALSE, same range lablled as - and + can overlap with each other
  if (class(method)!="character"){
    stop("method must be a character ('intersection' or 'union' or NULL)...")
  }

  if (length(selectedlist)>1){
    if(method !="union" & method!= "intersection"){
      stop("method must be 'union' or 'intersection' if more than two input GR given...")
    }    
  }

  #check wether selectedfix properly set (put it only when input GR list is only 1)
  if(length(selectedlist)==1){
    if(is.null(selectedfix)){
      stop("'selectedfix' must be the Fix of the only GR given in input...")
    }
  }

  if(length(selectedlist)>1){
    if(!is.null(selectedfix)){
      stop("'selectedfix' must be NULL if multiple GR given in input...")
    }    
  }

  if (!is.null(selectedfix)){
    if(class(selectedfix)!="GRanges"){
      stop("'selectedfix' must be of class GRanges...")
    }  

    if(length(selectedfix)!=length(selectedlist[[1]])){
      stop("'selectedfix' must have same length of the unique GR given in input...")
    }
  }  
  

  if (criterion1!="stringent" & criterion1!="permissive" & criterion1!="allofthem"){
    stop("criterion1 must be 'stringent' or 'permissive'...")
  }

  if (criterion2!="intersection" & criterion2!="union"){
    stop("criterion2 must be 'intersection' or 'union'...")
  }

  if (class(bamlist)!="list"){
    stop("'bamlist' must be a list of baseCoverage of BAM files")
  }


  if (!is.null(overlaplist)){
    if (class(overlaplist[[1]])!="GRanges"){
      stop("overlaplist does not contain all GRanges...")
    }
  }

  if (!is.null(notoverlaplist)){
    if (class(notoverlaplist[[1]])!="GRanges"){
      stop("overlaplist does not contain all GRanges...")
    }
  }

  if(class(minbp)!="numeric" & class(minbp)!="integer"){
    stop("'minbp' parameter must be a number...")
  }

  #find starting point (either single range or union/intersection of GRanges, if more than 1 GR provided)
  if (length(selectedlist)>1){
    #will be union or intersection of multiple ROIs... so BAMlist is cleared
    finalbam=list()
    if (method=="intersection"){
      startingpoint=selectedlist[[1]]
      for (k in 2:length(selectedlist)){
        if (length(startingpoint)==0){
          break
        }
        startingpoint=intersect(startingpoint,selectedlist[[k]])
      }


      #startingpoint=Reduce(intersect,selectedlist)
    }else{
      startingpoint=Reduce(union,selectedlist)
    } 
    #reconstruct fix according to the middle point of the new,fused, GR:
    #now selectedfix should not be NULL in either cases (startingpoint==1 or >1)
    selectedfix=resize(startingpoint,width=1,fix="center") 
  }else{
    #otherwise only one input
    startingpoint=selectedlist[[1]]
    #...and BAM file list is not altered
    finalbam=bamlist
  }

  #starting point strand: consider or not for the overlap? (strandSpecific?)
  #if starting point is already * this will not change anything
  tokeep_strand=strand(startingpoint)
  if(!strandSpecific){
    #make * everywhere, we want all overlaps, regardless of the strand:
    #startingpointprovv=startingpoint
    strand(startingpoint)="*"
  }else{
    #startingpointprovv=startingpoint
  }

  if (!is.null(overlaplist)){
    if (criterion1=="stringent"){
      #find GR for positive overlap:
      #must overlap with ALL, so intersect
      overlapwith=Reduce(intersect,overlaplist)  
      positiveov= countOverlaps(startingpoint,overlapwith,minoverlap=minbp)
      startingpoint=startingpoint[positiveov>0]  
      if (length(finalbam)>0){
        #keep elements of bam list, according to the calculated overlap
        finalbam=lapply(finalbam,function(k) {k[positiveov>0]})
      }
      selectedfix=selectedfix[positiveov>0]
      tokeep_strand=tokeep_strand[positiveov>0]
    }

    if (criterion1=="permissive"){
      for (i in 1:length(overlaplist)){
        ov=countOverlaps(startingpoint,overlaplist[[i]],minoverlap=minbp)
        startingpoint=startingpoint[ov>0]
        selectedfix=selectedfix[ov>0]
        tokeep_strand=tokeep_strand[ov>0]
        finalbam=lapply(finalbam,function(k) {k[ov>0]})
      }
    }

    if (criterion1=="allofthem"){
      overlapwith=Reduce(union,overlaplist)
      positiveov= countOverlaps(startingpoint,overlapwith,minoverlap=minbp)
      startingpoint=startingpoint[positiveov>0]
      selectedfix=selectedfix[positiveov>0]
      tokeep_strand=tokeep_strand[positiveov>0]
      finalbam=lapply(finalbam,function(k) {k[positiveov>0]})
    }

  }

  #starting point strand: consider or not for the overlap? (strandSpecific?)
  #if starting point is already * this will not change anything
  if(!strandSpecific){
    #make * everywhere, we want all overlaps, regardless of the strand:
    #startingpointprovv=startingpoint
    strand(startingpoint)="*"
  }else{
    #startingpointprovv=startingpoint
  }


  if (!is.null(notoverlaplist)){
    if (criterion2=="union"){
      #find GR for negative overlap:
      #I don't want it to overlap with ANY bp of notoverlaplist
      notoverlapwith=Reduce(union,notoverlaplist)   
      negativeov= countOverlaps(startingpoint,notoverlapwith,minoverlap=minbp)
      startingpoint=startingpoint[negativeov==0]
      selectedfix=selectedfix[negativeov==0]
      tokeep_strand=tokeep_strand[negativeov==0]
      finalbam=lapply(finalbam,function(k) {k[negativeov==0]})    
    }else{
      notoverlapwith=Reduce(intersect,notoverlaplist)   
      negativeov= countOverlaps(startingpoint,notoverlapwith,minoverlap=minbp)
      startingpoint=startingpoint[negativeov==0]  
      selectedfix=selectedfix[negativeov==0]
      tokeep_strand=tokeep_strand[negativeov==0]
      finalbam=lapply(finalbam,function(k) {k[negativeov==0]})          
    }
 
  }
  strand(startingpoint)=tokeep_strand

  return(list(startingpoint,finalbam,selectedfix))
}




#this function calculates the partial correlation between 2 selected features
#given a matrix
plotpcor<-function (mat,idx) {
  #mat: matrix/df with m features and n observations
  #idx: vector numeric of length 2, for indexes of cols in the matrix to consider
  #   for the partial correlation. The other elements will be considered as the confounding factors
  
  if (class(mat)!="matrix" & class(mat)!="data.frame"){
    stop("mat must be a matrix or a data.frame...")
  }
  if (length(idx) != 2){
    stop("idx must be a numeric vector of length 2...")
  }
  if (class(idx)!= "numeric" & class(idx)!= "integer"){
    stop("idx is not numeric...")
  }else{
    if(max(idx)>ncol(mat)){
      stop("max idx cannot be greater than number of mat columns...")
    }
  }
  # if (idx[1]==idx[2]){
  #   stop("two indexes are the same!...")
  # }

  #transform mat in data.frame
  mat=as.data.frame(mat)
  name_x=colnames(mat)[idx[2]]
  name_y=colnames(mat)[idx[1]]
  colnames(mat)[idx[1]]="y"
  colnames(mat)[idx[2]]="x"
  #exclude other variables
  others_lessx=mat[,-c(idx[2])]
  others_lessy=mat[,-c(idx[1])]
  y=mat[,idx[1]]
  x=mat[,idx[2]]
  modelx=lm(x~.,others_lessy)
  modely=lm(y~.,others_lessx)
  resx=modelx$residuals
  resy=modely$residuals

  return(list(resx,resy))
}



#clean chromosomes (_random, _alt, ChrUn...)
#this function will return the same starting GRange without _random, _alt, ChrUn

# cleanChromosomes<-function(range){
#   #range is the GenomicRange in input
#   if (class(range)!="GRanges"){
#     stop("'range' must be of class GRange...")
#   }
#   chrs=as.character(seqnames(range))
#   pos_random=grepl("_random$",chrs)
#   pos_alt=grepl("_alt$",chrs)
#   pos_fix=grepl("_fix$",chrs)
#   pos_ChrUn=grepl("^chrUn_",chrs)
#   tokeep=!(pos_random | pos_alt |pos_ChrUn |pos_fix)
#   return(range[tokeep])
# }


#translation nomenclature GRanges: use them in specific sections:
#when open a BED/GTF/GFF, clean chromosomes and keep standard. Convert to UCSC
# (promoters, BSgenome and other DBs work using UCSC)

#when associating WIG/BAM. If all are zeroes, it means that BAM/WIG are in ENCODE format.
# in that case, temporary convert to ENCODE the ROI to be associated! Here the lentgh
# should be the same, because chromosomes have been cleaned when opening the input from
# the beginning
#output can be smaller than the input
convertNomenclatureGR <-function(range, to="UCSC") {
  #range is the GenomicRange in input
  if (class(range)!="GRanges"){
    stop("'range' must be of class GRange...")
  }
  if(to !="UCSC" & to != "NCBI"){
  	stop("'to' must be either 'UCSC' or 'NCBI'...")
  }
  #remove non-standard chromosomes:
  gr=keepStandardChromosomes(range,pruning.mode="coarse")
  newStyle <- mapSeqlevels(seqlevels(gr), to)
  gr <- renameSeqlevels(gr, newStyle)
  return(gr)
}


# #verify that the output of coverage functions from ROIs (can be either a list or a matrix,
# #if ranges are fixed length) are all zeroes
# verifyzerocov<-function(covresult) {
#   if (class(covresult)=="matrix"){
#     if(all(covresult==0)){
#       return(TRUE)
#     }else{
#       return(FALSE)
#     }
#   }else if (class(covresult)=="list"){
#     #extremely RAM expensive and inefficient
#     covresult2=unlist(covresult)
#     # function "all" extremely RAM expensive and inefficient
#     if(all(covresult2==0)){
#       return(TRUE)
#     }else{
#       return(FALSE)
#     }
#   }else{
#     stop("'covresult' must be a matrix or a list...")
#   }
# }


#extremely efficient implementation of zero check of coverages in CPP. No more RAM peaks
#much faster
verifyzerocov<-cxxfunction(signature(covresult='List'), plugin='Rcpp', body = '  
     Rcpp::List xlist(covresult); 

     int n = xlist.size(); 
     //bool allzeros=true;
     //loop through all ranges in list
     int tempsum=0;
     Rcpp::LogicalVector ALLzeros(1,true);
     for(int i=0; i<n; i++) {  
       //printf("LOOP %i NEW ",i);   
         Rcpp::NumericVector y(xlist[i]); 
         int rangelen=y.size();
         // start from 1-based position (maybe Rcpp is 1-based)
         for(int k=1; k<= rangelen; k++){
          tempsum+=y[k];
         }

         if(tempsum>0){
          ALLzeros(0)=false;
          break;
         }
     }
     
     //ALLzeros(0)=allzeros;
     //return boolean. "true" if 
     return (ALLzeros); 
')




#from a GRange, extract ranges with a specific pattern
#if bothstrands (TRUE), each range will be considered in both directions.
#sequence CATTCC will be looked at ------> and <------- in the same range (as if the range is *)
# => CATTCC can be read in <---- and ----->
#example:
#    -----------CATTCC---->
#   or
#    <--CCTTAC-------------
#if bothstrand=FALSE and in Subject we have the strand information, we will look only at the current strand:
#example (for a range labelled as strand +):
#   (+) --------CATTCC------->
#   but not:
# [[(-) <--CCTTAC-------------]] <<- NO!
#because THAT range had (+) strand, therefore only + direction was invstigated for the pattern

extractPattern<-function(Subject,BSgenomeDB,pattern,bothstrands=TRUE){
  
  #Subject is the range from which patterns shuld be xtracted. If null,
  #   the pattern will be extracted from the entire genome
  #BSgenomeDB is the BSgenome database for the sequence (UCSC) with the correct genome assembly (for example, mm9/mm10), 
  # (BSgenome.Mmusculus.UCSC.mm9)
  #pattern is a string with character corresponding to the pattern to search, IUPAC codes available
  #bothstrands: if Subject has strands specified, search the motif only in the strand specified
  # otherwise, search the patterns in + strand even if that range is in the - strand
  
  ###############################################################
  ###############################################################
  ###############################################################
  #checks
  ###############################################################
  ###############################################################
  ###############################################################

  if(!is.null(Subject)){
    if(class(Subject)!="GRanges"){
      stop("'Subject' must be a GRanges object or NULL, for entire genome pattern search")
    }    
  }

  if(class(BSgenomeDB)!="BSgenome"){
    stop("'BSgenomeDB' should be a BSgenome database")
  }
  if(class(pattern)!="character"){
    stop("'pattern' must be a string representing the pattern to search")
  }
  splitted_pattern=strsplit(pattern,split="")[[1]]
  check=splitted_pattern=="." | splitted_pattern=="-" |
    splitted_pattern=="A" | splitted_pattern=="C" | splitted_pattern=="G" | splitted_pattern=="T" |
    splitted_pattern=="U" | splitted_pattern=="R" | splitted_pattern=="Y" | splitted_pattern=="S" |
    splitted_pattern=="W" | splitted_pattern=="K" | splitted_pattern=="M" | splitted_pattern=="B" |
    splitted_pattern=="D" | splitted_pattern=="H" | splitted_pattern=="V" | splitted_pattern=="N" 
  if(!all(check)){
    stop("some not valid letters (IUPAC code) are present in 'pattern'")
  }

  ###############################################################
  ###############################################################
  ###############################################################
  strToSearch=DNAString(pattern)

  
  if(!is.null(Subject)){
    #if strand *, + and -, if + only +, if - only -
    #sometimes, resulting motifs spans contiguous regions. Treat them separately
    
    #clean Subjcet from strange chromosome range (_alt,_fix_random.....)
    #Subject=cleanChromosomes(Subject)
    Subject=keepStandardChromosomes(Subject,pruning.mode="coarse")
    #########
    #should extract seq only for seqnames of GRange (Subject) found in DB
    ourseq=names(table(seqnames(Subject)))
    totalseq=seqnames(BSgenomeDB)
    tosave=ourseq[ourseq%in%totalseq]
    Subject=Subject[seqnames(Subject)%in%tosave]
    if (length(Subject)==0){
      return(NULL)
    }
    #if (length(Subject)==0){
    # 
    #}
    #########
    if( unique(as.character(strand(Subject)))[1]=="*" | bothstrands){
      #check both strands
      pos_pos=rep(TRUE,length(Subject))
      pos_neg=pos_pos   
      strand(Subject)="*"  
    }else{
      #is strand-specific and we have strand information => analyze only strand orinted patterns
      pos_pos=as.character(strand(Subject))=="+"
      pos_neg=as.character(strand(Subject))=="-"
    }


    
    chromosome_names=as.character(seqnames(Subject))

    #to avoid errors, keep only the seqnames in common between the two (BSgenomeDB,Subject)
    #For hg38, Subject seqnames are much more!

    yseq=getSeq(BSgenomeDB,Subject) # if - strand, reverse complement is done by default
    starts=start(Subject)
    ends=end(Subject)

    #split strand
    yseq_positive=yseq[pos_pos]
    starts_positive=starts[pos_pos]
    ends_positive=ends[pos_pos]
    chromosome_names_positive=chromosome_names[pos_pos]
    yseq_negative=yseq[pos_neg]
    #because if one of the two, strand was * and sequence was positive. We have to reverse complement.
    #otherwise, getting seq from - strand is already reversedcomplement
    if(unique(as.character(strand(Subject)))[1]=="*" | bothstrands){
      yseq_negative=reverseComplement(yseq_negative)
    }
    starts_negative=starts[pos_neg]
    ends_negative=ends[pos_neg]
    chromosome_names_negative=chromosome_names[pos_neg]


    #plus
    mtcp=vmatchPattern(strToSearch,yseq_positive,with.indels=TRUE, fixed="subject")
    stt_p=startIndex(mtcp)
    edd_p=endIndex(mtcp)  
    notvalid_pos_p=sapply(stt_p,is.null)
    chr_toput_p=chromosome_names_positive[!notvalid_pos_p]
    stt_p=stt_p[!notvalid_pos_p]
    edd_p=edd_p[!notvalid_pos_p]
    lengths=sapply(stt_p,length)
    chr_toput=rep(chr_toput_p,times=lengths)
    stt_p=unlist(stt_p)
    edd_p=unlist(edd_p)
    #start/end of + must be summed to starts of original Subject (offset)
    starts_global_p=starts_positive[!notvalid_pos_p]
    starts_global_p=rep(starts_global_p,times=lengths)
    newstarts_p=starts_global_p+stt_p-1
    newends_p=starts_global_p+edd_p-1
    #create range with + strand
    patternrange_p=GRanges(Rle(chr_toput),IRanges(newstarts_p,newends_p),strand=rep("+",length(newstarts_p)))

    #minus
    mtcm=vmatchPattern(strToSearch,yseq_negative,with.indels=TRUE, fixed="subject")
    stt_m=startIndex(mtcm)
    edd_m=endIndex(mtcm)
    notvalid_pos_m=sapply(stt_m,is.null)
    chr_toput_m=chromosome_names_negative[!notvalid_pos_m]
    stt_m=stt_m[!notvalid_pos_m]
    edd_m=edd_m[!notvalid_pos_m]
    lengths=sapply(stt_m,length)
    chr_toput=rep(chr_toput_m,times=lengths)
    stt_m=unlist(stt_m)
    edd_m=unlist(edd_m)
    #start/end of - must be the end of the Subject - start and -end of the motif
    ends_global_m=ends_negative[!notvalid_pos_m]
    ends_global_m=rep(ends_global_m,times=lengths)
    newstarts_m=ends_global_m-edd_m+1
    newends_m=ends_global_m-stt_m+1
    #create range with - strand
    patternrange_m=GRanges(Rle(chr_toput),IRanges(newstarts_m,newends_m),strand=rep("-",length(newstarts_m)))

    #join ranges
    patternrange=c(patternrange_p,patternrange_m)
  }else{
    #here, Subject is null, so pattern will be searched through all the genome.
    #WARNING: time consuming! (several seconds)
    patternrange=vmatchPattern(strToSearch,BSgenomeDB,with.indels=TRUE, fixed="subject")
  }

  

  #if no pattern, length of this range is ==0
  return(patternrange)

}

###############################################################
###############################################################
###############################################################
#GSEA functions (read gmt file of MSigDB and perform GOs)
###############################################################
###############################################################
###############################################################


readGMT<-function(fileName) {
  if(!file.exists(fileName)){
    stop("'fileName' does not exist")
  }
  gmt=readLines(fileName)
  #split each elment by "\t". Take the first (term) and from the 3rd (genes)
  splitted=sapply(gmt,strsplit,split="\t")
  number=unname(sapply(splitted,length))
  genes=lapply(1:length(splitted),function(i){return(splitted[[i]][3:number[i]])})
  genes=unlist(genes,use.names=FALSE)
  number=number-2
  terms=rep(unname(sapply(splitted,"[[",1)),times=number)
  df=data.frame(ont=terms,gene=genes)
  df=unique(df)
  return(df)
}



#the code for this function was inspired from enricher_internal function in DOSE package 
#(author Guangchuang Yu, http://guangchuangyu.github.io)
GOcalc<-function(gene,terms,minsize=10,maxsize=500,padj_method="BH") {
  #gene: list of genes
  #term: a data frame with term->gene associations, resulting from readGMT function
  #       or a combination of different catgories
  #minsize and maxsize: filter of size of genesets to consider
  #padj_method: method to calculate the padjusted from hypergeometric test
  if(class(gene)!="character"){
    stop("'gene' must be a character vector with gene symbols")
  }

  if(minsize<0 | maxsize<0){
    stop("'minsize' and 'maxsize' must be > 0")
  }
  if(!(padj_method%in%c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"))){
    stop("wrong padj method given")
  }

  gene=unique(gene)
  gene=gene[!is.na(gene)]

  if(length(gene)==0){
    stop("gens in query is empty...")
  }

  #take only those queried:
  pos=which(as.character(terms$gene)%in%gene)

  if(sum(pos)==0){
    #no genes provided matches those in genesets
    return(NULL)
  }

  
  queried_terms=terms[pos,]
  keyval=split(as.character(queried_terms[,2]),as.character(queried_terms[,1]))
  #retrieve universe (all genes in the term GMT provided)
  #if universe is custom, put here the universe as INTERSECTION.
  #in case, intersect also with keyval
  allgenes=unique(as.character(terms$gene))
  # if(!is.null(background)){
  #   allgenes=intersect(allgenes,background)
  # }
  
  #take all genes of terms specified found with at last one of our genes:
  allkeyval=split(as.character(terms[,2]),as.character(terms[,1]))
  pos=names(allkeyval)%in%names(keyval)
  allkeyval=allkeyval[pos]
  
  #filter keyval and allkeyval for min and max
  sizes=sapply(allkeyval,length)
  idx= sizes>=minsize & sizes<=maxsize
  #if no idx, no geneset satisfies criteria
  if(sum(idx)==0){
    return(NULL)
  }
  
  keyval=keyval[idx]
  allkeyval=allkeyval[idx]

  IDnames=names(keyval)
  
  #do hypergeometric test
  query=sapply(keyval, length)
  totalBG=sapply(allkeyval, length)
  genes_selected=sum(gene%in%as.character(terms$gene))
  
  #calculate pvalues with hypergeometric model
  HyperGeomTests=sapply(1:length(query),function(i){
                                  phyper(q=query[i]-1, m=totalBG[i], n=length(allgenes)-totalBG[i], k=genes_selected, lower.tail = FALSE, log.p = FALSE)
                                }
                        )
  
  tocalc=lapply(1:length(query),  function(i){
                lst=list(unname(query[i][1]),unname(genes_selected),unname(totalBG[i]),unname(length(allgenes)))
                names(lst)=c("overlap","query","geneset","universe")
                return(lst)
  })
  tocalc=do.call("rbind",tocalc)
  

  ## gene ratio: selected genes / genes for each geneset
  ## background ratio: ratio between all genes in a geneset (regardless the query) and the universe (union of all genes of all genesets used)

  gene_ratio = unlist(tocalc[,1])/unlist(tocalc[,2])
  background_ratio= unlist(tocalc[,3])/unlist(tocalc[,4])
  
  HyperGeomTests_padj=p.adjust(HyperGeomTests, method=padj_method)

  geneID = sapply(keyval, function(i) paste(i, collapse="/"))
  
  finaldf=data.frame(Term=IDnames,
                    Genes=geneID,
                    overlap=unlist(tocalc[,1]),
                    query=unlist(tocalc[,2]),
                    geneset=unlist(tocalc[,3]),
                    universe=unlist(tocalc[,4]),
                    gene_ratio=gene_ratio,
                    background_ratio=background_ratio,
                    pval=HyperGeomTests,
                    padj=HyperGeomTests_padj)
  #rank according to padj:
  ord=order(finaldf$padj)
  finaldf=finaldf[ord,]
  
  if(nrow(finaldf)==0){
    return(NULL)
  }
  return(finaldf)

}



#function for filter matrix of results from GO (padj, generatio....)
filterGOres<-function(GOres,padjthresh,generatiothresh,topN) {
  #inputs: the result matrix from GO analysis (serverROI) with thresholds from input GUI
  #filter padj (at least one ROI must be <thresh)

  mat=GOres
  pos_tokeep=rep(TRUE,nrow(mat))

  partdfpadj=mat[,grepl("_padj$",colnames(mat)),drop=FALSE]
  
  if(min(dim(mat))>0){
    tokeeppadj=apply(partdfpadj,1,function(i){any(i<padjthresh)})
    pos_tokeep=pos_tokeep&tokeeppadj
  }

  #filter gene ratio
  if(min(dim(mat))>0){
    partdfgeneratio=mat[,grepl("_gene_ratio$",colnames(mat)),drop=FALSE]
    tokeepgeneratio=apply(partdfgeneratio,1,function(i){any(i>=generatiothresh)})
    pos_tokeep=pos_tokeep&tokeepgeneratio     
  }
  


  #filter topN elements
  if(min(dim(mat))>0){
    #filter with previous filters:
    partdfpadjfilt=partdfpadj[pos_tokeep,]
    minvals=apply(partdfpadj,1,min)
    minvalsfilt=minvals[pos_tokeep]
    #keep the first N positions of the already filtered matrix (top 10 for example of the matrix
    #already filtered for the other parameters)
    minvalsfilt=sort(minvalsfilt)
    if(length(minvalsfilt)>topN){
      #in this case we have more elements than topN => filter
      val=minvalsfilt[topN]
      tokeeptopN=minvals<=val
      pos_tokeep=pos_tokeep&tokeeptopN

    }else{

    }
  }

  return(pos_tokeep)
}






#DEFINE CLASS ROI:
# this class will have 
# - a range derived from BED files or from prvious ROIs
# - a name
# - a list of BAM files associated (list, that derives from GRbaseCoverage2 function) - to be considered obsolete
# - a fixed point: this point is a range that represents a position inside the ranges of the ROI
#   that could be the midpoint (for "general" ROIs), the TSS (for promoters), the summit etc...
#   this is generated automatically when creating new ROIs. It becomes the summit when a user
#   will create a ROI summit from another ROI
# - annotation: a data frame rpresenting the annotation of the ROI (associated gene IDs/symbols...)
# - a flag: this indicates which kind of ROI a user is working on. This can be "general" ROI, a transcript, a TSS....
# - a source: this is a list of steps (character vectors) that represent how the current ROI had been generated

setClass("RegionOfInterest",
  slots=c(  name = "character", 
            range = "GenomicRanges",
            fixed="GenomicRanges",
            BAMlist="list",
            annotation="data.frame",
            flag="character",
            source="list"),
  prototype=list(range=GRanges(Rle("chr1"),IRanges(0,0 )),annotation=data.frame(),
                  fixed=GRanges(Rle("chr1"),IRanges(0,0 )),flag="normalFlag",source=list() ) #for default values

  #validity=function(object) if check for validity during construction
)




setGeneric(name="getName",def=function(object) {standardGeneric("getName")} )
setMethod(f="getName",signature="RegionOfInterest",definition=function(object){return(object@name)})
setGeneric(name="getRange",def=function(object) {standardGeneric("getRange")} )
setMethod(f="getRange",signature="RegionOfInterest",definition=function(object){return(object@range)})
setGeneric(name="getBAMlist",def=function(object) {standardGeneric("getBAMlist")} )
setMethod(f="getBAMlist",signature="RegionOfInterest",definition=function(object){return(object@BAMlist)})
setGeneric(name="getLength",def=function(object) {standardGeneric("getLength")} )
setMethod(f="getLength",signature="RegionOfInterest",definition=function(object){return(length(object@range))})
setGeneric(name="getWidth",def=function(object) {standardGeneric("getWidth")} )
setMethod(f="getWidth",signature="RegionOfInterest",definition=function(object){return(width(object@range))})
setGeneric(name="getAnnotation",def=function(object) {standardGeneric("getAnnotation")} )
setMethod(f="getAnnotation",signature="RegionOfInterest",definition=function(object){return(length(object@annotation))})
setGeneric(name="getFixed",def=function(object) {standardGeneric("getFixed")} )
setMethod(f="getFixed",signature="RegionOfInterest",definition=function(object){return(object@fixed)})
setGeneric(name="getFlag",def=function(object) {standardGeneric("getFlag")} )
setMethod(f="getFlag",signature="RegionOfInterest",definition=function(object){return(object@flag)})
setGeneric(name="getSource",def=function(object) {standardGeneric("getSource")} )
setMethod(f="getSource",signature="RegionOfInterest",definition=function(object){return(object@source)})





#set functions must be used, to change the object in: object=set...(object,...)
#setname
setGeneric(name="setName",def=function(object,name) {standardGeneric("setName")} )
setMethod(f="setName",signature="RegionOfInterest",def=function(object,name){object@name <- name;return(object)})
#setrange
setGeneric(name="setRange",def=function(object,range) {standardGeneric("setRange")} )
setMethod(f="setRange",signature="RegionOfInterest",def=function(object,range){object@range <- range;return(object)})

#setFix
setGeneric(name="setFix",def=function(object,range) {standardGeneric("setFix")} )
setMethod(f="setFix",signature="RegionOfInterest",def=function(object,range){
  if(!all(width(range))==1){
    stop("Fix must be range of length 0...")
  }
  if (class(object)!="RegionOfInterest"){
    stop("'object' must be of class 'RegionOfInterest'...")
  }
  if (class(range)!="GRanges"){
    stop("'range' must be of class 'GenomicRanges'...")
  }
  object@fixed <- range ;return(object)
})


#reset the BAMlist
setGeneric(name="resetBAMlist",def=function(object) {standardGeneric("resetBAMlist")} )
setMethod(f="resetBAMlist",signature="RegionOfInterest",def=function(object){object@BAMlist <- list(); return(object)})

#set the BAMlist
setGeneric(name="setBAMlist",def=function(object,bamlist) {standardGeneric("setBAMlist")} )
setMethod(f="setBAMlist",signature="RegionOfInterest",def=function(object,bamlist){object@BAMlist <- bamlist; return(object)})

##obsolete, we use the function directly 
#add new list (using GRbasecoverage2 function) in the BAMlist of the range
setGeneric(name="cover",def=function(Object,signalfile,signalfileNorm,signalControl,signalControlSpike) {standardGeneric("cover")} )
setMethod(f="cover",signature="RegionOfInterest",def=function(Object,signalfile,signalfileNorm,signalControl,signalControlSpike){
  #name is the basename(signalfile). signalfile is the complete path of the enrichment file
  if (!file.exists(signalfile)){
    stop("signalfile file doesn't exist...")
  }

  if (class(Object)!="RegionOfInterest"){
    stop("'Object' must be of class 'RegionOfInterest'...")
  }
  rang=getRange(Object)
  cov=GRbaseCoverage2(Object=rang, signalfile=signalfile,signalfileNorm=signalfileNorm,signalControl=signalControl,signalControlSpike=signalControlSpike, multiplFactor=1e+06)
  return(list(cov[[1]],cov[[2]]))

})


#if strand info is present, invert negative strand (for profiles and heatmaps)
#to see the asimmetry of profile. If strand is *, keep as +
setGeneric(name="unifyStrand",def=function(object) {standardGeneric("unifyStrand")} )
setMethod(f="unifyStrand",signature="RegionOfInterest",def=function(object){
  if (class(object)!="RegionOfInterest"){
    stop("'object' must be of class 'RegionOfInterest'...")
  }
  
  rang=getRange(object)
  pos_negative= as.character(strand(rang))=="-"
  #get BAMs of the RegionOfInterest
  bamlist=getBAMlist(object)
  #invert range in pos_negative positions
  bamlist2=list()
  
  if(length(bamlist)>0 & any(pos_negative)){

    for(i in 1:length(bamlist)){
      element_i=bamlist[[i]]
      cov_neg=element_i[pos_negative]
      cov_pos=element_i[!pos_negative]
      #invert. VERY SLOW. 
      cov_adj=lapply(cov_neg,rev)
      element_i[pos_negative]=cov_adj
      bamlist2[[i]]=element_i
    }  

    names(bamlist2)=names(bamlist)
    object@BAMlist <- bamlist2  
  }
   
  return(object)

})



#separate + and - strand if present. * is considered as +
#useful for the digital heatmap
setGeneric(name="splitStrand",def=function(object) {standardGeneric("splitStrand")} )
setMethod(f="splitStrand",signature="RegionOfInterest",def=function(object){
  if (class(object)!="RegionOfInterest"){
    stop("'object' must be of class 'RegionOfInterest'...")
  }
  
  rang=getRange(object)
  pos_negative= as.character(strand(rang))=="-"
  rang_neg=rang[pos_negative]
  rang_pos=rang[!pos_negative]
  object@range <- c(rang_pos,rang_neg)  
 
  return(object)

})




#use it if duplicated ROI to remove (for example promoters of alternative transcripts, but duplicated)
setGeneric(name="uniqueROI",def=function(object) {standardGeneric("uniqueROI")} )
setMethod(f="uniqueROI",signature="RegionOfInterest",def=function(object){
  if (class(object)!="RegionOfInterest"){
    stop("'object' must be of class 'RegionOfInterest'...")
  }
  
  rang=getRange(object)
  fix=getFixed(object)
  pos_dup=!duplicated(rang)
  newrange=rang[pos_dup]
  newfix=fix[pos_dup]

  bamlist=getBAMlist(object)
  #invert range in pos_negative positions
  bamlist2=list()

  if(length(bamlist)>0){
    for(i in 1:length(bamlist)){
    element_i=bamlist[[i]]
    #remove duplicated positions
    bamlist2[[i]]=element_i[pos_dup]
    }
    names(bamlist2)=names(bamlist)
    object@BAMlist<-bamlist2  
  }
  object@range <- newrange
  object@fixed<-newfix
  
  return(object)

})





#convert from/to UCSC/NCBI nomenclature for GRanges of a ROI
setGeneric(name="convertNomenclatureROI",def=function(Object,To) {standardGeneric("convertNomenclatureROI")} )
setMethod(f="convertNomenclatureROI",signature="RegionOfInterest",def=function(Object,To){
  if(To !="UCSC" & To != "NCBI"){
  	stop("'To' must be either 'UCSC' or 'NCBI'...")
  }
  if (class(Object)!="RegionOfInterest"){
    stop("'Object' must be of class 'RegionOfInterest'...")
  }
  rang=getRange(object)
  #call the convertNomenclatureGR function to convert GRange nomenclature
  newgr=convertNomenclatureGR(range=rang,to=To)
  Object@range <- newgr
  return(Object)
})



#annotate, using transcripts from a specific txdb database. TO BE IMPLEMENTED
#modify range of the object
setGeneric(name="annotate",def=function(object,transcripts) {standardGeneric("annotate")} )
setMethod(f="annotate",signature="RegionOfInterest",def=function(object,transcripts){})