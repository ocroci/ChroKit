
# Define global shiny server routine
shinyServer(function(input, output,session) {
  #############################################################################
  #############################################################################
  ### this block of code is executed only when session is closed. Useful if ChroKit is put inside a cluster

  # session$onSessionEnded(function() {

  # })
  #############################################################################
  #############################################################################

  ### these lines are needed for surfing the remote filsystem (if different from the machine
  ### in which shiny is running). These functions are contained in shinyFiles library
  #set directory in which take the files. If files taken, input$file will be set to filename
  shinyFileChoose(input, id='file', roots=c(wd=rootdir), session=session, restrictions=system.file(package='base'))
  #for searching the filesystem for BAM files
  shinyFileChoose(input, id='fileBAM', roots=c(wd=rootdir), session=session, restrictions=system.file(package='base'),
            ,filetypes=c('bw','BW','bigWig','bigwig','bam','BAM'))
  #for loading the session rds file
  shinyFileChoose(input, id='loadenv',  roots=c(wd=rootdir),session=session, restrictions=system.file(package='base'),filetypes=c('rds'))
  #for loading Genelist file
  shinyFileChoose(input, id='fileGENELISTS', roots=c(wd=rootdir), session=session, restrictions=system.file(package='base'))
  #save button for rda
  shinyFileSave(input, id="saveWork", roots=c(wd=rootsavedir), session=session)
  
  ################################################################################
  ################################################################################
  ################################################################################
  # REACTIVE VALUES DEFINITIONS
  ################################################################################
  ################################################################################
  ################################################################################

  ### thesee variables, when defined, are keep during the entire working session
  ### and the are not eliminated inside each single reactive block


  #variables regarding bed files and GRanges derived
  BEDvariables<-reactiveValues( tempBED=NULL,tempBEDname=NULL,opened=FALSE,sfn=NULL,completepath=NULL) # currentBED?
       
  # #variables regarding Granges (derived from BEDs)
  # GRvariables<- reactiveValues(listGR=list(),names=NULL,selected=NULL,selectedname=NULL,
  #                               listfordensity=NULL,colorsfordensity=NULL,listselected=NULL,listselectednames=NULL,changed=1)

  #variables regarding ROIs (built at runtime). primary ROIs are GRvariables$listGR!!
  #to be filled with RegionOfInterest objects. names are inside the object
  ROIvariables<-reactiveValues(listROI=list(),selected=NULL,selectedname=NULL,listfordensity=NULL,colorsfordensity=NULL,listselected=NULL,listselectednames=NULL,changed=1)
  
  #variables regarding bam files. all paths to BAMfiles
  BAMvariables<-reactiveValues(listBAM=c())

  #variables regarding databases. Databases of TXDB, dbSYMBOL, available (pre-calculated), downloadable (pre-calculated) or actual
  DATABASEvariables<-reactiveValues(currentTXDB=NULL,currentORG=NULL,currentDICTsymbol2id=NULL,currentDICTid2symbol=NULL,loadedTXDB=NULL,loadedORG=NULL,
                                    currentASSEMBLY=NULL,
                                    availASSEMBLIES=availASSEMBLIES,
                                    missingASSEMBLIES=missingASSEMBLIES)

  #variables regarding the heatmap (matrix, if put text or not...)
  heatvariables<-reactiveValues(matlist=NULL,completematrixes=NULL,binsAnalogHeat=NULL,BAMsForAnalogHeat=NULL,
                                sampleRandomAnalogHeat=NULL,ROIsForAnalogHeat=NULL,widthfix=NULL)

  #matrix for correlation/partial correlation to use in the observer
  corvariables<-reactiveValues(portionlist_boxes=NULL,is.density=FALSE)

  #variables regarding log of all actions (what the user did since the beginning)
  logvariables<-reactiveValues(msg=list(),temporary_GMTstorage=temporary_GMTstorage,temporary=list())

  #temporary variables for plotting
  toplot<-reactiveValues(viewDistributionPieSingleEval=list(),cmp=list(),digital=list(),analogic=list(),
                          gadgetanalogic=list(),profileAndBoxes=list(),viewROI=list(),dynamics=list(),GOvariables=list())

  #temporary variables for saving
  tosave<-reactiveValues(datatableROI=NULL,genelistROI=NULL,genelistROIwindow=NULL)
  #variables to store raw read counts of base coverage (to be stored as int 4 bytes) and corresponding normalization factors
  Enrichlist<-reactiveValues(rawcoverage=NULL,decryptkey=NULL,normfactlist=NULL)

  ################################################################################
  ################################################################################
  ################################################################################
  #BED FILE SELECTION 
  ################################################################################
  ################################################################################
  ################################################################################


  source(file.path("servers","serverBED.R"),local=TRUE)$value


  ################################################################################
  ################################################################################
  ################################################################################
  #BAM FILE SELECTION 
  ################################################################################
  ################################################################################
  ################################################################################


  source(file.path("servers","serverBAM.R"),local=TRUE)$value


  ################################################################################
  ################################################################################
  ################################################################################
  #LIBRARY (txdb) SELECTION
  ################################################################################
  ################################################################################
  ################################################################################
  source(file.path("servers","serverTXDB.R"),local=TRUE)$value

  ################################################################################
  ################################################################################
  ################################################################################
  #GENELISTS SELECTION 
  ################################################################################
  ################################################################################
  ################################################################################
  #source(file.path("servers","serverGENELISTS.R"),local=TRUE)$value


  ################################################################################
  ################################################################################
  ################################################################################
  # ROI management (from GRanges)
  ################################################################################
  ################################################################################
  ################################################################################

  source(file.path("servers","serverROI.R"),local=TRUE)$value

  ################################################################################
  ################################################################################
  ################################################################################
  # GENOMICS (deep analyses)
  ################################################################################
  ################################################################################
  ################################################################################

  source(file.path("servers","serverGENOMICS.R"),local=TRUE)$value



  ################################################################################
  ################################################################################
  ################################################################################
  # SAVE/LOAD (to save and load work)
  ################################################################################
  ################################################################################
  ################################################################################

  source(file.path("servers","serverSAVELOAD.R"),local=TRUE)$value







  ################################################################################
  ################################################################################
  ################################################################################
  # UPDATE BED UI graphics
  ################################################################################
  ################################################################################
  ################################################################################
  source(file.path("servers","updateBEDui.R"),local=TRUE)$value
  ################################################################################
  ################################################################################
  ################################################################################
  # UPDATE BAM UI graphics
  ################################################################################
  ################################################################################
  ################################################################################
  source(file.path("servers","updateBAMui.R"),local=TRUE)$value
  ################################################################################
  ################################################################################
  ################################################################################
  # UPDATE TXDB UI graphics
  ################################################################################
  ################################################################################
  ################################################################################
  source(file.path("servers","updateTXDBui.R"),local=TRUE)$value
  ################################################################################
  ################################################################################
  ################################################################################
  # UPDATE MANIPULATEROI UI graphics
  ################################################################################
  ################################################################################
  ################################################################################
  source(file.path("servers","updateMANIPULATEROIui.R"),local=TRUE)$value




  ################################################################################
  ################################################################################
  ################################################################################
  # UPDATE ROI UI tab 
  ################################################################################
  ################################################################################
  ################################################################################
  source(file.path("servers","updateROIui.R"),local=TRUE)$value

  ################################################################################
  ################################################################################
  ################################################################################
  # UPDATE GENOMICS UI tab 
  ################################################################################
  ################################################################################
  ################################################################################
  source(file.path("servers","updateGENOMICSui.R"),local=TRUE)$value



  ################################################################################
  ################################################################################
  ################################################################################
  # UPDATE SAVE/LOAD UI tab 
  ################################################################################
  ################################################################################
  ################################################################################
  #source(file.path("servers","updateSAVELOADui.R"),local=TRUE)$value


  ################################################################################
  # update logs 
  ################################################################################

  #show log of the actions made since the beginning of the program:
  output$showlogs<-renderText({
    paste(rev(logvariables$msg),collapse="\n")
  })

  output$showcurrentASSEMBLY<-renderText({
    if(!is.null(DATABASEvariables$currentASSEMBLY)){
      if(DATABASEvariables$currentASSEMBLY!=FALSE){
        paste("Current assembly: <b>",DATABASEvariables$currentASSEMBLY,"</b>",sep="")
      }else{
        "No genome assembly loaded..."
      }      
    }else{
      "No genome assembly loaded..."
    }
  })


  # #every time the ROI list is changing or the enrichments are changins, show the current memory used
  # #using gc
  # output$showRAMusageGC<-renderText({
  #   ROIvariables$listROI
  #   Enrichlist$rawcoverage
  #   mem=gc()
  #   #extract value of vectors used in Mb
  #   val=unname(mem[,2][2])
  #   if(length(val)>0){
  #     paste("<p style='font-size:20px'>RAM usage (Mb): <br><b>",val,"</b></p>",sep="")
  #   }else{
  #     paste("")
  #   }

  # })


  #every time the ROI list is changing or the enrichments are changins, show the current memory used
  #using gc
  output$showRAMusageGC<-renderText({
    ROIvariables$listROI
    Enrichlist$rawcoverage
    toplot$analogic$matrixes_processed
    toplot$digital$matrixes_processed
    mem=gc()
    #extract value of vectors used in Mb
    val=unname(mem[,2][2])
    val=round(mem_used()/1000000,1)
    if(length(val)>0){
      paste("<p style='font-size:20px'>RAM usage (Mb): <br><b>",val,"</b></p>",sep="")
    }else{
      paste("")
    }

  })



})


