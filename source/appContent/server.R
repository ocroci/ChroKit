
# Define global shiny server routine
shinyServer(function(input, output,session) {
  #############################################################################
  #############################################################################
  ### this block of code is executed only when USER variable (argv [2]) is set
  ### when the session is closed.
  ### in the architecture in which multiple user can use the program has been implemented
  ### for example, a PBS cluster. If the application is run normally (without USER)
  ### argument, this piece of code is not executed
  # #exit function when session is closed (execute the script closingSession.sh)
  session$onSessionEnded(function() {
    if (!is.null(USER)){
      #print(paste("bash /hpcnfs/data/BA/public_html/ocroci/closingSession.sh",USER))
      system(paste("bash /hpcnfs/data/BA/public_html/ocroci/closingSession.sh",USER))
    }
  })
  #############################################################################
  #############################################################################

  ### these lines are needed for surfing the remote filsystem (if different from the machine
  ### in which shiny is running). These functions are contained in shinyFiles library
  #set directory in which take the files. If files taken, input$file will be set to filename
  shinyFileChoose(input, id='file', roots=c(wd=rootdir), session=session, restrictions=system.file(package='base'))
  #for searching the filesystem for BAM files
  shinyFileChoose(input, id='fileBAM', roots=c(wd=rootdir), session=session, restrictions=system.file(package='base'))
  #for loading the session rds file
  shinyFileChoose(input, id='loadenv',  roots=c(wd=rootdir),session=session, restrictions=system.file(package='base'))
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
  logvariables<-reactiveValues(msg=list(),temporary_GMTstorage=temporary_GMTstorage)

  #temporary variables for plotting
  toplot<-reactiveValues(viewDistributionPieSingleEval=list(),cmp=list(),digital=list(),analogic=list(),
                          gadgetanalogic=list(),profileAndBoxes=list(),viewROI=list(),dynamics=list(),GOvariables=list())

  #temporary variables for saving
  tosave<-reactiveValues(datatableROI=NULL,genelistROI=NULL,genelistROIwindow=NULL)



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
  # UPDATE GENELISTS UI graphics
  ################################################################################
  ################################################################################
  ################################################################################
  #source(file.path("servers","updateGENELISTSui.R"),local=TRUE)$value




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
  source(file.path("servers","updateSAVELOADui.R"),local=TRUE)$value


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

})















# ################################################################################
# ################################################################################
# ################################################################################
# #HEATMAP+ CLUSTERING
# ################################################################################
# ################################################################################
# ################################################################################
#   observe({
#     set.seed(123)
#     x=sample(1:100,(input$cols*input$rows),replace=TRUE)
#     heatm=matrix(x,nrow=input$cols,ncol=input$rows)
#     variables$oldheat=heatm
#     variables$heat=heatm
#   })

#   #create heatmap (or load/calculate from data) if col changes
#   observe({
#     #heatmap
#     # set.seed(123)
#     # x=sample(1:100,(input$cols*input$rows),replace=TRUE)
#     # heat=matrix(x,nrow=input$cols,ncol=input$rows)
#     observeEvent(input$confirmDim,{
#       heat=variables$oldheat

#       #if a button is pressed, use input$rows and cols, but ISOLATED, so that is not instantly changed
#       textrows=1:input$rows
#       textcols=1:input$cols
#       rowdend=NULL
#       coldend=NULL


#       #if checkbox cluster row is checked, cluster for rows, and modif. vars
#       if (input$clusterrow & input$rows>=2 & input$cols>=2){
#           #hcl <- hclust(dist(t(heat),method=distmethod),method=clustmethod)
#           hcl <- variables$oldheat %>% t %>% dist(method=input$distmethod) %>% hclust(method = input$clustmethod) 
#            #now update with new order:
#           ord= hcl$order
#           heat=variables$oldheat[,ord] 
#           textrows=textrows[ord]
#           rowdend=hcl        
#       }

#       #if checkbox cluster row is checked, cluster for cols,and modif. vars
#       if (input$clustercols & input$rows>=2 & input$cols>=2){
#           #hcl <- hclust(dist(heat,method=distmethod),method=clustmethod )
#           hcl <- variables$oldheat %>% dist(method=input$distmethod) %>% hclust(method = input$clustmethod) 
#           #now update with new order:
#           ord= hcl$order
#           heat=variables$oldheat[ord,]
#           textcols=textcols[ord]
#           coldend=hcl        
#       }

#       #update reactive values (heat matrix, text cols, dendrograms (NULL or not))
#       variables$heat=heat
#       variables$textcols=textcols
#       variables$coldend=coldend
#       variables$rowdend=rowdend
#     })



#   })




# ################################################################################
# ################################################################################
# ################################################################################
# #PLOTS/TEXTS
# ################################################################################
# ################################################################################
# ################################################################################



#   #modify slider input of number of clusters, if heatmap dim are modified:
#   observe({
#      updateSliderInput(session, "BEDnumberclusterscols", value = round(input$cols /3+1),
#           min = 1, max = input$cols )
#      updateSliderInput(session, "BEDnumberclusterrow", value = round(input$rows/3+1),
#           min = 1, max = input$rows )
#   })



#   #output row cluster
#   output$rowcluster <- renderPlot({
#       par(mar = rep(0, 4),oma=rep(0, 4),xpd=NA)
#       #if there is clustering, plot dendr., otherwise empty plot
#       if (!is.null(variables$rowdend)){
#         if (ncol(variables$heat)<300){
#           #calculate which elements belong to which cluster, given by the desired number of clusters in input
#           clusterarray = cutree(variables$rowdend, input$BEDnumberclusterrow)
#           clusternum=length(unique(clusterarray))
#           #plot(as.dendrogram(variables$rowdend),horiz=TRUE,xlab="", sub="",axes=FALSE,ann=FALSE,frame.plot=FALSE,yaxs="i",xaxs="i")         
#           variables$rowdend %>% as.dendrogram %>% set("branches_k_color",k=clusternum) %>% plot  (horiz=TRUE,xlab="", sub="",axes=FALSE,ann=FALSE,frame.plot=FALSE,yaxs="i",xaxs="i")    
#         }else{
#           plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#         }
#       }else{
#         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#       }   
#   })

#   #output col cluster
#   output$colcluster <- renderPlot({
#       par(mar = rep(0, 4),oma=rep(0, 4),xpd=NA)
#       #if there is clustering, plot dendr., otherwise empty plot
#       if (!is.null(variables$coldend)){
#         if (nrow(variables$heat)<300){
#           #calculate which elements belong to which cluster, given by the desired number of clusters in input
#           clusterarray = cutree(variables$coldend, input$BEDnumberclusterscols)
#           clusternum=length(unique(clusterarray))
#           #plot(as.dendrogram(variables$coldend),horiz=FALSE,xlab="", sub="",axes=FALSE,ann=FALSE,frame.plot=FALSE,yaxs="i",xaxs="i")            
#           variables$coldend %>% as.dendrogram %>% set("branches_k_color",k=clusternum) %>% plot (horiz=FALSE,xlab="", sub="",axes=FALSE,ann=FALSE,frame.plot=FALSE,yaxs="i",xaxs="i")          
#         }else{
#           plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')         
#         }
#       }else{  
#         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#       }   
#   })


#   #output noplot
#   output$noplot <- renderPlot({
#     emptyplot()
#   })
#   #output noplot2
#   output$noplot2 <- renderPlot({
#     emptyplot()
#   })
#   #output noplot3
#   output$noplot3 <- renderPlot({
#     emptyplot()
#   })
#   #output noplot4
#   output$noplot4 <- renderPlot({
#     emptyplot()
#   })

#   #output heatmap image
#   output$heatmap <- renderPlot({
#       par(mar = rep(0, 4))
#       image(variables$heat)
#   })

#   #output text plot
#   output$testo <- renderPlot({
#       par(mar = rep(0, 4))
#       plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i")
#       #devo partire dalla coordinata x della metà della prima cella alla meta' dell'ultima
#       halfcellwidth=(1/input$cols)/2
#       #set cex according to number of columns
#       cextext=40/input$cols
#       if (cextext>3){
#         cextext=3
#       }
#       text(x=seq(0+halfcellwidth,1-halfcellwidth,length.out=input$cols),labels=variables$textcols,y=rep(0.5,input$cols),srt=90,cex=cextext  )
        
#   })







# ###############################################################################
# ###############################################################################
# ###############################################################################
# #COORDINATES
# ###############################################################################
# ###############################################################################
# ###############################################################################


# #calc coordinates of click
#   output$coordheat<-renderText({
#    paste0("x=", input$heatmap_click$x, "\ny=", input$heatmap_click$y)
#    ny=round(input$rows*input$heatmap_click$y)
#    nx=round(input$cols*input$heatmap_click$x)
#    paste0("cella x=", nx, "\ncella y=", ny)

#   })

# #calc coordinates of click in dendrogram row
#   output$coorddend<-renderText({
#    paste0("x=", input$rowdendrogram_click$x, "\ny=", input$rowdendrogram_click$y)
#   # ny=round(input$rows*input$heatmap_click$y)
#   # nx=round(input$cols*input$heatmap_click$x)
#   # paste0("cella dend x=", nx, "\ncella dend y=", ny)

#   })


# #calc coordinates & brush
# output$coordbrushheat<-renderText({
#    nymin=round(input$rows*input$heatmap_brush$ymin)
#    nxmin=round(input$cols*input$heatmap_brush$xmin)
#    nymax=round(input$rows*input$heatmap_brush$ymax)
#    nxmax=round(input$cols*input$heatmap_brush$xmax)
#     #adjust if <0 and > max col|rows
#    if (length(nymin)>0){
#        if (nymin<0){
#         nymin=0
#        }
#    }
#    if (length(nxmin)>0){
#        if (nxmin<0){
#          nxmin=0
#        }
#    }
#    if (length(nymax)>0){
#        if (nymax>input$rows){
#          nymax=input$rows
#        }
#    }
#    if (length(nxmax)>0){
#       if (nxmax>input$cols){
#         nxmax=input$cols
#       }
#    }
#    paste0("cella x1=", nxmin, "\ncella y1=", nymin,"\ncella x2=", nxmax, "\ncella y2=",nymax)
#   })


# })
