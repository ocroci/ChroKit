######################################################################
######################################################################
######################################################################
######################################################################
#SINGLE evaluation update ui (std plots)
######################################################################
######################################################################
######################################################################
######################################################################

#update ROI possible selection
observe({
  
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  historylist=as.list(nomi)
  lens=unlist(lapply(ROIvariables$listROI,getLength))
  lens2=paste("(",lens,")",sep="")
  if(length(nomi)>0){
    names(historylist)=paste(nomi,lens2)
  }else{
    names(historylist)=paste(nomi,lens)
  }
  getwdth=lapply(ROIvariables$listROI,getWidth) 
  getwdth=unlist(lapply(getwdth, function(k) {table(!duplicated(k))["TRUE"]==1}))
  #correction on NAs. Maybe fix th error :"NAs are not allowed in subscripted assignments"
  if(!is.null(getwdth)){
    getwdth[is.na(getwdth)]=FALSE
  }
  names(historylist)[getwdth]=paste(names(historylist)[getwdth],"(fixed size)")
  updateSelectInput(session,inputId="ROIchooseSingleEval",label=NULL,
                                    choices = historylist)   
  updateSelectInput(session,inputId="ROI1chooseCmp",label=NULL,
                                    choices = historylist)   
  updateCheckboxGroupInput(session,inputId="ROImaster",label=NULL,
                                    choices = historylist) 
  updateCheckboxGroupInput(session,inputId="ROIsForAnalogHeat",label=NULL,
                                    choices = historylist) 
  updateCheckboxGroupInput(session,inputId="ROIsForProfilesAndBox",label=NULL,
                                    choices = historylist)       
})

#update possible BAM selection for ROI
observe({
  
    input$ROIchooseSingleEval
    if (length(ROIvariables$listROI)>0){
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      pos=match(input$ROIchooseSingleEval,nomi)
      roi=ROIvariables$listROI[[pos]]
      if (!is.null(roi)){
        getbam=names(Enrichlist$rawcoverage[[pos]])
        #getbam=names(getBAMlist(roi))
        if (is.null(getbam)){
          updateSelectInput(session,"BAMchooseSingleEval",NULL,choices=character(0))
        }else{
          updateSelectInput(session,"BAMchooseSingleEval",NULL,choices=getbam)
        }
        
      }
    }

})




######################################################################
######################################################################
######################################################################
######################################################################
# comparison evaluation update ui (std cmp), with enrichments if available
######################################################################
######################################################################
######################################################################
######################################################################


#update ROI2 possible selection. Try to define according to ROI1...
observe({
  
      ROIvariables$listROI
      if (length(ROIvariables$listROI)>=2){
        nomi=unlist(lapply(ROIvariables$listROI,getName))

        pos=match(input$ROI1chooseCmp,nomi)
        if (!is.na(pos)){
          list2=ROIvariables$listROI[-pos]
          nomi2=nomi[- pos]
        }else{
          list2=ROIvariables$listROI
          nomi2=nomi
        }
        
        historylist=as.list(nomi2)
        lens=unlist(lapply(list2,getLength))
        lens2=paste("(",lens,")",sep="")
        if(length(nomi2)>0){
          names(historylist)=paste(nomi2,lens2)
        }else{
          names(historylist)=paste(nomi2,lens)
        }
        getwdth=lapply(list2,getWidth) 
        getwdth=unlist(lapply(getwdth, function(k) {table(!duplicated(k))["TRUE"]==1}))
        #correction on NAs. Maybe fix th error :"NAs are not allowed in subscripted assignments"
        if(!is.null(getwdth)){
          getwdth[is.na(getwdth)]=FALSE
        }
        names(historylist)[getwdth]=paste(names(historylist)[getwdth],"(fixed size)")
        updateSelectInput(session,inputId="ROI2chooseCmp",label=NULL,
                                    choices = historylist)    
      }else{
        nomi=NULL
        historylist=as.list(nomi)
        lens=NULL
        names(historylist)=paste(nomi,lens)
        updateSelectInput(session,inputId="ROI2chooseCmp",label=NULL,
                                    choices = historylist)
      }

})

#update BAM files available (take info form input$ROI1chooseCmp and input$ROI2chooseCmp)

observe({
  #input$ROI1chooseCmp

  if (length(ROIvariables$listROI)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos1=match(input$ROI1chooseCmp,nomi)
    roi1=ROIvariables$listROI[[pos1]]
    pos2=match(input$ROI2chooseCmp,nomi)
    roi2=ROIvariables$listROI[[pos2]]

    rawvals1=Enrichlist$rawcoverage[[pos1]]
    rawvals2=Enrichlist$rawcoverage[[pos2]]
    if (!is.null(roi1)){
      getbam1=names(rawvals1)
      
      if (is.null(getbam1)){
        output$BAMmenuchooseCmp1<-renderUI({NULL})
        #updateSelectInput(session,inputId="BAM1chooseCmp",label=NULL,choices=character(0))
      }else{
        output$BAMmenuchooseCmp1<-renderUI({ 
          list(HTML("<b>Choose enrichment-1:</b>"),
                selectInput(inputId="BAM1chooseCmp", label=NULL,choices=getbam1) ) })
        #updateSelectInput(session,inputId="BAM1chooseCmp",label=NULL,choices=getbam1)
      }
      
    }else{
      output$BAMmenuchooseCmp1<-renderUI({NULL})
      #updateSelectInput(session,inputId="BAM1chooseCmp",label=NULL,choices=character(0))
    }


    if (!is.null(roi2)){
      getbam2=names(rawvals2)
      
      if (is.null(getbam2)){   
        output$BAMmenuchooseCmp2<-renderUI({NULL})    
        #updateSelectInput(session,inputId="BAM2chooseCmp",label=NULL,choices=character(0))
      }else{
        output$BAMmenuchooseCmp2<-renderUI({ 
          list(HTML("<b>Choose enrichment-2:</b>"),
                selectInput(inputId="BAM2chooseCmp", label=NULL,choices=getbam2) ) })
        #updateSelectInput(session,inputId="BAM2chooseCmp",label=NULL,choices=getbam2)
      }
      
    }else{
      output$BAMmenuchooseCmp2<-renderUI({NULL})
      #updateSelectInput(session,inputId="BAM2chooseCmp",label=NULL,choices=character(0))
    }


  }else{
    output$BAMmenuchooseCmp1<-renderUI({NULL})
    output$BAMmenuchooseCmp2<-renderUI({NULL})
  }

})









######################################################################
######################################################################
######################################################################
######################################################################
# update Digital Heatmap lists and options
######################################################################
######################################################################
######################################################################
######################################################################



#update ROI available to view (even the master ROI itself)
observe({
  input$ROImaster
  
  if (length(ROIvariables$listROI)>0 & length(input$ROImaster)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    historylist=as.list(nomi)
    lens=unlist(lapply(ROIvariables$listROI,getLength))
    lens2=paste("(",lens,")",sep="")
    if(length(nomi)>0){
      names(historylist)=paste(nomi,lens2)
    }else{
      names(historylist)=paste(nomi,lens)
    }
    getwdth=lapply(ROIvariables$listROI,getWidth) 
    getwdth=unlist(lapply(getwdth, function(k) {table(!duplicated(k))["TRUE"]==1}))
    names(historylist)[getwdth]=paste(names(historylist)[getwdth],"(fixed size)")    
    updateCheckboxGroupInput(session,"ROIsForDigitalHeat",NULL,choices=historylist)

  }else{
    updateCheckboxGroupInput(session,"ROIsForDigitalHeat",NULL,choices=character(0))
  }
})



#update list ROIs to reorder for digital heatmap
#reorder ROI
observeEvent(input$ROIsForDigitalHeat,{

  #check if ROI are present, valid
  if(is.null(ROIvariables$listROI) | length(ROIvariables$listROI)<1 | length(input$ROIsForDigitalHeat)<1){
    output$reorderROImenuDigitalHeat<-renderUI({NULL})
    return()
  }

  lista=list()
  nomi=input$ROIsForDigitalHeat
  choicelist=as.list(1:length(nomi))
  names(choicelist)=as.character(1:length(nomi))
  for (i in 1:length(nomi)){
    lista[[i]]=fluidRow(column(4,
                            selectInput(inputId = paste("reorderROIdigitalHeat",i,sep=""), label = NULL, 
                          choices = choicelist,selected=i)),
            column(8,HTML(nomi[i])))
  }
     
  output$reorderROImenuDigitalHeat<-renderUI({   
    list(
      HTML("<b>ROIs ordering</b>"),
      wellPanel(id = "logPanelBAM",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
        fluidPage(
          lista
        )
      )
    ) 

  })

},ignoreNULL=FALSE)




#update ROI to cluster digital heatmap
observe({
  input$ROIsForDigitalHeat
  
  if(length(input$ROIsForDigitalHeat)>0){
    wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 100px; max-width: 300px; background-color: #ffffff;",
      updateCheckboxGroupInput(session,"ROIforClusteringDigitalHeat",NULL,choices=input$ROIsForDigitalHeat)
    )
  }else{
    updateCheckboxGroupInput(session,"ROIforClusteringDigitalHeat",NULL,choices=character(0))
  }

})


#update random sample to choose for digital heatmap:
observe({
  input$ROImaster
  
  if(length(ROIvariables$listROI)>0){
    if (length(input$ROImaster)>0){
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      pos=match(input$ROImaster,nomi)
      #all rois selected
      roi=ROIvariables$listROI[pos]
      bigrange=list()
      for(i in 1:length(roi)){ 
        #use granges() to remove eventually gene ID if TSS or transcripts 
        if(!is.null(roi[[i]]))  {
          bigrange[[i]]=granges(getRange(roi[[i]]) )
        }
        
      }
      finalrange=suppressWarnings(Reduce(c,bigrange))
      lbr=length(finalrange)
      if (lbr>0){
        minim=1
      }else{
        minim=0
      }
      if (lbr>2000){
        valshown=2000
      }else{
        valshown=lbr
      }
      updateNumericInput(session,inputId = 'sampleRandomDigitalHeat',label=NULL,min = minim, max = lbr, step = 1000,value=valshown)
    }else{
      updateNumericInput(session,inputId = 'sampleRandomDigitalHeat',label=NULL,min = 0, max = 0, step = 1000,value=0)
              
    }    
  }else{
    updateNumericInput(session,inputId = 'sampleRandomDigitalHeat',label=NULL,min = 0, max = 0, step = 1000,value=0)
  }

})


##observer for type of clustering and clustering parameters
#for the type of clustering, react to "ROIs for cluster" (input$ROIforClusteringDigitalHeat)
observe({
  
  input$ROIforClusteringDigitalHeat
  #if something inside, show the choise between hclust and kmeans
  #otherwise, hide (change in no selection or not selected from beginning)
  if(length(input$ROIforClusteringDigitalHeat)>0){
    output$clustertypeDigitalHeat<-renderUI({
      radioButtons("clusterTypeDigitalHeat",label="Clustering type",choiceNames=c("K-means","hierarchical"),choiceValues=c("kmean","hierarchical"))
    }) 
    output$clusternumbershowDigitalHeat<-renderUI({
      #must implement a check in serverGENOMICS to correct if this number <=0 or > length roi shown
      numericInput(inputId = 'clustnumDigitalHeat',label="Cluster number",min = 1, max = 433, step = 1,value=4)
    }) 
  }else{
    output$clustertypeDigitalHeat<-renderUI({NULL})
    output$clusternumbershowDigitalHeat<-renderUI({NULL})
  }
})

#hidden, advanced options for clustering: reactive to input$ROIforClusteringDigitalHeat (if any)
#and input$clusterTypeDigitalHeat ("hierarchical" or K-means)


observe({
  
  input$ROIforClusteringDigitalHeat
  input$clusterTypeDigitalHeat
  if(length(input$ROIforClusteringDigitalHeat)>0 & !is.null(input$clusterTypeDigitalHeat)){
    #if ==0, do not execute, because radiobutton dosn't exist!
    if(input$clusterTypeDigitalHeat=="hierarchical"){
      output$clusterHDistMethodDigitalHeat<-renderUI({
        selectInput("distmethodDigitalHeat",label="Distance method:",c("Euclidean"="euclidean",
                                                                  "Manhattan"="manhattan",
                                                                  "Canberra"="canberra",
                                                                  "Minkowski"="minkowski"))
      })
      output$clusterHClustMethodDigitalHeat<-renderUI({
        selectInput("clustmethodDigitalHeat",label="Clustering method:",c("Average"="average",
                                                                   "Complete"="complete",
                                                                   "Median"="median",
                                                                   "Centroid"="centroid"))
      })
      output$clusterKstartsDigitalHeat<-renderUI({NULL})
      output$clusterKiterationsDigitalHeat<-renderUI({NULL})
    }else{
      output$clusterKstartsDigitalHeat<-renderUI({
        numericInput(inputId = 'clustrandomstartsDigitalHeat',label="Num. starting points",min = 1, max = 40, step = 1,value=5)
      })
      output$clusterKiterationsDigitalHeat<-renderUI({
        numericInput(inputId = 'clustnumiterationsDigitalHeat',label="Num. iterations",min = 1, max = 40, step = 1,value=10)  
      })
      output$clusterHDistMethodDigitalHeat<-renderUI({NULL})
      output$clusterHClustMethodDigitalHeat<-renderUI({NULL})
    }
  
  }else{
    output$clusterKstartsDigitalHeat<-renderUI({NULL})
    output$clusterKiterationsDigitalHeat<-renderUI({NULL})
    output$clusterHDistMethodDigitalHeat<-renderUI({NULL})
    output$clusterHClustMethodDigitalHeat<-renderUI({NULL})
  }
})



#observer for changing the color scheme:


# toListenDigital_colors <- reactive({
#     list(input$optioncolorsforDigitalHeat,input$confirmUpdateDigitalHeat1,input$confirmUpdateDigitalHeat2)
# })

observeEvent(input$optioncolorsforDigitalHeat,{
  optionDigital=input$optioncolorsforDigitalHeat

  if(!is.null(optionDigital)){
    #if chosen global, show select color input
    if(input$optioncolorsforDigitalHeat=="global"){
      output$showcolorsDigitalheat<-renderUI({
        selectInput("colorCustomDigitalHeat_global",label="Choose global color:",c("white/red"="white_red4",
                                                          "white/blue"="white_blue",
                                                          "white/green"= "white_green"))
      })
      
    }else if(optionDigital=="custom"){
      #if BAM files in heatvariables is >0, loop and show different color blocks,
      #otherwise NULL. Variables are color1, color2, ... colorn

      bams=toplot$digital$ROIsForDigitalHeat
      lista=list()

      if(length(bams)>0){
        output$showcolorsDigitalheat<-renderUI({
          for (i in 1:length(bams)){
            lista[[i]]=fluidRow(
                                column(12,
                                      colorSelectorInput(inputId=paste("colorCustomDigitalHeat",i,sep=""),label=bams[i],choices=ColsArray,
                                                    selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                )
                        )

          }
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
            lista
          )
        }) 
              
      }else{
        #bams are NULL. show nothing
        output$showcolorsDigitalheat<-renderUI({NULL})
      }
    }else{
      output$showcolorsDigitalheat<-renderUI({NULL})
    }    
  }else{
    output$showcolorsDigitalheat<-renderUI({NULL})
  }

})





#observer for chosen color in Digital heat.
observe({
  optionMenu=isolate(input$optioncolorsforDigitalHeat)
  input$colorCustomDigitalHeat1
  input$colorCustomDigitalHeat_global

  if(!is.null(optionMenu)){
    if(optionMenu=="global"){
      if(!is.null(input$colorCustomDigitalHeat_global)){
        toplot$digital$colorpalettes=input$colorCustomDigitalHeat_global
      }else{
        toplot$digital$colorpalettes="white_red4"
      }

    }else if(optionMenu=="custom"){

      bams=isolate(toplot$digital$ROIsForDigitalHeat)
      tosearch=paste("colorCustomDigitalHeat",1:length(bams),sep="")
      # #find all possible colors selected.
      # inputs=names(isolate(input))
      # stringval=grep("colorCustomDigitalHeat\\d+$",inputs,value=TRUE)


      if(length(tosearch)>0){
        listinputcols=list()

        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }


        listinputcols=unlist(listinputcols)
        if(length(listinputcols)==length(tosearch)){
          toplot$digital$colorpalettes=listinputcols
        }else{
          toplot$digital$colorpalettes="white_red4"
        }
      }else{
        toplot$digital$colorpalettes="white_red4"
      }
    }else{
      toplot$digital$colorpalettes="white_red4"
    }
  }else{
    toplot$digital$colorpalettes="white_red4"
  }
})









#observer for the number of bins? (check if nbins > width of some range)
observe({
  
  input$ROIsForDigitalHeat
  if (length(ROIvariables$listROI)>0 & length(input$ROIsForDigitalHeat)>0){

    #adapt numbr of bins. bins must be <= length of the smallest of range
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForDigitalHeat,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    if(all(sapply(roi,class)=="RegionOfInterest")){
      maxtoshow=suppressWarnings(min(unlist(lapply(roi,checkMaxBins))))
      if(maxtoshow<5){
        valuetoshow=maxtoshow
      }else{
        valuetoshow=5
      }
      updateNumericInput(session,"binsDigitalHeat",label=NULL,min = 1, max = maxtoshow,value=valuetoshow,step = 1)      
    }else{
      updateNumericInput(session,"binsDigitalHeat",label=NULL,min = 1, max = 500,value=5,step = 1)
    }
  }else{
    updateNumericInput(session,"binsDigitalHeat",label=NULL,min = 1, max = 500,value=5,step = 1)
  }
})




#react to clicks in the clustering image on the left of the heatmap
#(if n ROI==1 and clustering=TRUE and ROI for clustering ==TRUE)
#create "nw ROI" button as in the analog heatmap
observeEvent(input$rowdendrogram_click_Digital,{
  #check coordinates. If not in the last (and only the last) column (use x coord), only if
  #clustering is activated, do the stuff, otherwise do nothing

  x=round(as.numeric(input$rowdendrogram_click_Digital$x))
  y=as.integer(input$rowdendrogram_click_Digital$y)
  y=y+1 #because is 0-based

  totROIsamples=length(toplot$digital$ROIsForDigitalHeat)
  totrow=toplot$digital$sampleRandomDigitalHeat
  #essential checks:
  if(length(y)>0 & length(totROIsamples)>0 & !is.null(toplot$digital$matrixes_processed)){
    #check if lateral image of clustering colors has been drown
    if(!is.null(toplot$digital$clusternumbermatrix)){
      #from here, heatmap exists, clustering has been made, so we can proceed:
      set.seed(123)
      #numbers go from top of the heatmap (1,2,3...) to the bottom of th heatmap (...,4,5,6)
      #find which cluster was selected
      value=toplot$digital$clusternumbermatrix[,y]  
      clusterselected=which(rev(toplot$digital$clusternumbermatrix==value))
      ybottom=min(clusterselected)
      ytop=max(clusterselected)
      toplot$gadgetdigital$ytop=ytop
      toplot$gadgetdigital$ybottom=ybottom

      ##create the text input and button
      output$newROIfromDigitalHeat_out<-renderUI({
        list(textInput("newROIfromDigitalHeat",label="New ROI from selection",value="",placeholder = "type new ROI name here"),
          actionButton("confirmImportROIfromDigitalHeat", "Import ROI"))
      })

      ##create the text to say the cluster number and the number of intervals selected
      output$textselectedelementsDigitalHeat<-renderText({
        paste("Cluster selected: ",value," (",(ytop-ybottom)+1," regions)",sep="")
      })      

      ## use this block if you want to update the plots in the digital or do something else
      ## ROI is =1 by definition (otherwise toplot$digital$clusternumbermatrix would have been NULL)
      # mat=toplot$digital$matrixes_processed
      # masterROIname=toplot$digital$ROImaster
      
      # rownames(mat)=rep(masterROIname,nrow(mat))
      # matjoin_slice=mat[clusterselected,,drop=FALSE]
      # ROInametoput=rep(toplot$digital$ROIsForDigitalHeat,each=toplot$digital$binsDigitalHeat)
      # colnames(matjoin_slice)=ROInametoput
      

      #create new sliced ROI and put in temporary variable toplot$digital$extractedROI
      #otherwise NULL and no button for NEW ROI
      #we ned the variable containing the range/fix/BAM/... subselected and then r-ordered
      #based on clustering
      #use:
      # toplot$digital$range_tokeep
      # toplot$digital$bams_tokeep
      # toplot$digital$fix_tokeep
      #name, source, kind are derived from the old ones




    }else{
      #clustering not valid in some way, remove new ROI button and do nothing
      output$newROIfromDigitalHeat_out<-renderUI({NULL})
      output$textselectedelementsDigitalHeat<-renderText({NULL})
    }
  }else{
    #remove new ROI button and do nothing
    output$newROIfromDigitalHeat_out<-renderUI({NULL})
    output$textselectedelementsDigitalHeat<-renderText({NULL})
  }

})









######################################################################
######################################################################
######################################################################
######################################################################
# update lists in Analogic Heatmap. ROIs available, intersection of BAM available in all the ROI
#selected and list for choose for which bam to cluster, according to which
#BAM files were selected
######################################################################
######################################################################
######################################################################
######################################################################





#update list of intersection if BAMs available for ROIs selected
observe({
  
  input$ROIsForAnalogHeat

  if (length(ROIvariables$listROI)>0 & length(input$ROIsForAnalogHeat)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForAnalogHeat,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    rawvals=Enrichlist$rawcoverage[pos]
    if (!is.null(roi)){
      #for loop, intersection of bam files
      getbam_current=list()
      for (i in 1:length(roi)){
        if(!is.null(roi[[i]])){

          nms=names(rawvals[[i]])
          if (length(nms)>0){
            getbam_current[[i]]=nms
          }else{
            getbam_current[[i]]=NA
          }           
        }

      }
      finalBAMs=Reduce(intersect,getbam_current)
      if (length(finalBAMs)==1){
        if(is.na(finalBAMs)){
          finalBAMs=character(0)
        }        
      }
      #if at least one BAM in common found, put the chioces
      if (length(finalBAMs)>0){
        updateCheckboxGroupInput(session,"BAMsForAnalogHeat",NULL,choices=finalBAMs)
      }else{
        updateCheckboxGroupInput(session,"BAMsForAnalogHeat",NULL,choices=character(0))
      }
 
    }
  }else{
    updateCheckboxGroupInput(session,"BAMsForAnalogHeat",NULL,choices=character(0))
  }
})



#here, choose ordering of BAM files that will appear in the menu.
#should depend on BAMs selected for analogic heatmap


observeEvent(input$BAMsForAnalogHeat,{ 

  #checks: there must be some ROI
  if(is.null(ROIvariables$listROI) | length(ROIvariables$listROI)<1 | length(input$ROIsForAnalogHeat)<1){
    output$reorderBAMmenuAnalogHeat<-renderUI({NULL})
    return()
  }
  #checks: some BAMs must be selected (>0)
  if(length(input$BAMsForAnalogHeat)<1 ){
    output$reorderBAMmenuAnalogHeat<-renderUI({NULL})
    return()    
  }
    
  getbam=input$BAMsForAnalogHeat
  choicelist=as.list(1:length(getbam))
  names(choicelist)=as.character(1:length(getbam))
  lista=list()
  for (i in 1:length(getbam)){
    lista[[i]]=fluidRow(column(4,          
                                  selectInput(inputId = paste("reorderBAManalogHeat",i,sep=""), label = NULL, 
                                              choices = choicelist,selected=i)                                
                              ),
                        column(8,
                                HTML(getbam[i])
                              )
                        )
  }  

  output$reorderBAMmenuAnalogHeat<-renderUI({   
    list(
      HTML("<b>Enrichment order</b>"),
      wellPanel(id = "logPanelBAM",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
        fluidPage(
          lista
        )
      )
    ) 

  })   
}, ignoreNULL = FALSE)






#update ranking/clustering UI
observe({
  
  input$chooseOrderingAnalogHeat
  input$BAMsForAnalogHeat
  input$ROIsForAnalogHeat

  if(length(input$ROIsForAnalogHeat)==1){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForAnalogHeat,nomi)
    #all rois selected
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      numranges=length(getRange(roi))
    }
  }else{
    numranges=1
  }

  #if clustering is selected 
  if(input$chooseOrderingAnalogHeat=="clustering"){
    #observer for number of clusters. uiOutput only if cluster is chosen (not ranking), roi is >0
    #and only ONE ROI is selected for the heatmap
    output$orderingAnalogHeat<-renderUI({
      if(length(input$BAMsForAnalogHeat)>0 ){
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 100px; max-width: 300px; background-color: #ffffff;",
          checkboxGroupInput("BAMsForClusteringAnalogHeat",NULL,choices=input$BAMsForAnalogHeat)
        )
      }else{
        checkboxGroupInput("BAMsForClusteringAnalogHeat",NULL,choices=character(0)) 
      }

    })

    #block whther to show the choice for the cluster number
    if(length(input$ROIsForAnalogHeat)>0){
      if(length(input$BAMsForAnalogHeat)>0 ){
        output$clusternumbershowAnalogHeat<-renderUI({
          #numericInput("clustnumAnalogHeat",NULL)
          numericInput(inputId = 'clustnumAnalogHeat',label="Cluster number",min = 1, max = numranges, step = 1,value=4)
        }) 
        output$clustertypeAnalogHeat<-renderUI({
          radioButtons("clusterTypeAnalogHeat",label="Clustering type",choiceNames=c("K-means","hierarchical"),choiceValues=c("kmean","hierarchical"))
        })
      }else{
        output$clusternumbershowAnalogHeat<-renderUI({
          checkboxGroupInput("clustnumAnalogHeat",NULL,choices=character(0))
          #numericInput(inputId = 'clustnumAnalogHeat',label="Cluster number",min = 1, max = numranges, step = 1,value=4)
        })
        output$clustertypeAnalogHeat<-renderUI({
          NULL
        })
      }

    }else{
      output$clusternumbershowAnalogHeat<-renderUI({
        checkboxGroupInput("clustnumAnalogHeat",NULL,choices=character(0))
        #numericInput(inputId = 'clustnumAnalogHeat',label="Cluster number",min = 1, max = numranges, step = 1,value=4)
      })

      output$clustertypeAnalogHeat<-renderUI({
        NULL
      })
    }



    
  #if cluster is not selected (ranking selected)
  }else{
    output$orderingAnalogHeat<-renderUI({
      if(length(input$BAMsForAnalogHeat)>0){
          selectInput("BAMsForRankingAnalogHeat",label=NULL,choices=input$BAMsForAnalogHeat)
      }else{
          checkboxGroupInput("BAMsForClusteringAnalogHeat",NULL,choices=character(0))
      }
    })

    output$clusternumbershowAnalogHeat<-renderUI({
      checkboxGroupInput("clustnumAnalogHeat",NULL,choices=character(0))
      #numericInput(inputId = 'clustnumAnalogHeat',label="Cluster number",min = 1, max = numranges, step = 1,value=4)
    })
    output$clustertypeAnalogHeat<-renderUI({
      NULL
    })

  }
})



#reactive to the kind of clustering (hclust or kmean) and the selected bam for clustering
#only if BAM for clustering is not empty
#in this case, modify clusterin options accordingly
observe({
  
  input$chooseOrderingAnalogHeat
  input$clusterTypeAnalogHeat
  input$BAMsForClusteringAnalogHeat

  #if "ranking" instead of "clustering", do not show options for clustering!
  if(input$chooseOrderingAnalogHeat=="clustering"){
    if(length(input$BAMsForClusteringAnalogHeat)>0 & !is.null(input$clusterTypeAnalogHeat)){

      #NOW define which parameters on which type of clustering:
      if(input$clusterTypeAnalogHeat=="hierarchical"){
        output$clusterHDistMethodAnalogHeat<-renderUI({
          selectInput("distmethodAnalogHeat",label="Distance method:",c("Euclidean"="euclidean",
                                                                    "Manhattan"="manhattan",
                                                                    "Canberra"="canberra",
                                                                    "Minkowski"="minkowski"))
        })
        output$clusterHClustMethodAnalogHeat<-renderUI({
          selectInput("clustmethodAnalogHeat",label="Clustering method:",c("Average"="average",
                                                                     "Complete"="complete",
                                                                     "Median"="median",
                                                                     "Centroid"="centroid"))
        })
        output$clusterKstartsAnalogHeat<-renderUI({NULL})
        output$clusterKiterationsAnalogHeat<-renderUI({NULL})

      }else{
        output$clusterKstartsAnalogHeat<-renderUI({
          numericInput(inputId = 'clustrandomstartsAnalogHeat',label="Num. starting points",min = 1, max = 40, step = 1,value=5)
        })
        output$clusterKiterationsAnalogHeat<-renderUI({
          numericInput(inputId = 'clustnumiterationsAnalogHeat',label="Num. iterations",min = 1, max = 40, step = 1,value=10)  
        })
        output$clusterHDistMethodAnalogHeat<-renderUI({NULL})
        output$clusterHClustMethodAnalogHeat<-renderUI({NULL})

      }

    }else{
      output$clusterKstartsAnalogHeat<-renderUI({NULL})
      output$clusterKiterationsAnalogHeat<-renderUI({NULL})
      output$clusterHDistMethodAnalogHeat<-renderUI({NULL})
      output$clusterHClustMethodAnalogHeat<-renderUI({NULL})
    }    
  }else{
    output$clusterKstartsAnalogHeat<-renderUI({NULL})
    output$clusterKiterationsAnalogHeat<-renderUI({NULL})
    output$clusterHDistMethodAnalogHeat<-renderUI({NULL})
    output$clusterHClustMethodAnalogHeat<-renderUI({NULL})    
  }
  
})



#observer for the radiobutton of wich type of coloring scheme is selected for the analog heatmap
#observe({input$optioncolorsforAnalogHeat})
observe({
  input$optioncolorsforAnalogHeat
  if(!is.null(input$optioncolorsforAnalogHeat)){
    #if chosen global, show select color input
    if(input$optioncolorsforAnalogHeat=="global"){
      output$showcolorsheat<-renderUI({
        selectInput("colorCustomAnalogHeat_global",label="Choose global color:",c("white/red"="white_red4",
                                                          "white/blue"="white_blue",
                                                          "rainbow"="rainbow",
                                                          "blue/white/red"="blue_white_red",
                                                          "exponential blue"="white_white_white_blue_blue4",
                                                          "gray scale"="white_grey90_grey80_grey70_grey50_black"))
      })
    }else if(input$optioncolorsforAnalogHeat=="custom"){
      #if BAM files in heatvariables is >0, loop and show different color blocks,
      #otherwise NULL. Variables are color1, color2, ... colorn
      bams=heatvariables$BAMsForAnalogHeat
      lista=list()

      if(length(bams)>0){
        output$showcolorsheat<-renderUI({
          for (i in 1:length(bams)){
            lista[[i]]=fluidRow(
                                column(12,
                                      colorSelectorInput(inputId=paste("colorCustomAnalogHeat",i,sep=""),label=bams[i],choices=ColsArray,
                                                    selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                )
                        )

          }
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
            lista
          )
        })        
      }else{
        #bams are NULL. show nothing
        output$showcolorsheat<-renderUI({NULL})
      }
    }else{
      output$showcolorsheat<-renderUI({NULL})
    }    
  }else{
    output$showcolorsheat<-renderUI({NULL})
  }

})





#observer for chosen color in analog heat.
observe({
  optionMenu=isolate(input$optioncolorsforAnalogHeat)
  input$colorCustomAnalogHeat1
  input$colorCustomAnalogHeat_global
  if(!is.null(optionMenu)){
    if(optionMenu=="global"){
      if(!is.null(input$colorCustomAnalogHeat_global)){
        toplot$analogic$colorpalettes=input$colorCustomAnalogHeat_global
      }else{
        toplot$analogic$colorpalettes="white_red4"
      }

    }else if(optionMenu=="custom"){

      bams=isolate(heatvariables$BAMsForAnalogHeat)
      tosearch=paste("colorCustomAnalogHeat",1:length(bams),sep="")
      # #find all possible colors selected.
      # inputs=names(isolate(input))
      # stringval=grep("colorCustomDigitalHeat\\d+$",inputs,value=TRUE)

      
      if(length(tosearch)>0){
        listinputcols=list()

        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }


        listinputcols=unlist(listinputcols)
        if(length(listinputcols)==length(tosearch)){
          toplot$analogic$colorpalettes=listinputcols
        }else{
          toplot$analogic$colorpalettes="white_red4"
        }
      }else{
        toplot$analogic$colorpalettes="white_red4"
      }
    }else{
      toplot$analogic$colorpalettes="white_red4"
    }
  }else{
    toplot$analogic$colorpalettes="white_red4"
  }
})








#observer for the number of bins? (check if nbins > width of some range)
observe({
  
  input$ROIsForAnalogHeat
  
  if (length(ROIvariables$listROI)>0 & length(input$ROIsForAnalogHeat)>0){
    #adapt numbr of bins. bins must be <= length of the smallest of range
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForAnalogHeat,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    checknulls=any(sapply(roi,is.null))
    if(!checknulls){
      rangeSelected=lapply(roi,getRange)
      rangeSelected=lapply(rangeSelected,granges)
      rangeSelected=suppressWarnings(Reduce(c,rangeSelected))
      maxtoshow=min(width(rangeSelected))
      if(maxtoshow<50){
        valuetoshow=maxtoshow
      }else{
        valuetoshow=50
      }
      updateNumericInput(session,"binsAnalogHeat",label=NULL,min = 1, max = maxtoshow,value=valuetoshow,step = 1)

    }else{
      updateNumericInput(session,"binsAnalogHeat",label=NULL,min = 1, max = 500,value=50,step = 1)
    }

  }else{
    updateNumericInput(session,"binsAnalogHeat",label=NULL,min = 1, max = 500,value=50,step = 1)
  }
})





#update possible random sample, from 1000 to length(sum(ROIs)), with step=1000
#depends only on the ROI selected (input$ROIsForAnalogHeat)
observe({
  
  input$ROIsForAnalogHeat
  if (length(input$ROIsForAnalogHeat)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForAnalogHeat,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    bigrange=list()
    for(i in 1:length(roi)){ 
      #use granges() to remove eventually gene ID if TSS or transcripts   
      if(!is.null(roi[[i]])) {
        bigrange[[i]]=granges(getRange(roi[[i]]) )
      }
      
    }
    
    finalrange=suppressWarnings(Reduce(c,bigrange))
    lbr=length(finalrange)

    if (lbr>0){
      minim=1
    }else{
      minim=0
    }
    if (lbr>2000){
      valshown=2000
    }else{
      valshown=lbr
    }

    updateNumericInput(session,inputId = 'sampleRandomAnalogHeat',label="Random sample of:",min = minim, max = lbr, step = 1000,value=valshown)
  }else{
    updateNumericInput(session,inputId = 'sampleRandomAnalogHeat',label="Random sample of:",min = 0, max = 0, step = 1000,value=0)
            
  }
})






#update x,y heatmap brushed coordinates
#using observeEvent, because with brush is execeuted only once and not twice
#but is like a normal "observe"
observeEvent(input$heatmap_brush,{
  set.seed(123)
  
  xleft=round(as.numeric(input$heatmap_brush$xmin))
  xright=round(as.numeric(input$heatmap_brush$xmax))
  ybottom=as.integer(input$heatmap_brush$ymin)
  ytop=as.integer(input$heatmap_brush$ymax)
  #coord is 0-based, so:
  ybottom=ybottom+1
  ytop=ytop+1



  if (length(xleft)>0){
    valueleft=xleft%/%heatvariables$binsAnalogHeat
    if (xleft%%heatvariables$binsAnalogHeat>0){
      xleft=valueleft+1
    }else{
      xleft=valueleft
    }
  }

  if(length(xright)>0){
    valueright=xright%/%heatvariables$binsAnalogHeat
    if (xright%%heatvariables$binsAnalogHeat>0){
      xright=valueright+1
    }else{
      xright=valueright
    }
  }

  totbamsamples=length(heatvariables$BAMsForAnalogHeat)
  totrow=heatvariables$sampleRandomAnalogHeat
  #if selected outside the max cols/rows or negative (white part of heatmap, if present)
  #adjust if <0 and > max col|rows
  if (length(ybottom)>0){
      if (ybottom<0){
       ybottom=0
      }
  }

  # if (length(ytop)>0){
  #     if (ytop>=totrow){
  #       ytop=totrow-1
  #     }
  # }

  toplot$gadgetanalogic$xright=xright
  toplot$gadgetanalogic$ytop=ytop
  toplot$gadgetanalogic$xleft=xleft
  toplot$gadgetanalogic$ybottom=ybottom
  toplot$gadgetanalogic$totbamsamples=totbamsamples


  #here plot profiles, boxplots, and prepare for t-SNE and PCA as well as gene list ID xtraction (
  #only when present)
  if(length(xright)>0 & length(ytop)>0 & length(xleft)>0 & length(ybottom)>0 & length(totbamsamples)>0 & !is.null(heatvariables$matlist)){
    #extract ROI and BAM from current heatvariables$matlist
    if( (xright-xleft)>=0 & (ytop-ybottom)>=0 & xright<=totbamsamples){
      
      # if (!is.null(heatvariables$completematrixes)){
      #   matlist=heatvariables$completematrixes
      # }else{
      #   matlist=heatvariables$matlist
      # }
      matlist=heatvariables$matlist
      nrowlist=lapply(matlist,nrow)
      ncollist=lapply(matlist,ncol)
      #which ROI and which pieace of ROI?
      ROInames=heatvariables$ROIsForAnalogHeat

      #put the correct label to matlist 
      for(i in 1:length(matlist)){
        rownames(matlist[[i]])=rep(heatvariables$ROIsForAnalogHeat[i],nrow(matlist[[i]]))
      }
      

      #join the matrixes
      matjoin=do.call(rbind,matlist)
      #extract the part of matrix, according to y coordinates chosen in input:
      matjoin_slice=matjoin[(totrow-ytop+1):(totrow-ybottom+1),,drop=FALSE]



      #conserve labels and split matrix without labels, according to labels
      labs=rownames(matjoin_slice)

      #use toplot$analogic$finalrange and toplot$analogic$finalbam to extract, with 
      #coordinates, the ROIs and BAM associated to them 



      #else, continue with binning
      matjoin_slice=as.data.frame(matjoin_slice)
      matjoin_slice_split=split(matjoin_slice,factor(labs,levels=unique(labs)))
      labs2=split(labs,factor(labs,levels=unique(labs)))
      matjoin_slice_split=lapply(matjoin_slice_split,as.matrix)
      matjoin_slice_split=lapply(1:length(matjoin_slice_split),function(i) { mat=matjoin_slice_split[[i]]; rownames(mat) =labs2[[i]];return(mat) })
      names(matjoin_slice_split)=unique(labs)

      #if nROI==1 and clustering, split mat by clustering     
      bycluster=FALSE 
      if(isolate(input$chooseOrderingAnalogHeat)=="clustering" & length(heatvariables$ROIsForAnalogHeat)==1){
        if(!is.na(heatvariables$clustnumAnalogHeat) & heatvariables$clustnumAnalogHeat<=433 & heatvariables$clustnumAnalogHeat>0){
          #split matjoin_slice according to clusters as if they were different ROIs
          bycluster=TRUE
          numclust=as.numeric(toplot$analogic$clusternumbermatrix)
          labs=rev(numclust)[(totrow-ytop+1):(totrow-ybottom+1)]
          matjoin_slice_split=split(matjoin_slice,factor(labs,levels=unique(labs)))
          labs2=split(labs,factor(labs,levels=unique(labs)))
          matjoin_slice_split=lapply(matjoin_slice_split,as.matrix)
          matjoin_slice_split=lapply(1:length(matjoin_slice_split),function(i) { mat=matjoin_slice_split[[i]]; rownames(mat) =labs2[[i]];return(mat) })
          names(matjoin_slice_split)=paste("cluster_",unique(labs),sep="")
        }
      }



      #length for each ROI for each selection:
      lengths_selections=sapply(matjoin_slice_split,nrow)

      #we now have the correct list with correct labels
      #truncate BAM according to xleft; xright
      portionlist=list()

      counter=1
      for(i in 1:length(matjoin_slice_split)){
        currentROI=matjoin_slice_split[[i]]
        #extract BAM selected for the selectd ROI
        toadd_number_selection=lengths_selections[[i]]
        bamselected=heatvariables$BAMsForAnalogHeat
        BAMnametoput=rep(bamselected,each=heatvariables$binsAnalogHeat)
        colnames(currentROI)=BAMnametoput
        #split, this time for columns (BAM files)
        BAMSforthisROI1=lapply(lapply(split( as.data.frame(t(currentROI)) ,factor(BAMnametoput,levels=unique(BAMnametoput))),as.matrix),t)
        BAMSforthisROI=BAMSforthisROI1[xleft:xright]
        #BAMSforthisROI=lapply(BAMSforthisROI1,function(z) { return( z[xleft:xright] ) })
        names(BAMSforthisROI)=paste(names(matjoin_slice_split)[i],"(",toadd_number_selection,"); ",bamselected[xleft:xright],sep="")
        for(k in 1:length(BAMSforthisROI)){
          portionlist[[counter]]=BAMSforthisROI[[k]]
          names(portionlist)[counter]=names(BAMSforthisROI)[k]
          counter=counter+1
        }
      } 

      

      #apply mean on the matrix/median
      portionlist_profile=lapply(1:length(portionlist),function(i){
        mat=portionlist[[i]]
        mat_profile=apply(mat,2,mean)
      })
      portionlist_boxes=lapply(1:length(portionlist),function(i){
        mat=portionlist[[i]]
        mat_boxes=apply(mat,1,sum)
      })

      names(portionlist_profile)=names(portionlist_boxes)=names(portionlist)

      n=length(portionlist_profile)
      #ROIvariables$colorsfordensity <- distinctColorPalette(n)
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      color_distinct <-sample(col_vector, n)


      toplot$gadgetanalogic$Log2BoxAnalogHeat=input$Log2BoxAnalogHeat
      toplot$gadgetanalogic$portionlist_profile=portionlist_profile
      toplot$gadgetanalogic$color_distinct=color_distinct
      toplot$gadgetanalogic$portionlist_boxes=portionlist_boxes
      toplot$gadgetanalogic$lengths_selections=lengths_selections


      roinumber=length(matjoin_slice_split)
      bamnumber=length(portionlist_boxes)/roinumber

      toplot$gadgetanalogic$roinumber=roinumber
      toplot$gadgetanalogic$bamnumber=bamnumber
      toplot$gadgetanalogic$roiname=names(matjoin_slice_split)
      toplot$gadgetanalogic$bamname=bamselected[xleft:xright]

      #output of the number of selected elements within brushed area:
      output$textselectedelementsAnalogHeat<-renderText({
        paste("Intervals selected: ",(ytop-ybottom)+1,sep="")
      })


      ##create the text input and button
      output$newROIfromAnalogHeat_out<-renderUI({
        list(textInput("newROIfromAnalogHeat",label="New ROI from selection",value="",placeholder = "type new ROI name here"),
          actionButton("confirmImportROIfromAnalogHeat", "Import ROI"))
      })

      #PLOTS

      #here, legends are ROIs or clusters, depending on bycluster variable. If TRUE, clsuter 1,2,3,
      #otherwise, ROI1,2,3 and in the legend we need the explanation: ROI1= blablabla
      #plot the profiles of regions selected
      output$profileAnalogHeat<-renderPlot({
        if(input$Log2BoxAnalogHeat){
          portionlist_profile2=lapply(portionlist_profile,log2)
          yl="Log2 read density (rpm/bp)"
        }else{
          portionlist_profile2=portionlist_profile
          yl="Read density (rpm/bp)"
        }
        maxval=max(unlist(lapply(portionlist_profile2,max)))
        minval=min(unlist(lapply(portionlist_profile2,min)))

        par(mar=c(4,4,2,2))
        plot(1, type="n", xlab="",xaxt="n", ylab=yl, xlim=c(1, length(portionlist_profile2[[1]])), ylim=c(minval, maxval))
        axis(1,at=c(1,length(portionlist_profile2[[1]])/2 +0.5,length(portionlist_profile2[[1]])),labels=c("start","center","end"))

        for(i in 1:length(portionlist_profile2)){
          lines(portionlist_profile2[[i]],lwd=2,col=color_distinct[i])
        }
        legend("topright",legend=names(portionlist_profile),col=color_distinct,lty=rep(1,length(color_distinct)),cex=0.7,bg="transparent")   
      })




      #plot boxplots of regions selected. Use portionlist_boxes
      output$boxplotByBAMAnalogHeat<-renderPlot({
        #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
        #invert the order
        newlist=list()
        newcols=c()
        newnames=c()
        #if grouping colors, adjust legend and colors by groups
        if(!input$GroupColorsAnalogHeat){
          for(i in 1:bamnumber){
            pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
            newlist=c(newlist,portionlist_boxes[pos_inverted])
            newcols=c(newcols,color_distinct[pos_inverted])
            newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
          }          
        }else{
          #in this case, take n° ROI colors, repeated for n° BAMs
          #in the legend, put n° ROI colors and tell which ROI with number
          #in xlab, put BAMs
          for(i in 1:bamnumber){
            pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
            newlist=c(newlist,portionlist_boxes[pos_inverted])
          }
          #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
          newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
          newcols=rep(color_distinct[1:toplot$gadgetanalogic$roinumber],toplot$gadgetanalogic$bamnumber)
        }


        factor_add=rep(0:(bamnumber-1),each=roinumber)
        addingfactor=1:length(newlist)+factor_add

        par(mar=c(14,4,1,1))
        if(input$Log2BoxAnalogHeat){
          newlist2=lapply(newlist,log2)
          yl="Log2 read density (rpm/bp)"
        }else{
          newlist2=newlist
          yl="Read density (rpm/bp)"
        }

        suppressWarnings(boxplot(newlist2,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
        ats=c()
        for(i in 1:bamnumber){
          window=addingfactor[(((i-1)*roinumber)+1):(i*roinumber)]
          currentvalue=(window[length(window)]-window[1])/2
          currentvalue=window[1]+currentvalue
          ats=c(ats,currentvalue)
        }

        axis(1,at=ats,label=bamselected[xleft:xright],las=2)
        legend("topright",legend=newnames,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))
     
      })

      output$boxplotByROIAnalogHeat<-renderPlot({
        factor_add=rep(0:(roinumber-1),each=bamnumber)
        addingfactor=1:length(portionlist_boxes)+factor_add
        if(input$Log2BoxAnalogHeat){
          portionlist_boxes2=lapply(portionlist_boxes,log2)
          yl="Log2 read density (rpm/bp)"
        }else{
          portionlist_boxes2=portionlist_boxes
          yl="Read density (rpm/bp)"
        }

        isgrouped=input$GroupColorsAnalogHeat
        #if grouping colors, adjust legend and colors by groups
        if(!isgrouped){
          newcols=color_distinct
          newnames=names(portionlist_boxes)
        }else{
          #in this case, take n° ROI colors, repeated for n° BAMs
          #in the legend, put n° ROI colors and tell which ROI with number
          #in xlab, put BAMs
          #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
          newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
          newcols=rep(color_distinct[1:toplot$gadgetanalogic$bamnumber],toplot$gadgetanalogic$roinumber)
        }

        par(mar=c(14,4,1,1))
        suppressWarnings(boxplot(portionlist_boxes2,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
        ats=c()
        for(i in 1:roinumber){
          window=addingfactor[(((i-1)*bamnumber)+1):(i*bamnumber)]
          currentvalue=(window[length(window)]-window[1])/2
          currentvalue=window[1]+currentvalue
          ats=c(ats,currentvalue)
        }

        if(!isgrouped){
          axis(1,at=ats,label=toplot$gadgetanalogic$roiname,las=2)
          legend("topright",legend=newnames,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))          
        }else{
          axis(1,at=ats,label=newnames,las=2)
          legend("topright",legend=toplot$gadgetanalogic$bamname,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))
        }

      })  





      #buttons for download PDF of profile, boxplot by ROI, boxplot by BAM
      output$saveprofileAnalogHeat=renderUI({downloadButton('saveprofileAnalogHeatbutton', 'Get PDF')})
      output$saveboxplotByROIAnalogHeat=renderUI({downloadButton('saveboxplotByROIAnalogHeatbutton', 'Get PDF')})
      output$saveboxplotByBAMAnalogHeat=renderUI({downloadButton('saveboxplotByBAMAnalogHeatbutton', 'Get PDF')})

      #plot correlation heatmap (and partial correlation?)
      if(roinumber==1 & length(portionlist_boxes)>=2 & (ytop-ybottom)>1){
        if(input$Log2BoxAnalogHeat){
          portionlist_boxescors=lapply(portionlist_boxes,log2)
        }else{
          portionlist_boxescors=portionlist_boxes
        }
        mat=do.call(cbind,portionlist_boxescors)
        
        #problem: when playing with heatmap, Error in colnames<-: length of 'dimnames' [2] not equal to array extent.
        #maybe a simple trycatch would do the job
        tryCatch({
          colnames(mat)=bamselected[xleft:xright]
        },
        warning = function( w ){
          colnames(mat)=1:ncol(mat)
        },
        error = function( err ){
          colnames(mat)=1:ncol(mat)
        }
        )
        
        mat[is.infinite(mat) &mat<0 ]=0
        correlation_total=cor(mat)
        trasp_cor=t(correlation_total)
        brk=c( seq( -1 , 1,0.01))
        my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

        toplot$gadgetanalogic$trasp_cor=trasp_cor
        toplot$gadgetanalogic$my_palette=my_palette
        toplot$gadgetanalogic$brk=brk
        toplot$gadgetanalogic$correlation_total=correlation_total
        

        output$corAnalogHeat<-renderPlot({
          par(mar=c(12,12,1,1),xpd=TRUE)
          image(0:nrow(trasp_cor), 0:ncol(trasp_cor),trasp_cor[,ncol(trasp_cor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
          axis( 2, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= rev(colnames( trasp_cor )), las= 2 )
          axis( 1, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= colnames( trasp_cor ), las= 2 )
          for (x in (nrow(correlation_total)-1+0.5):0.5  )
            for (y in 0.5: ((ncol(correlation_total)-1+0.5)   ))
              text(y,x, round(correlation_total[ncol(correlation_total)-x+0.5,y+0.5],2),col="blue")
        })
        output$savecorAnalogHeat=renderUI({downloadButton('savecorAnalogHeatbutton', 'Get PDF')})

        if (length(portionlist_boxescors)>2 & (ytop-ybottom)>1){
          #if number of BAMs is >2, calculate the partial correlation too
          correlation_partial=pcor(mat)$estimate
          #warning: The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero
          colnames(correlation_partial)=rownames(correlation_partial)=colnames(mat)
          trasp_pcor=t(correlation_partial)
          toplot$gadgetanalogic$trasp_pcor=trasp_pcor
          toplot$gadgetanalogic$correlation_partial=correlation_partial
          output$pcorAnalogHeat<-renderPlot({
            par(mar=c(12,12,1,1),xpd=TRUE)
            image(0:nrow(trasp_pcor), 0:ncol(trasp_pcor),trasp_pcor[,ncol(trasp_pcor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
            axis( 2, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= rev(colnames( trasp_pcor )), las= 2 )
            axis( 1, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= colnames( trasp_pcor ), las= 2 )
            for (x in (nrow(correlation_partial)-1+0.5):0.5  )
              for (y in 0.5: ((ncol(correlation_partial)-1+0.5)   ))
                text(y,x, round(correlation_partial[ncol(correlation_partial)-x+0.5,y+0.5],2),col="blue")
          }) 
          output$savepcorAnalogHeat=renderUI({downloadButton('savepcorAnalogHeatbutton', 'Get PDF')})         
        }else{
          #number of BAM==2 => no pcor
          output$pcorAnalogHeat<-renderPlot({NULL})
          output$savepcorAnalogHeat=renderUI({NULL})
        }

      }else{
        output$corAnalogHeat<-renderPlot({NULL})
        output$pcorAnalogHeat<-renderPlot({NULL})
        output$savecorAnalogHeat=renderUI({NULL})
        output$savepcorAnalogHeat=renderUI({NULL})
      }
    }else{
      output$textselectedelementsAnalogHeat<-renderText({NULL})
      output$newROIfromAnalogHeat_out<-renderUI({NULL})
      #output$confirmImportROIfromAnalogHeat<-renderUI({NULL})
      output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
      output$boxplotByROIAnalogHeat<-renderPlot({NULL})
      output$profileAnalogHeat<-renderPlot({NULL})
      output$saveprofileAnalogHeat=renderUI({NULL})
      output$saveboxplotByROIAnalogHeat=renderUI({NULL})
      output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
      output$corAnalogHeat<-renderPlot({NULL})
      output$pcorAnalogHeat<-renderPlot({NULL})
      output$savecorAnalogHeat=renderUI({NULL})
      output$savepcorAnalogHeat=renderUI({NULL})
    } 
  }else{
    #output of the number of selected elements within brushed area:
    output$textselectedelementsAnalogHeat<-renderText({NULL})
    output$newROIfromAnalogHeat_out<-renderUI({NULL})
    #output$confirmImportROIfromAnalogHeat<-renderUI({NULL})
    output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
    output$boxplotByROIAnalogHeat<-renderPlot({NULL})
    output$profileAnalogHeat<-renderPlot({NULL})
    output$saveprofileAnalogHeat=renderUI({NULL})
    output$saveboxplotByROIAnalogHeat=renderUI({NULL})
    output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
    output$corAnalogHeat<-renderPlot({NULL})
    output$pcorAnalogHeat<-renderPlot({NULL})
    output$savecorAnalogHeat=renderUI({NULL})
    output$savepcorAnalogHeat=renderUI({NULL})
  }
 # paste0("cella x1=", xleft, "\ncella y1=", ybottom,"\ncella x2=", xright, "\ncella y2=",ytop)
})





##respond to heatmap$heatmap_click (click the clustering rectangle)
##int this case:
#-clustering must be present (image clustering done with the colors)
#-clustering shown as a column on the right (maybe optional)
#-the xright should be exactly at the number of BAM (xright==totbamsamples+1)
#-if everything is satisfied, clear brush selection with session$resetBrush("plotBrush") (shiny v. >= 0.14, 2016)
#-retrieve ymin,ymax of that cluster
#-perform plots as above
observeEvent(input$rowdendrogram_click_Analog,{
  #check coordinates. If not in the last (and only the last) column (use x coord), only if
  #clustering is activated, do the stuff, otherwise do nothing

  x=round(as.numeric(input$rowdendrogram_click_Analog$x))
  y=as.integer(input$rowdendrogram_click_Analog$y)
  y=y+1 #because is 0-based

  totbamsamples=length(heatvariables$BAMsForAnalogHeat)
  totrow=heatvariables$sampleRandomAnalogHeat

  #essential checks: everything must exist
  if(length(y)>0 & length(totbamsamples)>0 & !is.null(heatvariables$matlist)){
    #now check if x is exactly in the totbamsamples+1 column (where cluster blocks are drown)
    #and clustering must exist
    #& length(isolate(input$BAMsForClusteringAnalogHeat))>0 & length(isolate(input$ROIsForAnalogHeat)) ==1
    #are FALSE, but becayse every input field is reset. For this, keep all
    if(isolate(input$chooseOrderingAnalogHeat)=="clustering" & length(heatvariables$ROIsForAnalogHeat)==1){
      if(!is.na(heatvariables$clustnumAnalogHeat) & heatvariables$clustnumAnalogHeat<=433 & heatvariables$clustnumAnalogHeat>0){
        #selection is valid. So, clear the brush
        session$resetBrush("heatmap_brush")
        
        set.seed(123)
        #numbers go from top of the heatmap (1,2,3...) to the bottom of th heatmap (...,4,5,6)
        #find which cluster was selected
        value=toplot$analogic$clusternumbermatrix[,y]        
        clusterselected=which(toplot$analogic$clusternumbermatrix==value)

        ybottom=min(clusterselected)
        ytop=max(clusterselected)

        toplot$gadgetanalogic$ytop=ytop
        toplot$gadgetanalogic$ybottom=ybottom

        #take the matrix. ROI is ==1 by definition, so no problems of splitting
        mat=heatvariables$matlist[[1]]
        ROIname=heatvariables$ROIsForAnalogHeat
        rownames(mat)=rep(ROIname,nrow(mat))


        #extract the part of matrix, according to cluster number chosen in input:
        matjoin_slice=mat[nrow(mat):1,]
        matjoin_slice=matjoin_slice[clusterselected,,drop=FALSE]

        
        #no slicing on x axis (cluster takes all bams)
        bamselected=heatvariables$BAMsForAnalogHeat
        BAMnametoput=rep(bamselected,each=heatvariables$binsAnalogHeat)
        colnames(matjoin_slice)=BAMnametoput
        BAMSforthisROI=lapply(lapply(split( as.data.frame(t(matjoin_slice)) ,factor(BAMnametoput,levels=unique(BAMnametoput))),as.matrix),t)
        
        #we have to end up with "portionlist" variable (as in the brush)
        portionlist=BAMSforthisROI
        names(portionlist)=paste(ROIname,"(",nrow(matjoin_slice),"); ",bamselected,sep="")
        

        #apply mean on the matrix/median
        portionlist_profile=lapply(1:length(portionlist),function(i){
          mat=portionlist[[i]]
          mat_profile=apply(mat,2,mean)
        })
        portionlist_boxes=lapply(1:length(portionlist),function(i){
          mat=portionlist[[i]]
          mat_boxes=apply(mat,1,sum)
        })

        names(portionlist_profile)=names(portionlist_boxes)=names(portionlist)

        n=length(portionlist_profile)
        #ROIvariables$colorsfordensity <- distinctColorPalette(n)
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        color_distinct <-sample(col_vector, n)


        toplot$gadgetanalogic$Log2BoxAnalogHeat=input$Log2BoxAnalogHeat
        toplot$gadgetanalogic$portionlist_profile=portionlist_profile
        toplot$gadgetanalogic$color_distinct=color_distinct
        toplot$gadgetanalogic$portionlist_boxes=portionlist_boxes


        roinumber=1#length(matjoin_slice_split)
        bamnumber=length(bamselected)#length(portionlist_boxes)/roinumber

        toplot$gadgetanalogic$roinumber=roinumber
        toplot$gadgetanalogic$bamnumber=bamnumber
        toplot$gadgetanalogic$roiname=ROIname
        toplot$gadgetanalogic$bamname=bamselected

        #output of the number of selected elements within brushed area:
        output$textselectedelementsAnalogHeat<-renderText({
          paste("Cluster selected: ",value," (",(ytop-ybottom)+1," regions)",sep="")
        })


        ##create the text input and button
        output$newROIfromAnalogHeat_out<-renderUI({
          list(textInput("newROIfromAnalogHeat",label="New ROI from selection",value="",placeholder = "type new ROI name here"),
            actionButton("confirmImportROIfromAnalogHeat", "Import ROI"))
        })

        #PLOTS
        #plot the profiles of regions selected
        output$profileAnalogHeat<-renderPlot({
          if(input$Log2BoxAnalogHeat){
            portionlist_profile2=lapply(portionlist_profile,log2)
            yl="Log2 read density (rpm/bp)"
          }else{
            portionlist_profile2=portionlist_profile
            yl="Read density (rpm/bp)"
          }
          maxval=max(unlist(lapply(portionlist_profile2,max)))
          minval=min(unlist(lapply(portionlist_profile2,min)))

          par(mar=c(4,4,2,2))
          plot(1, type="n", xlab="",xaxt="n", ylab=yl, xlim=c(1, length(portionlist_profile2[[1]])), ylim=c(minval, maxval))
          axis(1,at=c(1,length(portionlist_profile2[[1]])/2 +0.5,length(portionlist_profile2[[1]])),labels=c("start","center","end"))

          for(i in 1:length(portionlist_profile2)){
            lines(portionlist_profile2[[i]],lwd=2,col=color_distinct[i])
          }
          legend("topright",legend=names(portionlist_profile),col=color_distinct,lty=rep(1,length(color_distinct)),cex=0.7,bg="transparent")   
        })




        #plot boxplots of regions selected. Use portionlist_boxes
        output$boxplotByBAMAnalogHeat<-renderPlot({
          #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
          #invert the order
          newlist=list()
          newcols=c()
          newnames=c()
          #if grouping colors, adjust legend and colors by groups
          if(!input$GroupColorsAnalogHeat){
            for(i in 1:bamnumber){
              pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
              newlist=c(newlist,portionlist_boxes[pos_inverted])
              newcols=c(newcols,color_distinct[pos_inverted])
              newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
            }          
          }else{
            #in this case, take n° ROI colors, repeated for n° BAMs
            #in the legend, put n° ROI colors and tell which ROI with number
            #in xlab, put BAMs
            for(i in 1:bamnumber){
              pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
              newlist=c(newlist,portionlist_boxes[pos_inverted])
            }
            #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
            newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
            newcols=rep(color_distinct[1:toplot$gadgetanalogic$roinumber],toplot$gadgetanalogic$bamnumber)
          }


          factor_add=rep(0:(bamnumber-1),each=roinumber)
          addingfactor=1:length(newlist)+factor_add

          par(mar=c(14,4,1,1))
          if(input$Log2BoxAnalogHeat){
            newlist2=lapply(newlist,log2)
            yl="Log2 read density (rpm/bp)"
          }else{
            newlist2=newlist
            yl="Read density (rpm/bp)"
          }

          suppressWarnings(boxplot(newlist2,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
          ats=c()
          for(i in 1:bamnumber){
            window=addingfactor[(((i-1)*roinumber)+1):(i*roinumber)]
            currentvalue=(window[length(window)]-window[1])/2
            currentvalue=window[1]+currentvalue
            ats=c(ats,currentvalue)
          }

          axis(1,at=ats,label=bamselected,las=2)
          legend("topright",legend=newnames,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))
       
        })

        output$boxplotByROIAnalogHeat<-renderPlot({
          factor_add=rep(0:(roinumber-1),each=bamnumber)
          addingfactor=1:length(portionlist_boxes)+factor_add
          if(input$Log2BoxAnalogHeat){
            portionlist_boxes2=lapply(portionlist_boxes,log2)
            yl="Log2 read density (rpm/bp)"
          }else{
            portionlist_boxes2=portionlist_boxes
            yl="Read density (rpm/bp)"
          }

          isgrouped=input$GroupColorsAnalogHeat
          #if grouping colors, adjust legend and colors by groups
          if(!isgrouped){
            newcols=color_distinct
            newnames=names(portionlist_boxes)
          }else{
            #in this case, take n° ROI colors, repeated for n° BAMs
            #in the legend, put n° ROI colors and tell which ROI with number
            #in xlab, put BAMs
            #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
            newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
            newcols=rep(color_distinct[1:toplot$gadgetanalogic$bamnumber],toplot$gadgetanalogic$roinumber)
          }

          par(mar=c(14,4,1,1))
          suppressWarnings(boxplot(portionlist_boxes2,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
          ats=c()
          for(i in 1:roinumber){
            window=addingfactor[(((i-1)*bamnumber)+1):(i*bamnumber)]
            currentvalue=(window[length(window)]-window[1])/2
            currentvalue=window[1]+currentvalue
            ats=c(ats,currentvalue)
          }

          if(!isgrouped){
            axis(1,at=ats,label=toplot$gadgetanalogic$roiname,las=2)
            legend("topright",legend=newnames,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))          
          }else{
            axis(1,at=ats,label=newnames,las=2)
            legend("topright",legend=toplot$gadgetanalogic$bamname,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))
          }

        })  

        #buttons for download PDF of profile, boxplot by ROI, boxplot by BAM
        output$saveprofileAnalogHeat=renderUI({downloadButton('saveprofileAnalogHeatbutton', 'Get PDF')})
        output$saveboxplotByROIAnalogHeat=renderUI({downloadButton('saveboxplotByROIAnalogHeatbutton', 'Get PDF')})
        output$saveboxplotByBAMAnalogHeat=renderUI({downloadButton('saveboxplotByBAMAnalogHeatbutton', 'Get PDF')})
    
        #plot correlation heatmap (and partial correlation?)
        if(roinumber==1 & length(portionlist_boxes)>=2& length(value)>1){
          if(input$Log2BoxAnalogHeat){
            portionlist_boxescors=lapply(portionlist_boxes,log2)
          }else{
            portionlist_boxescors=portionlist_boxes
          }
          mat=do.call(cbind,portionlist_boxescors)
          
          #problem: when playing with heatmap, Error in colnames<-: length of 'dimnames' [2] not equal to array extent.
          #maybe a simple trycatch would do the job
          tryCatch({
            colnames(mat)=bamselected
          },
          warning = function( w ){
            colnames(mat)=1:ncol(mat)
          },
          error = function( err ){
            colnames(mat)=1:ncol(mat)
          }
          )
          
          mat[is.infinite(mat) &mat<0 ]=0
          correlation_total=cor(mat)
          trasp_cor=t(correlation_total)
          brk=c( seq( -1 , 1,0.01))
          my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

          toplot$gadgetanalogic$trasp_cor=trasp_cor
          toplot$gadgetanalogic$my_palette=my_palette
          toplot$gadgetanalogic$brk=brk
          toplot$gadgetanalogic$correlation_total=correlation_total
          

          output$corAnalogHeat<-renderPlot({
            par(mar=c(12,12,1,1),xpd=TRUE)
            image(0:nrow(trasp_cor), 0:ncol(trasp_cor),trasp_cor[,ncol(trasp_cor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
            axis( 2, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= rev(colnames( trasp_cor )), las= 2 )
            axis( 1, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= colnames( trasp_cor ), las= 2 )
            for (x in (nrow(correlation_total)-1+0.5):0.5  )
              for (y in 0.5: ((ncol(correlation_total)-1+0.5)   ))
                text(y,x, round(correlation_total[ncol(correlation_total)-x+0.5,y+0.5],2),col="blue")
          })
          output$savecorAnalogHeat=renderUI({downloadButton('savecorAnalogHeatbutton', 'Get PDF')})


          if (length(portionlist_boxescors)>2 & length(value)>1){
            #if number of BAMs is >2, calculate the partial correlation too
            correlation_partial=pcor(mat)$estimate
            #warning: The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero
            colnames(correlation_partial)=rownames(correlation_partial)=colnames(mat)
            trasp_pcor=t(correlation_partial)

            toplot$gadgetanalogic$trasp_pcor=trasp_pcor
            toplot$gadgetanalogic$correlation_partial=correlation_partial

            output$pcorAnalogHeat<-renderPlot({
              par(mar=c(12,12,1,1),xpd=TRUE)
              image(0:nrow(trasp_pcor), 0:ncol(trasp_pcor),trasp_pcor[,ncol(trasp_pcor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
              axis( 2, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= rev(colnames( trasp_pcor )), las= 2 )
              axis( 1, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= colnames( trasp_pcor ), las= 2 )
              for (x in (nrow(correlation_partial)-1+0.5):0.5  )
                for (y in 0.5: ((ncol(correlation_partial)-1+0.5)   ))
                  text(y,x, round(correlation_partial[ncol(correlation_partial)-x+0.5,y+0.5],2),col="blue")
            }) 


            output$savepcorAnalogHeat=renderUI({downloadButton('savepcorAnalogHeatbutton', 'Get PDF')})         
          
          }else{
            #number of BAM==2 => no pcor
            output$pcorAnalogHeat<-renderPlot({NULL})
            output$savepcorAnalogHeat=renderUI({NULL})
          }

        }else{
          output$corAnalogHeat<-renderPlot({NULL})
          output$pcorAnalogHeat<-renderPlot({NULL})
          output$savecorAnalogHeat=renderUI({NULL})
          output$savepcorAnalogHeat=renderUI({NULL})
        }
      }else{
        output$textselectedelementsAnalogHeat<-renderText({NULL})
        output$newROIfromAnalogHeat_out<-renderUI({NULL})
        #output$confirmImportROIfromAnalogHeat<-renderUI({NULL})
        output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
        output$boxplotByROIAnalogHeat<-renderPlot({NULL})
        output$profileAnalogHeat<-renderPlot({NULL})
        output$saveprofileAnalogHeat=renderUI({NULL})
        output$saveboxplotByROIAnalogHeat=renderUI({NULL})
        output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
        output$corAnalogHeat<-renderPlot({NULL})
        output$pcorAnalogHeat<-renderPlot({NULL})
        output$savecorAnalogHeat=renderUI({NULL})
        output$savepcorAnalogHeat=renderUI({NULL})        
      }  
    } 
  }


})




######################################################################
######################################################################
######################################################################
######################################################################
# update lists in Profiles & boxplots tab. Choose ROIs and BAMs in common to these ROIs 
######################################################################
######################################################################
######################################################################
######################################################################




#update list of intersection if BAMs available for ROIs selected
observe({
  
  input$ROIsForProfilesAndBox
  if (length(ROIvariables$listROI)>0 & length(input$ROIsForProfilesAndBox)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForProfilesAndBox,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    rawvals=Enrichlist$rawcoverage[pos]
    if (!is.null(roi)){
      #for loop, intersection of bam files
      getbam_current=list()
      for (i in 1:length(roi)){
        if(!is.null(roi[[i]])){
          nms=names(rawvals[[i]])
          if (length(nms)>0){
            getbam_current[[i]]=nms
          }else{
            getbam_current[[i]]=NA
          }           
        }
      }
      finalBAMs=Reduce(intersect,getbam_current)
      if (length(finalBAMs)==1){
        if(is.na(finalBAMs)){
          finalBAMs=character(0)
        }        
      }
      #if at least one BAM in common found, put the chioces
      if (length(finalBAMs)>0){
        updateCheckboxGroupInput(session,"BAMsForProfilesAndBox",NULL,choices=finalBAMs)
      }else{
        updateCheckboxGroupInput(session,"BAMsForProfilesAndBox",NULL,choices=character(0))
      }
 
    }
  }else{
    updateCheckboxGroupInput(session,"BAMsForProfilesAndBox",NULL,choices=character(0))
  }
})


#observer for the number of bins? (check if nbins > width of some range)
observe({
  
  input$ROIsForProfilesAndBox
  if (length(ROIvariables$listROI)>0 & length(input$ROIsForProfilesAndBox)>0){
    #adapt numbr of bins. bins must be <= length of the smallest of range
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForProfilesAndBox,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    
    counter=0
    for(i in roi){
      if(is.null(i)){
        counter=1
      }
    }
    if(counter==0){
      rangeSelected=lapply(roi,getRange)
      rangeSelected=lapply(rangeSelected,granges)
      rangeSelected=suppressWarnings(Reduce(c,rangeSelected))
      maxtoshow=min(width(rangeSelected))
      if(maxtoshow<50){
        valuetoshow=maxtoshow
      }else{
        valuetoshow=50
      }
      updateNumericInput(session,"binsProfilesAndBox",label=NULL,min = 1, max = maxtoshow,value=valuetoshow,step = 1)
      
    }



  }else{
    updateNumericInput(session,"binsProfilesAndBox",label=NULL,min = 1, max = 500,value=50,step = 1)
  }
})



toListenCor <- reactive({
    list(input$cor_click$y,input$cor_click$x,input$Log2BoxProfilesAndBox)
})


#observer for click the correlation/partial correlation heatmap
observeEvent(toListenCor(),{
  
  set.seed(123)
  if(!is.null(input$cor_click$x) & !is.null(input$cor_click$y) ){
    toplot$profileAndBoxes$x=floor(as.numeric(input$cor_click$x))+1
    toplot$profileAndBoxes$y=floor(as.numeric(input$cor_click$y))+1
        
  }
  portionlist_boxes=corvariables$portionlist_boxes 
  if(!is.null(portionlist_boxes) & length(toplot$profileAndBoxes$x)>0 & length(toplot$profileAndBoxes$y) >0){

    array_x=portionlist_boxes[[toplot$profileAndBoxes$x]]
    array_y=portionlist_boxes[[length(portionlist_boxes)-toplot$profileAndBoxes$y+1]]
    name_x=names(portionlist_boxes)[toplot$profileAndBoxes$x]
    name_y=names(portionlist_boxes)[length(portionlist_boxes)-toplot$profileAndBoxes$y+1]
    #take a random sample (2000) of the total number of points
    if(length(array_x)>2000){
      smpl=sort(sample(length(array_x),2000,replace=FALSE))
      array_x=array_x[smpl]
      array_y=array_y[smpl]
    }
    if(input$Log2BoxProfilesAndBox){
      array_x=log2(array_x)
      array_y=log2(array_y)
      name_x=paste("log2",name_x)
      name_y=paste("log2",name_y)
    }
    if(corvariables$is.density){
      name_x=paste(name_x,"read density (rpm/bp)")
      name_y=paste(name_y,"read density (rpm/bp)")
    }else{
      name_x=paste(name_x,"normalized reads (rpm)")
      name_y=paste(name_y,"normalized reads (rpm)")
    }
    #remove -Inf, because log2(0)=-Inf
    array_x[is.infinite(array_x) &array_x<0 ]=0
    array_y[is.infinite(array_y) &array_y<0 ]=0

    #find correlation values (both pearson and spearman)
    corP=round(cor(array_x,array_y,method="pearson"),2)
    corS=round(cor(array_x,array_y,method="spearman"),2)
    corlegend=c(paste("cor Pearson:",corP),paste("cor Spearman:",corS))

    # toplot$profileAndBoxes$logicscatter="ON"
    toplot$profileAndBoxes$array_x=array_x
    toplot$profileAndBoxes$array_y=array_y
    toplot$profileAndBoxes$name_x=name_x
    toplot$profileAndBoxes$name_y=name_y
    toplot$profileAndBoxes$corlegend=corlegend
    toplot$profileAndBoxes$px=NULL
    toplot$profileAndBoxes$py=NULL

    output$scatterProfilesAndBox<-renderPlot({
      my.cols <- rev(brewer.pal(5, "RdYlBu"))
      z <- kde2d(array_x, array_y, n=50)
      plot(array_x,array_y,xlab=name_x,ylab=name_y,main="Correlation (random 2000)")
      contour(z, drawlabels=FALSE, nlevels=5, col=my.cols, add=TRUE, lwd=2)
      legend("bottomright",legend=corlegend)
    })

    output$savescatterProfilesAndBox=renderUI({downloadButton('savescatterProfilesAndBoxbutton', 'Get PDF')})

  }else{
    # toplot$profileAndBoxes$logicscatter=NULL
    output$scatterProfilesAndBox<-renderPlot({NULL})
    output$savescatterProfilesAndBox=renderUI({NULL})
  }
})








#observer for the partial correlation click

toListenParcor <- reactive({
    list(input$pcor_click$y,input$pcor_click$x,input$Log2BoxProfilesAndBox)
})

observeEvent(toListenParcor(),{
  set.seed(123)
  
  if(!is.null(input$pcor_click$x) & !is.null(input$pcor_click$y) ){
    toplot$profileAndBoxes$px=floor(as.numeric(input$pcor_click$x))+1
    toplot$profileAndBoxes$py=floor(as.numeric(input$pcor_click$y))+1     
  }
  portionlist_boxes=corvariables$portionlist_boxes 
  if(!is.null(portionlist_boxes) & length(toplot$profileAndBoxes$px)>0 & length(toplot$profileAndBoxes$py) >0 &length(portionlist_boxes)>2 ){

    mat=do.call(cbind,portionlist_boxes)
    #remove -Inf, because log2(0)=-Inf
    mat[is.infinite(mat) &mat<0 ]=0
    if(input$Log2BoxProfilesAndBox){
      mat=log2(mat)
    }
    #take a random sample (2000) of the total number of points
    if(nrow(mat)>2000){
      smpl=sort(sample(nrow(mat),2000,replace=FALSE))
      mat=mat[smpl,]
    }
    #select idx of the matrix
    name_x=names(portionlist_boxes)[toplot$profileAndBoxes$px]
    name_y=names(portionlist_boxes)[length(portionlist_boxes)-toplot$profileAndBoxes$py+1]
    mat[is.infinite(mat) &mat<0 ]=0
    idx1=which(name_x==colnames(mat))
    idx2=which(name_y==colnames(mat))

    model=plotpcor(mat,idx=c(idx1,idx2))

    array_x=model[[1]]
    array_y=model[[2]]

    if(corvariables$is.density){
      name_x=paste(name_x,"read density (rpm/bp)")
      name_y=paste(name_y,"read density (rpm/bp)")
    }else{
      name_x=paste(name_x,"normalized reads (rpm)")
      name_y=paste(name_y,"normalized reads (rpm)")
    }

    toplot$profileAndBoxes$array_x=array_x
    toplot$profileAndBoxes$array_y=array_y
    toplot$profileAndBoxes$name_x=name_x
    toplot$profileAndBoxes$name_y=name_y
    toplot$profileAndBoxes$x=NULL
    toplot$profileAndBoxes$y=NULL

    output$scatterProfilesAndBox<-renderPlot({
      my.cols <- rev(brewer.pal(5, "RdYlBu"))
      z <- kde2d(array_x, array_y, n=50)
      plot(array_x,array_y,xlab=name_x,ylab=name_y,main="Correlation (random 2000)")
      contour(z, drawlabels=FALSE, nlevels=5, col=my.cols, add=TRUE, lwd=2)
    })

    output$savescatterProfilesAndBox=renderUI({downloadButton('savescatterProfilesAndBoxbutton', 'Get PDF')})

  }else{
    # toplot$profileAndBoxes$logicscatter=NULL
    output$scatterProfilesAndBox<-renderPlot({NULL})
    output$savescatterProfilesAndBox=renderUI({NULL})
  }
})








###############################################################################
###############################################################################
###############################################################################
# observers for the dynamics
###############################################################################
###############################################################################
###############################################################################

#observer for the genelists available
observe({
  
  ROIvariables$listROI
  #read ROIs, and put only genelists, if both promoters/transcripts/TES of that geneset are found
  nomi=unlist(lapply(ROIvariables$listROI,getName))

  gr1=grepl("promoters_genelist_",nomi)
  gr2=grepl("transcripts_genelist_",nomi)
  gr3=grepl("TES_genelist_",nomi)
  nomi_promoters=nomi[gr1]
  nomi_transcripts=nomi[gr2]
  nomi_TES=nomi[gr3]

  #find common genelists
  genelist_promoters=substr(nomi_promoters,20,nchar(nomi_promoters))
  genelist_transcripts=substr(nomi_transcripts,22,nchar(nomi_transcripts))
  genelist_TES=substr(nomi_TES,14,nchar(nomi_TES))

  common=intersect(genelist_promoters,genelist_transcripts)
  common=intersect(common,genelist_TES)
  if(length(common)>0){
    #find length
    idx=match(paste("promoters_genelist_",common,sep=""),nomi)

    rois=ROIvariables$listROI[idx]
    
    lens=unlist(lapply(rois,getLength)) 
    #de-comment this if you want to add "genelist " before the name of each gene list
    #historylist=as.list(paste("genelist",common))
    historylist=as.list(common)
    names(historylist)=paste(historylist,lens)
    updateCheckboxGroupInput(session,inputId="genelistsforDynamics",label=NULL,
                                      choices = historylist)     
  }else{
    updateCheckboxGroupInput(session,inputId="genelistsforDynamics",label=NULL,
                                      choices = character(0))  
  }
   
})



#observer for the BAM files common to the genelists selected
observe({
  input$genelistsforDynamics
  
  if(length(ROIvariables$listROI)>0&length(input$genelistsforDynamics)>0){
    #if we are here, promoters, transcripts and TES are available for the same genelist
    #for each gene, find promoters, transcripts and TES and the common BAM files associated
    commonBAMs=as.list(rep(NA,length(input$genelistsforDynamics)))
    commonBAMs=lapply(commonBAMs,function(x) {return(NULL)})
    nomi=unlist(lapply(ROIvariables$listROI,getName))

    triadsROItoshow=c()


    for(i in 1:length(input$genelistsforDynamics)){
      #change this if in the menu "genelist " is before the actual name of the genelist
      #selected=strsplit(input$genelistsforDynamics[i],split=" ")[[1]][2]
      selected=input$genelistsforDynamics[i]
      pos_promoters= which(paste("promoters_genelist_",selected,sep="")==nomi)
      pos_transcripts= which(paste("transcripts_genelist_",selected,sep="")==nomi)
      pos_TES= which(paste("TES_genelist_",selected,sep="")==nomi)
      triadsROItoshow=c(triadsROItoshow,paste("promoters_genelist_",selected,sep=""),
                                        paste("transcripts_genelist_",selected,sep=""),
                                        paste("TES_genelist_",selected,sep=""))
      if (length(pos_promoters)>0 & length(pos_transcripts)>0 & length(pos_TES)){
        promoters_roi=ROIvariables$listROI[[pos_promoters]]
        transcripts_roi=ROIvariables$listROI[[pos_transcripts]]
        TES_roi=ROIvariables$listROI[[pos_TES]]

        #get BAM for this genelist, checking if the BAM is associated to all 
        #promoters,transcripts,TES
        rawvals_promoters=Enrichlist$rawcoverage[[pos_promoters]]  
        rawvals_transcripts=Enrichlist$rawcoverage[[pos_transcripts]]   
        rawvals_TES=Enrichlist$rawcoverage[[pos_TES]]
              
        BAM_promoters=names(rawvals_promoters)
        BAM_transcripts=names(rawvals_transcripts)      
        BAM_TES=names(rawvals_TES) 
        commonBAMgenelist=intersect(BAM_promoters,BAM_transcripts) 
        commonBAMgenelist=intersect(commonBAMgenelist,BAM_TES)
        if(!is.null(commonBAMgenelist)){
          commonBAMs[[i]]=commonBAMgenelist
        }        
      }else{
        commonBAMs[[i]]=NULL
      }
    }
    #show triads ROI for selected genelists
    output$showROItriadGeneList<-renderUI({
      list(
        HTML("<b>ROIs associated to selected gene lists:</b>"),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
          HTML(triadsROItoshow,collapse="<br>")
        )
      )
      
    }) 

    common=Reduce(intersect, commonBAMs)

    if(length(common)!=0){
      historylist=as.list(common)
      names(historylist)=common
      #then, find BAM files common to the genelists selected:    
      updateCheckboxGroupInput(session,inputId="BAMforDynamics",label=NULL,
                                        choices = historylist)  
    }else{
      updateCheckboxGroupInput(session,inputId="BAMforDynamics",label=NULL,
                                        choices = character(0))        
    }


  }else{
    updateCheckboxGroupInput(session,inputId="BAMforDynamics",label=NULL,
                                        choices = character(0))  
    #print nothing for ROI triads
    output$showROItriadGeneList<-renderUI({NULL})                                      
  }


})

