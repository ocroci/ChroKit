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
          output$BAMmenuchoose_singleeval<-renderUI({NULL})
        }else{
          output$BAMmenuchoose_singleeval<-renderUI({ 
          list(HTML("<b>2) Select associated enrichment:</b>"),
                selectInput(inputId="BAMchooseSingleEval", label=NULL,choices=getbam) ) })

        }
        
      }
    }

})


#observer for pie and bar single eval colors
observe({
  #respond to color change
  input$chooseColorPaletteSingleEval
  #set variables inside plot area for PDFs
  if(isvalid(input$chooseColorPaletteSingleEval)){
    toplot$viewDistributionPieSingleEval$Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]
  }
})

#observer for density plot single eval colors
observe({
  input$chooseColorPaletteSingleEval_density
  #set variables inside plot area for PDFs
  if(isvalid(input$chooseColorPaletteSingleEval_density)) {
    toplot$viewDistributionPieSingleEval$Colors_density=strsplit(input$chooseColorPaletteSingleEval_density,split="_")[[1]]
  }
})

#observer for box input fields (normalization and color)
observe({
  input$chooseColorPaletteSingleEval_box
  if(isvalid(input$chooseColorPaletteSingleEval_box)) {
    toplot$viewDistributionPieSingleEval$Colors_box=strsplit(input$chooseColorPaletteSingleEval_box,split="_")[[1]]
  }
})  
observe({
  input$chooseNormalizationSingleEval_box
  if(isvalid(input$chooseNormalizationSingleEval_box)) {
    toplot$viewDistributionPieSingleEval$normalization_box=strsplit(input$chooseNormalizationSingleEval_box,split="_")[[1]]
  }
})


#observer for profile input fields (profile to plot, color,ylab)
observe({
  input$chooseColorPaletteSingleEval_profile
  if(isvalid(input$chooseColorPaletteSingleEval_profile)) {
    toplot$viewDistributionPieSingleEval$Colors_profile=strsplit(input$chooseColorPaletteSingleEval_profile,split="_")[[1]]
  }
})

observe({
  input$chooseNormalizationSingleEval_profile
  if(isvalid(input$chooseNormalizationSingleEval_profile) & length(toplot$viewDistributionPieSingleEval$matrixes)>0) {
      normalization=input$chooseNormalizationSingleEval_profile
      matrixes=toplot$viewDistributionPieSingleEval$matrixes
      lengths_sampled=toplot$viewDistributionPieSingleEval$lengths_sampled
      if(normalization=="totread"){
        ylab="Total reads"
      }else{
        ylab="Read density (reads/bp)"
        #divide for length of ranges
        for (i in 1:length(matrixes)){
          if (length(lengths_sampled[[i]])>0){
            matrixes[[i]]=matrixes[[i]]/lengths_sampled[[i]]
          }
        }
      }
    profile_to_plot=lapply(matrixes,function(i) {apply(i,2,mean)})
    names(profile_to_plot)=c("all","promoters","genebody","intergenic")
    toplot$viewDistributionPieSingleEval$profile_to_plot=profile_to_plot
    toplot$viewDistributionPieSingleEval$ylabprofile=ylab

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

    rawvals1=Enrichlist$decryptkey[[pos1]]
    rawvals2=Enrichlist$decryptkey[[pos2]]

    if (!is.null(roi1)){
      getbam1=names(rawvals1)
      
      if (is.null(getbam1)){
        output$BAMmenuchooseCmp1<-renderUI({NULL})
        #updateSelectInput(session,inputId="BAM1chooseCmp",label=NULL,choices=character(0))
      }else{
        output$BAMmenuchooseCmp1<-renderUI({ 
          list( list(HTML("<b>Choose enrichment of ROI-1:</b>"),htmlhelp("","help_pairwiseOverlaps_parameters_enrich1")),
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
          list(list(HTML("<b>Choose enrichment of ROI-2:</b>"),htmlhelp("","help_pairwiseOverlaps_parameters_enrich2")),
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
    output$ROIsForDigitalHeat_menu<-renderUI({
      list(
        list(HTML("<b>ROIs to view:</b>"),htmlhelp("","help_digitalHeatmap_parameters_ROItoview")),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
          checkboxGroupInput("ROIsForDigitalHeat",NULL,choices=historylist)
        )
      )
    })

  }else{
    output$ROIsForDigitalHeat_menu<-renderUI({checkboxGroupInput("ROIsForDigitalHeat",NULL,choices=NULL)})
  }
})


#update for final button
observe({
  if (!isvalid(input$ROImaster)|!isvalid(input$ROIsForDigitalHeat)){
    output$show_confirmUpdateDigitalHeat1<-renderUI({NULL})
    return()
  }
  output$show_confirmUpdateDigitalHeat1<-renderUI({actionButton("confirmUpdateDigitalHeat1", "Update plot")})
})



#update list ROIs to reorder for digital heatmap
#reorder ROI
observe({
  input$ROIsForDigitalHeat
  #check if ROI are present, valid
  if(is.null(ROIvariables$listROI) | length(ROIvariables$listROI)<1 | length(input$ROIsForDigitalHeat)<1|
          length(input$ROImaster)<1 | is.null(input$ROImaster)){
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
      list(HTML("<b>ROIs ordering</b>"),htmlhelp("","help_digitalHeatmap_parameters_ROIordering")),
      wellPanel(id = "logPanelBAM",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
        fluidPage(
          lista
        )
      )
    ) 

  })

})




#update ROI to cluster digital heatmap
observe({
  input$ROIsForDigitalHeat
  
  if(length(input$ROIsForDigitalHeat)>0 & !is.null(input$ROImaster) & length(input$ROImaster)>0){
    output$ROIsForClusterDigital_menu<-renderUI ({
      list(
        list(HTML("<b>ROIs for cluster:</b>"),htmlhelp("","help_digitalHeatmap_parameters_ROIforcluster")),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 100px; max-width: 300px; background-color: #ffffff;",
          checkboxGroupInput("ROIforClusteringDigitalHeat",NULL,choices=input$ROIsForDigitalHeat)
        )
      )

    })

  
  }else{
    output$ROIsForClusterDigital_menu<-renderUI ({NULL})
  }

})


#update random sample to choose for digital heatmap:
observe({
  input$ROImaster
  input$ROIsForDigitalHeat
  #to show random sample, must choose master ROI and at least one ROI to show
  if(length(ROIvariables$listROI)>0 & isvalid(input$ROImaster) & isvalid(input$ROIsForDigitalHeat)){
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

      output$sampleRandomROIDigital<-renderUI({
        list(
          list(HTML("<b>Random sample of genomic ranges to show:</b>"),htmlhelp("","help_digitalHeatmap_parameters_randomsample")),
          numericInput(inputId = 'sampleRandomDigitalHeat',label=NULL,min = minim, max = lbr, step = 1000,value=valshown)
        )
      })

    }else{
      output$sampleRandomROIDigital<-renderUI({NULL})
              
    }    
  }else{
    output$sampleRandomROIDigital<-renderUI({NULL})
  }

})


##observer for type of clustering and clustering parameters
#for the type of clustering, react to "ROIs for cluster" (input$ROIforClusteringDigitalHeat)
observe({
  
  input$ROIforClusteringDigitalHeat
  #if something inside, show the choise between hclust and kmeans
  #otherwise, hide (change in no selection or not selected from beginning)
  if(length(input$ROIforClusteringDigitalHeat)>0 & length(input$ROIsForDigitalHeat)>0 & length(input$ROImaster)>0){
    output$clustertypeDigitalHeat<-renderUI({
      radioButtons("clusterTypeDigitalHeat",label="Clustering type",
                    choiceNames=list(
                      htmlhelp("K-means","help_digitalHeatmap_parameters_clusterKmeans"),
                      htmlhelp("hierarchical","help_digitalHeatmap_parameters_clusterHierarchical")),
                    choiceValues=list("kmean","hierarchical")
                  )
    }) 
    output$clusternumbershowDigitalHeat<-renderUI({
      #must implement a check in serverGENOMICS to correct if this number <=0 or > length roi shown
      numericInput(inputId = 'clustnumDigitalHeat',label=list("Cluster number",htmlhelp("","help_digitalHeatmap_parameters_clusternumber")),min = 1, max = 433, step = 1,value=4)
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
  if(length(input$ROIforClusteringDigitalHeat)>0 & !is.null(input$clusterTypeDigitalHeat) 
        & length(input$ROIsForDigitalHeat)>0 & length(input$ROImaster)>0){
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
      toplot$digital$colorpalettes="white_red4"
      
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
      # #find all possible colors selected.
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


#observer for strand-specific overlap (appear only when have ROI to show selected)
observe({
  input$ROIsForDigitalHeat
  input$ROImaster
  if (!isvalid(input$ROIsForDigitalHeat)|!isvalid(input$ROImaster)){
    output$showStrandSpecOverlap<-renderUI({NULL})
    return()
  }
  output$showStrandSpecOverlap<-renderUI({
    checkboxInput("StrandSpecOverlap", label=list("Strand-specific overlaps",htmlhelp("","help_digitalHeatmap_parameters_strandspecific")),value = FALSE, width = NULL)
  })
})






#observer for the number of bins? (check if nbins > width of some range)
observe({
  
  input$ROIsForDigitalHeat
  input$ROImaster
  if (length(ROIvariables$listROI)>0 & length(input$ROIsForDigitalHeat)>0 & isvalid(input$ROImaster)){
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

      output$showbinsDigitalHeat<-renderUI({
        list(
          list(HTML("<b>Number of bins:</b>"),htmlhelp("","help_digitalHeatmap_parameters_nbins")),
          numericInput("binsDigitalHeat",label=NULL,min = 1, max = maxtoshow,value=valuetoshow,step = 1)  
        )

      })
          

    }else{
      output$showbinsDigitalHeat<-renderUI({NULL})
    }
  }else{
    output$showbinsDigitalHeat<-renderUI({NULL})
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
        output$showBAMsForAnalogHeat<-renderUI({
          list(
            list(HTML("<b>Select enrichments to show:</b>"),htmlhelp("","help_analogicHeatmap_parameters_enrichments")),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("BAMsForAnalogHeat",NULL,choices=finalBAMs)
            )
          )
        })
      }else{
        output$showBAMsForAnalogHeat<-renderUI({
          list(
            HTML("<font color='red'>Some of the ROI(s) selected do not have enrichment associated.\nGo to 'ROI preparation' to associate enrichments to ROIs</font>"),
            checkboxGroupInput("BAMsForAnalogHeat",NULL,choices=NULL)
          )
        })
        #output$showBAMsForAnalogHeat<-renderUI({NULL})
        output$reorderBAMmenuAnalogHeat<-renderUI({NULL})
        output$showrankingmethod<-renderUI({NULL})
        output$clustertypeAnalogHeat<-renderUI({NULL})
        output$orderingAnalogHeat<-renderUI({NULL})
        output$clusternumbershowAnalogHeat<-renderUI({NULL})
        output$clusterKstartsAnalogHeat<-renderUI({NULL})
        output$clusterKiterationsAnalogHeat<-renderUI({NULL})
        output$clusterHDistMethodAnalogHeat<-renderUI({NULL})
        output$clusterHClustMethodAnalogHeat<-renderUI({NULL})
        output$showbinsAnalogHeat<-renderUI({NULL})
        output$showsampleRandomAnalogHeat<-renderUI({NULL})
      }
 
    }
  }else{
    output$showBAMsForAnalogHeat<-renderUI({checkboxGroupInput("BAMsForAnalogHeat",NULL,choices=NULL)})
    #output$showBAMsForAnalogHeat<-renderUI({NULL})
    output$reorderBAMmenuAnalogHeat<-renderUI({NULL})
    output$showrankingmethod<-renderUI({NULL})
    output$clustertypeAnalogHeat<-renderUI({NULL})
    output$orderingAnalogHeat<-renderUI({NULL})
    output$clusternumbershowAnalogHeat<-renderUI({NULL})
    output$clusterKstartsAnalogHeat<-renderUI({NULL})
    output$clusterKiterationsAnalogHeat<-renderUI({NULL})
    output$clusterHDistMethodAnalogHeat<-renderUI({NULL})
    output$clusterHClustMethodAnalogHeat<-renderUI({NULL})
    output$showbinsAnalogHeat<-renderUI({NULL})
    output$showsampleRandomAnalogHeat<-renderUI({NULL})
  }
})



#here, choose ordering of BAM files that will appear in the menu.
#should depend on BAMs selected for analogic heatmap


observe({ 
  input$ROIsForAnalogHeat
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
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  pos=match(input$ROIsForAnalogHeat,nomi)
  #all rois selected
  roi=ROIvariables$listROI[pos]
  rawvals=Enrichlist$rawcoverage[pos]    

  #verify that this roi have selected BAM files. input could be remained from previous analysis




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
      list(HTML("<b>Enrichment order</b>"),htmlhelp("","help_analogicHeatmap_parameters_enrichmentorder")),
      wellPanel(id = "logPanelBAM",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
        fluidPage(
          lista
        )
      )
    ) 

  })   
})



#appear menu of choosing order (rank or cluster) only if some BAM selected
observe({
  #checks: there must be some ROI
  if(is.null(ROIvariables$listROI) | length(ROIvariables$listROI)<1 | length(input$ROIsForAnalogHeat)<1){
    output$showrankingmethod<-renderUI({NULL})
    return()
  }

  #checks: some BAMs must be selected (>0)
  if(length(input$BAMsForAnalogHeat)<1 ){
    output$showrankingmethod<-renderUI({NULL})
    return()    
  }

  output$showrankingmethod<-renderUI({
    list(
      HTML("<b>Clustering/ranking</b>"),
      radioButtons("chooseOrderingAnalogHeat",label=NULL,
                            choiceNames=list(
                              htmlhelp("Ranking","help_analogicHeatmap_parameters_ranking"),
                              htmlhelp("Clustering","help_analogicHeatmap_parameters_clustering")
                            ),
                            choiceValues=list("ranking","clustering")
                  )
    )
  })

})




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


  if (!isvalid(input$chooseOrderingAnalogHeat)|!isvalid(input$ROIsForAnalogHeat)){
    output$clustertypeAnalogHeat<-renderUI({NULL})
    output$orderingAnalogHeat<-renderUI({NULL})
    output$clusternumbershowAnalogHeat<-renderUI({NULL})
    return()
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
          numericInput(inputId = 'clustnumAnalogHeat',label=list("Cluster number",htmlhelp("","help_analogicHeatmap_parameters_clusternumber")),min = 1, max = numranges, step = 1,value=4)
        }) 
        output$clustertypeAnalogHeat<-renderUI({
          radioButtons("clusterTypeAnalogHeat",label="Clustering type",
            choiceNames=list(
              htmlhelp("K-means","help_analogicHeatmap_parameters_clusterKmeans"),
              htmlhelp("hierarchical","help_analogicHeatmap_parameters_clusterHierarchical")
            ),
            choiceValues=c("kmean","hierarchical")
          )
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
          checkboxGroupInput("BAMsForClusteringAnalogHeat",NULL,choices=NULL)
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


  if (!isvalid(input$chooseOrderingAnalogHeat)|!isvalid(input$ROIsForAnalogHeat)){
    output$clusterKstartsAnalogHeat<-renderUI({NULL})
    output$clusterKiterationsAnalogHeat<-renderUI({NULL})
    output$clusterHDistMethodAnalogHeat<-renderUI({NULL})
    output$clusterHClustMethodAnalogHeat<-renderUI({NULL})  
    return()    
  }

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
observeEvent(input$optioncolorsforAnalogHeat,{

  if (!isvalid(input$optioncolorsforAnalogHeat)){
    output$showcolorsheat<-renderUI({NULL})
    return()
  }

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


})





#observer for chosen color in analog heat.
observe({
  optionMenu=input$optioncolorsforAnalogHeat
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
      ##find all possible colors selected.
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
  #if no BAM selected, useless to show the bins
  if (!isvalid(input$BAMsForAnalogHeat)){
    output$showbinsAnalogHeat<-renderUI({NULL})
    return()
  }
  
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

      output$showbinsAnalogHeat<-renderUI({
        list(
          list(HTML("<b>Number of bins:</b>"),htmlhelp("","help_analogicHeatmap_parameters_nbins")),
          numericInput("binsAnalogHeat",label=NULL,min = 1, max = maxtoshow,value=valuetoshow,step = 1)
        )
        
      })
      

    }else{
      output$showbinsAnalogHeat<-renderUI({NULL})
    }

  }else{
    output$showbinsAnalogHeat<-renderUI({NULL})
    #updateNumericInput(session,"binsAnalogHeat",label=NULL,min = 1, max = 500,value=50,step = 1)
  }
})





#update possible random sample, from 1000 to length(sum(ROIs)), with step=1000
#depends only on the ROI selected (input$ROIsForAnalogHeat)
observe({
  
  input$ROIsForAnalogHeat
  #if no BAM selected, useless to show random sample to choose

  if (!isvalid(input$BAMsForAnalogHeat)){
    output$showsampleRandomAnalogHeat<-renderUI({NULL})
    return()
  }
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

    output$showsampleRandomAnalogHeat<-renderUI({
      numericInput(inputId = 'sampleRandomAnalogHeat',label=list("Random sample of:",htmlhelp("","help_analogicHeatmap_parameters_subsample")),min = minim, max = lbr, step = 1000,value=valshown)
    })
    
  }else{
    output$showsampleRandomAnalogHeat<-renderUI({NULL})
            
  }
})


#observer for button: it will appear only if valid ROI/BAM selected
#and, if clustering, some ROI have been selected for clustering
observe({
  if (!isvalid(input$BAMsForAnalogHeat)|!isvalid(input$ROIsForAnalogHeat)|!isvalid(input$chooseOrderingAnalogHeat)){
    output$show_confirmUpdateAnalogHeat<-renderUI({NULL})
    return()
  }  
  if(input$chooseOrderingAnalogHeat=="clustering" & length(input$BAMsForClusteringAnalogHeat)==0){
    output$show_confirmUpdateAnalogHeat<-renderUI({NULL})
    return()
  }
  output$show_confirmUpdateAnalogHeat<-renderUI({actionButton("confirmUpdateAnalogHeat", "Update plot")})

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
      # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      # color_distinct <-sample(col_vector, n)

      color_distinct=colors_list[1:n]
      toplot$gadgetanalogic$color_distinct=color_distinct
      toplot$gadgetanalogic$Log2BoxAnalogHeat=input$Log2BoxAnalogHeat
      toplot$gadgetanalogic$portionlist_profile=portionlist_profile
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


      #options specific to profile

      output$showprofileAnalogHeat_logOptions<-renderUI({
        checkboxInput("isLog2_profileAnalogHeat", label="log2",value = FALSE, width = NULL)
      })

      output$showprofileAnalogHeat_colorschemeOptions<-renderUI({
        radioButtons("colorscheme_profileAnalogHeat",label="Choose color scheme:",choices=c(
                                                  "Random colors"="random",
                                                  "Custom colors"="custom"
                                                        ),selected="random")       
      })
      #if custom, menu will be updated in an observer


      #here, legends are ROIs or clusters, depending on bycluster variable. If TRUE, clsuter 1,2,3,
      #otherwise, ROI1,2,3 and in the legend we need the explanation: ROI1= blablabla
      #plot the profiles of regions selected
      output$profileAnalogHeat<-renderPlot({
        #check specific inputs
        if(isvalid(input$isLog2_profileAnalogHeat)){
          islog=input$isLog2_profileAnalogHeat
        }else{  
          islog=FALSE
        }

        if (isvalid(input$colorscheme_profileAnalogHeat)){
          if(input$colorscheme_profileAnalogHeat=="custom"){
            #only the first if custom
            if (isvalid(input$colorCustomAnalogHeat_profile1)){
              #extract colors based on choice (loop through them):
              tosearch=paste("colorCustomAnalogHeat_profile",1:length(portionlist_profile),sep="")              
              if(length(tosearch)>0){
                listinputcols=list()
                for(i in 1:length(tosearch)){
                  listinputcols[[i]]=input[[tosearch[i]]]
                }
                listinputcols=unlist(listinputcols)
                if(length(listinputcols)==length(tosearch)){
                  colors=listinputcols
                }else{
                  colors=color_distinct
                }
              }else{
                colors=color_distinct
              }
            }else{
              colors=color_distinct
            }
          }else{
            colors=color_distinct
          }
        }else{
          colors=color_distinct
        }
        plot_analog_profile(islog2=islog,portionlist_profile=portionlist_profile,colors=colors)
      })




      #options specific to boxplots
      output$showboxAnalogHeat_logOptions<-renderUI({
        checkboxInput("isLog2_boxAnalogHeat", label="log2",value = FALSE, width = NULL)
      })

      output$showboxAnalogHeat_colorschemeOptions<-renderUI({
        radioButtons("colorscheme_boxAnalogHeat",label="Choose color scheme:",choices=c(
                                                  "Random colors"="random",
                                                  "Custom colors"="custom"
                                                        ),selected="random")       
      })




      #plot boxplots of regions selected. Use portionlist_boxes
      output$boxplotByBAMAnalogHeat<-renderPlot({
          if (isvalid(input$isLog2_boxAnalogHeat)){
            islog=input$isLog2_boxAnalogHeat
          }else{
            islog=FALSE
          }
          if (isvalid(input$GroupColorsAnalogHeat_box)){
            grouped=input$GroupColorsAnalogHeat_box
          }else{
            grouped=FALSE
          }
          #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
          #invert the order
          newlist=list()
          newcols=c()
          newnames=c()
          if (isvalid(input$colorscheme_boxAnalogHeat)){
            #if random, check if group or not
            if (input$colorscheme_boxAnalogHeat=="random"){
              #if grouping colors, adjust legend and colors by groups
              if(!grouped){
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
            }else{
              #extract all colors,because color scheem is custom
              tosearch=paste("colorCustomAnalogHeat_box",1:length(portionlist_boxes),sep="")              
              if(length(tosearch)>0){
                listinputcols=list()
                for(i in 1:length(tosearch)){
                  listinputcols[[i]]=input[[tosearch[i]]]
                }
                listinputcols=unlist(listinputcols)
                if(length(listinputcols)==length(tosearch)){
                  newcols=listinputcols
                  for(i in 1:bamnumber){
                    pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                    newlist=c(newlist,portionlist_boxes[pos_inverted])
                    newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                  }
                }else{
                  for(i in 1:bamnumber){
                    pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                    newlist=c(newlist,portionlist_boxes[pos_inverted])
                    newcols=c(newcols,color_distinct[pos_inverted])
                    newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                  }
                }
              }else{
                for(i in 1:bamnumber){
                  pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                  newlist=c(newlist,portionlist_boxes[pos_inverted])
                  newcols=c(newcols,color_distinct[pos_inverted])
                  newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                }
              }
            }            
          }else{
            #here option is not valid=> default, random colors
            for(i in 1:bamnumber){
              pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
              newlist=c(newlist,portionlist_boxes[pos_inverted])
              newcols=c(newcols,color_distinct[pos_inverted])
              newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
            } 

          }

          plot_analog_boxByBAM(materialtoplot=newlist,roinumber=roinumber,newcols=newcols,newnames=newnames,
                                bamname=toplot$gadgetanalogic$bamname,islog=islog,colors=newcols)

      })



      output$boxplotByROIAnalogHeat<-renderPlot({

          if (isvalid(input$isLog2_boxAnalogHeat)){
            islog=input$isLog2_boxAnalogHeat
          }else{
            islog=FALSE
          }
          if (isvalid(input$GroupColorsAnalogHeat_box)){
            grouped=input$GroupColorsAnalogHeat_box
          }else{
            grouped=FALSE
          }
          #if grouping colors, adjust legend and colors by groups
          if (isvalid(input$colorscheme_boxAnalogHeat)){

            if (input$colorscheme_boxAnalogHeat=="random"){
              if(!grouped){
                newcols=color_distinct
                newnames=names(portionlist_boxes)
              }else{
                newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
                newcols=rep(color_distinct[1:toplot$gadgetanalogic$bamnumber],toplot$gadgetanalogic$roinumber)
              }

            }else{
              grouped=FALSE
              #here is valid color scheme but not random: extract single colors
              tosearch=paste("colorCustomAnalogHeat_box",1:length(portionlist_boxes),sep="")              
              if(length(tosearch)>0){
                listinputcols=list()
                for(i in 1:length(tosearch)){
                  listinputcols[[i]]=input[[tosearch[i]]]
                }
                listinputcols=unlist(listinputcols)
                if(length(listinputcols)==length(tosearch)){
                  newcols=listinputcols
                  newnames=names(portionlist_boxes)
                }else{
                  newcols=color_distinct
                  newnames=names(portionlist_boxes)
                }
              }else{
                newcols=color_distinct
                newnames=names(portionlist_boxes)
              }
            }
          }else{
            grouped=FALSE
            newcols=color_distinct
            newnames=names(portionlist_boxes)
          }

          plot_analog_boxByROI(materialtoplot=portionlist_boxes,roiname=toplot$gadgetanalogic$roiname,
                              bamname=toplot$gadgetanalogic$bamname,newnames=newnames,islog=islog,
                              isgrouped=grouped,colors=newcols)

      })  





      #buttons for download PDF of profile, boxplot by ROI, boxplot by BAM
      output$saveprofileAnalogHeat=renderUI({downloadButton('saveprofileAnalogHeatbutton', 'Get PDF')})
      output$saveboxplotByROIAnalogHeat=renderUI({downloadButton('saveboxplotByROIAnalogHeatbutton', 'Get PDF')})
      output$saveboxplotByBAMAnalogHeat=renderUI({downloadButton('saveboxplotByBAMAnalogHeatbutton', 'Get PDF')})

      #plot correlation heatmap (and partial correlation?)
      # if(roinumber==1 & length(portionlist_boxes)>=2 & (ytop-ybottom)>1){
      #   if(input$Log2BoxAnalogHeat){
      #     portionlist_boxescors=lapply(portionlist_boxes,log2)
      #   }else{
      #     portionlist_boxescors=portionlist_boxes
      #   }
      #   mat=do.call(cbind,portionlist_boxescors)
        
      #   #problem: when playing with heatmap, Error in colnames<-: length of 'dimnames' [2] not equal to array extent.
      #   #maybe a simple trycatch would do the job
      #   tryCatch({
      #     colnames(mat)=bamselected[xleft:xright]
      #   },
      #   warning = function( w ){
      #     colnames(mat)=1:ncol(mat)
      #   },
      #   error = function( err ){
      #     colnames(mat)=1:ncol(mat)
      #   }
      #   )
        
      #   mat[is.infinite(mat) &mat<0 ]=0
      #   correlation_total=cor(mat)
      #   trasp_cor=t(correlation_total)
      #   brk=c( seq( -1 , 1,0.01))
      #   my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

      #   toplot$gadgetanalogic$trasp_cor=trasp_cor
      #   toplot$gadgetanalogic$my_palette=my_palette
      #   toplot$gadgetanalogic$brk=brk
      #   toplot$gadgetanalogic$correlation_total=correlation_total
        

      #   output$corAnalogHeat<-renderPlot({
      #     par(mar=c(12,12,1,1),xpd=TRUE)
      #     image(0:nrow(trasp_cor), 0:ncol(trasp_cor),trasp_cor[,ncol(trasp_cor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
      #     axis( 2, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= rev(colnames( trasp_cor )), las= 2 )
      #     axis( 1, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= colnames( trasp_cor ), las= 2 )
      #     for (x in (nrow(correlation_total)-1+0.5):0.5  )
      #       for (y in 0.5: ((ncol(correlation_total)-1+0.5)   ))
      #         text(y,x, round(correlation_total[ncol(correlation_total)-x+0.5,y+0.5],2),col="blue")
      #   })
      #   output$savecorAnalogHeat=renderUI({downloadButton('savecorAnalogHeatbutton', 'Get PDF')})

      #   if (length(portionlist_boxescors)>2 & (ytop-ybottom)>1){
      #     #if number of BAMs is >2, calculate the partial correlation too
      #     correlation_partial=pcor(mat)$estimate
      #     #warning: The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero
      #     colnames(correlation_partial)=rownames(correlation_partial)=colnames(mat)
      #     trasp_pcor=t(correlation_partial)
      #     toplot$gadgetanalogic$trasp_pcor=trasp_pcor
      #     toplot$gadgetanalogic$correlation_partial=correlation_partial
      #     output$pcorAnalogHeat<-renderPlot({
      #       par(mar=c(12,12,1,1),xpd=TRUE)
      #       image(0:nrow(trasp_pcor), 0:ncol(trasp_pcor),trasp_pcor[,ncol(trasp_pcor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
      #       axis( 2, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= rev(colnames( trasp_pcor )), las= 2 )
      #       axis( 1, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= colnames( trasp_pcor ), las= 2 )
      #       for (x in (nrow(correlation_partial)-1+0.5):0.5  )
      #         for (y in 0.5: ((ncol(correlation_partial)-1+0.5)   ))
      #           text(y,x, round(correlation_partial[ncol(correlation_partial)-x+0.5,y+0.5],2),col="blue")
      #     }) 
      #     output$savepcorAnalogHeat=renderUI({downloadButton('savepcorAnalogHeatbutton', 'Get PDF')})         
      #   }else{
      #     #number of BAM==2 => no pcor
      #     output$pcorAnalogHeat<-renderPlot({NULL})
      #     output$savepcorAnalogHeat=renderUI({NULL})
      #   }

      # }else{
      #   output$corAnalogHeat<-renderPlot({NULL})
      #   output$pcorAnalogHeat<-renderPlot({NULL})
      #   output$savecorAnalogHeat=renderUI({NULL})
      #   output$savepcorAnalogHeat=renderUI({NULL})
      # }
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
      output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
      output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
      output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
      output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
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
    output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
    output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
    output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
    output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
    output$saveboxplotByROIAnalogHeat=renderUI({NULL})
    output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
    output$corAnalogHeat<-renderPlot({NULL})
    output$pcorAnalogHeat<-renderPlot({NULL})
    output$savecorAnalogHeat=renderUI({NULL})
    output$savepcorAnalogHeat=renderUI({NULL})
  }
 # paste0("cella x1=", xleft, "\ncella y1=", ybottom,"\ncella x2=", xright, "\ncella y2=",ytop)
})









#respond to color scheme of profile in Analog heat
observe({
  colorscheme=input$colorscheme_profileAnalogHeat
  if (isvalid(colorscheme)){
    #here is something. appear menu only if custom
    if(colorscheme=="custom"){

      bams=names(isolate(toplot$gadgetanalogic$portionlist_profile))
      lista=list()
      #bams are the names of the lines : ROI; enrichment
      if(length(bams)>0){
        output$showprofileAnalogHeat_colorlistOptions<-renderUI({
          for (i in 1:length(bams)){
            lista[[i]]=fluidRow(
                                column(12,
                                      colorSelectorInput(inputId=paste("colorCustomAnalogHeat_profile",i,sep=""),label=bams[i],choices=ColsArray,
                                                    selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                )
                        )
          }
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; background-color: #ffffff;",
            lista
          )
        })  
      }else{
        output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
      }
    }else{
      output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
    }
  }else{  
    #is null, => no menu
    output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
  }
})



#respond to color scheme of box in Analog heat
observe({
  colorscheme=input$colorscheme_boxAnalogHeat
  if (isvalid(colorscheme)){
    #here is something. appear menu only if custom
    if(colorscheme=="custom"){
      #if custom, cannot group colors
      output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
      bams=names(isolate(toplot$gadgetanalogic$portionlist_boxes))
      lista=list()
      #bams are the names of the lines : ROI; enrichment
      if(length(bams)>0){
        output$showboxAnalogHeat_colorlistOptions<-renderUI({
          for (i in 1:length(bams)){
            lista[[i]]=fluidRow(
                                column(12,
                                      colorSelectorInput(inputId=paste("colorCustomAnalogHeat_box",i,sep=""),label=bams[i],choices=ColsArray,
                                                    selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                )
                        )
          }
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; background-color: #ffffff;",
            lista
          )
        })  
      }else{
        output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
      }
    }else{
      output$showboxAnalogHeat_groupcolOptions<-renderUI({
        checkboxInput("GroupColorsAnalogHeat_box", label="Group colors",value = FALSE, width = NULL)
      })
      output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
    }
  }else{  
    #is null, => no menu
    output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
    output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
  }
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
    if(isolate(toplot$analogic$chooseOrderingAnalogHeat)=="clustering" & length(heatvariables$ROIsForAnalogHeat)==1){
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
        # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        # color_distinct <-sample(col_vector, n)

        color_distinct=colors_list[1:n]

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
        #options specific to profile
        output$showprofileAnalogHeat_logOptions<-renderUI({
          checkboxInput("isLog2_profileAnalogHeat", label="log2",value = FALSE, width = NULL)
        })

        output$showprofileAnalogHeat_colorschemeOptions<-renderUI({
          radioButtons("colorscheme_profileAnalogHeat",label="Choose color scheme:",choices=c(
                                                    "Random colors"="random",
                                                    "Custom colors"="custom"
                                                          ),selected="random")       
        })



        #if custom, menu will be updated in an observer
        #plot the profiles of regions selected
        output$profileAnalogHeat<-renderPlot({
          #check specific inputs
          if(isvalid(input$isLog2_profileAnalogHeat)){
            islog=input$isLog2_profileAnalogHeat
          }else{  
            islog=FALSE
          }

          if (isvalid(input$colorscheme_profileAnalogHeat)){
            if(input$colorscheme_profileAnalogHeat=="custom"){
              #only the first if custom
              if (isvalid(input$colorCustomAnalogHeat_profile1)){
                #extract colors based on choice (loop through them):
                tosearch=paste("colorCustomAnalogHeat_profile",1:length(portionlist_profile),sep="")              
                if(length(tosearch)>0){
                  listinputcols=list()
                  for(i in 1:length(tosearch)){
                    listinputcols[[i]]=input[[tosearch[i]]]
                  }
                  listinputcols=unlist(listinputcols)
                  if(length(listinputcols)==length(tosearch)){
                    colors=listinputcols
                  }else{
                    colors=color_distinct
                  }
                }else{
                  colors=color_distinct
                }
              }else{
                colors=color_distinct
              }
            }else{
              colors=color_distinct
            }
          }else{
            colors=color_distinct
          }
          plot_analog_profile(islog2=islog,portionlist_profile=portionlist_profile,colors=colors)
  
        })


        #options specific to boxplots
        output$showboxAnalogHeat_logOptions<-renderUI({
          checkboxInput("isLog2_boxAnalogHeat", label="log2",value = FALSE, width = NULL)
        })

        output$showboxAnalogHeat_colorschemeOptions<-renderUI({
          radioButtons("colorscheme_boxAnalogHeat",label="Choose color scheme:",choices=c(
                                                    "Random colors"="random",
                                                    "Custom colors"="custom"
                                                          ),selected="random")       
        })



        #plot boxplots of regions selected. Use portionlist_boxes
        output$boxplotByBAMAnalogHeat<-renderPlot({
          if (isvalid(input$isLog2_boxAnalogHeat)){
            islog=input$isLog2_boxAnalogHeat
          }else{
            islog=FALSE
          }
          if (isvalid(input$GroupColorsAnalogHeat_box)){
            grouped=input$GroupColorsAnalogHeat_box
          }else{
            grouped=FALSE
          }
          #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
          #invert the order
          newlist=list()
          newcols=c()
          newnames=c()
          if (isvalid(input$colorscheme_boxAnalogHeat)){
            #if random, check if group or not
            if (input$colorscheme_boxAnalogHeat=="random"){
              #if grouping colors, adjust legend and colors by groups
              if(!grouped){
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
            }else{
              #extract all colors,because color scheem is custom
              tosearch=paste("colorCustomAnalogHeat_box",1:length(portionlist_boxes),sep="")              
              if(length(tosearch)>0){
                listinputcols=list()
                for(i in 1:length(tosearch)){
                  listinputcols[[i]]=input[[tosearch[i]]]
                }
                listinputcols=unlist(listinputcols)
                if(length(listinputcols)==length(tosearch)){
                  newcols=listinputcols
                  for(i in 1:bamnumber){
                    pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                    newlist=c(newlist,portionlist_boxes[pos_inverted])
                    newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                  }
                }else{
                  for(i in 1:bamnumber){
                    pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                    newlist=c(newlist,portionlist_boxes[pos_inverted])
                    newcols=c(newcols,color_distinct[pos_inverted])
                    newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                  }
                }
              }else{
                for(i in 1:bamnumber){
                  pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                  newlist=c(newlist,portionlist_boxes[pos_inverted])
                  newcols=c(newcols,color_distinct[pos_inverted])
                  newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                }
              }
            }            
          }else{
            #here option is not valid=> default, random colors
            for(i in 1:bamnumber){
              pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
              newlist=c(newlist,portionlist_boxes[pos_inverted])
              newcols=c(newcols,color_distinct[pos_inverted])
              newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
            } 

          }

          plot_analog_boxByBAM(materialtoplot=newlist,roinumber=roinumber,newcols=newcols,newnames=newnames,
                                bamname=toplot$gadgetanalogic$bamname,islog=islog,colors=newcols)

        })




        output$boxplotByROIAnalogHeat<-renderPlot({
          if (isvalid(input$isLog2_boxAnalogHeat)){
            islog=input$isLog2_boxAnalogHeat
          }else{
            islog=FALSE
          }
          if (isvalid(input$GroupColorsAnalogHeat_box)){
            grouped=input$GroupColorsAnalogHeat_box
          }else{
            grouped=FALSE
          }

          #if grouping colors, adjust legend and colors by groups
          if (isvalid(input$colorscheme_boxAnalogHeat)){

            if (input$colorscheme_boxAnalogHeat=="random"){
              if(!grouped){
                newcols=color_distinct
                newnames=names(portionlist_boxes)
              }else{
                newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
                newcols=rep(color_distinct[1:toplot$gadgetanalogic$bamnumber],toplot$gadgetanalogic$roinumber)
              }

            }else{
              grouped=FALSE
              #here is valid color scheme but not random: extract single colors
              tosearch=paste("colorCustomAnalogHeat_box",1:length(portionlist_boxes),sep="")              
              if(length(tosearch)>0){
                listinputcols=list()
                for(i in 1:length(tosearch)){
                  listinputcols[[i]]=input[[tosearch[i]]]
                }
                listinputcols=unlist(listinputcols)
                if(length(listinputcols)==length(tosearch)){
                  newcols=listinputcols
                  newnames=names(portionlist_boxes)
                }else{
                  newcols=color_distinct
                  newnames=names(portionlist_boxes)
                }
              }else{
                newcols=color_distinct
                newnames=names(portionlist_boxes)
              }
            }
          }else{
            grouped=FALSE
            newcols=color_distinct
            newnames=names(portionlist_boxes)
          }
          
          plot_analog_boxByROI(materialtoplot=portionlist_boxes,roiname=ROIname,
                                bamname=toplot$gadgetanalogic$bamname,newnames=newnames,islog=islog,
                                isgrouped=grouped,colors=newcols)
        })  





        #buttons for download PDF of profile, boxplot by ROI, boxplot by BAM
        output$saveprofileAnalogHeat=renderUI({downloadButton('saveprofileAnalogHeatbutton', 'Get PDF')})
        output$saveboxplotByROIAnalogHeat=renderUI({downloadButton('saveboxplotByROIAnalogHeatbutton', 'Get PDF')})
        output$saveboxplotByBAMAnalogHeat=renderUI({downloadButton('saveboxplotByBAMAnalogHeatbutton', 'Get PDF')})
    
        #plot correlation heatmap (and partial correlation?)
        # if(roinumber==1 & length(portionlist_boxes)>=2& length(value)>1){
        #   if(input$isLog2_boxAnalogHeat){
        #     portionlist_boxescors=lapply(portionlist_boxes,log2)
        #   }else{
        #     portionlist_boxescors=portionlist_boxes
        #   }
        #   mat=do.call(cbind,portionlist_boxescors)
          
        #   #problem: when playing with heatmap, Error in colnames<-: length of 'dimnames' [2] not equal to array extent.
        #   #maybe a simple trycatch would do the job
        #   tryCatch({
        #     colnames(mat)=bamselected
        #   },
        #   warning = function( w ){
        #     colnames(mat)=1:ncol(mat)
        #   },
        #   error = function( err ){
        #     colnames(mat)=1:ncol(mat)
        #   }
        #   )
          
        #   mat[is.infinite(mat) &mat<0 ]=0
        #   correlation_total=cor(mat)
        #   trasp_cor=t(correlation_total)
        #   brk=c( seq( -1 , 1,0.01))
        #   my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

        #   toplot$gadgetanalogic$trasp_cor=trasp_cor
        #   toplot$gadgetanalogic$my_palette=my_palette
        #   toplot$gadgetanalogic$brk=brk
        #   toplot$gadgetanalogic$correlation_total=correlation_total
          

        #   output$corAnalogHeat<-renderPlot({
        #     par(mar=c(12,12,1,1),xpd=TRUE)
        #     image(0:nrow(trasp_cor), 0:ncol(trasp_cor),trasp_cor[,ncol(trasp_cor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
        #     axis( 2, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= rev(colnames( trasp_cor )), las= 2 )
        #     axis( 1, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= colnames( trasp_cor ), las= 2 )
        #     for (x in (nrow(correlation_total)-1+0.5):0.5  )
        #       for (y in 0.5: ((ncol(correlation_total)-1+0.5)   ))
        #         text(y,x, round(correlation_total[ncol(correlation_total)-x+0.5,y+0.5],2),col="blue")
        #   })
        #   output$savecorAnalogHeat=renderUI({downloadButton('savecorAnalogHeatbutton', 'Get PDF')})


        #   if (length(portionlist_boxescors)>2 & length(value)>1){
        #     #if number of BAMs is >2, calculate the partial correlation too
        #     correlation_partial=pcor(mat)$estimate
        #     #warning: The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero
        #     colnames(correlation_partial)=rownames(correlation_partial)=colnames(mat)
        #     trasp_pcor=t(correlation_partial)

        #     toplot$gadgetanalogic$trasp_pcor=trasp_pcor
        #     toplot$gadgetanalogic$correlation_partial=correlation_partial

        #     output$pcorAnalogHeat<-renderPlot({
        #       par(mar=c(12,12,1,1),xpd=TRUE)
        #       image(0:nrow(trasp_pcor), 0:ncol(trasp_pcor),trasp_pcor[,ncol(trasp_pcor):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
        #       axis( 2, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= rev(colnames( trasp_pcor )), las= 2 )
        #       axis( 1, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= colnames( trasp_pcor ), las= 2 )
        #       for (x in (nrow(correlation_partial)-1+0.5):0.5  )
        #         for (y in 0.5: ((ncol(correlation_partial)-1+0.5)   ))
        #           text(y,x, round(correlation_partial[ncol(correlation_partial)-x+0.5,y+0.5],2),col="blue")
        #     }) 


        #     output$savepcorAnalogHeat=renderUI({downloadButton('savepcorAnalogHeatbutton', 'Get PDF')})         
          
        #   }else{
        #     #number of BAM==2 => no pcor
        #     output$pcorAnalogHeat<-renderPlot({NULL})
        #     output$savepcorAnalogHeat=renderUI({NULL})
        #   }

        # }else{
        #   output$corAnalogHeat<-renderPlot({NULL})
        #   output$pcorAnalogHeat<-renderPlot({NULL})
        #   output$savecorAnalogHeat=renderUI({NULL})
        #   output$savepcorAnalogHeat=renderUI({NULL})
        # }
      }else{
        output$textselectedelementsAnalogHeat<-renderText({NULL})
        output$newROIfromAnalogHeat_out<-renderUI({NULL})
        #output$confirmImportROIfromAnalogHeat<-renderUI({NULL})
        output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
        output$boxplotByROIAnalogHeat<-renderPlot({NULL})
        output$profileAnalogHeat<-renderPlot({NULL})
        output$saveprofileAnalogHeat=renderUI({NULL})
        output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
  output$showprofileAnalogHeat_logOptions<-renderUI({NULL})

      output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
      output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
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
        output$showBAMsforProfilesAndBox<-renderUI({
          list(
            list(HTML("<b>Select enrichments to show:</b>"),htmlhelp("","help_enrichmentInRois_parameters_enrichments")),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("BAMsForProfilesAndBox",NULL,choices=finalBAMs)
            )
          )
        })
        
      }else{
        output$showBAMsforProfilesAndBox<-renderUI({
          list(
            HTML("<font color='red'>Some of the ROI(s) selected do not have enrichment associated.\nGo to 'ROI preparation' to associate enrichments to ROIs</font>"),
            checkboxGroupInput("BAMsForProfilesAndBox",NULL,choices=NULL)
          )
          
        })
      }
    }
  }else{
    output$showBAMsforProfilesAndBox<-renderUI({checkboxGroupInput("BAMsForProfilesAndBox",NULL,choices=NULL)})
  }
})




#observer for the number of bins? (check if nbins > width of some range)
observe({
  
  input$ROIsForProfilesAndBox
  if (!isvalid(input$BAMsForProfilesAndBox)){
    output$showbinsforProfilesAndBox<-renderUI({NULL})
    return()
  }
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
      output$showbinsforProfilesAndBox<-renderUI({
        list(
          list(HTML("<b>Number of bins:</b>"),htmlhelp("","help_enrichmentInRois_parameters_bins")),
          numericInput(inputId = 'binsProfilesAndBox',label=NULL,min = 1, max = maxtoshow,value=valuetoshow,step = 1)
        )
      })    
    }
  }else{
    output$showbinsforProfilesAndBox<-renderUI({NULL})
  }
})



#observer for button
observe({
  if (!isvalid(input$ROIsForProfilesAndBox)|!isvalid(input$BAMsForProfilesAndBox)){
    output$show_confirmUpdateProfilesAndBox<-renderUI({NULL})
    return()
  }
  output$show_confirmUpdateProfilesAndBox<-renderUI({
    actionButton("confirmUpdateProfilesAndBox", "Update plot")
  })
})






#respond to color scheme of profile in profiles&box, profile plot
observe({
  colorscheme=input$colorscheme_profileProfileAndBox
  if (isvalid(colorscheme)){
    #here is something. appear menu only if custom
    if(colorscheme=="custom"){
      bams=names(isolate(toplot$profileAndBoxes$portionlist_profile))
      lista=list()
      #bams are the names of the lines : ROI; enrichment
      if(length(bams)>0){
        output$showprofileProfileAndBox_colorlistOptions<-renderUI({
          for (i in 1:length(bams)){
            lista[[i]]=fluidRow(
                                column(12,
                                      colorSelectorInput(inputId=paste("colorCustomProfileAndBox_profile",i,sep=""),label=bams[i],choices=ColsArray,
                                                    selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                )
                        )
          }
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; background-color: #ffffff;",
            lista
          )
        })  
      }else{
        output$showprofileProfileAndBox_colorlistOptions<-renderUI({NULL})
      }
    }else{
      output$showprofileProfileAndBox_colorlistOptions<-renderUI({NULL})
    }
  }else{  
    #is null, => no menu
    output$showprofileProfileAndBox_colorlistOptions<-renderUI({NULL})
  }
})






#respond to color scheme of box in profiles&box
observe({
  colorscheme=input$colorscheme_boxProfileAndBox
  if (isvalid(colorscheme)){
    #here is something. appear menu only if custom
    if(colorscheme=="custom"){
      #if custom, cannot group colors
      output$showBoxProfileAndBox_groupcolOptions<-renderUI({NULL})
      bams=names(isolate(toplot$profileAndBoxes$portionlist_boxes))
      lista=list()
      #bams are the names of the lines : ROI; enrichment
      if(length(bams)>0){
        output$showBoxProfileAndBox_colorlistOptions<-renderUI({
          for (i in 1:length(bams)){
            lista[[i]]=fluidRow(
                                column(12,
                                      colorSelectorInput(inputId=paste("colorCustomProfileAndBox_box",i,sep=""),label=bams[i],choices=ColsArray,
                                                    selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                )
                        )
          }
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; background-color: #ffffff;",
            lista
          )
        })  
      }else{
        output$showBoxProfileAndBox_colorlistOptions<-renderUI({NULL})
      }
    }else{
      output$showBoxProfileAndBox_groupcolOptions<-renderUI({
        checkboxInput("GroupColorsProfileAndBox_box", label="Group colors",value = FALSE, width = NULL)
      })
      output$showBoxProfileAndBox_colorlistOptions<-renderUI({NULL})
    }
  }else{  
    #is null, => no menu
    output$showBoxProfileAndBox_groupcolOptions<-renderUI({NULL})
    output$showBoxProfileAndBox_colorlistOptions<-renderUI({NULL})
  }
})





#observer for log2 or not depending on the cor method (if spearman, do not show log2)
observe({
  cormethod=input$corMethodProfilesAndBox_Cor
  if (isvalid(cormethod)){
    if (cormethod=="pearson"){
      output$showCorProfileAndBox_logOptions<-renderUI({
        checkboxInput("isLog2_CorProfileAndBox", label="log2",value = FALSE, width = NULL)
      })
    }else{
      #log2 option must be NULL, because cor method is spearman
      output$showCorProfileAndBox_logOptions<-renderUI({NULL})
    }
  }else{
    #log2 option must be NULL
    output$showCorProfileAndBox_logOptions<-renderUI({NULL})
  }
})




toListenCor <- reactive({
    list(input$cor_click$y,input$cor_click$x,input$isLog2_CorProfileAndBox,input$normalization_CorProfileAndBox)
})


#observer for click the correlation/partial correlation heatmap
observeEvent(toListenCor(),{
  set.seed(123)
  if(!is.null(input$cor_click$x) & !is.null(input$cor_click$y) ){
    toplot$profileAndBoxes$x=floor(as.numeric(input$cor_click$x))+1
    toplot$profileAndBoxes$y=floor(as.numeric(input$cor_click$y))+1  
  }
  portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
  if(!is.null(portionlist_boxes) & length(toplot$profileAndBoxes$x)>0 & length(toplot$profileAndBoxes$y) >0){
    #Snorm
    if(input$normalization_CorProfileAndBox=="readdensity"){
      for (i in 1:length(portionlist_boxes)){
        portionlist_boxes[[i]]=portionlist_boxes[[i]]/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
      }
    } 
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
    if(input$isLog2_CorProfileAndBox){
      array_x=log2(array_x)
      array_y=log2(array_y)
      name_x=paste("log2",name_x)
      name_y=paste("log2",name_y)
    }
    if(input$normalization_CorProfileAndBox=="readdensity"){
      name_x=paste(name_x,"Read density (reads/bp)")
      name_y=paste(name_y,"Read density (reads/bp)")
    }else{
      name_x=paste(name_x,"Total reads")
      name_y=paste(name_y,"Total reads")
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
    list(input$pcor_click$y,input$pcor_click$x,input$isLog2_CorProfileAndBox,input$normalization_CorProfileAndBox)
})

observeEvent(toListenParcor(),{
  set.seed(123)
  if(!is.null(input$pcor_click$x) & !is.null(input$pcor_click$y) ){
    toplot$profileAndBoxes$px=floor(as.numeric(input$pcor_click$x))+1
    toplot$profileAndBoxes$py=floor(as.numeric(input$pcor_click$y))+1     
  }
  portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes

  if(!is.null(portionlist_boxes) & length(toplot$profileAndBoxes$px)>0 & length(toplot$profileAndBoxes$py) >0 &length(portionlist_boxes)>2 ){
    #Snorm
    if(input$normalization_CorProfileAndBox=="readdensity"){
      for (i in 1:length(portionlist_boxes)){
        portionlist_boxes[[i]]=portionlist_boxes[[i]]/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
      }
    } 

    mat=do.call(cbind,portionlist_boxes)
    #remove -Inf, because log2(0)=-Inf
    mat[is.infinite(mat) &mat<0 ]=0
    if(input$isLog2_CorProfileAndBox){
      mat=log2(mat)
    }
    #take a random sample (2000) of the total number of points
    if(nrow(mat)>2000){
      smpl=sort(sample(nrow(mat),2000,replace=FALSE))
      mat=mat[smpl,]
    }
    colnames(mat)=names(portionlist_boxes)
    #select idx of the matrix
    name_x=names(portionlist_boxes)[toplot$profileAndBoxes$px]
    name_y=names(portionlist_boxes)[length(portionlist_boxes)-toplot$profileAndBoxes$py+1]
    mat[is.infinite(mat) &mat<0 ]=0
    idx1=which(name_x==colnames(mat))
    idx2=which(name_y==colnames(mat))

    model=plotpcor(mat,idx=c(idx1,idx2))

    array_x=model[[1]]
    array_y=model[[2]]

    if(input$normalization_CorProfileAndBox=="readdensity"){
      name_x=paste(name_x,"Read density (reads/bp)")
      name_y=paste(name_y,"Read density (reads/bp)")
    }else{
      name_x=paste(name_x,"Total reads")
      name_y=paste(name_y,"Total reads")
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

    output$show_genelistsforDynamics<-renderUI({
      list(
        list(HTML("<b>Select gene list(s):</b>"),htmlhelp("","help_dynamicsOnGenes_parameters_genelist")),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
          checkboxGroupInput("genelistsforDynamics",label=NULL,choices = historylist)
        )
      )
    })
 
  }else{
    output$show_genelistsforDynamics<-renderUI({HTML("<font color='red'>No genelists available. To import genelists, go to 'ROIs' -> 'Get promoters, transcripts, TES coordinates of a list of genes'</font>")})  
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
        list(HTML("<b>ROIs constituting the selected gene list(s):</b>"),htmlhelp("","help_dynamicsOnGenes_parameters_ROIassociated")),
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
      output$show_BAMforDynamics<-renderUI({
        list(
          list(HTML("<b>Select enrichments to show:</b>"),htmlhelp("","help_dynamicsOnGenes_parameters_ernichments")),
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
            checkboxGroupInput(inputId="BAMforDynamics",label=NULL,choices = historylist)
          )
        )
      })

    }else{
      output$show_BAMforDynamics<-renderUI({
        list(
          HTML("<font color='red'>Some of genelist(s) selected do not have enrichment associated.\nGo to 'ROI preparation' to associate enrichments to genelists</font>"),
          checkboxGroupInput(inputId="BAMforDynamics",label=NULL,choices = NULL)
        )    
      })        
    }


  }else{
    output$show_BAMforDynamics<-renderUI({checkboxGroupInput(inputId="BAMforDynamics",label=NULL,choices = NULL)})      
    #print nothing for ROI triads
    output$showROItriadGeneList<-renderUI({NULL})
  }

})


#observer for button
observe({
  if(!isvalid(input$genelistsforDynamics)|!isvalid(input$BAMforDynamics)){
    output$show_plotDynamics<-renderUI({NULL})
    return()
  }
  output$show_plotDynamics<-renderUI({actionButton("plotDynamics","Update plot")})
})



#observer for number of bins. Depends on ROI selected/existent and enrichment selected
observe({
  input$genelistsforDynamics
  input$BAMforDynamics
  if (!isvalid(ROIvariables$listROI)|!isvalid(input$genelistsforDynamics)|!isvalid(input$BAMforDynamics)){
    output$show_binsforDynamics<-renderUI({NULL})
    return()
  }
  output$show_binsforDynamics<-renderUI({
    list(
        list(HTML("<b>Number of bins:</b>"),htmlhelp("","help_dynamicsOnGenes_parameters_nbins")),
        numericInput(inputId = 'binsforDynamics',label=NULL,min = 1, max = 300, step = 1,value=100)
    )
  })

})





#respond to color scheme in dynamics
observe({
  colorscheme=input$colorschemeDynamics
  if (isvalid(colorscheme)){
    #here is something. appear menu only if custom
    if(colorscheme=="custom"){

      bams=isolate(toplot$dynamics$totalnames)
      lista=list()
      #bams are the names of the lines : ROI; enrichment
      if(length(bams)>0){
        output$show_colorsDynamics<-renderUI({
          for (i in 1:length(bams)){
            lista[[i]]=fluidRow(
                                column(12,
                                      colorSelectorInput(inputId=paste("colorCustomDynamics",i,sep=""),label=bams[i],choices=ColsArray,
                                                    selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                )
                        )
          }
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; background-color: #ffffff;",
            lista
          )
        })  
      }else{
        output$show_colorsDynamics<-renderUI({NULL})
      }
    }else{
      output$show_colorsDynamics<-renderUI({NULL})
    }
  }else{  
    #is null, => no menu
    output$show_colorsDynamics<-renderUI({NULL})
  }
})














######################################################################################
######################################################################################
########################################################################
#GO analyses
########################################################################
######################################################################################
######################################################################################

#observer for choices (from ROI or from a custom genelist?)
observe({

  input$chooseSourceGO
  if(input$chooseSourceGO=="fromGeneList"){
    #create text area input to put gene symbols/IDs

    output$viewSelectGenesGO<-renderUI({
      list(
        HTML("<b>Paste genes here (1 gene per line):</b>"),
        textAreaInput("pastedGenesGO",NULL,value="",height=250)
      )
      
    })
  }else{
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    if(isvalid(nomi)){
      if("promoters"%in%nomi){
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
        #correction on NAs. Maybe fix the error :"NAs are not allowed in subscripted assignments"
        if(!is.null(getwdth)){
          getwdth[is.na(getwdth)]=FALSE
        }
        names(historylist)[getwdth]=paste(names(historylist)[getwdth],"(fixed size)")
        output$viewSelectGenesGO<-renderUI({
          list(
            HTML("<b>Select ROI(s):</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("selectROIGO",label=NULL,choices=historylist)
            )          
          )
        })      
      
      }else{
        #DB not present, ROI not annotated => cannot extract genes from ROIs => warning
        output$viewSelectGenesGO<-renderUI({
          HTML("<font color='red'>Need a genome assembly. Choose the appropriate database from the 'Assembly' section</font><br>")
        })
      }      
    }else{
      #no ROI present
      output$viewSelectGenesGO<-renderUI({
        HTML("<font color='red'>No ROIs present. To import ROIs, go to 'ROIs' or load a session file in 'Load session file'</font><br>")
      })      
    }

  }

})


#observer for signatures (only if some ROI or genes have been typed or selected)
observe({
  input$chooseSourceGO
  input$pastedGenesGO
  input$selectROIGO
  if (input$chooseSourceGO=="fromROI" & length(input$selectROIGO)==0){
    output$show_selectedGenesetsGO<-renderUI({checkboxGroupInput(inputId="selectedGenesetsGO",label=NULL,choices=NULL)})
    return()    
  }

  if (input$chooseSourceGO=="fromGeneList" & !isvalid(input$pastedGenesGO)){
    output$show_selectedGenesetsGO<-renderUI({checkboxGroupInput(inputId="selectedGenesetsGO",label=NULL,choices=NULL)})
    return()    
  }

  output$show_selectedGenesetsGO<-renderUI({
    list(
      list(HTML("<b>Select signature(s):</b>"),htmlhelp("","help_goAnalysis_parameters_signatures")),
      wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
        checkboxGroupInput(inputId="selectedGenesetsGO",label=NULL,choices=names(GenesetsGMT))
      )
    )
  })

})


#observer for button
observe({
  input$selectROIGO
  input$pastedGenesGO
  if (!isvalid(input$chooseSourceGO)|!isvalid(input$selectedGenesetsGO)){
    output$show_doTheGO<-renderUI({NULL})
    return()
  }
  output$show_doTheGO<-renderUI({actionButton("doTheGO","GO!")})
})




#observer for additional parameters. If from geneset, if textArea filled, ask if symbls or other IDs
#if ROI, ask if genes inside a window or the nearest
observe({
  #must react to the initial choice
  input$chooseSourceGO
  #must react to area or ROI selected
  input$pastedGenesGO
  input$selectROIGO
  if(input$chooseSourceGO=="fromROI"){
    #choice is ROI, now look at the ROIs selected. If >1, show the additional input
    if(isvalid(input$selectROIGO)){
      output$additionalparametersGO<-renderUI({
        list(
          radioButtons("chooseCriteriaROIGO",label="Choose which annotated genes" ,
                                   choiceNames=list(
                                    htmlhelp("Nearest genes","help_goAnalysis_parameters_nearestgenes"),
                                    htmlhelp("Genes inside window","help_goAnalysis_parameters_genewindow")),
                                   choiceValues=list("nearestGene","windowGene"),
                                   selected="nearestGene"
                      )
        )
      })
    }else{
      output$additionalparametersGO<-renderUI({NULL})
    }
  }else{
    #from custom list of genes. If textArea is not empty, check database. If present, 
    #show options for various IDs

    if(isvalid(input$pastedGenesGO)){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
      if("promoters"%in%nomi){
        output$additionalparametersGO<-renderUI({
          radioButtons("chooseIDgeneGO",label=list("Which kinkd of identifiers?",htmlhelp("","help_goAnalysis_parameters_kindofID")),choices=c(
                                                "Symbols"="symbol",
                                                "ENTREZ IDs"="entrez",
                                                "ENSEMBL IDs"="ensembl",         
                                                "RefSeq IDs"="refseq"
                                          ),selected="symbol")     
        })
      }else{
        #must be symbols
        output$additionalparametersGO<-renderUI({HTML("Without a genome assembly, only symbols allowed. For other identifiers, go to 'Databases' section<br>")})
      }      
    }else{
      #empty text area...
      output$additionalparametersGO<-renderUI({NULL})
    }

  }
})


#observer for showing the genomic window from ROI to take genes
observe({
  input$chooseCriteriaROIGO
  #if everything is ok, show, otherwise hide
  if(input$chooseSourceGO=="fromROI" & isvalid(input$selectROIGO) & isvalid(input$chooseCriteriaROIGO)){
    if(input$chooseCriteriaROIGO=="windowGene"){
      output$chooseWindowROIGO<-renderUI({
        numericInput(inputId = 'WindowROIGO',label="Choose genomic window:",min = 1, max = 20000, step = 50,value=1000) 
      })      
    }else{
      output$chooseWindowROIGO<-renderUI({NULL})
    }
  }else{
    output$chooseWindowROIGO<-renderUI({NULL})
  }
})





### better to keep fixde in the UI.
### so heatmap is not reactive to changes in ROI/genelist input
### but only in the way we rank

# #here, observer for clustering/ranking pf padj matrix (result of GO)
observe({
  #menu of ranking/clustering should be put only when both genelist (from ROI, from
  #single list is only ranked) and at least one geneset is selected and at least 2 ROIs selected (otherwise barplot...)
  if(input$chooseSourceGO=="fromROI" & isvalid(input$selectROIGO) & length(input$selectROIGO)>1 & isvalid(input$chooseCriteriaROIGO) & isvalid(input$selectedGenesetsGO) & length(input$selectedGenesetsGO)>0){
    #ok, here you can show the options "ranking" (according to best padj) or "clustering"
    output$chooseOrderingGO_widget<-renderUI({
      radioButtons("chooseOrderingGO",label=list("How to order results",htmlhelp("","help_goAnalysis_parameters_orderresults")),choices=c(
                                            "Ranking best padj"="ranking",
                                            "Custering"="clustering"
                                      ),selected="ranking")     
    })    
  }else{
    #else, do not show the radiobutton menu of the ranking
    output$chooseOrderingGO_widget<-renderUI({NULL})
  }

})



#here, observer for cluster type if "clustering" has been selected
observe({
  if (input$chooseSourceGO=="fromROI" & isvalid(input$selectROIGO) & length(input$selectROIGO)>1 & isvalid(input$chooseCriteriaROIGO) 
          & isvalid(input$selectedGenesetsGO) & length(input$selectedGenesetsGO)>0&isvalid(input$chooseOrderingGO)){
    #same conditions as before, PLUS the "clustering" choice
    if(input$chooseOrderingGO=="clustering"){
      output$clustertypeGO_widget<-renderUI({
        radioButtons("clustertypeGO","Cluster type",choices=c(
                                              "K-means"="kmean",
                                              "hierarchical"="hierarchical"
                                        ),selected="kmean")       
      })
    }else{
      #if we are here it means that "ranking" or no menu for ranking/clustering appeared
      output$clustertypeGO_widget<-renderUI({NULL})
    }
  }else{
    #if we are here, some elements are not valid in the chain upstream
    output$clustertypeGO_widget<-renderUI({NULL})
  }
})






#cluster number and K parameter if "kmeans" selected from the radiobutton,
#otherwise hierarchical clustering parameters
observe({
  if (input$chooseSourceGO=="fromROI" & isvalid(input$selectROIGO) & length(input$selectROIGO)>1 & isvalid(input$chooseCriteriaROIGO) 
          & isvalid(input$selectedGenesetsGO) & length(input$selectedGenesetsGO)>0&isvalid(input$chooseOrderingGO)){
 


    if(input$chooseOrderingGO=="clustering" & isvalid(input$clustertypeGO)){
      if(input$clustertypeGO=="kmean"){
        #cluster number, kstarts and kiterations. Hclust parameters are NULL
        output$clusternumbershowGO_widget<-renderUI({
          numericInput(inputId = 'clustnumGO',label="Cluster number",min = 1, max = 100, step = 1,value=4)
        })
        output$clusterKstartsGO_widget<-renderUI({
          numericInput(inputId = 'clustrandomstartsGO',label="Num. starting points",min = 1, max = 40, step = 1,value=5)
        })
        output$clusterKiterationsGO_widget<-renderUI({
          numericInput(inputId = 'clustnumiterationsGO',label="Num. iterations",min = 1, max = 40, step = 1,value=10)  
        })    
        output$clusterHDistMethodGO_widget<-renderUI({NULL})    
        output$clusterHClustMethodGO_widget<-renderUI({NULL})
      }else if(input$clustertypeGO=="hierarchical"){
        #distmethod, clustmethod. Kmeans parameters are NULL
        output$clusterHDistMethodGO_widget<-renderUI({
          selectInput("distmethodGO",label="Distance method:",c("Euclidean"="euclidean",
                                                                    "Manhattan"="manhattan",
                                                                    "Canberra"="canberra",
                                                                    "Minkowski"="minkowski"))
        })
        output$clusterHClustMethodGO_widget<-renderUI({
          selectInput("clustmethodGO",label="Clustering method:",c("Average"="average",
                                                                     "Complete"="complete",
                                                                     "Median"="median",
                                                                     "Centroid"="centroid"))
        })
        output$clusternumbershowGO_widget<-renderUI({NULL})
        output$clusterKstartsGO_widget<-renderUI({NULL})
        output$clusterKiterationsGO_widget<-renderUI({NULL})
      }

    }else{
      #no clustering parameters!
      output$clusterHDistMethodGO_widget<-renderUI({NULL})    
      output$clusterHClustMethodGO_widget<-renderUI({NULL})
      output$clusternumbershowGO_widget<-renderUI({NULL})
      output$clusterKstartsGO_widget<-renderUI({NULL})
      output$clusterKiterationsGO_widget<-renderUI({NULL})
    }
  }else{
    #some parameters are not valid....
    output$clusterHDistMethodGO_widget<-renderUI({NULL})    
    output$clusterHClustMethodGO_widget<-renderUI({NULL})
    output$clusternumbershowGO_widget<-renderUI({NULL})
    output$clusterKstartsGO_widget<-renderUI({NULL})
    output$clusterKiterationsGO_widget<-renderUI({NULL})
  }
})



#observer for min/max size (some signature must be selected)
observe({
  input$selectedGenesetsGO
  if (!isvalid(input$selectedGenesetsGO) | length(input$selectedGenesetsGO)==0){
    output$show_minmaxSizeGO<-renderUI({NULL})
    return()
  }

  if (input$chooseSourceGO=="fromROI" & length(input$selectROIGO)==0){
    output$show_minmaxSizeGO<-renderUI({NULL})
    return()    
  }
  if (input$chooseSourceGO=="fromGeneList" & length(input$pastedGenesGO)==0){
    output$show_minmaxSizeGO<-renderUI({NULL})
    return()    
  }

  output$show_minmaxSizeGO<-renderUI({
    list(
      list(HTML("<b>Min signature size:</b>"),htmlhelp("","help_goAnalysis_parameters_minsize")),
      numericInput(inputId = 'minSizeGO',label=NULL,min = 1, step = 5,value=15),
      list(HTML("<b>Max signature size:</b>"),htmlhelp("","help_goAnalysis_parameters_maxsize")),
      numericInput(inputId = 'maxSizeGO',label=NULL,min = 1, step = 5,value=500)
    )
  })

})





#observer for GO plot interactivity (click/hover mouse)
observeEvent(input$GO_click,{
  x=input$GO_click$x
  y=input$GO_click$y

  #extract from cordinates the correct genes, GO and ROI from heatmap (if heatmap)
  mat=toplot$GOvariables$completemat

  if(!is.null(mat)){
    x=ceiling(x)
    y=ceiling(y)

    if(x==0){
      x=1
    }
    if(y==0){
      y=1
    }
    if(y==11& (ncol(mat)/3) <=30){
      y=10
    }
    if(x==31 & nrow(mat)<=30){
      x=30
    }
    
    #extract name of ROIs/genelist from
    extr=grep("_Genes$",colnames(mat),value=TRUE)
    blockNames=unlist(strsplit(extr,split="_Genes"))
    termNames=rownames(mat)

    if(x<=length(termNames)){

      #retrieve all genes for a term, if wanted by the user:
      termsdt=toplot$GOvariables$termsdt
      #                   ont  gene
      # BIOCARTA_RELA_PATHWAY IKBKG
      # BIOCARTA_RELA_PATHWAY  CHUK
      # BIOCARTA_RELA_PATHWAY EP300
      # BIOCARTA_RELA_PATHWAY  RELA
      # BIOCARTA_RELA_PATHWAY   TNF
      # BIOCARTA_RELA_PATHWAY IKBKB
      #only if we are capturing something internal to the heatmap
      if(length(blockNames)==1){
        #barplot. Only one single gene list
        Term=rownames(mat)[x]
        posGene=grepl("_Genes$",colnames(mat))
        
        #if genes in the overlap:

        genes=mat[x,posGene]
        #format the genes in the overlaps to show
        if(length(genes)==1){
          if(is.na(genes)){
            genes="No genes in overlap"
          }else{
            genes=as.list(paste(strsplit(genes,split="/")[[1]],"<br>",sep=""))
            genes=paste(genes,collapse="\n")          
          }
        }else{
          genes=as.list(paste(strsplit(genes,split="/")[[1]],"<br>",sep=""))
          genes=paste(genes,collapse="\n")
        }


        #if genes in the entire term (use toplot$GOvariables$termsdt if not null to find genes)


        output$showTermClicked<-renderUI({
          list(
            HTML("<b>Term:</b><br>"),
            wellPanel(id = "logPanel",style = "overflow-x:scroll; background-color: #ffffff;",
              HTML(paste(Term,sep=""))
            )
          )
        })

        output$GenesClicked<-renderText({genes})
        output$showGenesClicked<-renderUI({
          list(
            HTML("<b>Common genes in Term:</b><br>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; max-height: 250px; background-color: #ffffff;",
              htmlOutput("GenesClicked")
            ) 
          )
        })

      }else{
        if(y<=length(blockNames)){
          #multiple gene lists: is a heatmap
          #x,y: detect which ROI and which ontological Term
          #y is reversed
          nameblockrev=rev(blockNames)
          current=nameblockrev[y]
          colname_gene=paste(current,"_Genes",sep="")
          colname_geneRatio=paste(current,"_gene_ratio",sep="")
          colname_padj=paste(current,"_padj",sep="")
          posGene=which(colnames(mat)==colname_gene)
          posgeneRatio=which(colnames(mat)==colname_geneRatio)
          pospadj=which(colnames(mat)==colname_padj)
          genes=mat[x,posGene]
          #format the genes in the overlaps to show
          if(length(genes)==1){
            if(is.na(genes)){
              genes="No genes in overlap"
            }else{
              genes=as.list(paste(strsplit(genes,split="/")[[1]],"<br>",sep=""))
              genes=paste(genes,collapse="\n")          
            }
          }else{
            genes=as.list(paste(strsplit(genes,split="/")[[1]],"<br>",sep=""))
            genes=paste(genes,collapse="\n")
          }

          geneRatio=mat[x,posgeneRatio]
          padj=mat[x,pospadj]
          Term=rownames(mat)[x]

          #output Term clicked and the ROI clicked
          output$showTermClicked<-renderUI({
            list(
              HTML("<b>Term:</b><br>"),
              wellPanel(id = "logPanel",style = "overflow-x:scroll; background-color: #ffffff;",
                HTML(paste(Term,sep=""))
              ),
              HTML("<b>ROI:</b><br>"),
              wellPanel(id = "logPanel",style = "overflow-x:scroll; background-color: #ffffff;",
                HTML(current)
              )
            )
          })

          #output window with genes
          output$GenesClicked<-renderText({genes})
          output$showGenesClicked<-renderUI({
            list(
              HTML("<b>Common genes in Term:</b><br>"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; max-height: 250px; background-color: #ffffff;",
                htmlOutput("GenesClicked")
              ) 
            )
          })        
        }else{
          output$GenesClicked<-renderText({NULL})
          output$showGenesClicked<-renderUI({NULL})    
          output$showTermClicked<-renderUI({NULL})        
        }
        

      }    
    }else{
      output$GenesClicked<-renderText({NULL})
      output$showGenesClicked<-renderUI({NULL})    
      output$showTermClicked<-renderUI({NULL})
    }

  }else{
    output$GenesClicked<-renderText({NULL})
    output$showGenesClicked<-renderUI({NULL})    
    output$showTermClicked<-renderUI({NULL})    
  }
  


})
