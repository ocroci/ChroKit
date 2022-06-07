
  ######################################################################################
  ######################################################################################
  ######################################################################################
  #CREATE ROI 
  ######################################################################################
  ######################################################################################
  ######################################################################################




#observer for all the ROI menus, if ROIlist changes
observe({

  #historylist is simply the current primary ROIs opened (ROIvariables$names)
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
  updateCheckboxGroupInput(session,inputId="selectedROIs",label=NULL,
                                    choices = historylist)
  # updateCheckboxGroupInput(session,inputId="overlapROIs",label=NULL,
  #                                   choices = historylist)  
  # updateCheckboxGroupInput(session,inputId="notoverlapROIs",label=NULL,
  #                                   choices = historylist)  

  #those 2 have been moved in uiBED, but kept here
  updateCheckboxGroupInput(session,inputId="selectedCustomROItoRemove",label=NULL,
                                  choices = historylist)  
  updateSelectInput(session,inputId="selectedCustomROItoRename",label=NULL,
                                  choices = historylist) 
  
      
  updateCheckboxGroupInput(session,"confirmviewROI", NULL,choices=historylist) 
  updateSelectInput(session,"listgetROI",NULL,choices=historylist)
  updateCheckboxGroupInput(session,"selectROItoBAMassociate", NULL,choices=historylist)
  updateSelectInput(session,inputId="selectROIforBAMrename",label="Select ROI:",
                              choices = historylist) 
  updateSelectInput(session,"selectROItoresize", "ROI to resize:",choices=historylist)
  updateSelectInput(session,"selectROItoCenterSummit", "ROI to center on summit:",choices=historylist)
  updateSelectInput(session,"selectROItoFilter", "ROI to filter:",choices=historylist)
  updateSelectInput(session,"selectROItoFilterWIDTH", "ROI to filter:",choices=historylist)
  updateSelectInput(session,"selectROItoSample", "ROI to sample:",choices=historylist)
  updateSelectInput(session,"selectROIpredefPipeline", "Select ROI for preparation:",choices=historylist)
  updateSelectInput(session,"selectROItoExtractPattern","ROI to extract pattern from:",choices=historylist)

  #?
  

  #show current updated ROI list also in coordinate files section
  output$showBEDfiles<-renderText({paste(historylist,collapse="<br>")})


})



######################################################################################
######################################################################################
######################################################################################
#create ROI from overlaps
######################################################################################
######################################################################################
######################################################################################

#react to ROI selected and put radiobutton for union/intersection
observe({
  input$selectedROIs
  if(length(input$selectedROIs)>1){
    output$RadiobuttonStartingROI<-renderUI({
      radioButtons(inputId="choiceROI",label="Criteria for building the aggregated reference ROI:" ,
                                  choices=c("Intersection"="intersection",
                                            "Union"="union"),
                                  selected="union"
                  )       
    })

  }else{
    output$RadiobuttonStartingROI<-renderUI({NULL})
  }
})


#react to ROI selected to put overlapping and not overlapping lists
observe({
  input$selectedROIs

  #historylist is simply the current primary ROIs opened (ROIvariables$names)
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  historylist=as.list(nomi)
  lens=unlist(lapply(ROIvariables$listROI,getLength))
  lens2=paste("(",lens,")",sep="")
  if(length(nomi)>0){names(historylist)=paste(nomi,lens2)} else { names(historylist)=paste(nomi,lens)}
  getwdth=lapply(ROIvariables$listROI,getWidth) 
  getwdth=unlist(lapply(getwdth, function(k) {table(!duplicated(k))["TRUE"]==1}))
  #correction on NAs. Maybe fix th error :"NAs are not allowed in subscripted assignments"
  if(!is.null(getwdth)){
    getwdth[is.na(getwdth)]=FALSE
  }
  names(historylist)[getwdth]=paste(names(historylist)[getwdth],"(fixed size)")

  if(length(input$selectedROIs)>0){
    output$columnOverlap<-renderUI({
      list( HTML("Select ROI(s):"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput(inputId="overlapROIs",label=NULL,choices = historylist)
            )
      )
    })
    output$columnNotOverlap<-renderUI({
      list(
            HTML("Select ROI(s):"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput(inputId="notoverlapROIs",label=NULL,choices = historylist)
            )
      )
    })
  }else{
    output$columnOverlap<-renderUI({HTML(paste("(select at least one reference ROI)",sep=""))})
    output$columnNotOverlap<-renderUI({HTML(paste("(select at least one reference ROI)",sep=""))})
  }
})


#react to ROI selected overlap to put radiobutton for union/intersection for contrast overlapping ROI
observe ({
  input$overlapROIs

  if(length(input$overlapROIs)>1){
    output$overlapcriteria<-renderUI({
      radioButtons("choiceoverlapROI",label="Criteria for building the aggregated contrast ROI:" ,
                              choices=c("Intersection of the selected ROIs"="stringent",
                                        #"Overlap of the contrast ROIs"="permissive",
                                        "Union of the selected ROIs"="allofthem"),
                              selected="allofthem"
      )        
    })

  }else{
    output$overlapcriteria<-renderUI({NULL})
  }
})

#react to ROI selected notOverlap to put radiobutton for union/intersection for contrast not overlapping ROI
observe ({
  input$overlapROIs

  if(length(input$notoverlapROIs)>1){
    output$notoverlapcriteria<-renderUI({
      radioButtons("choicenotoverlapROI",label="Criteria for building the aggregated contrast ROI:" ,
                          choices=c("Intersection of the selected ROIs"="intersection",
                                                  "Union of the selected ROIs"="union"),
                          selected="union"
      ) 
    })
  }else{
    output$notoverlapcriteria<-renderUI({NULL})
  }
})

######################################################################################
######################################################################################
######################################################################################
#resize UI adaptation
######################################################################################
######################################################################################
######################################################################################

#if at least one ROI selected in input, make menus with radiobuttons and values to resize
observe({
  input$chooseResizeType
  if (input$chooseResizeType=="fixedVal"){
      updateRadioButtons(session,"choosePointResize",label="Choose fixed point for resize:",choices=c(
                              "From midpoint/TSS"="fromMid",
                              "From starts"="fromStart",
                              "From ends"="fromEnd"
                            ),selected="fromMid"
      )
  }else{
        updateRadioButtons(session,"choosePointResize",label="Choose if increment/decrement width:",choices=c(
                                "Increment width from midpoint/TSS"="increment",
                                "Decrement width from midpoint/TSS"="decrement"
                              ),selected="increment"
        ) 

  }
})





observe({
  choice=input$choosePointResize
  if (length(choice)>0){
    if(choice=="fromMid" |choice=="fromStart" | choice=="fromEnd"){

      if(choice=="fromMid"){label_fix="<h4>Range center</h4>"}
      if(choice=="fromStart"){label_fix="<h4>Range start</h4>"}
      if(choice=="fromEnd"){label_fix="<h4>Range end</h4>"}
      output$showFixedMenuResize<-renderUI({
        list(
          fluidRow(
            column(width=4,HTML("<b>Upstream</b>")),
            column(width=4),
            column(width=4,HTML("<b>Downstream</b>"))
          ),
          fluidRow(
            column(width=4,
              numericInput(inputId = 'sliderUpstreamROI',label=NULL,min = 0, max = 20000, step = 50,value=1000) 
            ),
            column(width=4,
              HTML(label_fix)
            ),
            column(width=4,
              numericInput(inputId = 'sliderDownstreamROI',label=NULL,min = 0, max = 20000, step = 50,value=1000) 
            )
          )
        )
      }) 
    }else if (choice=="increment" | choice=="decrement"){
      output$showFixedMenuResize<-renderUI({
        sliderInput('chooseWidthPercResize',label="Select the % of range width for resize",min = 0, max = 100, value = 0.3,step=1)
      })
    }


  }else{
    output$showFixedMenuResize<-renderUI({NULL})
  }
  

})




######################################################################################
######################################################################################
######################################################################################
#center on summit ROI
######################################################################################
######################################################################################
######################################################################################



#observer for center on summit BAM
observe({
  
  ROIvariables$listROI
  input$selectROItoCenterSummit
  if(length(ROIvariables$listROI)>0 & length(input$selectROItoCenterSummit)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoCenterSummit,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      getbam=names(Enrichlist$rawcoverage[[pos]])
      if (!is.null(getbam)){
        updateSelectInput(session,"selectBAMtoCenterSummit", "Enrichment to use for summit:",choices=getbam)
      }else{
        updateSelectInput(session,"selectBAMtoCenterSummit", "Enrichment to use for summit:",choices=character(0))
      }
    }else{
      updateSelectInput(session,"selectBAMtoCenterSummit", "Enrichment to use for summit:",choices=character(0))
    }
  }else{
    updateSelectInput(session,"selectBAMtoCenterSummit", "Enrichment to use for summit:",choices=character(0))
  }
})


######################################################################################
######################################################################################
######################################################################################
#Filter ROI ENRICHMENT
######################################################################################
######################################################################################
######################################################################################



#observer for filter BAM
observe({
  
  ROIvariables$listROI
  input$selectROItoFilter
  if(length(ROIvariables$listROI)>0 & length(input$selectROItoFilter)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoFilter,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      getbam=names(Enrichlist$rawcoverage[[pos]])
      if (!is.null(getbam)){
        updateSelectInput(session,"selectBAMtoFilter", "Enrichment to use for filtering:",choices=getbam)
      }else{
        updateSelectInput(session,"selectBAMtoFilter", "Enrichment to use for filtering:",choices=character(0))
      }
    }else{
      updateSelectInput(session,"selectBAMtoFilter", "Enrichment to use for filtering:",choices=character(0))
    }
  }else{
    updateSelectInput(session,"selectBAMtoFilter", "Enrichment to use for filtering:",choices=character(0))
  }
})




#need observer for enrichment thresholds reactive to each other
#if ROI or BAM changes, start from 0.9 quantile
observe({
  
  input$selectBAMtoFilter
  if(length(input$selectBAMtoFilter)>0){
    #if a BAM is available, take the bam from the roi selected, and put 0.9 as quantile, and the value corresponding to that 
    #in sum of base coverage for each range
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoFilter,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      rawcovs=Enrichlist$rawcoverage[[pos]]
      decrkeys=Enrichlist$decryptkey[[pos]]
      normfacts=Enrichlist$normfactlist[[pos]]
      getbam=names(rawcovs)
      pos2=match(input$selectBAMtoFilter,getbam)
      bamselected_tmp=rawcovs[[pos2]]
      decrkeyselected=decrkeys[[pos2]]
      normfactselected=normfacts[[pos2]]
      if (!is.null(bamselected_tmp)){
        #bamselected is the baseCoverageOutput selected
        # bamselected=decryptcov( list(bamselected_tmp,decrkeyselected),chunk=length(bamselected_tmp))
        # sums=unlist(lapply(bamselected,sum))
        # sums=sums*normfactselected
        sums=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bamselected_tmp,Nbins=1,Snorm=FALSE,key=decrkeyselected,norm_factor=normfactselected)

        maxx=max(sapply(sums,max))
        quant1=round(quantile(sums,0.9))
        quant2=round(quantile(sums,0))
        fract=round(maxx/100)
        updateSliderInput(session,'quantileThreshFilter',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
        updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
        updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2))
      }else{
        updateSliderInput(session,'quantileThreshFilter',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
        updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = 0, step = 50,value=0)
        updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = 0, step = 50,value=0)
      }
    }else{
      updateSliderInput(session,'quantileThreshFilter',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
      updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = 0, step = 50,value=0)
      updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = 0, step = 50,value=0)
    }
  }else{
    updateSliderInput(session,'quantileThreshFilter',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
    updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = 0, step = 50,value=0)
    updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = 0, step = 50,value=0)
  }
})


#observer for adapt absolute if quantile changes
observe({
  
  input$quantileThreshFilter
  if(length(input$selectBAMtoFilter)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoFilter,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      rawcovs=Enrichlist$rawcoverage[[pos]]
      decrkeys=Enrichlist$decryptkey[[pos]]
      normfacts=Enrichlist$normfactlist[[pos]]

      getbam=names(rawcovs)
      pos2=match(input$selectBAMtoFilter,getbam)
      bamselected_tmp=rawcovs[[pos2]]
      decrkeyselected=decrkeys[[pos2]]
      normfactselected=normfacts[[pos2]]
      if (!is.null(bamselected_tmp)){
        #bamselected is the baseCoverageOutput selected
        # bamselected=decryptcov( list(bamselected_tmp,decrkeyselected),chunk=length(bamselected_tmp))
        # sums=unlist(lapply(bamselected,sum))
        # sums=sums*normfactselected
        sums=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bamselected_tmp,Nbins=1,Snorm=FALSE,key=decrkeyselected,norm_factor=normfactselected)

        maxx=max(sapply(sums,max))
        quant1=round(quantile(sums,input$quantileThreshFilter[1]))
        quant2=round(quantile(sums,input$quantileThreshFilter[2]))
        fract=round(maxx/100)
        updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
        updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2))
      }else{
        updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = 0, step = 50,value=0)
        updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = 0, step = 50,value=0)
      }
    }else{
      updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = 0, step = 50,value=0)
      updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = 0, step = 50,value=0)
    }
  }else{
    updateNumericInput(session,inputId = 'absoluteFilter1',label=NULL,min = 0, max = 0, step = 50,value=0)
    updateNumericInput(session,inputId = 'absoluteFilter2',label=NULL,min = 0, max = 0, step = 50,value=0)
  }
})



######################################################################################
######################################################################################
######################################################################################
#Filter ROI WIDTH
######################################################################################
######################################################################################
######################################################################################

#observer of the input$selectROItoFilterWIDTH
observe({
  
  input$selectROItoFilterWIDTH
  if(length(input$selectROItoFilterWIDTH)>0){
    #if a BAM is available, take the bam from the roi selected, and put 0.9 as quantile, and the value corresponding to that 
    #in sum of base coverage for each range
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoFilterWIDTH,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){

      #bamselected is the baseCoverageOutput selected
      widths=width(getRange(roi))
      maxx=max(widths)
      quant1=round(quantile(widths,0.9))
      quant2=round(quantile(widths,0))
      fract=round(maxx/100)
      updateSliderInput(session,'quantileThreshFilterWIDTH',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
      updateNumericInput(session,inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
      updateNumericInput(session,inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2))

    }else{
      updateSliderInput(session,'quantileThreshFilterWIDTH',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
      updateNumericInput(session,inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
      updateNumericInput(session,inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
    }
  }else{
    updateSliderInput(session,'quantileThreshFilterWIDTH',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
    updateNumericInput(session,inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
    updateNumericInput(session,inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
  }
})





#observer for adapt absolute if quantile changes
observe({
  
  input$quantileThreshFilterWIDTH
  if(length(input$selectROItoFilterWIDTH)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoFilterWIDTH,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      widths=width(getRange(roi))
      maxx=max(widths)
      quant1=round(quantile(widths,input$quantileThreshFilterWIDTH[1]))
      quant2=round(quantile(widths,input$quantileThreshFilterWIDTH[2]))
      fract=round(maxx/100)
      updateNumericInput(session,inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
      updateNumericInput(session,inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2))

    }else{
      updateNumericInput(session,inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
      updateNumericInput(session,inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
    }
  }else{
    updateNumericInput(session,inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
    updateNumericInput(session,inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
  }
})

######################################################################################
######################################################################################
# Sample ROI
######################################################################################
######################################################################################
observe({
  
  input$selectROItoSample
  if(length(input$selectROItoSample)>0){
    #if a BAM is available, take the bam from the roi selected, and put 0.9 as quantile, and the value corresponding to that 
    #in sum of base coverage for each range
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoSample,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){

      #widths=width(getRange(roi))
      maxx=length(getRange(roi))
      quant=round(quantile(1:maxx,0.3))
      fract=round(maxx/100)

      updateSliderInput(session,'quantileThreshSample',label="Fraction to keep:",min = 0, max = 1, value = 0.3,step=0.05)
      updateNumericInput(session,inputId = 'numberSample',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant))

    }else{
      updateSliderInput(session,'quantileThreshSample',label="Fraction to keep:",min = 0, max = 1, value = 0.3,step=0.05)
      updateNumericInput(session,inputId = 'numberSample',label=NULL,min = 0, max = 0, step = 50,value=0)
    }
  }else{
    updateSliderInput(session,'quantileThreshSample',label="Fraction to keep:",min = 0, max = 1, value = 0.3,step=0.05)
    updateNumericInput(session,inputId = 'numberSample',label=NULL,min = 0, max = 0, step = 50,value=0)
  }
})

#observer for adapt absolute if quantile changes
observe({
  
  input$quantileThreshSample
  if(length(input$selectROItoSample)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoSample,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      maxx=length(getRange(roi))
      quant=round(quantile(1:maxx,input$quantileThreshSample))
      fract=round(maxx/100)
      updateNumericInput(session,inputId = 'numberSample',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant))
      
    }else{
      updateNumericInput(session,inputId = 'numberSample',label=NULL,min = 0, max = 0, step = 50,value=0)
    }
  }else{
    updateNumericInput(session,inputId = 'numberSample',label=NULL,min = 0, max = 0, step = 50,value=0)
  }
})


######################################################################################
#Pattern extraction from ROI
######################################################################################

#observer to decide if motif coms from a ROI (default) or the entire genome (slow)
# observe({
#   input$choiceWherePattern
#   nomi=unlist(lapply(ROIvariables$listROI,getName))
#   historylist=as.list(nomi)
#   lens=unlist(lapply(ROIvariables$listROI,getLength))
#   lens2=paste("(",lens,")",sep="")
#   if(length(nomi)>0){
#     names(historylist)=paste(nomi,lens2)
#   }else{
#     names(historylist)=paste(nomi,lens)
#   }
#   getwdth=lapply(ROIvariables$listROI,getWidth) 
#   getwdth=unlist(lapply(getwdth, function(k) {table(!duplicated(k))["TRUE"]==1}))
#   #correction on NAs. Maybe fix th error :"NAs are not allowed in subscripted assignments"
#   if(!is.null(getwdth)){
#     getwdth[is.na(getwdth)]=FALSE
#   }
#   names(historylist)[getwdth]=paste(names(historylist)[getwdth],"(fixed size)")
    

# })


#find current BSgenome in the installed packages. If not found, put a button for installation.
#look DATABASEvariables$currentASSEMBLY. If not null, search BSgenome in the package. If 
#present, import, put in DBvariables$BSgenomeDB and do nothing (maybe print what we are using), 
#otherwise, downloadbutton 
observe({
  #ROIvariables$listROI
  DATABASEvariables$currentASSEMBLY
  #use an old name for retro-compatibility
  DATABASEvariables$currentORG
  if(length(DATABASEvariables$currentASSEMBLY)>0){
    #calculate BSgenome string
    #BSgenome.Xyyyyy.UCSC.(genomeass)
    asm=DATABASEvariables$currentASSEMBLY
    avail_spl=strsplit(all_avail_assemblies,split="\\.")
    org=sapply(avail_spl,"[[",2)
    asms=sapply(avail_spl,"[[",4)
    pos=match(asm,asms)
    BSstring=paste("BSgenome.",org[pos],".UCSC.",asm,sep="")
    x=rownames(installed.packages())
    pos_pkg=match(BSstring,x)
    if(!is.na(pos_pkg)){
      #BSgenome found. Import library and do nothing
      library(BSstring,character.only=TRUE)
      output$showWarningBSgenome<-renderUI({NULL})
      output$showWarningBSgenome2<-renderUI({NULL})
    }else{
      #button for download the package
      output$showWarningBSgenome<-renderUI({
        list(
        HTML(paste("<font color='red'>Need ",BSstring," package to be installed. Install it now?</font><br>",sep="")),
        actionButton("DownloadBSgenome","Download")
        )
      })

      output$showWarningBSgenome2<-renderUI({
        list(
        HTML(paste("<font color='red'>Need ",BSstring," package to be installed. Install it now?</font><br>",sep="")),
        actionButton("DownloadBSgenome2","Download")
        )
      })
    }


  }else{
    #warning message: I need database of a genome assembly
    output$showWarningBSgenome<-renderUI({HTML("<font color='red'>Warning: choose a genome assembly from 'Databases' section</font>")})
    output$showWarningBSgenome2<-renderUI({HTML("<font color='red'>Warning: choose a genome assembly from 'Databases' section</font>")})
  }

})



#observer for options of extracting pattern from ROI:
observe({
  input$choiceWherePattern
  input$selectROItoExtractPattern
  if(isvalid(input$selectROItoExtractPattern)>0){

    #find strand. If * => only print text, if + and - , give the choice
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    if(isvalid(nomi)){
      pos=match(input$selectROItoExtractPattern,nomi)


      if(!is.null(ROIvariables$listROI[[pos]])){
        range_roi=getRange(ROIvariables$listROI[[pos]])
        strand_roi=unique(as.character(strand(range_roi)))

        if("*" %in% strand_roi){
          #pattern search on both strands, because strand not defined
          output$showStrandOptsPattern<-renderUI({HTML("Strands not defined: pattern search on both strands<br>")})
        }else{
          #stranded (+ and -) => offer the choice
          output$showStrandOptsPattern<-renderUI({
            list(
            HTML("<b>Strand selection:</b><br>"),
            radioButtons("strandOptsPattern",label=NULL ,
                                     choices=c("Both strands"="bothStrands",
                                               "Strand-specific"="strandSpecific"),
                                     selected="strandSpecific"
            )
            )
          })
        }          
      }else{
        output$showStrandOptsPattern<-renderUI({NULL})
      }
      
    }else{
      output$showStrandOptsPattern<-renderUI({NULL})
    }

      
  }else{
    #NULL
    output$showStrandOptsPattern<-renderUI({NULL})
  }    
  
})




#missing adapting quantile sliderINput for filtering, following absoluteFilter modification

  ######################################################################################
  ######################################################################################
  ######################################################################################
  #VIEW ROI statistics
  ######################################################################################
  ######################################################################################
  ######################################################################################



  #DENSITY plot for ROI width. Respond to changes in selection
  observeEvent(input$updatechoiceROI,{

      #plot density plot of width of peak
      output$viewROIwidth<-renderPlot ( {
        par(mar=c(4,4,4,4))
        if (!is.null(ROIvariables$listfordensity) & length(ROIvariables$listfordensity)>0){
          coords=lapply(ROIvariables$listfordensity,function(k) {return(list(k$x,k$y))})
          coordsx=unlist(lapply(coords,"[[",1))
          coordsy=unlist(lapply(coords,"[[",2))
          xmax=max(coordsx)
          xmin=min(coordsx)
          ymax=max(coordsy)
          ymin=min(coordsy)
          leg=unlist(ROIvariables$listselectednames)
          plot(1, type="n", xlab="log2 width", ylab="density", xlim=c(xmin, xmax), ylim=c(ymin, ymax),main="frequency plot width")
          for (i in 1:length(ROIvariables$listfordensity)){
            lines(ROIvariables$listfordensity[[i]],col=ROIvariables$colorsfordensity[[i]],lwd=3)
          }
          legend("topright",leg,col=unlist(ROIvariables$colorsfordensity),lty=1,lwd=3,bty = "n")

        }else{
          plot(1, type="n", xlab="log2 width", ylab="density",main="no ROI selected")
        }
      })
  })


  observe({
    
        input$confirmviewROI
        output$viewROIstat<-renderText(
          if(length(input$confirmviewROI)==1){
            nomi=unlist(lapply(ROIvariables$listROI,getName))
            pos=match(input$confirmviewROI,nomi)
            roi=ROIvariables$listROI[[pos]]
            if(!is.null(roi)){
              selectedrange=getRange(roi)
            
              wdth=width(selectedrange)
              x=quantile(wdth,input$quantileROIwidth)
              y=length(selectedrange[wdth>x])
              paste("Width at that quantile: <b>",round(x,0),"</b><br>Number peaks with greater width: <b>",y,"</b>",sep="")
            }else{
              paste("You have selected a non existent ROI...")
            }

          }else{
            paste("You have to select ONE ROI...")
          }
        )

  })



  #BARPLOT for ROI number
  observeEvent(ROIvariables$changed,{
      
      if (length(ROIvariables$listselected)>0 & length(ROIvariables$listROI)>0 ){

        output$viewpeaksnumberROI<-renderPlot({
          par(mar=c(15,5,1,3))
          peaknumberlist=unlist(lapply(ROIvariables$listselected,length))
          names(peaknumberlist)=ROIvariables$listselectednames
          barplot (peaknumberlist,las=2,col=unlist(ROIvariables$colorsfordensity))
          mtext(side=2, text="Peaks number", line=4)
          
        })
      }else{
        output$viewpeaksnumberROI<-renderPlot({NULL})
      }
  })



  observe({
    
        input$confirmviewROI
        output$viewROIsource<-renderText(
          if(length(input$confirmviewROI)==1){
           
            nomi=unlist(lapply(ROIvariables$listROI,getName))
            pos=match(input$confirmviewROI,nomi)
            roi=ROIvariables$listROI[[pos]]
            logsource=getSource(roi)
            paste(unlist(logsource),collapse="<br>")
          }else{
            paste("No ROI selected or multiple selection...")
          }
        )

  })


######################################################################################
######################################################################################
######################################################################################
#GET ROI tab (update checkbox for the features of the ROI selected)
#and clean the data table if other ROI is selected 
######################################################################################
######################################################################################
######################################################################################
observe({
  
  input$listgetROI
  ROIcompleteList=ROIvariables$listROI
  #find and show all possible fatures of the ROI in CheckboxGroupInput, interactively
  #if a ROI is selected

  if(  !(!is.null(input$listgetROI)&length(input$listgetROI)>0 & input$listgetROI!="" & length(ROIcompleteList)>0 )   )  {
    output$showROIoptionsToGET<-renderUI({
      paste("No ROI available...")
    }) 
    output$showROIoptionsToGET<-renderUI({
      paste("No ROI available...")
    }) 
    output$previewROItodownload<-renderUI({
      paste("Select a ROI...")
    })
    output$previewROItodownloadbutton<-renderUI({NULL})
    tosave$datatableROI=NULL

  }else{
    nomi=unlist(lapply(ROIcompleteList,getName))
    pos=match(input$listgetROI,nomi)
    roi=ROIcompleteList[[pos]]
    if(!is.null(roi)){
      listGUI=list()
      emd=as.data.frame(elementMetadata(getRange(roi)))
      if (input$choosegetROItype=="eachRange"){
        #chr, start, end and strand are by definition, present
        toadd_range=c("chr; start; end; strand"="ranges")
        listGUI=c(listGUI,list(checkboxGroupInput("ROIoptionsToViewRANGE","View Ranges:",choices=toadd_range)))
        if(ncol(emd)!=0){
          #add metadata cols (for example, the annotation) in the options to provide
          toadd_emd=colnames(emd)
          names(toadd_emd)=toadd_emd
          listGUI=c(listGUI,list(checkboxGroupInput("ROIoptionsToViewMETADATA","Select annotations to show:",choices=toadd_emd)))
        }else{
          listGUI=c(listGUI,list(renderText("No annotation found for the selected ROI.")))
        }
        ####newenrichimplementation####
        rawvals=Enrichlist$rawcoverage[[pos]]
        ###############################

        bams=names(rawvals)
        if(length(bams)>0){
          #add bam file enrichment to the options to provide
          toadd_bams=bams
          names(toadd_bams)=paste(bams,"enrichment")
          listGUI=c(listGUI, list(
              HTML("<b>Enrichments:</b>"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; background-color: #ffffff;",
                checkboxGroupInput("ROIoptionsToViewENRICHMENTS",NULL,choices=toadd_bams)
              )))

        }else{
          listGUI=c(listGUI,list(renderText("No enrichments associated to the selected ROI.")) )
        }
        listGUI=c(listGUI,list(HTML("<br>"),actionButton("showdataframeROI", "Preview ROI")))

      }else if (input$choosegetROItype=="genesWindow"){

        if(ncol(emd)!=0){
          toadd_emd=colnames(emd)
          #show the window and the button for the window annotation
          if("gene_id"%in%toadd_emd & "symbol"%in%toadd_emd & "ensembl_id" %in%toadd_emd & "refSeq_id"%in%toadd_emd ){
            listGUI=list(HTML('<hr size="10">'),
              HTML("Get annotated genes within a genomic window:"),
              #show all IDs (excluding distancefromTSS)
              radioButtons("IDorsymbolwindow",label=NULL ,
                                       choices=c("gene_id"="gene_id",
                                                 "symbol"="symbol",
                                                 "ensembl_id"="ensembl_id",
                                                 "refSeq_id"="refSeq_id"),
                                       selected="gene_id"
              ),
              HTML("<br>"),
              HTML("Select the genomic window (bp):"),
              numericInput(inputId = 'windowAnnotateGenes',label=NULL,min = 100, max = 200000, step = 1000,value=20000),
              actionButton("showgenelistWindowROI", "Preview genes in window")     
            )   
          }else{
            listGUI<-list(paste("No annotation found for the selected ROI. Please import a genome assembly first."))
            output$previewROItodownload<-renderUI({
              paste("")
            })
            output$previewROItodownloadbutton<-renderUI({NULL})
            tosave$datatableROI=NULL           
          }
        }else{
          listGUI<-list(paste("No annotation found for the selected ROI. Please import a genome assembly first."))
          output$previewROItodownload<-renderUI({
            paste("")
          })
          output$previewROItodownloadbutton<-renderUI({NULL})
          tosave$datatableROI=NULL                
        }
        

      }else{
        #here edit/download ROI notes
        #get value of source of the selected ROI:
        text=getSource(roi)
        windowNotesToShow<-list(textAreaInput(inputId="ROInotes",label="View and edit the notes of the selected ROI:",value=text,cols=60,rows=10),
                      actionButton("saveNotesROI", "save notes"),
                      downloadButton('downloadNotesROI', 'Download notes'))
        output$previewROItodownload<-renderUI({windowNotesToShow})
        output$previewROItodownloadbutton<-renderUI({NULL})
        listGUI<-list(paste(""))
      }  







      output$showROIoptionsToGET<-renderUI({listGUI})





    }else{
      output$showROIoptionsToGET<-renderUI({
        paste("No ROI available...")
      }) 
      output$previewROItodownload<-renderUI({
        paste("Select a ROI...")
      })
      output$previewROItodownloadbutton<-renderUI({NULL})
      tosave$datatableROI=NULL      
    }    
  
  }


})


######################################################################################
######################################################################################
######################################################################################
#ASSOCIATE BAM tab update
######################################################################################
######################################################################################
######################################################################################





#reactive through the selectd ROI to cover with enrichments, showing
#the union of enrichments missing and the union of available enrichments
observe({
  
  input$selectROItoBAMassociate
  if(isvalid(input$selectROItoBAMassociate) & length(ROIvariables$listROI)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoBAMassociate,nomi)
    rois=ROIvariables$listROI[pos] 
    allBAMavailable=c()
    allBAMmissing=c()
    bamblock=Enrichlist$rawcoverage[pos]
    #in the bampresent=names(getBAMlist(rois[[i]])) can raise an error, not always
    #it tells that unable to find an inherited method for function ‘getBAMlist’ for signature ‘"NULL"’
    for (i in 1:length(rois)){

      bampresent=names(bamblock[[i]])

      allbams=names(BAMvariables$listBAM)
      remainingbams=setdiff(allbams,bampresent)
      allBAMmissing=c(allBAMmissing,remainingbams)
      allBAMavailable=c(allBAMavailable,bampresent)
    }
    #define the union of them and show in the menu
    allBAMmissing=unique(allBAMmissing)
    allBAMavailable=unique(allBAMavailable)
    if(is.null(allBAMmissing)){
      allBAMmissing=character(0)
    }
    if(is.null(allBAMavailable)){
      allBAMavailable=character(0)
    }

    updateCheckboxGroupInput(session,"selectBAMtoassociate",NULL,choices=allBAMmissing)
    updateCheckboxGroupInput(session,"selectBAMtoDeassociate",NULL,choices=allBAMavailable)

  }else{
    #null input fields
    updateCheckboxGroupInput(session,"selectBAMtoDeassociate",NULL,choices=character(0))
    updateCheckboxGroupInput(session,"selectBAMtoassociate",NULL,choices=character(0))
  }


})


#observer for the menu to choose how to normalize (react to BAM selected for association)
observe({
  input$selectBAMtoassociate

  listbams=isolate(BAMvariables$listBAM)
  pos2=match(input$selectBAMtoassociate,names(listbams))
  paths=listbams[pos2]

  #if selected are all wig files, do not show normalization menu
  if (all(substring(paths,nchar(paths)-2,nchar(paths))==".bw" | tolower(substring(paths,nchar(paths)-6,nchar(paths)))==".bigwig")){
    output$radioForNorm<-renderUI({NULL})
  }else{
    if(length(input$selectBAMtoassociate)>0 & length(BAMvariables$listBAM)<3){
      output$radioForNorm<-renderUI({
        radioButtons("selectMethodForNorm",label="Normalization method:" ,
                                     choices=c("no normalization"="nonorm",
                                               "library-size"="librarysize",
                                               "custom normalizer"="customnorm"),
                                     selected="librarysize"
                        )
      })
    #spike in normalization only if we have more than 3 bam file (current and input+spikeInput)
    }else if (length(input$selectBAMtoassociate)>0 & length(BAMvariables$listBAM)>=3){
      output$radioForNorm<-renderUI({
        radioButtons("selectMethodForNorm",label="Normalization method:" ,
                                     choices=c("no normalization"="nonorm",
                                               "library-size"="librarysize",
                                               "custom normalizer"="customnorm",
                                               "spike-in"="spikein"),
                                     selected="librarysize"
                        )
      })
    }else{
      #null menu
      output$radioForNorm<-renderUI({NULL})
    }    
  }
  
})

#observer for normalizer selection(react to radiobutton, list only BAM available (not wig if possible))
observe({
  input$selectMethodForNorm
  if(length(input$selectMethodForNorm)>0){
    
    if(length(input$selectBAMtoassociate)>0){
      #find, among all enrichment files, the BAM files
      allbams=BAMvariables$listBAM
      #find .bam at the end of file path
      pos=grepl("\\.bam$",allbams)
      namestoshow=names(allbams)[pos]


      if(input$selectMethodForNorm=="customnorm"){
        output$menuForNorm<-renderUI({
          selectInput("customNormalizer",choices=namestoshow,label="Choose the normalizer reads:",NULL)
        })              
      }else if (input$selectMethodForNorm=="spikein"){
        #3 menus with bam file
        output$menuForNorm<-renderUI({ list (
          selectInput("customNormalizer",choices=namestoshow,label="Choose the spike-in reads:",NULL),
          selectInput("ctrlNormalizer",choices=namestoshow,label="Choose the control reads:",NULL),
          selectInput("ctrlSpikeinNormalizer",choices=namestoshow,label="Choose the control spike-in reads:",NULL)
          )

        })
        
      }else{
        output$menuForNorm<-renderUI({NULL})
      }
      
    }else{
      output$menuForNorm<-renderUI({NULL})
    }
  }else{
    output$menuForNorm<-renderUI({NULL})
  }
})




  ######################################################################################
  ######################################################################################
  ######################################################################################
  #Rename/reorder BAM tab 
  ######################################################################################
  ######################################################################################
  ######################################################################################


  #update selectInput of BAMs associated to selected ROI




#based on the radiobutton, rename or reorder menu
observe({
  input$selectRenameOrReorder
  if(!is.null(ROIvariables$listROI) & length(ROIvariables$listROI)>=1){
    if(input$selectRenameOrReorder=="rename"){
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      pos=match(input$selectROIforBAMrename,nomi)
      roi=ROIvariables$listROI[[pos]]

      if (!is.null(roi) & length(BAMvariables$listBAM)>0){

        getbam=names(Enrichlist$rawcoverage[[pos]])
        if (length(getbam)!=0){
          output$RenameReorder<-renderUI({
            list(
              HTML("<h3>Rename enrichments</h3><br>"),
              selectInput("selectedBAMtoRename",label="Select enrichment to rename:",choices=getbam),
              textInput("newBAMname","Select new name:",placeholder="type new enrichment name here",value=""),
              actionButton("renameBAM", "Rename")
            )
          })
        }else{
          output$RenameReorder<-renderUI({ HTML("No enrichments associated to the selected ROI...")})
        }

      }else{
        output$RenameReorder<-renderUI({ HTML("No ROI...")})
      }
    }else{
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      pos=match(input$selectROIforBAMrename,nomi)
      roi=ROIvariables$listROI[[pos]]
      if (!is.null(roi) ){
        getbam=names(Enrichlist$rawcoverage[[pos]])
        if(length(getbam)>1){
          choicelist=as.list(1:length(getbam))
          names(choicelist)=as.character(1:length(getbam))
          lista=list()
          for (i in 1:length(getbam)){
            lista[[i]]=fluidRow(column(3,
                                      selectInput(inputId = paste("reorderoptionBAM",i,sep=""), label = NULL, 
                                   choices = choicelist,selected=i)),
                      column(4,HTML(getbam[i])))
          }

          output$RenameReorder<-renderUI({

            list(
              HTML("<h3>Reorder enrichments</h3><br>"),
              wellPanel(id = "logPanelBAM",style = "overflow-y:scroll; max-height: 400px",
                fluidPage(
                  lista
                )
              ),
              actionButton("reorderBAM","Reorder!")
            )
          })





        }else{
          output$RenameReorder<-renderUI({ HTML("At least 2 enrichments must be associated to the selected ROI...")})
        }
      }else{
        output$RenameReorder<-renderUI({ HTML("No ROI...")})
      }

    }    
  }else{
    output$RenameReorder<-renderUI({NULL})
  }

})



  #update text field for renaming the BAM to blank for any change in ROIlist
  # observe({
    
  #   nomi=unlist(lapply(ROIvariables$listROI,getName))
  #   updateTextInput(session,inputId="newBAMname",label="New enrichment name:",value="")    
  # })



  # #reorder BAM for the ROI. Must respond to input$selectROIforBAMrename, the same!
  # observeEvent(ROIvariables$listROI,{
    
  #   output$dinamicBAM<-renderUI({
  #     if (!is.null(ROIvariables$listROI) & length(ROIvariables$listROI)>=1){
        

  #       if (!is.null(roi) ){
  #         getbam=names(Enrichlist$rawcoverage[[pos]])
  #         if (length(getbam)>0){
  #           choicelist=as.list(1:length(getbam))
  #           names(choicelist)=as.character(1:length(getbam))
  #           for (i in 1:length(getbam)){
  #             lista[[i]]=fluidRow(column(3,
  #                                     selectInput(inputId = paste("reorderoptionBAM",i,sep=""), label = NULL, 
  #                                  choices = choicelist,selected=i)),
  #                     column(4,HTML(getbam[i])))
  #           }
  #           return(lista)  
  #         }else{
  #           return(HTML("No enrichment..."))
  #         }

  #       }else{
  #         return(HTML("No enrichment..."))
  #       }     
  #     }else{
  #       return(HTML("No enrichment..."))
  #     }

  #   })   
  # })




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
            HTML("<b>Select ROI:</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("selectROIGO",label=NULL,choices=historylist)
            )          
          )
        })      
      
      }else{
        #DB not present, ROI not annotated => cannot extract genes from ROIs => warning
        output$viewSelectGenesGO<-renderUI({
          HTML("<font color='red'>Need a genome assembly. Choose the appropriate database from the 'Databases' section</font><br>")
        })
      }      
    }else{
      #no ROI present
      output$viewSelectGenesGO<-renderUI({
        HTML("<font color='red'>No ROI present</font><br>")
      })      
    }

  }

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
                                   choices=c("Nearest genes"="nearestGene",
                                             "Genes inside window"="windowGene"),
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
          radioButtons("chooseIDgeneGO","Which kinkd of identifiers?",choices=c(
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
      radioButtons("chooseOrderingGO","How to order results",choices=c(
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



######################################################################################
######################################################################################
########################################################################
#Predefined Pipeline for ROI preparation for analogic heatmap
########################################################################
######################################################################################
######################################################################################




#react to slider of fraction of ranges to keep in predefined pipeline for ROI preparation
observe({
  input$quantileThreshPredefPipeline
  if(length(input$selectROIpredefPipeline)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROIpredefPipeline,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      maxx=length(getRange(roi))
      quant=round(quantile(1:maxx,input$quantileThreshPredefPipeline))
      fract=round(maxx/100)
        
      if(unname(quant)>50000){
        warn="WARNING: genomic intervals will be > 50000"
      }else{
        warn=""
      }
      msg= paste(paste("Genomic intervals to keep: ",unname(quant),"/",maxx,sep=""),
                  "<br>","<font color='red'>",warn,"</font>",sep="")
      #print(msg)         
      output$textNumRangesPredefPipeline<-renderText({msg})
 
    }else{
      output$textNumRangesPredefPipeline<-renderText({NULL})
    }
  }else{
    output$textNumRangesPredefPipeline<-renderText({NULL})
  }
})







#observer for BAM association in ROI heatmap preparation
observe({
  ROIvariables$listROI
  input$selectROIpredefPipeline


    if(length(ROIvariables$listROI)>0 & length(input$selectROIpredefPipeline)>0){
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      pos=match(input$selectROIpredefPipeline,nomi)
      roi=ROIvariables$listROI[[pos]]
      if (!is.null(roi)){
        #look into opened enrichment files
        allbams=names(BAMvariables$listBAM)
        if(!is.null(allbams) & length(allbams)>0){
          output$menuEnrichPredefPipeline<-renderUI({
            list(
              HTML("<b>Select enrichment to associate to the new ROI:</b>"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 400px; background-color: #ffffff;",
                    checkboxGroupInput(inputId="enrichAllPredefPipeline",label=NULL,choices=allbams)
              )
            )
          }) 


          if (input$choiceSummitPredefPipeline=="Yes"){
            output$menuSummitPredefPipeline<-renderUI({
              list(
                HTML("<b>Select enrichment to use for summit:</b>"),
                selectInput("selectBAMsummitPredefPipeline",label=NULL,choices=allbams)
              )
            })  
          }else{
            output$menuSummitPredefPipeline<-renderUI({NULL})
          }
 



        }else{
          output$menuEnrichPredefPipeline<-renderUI({
            paste("No enrichment files opened")
          })
          output$menuSummitPredefPipeline<-renderUI({
            paste("No enrichment files opened")
          })
        }
      }else{
        output$menuEnrichPredefPipeline<-renderUI({NULL})
        output$menuSummitPredefPipeline<-renderUI({NULL})
      }
    }else{
      output$menuEnrichPredefPipeline<-renderUI({NULL})
      output$menuSummitPredefPipeline<-renderUI({NULL})
    }


})






