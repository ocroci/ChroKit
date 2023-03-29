#observer to show ROI menus based on existing ROIs
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
  if (length(historylist)>0){
    output$show_selectROIpredefPipeline<-renderUI({
        selectInput("selectROIpredefPipeline", label="Select ROI to prepare for heatmaps:",choices=historylist)
    })
    output$show_selectedROIs<-renderUI({
      wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
        checkboxGroupInput(inputId="selectedROIs",label=NULL,choices = historylist)
      )
    })

    output$show_selectROItoBAMassociate<-renderUI({
      list(
        HTML("<b>Choose ROI(s):</b><br>"),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
          checkboxGroupInput("selectROItoBAMassociate",label=NULL,choices = historylist)
        )
      )
    })
    output$show_selectROIforBAMrename<-renderUI({
      list(
        selectInput("selectROIforBAMrename",label="Select ROI for which rename enrichments:",choices=historylist)
      )
    })
    output$show_selectROIforBAMreorder<-renderUI({
      list(
        selectInput("selectROIforBAMreorder",label="Select ROI for which reorder enrichments:",choices=historylist)
      )
    })

    output$show_selectROIultraeasy<-renderUI({
      list(
        selectInput("selectROIultraeasy",label="Select ROI to prepare for analyses:",choices=historylist)
      )
    })


    output$show_selectROItoresize<-renderUI({
      list(
        selectInput("selectROItoresize",label="Select ROI to resize:",choices =historylist),
        radioButtons("chooseResizeType",label=NULL,choices=c(
                                "Resize using fixed value"="fixedVal",
                              "Resize using percentage of range width"="percentageVal"
                            ),selected="fixedVal")
      )
    })

    output$show_selectROItoFilterWIDTH<-renderUI({
      list(
        selectInput("selectROItoFilterWIDTH",label="Select ROI to filter:",choices =historylist)
      )
    })


    output$show_selectROItoSample<-renderUI({
      selectInput("selectROItoSample",label="Select ROI to sample:",choices =historylist)
    })

    output$show_selectROItoExtractPattern<-renderUI({
      selectInput("selectROItoExtractPattern","ROI to extract pattern from:",choices=historylist)
    })

    output$show_selectROItoCenterSummit<-renderUI({
      selectInput("selectROItoCenterSummit",label="Select ROI to center on summit:",choices=historylist)
    })
    
    output$show_selectROItoFilter<-renderUI({
      selectInput("selectROItoFilter",label="Select ROI to filter:",choices=historylist)
    })


    output$show_selectROItoRename<-renderUI({
      list(
        selectInput(inputId="selectedCustomROItoRename",label="ROI to rename:",choices = historylist),
        textInput("newfilenameROI","New ROI name:",placeholder="type new ROI name here",value=""),
        actionButton("renameROI", "Rename")
      )
    })


    output$show_selectROItoEditNotes<-renderUI({
      selectInput("selectROItoEditNotes",label="Select ROI for which notes should be edited:",choices=historylist)
    })

  }else{
    output$show_selectedROIs<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROIpredefPipeline<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROIultraeasy<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoBAMassociate<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROIforBAMrename<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROIforBAMreorder<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoresize<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoFilterWIDTH<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoSample<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoExtractPattern<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoCenterSummit<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoFilter<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoRename<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    output$show_selectROItoEditNotes<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
  }

})


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#observers for preparing ROI basic (ultraeasy)
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#observer for enrichment menu
observe({
  if (!isvalid(input$selectROIultraeasy)){
    output$show_selectBAMultraeasy<-renderUI({checkboxGroupInput(inputId="selectBAMultraeasy",label=NULL,choices=NULL)})
    return()
  }
  allbams=names(BAMvariables$listBAM)   
  if(!isvalid(allbams)){
    output$show_selectBAMultraeasy<-renderUI({
      list(
        HTML("<font color='red'>No enrichment files found. Need to import enrichment file(s)</font>"),
        checkboxGroupInput(inputId="selectBAMultraeasy",label=NULL,choices=NULL)
      ) 
    })
    return()    
  }

  nomi=unlist(lapply(ROIvariables$listROI,getName))
  pos=match(input$selectROIultraeasy,nomi)
  roi=ROIvariables$listROI[[pos]] 
  allBAMmissing=c()

  bamblock=Enrichlist$rawcoverage[[pos]]
  #in the bampresent=names(getBAMlist(rois[[i]])) can raise an error, not always
  #it tells that unable to find an inherited method for function ‘getBAMlist’ for signature ‘"NULL"’
  bampresent=names(bamblock)
  remainingbams=setdiff(allbams,bampresent)
  allBAMmissing=c(allBAMmissing,remainingbams)

  if(is.null(allBAMmissing)){
    allBAMmissing=character(0)
  }

  if (!isvalid(allBAMmissing)){
    output$show_selectBAMultraeasy<-renderUI({
      list(
        HTML("<font color='red'>All available enrichment files have been associated to selected ROI</font>"),
        checkboxGroupInput(inputId="selectBAMultraeasy",label=NULL,choices=NULL)
      ) 
    })
    return()
  }


  output$show_selectBAMultraeasy<-renderUI({
    list(
      HTML("<br>"),
      list(HTML("<b>Select enrichment(s) to associate to the selected ROI:</b>"),htmlhelp("","help_ultraeasy_enrichmentAssoc")),
      wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 400px; background-color: #ffffff;",
            checkboxGroupInput(inputId="selectBAMultraeasy",label=NULL,choices=allBAMmissing)
      )
    )
  }) 




})



#observer for button ultraeasy ROI prep
observe({
  if (!isvalid(input$selectROIultraeasy)|!isvalid(input$selectBAMultraeasy)){
    output$show_Buttonultraeasy<-renderUI({NULL})
    return()
  }
  output$show_Buttonultraeasy<-renderUI({actionButton("Buttonultraeasy","Prepare ROI")})
})




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#observers for prepare for heatmap (EASY)
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



#observer to show or not the quantile bar (%genomic ranges to show) and radio for summit based on selection of a ROI
observe({
  input$selectROIpredefPipeline
  if (isvalid(input$selectROIpredefPipeline)){
    output$show_quantileThreshPredefPipeline<-renderUI({
      list(
        HTML("<br>"),
        #subsample %
        sliderInput('quantileThreshPredefPipeline',label=list("% genomic ranges to keep:",htmlhelp("","help_prepheat_subsample")),min = 0, max = 1, value = 0.3,step=0.05)
      )
    })
    output$show_radioButtonPredefSummit<-renderUI({
      list(
        HTML("<br>"),
        #subsample %
        radioButtons("choiceSummitPredefPipeline",label=list("Do you want to center on summit?",htmlhelp("","help_prepheat_summit")) ,
                                      choices=c("Yes"="Yes",
                                                "No"="No"),
                                      selected="No"
            )
      )
    })

    output$show_updownstreamPredefPipeline<-renderUI({
      list(
        HTML("<br>"),
        list(HTML("<b>Resize:</b>"),htmlhelp("","help_prepheat_resize"),HTML("<br>")),
        #set intervals from midpoint/summit
        fluidRow(
          column(width=5,HTML("<b>Upstream (bp from center)</b>")),

          column(width=5,HTML("<b>Downstream (bp from center)</b>")),
          column(width=2)
        ),
        fluidRow(
          column(width=5,
            numericInput(inputId = 'sliderUpstreamPredefPipeline',label=NULL,min = 0, max = 20000, step = 50,value=2000) 
          ),

          column(width=5,
            numericInput(inputId = 'sliderDownstreamPredefPipeline',label=NULL,min = 0, max = 20000, step = 50,value=2000) 
          ),
          column(width=2)
        )
      )
    })


    output$show_ROInamePredefPipeline<-renderUI({
      list(
        HTML("<br>"),
        HTML("<b>Name of new ROI:</b>"),
        fluidRow(
          column(width=4,
            textInput("ROInamePredefPipeline",label=NULL,placeholder="type new ROI name here")
          ),
          column(width=4,
            actionButton("PrepareROIpredefPipeline","Create ROI")
          ),
          column(width=4)
        )

      )
    })




  }else{
    output$show_quantileThreshPredefPipeline<-renderUI({NULL})
    output$show_radioButtonPredefSummit<-renderUI({NULL})
    output$show_updownstreamPredefPipeline<-renderUI({NULL})
  }

})







#observer for showing the number of ranges selected with the bar
observe({
  input$quantileThreshPredefPipeline
  if(isvalid(input$quantileThreshPredefPipeline)&isvalid(input$selectROIpredefPipeline)){
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
      output$textNumRangesPredefPipeline<-renderUI({NULL})
    }
  }else{
    output$textNumRangesPredefPipeline<-renderUI({NULL})
  }
})





#observer for BAM association in ROI heatmap preparation
observe({
  ROIvariables$listROI
  input$selectROIpredefPipeline

    if(length(ROIvariables$listROI)>0 & isvalid(input$selectROIpredefPipeline)>0){
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      pos=match(input$selectROIpredefPipeline,nomi)
      roi=ROIvariables$listROI[[pos]]
      if (!is.null(roi)){
        #look into opened enrichment files
        allbams=names(BAMvariables$listBAM)
        if(!is.null(allbams) & length(allbams)>0){
          output$menuEnrichPredefPipeline<-renderUI({
            list(
              HTML("<br>"),
              list(HTML("<b>Select enrichment to associate to the new ROI:</b>"),htmlhelp("","help_prepheat_enrichmentAssoc")),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 400px; background-color: #ffffff;",
                    checkboxGroupInput(inputId="enrichAllPredefPipeline",label=NULL,choices=allbams)
              )
            )
          }) 


          if (isvalid(input$choiceSummitPredefPipeline)){
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






##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#observers for prepare for genelist (EASY)
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


#here update genelist menu
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

    output$show_selectROIGenelistPipeline<-renderUI({
      list(
        HTML("<b>Select gene list(s) to prepare for metagene profiles:</b>"),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
          checkboxGroupInput("selectROIGenelistPipeline",label=NULL,choices = historylist)
        )
      )
    })
 
  }else{
    output$show_selectROIGenelistPipeline<-renderUI({HTML("<font color='red'>No genelists available. Go to 'ROI' section to import genelists</font>")})  
  }
   
})



#observer for the BAM files common to the genelists selected
observe({
  input$selectROIGenelistPipeline
  if(length(ROIvariables$listROI)>0&length(input$selectROIGenelistPipeline)>0){
    #if we are here, promoters, transcripts and TES are available for the same genelist
    #for each gene, find promoters, transcripts and TES and the common BAM files associated
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    triadsROItoshow=c()
    allBAMmissings=list()
    for(i in 1:length(input$selectROIGenelistPipeline)){
      #change this if in the menu "genelist " is before the actual name of the genelist
      #selected=strsplit(input$selectROIGenelistPipeline[i],split=" ")[[1]][2]
      selected=input$selectROIGenelistPipeline[i]
      pos_promoters= which(paste("promoters_genelist_",selected,sep="")==nomi)
      pos_transcripts= which(paste("transcripts_genelist_",selected,sep="")==nomi)
      pos_TES= which(paste("TES_genelist_",selected,sep="")==nomi)
      triadsROItoshow=c(triadsROItoshow,paste("promoters_genelist_",selected,sep=""),
                                        paste("transcripts_genelist_",selected,sep=""),
                                        paste("TES_genelist_",selected,sep=""))
      pos=c(pos_promoters,pos_transcripts,pos_TES)
      if (length(pos_promoters)>0 & length(pos_transcripts)>0 & length(pos_TES)>0){
        rois=ROIvariables$listROI[pos] 
        allBAMmissing=c()
        bamblock=Enrichlist$rawcoverage[pos]
        #in the bampresent=names(getBAMlist(rois[[i]])) can raise an error, not always
        #it tells that unable to find an inherited method for function ‘getBAMlist’ for signature ‘"NULL"’
        for (k in 1:length(rois)){
          bampresent=names(bamblock[[k]])
          allbams=names(BAMvariables$listBAM)
          remainingbams=setdiff(allbams,bampresent)
          allBAMmissing=c(allBAMmissing,remainingbams)
        }
        #define the union of them and show in the menu
        allBAMmissing=unique(allBAMmissing)
        if(is.null(allBAMmissing)){
          allBAMmissing=character(0)
        }
       
      }else{
        allBAMmissing=character(0)
      }

      allBAMmissings[[i]]=allBAMmissing
    }
    #show triads ROI for selected genelists
    output$showROItriadGeneListPipeline<-renderUI({
      list(
        list(HTML("<b>ROIs constituting the selected gene list(s):</b>"),htmlhelp("","help_prepgenelist_ROIassociated")),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
          HTML(triadsROItoshow,collapse="<br>")
        )
      )
      
    }) 

    allBAMavailableToAssoc=Reduce(union, allBAMmissings)

    if(length(allBAMavailableToAssoc)!=0){
      historylist=as.list(allBAMavailableToAssoc)
      names(historylist)=allBAMavailableToAssoc
      #then, find BAM files allBAMavailableToAssoc to the genelists selected: 
      output$show_BAMforGenelistPipeline<-renderUI({
        list(
          list(HTML("<b>Select enrichments to associate:</b>"),htmlhelp("","help_prepgenelist_enrichmentAssoc")),
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
            checkboxGroupInput(inputId="BAMforGenelistPipeline",label=NULL,choices = historylist)
          )
        )
      })


    }else{
      output$show_BAMforGenelistPipeline<-renderUI({NULL})            
    }
  }else{
    output$show_BAMforGenelistPipeline<-renderUI({NULL})      
    #print nothing for ROI triads
    output$showROItriadGeneListPipeline<-renderUI({NULL})
  }
})


#observer for button of genelist pipeline
observe({
  input$selectROIGenelistPipeline
  input$BAMforGenelistPipeline
  if (!isvalid(input$selectROIGenelistPipeline) ){
    output$show_buttonGenelistPipeline<-renderUI({NULL})
    return()
  }
  if (isvalid(input$selectROIGenelistPipeline)& !isvalid(input$BAMforGenelistPipeline)){
    output$show_buttonGenelistPipeline<-renderUI({NULL})
    return()
  }


  nomi=unlist(lapply(ROIvariables$listROI,getName))
  triadsROItoshow=c()
  allBAMmissings=list()
  for(i in 1:length(input$selectROIGenelistPipeline)){
    #change this if in the menu "genelist " is before the actual name of the genelist
    #selected=strsplit(input$selectROIGenelistPipeline[i],split=" ")[[1]][2]
    selected=input$selectROIGenelistPipeline[i]
    pos_promoters= which(paste("promoters_genelist_",selected,sep="")==nomi)
    pos_transcripts= which(paste("transcripts_genelist_",selected,sep="")==nomi)
    pos_TES= which(paste("TES_genelist_",selected,sep="")==nomi)
    triadsROItoshow=c(triadsROItoshow,paste("promoters_genelist_",selected,sep=""),
                                      paste("transcripts_genelist_",selected,sep=""),
                                      paste("TES_genelist_",selected,sep=""))
    pos=c(pos_promoters,pos_transcripts,pos_TES)
    if (length(pos_promoters)>0 & length(pos_transcripts)>0 & length(pos_TES)>0){
      rois=ROIvariables$listROI[pos] 
      allBAMmissing=c()
      bamblock=Enrichlist$rawcoverage[pos]
      #in the bampresent=names(getBAMlist(rois[[i]])) can raise an error, not always
      #it tells that unable to find an inherited method for function ‘getBAMlist’ for signature ‘"NULL"’
      for (k in 1:length(rois)){
        bampresent=names(bamblock[[k]])
        allbams=names(BAMvariables$listBAM)
        remainingbams=setdiff(allbams,bampresent)
        allBAMmissing=c(allBAMmissing,remainingbams)
      }
      #define the union of them and show in the menu
      allBAMmissing=unique(allBAMmissing)
      if(is.null(allBAMmissing)){
        allBAMmissing=character(0)
      }
     
    }else{
      allBAMmissing=character(0)
    }

    allBAMmissings[[i]]=allBAMmissing
  }
  allBAMavailableToAssoc=Reduce(union, allBAMmissings)

  if(length(allBAMavailableToAssoc)==0){
    output$show_buttonGenelistPipeline<-renderUI({HTML("<font color='red'>All available enrichments already associated to selected genelist(s)</font>")})
    return()    
  }

  output$show_buttonGenelistPipeline<-renderUI({
    actionButton("PrepareROIgenelistPipeline","Prepare genelist")
  })

})












##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#observers for manual ROI management
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#observers to help buttons
observeEvent(input$help_roimanual_overlaps_selectref, {boxHelpServer(help_roimanual_overlaps_selectref)})
observeEvent(input$help_roimanual_overlaps_intersectionref, {boxHelpServer(help_roimanual_overlaps_intersectionref)})
observeEvent(input$help_roimanual_overlaps_unionref, {boxHelpServer(help_roimanual_overlaps_unionref)})
observeEvent(input$help_roimanual_overlaps_selectoverlapwith, {boxHelpServer(help_roimanual_overlaps_selectoverlapwith)})
observeEvent(input$help_roimanual_overlaps_intersectionoverlapwith, {boxHelpServer(help_roimanual_overlaps_intersectionoverlapwith)})
observeEvent(input$help_roimanual_overlaps_unionoverlapwith, {boxHelpServer(help_roimanual_overlaps_unionoverlapwith)})
observeEvent(input$help_roimanual_overlaps_selectnotoverlapwith, {boxHelpServer(help_roimanual_overlaps_selectnotoverlapwith)})
observeEvent(input$help_roimanual_overlaps_intersectionnotoverlapwith, {boxHelpServer(help_roimanual_overlaps_intersectionnotoverlapwith)})
observeEvent(input$help_roimanual_overlaps_unionnotoverlapwith, {boxHelpServer(help_roimanual_overlaps_unionnotoverlapwith)})
observeEvent(input$help_roimanual_overlaps_minimumbp, {boxHelpServer(help_roimanual_overlaps_minimumbp)})
observeEvent(input$help_roimanual_overlaps_strandspecific, {boxHelpServer(help_roimanual_overlaps_strandspecific)})





##########################################################################################
# overlaps
##########################################################################################

#update for options of overlap (min bp, name of new ROI, button) only after at least one selected ROI
observe({
  input$selectedROIs
  if(length(input$selectedROIs)<1){
    output$show_optionsOverlapButton<-renderUI({NULL})
    return()
  }else{
    output$show_optionsOverlapButton<-renderUI({
      list(
        list(HTML("<b>Minimum number of bp to consider for overlaps:</b>"),htmlhelp("","help_roimanual_overlaps_minimumbp")),
        numericInput( inputId = 'minOverlapNEWROI',label=NULL,min = 1, step = 5,value=1),
        checkboxInput("StrandSpecOverlapNEWROI", label=list("Strand-specific overlaps",htmlhelp("","help_roimanual_overlaps_strandspecific")),value = FALSE, width = NULL),

        HTML("<b>Name of new ROI:</b>"),
        fluidRow(
          column(width=8,
            textInput("ROIname",label=NULL,placeholder="type new ROI name here")
          ),
          column(width=4,
            actionButton("maketheROI", "Build ROI") 
          )
        )
      )
    })

  }
})





#update for radio menu for union/intersection of the main ROI choice. activate only if more than 1 selected
observe({
  input$selectedROIs
  if(length(input$selectedROIs)>1){
    output$RadiobuttonStartingROI<-renderUI({
      radioButtons(inputId="choiceROI",label="Criteria for building the aggregated reference ROI:" ,
                                  choiceNames=list(
                                    htmlhelp("Intersection","help_roimanual_overlaps_intersectionref"),
                                    htmlhelp("Union","help_roimanual_overlaps_unionref")
                                  ),
                                  choiceValues=list(
                                    "intersection",
                                    "union"
                                  ),
                                  selected="union"
                  )       
    })
  }else{
    output$RadiobuttonStartingROI<-renderUI({NULL})
  }
})


#react to ROI selected of overlaps to put overlapping and not overlapping lists
observe({
  input$selectedROIs
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
      list( list(HTML("Select ROI(s):"),htmlhelp("","help_roimanual_overlaps_selectoverlapwith")),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput(inputId="overlapROIs",label=NULL,choices = historylist)
            )
      )
    })
    output$columnNotOverlap<-renderUI({
      list(
            list(HTML("Select ROI(s):"),htmlhelp("","help_roimanual_overlaps_selectnotoverlapwith")),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput(inputId="notoverlapROIs",label=NULL,choices = historylist)
            )
      )
    })
  }else{
    output$columnOverlap<-renderUI({HTML(paste("(need at least one reference ROI selected)",sep=""))})
    output$columnNotOverlap<-renderUI({HTML(paste("(need at least one reference ROI selected)",sep=""))})
  }
})



#observer for union/intersection criteria of overlapping and not overlapping ROIs in overlaps
observe ({
  input$overlapROIs
  if(length(input$overlapROIs)>1){
    output$overlapcriteria<-renderUI({
      radioButtons("choiceoverlapROI",label="Criteria for building the aggregated contrast ROI:" ,
                              choiceNames=list(
                                htmlhelp("Intersection of the selected ROIs","help_roimanual_overlaps_intersectionoverlapwith"),
                                htmlhelp("Union of the selected ROIs","help_roimanual_overlaps_unionoverlapwith")
                              ),
                              choiceValues=list(
                                "stringent",
                                "allofthem"
                              ),
                              selected="allofthem"
      )        
    })
  }else{
    output$overlapcriteria<-renderUI({NULL})
  }
})
observe ({
  input$notoverlapROIs
  if(length(input$notoverlapROIs)>1){
    output$notoverlapcriteria<-renderUI({
      radioButtons("choicenotoverlapROI",label="Criteria for building the aggregated contrast ROI:" ,
                          choiceNames=list(
                            htmlhelp("Intersection of the selected ROIs","help_roimanual_overlaps_intersectionnotoverlapwith"),
                            htmlhelp("Union of the selected ROIs","help_roimanual_overlaps_unionnotoverlapwith")
                          ),
                          choiceValues=list(
                            "intersection",
                            "union"
                          ),
                          selected="union"
      ) 
    })
  }else{
    output$notoverlapcriteria<-renderUI({NULL})
  }
})




##########################################################################################
# associate enrichments
##########################################################################################

#observers for help buttons
observeEvent(input$help_roimanual_enrichmentAssoc_enrichments, {boxHelpServer(help_roimanual_enrichmentAssoc_enrichments)})
observeEvent(input$help_roimanual_enrichmentAssoc_NOnorm, {boxHelpServer(help_roimanual_enrichmentAssoc_NOnorm)})
observeEvent(input$help_roimanual_enrichmentAssoc_librarysize, {boxHelpServer(help_roimanual_enrichmentAssoc_librarysize)})
observeEvent(input$help_roimanual_enrichmentAssoc_customnorm, {boxHelpServer(help_roimanual_enrichmentAssoc_customnorm)})
observeEvent(input$help_roimanual_enrichmentAssoc_spikein, {boxHelpServer(help_roimanual_enrichmentAssoc_spikein)})
observeEvent(input$help_roimanual_enrichmentAssoc_numbercores, {boxHelpServer(help_roimanual_enrichmentAssoc_numbercores)})

observeEvent(input$help_roimanual_reorderenrich_index, {boxHelpServer(help_roimanual_reorderenrich_index)})




#observer for showing enrichments to associate or deassociate, based on the currently selected ROIs
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


    if (!isvalid(input$Assoc_or_deassoc)){
      output$show_selectBAMtoDeassociate<-renderUI({NULL})
      output$show_selectBAMtoassociate<-renderUI({checkboxGroupInput("selectBAMtoassociate",NULL,choices=NULL)})
      output$show_coresAndButtonAssociate<-renderUI({NULL}) 
      return()     
    }

    if (input$Assoc_or_deassoc=="associate"){
      if (length(allBAMmissing)>0){
        output$show_selectBAMtoassociate<-renderUI({
          list(
            #HTML("<h3>Associate enrichment(s) to ROI(s)</h3><br>"),
            list(HTML("<b>Select enrichment(s) to associate to ROI(s):</b>"),htmlhelp("","help_roimanual_enrichmentAssoc_enrichments"),HTML("<br>")),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 180px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("selectBAMtoassociate",NULL,choices=allBAMmissing)
            )
          )
        })
        output$show_coresAndButtonAssociate<-renderUI({
          list(
            list(HTML("<b>Number of cores:</b>"),htmlhelp("","help_roimanual_enrichmentAssoc_numbercores")),
            numericInput(inputId = 'coresCoverage',label=NULL,min = 1, max = nc, step = 1,value=1),
            HTML("<i>WARNING</i>: time consuming (some minutes)<br>"),
            actionButton("confirmBAMassociate", "Associate!")
          )
        }) 
        output$show_selectBAMtoDeassociate<-renderUI({NULL})        
      }else{
        if (length(allBAMavailable)>0){
          #ROI selected, but all enrichments available already associated to the selected ROIs
          output$show_selectBAMtoassociate<-renderUI({
            list(
              HTML("<font color='red'>All available enrichments already associated to selected ROI(s)</font>"),
              checkboxGroupInput("selectBAMtoassociate",NULL,choices=NULL)
            )   
          })
          output$show_coresAndButtonAssociate<-renderUI({NULL}) 
          output$show_selectBAMtoDeassociate<-renderUI({NULL})         
        }else{
          output$show_selectBAMtoassociate<-renderUI({            
            list(
              HTML("<font color='red'>All available enrichments already associated to selected ROI(s)</font>"),
              checkboxGroupInput("selectBAMtoassociate",NULL,choices=NULL)
            ) 
          })
          output$show_coresAndButtonAssociate<-renderUI({NULL})  
          output$show_selectBAMtoDeassociate<-renderUI({NULL})           
        }

      }
     
    }else{
      #here, de-association selected
      if (length(allBAMavailable)>0){
        output$show_selectBAMtoDeassociate<-renderUI({
          list(
            #HTML("<h3>Enrichments associated to ROI(s)</h3><br>"),
            HTML("<b>Select enrichment(s) to remove from selected ROI(s):</b><br>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 180px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("selectBAMtoDeassociate",NULL,choices=allBAMavailable)
            ),
            actionButton("confirmBAMDeassociate", "Remove")
          )
        })
        output$show_selectBAMtoassociate<-renderUI({NULL})
        output$show_coresAndButtonAssociate<-renderUI({NULL})
      }else{
        if (length(allBAMmissing)>0){
          output$show_selectBAMtoDeassociate<-renderUI({HTML("<font color='red'>No enrichments associated to the selected ROI(s)</font>")})
          output$show_selectBAMtoassociate<-renderUI({NULL})
          output$show_coresAndButtonAssociate<-renderUI({NULL})
        }else{
          #bam available ==0 and bam missing=0 => no enrichment at all!
          output$show_selectBAMtoDeassociate<-renderUI({HTML("<font color='red'>No enrichments available</font>")})
          output$show_selectBAMtoassociate<-renderUI({NULL})
          output$show_coresAndButtonAssociate<-renderUI({NULL})
        }
      }
    }

  }else{
    output$show_selectBAMtoDeassociate<-renderUI({NULL})
    output$show_selectBAMtoassociate<-renderUI({NULL})
    output$show_coresAndButtonAssociate<-renderUI({NULL})
  }
})


#observer for choice to assoc or deassoc enrichments
observe({
  input$selectROItoBAMassociate
  if (isvalid(input$selectROItoBAMassociate)){
    output$show_choiceAssocDeassoc<-renderUI({
      radioButtons("Assoc_or_deassoc",label="Associate or remove from selected ROI(s)?" ,
                                       choices=c("Associate enrichment(s)"="associate",
                                                 "Remove enrichment(s)"="deassociate"),
                                       selected=character(0)
      )
    })
  }else{
    output$show_choiceAssocDeassoc<-renderUI({NULL})  
  }
  
})



#observer for normalization menu (spike, no norm, library size...)
observe({
  input$selectROItoBAMassociate
  input$Assoc_or_deassoc
  if (!isvalid(input$selectROItoBAMassociate)){
    output$radioForNorm<-renderUI({NULL})
    return()    
  }
  if (!isvalid(input$Assoc_or_deassoc)){
    output$radioForNorm<-renderUI({NULL})
    return()
  }
  if (input$Assoc_or_deassoc=="deassociate"){
    output$radioForNorm<-renderUI({NULL})
    return()
  }
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
                                     choiceNames=list(
                                        htmlhelp("no normalization","help_roimanual_enrichmentAssoc_NOnorm"),
                                        htmlhelp("library-size","help_roimanual_enrichmentAssoc_librarysize"),
                                        htmlhelp("custom normalizer","help_roimanual_enrichmentAssoc_customnorm")
                                     ),
                                     choiceValues=list(
                                      "nonorm",
                                      "librarysize",
                                      "customnorm"
                                     ),selected="librarysize"
                        )
      })
    #spike in normalization only if we have more than 3 bam file (current and input+spikeInput)
    }else if (length(input$selectBAMtoassociate)>0 & length(BAMvariables$listBAM)>=3){
      output$radioForNorm<-renderUI({
        radioButtons("selectMethodForNorm",label="Normalization method:" ,
                                     choiceNames=list(
                                        htmlhelp("no normalization","help_roimanual_enrichmentAssoc_NOnorm"),
                                        htmlhelp("library-size","help_roimanual_enrichmentAssoc_librarysize"),
                                        htmlhelp("custom normalizer","help_roimanual_enrichmentAssoc_customnorm"),
                                        htmlhelp("spike-in","help_roimanual_enrichmentAssoc_spikein")
                                     ),
                                     choiceValues=list(
                                      "nonorm",
                                      "librarysize",
                                      "customnorm",
                                      "spikein"
                                     ),selected="librarysize"
                        )
      })
    }else{
      #null menu
      output$radioForNorm<-renderUI({NULL})
    }    
  }
})


#observer for enrichment after choosing normalization:
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






#based on the radiobutton of BAM rename, edit rename enrichment menu
observe({
  input$selectROIforBAMrename
  if(!isvalid(ROIvariables$listROI) | !isvalid(input$selectROIforBAMrename) ){
    output$show_BAMrename<-renderUI({NULL})
    return()
  }
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  pos=match(input$selectROIforBAMrename,nomi)
  roi=ROIvariables$listROI[[pos]]
  if (!is.null(roi) & length(BAMvariables$listBAM)>0){
    getbam=names(Enrichlist$rawcoverage[[pos]])
    if (length(getbam)!=0){
      output$show_BAMrename<-renderUI({
        list(
          selectInput("selectedBAMtoRename",label="Select enrichment to rename:",choices=getbam),
          textInput("newBAMname","Select new name:",placeholder="type new enrichment name here",value=""),
          actionButton("renameBAM", "Rename")
        )
      })
    }else{
      output$show_BAMrename<-renderUI({ HTML("No enrichments associated to the selected ROI...")})
    }
  }else{
    output$show_BAMrename<-renderUI({ HTML("No ROI...")})
  }
})


#based on the radiobutton of BAM reorder, edit reorder enrichment menu
observe({
  input$selectROIforBAMreorder
  if(!isvalid(ROIvariables$listROI) | !isvalid(input$selectROIforBAMreorder) ){
    output$show_BAMreorder<-renderUI({NULL})
    return()
  }
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  pos=match(input$selectROIforBAMreorder,nomi)
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

      output$show_BAMreorder<-renderUI({
        list(
          list(HTML("<b>Change the index for the new enrichments order</b>"),htmlhelp("","help_roimanual_reorderenrich_index")),
          wellPanel(id = "logPanelBAM",style = "overflow-y:scroll; max-height: 400px",
            fluidPage(
              lista
            )
          ),
          actionButton("reorderBAM","Reorder!")
        )
      })

    }else{
      output$show_BAMreorder<-renderUI({ HTML("At least 2 enrichments must be associated to the selected ROI...")})
    }
  }else{
    output$show_BAMreorder<-renderUI({ HTML("No ROI...")})
  }


})






##########################################################################################
# resize ROI
##########################################################################################


#observer for name of new ROI and button
observe({
  if (!isvalid(input$selectROItoresize)){
    output$show_namebuttonres<-renderUI({NULL})
    return()
  }
  output$show_namebuttonres<-renderUI({
    list(

      HTML("<br>"),
      HTML("<b>Name of new ROI:</b>"),
      fluidRow(
        column(width=6,
          textInput("ROInameResize",label=NULL,placeholder="type new ROI name here")
        ),
        column(width=6,
          actionButton("resizeROI","Create ROI")
        )
      )

    )
  })
})

#observer for type of resize
observe({
  input$chooseResizeType
  if (!isvalid(input$selectROItoresize)){
    output$show_choosePointResize<-renderUI({NULL})
    return()
  }
  ROI=isolate(input$selectROItoresize)
  nomi=unlist(lapply(isolate(ROIvariables$listROI),getName))
  pos=match(ROI,nomi)
  roi=isolate(ROIvariables$listROI)[[pos]] 

  nameroi=getName(roi)
  flag=getFlag(roi) #transcriptFlag, Pattern, normalFlag

  if (flag=="Pattern"){
    toput="pattern"
  }else if (flag=="promoterFlag"){
    toput="TSS"
  }else{
    toput="midpoint"
  }
  if (!isvalid(input$chooseResizeType)){
    output$show_choosePointResize<-renderUI({NULL})
    return()
  }
  if (input$chooseResizeType=="fixedVal"){
    output$show_choosePointResize<-renderUI({
      radioButtons("choosePointResize",label="Choose fixed point for resize:",
                      choiceNames=list(
                              paste("From",toput),
                              "From starts",
                              "From ends"
                            ),choiceValues=list("fromMid","fromStart","fromEnd"),selected="fromMid"
      )
    })

  }else{
    output$show_choosePointResize<-renderUI({
      radioButtons("choosePointResize",label="Choose if increment/decrement width:",choices=c(
                              "Increment width from midpoint/TSS"="increment",
                              "Decrement width from midpoint/TSS"="decrement"
                            ),selected="increment"
      ) 
    })
  }
})


#observer for intervals or % intervals for resize
observe({
  choice=input$choosePointResize
  if (length(choice)>0){
    if(choice=="fromMid" |choice=="fromStart" | choice=="fromEnd"){

      if(choice=="fromMid"){label_fix="Range center"}
      if(choice=="fromStart"){label_fix="Range start"}
      if(choice=="fromEnd"){label_fix="Range end"}

      output$showFixedMenuResize<-renderUI({
        list(
          fluidRow(
            column(width=5,HTML( paste("<b>Upstream (bp from",label_fix,")</b>"))),

            column(width=5,HTML( paste("<b>Downstream (bp from",label_fix,")</b>")))
          ),
          fluidRow(
            column(width=5,
              numericInput(inputId = 'sliderUpstreamROI',label=NULL,min = 0, max = 20000, step = 50,value=1000) 
            ),

            column(width=5,
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



##########################################################################################
# filter for width ROI
##########################################################################################

#observer for help
observeEvent(input$help_roimanual_filterwidth_min, {boxHelpServer(help_roimanual_filterwidth_min)})
observeEvent(input$help_roimanual_filterwidth_max, {boxHelpServer(help_roimanual_filterwidth_max)})


#observer for name of new ROI and button
observe({
  input$selectROItoFilterWIDTH
  if (!isvalid(input$selectROItoFilterWIDTH)){
    output$show_namebuttonwidth<-renderUI({NULL})
    return()
  }
  output$show_namebuttonwidth<-renderUI({
    list(

      HTML("<br>"),
      HTML("<b>Name of new ROI:</b>"),
      fluidRow(
        column(width=6,
          textInput("ROInameFilterWIDTH",label=NULL,placeholder="type new ROI name here")
        ),
        column(width=6,
          actionButton("FilterROIWIDTH","Create ROI")
        )
      )
    )
  })
})





#observer to change slider input/ quantile and absolute width thresholds
observe({
  input$selectROItoFilterWIDTH
  if(length(input$selectROItoFilterWIDTH)>0){
    #if a BAM is available, take the bam from the roi selected, and put 0.9 as quantile, and the value corresponding to that 
    #in sum of base coverage for each range
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoFilterWIDTH,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      widths=width(getRange(roi))
      maxx=max(widths)
      quant1=round(quantile(widths,0.9))
      quant2=round(quantile(widths,0))
      fract=round(maxx/100)
      output$show_quantileThreshFilterWIDTH<-renderUI({
        sliderInput('quantileThreshFilterWIDTH',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
      })
      
      output$show_absoluteFiltersWIDTH<-renderUI({
        list(
          fluidRow(
            column(width=4,list(HTML("<b>MIN width (bp)</b>"),htmlhelp("","help_roimanual_filterwidth_min"))  ),
            column(width=4),
            column(width=4,list(HTML("<b>MAX width (bp)</b>"),htmlhelp("","help_roimanual_filterwidth_max"))   )
          ),
          fluidRow(
            column(width=4,
              numericInput(inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
            ),
            column(width=4,
              HTML("<h4>      &lt; --- &gt; </h4>")
            ),
            column(width=4,
              numericInput(inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2)) 
            )
          )
        )
      })


    }else{
      output$show_quantileThreshFilterWIDTH<-renderUI({NULL})
      output$show_absoluteFiltersWIDTH<-renderUI({NULL})
    }
  }else{
    output$show_quantileThreshFilterWIDTH<-renderUI({NULL})
    output$show_absoluteFiltersWIDTH<-renderUI({NULL})
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


      output$show_absoluteFiltersWIDTH<-renderUI({
        list(
          fluidRow(
            column(width=4,list(HTML("<b>MIN width (bp)</b>"),htmlhelp("","help_roimanual_filterwidth_min"))),
            column(width=4),
            column(width=4,list(HTML("<b>MAX width (bp)</b>"),htmlhelp("","help_roimanual_filterwidth_max")))
          ),
          fluidRow(
            column(width=4,
              numericInput(inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
            ),
            column(width=4,
              HTML("<h4>      &lt; --- &gt; </h4>")
            ),
            column(width=4,
              numericInput(inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2)) 
            )
          )
        )
      })

    }else{
      output$show_absoluteFiltersWIDTH<-renderUI({NULL})
    }
  }else{
    output$show_absoluteFiltersWIDTH<-renderUI({NULL})
  }
})



##########################################################################################
# random sample ROI
##########################################################################################

#observer for name of new ROI and button
observe({
  input$selectROItoSample
  if (!isvalid(input$selectROItoSample)){
    output$show_namebuttonsample<-renderUI({NULL})
    return()
  }
  output$show_namebuttonsample<-renderUI({
    list(
      HTML("<br>"),
      HTML("<b>Name of new ROI:</b>"),
      fluidRow(
        column(width=6,
          textInput("ROInameSample",label=NULL,placeholder="type new ROI name here")
        ),
        column(width=6,
          actionButton("SampleROI","Create ROI")
        )
      )
    )
  })
})



#observer for quantile bar and absolute value for number of sample given selected ROIs
observe({
  input$selectROItoSample
  if(length(input$selectROItoSample)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoSample,nomi)
    roi=ROIvariables$listROI[[pos]]
    if (!is.null(roi)){
      maxx=length(getRange(roi))
      quant=round(quantile(1:maxx,0.3))
      fract=round(maxx/100)
      output$show_quantileThreshSample<-renderUI({
        sliderInput('quantileThreshSample',label="Fraction of genomic ranges to keep:",min = 0, max = 1, value = 0.3,step=0.05)
      })
      output$show_numberSample<-renderUI({
        numericInput(inputId = 'numberSample',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant))
      })

    }else{
      output$show_quantileThreshSample<-renderUI({NULL})
      output$show_numberSample<-renderUI({NULL})
    }
  }else{
    output$show_quantileThreshSample<-renderUI({NULL})
    output$show_numberSample<-renderUI({NULL})
  }
})


#observer for absolute value for number of sample given quantile bar
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
      output$show_numberSample<-renderUI({
        numericInput(inputId = 'numberSample',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant))
      })  
      
    }else{
      output$show_numberSample<-renderUI({NULL})
    }
  }else{
    output$show_numberSample<-renderUI({NULL})
  }
})




##########################################################################################
# extract sequence pattern from ROI
##########################################################################################

#observers for help buttons
observeEvent(input$help_roimanual_pattern_IUPAC, {boxHelpServer(help_roimanual_pattern_IUPAC)})
observeEvent(input$help_roimanual_pattern_bothstrands, {boxHelpServer(help_roimanual_pattern_bothstrands)})
observeEvent(input$help_roimanual_pattern_strandspecific, {boxHelpServer(help_roimanual_pattern_strandspecific)})



#observer for strand option of extracting pattern from ROI:
observe({
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
          output$showStrandOptsPattern<-renderUI({HTML("ROI is not strand-specific: pattern will be searched on both strands<br>")})
        }else{
          #stranded (+ and -) => offer the choice
          output$showStrandOptsPattern<-renderUI({
            list(
              HTML("<b>Strand selection:</b><br>"),
              radioButtons("strandOptsPattern",label=NULL ,
                                    choiceNames=list(
                                      htmlhelp("Both strands","help_roimanual_pattern_bothstrands"),
                                      htmlhelp("Strand-specific","help_roimanual_pattern_strandspecific")
                                    ),
                                    choiceValues=list(
                                      "bothStrands",
                                      "strandSpecific"
                                    ),selected="strandSpecific"
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



#find current BSgenome in the installed packages. If not found, put a button for installation.
#look DATABASEvariables$currentASSEMBLY. If not null, search BSgenome in the package. If 
#present, import, put in DBvariables$BSgenomeDB and do nothing (maybe print what we are using), 
#otherwise, downloadbutton 
observe({
  input$selectROItoExtractPattern
  if (!isvalid(input$selectROItoExtractPattern)){
    output$showWarningBSgenome<-renderUI({NULL})
    output$show_namebuttonpattern<-renderUI({NULL})
    output$show_PatternToSearch<-renderUI({NULL})
    return()
  }
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

    if (is.na(pos)){
      output$showWarningBSgenome<-renderUI({HTML("<font color='red'>BSegnome not available. It seems that you are using a non-standard genome</font>")})
      output$show_PatternToSearch<-renderUI({NULL})
      output$showStrandOptsPattern<-renderUI({NULL})
      output$show_namebuttonpattern<-renderUI({NULL})      
      return()
    }

    BSstring=paste("BSgenome.",org[pos],".UCSC.",asm,sep="")
    x=rownames(installed.packages())
    pos_pkg=match(BSstring,x)
    if(!is.na(pos_pkg)){
      #BSgenome found. Import library and do nothing
      library(BSstring,character.only=TRUE)
      output$showWarningBSgenome<-renderUI({NULL})
      output$show_namebuttonpattern<-renderUI({
        list(
          HTML("<br>"),
          HTML("<b>Name of new ROI:</b>"),
          fluidRow(
            column(width=6,
              textInput("ROInamePattern",label=NULL,placeholder="type new ROI name here")
            ),
            column(width=6,
              actionButton("ExtractPatternROI","Create ROI")
            )
          )
        )
      })

      output$show_PatternToSearch<-renderUI({
        textInput("PatternToSearch",label=list("Select pattern (IUPAC nomenclature)",htmlhelp("","help_roimanual_pattern_IUPAC")),placeholder="ATCNYGG")
      })
    }else{
      #button for download the package
      output$showWarningBSgenome<-renderUI({
        list(
        HTML(paste("<font color='red'>Need ",BSstring," package to be installed. Install it now?</font><br>",sep="")),
        actionButton("DownloadBSgenome","Download")
        )
      })
      output$show_PatternToSearch<-renderUI({NULL})
      output$showStrandOptsPattern<-renderUI({NULL})
      output$show_namebuttonpattern<-renderUI({NULL})

      #here should nullify all other options!

    }

  }else{
    #warning message: I need database of a genome assembly
    output$showWarningBSgenome<-renderUI({HTML("<font color='red'>Warning: choose a genome assembly from 'Databases' section</font>")})
    #here should nullify all other options!
    output$show_PatternToSearch<-renderUI({NULL})
    output$showStrandOptsPattern<-renderUI({NULL})
    output$show_namebuttonpattern<-renderUI({NULL})
  }

})











##########################################################################################
# center ROI on summit
##########################################################################################

#observer for help button
observeEvent(input$help_roimanual_summit_enrichmenttouse, {boxHelpServer(help_roimanual_summit_enrichmenttouse)})



#observer for updating enrichment list for summit and text/button
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
        output$show_selectBAMtoCenterSummit<-renderUI({
          selectInput("selectBAMtoCenterSummit", label=list("Enrichment to use for summit:",htmlhelp("","help_roimanual_summit_enrichmenttouse")),choices=getbam)
        })

        output$show_namebuttonsummit<-renderUI({
          list(
            HTML("<br>"),
            HTML("<b>Name of new ROI:</b>"),
            fluidRow(
              column(width=6,
                textInput("ROInameSummit",label=NULL,placeholder="type new ROI name here")
              ),
              column(width=6,
                actionButton("SummitROI","Create ROI")
              )
            )
          )
        })
      }else{
        output$show_selectBAMtoCenterSummit<-renderUI({HTML("<font color='red'>At least one enrichment must be assciated to the selected ROI</font>")})
        output$show_namebuttonsummit<-renderUI({NULL})
      }
    }else{
      output$show_selectBAMtoCenterSummit<-renderUI({NULL})
      output$show_namebuttonsummit<-renderUI({NULL})
    }
  }else{
    output$show_selectBAMtoCenterSummit<-renderUI({NULL})
    output$show_namebuttonsummit<-renderUI({NULL})
  }
})




##########################################################################################
# filter ROI on enrichments
##########################################################################################

#observers for helps
observeEvent(input$help_roimanual_filterenrich_min, {boxHelpServer(help_roimanual_filterenrich_min)})
observeEvent(input$help_roimanual_filterenrich_max, {boxHelpServer(help_roimanual_filterenrich_max)})


#observer for update enrichment list in filter for enrichment
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
        output$show_selectBAMtoFilter<-renderUI({
          selectInput("selectBAMtoFilter", label="Enrichment to use for filtering:",choices=getbam)
        })
        
      }else{
        output$show_selectBAMtoFilter<-renderUI({HTML("<font color='red'>At least one enrichment must be assciated to the selected ROI</font>")})
      }
    }else{
      output$show_selectBAMtoFilter<-renderUI({NULL})
    }
  }else{
    output$show_selectBAMtoFilter<-renderUI({NULL})
  }
})



#need observer for enrichment thresholds reactive to each other
observe({
  
  input$selectBAMtoFilter
  if (!isvalid(isolate(input$selectROItoFilter))){
    output$show_quantileThreshFilter<-renderUI({NULL })
    output$show_absoluteFilters<-renderUI({NULL })
    return()
  }
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
        sums=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bamselected_tmp,Nbins=1,Snorm=FALSE,key=decrkeyselected,norm_factor=normfactselected)
        maxx=max(sapply(sums,max))
        quant1=round(quantile(sums,0.9))
        quant2=round(quantile(sums,0))
        fract=round(maxx/100)


        output$show_quantileThreshFilter<-renderUI({
          sliderInput('quantileThreshFilter',label="Select quantiles:",min = 0, max = 1, value = c(0,0.9),step=0.002)
        })

        output$show_absoluteFilters<-renderUI({
          list(
            fluidRow(
              column(width=4,list(HTML("<b>MIN enrichment</b>"),htmlhelp("","help_roimanual_filterenrich_min"))  ),
              column(width=4),
              column(width=4,list(HTML("<b>MAX enrichment</b>"),htmlhelp("","help_roimanual_filterenrich_max"))  )
            ),
            fluidRow(
              column(width=4,
                numericInput(inputId = 'absoluteFilter1',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
              ),
              column(width=4,
                HTML("<h4>      &lt; --- &gt; </h4>")
              ),
              column(width=4,
                numericInput(inputId = 'absoluteFilter2',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2))
              )
            )
          )
        })

        output$show_namebuttonFilterenrich<-renderUI({
          list(
            HTML("<br>"),
            HTML("<b>Name of new ROI:</b>"),
            fluidRow(
              column(width=6,
                textInput("ROInameFilter",label=NULL,placeholder="type new ROI name here")
              ),
              column(width=6,
                actionButton("FilterROI","Create ROI")
              )
            )
          )
        })        
        
      }else{
        output$show_quantileThreshFilter<-renderUI({NULL })
        output$show_absoluteFilters<-renderUI({NULL })
        output$show_namebuttonFilterenrich<-renderUI({NULL})
      }
    }else{
      output$show_quantileThreshFilter<-renderUI({NULL })
      output$show_absoluteFilters<-renderUI({NULL })
      output$show_namebuttonFilterenrich<-renderUI({NULL})
    }
  }else{
    output$show_quantileThreshFilter<-renderUI({NULL })
    output$show_absoluteFilters<-renderUI({NULL })
    output$show_namebuttonFilterenrich<-renderUI({NULL})
  }
})



#observer for adapt absolute value thresholds if quantile changes
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
        sums=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bamselected_tmp,Nbins=1,Snorm=FALSE,key=decrkeyselected,norm_factor=normfactselected)
        maxx=max(sapply(sums,max))
        quant1=round(quantile(sums,input$quantileThreshFilter[1]))
        quant2=round(quantile(sums,input$quantileThreshFilter[2]))
        fract=round(maxx/100)
        output$show_absoluteFilters<-renderUI({
          list(
            fluidRow(
              column(width=4,list(HTML("<b>MIN enrichment</b>"),htmlhelp("","help_roimanual_filterenrich_min"))),
              column(width=4),
              column(width=4,list(HTML("<b>MAX enrichment</b>"),htmlhelp("","help_roimanual_filterenrich_max")))
            ),
            fluidRow(
              column(width=4,
                numericInput(inputId = 'absoluteFilter1',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant1))
              ),
              column(width=4,
                HTML("<h4>      &lt; --- &gt; </h4>")
              ),
              column(width=4,
                numericInput(inputId = 'absoluteFilter2',label=NULL,min = 0, max = maxx, step = fract,value=unname(quant2))
              )
            )
          )
        })
      }else{
        output$show_absoluteFilters<-renderUI({NULL })
      }
    }else{
      output$show_absoluteFilters<-renderUI({NULL })
    }
  }else{
    output$show_absoluteFilters<-renderUI({NULL })
  }
})






######################################################################################
#reorder ROIs
######################################################################################

#observer for help
observeEvent(input$help_roimanual_reorderroi_index, {boxHelpServer(help_roimanual_reorderroi_index)})


#observer for ROI menu with all indexes for reorder
observe({
  if (!is.null(ROIvariables$listROI) & length(ROIvariables$listROI)>=1){
    lista=list()
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    choicelist=as.list(1:length(nomi))
    names(choicelist)=as.character(1:length(nomi))
    for (i in 1:length(nomi)){
      lista[[i]]=fluidRow(column(3,
                              selectInput(inputId = paste("reorderoptionROI",i,sep=""), label = NULL, 
                            choices = choicelist,selected=i)),
              column(4,HTML(nomi[i])))
    }
    output$dinamicROI<-renderUI({
      list(
        list(HTML("<b>Change the index for the new ROIs order</b>"),htmlhelp("","help_roimanual_reorderroi_index")),
        wellPanel(id = "logPanelROI",style = "overflow-y:scroll; max-height: 300px",
          fluidPage(
            lista
          )
        ),
        #button to confirm to reorder
        actionButton("reorderROI","Reorder!")
      )
      
    })
  }else{
    output$dinamicROI<-renderUI({HTML("<font color='red'>Need to have at least one ROI</font>")})
    
  }

}) 



######################################################################################
#view/edit/save notes for a ROI
######################################################################################  


#observer for updating text of notes for selected ROI (input$selectROItoEditNotes)
observe({
  input$selectROItoEditNotes
  if (!isvalid(ROIvariables$listROI) | !isvalid(input$selectROItoEditNotes)){
    output$ViewNotesToEdit<-renderUI({NULL})
    return()
  }

  nomi=unlist(lapply(ROIvariables$listROI,getName))
  pos=match(input$selectROItoEditNotes,nomi)
  roi=ROIvariables$listROI[[pos]]
  text=getSource(roi)
  windowNotesToShow<-list(textAreaInput(inputId="ROInotes",label="View and edit the notes of the selected ROI:",value=text,cols=60,rows=10),
                  actionButton("saveNotesROI", "save notes"),
                  downloadButton('downloadNotesROI', 'Download notes'))
  output$ViewNotesToEdit<-renderUI({windowNotesToShow})

})



##########################################################################################
##########################################################################################
##########################################################################################
# manual manipulation menu
##########################################################################################
##########################################################################################
##########################################################################################

#main observer for manual ROI page
observe({
  input$chooseManualManipulation
  if (!isvalid(input$chooseManualManipulation)){
    output$show_manualROImanipulationPage<-renderUI({NULL})
    return()
  }

  if (input$chooseManualManipulation=="overlaps"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=9,
          fluidRow(
            column(width=4,
              HTML("<h4><b>Choose the reference ROI:</b></h4>"),
              htmlhelp("Select ROI(s):","help_roimanual_overlaps_selectref"),
              uiOutput("show_selectedROIs"),
              uiOutput("RadiobuttonStartingROI")          
            ),
            column(width=4,
              HTML("<h4><b>...that overlaps with:</b></h4>"),  
              uiOutput("columnOverlap"),
              uiOutput("overlapcriteria") 
            ),
            column(width=4,
              HTML("<h4><b>...that doesn't overlap with:</b></h4>"),
              uiOutput("columnNotOverlap"),
              uiOutput("notoverlapcriteria")
            )
          )  
        ),
        column(width=3,
          uiOutput("show_optionsOverlapButton")
        )

      )

    })
  }else if (input$chooseManualManipulation=="associateEnrichments"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=6,
          uiOutput("show_selectROItoBAMassociate")
        ),
        column(width=6,
          uiOutput("show_choiceAssocDeassoc"),
          uiOutput("show_selectBAMtoassociate"),
          uiOutput("radioForNorm"),
          uiOutput("menuForNorm"),
          uiOutput("show_coresAndButtonAssociate"),
          uiOutput("show_selectBAMtoDeassociate")                          
        )
      )
    })
  }else if (input$chooseManualManipulation=="renameEnrichments"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=5,
          uiOutput("show_selectROIforBAMrename"),
          uiOutput("show_BAMrename")
        ),
        column(width=7)
      )
    })
  }else if (input$chooseManualManipulation=="reorderEnrichments"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=5,
          uiOutput("show_selectROIforBAMreorder"),
          uiOutput("show_BAMreorder")
        ),
        column(width=7)
      )
    })

  }else if (input$chooseManualManipulation=="resize"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=4,
          uiOutput("show_selectROItoresize"),
          uiOutput("show_choosePointResize"),
          uiOutput("showFixedMenuResize"),
          uiOutput("show_namebuttonres")
        ),
        column(width=8)
      )
    })
  }else if (input$chooseManualManipulation=="filterWidth"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=4,
          uiOutput("show_selectROItoFilterWIDTH"),
          uiOutput("show_quantileThreshFilterWIDTH"),
          uiOutput("show_absoluteFiltersWIDTH"),
          uiOutput("show_namebuttonwidth")
        ),
        column(width=6)
      )


    })
  }else if (input$chooseManualManipulation=="sample"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=4,
          uiOutput("show_selectROItoSample"),
          uiOutput("show_quantileThreshSample"),
          uiOutput("show_numberSample"),
          uiOutput("show_namebuttonsample")
        )
      )
    })
  }else if (input$chooseManualManipulation=="pattern"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=4,
          uiOutput("show_selectROItoExtractPattern"),
          #warning in case BSgenome DB not present
          uiOutput("showWarningBSgenome"),
          #motif text input
          uiOutput("show_PatternToSearch"),
          #menu or text to choose if both strands or strand in range
          uiOutput("showStrandOptsPattern"),
          #name of the new motif
          uiOutput("show_namebuttonpattern")
        )
      )

    })
  }else if (input$chooseManualManipulation=="summit"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=4,
          uiOutput("show_selectROItoCenterSummit"),
          uiOutput("show_selectBAMtoCenterSummit"),
          uiOutput("show_namebuttonsummit")
        )
      )
    })
  }else if (input$chooseManualManipulation=="filterEnrichment"){
    output$show_manualROImanipulationPage<-renderUI({
      fluidRow(
        column(width=4,
          uiOutput("show_selectROItoFilter"),
          uiOutput("show_selectBAMtoFilter"),
          uiOutput("show_quantileThreshFilter"),
          uiOutput("show_absoluteFilters"),
          uiOutput("show_namebuttonFilterenrich")
        )
      )
    })
  }else if (input$chooseManualManipulation=="renameRoi"){
    output$show_manualROImanipulationPage<-renderUI({
      column(width=4,
        uiOutput("show_selectROItoRename")
      )
    })
  }else if (input$chooseManualManipulation=="reorderRoi"){
    output$show_manualROImanipulationPage<-renderUI({
      column(width=4,
        uiOutput("dinamicROI")
      )  
    })
  }else if (input$chooseManualManipulation=="editNotes"){
    output$show_manualROImanipulationPage<-renderUI({
      column(width=4,
        uiOutput("show_selectROItoEditNotes"),
        uiOutput("ViewNotesToEdit")
      )
    })
  }


})






##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#main observer for ROI manipulation
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

###react to help buttons for box   
observeEvent(input$msg_ROImanipulation, {
  boxHelpServer(msg_ROImanipulation)
})
###react to help buttons for all choices
observeEvent(input$msg_prepare_ultraeasy, {
  boxHelpServer(msg_prepare_ultraeasy)
})
observeEvent(input$msg_prepare_forheat, {
  boxHelpServer(msg_prepare_forheat)
})
observeEvent(input$msg_prepare_formetagene, {
  boxHelpServer(msg_prepare_formetagene)
})
observeEvent(input$msg_prepare_manual, {
  boxHelpServer(msg_prepare_manual)
})

#enrichment association menu in ultra-easy
observeEvent(input$help_ultraeasy_enrichmentAssoc, {
  boxHelpServer(help_ultraeasy_enrichmentAssoc)
})

#options in prepare for heatmaps
observeEvent(input$help_prepheat_subsample, {
  boxHelpServer(help_prepheat_subsample)
})
observeEvent(input$help_prepheat_summit, {
  boxHelpServer(help_prepheat_summit)
})
observeEvent(input$help_prepheat_resize, {
  boxHelpServer(help_prepheat_resize)
})
observeEvent(input$help_prepheat_enrichmentAssoc, {
  boxHelpServer(help_prepheat_enrichmentAssoc)
})

#options in prepare for genelists
observeEvent(input$help_prepgenelist_ROIassociated, {
  boxHelpServer(help_prepgenelist_ROIassociated)
})
observeEvent(input$help_prepgenelist_enrichmentAssoc, {
  boxHelpServer(help_prepgenelist_enrichmentAssoc)
})


#options in prepare manual
observeEvent(input$help_roimanual_overlaps, {boxHelpServer(help_roimanual_overlaps)})

observeEvent(input$help_roimanual_resize, {boxHelpServer(help_roimanual_resize)})
observeEvent(input$help_roimanual_pattern, {boxHelpServer(help_roimanual_pattern)})
observeEvent(input$help_roimanual_summit, {boxHelpServer(help_roimanual_summit)})
observeEvent(input$help_roimanual_enrichmentAssoc, {boxHelpServer(help_roimanual_enrichmentAssoc)})
observeEvent(input$help_roimanual_subsample, {boxHelpServer(help_roimanual_subsample)})
observeEvent(input$help_roimanual_filterwidth, {boxHelpServer(help_roimanual_filterwidth)})
observeEvent(input$help_roimanual_filterenrich, {boxHelpServer(help_roimanual_filterenrich)})
observeEvent(input$help_roimanual_renameroi, {boxHelpServer(help_roimanual_renameroi)})
observeEvent(input$help_roimanual_editnotes, {boxHelpServer(help_roimanual_editnotes)})
observeEvent(input$help_roimanual_renameenrich, {boxHelpServer(help_roimanual_renameenrich)})
observeEvent(input$help_roimanual_reorderroi, {boxHelpServer(help_roimanual_reorderroi)})
observeEvent(input$help_roimanual_reorderenrich, {boxHelpServer(help_roimanual_reorderenrich)})
















#observer for initial choice of the pipeline
observe({
		input$choose_ROImanipulation_pipeline
    logvariables$temporary_ROI_present$ROIexisting
		if (!isvalid(input$choose_ROImanipulation_pipeline)){
			output$show_ROImanipulationUI<-renderUI({NULL})
			return()
		}

    if (input$choose_ROImanipulation_pipeline=="prepare_easy"){
      output$show_ROImanipulationUI<-renderUI({
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title="Prepare ROI basic (ULTRA-EASY)",
          fluidRow(
            column(width=6,
              uiOutput("show_selectROIultraeasy"),
              uiOutput("show_selectBAMultraeasy"),
              uiOutput("show_Buttonultraeasy")
            ),
            column(width=6)
          )
        )
      })

		}else if (input$choose_ROImanipulation_pipeline=="prepare_heat"){
			#GUI is the "prepare for heatmap"
      output$show_ROImanipulationUI<-renderUI({
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title="Prepare for heatmaps (EASY)",
          fluidRow(
            column(width=6,
              uiOutput("show_selectROIpredefPipeline"),
              uiOutput("show_quantileThreshPredefPipeline"),
              uiOutput("textNumRangesPredefPipeline"),
              uiOutput("show_radioButtonPredefSummit"),
              uiOutput("menuSummitPredefPipeline"),
              uiOutput ("show_updownstreamPredefPipeline"),
              uiOutput("menuEnrichPredefPipeline"),
              uiOutput("show_ROInamePredefPipeline")
            ),
            column(width=6)
          )
        )
      })   



		}else if (input$choose_ROImanipulation_pipeline=="prepare_metagene"){
      output$show_ROImanipulationUI<-renderUI({
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title="Prepare for metagene profiles (EASY)",
          fluidRow(
            column(width=6,
              uiOutput("show_selectROIGenelistPipeline"),
              uiOutput("showROItriadGeneListPipeline"),
              uiOutput("show_BAMforGenelistPipeline"),
              uiOutput("show_buttonGenelistPipeline")
            ),
            column(width=6)
          )
        )
      })  
















		}else if (input$choose_ROImanipulation_pipeline=="prepare_manual"){
      output$show_ROImanipulationUI<-renderUI({
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title="Manual manipulation (ADVANCED)",
          fluidRow(
            column(width=2,style='border-right: 1px solid black',
              #fixed radiobutton with all possible manipulations (overlaps, resize....)
              HTML("<font color='grey'>&emsp;&nbsp;&nbsp;Modify coordinates</font><br>"),
              radioButtons("chooseManualManipulation",label=NULL,choiceNames=list(
                              htmlhelp("Overlaps between ROIs","help_roimanual_overlaps"),
                              htmlhelp("Resize ROI","help_roimanual_resize"),
                              htmlhelp("Extract sequence pattern","help_roimanual_pattern"),
                              list(htmlhelp("Center on summit","help_roimanual_summit"),HTML("<br><br><font color='grey'>Enrichment association</font>")),
                              list(htmlhelp("Associate/remove enrichments","help_roimanual_enrichmentAssoc"),HTML("<br><br><font color='grey'>Filtering</font>")),
                              htmlhelp("Random sample","help_roimanual_subsample"),
                              htmlhelp("Filter for range width","help_roimanual_filterwidth"),
                              list(htmlhelp("Filter for enrichment","help_roimanual_filterenrich"),HTML("<br><br><font color='grey'>Editing</font>")),
                              htmlhelp("Rename a ROI","help_roimanual_renameroi"),
                              htmlhelp("Edit ROI notes","help_roimanual_editnotes"),
                              htmlhelp("Rename enrichments","help_roimanual_renameenrich"),
                              htmlhelp("Reorder ROIs","help_roimanual_reorderroi"),
                              htmlhelp("Reorder enrichments" ,"help_roimanual_reorderenrich")
                            ),choiceValues=list(
                              "overlaps",
                              "resize",
                              "pattern",
                              "summit",
                              "associateEnrichments",
                              "sample",
                              "filterWidth",
                              "filterEnrichment",
                              "renameRoi",
                              "editNotes",
                              "renameEnrichments",
                              "reorderRoi",
                              "reorderEnrichments"
                            ),selected=character(0))
            ),
            column(width=10,
              #different outputs if different choice 
              uiOutput("show_manualROImanipulationPage")
            )
          )


        )    
      })
		}



})