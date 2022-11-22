
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

  #those 2 have been moved in uiBED, but kept here
  updateCheckboxGroupInput(session,inputId="selectedCustomROItoRemove",label=NULL,
                                  choices = historylist)  
  
      
  updateCheckboxGroupInput(session,"confirmviewROI", label=NULL,choices=historylist) 
  updateSelectInput(session,"listgetROI",NULL,choices=historylist)
 
  #show current updated ROI list also in coordinate files section
  output$showBEDfiles<-renderText({paste(historylist,collapse="<br>")})


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
              HTML("Select gene identifier to show:"),
              #show all IDs (excluding distancefromTSS)
              radioButtons("IDorsymbolwindow",label=NULL ,
                                       choices=c("gene_id"="gene_id",
                                                 "symbol"="symbol",
                                                 "ensembl_id"="ensembl_id",
                                                 "refSeq_id"="refSeq_id"),
                                       selected="gene_id"
              ),
              HTML("<br>"),
              htmlhelp("Select the genomic window (bp):","help_BED_getroi_windowvalue"),
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
        

      }#else{
      #   #here edit/download ROI notes
      #   #get value of source of the selected ROI:
      #   text=getSource(roi)
      #   windowNotesToShow<-list(textAreaInput(inputId="ROInotes",label="View and edit the notes of the selected ROI:",value=text,cols=60,rows=10),
      #                 actionButton("saveNotesROI", "save notes"),
      #                 downloadButton('downloadNotesROI', 'Download notes'))
      #   output$previewROItodownload<-renderUI({windowNotesToShow})
      #   output$previewROItodownloadbutton<-renderUI({NULL})
      #   listGUI<-list(paste(""))
      # }  







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



