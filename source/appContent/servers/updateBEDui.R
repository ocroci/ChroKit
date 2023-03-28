################################################################################
################################################################################
################################################################################
################################################################################
# update UI for coordinate files section
################################################################################
################################################################################
################################################################################
################################################################################
#show head selected files, between file choosing process and file confirmation/new GRange
#use rendereTable instead of renderDataTable according to your needs
output$fileHead <- renderDataTable({

  if (!is.null(BEDvariables$tempBED)& BEDvariables$opened){
    as.matrix(BEDvariables$tempBED)

  }else{
    NULL
    #as.matrix(data.frame("BED_file_now_not_selected_or_not_valid..."))
  }
  
},options=list(pagingType = "simple",pageLength = 5,lengthChange=FALSE,searching=FALSE))



#show current file opened
observe({
  BEDvariables$sfn
  if (!isvalid(BEDvariables$sfn)){
    output$showcurrentfile<-renderText({NULL})
    output$show_filepreviewtitle<-renderUI({NULL})
    return()
  }
  output$showcurrentfile<-renderText({
    paste("File content: <b>",BEDvariables$sfn,"</b>",sep="")
  })  

  output$show_filepreviewtitle<-renderUI({
    list(HTML("<h3>File preview"),htmlhelp("","help_BED_filepreview"),HTML("</h3>"))
  })

})






#observer for missing BSgenome
observe({
  ROIvariables$listROI
  if (!isvalid(ROIvariables$listROI)){
    output$showWarningBSgenome2<-renderUI({NULL})
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
    BSstring=paste("BSgenome.",org[pos],".UCSC.",asm,sep="")
    x=rownames(installed.packages())
    pos_pkg=match(BSstring,x)
    if(!is.na(pos_pkg)){
      #BSgenome found. Import library and do nothing
      library(BSstring,character.only=TRUE)
      output$showWarningBSgenome2<-renderUI({NULL})
    }else{
      #button for download the package
      output$showWarningBSgenome2<-renderUI({
        list(
        HTML(paste("<font color='red'>Need ",BSstring," package to be installed. Install it now?</font><br>",sep="")),
        actionButton("DownloadBSgenome2","Download")
        )
      })
    }

  }else{
    #warning message: I need database of a genome assembly
    output$showWarningBSgenome2<-renderUI({HTML("<font color='red'>Warning: choose a genome assembly from 'Databases' section</font>")})
    #here should nullify all other options!
  }

})




#react to main radiobutton (from where ROI?)
observe({
  if(!is.null(input$importROImainchoice)){
    if(input$importROImainchoice=="fromfile"){
      output$importROIwindowToShow<-renderUI({
        list(
          column(width=4,
            HTML("<h3>Open file:</h3><br>"),
            radioButtons("loadBEDsource",NULL,choices=c(
                                                      "Choose file from filesystem"="filesystem",
                                                      "Manually type the path of the file"="path"
                                                            ),selected="filesystem"),
            uiOutput("loadBEDsource"),
            HTML("<br>"),
            HTML('<hr size=3>'),
            HTML("<h4>Parameters:</h4>"),
            checkboxInput("readheader",label=list("Header",htmlhelp("","help_BED_headeroption")),value=TRUE),
            numericInput(inputId = 'skiplines',label=list("Lines to skip:",htmlhelp("","help_BED_linesskip")),min = 0, step = 1,value=0)
            
          ),
          column(width=8,
            uiOutput("show_filepreviewtitle"),
            HTML("<br>"),
            htmlOutput("showcurrentfile"),
            dataTableOutput("fileHead"),
            HTML("<br><br>"),
            #open button and cancel button
            fluidRow(
              column(4,uiOutput('cancelfilebutton')),
              column(2,uiOutput('openfilebutton'))
            )          
          )

        )

      })
    }else if(input$importROImainchoice=="fromgenelist"){

      #check existence of TXDB database
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      if("promoters"%in% nomi & "transcripts" %in% nomi & "TES" %in% nomi & length(DATABASEvariables$currentASSEMBLY)>0 ){

        output$importROIwindowToShow<-renderUI({
          list(
            column(width=6,  
              HTML("<h3>Genes to import</h3>"),

              radioButtons("loadGenelistsource",NULL,choices=c(
                                                        "Paste IDs/symbols"="paste",
                                                        "Choose gene list from filesystem"="filesystem",
                                                        "Manually type the path of the file"="path"
                                                              ),selected="paste"),
              uiOutput("loadGenelistsource")
            ),
            column(width=6,
              HTML("<h3>Parameters</h3>"),
              HTML("<br>"),
              radioButtons("symbolORid",label=list("What kind of identifiers are you importing?",htmlhelp("","help_BED_kindofID")),choices=c(
                                                        "ENTREZ IDs"="entrez",
                                                        "ENSEMBL IDs"="ensembl",
                                                        "Symbols"="symbol",
                                                        "RefSeq IDs"="refseq"
                                                              ),selected="symbol"),
              HTML("<br>"),
              list(HTML("<b>Max length for transcripts:</b>"),htmlhelp("","help_BED_maxtranscriptlen")),
              numericInput(inputId = 'thresholdTranscripts',label=NULL,min = 0, step = 100000,value=200000)       
            )
          )
        })
      }else{
        output$importROIwindowToShow<-renderUI({HTML("<font color='red'>You need to select a genome Assembly to import gene lists. Go to 'Assembly' to select the correct genome assembly.</font>")})
      }



    }else{
            #check existence of TXDB database
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      if("promoters"%in% nomi & "transcripts" %in% nomi & "TES" %in% nomi & length(DATABASEvariables$currentASSEMBLY)>0 ){

        output$importROIwindowToShow<-renderUI({
          list(
            #warning in case BSgenome DB not present (copy of what seen in extract pattern from modifyROI)
            uiOutput("showWarningBSgenome2"),
            #motif text input (copy of that in motifyROI)
            textInput("PatternToSearch2",label=list("Select pattern (IUPAC nomenclature)",htmlhelp("","help_BED_IUPACpattern")),placeholder="ATCNYGG"),
            #new ROI name
            textInput("ROInamePattern2",label="Name of the ROI",placeholder="type new ROI name here"),
            actionButton("ExtractPatternROI2","Create ROI") 
          )
        })
      }else{
        output$importROIwindowToShow<-renderUI({HTML("<font color='red'>You need to select a genome Assembly to extract sequence patterns. Go to 'Assembly' to select the correct genome assembly.</font>")})
      }

    }
  }else{
    output$importROIwindowToShow<-renderUI({NULL})
  }
})




#react to radiobutton, whether to choose from filesystem the coordinate file
observe({
  if (!is.null(input$loadBEDsource)){
    if(input$loadBEDsource=="filesystem"){
      output$loadBEDsource<-renderUI({
        shinyFilesButton('file', 'Select a file', 'Please select a file', FALSE)
      })
    }else{
      output$loadBEDsource<-renderUI({
        list(
          textInput("BEDfrompath",NULL,value="",placeholder = "/path/to/BEDorGTF"),
          actionButton("confirmImportBEDfrompath", "Open file")
        )
      })
    }
  }else{
    output$loadBEDsource<-renderUI({NULL})
  }
})


#react to radiobutton, whether to choose from filesystem or other source the genelist file
observe({
  if (!is.null(input$loadGenelistsource)){

    if(input$loadGenelistsource=="paste"){
      output$loadGenelistsource<-renderUI({
        list(
          HTML("<b>Paste IDs/symbols here:</b><br>"),
          textAreaInput("pastedGENELISTS",NULL,value="",height=150),
          textInput("nameGENELISTS",NULL,placeholder="new genelist name",value=""),
          actionButton("createGENELISTSfrompaste", "Import")  
        )
      })
    }else if (input$loadGenelistsource=="filesystem"){
      output$loadGenelistsource<-renderUI({
        shinyFilesButton('fileGENELISTS', 'Open gene list text file', 'Please select a txt file', FALSE)
      })
    }else{
      output$loadGenelistsource<-renderUI({
        list(
          textInput("GENELISTSfrompath",NULL,value=NULL,placeholder = "/path/to/geneList.txt"),
          actionButton("createGENELISTSfrompath", "Open gene list")
        )
      })
    }


  }else{
    output$loadGenelistsource<-renderUI({NULL})
  }

})




#create/hide "Open it!" button and cancel button for opening coord. file, only when file is opened
#or changed but always valid and opened
observe({
  if(!is.null(BEDvariables$tempBED) & !is.null(BEDvariables$tempBEDname)){
    output$openfilebutton<-renderUI({
      actionButton("confirmation", "Confirm and import as ROI")
    })
    output$cancelfilebutton<-renderUI({
      actionButton("cancellation", "Cancel")
    })      
  }else{
    output$openfilebutton<-renderText({
      c("")
    })
    output$cancelfilebutton<-renderUI({
      c("")
    })
  } 
})



##react to lines to skip (temp already there) and header (TRUE/FALSE), keeping the same file name
##until cancel
observe({
  input$readheader
  input$skiplines
  #now check: if file selected (sfn and tempbed must not be null or false), re-read the BED/GTF
  #otherwise skip (withOUT error popup).
  if(!is.null(isolate(BEDvariables$tempBED)) & !is.null(isolate(BEDvariables$sfn)) & !is.null(BEDvariables$completepath) & isvalid(input$skiplines)){
    tryCatch({
      BEDvariables$tempBED=readBEDGTFF(BEDvariables$completepath,Header=input$readheader,Skip=input$skiplines)
    },warning = function( w ){
      #if read file is not success, do not throw the file! simply keep the same as before with popup!
      sendSweetAlert(
        session = session,
        title = "Error in opening file",
        text = "Try to change 'Header' or 'lines to skip' parameters",
        type = "error"
      ) 
    },error = function( err ){
      sendSweetAlert(
        session = session,
        title = "Error in opening file",
        text = "Try to change 'Header' or 'lines to skip' parameters",
        type = "error"
      )       
    })
    #this way, if new header/lines to skip are not ok, keep the previous that were valid

  }else{

  }
})









  ######################################################################################
  #VIEW ROI statistics
  ######################################################################################
  ######################################################################################
  ######################################################################################

#observe for get roi and view roi boxes. NULL if no ROI
observe({
  ROIvariables$listROI
  if (!isvalid(ROIvariables$listROI)){
    removeTab(inputId = "newROItabset",target="managingROItabPanel")
    removeTab(inputId = "newROItabset",target="downloadROItabPanel")
    return()    
  }


  appendTab(inputId = "newROItabset",
    tabPanel("Managing ROIs",value="managingROItabPanel",

      fluidRow(
        column(width=9,style='padding:0px;',
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_quickviewROIs",title="Quick ROI preview"),
            fluidRow(
              column(width=4,
                uiOutput("show_confirmviewROI")
                # HTML("<b>Select ROI to view:</b>"),
                # wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 300px; max-width: 300px; background-color: #ffffff;",
                #   checkboxGroupInput("confirmviewROI", label=NULL,choices=NULL)
                # )
              ),

              column(width=8,
                uiOutput("show_chooseROIvisualiz"),
                plotOutput('viewROImaterial'),
                htmlOutput("saveviewpeaksROImaterial")


              )
            )
          )
        ),
        column(width=3,style='padding:0px;',
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_deleteRois_deleteRois",title="Loaded ROIs"),

            
            HTML("<b>Available ROIs:</b>"),
            uiOutput("show_selectedCustomROItoRemove"),

            HTML("Select ROIs to delete<br><br>"),
            actionButton("deleteROI", "Delete")
          ) 
        )
      )
    )
  )



  appendTab(inputId = "newROItabset",
    tabPanel("Download ROIs",value="downloadROItabPanel",
      fluidRow(
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_getRois_BOX",title="Download ROI"),
          fluidRow(
            column(width=4,
              uiOutput("show_choosegetROImenu"),
              uiOutput("showROIoptionsToGET")    
            ),
            column(width=8, 
              htmlOutput("previewROItodownload"),
              htmlOutput("previewROItodownloadbutton")
            )
          )
        )
      )

    )
  )



})






























#observer for options to view (width distribtion or number of ranges)
observe({
  input$confirmviewROI
  if (!isvalid(input$confirmviewROI)){
    output$show_chooseROIvisualiz<-renderUI({NULL})
    return()
  }

  output$show_chooseROIvisualiz<-renderUI({ radioButtons("chooseROIvisualiz",label=list(HTML("Select how to view ROIs info"),htmlhelp("","help_BED_viewoptions")),
                                choiceNames=list(  
                                  "Distribution of ranges width",
                                  "Number of ranges"
                                ),
                                choiceValues=list(
                                  "ROIwidth",
                                  "ROIintervals"
                                ),selected="ROIwidth",inline=TRUE)
  })

})



######################################################################################
#GET ROI box
######################################################################################






  # observe({
    
  #       input$confirmviewROI
  #       output$viewROIstat<-renderText(
  #         if(length(input$confirmviewROI)==1){
  #           nomi=unlist(lapply(ROIvariables$listROI,getName))
  #           pos=match(input$confirmviewROI,nomi)
  #           roi=ROIvariables$listROI[[pos]]
  #           if(!is.null(roi)){
  #             selectedrange=getRange(roi)
            
  #             wdth=width(selectedrange)
  #             x=quantile(wdth,input$quantileROIwidth)
  #             y=length(selectedrange[wdth>x])
  #             paste("Width at that quantile: <b>",round(x,0),"</b><br>Number peaks with greater width: <b>",y,"</b>",sep="")
  #           }else{
  #             paste("You have selected a non existent ROI...")
  #           }

  #         }else{
  #           paste("You have to select ONE ROI...")
  #         }
  #       )

  # })