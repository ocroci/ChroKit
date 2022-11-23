#react to radiobutton, whether to choose from filesystem the enrichment file
observe({
  if (!is.null(input$loadEnrichmentsource)){
    if(input$loadEnrichmentsource=="filesystem"){
      output$loadEnrichmentsource<-renderUI({
        shinyFilesButton('fileBAM', 'Select a file', 'Please select an enrichment file', FALSE,multiple=TRUE)
      })
    }else{
      output$loadEnrichmentsource<-renderUI({
        list(
          textInput("BAMfrompath",NULL,value=NULL,placeholder = "/path/to/BAM.bam or bigWig.bw"),
          actionButton("confirmImportBAMfrompath", "Open file")
        )
      })
    }
  }else{
    output$loadEnrichmentsource<-renderUI({NULL})
  }
})





#show BAM files to delete
 observe({
 	if (length(BAMvariables$listBAM)>0){
    historylist=unlist(BAMvariables$listBAM)
    bampaths=names(BAMvariables$listBAM)
    finalname=paste(bampaths,"(",historylist,")")
    historylist=as.list(names(BAMvariables$listBAM))
    names(historylist)=finalname
 	}else{
 		historylist=list()
 	} 

  output$show_selectedBAMtoRemove<-renderUI({
    wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 450px; background-color: #ffffff;",
      checkboxGroupInput(inputId="selectedBAMtoRemove",label="",
                                choices = historylist) 
    )
  })


})


#show BAM files to rename. Show only the name of list, not the basename of the file
observe({
  BAMvariables$listBAM
  if(length(BAMvariables$listBAM)>0){
    historylist=as.list(names(BAMvariables$listBAM))
    names(historylist)=names(BAMvariables$listBAM)
  }else{
    historylist=list()
  } 

  output$show_selectedBAMfiletoRename<-renderUI({
    selectInput(inputId="selectedBAMfiletoRename",label="Enrichment to rename:",
                                choices = historylist) 
  })

})



#observer for boxes
observe({
  BAMvariables$listBAM
  if(!isvalid(BAMvariables$listBAM)){
    output$show_boxremovebam<-renderUI({NULL})
    output$show_boxrenamebam<-renderUI({NULL})
    return()
  }

  output$show_boxremovebam<-renderUI({
    box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
      title=boxHelp(ID="msg_enrichmentFiles_deleteEnrichment",title="Opened enrichment files"),
      HTML("<b>Available enrichment files:</b>"),
      uiOutput("show_selectedBAMtoRemove"),
      HTML("Select file references to delete<br>"),
      HTML("<br><br>"),
      actionButton("deleteBAM", "Delete")
    )
  })

  output$show_boxrenamebam<-renderUI({
    box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
      title=boxHelp(ID="msg_enrichmentFiles_renameEnrichment",title="Rename enrichment file"),

      uiOutput("show_selectedBAMfiletoRename"),
      
      HTML("<br><br>"),
      textInput("newfilenameBAMfile","New enrichment name:",placeholder="type new enrichment name here",value=""),
      HTML("<br><br>"),
      actionButton("renameBAMfile", "Rename")
    )     
  })

})
















