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

# #show BED file name history
# output$showBEDfiles<-renderText({
#   paste(as.list(BEDvariables$BEDfilehistory),collapse="<br>")
# })


# #show current problem or not in opening the file, if any over file preview
# output$showproblem<-renderText({
#   paste(logvariables$currentproblem,collapse="\n")
# })

#show current file opened
output$showcurrentfile<-renderText({
  paste("File content: <b>",BEDvariables$sfn,"</b>",sep="")
})

# #update checkbox input of the files to be deleted, if BEDfilehistory changes
# observe({
#     historylist=as.list(GRvariables$names)
#     names(historylist)=GRvariables$names
#     updateCheckboxGroupInput(session,inputId="selectedBEDtoRemove",label="Coordinates to delete:",
#                                 choices = historylist)    
# })





#create/hide "Open it!" button and cancel button for opening coord. file, only when file is opened
#or changed but always valid and opened
observe({
  if(!is.null(BEDvariables$tempBED) & !is.null(BEDvariables$tempBEDname)){
    output$openfilebutton<-renderUI({
      actionButton("confirmation", "Import as ROI")
    })
    output$cancelfilebutton<-renderUI({
      actionButton("cancellation", "Cancel")
    })      
  }else{
    output$openfilebutton<-renderText({
      c("select a file...")
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
######################################################################################
######################################################################################
#RENAME ROI
######################################################################################
######################################################################################
######################################################################################


#update text field for renaming the ROI, if ROIvariables$names changed
observe({
  
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  updateTextInput(session,inputId="newfilenameROI",label="New ROI name:",value="")    
})

######################################################################################
######################################################################################
######################################################################################
#reorder ROI
######################################################################################
######################################################################################
######################################################################################

#reorder ROI
observeEvent(ROIvariables$listROI,{
  
  output$dinamicROI<-renderUI({
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
      return(lista)       
    }else{
      return(HTML("No ROI..."))
    }

  })   
})
