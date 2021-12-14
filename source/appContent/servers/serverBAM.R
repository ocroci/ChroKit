
###react to help buttons:
#import enrichments
observeEvent(input$msg_enrichmentFiles_importEnrichment, {
  boxHelpServer(msg_enrichmentFiles_importEnrichment)
})
#delete enrichments
observeEvent(input$msg_enrichmentFiles_deleteEnrichment, {
  boxHelpServer(msg_enrichmentFiles_deleteEnrichment)
})
#rename enrichments
observeEvent(input$msg_enrichmentFiles_renameEnrichment, {
  boxHelpServer(msg_enrichmentFiles_renameEnrichment)
})




#open BAM file paths
observeEvent(input$fileBAM, {
  bamchoice=shinyFilePath(input$fileBAM)
  #if the file is not NULL 
  if (!is.null(bamchoice)& is.na(suppressWarnings(as.integer(bamchoice)))){
    toopen=basename(bamchoice)
    if (length(BAMvariables$listBAM)>0){
      alreadyopened=sapply(BAMvariables$listBAM,basename)     
    }else{
      alreadyopened=NULL
    }
    #file must end with ".bam"
    #if file is not in the history of the BAM files already opened
    if (!toopen %in% alreadyopened){

      #if is a BAM file (extension ends with .bam)
      if ( substring(toopen,nchar(toopen)-3,nchar(toopen))==".bam"){
        #find if bam file has its own index (.bam.bai)
        baifile=paste(bamchoice,".bai",sep="")
        if (file.exists(baifile)){ 
          #add bam file path to the list!
          BAMvariables$listBAM[length(BAMvariables$listBAM)+1]<-bamchoice
          names(BAMvariables$listBAM)[length(BAMvariables$listBAM)]<-toopen
          logvariables$msg[[length(logvariables$msg)+1]]= paste("Added ",toopen," BAM file <br>",sep="")
          print(paste("Added",toopen,"BAM file"))
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">BAM doesn\'t have its own index (.bai). You can create it with samtools index <bam file>...<br></font>'
          sendSweetAlert(
            session = session,
            title = "bam index (.bai) not found",
            text = "BAM doesn\'t have its own index (.bai) in the same directory. You can create it with samtools index <bam file>",
            type = "error"
          )         
        }
      }

      #if is a bigWig file (extension ends with .bw/.bigWig)
      else if (substring(toopen,nchar(toopen)-2,nchar(toopen))==".bw" | tolower(substring(toopen,nchar(toopen)-6,nchar(toopen)))==".bigwig"){
        if (.Platform$OS.type=="unix"){
          BAMvariables$listBAM[length(BAMvariables$listBAM)+1]<-bamchoice
          names(BAMvariables$listBAM)[length(BAMvariables$listBAM)]<-toopen
          logvariables$msg[[length(logvariables$msg)+1]]= paste("Added ",toopen," WIG file <br>",sep="")
          print(paste("Added",toopen,"WIG file"))          
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Cannot use bigWig on Windows OS...<br></font>'
          sendSweetAlert(
            session = session,
            title = "Windows OS detected",
            text = "Cannot use bigWig enrichments on Windows OS, use bam files instead",
            type = "error"
          )  
        }
      #else, not bam neither bigWig
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">File selected doesn\'t end with \'.bam\' or \'.bw/.bigWig\' extension...<br></font>'
        sendSweetAlert(
          session = session,
          title = "Bad file format",
          text = "Selected file name doesn\'t end with \'.bam\' or \'.bw/.bigWig\' extension",
          type = "error"
        )       
      }
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">File selected \'',toopen,'\' already opened...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "File already opened",
        text = paste('File selected \'',toopen,'\' already opened',sep=""),
        type = "error"
      )     
    }
  }
})



#observe key pressed for importing the BAM file from the link (input$BAMfrompath)
observeEvent(input$confirmImportBAMfrompath, {

  if(length(input$BAMfrompath)>0 ){
    if(input$BAMfrompath!=""){
      toopen=basename(input$BAMfrompath)
      if (length(BAMvariables$listBAM)>0){
        alreadyopened=sapply(BAMvariables$listBAM,basename)       
      }else{
        alreadyopened=NULL
      }

      if(file.exists(input$BAMfrompath)){
        if (!toopen %in% alreadyopened){

          #if is a BAM file (extension ends with .bam)
          if (substring(toopen,nchar(toopen)-3,nchar(toopen))==".bam"){
            baifile=paste(input$BAMfrompath,".bai",sep="")
            if(file.exists(baifile)){
              BAMvariables$listBAM[length(BAMvariables$listBAM)+1]<-input$BAMfrompath
              names(BAMvariables$listBAM)[length(BAMvariables$listBAM)]<-toopen
              logvariables$msg[[length(logvariables$msg)+1]]= paste("Added ",toopen," BAM file <br>",sep="")
              print(paste("Added",toopen,"BAM file"))
              #clean the input field, rady for a new BAM file path
              #updateTextInput(session, inputId="BAMfrompath", label = NULL, value = "")
            }else{
               #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">BAM doesn\'t have its own index (.bai). You can create it with samtools index <bam file>...<br></font>'
              sendSweetAlert(
                session = session,
                title = "bam index (.bai) not found",
                text = "BAM doesn\'t have its own index (.bai) in the same directory. You can create it with samtools index <bam file>",
                type = "error"
              )                 
               #updateTextInput(session, inputId="BAMfrompath", label = NULL, value = "")
            }
          }

          #if is a bigWig file (extension ends with .bw/.bigWig)
          else if(substring(toopen,nchar(toopen)-2,nchar(toopen))==".bw" | tolower(substring(toopen,nchar(toopen)-6,nchar(toopen)))==".bigwig"){
            if (.Platform$OS.type=="unix"){
              BAMvariables$listBAM[length(BAMvariables$listBAM)+1]<-input$BAMfrompath
              names(BAMvariables$listBAM)[length(BAMvariables$listBAM)]<-toopen
              logvariables$msg[[length(logvariables$msg)+1]]= paste("Added ",toopen," WIG file <br>",sep="")
              print(paste("Added",toopen,"WIG file"))
            }else{
              #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Cannot use bigWig on Windows OS...<br></font>'
              sendSweetAlert(
                session = session,
                title = "Windows OS detected",
                text = "Cannot use bigWig enrichments on Windows OS, use bam files instead",
                type = "error"
              )            
            }


          }

          #else, not bam neither WIG file (extension)
          else{
            #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">File selected doesn\'t end with \'.bam\' or \'.bw/.bigWig\' extension...<br></font>'
            sendSweetAlert(
              session = session,
              title = "Bad file format",
              text = "Selected file name doesn\'t end with \'.bam\' or \'.bw/.bigWig\' extension",
              type = "error"
            )             
            #updateTextInput(session, inputId="BAMfrompath", label = NULL, value = "")
          }
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">File selected \'',toopen,'\' already opened...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "File already opened",
            text = paste('File selected \'',toopen,'\' already opened',sep=""),
            type = "error"
          )           
          #updateTextInput(session, inputId="BAMfrompath", label = NULL, value = "")
        }      
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">File selected \'',toopen,'\' doesn t exist or wrong file path...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "File not found",
          text = paste("File selected \'",toopen,"\' doesn't exist. Check if the path and name of the file are correct (also avoid spaces)",sep=""),
          type = "error"
        )        
        #updateTextInput(session, inputId="BAMfrompath", label = NULL, value = "")
      }      
    }

  }
})




#delete BAM file paths
observeEvent(input$deleteBAM,{
  #check if at least one BAM is selected from the checkbutton list
  if (length(input$selectedBAMtoRemove)>0){
    #match with the names of the BAMs eventually modified
    pos=match(input$selectedBAMtoRemove,names(BAMvariables$listBAM))
    #now update GRvariables$names (renamed names of GR/BAMs), GRvariables$listGR (the GR object), BAMvariables$BAMfilehistory (the history of BAM files opened)
    BAMvariables$listBAM= BAMvariables$listBAM[-c(pos)]
    for (i in input$selectedBAMtoRemove){
      logvariables$msg[[length(logvariables$msg)+1]]= paste("Removed ",i," enrichment file<br>",sep="")
      print(paste("Removed",i,"enrichment file"))
    }
  }else{
    #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">No enrichment files selected to remove...<br></font>'
    sendSweetAlert(
      session = session,
      title = "No enrichment selected",
      text = "Select at least a enrichment file to remove",
      type = "error"
    )
    
  }
})






#rename BAM files
observeEvent(input$renameBAMfile,{
  #check if one BAM file is selected
  if (nchar(input$selectedBAMfiletoRename)>=1) {
    #find position in BAM list
    totalbams=names(BAMvariables$listBAM)  
    pos=match(input$selectedBAMfiletoRename,totalbams)
    #check if name does not match with another existing BAM
    if (!input$newfilenameBAMfile %in% totalbams & length(pos)>0 ){
      oldname=totalbams[pos]
      names(BAMvariables$listBAM)[pos]=input$newfilenameBAMfile
      logvariables$msg[[length(logvariables$msg)+1]]= paste('Renamed ',oldname,' enrichment file in ',input$newfilenameBAMfile,'<br></font>',sep="")
      print(paste("Renamed",oldname,"enrichment file in",input$newfilenameBAMfile))
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newfilenameBAMfile,' name already exist in enrichment files list...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "File name exists",
        text = paste("\'",input$newfilenameBAMfile,"\' name already exists in the list of enrichment files",sep=""),
        type = "error"
      )     
    }
  }

})


