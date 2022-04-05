################################################################################
# FILE MANAGEMENT 
################################################################################


###react to help buttons:
#choose coordinate box
observeEvent(input$msg_coordinateFiles_chooseCoordinates, {
  boxHelpServer(msg_coordinateFiles_chooseCoordinates)
})
#preview window of coordinates
observeEvent(input$msg_coordinateFiles_filePreview, {
  boxHelpServer(msg_coordinateFiles_filePreview)
})
#available ROIs list
observeEvent(input$msg_coordinateFiles_roisAvailable, {
  boxHelpServer(msg_coordinateFiles_roisAvailable)
})



observeEvent(input$file, {
  temp=shinyFileName(input$file)

  #if the file is not NULL and not an integer (file yet not chosen)
  if (!is.null(temp) & !inherits(temp,"integer")   ){
    
    BEDvariables$completepath=shinyFilePath(input$file)
    BEDvariables$opened=TRUE
    BEDvariables$sfn=temp
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    #if file is not in the history of the BED files already opened
    if ( !BEDvariables$sfn %in% nomi){
      #set current BED file opened and its name, with our customed functions
      tryCatch({
          BEDvariables$tempBED=readBEDGTFF(BEDvariables$completepath,Header=input$readheader,Skip=input$skiplines)
          BEDvariables$tempBEDname=BEDvariables$sfn
       },
       warning = function( w ){

         BEDvariables$tempBED=NULL
         BEDvariables$tempBEDname=NULL
         BEDvariables$completepath=NULL
         #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file... try changin "Header" or "lines to skip"<br></font>'
        sendSweetAlert(
          session = session,
          title = "Cannot open file",
          text = "Make sure the file is in the correct format (BED or GTF/GFF) or try to change 'Header' or 'lines to skip' parameters",
          type = "error"
        )           
         BEDvariables$opened=FALSE
         BEDvariables$sfn=NULL
       },
       error = function( err ){

         BEDvariables$tempBED=NULL
         BEDvariables$tempBEDname=NULL
         BEDvariables$completepath=NULL
         #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file... try changin "Header" or "lines to skip"<br></font>'
        sendSweetAlert(
          session = session,
          title = "Cannot open file",
          text = "Make sure the file is in the correct format (BED or GTF/GFF) or try to change 'Header' or 'lines to skip' parameters",
          type = "error"
        )           
         BEDvariables$opened=FALSE
         BEDvariables$sfn=NULL
       })

    }else{
      BEDvariables$tempBED=NULL
      BEDvariables$tempBEDname=NULL
      BEDvariables$sfn=NULL
      BEDvariables$completepath=NULL
      #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Coordinates file already opened...<br></font>'
      sendSweetAlert(
        session = session,
        title = "File already opened",
        text = paste("File \'",temp,"\' had been already opened and is already in the list of ROIs",sep=""),
        type = "error"
      )        
      BEDvariables$opened=FALSE
    }
  }
},ignoreInit=TRUE)


#observe key pressed for importing the BED file from the link (input$BEDfrompath)
# toListenpath <- reactive({
#   list(input$confirmImportBEDfrompath)
# })
observeEvent(input$confirmImportBEDfrompath, {

  if(length(input$BEDfrompath)>0 ){

    if(input$BEDfrompath!=""){
      temp=basename(input$BEDfrompath)
      nomi=unlist(lapply(ROIvariables$listROI,getName))

      if(file.exists(input$BEDfrompath)){
        BEDvariables$completepath=input$BEDfrompath
        BEDvariables$sfn=temp
        if (!BEDvariables$sfn %in% nomi){
          BEDvariables$opened=TRUE
          tryCatch({
            BEDvariables$tempBED=readBEDGTFF(BEDvariables$completepath,Header=input$readheader,Skip=input$skiplines)
            BEDvariables$tempBEDname=BEDvariables$sfn
            
          },
          warning = function( w ){

            BEDvariables$tempBED=NULL
            BEDvariables$tempBEDname=NULL
            BEDvariables$completepath=NULL
            #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file... try changin "Header" or "lines to skip"<br></font>'
            sendSweetAlert(
              session = session,
              title = "Cannot open file",
              text = "Make sure the file is in the correct format (BED or GTF/GFF) or try to change 'Header' or 'lines to skip' parameters",
              type = "error"
            )  
            BEDvariables$opened=FALSE
            BEDvariables$sfn=NULL
            #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
          },
          error = function( err ){

            BEDvariables$tempBED=NULL
            BEDvariables$tempBEDname=NULL
            BEDvariables$completepath=NULL
            #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file... try changin "Header" or "lines to skip"<br></font>'
            sendSweetAlert(
              session = session,
              title = "Cannot open file",
              text = "Make sure the file is in the correct format (BED or GTF/GFF) or try to change 'Header' or 'lines to skip' parameters",
              type = "error"
            )            

            BEDvariables$opened=FALSE
            BEDvariables$sfn=NULL
            #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
          })


        }else{
          BEDvariables$tempBED=NULL
          BEDvariables$tempBEDname=NULL
          BEDvariables$sfn=NULL
          BEDvariables$completepath=NULL
          #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Coordinates file already opened...<br></font>'
          
          sendSweetAlert(
            session = session,
            title = "File already opened",
            text = paste("File \'",temp,"\' had been already opened and is already in the list of ROIs",sep=""),
            type = "error"
          )

          BEDvariables$opened=FALSE  
          #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")       
        }

      }else{
        BEDvariables$tempBED=NULL
        BEDvariables$tempBEDname=NULL
        BEDvariables$sfn=NULL
        BEDvariables$opened=FALSE
        BEDvariables$completepath=NULL
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">File selected \'',temp,'\' doesn t exist or wrong file path...<br></font>',sep="")
        #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
        sendSweetAlert(
          session = session,
          title = "File not found",
          text = paste("File selected \'",temp,"\' doesn't exist. Check if the path and name of the file are correct (also avoid spaces)",sep=""),
          type = "error"
        )
      }      
    }else{
      #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
      sendSweetAlert(
        session = session,
        title = "Empty file path",
        text = "TO open a coordinate file from path (either in BED or GTF/GFF format), plase typ the correct path to the fil in the field",
        type = "error"
      )
    }

  }else{
    #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
    sendSweetAlert(
      session = session,
      title = "Empty file path",
      text = "TO open a coordinate file from path (either in BED or GTF/GFF format), plase typ the correct path to the fil in the field",
      type = "error"
    )
  }
},ignoreInit=TRUE)












#READ bed file
#verify the confirmation of the BED file ("open it!" button). If ok, transform to GRanges and add to GRlist
observeEvent(input$confirmation, {
  
  #if confirm button, read GRanges and create the new ROI
  if (!is.null(BEDvariables$tempBEDname)){
    #check if current opened file is not NULL and have sufficient columns
    if (!is.null(BEDvariables$tempBED) & !ncol(BEDvariables$tempBED)<3){
      # it's all ok, so from the data.frame of the file, obtain GRanges object
      #do a trycatch, otherwise fatal error
      tryCatch({
        if(ncol(BEDvariables$tempBED)==3){
          range=GRanges(Rle(BEDvariables$tempBED[,1]),IRanges(BEDvariables$tempBED[,2],BEDvariables$tempBED[,3]))
        }else{
          range=GRanges(Rle(BEDvariables$tempBED[,1]),IRanges(BEDvariables$tempBED[,2],BEDvariables$tempBED[,3]),strand(BEDvariables$tempBED[,4]))
        }

        #here, force UCSC nomenclature (chr1, chr2,...chrM.... etc)
        range=convertNomenclatureGR(range,to="UCSC")

        ##this ROI is not a subset, but is a new ROI. Re-annotate if DB present (not FALSE, not null)
        if(length(DATABASEvariables$currentASSEMBLY)>0){
          if(DATABASEvariables$currentASSEMBLY!=FALSE){
            #we have a database. So extract the fix of promoters and apply the distanceFromTSS3 function
            #for this newly created range
            nomi=unlist(lapply(ROIvariables$listROI,getName))
            pos_promo=match("promoters",nomi)
            promo=ROIvariables$listROI[[pos_promo]]
            fix_promoters=getFixed(promo)
            annotatedpart=suppressWarnings(distanceFromTSS3(Object=range,Tss=fix_promoters,criterion="midpoint"))
            elementMetadata(range)=annotatedpart
          }
        }
        

        # log the action just done:
        logvariables$msg[[length(logvariables$msg)+1]]= paste("Added ",BEDvariables$tempBEDname," new ROI<br>",sep="")
        print(paste("Added",BEDvariables$tempBEDname,"new ROI"))
        ####newenrichimplementation####
        #new ROI imported has no list of enrichments at the beginning! => initialization
        Enrichlist$rawcoverage[[BEDvariables$tempBEDname]]=list()
        Enrichlist$normfactlist[[BEDvariables$tempBEDname]]=list()
        ################################
        #temporary BED and BED file name must be returned to NULL for the next open
        ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                      name=BEDvariables$tempBEDname,
                                      range=range,
                                      fixed=resize(range,width=1,fix="center"),
                                      flag="normalFlag",
                                      source=list(BEDvariables$tempBEDname)
                                      )  

        BEDvariables$tempBED=NULL
        BEDvariables$tempBEDname=NULL
        BEDvariables$opened=FALSE
        BEDvariables$sfn=NULL

        #alert the user that the file has been correctly import as ROI
        sendSweetAlert(
          session = session,
          title = "ROI imported!",
          text = paste("The coordinate file '",BEDvariables$tempBEDname,"' has been imported as a new ROI",sep=""),
          type = "success"
        )         
      },
      warning = function( w ){
        BEDvariables$sfn=NULL
        #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in reading the file (check the file format)...<br></font>'
        sendSweetAlert(
          session = session,
          title = "Unrecognised format",
          text = paste("Problems in reading the file (check the file format)",sep=""),
          type = "error"
        )            
        #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
      },
      error = function( err ){
        #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in reading the file (check the file format)...<br></font>'
        sendSweetAlert(
          session = session,
          title = "Unrecognised format",
          text = paste("Problems in reading the file (check the file format)",sep=""),
          type = "error"
        )            
        BEDvariables$sfn=NULL
        #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
      })

    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Coordinates file not added because columns <3 or data.frame==NULL...<br></font>'
      sendSweetAlert(
        session = session,
        title = "Bad format",
        text = paste("Coordinates file not added because n columns <3 (chr, start, end not found) or the dataframe is empty (NULL)",sep=""),
        type = "error"
      )           
      #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
      
    }
 
  }else{
    BEDvariables$sfn=NULL
    #updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")
    #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Coordinates file not added because name of the file is NULL...<br></font>'
    sendSweetAlert(
      session = session,
      title = "Bad file name",
      text = paste("Coordinates file not added because name of the file is NULL",sep=""),
      type = "error"
    )
  }
},ignoreInit=TRUE)


#cancel file preview:
observeEvent(input$cancellation, {
            BEDvariables$tempBED=NULL
            BEDvariables$tempBEDname=NULL
            BEDvariables$opened=FALSE
            BEDvariables$sfn=NULL
            updateTextInput(session, inputId="BEDfrompath", label = NULL, value = "")  
},ignoreInit=TRUE)































##########################################################
##########################################################
##########################################################
#GENELISTS
##########################################################
##########################################################
##########################################################






###react to help buttons:
#import genelists 
observeEvent(input$msg_genelists_importGenelist, {
  boxHelpServer(msg_genelists_importGenelist)
})
#parameters for genelists
observeEvent(input$msg_genelists_parametersGenelist, {
  boxHelpServer(msg_genelists_parametersGenelist)
})
#opened genelists
observeEvent(input$msg_genelists_openedGenelist, {
  boxHelpServer(msg_genelists_openedGenelist)
})



#open gene list file
observeEvent(input$fileGENELISTS, {
  genelistchoice=shinyFilePath(input$fileGENELISTS)

  #if the file is not NULL & is != 1 (sometimes is ==1 if button clicked, but not valid...)
  if (!is.null(genelistchoice) & file.exists(genelistchoice) ){
    toopen=basename(genelistchoice)


    #if file is not in the history of the BAM files already opened
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    regex=paste(c("promoters_genelist_","transcripts_genelist_","TES_genelist_"),toopen,sep="")
    pos=match(regex,nomi)
    #if there isn't something related to this genelist
    if(!any(!is.na(pos))){
      if("promoters"%in% nomi & "transcripts" %in% nomi & "TES" %in% nomi & length(DATABASEvariables$currentASSEMBLY)>0 ){

        #promoters, transcripts and TES are in the ROI list. 
        #match ID/symbols using the dictionary
        wherepromoters=match("promoters",nomi)
        roi=ROIvariables$listROI[[wherepromoters]]
        range=getRange(roi)

        oldSource_promoters=getSource(ROIvariables$listROI[[wherepromoters]])
        newfix_promoters=getFixed(ROIvariables$listROI[[wherepromoters]])
        
        wheretranscripts=match("transcripts",nomi)
        oldSource_transcripts=getSource(ROIvariables$listROI[[wheretranscripts]])
        newfix_transcripts=getFixed(ROIvariables$listROI[[wheretranscripts]])
        range_transcripts=getRange(ROIvariables$listROI[[wheretranscripts]])
        
        whereTES=match("TES",nomi)
        oldSource_TES=getSource(ROIvariables$listROI[[whereTES]])
        newfix_TES=getFixed(ROIvariables$listROI[[whereTES]])
        range_TES=getRange(ROIvariables$listROI[[whereTES]])
        #if symbolORid == IDs, use ids provided
        #else, use DATABASEvariables$currentDICTsymbol2id dictionary to convert into ID
        tryCatch({
          genelist=as.character(read.table(genelistchoice,header=FALSE)[,1])
          uniquegenelist=unique(genelist)
          lostinunique=length(genelist)-length(uniquegenelist)
        },
        warning = function( w ){
          #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file...<br></font>'
          sendSweetAlert(
            session = session,
            title = "Cannot open file",
            text = "File must be in a simple text format (txt, not formatted)",
            type = "error"
          )           
        },
        error = function( err ){
          #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file...<br></font>'
          sendSweetAlert(
            session = session,
            title = "Cannot open file",
            text = "File must be in a simple text format (txt, not formatted)",
            type = "error"
          )        
        })  

        #findPositionFromGene() function
        pos_match=findPositionFromGene(genelist=uniquegenelist,annotatedrange=range,kindofID=input$symbolORid,thresh=input$thresholdTranscripts)

        lostnotfound=pos_match$notfound
        losttoolarge=pos_match$losttoolarge
        lostisoformstoolarge=pos_match$lostisoformstoolarge


        if (length(pos_match$totake)>0 & !all(is.na(pos_match$totake))){
          
          newrange_promoters=range[pos_match$totake]
          newrange_transcripts=range_transcripts[pos_match$totake]
          newrange_TES=range_TES[pos_match$totake]

          newfix_promoters=newfix_promoters[pos_match$totake]
          newfix_transcripts=newfix_transcripts[pos_match$totake]
          newfix_TES=newfix_TES[pos_match$totake]
          

          #build the message
          mainmsg=paste("in ",toopen," genelist (",length(genelist)," original ",input$symbolORid,", ",lostinunique+length(lostnotfound)+length(losttoolarge)+length(lostisoformstoolarge)," lost ",sep="")
          if(lostinunique>0){
            mainmsg=paste(mainmsg,", ",lostinunique," non-unique",sep="")
          }
          if(length(lostnotfound)>0){
            mainmsg=paste(mainmsg,", ",paste(lostnotfound,collapse=";")," not found",sep="")
          }
          if(length(losttoolarge)>0){
            mainmsg=paste(mainmsg,", ",paste(losttoolarge,collapse=";")," has width > ",input$thresholdTranscripts,sep="")
          }
          if(length(lostisoformstoolarge)>0){
            mainmsg=paste(mainmsg,", some isoforms of ",paste(lostisoformstoolarge,collapse=";")," has width > ",input$thresholdTranscripts,sep="")
          }
          mainmsg=paste(mainmsg,")",sep="")

          newSource=paste("taken those ",mainmsg,sep="")
          newSource_promoters=c(oldSource_promoters,list(newSource))
          newSource_transcripts=c(oldSource_transcripts,list(newSource))
          newSource_TES=c(oldSource_TES,list(newSource))

          #create new promoters, transcripts and TES for this gene list
          ####newenrichimplementation####
          #new ROI imported has no list of enrichments at the beginning! => initialization
          Enrichlist$rawcoverage[[paste("promoters_genelist_",toopen,sep="")]]=list()
          Enrichlist$normfactlist[[paste("promoters_genelist_",toopen,sep="")]]=list()
          ################################
          ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                    name=paste("promoters_genelist_",toopen,sep=""),
                                                                    range=newrange_promoters,
                                                                    fixed=newfix_promoters,
                                                                    flag="promoterFlag",
                                                                    source=newSource_promoters) 
          ####newenrichimplementation####
          #new ROI imported has no list of enrichments at the beginning! => initialization
          Enrichlist$rawcoverage[[paste("transcripts_genelist_",toopen,sep="")]]=list()
          Enrichlist$normfactlist[[paste("transcripts_genelist_",toopen,sep="")]]=list()
          ################################
          ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                    name=paste("transcripts_genelist_",toopen,sep=""),
                                                                    range=newrange_transcripts,
                                                                    fixed=newfix_transcripts,
                                                                    flag="transcriptFlag",
                                                                    source=newSource_transcripts) 
          ####newenrichimplementation####
          #new ROI imported has no list of enrichments at the beginning! => initialization
          Enrichlist$rawcoverage[[paste("TES_genelist_",toopen,sep="")]]=list()
          Enrichlist$normfactlist[[paste("TES_genelist_",toopen,sep="")]]=list()
          ################################          
          ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                    name=paste("TES_genelist_",toopen,sep=""),
                                                                    range=newrange_TES,
                                                                    fixed=newfix_TES,
                                                                    flag="TESFlag",
                                                                    source=newSource_TES) 


          logvariables$msg[[length(logvariables$msg)+1]]= paste('Created promoters, transcripts and TES ',mainmsg,", ",length(newrange_promoters)," obtained<br>",sep="")  
          print(paste('Created promoters, transcripts and TES for ',toopen ,' gene list',sep=""))

        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Genes not valid: associated promoters/transcripts/TES not found, maybe you put wrong kind of identifiers...<br></font>'
          sendSweetAlert(
            session = session,
            title = "Cannot open file",
            text = "Annotated genomic elements not found for the provided genelist. Maybe you put wrong kind of identifiers (e.g. you have Symbols in the file, but you selected ENTREZ ID in the menu)",
            type = "error"
          )            
        }

      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">I need both promoters, transcripts and TES. Ask to the database...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Annotated elements not found",
          text = "Promoters, transcripts and TES of a specific genome assembly not found, but you need them: go to 'Database' section and choose a genome assembly",
          type = "error"
        )       
      }

    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Some ROI still exist from this \'',toopen,'\' genelist file. Remove them before...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Existing genelist found",
        text = paste('Found either promoters, transcripts or TES for \'',toopen,'\' genelist. Remove them to import the genelist with the same name',sep=""),
        type = "error"
      )     
    }

  }
})





#observe key pressed for importing the GENELIST from path
observeEvent(input$createGENELISTSfrompath, {

  if(length(input$GENELISTSfrompath)>0 ){

    if(input$GENELISTSfrompath!=""){
      toopen=basename(input$GENELISTSfrompath)
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      regex=paste(c("promoters_genelist_","transcripts_genelist_","TES_genelist_"),toopen,sep="")
      pos=match(regex,nomi)
      #if there isn't something related to this genelist
      if(!any(!is.na(pos))){
    

        if(file.exists(input$GENELISTSfrompath)){


          #same block as before, to create promoters, transcripts and TES from genelist
          if("promoters"%in% nomi & "transcripts" %in% nomi & "TES" %in% nomi & length(DATABASEvariables$currentASSEMBLY)>0 ){


            #promoters, transcripts and TES are in the ROI list. 
            #match ID/symbols using the dictionary
            wherepromoters=match("promoters",nomi)
            roi=ROIvariables$listROI[[wherepromoters]]
            range=getRange(roi)

            oldSource_promoters=getSource(ROIvariables$listROI[[wherepromoters]])
            newfix_promoters=getFixed(ROIvariables$listROI[[wherepromoters]])
            
            wheretranscripts=match("transcripts",nomi)
            oldSource_transcripts=getSource(ROIvariables$listROI[[wheretranscripts]])
            newfix_transcripts=getFixed(ROIvariables$listROI[[wheretranscripts]])
            range_transcripts=getRange(ROIvariables$listROI[[wheretranscripts]])
            
            whereTES=match("TES",nomi)
            oldSource_TES=getSource(ROIvariables$listROI[[whereTES]])
            newfix_TES=getFixed(ROIvariables$listROI[[whereTES]])
            range_TES=getRange(ROIvariables$listROI[[whereTES]])
            #if symbolORid == IDs, use ids provided
            #else, use DATABASEvariables$currentDICTsymbol2id dictionary to convert into ID
            tryCatch({

              genelist=as.character(read.table(input$GENELISTSfrompath,header=FALSE)[,1])
              uniquegenelist=unique(genelist)
              lostinunique=length(genelist)-length(uniquegenelist)
            },
            warning = function( w ){
              #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file...<br></font>'
              sendSweetAlert(
                session = session,
                title = "Cannot open file",
                text = "File must be in a simple text format (txt, not formatted)",
                type = "error"
              )               
            },
            error = function( err ){
              #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the file...<br></font>'
              sendSweetAlert(
                session = session,
                title = "Cannot open file",
                text = "File must be in a simple text format (txt, not formatted)",
                type = "error"
              )               
            })
            #findPositionFromGene() function
            pos_match=findPositionFromGene(genelist=uniquegenelist,annotatedrange=range,kindofID=input$symbolORid,thresh=input$thresholdTranscripts)

            lostnotfound=pos_match$notfound
            losttoolarge=pos_match$losttoolarge
            lostisoformstoolarge=pos_match$lostisoformstoolarge


            if (length(pos_match$totake)>0 & !all(is.na(pos_match$totake))){
              
              newrange_promoters=range[pos_match$totake]
              newrange_transcripts=range_transcripts[pos_match$totake]
              newrange_TES=range_TES[pos_match$totake]

              newfix_promoters=newfix_promoters[pos_match$totake]
              newfix_transcripts=newfix_transcripts[pos_match$totake]
              newfix_TES=newfix_TES[pos_match$totake]
              

              #build the message
              mainmsg=paste("in ",toopen," genelist (",length(genelist)," original ",input$symbolORid,", ",lostinunique+length(lostnotfound)+length(losttoolarge)+length(lostisoformstoolarge)," lost ",sep="")
              if(lostinunique>0){
                mainmsg=paste(mainmsg,", ",lostinunique," non-unique",sep="")
              }
              if(length(lostnotfound)>0){
                mainmsg=paste(mainmsg,", ",paste(lostnotfound,collapse=";")," not found",sep="")
              }
              if(length(losttoolarge)>0){
                mainmsg=paste(mainmsg,", ",paste(losttoolarge,collapse=";")," has width > ",input$thresholdTranscripts,sep="")
              }
              if(length(lostisoformstoolarge)>0){
                mainmsg=paste(mainmsg,", some isoforms of ",paste(lostisoformstoolarge,collapse=";")," has width > ",input$thresholdTranscripts,sep="")
              }
              mainmsg=paste(mainmsg,")",sep="")

              newSource=paste("taken those ",mainmsg,sep="")
              newSource_promoters=c(oldSource_promoters,list(newSource))
              newSource_transcripts=c(oldSource_transcripts,list(newSource))
              newSource_TES=c(oldSource_TES,list(newSource))

              #create new promoters, transcripts and TES for this gene list
              ####newenrichimplementation####
              #new ROI imported has no list of enrichments at the beginning! => initialization
              Enrichlist$rawcoverage[[paste("promoters_genelist_",toopen,sep="")]]=list()
              Enrichlist$normfactlist[[paste("promoters_genelist_",toopen,sep="")]]=list()
              ################################ 
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                        name=paste("promoters_genelist_",toopen,sep=""),
                                                                        range=newrange_promoters,
                                                                        fixed=newfix_promoters,
                                                                        flag="promoterFlag",
                                                                        source=newSource_promoters) 
              ####newenrichimplementation####
              #new ROI imported has no list of enrichments at the beginning! => initialization
              Enrichlist$rawcoverage[[paste("transcripts_genelist_",toopen,sep="")]]=list()
              Enrichlist$normfactlist[[paste("transcripts_genelist_",toopen,sep="")]]=list()
              ################################ 
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                        name=paste("transcripts_genelist_",toopen,sep=""),
                                                                        range=newrange_transcripts,
                                                                        fixed=newfix_transcripts,
                                                                        flag="transcriptFlag",
                                                                        source=newSource_transcripts) 
              ####newenrichimplementation####
              #new ROI imported has no list of enrichments at the beginning! => initialization
              Enrichlist$rawcoverage[[paste("TES_genelist_",toopen,sep="")]]=list()
              Enrichlist$normfactlist[[paste("TES_genelist_",toopen,sep="")]]=list()
              ################################               
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                        name=paste("TES_genelist_",toopen,sep=""),
                                                                        range=newrange_TES,
                                                                        fixed=newfix_TES,
                                                                        flag="TESFlag",
                                                                        source=newSource_TES) 

              logvariables$msg[[length(logvariables$msg)+1]]= paste('Created promoters, transcripts and TES ',mainmsg,", ",length(newrange_promoters)," obtained<br>",sep="")  
              print(paste('Created promoters, transcripts and TES for ',toopen ,' gene list',sep=""))
              updateTextInput(session, inputId="GENELISTSfrompath", label = NULL, value = "") 
            }else{
              #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Genes not valid: associated promoters/transcripts/TES not found, maybe you put wrong kind of identifiers...<br></font>'
              sendSweetAlert(
                session = session,
                title = "Cannot open file",
                text = "Annotated genomic elements not found for the provided genelist. Maybe you put wrong kind of identifiers (e.g. you have Symbols in the file, but you selected ENTREZ ID in the menu)",
                type = "error"
              )               
            }

     
          }else{
            #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">I need both promoters, transcripts and TES. Ask to the database...<br></font>',sep="")
            sendSweetAlert(
              session = session,
              title = "Annotated elements not found",
              text = "Promoters, transcripts and TES of a specific genome assembly not found, but you need them: go to 'Database' section and choose a genome assembly",
              type = "error"
            )             
          }

    
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">File selected \'',toopen,'\' doesn t exist or wrong file path...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "File not found",
            text = paste("File selected \'",toopen,"\' doesn't exist. Check if the path and name of the file are correct (also avoid spaces)",sep=""),
            type = "error"
          )          
        }   
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Some ROI still exist from this \'',toopen,'\' genelist file. Remove them before...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Existing genelist found",
          text = paste('Found either promoters, transcripts or TES for \'',toopen,'\' genelist. Remove them to import the genelist with the same name',sep=""),
          type = "error"
        )        
        updateTextInput(session, inputId="GENELISTSfrompath", label = NULL, value = "")        
      }
    }

  }
})







#observer for the text area pasted genes
observeEvent(input$createGENELISTSfrompaste,{
  
  if(input$pastedGENELISTS!=""){
    if(input$nameGENELISTS!=""){
      toopen=input$nameGENELISTS
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      regex=paste(c("promoters_genelist_","transcripts_genelist_","TES_genelist_"),toopen,sep="")
      pos=match(regex,nomi)
      #if there isn't something related to this genelist
      if(!any(!is.na(pos))){
        #if we have promoters and transcripts and TES taken from database
        if("promoters"%in% nomi & "transcripts" %in% nomi & "TES" %in% nomi & length(DATABASEvariables$currentASSEMBLY)>0 ){


          #find promoters and all their features
          wherepromoters=match("promoters",nomi)
          roi=ROIvariables$listROI[[wherepromoters]]
          range=getRange(roi)
          oldSource_promoters=getSource(ROIvariables$listROI[[wherepromoters]])
          newfix_promoters=getFixed(ROIvariables$listROI[[wherepromoters]])
          #find transcripts and all their features
          wheretranscripts=match("transcripts",nomi)
          oldSource_transcripts=getSource(ROIvariables$listROI[[wheretranscripts]])
          newfix_transcripts=getFixed(ROIvariables$listROI[[wheretranscripts]])
          range_transcripts=getRange(ROIvariables$listROI[[wheretranscripts]])
          #find TES and all their features
          whereTES=match("TES",nomi)
          oldSource_TES=getSource(ROIvariables$listROI[[whereTES]])
          newfix_TES=getFixed(ROIvariables$listROI[[whereTES]])
          range_TES=getRange(ROIvariables$listROI[[whereTES]])

          #read the pasted list, \n separated
          genelist=strsplit(input$pastedGENELISTS,split="\n")[[1]]
          uniquegenelist=unique(genelist)
          lostinunique=length(genelist)-length(uniquegenelist)
          #findPositionFromGene() function
          pos_match=findPositionFromGene(genelist=uniquegenelist,annotatedrange=range,kindofID=input$symbolORid,thresh=input$thresholdTranscripts)

          lostnotfound=pos_match$notfound
          losttoolarge=pos_match$losttoolarge
          lostisoformstoolarge=pos_match$lostisoformstoolarge


          if (length(pos_match$totake)>0 & !all(is.na(pos_match$totake))){
            
            newrange_promoters=range[pos_match$totake]
            newrange_transcripts=range_transcripts[pos_match$totake]
            newrange_TES=range_TES[pos_match$totake]

            newfix_promoters=newfix_promoters[pos_match$totake]
            newfix_transcripts=newfix_transcripts[pos_match$totake]
            newfix_TES=newfix_TES[pos_match$totake]
            

            #build the message
            mainmsg=paste("in ",toopen," genelist (",length(genelist)," original ",input$symbolORid,", ",lostinunique+length(lostnotfound)+length(losttoolarge)+length(lostisoformstoolarge)," lost ",sep="")
            if(lostinunique>0){
              mainmsg=paste(mainmsg,", ",lostinunique," non-unique",sep="")
            }
            if(length(lostnotfound)>0){
              mainmsg=paste(mainmsg,", ",paste(lostnotfound,collapse=";")," not found",sep="")
            }
            if(length(losttoolarge)>0){
              mainmsg=paste(mainmsg,", ",paste(losttoolarge,collapse=";")," has width > ",input$thresholdTranscripts,sep="")
            }
            if(length(lostisoformstoolarge)>0){
              mainmsg=paste(mainmsg,", some isoforms of ",paste(lostisoformstoolarge,collapse=";")," has width > ",input$thresholdTranscripts,sep="")
            }
            mainmsg=paste(mainmsg,")",sep="")

            newSource=paste("taken those ",mainmsg,sep="")
            newSource_promoters=c(oldSource_promoters,list(newSource))
            newSource_transcripts=c(oldSource_transcripts,list(newSource))
            newSource_TES=c(oldSource_TES,list(newSource))

            #create new promoters, transcripts and TES for this gene list
            ####newenrichimplementation####
            #new ROI imported has no list of enrichments at the beginning! => initialization
            Enrichlist$rawcoverage[[paste("promoters_genelist_",toopen,sep="")]]=list()
            Enrichlist$normfactlist[[paste("promoters_genelist_",toopen,sep="")]]=list()
            ################################ 
            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                      name=paste("promoters_genelist_",toopen,sep=""),
                                                                      range=newrange_promoters,
                                                                      fixed=newfix_promoters,
                                                                      flag="promoterFlag",
                                                                      source=newSource_promoters) 
            ####newenrichimplementation####
            #new ROI imported has no list of enrichments at the beginning! => initialization
            Enrichlist$rawcoverage[[paste("transcripts_genelist_",toopen,sep="")]]=list()
            Enrichlist$normfactlist[[paste("transcripts_genelist_",toopen,sep="")]]=list()
            ################################
            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                      name=paste("transcripts_genelist_",toopen,sep=""),
                                                                      range=newrange_transcripts,
                                                                      fixed=newfix_transcripts,
                                                                      flag="transcriptFlag",
                                                                      source=newSource_transcripts) 
            ####newenrichimplementation####
            #new ROI imported has no list of enrichments at the beginning! => initialization
            Enrichlist$rawcoverage[[paste("TES_genelist_",toopen,sep="")]]=list()
            Enrichlist$normfactlist[[paste("TES_genelist_",toopen,sep="")]]=list()
            ################################             
            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                                      name=paste("TES_genelist_",toopen,sep=""),
                                                                      range=newrange_TES,
                                                                      fixed=newfix_TES,
                                                                      flag="TESFlag",
                                                                      source=newSource_TES) 


            logvariables$msg[[length(logvariables$msg)+1]]= paste('Created promoters, transcripts and TES ',mainmsg,", ",length(newrange_promoters)," obtained<br>",sep="")  
            print(paste('Created promoters, transcripts and TES for ',toopen ,' gene list',sep=""))
            updateTextInput(session, inputId="nameGENELISTS", label = NULL, value = "")
            updateTextAreaInput(session,inputId="pastedGENELISTS",value="",label=NULL)
          }else{
            #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Genes not valid: associated promoters/transcripts/TES not found, maybe you put wrong kind of identifiers...<br></font>'
            sendSweetAlert(
              session = session,
              title = "Cannot open file",
              text = "Annotated genomic elements not found for the provided genelist. Maybe you put wrong kind of identifiers (e.g. you have Symbols in the file, but you selected ENTREZ ID in the menu)",
              type = "error"
            )            
          }


        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">I need both promoters, transcripts and TES. Ask to the database...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "Annotated elements not found",
            text = "Promoters, transcripts and TES of a specific genome assembly not found, but you need them: go to 'Database' section and choose a genome assembly",
            type = "error"
          )         
        }
      
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Some ROI still exist from this \'',toopen,'\' genelist. Remove them before, or choose another name...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Existing genelist found",
          text = paste('Found either promoters, transcripts or TES for \'',toopen,'\' genelist. Remove them to import the genelist with the same name',sep=""),
          type = "error"
        )       
      }
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">File name is empty...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Empty gene list name",
        text = "Put a name for the new gene list in 'genelist name' field",
        type = "error"
      )      
    }
  }else{
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Put some genes (IDs or symbols) in the text area...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Empty list",
      text = 'Put some genes (IDs or symbols) in the text area',
      type = "error"
    )   
  }
  #check if both teextarea AND filename are null or without length (pastedGENELISTS & nameGENELISTS)
})














##############################################################
# delete,rename,reorder ROIs
##############################################################





###react to help buttons:
#delete
observeEvent(input$msg_deleteRois_deleteRois, {
  boxHelpServer(msg_deleteRois_deleteRois)
})
#rename
observeEvent(input$msg_deleteRois_renameRois, {
  boxHelpServer(msg_deleteRois_renameRois)
})
#reorder
observeEvent(input$msg_deleteRois_reorderRois, {
  boxHelpServer(msg_deleteRois_reorderRois)
})

#remove ROI after pushing button:
observeEvent(input$deleteROI,{

  #check if at least one ROI is selected from the checkbutton list
  if (length(input$selectedCustomROItoRemove)>0){
    ROInametoRemove=input$selectedCustomROItoRemove
    #if "promoters","transcripts","TES" are being deleted,
    #make the DATABASEvariables$currentASSEMBLY to FALSE, because we cannot 
    #consider anymore the use of annotated elements
    #"promoters","transcripts","TES" must be deleted simoultaneously,
    #they are special ROIs
    thereisprom="promoters"%in%ROInametoRemove
    thereistransc="transcripts"%in%ROInametoRemove
    thereistes="TES"%in%ROInametoRemove
    anyofannotation=c(thereisprom,thereistransc,thereistes)
    if(any(anyofannotation) & !all(anyofannotation) ){
      #any is present, but not all together; cannot remove ROIs
      if(thereisprom){
        ROInametoRemove=ROInametoRemove[-which(ROInametoRemove=="promoters")]
      }
      if(thereistransc){
        ROInametoRemove=ROInametoRemove[-which(ROInametoRemove=="transcripts")]
      }
      if(thereistes){
        ROInametoRemove=ROInametoRemove[-which(ROInametoRemove=="TES")]
      }
    }else if(all(anyofannotation)){
      #all promoters, transcripts, TES present. We can remove all and set 
      #current DB assembly to FALSE, because we removed anything
      DATABASEvariables$currentASSEMBLY=FALSE
    }

    if(length(ROInametoRemove)>0){
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      #match with the names of the BEDs eventually modified
      pos=match(ROInametoRemove,nomi)

      ####newenrichimplementation####
      Enrichlist$rawcoverage[pos]<-NULL
      Enrichlist$normfactlist[pos]<-NULL
      ################################

      # #now update nomi, ROIvariables$listROI (the GR object)
      # for(i in pos){
      #   ROIvariables$listROI[[i]]=setBAMlist(ROIvariables$listROI[[i]],list())
      # }

      ROIvariables$listROI[pos]<-NULL

      for (i in ROInametoRemove){
        logvariables$msg[[length(logvariables$msg)+1]]= paste("Removed ROI ",i,"<br>",sep="")
        print(paste("Removed ROI",i))
      }

      #calling garbage collector can be very slow...
      print("forcing gc...")
      print(gc())   
   
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">You have to remove both promoters,transcripts,TES for the given genome assembly...<br></font>'
      sendSweetAlert(
        session = session,
        title = "Removing annotated elements",
        text = "You are trying to remove either 'promoters','transcripts' or 'TES' of a specific genome assembly. You have to remove all of them at once",
        type = "error"
      )    
    }
  }else{
    #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">No ROI selected to remove...<br></font>'
    sendSweetAlert(
      session = session,
      title = "No ROI selected",
      text = "You have to select at least one ROI to remove",
      type = "error"
    )  
  }
},ignoreInit=TRUE)


#rename ROI
observeEvent(input$renameROI,{
  #check if list ROI is not empty or NULL
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1){
  	#check if name is set
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$newfilenameROI)>=1 & !(any(input$ROInameResize == nottobe))) {
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$newfilenameROI) &!grepl(nottobe2[2],input$newfilenameROI) & !grepl(nottobe2[3],input$newfilenameROI)){
        #check if name does not match with another existing BED
        if (!input$newfilenameROI %in% nomi ){
          pos=match(input$selectedCustomROItoRename, nomi)
          oldname=nomi[pos]
          #also name in the object of the class
          ####newenrichimplementation####
          names(Enrichlist$rawcoverage)[pos]<-input$newfilenameROI
          names(Enrichlist$normfactlist)[pos]<-input$newfilenameROI
          ###############################
          ROIvariables$listROI[[pos]]=setName(ROIvariables$listROI[[pos]],input$newfilenameROI)
          #log the change:
          logvariables$msg[[length(logvariables$msg)+1]]= paste('Renamed ',oldname,' ROI in ',input$newfilenameROI,'<br></font>',sep="")
          print(paste('Renamed',oldname,'ROI in',input$newfilenameROI))
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newfilenameROI,' name already exist in ROIs...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "ROI already present",
            text = "A ROI with the same name already exists; use a different name",
            type = "error"
          )        
        }        
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newfilenameROI,' cannot use this kind of name ...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Bad ROI name",
          text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
          type = "error"
        )      
      }
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newfilenameROI,' is not a valid name...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Bad ROI name",
        text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
        type = "error"
      )    
    }    
  }

},ignoreInit=TRUE)




#reorder ROI
observeEvent(input$reorderROI,{
  #check new order, if all numbers of 1:length(GRvariables$names) are in input$reorderoption...
  if (!is.null(ROIvariables$listROI) & length(ROIvariables$listROI)>1){
  	nomi=unlist(lapply(ROIvariables$listROI,getName))
    allnumbers=as.character(1:length(nomi))
    #stringval=grep("reorderoptionROI",names(input),value=TRUE)
    listprovv=list()
    for (i in 1:length(ROIvariables$listROI)){
      stringval=grep(paste("reorderoptionROI",i,"$",sep=""),names(input),value=TRUE)
      listprovv[[i]]=input[[ stringval ]]
    }  
    listprovv=as.numeric(unlist(listprovv))
    #check if new order comprises all the possible positions
    if (identical(unique(sort(listprovv)),unique(sort(as.numeric(allnumbers)))) ){
      #ord=order(listprovv)
      ord=listprovv
      #reorder nomi (character)
      #nomi=nomi[ord]
      #reorder ROIvariables$listROI classes (list)

      ####newenrichimplementation####
      #in theory useless, but better to keep the order also in the list of enrichments
      Enrichlist$rawcoverage=Enrichlist$rawcoverage[order(ord)]
      Enrichlist$normfactlist=Enrichlist$normfactlist[order(ord)]
      ###############################
      ROIvariables$listROI=ROIvariables$listROI[order(ord)]
      if (!identical(unique(listprovv),unique(allnumbers))){
        logvariables$msg[[length(logvariables$msg)+1]]= paste('ROI files reordered...<br>',sep="")
        print("ROI files reordered")
      }  

    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red"> There arent\' all possible position in new ranking...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Ordering problem",
        text = "You didn't put all the ranking positions correct in the new ordering",
        type = "error"
      )    
    }    
  }

},ignoreInit=TRUE)









