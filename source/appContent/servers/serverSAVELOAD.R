


observeEvent(input$saveWork,{
	#for SAVE environment
	volumes=c("wd"=rootsavedir)
  	
	fileinfo <- parseSavePath(volumes, input$saveWork)
	#fileinfo$datapath is the complete path + file name
	if (nrow(fileinfo) > 0) {
		#put all persisten variables together
		ROIvariables_listROI=ROIvariables$listROI
		ROIvariables_selected=ROIvariables$selected
		ROIvariables_selectedname=ROIvariables$selectedname
		ROIvariables_listfordensity=ROIvariables$listfordensity
		ROIvariables_colorsfordensity=ROIvariables$colorsfordensity
		ROIvariables_listselected=ROIvariables$listselected
		ROIvariables_listselectednames=ROIvariables$listselectednames
		ROIvariables_changed=ROIvariables$changed


		BEDvariables_tempBED=BEDvariables$tempBED
		BEDvariables_tempBEDname=BEDvariables$tempBEDname
		BEDvariables_opened=BEDvariables$opened
		BEDvariables_sfn=BEDvariables$sfn


		DATABASEvariables_currentASSEMBLY=DATABASEvariables$currentASSEMBLY
		#DATABASEvariables_availASSEMBLIES=DATABASEvariables$availASSEMBLIES
		#DATABASEvariables_missingASSEMBLIES=DATABASEvariables$missingASSEMBLIES


		BAMvariables_listBAM=BAMvariables$listBAM
		
		logvariables_msg=logvariables$msg


		listOfPeakTimeVariables=list(
			ROIvariables_listROI,
			ROIvariables_selected,
			ROIvariables_selectedname,
			ROIvariables_listfordensity,
			ROIvariables_colorsfordensity,
			ROIvariables_listselected,
			ROIvariables_listselectednames,
			ROIvariables_changed,
			BEDvariables_tempBED,
			BEDvariables_tempBEDname,
			BEDvariables_opened,
			BEDvariables_sfn,
			DATABASEvariables_currentASSEMBLY,
			#DATABASEvariables_availASSEMBLIES,
			#DATABASEvariables_missingASSEMBLIES,
			BAMvariables_listBAM,
			logvariables_msg
		)


		names(listOfPeakTimeVariables)=c(

			"ROIvariables_listROI",
			"ROIvariables_selected",
			"ROIvariables_selectedname",
			"ROIvariables_listfordensity",
			"ROIvariables_colorsfordensity",
			"ROIvariables_listselected",
			"ROIvariables_listselectednames",
			"ROIvariables_changed",
			"BEDvariables_tempBED",
			"BEDvariables_tempBEDname",
			"BEDvariables_opened",
			"BEDvariables_sfn",
			"DATABASEvariables_currentASSEMBLY",
			#"DATABASEvariables_availASSEMBLIES",
			#"DATABASEvariables_missingASSEMBLIES",
			"BAMvariables_listBAM",
			"logvariables_msg")
      	#logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="blue">Wainting for the session to be saved...<br></font>',sep="")
      	
        tryCatch({
	      	saveRDS(listOfPeakTimeVariables, as.character(fileinfo$datapath) )
	      	logvariables$msg[[length(logvariables$msg)+1]]= paste('Session ',basename(as.character(fileinfo$datapath)),' saved<br>',sep="")
	      	print(paste('Session ',as.character(fileinfo$datapath),' saved',sep=""))
	        sendSweetAlert(
	          session = session,
	          title = "Session saved!",
	          text = paste("The session '",basename(as.character(fileinfo$datapath)),"' has been saved",sep=""),
	          type = "success"
	        ) 
        
        },
        warning = function( w ){
	        sendSweetAlert(
	          session = session,
	          title = "Cannot save the session",
	          text = "Cannot save the session in the specified directory. Maybe you don't have the write permission.",
	          type = "error"
	        ) 
        },
        error = function( err ){
	        sendSweetAlert(
	          session = session,
	          title = "Cannot save the session",
	          text = "Cannot save the session in the specified directory. Maybe you don't have the write permission.",
	          type = "error"
	        ) 
        })
      
    }
},ignoreInit=TRUE)










#load rds file
observeEvent(input$loadenv, {
    sessiontoload=shinyFileName(input$loadenv)
    #if the file is not NULL 

    if (!is.null(sessiontoload) & length(sessiontoload)>0 & isvalid(sessiontoload) & class(sessiontoload)!="integer"){
        #try to open RDS file session
        tryCatch({

        	#purge the memory to import the new session
        	print("forcing gc for new session...")
    		print(gc())

    	
            totalENV=readRDS(shinyFilePath(input$loadenv))
            #load the ROIvariables
           	#x=totalENV$ROIvariables
            ROIvariables$listROI=totalENV$ROIvariables_listROI
            ROIvariables$selected=totalENV$ROIvariables_selected
            ROIvariables$selectedname=totalENV$ROIvariables_selectedname
            ROIvariables$listfordensity=totalENV$ROIvariables_listfordensity
            ROIvariables$colorsfordensity=totalENV$ROIvariables_colorsfordensity
            ROIvariables$listselected=totalENV$ROIvariables_listselected
            ROIvariables$listselectednames=totalENV$ROIvariables_listselectednames
            ROIvariables$changed=totalENV$ROIvariables_changed

            #load the BEDvariables
            #y=totalENV$BEDvariables

            BEDvariables$tempBED=totalENV$BEDvariables_tempBED
            BEDvariables$tempBEDname=totalENV$BEDvariables_tempBEDname
            BEDvariables$opened=totalENV$BEDvariables_opened
            BEDvariables$sfn=totalENV$BEDvariables_sfn


			#load DB variables, 
			DATABASEvariables$currentASSEMBLY=totalENV$DATABASEvariables_currentASSEMBLY
			#DATABASEvariables$availASSEMBLIES=totalENV$DATABASEvariables_availASSEMBLIES
			#DATABASEvariables$missingASSEMBLIES=totalENV$DATABASEvariables_missingASSEMBLIES

			#load the BAMvariables
			#bam=totalENV$BAMvariables
			BAMvariables$listBAM=totalENV$BAMvariables_listBAM
			#DATABASE variables not loaded! libraries must be re-loaded!
			#load the heatvariables
			#load the logvariables
			#logs=totalENV$logvariables
			logvariables$msg=totalENV$logvariables_msg
			#load the clustvariables

            logvariables$msg[[length(logvariables$msg)+1]]= paste('Session ',sessiontoload,' loaded<br>',sep="")
            print(paste('Session ',shinyFilePath(input$loadenv),' loaded',sep="" ))
            sendSweetAlert(
	          session = session,
	          title = "Session loaded!",
	          text = paste("The session '",sessiontoload,"' has been loaded",sep=""),
	          type = "success"
	        ) 
            
          },
          warning = function( w ){
          	#logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the rds... maybe is not the rds of your environment...<br></font>'
		      sendSweetAlert(
		        session = session,
		        title = "Cannot open file",
		        text = "Error in opening the file: maybe is not the *.rds file of a saved environment",
		        type = "error"
		      )           
          },
          error = function( err ){
			#logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Problems in opening the rds... maybe is not the rds of your environment...<br></font>'
		      sendSweetAlert(
		        session = session,
		        title = "Cannot open file",
		        text = "Error in opening the file: maybe is not the *.rds file of a saved environment",
		        type = "error"
		      )              
          })

    }else{
	    #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">No file selected...<br></font>'
	      # sendSweetAlert(
	      #   session = session,
	      #   title = "No files selected",
	      #   text = "Select an rds file from a saved session to load",
	      #   type = "error"
	      # )      
    }
   
},ignoreInit=TRUE)


