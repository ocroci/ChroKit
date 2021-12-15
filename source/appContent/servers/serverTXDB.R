




###react to help buttons:
#extract annotated elements from database
observeEvent(input$msg_databases_extractAnnotatedElements, {
  boxHelpServer(msg_databases_extractAnnotatedElements)
})
#download database
observeEvent(input$msg_databases_downloadDatabases, {
  boxHelpServer(msg_databases_downloadDatabases)
})







##extract annotation from an available assembly (txdb+org databases):
observeEvent(input$confirmASSEMBLYforuse,{
	#check if input is not NULL
	if (!is.null(input$searchASSEMBLYforuse) & input$searchASSEMBLYforuse!="" & length(input$searchASSEMBLYforuse)>0){
		#find if a previous different DB was used
		oldASSEMBLY=DATABASEvariables$currentASSEMBLY
		if(is.null(oldASSEMBLY)){
			oldASSEMBLY=FALSE
		}
		#extract TSS,transcripts and TES from the selected assembly; needs upstream/downstream and other paramters
		logvariables$msg[[length(logvariables$msg)+1]]= paste('Defining the promoters, genebodies, TES for ',input$searchASSEMBLYforuse,' assembly databsase<br>',sep="")	

		#defining transcripts with all annotation using extractfromDB function	
		transc=extractFromDB(assembly=input$searchASSEMBLYforuse,avail_assemblies=all_avail_assemblies)$transcripts
		
		#remove not conventional chromosomes
		transc=convertNomenclatureGR(transc,to="UCSC")

		#define promoters and TES from the transcripts
		promo=suppressWarnings(resize(transc,width=1,fix="start"))
		fixed_promo=promo
		tes=suppressWarnings(resize(transc,width=1,fix="end"))
		fixed_tes=tes
	
		#for asimmetric resize, we have upstreamTSS, downstreamTSS; upstreamTES, downstreamTES
		# input$absoluteFilterUpstreamTSS
		# input$absoluteFilterDownstreamTSS
		# input$absoluteFilterUpstreamTES
		# input$absoluteFilterDownstreamTES

		#for TSS:
		#find negative positions: 
		pos_negative=strand(promo)=="-"
		#find positive position:  
		pos_positive=!pos_negative

		#change sign of upstreams if "-" in the value from UI:
		# upstreamTSS=-(input$upstreamTSS)
		# upstreamTES=-(input$upstreamTES)
		#otherwise, do not change
		upstreamTSS=input$absoluteFilterUpstreamTSS
		upstreamTES=input$absoluteFilterUpstreamTES

		#filter transc, promo,tes for those that will have width <0 with current up/downstream parameters 
		negative_positivestrand=start(promo)-upstreamTSS<0 | start(tes)-upstreamTES <0
		negative_negativestrand=start(promo)-input$absoluteFilterDownstreamTSS <0|start(tes)-input$absoluteFilterDownstreamTES<0

		toremove=negative_positivestrand&pos_positive | negative_negativestrand&pos_negative
		
		#filter everything 
		transc=transc[!toremove]
		promo=promo[!toremove]
		tes=tes[!toremove]
		fixed_promo=fixed_promo[!toremove]
		fixed_tes=fixed_tes[!toremove]
		pos_positive=pos_positive[!toremove]
		pos_negative=pos_negative[!toremove]

		#this code gives warnings, because maybe some start/end exceed the length
		#of the chromosome
		options(warn=-1)
		start(promo[pos_positive])=start(promo[pos_positive])-upstreamTSS
		end(promo[pos_positive])=end(promo[pos_positive])+input$absoluteFilterDownstreamTSS
		#the same, but opposite, for negative strand positions:
		start(promo[pos_negative])=start(promo[pos_negative])-input$absoluteFilterDownstreamTSS
		end(promo[pos_negative])=end(promo[pos_negative])+upstreamTSS

		#same for TES
		start(tes[pos_positive])=start(tes[pos_positive])-upstreamTES
		end(tes[pos_positive])=end(tes[pos_positive])+input$absoluteFilterDownstreamTES
		#the same, but opposite, for negative strand positions:
		start(tes[pos_negative])=start(tes[pos_negative])-input$absoluteFilterDownstreamTES
		end(tes[pos_negative])=end(tes[pos_negative])+upstreamTES

		options(warn=0)
		#now keep only those with annotated gene IDS or symbols, if selected by the user
		#keep ALL isoforms for a specific gene
		
		toadd=""
		#prepare the array. by default, without filtering, everything is kept
		tokeep=rep(TRUE,length(transc))
		#verify that the user selected at least one filter
		if (!is.null(input$includeOnly)){
			if("ENTREZ"%in%input$includeOnly){
				colstrings=as.character(transc$gene_id)
				pos=!is.na(colstrings)
				tokeep=tokeep&pos
			}
			if("SYMBOL"%in%input$includeOnly){
				colstrings=as.character(transc$symbol)
				pos=!is.na(colstrings)
				tokeep=tokeep&pos
			}
			if("ENSEMBL"%in%input$includeOnly){
				colstrings=as.character(transc$ensembl_id)
				pos=!is.na(colstrings)
				tokeep=tokeep&pos
			}
			if("RefSeq"%in%input$includeOnly){
				colstrings=as.character(transc$refSeq_id)
				pos=!is.na(colstrings)
				tokeep=tokeep&pos				
			}
			#keep only those that survives (tokeep) the filtering steps for the existence of IDs/symbols
			toremove=table(tokeep)["FALSE"]
			transc=transc[tokeep]
			promo=promo[tokeep]
			tes=tes[tokeep]
			fixed_promo=fixed_promo[tokeep]
			fixed_tes=fixed_tes[tokeep]
			logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Lost ',toremove,' transcripts without desired IDs or symbols...<br></font>',sep="")
			toadd="(only with gene IDs)"
		}
		
		#add these new three special ROIs
		nomi=unlist(lapply(ROIvariables$listROI,getName))

		pos_promoters=!is.na(match(nomi,"promoters"))
		pos_transcripts=!is.na(match(nomi,"transcripts"))
		pos_TES=!is.na(match(nomi,"TES"))

        ####newenrichimplementation####
        #new ROI imported has no list of enrichments at the beginning! => initialization
        Enrichlist$rawcoverage[["promoters"]]=list()
        Enrichlist$normfactlist[["promoters"]]=list()
        ################################  
		newROIpromoter=new("RegionOfInterest",
							name="promoters",
							range=promo,
							fixed=fixed_promo,
							flag="promoterFlag",
							source=list(paste("promoters",toadd,"(with",upstreamTSS,"bp upstream and",input$absoluteFilterDownstreamTSS,"downstream TSS)")))
        ####newenrichimplementation####
        #new ROI imported has no list of enrichments at the beginning! => initialization
        Enrichlist$rawcoverage[["transcripts"]]=list()
        Enrichlist$normfactlist[["transcripts"]]=list()
        ################################		
		newROItranscripts=new("RegionOfInterest",
							name="transcripts",
							range=transc,
							fixed=resize(transc,width=1,fix="start"),
							flag="transcriptFlag",
							source=list(paste("transcripts",toadd)))
        ####newenrichimplementation####
        #new ROI imported has no list of enrichments at the beginning! => initialization
        Enrichlist$rawcoverage[["TES"]]=list()
        Enrichlist$normfactlist[["TES"]]=list()
        ################################
		newROItes=new("RegionOfInterest",
							name="TES",
							range=tes,
							fixed=fixed_tes,
							flag="TESFlag",
							source=list(paste("TES",toadd,"(with",upstreamTES,"bp upstream and",input$absoluteFilterDownstreamTES,"downstream TES)")))

		#if found matched name, substitute the special ROI, otherwise, create it
		if(any(pos_promoters)){
			pos_promoters=which(pos_promoters)
			ROIvariables$listROI[[pos_promoters]]=newROIpromoter
		}else{
			ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=newROIpromoter		
		}

		if (any(pos_transcripts)){
			pos_transcripts=which(pos_transcripts)
			ROIvariables$listROI[[pos_transcripts]]=newROItranscripts
		}else{
			ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=	newROItranscripts		
		}

		if (any(pos_TES)){
			pos_TES=which(pos_TES)
			ROIvariables$listROI[[pos_TES]]=newROItes
		}else{
			ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=newROItes				
		}

		#log the database currently used
		DATABASEvariables$currentASSEMBLY=input$searchASSEMBLYforuse
		logvariables$msg[[length(logvariables$msg)+1]]= paste('Set the promoters, genebodies, TES according to ',input$searchASSEMBLYforuse,' databsase<br>',sep="")
		
		############################################################
		############################################################
		#ANNOTATION ROIs
		#if db changed or the requirement t have entrez or other IDs
		#makes the annotation different

		#re-annotate each ROI
		
		
		#get all the ROIs (except promoters, transcripts, TES and derivatives)
		#the first ones are already annotated by definition, the others were from genes whose coords are from previous assembly/organism
		nomi=unlist(lapply(ROIvariables$listROI,getName))
		toavoid=c("promoters","transcripts","TES","promoters_genelist_","transcripts_genelist_","TES_genelist_")
		pos=!nomi%in%toavoid
		rois=ROIvariables$listROI[pos]
		nomi_toprocess=nomi[pos]
		fix_promoters=fixed_promo

		#execute annotation if we have at least one ROI
		if (length(rois)>0){
			###keep it criterion="midpoint" for now...
			#if nc==1 or windows system(nc==1), do it in single core
			if(nc==1){
				print("Re-annotating existing ROIs in single core...")
				annotatedrois=lapply(1:length(rois),function(i){
					currentroi=rois[[i]]
					currentrange=getRange(currentroi)
					annotatedpart=suppressWarnings(distanceFromTSS3(Object=currentrange,Tss=fix_promoters,criterion="midpoint"))
					#use setRange function to set the annotated range of the ROI
					elementMetadata(currentrange)=annotatedpart
					currentroi=setRange(currentroi,currentrange)
					logvariables$msg[[length(logvariables$msg)+1]]= paste('Re-annotating ',nomi_toprocess[i],' using ',input$searchASSEMBLYforuse,' genome assembly<br>',sep="")
					return(currentroi)
				})
				ROIvariables$listROI[pos]=annotatedrois
			#otherwise, if nc>1, do it in parallel, keeping the warnings/error if the case
			}else{

				tryCatch({
					print("Re-annotating existing ROIs in multi core...")
					annotatedrois=mclapply(1:length(rois),function(i){
						currentroi=rois[[i]]
						if (getFlag(currentroi)=="normalFlag"){
							currentrange=getRange(currentroi)
							annotatedpart=suppressWarnings(distanceFromTSS3(Object=currentrange,Tss=fix_promoters,criterion="midpoint"))
							#use setRange function to set the annotated range of the ROI
							elementMetadata(currentrange)=annotatedpart
							currentroi=setRange(currentroi,currentrange)							
						}else{
							print(paste("ROI ",getName(currentroi)," not annotated because of the flag originally from another database"))
						}

						return(currentroi)
					},mc.cores=nc)
					ROIvariables$listROI[pos]=annotatedrois
					logvariables$msg[[length(logvariables$msg)+1]]= paste('Re-annotating ROIs using ',input$searchASSEMBLYforuse,' genome assembly<br>',sep="")
				},
				warning=function(w){
					print("warning: problems in annotating ROIs...")
			      sendSweetAlert(
			        session = session,
			        title = "Problems in (re)annotating ROIs",
			        text = "Cannot (re)annotate ROIs. Try to use 1 single core, or maybe the genome is incompatible with existing ROIs",
			        type = "error"
			      ) 
				},
				error=function(err){
					print("error: problems in annotating ROIs...")

			      sendSweetAlert(
			        session = session,
			        title = "Problems in (re)annotating ROIs",
			        text = "Cannot (re)annotate ROIs. Try to use 1 single core, or maybe the genome is incompatible with existing ROIs",
			        type = "error"
			      ) 
				})				
			}
		}
		
	}
},ignoreInit=TRUE)



##download the missing database (using the function)
#and re-set the available and missing databases, to show to the user the correct list
observeEvent(input$confirmASSEMBLYfordownload,{
	#check if input is valid:
	if (!is.null(input$searchASSEMBLYfordownload) & input$searchASSEMBLYfordownload!="" & length(input$searchASSEMBLYfordownload)>0){
		#from assembly string (ex:"mm9") to download of the correct database from bioconductor
		#and upadte of DATABASEvariables$avail/missing database
		# tryCatch({
		if (checkBiocConnection()){
			downloadDB(assembly=input$searchASSEMBLYfordownload,avail_assemblies=all_avail_assemblies)
	    	#alert the user that the the assembly has been correctly downloaded
	       	updateExistingDB=getExistingDB(avail_assemblies=all_avail_assemblies)
			DATABASEvariables$availASSEMBLIES=updateExistingDB$assemblies_we_have
			DATABASEvariables$missingASSEMBLIES=updateExistingDB$assemblies_we_donthave
			logvariables$msg[[length(logvariables$msg)+1]]= paste('Downloaded databases for ',input$searchASSEMBLYfordownload,' assembly<br>',sep="")
	        sendSweetAlert(
	          session = session,
	          title = "Genome assembly downloaded!",
	          text = paste("The ",input$searchASSEMBLYfordownload," genome assembly has been correctly downloaded!",sep=""),
	          type = "success"
	        ) 
	    }else{
	    	sendSweetAlert(
		        session = session,
		        title = "Connection problems",
		        text = "Problems in connecting to bioconductor site",
		        type = "error"
		    )
		    return()
	    }

	   #  },
	   #  warning = function( w ){
		  # sendSweetAlert(
	   #        session = session,
	   #        title = "Problems in downloading database",
	   #        text = "Check your internet connection, or change bioCversion (for example, for R 3.6, the bioCversion 3.9 is needed)",
	   #        type = "error"
	   #    )
	   #    return()
	   #  },
	   #  error = function( err ){
	   #      sendSweetAlert(
	   #        session = session,
	   #        title = "Problems in downloading database",
	   #        text = "Check your internet connection, or change bioCversion (for example, for R 3.6, the bioCversion 3.9 is needed)",
	   #        type = "error"
	   #      )
	   #     	return()
	   #  })
		

	}
},ignoreInit=TRUE)


