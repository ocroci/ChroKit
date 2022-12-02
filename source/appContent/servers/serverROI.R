


###react to help buttons:
#options
observeEvent(input$msg_newRois_options, {
  boxHelpServer(msg_newRois_options)
})

#ROI combination
observeEvent(input$msg_newRois_ROIcombination, {
  boxHelpServer(msg_newRois_ROIcombination)
})





#given all the inputs for ROIs (with controls, append to ROIlist with its names)
#look at the button to make the ROI
observeEvent(input$maketheROI, {
  #check if at least one bed is selected from the checkbutton list (primary or secondary)
  if (length(input$selectedROIs)>0 & !is.na(input$minOverlapNEWROI)>0){

    nomi=unlist(lapply(ROIvariables$listROI,getName))
    #check if name of the ROI has been set and is not equal to created ROIs or primary ROIs (GRvariables$names):
    if (!input$ROIname %in% nomi) {

      if (input$ROIname!="promoters" & input$ROIname!="transcripts" & input$ROIname!="TES" &
          !grepl("promoters_genelist_",input$ROIname) &!grepl("transcripts_genelist_",input$ROIname)&!grepl("TES_genelist_",input$ROIname)) {
        #now  if passed controls, put in a list GRanges for creating ROIs, based on names

        #selected ranges
        pos=match(input$selectedROIs,nomi)
        selectedranges=lapply(ROIvariables$listROI[pos],getRange) 
        #if only one GR selected at the beginning, preserve its fixed point
        if(length(pos)==1){
          selectedfix=getFixed(ROIvariables$listROI[[pos]])
        }else{
          selectedfix=NULL
        }
        
        #if only one selected, keep the flag of old ROI, otherwise use "normalFlag"
        #and set the names of the source of the new ROI accordingly
        if(length(pos)==1){
          newflag=getFlag(ROIvariables$listROI[[pos]])
          toadd=paste("From [",unlist(getSource(ROIvariables$listROI[[pos]])),"]" )
        }else{
          #maybe should put the sources of all the ROIS in intersection or union
          newflag="normalFlag"
          toadd=paste("From",input$choiceROI,"of [",paste(input$selectedROIs,collapse=", "),"]",collapse=" ")        
        }

        #overlap with ranges
        if (length(input$overlapROIs)>0){
          pos2=match(input$overlapROIs,nomi)
          overlapregions=lapply(ROIvariables$listROI[pos2],getRange)

          if(length(input$overlapROIs)>1){
            if(input$choiceoverlapROI=="stringent"){
              toadd=paste(toadd,"that overlaps (minoverlap=",input$minOverlapNEWROI,"bp) with the intersection of [",paste(input$overlapROIs,collapse=" AND "),"]",collapse=" ")
            }
            if(input$choiceoverlapROI=="permissive"){
              toadd=paste(toadd,"that overlaps (minoverlap=",input$minOverlapNEWROI,"bp) with both [",paste(input$overlapROIs,collapse=", "),"]",collapse=" ")
            }
            if(input$choiceoverlapROI=="allofthem"){
              toadd=paste(toadd,"that overlaps (minoverlap=",input$minOverlapNEWROI,"bp) with any of [",paste(input$overlapROIs,collapse=", "),"]",collapse=" ")
            }                  
          }else{
            toadd=paste(toadd,"that overlaps with ",paste(input$overlapROIs,collapse=" "),collapse=" ")                
          }
 
          
        }else{
          overlapregions=NULL
        }

        #not overlap with ranges
        if (length(input$notoverlapROIs)>0){
          pos3=match(input$notoverlapROIs,nomi)
          notoverlapregions=lapply(ROIvariables$listROI[pos3],getRange) 

          if(length(input$notoverlapROIs)>1){
            if(input$choicenotoverlapROI=="intersection"){
              toadd=paste(toadd,"and do not overlap (minoverlap=",input$minOverlapNEWROI,"bp) with the intersection of [",paste(input$notoverlapROIs,collapse=" AND "),"]",collapse=" ")
            }
            if(input$choicenotoverlapROI=="union"){
              toadd=paste(toadd,"and do not overlap (minoverlap=",input$minOverlapNEWROI,"bp) with the union of [",paste(input$notoverlapROIs,collapse=", "),"]",collapse=" ")
            }
          }else{
            toadd=paste(toadd,"and do not overlap (minoverlap=",input$minOverlapNEWROI,"bp) with any of [",paste(input$notoverlapROIs,collapse=", "),"]",collapse=" ")
          }

               
        }else{
          notoverlapregions=NULL
        }   

        #if one ROI selected, keep the BAM files associated
        if (length(input$selectedROIs)==1){
          #selectedbam=getBAMlist(ROIvariables$listROI[[pos]])
          ####newenrichimplementation####
          selectedbam=Enrichlist$rawcoverage[[pos]]
          selecteddecryptkey=Enrichlist$decryptkey[[pos]]
          selectednormfact=Enrichlist$normfactlist[[pos]]
          ################################
        }else{
          selectedbam=list()
          selecteddecryptkey=list()
          ####newenrichimplementation####
          selectednormfact=list()
          ################################
        }

        if (nchar(input$ROIname)>=1){
          #if null (no multiple ROI selected or contrast), attrib. a random value, does not matter
          if(is.null(input$choiceROI)) {METHOD="union"} else{METHOD=input$choiceROI}
          if(is.null(input$choiceoverlapROI)) {CRITERION1="permissive"} else{CRITERION1=input$choiceoverlapROI}
          if(is.null(input$choicenotoverlapROI)) {CRITERION2="union"} else {CRITERION2=input$choicenotoverlapROI}

          roigeneration=suppressWarnings(generateROI(selectedranges,selectedfix,overlapregions,notoverlapregions,
                        method=METHOD,criterion1=CRITERION1, criterion2=CRITERION2,bamlist=selectedbam,decryptkeylist=selecteddecryptkey,minbp=input$minOverlapNEWROI,strandSpecific=input$StrandSpecOverlapNEWROI))

          ROI=roigeneration[[1]]

          if(length(ROI)>0){
            BAMlist=roigeneration[[2]]
            decryptkeylist=roigeneration[[3]]
            newfix=roigeneration[[4]]

            #if union or intersection (input$choiceROI) // if selectedranges >0,
            #reannotate the ROI because it changed
            if(length(selectedranges)>1){
              if(length(DATABASEvariables$currentASSEMBLY)>0){
                if(DATABASEvariables$currentASSEMBLY!=FALSE){
                  #we have a database. So extract the fix of promoters and apply the distanceFromTSS3 function
                  #for this newly created range
                  nomi=unlist(lapply(ROIvariables$listROI,getName))
                  pos_promo=match("promoters",nomi)
                  promo=ROIvariables$listROI[[pos_promo]]
                  fix_promoters=getFixed(promo)
                  annotatedpart=suppressWarnings(distanceFromTSS3(Object=ROI,Tss=fix_promoters,criterion="midpoint"))
                  elementMetadata(ROI)=annotatedpart
                }
              }              
            }
            ####newenrichimplementation####
            Enrichlist$rawcoverage[[input$ROIname]]=BAMlist
            Enrichlist$decryptkey[[input$ROIname]]=decryptkeylist
            Enrichlist$normfactlist[[input$ROIname]]=selectednormfact
            ################################
            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                        name=input$ROIname,
                                        range=ROI,
                                        fixed=newfix,
                                        BAMlist=list(),
                                        flag=newflag,
                                        source=list(toadd))     
            logvariables$msg[[length(logvariables$msg)+1]]= paste("Created ROI named ",input$ROIname,"<br>",sep="")      
            print(paste("Created ROI named",input$ROIname))       


            }else{
              #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">new ROI has length ==0<br></font>'
              sendSweetAlert(
                session = session,
                title = "Empty ROI",
                text = "New ROI not created, because has length = 0. No genomic regions are present with the overlaps requested",
                type = "error"
              )             
            }



        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">ROI name selected has length ==0<br></font>'
          sendSweetAlert(
            session = session,
            title = "Empty ROI name",
            text = "Put a name for the new ROI in 'Name of the ROI' field",
            type = "error"
          )        
        }

      } else{
        #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">New ROI cannot have names: \'promoters\' , \'transcripts\', \'TES\'<br></font>'
        sendSweetAlert(
          session = session,
          title = "Bad ROI name",
          text = "New ROI cannot be named 'promoters','transcripts','TES','promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
          type = "error"
        )       
      }
      
    } else{
      #error/warning message in logs, in red color
      #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Already exists another ROI with the same name...<br></font>'
      sendSweetAlert(
        session = session,
        title = "ROI already present",
        text = "A ROI with the same name already exists; use a different name",
        type = "error"
      )    
    }
  }else{
    #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">No coordinates files selected or minimum bp defined...<br></font>'
    sendSweetAlert(
      session = session,
      title = "Missing parameters",
      text = "You didn't select any ROI to start from or the minimum number of bp to consider for the overlap",
      type = "error"
    )  
  }
},ignoreInit=TRUE)









###react to help buttons:
#resize
observeEvent(input$msg_modifyRois_resize, {
  boxHelpServer(msg_modifyRois_resize)
})
#center on summit
observeEvent(input$msg_modifyRois_summit, {
  boxHelpServer(msg_modifyRois_summit)
})
#random subsample
observeEvent(input$msg_modifyRois_sample, {
  boxHelpServer(msg_modifyRois_sample)
})
#filter for width
observeEvent(input$msg_modifyRois_width, {
  boxHelpServer(msg_modifyRois_width)
})
#filter for enrichment
observeEvent(input$msg_modifyRois_enrichment, {
  boxHelpServer(msg_modifyRois_enrichment)
})
#pattern search
observeEvent(input$msg_modifyRois_pattern, {
  boxHelpServer(msg_modifyRois_pattern)
})

#resize ROI box
observeEvent(input$resizeROI,{
  #check if ROI list is not null
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1 & isvalid(input$sliderUpstreamROI) & isvalid(input$sliderDownstreamROI)) {
    #check if name is set
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$ROInameResize)>=1 & !(any(input$ROInameResize == nottobe))) {
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      
      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$ROInameResize) &!grepl(nottobe2[2],input$ROInameResize) & !grepl(nottobe2[3],input$ROInameResize)){
        #check if name does not match with another existing BED
        if (!input$ROInameResize %in% nomi ){
          if (!is.null(input$selectROItoresize)){
            if (input$selectROItoresize!="transcripts"){
              #resize with input$sliderUpstreamROI and input$sliderDownstreamROI
              #check if new size is greater. 
              #if both down and up are smaller, resize range
              #enrichments associated cannot be preserved due to LZ4

              pos=match(input$selectROItoresize,nomi)
              roi=getRange(ROIvariables$listROI[[pos]])
              oldSource=getSource(ROIvariables$listROI[[pos]])
              fix=start(getFixed(ROIvariables$listROI[[pos]]))   
              #split positive or undetermined strand from negative strand
              pos_positive=as.logical(!strand(roi)=="-")
              pos_negative=as.logical(strand(roi)=="-")

              if (input$chooseResizeType=="fixedVal"){
                if (input$choosePointResize=="fromMid"){
                  #remove ranges that would have negative starts with new width in positive strand
                  idx=(fix-input$sliderUpstreamROI)>0
                  idx_positive_and_good=pos_positive & idx
                  #negative strand
                  idx=(fix-input$sliderDownstreamROI)>0
                  idx_negative_and_good=pos_negative & idx
                  idx_positive_or_negative_good=idx_positive_and_good|idx_negative_and_good
                  roi=roi[idx_positive_or_negative_good]
                  fix=fix[idx_positive_or_negative_good]

                  pos_positive=as.logical(!strand(roi)=="-")
                  pos_negative=as.logical(strand(roi)=="-")
                  roi_positive=roi[pos_positive]
                  roi_negative=roi[pos_negative]
                  fix_positive=fix[pos_positive]
                  fix_negative=fix[pos_negative]

                  #resize positive strand:
                  oldstarts_positive=start(roi_positive)
                  oldends_positive=end(roi_positive)
                  oldstarts_negative=start(roi_negative)
                  oldends_negative=end(roi_negative)
                  
                  newstart=input$sliderUpstreamROI
                  newend=input$sliderDownstreamROI

                  toadd=paste("resized to window [-",input$sliderUpstreamROI,"; +",input$sliderDownstreamROI,"] from the centre",sep="")

                  #resize plus strand
                  start(roi_positive)=fix_positive-newstart
                  end(roi_positive)=fix_positive+newend
                  #resize negative strand
                  start(roi_negative)=fix_negative-newend
                  end(roi_negative)=fix_negative+newstart
                  #join in the original roi
                  roi[pos_positive]=roi_positive
                  roi[pos_negative]=roi_negative
                  newfix=getFixed(ROIvariables$listROI[[pos]])[idx_positive_or_negative_good] 
                  newrange=roi

                  if(length(unique(fix_positive-oldstarts_positive))==1 & length(unique(oldends_positive-fix_positive))==1){
                    logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI [-',(fix_positive-oldstarts_positive)[1],'; +',(oldends_positive-fix_positive)[1],']<br>',sep="")
                    print(paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI [-',(fix_positive-oldstarts_positive)[1],'; +',(oldends_positive-fix_positive)[1],']',sep=""))
                  }else{
                    logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI with original widths<br>',sep="")
                    print(paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI with original widths',sep=""))
                  }

                }else if (input$choosePointResize=="fromStart"){

                  #structure:

                  #             start          end
                  #---------------|-------------|------ +

                  #        upstr.  downstr.
                  #-------|-------*----------|--------- + (* is the new fixed point) 


                  #             start          end
                  #---------------|-------------|------------- -

                  #                      downstr.   upstr.
                  #---------------------|-------*----------|-- - (* is the new fixed point) 

                  #remove ranges that would have negative starts with new width in positive strand
                  idx=(start(roi)-input$sliderUpstreamROI)>0
                  idx_positive_and_good=pos_positive & idx
                  #negative strand
                  idx=(end(roi)-input$sliderDownstreamROI)>0
                  idx_negative_and_good=pos_negative & idx
                  idx_positive_or_negative_good=idx_positive_and_good|idx_negative_and_good
                  roi=roi[idx_positive_or_negative_good]
                  pos_positive=as.logical(!strand(roi)=="-")
                  pos_negative=as.logical(strand(roi)=="-")
                  roi_positive=roi[pos_positive]
                  roi_negative=roi[pos_negative]
                  #resize positive strand:
                  oldstarts_positive=start(roi_positive)
                  oldends_positive=end(roi_positive)
                  oldstarts_negative=start(roi_negative)
                  oldends_negative=end(roi_negative)
                  
                  newstart=input$sliderUpstreamROI
                  newend=input$sliderDownstreamROI

                  toadd=paste("resized to window [-",input$sliderUpstreamROI,"; +",input$sliderDownstreamROI,"] from the start",sep="")

                  #resize plus strand
                  start(roi_positive)=oldstarts_positive-newstart
                  end(roi_positive)=oldstarts_positive+newend
                  #resize negative strand
                  start(roi_negative)=oldends_negative-newend
                  end(roi_negative)=oldends_negative+newstart
                  #join in the original roi
                  roi[pos_positive]=roi_positive
                  roi[pos_negative]=roi_negative

                  newfix=roi
                  start(newfix)[pos_positive]=end(newfix)[pos_positive]=oldstarts_positive
                  start(newfix)[pos_negative]=end(newfix)[pos_negative]=oldends_negative
                  newrange=roi
                  logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI (resized ',newstart,'bp upstream and ',newend,' bp downstream from start) from ',input$selectROItoresize,' ROI<br>',sep="")
                  print(paste('Created ',input$ROInameResize,' ROI (resized ',newstart,'bp upstream and ',newend,' bp from start) from ',input$selectROItoresize,' ROI',sep=""))             


                
                }else if (input$choosePointResize=="fromEnd"){

                  #structure:

                  #             start          end
                  #---------------|-------------|------ -

                  #        downstr. upstr.
                  #-------|-------*----------|--------- - (* is the new fixed point) 


                  #             start          end
                  #---------------|-------------|------------- +

                  #                      upstr.   downstr.
                  #---------------------|-------*----------|-- + (* is the new fixed point) 


                 #remove ranges that would have negative starts with new width in positive strand
                  idx=(end(roi)-input$sliderUpstreamROI)>0
                  idx_positive_and_good=pos_positive & idx
                  #negative strand
                  idx=(start(roi)-input$sliderDownstreamROI)>0
                  idx_negative_and_good=pos_negative & idx
                  idx_positive_or_negative_good=idx_positive_and_good|idx_negative_and_good
                  roi=roi[idx_positive_or_negative_good]
                  pos_positive=as.logical(!strand(roi)=="-")
                  pos_negative=as.logical(strand(roi)=="-")
                  roi_positive=roi[pos_positive]
                  roi_negative=roi[pos_negative]
                  #resize positive strand:
                  oldstarts_positive=start(roi_positive)
                  oldends_positive=end(roi_positive)
                  oldstarts_negative=start(roi_negative)
                  oldends_negative=end(roi_negative)
                  
                  newstart=input$sliderUpstreamROI
                  newend=input$sliderDownstreamROI

                  toadd=paste("resized to window [-",input$sliderUpstreamROI,"; +",input$sliderDownstreamROI,"] from the end",sep="")

                  #resize plus strand
                  start(roi_positive)=oldends_positive-newstart
                  end(roi_positive)=oldends_positive+newend
                  #resize negative strand
                  start(roi_negative)=oldstarts_negative-newend
                  end(roi_negative)=oldstarts_negative+newstart
                  #join in the original roi
                  roi[pos_positive]=roi_positive
                  roi[pos_negative]=roi_negative

                  newfix=roi
                  start(newfix)[pos_positive]=end(newfix)[pos_positive]=oldends_positive
                  start(newfix)[pos_negative]=end(newfix)[pos_negative]=oldstarts_negative
                  newrange=roi
                  logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI (resized ',newstart,'bp upstream and ',newend,' bp downstream from end) from ',input$selectROItoresize,' ROI<br>',sep="")
                  print(paste('Created ',input$ROInameResize,' ROI (resized ',newstart,'bp upstream and ',newend,' bp downstream from end) from ',input$selectROItoresize,' ROI',sep=""))             


                }
              }else{
                valueperc=input$chooseWidthPercResize
                value=floor(width(roi)* (valueperc/100) /2)
                if(input$choosePointResize=="increment"){
                  #independently on the strand
                  idx=(start(roi)-value)>0
                  roi=roi[idx]
                  fix=fix[idx] 
                  value=value[idx]
                  newfix=getFixed(ROIvariables$listROI[[pos]])[idx] 
                  #resize 
                  start(roi)=start(roi)-value
                  end(roi)=end(roi)+value
                  newrange=roi 
                  toadd=paste("Incremented width (",input$chooseWidthPercResize,"% from the centre)",sep="")                 
                  logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI (resized, incremented ',valueperc,' % from center) from ',input$selectROItoresize,' ROI<br>',sep="")
                  print(paste('Created ',input$ROInameResize,' ROI (resized, incremented ',valueperc,' % from center) from ',input$selectROItoresize,' ROI',sep=""))
                }else{
                  start(roi)=start(roi)+value
                  end(roi)=end(roi)-value    
                  newfix=getFixed(ROIvariables$listROI[[pos]])
                  newrange=roi 
                  toadd=paste("Decremented width (",input$chooseWidthPercResize,"% from the centre)",sep="")                 
                  logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI (resized, decremented ',valueperc,' % from center) from ',input$selectROItoresize,' ROI<br>',sep="")
                  print(paste('Created ',input$ROInameResize,' ROI (resized, decremented ',valueperc,' % from center) from ',input$selectROItoresize,' ROI',sep=""))             
                }


              }





              newflag=getFlag(ROIvariables$listROI[[pos]])
              newSource=c(oldSource,list(toadd))
              ####newenrichimplementation####
              Enrichlist$rawcoverage[[input$ROInameResize]]=list()
              Enrichlist$decryptkey[[input$ROInameResize]]=list()
              Enrichlist$normfactlist[[input$ROInameResize]]=list()
              ################################
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                            name=input$ROInameResize,
                                            range=newrange,
                                            fixed=newfix,
                                            BAMlist=list(),
                                            flag=newflag,
                                            source=newSource) 



            }else{
              #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Are you trying to resize transcripts around midpoint? You cannot.<br></font>',sep="")
              sendSweetAlert(
                session = session,
                title = "Bad choice",
                text = "You cannot resize transcripts",
                type = "error"
              )            
            }

          }
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameResize,' name already present in ROIs...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "ROI already present",
            text = "A ROI with the same name already exists; use a different name",
            type = "error"
          )        
        }

      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameResize,' is not a valid kind of name...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Bad ROI name",
          text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
          type = "error"
        )      
      }
    }else{
      #name is TSS or transcripts or is empty
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameResize,' is not a valid name...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Bad ROI name",
        text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
        type = "error"
      )    
    }
  }else{
    #ROI list is null or up/downstream not set
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Check input fields...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Missing parameters",
      text = "Up or Down stream bp are missing from the input fields",
      type = "error"
    )  
  }
},ignoreInit=TRUE)




#center on summit box
observeEvent(input$SummitROI,{
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1 & length(input$selectROItoCenterSummit)>0){
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$ROInameSummit)>=1 & !(any(input$ROInameSummit == nottobe))) {
      

      nomi=unlist(lapply(ROIvariables$listROI,getName))
        
      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$ROInameSummit) &!grepl(nottobe2[2],input$ROInameSummit) & !grepl(nottobe2[3],input$ROInameSummit)){


        if (!input$ROInameSummit %in% nomi){

          #new roi is the result of summitFromBaseCoverage function. BAMlist is removed
          #new fix is the same as the new range, flag is reset (summit on promoters is not "promoter" anymore)
          #name is the input$ROInameSummit name
          #let also the transcripts to be centered on the summit

          #extract the dsired ROI.
          nomi=unlist(lapply(ROIvariables$listROI,getName))
          pos=match(input$selectROItoCenterSummit,nomi)
          roi=ROIvariables$listROI[[pos]]
          ####newenrichimplementation####
          #we have to find maximum enrichment point. Raw cov w/o normalization is ok
          bam=Enrichlist$rawcoverage[[pos]]
          key=Enrichlist$decryptkey[[pos]]
          #normvals=Enrichlist$normfactlist[[pos]]
          ################################
          oldSource=getSource(roi)
          #getbam=names(getBAMlist(roi))
          getbam=names(bam)
          #if there is at least one BAM
          if (!is.null(getbam)){
            #bam=getBAMlist(roi)
            pos2=match(input$selectBAMtoCenterSummit,getbam)
            bamselected_tmp=bam[[pos2]]
            keyselected=key[[pos2]]
            #decrypt to integer for calculating summit
            #bamselected=decryptcov( list(bamselected_tmp,keyselected),chunk=length(bamselected_tmp))
            toadd=paste("Centered on summit using",input$selectBAMtoCenterSummit,"enrichment (range width=1)")
            #use the function to select summit:
            newrange=summitFromBaseCoverage(Object=getRange(roi),baseCoverageOutput=bamselected_tmp,keys=keyselected)

            ##this ROI is not a subset, but is a new ROI. Re-annotate if DB present (not FALSE, not null)
            if(length(DATABASEvariables$currentASSEMBLY)>0){
              if(DATABASEvariables$currentASSEMBLY!=FALSE){
                #we have a database. So extract the fix of promoters and apply the distanceFromTSS3 function
                #for this newly created range
                nomi=unlist(lapply(ROIvariables$listROI,getName))
                pos_promo=match("promoters",nomi)
                promo=ROIvariables$listROI[[pos_promo]]
                fix_promoters=getFixed(promo)
                annotatedpart=suppressWarnings(distanceFromTSS3(Object=newrange,Tss=fix_promoters,criterion="midpoint"))
                elementMetadata(newrange)=annotatedpart
              }
            }
            ####newenrichimplementation####
            Enrichlist$rawcoverage[[input$ROInameSummit]]=list()
            Enrichlist$decryptkey[[input$ROInameSummit]]=list()
            Enrichlist$normfactlist[[input$ROInameSummit]]=list()
            ################################
            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                            name=input$ROInameSummit,
                                            range=newrange,
                                            fixed=newrange,
                                            BAMlist=list(),
                                            flag="normalFlag",
                                            source=c(oldSource,list(toadd))) 
            logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameSummit,' ROI as summit of ',input$selectROItoCenterSummit,' ROI using ',input$selectBAMtoCenterSummit,' bam<br>',sep="")
            print(paste('Created ',input$ROInameSummit,' ROI as summit of ',input$selectROItoCenterSummit,' ROI using ',input$selectBAMtoCenterSummit,' bam',sep=""))
          }

        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameSummit,' name already exist in ROIs...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "ROI already present",
            text = "A ROI with the same name already exists; use a different name",
            type = "error"
          )        
        }

      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameSummit,' is not a valid kind of name...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Bad ROI name",
          text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
          type = "error"
        )       
      }

    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameSummit,' is not a valid name...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Bad ROI name",
        text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
        type = "error"
      )    
    }
  }else{
    #ROI not present
    sendSweetAlert(
      session = session,
      title = "No ROI selected",
      text = "Select one ROI",
      type = "error"
    )
  }

},ignoreInit=TRUE)


#Filter the ROI button for enrichment 
observeEvent(input$FilterROI,{
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1 & length(input$selectROItoFilter)>0 & isvalid(input$absoluteFilter1) & isvalid(input$absoluteFilter2)){
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$ROInameFilter)>=1 & !(any(input$ROInameFilter == nottobe))) {
      nomi=unlist(lapply(ROIvariables$listROI,getName))


      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$ROInameFilter) &!grepl(nottobe2[2],input$ROInameFilter) & !grepl(nottobe2[3],input$ROInameFilter)){

        if (!input$ROInameFilter %in% nomi ){
          nomi=unlist(lapply(ROIvariables$listROI,getName))
          pos=match(input$selectROItoFilter,nomi)
          roi=ROIvariables$listROI[[pos]]
          oldSource=getSource(roi)
          #getbam=names(getBAMlist(roi))
          ####newenrichimplementation####
          #we have to find maximum enrichment point. Raw cov w/o normalization is ok
          bam=Enrichlist$rawcoverage[[pos]]
          keyvals=Enrichlist$decryptkey[[pos]]
          normvals=Enrichlist$normfactlist[[pos]]
          ################################
          getbam=names(bam)
          if (!is.null(getbam)){

            #filter based not normalized coverage should be ok (it's a fraction, not abs number)
            #bam=getBAMlist(roi)
            pos2=match(input$selectBAMtoFilter,getbam)
            bamselected_tmp=bam[[pos2]]  
            keyselected=keyvals[[pos2]]
            norm_selected=normvals[[pos2]]
            
            sums=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bamselected_tmp,Nbins=1,Snorm=FALSE,key=keyselected,norm_factor=norm_selected)
            #bamselected=decryptcov( list(bamselected_tmp,keyselected),chunk=length(bamselected_tmp))
            #sums=unlist(lapply(bamselected,sum))
            #sums=sums*norm_selected
            
            #find the position in which sums is > than input$absoluteFilter1
            # and < input$absoluteFilter2
            #then make the new roi, with input$ROInameFilter name, newrange is old range in 
            #position in which sums > input$absoluteFilter; new fix is the old fix, new flag is the old flag

            selectedPositions1=sums>input$absoluteFilter1
            selectedPositions2=sums<input$absoluteFilter2
            
            excludedlow=round(table(selectedPositions1)["FALSE"]/length(sums)*100,1)
            excludedhigh=round(table(selectedPositions2)["FALSE"]/length(sums)*100,1)
            
            selectedPositions=selectedPositions1 & selectedPositions2
            tabForMessage=table(selectedPositions)
            newrange=getRange(roi)[selectedPositions]


            if(length(newrange)>0){
              newBAMlist=list()
              newkeylist=list()
              for(i in 1:length(bam)){
                newBAMlist[[i]]=bam[[i]][selectedPositions]
                newkeylist[[i]]=keyvals[[i]][selectedPositions]
              }
              names(newBAMlist)=names(newkeylist)=getbam
              newfix=getFixed(roi)[selectedPositions]
              newflag=getFlag(roi)
              perc=round(tabForMessage["TRUE"]/length(selectedPositions)*100,2)
              toadd=paste("Filtered on",input$selectBAMtoFilter,"BAM file enrichment (kept ",tabForMessage["TRUE"],"/",length(selectedPositions),"ranges,",perc,"%), excluding",excludedlow,"% low and",excludedhigh,"% high")
              newSource=c(oldSource,list(toadd))

              ####newenrichimplementation####
              Enrichlist$rawcoverage[[input$ROInameFilter]]=newBAMlist
              Enrichlist$decryptkey[[input$ROInameFilter]]=newkeylist
              Enrichlist$normfactlist[[input$ROInameFilter]]=normvals
              ################################ 
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                              name=input$ROInameFilter,
                                              range=newrange,
                                              fixed=newfix,
                                              BAMlist=list(),
                                              flag=newflag,
                                              source=newSource)
              logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameFilter,' ROI after filtering of ',input$selectROItoFilter,' ROI using ',input$selectBAMtoFilter,' bam (kept ',tabForMessage["TRUE"],'/',length(selectedPositions),' ranges)<br>',sep="")    
              print(paste('Created ',input$ROInameFilter,' ROI after filtering of ',input$selectROItoFilter,' ROI using ',input$selectBAMtoFilter,' bam (kept ',tabForMessage["TRUE"],'/',length(selectedPositions),' ranges)',sep=""))          
            }else{
              #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$selectROItoFilter,' filtering produced a ROI with length=0...<br></font>',sep="")
              sendSweetAlert(
                session = session,
                title = "Empty ROI produced",
                text = paste("Filtering '",input$selectROItoFilter,"' produced a ROI with length=0",sep=""),
                type = "error"
              )            

            }
          }
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameFilter,' name already exist in ROIs...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "ROI already present",
            text = "A ROI with the same name already exists; use a different name",
            type = "error"
          )        
        }
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameFilter,' is not a valid kind of name...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Bad ROI name",
          text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
          type = "error"
        )
      }
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameFilter,' is not a valid name...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Bad ROI name",
        text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
        type = "error"
      )    
    }
  }else{
    #ROI not present or missing input fields
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Check input fields...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Missing parameters",
      text = "Filtering thresholds not specified in the input fields",
      type = "error"
    )  
  }  
},ignoreInit=TRUE)






#Filter the ROI button for width
observeEvent(input$FilterROIWIDTH,{
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1 & length(input$selectROItoFilterWIDTH)>0 & isvalid(input$absoluteFilter1WIDTH) & isvalid(input$absoluteFilter2WIDTH)){
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$ROInameFilterWIDTH)>=1 & !(any(input$ROInameFilterWIDTH == nottobe))) {
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$ROInameFilterWIDTH) &!grepl(nottobe2[2],input$ROInameFilterWIDTH) & !grepl(nottobe2[3],input$ROInameFilterWIDTH)){

        if (!input$ROInameFilterWIDTH %in% nomi ){
          nomi=unlist(lapply(ROIvariables$listROI,getName))
          pos=match(input$selectROItoFilterWIDTH,nomi)
          roi=ROIvariables$listROI[[pos]]
          oldSource=getSource(roi)
          
          #check if ROI is "fixed size". In this case, we cannot filtered based on width
          getwdth=getWidth(roi)
          isdup=table(!duplicated(getwdth))["TRUE"]==1
          
          if(!isdup){
            sums=width(getRange(roi))
            #find the position in which sums is > than input$absoluteFilter1
            # and < input$absoluteFilter2
            #then make the new roi, with input$ROInameFilterWIDTH name, newrange is old range in 
            #position in which sums > input$absoluteFilter; new fix is the old fix, new flag is the old flag
            selectedPositions1=sums>input$absoluteFilter1WIDTH
            selectedPositions2=sums<input$absoluteFilter2WIDTH

            excludedlow=round(table(selectedPositions1)["FALSE"]/length(sums)*100,1)
            excludedhigh=round(table(selectedPositions2)["FALSE"]/length(sums)*100,1)

            selectedPositions=selectedPositions1 & selectedPositions2
            tabForMessage=table(selectedPositions)
            newrange=getRange(roi)[selectedPositions]

            if(length(newrange)>0){
              #bams=getBAMlist(roi)

              ####newenrichimplementation####
              #we have to find maximum enrichment point. Raw cov w/o normalization is ok
              bams=Enrichlist$rawcoverage[[pos]]
              keyvals=Enrichlist$decryptkey[[pos]]
              normvals=Enrichlist$normfactlist[[pos]]
              ################################

              if(length(bams)>0){
                newBAMlist=list()
                newkeylist=list()
                for(i in 1:length(bams)){
                  newBAMlist[[i]]=bams[[i]][selectedPositions]
                  newkeylist[[i]]=keyvals[[i]][selectedPositions]
                }
                names(newBAMlist)=names(newkeylist)=names(bams)            
              }else{
                newBAMlist=list()
                newkeylist=list()
              }
              
              newfix=getFixed(roi)[selectedPositions]
              newflag=getFlag(roi)
              perc=round(tabForMessage["TRUE"]/length(selectedPositions)*100,2)
              toadd=paste("Filtered on width (kept ",tabForMessage["TRUE"],"/",length(selectedPositions),"ranges,",perc,"%), excluding",excludedlow,"% low and",excludedhigh,"% high")
              newSource=c(oldSource,list(toadd))

              ####newenrichimplementation####
              Enrichlist$rawcoverage[[input$ROInameFilterWIDTH]]=newBAMlist
              Enrichlist$decryptkey[[input$ROInameFilterWIDTH]]=newkeylist
              Enrichlist$normfactlist[[input$ROInameFilterWIDTH]]=normvals
              ################################ 

              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                              name=input$ROInameFilterWIDTH,
                                              range=newrange,
                                              fixed=newfix,
                                              BAMlist=list(),
                                              flag=newflag,
                                              source=newSource)
              logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameFilterWIDTH,' ROI after filtering of ',input$selectROItoFilter,' ROI by length (kept ',tabForMessage["TRUE"],'/',length(selectedPositions),' ranges)<br>',sep="")   
              print(paste('Created ',input$ROInameFilterWIDTH,' ROI after filtering of ',input$selectROItoFilter,' ROI by length (kept ',tabForMessage["TRUE"],'/',length(selectedPositions),' ranges)',sep=""))    

            }else{
              #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$selectROItoFilterWIDTH,' filtering produced a ROI with length=0...<br></font>',sep="")
              sendSweetAlert(
                session = session,
                title = "Empty ROI produced",
                text = paste("Filtering '",input$selectROItoFilterWIDTH,"' produced a ROI with length=0",sep=""),
                type = "error"
              )            
            }

          }else{
            #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$selectROItoFilterWIDTH,' cannot be selected by width, has a fixed size...<br></font>',sep="")
            sendSweetAlert(
              session = session,
              title = "Cannot resize",
              text = paste("Cannot filter '",input$selectROItoFilterWIDTH,"' for width, because it has fixed width",sep=""),
              type = "error"
            )          
          }

        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameFilterWIDTH,' is not a valid name...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "ROI already present",
            text = "A ROI with the same name already exists; use a different name",
            type = "error"
          )        

        }

      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameFilterWIDTH,' is not a valid kind of name...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Bad ROI name",
          text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
          type = "error"
        )

      }
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameFilterWIDTH,' is not a valid name...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Bad ROI name",
        text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
        type = "error"
      )

    }
  }else{
    #ROI not present or Width input filds not set
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Check input fields...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Missing parameters",
      text = "Filtering thresholds not specified in the input fields",
      type = "error"
    )  
  }  
},ignoreInit=TRUE)


#Filter the ROI button for random sampling
observeEvent(input$SampleROI,{
  set.seed(123)
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1 & length(input$selectROItoSample)>0 &isvalid(input$numberSample)){
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$ROInameSample)>=1 & !(any(input$ROInameSample == nottobe))) {
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$ROInameSample) &!grepl(nottobe2[2],input$ROInameSample) & !grepl(nottobe2[3],input$ROInameSample)){

        if (!input$ROInameSample %in% nomi ){
          nomi=unlist(lapply(ROIvariables$listROI,getName))
          pos=match(input$selectROItoSample,nomi)
          roi=ROIvariables$listROI[[pos]]
          oldSource=getSource(roi)

          #sample ranges and BAM associated (with the same order)
          #how many? look at the input$numberSample

          rangefromROI=getRange(roi)
          selectedPositions=sort(sample(length(rangefromROI),input$numberSample,replace=FALSE))
          newrange=rangefromROI[selectedPositions]

          if(length(newrange)>0){
            #bams=getBAMlist(roi)
            ####newenrichimplementation####
            #we have to find maximum enrichment point. Raw cov w/o normalization is ok
            bams=Enrichlist$rawcoverage[[pos]]
            keys=Enrichlist$decryptkey[[pos]]
            normvals=Enrichlist$normfactlist[[pos]]
            ################################
            if(length(bams)>0){
              newBAMlist=list()
              newkeylist=list()
              for(i in 1:length(bams)){
                newBAMlist[[i]]=bams[[i]][selectedPositions]
                newkeylist[[i]]=keys[[i]][selectedPositions]
              }
              names(newBAMlist)=names(newkeylist)=names(bams)            
            }else{
              newBAMlist=list()
              newkeylist=list()
            }
            
            newfix=getFixed(roi)[selectedPositions]
            newflag=getFlag(roi)
            perc=round(input$numberSample/length(rangefromROI)*100,2)

            toadd=paste("ROI sampled (kept ",input$numberSample,"/",length(rangefromROI),"ranges,",perc,"%)")
            newSource=c(oldSource,list(toadd))
            ####newenrichimplementation####
            Enrichlist$rawcoverage[[input$ROInameSample]]=newBAMlist
            Enrichlist$decryptkey[[input$ROInameSample]]=newkeylist
            Enrichlist$normfactlist[[input$ROInameSample]]=normvals
            ################################ 
            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                            name=input$ROInameSample,
                                            range=newrange,
                                            fixed=newfix,
                                            BAMlist=list(),
                                            flag=newflag,
                                            source=newSource)

            logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameSample,' ROI after sampling of ',input$selectROItoSample,' (kept ',perc,'%, ',input$numberSample,'/',length(rangefromROI),' ranges)<br>',sep="")  
            print(paste('Created ',input$ROInameSample,' ROI after sampling of ',input$selectROItoSample,' (kept ',perc,'%, ',input$numberSample,'/',length(rangefromROI),' ranges)',sep=""))          
          }else{
            #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$selectROItoSample,' sampling produced a ROI with length=0...<br></font>',sep="")
            sendSweetAlert(
              session = session,
              title = "Empty ROI produced",
              text = paste("Random sampling '",input$selectROItoSample,"' produced a ROI with length=0",sep=""),
              type = "error"
            )          
          }

        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameSample,' is not a valid name...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "ROI already present",
            text = "A ROI with the same name already exists; use a different name",
            type = "error"
          )
        }

      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameSample,' is not a valid kind of name...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Bad ROI name",
          text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
          type = "error"
        )

      }
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$ROInameSample,' is not a valid name...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Bad ROI name",
        text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
        type = "error"
      )

    }
  }else{
    #ROI not present or missing sample field
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Check the input field...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Missing parameters",
      text = "The number of random ranges to select from the ROI is missing from the input field",
      type = "error"
    )  
  }  
},ignoreInit=TRUE)


# toListenDOWNLOADbsgenome <- reactive({
#     list(input$DownloadBSgenome,input$DownloadBSgenome2)
# })
#download BSgenome selected; when done, set the variable ($currentORG and import library)
observeEvent(input$DownloadBSgenome,{
  asm=DATABASEvariables$currentASSEMBLY
  avail_spl=strsplit(all_avail_assemblies,split="\\.")
  org=sapply(avail_spl,"[[",2)
  asms=sapply(avail_spl,"[[",4)
  pos=match(asm,asms)
  BSstring=paste("BSgenome.",org[pos],".UCSC.",asm,sep="")

  x=rownames(installed.packages())
  #install with correct R version from bioC
  Rversion=strsplit(version$version.string,split=" ")[[1]][3]
  Rversion=strsplit(Rversion,split="\\.")[[1]]
  Rversion_main=as.numeric(Rversion[1])
  Rversion_submain=as.numeric(Rversion[2])
  if(Rversion_main==3){
    if(Rversion_submain>=5){
      R35=TRUE
    }else{
      R35=FALSE
    }
  }else if(Rversion_main>3){
    R35=TRUE
  }else{
    R35=FALSE
  }

  #BiocManager if R version> 3.5
  if(! ("BiocManager" %in% x) & R35){
    print("Installing BiocManager package for R > 3.5.0 ...")
    install.packages("BiocManager")
  }else{
    #print("BiocManager package already installed...")
  }

  #install the correct BSgenome package

  # tryCatch({
  if (checkBiocConnection()){
    if(R35){
      print(paste("updating bioconductor packages to the version",bioCversion))
      BiocManager::install(version = bioCversion,ask=FALSE,force=TRUE)
      print(paste("Downloading",BSstring,"..."))
      BiocManager::install(BSstring, version = bioCversion,ask=FALSE,force=TRUE,type="source")
    }else{
      print(paste("Downloading",BSstring,"..."))
      biocLite(BSstring,suppressUpdates=TRUE,ask=FALSE)
    }

    #update reactive variable
    DATABASEvariables$currentORG=BSstring
    #import library
    library(BSstring,character.only=TRUE)

    #alert the user the package has been installed
    sendSweetAlert(
      session = session,
      title = paste(BSstring, "installed!"),
      text = paste("The package '",BSstring,"' has been installed and imported for pattern extraction from ROIs",sep=""),
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


},ignoreInit=TRUE)







#observer for search pattern button.
observeEvent(input$ExtractPatternROI,{
  #checks:
  #-ROI must exist
  #-database must be active
  #-BSgenome of correct assembly must be present
  #-pattern should not be empty and valid (check validity in the function)
  #-input$strandOptsPattern if NULL, bothstrands is TRUE, otherwise look
  #-new ROI name must not exist, and be not promoter,...
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1 ){
    #now check if assembly present
    if(length(DATABASEvariables$currentASSEMBLY)>0){
      #now check if BSgenome of that assembly is present
      asm=DATABASEvariables$currentASSEMBLY
      avail_spl=strsplit(all_avail_assemblies,split="\\.")
      org=sapply(avail_spl,"[[",2)
      asms=sapply(avail_spl,"[[",4)
      pos=match(asm,asms)
      BSstring=paste("BSgenome.",org[pos],".UCSC.",asm,sep="")
      x=rownames(installed.packages())
      pos_pkg=match(BSstring,x)
      if(!is.na(pos_pkg)){
        #now check if pattern is not empty (and maybe with length 2<x<n )
        if(isvalid(input$PatternToSearch)>0){
          #now check if pattern is valid
          pattern=toupper(input$PatternToSearch)
          pattern_splitted=strsplit(pattern,split="")[[1]]

          check=pattern_splitted=="." | pattern_splitted=="-" |
                pattern_splitted=="A" | pattern_splitted=="C" | pattern_splitted=="G" | pattern_splitted=="T" |
                pattern_splitted=="U" | pattern_splitted=="R" | pattern_splitted=="Y" | pattern_splitted=="S" |
                pattern_splitted=="W" | pattern_splitted=="K" | pattern_splitted=="M" | pattern_splitted=="B" |
                pattern_splitted=="D" | pattern_splitted=="H" | pattern_splitted=="V" | pattern_splitted=="N"
          if(all(check)){
            #check if, if fromGenome selected, pattern shuld be >=4 bp

            
            #check if ROI name is ok
            nottobe=c("promoters","transcripts","TES")
            if(!(any(input$ROInamePattern == nottobe)) & nchar(input$ROInamePattern)>0){
              #check if ROI name is *****_genelist_....
              nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
              if(!grepl(nottobe2[1],input$ROInamePattern) &!grepl(nottobe2[2],input$ROInamePattern) & !grepl(nottobe2[3],input$ROInamePattern)){
                #now check if ROI name already present
                nomi=unlist(lapply(ROIvariables$listROI,getName))
                if(!input$ROInamePattern %in% nomi){
                  #execute. Extract range with pattern (Flag="pattern").
                  #BAM file is NULL (no sense to keep enrichment of 6-10 bp...)
                  #shuld re-annotate ranges (DB is present by definition at this point)
                  #fix is the center of new range (center of the pattern)
                  #source: add "extracted CNNTACGT from *** ROI (both/s-specific)"

                  if(isvalid(input$strandOptsPattern)){
                    if(input$strandOptsPattern=="bothStrands"){
                      both=TRUE
                    }else{
                      both=FALSE
                    }
                  }else{
                    both=TRUE
                  }

                  toadd=paste("extracting of ",pattern," pattern",sep="")
                  if(both){
                    toadd=paste(toadd," searched in both strands",sep="")
                  }else{
                    toadd=paste(toadd," searched in strand-specific way",sep="")
                  }


                  pos=match(input$selectROItoExtractPattern,nomi)
                  roi=ROIvariables$listROI[[pos]]
                  roi_range=getRange(roi)
                  oldSource=getSource(roi)
                  newSource=c(oldSource,list(toadd))
                  #apply the function extractPattern:
                  extracted=suppressWarnings(extractPattern(Subject=roi_range,BSgenomeDB=eval(parse(text=BSstring)),pattern=pattern,bothstrands=both))

                  if(length(extracted)>0){
                    #find Fix of new ranges (midpoint)
                    newFix=resize(extracted,fix="center",width=1)
                    #annotate this new ROI. By definition at this point we have the database
                    nomi=unlist(lapply(ROIvariables$listROI,getName))
                    pos_promo=match("promoters",nomi)
                    promo=ROIvariables$listROI[[pos_promo]]
                    fix_promoters=getFixed(promo)
                    annotatedpart=suppressWarnings(distanceFromTSS3(Object=extracted,Tss=fix_promoters,criterion="midpoint"))
                    elementMetadata(extracted)=annotatedpart
                    ####newenrichimplementation####
                    #new ROI imported has no list of enrichments at the beginning! => initialization
                    Enrichlist$rawcoverage[[input$ROInamePattern]]=list()
                    Enrichlist$decryptkey[[input$ROInamePattern]]=list()
                    Enrichlist$normfactlist[[input$ROInamePattern]]=list()
                    ################################
                    ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                    name=input$ROInamePattern,
                                                    range=extracted,
                                                    fixed=newFix,
                                                    BAMlist=list(),
                                                    flag="Pattern",
                                                    source=newSource)

                    print(paste(input$ROInamePattern," ROI created, ",toadd," from ",input$selectROItoExtractPattern," ROI",sep=""))
            
                  }else if (length(extracted)==0 & !is.null(extracted)){
                    sendSweetAlert(
                      session = session,
                      title = "Empty ROI produced",
                      text = paste("pattern '",pattern,"' is not present in any range of ROI ",input$selectROItoExtractPattern,sep=""),
                      type = "error"
                    ) 
                  #if it's NULL, it means that any chromosome format is not recognised in BSgenome DB                    
                  }else if (length(extracted)==0 & is.null(extracted)){
                    sendSweetAlert(
                      session = session,
                      title = "chromosome names format not ok",
                      text = paste("Chromosome names of the ROI does not match those within BSgenome database. Check the format of the
                                  ROI, using 'get ROI' tab, looking at the correct chromosome names",sep=""),
                      type = "error"
                    )  
                  }

                }else{
                  sendSweetAlert(
                    session = session,
                    title = "ROI already present",
                    text = "A ROI with the same name already exists; use a different name",
                    type = "error"
                  )                  
                }
              }else{
                sendSweetAlert(
                  session = session,
                  title = "Bad ROI name",
                  text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
                  type = "error"
                )                
              }
            }else{
              sendSweetAlert(
                session = session,
                title = "Bad ROI name",
                text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
                type = "error"
              )             
            }
          
    
          }else{
            sendSweetAlert(
              session = session,
              title = "Pattern not valid",
              text = "The sequence pattern is not valid. Yu have to be consistent with IUPAC nomenclature (available letters are: A,C,G,T,U,R,Y,S,W,K,M,B,D,H,V,N,.,-)",
              type = "error"
            )             
          }

        }else{
          sendSweetAlert(
            session = session,
            title = "Missing Pattern",
            text = "You must choose a pattern to search (a sequence string)",
            type = "error"
          )           
        }

      }else{
        sendSweetAlert(
          session = session,
          title = "Missing BSgenome",
          text = paste("You must have ",BSstring," library for sequences. Download/Import it using the button"),
          type = "error"
        )         
      }
    }else{
      sendSweetAlert(
        session = session,
        title = "Missing database",
        text = "You must select a genome assembly. Go to 'Assembly' section to import a genome assembly",
        type = "error"
      )       
    }

  }else{
    sendSweetAlert(
      session = session,
      title = "Missing Elements",
      text = "It seems that no ROI and no annotation is present in the current session",
      type = "error"
    )      
  }
})












####################################################################################
####################################################################################
#view ROI tab
####################################################################################
####################################################################################


###react to help buttons:
#box help
observeEvent(input$msg_quickviewROIs, {
  boxHelpServer(msg_quickviewROIs)
})
#help options for view ROI
observeEvent(input$help_BED_viewoptions, {
  boxHelpServer(help_BED_viewoptions)
})








#view ROI
# view stat (width, number) ranges
observe({
  set.seed(123)
  
  if (!isvalid(input$chooseROIvisualiz)|!isvalid(input$confirmviewROI)){
    output$viewROImaterial<-renderPlot ({NULL})
    output$saveviewpeaksROImaterial<-renderUI ({NULL})
    return()
  }

  nomi=unlist(lapply(isolate(ROIvariables$listROI),getName))
  selection=list()
  for (i in 1:length(input$confirmviewROI)){
      pos2=match(input$confirmviewROI[i],nomi)
      selection[[i]]=getRange(isolate(ROIvariables$listROI)[[pos2]])
  }
  lun=unlist(lapply(selection,length))
  #continue only if all ROIs selected have length >0
  if (!all(lun>0)){
    output$viewROImaterial<-renderPlot({NULL})
    output$saveviewpeaksROImaterial=renderUI({NULL})
    sendSweetAlert(
      session = session,
      title = "Missing ROI",
      text = "Select at least one ROI",
      type = "error"
    )
    return()
  }


  #here we have >=1 selected ROI => plot in any case density (with different colors) and barplot for peak number
  #have to set the densities and colors. then select the max(x) and max (y) for the blank plot
  listfordensity=list()
  listselected=list()
  listselectednames=list()
  for (i in 1:length(input$confirmviewROI)){
    pos=match(input$confirmviewROI[i],nomi)
    listselected[[i]]=getRange(isolate(ROIvariables$listROI)[[pos]])
    lungh=length(listselected[[i]])
    listselectednames[[i]]=paste(nomi[pos]," (",lungh,")",sep="")
    x=listselected[[i]]
    listfordensity[[i]]=density(log2(width( x  )))
  }

  n=length(input$confirmviewROI)
  colorsfordensity=colors_list[1:n]
  ROIvariables$colorsfordensity<-colorsfordensity
  ROIvariables$listfordensity=listfordensity
  ROIvariables$listselected=listselected
  ROIvariables$listselectednames=listselectednames
  if (input$chooseROIvisualiz=="ROIwidth"){
        
    if (!is.null(listfordensity) & length(listfordensity)>0){
      coords=lapply(listfordensity,function(k) {return(list(k$x,k$y))})
      coordsx=unlist(lapply(coords,"[[",1))
      coordsy=unlist(lapply(coords,"[[",2))
      xmax=max(coordsx)
      xmin=min(coordsx)
      ymax=max(coordsy)
      ymin=min(coordsy)
      leg=unlist(listselectednames)

      output$viewROImaterial<-renderPlot ( {
        par(mar=c(4,4,4,1))
        plot(1, type="n", xlab="log2 width", ylab="probability", xlim=c(xmin, xmax), ylim=c(ymin, ymax),main="frequency plot width")
        for (i in 1:length(listfordensity)){
          lines(listfordensity[[i]],col=colorsfordensity[[i]],lwd=3)
        }
        legend("topright",leg,col=unlist(colorsfordensity),lty=1,lwd=3,bty = "n")
      })
    }else{
      output$viewROImaterial<-renderPlot ( {
        NULL
        #plot(1, type="n", xlab="log2 width", ylab="density",main="no ROI selected")
      })
    }

    
    output$saveviewpeaksROImaterial=renderUI({downloadButton('saveviewpeakswidthROIbutton', 'Get PDF')})

  }else{
    output$viewROImaterial<-renderPlot({
      par(mar=c(15,5,1,1))
      peaknumberlist=unlist(lapply(listselected,length))
      names(peaknumberlist)=listselectednames
      barplot (peaknumberlist,las=2,col=unlist(colorsfordensity))
      mtext(side=2, text="Intervals number", line=4)
      
    })
    #download button for PDF peak number 
    output$saveviewpeaksROImaterial=renderUI({downloadButton('saveviewpeaksnumberROIbutton', 'Get PDF')})        
  }


})




#observer for PDF download of peak number
output$saveviewpeaksnumberROIbutton<- downloadHandler(
  filename=function() {
      paste('Barplot_interval_number.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    par(mar=c(15,5,1,3))
    peaknumberlist=unlist(lapply(ROIvariables$listselected,length))
    names(peaknumberlist)=ROIvariables$listselectednames
    barplot (peaknumberlist,las=2,col=unlist(ROIvariables$colorsfordensity))
    mtext(side=2, text="Number of intervals", line=4)
    dev.off()
  } 
)

output$saveviewpeakswidthROIbutton<-downloadHandler(
  filename=function() {
      paste('Distribution_interval_width.pdf', sep='')
  },
  content=function(file) {
    #ROIcariables$listfordensity is not null, because the call to pdf button is inside
    #the "if" above here
    pdf(file)
    par(mar=c(5,5,4,4))
    coords=lapply(ROIvariables$listfordensity,function(k) {return(list(k$x,k$y))})
    coordsx=unlist(lapply(coords,"[[",1))
    coordsy=unlist(lapply(coords,"[[",2))
    xmax=max(coordsx)
    xmin=min(coordsx)
    ymax=max(coordsy)
    ymin=min(coordsy)
    leg=unlist(ROIvariables$listselectednames)
    plot(1, type="n", xlab="log2 width", ylab="probability", xlim=c(xmin, xmax), ylim=c(ymin, ymax),main="frequency plot width")
    for (i in 1:length(ROIvariables$listfordensity)){
      lines(ROIvariables$listfordensity[[i]],col=ROIvariables$colorsfordensity[[i]],lwd=3)
    }
    legend("topright",leg,col=unlist(ROIvariables$colorsfordensity),lty=1,lwd=3,bty = "n")

    dev.off()
  } 


)








##############################################################################
##############################################################################
##############################################################################
# observers for buttons for preview of the ROIs tables or genes associated
# and get ROI table or genelist
##############################################################################
##############################################################################
##############################################################################

###react to help buttons:
#ROI selection
observeEvent(input$msg_getRois_BOX, {
  boxHelpServer(msg_getRois_BOX)
})

#help each range option
observeEvent(input$help_BED_getroi_eachGR, {
  boxHelpServer(help_BED_getroi_eachGR)
})
#help genomic window option
observeEvent(input$help_BED_getroi_genomicWindow, {
  boxHelpServer(help_BED_getroi_genomicWindow)
})
#help genomic window value
observeEvent(input$help_BED_getroi_windowvalue, {
  boxHelpServer(help_BED_getroi_windowvalue)
})




observeEvent(input$showdataframeROI,{
  if(length(ROIvariables$listROI)>0 & length(input$listgetROI)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$listgetROI,nomi)
    ROI=ROIvariables$listROI[[pos]]
    range=getRange(ROI)
    ####newenrichimplementation####
    #we have to find maximum enrichment point. Raw cov w/o normalization is ok
    bams=Enrichlist$rawcoverage[[pos]]
    keys=Enrichlist$decryptkey[[pos]]
    normvals=Enrichlist$normfactlist[[pos]]
    ################################
    #build the dataframe (df) to show according to the various checks in the input
    #contained in the input$ROIoptionsToView (set in the updateROI.R)
    
    df=data.frame(matrix(nrow=length(range), ncol=0))
    if(!is.null(input$ROIoptionsToViewRANGE)){
      #append range to the dataframe:
      df2=data.frame(chr=seqnames(range),start=start(range),end=end(range),strand=strand(range))
      df=cbind(df,df2)
    }else{
      #no range
    }

    if(!is.null(input$ROIoptionsToViewMETADATA)){
      #append selected annotations to dataframe
      em=as.data.frame(elementMetadata(range))
      pos=match(input$ROIoptionsToViewMETADATA,colnames(em))
      df2=em[,pos,drop=FALSE]
      colnames(df2)=colnames(em)[pos]
      df=cbind(df,df2)
    }else{
      #no metadata
    }
    
    if(!is.null(input$ROIoptionsToViewENRICHMENTS)&length(bams)>0){
      #append selected enrichments to dataframe

      #bams=getBAMlist(ROI)
      bamnames=names(bams)
      pos=match(input$ROIoptionsToViewENRICHMENTS,bamnames)
      bams_selected_tmp=bams[pos]
      keys_selected=keys[pos]
      #summed values must be normalized for the normalization factor
      normvals_selected=normvals[pos]
      if(nc==1){
        megalist=lapply(1:length(bams_selected_tmp),function(i) {

          summ=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bams_selected_tmp[[i]],Nbins=1,Snorm=FALSE,key=keys_selected[[i]],norm_factor=normvals_selected[[i]])

          return(summ)          
        })
      }else{
        decision=0
        tryCatch({
          megalist=mclapply(1:length(bams_selected_tmp),function(i) {
            summ=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bams_selected_tmp[[i]],Nbins=1,Snorm=FALSE,key=keys_selected[[i]],norm_factor=normvals_selected[[i]])
            return(summ)          
          },mc.cores=nc)

          decision=1  
        },
        warning = function( w ){
          print("Warning: using single core for memory limit...")
        },
        error = function( err ){
          print("Warning: using single core for memory limit...")
        })

        if(decision==0){
          megalist=lapply(1:length(bams_selected_tmp),function(i) {
            summ=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bams_selected_tmp[[i]],Nbins=1,Snorm=FALSE,key=keys_selected[[i]],norm_factor=normvals_selected[[i]])
            return(summ)          
          })            
        }
      }

      #append to dataframe the sum of the enrichments
      for (i in 1:length(megalist)){
        name=names(bams_selected_tmp)[i]
        df[,ncol(df)+1]=as.character(megalist[[i]])
        colnames(df)[ncol(df)]=paste("enrichment_",name,sep="")
      } 


    }else{
      #no enrichment
    }

    if(ncol(df)>0){
      print(paste("Viewing the ROI",input$listgetROI))
      names(df)=gsub(" ","_",names(df))
      tosave$datatableROI<-df

      if(ncol(df)==1){
        #one single element: show txt instead of datatable
        output$previewROItodownload<-renderUI({
          textAreaInput("genelistROIhead",label=NULL,value=paste(as.character(df[,1]),collapse="\n"),height=400)
        }) 
        output$previewROItodownloadbutton<-renderUI({
          downloadButton('saveROIdatatable', 'Download')
          #shinySaveButton("saveROIdatatable","Download",NULL,filetype=list(rds="xls"))
        }) 
      }else{
        #put in previewROItodownload the table

        output$tableROIhead <- renderDataTable({
          as.matrix(df)
        },options=list(pagingType = "simple",pageLength = 10,lengthChange=FALSE,searching=TRUE,autowidth=TRUE))

        output$previewROItodownload<-renderUI({
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 650px; max-width: 1000px; background-color: #ffffff;",
            dataTableOutput("tableROIhead")
          )
          
          #shinySaveButton("saveROIdatatable","Download",NULL,filetype=list(rds="xls"))
        })
        output$previewROItodownloadbutton<-renderUI({
          downloadButton('saveROIdatatable', 'Download')
          #shinySaveButton("saveROIdatatable","Download",NULL,filetype=list(rds="xls"))
        })        
      }

            
    }else{
      #no data to show for this ROI
      output$previewROItodownload<-renderUI({
        paste("Select at least one parameter of the ROI...")
      })
      output$previewROItodownloadbutton<-renderUI({
        NULL
      })
      tosave$datatableROI=NULL 
      sendSweetAlert(
        session = session,
        title = "Missing data",
        text = "Select at least one parameter to show for this ROI",
        type = "error"
      )     
    }



  }else{
    output$previewROItodownload<-renderUI({
      paste("Select a ROI...")
    })
    output$previewROItodownloadbutton<-renderUI({
      NULL
    })
    tosave$datatableROI=NULL
    sendSweetAlert(
      session = session,
      title = "Missing ROI",
      text = "Select a ROI",
      type = "error"
    ) 
  }
},ignoreInit=TRUE)




#observer for the annotation of ROI in genomic window
observeEvent(input$showgenelistWindowROI,{
  #check on the existence of ROIs
  if(length(ROIvariables$listROI)>0 & length(input$listgetROI)>0){
    #check on the valid window selected
    if (input$windowAnnotateGenes>0 & input$windowAnnotateGenes<200000 & isvalid(input$windowAnnotateGenes)){
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      
      #check if promoters exist:
      if("promoters" %in% nomi){

        pos=match(input$listgetROI,nomi)
        ROI=ROIvariables$listROI[[pos]]
        range=getRange(ROI)
        pospromo=match("promoters",nomi)
        ROIpromo=ROIvariables$listROI[[pospromo]]
        fixedpromo=getFixed(ROIpromo)
        #given the range, center on the midpoint, take promoters and 
        range=resize(range,width=1,fix="center")
        #extract fixed point of promoters and resize of double the input window:
        #for example, if 10k window from the range midpoint is equal to do:
        #range midpoints:   -                      -         
        #                         10k  10k      
        #promoters:             |----|----|       |----|----|
        promo=suppressWarnings(unique(resize(fixedpromo,width=input$windowAnnotateGenes*2,fix="center")))
        ov=suppressWarnings(countOverlaps(promo,range))
        #WARNING: more than one midpoint of grange can be associated with the same
        # promoter!!
        #extract gene id of overlapping promoters with that window
        promo_selected=promo[ov>0]

        emd=as.data.frame(elementMetadata(promo_selected))
        pos=match(input$IDorsymbolwindow,colnames(emd))

        toshow=emd[,pos,drop=FALSE]
        toshow=unique(toshow)
        colnames(toshow)=NULL
        tosave$genelistROIwindow<-toshow


        output$previewROItodownload<-renderUI({
          textAreaInput("windowROIhead",label=NULL,value=paste(as.character(toshow[,1]),collapse="\n"),height=400)
        })   
        output$previewROItodownloadbutton<-renderUI({
          downloadButton('saveROIwindow', 'Download')
          #shinySaveButton("saveROIgenelist","Download",NULL,filetype=list(rds="txt"))
        }) 
        print(paste("Viewing genes annotated inside window to ROI",input$listgetROI))

      }else{
        output$previewROItodownload<-renderUI({
          paste("require promoters (choose an assembly for this)")
        })  
        output$previewROItodownloadbutton<-renderUI({
          NULL
        }) 
        tosave$genelistROIwindow<-NULL       
        sendSweetAlert(
          session = session,
          title = "Annotated elements not found",
          text = "Promoters of a specific genome assembly not found, but you need them: go to 'Assembly' section and choose a genome assembly",
          type = "error"
        ) 

      }
    }else{
      output$previewROItodownload<-renderUI({
        paste("Window not valid (==0 or > 200000 or missing)...")
      })  
      output$previewROItodownloadbutton<-renderUI({
        NULL
      }) 
      tosave$genelistROIwindow<-NULL
      sendSweetAlert(
        session = session,
        title = "Window not valid",
        text = "The genomic window must be between 0 and 200000",
        type = "error"
      )

    }
  }else{
    output$previewROItodownload<-renderUI({
      paste("Select a ROI...")
    }) 
    output$previewROItodownloadbutton<-renderUI({
      NULL
    }) 
    tosave$genelistROIwindow<-NULL
    sendSweetAlert(
      session = session,
      title = "Missing ROI",
      text = "Select a ROI",
      type = "error"
    ) 
  }
},ignoreInit=TRUE)






#view /download notes for ROI is in MANIPULATE roi advanced menu

#here, add view/edit notes of a ROI (the source)
observeEvent(input$saveNotesROI,{
  #check on the existence of ROIs
  if(length(ROIvariables$listROI)>0 & length(input$selectROItoEditNotes)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoEditNotes,nomi)
    ROI=ROIvariables$listROI[[pos]]
    #get value from input$ROInotes (text) and put as source for the ROI named by input$selectROItoEditNotes
    newsource=input$ROInotes
    ROI@source=as.list(newsource)
    ROIvariables$listROI[[pos]]=ROI
  }else{
    sendSweetAlert(
      session = session,
      title = "ROI not selected",
      text = "You didn't select a ROI",
      type = "error"
    )    
  }
},ignoreInit=TRUE)

#here, download text file of the notes of the ROI (source)
output$downloadNotesROI<- downloadHandler(
  filename=function() {paste('ROInotes.txt', sep='')},
  content=function(file) {
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoEditNotes,nomi)
    ROI=ROIvariables$listROI[[pos]]
    source_ROI=getSource(ROI)
    #download the source of the ROI saved, NOT the text actually written on the textArea    
    print("Downloading ROI source")
    fileConn<-file(file)

    writeLines(unlist(source_ROI), fileConn)
    close(fileConn)
  }
)






#observer for download of peaks
output$saveROIdatatable<- downloadHandler(
  filename=function() {paste('ROI.xls', sep='')},
  content=function(file) {
    print("Downloading ROI table")
    write.table(tosave$datatableROI,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  }
)
#observer for download of gene list
output$saveROIgenelist<- downloadHandler(
  filename=function() {paste('genelist.txt', sep='')},
  content=function(file) {
    print("Downloading nearest genes annotated")
    write.table(tosave$genelistROI,file=file,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE   ) 
  }
)
#observer for download of gene list from window
output$saveROIwindow<- downloadHandler(
  filename=function() {paste('genelist_window.txt', sep='')},
  content=function(file) {
    print("Downloading genes annotated in window")
    write.table(tosave$genelistROIwindow,file=file,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE   ) 
  }
)



####################################
#associate enrichments
####################################

###react to help buttons:
#associate/remove enrichments
observeEvent(input$msg_associateEnrichments_associateRemove, {
  boxHelpServer(msg_associateEnrichments_associateRemove)
})
#rename/remove enrichments
observeEvent(input$msg_renameEnrichments_renameOrderEnrichments, {
  boxHelpServer(msg_renameEnrichments_renameOrderEnrichments)
})





#respond to click in Deassociate BAM section (input$confirmBAMDeassociate), input$selectBAMtoDeassociate
observeEvent(input$confirmBAMDeassociate,{
  if (length(ROIvariables$listROI)>0 & length(input$selectBAMtoDeassociate)>0 & isvalid(input$selectROItoBAMassociate)){
    #check on length of Range is nonsense, because if length range ==0, there aren't bam files assoc.
    #get names of BAMlist in ROI selected 
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoBAMassociate,nomi)
    #ROIs=ROIvariables$listROI[pos]

    #for each ROI selected, get BAMlist and remove, if present, the selected enrichment files
    for(i in 1:length(pos)){
      #bams=getBAMlist(ROIvariables$listROI[[pos[i]]])

      ####newenrichimplementation####
      tempbam=Enrichlist$rawcoverage
      tempkey=Enrichlist$decryptkey
      tempnorm=Enrichlist$normfactlist
      bams=tempbam[[pos[i]]]
      keys=tempkey[[pos[i]]]
      normvals=tempnorm[[pos[i]]]
      ###############################  

      bamnames=names(bams)  
      pos2=match(input$selectBAMtoDeassociate,bamnames)
      #remove NA (not all selected BAMs are present in all ROIs)
      pos2=pos2[!is.na(pos2)]

      ####newenrichimplementation####
      bams[pos2]<-NULL
      keys[pos2]<-NULL
      normvals[pos2]<-NULL
      if(is.null(bams)){
        bams=list()
        keys=list()
        normvals=list()
      }
      tempbam[[pos[i]]]=bams
      tempkey[[pos[i]]]=keys
      tempnorm[[pos[i]]]=normvals
      ###############################

      #if NULL, recreate empty BAMlist
      # if (is.null(bams)){
      #   bams=list()
      # } 
      ####newenrichimplementation####
      Enrichlist$rawcoverage<-NULL
      Enrichlist$decryptkey<-NULL
      Enrichlist$normfactlist<-NULL
      gc()

      Enrichlist$rawcoverage=tempbam
      Enrichlist$decryptkey=tempkey
      Enrichlist$normfactlist=tempnorm


      ###############################  
      #ROIvariables$listROI[[pos[i]]]=setBAMlist(ROIvariables$listROI[[pos[i]]],bams)   
      for (i in input$selectBAMtoDeassociate){
        print(paste("Removed",i,"from ROI",input$selectROItoBAMassociate))
      }
      
      print("forcing gc...")
      print(gc())

    }
  }else{
    sendSweetAlert(
      session = session,
      title = "Missing selection",
      text = "Select at least one ROI and an enrichment file to remove from it",
      type = "error"
    )     
  }
},ignoreInit=TRUE)
    





#react to button "input$confirmBAMassociate" to associate enrichment files:
observeEvent(input$confirmBAMassociate,{
  #input$selectROItoBAMassociate
  #input$selectBAMtoassociate
  if (length(ROIvariables$listROI)>0 & length(input$selectBAMtoassociate)>0 & length(input$selectROItoBAMassociate)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROItoBAMassociate,nomi)
    ROIs=ROIvariables$listROI[pos] 
    totallist=list()   
     
    #now for each of these ROIs selected,find BAM (among all those selected) to associate with:
    for(i in 1:length(ROIs)){
      ####newenrichimplementation####
      enrichpresent=names(Enrichlist$rawcoverage[pos][[i]])
      ###############################     
      allbams=names(BAMvariables$listBAM)
      remaining=setdiff(allbams,enrichpresent)

      #if no BAM associatable with that ROI, simply skip
      if(length(remaining)>0){
        #find the requested enrichment over the remaining enrichments
        remaining_selected=remaining[match(input$selectBAMtoassociate,remaining)]
        #remove those that are NA:
        remaining_selected=remaining_selected[!is.na(remaining_selected)]
        pos2=match(remaining_selected,names(BAMvariables$listBAM))
        paths=BAMvariables$listBAM[pos2]
        partlist=rep(list(ROIs[[i]]),length(paths))
        names(partlist)=paths
        totallist=c(totallist,partlist)
      }else{

      }


    }



    #HERE the check for feasibility considering the RAM available. Sum width of all genomic ranges * number of coverage
    ###########################################################################################
    totalBPtobecoveraged=sum(sapply(lapply(lapply(totallist,getRange),width),sum))
    #here, let's suppose an average compression of 2x from the total number of bp (1-byte encryption) (assumption!)
    totalRAMtobeused=totalBPtobecoveraged/2

    avRAM=availRAM()
    if (totalRAMtobeused/1000000>avRAM){
      #not enough RAM for final storage: exit with message
      print(paste("available RAM:",avRAM,"not enough. Predict to occupy about",totalRAMtobeused/1000000))
      sendSweetAlert(
        session = session,
        title = "Not enough memory",
        text = "The RAM available seems not enough for the coverage. Try to reduce the number of ROI or the number of enrichment files for the coverage, or remove enrichment files aready associated.",
        type = "error"
      )          
      return()   
    }
    ###########################################################################################


    if(!is.null(input$selectMethodForNorm)){
	    if(input$selectMethodForNorm=="librarysize"){
	      Normmethod="lib"
	      normalizer=NULL
	      normCtrl=NULL
	      normCtrlSpike=NULL      
	    }else if(input$selectMethodForNorm=="customnorm"){
	      if(length(input$customNormalizer)>0){
	        Normmethod="custom"
	        normalizer=BAMvariables$listBAM[match(input$customNormalizer,names(BAMvariables$listBAM))][[1]]
	        normCtrl=NULL
	        normCtrlSpike=NULL
	      }
	    }else if(input$selectMethodForNorm=="spikein") {
	      if(length(input$customNormalizer)>0 & length(input$ctrlNormalizer)>0 & input$ctrlSpikeinNormalizer>0){
	        Normmethod="spikein"
	        normalizer=BAMvariables$listBAM[match(input$customNormalizer,names(BAMvariables$listBAM))][[1]]
	        normCtrl=BAMvariables$listBAM[match(input$ctrlNormalizer,names(BAMvariables$listBAM))][[1]]
	        normCtrlSpike=BAMvariables$listBAM[match(input$ctrlSpikeinNormalizer,names(BAMvariables$listBAM))][[1]]
	      }
	    }else if(input$selectMethodForNorm=="nonorm"){
	      Normmethod="no"
	      normalizer=NULL
	      normCtrl=NULL
	      normCtrlSpike=NULL
	    }else{
	      print("something wrong, no normalization existing...")
	    }		    	
    }else{
    	#no norm. menu, => wig file => no normalization
	    Normmethod="no"
	    normalizer=NULL
	    normCtrl=NULL
	    normCtrlSpike=NULL    	
    }




    ##here use totallist to parallelize the code. Each iteration has path to enrichment and a ROI
    
    tryCatch({

      #parallelize only if system RAM is very high
      if(input$coresCoverage==1 | length(totallist)==1){
        print ("Computation in single core")
        finallist=lapply(1:length(totallist),function(i) {
          if(Normmethod=="lib"){
            tonorm=names(totallist)[i]
          }else{
          	#if no normaliz., custom or spikein, normalizers were defined as normalizer,normCtrl,normCtrlSpikein
            tonorm=normalizer
          }          
          #here put check if transcripts:
          if(getFlag(totallist[[i]])!="transcriptFlag"){
            rangetocov=getRange(totallist[[i]])
            #here use chunksforassociate function to determine the number of chunks
            nchunks=chunksforassociate(rangetocov)
            singlecover=GRbaseCoverage2(Object=rangetocov,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                      signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=nchunks)
            nfact=singlecover[[3]]
            decryptkey=singlecover[[2]]
            singlecover=singlecover[[1]]

            #if all(singlecover==0) -> change nomenclature
            if(verifyzerocov(singlecover)){
              print ("cov is 0s... converting nomelclature to NCBI for coverage...")
              temprange=convertNomenclatureGR(rangetocov,to="NCBI")
              #temproi=setRange(totallist[[i]],temprange)
              singlecover=GRbaseCoverage2(Object=temprange,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                      signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=nchunks)
              nfact=singlecover[[3]]
              decryptkey=singlecover[[2]]
              singlecover=singlecover[[1]]              
            }

            print (paste("coverage",names(totallist)[i]))
            return(list(singlecover,decryptkey,nfact))
          }else{
            rangetransc=getRange(totallist[[i]])
            thirtypercent=round((width(rangetransc)/10)*3)
            #try suppress warnings... BE CAREFUL... maybe it will trim some bases,
            #this could cause problems in the future...
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="start"))
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="end"))
            #here use chunksforassociate function to determine the number of chunks
            nchunks=chunksforassociate(rangetransc)
            singlecover=GRbaseCoverage2(Object=rangetransc, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
            												signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks=nchunks)
            nfact=singlecover[[3]]
            decryptkey=singlecover[[2]]
            singlecover=singlecover[[1]]
            
            #if all(singlecover==0) -> change nomenclature
            if(verifyzerocov(singlecover)){
              print ("cov is 0s... converting nomelclature to NCBI for coverage...")
              temprange=convertNomenclatureGR(rangetransc,to="NCBI")
              singlecover=GRbaseCoverage2(Object=temprange, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
            												signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks)
              nfact=singlecover[[3]]
              decryptkey=singlecover[[2]]
              singlecover=singlecover[[1]]
            }
            print (paste("coverage",names(totallist)[i]))
            return(list(singlecover,decryptkey,nfact))
          }
        })#,mc.cores=nc)  

      }else{
        if(input$coresCoverage>nc){
          newnc=nc
        }else{
          newnc=input$coresCoverage
        }
        print (paste("parallel computation with ",newnc," cores",sep=""))
        finallist=mclapply(1:length(totallist),function(i) {
          if(Normmethod=="lib"){
            tonorm=names(totallist)[i]
          }else{
          	#if no normaliz., custom or spikein, normalizers were defined as normalizer,normCtrl,normCtrlSpikein
            tonorm=normalizer
          }  
       
          #here put check if transcripts:
          if(getFlag(totallist[[i]])!="transcriptFlag"){
            rangetocov=getRange(totallist[[i]])
            singlecover=GRbaseCoverage2(Object=rangetocov,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                  signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=1)
            nfact=singlecover[[3]]
            decryptkey=singlecover[[2]]
            singlecover=singlecover[[1]]            
            #if all(singlecover==0) -> change nomenclature
            if(verifyzerocov(singlecover)){
              print ("cov is 0s... converting nomelclature to NCBI for coverage...")
              temprange=convertNomenclatureGR(rangetocov,to="NCBI")
              #temproi=setRange(totallist[[i]],temprange)
              singlecover=GRbaseCoverage2(Object=temprange,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                      signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=1)
              nfact=singlecover[[3]]
              decryptkey=singlecover[[2]]
              singlecover=singlecover[[1]]            
            }

            print (paste("coverage",names(totallist)[i]))
            return(list(singlecover,decryptkey,nfact))
          }else{
            rangetransc=getRange(totallist[[i]])
            thirtypercent=round((width(rangetransc)/10)*3)
            #try suppress warnings... BE CAREFUL... maybe it will trim some bases,
            #this could cause problems in the future...
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="start"))
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="end"))
            singlecover=GRbaseCoverage2(Object=rangetransc, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks=1)
            nfact=singlecover[[3]]
            decryptkey=singlecover[[2]]
            singlecover=singlecover[[1]]            
            #if all(singlecover==0) -> change nomenclature
            if(verifyzerocov(singlecover)){
              print ("cov is 0s... converting nomelclature to NCBI for coverage...")
              temprange=convertNomenclatureGR(rangetransc,to="NCBI")
              singlecover=GRbaseCoverage2(Object=temprange, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
            												signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks=1)
              nfact=singlecover[[3]]
              decryptkey=singlecover[[2]]
              singlecover=singlecover[[1]]
            }
            print (paste("coverage",names(totallist)[i]))
            return(list(singlecover,decryptkey,nfact))
          }
        },mc.cores=newnc)  

      }



      if(Normmethod=="lib"){
        msg="library-size normalized"
      }else if (Normmethod=="custom"){
        msg=paste("normalized for ",input$customNormalizer,sep="")
      }else if (Normmethod=="no"){
      	msg="not normalized"
      }else if (Normmethod=="spikein"){
      	msg=paste("normalized for ",input$customNormalizer," as sample spike-in, ",
      		input$ctrlNormalizer," as control and ",input$ctrlSpikeinNormalizer," as control spike-in",sep="")
      }

      
      ##now attrib to each ROI the correct enrichments just calculated
      nameROIs=sapply(totallist,getName)
      finallist2=lapply(finallist,"[[",1)
      finaldecryptkeylist=lapply(finallist,"[[",2)
      normfinallist=lapply(finallist,"[[",3)

      splittedresultlist=split(finallist2,nameROIs)
      splittedkeylist=split(finaldecryptkeylist,nameROIs)
      splittednormlist=split(normfinallist,nameROIs)
      rm(finallist2)
      splittedtotallist=split(totallist,nameROIs)
      #"for each ROI selected for enrich. computation..."
      for(i in 1:length(splittedresultlist)){
        
        #previouslist=getBAMlist(splittedtotallist[[i]][[1]])
        ####newenrichimplementation####
        nme=getName(splittedtotallist[[i]][[1]])
        roilocation=match(nme,names(Enrichlist$rawcoverage))
        previousrawvals=Enrichlist$rawcoverage[[roilocation]]
        previousdecryptkey=Enrichlist$decryptkey[[roilocation]]
        previousnormvals=Enrichlist$normfactlist[[roilocation]]
        ###############################

        realnames=names(splittedtotallist[[i]])
        #get from bamlist the name
        pos=match(realnames,BAMvariables$listBAM)
        toname=names(BAMvariables$listBAM)[pos]
        names(splittedresultlist[[i]])=names(splittedkeylist[[i]])=names(splittednormlist[[i]])=toname
        #here, check if some finallist2 was NULL.
        notnulls=!(sapply(splittedresultlist[[i]],is.null))
        #print message of something is NULL
        if(sum(notnulls)<length(toname)){
          sendSweetAlert(
            session = session,
            title = "Some problems",
            text = paste("Some enrichment files gave 'NULL' (",toname[!notnulls],") and were not associated. Please, try again with them",sep=""),
            type="warning"
          )          
        }


        splittedresultlist_notnull=splittedresultlist[[i]][notnulls]
        splittedkeylist_notnull=splittedkeylist[[i]][notnulls]
        splittednormlist_notnull=splittednormlist[[i]][notnulls]
        #finalfinallist=c(previouslist,splittedresultlist_notnull)


        ####newenrichimplementation####
        #to be modified
        finalrawvals=c(previousrawvals,splittedresultlist_notnull)
        finaldecryptkeylist=c(previousdecryptkey,splittedkeylist_notnull)
        finalnormvals=c(previousnormvals,splittednormlist_notnull)
        
        
        rm(previousrawvals)
        rm(previousnormvals)
        rm(previousdecryptkey)
        ###############################
        #which ROI?
        posROI=match(names(splittedresultlist)[i],nomi)

        # ROIvariables$listROI[[posROI]]=setBAMlist(ROIvariables$listROI[[posROI]],finalfinallist)
        ####newenrichimplementation####
        #fill enrichlist
        Enrichlist$rawcoverage[[posROI]]=finalrawvals
        Enrichlist$decryptkey[[posROI]]=finaldecryptkeylist
        Enrichlist$normfactlist[[posROI]]=finalnormvals
        ################################


        rm(finalrawvals)
        rm(finalnormvals)
        rm(finaldecryptkeylist)
        #rm(finalfinallist)
        rm(splittedresultlist_notnull)
        rm(splittedkeylist_notnull)
        rm(splittednormlist_notnull)
        
        for(k in toname){
          logvariables$msg[[length(logvariables$msg)+1]]=paste('Associated ',k,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,'<br>',sep="")
        }
        print(paste('Associated ',toname,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,sep=""))
      }
      rm(splittedresultlist)

      print("forcing gc...")
      print(gc())  
      #alert the user that the block of files has been associated corectly to ROIs
      sendSweetAlert(
        session = session,
        title = "Enrichments associated!",
        text = paste("Selected enrichment files has been associated to selected ROIs",sep=""),
        type = "success"
      ) 
      return()

    },warning = function( w ){
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Warning: Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Problems in association",
        #text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
        text=p(style="text-align: left;","Problems in association. This is potentially due to: ",
        tags$li(style="text-align: left;",style="text-align: left;","Enrichment file not found. Check the correct path."),	
        tags$li(style="text-align: left;",style="text-align: left;","The enrichment file is corrupted"),
        tags$li(style="text-align: left;",style="text-align: left;","Memory not sufficient. Options: 1) Try to set the number of cores=1 2) Increase the RAM 3) Save the working session, exit from R, re-open the program and load the working session again (should solve memory fragmentation)")), 
        type = "error"
      )
    },error = function( err ){
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Problems in association",
        #text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
        text=p(style="text-align: left;","Problems in association. This is potentially due to: ",
        tags$li(style="text-align: left;","Enrichment file not found. Check the correct path."),	
        tags$li(style="text-align: left;","The enrichment file is corrupted"),
        tags$li(style="text-align: left;","Memory not sufficient. Options: 1) Try to set the number of cores=1 2) Increase the RAM 3) Save the working session, exit from R, re-open the program and load the working session again (should solve memory fragmentation)")), 
        type = "error"
      )          
    }) 




  }else{
    sendSweetAlert(
      session = session,
      title = "Missing selection",
      text = "Select at least one ROI and an enrichment file to associate",
      type = "error"
    )  
  }


})








###react to help buttons:
#rename enrichments
observeEvent(input$msg_renameEnrichments_renameEnrichments, {
  boxHelpServer(msg_renameEnrichments_renameEnrichments)
})
#reorder enrichments
observeEvent(input$msg_renameEnrichments_reorderEnrichments, {
  boxHelpServer(msg_renameEnrichments_reorderEnrichments)
})




#observe to rename a BAM of a specific ROI
observeEvent(input$renameBAM,{
  if (!is.null(ROIvariables$listROI)& length(ROIvariables$listROI)>=1){
    #check if name is set
    if (nchar(input$newBAMname)>=1) {
      #get the ROI in which rename the BAM file:
      nomi=unlist(lapply(ROIvariables$listROI,getName))
      pos=match(input$selectROIforBAMrename,nomi)
      roi=ROIvariables$listROI[[pos]]
      ####newenrichimplementation####
      rawvals=Enrichlist$rawcoverage[[pos]]
      keyvals=Enrichlist$decryptkey[[pos]]
      normvals=Enrichlist$normfactlist[[pos]]
      ###############################
      getbam=names(rawvals)  
      #check if name does not match with another existing BAM and at this ROI have at least one BAM
      if (!input$newBAMname %in% getbam & length(getbam)>0 ){

        for (i in 1:length(ROIvariables$listROI)){
          bamselected=Enrichlist$rawcoverage[[i]]
          getbam=names(bamselected)
          if(!input$newBAMname %in% getbam & length(getbam)>0){
            pos2=match(input$selectedBAMtoRename,getbam)
            #if BAM found in the current ROI in the loop
            if(!is.na(pos2)){
              oldname=getbam[pos2]
              names(bamselected)[pos2]=input$newBAMname
              names(Enrichlist$rawcoverage[[i]])=names(bamselected)
              names(Enrichlist$decryptkey[[i]])=names(bamselected)
              names(Enrichlist$normfactlist[[i]])=  names(bamselected)            
              #ROIvariables$listROI[[i]]=setBAMlist(ROIvariables$listROI[[i]],bamselected)
            }
          }
        }
        logvariables$msg[[length(logvariables$msg)+1]]= paste('Renamed ',oldname,' ROI in ',input$newBAMname,'<br></font>',sep="")
        print(paste('Renamed ',oldname,' ROI in ',input$newBAMname,sep=""))
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newBAMname,' name already exist in enrichments list for that ROI, or enrichment not found...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Name already present",
          text = paste("Name '",input$newBAMname,"'  already exists in enrichments list for that ROI, or enrichment not found",sep=""),
          type = "error"
        )      
      }
    }    
  }

},ignoreInit=TRUE)




#observe to reorder a BAM of a specific ROI
observeEvent(input$reorderBAM,{
  if (!is.null(ROIvariables$listROI) & length(ROIvariables$listROI)>0){
    #find ROI selected and if there are BAM files associated with it
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$selectROIforBAMrename,nomi)
    roi=ROIvariables$listROI[[pos]]
    ####newenrichimplementation####
    rawvals=Enrichlist$rawcoverage[[pos]]
    keyvals=Enrichlist$decryptkey[[pos]]
    normvals=Enrichlist$normfactlist[[pos]]
    ###############################
    getbam=names(rawvals)  

    if (length(getbam)>1 ){
      allnumbers=as.character(1:length(getbam))
      listprovv=list()
      for (i in 1:length(getbam)){
        stringval=grep(paste("reorderoptionBAM",i,"$",sep=""),names(input),value=TRUE)
        #extract what is contained inside each cell
        listprovv[[i]]=input[[ stringval ]]
      }
       
      listprovv=as.numeric(unlist(listprovv))
      #check if new order comprises all the possible positions
      if (identical(unique(sort(listprovv)),unique(sort(as.numeric(allnumbers)))) ){
        # ord=order(listprovv)
        ord=listprovv
        #reorder BAM 
        rawvals=rawvals[order(ord)]
        keyvals=keyvals[order(ord)]
        normvals=normvals[order(ord)]


        Enrichlist$rawcoverage[[pos]]=rawvals
        Enrichlist$decryptkey[[pos]]=keyvals
        Enrichlist$normfactlist[[pos]]=normvals
        #ROIvariables$listROI[[pos]]=setBAMlist(ROIvariables$listROI[[pos]],rawvals)
        
        if (!identical(unique(listprovv),unique(allnumbers))){
          logvariables$msg[[length(logvariables$msg)+1]]= paste('enrichment files reordered...<br>',sep="")
          print("enrichment files reordered")
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
    
    }#else: cannot reorder something with length ==1
      
  }
},ignoreInit=TRUE)


















##############################################################################
##############################################################################
##############################################################################
##############################################################################
# predefined pipeline for ROI preparation for heatmap
##############################################################################
##############################################################################
##############################################################################
##############################################################################

#ROI combination
observeEvent(input$msg_PredefPipeline_parameters, {
  boxHelpServer(msg_PredefPipeline_parameters)
})


observeEvent(input$PrepareROIpredefPipeline,{
  #here, the inputs are: 
  #input$selectROIpredefPipeline -> starting ROI to prepare
  #input$quantileThreshPredefPipeline -> random fraction of ranges to take
  #input$choiceSummitPredefPipeline -> if do summit center (yes/no)
  #   input$selectBAMsummitPredefPipeline -> enrichment file to use for summit, if previous was "yes"
  #input$sliderUpstreamPredefPipeline, input$sliderDownstreamPredefPipeline -> window resize
  #input$enrichAllPredefPipeline -> enrichment files to associate to prepared ROI
  #input$ROInamePredefPipeline -> name of the new ROI



  ##################################################################################
  #checks
  ##################################################################################
  if(!isvalid(input$selectROIpredefPipeline)){
    #no ROI, do nothing!    
    return()
  }

  #general criteria
  if(is.null(ROIvariables$listROI) | length(ROIvariables$listROI)<1 | length(input$selectROIpredefPipeline)==0 | !isvalid(input$sliderUpstreamPredefPipeline) | !isvalid(input$sliderDownstreamPredefPipeline)){
    sendSweetAlert(
      session = session,
      title = "Missing ROI",
      text = "You have to select at least one ROI",
      type = "error"
    )          
    return()
  }

  #check if valid name
  nottobe=c("promoters","transcripts","TES")
  if (nchar(input$ROInamePredefPipeline)<1 | any(input$ROInamePredefPipeline == nottobe)) {
    sendSweetAlert(
      session = session,
        title = "Bad ROI name",
        text = "File name for the new ROI is missing, or you are trying to use reserved names: 'promoters','transcripts' or 'TES'",
        type = "error"
    )          
    return()    
  }

  nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
  if(grepl(nottobe2[1],input$ROInamePredefPipeline) | grepl(nottobe2[2],input$ROInamePredefPipeline) | grepl(nottobe2[3],input$ROInamePredefPipeline)){
    sendSweetAlert(
      session = session,
      title = "Bad ROI name",
      text = "New ROI cannot be named 'promoters_genelist_*','transcripts_genelist_*' or 'TES_genelist_*', because are reserved names",
      type = "error"
    )  
    return()      
  }

  #check if name existing
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  if(input$ROInamePredefPipeline %in% nomi){
    sendSweetAlert(
      session = session,
      title = "ROI already present",
      text = "A ROI with the same name already exists; use a different name",
      type = "error"
    ) 
    return()     
  }

  #check on transcript type?
  if (input$selectROIpredefPipeline=="transcripts" | grepl("transcripts_genelist_",input$selectROIpredefPipeline)){
    sendSweetAlert(
      session = session,
      title = "Bad choice",
      text = "Transcripts cannot be used, because they cannot be resized",
      type = "error"
    )     
    return() 
  }




  #if summit, check if enrichment exist
  if(input$choiceSummitPredefPipeline=="Yes"){
    if(!isvalid(input$selectBAMsummitPredefPipeline) | length(input$selectBAMsummitPredefPipeline)<1){
      sendSweetAlert(
        session = session,
        title = "Bad choice",
        text = "Summit detection selected, but no enrichment found...",
        type = "error"
      )     
      return() 
    }

    #check existence of file and .bai extension if BAM
    pos2=match(input$selectBAMsummitPredefPipeline,names(BAMvariables$listBAM))
    paths_summit=BAMvariables$listBAM[[pos2]]
    if(!file.exists(paths_summit)){
      sendSweetAlert(
        session = session,
        title = "Enrichment not found",
        text = "Enrichment file for summit detection not found. Maybe the link to that file changed or file has been removed from filesystem... try to reopen it",
        type = "error"
      )     
      return()       
    }

    #if enrichment for summit is a BAM, check if .bai exists with the same name    
    if ( substring(input$selectBAMsummitPredefPipeline,nchar(input$selectBAMsummitPredefPipeline)-3,nchar(input$selectBAMsummitPredefPipeline))==".bam"){
      #find if bam file has its own index (.bam.bai)
      baifile=paste(paths_summit,".bai",sep="")
      if (!file.exists(baifile)){
        sendSweetAlert(
          session = session,
          title = "BAM index not found",
          text = "BAM index of the selected BAM file for summit detection not found.",
          type = "error"
        )     
        return()           
      }
    }
  }




  #if some enrichments selected, check existence. They are not mandatory for digital heatmap
  if(isvalid(input$enrichAllPredefPipeline) & length(input$enrichAllPredefPipeline)>=1 ){
    for(i in 1:length(input$enrichAllPredefPipeline)){
      pos2=match(input$enrichAllPredefPipeline[i],names(BAMvariables$listBAM))
      currentenrichpath=BAMvariables$listBAM[[pos2]]
      if(!file.exists(currentenrichpath)){
        sendSweetAlert(
          session = session,
          title = "Enrichment not found",
          text = paste("Enrichment file ",input$enrichAllPredefPipeline[i]," not found. Maybe the link to that file changed or file has been removed from filesystem... try to reopen it",sep="") ,
          type = "error"
        )     
        return()           
      }

      #check .bai file if BAM file selected
      if ( substring(input$enrichAllPredefPipeline[i],nchar(input$enrichAllPredefPipeline[i])-3,nchar(input$enrichAllPredefPipeline[i]))==".bam"){
        #find if bam file has its own index (.bam.bai)
        baifile=paste(currentenrichpath,".bai",sep="")
        if (!file.exists(baifile)){
          sendSweetAlert(
            session = session,
            title = "BAM index not found",
            text = paste("BAM index of the BAM file ",input$enrichAllPredefPipeline[i]," not found.",sep=""),
            type = "error"
          )     
          return()           
        }
      }


    }
  }


  ##################################################################################
  #random sample
  ##################################################################################
  
  pos=match(input$selectROIpredefPipeline,nomi)
  roi=ROIvariables$listROI[[pos]]
  oldSource=getSource(roi)
  rangefromROI=getRange(roi)
  newflag=getFlag(roi)

  #check flag: transcript flag cannot undergo this pipeline, they cannot be resized!!
  if(newflag=="transcriptFlag"){
    sendSweetAlert(
      session = session,
      title = "Bad choice",
      text = "Transcripts cannot be used, because they cannot be resized around their midpoint",
      type = "error"
    )     
    return()     
  }
  
  
  
  maxx=length(getRange(roi))
  quant=unname(round(quantile(1:maxx,input$quantileThreshPredefPipeline)))
  selectedPositions=sort(sample(length(rangefromROI),quant,replace=FALSE))
  newrange=rangefromROI[selectedPositions]

  if(length(newrange)==0){
    sendSweetAlert(
      session = session,
      title = "Empty ROI produced",
      text = paste("Random sampling '",input$selectROIpredefPipeline,"' produced a ROI with length=0",sep=""),
      type = "error"
    )  
    return()    
  }


  #bams=getBAMlist(roi)
  ####newenrichimplementation####
  bams=Enrichlist$rawcoverage[[pos]]
  keys=Enrichlist$decryptkey[[pos]]
  normvals=Enrichlist$normfactlist[[pos]]
  ################################
  if(length(bams)>0){
    newBAMlist=list()
    newkeylist=list()
    for(i in 1:length(bams)){
      newBAMlist[[i]]=bams[[i]][selectedPositions]
      newkeylist[[i]]=keys[[i]][selectedPositions]
    }
    names(newBAMlist)=names(newkeylist)=names(bams)            
  }else{
    newBAMlist=list()
    newkeylist=list()
  }
            
  newfix=getFixed(roi)[selectedPositions]

  perc=round(quant/length(rangefromROI)*100,2)
  toadd=paste("ROI sampled (kept ",quant,"/",length(rangefromROI),"ranges,",perc,"%)")
  newSource=c(oldSource,list(toadd))

  ####newenrichimplementation####
  Enrichlist$rawcoverage[[input$ROInamePredefPipeline]]=newBAMlist
  Enrichlist$decryptkey[[input$ROInamePredefPipeline]]=newkeylist
  Enrichlist$normfactlist[[input$ROInamePredefPipeline]]=normvals
  ################################  
  ROI_after_sample=new("RegionOfInterest",
                                            name=input$ROInamePredefPipeline,
                                            range=newrange,
                                            fixed=newfix,
                                            BAMlist=list(),
                                            flag=newflag,
                                            source=newSource)

  # logvariables$msg[[length(logvariables$msg)+1]]= paste('Sampling of ',input$selectROIpredefPipeline,' (kept ',perc,'%, ',quant,'/',length(rangefromROI),'ranges)<br>',sep="")  
  # print(paste('Sampling of ',input$selectROIpredefPipeline,' (kept ',perc,'%, ',quant,'/',length(rangefromROI),'ranges)',sep=""))          


  ##################################################################################
  #summit (if yes)
  ##################################################################################
  
	totalBPtobecoveraged=sum(width(getRange(ROI_after_sample)))
	#here, let's suppose an average compression of 2x from the total number of bp (1-byte encryption) (assumption!)
	totalRAMtobeused=totalBPtobecoveraged/2
	avRAM=availRAM()
	if (totalRAMtobeused/1000000>avRAM){
	  #not enough RAM for final storage: exit with message
	  print(paste("available RAM:",avRAM,"not enough. Predict to occupy about",totalRAMtobeused/1000000))
	  sendSweetAlert(
	    session = session,
	    title = "Not enough memory",
	    text = "The RAM available seems not enough for calculating the summit.",
	    type = "error"
	  )          
	  return()   
	}



  #with elements produced in previous block, do the summit
  #check if enrichment file of the summit already associated
  bams_names=names(newBAMlist)
  if(input$choiceSummitPredefPipeline=="Yes"){
    if(!input$selectBAMsummitPredefPipeline%in%bams_names){
      rangetocov=getRange(ROI_after_sample)
      nchunks=chunksforassociate(rangetocov)
      #associate BAM/WIG, no need to normalize, it's just summit. Useless for the same enrich. file
      singlecover=GRbaseCoverage2(Object=rangetocov,signalfile=paths_summit,signalfileNorm=NULL,
      						signalControl=NULL,signalControlSpike=NULL,multiplFactor=1e+06,nchunks=nchunks)
      #we don't need to keep any norm fact here, it's just the summit
      key=singlecover[[2]]
      singlecover=singlecover[[1]]
      
      if(verifyzerocov(singlecover)){
        print ("cov is 0s... converting nomelclature to NCBI for summit cov detection...")
        temprange=convertNomenclatureGR(rangetocov,to="NCBI")
        #temproi=setRange(ROI_after_sample,temprange)
        singlecover=GRbaseCoverage2(Object=temprange,signalfile=paths_summit,signalfileNorm=NULL,signalControl=NULL,signalControlSpike=NULL,multiplFactor=1e+06,nchunks=nchunks)
        key=singlecover[[2]]
        singlecover=singlecover[[1]]

      }
    }else{
      #do nothing, enrichment already in ROI. Use BAM already present to center on summit
      pos2=match(input$selectBAMsummitPredefPipeline,bams_names)
      singlecover=newBAMlist[[pos2]]
      key=newkeylist[[pos2]]
    }
    #newrange and singlecover must have same length (both derive from correctly subsampled ROI)
    newrange=summitFromBaseCoverage(Object=newrange,baseCoverageOutput=singlecover,keys=key)
    #newrange=summitFromBaseCoverage(Object=newrange,baseCoverageOutput=singlecover)
    oldSource=newSource
    toadd=paste("Centered on summit using",input$selectBAMsummitPredefPipeline,"enrichment (range width=1)")
    
    #re-annotate everything, because center on summit
    if(length(DATABASEvariables$currentASSEMBLY)>0){
      if(DATABASEvariables$currentASSEMBLY!=FALSE){
        #we have a database. So extract the fix of promoters and apply the distanceFromTSS3 function
        #for this newly created range
        nomi=unlist(lapply(ROIvariables$listROI,getName))
        pos_promo=match("promoters",nomi)
        promo=ROIvariables$listROI[[pos_promo]]
        fix_promoters=getFixed(promo)
        annotatedpart=suppressWarnings(distanceFromTSS3(Object=newrange,Tss=fix_promoters,criterion="midpoint"))
        elementMetadata(newrange)=annotatedpart
      }
    }    
    ####newenrichimplementation####
    Enrichlist$rawcoverage[[input$ROInamePredefPipeline]]=list()
    Enrichlist$decryptkey[[input$ROInamePredefPipeline]]=list()
    Enrichlist$normfactlist[[input$ROInamePredefPipeline]]=list()
    ################################ 
    ROI_after_summit=new("RegionOfInterest",
                                            name=input$ROInamePredefPipeline,
                                            range=newrange,
                                            fixed=newrange,
                                            BAMlist=list(),
                                            flag="normalFlag",
                                            source=c(oldSource,list(toadd))

    )

  }else{
    #after this passage, ROI won't change 
    ROI_after_summit=ROI_after_sample
  }


  ##################################################################################
  #resize
  ##################################################################################

  #get all starting info from ROI_after_summit
  roi=getRange(ROI_after_summit)
  oldSource=getSource(ROI_after_summit)
  fix=start(getFixed(ROI_after_summit))
  #bamlist=getBAMlist(ROI_after_summit)
  bamlist=Enrichlist$rawcoverage[[input$ROInamePredefPipeline]]
  keylist=Enrichlist$decryptkey[[input$ROInamePredefPipeline]]
  normfact=Enrichlist$normfactlist[[input$ROInamePredefPipeline]]
  #split positive or undetermined strand from negative strand
  pos_positive=as.logical(!strand(roi)=="-")
  pos_negative=as.logical(strand(roi)=="-")
  #remove ranges that would have negative starts with new width in positive strand
  idx=(fix-input$sliderUpstreamPredefPipeline)>0
  idx_positive_and_good=pos_positive & idx
  #negative strand
  idx=(fix-input$sliderDownstreamPredefPipeline)>0
  idx_negative_and_good=pos_negative & idx
  idx_positive_or_negative_good=idx_positive_and_good|idx_negative_and_good
  roi=roi[idx_positive_or_negative_good]
  fix=fix[idx_positive_or_negative_good]
  bamlist=lapply(bamlist,function(i) {i[idx_positive_or_negative_good]})
  keylist=lapply(keylist,function(i) {i[idx_positive_or_negative_good]})

  pos_positive=as.logical(!strand(roi)=="-")
  pos_negative=as.logical(strand(roi)=="-")
  roi_positive=roi[pos_positive]
  roi_negative=roi[pos_negative]
  fix_positive=fix[pos_positive]
  fix_negative=fix[pos_negative]
  bamlist_positive=lapply(bamlist,function(i) {i[pos_positive]})
  keylist_positive=lapply(keylist,function(i) {i[pos_positive]})
  bamlist_negative=lapply(bamlist,function(i) {i[pos_negative]})
  keylist_negative=lapply(keylist,function(i) {i[pos_negative]})

  #resize positivee strand:
  oldstarts_positive=start(roi_positive)
  oldends_positive=end(roi_positive)
  oldstarts_negative=start(roi_negative)
  oldends_negative=end(roi_negative)
  
  newstart=input$sliderUpstreamPredefPipeline
  newend=input$sliderDownstreamPredefPipeline

  toadd=paste("resized to window [-",input$sliderUpstreamPredefPipeline,"; +",input$sliderDownstreamPredefPipeline,"] from the centre",sep="")

  #resize plus strand
  start(roi_positive)=fix_positive-newstart
  end(roi_positive)=fix_positive+newend
  #resize negative strand
  start(roi_negative)=fix_negative-newend
  end(roi_negative)=fix_negative+newstart
  #join in the original roi
  roi[pos_positive]=roi_positive
  roi[pos_negative]=roi_negative
  newfix=getFixed(ROI_after_summit)[idx_positive_or_negative_good]
  newflag=getFlag(ROI_after_summit)
  #ROI_after_summit=setRange(ROI_after_summit,roi)
  #ROI_after_summit=setFix(ROI_after_summit,newfix)
  newrange=roi

  newBAMlist=list()
  newkeylist=list()
  normfact=list()
  #}

  newSource=c(oldSource,list(toadd))
  ####newenrichimplementation####
  Enrichlist$rawcoverage[[input$ROInamePredefPipeline]]=newBAMlist
  Enrichlist$decryptkey[[input$ROInamePredefPipeline]]=newkeylist
  Enrichlist$normfactlist[[input$ROInamePredefPipeline]]=normfact
  ################################   
  ROI_after_resize=new("RegionOfInterest",
                                name=input$ROInamePredefPipeline,
                                range=newrange,
                                fixed=newfix,
                                BAMlist=list(),
                                flag=newflag,
                                source=newSource) 



  ##################################################################################
  #enrichment association
  ##################################################################################
  #ncores=input$coresPredefPipeline
  ncores=1
  enrichments_selected=input$enrichAllPredefPipeline
  #enrichments_alreadyAssociated=getBAMlist(ROI_after_resize)
  enrichments_alreadyAssociated=Enrichlist$rawcoverage[[input$ROInamePredefPipeline]]
  previousdecryptkey=Enrichlist$decryptkey[[input$ROInamePredefPipeline]]
  previousnormvals =Enrichlist$normfactlist[[input$ROInamePredefPipeline]]
  
  #determine which enrichments are already associated
  enrichments_toAssociate=setdiff(enrichments_selected,names(enrichments_alreadyAssociated))
  pos2=match(enrichments_toAssociate,names(BAMvariables$listBAM))
  paths=BAMvariables$listBAM[pos2]  



	totalBPtobecoveraged=sum(width(getRange(ROI_after_resize)))
	#here, let's suppose an average compression of 2x from the total number of bp (1-byte encryption) (assumption!)
	totalRAMtobeused=totalBPtobecoveraged/2
	avRAM=availRAM()
	if (totalRAMtobeused/1000000>avRAM){
	  #not enough RAM for final storage: exit with message
	  print(paste("available RAM:",avRAM,"not enough. Predict to occupy about",totalRAMtobeused/1000000))
	  sendSweetAlert(
	    session = session,
	    title = "Not enough memory",
	    text = "The RAM available seems not enough for calculating enrichments.",
	    type = "error"
	  )          
	  return()   
	}

  #depending if only one operation, lapply or mclapply with n cores
  #by default, simple library-size normalization (except WIGs)
  #IMPORTANT: NEVER TANSCRIPT FLAG
  if(length(enrichments_selected)>=1){
    if(ncores==1 | length(enrichments_toAssociate)==1 | length(paths)==1){
      #single core (lapply)
      finallist=lapply(1:length(paths),function(i) {
        rangetocov=getRange(ROI_after_resize)
        nchunks=chunksforassociate(rangetocov)
        singlecover=GRbaseCoverage2(Object=rangetocov,signalfile=paths[[i]],signalfileNorm=paths[[i]],signalControl=NULL,signalControlSpike=NULL,nchunks=nchunks)
        normfacts=singlecover[[3]]
        decryptkey=singlecover[[2]]
        singlecover=singlecover[[1]]
        print(paste("Associating ",paths[[i]]," in single core",sep=""))
        if(verifyzerocov(singlecover)){
        	print ("cov is 0s... converting nomelclature to NCBI for cov...")
        	temprange=convertNomenclatureGR(rangetocov,to="NCBI")
        	#temproi=setRange(ROI_after_resize,temprange)
        	singlecover=GRbaseCoverage2(Object=temprange,signalfile=paths[[i]],signalfileNorm=paths[[i]],signalControl=NULL,signalControlSpike=NULL,nchunks=nchunks)
          normfacts=singlecover[[3]]
          decryptkey=singlecover[[2]]
          singlecover=singlecover[[1]]
        }
        return(list(singlecover,decryptkey,normfacts))
      })

    }else{
      #multicore (mclapply)
      finallist=mclapply(1:length(paths),function(i) {
        rangetocov=getRange(ROI_after_resize)
        singlecover=GRbaseCoverage2(Object=rangetocov,signalfile=paths[[i]],signalfileNorm=paths[[i]],signalControl=NULL,signalControlSpike=NULL,nchunks=1)
        normfacts=singlecover[[3]]
        decryptkey=singlecover[[2]]
        singlecover=singlecover[[1]]        
        print(paste("Associating ",paths[[i]]," in multi core",sep=""))
        if(verifyzerocov(singlecover)){
        	print ("cov is 0s... converting nomelclature to NCBI for cov...")
        	temprange=convertNomenclatureGR(rangetocov,to="NCBI")
        	#temproi=setRange(ROI_after_resize,temprange)
        	singlecover=GRbaseCoverage2(Object=temprange,signalfile=paths[[i]],signalfileNorm=paths[[i]],signalControl=NULL,signalControlSpike=NULL,nchunks=1)
          normfacts=singlecover[[3]]
          decryptkey=singlecover[[2]]
          singlecover=singlecover[[1]]        
        }
        return(list(singlecover,decryptkey,normfacts))
      },mc.cores=ncores)
    }

    #put names to new enrichments
    finallist2=lapply(finallist,"[[",1)
    decryptkeyfinallist=lapply(finallist,"[[",2)
    normfinallist=lapply(finallist,"[[",3)
    names(finallist2)=names(decryptkeyfinallist)=names(normfinallist)=names(paths)
    finalfinallist=c(enrichments_alreadyAssociated,finallist2)
    finaldecryptkeylist=c(previousdecryptkey,decryptkeyfinallist)
    finalnormvals=c(previousnormvals,normfinallist)

    ####newenrichimplementation####
    #fill enrichlist
    Enrichlist$rawcoverage[[input$ROInamePredefPipeline]]=finalfinallist
    Enrichlist$decryptkey[[input$ROInamePredefPipeline]]=finaldecryptkeylist
    Enrichlist$normfactlist[[input$ROInamePredefPipeline]]=finalnormvals

    ROI_after_associate=ROI_after_resize
    #ROI_after_associate=setBAMlist(ROI_after_resize,finalfinallist)
  }else{
    ROI_after_associate=ROI_after_resize
  }

  #logvariables$msg[[length(logvariables$msg)+1]]=paste('Associated ',k,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,'<br>',sep="")

  
  ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=ROI_after_associate
  #here, single big log message, summary of previous operations and append the new ROI prepared

  logvariables$msg[[length(logvariables$msg)+1]]= paste("Prepared ROI ",input$ROInamePredefPipeline," using the default heatmap pipeline starting from ",input$selectROIpredefPipeline," ROI: ",
        'sampling of ',input$selectROIpredefPipeline,' (kept ',perc,'%, ',quant,'/',length(rangefromROI),' ranges);<br>',
        "centered on summit using ",input$selectBAMsummitPredefPipeline," enrichment;<br>",
        'resized at [-',input$sliderUpstreamPredefPipeline,'; +',input$sliderDownstreamPredefPipeline,'] intervals from center.',sep="") 

  print(paste("Prepared ROI ",input$ROInamePredefPipeline," using the default heatmap pipeline starting from ",input$selectROIpredefPipeline," ROI: ",
        'sampling of ',input$selectROIpredefPipeline,' (kept ',perc,'%, ',quant,'/',length(rangefromROI),' ranges); ',
        "centered on summit using ",input$selectBAMsummitPredefPipeline," enrichment; ",
        'resized at [-',input$sliderUpstreamPredefPipeline,'; +',input$sliderDownstreamPredefPipeline,'] intervals from center',sep=""))
  
})
























##############################################################################
##############################################################################
##############################################################################
##############################################################################
# predefined pipeline for genelist preparation for metagene profile
##############################################################################
##############################################################################
##############################################################################
##############################################################################

#ROI combination
observeEvent(input$msg_GenelistPipeline_parameters, {
  boxHelpServer(msg_GenelistPipeline_parameters)
})


observeEvent(input$PrepareROIgenelistPipeline,{
  #here, the inputs are: 
  #input$selectROIGenelistPipeline -> starting genelists (triad of ROIs) to prepare
  #input$BAMforGenelistPipeline -> enrichments to asociate to genelists
  #input$PrepareROIgenelistPipeline -> button to execute



  ##################################################################################
  #checks (useless if great UI)
  ##################################################################################
  if(!isvalid(input$selectROIGenelistPipeline)){
    #no ROI, do nothing!    
    return()
  }

  #general criteria
  if(is.null(ROIvariables$listROI) | length(ROIvariables$listROI)<1 | length(input$selectROIGenelistPipeline)==0 ){
    sendSweetAlert(
      session = session,
      title = "Missing ROI",
      text = "You have to select at least one ROI",
      type = "error"
    )          
    return()
  }


  #check if enrichment selected (useless if great UI)
  if(!isvalid(input$BAMforGenelistPipeline)){
    sendSweetAlert(
      session = session,
      title = "Missing enrichment",
      text = "You have to select at least one enrichment",
      type = "error"
    )          
    return()
  }


  #if some enrichments selected, check existence. 
  if(isvalid(input$BAMforGenelistPipeline) & length(input$BAMforGenelistPipeline)>=1 ){
    for(i in 1:length(input$BAMforGenelistPipeline)){
      pos2=match(input$BAMforGenelistPipeline[i],names(BAMvariables$listBAM))
      currentenrichpath=BAMvariables$listBAM[[pos2]]
      if(!file.exists(currentenrichpath)){
        sendSweetAlert(
          session = session,
          title = "Enrichment not found",
          text = paste("Enrichment file ",input$BAMforGenelistPipeline[i]," not found. Maybe the link to that file changed or file has been removed from filesystem... try to reopen it",sep="") ,
          type = "error"
        )     
        return()           
      }

      #check .bai file if BAM file selected
      if ( substring(input$BAMforGenelistPipeline[i],nchar(input$BAMforGenelistPipeline[i])-3,nchar(input$BAMforGenelistPipeline[i]))==".bam"){
        #find if bam file has its own index (.bam.bai)
        baifile=paste(currentenrichpath,".bai",sep="")
        if (!file.exists(baifile)){
          sendSweetAlert(
            session = session,
            title = "BAM index not found",
            text = paste("BAM index of the BAM file ",input$BAMforGenelistPipeline[i]," not found.",sep=""),
            type = "error"
          )     
          return()           
        }
      }

    }
  }







  ######################################################################
  # execute
  ######################################################################


  nomi=unlist(lapply(ROIvariables$listROI,getName))
  pos=list()
  for(i in 1:length(input$selectROIGenelistPipeline)){
    selected=input$selectROIGenelistPipeline[i]
    pos_promoters= which(paste("promoters_genelist_",selected,sep="")==nomi)
    pos_transcripts= which(paste("transcripts_genelist_",selected,sep="")==nomi)
    pos_TES= which(paste("TES_genelist_",selected,sep="")==nomi)
    pos[[i]]=c(pos_promoters,pos_transcripts,pos_TES)
  }
  pos=Reduce(union, pos)
  ROIs=ROIvariables$listROI[pos] 
  totallist=list()   
  
  #now for each of these ROIs selected,find BAM (among all those selected) to associate with:
  for(i in 1:length(ROIs)){
    ####newenrichimplementation####
    enrichpresent=names(Enrichlist$rawcoverage[pos][[i]])
    ###############################     
    allbams=names(BAMvariables$listBAM)
    remaining=setdiff(allbams,enrichpresent)

    #if no BAM associatable with that ROI, simply skip
    if(length(remaining)>0){
      #find the requested enrichment over the remaining enrichments
      remaining_selected=remaining[match(input$BAMforGenelistPipeline,remaining)]
      #remove those that are NA:
      remaining_selected=remaining_selected[!is.na(remaining_selected)]
      pos2=match(remaining_selected,names(BAMvariables$listBAM))
      paths=BAMvariables$listBAM[pos2]
      partlist=rep(list(ROIs[[i]]),length(paths))
      names(partlist)=paths
      totallist=c(totallist,partlist)
    }else{

    }
  }

  
  #HERE the check for feasibility considering the RAM available. Sum width of all genomic ranges * number of coverage
  ###########################################################################################
  totalBPtobecoveraged=sum(sapply(lapply(lapply(totallist,getRange),width),sum))
  #here, let's suppose an average compression of 2x from the total number of bp (1-byte encryption) (assumption!)
  totalRAMtobeused=totalBPtobecoveraged/2

  avRAM=availRAM()
  if (totalRAMtobeused/1000000>avRAM){
    #not enough RAM for final storage: exit with message
    print(paste("available RAM:",avRAM,"not enough. Predict to occupy about",totalRAMtobeused/1000000))
    sendSweetAlert(
      session = session,
      title = "Not enough memory",
      text = "The RAM available seems not enough for the coverage. Try to reduce the number of ROI or the number of enrichment files for the coverage, or remove enrichment files aready associated.",
      type = "error"
    )          
    return()   
  }
  ###########################################################################################


  #######################
  # normalization, only lib size for pipeline for genelists
  # if you want custom normalization, associate enrichments manually in the advanced ROI manaegment menu
  #######################

  Normmethod="lib"
  normalizer=NULL
  normCtrl=NULL
  normCtrlSpike=NULL 


  #nc=1, because in pipeline for genelists we keep it simple

  ##here use totallist to parallelize the code. Each iteration has path to enrichment and a ROI
  tryCatch({

    print ("Associating enrichments to genelist(s) (single core)")
    finallist=lapply(1:length(totallist),function(i) {

      tonorm=names(totallist)[i]
     
      #here put check if transcripts:
      if(getFlag(totallist[[i]])!="transcriptFlag"){
        rangetocov=getRange(totallist[[i]])
        #here use chunksforassociate function to determine the number of chunks
        nchunks=chunksforassociate(rangetocov)
        singlecover=GRbaseCoverage2(Object=rangetocov,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                  signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=nchunks)
        nfact=singlecover[[3]]
        decryptkey=singlecover[[2]]
        singlecover=singlecover[[1]]

        #if all(singlecover==0) -> change nomenclature
        if(verifyzerocov(singlecover)){
          print ("cov is 0s... converting nomelclature to NCBI for coverage...")
          temprange=convertNomenclatureGR(rangetocov,to="NCBI")
          #temproi=setRange(totallist[[i]],temprange)
          singlecover=GRbaseCoverage2(Object=temprange,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                  signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=nchunks)
          nfact=singlecover[[3]]
          decryptkey=singlecover[[2]]
          singlecover=singlecover[[1]]              
        }

        print (paste("coverage",names(totallist)[i]))
        return(list(singlecover,decryptkey,nfact))
      }else{
        rangetransc=getRange(totallist[[i]])
        thirtypercent=round((width(rangetransc)/10)*3)
        #try suppress warnings... BE CAREFUL... maybe it will trim some bases,
        #this could cause problems in the future...
        rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="start"))
        rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="end"))
        #here use chunksforassociate function to determine the number of chunks
        nchunks=chunksforassociate(rangetransc)
        singlecover=GRbaseCoverage2(Object=rangetransc, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
                                signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks=nchunks)
        nfact=singlecover[[3]]
        decryptkey=singlecover[[2]]
        singlecover=singlecover[[1]]
        
        #if all(singlecover==0) -> change nomenclature
        if(verifyzerocov(singlecover)){
          print ("cov is 0s... converting nomelclature to NCBI for coverage...")
          temprange=convertNomenclatureGR(rangetransc,to="NCBI")
          singlecover=GRbaseCoverage2(Object=temprange, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
                                signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks)
          nfact=singlecover[[3]]
          decryptkey=singlecover[[2]]
          singlecover=singlecover[[1]]
        }
        print (paste("coverage",names(totallist)[i]))
        return(list(singlecover,decryptkey,nfact))
      }
    })  




    msg="library-size normalized"

    
    ##now attrib to each ROI the correct enrichments just calculated
    nameROIs=sapply(totallist,getName)
    finallist2=lapply(finallist,"[[",1)
    finaldecryptkeylist=lapply(finallist,"[[",2)
    normfinallist=lapply(finallist,"[[",3)

    splittedresultlist=split(finallist2,nameROIs)
    splittedkeylist=split(finaldecryptkeylist,nameROIs)
    splittednormlist=split(normfinallist,nameROIs)
    rm(finallist2)
    splittedtotallist=split(totallist,nameROIs)
    #"for each ROI selected for enrich. computation..."
    for(i in 1:length(splittedresultlist)){
      
      #previouslist=getBAMlist(splittedtotallist[[i]][[1]])
      ####newenrichimplementation####
      nme=getName(splittedtotallist[[i]][[1]])
      roilocation=match(nme,names(Enrichlist$rawcoverage))
      previousrawvals=Enrichlist$rawcoverage[[roilocation]]
      previousdecryptkey=Enrichlist$decryptkey[[roilocation]]
      previousnormvals=Enrichlist$normfactlist[[roilocation]]
      ###############################

      realnames=names(splittedtotallist[[i]])
      #get from bamlist the name
      pos=match(realnames,BAMvariables$listBAM)
      toname=names(BAMvariables$listBAM)[pos]
      names(splittedresultlist[[i]])=names(splittedkeylist[[i]])=names(splittednormlist[[i]])=toname
      #here, check if some finallist2 was NULL.
      notnulls=!(sapply(splittedresultlist[[i]],is.null))
      #print message of something is NULL
      if(sum(notnulls)<length(toname)){
        sendSweetAlert(
          session = session,
          title = "Some problems",
          text = paste("Some enrichment files gave 'NULL' (",toname[!notnulls],") and were not associated. Please, try again with them",sep=""),
          type="warning"
        )          
      }


      splittedresultlist_notnull=splittedresultlist[[i]][notnulls]
      splittedkeylist_notnull=splittedkeylist[[i]][notnulls]
      splittednormlist_notnull=splittednormlist[[i]][notnulls]
      #finalfinallist=c(previouslist,splittedresultlist_notnull)


      ####newenrichimplementation####
      #to be modified
      finalrawvals=c(previousrawvals,splittedresultlist_notnull)
      finaldecryptkeylist=c(previousdecryptkey,splittedkeylist_notnull)
      finalnormvals=c(previousnormvals,splittednormlist_notnull)
      
      
      rm(previousrawvals)
      rm(previousnormvals)
      rm(previousdecryptkey)
      ###############################
      #which ROI?
      posROI=match(names(splittedresultlist)[i],nomi)

      # ROIvariables$listROI[[posROI]]=setBAMlist(ROIvariables$listROI[[posROI]],finalfinallist)
      ####newenrichimplementation####
      #fill enrichlist
      Enrichlist$rawcoverage[[posROI]]=finalrawvals
      Enrichlist$decryptkey[[posROI]]=finaldecryptkeylist
      Enrichlist$normfactlist[[posROI]]=finalnormvals
      ################################


      rm(finalrawvals)
      rm(finalnormvals)
      rm(finaldecryptkeylist)
      #rm(finalfinallist)
      rm(splittedresultlist_notnull)
      rm(splittedkeylist_notnull)
      rm(splittednormlist_notnull)
      
      for(k in toname){
        logvariables$msg[[length(logvariables$msg)+1]]=paste('Associated ',k,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,'<br>',sep="")
      }
      print(paste('Associated ',toname,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,sep=""))
    }
    rm(splittedresultlist)

    print("forcing gc...")
    print(gc())  
    #alert the user that the block of files has been associated corectly to ROIs
    sendSweetAlert(
      session = session,
      title = "Enrichments associated!",
      text = paste("Enrichment file(s) associated to genelist(s)",sep=""),
      type = "success"
    ) 
    return()

  },warning = function( w ){
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Warning: Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Problems in association",
      #text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
      text=p(style="text-align: left;","Problems in association. This is potentially due to: ",
      tags$li(style="text-align: left;",style="text-align: left;","Enrichment file not found. Check the correct path."),  
      tags$li(style="text-align: left;",style="text-align: left;","The enrichment file is corrupted"),
      tags$li(style="text-align: left;",style="text-align: left;","Memory not sufficient. Options: 1) Try to set the number of cores=1 2) Increase the RAM 3) Save the working session, exit from R, re-open the program and load the working session again (should solve memory fragmentation)")), 
      type = "error"
    )
  },error = function( err ){
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Problems in association",
      #text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
      text=p(style="text-align: left;","Problems in association. This is potentially due to: ",
      tags$li(style="text-align: left;","Enrichment file not found. Check the correct path."),  
      tags$li(style="text-align: left;","The enrichment file is corrupted"),
      tags$li(style="text-align: left;","Memory not sufficient. Options: 1) Try to set the number of cores=1 2) Increase the RAM 3) Save the working session, exit from R, re-open the program and load the working session again (should solve memory fragmentation)")), 
      type = "error"
    )          
  }) 



})








observeEvent(input$Buttonultraeasy,{
  #checks should be useless
  if (!isvalid(input$selectROIultraeasy)| !isvalid(input$selectBAMultraeasy) | length(ROIvariables$listROI)==0){
    sendSweetAlert(
      session = session,
      title = "Input not valid",
      text = "Check the inputs (ROI and one or more enrichment files)",
      type = "error"
    )
    return()
  }


  #execute enrichment association
  nomi=unlist(lapply(ROIvariables$listROI,getName))
  pos=match(input$selectROIultraeasy,nomi)
  ROIs=ROIvariables$listROI[[pos]]
  totallist=list()   
     

  ####newenrichimplementation####
  enrichpresent=names(Enrichlist$rawcoverage[[pos]])
  ###############################     
  allbams=names(BAMvariables$listBAM)
  remaining=setdiff(allbams,enrichpresent)

  #if no BAM associatable with that ROI, simply skip
  if(length(remaining)>0){
    #find the requested enrichment over the remaining enrichments
    remaining_selected=remaining[match(input$selectBAMultraeasy,remaining)]
    #remove those that are NA:
    remaining_selected=remaining_selected[!is.na(remaining_selected)]
    pos2=match(remaining_selected,names(BAMvariables$listBAM))
    paths=BAMvariables$listBAM[pos2]
    partlist=rep(list(ROIs),length(paths))
    names(partlist)=paths
    totallist=c(totallist,partlist)
  }



  #HERE the check for feasibility considering the RAM available. Sum width of all genomic ranges * number of coverage
  ###########################################################################################
  totalBPtobecoveraged=sum(sapply(lapply(lapply(totallist,getRange),width),sum))
  #here, let's suppose an average compression of 2x from the total number of bp (1-byte encryption) (assumption!)
  totalRAMtobeused=totalBPtobecoveraged/2

  avRAM=availRAM()
  if (totalRAMtobeused/1000000>avRAM){
    #not enough RAM for final storage: exit with message
    print(paste("available RAM:",avRAM,"not enough. Predict to occupy about",totalRAMtobeused/1000000))
    sendSweetAlert(
      session = session,
      title = "Not enough memory",
      text = "The RAM available seems not enough for the coverage. Try to reduce the number of ROI or the number of enrichment files for the coverage, or remove enrichment files aready associated.",
      type = "error"
    )          
    return()   
  }
  ###########################################################################################



  Normmethod="lib"
  normalizer=NULL
  normCtrl=NULL
  normCtrlSpike=NULL      



  ##here use totallist to parallelize the code. Each iteration has path to enrichment and a ROI
  tryCatch({


    print ("Computation in single core")
    finallist=lapply(1:length(totallist),function(i) {
      if(Normmethod=="lib"){
        tonorm=names(totallist)[i]
      }else{
        #if no normaliz., custom or spikein, normalizers were defined as normalizer,normCtrl,normCtrlSpikein
        tonorm=normalizer
      }          
      #here put check if transcripts:
      if(getFlag(totallist[[i]])!="transcriptFlag"){
        rangetocov=getRange(totallist[[i]])
        #here use chunksforassociate function to determine the number of chunks
        nchunks=chunksforassociate(rangetocov)
        singlecover=GRbaseCoverage2(Object=rangetocov,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                  signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=nchunks)
        nfact=singlecover[[3]]
        decryptkey=singlecover[[2]]
        singlecover=singlecover[[1]]

        #if all(singlecover==0) -> change nomenclature
        if(verifyzerocov(singlecover)){
          print ("cov is 0s... converting nomelclature to NCBI for coverage...")
          temprange=convertNomenclatureGR(rangetocov,to="NCBI")
          #temproi=setRange(totallist[[i]],temprange)
          singlecover=GRbaseCoverage2(Object=temprange,signalfile=names(totallist)[i],signalfileNorm=tonorm,
                  signalControl=normCtrl,signalControlSpike=normCtrlSpike,multiplFactor=1e+06,nchunks=nchunks)
          nfact=singlecover[[3]]
          decryptkey=singlecover[[2]]
          singlecover=singlecover[[1]]              
        }

        print (paste("coverage",names(totallist)[i]))
        return(list(singlecover,decryptkey,nfact))
      }else{
        rangetransc=getRange(totallist[[i]])
        thirtypercent=round((width(rangetransc)/10)*3)
        #try suppress warnings... BE CAREFUL... maybe it will trim some bases,
        #this could cause problems in the future...
        rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="start"))
        rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="end"))
        #here use chunksforassociate function to determine the number of chunks
        nchunks=chunksforassociate(rangetransc)
        singlecover=GRbaseCoverage2(Object=rangetransc, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
                                signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks=nchunks)
        nfact=singlecover[[3]]
        decryptkey=singlecover[[2]]
        singlecover=singlecover[[1]]
        
        #if all(singlecover==0) -> change nomenclature
        if(verifyzerocov(singlecover)){
          print ("cov is 0s... converting nomelclature to NCBI for coverage...")
          temprange=convertNomenclatureGR(rangetransc,to="NCBI")
          singlecover=GRbaseCoverage2(Object=temprange, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
                                signalControlSpike=normCtrlSpike, multiplFactor=1e+06,nchunks)
          nfact=singlecover[[3]]
          decryptkey=singlecover[[2]]
          singlecover=singlecover[[1]]
        }
        print (paste("coverage",names(totallist)[i]))
        return(list(singlecover,decryptkey,nfact))
      }
    })#,mc.cores=nc)  





    msg="library-size normalized"


      
    ##now attrib to each ROI the correct enrichments just calculated
    nameROIs=sapply(totallist,getName)
    finallist2=lapply(finallist,"[[",1)
    finaldecryptkeylist=lapply(finallist,"[[",2)
    normfinallist=lapply(finallist,"[[",3)

    splittedresultlist=split(finallist2,nameROIs)
    splittedkeylist=split(finaldecryptkeylist,nameROIs)
    splittednormlist=split(normfinallist,nameROIs)
    rm(finallist2)
    splittedtotallist=split(totallist,nameROIs)
    #"for each ROI selected for enrich. computation..."
    for(i in 1:length(splittedresultlist)){
      
      #previouslist=getBAMlist(splittedtotallist[[i]][[1]])
      ####newenrichimplementation####
      nme=getName(splittedtotallist[[i]][[1]])
      roilocation=match(nme,names(Enrichlist$rawcoverage))
      previousrawvals=Enrichlist$rawcoverage[[roilocation]]
      previousdecryptkey=Enrichlist$decryptkey[[roilocation]]
      previousnormvals=Enrichlist$normfactlist[[roilocation]]
      ###############################

      realnames=names(splittedtotallist[[i]])
      #get from bamlist the name
      pos=match(realnames,BAMvariables$listBAM)
      toname=names(BAMvariables$listBAM)[pos]
      names(splittedresultlist[[i]])=names(splittedkeylist[[i]])=names(splittednormlist[[i]])=toname
      #here, check if some finallist2 was NULL.
      notnulls=!(sapply(splittedresultlist[[i]],is.null))
      #print message of something is NULL
      if(sum(notnulls)<length(toname)){
        sendSweetAlert(
          session = session,
          title = "Some problems",
          text = paste("Some enrichment files gave 'NULL' (",toname[!notnulls],") and were not associated. Please, try again with them",sep=""),
          type="warning"
        )          
      }


      splittedresultlist_notnull=splittedresultlist[[i]][notnulls]
      splittedkeylist_notnull=splittedkeylist[[i]][notnulls]
      splittednormlist_notnull=splittednormlist[[i]][notnulls]
      #finalfinallist=c(previouslist,splittedresultlist_notnull)


      ####newenrichimplementation####
      #to be modified
      finalrawvals=c(previousrawvals,splittedresultlist_notnull)
      finaldecryptkeylist=c(previousdecryptkey,splittedkeylist_notnull)
      finalnormvals=c(previousnormvals,splittednormlist_notnull)
      
      
      rm(previousrawvals)
      rm(previousnormvals)
      rm(previousdecryptkey)
      ###############################
      #which ROI?
      posROI=match(names(splittedresultlist)[i],nomi)

      # ROIvariables$listROI[[posROI]]=setBAMlist(ROIvariables$listROI[[posROI]],finalfinallist)
      ####newenrichimplementation####
      #fill enrichlist
      Enrichlist$rawcoverage[[posROI]]=finalrawvals
      Enrichlist$decryptkey[[posROI]]=finaldecryptkeylist
      Enrichlist$normfactlist[[posROI]]=finalnormvals
      ################################


      rm(finalrawvals)
      rm(finalnormvals)
      rm(finaldecryptkeylist)
      #rm(finalfinallist)
      rm(splittedresultlist_notnull)
      rm(splittedkeylist_notnull)
      rm(splittednormlist_notnull)
      
      for(k in toname){
        logvariables$msg[[length(logvariables$msg)+1]]=paste('Associated ',k,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,'<br>',sep="")
      }
      print(paste('Associated ',toname,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,sep=""))
    }



    rm(splittedresultlist)

    print("forcing gc...")
    print(gc())  
    #alert the user that the block of files has been associated corectly to ROIs
    sendSweetAlert(
      session = session,
      title = "Enrichments associated!",
      text = paste("Selected enrichment files has been associated to selected ROI",sep=""),
      type = "success"
    ) 
    return()
  },warning = function( w ){
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Warning: Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Problems in association",
      #text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
      text=p(style="text-align: left;","Problems in association. This is potentially due to: ",
      tags$li(style="text-align: left;",style="text-align: left;","Enrichment file not found. Check the correct path."),  
      tags$li(style="text-align: left;",style="text-align: left;","The enrichment file is corrupted"),
      tags$li(style="text-align: left;",style="text-align: left;","Memory not sufficient. Options: 1) Try to set the number of cores=1 2) Increase the RAM 3) Save the working session, exit from R, re-open the program and load the working session again (should solve memory fragmentation)")), 
      type = "error"
    )
  },error = function( err ){
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Problems in association",
      #text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
      text=p(style="text-align: left;","Problems in association. This is potentially due to: ",
      tags$li(style="text-align: left;","Enrichment file not found. Check the correct path."),  
      tags$li(style="text-align: left;","The enrichment file is corrupted"),
      tags$li(style="text-align: left;","Memory not sufficient. Options: 1) Try to set the number of cores=1 2) Increase the RAM 3) Save the working session, exit from R, re-open the program and load the working session again (should solve memory fragmentation)")), 
      type = "error"
    )          
  }) 


})










