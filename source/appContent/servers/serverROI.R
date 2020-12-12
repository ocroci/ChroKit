

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
          selectedbam=getBAMlist(ROIvariables$listROI[[pos]])
        }else{
          selectedbam=list()
        }

        if (nchar(input$ROIname)>=1){

          roigeneration=suppressWarnings(generateROI(selectedranges,selectedfix,overlapregions,notoverlapregions,
                        input$choiceROI,input$choiceoverlapROI, input$choicenotoverlapROI,bamlist=selectedbam,minbp=input$minOverlapNEWROI,strandSpecific=input$StrandSpecOverlapNEWROI))
         
          ROI=roigeneration[[1]]

          if(length(ROI)>0){
            BAMlist=roigeneration[[2]]
            newfix=roigeneration[[3]]

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


            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                        name=input$ROIname,
                                        range=ROI,
                                        fixed=newfix,
                                        BAMlist=BAMlist,
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
              #check if new size is greater. If so, destroy BAM files matrix list from ROI object
              # (be smart and re-calculate with bam files the missing parts)
              #if both down and up are smaller, resize range and BAMlist (if not null)

              #nomi=lapply(ROIvariables$listROI,getName)
              pos=match(input$selectROItoresize,nomi)
              roi=getRange(ROIvariables$listROI[[pos]])
              oldSource=getSource(ROIvariables$listROI[[pos]])
              fix=start(getFixed(ROIvariables$listROI[[pos]]))
              bamlist=getBAMlist(ROIvariables$listROI[[pos]])

              #split positive or undetermined strand from negative strand
              pos_positive=as.logical(!strand(roi)=="-")
              pos_negative=as.logical(strand(roi)=="-")
              #remove ranges that would have negative starts with new width in positive strand
              idx=(fix-input$sliderUpstreamROI)>0
              idx_positive_and_good=pos_positive & idx
              #negative strand
              idx=(fix-input$sliderDownstreamROI)>0
              idx_negative_and_good=pos_negative & idx
              idx_positive_or_negative_good=idx_positive_and_good|idx_negative_and_good
              roi=roi[idx_positive_or_negative_good]
              fix=fix[idx_positive_or_negative_good]
              bamlist=lapply(bamlist,function(i) {i[idx_positive_or_negative_good]})

              pos_positive=as.logical(!strand(roi)=="-")
              pos_negative=as.logical(strand(roi)=="-")
              roi_positive=roi[pos_positive]
              roi_negative=roi[pos_negative]
              fix_positive=fix[pos_positive]
              fix_negative=fix[pos_negative]
              bamlist_positive=lapply(bamlist,function(i) {i[pos_positive]})
              bamlist_negative=lapply(bamlist,function(i) {i[pos_negative]})

              #resize positivee strand:
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
              newflag=getFlag(ROIvariables$listROI[[pos]])
              #ROIvariables$listROI[[pos]]=setRange(ROIvariables$listROI[[pos]],roi)
              #ROIvariables$listROI[[pos]]=setFix(ROIvariables$listROI[[pos]],newfix)
              newrange=roi

              #to be checked when BAM files will be available
              if (length(bamlist)>0){
                #if not everything is less than the old size in any position, reset BAMlist!
                if (! (all((fix_positive-oldstarts_positive) >= newstart) & all((oldends_positive-fix_positive)>= newend) )
                          | !(all((fix_negative-oldstarts_negative) >= newend) & all((oldends_negative-fix_negative)>= newstart) )    ){
                  #reset BAMlist or improve only the delta margins of already existing matrixes
                  #ROIvariables$listROI[[pos]]=resetBAMlist(ROIvariables$listROI[[pos]])
                  logvariables$msg[[length(logvariables$msg)+1]]= paste('Removed ',input$selectROItoresize,' associated BAMs, because new size > old size...<br>',sep="")
                  newBAMlist=list()
                }else{
                  toadd=paste(toadd,"(kept associated BAMs)")
                  #split - and + for BAM resize 
                 
                  shiftleft_positive= (fix_positive-oldstarts_positive)- newstart
                  shiftright_positive=(oldends_positive-fix_positive)-newend


                  if (any(pos_positive)){
                    bamlistnew_positive=lapply(1:length(bamlist_positive),function(i) {
                      #cut also BAM file if some ranges have a resize that will lead to a negative start
                      x=bamlist_positive[[i]]
                      slicedbam=lapply(1:length(x), function(k) {
                                            return(x[[k]] [(shiftleft_positive[k]+1): (length(x[[k]])-shiftright_positive[k]) ])
                                            })
                      return(slicedbam)
                    } )

                  }

                  shiftleft_negative= (fix_negative-oldstarts_negative)- newend
                  shiftright_negative=(oldends_negative-fix_negative)-newstart


                  if (any(pos_negative)){
                    bamlistnew_negative=lapply(1:length(bamlist_negative),function(i) {
                      #cut also BAM file if some ranges have a resize that will lead to a negative start
                      x=bamlist_negative[[i]]
                      slicedbam=lapply(1:length(x), function(k) {
                                            return(x[[k]] [(shiftleft_negative[k]+1): (length(x[[k]])-shiftright_negative[k]) ])
                                            })
                      return(slicedbam)
                    } )                    
                  }



                  #join bamlistnew_positive and bamlistnew_negative back to bamlist
                  bamlistnew=lapply(1:length(bamlist),function(i) {
                    x=bamlist[[i]]
                    if (any(pos_positive)){
                      x[pos_positive]=bamlistnew_positive[[i]]
                    }
                    
                    if (any(pos_negative)){
                      x[pos_negative]=bamlistnew_negative[[i]]
                    }
                    return(x)
                  })

                  names(bamlistnew)=names(bamlist)
                  newBAMlist=bamlistnew
                  #ROIvariables$listROI[[pos]]=setBAMlist(ROIvariables$listROI[[pos]],bamlistnew)

                  #x is single list of baseCoverage (BAM)
                  logvariables$msg[[length(logvariables$msg)+1]]= paste('Conserved ',input$selectROItoresize,' associated BAMs, because new size < old size...<br>',sep="")

                }
              }else{
                newBAMlist=list()
              }

              newSource=c(oldSource,list(toadd))
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                            name=input$ROInameResize,
                                            range=newrange,
                                            fixed=newfix,
                                            BAMlist=newBAMlist,
                                            flag=newflag,
                                            source=newSource) 


              if(length(unique(fix_positive-oldstarts_positive))==1 & length(unique(oldends_positive-fix_positive))==1){
                logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI [-',(fix_positive-oldstarts_positive)[1],'; +',(oldends_positive-fix_positive)[1],']<br>',sep="")
                print(paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI [-',(fix_positive-oldstarts_positive)[1],'; +',(oldends_positive-fix_positive)[1],']',sep=""))
              }else{
                logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI with original widths<br>',sep="")
                print(paste('Created ',input$ROInameResize,' ROI [-',input$sliderUpstreamROI,'; +',input$sliderDownstreamROI,'] from ',input$selectROItoresize,' ROI with original widths',sep=""))
              }
            }else{
              #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Are you trying to resize transcripts around midpoint? You cannot.<br></font>',sep="")
              sendSweetAlert(
                session = session,
                title = "Bad choice",
                text = "You cannot resize transcripts using their midpoint",
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
          oldSource=getSource(roi)
          getbam=names(getBAMlist(roi))
          #if there is at least one BAM
          if (!is.null(getbam)){
            bam=getBAMlist(roi)
            pos2=match(input$selectBAMtoCenterSummit,getbam)
            bamselected=bam[[pos2]]
            toadd=paste("Centered on summit using",input$selectBAMtoCenterSummit,"enrichment (range width=1)")
            #use the function to select summit:
            newrange=summitFromBaseCoverage(Object=getRange(roi),baseCoverageOutput=bamselected)

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
          getbam=names(getBAMlist(roi))
          if (!is.null(getbam)){
            bam=getBAMlist(roi)
            pos2=match(input$selectBAMtoFilter,getbam)
            bamselected=bam[[pos2]]  
            sums=unlist(lapply(bamselected,sum))
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
              for(i in 1:length(bam)){
                newBAMlist[[i]]=bam[[i]][selectedPositions]
              }
              names(newBAMlist)=getbam
              newfix=getFixed(roi)[selectedPositions]
              newflag=getFlag(roi)
              perc=round(tabForMessage["TRUE"]/length(selectedPositions)*100,2)
              toadd=paste("Filtered on",input$selectBAMtoFilter,"BAM file enrichment (kept ",tabForMessage["TRUE"],"/",length(selectedPositions),"ranges,",perc,"%), excluding",excludedlow,"% low and",excludedhigh,"% high")
              newSource=c(oldSource,list(toadd))
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                              name=input$ROInameFilter,
                                              range=newrange,
                                              fixed=newfix,
                                              BAMlist=newBAMlist,
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
              bams=getBAMlist(roi)
              if(length(bams)>0){
                newBAMlist=list()
                for(i in 1:length(bams)){
                  newBAMlist[[i]]=bams[[i]][selectedPositions]
                }
                names(newBAMlist)=names(bams)            
              }else{
                newBAMlist=list()
              }
              
              newfix=getFixed(roi)[selectedPositions]
              newflag=getFlag(roi)
              perc=round(tabForMessage["TRUE"]/length(selectedPositions)*100,2)
              toadd=paste("Filtered on width (kept ",tabForMessage["TRUE"],"/",length(selectedPositions),"ranges,",perc,"%), excluding",excludedlow,"% low and",excludedhigh,"% high")
              newSource=c(oldSource,list(toadd))
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                              name=input$ROInameFilterWIDTH,
                                              range=newrange,
                                              fixed=newfix,
                                              BAMlist=newBAMlist,
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
            bams=getBAMlist(roi)
            if(length(bams)>0){
              newBAMlist=list()
              for(i in 1:length(bams)){
                newBAMlist[[i]]=bams[[i]][selectedPositions]
              }
              names(newBAMlist)=names(bams)            
            }else{
              newBAMlist=list()
            }
            
            newfix=getFixed(roi)[selectedPositions]
            newflag=getFlag(roi)
            perc=round(input$numberSample/length(rangefromROI)*100,2)

            toadd=paste("ROI sampled (kept ",input$numberSample,"/",length(rangefromROI),"ranges,",perc,"%)")
            newSource=c(oldSource,list(toadd))
            ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                            name=input$ROInameSample,
                                            range=newrange,
                                            fixed=newfix,
                                            BAMlist=newBAMlist,
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

  tryCatch({
	  if(R35){
	    print(paste("Downloading",BSstring,"..."))
	    BiocManager::install(BSstring, version = bioCversion,ask=FALSE)
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


  },
  warning = function( w ){
	sendSweetAlert(
        session = session,
        title = "Problems in downloading BSgenome",
        text = "Check your internet connection, or change bioCversion (for example, for R 3.6, the bioCversion 3.9 is needed)",
        type = "error"
    )
  },
  error = function( err ){
      sendSweetAlert(
        session = session,
        title = "Problems in downloading BSgenome",
        text = "Check your internet connection, or change bioCversion (for example, for R 3.6, the bioCversion 3.9 is needed)",
        type = "error"
      )
  })


  

})



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
            allowedpattern=TRUE
            if(input$choiceWherePattern=="fromGenome"){
              checkComplexity=pattern_splitted!="N" &pattern_splitted!="." &pattern_splitted!="-"
              if( sum(checkComplexity)<4){
                allowedpattern=FALSE
              }
            }
            if(allowedpattern){
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

                    if(input$choiceWherePattern=="fromGenome"){
                      both=TRUE
                    }

                    toadd=paste("extracting of ",pattern," pattern",sep="")
                    if(both){
                      toadd=paste(toadd," searched in both strands",sep="")
                    }else{
                      toadd=paste(toadd," searched in strand-specific way",sep="")
                    }

                    if(input$choiceWherePattern=="fromROI"){
                      pos=match(input$selectROItoExtractPattern,nomi)
                      roi=ROIvariables$listROI[[pos]]
                      roi_range=getRange(roi)
                      oldSource=getSource(roi)
                      newSource=c(oldSource,list(toadd))
                      #apply the function extractPattern:
                      extracted=suppressWarnings(extractPattern(Subject=roi_range,BSgenomeDB=eval(parse(text=BSstring)),pattern=pattern,bothstrands=both))
                      
                    }else{
                      newSource=list(toadd)
                      extracted=suppressWarnings(extractPattern(Subject=NULL,BSgenomeDB=eval(parse(text=BSstring)),pattern=pattern,bothstrands=both))
                    }



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
                      ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                                      name=input$ROInamePattern,
                                                      range=extracted,
                                                      fixed=newFix,
                                                      BAMlist=list(),
                                                      flag="Pattern",
                                                      source=newSource)
                      if(input$choiceWherePattern=="fromGenome"){
                        print(paste(input$ROInamePattern,"ROI created, ",toadd," from entire genome",sep=""))
                      }else{
                        print(paste(input$ROInamePattern,"ROI created, ",toadd," from ",input$selectROItoExtractPattern," ROI",sep=""))
                      }
                      
                      
                    }else{
                      sendSweetAlert(
                        session = session,
                        title = "Empty ROI produced",
                        text = paste("pattern '",pattern,"' is not present in any range of ROI ",input$selectROItoExtractPattern,sep=""),
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
                text = "The sequence pattern must be >= 4 not-N-letters when searching the entire genome. Try to add more letters (that are not 'N')",
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
        text = "You must select a genome assembly. Go to 'Database' section to import a genome assembly",
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
#select ROI to visualize
observeEvent(input$msg_viewRois_selectRoi, {
  boxHelpServer(msg_viewRois_selectRoi)
})
#visualization
observeEvent(input$msg_viewRois_visualization, {
  boxHelpServer(msg_viewRois_visualization)
})
#information
observeEvent(input$msg_viewRois_information, {
  boxHelpServer(msg_viewRois_information)
})

#view ROI
# view stat (width, number) ranges
observeEvent(input$updatechoiceROI,{
  set.seed(123)
  #response to checkbutton. Now check if there is at least one ROI selected and at least one ROI exists in 
  # the list of ROI

  if (length(input$confirmviewROI)>0 & length(ROIvariables$listROI)){
    #now, if length(input$confirmviewROI)==1 (one and only one selected, we have a selected ROI => sliderinput for quantile width)
    nomi=unlist(lapply(ROIvariables$listROI,getName))

    selection=list()
    for (i in 1:length(input$confirmviewROI)){
        pos2=match(input$confirmviewROI[i],nomi)
        selection[[i]]=getRange(ROIvariables$listROI[[pos2]])
    }
    lun=unlist(lapply(selection,length))

    #continue only if all ROIs selected have length >0
    if (all(lun>0)){
      if (length(input$confirmviewROI)==1){
        pos=match(input$confirmviewROI,nomi)   

        ROIvariables$selected=getRange(ROIvariables$listROI[[pos]])
        lungh=length(ROIvariables$selected)

        ROIvariables$selectedname=paste(nomi[pos]," (",lungh,")",sep="")
        
      }else{
        ROIvariables$selected=NULL
        ROIvariables$selectedname=NULL     
      }
      #here we have >=1 selected ROI => plot in any case density (with different colors) and barplot for peak number
      #have to set the densities and colors. then select the max(x) and max (y) for the blank plot
      ROIvariables$listfordensity=list()
      ROIvariables$listselected=list()
      ROIvariables$listselectednames=list()
      for (i in 1:length(input$confirmviewROI)){
        pos=match(input$confirmviewROI[i],nomi)
        ROIvariables$listselected[[i]]=getRange(ROIvariables$listROI[[pos]])
        lungh=length(ROIvariables$listselected[[i]])

        ROIvariables$listselectednames[[i]]=paste(nomi[pos]," (",lungh,")",sep="")
        x=ROIvariables$listselected[[i]]
        ROIvariables$listfordensity[[i]]=density(log2(width( x  )))

      }

      n=length(input$confirmviewROI)
      #ROIvariables$colorsfordensity <- distinctColorPalette(n)
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      ROIvariables$colorsfordensity <-sample(col_vector, n)

      ROIvariables$changed=ROIvariables$changed+1

      #download button for PDF peak number 
      output$saveviewpeaksnumberROI=renderUI({downloadButton('saveviewpeaksnumberROIbutton', 'Get PDF')})
      output$saveviewpeakswidthROI=renderUI({downloadButton('saveviewpeakswidthROIbutton', 'Get PDF')})

    }else{
      ROIvariables$colorsfordensity=NULL
      ROIvariables$listfordensity=NULL
      ROIvariables$listselected=NULL
      ROIvariables$listselectednames=NULL
      ROIvariables$selected=NULL
      ROIvariables$selectedname=NULL 
      output$saveviewpeaksnumberROI<-renderUI({NULL})
      output$saveviewpeakswidthROI<-renderUI({NULL})
      ROIvariables$changed=ROIvariables$changed+1
      logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red"> All GRanges selected must have length >0<br></font>',sep="")
    }


  }else{
    ROIvariables$colorsfordensity=NULL
    ROIvariables$listfordensity=NULL
    ROIvariables$listselected=NULL
    ROIvariables$listselectednames=NULL
    ROIvariables$selected=NULL
    ROIvariables$selectedname=NULL
    output$saveviewpeaksnumberROI<-renderUI({NULL}) 
    output$saveviewpeakswidthROI<-renderUI({NULL})
    ROIvariables$changed=ROIvariables$changed+1
    sendSweetAlert(
      session = session,
      title = "Missing ROI",
      text = "Select at least one ROI",
      type = "error"
    )
  }

},ignoreInit=TRUE)




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
    plot(1, type="n", xlab="log2 width", ylab="density", xlim=c(xmin, xmax), ylim=c(ymin, ymax),main="frequency plot width")
    for (i in 1:length(ROIvariables$listfordensity)){
      lines(ROIvariables$listfordensity[[i]],col=ROIvariables$colorsfordensity[[i]],lwd=2)
    }
    legend("topright",leg,col=unlist(ROIvariables$colorsfordensity),lty=1,lwd=2,bty = "n")

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
observeEvent(input$msg_getRois_roiSelection, {
  boxHelpServer(msg_getRois_roiSelection)
})
#table explore 
observeEvent(input$msg_getRois_preview, {
  boxHelpServer(msg_getRois_preview)
})


observeEvent(input$showdataframeROI,{
  if(length(ROIvariables$listROI)>0 & length(input$listgetROI)>0){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$listgetROI,nomi)
    ROI=ROIvariables$listROI[[pos]]
    range=getRange(ROI)

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
    
    if(!is.null(input$ROIoptionsToViewENRICHMENTS)){
      #append selected enrichments to dataframe
      bams=getBAMlist(ROI)
      bamnames=names(getBAMlist(ROI))
      pos=match(input$ROIoptionsToViewENRICHMENTS,bamnames)
      bams_selected=bams[pos]
      if(nc==1){
        megalist=lapply(1:length(bams_selected),function(i) {
          summ=sapply(bams_selected[[i]],sum)
          return(summ)          
        })
      }else{
        decision=0
        tryCatch({
          megalist=mclapply(1:length(bams_selected),function(i) {
            summ=sapply(bams_selected[[i]],sum)
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
          megalist=lapply(1:length(bams_selected),function(i) {
            summ=sapply(bams_selected[[i]],sum)
            return(summ)          
          })            
        }
      }

      #append to dataframe the sum of the enrichments
      for (i in 1:length(megalist)){
        name=names(bams_selected)[i]
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
          paste("require promoters (use a database for this)")
        })  
        output$previewROItodownloadbutton<-renderUI({
          NULL
        }) 
        tosave$genelistROIwindow<-NULL       
        sendSweetAlert(
          session = session,
          title = "Annotated elements not found",
          text = "Promoters of a specific genome assembly not found, but you need them: go to 'Databases' section and choose a genome assembly",
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
      bams=getBAMlist(ROIvariables$listROI[[pos[i]]])
      bamnames=names(bams)  
      pos2=match(input$selectBAMtoDeassociate,bamnames)
      #remove NA (not all selected BAMs are present in all ROIs)
      pos2=pos2[!is.na(pos2)]
      bams[pos2]<-NULL
      #if NULL, recreate empty BAMlist
      if (is.null(bams)){
        bams=list()
      }   
      ROIvariables$listROI[[pos[i]]]=setBAMlist(ROIvariables$listROI[[pos[i]]],bams)   
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
      bampresent=names(getBAMlist(ROIs[[i]]))
      allbams=names(BAMvariables$listBAM)
      remaining=setdiff(allbams,bampresent)
      #find the requested enrichment over the remaining enrichments
      remaining_selected=remaining[match(input$selectBAMtoassociate,remaining)]
      #if no BAM associatable with that ROI, simply skip
      if(!all(is.na(remaining_selected))){
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



    if(input$selectMethodForNorm=="librarysize"){
      Normmethod="lib"
      normalizer=NULL
      normCtrl=NULL
      normCtrlSpike=NULL      
    }else if(input$selectMethodForNorm=="customnorm"){
      if(length(input$customNormalizer)>0){
        Normmethod="custom"
        normalizer=BAMvariables$listBAM[match(input$customNormalizer,names(BAMvariables$listBAM))]
        normCtrl=NULL
        normCtrlSpike=NULL
      }
    }else if(input$selectMethodForNorm=="spikein") {
      if(length(input$customNormalizer)>0 & length(input$ctrlNormalizer)>0 & input$ctrlSpikeinNormalizer>0){
        Normmethod="spikein"
        normalizer=BAMvariables$listBAM[match(input$customNormalizer,names(BAMvariables$listBAM))]
        normCtrl=BAMvariables$listBAM[match(input$ctrlNormalizer,names(BAMvariables$listBAM))]
        normCtrlSpike=BAMvariables$listBAM[match(input$ctrlSpikeinNormalizer,names(BAMvariables$listBAM))]
      }
    }else if(input$selectMethodForNorm=="nonorm"){
      Normmethod="no"
      normalizer=NULL
      normCtrl=NULL
      normCtrlSpike=NULL
    }else{
      print("something wrong, no normalization existing...")
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
            singlecover=cover(Object=totallist[[i]],signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,signalControlSpike=normCtrlSpike)
            print (paste("coverage",names(totallist)[i]))
            return(singlecover)
          }else{
            rangetransc=getRange(totallist[[i]])
            thirtypercent=round((width(rangetransc)/10)*3)
            #try suppress warnings... BE CAREFUL... maybe it will trim some bases,
            #this could cause problems in the future...
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="start"))
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="end"))

            singlecover=GRbaseCoverage2(Object=rangetransc, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,
            												signalControlSpike=normCtrlSpike, multiplFactor=1e+06)
            print (paste("coverage",names(totallist)[i]))
            return(singlecover)
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
            singlecover=cover(Object=totallist[[i]],signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,signalControlSpike=normCtrlSpike)
            print (paste("coverage",names(totallist)[i]))
            return(singlecover)
          }else{
            rangetransc=getRange(totallist[[i]])
            thirtypercent=round((width(rangetransc)/10)*3)
            #try suppress warnings... BE CAREFUL... maybe it will trim some bases,
            #this could cause problems in the future...
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="start"))
            rangetransc=suppressWarnings(resize(rangetransc,width=thirtypercent+width(rangetransc),fix="end"))

            singlecover=GRbaseCoverage2(Object=rangetransc, signalfile=names(totallist)[i],signalfileNorm=tonorm,signalControl=normCtrl,signalControlSpike=normCtrlSpike, multiplFactor=1e+06)
            print (paste("coverage",names(totallist)[i]))
            return(singlecover)
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
      splittedresultlist=split(finallist,nameROIs)
      splittedtotallist=split(totallist,nameROIs)
      for(i in 1:length(splittedresultlist)){
        previouslist=getBAMlist(splittedtotallist[[i]][[1]])
        realnames=names(splittedtotallist[[i]])
        #get from bamlist the name
        pos=match(realnames,BAMvariables$listBAM)
        toname=names(BAMvariables$listBAM)[pos]
        names(splittedresultlist[[i]])=toname
        #here, check if some finallist was NULL.
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
        finalfinallist=c(previouslist,splittedresultlist_notnull)
        #which ROI?
        posROI=match(names(splittedresultlist)[i],nomi)
        ROIvariables$listROI[[posROI]]=setBAMlist(ROIvariables$listROI[[posROI]],finalfinallist)
        for(k in toname){
          logvariables$msg[[length(logvariables$msg)+1]]=paste('Associated ',k,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,'<br>',sep="")
        }
        print(paste('Associated ',toname,' enrichment to ',names(splittedresultlist)[i],' ROI, ',msg,sep=""))
      }


      print("forcing gc...")
      print(gc())  
      #alert the user that the block of files has been associated corectly to ROIs
      sendSweetAlert(
        session = session,
        title = "Enrichments associated!",
        text = paste("Selected enrichment files has been associated to selected ROIs",sep=""),
        type = "success"
      ) 


    },warning = function( w ){
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Warning: Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Problems in association",
        text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
        type = "error"
      )
    },error = function( err ){
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Problems in associating BAM or WIG. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Problems in association",
        text = "BAM or WIG files not associated. Corrupted file, or, if number cores>1, SEGFAULT. Try using less cores or even 1. Maybe enrichment files not found: be sure that correct file path are provided, and files are on the same machine of the running program",
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
      getbam=names(getBAMlist(roi))  
      #check if name does not match with another existing BAM and at this ROI have at least one BAM
      if (!input$newBAMname %in% getbam & length(getbam)>0 ){

        for (i in 1:length(ROIvariables$listROI)){
          getbam=names(getBAMlist(ROIvariables$listROI[[i]]))
          if(!input$newBAMname %in% getbam & length(getbam)>0){
            pos2=match(input$selectedBAMtoRename,getbam)
            #if BAM found in the current ROI in the loop
            if(!is.na(pos2)){
              bamselected=getBAMlist(ROIvariables$listROI[[i]])
              oldname=names(bamselected)[pos2]
              names(bamselected)[pos2]=input$newBAMname
              ROIvariables$listROI[[i]]=setBAMlist(ROIvariables$listROI[[i]],bamselected)
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
    getbam=names(getBAMlist(roi))  
    bamsl=getBAMlist(roi)
    
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
        bamsl=bamsl[order(ord)]

        ROIvariables$listROI[[pos]]=setBAMlist(ROIvariables$listROI[[pos]],bamsl)
        if (!identical(unique(listprovv),unique(allnumbers))){
          logvariables$msg[[length(logvariables$msg)+1]]= paste('BAM files reordered...<br>',sep="")
          print("BAM files reordered")
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
# GO gene ontology analysis
##############################################################################
##############################################################################
##############################################################################
##############################################################################


###react to help buttons:
#parameters
observeEvent(input$msg_goAnalysis_parameters, {
  boxHelpServer(msg_goAnalysis_parameters)
})
#GO plot
observeEvent(input$msg_goAnalysis_goPlot, {
  boxHelpServer(msg_goAnalysis_goPlot)
})

#GO table
observeEvent(input$msg_goAnalysis_goTable, {
  boxHelpServer(msg_goAnalysis_goTable)
})





toListenGO <- reactive({
    list(input$doTheGO)#,input$doTheGO2)
})

observeEvent(toListenGO(),{
  if(!isvalid(input$minSizeGO)|!isvalid(input$maxSizeGO)){
    #numeric input not valid
    sendSweetAlert(
      session = session,
      title = "Geneset size thresholds not valid",
      text = "Min or Max geneset size thresholds are not valid",
      type = "error"
    ) 
    return() 
  }
  #check if any geneset was selected!!
  if(isvalid(input$selectedGenesetsGO)){
    #if selected ROI:
    # if valid input$chooseCriteriaROIGO, ROI ok, otherwise ROI not seelcted or no DB or no ROI
    #if selected fromcutomGenelist
    # testarea must not be empty. If no DB, must be symbols
    if(input$chooseSourceGO=="fromROI"){

      if(!isvalid(input$selectROIGO)){
        sendSweetAlert(
          session = session,
          title = "No ROI selected",
          text = "You have to select at least one ROI",
          type = "error"
        )          
        return()
      }

      if(isvalid(input$chooseCriteriaROIGO)){
        #create list for gene lists in input (n= number of ROI selected)
        inputGeneList=as.list(rep(NA,length(input$selectROIGO)))
        names(inputGeneList)=input$selectROIGO
        nomi=unlist(lapply(ROIvariables$listROI,getName))
        pos=match(input$selectROIGO,nomi)
        ROI=ROIvariables$listROI[pos]


        if(input$chooseCriteriaROIGO=="windowGene"){
          if(isvalid(input$WindowROIGO)){
            print( paste("GO analysis in genomic window from", paste(input$selectROIGO,collapse="; ")) )

            pospromo=match("promoters",nomi)
            ROIpromo=ROIvariables$listROI[[pospromo]]
            fixedpromo=getFixed(ROIpromo)

            #open the window around the Fix of promoters (TSSs), for example for 10kb window:
            #                         10k  10k      
            #promoters:             |----|----|       |----|----|
            promo=suppressWarnings(unique(resize(fixedpromo,width=input$WindowROIGO*2,fix="center")))
            #if (!is.null(ROI)) {}
            for (i in 1:length(ROI)){
              #center the range in the Fixed point
              range=getFixed(ROI[[i]])

              ov=suppressWarnings(countOverlaps(promo,range))
              #WARNING: more than one midpoint of grange can be associated with the same
              # promoter!!
              #extract gene id of overlapping promoters with that window
              promo_selected=promo[ov>0]

              emd=as.data.frame(elementMetadata(promo_selected))
              symbols=unique(as.character(emd$symbol))
              inputGeneList[[i]]=symbols
            }

          }else{
            #window not valid
            sendSweetAlert(
              session = session,
              title = "Genomic window not valid",
              text = "The genomic window is empty or is not a number",
              type = "error"
            )
            return()
          }
          
        }else{
          #nearest gene, assign symbols in inputGeneList elements
          print( paste("GO analysis of nearest genes annotated to", paste(input$selectROIGO,collapse="; ")) )
          for (i in 1:length(ROI)){
            range=getRange(ROI[[i]])
            emd=as.data.frame(elementMetadata(range))
            symbols=unique(as.character(emd$symbol))
            inputGeneList[[i]]=symbols
          }

        }
      }else{
        #no ROI, no ROI selected or no database (ROI not annotated)
        sendSweetAlert(
          session = session,
          title = "No database or no ROI",
          text = "Be sure to select at least one ROI and that a genome assembly is active. Go to 'Databases' section for that",
          type = "error"
        )   
        return()
        
      }

    }else{
      #here, chosen custom gene list instead of ROI
      if(isvalid(input$pastedGenesGO)){
        print("GO analysis from a pasted set of genes")
        genelist=strsplit(input$pastedGenesGO,split="\n")[[1]]
        uniquegenelist=unique(genelist)
        lostinunique=length(genelist)-length(uniquegenelist)
        print (paste(lostinunique," genes lost because duplicated",sep=""))
        if(isvalid(input$chooseIDgeneGO)){

          nomi=unlist(lapply(ROIvariables$listROI,getName))
          pospromo=match("promoters",nomi)
          ROIpromo=ROIvariables$listROI[[pospromo]]
          rangepromo=getFixed(ROIpromo)
          emd=as.data.frame(elementMetadata(rangepromo))

          #here, take promoters and translate ***-> to symbols
          if(input$chooseIDgeneGO=="entrez"){
            #gene_id column
            pos=match(uniquegenelist,emd$gene_id)
          }else if(input$chooseIDgeneGO=="ensembl"){
            #ensembl_id column
            pos=match(uniquegenelist,emd$ensembl_id)
          }else if(input$chooseIDgeneGO=="refseq"){
            #refSeq_id column
            pos=match(uniquegenelist,emd$refSeq_id)
          }else{
            #match with itself, but same capital letters
            pos=match(toupper(uniquegenelist),toupper(emd$symbol))
          }
          translatedSymbols=as.character(emd[pos,]$symbol)

        }else{
          #only symbols shuld be provided
          translatedSymbols=uniquegenelist
        }


        translatedSymbols=translatedSymbols[!is.na(translatedSymbols)]
        if(length(translatedSymbols)<1){
          #no match found... problems
          sendSweetAlert(
            session = session,
            title = "Invalid genes",
            text = "Are you sure to have typed the correct gene IDs or symbols in the text area? Maybe you put spaces or the wrong kind of IDs",
            type = "error"
          )
          return()
        }



        inputGeneList=as.list(rep(NA,1))
        names(inputGeneList)="custom_gene_list"
        inputGeneList[[1]]=translatedSymbols

      }else{
        sendSweetAlert(
          session = session,
          title = "Empty/not valid genes",
          text = "I didn't find genes in the text area or those genes are not valid IDs/names",
          type = "error"
        )         
        return()
      }
      
    }

    #from here, we have the list of gene symbols to query
    #?parallel GO on elements of inputGeneList
    #read geneSets files and store them in logvariables$temporary_GMTstorage list
    #the elements of this list correspond to gene sets and are NA at the beginning
    posSets=match(input$selectedGenesetsGO,names(GenesetsGMT))
    SetsToRead=GenesetsGMT[posSets]
    allterms=as.list(rep(NA,length(input$selectedGenesetsGO)))
    names(allterms)=names(SetsToRead)
    
    #do it in parallel even if reading files?
    for(i in 1:length(SetsToRead)){
      if(!is.na(logvariables$temporary_GMTstorage[names(SetsToRead)[i]])){
        #take cached values from temporary file
        allterms[[names(SetsToRead)[i]]]<-logvariables$temporary_GMTstorage[[names(SetsToRead)[i]]]
      }else{
        #read the geneset GMT file ex-novo
        GMT=readGMT(SetsToRead[[i]])
        logvariables$temporary_GMTstorage[[names(SetsToRead)[i]]]=GMT
        allterms[[names(SetsToRead)[i]]]<-GMT
      }
    }

    #union of all genesets of all categories
    terms=do.call("rbind",allterms)
    dt=data.table(terms)
    dt=unique(dt)
    termsdt=as.data.frame(unique(dt))
    
    ##################################################
    #calculate how many genes are correctly found
    allgenes=toupper(Reduce("union",inputGeneList))
    identif=match(allgenes,as.character(termsdt$gene))
    ##################################################



    #now use the GO function for each roi or for genelist. Try in parallel
    minsizeGO=input$minSizeGO
    maxsizeGO=input$maxSizeGO
    #execeute the GO in parallel. If nc==1 (maybe windows OS) use lapply
    #be careful! if window gene, block of current range could be 0!!
    if(nc==1){
      reslist=lapply(1:length(inputGeneList),function(i){
        currentblock=inputGeneList[[i]]
        #transform uppercase 
        currentblock=toupper(currentblock)
        if (length(currentblock)>1){
          #do the GO with hypergeometric test
          res=GOcalc(gene=currentblock,terms=termsdt,minsize=minsizeGO,maxsize=maxsizeGO,padj_method="BH")
        }else{
          res=NULL
        }
        return(res)
      })
    }else{
      reslist=mclapply(1:length(inputGeneList),function(i){
        currentblock=inputGeneList[[i]]
        #transform uppercase 
        currentblock=toupper(currentblock)
        if (length(currentblock)>1){
          #do the GO with hypergeometric test
          res=GOcalc(gene=currentblock,terms=termsdt,minsize=minsizeGO,maxsize=maxsizeGO,padj_method="BH")
        }else{
          res=NULL
        }
        return(res)
      },mc.cores=nc)      
    }


    names(reslist)=names(inputGeneList)
    #union of all terms resulting from all ROIs (if ROIs). If null, padj is the max(1) for each term.
    #if all NULL, return problem with sweetalert
    unionterms=lapply(reslist,rownames)
    unionterms=unlist(unionterms)
    unionterms=unique(unionterms)
    
    #union of all blocks in single df. Terms not present: put NA as genes, 0 as ratio, 1 as padj
    #collector of all dfs
    final_list=as.list(rep(NA,length(reslist)))
    #names(final_list)=names(reslist)
    #for loop is enough efficient, because at most we have 10/20 ROIs...
    tofillNULL=c("NA",0,1)
    for(i in 1:length(reslist)){
      #for each ROI, fill entire unionterms, if NULL, fill with NA, 0, 1
      #extract info from result
      tempres=reslist[[i]]
      tempdf=data.frame(matrix(rep(tofillNULL,length(unionterms)),byrow=TRUE,ncol=3))
      rownames(tempdf)=unionterms
      tempdf[,1]=as.character(tempdf[,1])
      tempdf[,2]=as.numeric(tempdf[,2])
      tempdf[,3]=as.numeric(tempdf[,3])
      if(!is.null(tempres)){
        #fill if the term is present with real values
        partres=reslist[[i]][,c(2,7,10)]
        partres[,1]=as.character(partres[,1])
        colnames(partres)=colnames(tempdf)=paste(names(reslist)[i],c("Genes","gene_ratio","padj"),sep="_")
        pos=match(unionterms,rownames(reslist[[i]]))
        pos2=pos[!is.na(pos)]
        tempdf[!is.na(pos),]=partres[pos2,]

      }
      final_list[[i]]=tempdf
    }
    finaldf=do.call("cbind",final_list)



    ######## ORDERING THE MATRIX ######




    mat=finaldf
    #extract -log10 padj values of the matrix 
    #this block can remain. It is removing NA,infinite from matrix and creating padj matrix
    #the important thing is that the order keeps the same
    partdfpadj=mat[,grepl("_padj$",colnames(mat)),drop=FALSE]
    partdfpadj=-log10(partdfpadj)
    #remove inf and NA values
    posNA=apply(partdfpadj,1,function(i){any(is.na(i))})
    posInf=apply(partdfpadj,1,function(i){any(is.infinite(i))})
    partdfpadj=partdfpadj[!posNA & !posInf,,drop=FALSE]
    mat=mat[!posNA & !posInf,,drop=FALSE]

    ########################################################################
	  ########################################################################
	  ########################################################################
	  #in this part rank or cluster the matrix (put the desired ordering according to the inputs)
	  chooseOrderingGO=input$chooseOrderingGO
	  clustertypeGO=input$clustertypeGO
    


	  #if clustering is valid and selected,
	  #respond interactively to the kind of clustering,
	  #otherwise, rank


	  RANKING=TRUE
	  ## if number of genelists are <=2 , only ranking! clustering won't make sense
	  if(ncol(partdfpadj)>1 &nrow(partdfpadj)>2 & isvalid(chooseOrderingGO)){
	    if(chooseOrderingGO=="clustering"){
	      RANKING=FALSE
	    }else{
	      RANKING=TRUE          
	    }
	  }else{
	    RANKING=TRUE
	  }

	  
	  if(RANKING){
	    #ranking. We don't have a matrix, we will do simple barplot
	    # maxvals=apply(partdfpadj,1,max)
	    # ord=order(-maxvals)
	    print ("ranking")
		  maxvals=apply(partdfpadj,1,max)
		  ord=order(-maxvals)
	  }else{
  		clustnumGO=input$clustnumGO
  		clustrandomstartsGO=input$clustrandomstartsGO
  		clustnumiterationsGO=input$clustnumiterationsGO
  		distmethodGO=input$distmethodGO
  		clustmethodGO=input$clustmethodGO
	  	if (!isvalid(clustertypeGO)){
	  	  print ("clsuter not valid")
	  	  return()
	  	}
	  	if( !(isvalid(clustnumGO)&isvalid(clustrandomstartsGO)&isvalid(clustnumiterationsGO) ) &
	  		 !(isvalid(distmethodGO)&isvalid(clustmethodGO))){
	  	  print ("cluster parameters not valid")
	  	  return()	
	  	}

	  	provv1=apply(partdfpadj, 2, as.list)
	  	matlist=lapply(provv1,as.matrix)
			#here,clustering. Use Kmeans or hierarchical, according to what has been chosen
	    if(clustertypeGO=="kmean"){
	      #do kmean
	      print ("kmean clustering")
	      clustobj=clusterMatrixKmeans(matlist=matlist,clustinds=1:length(matlist),numberclusters=clustnumGO,
	      														startingpoints=clustrandomstartsGO,iter=clustnumiterationsGO)

	    }else{
	      #do hierarchical
	      print ("h clustering")
	      clustobj=clusterMatrix(matlist=matlist,clustinds=1:length(matlist),distmethod=distmethodGO,clustmethod=clustmethodGO)
	    	
	    } 
	  ord=clustobj$ord
	  }
	  ########################################################################
	  ########################################################################
	  ########################################################################
	  #whatever os the ordering (ranking, kmean, hclust..), order the matrix!
	  partdfpadj=partdfpadj[ord,,drop=FALSE]
	  mat=mat[ord,,drop=FALSE]


	  #save variables to be used later in other reactive contexts
    toplot$GOvariables$tempmat=mat
    toplot$GOvariables$temppadj=partdfpadj
    toplot$GOvariables$termsdt=termsdt

  }else{
    #no geneset selected
    sendSweetAlert(
      session = session,
      title = "No geneset selected",
      text = "You have to select at least one geneset from the menu",
      type = "error"
    )     
  }

},ignoreInit=TRUE)



##observer for GO plot/table from toplot$GOvariables$tempmat matrix
observe({

  #here the reactive variables of the observer. If they are not satisfied,
  #return immediately


  #observer using reactive parameters
  if(!is.null(toplot$GOvariables$tempmat)){
    mat=toplot$GOvariables$tempmat
    partdfpadj=toplot$GOvariables$temppadj
    padjthresh=10^(-(input$log10padjGO))
    generatiothresh=input$quantileGeneRatioGO
    topN=input$topNGO
    
    #here,filter the results according to parameters, and topN.
    #topN: keep the topN hits based on the best padj in all the ROIs queried.
    # padjthresh: save hit if any of ROI have a padj that is ok for this threshold
    #maybe one of the 2 can be removed if text of genesets in the plot is reduced
    pos_tokeep=filterGOres(GOres=mat,padjthresh=padjthresh,generatiothresh=generatiothresh,topN=topN)
    mat=mat[pos_tokeep,,drop=FALSE]
    partdfpadj=partdfpadj[pos_tokeep,,drop=FALSE]

    #reset selected genes from the previous analysis
    output$GenesClicked<-renderText({NULL})
    output$showGenesClicked<-renderUI({NULL}) 
    output$showTermClicked<-renderUI({NULL})
    
    #if real matrix, heatmap, otherwise barplot
    if(min(dim(mat))>0){
      
      #here we should filter also padj matrix... or recreate it from filtered matrix...
      toplot$GOvariables$partdfpadj=partdfpadj
      toplot$GOvariables$completemat=mat

      if(ncol(partdfpadj)>1){
        #multiplying factor for height(n ROIs,=>ncol)
        if(ncol(partdfpadj)<10){
          factormult=10/ncol(partdfpadj)
        }else{
          factormult=1
        }
        #multiplying factor for width(n ontologies,=>nrow)
        if(nrow(partdfpadj)<30){
          factormult2=30/nrow(partdfpadj)
        }else{
          factormult2=1
        }

        trasp=as.matrix(partdfpadj)


        output$plotOntology<-renderPlot({
          colorpal=input$colorScaleGO 
          scaleQuantile=input$scaleQuantileGO
          brk=201
          trasp[trasp>quantile(trasp,scaleQuantile)]=quantile(trasp,scaleQuantile)
          colorsplitted=strsplit(colorpal,split="_")[[1]]
          palette=colorRampPalette(colorsplitted)(n=brk-1)
          par(mar = rep(0, 4))
          image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", 
                ylab = "",col=palette,xlim=c(0,nrow(trasp)*factormult2),ylim=c(0,ncol(trasp)*factormult),useRaster=FALSE  )
          #lines for the columns
          for(csep in 1:nrow(trasp)){
            rect(xleft = csep , ybottom = rep(0,length(csep)), 
              xright = csep  + 0.01, 
              ytop = rep(ncol(trasp), csep), lty = 1, lwd = 1, 
              col = "black", border = "black"
            )
          }
          #lines for the rows
          for (csep2 in 1:ncol(trasp)){
            rect(xleft = rep(0,length(csep2))  , ybottom = csep2 , 
              xright = rep(nrow(trasp) , csep2), 
              ytop = csep2 +0.01, lty = 1, lwd = 1, 
              col = "black", border = "black")
          }        
        })

        #left names of ROIs
        output$plotMaterialLeft<-renderPlot({

          par(mar = rep(0, 4))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i",yaxs="i")
          #have to start from x coordinate of the half of first cell to half of the last
          numbersample=ncol(partdfpadj)
          if(numbersample<10){
            #divide by 10 and multiply by number sample
            newmax=(1/10)*numbersample
            halfcellwidthmin=(newmax/numbersample)/2
            halfcellwidthmax=newmax-halfcellwidthmin

          }else{
            halfcellwidthmin=(1/numbersample)/2
            halfcellwidthmax=1-halfcellwidthmin
          }
          cextext=0.8
          text(y=seq(halfcellwidthmin, halfcellwidthmax ,length.out=numbersample),labels=rev(colnames(partdfpadj)),x=rep(1,numbersample),cex=cextext,adj=1  )
        })


        #plot color legend (only if heatmap)
        output$colorScaleGO<-renderPlot({    
          scaleQuantile=input$scaleQuantileGO
          trasp[trasp>quantile(trasp,scaleQuantile)]=quantile(trasp,scaleQuantile)          
          brk=201
          palette_col=input$colorScaleGO 
            
            if(palette_col=="rainbow"){
              my_palette <- rev(colorRampPalette(brewer.pal(9, 'Spectral'))(n = brk-1))
            }else{
              colorsplitted=strsplit(palette_col,split="_")[[1]]
              my_palette <- colorRampPalette(colorsplitted)(n = brk-1 )
            }
            color.bar2(my_palette, min=round(0),max=round(max(trasp),2))                    

        })
        output$saveheatmapGO<-renderUI({downloadButton('saveheatmapGObutton', 'Get PDF')})
      }else{

        #left axis
        output$plotMaterialLeft<-renderPlot({
       
          if(nchar(round(max(partdfpadj[,1])))>2){
            maxval=max(partdfpadj[,1])
            labs=seq(0, maxval, round(max(partdfpadj[,1],1)/8))
          }else{
            maxval=round(max(partdfpadj[,1]),3)
            labs=seq(0, maxval , signif(max(partdfpadj[,1],1)/8,digits=2)   )
          }
          #labs=c(labs,signif(max(partdfpadj[,1]),digits=2))
          maxpoint=labs[length(labs)] /maxval
          ats=seq(0,maxpoint,(1/(length(labs))))

          ats[1]=ats[1]+0.01
          ats[length(ats)]=ats[length(ats)]-0.01
          par(mar = c(0,0,0,0))
          plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i", yaxs="i")

          axis(side=2, labels=labs,at=ats,pos=.98)
          text(x=0.80,y=0.5,las=2,label="-log10 padj",srt = 90)
        })

        output$plotOntology<-renderPlot({
          par(mar = c(0,0,0,0))
          if(nrow(partdfpadj)<30){
            factormult2=30/nrow(partdfpadj)
          }else{
            factormult2=1
          }
          
          barplot(partdfpadj[,1],names=rownames(partdfpadj),las=2,ylab="-log10 padj",ylim=c(0,max(partdfpadj[,1])),xlim=c(0,nrow(partdfpadj)*factormult2),xaxt = 'n', yaxt = 'n',xaxs="i",space=rep(0,nrow(partdfpadj)))
        })  

        #no color scale if barplot
        output$colorScaleGO<-renderPlot({NULL})
        output$saveheatmapGO<-renderUI({downloadButton('savebarplotGObutton', 'Get PDF')})
      }

      #prepare table to show and download
      toshow=mat  
      toshow$Term=rownames(mat)
      #rerder cols, Gene list must be at the end
      pospadj=which(grepl("_padj$",colnames(toshow)))
      posratio=which(grepl("_gene_ratio$",colnames(toshow)))
      posgene=which(grepl("_Genes$",colnames(toshow)))
      posterm=which(grepl("Term$",colnames(toshow)))
      newpos=c(posterm,pospadj,posratio,posgene)
      toshow=toshow[,newpos ]
      toplot$GOvariables$filtereddf=toshow


      output$tableOntology <- renderDataTable({
        toshow
      },options=list(pagingType = "simple",pageLength = 10,lengthChange=FALSE,searching=TRUE,autowidth=FALSE )   )

      output$textNameGO <- renderPlot({
        par(mar = rep(0, 4))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i")
        #have to start from x coordinate of the half of first cell to half of the last
        numbersample=nrow(partdfpadj)
        if(numbersample<30){
          #divide by 10 and multiply by number sample
          newmax=(1/30)*numbersample
          halfcellwidthmin=(newmax/numbersample)/2
          halfcellwidthmax=newmax-halfcellwidthmin
          cextext=0.8
        }else{
          halfcellwidthmin=(1/numbersample)/2
          halfcellwidthmax=1-halfcellwidthmin
          coeff=(nrow(partdfpadj)/30 )
          coeff2=(1-(30/nrow(partdfpadj)))/3
          cextext=0.8/coeff +coeff2
        }
        
        text(x=seq(halfcellwidthmin, halfcellwidthmax ,length.out=numbersample),labels=rownames(partdfpadj),y=rep(0.98,numbersample),srt=90,cex=cextext,adj=1  )     
      })

      output$tableGOdownloadButton<-renderUI({downloadButton('saveGOdataTable', 'Download data')})
      
    }else{
      output$plotOntology<-renderPlot({NULL})
      output$tableOntology <- renderDataTable({NULL})
      output$textNameGO <- renderPlot({NULL})
      output$tableGOdownloadButton<-renderUI({NULL})
      output$plotMaterialLeft<-renderPlot({NULL})
      output$colorScaleGO<-renderPlot({NULL})
      output$saveheatmapGO<-renderUI({NULL})
    }    
  }else{
    output$plotOntology<-renderPlot({NULL})
    output$tableOntology <- renderDataTable({NULL})   
    output$textNameGO <- renderPlot({NULL}) 
    output$tableGOdownloadButton<-renderUI({NULL})
    output$plotMaterialLeft<-renderPlot({NULL})
    output$colorScaleGO<-renderPlot({NULL})
    output$saveheatmapGO<-renderUI({NULL})
  }
  
})



#observer for download data button
output$saveGOdataTable<- downloadHandler(
  filename=function() {
      paste('GeneOntology_data.xls', sep='')
  },
  content=function(file) {

      write.table(toplot$GOvariables$filtereddf,file=file,row.names=FALSE,sep="\t",quote=FALSE,col.names=TRUE   ) 
  } 
  
)




output$saveheatmapGObutton<- downloadHandler(
  filename=function() {
      paste('Heatmap_GO.pdf', sep='')
  },
  content=function(file) {
    partdfpadj=toplot$GOvariables$partdfpadj    
    #multiplying factor for height(n ROIs,=>ncol)
    if(ncol(partdfpadj)<10){
      factormult=10/ncol(partdfpadj)
    }else{
      factormult=1
    }
    #multiplying factor for width(n ontologies,=>nrow)
    if(nrow(partdfpadj)<30){
      factormult2=30/nrow(partdfpadj)
      cextext=1
    }else{
      factormult2=1
      coeff=(nrow(partdfpadj)/30 )
      coeff2=(1-(30/nrow(partdfpadj)))/3
      cextext=1/coeff +coeff2
    }

    trasp=as.matrix(partdfpadj)
    colorpal=input$colorScaleGO 
    scaleQuantile=input$scaleQuantileGO
    brk=201
    trasp[trasp>quantile(trasp,scaleQuantile)]=quantile(trasp,scaleQuantile)
    colorsplitted=strsplit(colorpal,split="_")[[1]]
    palette=colorRampPalette(colorsplitted)(n=brk-1)

    pdf(file,width=10,height=8)
    block1=c(rep(1,9),2)
    block2=c(rep(3,9),4)
    mlay=matrix( c(block1,rep(block2,9)),ncol=10,nrow=10,byrow=TRUE)
    layout(mlay)
    par(mar=rep(0,4))
    plot.new()

    #draw color scale
    color.bar2(palette, min=round(0),max=round(max(trasp),2),margins=c(2,1.2,2,1.2))  

    par(mar = c(22,22,0,0))
    image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", 
          ylab = "",col=palette,xlim=c(0,nrow(trasp)*factormult2),ylim=c(0,ncol(trasp)*factormult),useRaster=FALSE  )
    #lines for the columns
    for(csep in 1:nrow(trasp)){
      rect(xleft = csep , ybottom = rep(0,length(csep)), 
        xright = csep  + 0.01, 
        ytop = rep(ncol(trasp), csep), lty = 1, lwd = 1, 
        col = "black", border = "black"
      )
    }
    #lines for the rows
    for (csep2 in 1:ncol(trasp)){
      rect(xleft = rep(0,length(csep2))  , ybottom = csep2 , 
        xright = rep(nrow(trasp) , csep2), 
        ytop = csep2 +0.01, lty = 1, lwd = 1, 
        col = "black", border = "black")
    }    
    #draw axes
    
    axis(1,at=(1:nrow(trasp)-0.5),labels=rownames(trasp),las= 2,tick=FALSE,cex.axis=cextext )
    axis(2,at=(1:ncol(trasp)-0.5),labels=rev(colnames(trasp)),las= 1,tick=FALSE )
    par(mar=rep(0,4))
    plot.new()
  

    dev.off()

  } 
)



output$savebarplotGObutton<- downloadHandler(
  filename=function() {
      paste('Barplot_GO.pdf', sep='')
  },
  content=function(file) {
    partdfpadj=toplot$GOvariables$partdfpadj 
    pdf(file,height=8,width=10)
    par(mar = c(15,6,0,0))
    if(nrow(partdfpadj)<30){
      factormult2=30/nrow(partdfpadj)
    }else{
      factormult2=1
    }
    
    barplot(partdfpadj[,1],names=rownames(partdfpadj),las=2,ylab="-log10 padj",xlim=c(0,nrow(partdfpadj)*factormult2),space=rep(0,nrow(partdfpadj)))

    dev.off()
  } 
)














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
  #input$coresPredefPipeline -> number of cores for enrichment predef.pipeline association
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
      text = "Transcripts cannot be used, because they cannot be resized around their midpoint",
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
    paths_summit=BAMvariables$listBAM[pos2]
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
      currentenrichpath=BAMvariables$listBAM[pos2]
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


  bams=getBAMlist(roi)
  if(length(bams)>0){
    newBAMlist=list()
    for(i in 1:length(bams)){
      newBAMlist[[i]]=bams[[i]][selectedPositions]
    }
    names(newBAMlist)=names(bams)            
  }else{
    newBAMlist=list()
  }
            
  newfix=getFixed(roi)[selectedPositions]
  
  perc=round(quant/length(rangefromROI)*100,2)
  toadd=paste("ROI sampled (kept ",quant,"/",length(rangefromROI),"ranges,",perc,"%)")
  newSource=c(oldSource,list(toadd))

  ROI_after_sample=new("RegionOfInterest",
                                            name=input$ROInamePredefPipeline,
                                            range=newrange,
                                            fixed=newfix,
                                            BAMlist=newBAMlist,
                                            flag=newflag,
                                            source=newSource)

  # logvariables$msg[[length(logvariables$msg)+1]]= paste('Sampling of ',input$selectROIpredefPipeline,' (kept ',perc,'%, ',quant,'/',length(rangefromROI),'ranges)<br>',sep="")  
  # print(paste('Sampling of ',input$selectROIpredefPipeline,' (kept ',perc,'%, ',quant,'/',length(rangefromROI),'ranges)',sep=""))          


  ##################################################################################
  #summit (if yes)
  ##################################################################################
  
  #with elements produced in previous block, do the summit
  #check if enrichment file of the summit already associated
  bams_names=names(newBAMlist)
  if(input$choiceSummitPredefPipeline=="Yes"){
    if(!input$selectBAMsummitPredefPipeline%in%bams_names){
      #associate BAM/WIG
      singlecover=cover(Object=ROI_after_sample,signalfile=paths_summit,signalfileNorm=paths_summit)
    }else{
      #do nothing, enrichment already in ROI. Use BAM already present to center on summit
      pos2=match(input$selectBAMsummitPredefPipeline,bams_names)
      singlecover=newBAMlist[[pos2]]
    }
    #newrange and singlecover must have same length (both derive from correctly subsampled ROI)
    newrange=summitFromBaseCoverage(Object=newrange,baseCoverageOutput=singlecover)
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
  bamlist=getBAMlist(ROI_after_summit)

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

  pos_positive=as.logical(!strand(roi)=="-")
  pos_negative=as.logical(strand(roi)=="-")
  roi_positive=roi[pos_positive]
  roi_negative=roi[pos_negative]
  fix_positive=fix[pos_positive]
  fix_negative=fix[pos_negative]
  bamlist_positive=lapply(bamlist,function(i) {i[pos_positive]})
  bamlist_negative=lapply(bamlist,function(i) {i[pos_negative]})

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

  #to be checked when BAM files will be available
  if (length(bamlist)>0){
    #if not everything is less than the old size in any position, reset BAMlist!
    if (! (all((fix_positive-oldstarts_positive) >= newstart) & all((oldends_positive-fix_positive)>= newend) )
              | !(all((fix_negative-oldstarts_negative) >= newend) & all((oldends_negative-fix_negative)>= newstart) )    ){
      #reset BAMlist or improve only the delta margins of already existing matrixes
      #ROI_after_summit=resetBAMlist(ROI_after_summit)
      logvariables$msg[[length(logvariables$msg)+1]]= paste('Removed associated BAMs, because new size > old size...<br>',sep="")
      newBAMlist=list()
    }else{
      toadd=paste(toadd,"(kept associated BAMs)")
      #split - and + for BAM resize 
      
      shiftleft_positive= (fix_positive-oldstarts_positive)- newstart
      shiftright_positive=(oldends_positive-fix_positive)-newend


      if (any(pos_positive)){
        bamlistnew_positive=lapply(1:length(bamlist_positive),function(i) {
          #cut also BAM file if some ranges have a resize that will lead to a negative start
          x=bamlist_positive[[i]]
          slicedbam=lapply(1:length(x), function(k) {
                                return(x[[k]] [(shiftleft_positive[k]+1): (length(x[[k]])-shiftright_positive[k]) ])
                                })
          return(slicedbam)
        } )

      }

      shiftleft_negative= (fix_negative-oldstarts_negative)- newend
      shiftright_negative=(oldends_negative-fix_negative)-newstart


      if (any(pos_negative)){
        bamlistnew_negative=lapply(1:length(bamlist_negative),function(i) {
          #cut also BAM file if some ranges have a resize that will lead to a negative start
          x=bamlist_negative[[i]]
          slicedbam=lapply(1:length(x), function(k) {
                                return(x[[k]] [(shiftleft_negative[k]+1): (length(x[[k]])-shiftright_negative[k]) ])
                                })
          return(slicedbam)
        } )                    
      }



      #join bamlistnew_positive and bamlistnew_negative back to bamlist
      bamlistnew=lapply(1:length(bamlist),function(i) {
        x=bamlist[[i]]
        if (any(pos_positive)){
          x[pos_positive]=bamlistnew_positive[[i]]
        }
        
        if (any(pos_negative)){
          x[pos_negative]=bamlistnew_negative[[i]]
        }
        return(x)
      })

      names(bamlistnew)=names(bamlist)
      newBAMlist=bamlistnew
      #ROI_after_summit=setBAMlist(ROI_after_summit,bamlistnew)

      #x is single list of baseCoverage (BAM)
      logvariables$msg[[length(logvariables$msg)+1]]= paste('Conserved associated BAMs, because new size < old size...<br>',sep="")

    }
  }else{
    newBAMlist=list()
  }

  newSource=c(oldSource,list(toadd))
  ROI_after_resize=new("RegionOfInterest",
                                name=input$ROInamePredefPipeline,
                                range=newrange,
                                fixed=newfix,
                                BAMlist=newBAMlist,
                                flag=newflag,
                                source=newSource) 

  # if(length(unique(fix_positive-oldstarts_positive))==1 & length(unique(oldends_positive-fix_positive))==1){
  #   logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInamePredefPipeline,' ROI [-',input$sliderUpstreamPredefPipeline,'; +',input$sliderDownstreamPredefPipeline,'] from ',input$selectROIpredefPipeline,' ROI [-',(fix_positive-oldstarts_positive)[1],'; +',(oldends_positive-fix_positive)[1],']<br>',sep="")
  #   print(paste('Created ',input$ROInamePredefPipeline,' ROI [-',input$sliderUpstreamPredefPipeline,'; +',input$sliderDownstreamPredefPipeline,'] from ',input$selectROIpredefPipeline,' ROI [-',(fix_positive-oldstarts_positive)[1],'; +',(oldends_positive-fix_positive)[1],']',sep=""))
  # }else{
  #   logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$ROInamePredefPipeline,' ROI [-',input$sliderUpstreamPredefPipeline,'; +',input$sliderDownstreamPredefPipeline,'] from ',input$selectROIpredefPipeline,' ROI with original widths<br>',sep="")
  #   print(paste('Created ',input$ROInamePredefPipeline,' ROI [-',input$sliderUpstreamPredefPipeline,'; +',input$sliderDownstreamPredefPipeline,'] from ',input$selectROIpredefPipeline,' ROI with original widths',sep=""))
  # }




  ##################################################################################
  #enrichment association
  ##################################################################################
  ncores=input$coresPredefPipeline
  enrichments_selected=input$enrichAllPredefPipeline
  enrichments_alreadyAssociated=getBAMlist(ROI_after_resize)
  
  #determine which enrichments are already associated
  enrichments_toAssociate=setdiff(enrichments_selected,names(enrichments_alreadyAssociated))
  pos2=match(enrichments_toAssociate,names(BAMvariables$listBAM))
  paths=BAMvariables$listBAM[pos2]  


  #depending if only one operation, lapply or mclapply with n cores
  #IMPORTANT: NEVER TANSCRIPT FLAG
  if(length(enrichments_selected)>=1){
    if(ncores==1 | length(enrichments_toAssociate)==1 | length(paths)==1){
      #single core (lapply)
      finallist=lapply(1:length(paths),function(i) {
        singlecover=cover(Object=ROI_after_resize,signalfile=paths[i],signalfileNorm=paths[i])
        print(paste("Associating ",paths[i]," in single core",sep=""))
        return(singlecover)
      })

    }else{
      #multicore (mclapply)
      finallist=mclapply(1:length(paths),function(i) {
        singlecover=cover(Object=ROI_after_resize,signalfile=paths[i],signalfileNorm=paths[i])
        print(paste("Associating ",paths[i]," in multi core",sep=""))
        return(singlecover)
      },mc.cores=ncores)
    }

    #put names to new enrichments
    names(finallist)=names(paths)
    finalfinallist=c(enrichments_alreadyAssociated,finallist)
    ROI_after_associate=setBAMlist(ROI_after_resize,finalfinallist)
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





