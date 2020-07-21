

#show BAM files to delete
 observe({
 	if (length(BAMvariables$listBAM)>0){
    historylist=unlist(BAMvariables$listBAM)
    bampaths=names(BAMvariables$listBAM)
    finalname=paste(bampaths,"(",historylist,")")
    historylist=as.list(names(BAMvariables$listBAM))
    names(historylist)=finalname
 		updateCheckboxGroupInput(session,inputId="selectedBAMtoRemove",label="",
                                  choices = historylist) 
 	}else{
 		historylist=list()
 		updateCheckboxGroupInput(session,inputId="selectedBAMtoRemove",label="",
                                  choices = historylist) 
 	} 
})


#show BAM files to rename. Show only the name of list, not the basename of the file
observe({
  BAMvariables$listBAM
  if(length(BAMvariables$listBAM)>0){
    historylist=as.list(names(BAMvariables$listBAM))
    names(historylist)=names(BAMvariables$listBAM)
    updateSelectInput(session,inputId="selectedBAMfiletoRename",label="Enrichment to rename:",
                                  choices = historylist) 
  }else{
    historylist=list()
    updateSelectInput(session,inputId="selectedBAMfiletoRename",label="Enrichment to rename:",
                                  choices = historylist) 
  } 
})
