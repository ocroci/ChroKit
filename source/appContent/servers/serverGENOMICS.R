##################################################################
##################################################################
##################################################################
# respond to plot Single evaluation
##################################################################
##################################################################
##################################################################

###react to help buttons:
#parameters
observeEvent(input$msg_singleEvaluation_parameters, {
  boxHelpServer(msg_singleEvaluation_parameters)
})
#distribution
observeEvent(input$msg_singleEvaluation_distribution, {
  boxHelpServer(msg_singleEvaluation_distribution)
})
#width distribution
observeEvent(input$msg_singleEvaluation_widthDistribution, {
  boxHelpServer(msg_singleEvaluation_widthDistribution)
})
#enrichment boxplot
observeEvent(input$msg_singleEvaluation_enrichmentBoxplot, {
  boxHelpServer(msg_singleEvaluation_enrichmentBoxplot)
})
#profile 
observeEvent(input$msg_singleEvaluation_peakProfile, {
  boxHelpServer(msg_singleEvaluation_peakProfile)
})



observeEvent(input$plotSingleEval,{
	set.seed(123)
	if (!is.null(input$ROIchooseSingleEval) & length(ROIvariables$listROI)>0) {
		nomi=unlist(lapply(ROIvariables$listROI,getName))
		
    	if ("promoters"%in% nomi & "transcripts"%in% nomi) {
      		pos=match(input$ROIchooseSingleEval,nomi)
      		roi=uniqueROI(ROIvariables$listROI[[pos]]) 

          roi=unifyStrand(roi)

      		if (!is.null(roi)){
      			ROInameSingleEval=input$ROIchooseSingleEval
      			#calculate overlaps and plot pie/barplot/density of width
            promoter_roi=ROIvariables$listROI[[which(nomi=="promoters")]]
      			promo=getRange(promoter_roi)
      			transc=getRange(ROIvariables$listROI[[which(nomi=="transcripts")]])
    			  range=getRange(roi)
    			  ov_range=countOverlaps(range,promo)
    			  range_promo=range[ov_range>0]
    			  range_ws=range[ov_range==0]
    			  ov_transcripts=countOverlaps(range_ws,transc)
    			  range_intra=range_ws[ov_transcripts>0]
    			  range_inter=range_ws[ov_transcripts==0]
    			  #print (paste(length(range_promo),length(range_intra),length(range_inter)))
    			  elements<-c(length(range_promo),length(range_intra),length(range_inter))
				    perc<-c(round(length(range_promo)/length(range),2)*100,round(length(range_intra)/length(range),2)*100,round(length(range_inter)/length(range),2)*100)
				    label<-c(paste(perc[1],"% (",elements[1],")",sep=''),paste(perc[2],"% (",elements[2],")",sep=''),paste(perc[3],"% (",elements[3],")",sep=''))
            label<-paste(c("Promoters:","Genebody:","Intergenic:"),label,sep="")
				    maintitle=paste("ROI: ",ROInameSingleEval," (",length(range)," total ranges)",sep='')

            toplot$viewDistributionPieSingleEval$Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]
            toplot$viewDistributionPieSingleEval$elements=elements
            toplot$viewDistributionPieSingleEval$maintitle=maintitle
            toplot$viewDistributionPieSingleEval$label=label
            toplot$viewDistributionPieSingleEval$nomi=nomi
            toplot$viewDistributionPieSingleEval$range=range
            toplot$viewDistributionPieSingleEval$range_promo=range_promo
            toplot$viewDistributionPieSingleEval$range_intra=range_intra
            toplot$viewDistributionPieSingleEval$range_inter=range_inter


            #pie plot
   				  output$viewDistributionPieSingleEval<-renderPlot({
              Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]
   					  #improve the lgend position
   					  par(mar=c(7,7,7,7))
   					  pie(elements,col=Colors[2:4],main=maintitle,cex.main=1.2,labels=label,cex=1)
   					  #legend("bottomright",legend=c("promoter","genebody","intergenic"),col=Colors,pch=rep(19,3) )
   				  })



            #barplot
            output$viewDistributionBarSingleEval<-renderPlot({
              Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]
              par(mar=c(10,6,5,5))
              names(elements)=label
              barplot(elements,col=Colors[2:4],main=maintitle,cex.main=1.2,las=2)
              title(ylab = "Interval number", cex.lab = 1,
                        line = 4)
            })
            #download button for PDF pie
            output$saveviewDistributionPieSingleEval=renderUI({downloadButton('saveviewDistributionPieSingleEvalbutton', 'Get PDF')})
            #download button for PDF barplot
            output$saveviewDistributionBarSingleEval=renderUI({downloadButton('saveviewDistributionBarSingleEvalbutton', 'Get PDF')})
            print(paste("Single evaluation of",input$ROIchooseSingleEval))
            #width distribution
            if (length(range)>2){
              
              quant=0.95
              Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]
              cuttedga=range[width(range)<quantile(width(range),quant)]
              cuttedprom=range_promo[width(range_promo)<quantile(width(range_promo),quant)]
              cuttedintra=range_intra[width(range_intra)<quantile(width(range_intra),quant)]
              cuttedinter=range_inter[width(range_inter)<quantile(width(range_inter),quant)]

              if (length(range_promo)>2&length(cuttedga)>2){

                output$widthDistributionSingleEval<-renderPlot( {
                  if (length(range_intra)>2){
                    if(length(range_inter)>2){
                      xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x,density(log2(width(cuttedinter)))$x), 
                              max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x,density(log2(width(cuttedinter)))$x)  )
                      ylim=c(0, max(max(density(log2(width(cuttedprom)))$y),max(density(log2(width(cuttedintra)))$y),max(density(log2(width(cuttedinter)))$y) ) )
                      plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 length",cex.lab=1,cex.axis=1,ylab="Frequency",ylim=ylim,xlim=xlim,
                        main=paste("Interval length (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=2,type="l",col=Colors[1])
                      lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=Colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=2)
                      lines(density(log2(width(cuttedintra)))$x,density(log2(width(cuttedintra)))$y,cex=1,col=Colors[3],ylim=ylim,xlim=xlim,lty=1,lwd=2)
                      lines(density(log2(width(cuttedinter)))$x,density(log2(width(cuttedinter)))$y,cex=1,col=Colors[4],ylim=ylim,xlim=xlim,lty=1,lwd=2)
                
                    }else{
                      xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x), 
                              max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x  ))
                      ylim=c(0, max(max(density(log2(width(cuttedprom)))$y),max(density(log2(width(cuttedintra)))$y) ) )
                      plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 length",cex.lab=1,cex.axis=1,ylab="Frequency",ylim=ylim,xlim=xlim,
                        main=paste("Interval length (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=2,type="l",col=Colors[1])
                      lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=Colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=2)
                      lines(density(log2(width(cuttedintra)))$x,density(log2(width(cuttedintra)))$y,cex=1,col=Colors[3],ylim=ylim,xlim=xlim,lty=1,lwd=2)
                
                    }
                    
                  }else{

                    xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x), 
                              max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x  ))
                    ylim=c(0, max(max(density(log2(width(cuttedprom)))$y) ) )
                    plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 length",cex.lab=1,cex.axis=1,ylab="Frequency",ylim=ylim,xlim=xlim,
                      main=paste("Interval length (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=2,type="l",col=Colors[1])
                    lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=Colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=2)
                  }
                  legend("topright" , c("All intervals","Promoter","Genebody","Intergenic") , col=Colors , lty=rep(1,4),lwd=3,cex=1)
                })
                #download button for width distribution
                output$savewidthDistributionSingleEval<-renderUI({downloadButton('savewidthDistributionSingleEvalbutton', 'Get PDF')})
              
              }else{
                #only total plot, but at least one between promoters, intra, inter must be >0 by definition
                output$widthDistributionSingleEval<-renderPlot({NULL})
                output$savewidthDistributionSingleEval<-renderUI({NULL})
              }

            }else{
              #range has length==0
              output$widthDistributionSingleEval<-renderPlot({NULL})
              output$savewidthDistributionSingleEval=renderUI({NULL})
            }



      			getbam=names(getBAMlist(roi))

            toplot$viewDistributionPieSingleEval$getbam=getbam
            toplot$viewDistributionPieSingleEval$chooseNormalizationSingleEval=input$chooseNormalizationSingleEval

      			if (!is.null(getbam)){
      				pos2=match(input$BAMchooseSingleEval,getbam)
  					  bam_orig=getBAMlist(roi)[[pos2]]
              #bam=unlist(lapply(bam,sum))
  					  #calculate enrichments for boxplots
  					  bam_promo_orig=bam_orig[ov_range>0]
  					  bam_ws_orig=bam_orig[ov_range==0]
  					  bam_intra_orig=bam_ws_orig[ov_transcripts>0]
  					  bam_inter_orig=bam_ws_orig[ov_transcripts==0]

              bam=unlist(lapply(bam_orig,sum))
              bam_promo=unlist(lapply(bam_promo_orig,sum))
              bam_intra=unlist(lapply(bam_intra_orig,sum))
              bam_inter=unlist(lapply(bam_inter_orig,sum))

              toplot$viewDistributionPieSingleEval$bam=bam
              toplot$viewDistributionPieSingleEval$bam_promo=bam_promo
              toplot$viewDistributionPieSingleEval$bam_intra=bam_intra
              toplot$viewDistributionPieSingleEval$bam_inter=bam_inter

              output$enrichmentBoxSingleEval<-renderPlot({
                par(mar=c(8,5,5,5))
                Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]
                if (input$chooseNormalizationSingleEval=="totread"){
                  suppressWarnings(boxplot(bam,bam_promo,bam_intra,bam_inter,notch=TRUE,outline=FALSE,varwidth=TRUE,col=Colors,xaxt="n",ylab="Normalized reads (rpm)"))
                }else{
                  bam=round( bam/length(range),3)
                  bam_promo=round( bam_promo/length(range_promo),3)
                  bam_intra=round( bam_intra/length(range_intra),3)
                  bam_inter=round( bam_inter/length(range_inter),3)
                  suppressWarnings(boxplot(bam,bam_promo,bam_intra,bam_inter,notch=TRUE,outline=FALSE,varwidth=TRUE,col=Colors,xaxt="n",ylab="Read density (rpm/bp)"))
                }
                axis(1,at=1:4,labels=c("All intervals","Promoter","Genebody","Intergenic"),las=2)

              })


              # bam_orig
              # bam_promo_orig
              # bam_intra_orig
              # bam_inter_orig
              listtoprofile=list(bam_orig,bam_promo_orig,bam_intra_orig,bam_inter_orig)
              

              for(i in 1:length(listtoprofile)){
                set.seed(123)
                smp=sample(1:length(listtoprofile[[i]]),floor(length(listtoprofile[[i]])/10),replace=FALSE)
                listtoprofile[[i]]=listtoprofile[[i]][smp]
              }

              #find matrix in bins for promo,intra,inter
              matrixes=list()
              if(input$chooseNormalizationSingleEval=="totread"){
                normmethod=FALSE
                ylab="Normalized reads (rpm)"
              }else{
                normmethod=TRUE
                ylab="Read density (rpm/bp)"
              }
              toplot$viewDistributionPieSingleEval$ylabprofile=ylab
              for(i in 1:length(listtoprofile)){
                if(length(listtoprofile[[i]])>0){
                  matrixes[[i]]=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=listtoprofile[[i]],Nbins=50,Snorm=normmethod)
                }else{
                  matrixes[[i]]=matrix(rep(0,50),nrow=1)
                }
                
              }

              #compact the matrix using mean
              profile_to_plot=lapply(matrixes,function(i) {apply(i,2,mean)})
              names(profile_to_plot)=c("all","promoters","genebody","intergenic")
              toplot$viewDistributionPieSingleEval$profile_to_plot=profile_to_plot
              output$TSSprofileSingleEval<-renderPlot({
                Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]

                mins=Reduce(min,profile_to_plot)
                maxs=Reduce(max,profile_to_plot)

                plot(0,type='n',xaxt="n",ylim=c(mins,maxs),xlim=c(0,50),ylab=ylab,xlab="Genomic Window",
                        main="Shape profile",cex.main=1.2,cex=1,xaxt="n",cex.lab=1.2,cex.axis=1)
                for(i in 1:length(profile_to_plot)){
                  lines(profile_to_plot[[i]],pch=".",lwd=2,col=Colors[i])
                }

                axis(1,at=seq(1,length(profile_to_plot[[1]]),length(profile_to_plot[[1]])/2-1),
                      labels=c( "start","midpoint","end" ),cex.axis=1)

              })


              output$saveboxdataSingleEval=renderUI({downloadButton('saveenrichmentBoxSingleEvaldata', 'Save data')})
              output$saveenrichmentBoxSingleEval=renderUI({downloadButton('saveenrichmentBoxSingleEvalbutton', 'Get PDF')})
              output$saveenrichmentProfileSingleEval=renderUI({downloadButton('saveenrichmentProfileSingleEvalbutton', 'Get PDF')})

      			}else{
      				#roi doesn't have bam...
              output$enrichmentBoxSingleEval<-renderPlot({NULL})
              output$TSSprofileSingleEval<-renderPlot({NULL})
              output$saveenrichmentBoxSingleEval=renderUI({NULL})
              output$saveenrichmentProfileSingleEval=renderUI({NULL})
      				#logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">You have to associate a BAM file to this ROI for the enrichments...<br></font>'
              toplot$viewDistributionPieSingleEval$bam=NULL
              toplot$viewDistributionPieSingleEval$bam_promo=NULL
              toplot$viewDistributionPieSingleEval$bam_intra=NULL
              toplot$viewDistributionPieSingleEval$bam_inter=NULL
              output$saveboxdataSingleEval=renderUI({NULL})
            }
      		}else{
            #not easy to be here. if roi is null, something should go wrong upstream...
            output$viewDistributionPieSingleEval<-renderPlot({NULL})
            output$viewDistributionBarSingleEval<-renderPlot({NULL})
            output$widthDistributionSingleEval<-renderPlot({NULL})
            output$enrichmentBoxSingleEval<-renderPlot({NULL})
            output$TSSprofileSingleEval<-renderPlot({NULL})
            output$saveviewDistributionPieSingleEval=renderUI({NULL})
            output$saveviewDistributionBarSingleEval=renderUI({NULL})
            output$savewidthDistributionSingleEval=renderUI({NULL})
            output$saveenrichmentBoxSingleEval=renderUI({NULL})
            output$saveenrichmentProfileSingleEval=renderUI({NULL})
      			#logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">ROI not found...<br></font>'
            toplot$viewDistributionPieSingleEval$Colors=NULL
            toplot$viewDistributionPieSingleEval$elements=NULL
            toplot$viewDistributionPieSingleEval$maintitle=NULL
            toplot$viewDistributionPieSingleEval$label=NULL
            toplot$viewDistributionPieSingleEval$nomi=NULL
            toplot$viewDistributionPieSingleEval$range=NULL
            toplot$viewDistributionPieSingleEval$range_promo=NULL
            toplot$viewDistributionPieSingleEval$range_intra=NULL
            toplot$viewDistributionPieSingleEval$range_inter=NULL
            toplot$viewDistributionPieSingleEval$getbam=NULL
            output$saveboxdataSingleEval=renderUI({NULL})
      			#roi is null...
      		}

    	}else{
        output$viewDistributionPieSingleEval<-renderPlot({NULL})
        output$viewDistributionBarSingleEval<-renderPlot({NULL})
        output$widthDistributionSingleEval<-renderPlot({NULL})
        output$enrichmentBoxSingleEval<-renderPlot({NULL})
        output$TSSprofileSingleEval<-renderPlot({NULL})
        output$saveviewDistributionPieSingleEval=renderUI({NULL})
        output$saveviewDistributionBarSingleEval=renderUI({NULL})
        output$savewidthDistributionSingleEval=renderUI({NULL})
        output$saveenrichmentBoxSingleEval=renderUI({NULL})
        output$saveenrichmentProfileSingleEval=renderUI({NULL})
    		#logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">promoters/transcripts not found... ask to database<br></font>'
        sendSweetAlert(
          session = session,
          title = "Annotated elements not found",
          text = "Promoters, transcripts and TES of a specific genome assembly not found, but you need them: go to 'Databases' section and choose a genome assembly",
          type = "error"
        )        
        toplot$viewDistributionPieSingleEval$Colors=NULL
        toplot$viewDistributionPieSingleEval$elements=NULL
        toplot$viewDistributionPieSingleEval$maintitle=NULL
        toplot$viewDistributionPieSingleEval$label=NULL
        toplot$viewDistributionPieSingleEval$nomi=NULL
        toplot$viewDistributionPieSingleEval$range=NULL
        toplot$viewDistributionPieSingleEval$range_promo=NULL
        toplot$viewDistributionPieSingleEval$range_intra=NULL
        toplot$viewDistributionPieSingleEval$range_inter=NULL
        toplot$viewDistributionPieSingleEval$getbam=NULL
        output$saveboxdataSingleEval=renderUI({NULL})
      		#print("promoters/transcripts not found... ask to database")
    	}

	}else{
		#roi not found 
    output$viewDistributionPieSingleEval<-renderPlot({NULL})
    output$viewDistributionBarSingleEval<-renderPlot({NULL})
    output$widthDistributionSingleEval<-renderPlot({NULL})
    output$enrichmentBoxSingleEval<-renderPlot({NULL})
    output$TSSprofileSingleEval<-renderPlot({NULL})
    output$saveviewDistributionPieSingleEval=renderUI({NULL})
    output$saveviewDistributionBarSingleEval=renderUI({NULL})
    output$savewidthDistributionSingleEval=renderUI({NULL})
    output$saveenrichmentBoxSingleEval=renderUI({NULL})
    output$saveenrichmentProfileSingleEval=renderUI({NULL})
    toplot$viewDistributionPieSingleEval$Colors=NULL
    toplot$viewDistributionPieSingleEval$elements=NULL
    toplot$viewDistributionPieSingleEval$maintitle=NULL
    toplot$viewDistributionPieSingleEval$label=NULL
    toplot$viewDistributionPieSingleEval$nomi=NULL
    toplot$viewDistributionPieSingleEval$range=NULL
    toplot$viewDistributionPieSingleEval$range_promo=NULL
    toplot$viewDistributionPieSingleEval$range_intra=NULL
    toplot$viewDistributionPieSingleEval$range_inter=NULL
    toplot$viewDistributionPieSingleEval$getbam=NULL
    output$saveboxdataSingleEval=renderUI({NULL})
	}

},ignoreInit=TRUE)






output$saveviewDistributionPieSingleEvalbutton<- downloadHandler(
  filename=function() {
      paste('Pie_location_single.pdf', sep='')
  },
  content=function(file) {
      pdf(file)
      par(mar=c(7,7,7,7))
      pie(toplot$viewDistributionPieSingleEval$elements,col=toplot$viewDistributionPieSingleEval$Colors[2:4],
            main=toplot$viewDistributionPieSingleEval$maintitle,cex.main=1.2,labels=toplot$viewDistributionPieSingleEval$label,cex=1)
      dev.off()
  } 
)


#respond for PDF button for barplot
output$saveviewDistributionBarSingleEvalbutton<- downloadHandler(
  filename=function() {
      paste('Bar_location_single.pdf', sep='')
  },
  content=function(file) {
      pdf(file)
      par(mar=c(10,6,5,5))
      names(toplot$viewDistributionPieSingleEval$elements)=toplot$viewDistributionPieSingleEval$label
      barplot(toplot$viewDistributionPieSingleEval$elements,col=toplot$viewDistributionPieSingleEval$Colors[2:4],main=toplot$viewDistributionPieSingleEval$maintitle,cex.main=1.2,las=2)
      title(ylab = "Interval number", cex.lab = 1,line = 4)
      dev.off()
  } 
)



#button for downloading PDF of width distribution
output$savewidthDistributionSingleEvalbutton<- downloadHandler(
  filename=function() {
      paste('Width_distribution_single.pdf', sep='')
  },
  content=function(file) {
        pdf(file)
        quant=0.95

        cuttedga=toplot$viewDistributionPieSingleEval$range[width(toplot$viewDistributionPieSingleEval$range)<quantile(width(toplot$viewDistributionPieSingleEval$range),quant)]
        cuttedprom=toplot$viewDistributionPieSingleEval$range_promo[width(toplot$viewDistributionPieSingleEval$range_promo)<quantile(width(toplot$viewDistributionPieSingleEval$range_promo),quant)]
        cuttedintra=toplot$viewDistributionPieSingleEval$range_intra[width(toplot$viewDistributionPieSingleEval$range_intra)<quantile(width(toplot$viewDistributionPieSingleEval$range_intra),quant)]
        cuttedinter=toplot$viewDistributionPieSingleEval$range_inter[width(toplot$viewDistributionPieSingleEval$range_inter)<quantile(width(toplot$viewDistributionPieSingleEval$range_inter),quant)]

        if (length(toplot$viewDistributionPieSingleEval$range_promo)>2&length(cuttedga)>2){
          if (length(toplot$viewDistributionPieSingleEval$range_intra)>2){
            if(length(toplot$viewDistributionPieSingleEval$range_inter)>2){
              xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x,density(log2(width(cuttedinter)))$x), 
                      max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x,density(log2(width(cuttedinter)))$x)  )
              ylim=c(0, max(max(density(log2(width(cuttedprom)))$y),max(density(log2(width(cuttedintra)))$y),max(density(log2(width(cuttedinter)))$y) ) )
              plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 length",cex.lab=1,cex.axis=1,ylab="Frequency",ylim=ylim,xlim=xlim,
                main=paste("Interval length (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=2,type="l",col=toplot$viewDistributionPieSingleEval$Colors[1])
              lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=toplot$viewDistributionPieSingleEval$Colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=2)
              lines(density(log2(width(cuttedintra)))$x,density(log2(width(cuttedintra)))$y,cex=1,col=toplot$viewDistributionPieSingleEval$Colors[3],ylim=ylim,xlim=xlim,lty=1,lwd=2)
              lines(density(log2(width(cuttedinter)))$x,density(log2(width(cuttedinter)))$y,cex=1,col=toplot$viewDistributionPieSingleEval$Colors[4],ylim=ylim,xlim=xlim,lty=1,lwd=2)
              legend("topright" , c("All ranges","Promoter","Genebody","Intergenic") , col=toplot$viewDistributionPieSingleEval$Colors , lty=rep(1,4),lwd=3,cex=1)
              dev.off()
            }else{
              xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x), 
                      max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x  ))
              ylim=c(0, max(max(density(log2(width(cuttedprom)))$y),max(density(log2(width(cuttedintra)))$y) ) )
              plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 length",cex.lab=1,cex.axis=1,ylab="Frequency",ylim=ylim,xlim=xlim,
                main=paste("Interval length (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=2,type="l",col=toplot$viewDistributionPieSingleEval$Colors[1])
              lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=toplot$viewDistributionPieSingleEval$Colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=2)
              lines(density(log2(width(cuttedintra)))$x,density(log2(width(cuttedintra)))$y,cex=1,col=toplot$viewDistributionPieSingleEval$Colors[3],ylim=ylim,xlim=xlim,lty=1,lwd=2)
              legend("topright" , c("All ranges","Promoter","Genebody","Intergenic") , col=toplot$viewDistributionPieSingleEval$Colors , lty=rep(1,4),lwd=3,cex=1)
              dev.off()
            }
            
          }else{

            xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x), 
                      max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x  ))
            ylim=c(0, max(max(density(log2(width(cuttedprom)))$y) ) )
            plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 length",cex.lab=1,cex.axis=1,ylab="Frequency",ylim=ylim,xlim=xlim,
              main=paste("Interval length (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=2,type="l",col=toplot$viewDistributionPieSingleEval$Colors[1])
            lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=toplot$viewDistributionPieSingleEval$Colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=2)
            legend("topright" , c("All ranges","Promoter","Genebody","Intergenic") , col=toplot$viewDistributionPieSingleEval$Colors , lty=rep(1,4),lwd=3,cex=1)
            dev.off()
          }      
        }else{
          #only total plot, but at least one between promoters, intra, inter must be >0 by definition
          dev.off()
        }       
  } 
)



#button for saving boxplot of SingleEval
output$saveenrichmentBoxSingleEvalbutton<- downloadHandler(
  filename=function() {
      paste('Boxplot_single.pdf', sep='')
  },
  content=function(file) {
      pdf(file)
      par(mar=c(8,5,5,5))
      if (input$chooseNormalizationSingleEval=="totread"){
        suppressWarnings(boxplot(toplot$viewDistributionPieSingleEval$bam,toplot$viewDistributionPieSingleEval$bam_promo,toplot$viewDistributionPieSingleEval$bam_intra,toplot$viewDistributionPieSingleEval$bam_inter,notch=TRUE,outline=FALSE,varwidth=TRUE,col=toplot$viewDistributionPieSingleEval$Colors,xaxt="n",ylab="Normalized reads (rpm)"))
      }else{
        bam=round( toplot$viewDistributionPieSingleEval$bam/length(toplot$viewDistributionPieSingleEval$range),3)
        bam_promo=round( toplot$viewDistributionPieSingleEval$bam_promo/length(toplot$viewDistributionPieSingleEval$range_promo),3)
        bam_intra=round( toplot$viewDistributionPieSingleEval$bam_intra/length(toplot$viewDistributionPieSingleEval$range_intra),3)
        bam_inter=round( toplot$viewDistributionPieSingleEval$bam_inter/length(toplot$viewDistributionPieSingleEval$range_inter),3)
        suppressWarnings(boxplot(bam,bam_promo,bam_intra,bam_inter,notch=TRUE,outline=FALSE,varwidth=TRUE,col=toplot$viewDistributionPieSingleEval$Colors,xaxt="n",ylab="Read density (rpm/bp)"))
      }
      axis(1,at=1:4,labels=c("All ranges","Promoter","Genebody","Intergenic"),las=2)
      dev.off()
  } 
)

#observer for download of xls table of the boxplot data of SingleEval
output$saveenrichmentBoxSingleEvaldata<- downloadHandler(
  filename=function() {
      paste('boxplot_data.xls', sep='')
  },
  content=function(file) {
      if(input$chooseNormalizationSingleEval=="totread"){
        maxval=max(length(toplot$viewDistributionPieSingleEval$bam),
                  length(toplot$viewDistributionPieSingleEval$bam_promo),
                  length(toplot$viewDistributionPieSingleEval$bam_intra),
                  length(toplot$viewDistributionPieSingleEval$bam_inter))
        arr=matrix(rep("",maxval*4),ncol=4)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam),1]=toplot$viewDistributionPieSingleEval$bam
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_promo),2]=toplot$viewDistributionPieSingleEval$bam_promo
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_intra),3]=toplot$viewDistributionPieSingleEval$bam_intra
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_inter),4]=toplot$viewDistributionPieSingleEval$bam_inter
        
      }else{
        maxval=max(length(toplot$viewDistributionPieSingleEval$bam),
                  length(toplot$viewDistributionPieSingleEval$bam_promo),
                  length(toplot$viewDistributionPieSingleEval$bam_intra),
                  length(toplot$viewDistributionPieSingleEval$bam_inter))
        arr=matrix(rep("",maxval*4),ncol=4)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam),1]=toplot$viewDistributionPieSingleEval$bam/length(toplot$viewDistributionPieSingleEval$range)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_promo),2]=toplot$viewDistributionPieSingleEval$bam_promo/length(toplot$viewDistributionPieSingleEval$range_promo)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_intra),3]=toplot$viewDistributionPieSingleEval$bam_intra/length(toplot$viewDistributionPieSingleEval$range_intra)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_inter),4]=toplot$viewDistributionPieSingleEval$bam_inter/length(toplot$viewDistributionPieSingleEval$range_inter)
      }
      colnames(arr)=c("all","promoters","genebody","intergenic")
      write.table(arr,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  } 
  
)

output$saveenrichmentProfileSingleEvalbutton<- downloadHandler(
  filename=function() {
      paste('Profile_single.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    par(mar=c(6,6,3,3))
    Colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]]
    mins=Reduce(min,toplot$viewDistributionPieSingleEval$profile_to_plot)
    maxs=Reduce(max,toplot$viewDistributionPieSingleEval$profile_to_plot)

    plot(0,type='n',xaxt="n",ylim=c(mins,maxs),xlim=c(0,50),ylab=toplot$viewDistributionPieSingleEval$ylabprofile,xlab="Genomic Window",
                        main="Shape profile",cex.main=1.2,cex=1,xaxt="n",cex.lab=1.2,cex.axis=1)
    for(i in 1:length(toplot$viewDistributionPieSingleEval$profile_to_plot)){
      lines(toplot$viewDistributionPieSingleEval$profile_to_plot[[i]],pch=".",lwd=2,col=Colors[i])
    }

    axis(1,at=seq(1,length(toplot$viewDistributionPieSingleEval$profile_to_plot[[1]]),length(toplot$viewDistributionPieSingleEval$profile_to_plot[[1]])/2-1),
                      labels=c( "start","midpoint","end" ),cex.axis=1)
    dev.off()
  } 
)



##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
# respond to plot Pairwise overlaps
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################

###react to help buttons:
#parameters
observeEvent(input$msg_pairwiseOverlaps_parameters, {
  boxHelpServer(msg_pairwiseOverlaps_parameters)
})
#overlap
observeEvent(input$msg_pairwiseOverlaps_overlap, {
  boxHelpServer(msg_pairwiseOverlaps_overlap)
})
#overlap and enrichment
observeEvent(input$msg_pairwiseOverlaps_overlapAndEnrichment, {
  boxHelpServer(msg_pairwiseOverlaps_overlapAndEnrichment)
})


observeEvent(input$plotCmp,{
  #verify to have at least input$ROI2chooseCmp and input$ROI1chooseCmp 
  if (length(ROIvariables$listROI)>0&length(input$ROI2chooseCmp)>0 & length(input$ROI1chooseCmp)>0 & !is.na(input$minOverlapCmp) ){
    n1=input$ROI1chooseCmp
    n2=input$ROI2chooseCmp

    if(is.null(input$BAM1chooseCmp)){
      BAM1_choose=""
    }else{
      BAM1_choose=input$BAM1chooseCmp
    }
    if(is.null(input$BAM2chooseCmp)){
      BAM2_choose=""
    }else{
      BAM2_choose=input$BAM2chooseCmp
    }

    #find the corresponding ROI objects from their names, as usual
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos1=match(input$ROI1chooseCmp,nomi)
    roi1=uniqueROI(ROIvariables$listROI[[pos1]])
    pos2=match(input$ROI2chooseCmp,nomi)
    roi2=uniqueROI(ROIvariables$listROI[[pos2]])
    #get Range and find the overlaps, with minimum bp set by the user
    range1=getRange(roi1)
    range2=getRange(roi2)
    ov=findOverlaps(range1,range2,minoverlap=input$minOverlapCmp)
    toplot$cmp$ov=ov
    sh=subjectHits(ov)
    qh=queryHits(ov)
    #find common ranges for later. These are AREAS of overlap.
    #this is for putting a single value in Venn diagram:
    #     ----    ----
    #    -----------
    # is 1 area of overlap
    common=union( unique(range1[qh]) ,unique(range2[sh]))
    totalunion=union(range1,range2)
    if (length(sh)==0&length(qh)==0){
      #if no overlap at all, ranges are only "exclusive"
      range1_ov=NULL
      range2_ov=NULL
      range1_notov=range1
      range2_notov=range2
      perc1=0
      perc2=0
      jaccard= 0
    }else{
      #otherwise, calculate how many overlaps. Ranges in 1 that overlaps with 2,
      #ranges of 2 that overlaps with 1, ranges exclusively in 1 and excl. in 2
      #and find numbers and fractions.
      range1_ov=unique(range1[qh])
      range2_ov=unique(range2[sh])
      range1_notov=range1[-qh]
      range2_notov=range2[-sh]    
      perc1=  round( (length(range1_ov)/length(range1))*100 )
      perc2=  round( (length(range2_ov)/length(range2))*100 )
      jaccard=round(  length(common)  /length(totalunion),2)
    }


    #make the matrix for barplot and barplot of the overlap (the most correct representation)
    matbar=as.matrix(data.frame(n1=c(length(range1_ov),length(range1_notov),0,0),
                      n2=c(0,0,length(range2_ov),length(range2_notov))))
    colnames(matbar)=paste(c("ROI 1","ROI 2")," (",c(length(range1_ov),length(range2_ov)),"; ",c(perc1,perc2),"%)",sep="")
    
    toplot$cmp$Colors=strsplit(input$chooseColorPaletteCmp,split="_")[[1]]
    toplot$cmp$matbar=matbar
    toplot$cmp$jaccard=jaccard
    toplot$cmp$n1=n1
    toplot$cmp$n2=n2


    output$viewBarplotCmp<-renderPlot({
      Colors=strsplit(input$chooseColorPaletteCmp,split="_")[[1]]
      par(xpd=TRUE,mar=c(14,6,5,5))
      barplot(matbar,col=c(Colors[3],Colors[1],Colors[4],Colors[2]),main=paste("Jaccard idx:",jaccard),cex.main=1.2,las=2)
      title(ylab = "Interval number", cex.lab = 1,line = 4)
      legend("bottomleft",inset=c(0,-1.5),legend=c(paste(n1,"alone"),
                                                      paste(n1,"overlapping"),
                                                      paste(n2,"overlapping"),
                                                      paste(n2,"alone")),pch=19,
              col=c(Colors[1],Colors[3],Colors[4],Colors[2]))
    })

    #slow function...
    #for Venn Diagram, overlap number is number of overlapping AREAS (see above)
    #Venn diagram is an approximation of the overlap representation.
    #Good when ranges have similar properties or similar distribution of widths
    #Bad when globally they have different widths
    area1=length(range1_notov)
    area2=length(range2_notov)
    toplot$cmp$area1=area1
    toplot$cmp$area2=area2
    toplot$cmp$common=common
    output$viewVennCmp<-renderPlot({
      Colors=strsplit(input$chooseColorPaletteCmp,split="_")[[1]]
      griglia<-draw.pairwise.venn(area1=area1+length(common), area2=area2+length(common), 
                cross.area=length(common),c(n1,n2),cex=1,ext.dist=c(.01,-0.04),
                ext.line.lwd=0,ext.length=0,ext.pos=90,alpha=c(0.1,0.1),col=c(Colors[1],Colors[2]),fill=c(Colors[1],Colors[2]),
                label.col=c(Colors[1],Colors[3],Colors[4]),cat.pos=c(310,50),cat.dist=c(0.1,0.1),lwd=c(2,2),cat.cex=c(1.2,1.2),
                cat.col=c(Colors[1],Colors[4]),margin=0.3,main="Intervals overlap",main.cex=1.2,main.col="black",main.pos=c(1,1),filename=NULL)

    })
    #download button for PDF Barplot of overlap
    output$saveviewBarplotCmp=renderUI({downloadButton('saveviewBarplotCmpbutton', 'Get PDF')})
    #download button for Venn
    output$saveviewVennCmp=renderUI({downloadButton('saveviewVennCmpbutton', 'Get PDF')})
    print(paste("Pairwise comparison between",input$ROI1chooseCmp,"and",input$ROI2chooseCmp))


    #now, check BAMlists of the 2 ROIs selected
    getbam1=names(getBAMlist(roi1))  
    getbam2=names(getBAMlist(roi2))   

    toplot$cmp$getbam1=getbam1
    toplot$cmp$getbam2=getbam2
    toplot$cmp$BAM1chooseCmp=BAM1_choose
    toplot$cmp$BAM2chooseCmp=BAM2_choose

    #if there is at least one overlap and at least one BAMfile has been selectd by the user
    #(meaning that is associated with at least one of the 2 ROIs selected), proceed
    if ( (nchar(BAM1_choose)>0 |nchar(BAM2_choose)>0)  & (!is.null(getbam1) | !is.null(getbam2)) & length(ov)>0){
      ov_range1=findOverlaps(common,range1,minoverlap=input$minOverlapCmp)
      ov_range2=findOverlaps(common,range2,minoverlap=input$minOverlapCmp)

      if(input$chooseNormalizationCmp=="readdensity"){
        lab_axis="Read density (rpm/bp)"
      }else{
        lab_axis="Normalized reads (rpm)"
      }
      #if BAM1 is present and selected, proceed
      if(nchar(BAM1_choose)>0){
        pos_bam1_1=match(BAM1_choose,getbam1)
        bam1_1=getBAMlist(roi1)[[pos_bam1_1]]
        bam1_1=unlist(lapply(bam1_1,sum))
        if(!is.null(bam1_1)){
          range1only_bam1=bam1_1[-qh]
          if (input$chooseNormalizationCmp=="readdensity"){
            range1only_bam1=range1only_bam1/width(range1_notov)
          }
          #extract BAM1 of range1, in positions of range1 that overlaps with range2 
          bam1_common=bam1_1[subjectHits(ov_range1)]
          if (input$chooseNormalizationCmp=="readdensity"){
            width_range1_common=width(range1[subjectHits(ov_range1)])
            bam1_common=bam1_common/width_range1_common
            #here, we can have that multiple range1 overlap with a single range 2 in "common" regions.
            #in this case, for each overlapping AREA, find the MEAN of signals of BAM1 of ranges1
            #that overlap with this single range2:
            #range1: ------ ---------           ---
            #range2:  ----------------------------
            #the average of signals of BAM1 calculated on range1 is the signal that
            #will appear as "common regions"
            common_bam1=tapply(bam1_common,queryHits(ov_range1),mean)
          }else{
            common_bam1=tapply(bam1_common,queryHits(ov_range1),mean)
          }

          #if BAM1 is associated also with range2, we can have the reads of BAM1
          #but in the exclusive regions of ROI2!!
          pos_bam1_2=match(BAM1_choose,getbam2)
          if(!is.na(pos_bam1_2) & length(pos_bam1_2)>0){
            bam1_2=getBAMlist(roi2)[[pos_bam1_2]]
            bam1_2=unlist(lapply(bam1_2,sum)) 
            range2only_bam1=bam1_2[-sh]       
            if(input$chooseNormalizationCmp=="readdensity") {
              range2only_bam1=range2only_bam1/width(range2_notov)
            }
          }else{
            #bam2 for range2, but not bam1 for range2
            bam1_2=NA
            range2only_bam1=NA
          }          
        }else{
          common_bam1=NA
          range1only_bam1=NA
          range2only_bam1=NA          
        }
        
      }else{
        common_bam1=NA
        range1only_bam1=NA
        range2only_bam1=NA
      }

      #if BAM2 is present and selected, proceed
      if(nchar(BAM2_choose)>0){
        pos_bam2_2=match(BAM2_choose,getbam2)
        bam2_2=getBAMlist(roi2)[[pos_bam2_2]]
        bam2_2=unlist(lapply(bam2_2,sum))

        if(!is.null(bam2_2)){
          range2only_bam2=bam2_2[-sh]
          if(input$chooseNormalizationCmp=="readdensity"){
            range2only_bam2=range2only_bam2/width(range2_notov)
          }
          bam2_common=bam2_2[subjectHits(ov_range2)]
          if (input$chooseNormalizationCmp=="readdensity"){
            width_range2_common=width(range2[subjectHits(ov_range2)])
            bam2_common=bam2_common/width_range2_common
            common_bam2=tapply(bam2_common,queryHits(ov_range2),mean)
          }else{
            common_bam2=tapply(bam2_common,queryHits(ov_range2),mean)
          }
          #if BAM2 is associated also with ROI 1,
          pos_bam2_1=match(BAM2_choose,getbam1)
          if(!is.na(pos_bam2_1) & length(pos_bam2_1)>0){
            bam2_1=getBAMlist(roi1)[[pos_bam2_1]]
            bam2_1=unlist(lapply(bam2_1,sum))
            range1only_bam2=bam2_1[-qh]
            if(input$chooseNormalizationCmp=="readdensity"){
              range1only_bam2=range1only_bam2/width(range1_notov)
            }
          } else{
            #bam1 for range1, but not bam2 for range1
            bam2_1=NA
            range1only_bam2=NA
          }          
        }else{
          range2only_bam2=NA
          range1only_bam2=NA
          common_bam2=NA          
        }
         
      }else{
        #no bam2 for range2. => no bam2 for range2 and range1
        range2only_bam2=NA
        range1only_bam2=NA
        common_bam2=NA
      }

      toplot$cmp$range2only_bam1=range2only_bam1
      toplot$cmp$range1only_bam1=range1only_bam1
      toplot$cmp$common_bam1=common_bam1
      toplot$cmp$range1only_bam2=range1only_bam2
      toplot$cmp$range2only_bam2=range2only_bam2
      toplot$cmp$common_bam2=common_bam2
      toplot$cmp$lab_axis=lab_axis

      toplot$cmp$BAM1chooseCmp=BAM1_choose
      toplot$cmp$BAM2chooseCmp=BAM2_choose

      #make the boxplot with data available. If both bam files associated with
      #both selected ROIs, the boxplot will be complete
      #as common regions, the signal of BAM1 is the BAM1 signal of range1 inside
      #common areas (the average for each single overlapping area)
      #the same for BAM2: is the average of BAM2 signals inside ranges2 in common areas
      output$viewBoxplotCmp<-renderPlot({
        Colors=strsplit(input$chooseColorPaletteCmp,split="_")[[1]]
        par(xpd=TRUE,mar=c(16,4,2,2))
        if(input$islogCmp){
          arraytoplot=list(log2(range2only_bam1),log2(range1only_bam1),log2(common_bam1),
                log2(range1only_bam2),log2(range2only_bam2),log2(common_bam2))
          laby=paste("log2",lab_axis)
        }else{
          arraytoplot=list(range2only_bam1,range1only_bam1,common_bam1,
                range1only_bam2,range2only_bam2,common_bam2)
          laby=lab_axis
        }

        suppressWarnings(boxplot(arraytoplot,col=c(Colors[2],Colors[1],Colors[3],Colors[1],Colors[2],Colors[4]),
                outline=FALSE,notch=TRUE,varwidth=TRUE,at =c(1,2,3, 5,6,7),xaxt="n",ylab=laby))

        axis(1,at=c(2,6),labels=c("enrich. 1","enrich. 2"),las=2,cex=1)
        legend("bottomleft",inset=c(0,-1),legend=c(paste(n1,"alone"),
                                                      paste(n1,"overlapping"),
                                                      paste(n2,"overlapping"),
                                                      paste(n2,"alone")),pch=19,
              col=c(Colors[1],Colors[3],Colors[4],Colors[2]))
        #another legend for complete name of enrichments
        legend("bottomleft",inset=c(0,-1.5),legend=c(paste("enrich. 1: ",BAM1_choose,sep=""),paste("enrich. 2: ",BAM2_choose,sep="")),box.col = "black",bg = "white")
      })
      #pdf button for boxplot
      output$saveviewBoxplotCmp=renderUI({downloadButton('saveviewBoxplotCmpbutton', 'Get PDF')})
      #now the download data button will appear
      output$saveboxdataCmp=renderUI({downloadButton('saveenrichmentBoxCmpdata', 'Save data')})



      quantiles=seq(0,0.99,0.05) #will be x coordinates for the final plot
      #here, at least bam1_1 or bam2_2 must exist by definition

      if(exists("bam1_1")& !exists("bam2_2")){
        #only BAM1 of ROI1
        fracts=lapply(quantiles,test_ov,ov=ov,bam1_1=bam1_1,bam2_2=NULL)
      }
      if(exists("bam2_2")& !exists("bam1_1")){
        #only BAM2 of ROI2
        fracts=lapply(quantiles,test_ov,ov=ov,bam1_1=NULL,bam2_2=bam2_2)
      }
      if(exists("bam1_1")& exists("bam2_2")){
        #both BAM1 of ROI1 and BAM2 of ROI2
        fracts=lapply(quantiles,test_ov,ov=ov,bam1_1=bam1_1,bam2_2=bam2_2)
      }

      #extract values for the calibration overlap plot:
      only1=unlist(lapply(fracts,"[[",1))
      only2=unlist(lapply(fracts,"[[",2))
      combined1=unlist(lapply(fracts,"[[",3))
      combined2=unlist(lapply(fracts,"[[",4))

      toplot$cmp$quantiles=quantiles
      toplot$cmp$only1=only1
      toplot$cmp$only2=only2
      toplot$cmp$combined1=combined1
      toplot$cmp$combined2=combined2

      #plot the calibration of overlap
      output$viewCalibrationCmp<-renderPlot({
        Colors=strsplit(input$chooseColorPaletteCmp,split="_")[[1]]
        par(xpd=TRUE,mar=c(13,5,2,2))
        plot(1, type="n", xlab="Enrichment threshold (quantile)", ylab="fraction of overlap", xlim=c(0, 1), ylim=c(0, 1),xaxt="n")
        lines(quantiles,only1,type="p",col=Colors[1])
        lines(quantiles,only2,type="p",col=Colors[2])
        lines(quantiles,combined1,type="p",col=Colors[3])
        lines(quantiles,combined2,type="p",col=Colors[4])
        if(exists("bam1_1")& !exists("bam2_2")){
          legend("bottomleft",inset=c(0,-0.8),legend=paste(n1,"overlap"),col=Colors[1],bg="transparent")
        }else if(exists("bam2_2")& !exists("bam1_1")){
          legend("bottomleft",inset=c(0,-0.8),legend=paste(n2,"overlap"),col=Colors[2],bg="transparent")
        }else if(exists("bam1_1")& exists("bam2_2")){
          legend("bottomleft",inset=c(0,-0.8),
                  legend=c(paste(n1,"overlap"),paste(n2,"overlap"),paste(n1,"overlap combined"),paste(n2,"overlap combined")),
                      col=c(Colors[1],Colors[2],Colors[3],Colors[4]),pch=19,bg="transparent" )
        }
        axis(1,quantiles)
      })
      #button for PDF of the calibration overlap plot
      output$saveviewCalibrationCmp=renderUI({downloadButton('saveviewCalibrationCmpbutton', 'Get PDF')})


      #common_ regions are 2 slightly different subsets:
      #range1:       -------     -----        -----
      #range2:          ------------------------

      #common_bam1:  -------     -----        ----- (avg. BAM1 signal)
      #common_bam2:     ------------------------    (avg. BAM2 signal)
      #for each AREA of overlap, common_bam1 is average signal of 1 calculated on 
      #range1 regions, while common_bam2 is the avg. signal of 2 calculated on range2 regions.
      #=> the number of common_bam1 == number common_bam2 == overlapping AREAS
      #=> they can be treated as ONE subset of regions; this number is equal 
      #to the number written in the intersection of the Venn diagram


      #scatterplot only if both ranges have both bam files
      #if both bam files were associated to all the subsets (3) => 2 BAM files for each subset (region1 only, rgion2 only, cmmon regions)
      if(all(!is.na(range2only_bam1)) & all(!is.na(range1only_bam1)) & all(!is.na(range2only_bam2)) & all(!is.na(range1only_bam2)) & all(!is.na(common_bam2)) & all(!is.na(common_bam1)) ){
         
        df=as.data.frame(matrix( rep(NA,(length(range2only_bam2)+length(range1only_bam1)+length(common_bam1))*2),ncol=2 ))
        labs=c(rep("2only",length(range2only_bam2)), rep("1only",length(range1only_bam1)),rep("common",length(common_bam1)) )
        df[,3]=labs
        colnames(df)=c("bam1","bam2","label")

        #fill df with enrichment values with the proper belonging subset
        if(length(range2only_bam1)>0){
          df[1:length(range2only_bam2),1]=range2only_bam1
          df[1:length(range2only_bam2),2]=range2only_bam2          
        }
        if(length(range1only_bam1)>0){
          df[(length(range2only_bam2)+1):(length(range2only_bam2)+length(range1only_bam2)),1]=range1only_bam1
          df[(length(range2only_bam2)+1):(length(range2only_bam2)+length(range1only_bam2)),2]=range1only_bam2          
        }
        if(length(common_bam1)>0){
          df[(nrow(df)-length(common_bam1)+1):nrow(df),1]=common_bam1
          df[(nrow(df)-length(common_bam2)+1):nrow(df),2]=common_bam2          
        }


        #remove 0, negative value and NAs for logarithm
        arraylogical= df[,1]<=0 | df[,2]<=0 | is.na(df[,1]==0) | is.na(df[,2]==0) 
        dfprovv=df[!arraylogical,]

        cor_all=round(cor(log2(dfprovv[,1]),log2(dfprovv[,2])),2)
        cor_common=round(cor(log2(dfprovv[dfprovv[,3]=="common",1]),log2(dfprovv[dfprovv[,3]=="common",2])),2)

        toplot$cmp$completedf=df
        #random sample df if more than 10.000
        if (nrow(df)>5000){
          set.seed(123) # set the seed, so the result will be the same every run with same parameters
          df=df[sort(sample(nrow(df), 5000,replace=FALSE)), ]
        }

        toplot$cmp$df=df
        toplot$cmp$cor_common=cor_common
        toplot$cmp$cor_all=cor_all
        

        #here try to put the checkbox input with the choice of what to show in the scatterplot
        listscatterchoice=list("exclusive ROI-1 ranges"="1only","exclusive ROI-2 ranges"="2only","common ranges"="common")
        output$showScatterChoice<-renderUI({checkboxGroupInput("scatterplotChoice",label="In scatterplot, show only:",choices=listscatterchoice,selected=unlist(listscatterchoice) )})

        output$viewScatterplotCmp<-renderPlot({
          Colors=strsplit(input$chooseColorPaletteCmp,split="_")[[1]]

          #here, based on what you want to show, remove from df the rows "1only","2only" or "common"
          subsetToShow=input$scatterplotChoice

          if(input$islogCmp){
            toplot1=log2(df[,1])
            toplot2=log2(df[,2])
            #set -inf values to min
            toplot1[is.infinite(toplot1)]=min(toplot1[!is.infinite(toplot1)])
            toplot2[is.infinite(toplot2)]=min(toplot2[!is.infinite(toplot2)])
            xlims=c(min(toplot1),max(toplot1))
            ylims=c(min(toplot2),max(toplot2))
            laby=paste("log2",toplot$cmp$BAM2chooseCmp,lab_axis)
            labx=paste("log2",toplot$cmp$BAM1chooseCmp,lab_axis)
          }else{
            toplot1=df[,1]
            toplot2=df[,2]
            xlims=c(0,max(toplot1))
            ylims=c(0,max(toplot2))
            laby=paste(toplot$cmp$BAM2chooseCmp,lab_axis)
            labx=paste(toplot$cmp$BAM1chooseCmp,lab_axis)
          }
          #here, define a priori the x and ylim:

          pos=df[,3] %in% subsetToShow
          df2=df[pos,]
          toplot1=toplot1[pos]
          toplot2=toplot2[pos]

          par(xpd=TRUE,mar=c(15,5,2,2))
          plot(toplot1,toplot2,col=ifelse(df2[,3]=="1only",Colors[1],ifelse(df2[,3]=="2only",Colors[2],"black")), 
                            xlab=labx,ylab=laby,xlim=xlims,ylim=ylims)
          legend("bottomleft",inset=c(0,-0.75),legend=c(paste("cor common intervals: ",cor_common,sep=""),paste("cor all intervals: ",cor_all,sep="")),box.col = "black",bg = "white")
          legend("bottomleft",inset=c(0,-1.2),legend=c(paste(n2,"intervals"),paste(n1,"intervals"),"common"),box.col = "black",bg = "white",pch=19,col=c(Colors[2],Colors[1],"black"))
        })
        #button for PDF of the scatterplot
        output$saveviewScatterplotCmp=renderUI({downloadButton('saveviewScatterplotCmpbutton', 'Get PDF')})
        #button to download data about the scatterplot
        output$saveScatterdataCmp=renderUI({downloadButton('saveenrichmentScatterCmpdata', 'Save data')})

      }else{
        output$viewScatterplotCmp<-renderPlot({NULL})
        output$showScatterChoice<-renderUI({NULL})
        output$saveScatterdataCmp=renderUI({NULL})
        output$saveviewScatterplotCmp=renderUI({NULL})
      }
        
    }else{
      #plot enrichment stuff (box & scatter) are NULL
      output$viewCalibrationCmp<-renderPlot({NULL})
      output$viewBoxplotCmp<-renderPlot({NULL})
      output$viewScatterplotCmp<-renderPlot({NULL})
      output$showScatterChoice<-renderUI({NULL})
      output$saveboxdataCmp=renderUI({NULL})
      output$saveviewBoxplotCmp=renderUI({NULL})
      output$saveScatterdataCmp=renderUI({NULL})
      output$saveviewScatterplotCmp=renderUI({NULL})
    }
  }else{
    #all plots are NULL
    output$viewCalibrationCmp<-renderPlot({NULL})
    output$viewVennCmp<-renderPlot({NULL})
    output$viewBarplotCmp<-renderPlot({NULL})
    output$saveviewBarplotCmp=renderUI({NULL})
    output$saveviewVennCmp=renderUI({NULL})
    output$viewBoxplotCmp<-renderPlot({NULL})
    output$saveviewBoxplotCmp=renderUI({NULL})
    output$saveboxdataCmp=renderUI({NULL})
    output$viewScatterplotCmp<-renderPlot({NULL})
    output$showScatterChoice<-renderUI({NULL})
    output$showScatterChoice<-renderUI({NULL})
    output$saveScatterdataCmp=renderUI({NULL})
    output$saveviewScatterplotCmp=renderUI({NULL})
    sendSweetAlert(
      session = session,
      title = "ROIs not selected",
      text = "For this analysis you need select two ROIs (ROI1 and ROI2)",
      type = "error"
    )
  }
},ignoreInit=TRUE)





#observer to get PDF button for barplot cmp
output$saveviewBarplotCmpbutton<- downloadHandler(
  filename=function() {
      paste('Barplot_overlap.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    par(xpd=TRUE,mar=c(16,6,5,5))
    barplot(toplot$cmp$matbar,col=c(toplot$cmp$Colors[3],toplot$cmp$Colors[1],toplot$cmp$Colors[4],toplot$cmp$Colors[2]),main=paste("Jaccard idx:",toplot$cmp$jaccard),cex.main=1.2,las=2)
    title(ylab = "Interval number", cex.lab = 1,line = 4)
    legend("bottomleft",inset=c(0,-1),legend=c(paste(toplot$cmp$n1,"alone"),
                                                      paste(toplot$cmp$n1,"overlapping"),
                                                      paste(toplot$cmp$n2,"overlapping"),
                                                      paste(toplot$cmp$n2,"alone")),pch=19,
            col=c(toplot$cmp$Colors[1],toplot$cmp$Colors[3],toplot$cmp$Colors[4],toplot$cmp$Colors[2]))
    dev.off()
  } 
)





#observer to get PDF button for venn cmp
output$saveviewVennCmpbutton<- downloadHandler(
  filename=function() {
      paste('Venn_overlap.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    griglia<-draw.pairwise.venn(area1=toplot$cmp$area1+length(toplot$cmp$common), area2=toplot$cmp$area2+length(toplot$cmp$common), 
                cross.area=length(toplot$cmp$common),c(toplot$cmp$n1,toplot$cmp$n2),cex=1,ext.dist=c(.01,-0.04),
                ext.line.lwd=0,ext.length=0,ext.pos=90,alpha=c(0.1,0.1),col=c(toplot$cmp$Colors[1],toplot$cmp$Colors[2]),fill=c(toplot$cmp$Colors[1],toplot$cmp$Colors[2]),
                label.col=c(toplot$cmp$Colors[1],toplot$cmp$Colors[3],toplot$cmp$Colors[4]),cat.pos=c(310,50),cat.dist=c(0.1,0.1),lwd=c(2,2),cat.cex=c(1.2,1.2),
                cat.col=c(toplot$cmp$Colors[1],toplot$cmp$Colors[4]),margin=0.3,main="Intervals overlap",main.cex=1.2,main.col="black",main.pos=c(1,1),filename=NULL)

    dev.off()
  } 
)



#observer to get PDF button for boxplot cmp
output$saveviewBoxplotCmpbutton<- downloadHandler(
  filename=function() {
      paste('Boxplot_enrichment_overlap.pdf', sep='')
  },
  content=function(file) {
        if(input$islogCmp){
          arraytoplot=list(log2(toplot$cmp$range2only_bam1),log2(toplot$cmp$range1only_bam1),log2(toplot$cmp$common_bam1),
                log2(toplot$cmp$range1only_bam2),log2(toplot$cmp$range2only_bam2),log2(toplot$cmp$common_bam2))
          laby=paste("log2",toplot$cmp$lab_axis)
        }else{
          arraytoplot=list(toplot$cmp$range2only_bam1,toplot$cmp$range1only_bam1,toplot$cmp$common_bam1,
                toplot$cmp$range1only_bam2,toplot$cmp$range2only_bam2,toplot$cmp$common_bam2)
          laby=toplot$cmp$lab_axis
        }
        pdf(file)
        par(xpd=TRUE,mar=c(16,4,2,2))

        suppressWarnings(boxplot(arraytoplot,col=c(toplot$cmp$Colors[2],toplot$cmp$Colors[1],toplot$cmp$Colors[3],
                                                                  toplot$cmp$Colors[1],toplot$cmp$Colors[2],toplot$cmp$Colors[4]),
                outline=FALSE,notch=TRUE,varwidth=TRUE,at =c(1,2,3, 5,6,7),xaxt="n",ylab=laby))

        axis(1,at=c(2,6),labels=c("enrich. 1","enrich. 2"),las=2,cex=1)
        legend("bottomleft",inset=c(0,-0.6),legend=c(paste(toplot$cmp$n1,"alone"),
                                                      paste(toplot$cmp$n1,"overlapping"),
                                                      paste(toplot$cmp$n2,"overlapping"),
                                                      paste(toplot$cmp$n2,"alone")),pch=19,
              col=c(toplot$cmp$Colors[1],toplot$cmp$Colors[3],toplot$cmp$Colors[4],toplot$cmp$Colors[2]))
        legend("bottomleft",inset=c(0,-0.9),legend=c(paste("enrich. 1: ",toplot$cmp$BAM1chooseCmp,sep=""),paste("enrich. 2: ",toplot$cmp$BAM2chooseCmp,sep="")),box.col = "black",bg = "white")
        dev.off()
  } 
)




#observer to get PDF button for overlap calibration plot
output$saveviewCalibrationCmpbutton<- downloadHandler(
  filename=function() {
      paste('Plot_calibration_overlap.pdf', sep='')
  },
  content=function(file) {
        
        pdf(file)
        par(xpd=TRUE,mar=c(13,4,2,2))
        plot(1, type="n", xlab="Enrichment threshold (quantile)", ylab="fraction of overlap", xlim=c(0, 1), ylim=c(0, 1),xaxt="n")
        lines(toplot$cmp$quantiles,toplot$cmp$only1,type="p",col=toplot$cmp$Colors[1])
        lines(toplot$cmp$quantiles,toplot$cmp$only2,type="p",col=toplot$cmp$Colors[2])
        lines(toplot$cmp$quantiles,toplot$cmp$combined1,type="p",col=toplot$cmp$Colors[3])
        lines(toplot$cmp$quantiles,toplot$cmp$combined2,type="p",col=toplot$cmp$Colors[4])
        #print legnd according to the existence of bam1, bam2 or both
        if(!is.null(toplot$cmp$only1)& is.null(toplot$cmp$only2)){
          legend("bottomleft",inset=c(0,-0.6),legend=paste(toplot$cmp$n1,"overlap"),col=toplot$cmp$Colors[1],bg="transparent")
        }else if(is.null(toplot$cmp$only1)& !is.null(toplot$cmp$only2)){
          legend("bottomleft",inset=c(0,-0.6),legend=paste(toplot$cmp$n2,"overlap"),col=toplot$cmp$Colors[2],bg="transparent")
        }else if(!is.null(toplot$cmp$only1)& !is.null(toplot$cmp$only2)){
          legend("bottomleft",inset=c(0,-0.6),
                  legend=c(paste(toplot$cmp$n1,"overlap"),paste(toplot$cmp$n2,"overlap"),paste(toplot$cmp$n1,"overlap combined"),paste(toplot$cmp$n2,"overlap combined")),
                      col=c(toplot$cmp$Colors[1],toplot$cmp$Colors[2],toplot$cmp$Colors[3],toplot$cmp$Colors[4]),pch=19,bg="transparent" )
        }
        axis(1,toplot$cmp$quantiles)

        dev.off()
  } 
)


#observer for download data of boxplot
output$saveenrichmentBoxCmpdata<- downloadHandler(
  filename=function() {
      paste('boxplot_data.xls', sep='')
  },
  content=function(file) {
      if(input$islogCmp){
        arraytoplot=list(log2(toplot$cmp$range2only_bam1),log2(toplot$cmp$range1only_bam1),log2(toplot$cmp$common_bam1),
                log2(toplot$cmp$range1only_bam2),log2(toplot$cmp$range2only_bam2),log2(toplot$cmp$common_bam2))
      }else{
        arraytoplot=list(toplot$cmp$range2only_bam1,toplot$cmp$range1only_bam1,toplot$cmp$common_bam1,
                toplot$cmp$range1only_bam2,toplot$cmp$range2only_bam2,toplot$cmp$common_bam2)
      }
      maxval=max(length(arraytoplot[[1]]),
                  length(arraytoplot[[2]]),
                  length(arraytoplot[[3]]),
                  length(arraytoplot[[4]]),
                  length(arraytoplot[[5]]),
                  length(arraytoplot[[6]])     )
      arr=matrix(rep("",maxval*6),ncol=6)
      arr[1:length(arraytoplot[[1]]),1]=arraytoplot[[1]]
      arr[1:length(arraytoplot[[2]]),2]=arraytoplot[[2]]
      arr[1:length(arraytoplot[[3]]),3]=arraytoplot[[3]]
      arr[1:length(arraytoplot[[4]]),4]=arraytoplot[[4]]
      arr[1:length(arraytoplot[[5]]),5]=arraytoplot[[5]]
      arr[1:length(arraytoplot[[6]]),6]=arraytoplot[[6]]        

      colnames(arr)=c(paste(toplot$cmp$n2,"_intervals_only_",toplot$cmp$BAM1chooseCmp,"_enrichment",sep=""),
                      paste(toplot$cmp$n1,"_intervals_only_",toplot$cmp$BAM1chooseCmp,"_enrichment",sep=""),
                      paste("common_intervals_",toplot$cmp$BAM1chooseCmp,"_enrichment",sep=""),
                      paste(toplot$cmp$n1,"_intervals_only_",toplot$cmp$BAM2chooseCmp,"_enrichment",sep=""),
                      paste(toplot$cmp$n2,"_intervals_only_",toplot$cmp$BAM2chooseCmp,"_enrichment",sep=""),
                      paste("common_intervals_",toplot$cmp$BAM2chooseCmp,"_enrichment",sep=""))
      colnames(arr)=gsub(" ","_",colnames(arr))
      write.table(arr,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  } 
  
)




#observer to get PDF button for scatterplot cmp
output$saveviewScatterplotCmpbutton<- downloadHandler(
  filename=function() {
      paste('Scatterplot_overlap.pdf', sep='')
  },
  content=function(file) {

          df=toplot$cmp$df
          subsetToShow=input$scatterplotChoice

          if(input$islogCmp){
            toplot1=log2(df[,1])
            toplot2=log2(df[,2])
            #set -inf values to min
            toplot1[is.infinite(toplot1)]=min(toplot1[!is.infinite(toplot1)])
            toplot2[is.infinite(toplot2)]=min(toplot2[!is.infinite(toplot2)])
            xlims=c(min(toplot1),max(toplot1))
            ylims=c(min(toplot2),max(toplot2))
            laby=paste("log2",toplot$cmp$BAM2chooseCmp,toplot$cmp$lab_axis)
            labx=paste("log2",toplot$cmp$BAM1chooseCmp,toplot$cmp$lab_axis)
          }else{
            toplot1=df[,1]
            toplot2=df[,2]
            xlims=c(0,max(toplot1))
            ylims=c(0,max(toplot2))
            laby=paste(toplot$cmp$BAM2chooseCmp,toplot$cmp$lab_axis)
            labx=paste(toplot$cmp$BAM1chooseCmp,toplot$cmp$lab_axis)
          }
          #here, define a priori the x and ylim:

          pos=df[,3] %in% subsetToShow
          df2=df[pos,]
          toplot1=toplot1[pos]
          toplot2=toplot2[pos]


        pdf(file)
        par(xpd=TRUE,mar=c(15,5,2,2))
        plot(toplot1,toplot2,col=ifelse(df2[,3]=="1only",toplot$cmp$Colors[1],ifelse(df2[,3]=="2only",toplot$cmp$Colors[2],"black")), xlab=labx,ylab=laby,
                            xlim=xlims,ylim=ylims)
        legend("bottomleft",inset=c(0,-0.5),legend=c(paste("cor common intervals: ",toplot$cmp$cor_common,sep=""),paste("cor all intervals: ",toplot$cmp$cor_all,sep="")),box.col = "black",bg = "white")
        legend("bottomleft",inset=c(0,-0.8),legend=c(paste(toplot$cmp$n2,"intervals"),paste(toplot$cmp$n1,"intervals"),"common"),box.col = "black",bg = "white",pch=19,col=c(toplot$cmp$Colors[2],toplot$cmp$Colors[1],"black"))
        dev.off()
  } 
)



#respond to download data about scatterplot if present:
output$saveenrichmentScatterCmpdata<- downloadHandler(
  filename=function() {
      paste('scatter_data.xls', sep='')
  },
  content=function(file) {

      if(input$islogCmp){
        toplot$cmp$completedf[,1]=log2(toplot$cmp$completedf[,1])
        toplot$cmp$completedf[,2]=log2(toplot$cmp$completedf[,2])
      }

      # laby=paste(toplot$cmp$BAM2chooseCmp,toplot$cmp$lab_axis)
      # labx=paste(toplot$cmp$BAM1chooseCmp,toplot$cmp$lab_axis)
      df=toplot$cmp$completedf
      df[df[,3]=="2only",][,3]=toplot$cmp$n2
      df[df[,3]=="1only",][,3]=toplot$cmp$n1

      colnames(df)=c(toplot$cmp$BAM1chooseCmp,toplot$cmp$BAM2chooseCmp,"label")
      colnames(df)=gsub(" ","_",colnames(df))
      write.table(df,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  } 
  
)














##################################################################
##################################################################
##################################################################
# respond to plot Digital Heatmap
##################################################################
##################################################################
##################################################################



###react to help buttons:
#parameters
observeEvent(input$msg_digitalHeatmap_parameters, {
  boxHelpServer(msg_digitalHeatmap_parameters)
})
#digital heatmap
observeEvent(input$msg_digitalHeatmap_heatmap, {
  boxHelpServer(msg_digitalHeatmap_heatmap)
})
#jaccard idx 
observeEvent(input$msg_digitalHeatmap_jaccardIdx, {
  boxHelpServer(msg_digitalHeatmap_jaccardIdx)
})
#overlap frequency bias
observeEvent(input$msg_digitalHeatmap_overlapBias, {
  boxHelpServer(msg_digitalHeatmap_overlapBias)
})






toListenDigitalHeat <- reactive({
    list(input$confirmUpdateDigitalHeat1,input$confirmUpdateDigitalHeat2)
})

observeEvent(toListenDigitalHeat(),{
  set.seed(123)



  #check if at least one ROI in ROIsForDigitalHeat is selected
  if(length(ROIvariables$listROI)>0 & length(input$ROIsForDigitalHeat)>0 & 
              length(input$ROImaster)>0 & input$binsDigitalHeat>0 & 
              input$sampleRandomDigitalHeat>0 & isvalid(input$binsDigitalHeat) &
              isvalid(input$sampleRandomDigitalHeat)) {

    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROImaster,nomi)
    #all master rois selected
    roi=ROIvariables$listROI[pos]
    maxbinstohave=min(unlist(lapply(roi,checkMaxBins)))
    
    toplot$digital$maxbinstohave=maxbinstohave
    toplot$digital$ROIsForDigitalHeat=input$ROIsForDigitalHeat
    strandSpecificity=input$StrandSpecOverlap
    

    ###########################################################################################
    #here, give new ordering for the ROIs to show
    temp_names=toplot$digital$ROIsForDigitalHeat
    allnumbers=as.character(1:length(temp_names))
    listprovv=list()
    for (i in 1:length(temp_names)){
      stringval=grep(paste("reorderROIdigitalHeat",i,"$",sep=""),names(input),value=TRUE)
      listprovv[[i]]=input[[ stringval ]]
    }  
    listprovv=as.numeric(unlist(listprovv))
    #check if new order comprises all the possible positions
    if (identical(unique(sort(listprovv)),unique(sort(as.numeric(allnumbers)))) ){
      toplot$digital$ROIsForDigitalHeat=temp_names[order(listprovv)]
    }else{
      sendSweetAlert(
        session = session,
        title = "Ordering problem",
        text = "You didn't put all the ranking positions correct in the new ordering",
        type = "error"
      )   
      return() 
    }    
    ###########################################################################################


    if(input$binsDigitalHeat<=maxbinstohave){
      nameROI=unlist(lapply(roi,getName))
      nbin=input$binsDigitalHeat
      bigrange=list()
      bigbamlist=list()
      fixes=list()
      ranges_forclick=list()
      #loop through all master ROIs selected
      for(i in 1:length(roi)){
        unifROI=uniqueROI(roi[[i]])
        rangeroi=getRange(unifROI)
        bigrange[[i]]=granges(rangeroi)
        bigrange[[i]]$label=nameROI[i]

        ranges_forclick[[i]]=rangeroi
        bigbamlist[[i]]=getBAMlist(unifROI)
        fixes[[i]]=granges(getFixed(unifROI))

      }

      finalrange=Reduce(c,bigrange)

      #reduce number if sample>max length
      if(input$sampleRandomDigitalHeat>length(finalrange)){
        samplerandom=length(finalrange)
      }else{
        samplerandom=input$sampleRandomDigitalHeat
      }
      pos2=match(toplot$digital$ROIsForDigitalHeat,nomi)
      roiTooverlap=ROIvariables$listROI[pos2]
      roiTooverlap_range_list=lapply(roiTooverlap,getRange)
      numbersamples=length(roiTooverlap_range_list)
      names_samples=toplot$digital$ROIsForDigitalHeat

      #sample the matrix according to samplerandom:
      Digital_sample_pos=sort(sample(1:length(finalrange),samplerandom,replace=FALSE))
      
      #save the subselected range in the temporary variable
      finalrange_subsample=finalrange[Digital_sample_pos]
      
      #retrieve the labels of ranges (which ROI belongs to) after subsampling
      labelstosplit=finalrange[Digital_sample_pos]$label
      
      
      toplot$digital$ROImaster=input$ROImaster
      toplot$digital$binsDigitalHeat=input$binsDigitalHeat
      toplot$digital$sampleRandomDigitalHeat=input$sampleRandomDigitalHeat

      finalrange_splitted=split(finalrange,finalrange$label)

      if(length(finalrange_splitted)==1){
        #find positive/negative positions to flip matrix rows of overlapping bin

        #if only one label. Then, compute the overlap for everything, and then sample
        matlistTotal=as.list(rep(NA,length(finalrange_splitted)))
        
        #this loop make no sense
        for(i in 1:length(finalrange_splitted)){
          pos_negative=as.character(strand(finalrange_splitted[[i]]))=="-"
          pos_positive=as.character(strand(finalrange_splitted[[i]]))!="-"
          if(nc==1){
            listSingleROI=lapply(1:length(roiTooverlap),function(k) {
                  mat=countOverlapsInBins(finalrange_splitted[[i]],roiTooverlap_range_list[[k]],nbins=nbin,strandspecific=strandSpecificity)
                  if(any(pos_negative)){
                    if (table(pos_negative)["TRUE"]>1){
                      if(nbin==1){
                        #if nbin==1, no need to invert according to strand...               
                      }else{
                        mat[pos_negative,]<- mat[pos_negative,][ , ncol(mat[pos_negative,]):1]
                      }                      
                    }
                  }
                  return(mat)
                })            
          }else{

            decision=0
            tryCatch({
              listSingleROI=mclapply(1:length(roiTooverlap),function(k) {
                  mat=suppressWarnings(countOverlapsInBins(finalrange_splitted[[i]],roiTooverlap_range_list[[k]],nbins=nbin,strandspecific=strandSpecificity))
                  if(any(pos_negative)){
                    if (table(pos_negative)["TRUE"]>1){
                      if(nbin==1){
                        #if nbin==1, no need to invert according to strand...               
                      }else{
                        mat[pos_negative,]<- mat[pos_negative,][ , ncol(mat[pos_negative,]):1]
                      }                      
                    }
                  }
                  return(mat)
                },mc.cores=nc)  
              decision=1  
            },
            warning = function( w ){
              print("Warning: using single core for memory limit in overlaps single master ROI...")
            },
            error = function( err ){
              print("Warning: using single core for memory limit in overlaps single master ROI...")
            })
            if(decision==0){
              listSingleROI=lapply(1:length(roiTooverlap),function(k) {
                  mat=countOverlapsInBins(finalrange_splitted[[i]],roiTooverlap_range_list[[k]],nbins=nbin,strandspecific=strandSpecificity)
                  if(any(pos_negative)){
                    if (table(pos_negative)["TRUE"]>1){
                      if(nbin==1){
                        #if nbin==1, no need to invert according to strand...               
                      }else{
                        mat[pos_negative,]<- mat[pos_negative,][ , ncol(mat[pos_negative,]):1]
                      }                      
                    }
                  }
                  return(mat)
              })   
            }
          
          }

          names(listSingleROI)=toplot$digital$ROIsForDigitalHeat
          #for each roi master selected, loop countOverlapsInBins with all input$
          matlistTotal[[i]]=listSingleROI
        }

        matlistDigital=matlistTotal
        for(i in 1:length(matlistTotal)){
          for(k in 1:length(matlistTotal[[i]])){
            if(nrow(matlistTotal[[i]][[k]])==1){
              matlistDigital[[i]][[k]]=matlistTotal[[i]][[k]]
            }else{
              matlistDigital[[i]][[k]]=matlistTotal[[i]][[k]][Digital_sample_pos,,drop=FALSE]
            }
          }
        }

      }else{
        matlistTotal=NULL
        #more labels, no venn diagrams, multicore and split at the beginning
        finalrange_sample_join=finalrange[Digital_sample_pos]
        #xtract only the rows in the Digital_sample_pos random positions
        finalrange_sample=split(finalrange_sample_join,factor(labelstosplit,levels=unique(labelstosplit)))
        #create matrixes of countOverlapsInBins
        matlistDigital=as.list(rep(NA,length(finalrange_sample)))
        for(i in 1:length(finalrange_sample)){
          pos_negative=as.character(strand(finalrange_sample[[i]]))=="-"
          pos_positive=as.character(strand(finalrange_sample[[i]]))!="-"

          if(nc==1){
            listSingleROI=lapply(1:length(roiTooverlap),function(k) {
              mat=countOverlapsInBins(finalrange_sample[[i]],roiTooverlap_range_list[[k]],nbins=nbin,strandspecific=strandSpecificity)
              if(any(pos_negative)){
                if (table(pos_negative)["TRUE"]>1){
                  if(nbin==1){
                    #if nbin==1, no need to invert according to strand...               
                  }else{
                    mat[pos_negative,]<- mat[pos_negative,][ , ncol(mat[pos_negative,]):1]
                  }                      
                }
              }
              return(mat)
            })
          }else{
            decision=0
            tryCatch({
              listSingleROI=mclapply(1:length(roiTooverlap),function(k) {
                mat=suppressWarnings(countOverlapsInBins(finalrange_sample[[i]],roiTooverlap_range_list[[k]],nbins=nbin,strandspecific=strandSpecificity))
                if(any(pos_negative)){
                  if (table(pos_negative)["TRUE"]>1){
                    if(nbin==1){
                      #if nbin==1, no need to invert according to strand...               
                    }else{
                      mat[pos_negative,]<- mat[pos_negative,][ , ncol(mat[pos_negative,]):1]
                    }                      
                  }
                }
                return(mat)
             },mc.cores=nc) 
              decision=1
            },
            warning = function( w ){
              print("Warning: using single core for memory limit in overlaps...")
            },
            error = function( err ){
              print("Warning: using single core for memory limit in overlaps...")
            })
            if(decision==0){
              listSingleROI=lapply(1:length(roiTooverlap),function(k) {
                mat=countOverlapsInBins(finalrange_sample[[i]],roiTooverlap_range_list[[k]],nbins=nbin,strandspecific=strandSpecificity)
                if(any(pos_negative)){
                  if (table(pos_negative)["TRUE"]>1){
                    if(nbin==1){
                      #if nbin==1, no need to invert according to strand...               
                    }else{
                      mat[pos_negative,]<- mat[pos_negative,][ , ncol(mat[pos_negative,]):1]
                    }                      
                  }
                }
                return(mat)
              })
            }

            
         
          }

          names(listSingleROI)=toplot$digital$ROIsForDigitalHeat
          #for each roi master selected, loop countOverlapsInBins with all input$
          matlistDigital[[i]]=listSingleROI
        }
      }


      names(matlistDigital)=input$ROImaster

      ###prepare for clustering
      clusterTypeDigitalHeat=input$clusterTypeDigitalHeat
      distmethodDigitalHeat=input$distmethodDigitalHeat
      clustmethodDigitalHeat=input$clustmethodDigitalHeat
      #cluster and collapse second order matrixes. Should be parallel
      pos_cluster=match(input$ROIforClusteringDigitalHeat,toplot$digital$ROIsForDigitalHeat)
      zeros=c()
      for(i in 1:length(matlistDigital)){
        toinspect=matlistDigital[[i]]
        zeros=c(zeros,nrow(toinspect[[1]]))
      }

      #if clustering position is not found (no clustering) => no clustering
      #need "hierarchical" to have a clustering method (in this case, doesn't mean anything)
      #also dist/clust method are fake characters
      if(length(pos_cluster)==0 |  (any(zeros==1))){
        pos_cluster=NULL
        clusterTypeDigitalHeat="hierarchical"
        distmethodDigitalHeat="fake"
        clustmethodDigitalHeat="fake"
      }

      #if parameters of the clustering (numeric inputs) are not set (NA),
      #or not valid, in this case, no clustering. 
      #1- for any clustering, the number must be set.
      #2- If kmeans is without valid
      #advanced numeric parameters, set hclust with fake (no) clustering (see above...)

      clustnumDigitalHeat=input$clustnumDigitalHeat
      clustrandomstartsDigitalHeat=input$clustrandomstartsDigitalHeat
      clustnumiterationsDigitalHeat=input$clustnumiterationsDigitalHeat
      if(!isvalid(clustnumDigitalHeat)){
        pos_cluster=NULL
        clusterTypeDigitalHeat="hierarchical"
        distmethodDigitalHeat="fake"
        clustmethodDigitalHeat="fake"        
      }
      if(clusterTypeDigitalHeat=="kmean"){
        if(!isvalid(clustrandomstartsDigitalHeat) | !isvalid(clustnumiterationsDigitalHeat)){
          pos_cluster=NULL
          clusterTypeDigitalHeat="hierarchical"
          distmethodDigitalHeat="fake"
          clustmethodDigitalHeat="fake"           
        }
      }


      ##CUSTERING
      if(nc==1){
        finalmats_complete=lapply(1:length(matlistDigital),function(i) {
          matlist=matlistDigital[[i]]
          #be sure that matlist is a list of matrixes (if nbin==1, can be transformed in "numeric")
          matlist=lapply(matlist,as.matrix)

          if (clusterTypeDigitalHeat=="hierarchical"){
            joinmat=clusterMatrix(matlist=matlist,distmethod=distmethodDigitalHeat,clustmethod=clustmethodDigitalHeat,clustinds=pos_cluster)
          
          #if K-means selected            
          }else{
            startingpoints=clustrandomstartsDigitalHeat
            iter=clustnumiterationsDigitalHeat
            if (is.null(startingpoints)){
              startingpoints=1
            }else{
              if(startingpoints<=0){
                startingpoints=1
              }
              if(startingpoints>40){
                startingpoints=40
              }
            }
            #iterations for K-means clustering
            if(is.null(iter)){
              iter=1
            }else{
              if(iter<=0){
                iter=1
              }
              if(iter>40){
                iter=40
              }                
            }
            joinmat=clusterMatrixKmeans(matlist=matlist,clustinds=pos_cluster,numberclusters=clustnumDigitalHeat,startingpoints=startingpoints,iter=iter)
          }
          return(joinmat)
        })

      }else{
        decision=0
        tryCatch({
          finalmats_complete=mclapply(1:length(matlistDigital),function(i) {
            matlist=matlistDigital[[i]]
            #be sure that matlist is a list of matrixes (if nbin==1, can be transformed in "numeric")
            matlist=lapply(matlist,as.matrix)
            
            if (clusterTypeDigitalHeat=="hierarchical"){
              joinmat=clusterMatrix(matlist=matlist,distmethod=distmethodDigitalHeat,clustmethod=clustmethodDigitalHeat,clustinds=pos_cluster)
            
            #if K-means selected            
            }else{
              startingpoints=clustrandomstartsDigitalHeat
              iter=clustnumiterationsDigitalHeat
              if (is.null(startingpoints)){
                startingpoints=1
              }else{
                if(startingpoints<=0){
                  startingpoints=1
                }
                if(startingpoints>40){
                  startingpoints=40
                }
              }
              #iterations for K-means clustering
              if(is.null(iter)){
                iter=1
              }else{
                if(iter<=0){
                  iter=1
                }
                if(iter>40){
                  iter=40
                }                
              }
              joinmat=clusterMatrixKmeans(matlist=matlist,clustinds=pos_cluster,numberclusters=clustnumDigitalHeat,startingpoints=startingpoints,iter=iter)
            }
            return(joinmat)
          },mc.cores=nc)
          decision=1
        },
        warning = function( w ){
          print("Warning: using single core for memory limit in clustering...")
        },
        error = function( err ){
          print("Warning: using single core for memory limit in clustering...")
        })
        if(decision==0){
          finalmats_complete=lapply(1:length(matlistDigital),function(i) {
            matlist=matlistDigital[[i]]
            #be sure that matlist is a list of matrixes (if nbin==1, can be transformed in "numeric")
            matlist=lapply(matlist,as.matrix)
            
            if (clusterTypeDigitalHeat=="hierarchical"){
              joinmat=clusterMatrix(matlist=matlist,distmethod=distmethodDigitalHeat,clustmethod=clustmethodDigitalHeat,clustinds=pos_cluster)
            
            #if K-means selected            
            }else{
              startingpoints=clustrandomstartsDigitalHeat
              iter=clustnumiterationsDigitalHeat
              if (is.null(startingpoints)){
                startingpoints=1
              }else{
                if(startingpoints<=0){
                  startingpoints=1
                }
                if(startingpoints>40){
                  startingpoints=40
                }
              }
              #iterations for K-means clustering
              if(is.null(iter)){
                iter=1
              }else{
                if(iter<=0){
                  iter=1
                }
                if(iter>40){
                  iter=40
                }                
              }
              joinmat=clusterMatrixKmeans(matlist=matlist,clustinds=pos_cluster,numberclusters=clustnumDigitalHeat,startingpoints=startingpoints,iter=iter)
            }
            return(joinmat)
          })          
        }
      }
      #finalmats: extract matrix from finalmats_complete
      finalmats=lapply(finalmats_complete,function(element) {element$mat})

      names(finalmats)=nameROI
      matrixes_processed=do.call(rbind,finalmats)

      #if nROI==1, paste to names_samples the number and % of the regions overlapping
      #binding fraction CALCULATION
      if(length(finalmats)==1){
        for(i in 1:(length(roiTooverlap_range_list)) ){
          #take the piece of the matrix, columns corresponding to one overlap
          partialmat=matrixes_processed[,(((i-1)*nbin)+1):((i*nbin)) ]
          if (nbin>1){
            #calculate rows that have at least one 1
            logic=apply(partialmat,1,sum)>0
          }else{
            #number of bins ==1 => partialmat is a simple array
            logic=partialmat>0
          }

          overlap=table(logic)["TRUE"]
          if(is.na(overlap)){
            overlap=0
          }
          perc=round(overlap/length(logic)*100)
          names_samples[i]=paste(names_samples[i]," (",overlap,"/",length(logic),";",perc,"%)",sep="")            

        }
      }



      toplot$digital$numbersamples=numbersamples
      #toplot$digital$colorsDigitalHeat=input$colorsDigitalHeat
      toplot$digital$matrixes_processed=matrixes_processed
      toplot$digital$roiTooverlap_range_list=roiTooverlap_range_list
      toplot$digital$labelstosplit=labelstosplit
      toplot$digital$nbin=nbin
      toplot$digital$names_samples=names_samples
      toplot$digital$matlistTotal=matlistTotal
      toplot$digital$clustnumDigitalHeat= clustnumDigitalHeat
      toplot$digital$clustering=finalmats_complete[[1]]
      toplot$digital$clusterTypeDigitalHeat=clusterTypeDigitalHeat



      #reset the color menu:
      if(input$optioncolorsforDigitalHeat=="custom"){
        bams=toplot$digital$ROIsForDigitalHeat
        lista=list()

        if(length(bams)>0){
          output$showcolorsDigitalheat<-renderUI({
            for (i in 1:length(bams)){
              lista[[i]]=fluidRow(
                                  column(12,
                                        colorSelectorInput(inputId=paste("colorCustomDigitalHeat",i,sep=""),label=bams[i],choices=ColsArray,
                                                      selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                  )
                          )

            }
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
              lista
            )
          }) 
                
        }else{
          #bams are NULL. show nothing
          output$showcolorsDigitalheat<-renderUI({NULL})
        }
      }



      #clear button for the new roi: we have another analysis
      output$newROIfromDigitalHeat_out<-renderUI({NULL})
      output$textselectedelementsDigitalHeat<-renderText({NULL})

      #displays the number of random elements used for the heatmap. If all pssible selected, green
      output$textfractionelementsDigitalHeat<-renderText({
        if(samplerandom==length(finalrange)){
          color="'green'"
        }else{
          color="'red'"
        }
        paste("Intervals displayed: <font color=",color,">",samplerandom,"/",length(finalrange),"</font>",sep="")
      })

      output$heatmapDigital<-renderPlot({

        nbin=toplot$digital$nbin
        #find xlim according to how many numbersamples if numbersamples<10
        if(numbersamples<10){
          factormult=10/numbersamples
        }else{
          factormult=1
        }
        
        matProc=toplot$digital$matrixes_processed
        #color settings
        palette_col=toplot$digital$colorpalettes
        
        if(isolate(input$optioncolorsforDigitalHeat)=="global"){

          colorsplitted=strsplit(palette_col,split="_")[[1]]
          palette=colorRampPalette(colorsplitted)(n=2)

          
        }else{
          if(length(palette_col)==1){
            if(palette_col=="white_red4"){
              colorsplitted=strsplit(palette_col,split="_")[[1]]
              palette=colorRampPalette(colorsplitted)(n=2)  
            }else{
              palette=c()
              for(i in 1:length(palette_col)){
                palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n=2))
              }
            }
                      
          }else{
            #if input$optioncolorsforDigitalHeat is custom, palette should have length==nbams
            palette=c()
            for(i in 1:length(palette_col)){
              palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n=2))
            }  
            #if different colors, we need to change the order of magnitude 
            #for each bam

            matProc[matProc==1]=0.999
            matProc[matProc==0]=0.001
            for(i in 1:length(toplot$digital$ROIsForDigitalHeat)){
                piecematrix=matProc[, ((i-1)*nbin+1):(i*nbin),drop=FALSE ]
                asvec=unique(as.vector(piecematrix))
                if(length(asvec)==1){


                  #all elments are the same, it can break the color palette, change one
                  if(asvec==0.999){
                    piecematrix[1,1]=0.001
                  }else if(asvec==0.001){
                    piecematrix[1,1]=0.999
                  }
                }
                piecematrix=piecematrix+(i-1)
                matProc[, ((i-1)*nbin+1):(i*nbin) ]=piecematrix
                #palettes[i]=colorRampPalette(colorsplitted)(n=brk)
            }           
          }

        }



        trasp=t(matProc)
        if(max(dim(trasp))>30000){
          raster=FALSE
        }else{
          raster=TRUE
        }
        par(mar = rep(0, 4))
        image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1],axes = FALSE, xlab = "", ylab = "",col=palette
        #REMOVE x,y lim if cordinates don't match
              ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp)),useRaster=raster  
              )

        counter=seq(nbin,nbin*length(roiTooverlap_range_list),nbin)

        for (csep in counter){

          rect(xleft = csep , ybottom = rep(0,length(csep)), 
            xright = csep  + 0.01, 
            ytop = rep(ncol(trasp), csep), lty = 1, lwd = 1, 
            col = "black", border = "black")
        }


        #how many for each ROI? then revert, because the matrix is inverted in "image"
        #but the matrix original is correct, with ROI order and columns=bins; rows=ranges
        correctOrder=table(labelstosplit)[order(order(unique(labelstosplit)))]
        counter2=cumsum(rev(as.integer(correctOrder)))

        for (csep2 in counter2){
          rect(xleft = rep(0,length(csep2))  , ybottom = csep2 , 
            xright = rep(nrow(trasp) , csep2), 
            ytop = csep2 +0.01, lty = 1, lwd = 1, 
            col = "black", border = "black")
        }
      })

      output$saveheatmapDigital=renderUI({downloadButton('saveheatmapDigitalbutton', 'Get PDF')})
      
      print("Drawn digital heatmap")
      #here put code for colored cluster column at the left of the digital heatmap
      #conditions: must be "cluster" (not "ranked") and must have number ROI ==1
      #otherwise plot NULL. Less bins, more balanced clusters.
      set.seed(123)
      color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      
      toplot$digital$clusternumbermatrix=NULL
      if(length(input$ROImaster) ==1 & length(input$ROIforClusteringDigitalHeat)>0 ){
        color_distinct_cluster=sample(color, toplot$digital$clustnumDigitalHeat)
        #for this line of code, thanks to
        # http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
        # Jelena-bioinf / Megatron
        

        if(toplot$digital$clustnumDigitalHeat<=433 & toplot$digital$clustnumDigitalHeat>0){        
          if(!is.null(toplot$digital$clustering$clustobj)){

            #if custering will be drown, it means that we have 1 master ROI and clustering is ok.
            #therefore, keep info for subsequent subselection based on the click on cluster            
            #remember that this ROI must be re-annotated
            #1-define the starting "material", obtained at the beginning of the heat button:
            range_tokeep=ranges_forclick[[1]]
            bams_tokeep=bigbamlist[[1]]
            fix_tokeep=fixes[[1]]
            #2-subsample them (1 single ROI)
            range_tokeep=range_tokeep[Digital_sample_pos]
            bams_tokeep=lapply(bams_tokeep,function(bamob){bamob[Digital_sample_pos]})
            fix_tokeep=fix_tokeep[Digital_sample_pos]
            #3-re-ordering based on clustering (toplot$digital$clustering$clustobj$ord)
            range_tokeep=range_tokeep[toplot$digital$clustering$ord]
            bams_tokeep=lapply(bams_tokeep,function(bamob){bamob[toplot$digital$clustering$ord]})
            fix_tokeep=fix_tokeep[toplot$digital$clustering$ord]
            #4-save in temporary variables to be used by clicking the clsuter of digital
            toplot$digital$range_tokeep=range_tokeep
            toplot$digital$bams_tokeep=bams_tokeep
            toplot$digital$fix_tokeep=fix_tokeep


            if (clusterTypeDigitalHeat=="hierarchical"){
              # use  toplot$digital$clustering from clusterMatrix function
              # if number of masterROI==1, we have only one cluster object
              clusterarray=cutree( toplot$digital$clustering$clustobj,k=toplot$digital$clustnumDigitalHeat)
            }else{
              #in this case (kmeans) clusterarray is already defined by  toplot$digital$clustering$cluster
              clusterarray= toplot$digital$clustering$clustobject$cluster
            }

            #plot image of clusters 
            clusterarray=clusterarray[toplot$digital$clustering$ord]
            clusterarray=rev(clusterarray)
            matrix_like=matrix(clusterarray,nrow=1)
            toplot$digital$clusternumbermatrix=matrix_like
            toplot$digital$color_distinct_cluster=color_distinct_cluster

            output$clustersImageLeftDigital<-renderPlot({
              par(mar = c(0,0,0,0))
              image(0:nrow(matrix_like), 0:ncol(matrix_like),matrix_like,col=color_distinct_cluster,axes = FALSE, xlab = "", ylab = "",ylim=c(0,ncol(matrix_like)))               
            })

          }else{
            output$clustersImageLeftDigital<-renderPlot({NULL})
            toplot$digital$color_distinct_cluster=NULL
            toplot$digital$clusternumbermatrix=NULL
            toplot$digital$range_tokeep=NULL
            toplot$digital$bams_tokeep=NULL
            toplot$digital$fix_tokeep=NULL
          }
        }else{
          output$clustersImageLeftDigital<-renderPlot({NULL})
          toplot$digital$color_distinct_cluster=NULL
          toplot$digital$clusternumbermatrix=NULL
          toplot$digital$range_tokeep=NULL
          toplot$digital$bams_tokeep=NULL
          toplot$digital$fix_tokeep=NULL
        }

      }else{
        output$clustersImageLeftDigital<-renderPlot({NULL})
        toplot$digital$clusternumbermatrix=NULL
        toplot$digital$color_distinct_cluster=NULL
        toplot$digital$range_tokeep=NULL
        toplot$digital$bams_tokeep=NULL
        toplot$digital$fix_tokeep=NULL
      }



      #put names of the ROIs for the columns
      output$textNameDigitalHeat <- renderPlot({
        par(mar = rep(0, 4))
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i")
        #have to start from x coordinate of the half of first cell to half of the last
        
        if(numbersamples<10){
          #divide by 10 and multiply by number sample
          newmax=(1/10)*numbersamples
          halfcellwidthmin=(newmax/numbersamples)/2
          halfcellwidthmax=newmax-halfcellwidthmin

        }else{
          halfcellwidthmin=(1/numbersamples)/2
          halfcellwidthmax=1-halfcellwidthmin
        }
        cextext=0.8
        text(x=seq(halfcellwidthmin, halfcellwidthmax ,length.out=numbersamples),labels=names_samples,y=rep(0.5,numbersamples),srt=90,cex=cextext  )
              
      })


      #if nbins>1, do the plot of frequency to see if overlaps are concentrated
      #in paarticular places. This can be done even when matlistTotal is null,
      #so we have multiple master ROIs
      if(nbin>1){
        #set the colors (different for each single block).
        #take them from the original definition of the entire palette of 433 colors  
        #and save them temporary, to be consistent with the downloaded PDF of the plot
        ## calculate the bias of the position of the overlap (if more than 1 bin)
        ## matlistDigital is the list of the list of the matrix (first order: master ROI;
        ## second order: ROI to view). It is NOT ranked based on the clustering.
        ## For any reason, you have to rank each of these matrixes based on the 
        ## finalmats_complete[[n]]$ord order, for each master ROI
        ## scheme:
        ## 1 0 0
        ## 1 0 0
        ## 1 0 1
        ## 0 0 1
        ##
        ## 3 0 2 <- result for each ROI (matrix). For this specific ROI, 
        ## overlaps are more frequent at the left of the master ROI
        ## keep also the length of each subset, for the conversion frequency/percentage
        freq_overlap=as.list(rep(NA,length(matlistDigital) ))
        names(freq_overlap)=names(matlistDigital)
        length_subsets=c()
        for(i in 1:length(freq_overlap)){
          length_subsets[i]=nrow(matlistDigital[[i]][[1]])
          freq_overlap[[i]]=lapply(matlistDigital[[i]],function(k){apply(k,2,sum)})
        }
        names(length_subsets)=names(matlistDigital)


        
        toplot$digital$freq_overlap=freq_overlap
        toplot$digital$length_subsets=length_subsets
        color_distinct_bias=sample(color, length(freq_overlap)*length(freq_overlap[[1]]))
        toplot$digital$color_distinct_bias=color_distinct_bias

        output$frequencyOvDigitalHeat <- renderPlot({
          #interactively change freq/% if button is pressed:
          ovtoplot=freq_overlap
          if(input$FracToPercDigitalHeat){
            for(i in 1:length(freq_overlap)){
              ovtoplot[[i]]=lapply(freq_overlap[[i]],function(k){k/length_subsets[[i]]})
            }
            maxval=1
            ylabfreq="Fraction of overlaps"
          }else{
            maxval=max(length_subsets)
            ylabfreq="Number of overlaps"
          }
          #now plot (either % or absolute numbers)
          plot(0,type='n',xaxt="n",ylim=c(0,maxval),xlim=c(1,nbin),ylab=ylabfreq,xlab="Genomic Window",
                        main="Overlaps bias",cex.main=1.2,cex=1,xaxt="n",cex.lab=1.2,cex.axis=1)
          count=1
          legnames=c()
          for(i in 1:length(ovtoplot)){
            for(k in 1:length(ovtoplot[[i]])){
              #draw the lines/points
              lines(ovtoplot[[i]][[k]],col=color_distinct_bias[count])
              legnames[count]=paste(names(ovtoplot[[i]])[k],"in",names(ovtoplot)[i])
              count=count+1
            }
          }
          axis(1,at=seq(1,nbin,(nbin-1)/2),
                      labels=c( "start","midpoint","end" ),cex.axis=1)
          #legend for each line
          legend("topleft",legend=legnames,lty=rep(1,length(legnames)),col=color_distinct_bias,bg="transparent")
        })
        output$savefrequencyOvDigitalHeat=renderUI({downloadButton('savefrequencyOvDigitalHeatbutton', 'Get PDF')})
      }else{
        output$frequencyOvDigitalHeat <- renderPlot({NULL})
        output$savefrequencyOvDigitalHeat=renderUI({NULL})
      }



      #These plots refer to ALL ranges considered, not the random sample.
      #modify thresholds here if you want for example matrix even with 2 ROIs
      if (!is.null(matlistTotal)){
        matList=matlistTotal[[1]]

        if(length(matList)>1){
          #for each matrix of 0/1, find max (0 or 1 if overlap or not)
          mat=lapply(matList,function(i) {apply(i,1,max)})
          mat=do.call(cbind,mat)
          colnames(mat)=toplot$digital$ROIsForDigitalHeat


          #CHW-RUSKEY diagram: de-comment if you want it
          # if(length(matList)==400 | length(matList)==500){
          #   #plot chow-ruskey plot with Vennerable package

          #   new=table(apply(mat,1,paste0,collapse=""))
          #   venn2=Venn(Weight=new,SetNames=toplot$digital$ROIsForDigitalHeat)
          #   array=venn2@IndicatorWeight[,ncol(venn2@IndicatorWeight)]
          #   array[is.na(array)]=0
          #   venn2@IndicatorWeight[,ncol(venn2@IndicatorWeight)]=array

          #   toplot$digital$jaccardorvenn=venn2

          #   output$JaccardDigital<-renderPlot({
          #     plot(venn2,type="ChowRuskey")
          #   })
            
          # }

          # matfinal=matrix(rep(NA,(ncol(mat)*ncol(mat))),ncol=ncol(mat))
          matperc=matrix(rep(NA,(ncol(mat)*ncol(mat))),ncol=ncol(mat))
          #plot heatmap with overlap
          for(i in 1:ncol(mat)){
            for(k in i:ncol(mat)){
              localsum=nrow(mat)
              ov=sum(mat[,i] & mat[,k])
              #with currentsum and ov, we can calculate the % of co-binding in the ROI
              #then fisher test in 2 ways: global, using 250000 DNAse accessibile sites
              #(Dolfini et al., 2016) or the number of sites in the ROI
              globalsum=250000
              AnotB=sum(mat[,i]==1 & mat[,k]==0)
              BnotA=sum(mat[,k]==1 & mat[,i]==0)
              # #if(input$chooseOrderingDigitalHeat=="global"){
              #   matfisher=matrix(c(ov,BnotA,AnotB,globalsum-(ov+AnotB+BnotA)),byrow=TRUE,ncol=2)
              # #}else{
              #   #matfisher=matrix(c(ov,BnotA,AnotB,localsum-(ov+AnotB+BnotA)),byrow=TRUE,ncol=2)
              # #}

              # fish=fisher.test(matfisher)
              
              # if(fish$estimate>1){
              #   #significant overlap
              #   val= -log10(fish$p.value)
              # }else{
              #   #significant mutual exclusion
              #   val= log10(fish$p.value)
              # }
              # matfinal[i,k]=val
              matperc[i,k]=round(ov/(AnotB+BnotA+ov),2)

            }
          }

          # colnames(matfinal)=rownames(matfinal)=toplot$digital$ROIsForDigitalHeat
          colnames(matperc)=rownames(matperc)=toplot$digital$ROIsForDigitalHeat

          # matfinal[lower.tri(matfinal)] <- t(matfinal)[lower.tri(matfinal)]
          matperc[lower.tri(matperc)] <- t(matperc)[lower.tri(matperc)]


          # matfinal[matfinal>20]=20
          # matfinal[matfinal<(-20)]=-20

          matperc[is.na(matperc)]=0
          matperc[is.infinite(matperc)]=1


          toplot$digital$jaccardorvenn=matperc
          output$JaccardDigital<-renderPlot({
            absmax=max(abs(matperc))
            trasp=t(matperc)
            brk=c( seq( 0 , 1,0.01 ))
            my_palette <- colorRampPalette(c("white","blue"))(n = length(brk)-1 ) 
            strlengths=nchar(colnames( trasp ))
            if(max(strlengths)<=25){
              margin_12=13
            }else{
              margin_12=13+ ((max(strlengths)-25) /3)
            }
            if(margin_12>16){
              margin_12=16
            }

            par(mar=c(margin_12,margin_12,2,3),xpd=TRUE)
            image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk,main="Jaccard index")
            axis( 2, at=seq(0.5,ncol(trasp)+0.5-1,1 ), labels= rev(colnames( trasp )), las= 2 )
            axis( 1, at=seq(0.5,ncol(trasp)+0.5-1,1 ), labels= colnames( trasp ), las= 2 )
            for (x in (nrow(matperc)-1+0.5):0.5  )
              for (y in 0.5: ((ncol(matperc)-1+0.5)   ))
                text(y,x, round(matperc[ncol(matperc)-x+0.5,y+0.5],1),col="red")
          })


          output$saveJaccardDigital=renderUI({downloadButton('saveJaccardDigitalbutton', 'Get PDF')})

        }else{
          output$JaccardDigital<-renderPlot({NULL})
          output$saveJaccardDigital=renderUI({NULL})
          toplot$digital$matlistTotal=NULL
          #plot NULL, because only one ROI considered
        }
      }else{
        #plot NULL
        output$JaccardDigital<-renderPlot({NULL})
        output$saveJaccardDigital=renderUI({NULL})
        toplot$digital$matlistTotal=NULL
      }

    }else{
      #too many bins selected. Must reduce
      toplot$digital$color_distinct_cluster=NULL
      toplot$digital$ROIsForDigitalHeat=NULL
      toplot$digital$ROImaster=NULL
      toplot$digital$binsDigitalHeat=NULL
      toplot$digital$sampleRandomDigitalHeat=NULL
      output$textselectedelementsDigitalHeat<-renderText({NULL})
      output$heatmapDigital<-renderPlot({NULL})
      output$saveheatmapDigital=renderUI({NULL})
      output$textNameDigitalHeat<-renderPlot({NULL}) 
      output$JaccardDigital<-renderPlot({NULL})  
      output$saveJaccardDigital=renderUI({NULL})
      toplot$digital$matlistTotal=NULL 
      toplot$digital$clusternumbermatrix=NULL
      output$clustersImageLeftDigital<-renderPlot({NULL})
      toplot$digital$range_tokeep=NULL
      toplot$digital$bams_tokeep=NULL
      toplot$digital$fix_tokeep=NULL
      output$newROIfromDigitalHeat_out<-renderUI({NULL})
      output$textfractionelementsDigitalHeat<-renderText({NULL})
      sendSweetAlert(
        session = session,
        title = "Too many bins",
        text = "The selected number of bins exceeds the maximum allowed, because some ranges have length < nbins. Try with a lower number",
        type = "error"
      )
    }

  }else{
    toplot$digital$ROIsForDigitalHeat=NULL
    toplot$digital$color_distinct_cluster=NULL
    toplot$digital$ROImaster=NULL
    toplot$digital$binsDigitalHeat=NULL
    toplot$digital$sampleRandomDigitalHeat=NULL
    output$textselectedelementsDigitalHeat<-renderText({NULL})
    output$heatmapDigital<-renderPlot({NULL})
    output$saveheatmapDigital=renderUI({NULL})
    output$textNameDigitalHeat<-renderPlot({NULL})
    output$JaccardDigital<-renderPlot({NULL})
    output$saveJaccardDigital=renderUI({NULL})
    toplot$digital$matlistTotal=NULL
    toplot$digital$clusternumbermatrix=NULL
    output$clustersImageLeftDigital<-renderPlot({NULL})
    toplot$digital$range_tokeep=NULL
    toplot$digital$bams_tokeep=NULL
    toplot$digital$fix_tokeep=NULL
    output$newROIfromDigitalHeat_out<-renderUI({NULL})
    output$textfractionelementsDigitalHeat<-renderText({NULL})
    sendSweetAlert(
      session = session,
      title = "Wrong/missing parameters",
      text = "Some parameters are missing or wrong; check them and then click the button",
      type = "error"
    )
  }
},ignoreInit=TRUE)




#observer for button click for pdf save of digital heatmap
output$saveheatmapDigitalbutton<- downloadHandler(
  filename=function() {
      paste('Heatmap_digital.pdf', sep='')
  },
  content=function(file) {
      
      #find xlim according to how many numbersamples if numbersamples<8
      if(toplot$digital$numbersamples<8){
        factormult=10/toplot$digital$numbersamples
      }else{
        factormult=1
      }

      matrixes_processed=toplot$digital$matrixes_processed
      nbin=toplot$digital$nbin

      #modify the colors if required: colors are 2
      palette_col=toplot$digital$colorpalettes

      

      if(isolate(input$optioncolorsforDigitalHeat)=="global"){

        colorsplitted=strsplit(palette_col,split="_")[[1]]
        palette=colorRampPalette(colorsplitted)(n=2)

        
      }else{
        if(length(palette_col)==1){
          if(palette_col=="white_red4"){
            colorsplitted=strsplit(palette_col,split="_")[[1]]
            palette=colorRampPalette(colorsplitted)(n=2)  
          }else{
            palette=c()
            for(i in 1:length(palette_col)){
              palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n=2))
            }
          }
                    
        }else{
          #if input$optioncolorsforDigitalHeat is custom, palette should have length==nbams
          palette=c()
          for(i in 1:length(palette_col)){
            palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n=2))
          }  
          #if different colors, we need to change the order of magnitude 
          #for each bam

          matrixes_processed[matrixes_processed==1]=0.999
          matrixes_processed[matrixes_processed==0]=0.001
          for(i in 1:length(toplot$digital$ROIsForDigitalHeat)){
              piecematrix=matrixes_processed[, ((i-1)*nbin+1):(i*nbin),drop=FALSE ]
              asvec=unique(as.vector(piecematrix))
              if(length(asvec)==1){
                #all elments are the same, it can break the color palette, change one
                if(asvec==0.999){
                  piecematrix[1,1]=0.001
                }else if(asvec==0.001){
                  piecematrix[1,1]=0.999
                }
              }
              piecematrix=piecematrix+(i-1)
              matrixes_processed[, ((i-1)*nbin+1):(i*nbin) ]=piecematrix
              #palettes[i]=colorRampPalette(colorsplitted)(n=brk)
          }           
        }

      }



      intervals_vertical=table(toplot$digital$labelstosplit)[order(order(unique(toplot$digital$labelstosplit)))]
      #revert order because the "0" starts at the bottom, (last ROI)
      intervals_vertical=rev(intervals_vertical)
      cumulative=cumsum(intervals_vertical)
      tosubtract=intervals_vertical/2
      pos_vertical=cumulative-tosubtract

      #calculate % out of the total length of ROIs considered...
      percs_rois=round({intervals_vertical/sum(intervals_vertical)}*100,2)
      labels_vert=paste("(",unname(intervals_vertical),"; ",percs_rois,"%)",sep="")

      trasp=t(matrixes_processed)
      pdf(file)

      #if cluster image left present, plot it
      mlay=matrix(c(1,2),ncol=2,nrow=1)
      if(!is.null(toplot$digital$clusternumbermatrix)){
        layout(mlay,widths=c(40,60),heights=c(100,100))
      }else{
        #layout(mlay,widths=c(0,100),heights=c(100,100))
      }

      # par(mar = c(16,12,1,1))
      # image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
      # #REMOVE x,y lim if cordinates don't match
      #       ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp))   
      #       )

      if(!is.null(toplot$digital$clusternumbermatrix)){
        par(mar = c(12,12,0,0))
        image(0:nrow(toplot$digital$clusternumbermatrix), 0:ncol(toplot$digital$clusternumbermatrix),toplot$digital$clusternumbermatrix,col=toplot$digital$color_distinct_cluster,axes = FALSE, xlab = "", ylab = "",ylim=c(0,ncol(trasp)))
        axis( 2,at=pos_vertical,labels=paste(rev(toplot$digital$ROImaster),"\n",labels_vert),las=1,tick=FALSE)
        axis(1,0.5,labels="cluster",las= 2,tick=FALSE)
        par(mar = c(12,0,0,0))
        image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
        #REMOVE x,y lim if cordinates don't match
              ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp))   
              )            
      }else{
        par(mar = c(12,12,0,0))
        image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
        #REMOVE x,y lim if cordinates don't match
              ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp))   
              )
        axis( 2,at=pos_vertical,labels=paste(rev(toplot$digital$ROImaster),"\n",labels_vert),las=1,tick=FALSE)
      }





      counter=seq(toplot$digital$nbin,toplot$digital$nbin*length(toplot$digital$roiTooverlap_range_list),toplot$digital$nbin)
      for (csep in counter){
        rect(xleft = csep , ybottom = rep(0,length(csep)), 
          xright = csep  + 0.01, 
          ytop = rep(ncol(trasp), csep), lty = 1, lwd = 1, 
          col = "black", border = "black")
      }
      #how many for each ROI? then revert, because the matrix is inverted in "image"
      #but the matrix original is correct, with ROI order and columns=bins; rows=ranges
      correctOrder=table(toplot$digital$labelstosplit)[order(order(unique(toplot$digital$labelstosplit)))]
      counter2=cumsum(rev(as.integer(correctOrder)))
      for (csep2 in counter2){
        rect(xleft = rep(0,length(csep2))  , ybottom = csep2 , 
          xright = rep(nrow(trasp) , csep2), 
          ytop = csep2 +0.01, lty = 1, lwd = 1, 
          col = "black", border = "black")
      }

      axis( 1, at=seq(toplot$digital$nbin/2,length(toplot$digital$names_samples)*toplot$digital$nbin+toplot$digital$nbin/2-1,toplot$digital$nbin ), 
            labels= toplot$digital$names_samples, las= 2,tick=FALSE )


      # axis( 2,at=pos_vertical,labels=paste(rev(toplot$digital$ROImaster),"\n",labels_vert),las=1,tick=FALSE)
      
      dev.off()
  } 
)


output$savefrequencyOvDigitalHeatbutton<- downloadHandler(
  filename=function() {
    paste('Overlap_bias.pdf', sep='')
  },
  content=function(file) {
    ovtoplot=toplot$digital$freq_overlap
    if(input$FracToPercDigitalHeat){
      for(i in 1:length(toplot$digital$freq_overlap)){
        ovtoplot[[i]]=lapply(toplot$digital$freq_overlap[[i]],function(k){k/toplot$digital$length_subsets[[i]]})
      }
      maxval=1
      ylabfreq="Fraction of overlaps"
    }else{
      maxval=max(toplot$digital$length_subsets)
      ylabfreq="Number of overlaps"
    }
    pdf(file)
    #now plot (either % or absolute numbers)
    plot(0,type='n',xaxt="n",ylim=c(0,maxval),xlim=c(1,toplot$digital$nbin),ylab=ylabfreq,xlab="Genomic Window",
                  main="Overlaps bias",cex.main=1.2,cex=1,xaxt="n",cex.lab=1.2,cex.axis=1)
    count=1
    legnames=c()
    for(i in 1:length(ovtoplot)){
      for(k in 1:length(ovtoplot[[i]])){
        #draw the lines/points
        lines(ovtoplot[[i]][[k]],col=toplot$digital$color_distinct_bias[count])
        legnames[count]=paste(names(ovtoplot[[i]])[k],"in",names(ovtoplot)[i])
        count=count+1
      }
    }
    axis(1,at=seq(1,toplot$digital$nbin,(toplot$digital$nbin-1)/2),
                labels=c( "start","midpoint","end" ),cex.axis=1)
    #legend for each line
    legend("topleft",legend=legnames,lty=rep(1,length(legnames)),col=toplot$digital$color_distinct_bias,bg="transparent")    
    dev.off()
  }
)

#observer to plot multiple overlaps (venn or heatmap)
output$saveJaccardDigitalbutton<- downloadHandler(
  filename=function() {
      paste('Multiple_overlaps.pdf', sep='')
  },
  content=function(file) {
        matList=toplot$digital$matlistTotal[[1]]

        if(length(matList)>1){
          
          volumes=c("UserFolder"="/")
          fileinfo <- parseSavePath(volumes, input$saveJaccardDigital)
          pdf(file)

          if(length(matList)==2 | length(matList)==3){
            plot(toplot$digital$jaccardorvenn)  
          }
          # if(length(matList)==400 | length(matList)==500){
          #   plot(toplot$digital$jaccardorvenn,type="ChowRuskey")         
          # }
          if(length(matList)>3){
            absmax=max(abs(toplot$digital$jaccardorvenn))
            trasp=t(toplot$digital$jaccardorvenn)
            brk=c( seq( 0 , 1,0.01 ))
            my_palette <- colorRampPalette(c("white","blue"))(n = length(brk)-1 ) 
            strlengths=nchar(colnames( trasp ))
            if(max(strlengths)<=25){
              margin_12=13
            }else{
              margin_12=13+ ((max(strlengths)-25) /3)
            }
            if(margin_12>16){
              margin_12=16
            }
            par(mar=c(margin_12,margin_12,3,3),xpd=TRUE)
            image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk,main="Jaccard index")
            axis( 2, at=seq(0.5,ncol(trasp)+0.5-1,1 ), labels= rev(colnames( trasp )), las= 2 )
            axis( 1, at=seq(0.5,ncol(trasp)+0.5-1,1 ), labels= colnames( trasp ), las= 2 )
            for (x in (nrow(toplot$digital$jaccardorvenn)-1+0.5):0.5  )
              for (y in 0.5: ((ncol(toplot$digital$jaccardorvenn)-1+0.5)   ))
                text(y,x, round(toplot$digital$jaccardorvenn[ncol(toplot$digital$jaccardorvenn)-x+0.5,y+0.5],1),col="red")
          }
          dev.off()
        }
  } 
)




##respond to "new ROI" button of digital heatmap:
observeEvent(input$confirmImportROIfromDigitalHeat,{
  #essential checks: final range must exist, finalfix must exist
  if(length(toplot$digital$range_tokeep)>0 & length(toplot$digital$fix_tokeep)>0){
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$newROIfromDigitalHeat)>=1 & !(any(input$newROIfromDigitalHeat == nottobe))) {
      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$newROIfromDigitalHeat) &!grepl(nottobe2[2],input$newROIfromDigitalHeat) & !grepl(nottobe2[3],input$newROIfromDigitalHeat)){
        nomi=unlist(lapply(ROIvariables$listROI,getName))
        if (!input$newROIfromDigitalHeat %in% nomi){
          #totrow=length(toplot$digital$range_tokeep)

          ##implement a check for NULL bam/WIG files?
          ##extract range, fix, bamlist ordered corresponding to the cluster selected
          range_sel=toplot$digital$range_tokeep[toplot$gadgetdigital$ybottom:toplot$gadgetdigital$ytop]
          bams_sel=lapply(toplot$digital$bams_tokeep,function(bamob){bamob[toplot$gadgetdigital$ybottom:toplot$gadgetdigital$ytop]})
          fix_sel=toplot$digital$fix_tokeep[toplot$gadgetdigital$ybottom:toplot$gadgetdigital$ytop]

          # #take from listROI the toplot$digital$labelstosplit name
          label_selected=unique(toplot$digital$labelstosplit)
          nomitot=unlist(lapply(ROIvariables$listROI,getName))
          postot=match(unique(label_selected),nomitot)


          originalROI=ROIvariables$listROI[[postot]]
          oldSource=getSource(originalROI)
          flag=getFlag(originalROI)
          toadd=paste("extracted ",length(range_sel)," ranges from digital heatmap",sep="")
          newSource=c(oldSource,list(toadd))    

          ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                name=input$newROIfromDigitalHeat,
                                range=range_sel,
                                fixed=fix_sel,
                                BAMlist=bams_sel,
                                flag=flag,
                                source=newSource)      
          logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$newROIfromDigitalHeat,' ROI, extracting ',length(fix_sel),' elements from digital heatmap<br>',sep="")
          print(paste('Created ',input$newROIfromDigitalHeat,' ROI, extracting ',length(fix_sel),' elements from a cluster of digital heatmap',sep=""))
          sendSweetAlert(
            session = session,
            title = "ROI extracted!",
            text = paste('Created ',input$newROIfromDigitalHeat,' ROI, extracting ',length(fix_sel),' elements from a cluster of digital heatmap',sep=""),
            type = "success"
          )          
        }else{
          sendSweetAlert(
            session = session,
            title = "Bad name",
            text = paste("\'",input$newROIfromDigitalHeat,"'\' ROI is already present. Choose another name",sep=""),
            type = "error"
          )  
        }
      }else{
        sendSweetAlert(
          session = session,
          title = "Bad name",
          text = paste("\'",input$newROIfromDigitalHeat,"'\' name for the new ROI is not allowed (must not start with 'promoters_genelist*','transcripts_genelist*','TES_genelist*')",sep=""),
          type = "error"
        ) 
      }
    }else{
      sendSweetAlert(
        session = session,
        title = "Bad name",
        text = paste("\'",input$newROIfromDigitalHeat,"'\' name for the new ROI is not allowed (must not be 'promoters','transcripts','TES') or is missing",sep=""),
        type = "error"
      )  
    }
  }else{
    #but should not appear...
    sendSweetAlert(
      session = session,
      title = "Data not found",
      text = "I didn't find the data from the heatmap",
      type = "error"
    ) 
  }

})

















##################################################################
##################################################################
##################################################################
# respond to plot Analogic Heatmap
##################################################################
##################################################################
##################################################################

###react to help buttons:
#parameters
observeEvent(input$msg_analogicHeatmap_parameters, {
  boxHelpServer(msg_analogicHeatmap_parameters)
})
#heatmap
observeEvent(input$msg_analogicHeatmap_heatmap, {
  boxHelpServer(msg_analogicHeatmap_heatmap)
})
#profiles
observeEvent(input$msg_analogicHeatmap_profiles, {
  boxHelpServer(msg_analogicHeatmap_profiles)
})
#enrichments
observeEvent(input$msg_analogicHeatmap_enrichments, {
  boxHelpServer(msg_analogicHeatmap_enrichments)
})



toListenAnalogHeat <- reactive({
    list(input$confirmUpdateAnalogHeat,input$confirmUpdateAnalogHeat2)
})
observeEvent(toListenAnalogHeat(),{
  
  set.seed(123)
  #checks: ROIlist must be >0 and lelected ROI must be >0, otherwise plot NULL
  #        total length of union of ROI selected must be >0 (at least one range) ,otherwise plot NULL
  #        BAMlist selected must be >0, otherwise plot is NULL
  #        other parameters should be ok for situation, otherwise correct the update UI part accordingly 
        #if sample random > total length, take the total length

  samplerandom=input$sampleRandomAnalogHeat

  toplot$analogic$samplerandom=samplerandom
  toplot$analogic$ROIsForAnalogHeat=input$ROIsForAnalogHeat
  toplot$analogic$binsAnalogHeat=input$binsAnalogHeat
  toplot$analogic$BAMsForAnalogHeat=input$BAMsForAnalogHeat

  ################################################################################################
  ##here, if reorder BAM, act here to change the order of enrichments in toplot$analogic$BAMsForAnalogHeat variable
  getbam=toplot$analogic$BAMsForAnalogHeat
  allnumbers=as.character(1:length(getbam))
  listprovv=list()
  for (i in 1:length(getbam)){
    stringval=grep(paste("reorderBAManalogHeat",i,"$",sep=""),names(input),value=TRUE)
    #extract what is contained inside each cell
    listprovv[[i]]=input[[ stringval ]]
  }
  listprovv=as.numeric(unlist(listprovv))
  #check if new order comprises all the possible positions
  if (identical(unique(sort(listprovv)),unique(sort(as.numeric(allnumbers)))) ){
    #reorder BAM 
    toplot$analogic$BAMsForAnalogHeat=getbam[order(listprovv)]
  }else{
    sendSweetAlert(
      session = session,
      title = "Ordering problem",
      text = "You didn't put all the ranking positions correct in the new ordering",
      type = "error"
    ) 
    return()     
  } 
  ################################################################################################








  if (length(ROIvariables$listROI)>0 & length(input$ROIsForAnalogHeat)>0 & length(toplot$analogic$BAMsForAnalogHeat)>0 & input$binsAnalogHeat>0 & samplerandom>0 & isvalid(input$sampleRandomAnalogHeat) & isvalid(input$binsAnalogHeat)){

    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForAnalogHeat,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    maxbinstohave=min(unlist(lapply(roi,checkMaxBins)))

    toplot$analogic$maxbinstohave=maxbinstohave

    if(input$binsAnalogHeat<=maxbinstohave){
      nbintouse=input$binsAnalogHeat
      nameROI=unlist(lapply(roi,getName))
      bigrange=list()
      bigbamlist=list()
      fixes=list()
      BAMs=list()
      bamnames=list()
      completerange=list()
      for(i in 1:length(roi)){
        #may cause problems: slow and when re-build ROI from heatmap selection,
        #strands "-" are inverted, while should not be...
        #unifROI=unifyStrand(roi[[i]])
        unifROI=roi[[i]]
        unifROI=uniqueROI(unifROI)
        fixes[[i]]=granges(getFixed(unifROI))
        completerange[[i]]=getRange(unifROI)
        #we have to collapse => all ranges must have the same columns
        bigrange[[i]]=granges(completerange[[i]])
        bigrange[[i]]$label=nameROI[i]      
        totbams=getBAMlist(unifROI)
        BAMs[[i]]=totbams
        #order even if BAM reordered, should be fine
        bamnames[[i]]=names(totbams)
        pos2=match(toplot$analogic$BAMsForAnalogHeat,bamnames[[i]])
        bigbamlist[[i]]=totbams[pos2]
        
      }

      #collapse all selected ranges into a single, big range list
      #this is useful for random sampling, to keep the proportions of the lengths
      
      finalrange=Reduce(c,bigrange)
      # finalfixes=Reduce(c,fixes)
      
      #print("prepping strands/ROIs")
      if (length(finalrange)>0){
        #if sample is greater than the length, set it as the length
        if(samplerandom>length(finalrange)){
          samplerandom=length(finalrange)
        }

        #from now, bigrange is the list of ROI selected to put in heatmap
        
      

        #if parameters of the clustering (numeric inputs) are not set (NA),
        #or not valid, in this case, no clustering. 
        #1- for any clustering, the number must be set.
        #2- If kmeans is without valid
        #advanced numeric parameters, set hclust with fake (no) clustering (see above...)
        checkparsclust=TRUE

        if(input$chooseOrderingAnalogHeat=="clustering"){
          #check number cluster
          if(!isvalid(input$clustnumAnalogHeat)){
            checkparsclust=FALSE
          }
          #check advanced K-means parameters
          if(input$clusterTypeAnalogHeat=="kmean"){
            if(!isvalid(input$clustrandomstartsAnalogHeat) | !isvalid(input$clustnumiterationsAnalogHeat) |length(input$BAMsForClusteringAnalogHeat)==0){
              checkparsclust=FALSE
            }
          }
        }




        if(checkparsclust){
          #sample random according to the input$sampleRandomAnalogHeat (maximum:total length of sum of ROIs)
          #label of which ROI belongs to is preserved
          finalBAMs_sample_pos=sort(sample(1:length(finalrange),samplerandom,replace=FALSE))
          #retrieve the labels of ranges (which ROI belongs to) after subsampling
          finalrange_sampled=finalrange[finalBAMs_sample_pos]
          # fixes=finalfixes[finalBAMs_sample_pos]
          labelstosplit=finalrange_sampled$label


          criteria=TRUE
          if(!is.null(input$clustnumAnalogHeat) &input$chooseOrderingAnalogHeat=="clustering"){
            if(length(labelstosplit)>0){
              if(input$clustnumAnalogHeat>=min(table(labelstosplit))){
                criteria=FALSE
              }
            }
          }

          #check: cluster number must be < min of labelstosplit label,
          if(criteria ){
            #if sample is so strict that removes a poorly rpresentative roi,
            #stop execution, asking for more intervals
            if(length(roi)==length(unique(labelstosplit))){
              strand_sampled=split( as.factor(strand(finalrange_sampled)) ,factor(labelstosplit,levels=unique(labelstosplit)))
              positions_splitted=split( finalBAMs_sample_pos ,factor(labelstosplit,levels=unique(labelstosplit)))
              #print("preparing material for ROI extraction...")

              #########prepare for ROI xtraction from heatmap#############
              finalbam=list()
              for(i in 1:length(BAMs)){
                finalbam[[i]]=as.list(rep(NA,length(BAMs[[i]])))
              }
              cumulatesum=0
              #for each ROI
              for(i in 1:length(BAMs)){
                #for each BAM, subselect only for those selected in that label:
                #it's a RELATIVE position from the beginning of that label, so subtract what is previous
                for (k in 1:length(BAMs[[i]])){
                  finalbam[[i]][[k]]=BAMs[[i]][[k]][positions_splitted[[i]]-cumulatesum]
                }
                names(finalbam[[i]])=bamnames[[i]]
                #the same of ROIs (complete, if have different columns):
                completerange[[i]]=completerange[[i]][positions_splitted[[i]]-cumulatesum]
                fixes[[i]]=fixes[[i]][positions_splitted[[i]]-cumulatesum]
                cumulatesum=cumulatesum+length(BAMs[[i]][[1]])
              }
              names(finalbam)=input$ROIsForAnalogHeat
              #split also finalrange_sampled and fixes according to "labelstosplit"
              # finalrange_sampled_split=split( finalrange_sampled ,factor(labelstosplit,levels=unique(labelstosplit)))
              # fixes_sampled_split=split( fixes ,factor(labelstosplit,levels=unique(labelstosplit)))          
              ############################################################
              #print("prepared...")

              #combine lists from different ROIs (c)
              subselectedBAMlist=list()
              for(i in 1:length(bigbamlist[[1]])){
                provv=list()
                #for each ROI
                for (k in 1:length(bigbamlist)){
                  provv=c(provv,bigbamlist[[k]][[i]])
                }   
                #subselect "provv": provv is the big union of all ROI for a specific BAM
                provv=provv[finalBAMs_sample_pos]
                subselectedBAMlist[[i]]=provv
              }



              #print("subselect")

              #use only finalBAMs_sample_pos positions sample. This is because the "random" sample
              #number refers to the total number of rows plotted, and not for each range.
              #this means that we have to merge matrixes first and then subsample. This will keep 
              #the proportion between ROIs intact.


              #re-split the matlists according to the correct label (sampled with the same position)
              #ROIs must be in the higher order hierarchy of list 
              slicedbamlist=as.list(rep(NA,length(toplot$analogic$BAMsForAnalogHeat)))
              slicedbamlist=rep(list(slicedbamlist),length(input$ROIsForAnalogHeat))
              
              #print("prepping BAM pieces to process")
              for(i in 1:length(subselectedBAMlist)){
                slicedbamlist2=split(subselectedBAMlist[[i]],factor(labelstosplit,levels=unique(labelstosplit)))
                #for each ROI obtained after split
                for(k in 1:length(slicedbamlist2)){
                  slicedbamlist[[k]][[i]]=slicedbamlist2[[k]]
                }      
              }
              
              #transform the lists of baseCoverage output in binned matrixes lists
              #check if ROI have fixed size: if so, create the matrixes without bins, then bin
              #but matrix without bins is used for profiles

              nomi=unlist(lapply(ROIvariables$listROI,getName))
              pos=match(input$ROIsForAnalogHeat,nomi)
              roisSelected=ROIvariables$listROI[pos]
              getwdth=lapply(roisSelected,getWidth) 
              widths=Reduce(union,getwdth)

              if(length(slicedbamlist)<length(slicedbamlist[[1]])){
                #if BAMs>ROIs, paralllize at BAM level

                matlists=lapply(1:length(slicedbamlist),function(i) {
                  #in case, filter for transcript flag if mem problems

                  if(nc==1){
                    matlist=lapply(1:length(slicedbamlist[[i]]),function(k) {

                                        return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE))
                                      })
                  }else{

                    decision=0
                    tryCatch({
                      matlist=mclapply(1:length(slicedbamlist[[i]]),function(k) {
                                        return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE))
                                      },mc.cores=nc) 
                      #if multicore returned correctly, decision==1, following if won't be calculated..
                      decision=1
                    },
                    warning = function( w ){
                      print("warning: coverage in single core calculation, due to memory limit...")
                    },
                    error = function( err ){
                      print("coverage in single core calculation, due to memory limit...")
                    })

                    if(decision==0){
                      print("Detected problems, maybe memory limits...")
                      matlist=lapply(1:length(slicedbamlist[[i]]),function(k) {
                                        return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE))
                                      })
                    }
                  }

                  ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                  ###is a ROI. 
                  pos_toinvert=as.character(strand_sampled[[i]])=="-"
                  
                  #for k in matlist (k should be the bam)
                  for(k in 1:length(matlist)){
                    #for each of the BAM file in this current ROI;
                    matprovv=matlist[[k]]

                    if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                      matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
                    }
                    
                  }

                  return(matlist)    
                })
              }else{

                #if ROIs>BAMs, parallelize at ROI level
                if(nc==1){
                  matlists=lapply(1:length(slicedbamlist),function(i) {
                    #in case, filter for transcript flag if mem problems

                    matlist=lapply(1:length(slicedbamlist[[i]]),function(k) {
                                          return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE))
                                        }) 
                    ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                    ###is a ROI. 
                    pos_toinvert=as.character(strand_sampled[[i]])=="-"
                    for(k in 1:length(matlist)){
                      #for each of the BAM file in this current ROI;
                      matprovv=matlist[[k]]
                      if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                        matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
                      }
                    }
                    return(matlist)    
                  })                
                }else{
                  decision=0
                  tryCatch({
                    matlists=mclapply(1:length(slicedbamlist),function(i) {
                      #in case, filter for transcript flag if mem problems

                      matlist=lapply(1:length(slicedbamlist[[i]]),function(k) {
                                            return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE))
                                          }) 
                      ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                      ###is a ROI. 
                      pos_toinvert=as.character(strand_sampled[[i]])=="-"
                      for(k in 1:length(matlist)){
                        #for each of the BAM file in this current ROI;
                        matprovv=matlist[[k]]
                        if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                          matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
                        }
                      }
                      return(matlist)    
                    },mc.cores=nc)  
                    #if multicore returned correctly, decision==1, following if won't be calculated..
                    decision=1
                  },
                  warning = function( w ){
                    print("warning: coverage in single core calculation, due to memory limit...")
                  },
                  error = function( err ){
                    print("coverage in single core calculation, due to memory limit...")
                  })              
                  
                  if(decision==0){
                    print("Detected problems, maybe memory limits...")
                    matlists=lapply(1:length(slicedbamlist),function(i) {
                      #in case, filter for transcript flag if mem problems

                      matlist=lapply(1:length(slicedbamlist[[i]]),function(k) {
                                            return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE))
                                          }) 
                      ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                      ###is a ROI. 
                      pos_toinvert=as.character(strand_sampled[[i]])=="-"

                      #for k in matlist (k should be the bam)
                      for(k in 1:length(matlist)){
                        #for each of the BAM file in this current ROI;
                        matprovv=matlist[[k]]

                        if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                          matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
                        }
                        
                      }
                      return(matlist)    
                    }) 

                  }
                }
              }



              #return matlists in any case. If fixed size is the same for all, return complete matrix
              #print("matrixes from basecov done...")
              


              #transform the matrixlists using clustering or ordering. Keep the order
              #find cluster inds

              toorder_final=list()
              if(input$chooseOrderingAnalogHeat=="clustering"){
                heatvariables$clustnumAnalogHeat= input$clustnumAnalogHeat
                tempclustnum=input$clustnumAnalogHeat
                pos=match(input$BAMsForClusteringAnalogHeat,toplot$analogic$BAMsForAnalogHeat)
                if (length(pos)==0){
                  pos=NULL
                }

                #choose depending on hclust or kmeans (on kmeans, do a check on number
                #of clusters, between 2 and 433 included)
                if (input$clusterTypeAnalogHeat=="hierarchical"){
                  if(class(input$distmethodAnalogHeat)!="character" |class(input$clustmethodAnalogHeat)!="character"){
                    distmethod="euclidean"
                    clustmethod="average"
                  }else{
                    distmethod=input$distmethodAnalogHeat
                    clustmethod=input$clustmethodAnalogHeat
                  }
                  

                  if(nc==1){
                    clusteringblock=lapply(1:length(matlists),function(i){
                      clust=clusterMatrix(matlist=matlists[[i]],distmethod=distmethod,clustmethod=clustmethod,clustinds=pos)
                      return(clust)
                    })
                  }else{
                    decision=0
                    tryCatch({
                      clusteringblock=mclapply(1:length(matlists),function(i){
                        clust=clusterMatrix(matlist=matlists[[i]],distmethod=distmethod,clustmethod=clustmethod,clustinds=pos)
                        return(clust)
                      },mc.cores=nc) 
                      decision=1
                    },
                    warning = function( w ){
                      print("warning: clustering in single core calculation, due to memory limit...")
                    },
                    error = function( err ){
                      print("clustering in single core calculation, due to memory limit...")
                    })  
                    if(decision==0){
                      print("cannot exec parallel clustering...")
                      clusteringblock=lapply(1:length(matlists),function(i){
                        clust=clusterMatrix(matlist=matlists[[i]],distmethod=distmethod,clustmethod=clustmethod,clustinds=pos)
                        return(clust)
                      })
                    }               
                  }
                  finalmats=lapply(clusteringblock,function(cl){cl$mat})
                  toorder_final=lapply(clusteringblock,function(cl){cl$ord})

                #implements k-means clustering
                }else{
                  #if iterations or starting points are <=0, put 1 ... if they are > 40, put 40
                  startingpoints=input$clustrandomstartsAnalogHeat
                  iter=input$clustnumiterationsAnalogHeat
                  if (is.null(startingpoints)){
                    startingpoints=1
                  }else{
                    if(startingpoints<=0){
                      startingpoints=1
                    }
                    if(startingpoints>40){
                      startingpoints=40
                    }
                  }

                  if(is.null(iter)){
                    iter=1
                  }else{
                    if(iter<=0){
                      iter=1
                    }
                    if(iter>40){
                      iter=40
                    }                
                  }

                  #inside clusterMatrixKmeans implements check on the number of clusters (0<x<=433)
                  #lse return NULL cluster object and the ranking is the original (no particular clustering or ranking)
                  
                  if(nc==1){
                    clusteringblock=lapply(1:length(matlists),function(i){
                      clust=clusterMatrixKmeans(matlist=matlists[[i]],clustinds=pos,numberclusters=tempclustnum,startingpoints=startingpoints,iter=iter)
                      return(clust)
                    })
                  }else{
                    decision=0
                    tryCatch({
                      clusteringblock=mclapply(1:length(matlists),function(i){
                        clust=clusterMatrixKmeans(matlist=matlists[[i]],clustinds=pos,numberclusters=tempclustnum,startingpoints=startingpoints,iter=iter)
                        return(clust)
                      },mc.cores=nc) 
                      decision=1
                    },
                    warning = function( w ){
                      print("warning: clustering in single core calculation, due to memory limit...")
                    },
                    error = function( err ){
                      print("clustering in single core calculation, due to memory limit...")
                    })  
                    if(decision==0){
                      print("cannot exec parallel clustering...")
                      clusteringblock=lapply(1:length(matlists),function(i){
                        clust=clusterMatrixKmeans(matlist=matlists[[i]],clustinds=pos,numberclusters=tempclustnum,startingpoints=startingpoints,iter=iter)
                        return(clust)
                      })
                    }               
                  }
                  finalmats=lapply(clusteringblock,function(cl){cl$mat})
                  toorder_final=lapply(clusteringblock,function(cl){cl$ord})
                }

              }else{
                heatvariables$clustnumAnalogHeat=NA
                #ranking
                pos=match(input$BAMsForRankingAnalogHeat,toplot$analogic$BAMsForAnalogHeat)
                finalmats=list()
                for(i in 1:length(matlists)){
                  #rank the heatmap identified with "pos"
                  mat_capitain=matlists[[i]][[pos]]
                  values=apply(mat_capitain,1,sum)
                  ord_capitain=order(-values)
                  listmatprovv=lapply(matlists[[i]],function(k) {return(k[ord_capitain,,drop=FALSE])})
                  toorder_final[[i]]=ord_capitain
                  listmatprovv=do.call(cbind,listmatprovv)
                  finalmats[[i]]=listmatprovv
                }
              }
              
              #print("clustering done if selected...")
              #if completematrix is not NULL (fixed size, cluster it with the same order as returned by
              # finalmats!)

              #order all the material for ROI extraction, for each ROI depicted, one order
              for(i in 1:length(finalmats)){
                completerange[[i]]=completerange[[i]][toorder_final[[i]]]
                fixes[[i]]=fixes[[i]][toorder_final[[i]]]
                for(k in 1:length(finalbam[[i]])){
                  finalbam[[i]][[k]]=finalbam[[i]][[k]][toorder_final[[i]]]
                }
              }
              
              toplot$analogic$finalrange=completerange
              toplot$analogic$finalbam=finalbam
              toplot$analogic$finalfixes=fixes


              ########################################################################
              ########################################################################
              ########################################################################
              #finalmats are the final matrixes splitted by ROI
              #keep this in mind for retrieving data from interacting heatmap
              #need a reactive variable for keeping it. Then reactive to click brush
              #from the stored heatmap and retrieve data and plot and so on...
              ########################################################################
              ########################################################################
              ########################################################################
              heatvariables$matlist=finalmats
              # heatvariables$completematrixes=completematrixes

              heatvariables$ROIsForAnalogHeat=input$ROIsForAnalogHeat
              heatvariables$binsAnalogHeat=input$binsAnalogHeat
              heatvariables$BAMsForAnalogHeat=toplot$analogic$BAMsForAnalogHeat
              heatvariables$sampleRandomAnalogHeat=samplerandom

              # heatvariables$widthfix=widthfix
              #finally,combine the matrixes. Be VERY careful about the order.
              
              matrixes_processed=do.call(rbind,finalmats)
              lengthsfinalmats=sapply(finalmats,nrow)

              colnames(matrixes_processed)=rep(heatvariables$BAMsForAnalogHeat,each=toplot$analogic$binsAnalogHeat)
              rownames(matrixes_processed)=rep(heatvariables$ROIsForAnalogHeat,times=lengthsfinalmats)
              
              ###########################
              #finally, plot the heatmap (with image, see previous examples) to be reactive to x;y
              #respond to colors and input$quantileThreshAnalogHeat DIRECTLY (inside output$plot)
              #if too slow to plot, put these parameters outside output$...
              numbersamples=length(toplot$analogic$BAMsForAnalogHeat)
              numberbins=input$binsAnalogHeat


              toplot$analogic$numbersamples=numbersamples
              toplot$analogic$method=input$chooseQuantileMethodAnalogHeat
              toplot$analogic$matrixes_processed=matrixes_processed
              toplot$analogic$labelstosplit=labelstosplit
              toplot$analogic$clusterTypeAnalogHeat=input$clusterTypeAnalogHeat

              #here, put the code to generate clustering number arrays to plot
              #and set a variable "clustering" if all conditions are met, to draw the side colors
              isclusteringok= input$chooseOrderingAnalogHeat=="clustering" & length(input$ROIsForAnalogHeat) ==1 & length(input$BAMsForClusteringAnalogHeat)>0 &
                              !is.na(heatvariables$clustnumAnalogHeat) & heatvariables$clustnumAnalogHeat<=433 & heatvariables$clustnumAnalogHeat>0
              
              toplot$analogic$clusternumbermatrix=NULL
              if(isclusteringok){
                set.seed(123)
                #for this line of code, thanks to
                # http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
                # Jelena-bioinf / Megatron
                color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
                color_distinct_cluster=sample(color, heatvariables$clustnumAnalogHeat)
                toplot$analogic$color_distinct_cluster=color_distinct_cluster
                #if clustering is hierarchical
                if (toplot$analogic$clusterTypeAnalogHeat=="hierarchical"){
                  # use clust$clustobject from clusterMatrix function
                  clusterarray=cutree(clusteringblock[[1]]$clustobject,k=heatvariables$clustnumAnalogHeat)
                #if clustering is K-means
                }else{
                  #in this case (kmeans) clusterarray is already defined by clust$clustobject$cluster
                  clusterarray=clusteringblock[[1]]$clustobject$cluster
                } 
                clusterarray=clusterarray[clusteringblock[[1]]$ord]
                clusterarray=rev(clusterarray)
                matrix_like=matrix(clusterarray,nrow=1)
                toplot$analogic$clusternumbermatrix=matrix_like    

              }else{
                output$clustersImageLeft<-renderPlot({NULL})
              }



              #plot cluster color heatmap on the left;
              if(isclusteringok){
                output$clustersImageLeft<-renderPlot({
                  par(mar = c(0,0,0,0))
                  image(0:nrow(matrix_like), 0:ncol(matrix_like),matrix_like,col=color_distinct_cluster,axes = FALSE, xlab = "", ylab = "",ylim=c(0,ncol(matrix_like)))
                })              
              }else{
                output$clustersImageLeft<-renderPlot({NULL}) 
              }




              #save heatmap data xls button
              output$showsaveAnalogHeatdata=renderUI({downloadButton('saveAnalogHeatdata', 'Download matrix data')})
              #clear button for new ROI: we have another analysis
              output$newROIfromAnalogHeat_out<-renderUI({NULL})



              #reset the color menu:
              if(input$optioncolorsforAnalogHeat=="custom"){
                bams=heatvariables$BAMsForAnalogHeat
                lista=list()

                if(length(bams)>0){
                  output$showcolorsheat<-renderUI({
                    for (i in 1:length(bams)){
                      lista[[i]]=fluidRow(
                                          column(12,
                                                colorSelectorInput(inputId=paste("colorCustomAnalogHeat",i,sep=""),label=bams[i],choices=ColsArray,
                                                              selected=ColsArray[1],mode="radio",ncol=length(ColsArray))
                                          )
                                  )

                    }
                    wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
                      lista
                    )
                  })        
                }else{
                  #bams are NULL. show nothing
                  output$showcolorsheat<-renderUI({NULL})
                }
              }


              
              output$heatmapAnalog<-renderPlot({

                if(numbersamples<10){
                  factormult=10/numbersamples
                }else{
                  factormult=1
                }
                matProc_analogic=toplot$analogic$matrixes_processed

                #modify max values above quantile threshold. one for each BAM (not fixed value)
                quantThresh=input$quantileThreshAnalogHeat
                nbin=heatvariables$binsAnalogHeat
                method=input$chooseQuantileMethodAnalogHeat

                if(method=="eachBAM"){
                  for(i in 1:length(heatvariables$BAMsForAnalogHeat)){
                    piecematrix=matProc_analogic[, ((i-1)*nbin+1):(i*nbin),drop=FALSE ]
                    maxval=max(piecematrix)
                    quantval=quantile(piecematrix,quantThresh)
                    if(quantval>0){
                      piecematrix[piecematrix>quantval]=quantval
                      piecematrix=piecematrix/quantval
                    }else{
                      piecematrix[piecematrix>quantval]=maxval
                      piecematrix=piecematrix/maxval
                    }
                    
                    matProc_analogic[, ((i-1)*nbin+1):(i*nbin) ]=piecematrix
                  }
                }else{
                  maxval=max(matProc_analogic)
                  quantval=quantile(matProc_analogic,quantThresh)
                  if(quantval>0){
                    matProc_analogic[matProc_analogic>quantval]=quantval
                    matProc_analogic=matProc_analogic/quantval
                  }else{
                    matProc_analogic[matProc_analogic>quantval]=maxval
                    matProc_analogic=matProc_analogic/maxval
                  }
                  
                }

                #extract the value of color palette, resulting from a complex evaluation
                #of reactive inputs (see updateGENOMICS)
                palette_col=toplot$analogic$colorpalettes

                #for plotting only, each block of the matrix, for each bam, sum 1
                brk=201

                if(isolate(input$optioncolorsforAnalogHeat)=="global"){
                  #here, only one palette for everything
                  if(palette_col=="rainbow"){
                    palette <- rev(colorRampPalette(brewer.pal(9, 'Spectral'))(n =brk-1))
                  }else{
                    colorsplitted=strsplit(palette_col,split="_")[[1]]
                    palette=colorRampPalette(colorsplitted)(n=brk-1)
                  }
                  
                }else{


                  if(length(palette_col)==1 ){
                    if(palette_col=="white_red4"){
                      colorsplitted=strsplit(palette_col,split="_")[[1]]
                      palette=colorRampPalette(colorsplitted)(n =brk-1)  
                    }else{
                      palette=c()
                      for(i in 1:length(palette_col)){
                        palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n =brk-1))
                      }
                    }
                              
                  }else{
                    #if input$optioncolorsforAnalogHeat is custom, palette should have length==nbams
                    palette=c()
                    for(i in 1:length(palette_col)){
                      palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n=brk-1))
                    }

                    #if different colors, we need to change the order of magnitude 
                    #for each bam
                    matProc_analogic[matProc_analogic>0.99]=0.999
                    matProc_analogic[matProc_analogic<0.01]=0.001
                    for(i in 1:length(heatvariables$BAMsForAnalogHeat)){
                        piecematrix=matProc_analogic[, ((i-1)*nbin+1):(i*nbin),drop=FALSE ]
                        asvec=unique(as.vector(piecematrix))
                        if(length(asvec)==1){
                          #all elments are the same, it can break the color palette, change one
                          if(asvec==0.999){
                            piecematrix[1,1]=0.001
                          }else if(asvec==0.001){
                            piecematrix[1,1]=0.999
                          }
                        }
                        piecematrix=piecematrix+(i-1)
                        matProc_analogic[, ((i-1)*nbin+1):(i*nbin) ]=piecematrix
                        #palettes[i]=colorRampPalette(colorsplitted)(n=brk)
                    } 


                  }

                }



                trasp=t(matProc_analogic)
                if(class(trasp)!="matrix"){
                  trasp=matrix(trasp,ncol=1)
                }
                
                if(max(dim(trasp))>30000){
                  raster=FALSE
                }else{
                  raster=TRUE
                }


                par(mar = rep(0, 4))
                #here we are plotting with the parameters
                image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
                #REMOVE x,y lim if cordinates don't match
                      ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp)),useRaster=raster  
                      )


                counter=seq(numberbins,numberbins*numbersamples,numberbins)
                for (csep in counter){
                  rect(xleft = csep , ybottom = rep(0,length(csep)), 
                    xright = csep  + 0.01, 
                    ytop = rep(ncol(trasp), csep), lty = 1, lwd = 1, 
                    col = "black", border = "black")
                }


                #how many for each ROI? then revert, because the matrix is inverted in "image"
                #but the matrix original is correct, with ROI order and columns=bins; rows=ranges
                correctOrder=table(labelstosplit)[order(order(unique(labelstosplit)))]
                counter2=cumsum(rev(as.integer(correctOrder)))

                for (csep2 in counter2){
                  rect(xleft = rep(0,length(csep2))  , ybottom = csep2 , 
                    xright = rep(nrow(trasp) , csep2), 
                    ytop = csep2 +0.01, lty = 1, lwd = 1, 
                    col = "black", border = "black")
                }

              })


              print("Drawn analogic heatmap")
              #button for downloading heatmap analogic
              output$saveheatmapAnalog=renderUI({downloadButton('saveheatmapAnalogbutton', 'Get PDF')})






              #output of the number of heatmap elements as fraction. In red if not taken all:
              output$textfractionelementsAnalogHeat<-renderText({
                if(samplerandom==length(finalrange)){
                  color="'green'"
                }else{
                  color="'red'"
                }
                paste("Intervals displayed: <font color=",color,">",samplerandom,"/",length(finalrange),"</font>",sep="")
              })






              #name of columns for heatmap
              output$textNameAnalogHeat <- renderPlot({
                par(mar = rep(0, 4))
                plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i")
                #have to start from x coordinate of the half of first cell to half of the last
                numbersample=length(heatvariables$BAMsForAnalogHeat)
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
                text(x=seq(halfcellwidthmin, halfcellwidthmax ,length.out=numbersample),labels=heatvariables$BAMsForAnalogHeat,y=rep(0.5,numbersample),srt=90,cex=cextext  )
                      
              })







              #plot color scale of heatmap
              output$colorScaleAnalogHeat<-renderPlot({
                brk=201
                palette_col=toplot$analogic$colorpalettes
                if(isolate(input$optioncolorsforAnalogHeat)=="global"){
                  
                  if(palette_col=="rainbow"){
                    my_palette <- rev(colorRampPalette(brewer.pal(9, 'Spectral'))(n = brk-1))
                  }else{
                    colorsplitted=strsplit(palette_col,split="_")[[1]]
                    my_palette <- colorRampPalette(colorsplitted)(n = brk-1 )
                  }
                  color.bar(my_palette, min=0,max=1)                    
                }else{
                  NULL
                }

              })


            }else{
              heatvariables$BAMsForAnalogHeat=NULL
              output$heatmapAnalog<-renderPlot({NULL})
              output$textfractionelementsAnalogHeat<-renderText({NULL})
              output$colorScaleAnalogHeat<-renderPlot({NULL})
              output$textNameAnalogHeat <- renderPlot({NULL})
              toplot$analogic$finalrange=NULL
              toplot$analogic$finalbam=NULL
              toplot$analogic$finalfixes=NULL
              heatvariables$matlist=NULL   
              output$clustersImageLeft<-renderPlot({NULL})    
              output$saveheatmapAnalog=renderUI({NULL})
              output$showsaveAnalogHeatdata=renderUI({NULL})
              output$newROIfromAnalogHeat_out<-renderUI({NULL})
              output$textselectedelementsAnalogHeat<-renderText({NULL}) 
              #reset color scheme. Put it to global with reset parameters
              updateRadioButtons(session,inputId="optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("global color","custom colors"),choiceValues=c("global","custom"),selected="global")
              sendSweetAlert(
                session = session,
                title = "Too few clusters",
                text = "Increase the number of clusters",
                type = "error"
              ) 
            }
          


          }else{
            heatvariables$BAMsForAnalogHeat=NULL
            output$heatmapAnalog<-renderPlot({NULL})
            output$textfractionelementsAnalogHeat<-renderText({NULL})
            output$colorScaleAnalogHeat<-renderPlot({NULL})
            output$textNameAnalogHeat <- renderPlot({NULL})
            toplot$analogic$finalrange=NULL
            toplot$analogic$finalbam=NULL
            toplot$analogic$finalfixes=NULL
            heatvariables$matlist=NULL   
            output$clustersImageLeft<-renderPlot({NULL})    
            output$saveheatmapAnalog=renderUI({NULL})
            output$showsaveAnalogHeatdata=renderUI({NULL})
            output$newROIfromAnalogHeat_out<-renderUI({NULL})
            output$textselectedelementsAnalogHeat<-renderText({NULL})
            updateRadioButtons(session,inputId="optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("global color","custom colors"),choiceValues=c("global","custom"),selected="global")

            sendSweetAlert(
              session = session,
              title = "Too few samples",
              text = "Increase the random sample (advanced parameters tab), because some ROI with few ranges disappeared",
              type = "error"
            )             
          }

          



        }else{
          #some numeric parameters (clustering) are not set, thus NA
          print("WARNING: some parameters not set in analogic heatmap...")
          heatvariables$BAMsForAnalogHeat=NULL
          output$heatmapAnalog<-renderPlot({NULL})
          output$textfractionelementsAnalogHeat<-renderText({NULL})
          output$colorScaleAnalogHeat<-renderPlot({NULL})
          output$textNameAnalogHeat <- renderPlot({NULL})
          toplot$analogic$finalrange=NULL
          toplot$analogic$finalbam=NULL
          toplot$analogic$finalfixes=NULL
          heatvariables$matlist=NULL   
          output$clustersImageLeft<-renderPlot({NULL})    
          output$saveheatmapAnalog=renderUI({NULL})
          output$showsaveAnalogHeatdata=renderUI({NULL})
          output$newROIfromAnalogHeat_out<-renderUI({NULL})
          output$textselectedelementsAnalogHeat<-renderText({NULL})
          updateRadioButtons(session,inputId="optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("global color","custom colors"),choiceValues=c("global","custom"),selected="global")

          sendSweetAlert(
            session = session,
            title = "No clustering parameters",
            text = "Some paramters for the clustering are not set properly",
            type = "error"
          )             
        }

        



      }else{
        heatvariables$BAMsForAnalogHeat=NULL
        output$heatmapAnalog<-renderPlot({NULL})
        output$textfractionelementsAnalogHeat<-renderText({NULL})
        output$colorScaleAnalogHeat<-renderPlot({NULL})
        output$textNameAnalogHeat <- renderPlot({NULL})
        toplot$analogic$finalrange=NULL
        toplot$analogic$finalbam=NULL
        toplot$analogic$finalfixes=NULL
        heatvariables$matlist=NULL
        output$clustersImageLeft<-renderPlot({NULL})
        output$saveheatmapAnalog=renderUI({NULL})
        output$showsaveAnalogHeatdata=renderUI({NULL})
        output$newROIfromAnalogHeat_out<-renderUI({NULL})
        output$textselectedelementsAnalogHeat<-renderText({NULL})
        updateRadioButtons(session,inputId="optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("global color","custom colors"),choiceValues=c("global","custom"),selected="global")
        sendSweetAlert(
          session = session,
          title = "Empty range",
          text = "The final range has length = 0. Check what's wrong with the selected ranges",
          type = "error"
        ) 
      }
    }else{
      heatvariables$BAMsForAnalogHeat=NULL
      output$heatmapAnalog<-renderPlot({NULL})
      output$textfractionelementsAnalogHeat<-renderText({NULL})
      output$colorScaleAnalogHeat<-renderPlot({NULL})
      output$textNameAnalogHeat <- renderPlot({NULL})
      toplot$analogic$finalrange=NULL
      toplot$analogic$finalbam=NULL
      toplot$analogic$finalfixes=NULL
      heatvariables$matlist=NULL
      #logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">Number of bins is > than width of smallest range. Decrease number of bins or filter range based on width<br></font>'
      sendSweetAlert(
        session = session,
        title = "Too many bins",
        text = "The selected number of bins exceeds the maximum allowed, because some ranges have length < nbins. Try with a lower number",
        type = "error"
      )      
      output$clustersImageLeft<-renderPlot({NULL})
      output$saveheatmapAnalog=renderUI({NULL})
      output$showsaveAnalogHeatdata=renderUI({NULL})
      output$newROIfromAnalogHeat_out<-renderUI({NULL})
      output$textselectedelementsAnalogHeat<-renderText({NULL})
      updateRadioButtons(session,inputId="optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("global color","custom colors"),choiceValues=c("global","custom"),selected="global")
    }
  }else{
    heatvariables$BAMsForAnalogHeat=NULL
    output$heatmapAnalog<-renderPlot({NULL})
    output$textfractionelementsAnalogHeat<-renderText({NULL})
    output$colorScaleAnalogHeat<-renderPlot({NULL})
    output$textNameAnalogHeat <- renderPlot({NULL})
    toplot$analogic$finalrange=NULL
    toplot$analogic$finalbam=NULL
    toplot$analogic$finalfixes=NULL
    heatvariables$matlist=NULL
    output$clustersImageLeft<-renderPlot({NULL})
    output$saveheatmapAnalog=renderUI({NULL})
    output$showsaveAnalogHeatdata=renderUI({NULL})
    output$newROIfromAnalogHeat_out<-renderUI({NULL})
    output$textselectedelementsAnalogHeat<-renderText({NULL})
    updateRadioButtons(session,inputId="optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("global color","custom colors"),choiceValues=c("global","custom"),selected="global")
    sendSweetAlert(
      session = session,
      title = "Wrong/missing parameters",
      text = "Some parameters are missing or wrong; check them and then click the button",
      type = "error"
    )

  }
},ignoreInit=TRUE)










#reactive for download button PDF of analog heatmap
output$saveheatmapAnalogbutton<- downloadHandler(
  filename=function() {
      paste('Heatmap_analogic.pdf', sep='')
  },
  content=function(file) {
        
          #find xlim according to how many numbersamples if numbersamples<8
          if(toplot$analogic$numbersamples<10){
            factormult=10/toplot$analogic$numbersamples
          }else{
            factormult=1
          }
          #quantThresh should be reactive to the change in color saturation bar
          quantThresh=input$quantileThreshAnalogHeat
          method=input$chooseQuantileMethodAnalogHeat
          nbin=heatvariables$binsAnalogHeat
          #matrix temporary for each change in graphical parameter. picks data from original matrix
          #toplot$analogic$matrixes_processed should not be modified after the big matrix computation
          temp_matrix_processed=toplot$analogic$matrixes_processed
          if(method=="eachBAM"){
            for(i in 1:length(heatvariables$BAMsForAnalogHeat)){
              piecematrix=temp_matrix_processed[, ((i-1)*toplot$analogic$binsAnalogHeat+1):(i*toplot$analogic$binsAnalogHeat),drop=FALSE ]
              piecematrix[piecematrix>quantile(piecematrix,quantThresh)]=quantile(piecematrix,quantThresh)
              piecematrix=piecematrix/max(piecematrix)
              temp_matrix_processed[, ((i-1)*toplot$analogic$binsAnalogHeat+1):(i*toplot$analogic$binsAnalogHeat) ]=piecematrix
            }
          }else{
            temp_matrix_processed[temp_matrix_processed>quantile(temp_matrix_processed,quantThresh)]=quantile(temp_matrix_processed,quantThresh)
            temp_matrix_processed=temp_matrix_processed/max(temp_matrix_processed)
          }

          #get colors
          palette_col=toplot$analogic$colorpalettes
          brk=201


        

          if(isolate(input$optioncolorsforAnalogHeat)=="global"){
            #here, only one palette for everything
            if(palette_col=="rainbow"){
              palette <- rev(colorRampPalette(brewer.pal(9, 'Spectral'))(n =brk-1))
            }else{
              colorsplitted=strsplit(palette_col,split="_")[[1]]
              palette=colorRampPalette(colorsplitted)(n=brk-1)
            }
            
          }else{

            if(length(palette_col)==1 ){
              if(palette_col=="white_red4"){
                colorsplitted=strsplit(palette_col,split="_")[[1]]
                palette=colorRampPalette(colorsplitted)(n =brk-1)  
              }else{
                palette=c()
                for(i in 1:length(palette_col)){
                  palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n =brk-1))
                }
              }
                        
            }else{
              #if input$optioncolorsforAnalogHeat is custom, palette should have length==nbams
              palette=c()
              for(i in 1:length(palette_col)){
                palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n=brk-1))
              }

              #if different colors, we need to change the order of magnitude 
              #for each bam
              temp_matrix_processed[temp_matrix_processed>0.99]=0.999
              temp_matrix_processed[temp_matrix_processed<0.01]=0.001
              for(i in 1:length(heatvariables$BAMsForAnalogHeat)){
                  piecematrix=temp_matrix_processed[, ((i-1)*nbin+1):(i*nbin),drop=FALSE ]
                  asvec=unique(as.vector(piecematrix))
                  if(length(asvec)==1){
                    #all elments are the same, it can break the color palette, change one
                    if(asvec==0.999){
                      piecematrix[1,1]=0.001
                    }else if(asvec==0.001){
                      piecematrix[1,1]=0.999
                    }
                  }
                  piecematrix=piecematrix+(i-1)
                  temp_matrix_processed[, ((i-1)*nbin+1):(i*nbin) ]=piecematrix
                  #palettes[i]=colorRampPalette(colorsplitted)(n=brk)
              } 
            }
          }
                  

          


          trasp=t(temp_matrix_processed)

          #insert ROI analyzed as axis in position 2:
          #ROI names are heatvariables$ROIsForAnalogHeat
          #shuld write labels weighted for the size of each ROI
          #extract position of ROIs as rows (how many ranges for each ROI) to determine the 
          #vertical positions
          intervals_vertical=table(toplot$analogic$labelstosplit)[order(order(unique(toplot$analogic$labelstosplit)))]
          #revert order because the "0" starts at the bottom, (last ROI)
          intervals_vertical=rev(intervals_vertical)
          cumulative=cumsum(intervals_vertical)
          tosubtract=intervals_vertical/2
          pos_vertical=cumulative-tosubtract

          #calculate % out of the total length of ROIs considered...
          percs_rois=round({intervals_vertical/sum(intervals_vertical)}*100,2)
          labels_vert=paste("(",unname(intervals_vertical),"; ",percs_rois,"%)",sep="")


          pdf(file)

          #if cluster image left present, plot it
          mlay=matrix(c(1,2),ncol=2,nrow=1)
          if(!is.null(toplot$analogic$clusternumbermatrix)){
            layout(mlay,widths=c(40,60),heights=c(100,100))
          }else{
            #layout(mlay,widths=c(0,100),heights=c(100,100))
          }


          

          if(!is.null(toplot$analogic$clusternumbermatrix)){
            par(mar = c(12,12,0,0))
            image(0:nrow(toplot$analogic$clusternumbermatrix), 0:ncol(toplot$analogic$clusternumbermatrix),toplot$analogic$clusternumbermatrix,col=toplot$analogic$color_distinct_cluster,axes = FALSE, xlab = "", ylab = "",ylim=c(0,ncol(trasp)))
            axis( 2,at=pos_vertical,labels=paste(rev(heatvariables$ROIsForAnalogHeat),"\n",labels_vert),las=1,tick=FALSE)
            axis(1,0.5,labels="cluster",las= 2,tick=FALSE)
            par(mar = c(12,0,0,0))
            image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
            #REMOVE x,y lim if cordinates don't match
                  ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp))   
                  )            
          }else{
            par(mar = c(12,12,0,0))
            image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
            #REMOVE x,y lim if cordinates don't match
                  ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp))   
                  )
            axis( 2,at=pos_vertical,labels=paste(rev(heatvariables$ROIsForAnalogHeat),"\n",labels_vert),las=1,tick=FALSE)
            
          }



          #for each BAM (all bins/BAM), draw a vertical black line to split BAMs
          counter=seq(toplot$analogic$binsAnalogHeat,toplot$analogic$binsAnalogHeat*toplot$analogic$numbersamples,toplot$analogic$binsAnalogHeat)
          for (csep in counter){
            rect(xleft = csep , ybottom = rep(0,length(csep)), 
              xright = csep  + 0.01, 
              ytop = rep(ncol(trasp), csep), lty = 1, lwd = 1, 
              col = "black", border = "black")
          }


          #for each ROI, draw a horizontal black line to split ROIs (if >1)
          #how many for each ROI? then revert, because the matrix is inverted in "image"
          #but the matrix original is correct, with ROI order and columns=bins; rows=ranges
          correctOrder=table(toplot$analogic$labelstosplit)[order(order(unique(toplot$analogic$labelstosplit)))]
          counter2=cumsum(rev(as.integer(correctOrder)))

          for (csep2 in counter2){
            rect(xleft = rep(0,length(csep2))  , ybottom = csep2 , 
              xright = rep(nrow(trasp) , csep2), 
              ytop = csep2 +0.01, lty = 1, lwd = 1, 
              col = "black", border = "black")
          }

          
          #insert text as axes in the same graphical object as heatmap
          axis( 1, at=seq(toplot$analogic$binsAnalogHeat/2,length(heatvariables$BAMsForAnalogHeat)*toplot$analogic$binsAnalogHeat+toplot$analogic$binsAnalogHeat/2-1,toplot$analogic$binsAnalogHeat ), 
                labels= heatvariables$BAMsForAnalogHeat, las= 2,tick=FALSE )

          
          dev.off()
  } 
)



#save heatmap data (xls)
output$saveAnalogHeatdata<-downloadHandler(
  filename=function() {
      paste('heatmap_data.xls', sep='')
  },
  content=function(file) {
      write.table(toplot$analogic$matrixes_processed,file=file,row.names=TRUE,sep="\t",quote=FALSE,col.names=NA   ) 
  }   
)



############OBSERVERS for plots interactive upon analogic heatmap selection######




output$saveprofileAnalogHeatbutton<- downloadHandler(
  filename=function() {
      paste('Profile_analog_heatmap.pdf', sep='')
  },
  content=function(file) {
      pdf(file)
      
      if(input$Log2BoxAnalogHeat){
        portionlist_profile=lapply(toplot$gadgetanalogic$portionlist_profile,log2)
        yl="Log2 read density (rpm/bp)"
      }else{
        portionlist_profile=toplot$gadgetanalogic$portionlist_profile
        yl="Read density (rpm/bp)"
      }
      maxval=max(unlist(lapply(portionlist_profile,max)))
      minval=min(unlist(lapply(portionlist_profile,min)))

      par(mar=c(4,4,2,2))
      plot(1, type="n", xlab="",xaxt="n", ylab=yl, xlim=c(1, length(portionlist_profile[[1]])), ylim=c(minval, maxval))
      axis(1,at=c(1,length(portionlist_profile[[1]])/2 +0.5,length(portionlist_profile[[1]])),labels=c("start","center","end"))

      for(i in 1:length(portionlist_profile)){
        lines(portionlist_profile[[i]],lwd=2,col=toplot$gadgetanalogic$color_distinct[i])
      }
      legend("topright",legend=names(portionlist_profile),col=toplot$gadgetanalogic$color_distinct,lty=rep(1,length(toplot$gadgetanalogic$color_distinct)),cex=0.6,bg="transparent")   

      dev.off()
  } 
)




#observer for PDF to download boxplot by ROI
output$saveboxplotByROIAnalogHeatbutton<- downloadHandler(
  filename=function() {
      paste('BoxplotByROI_analog_heatmap.pdf', sep='')
  },
  content=function(file) {
      pdf(file)

      factor_add=rep(0:(toplot$gadgetanalogic$roinumber-1),each=toplot$gadgetanalogic$bamnumber)
      addingfactor=1:length(toplot$gadgetanalogic$portionlist_boxes)+factor_add
  
      if(input$Log2BoxAnalogHeat){
        portionlist_boxes=lapply(toplot$gadgetanalogic$portionlist_boxes,log2)
        yl="Log2 read density (rpm/bp)"
      }else{
        portionlist_boxes=toplot$gadgetanalogic$portionlist_boxes
        yl="Read density (rpm/bp)"
      }

      isgrouped=input$GroupColorsAnalogHeat
      
      #if grouping colors, adjust legend and colors by groups
      if(!isgrouped){
        newcols=toplot$gadgetanalogic$color_distinct
        newnames=names(portionlist_boxes)
      }else{
        #in this case, take n° ROI colors, repeated for n° BAMs
        #in the legend, put n° ROI colors and tell which ROI with number
        #in xlab, put BAMs
        #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
        newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
        newcols=rep(toplot$gadgetanalogic$color_distinct[1:toplot$gadgetanalogic$bamnumber],toplot$gadgetanalogic$roinumber)
      }


      par(mar=c(14,4,1,1))
      suppressWarnings(boxplot(portionlist_boxes,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
      ats=c()
      for(i in 1:toplot$gadgetanalogic$roinumber){
        window=addingfactor[(((i-1)*toplot$gadgetanalogic$bamnumber)+1):(i*toplot$gadgetanalogic$bamnumber)]
        currentvalue=(window[length(window)]-window[1])/2
        currentvalue=window[1]+currentvalue
        ats=c(ats,currentvalue)
      }

      if(!isgrouped){
        axis(1,at=ats,label=toplot$gadgetanalogic$roiname,las=2)
        legend("topright",legend=newnames,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))          
      }else{
        axis(1,at=ats,label=newnames,las=2)
        legend("topright",legend=toplot$gadgetanalogic$bamname,col=newcols,cex=0.7,bg="transparent",pch=rep(19,length(portionlist_boxes)))
      }

      # axis(1,at=ats,label=toplot$gadgetanalogic$roiname,las=2)
      # legend("right",inset=c(-0.7,0),legend=names(portionlist_boxes),col=toplot$gadgetanalogic$color_distinct,cex=0.6,bg="transparent",pch=rep(19,length(portionlist_boxes)))

      dev.off()
  } 
)





#observer for boxplot by BAM PDF download button
output$saveboxplotByBAMAnalogHeatbutton<- downloadHandler(
  filename=function() {
      paste('BoxplotByBAM_analog_heatmap.pdf', sep='')
  },
  content=function(file) {
      pdf(file)
      
      #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
      #invert the order
      newlist=list()
      newcols=c()
      newnames=c()
      #if grouping colors, adjust legend and colors by groups
      if(!input$GroupColorsAnalogHeat){
        for(i in 1:toplot$gadgetanalogic$bamnumber){
          pos_inverted=seq(i,length(toplot$gadgetanalogic$portionlist_boxes),toplot$gadgetanalogic$bamnumber)
          newlist=c(newlist,toplot$gadgetanalogic$portionlist_boxes[pos_inverted])
          newcols=c(newcols,toplot$gadgetanalogic$color_distinct[pos_inverted])
          newnames=c(newnames,names(toplot$gadgetanalogic$portionlist_boxes)[pos_inverted])
        }          
      }else{
        #in this case, take n° ROI colors, repeated for n° BAMs
        #in the legend, put n° ROI colors and tell which ROI with number
        #in xlab, put BAMs
        for(i in 1:toplot$gadgetanalogic$bamnumber){
          pos_inverted=seq(i,length(toplot$gadgetanalogic$portionlist_boxes),toplot$gadgetanalogic$bamnumber)
          newlist=c(newlist,toplot$gadgetanalogic$portionlist_boxes[pos_inverted])
        }
        #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(toplot$gadgetanalogic$portionlist_boxes)))
        newnames=unique(sapply(strsplit(names(toplot$gadgetanalogic$portionlist_boxes),split=";"),"[[",1))
        newcols=rep(toplot$gadgetanalogic$color_distinct[1:toplot$gadgetanalogic$roinumber],toplot$gadgetanalogic$bamnumber)
      }


      factor_add=rep(0:(toplot$gadgetanalogic$bamnumber-1),each=toplot$gadgetanalogic$roinumber)
      addingfactor=1:length(newlist)+factor_add
      par(mar=c(14,4,1,1))

      if(input$Log2BoxAnalogHeat){
        newlist=lapply(newlist,log2)
        yl="Log2 read density (rpm/bp)"
      }else{
        yl="Read density (rpm/bp)"
      }

      suppressWarnings(boxplot(newlist,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
      ats=c()
      for(i in 1:toplot$gadgetanalogic$bamnumber){
        window=addingfactor[(((i-1)*toplot$gadgetanalogic$roinumber)+1):(i*toplot$gadgetanalogic$roinumber)]
        currentvalue=(window[length(window)]-window[1])/2
        currentvalue=window[1]+currentvalue
        ats=c(ats,currentvalue)
      }

      axis(1,at=ats,label=toplot$gadgetanalogic$bamname,las=2)
      legend("topright",legend=newnames,col=newcols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$gadgetanalogic$portionlist_boxes)))
   
      dev.off()
  } 
)




#observer for PDF button for correlation heatmap
output$savecorAnalogHeatbutton<- downloadHandler(
  filename=function() {
      paste('Correlation_analog_heatmap.pdf', sep='')
  },
  content=function(file) {
        pdf(file)

        par(mar=c(12,12,1,1),xpd=TRUE)
        image(0:nrow(toplot$gadgetanalogic$trasp_cor), 0:ncol(toplot$gadgetanalogic$trasp_cor),toplot$gadgetanalogic$trasp_cor[,ncol(toplot$gadgetanalogic$trasp_cor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=toplot$gadgetanalogic$my_palette,breaks=toplot$gadgetanalogic$brk)
        axis( 2, at=seq(0.5,ncol(toplot$gadgetanalogic$trasp_cor)+0.5-1,1 ), labels= rev(colnames( toplot$gadgetanalogic$trasp_cor )), las= 2 )
        axis( 1, at=seq(0.5,ncol(toplot$gadgetanalogic$trasp_cor)+0.5-1,1 ), labels= colnames( toplot$gadgetanalogic$trasp_cor ), las= 2 )
        for (x in (nrow(toplot$gadgetanalogic$correlation_total)-1+0.5):0.5  )
          for (y in 0.5: ((ncol(toplot$gadgetanalogic$correlation_total)-1+0.5)   ))
            text(y,x, round(toplot$gadgetanalogic$correlation_total[ncol(toplot$gadgetanalogic$correlation_total)-x+0.5,y+0.5],2),col="blue")

        dev.off()
  } 
)







#observer for PDF button of partial correlation heatmap of analog heatmap
output$savepcorAnalogHeatbutton<- downloadHandler(
  filename=function() {
      paste('PartialCorrelation_analog_heatmap.pdf', sep='')
  },
  content=function(file) {
        pdf(file)

        par(mar=c(12,12,1,1),xpd=TRUE)
        image(0:nrow(toplot$gadgetanalogic$trasp_pcor), 0:ncol(toplot$gadgetanalogic$trasp_pcor),toplot$gadgetanalogic$trasp_pcor[,ncol(toplot$gadgetanalogic$trasp_pcor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=toplot$gadgetanalogic$my_palette,breaks=toplot$gadgetanalogic$brk)
        axis( 2, at=seq(0.5,ncol(toplot$gadgetanalogic$trasp_pcor)+0.5-1,1 ), labels= rev(colnames( toplot$gadgetanalogic$trasp_pcor )), las= 2 )
        axis( 1, at=seq(0.5,ncol(toplot$gadgetanalogic$trasp_pcor)+0.5-1,1 ), labels= colnames( toplot$gadgetanalogic$trasp_pcor ), las= 2 )
        for (x in (nrow(toplot$gadgetanalogic$correlation_partial)-1+0.5):0.5  )
          for (y in 0.5: ((ncol(toplot$gadgetanalogic$correlation_partial)-1+0.5)   ))
            text(y,x, round(toplot$gadgetanalogic$correlation_partial[ncol(toplot$gadgetanalogic$correlation_partial)-x+0.5,y+0.5],2),col="blue")
     
        dev.off()
  } 
)




#observer for the button "new ROI from selection in analogic heatmap", will use:
# toplot$gadgetanalogic$ytop
# toplot$gadgetanalogic$ybottom

# NOT important! if you fnd the button, x;y are already checked by the observer of x and y!
# if(length(xright)>0 & length(ytop)>0 & length(xleft)>0 & length(ybottom)>0 & length(totbamsamples)>0 & !is.null(heatvariables$matlist)){
# if( (xright-xleft)>=0 & (ytop-ybottom)>0 & xright<=totbamsamples){

# check if toplot$analogic$finalrange and toplot$analogic$finalbam are not NULL:
# if(!is.null(toplot$analogic$finalrange) & !is.null(toplot$analogic$finalbam)){
#....and are not made by different ROIs (labels):
# lab=length(unique(toplot$analogic$finalrange$label))
# if length(lab)==1 ...

# check if name of new ROI already exists or is in blacklist ("promoters"....)
# if(.....)

#if all ok, create new ROI with all the elements... log will be "selected from [ybottom;ytop] from analogic heatmap"
#the new name is input$newROIfromAnalogHeat

observeEvent(input$confirmImportROIfromAnalogHeat,{
  #checks on y,x for selected area on heatmap should be useless...
  if(!is.null(toplot$analogic$finalrange) & !is.null(toplot$analogic$finalbam) & !is.null(toplot$analogic$finalfixes)){
    #find labels in toplot$analogic$finalrange: should be 1 type only. If more than 1 (multiple ROI selected)
    #do not allow the new ROI
    #now check whether roi name is not length==0 or already exist or prohibited name
    nottobe=c("promoters","transcripts","TES")
    if (nchar(input$newROIfromAnalogHeat)>=1 & !(any(input$newROIfromAnalogHeat == nottobe))) {
      nottobe2=c("promoters_genelist_","transcripts_genelist_","TES_genelist_")
      if(!grepl(nottobe2[1],input$newROIfromAnalogHeat) &!grepl(nottobe2[2],input$newROIfromAnalogHeat) & !grepl(nottobe2[3],input$newROIfromAnalogHeat)){
        #now check if ROI name already exists:
        nomi=unlist(lapply(ROIvariables$listROI,getName))
        if (!input$newROIfromAnalogHeat %in% nomi){
          #find which label (roi name). If found, add source name of that roi, subset selected from hatmap
          
          totrow=length(toplot$analogic$labelstosplit)
          start=totrow-toplot$gadgetanalogic$ytop
          stop=totrow-toplot$gadgetanalogic$ybottom

          if(start<=0){
            start=1
          }
          if(stop>totrow){
            stop=totrow
          }

          label_selected=toplot$analogic$labelstosplit[start:stop]
          #extract positions of labels: if only one label, continue
          if(length(unique(label_selected))==1){
            #select the ROI from which the selection was made (for the fix, source, type and so on)
            nomi=names(toplot$analogic$finalbam)
            pos=match(unique(label_selected),nomi)   

            #SOMETIMES it gives this error/warning:
            #Warning: Error in [: only 0's may be mixed with negative subscripts
            #76: observeEventHandler [servers/serverGENOMICS.R#3347]'
   
            #temporary patch: if start==0, start=1
            start=start+1
            stop=stop+1
            #permanent patch should consider ytop, ybottom, that now are correct coordinates...

            #calculate the new start and stop inside the block
            beforelabel=toplot$analogic$labelstosplit[1:(start-1)]
            newstart=table(beforelabel==unique(label_selected))["TRUE"]+1
            newstop=newstart+(stop-start)

            correct_range=toplot$analogic$finalrange[[pos]][newstart:newstop]
            correct_fix=toplot$analogic$finalfixes[[pos]][newstart:newstop]
            correct_bams=toplot$analogic$finalbam[[pos]]
            newBAMlist=as.list(rep(NA,length(correct_bams)))
            for(i in 1:length(newBAMlist)){
              newBAMlist[[i]]=correct_bams[[i]][newstart:newstop]
            }
            names(newBAMlist)=names(correct_bams)
            nomitot=unlist(lapply(ROIvariables$listROI,getName))
            postot=match(unique(label_selected),nomitot)

            tryCatch({
              originalROI=ROIvariables$listROI[[postot]]
              oldSource=getSource(originalROI)
              flag=getFlag(originalROI)
              #take only the first element of the fix, because we have one ROI
              
              toadd=paste("extracted ",length(correct_range)," ranges from analogic heatmap",sep="")
              newSource=c(oldSource,list(toadd))
              #we should have all the elements... create ROI!
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                              name=input$newROIfromAnalogHeat,
                                              range=correct_range,
                                              fixed=correct_fix,
                                              BAMlist=newBAMlist,
                                              flag=flag,
                                              source=newSource)
              logvariables$msg[[length(logvariables$msg)+1]]= paste('Created ',input$newROIfromAnalogHeat,' ROI, extracting ',length(correct_fix),' elements from analogic heatmap<br>',sep="")
              print(paste('Created ',input$newROIfromAnalogHeat,' ROI, extracting ',length(correct_fix),' elements from analogic heatmap',sep=""))
              sendSweetAlert(
                session = session,
                title = "ROI extracted!",
                text = paste('Created ',input$newROIfromAnalogHeat,' ROI, extracting ',length(correct_fix),' elements from analogic heatmap',sep=""),
                type = "success"
              ) 
            },
            warning = function( w ){
              print("warning: cannot retrieve the ROI, maybe lost some of the original ROIs...")
            },
            error = function( err ){
              print("warning: cannot retrieve the ROI, maybe lost some of the original ROIs...")
            })

          }else{
            #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Multiple ROIs selected from heatmap...<br></font>',sep="")
            sendSweetAlert(
              session = session,
              title = "Multiple ROIs selected",
              text = "The genomic intervals (rows in the heatmap) of a single ROI must be selected from the heatmap",
              type = "error"
            )           
          }
        }else{
          #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newROIfromAnalogHeat,' ROI already exists...<br></font>',sep="")
          sendSweetAlert(
            session = session,
            title = "Bad name",
            text = paste("\'",input$newROIfromAnalogHeat,"'\' ROI is already present. Choose another name",sep=""),
            type = "error"
          )          
        }
      }else{
        #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newROIfromAnalogHeat,' name is not allowed...<br></font>',sep="")
        sendSweetAlert(
          session = session,
          title = "Bad name",
          text = paste("\'",input$newROIfromAnalogHeat,"'\' name for the new ROI is not allowed (must not start with 'promoters_genelist*','transcripts_genelist*','TES_genelist*')",sep=""),
          type = "error"
        )        
      }  
    }else{
      #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">',input$newROIfromAnalogHeat,' name is not allowed or missing name...<br></font>',sep="")
      sendSweetAlert(
        session = session,
        title = "Bad name",
        text = paste("\'",input$newROIfromAnalogHeat,"'\' name for the new ROI is not allowed (must not be 'promoters','transcripts','TES') or is missing",sep=""),
        type = "error"
      )     
    }

  }else{
    #logvariables$msg[[length(logvariables$msg)+1]]= paste('<font color="red">Heatmap data not found...<br></font>',sep="")
    sendSweetAlert(
      session = session,
      title = "Data not found",
      text = "I didn't find the data from the heatmap",
      type = "error"
    )  
  }
},ignoreInit=TRUE)






##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
# respond to plot Profiles and Boxplots, of ROIs and BAMs selected
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################


###react to help buttons:
#parameters
observeEvent(input$msg_enrichmentInRois_parameters, {
  boxHelpServer(msg_enrichmentInRois_parameters)
})
#profiles
observeEvent(input$msg_enrichmentInRois_profiles, {
  boxHelpServer(msg_enrichmentInRois_profiles)
})
#boxplots
observeEvent(input$msg_enrichmentInRois_boxplots, {
  boxHelpServer(msg_enrichmentInRois_boxplots)
})
#correlations
observeEvent(input$msg_enrichmentInRois_correlations, {
  boxHelpServer(msg_enrichmentInRois_correlations)
})
#scatterplots
observeEvent(input$msg_enrichmentInRois_scatterplot, {
  boxHelpServer(msg_enrichmentInRois_scatterplot)
})





observeEvent(input$confirmUpdateProfilesAndBox,{
  set.seed(123)

  #clear coordinates for scatterplots
  toplot$profileAndBoxes$x=NULL
  toplot$profileAndBoxes$y=NULL
  binstouse=input$binsProfilesAndBox

  toplot$profileAndBoxes$ROIsForProfilesAndBox=input$ROIsForProfilesAndBox
  toplot$profileAndBoxes$binsProfilesAndBox=input$binsProfilesAndBox
  toplot$profileAndBoxes$BAMsForProfilesAndBox=input$BAMsForProfilesAndBox

  if (length(ROIvariables$listROI)>0 & length(input$ROIsForProfilesAndBox)>0 & 
      input$binsProfilesAndBox>0 & length(input$BAMsForProfilesAndBox)>0 &isvalid(input$binsProfilesAndBox)){
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForProfilesAndBox,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]
    maxbinstohave=min(unlist(lapply(roi,checkMaxBins)))

    toplot$profileAndBoxes$maxbinstohave=maxbinstohave
    
    if(input$binsProfilesAndBox<maxbinstohave){
      nameROI=unlist(lapply(roi,getName))
      bigbamlist=list()
      strands=as.list(rep(NA,length(roi)))
      for(i in 1:length(roi)){
        unifROI=uniqueROI(roi[[i]])
        strands[[i]]=as.factor(strand(getRange(uniqueROI(roi[[i]]))))
        allbams=names(getBAMlist(unifROI) )
        pos2=match(input$BAMsForProfilesAndBox,allbams)
        #order even if BAM reordered, should be fine
        bigbamlist[[i]]=getBAMlist(unifROI)[pos2]
        names(bigbamlist[[i]])=input$BAMsForProfilesAndBox
      }
      names(bigbamlist)=nameROI
      
      #from bigbamlist (1st order: roi, 2nd order: BAM) do the binning from baseCOverage to matrix


      matlists=list()
      if(input$chooseNormalizationProfilesAndBox=="readdensity"){
        Snorm_logic=TRUE
      }else{
        Snorm_logic=FALSE
      }
    
      #print("doing matrix from basecoverage")

      if(nc==1){
        matlists=lapply(1:length(bigbamlist),function(i) {
          matlist=lapply(1:length(bigbamlist[[i]]),function(k) {
                      return(makeMatrixFrombaseCoverage(bigbamlist[[i]][[k]],Nbins=binstouse,Snorm=Snorm_logic))
                    }) 
          names(matlist)=names(bigbamlist[[i]])
          # ###HERE invert "-" strand. They are in strands list. each element of this list
          # ###is a ROI. 
          pos_toinvert=as.character(strands[[i]])=="-"
          #
          for(k in 1:length(matlist)){
            #for each of the BAM file in this current ROI;
            matprovv=matlist[[k]]
            if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
              matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
            }
          }
          return(matlist)
        })
      }else{

        decision=0
        tryCatch({
          #if number of ROIs is > than number of BAMs, parallelize ROIs
          if(length(bigbamlist)>length(bigbamlist[[i]])){
            matlists=mclapply(1:length(bigbamlist),function(i) {
              matlist=lapply(1:length(bigbamlist[[i]]),function(k) {
                          return(makeMatrixFrombaseCoverage(bigbamlist[[i]][[k]],Nbins=binstouse,Snorm=Snorm_logic))
                        }) 
              names(matlist)=names(bigbamlist[[i]])

              # ###HERE invert "-" strand. They are in strands list. each element of this list
              # ###is a ROI. 
              pos_toinvert=as.character(strands[[i]])=="-"
              #
              for(k in 1:length(matlist)){
                #for each of the BAM file in this current ROI;
                matprovv=matlist[[k]]
                if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                  matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
                }
              }

              return(matlist)
            },mc.cores=nc)            
          
          }else{
            matlists=lapply(1:length(bigbamlist),function(i) {
              matlist=mclapply(1:length(bigbamlist[[i]]),function(k) {
                          return(makeMatrixFrombaseCoverage(bigbamlist[[i]][[k]],Nbins=binstouse,Snorm=Snorm_logic))
                        },mc.cores=nc) 
              names(matlist)=names(bigbamlist[[i]])

              # ###HERE invert "-" strand. They are in strands list. each element of this list
              # ###is a ROI. 
              pos_toinvert=as.character(strands[[i]])=="-"
              #
              for(k in 1:length(matlist)){
                #for each of the BAM file in this current ROI;
                matprovv=matlist[[k]]
                if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                  matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
                }
              }

              return(matlist)
            })               
          }

          #if multicore returned correctly, decision==1, following if won't be calculated..
          decision=1
        },
        warning = function( w ){
          print("warning: single core calculation, due to memory limit...")
        },
        error = function( err ){
          print("single core calculation, due to memory limit...")
        })
        if(decision==0){
          matlists=lapply(1:length(bigbamlist),function(i) {
            matlist=lapply(1:length(bigbamlist[[i]]),function(k) {
                        return(makeMatrixFrombaseCoverage(bigbamlist[[i]][[k]],Nbins=binstouse,Snorm=Snorm_logic))
                      }) 
            names(matlist)=names(bigbamlist[[i]])

            # ###HERE invert "-" strand. They are in strands list. each element of this list
            # ###is a ROI. 
            pos_toinvert=as.character(strands[[i]])=="-"
            #
            for(k in 1:length(matlist)){
              #for each of the BAM file in this current ROI;
              matprovv=matlist[[k]]
              if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
              }
            }

            return(matlist)
          })
        }
      }
   
      names(matlists)=nameROI

      #print("reordering lists")
      #obtain list form list of lists
      portionlist=list()
      listforbox=list()
      for(i in 1:length(matlists)){
        provv=matlists[[i]]
        provv_box=bigbamlist[[i]]
        names(provv)=names(provv_box)=paste("ROI:",names(matlists)[i]," (",nrow(matlists[[i]][[1]]),"); ",names(matlists[[i]]),sep="" )
        portionlist=c(portionlist,provv)
        listforbox=c(listforbox,provv_box)

      }

      #print("collapse matrixes")
      #apply mean on the matrix/median
      portionlist_profile=lapply(1:length(portionlist),function(i){
        mat=portionlist[[i]]
        mat_profile=apply(mat,2,mean)
      })
      #apply cumulative enrichment for boxplots.
      #need to use the original bamlist, because in binning we loose the remaining
      #of the division
      portionlist_boxes=lapply(1:length(listforbox),function(i){
        block=listforbox[[i]]
        len=sapply(block,length)
        mat_boxes=sapply(block,sum)
        if(input$chooseNormalizationProfilesAndBox=="readdensity"){
          mat_boxes=mat_boxes/len
        }  
        return(mat_boxes) 
      })
      names(portionlist_profile)=names(portionlist_boxes)=names(portionlist)

      n=length(portionlist_profile)
      #ROIvariables$colorsfordensity <- distinctColorPalette(n)
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      color_distinct <-sample(col_vector, n)

      if(input$chooseNormalizationProfilesAndBox=="readdensity"){
        yl="Read density (rpm/bp)"
        corvariables$is.density=TRUE
      }else{
        yl="Normalized reads (rpm)"
        corvariables$is.density=FALSE
      }

      #kep the values for the scatterplot, when you select cells in the correlation heatmap
      corvariables$portionlist_boxes=portionlist_boxes
      names(portionlist_profile)=names(portionlist_boxes)=names(portionlist)
      toplot$profileAndBoxes$Log2BoxProfilesAndBox=input$Log2BoxProfilesAndBox
      toplot$profileAndBoxes$portionlist_profile=portionlist_profile
      toplot$profileAndBoxes$portionlist_boxes=portionlist_boxes
      toplot$profileAndBoxes$yl=yl
      toplot$profileAndBoxes$color_distinct=color_distinct



      print("Drawing profiles and boxplots")
      ##PLOT profiles and boxplots using the obtained matrix lists
      output$profileProfilesAndBox<-renderPlot({
        if(input$Log2BoxProfilesAndBox){
          portionlist_profile=lapply(1:length(portionlist_profile),function(i) {log2(portionlist_profile[[i]])})
          yl=paste("Log2",yl)
        }else{
        }
        #

        maxval=max(unlist(lapply(portionlist_profile,max)))
        if(min(unlist(lapply(portionlist_profile,min))) <0){
          minval=min(unlist(lapply(portionlist_profile,min)))
        }else{
          minval=0
        }
        par(mar=c(10,4,1,1),xpd=TRUE)
        plot(1, type="n", xlab="",xaxt="n", ylab=yl, xlim=c(1, length(portionlist_profile[[1]])), ylim=c(minval, maxval))
        axis(1,at=c(1,length(portionlist_profile[[1]])/2 +0.5,length(portionlist_profile[[1]])),labels=c("start","center","end"))
        for(i in 1:length(portionlist_profile)){
          lines(portionlist_profile[[i]],lwd=2,col=color_distinct[i])
        }
        legend("bottom",inset=c(0,-0.5),legend=names(portionlist),col=color_distinct,lty=rep(1,length(color_distinct)),cex=0.6,bg="transparent")   
      })

      #PDF download button for profiles
      output$saveprofileProfilesAndBox=renderUI({downloadButton('saveprofileProfilesAndBoxbutton', 'Get PDF')})

      roinumber=length(matlists)
      bamnumber=length(matlists[[1]])

      toplot$profileAndBoxes$roinumber=roinumber
      toplot$profileAndBoxes$bamnumber=bamnumber

      toplot$profileAndBoxes$bamname=names(bigbamlist[[1]])
      toplot$profileAndBoxes$roiname=names(bigbamlist)



      #plot boxplots of regions selected. Use portionlist_boxes
      output$boxByBAMProfilesAndBox<-renderPlot({
        #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
        #invert the order
        newlist=list()
        newcols=c()
        newnames=c()

        #if grouping colors, adjust legend and colors by groups
        if(!input$GroupColorsProfilesAndBox){
          for(i in 1:bamnumber){
            pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
            newlist=c(newlist,portionlist_boxes[pos_inverted])
            newcols=c(newcols,color_distinct[pos_inverted])
            newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
          }          
        }else{
          #in this case, take n° ROI colors, repeated for n° BAMs
          #in the legend, put n° ROI colors and tell which ROI with number
          #in xlab, put BAMs
          for(i in 1:bamnumber){
            pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
            newlist=c(newlist,portionlist_boxes[pos_inverted])
          }
          #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
          newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
          newcols=rep(color_distinct[1:toplot$profileAndBoxes$roinumber],toplot$profileAndBoxes$bamnumber)
        }

        factor_add=rep(0:(bamnumber-1),each=roinumber)
        addingfactor=1:length(newlist)+factor_add
        par(mar=c(10,4,1,1),xpd=TRUE)
        if(input$Log2BoxProfilesAndBox){
          newlist=lapply(1:length(newlist),function(i) {log2(newlist[[i]])})
          yl=paste("Log2",yl)
        }else{
        }
        suppressWarnings(boxplot(newlist,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
        ats=c()
        for(i in 1:bamnumber){
          window=addingfactor[(((i-1)*roinumber)+1):(i*roinumber)]
          currentvalue=(window[length(window)]-window[1])/2
          currentvalue=window[1]+currentvalue
          ats=c(ats,currentvalue)
        }
        axis(1,at=ats,label=names(bigbamlist[[1]]),las=2)
        legend("bottom",inset=c(0,-0.5),legend=newnames,col=newcols,cex=0.6,bg="transparent",pch=rep(19,length(portionlist_boxes)))
      })
 


      output$boxByROIProfilesAndBox<-renderPlot({
        factor_add=rep(0:(roinumber-1),each=bamnumber)
        addingfactor=1:length(portionlist_boxes)+factor_add
        if(input$Log2BoxProfilesAndBox){
          portionlist_boxes=lapply(1:length(portionlist_boxes),function(i) {log2(portionlist_boxes[[i]])})
          yl=paste("Log2",yl)
        }else{
        }
        isgrouped=input$GroupColorsProfilesAndBox
        
        #if grouping colors, adjust legend and colors by groups
        if(!isgrouped){
          newcols=color_distinct
          newnames=names(portionlist)
        }else{
          #in this case, take n° ROI colors, repeated for n° BAMs
          #in the legend, put n° ROI colors and tell which ROI with number
          #in xlab, put BAMs
          #Warning: Error in strsplit: non-character argument
          #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
          newnames=unique(sapply(strsplit(names(portionlist),split=";"),"[[",1))
          #remove "ROI:" from newnames
          newnames=sapply(strsplit(newnames,split="ROI:"),"[[",2)
          newcols=rep(color_distinct[1:toplot$profileAndBoxes$bamnumber],toplot$profileAndBoxes$roinumber)
        }

        par(mar=c(10,4,1,1),xpd=TRUE)
        suppressWarnings(boxplot(portionlist_boxes,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE,las=2))
        ats=c()
        for(i in 1:roinumber){
          window=addingfactor[(((i-1)*bamnumber)+1):(i*bamnumber)]
          currentvalue=(window[length(window)]-window[1])/2
          currentvalue=window[1]+currentvalue
          ats=c(ats,currentvalue)
        }

        if(!isgrouped){
          axis(1,at=ats,label=toplot$profileAndBoxes$roiname,las=2)
          legend("bottom",inset=c(0,-0.5),legend=newnames,col=newcols,cex=0.6,bg="transparent",pch=rep(19,length(portionlist_boxes)))          
        }else{
          axis(1,at=ats,label=newnames,las=2)
          legend("bottom",inset=c(0,-0.5),legend=toplot$profileAndBoxes$bamname,col=newcols,cex=0.6,bg="transparent",pch=rep(19,length(portionlist_boxes)))
        }


      })  




      #PDF download buttons for boxplots by ROI and by BAM
      output$saveboxByROIProfilesAndBox=renderUI({downloadButton('saveboxByROIProfilesAndBoxbutton', 'Get PDF')})
      output$saveboxByBAMProfilesAndBox=renderUI({downloadButton('saveboxByBAMProfilesAndBoxbutton', 'Get PDF')})

      #download data of boxplot by ROI
      output$saveboxdataProfANDbox=renderUI({downloadButton('saveenrichmentBoxProfANDboxdata', 'Download data')})

      if(bamnumber>1 & roinumber==1){

        toplot$profileAndBoxes$choosecorMethodProfilesAndBox=input$choosecorMethodProfilesAndBox

        output$corProfilesAndBox<-renderPlot({
          if(input$Log2BoxProfilesAndBox){
            portionlist_boxes=lapply(portionlist_boxes,log2)
          }
          mat=do.call(cbind,portionlist_boxes)
          #if log2, 0 will be -Inf. Correct them
          mat[is.infinite(mat) &mat<0 ]=0
          #colnames(mat)=bamselected[xleft:xright]
          colnames(mat)=names(matlists[[1]])
          #according to the input, use pearson or spearman
          correlation_total=cor(mat,method=input$choosecorMethodProfilesAndBox)    
          trasp_cor=t(correlation_total)
          brk=c( seq( -1 , 1,0.01))
          my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

          par(mar=c(12,12,1,1),xpd=TRUE)
          image(0:nrow(trasp_cor), 0:ncol(trasp_cor),trasp_cor[,ncol(trasp_cor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
          axis( 2, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= rev(colnames( trasp_cor )), las= 2 )
          axis( 1, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= colnames( trasp_cor ), las= 2 )
          for (x in (nrow(correlation_total)-1+0.5):0.5  )
            for (y in 0.5: ((ncol(correlation_total)-1+0.5)   ))
              text(y,x, round(correlation_total[ncol(correlation_total)-x+0.5,y+0.5],2),col="blue")
        })
        
        #PDF download button for cor matrix
        output$savecorProfilesAndBox=renderUI({downloadButton('savecorProfilesAndBoxbutton', 'Get PDF')})

        #at least 3 for partial correlation, otherwise problems with lm
        if(bamnumber>2){
          output$pcorProfilesAndBox<-renderPlot({
	          if(input$Log2BoxProfilesAndBox){
	            portionlist_boxes=lapply(portionlist_boxes,log2)
	          }
	          mat=do.call(cbind,portionlist_boxes)
	          #colnames(mat)=bamselected[xleft:xright]
	          colnames(mat)=names(matlists[[1]])
	          mat[is.infinite(mat) &mat<0 ]=0
	          #accroding to the input, use pearson or spearman
	          correlation_partial=pcor(mat,method=input$choosecorMethodProfilesAndBox)$estimate
	          colnames(correlation_partial)=rownames(correlation_partial)= colnames(mat)
	          #correction for likely a bug in pcor function

	          trasp_pcor=t(correlation_partial)
	          brk=c( seq( -1 , 1,0.01))
	          my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

	          par(mar=c(12,12,1,1),xpd=TRUE)
	          image(0:nrow(trasp_pcor), 0:ncol(trasp_pcor),trasp_pcor[,ncol(trasp_pcor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
	          axis( 2, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= rev(colnames( trasp_pcor )), las= 2 )
	          axis( 1, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= colnames( trasp_pcor ), las= 2 )
	          for (x in (nrow(correlation_partial)-1+0.5):0.5  )
	            for (y in 0.5: ((ncol(correlation_partial)-1+0.5)   ))
	              text(y,x, round(correlation_partial[ncol(correlation_partial)-x+0.5,y+0.5],2),col="blue")
          })
  
          #PDf download button for partial correlation
          output$savepcorProfilesAndBox=renderUI({downloadButton('savepcorProfilesAndBoxbutton', 'Get PDF')})


        }else{
        	output$pcorProfilesAndBox <-renderPlot({NULL})
          output$savepcorProfilesAndBox=renderUI({NULL})
	        toplot$profileAndBoxes$px=NULL
	        toplot$profileAndBoxes$py=NULL
        }


      }else{
        output$corProfilesAndBox <-renderPlot({NULL})
        output$savecorProfilesAndBox=renderUI({NULL}) 
        output$pcorProfilesAndBox <-renderPlot({NULL})
        output$savepcorProfilesAndBox=renderUI({NULL})
        corvariables$portionlist_boxes=NULL
        toplot$profileAndBoxes$x=NULL
        toplot$profileAndBoxes$y=NULL
        toplot$profileAndBoxes$px=NULL
        toplot$profileAndBoxes$py=NULL
      }
    }else{
      output$profileProfilesAndBox<-renderPlot({NULL})
      output$boxByROIProfilesAndBox<-renderPlot({NULL})
      output$saveboxdataProfANDbox=renderUI({NULL})
      output$boxByBAMProfilesAndBox<-renderPlot({NULL})   
      output$corProfilesAndBox <-renderPlot({NULL})
      output$savecorProfilesAndBox=renderUI({NULL})    
      output$pcorProfilesAndBox <-renderPlot({NULL})
      output$savepcorProfilesAndBox=renderUI({NULL})
      output$saveprofileProfilesAndBox=renderUI({NULL})
      output$saveboxByROIProfilesAndBox=renderUI({NULL})
      output$saveboxByBAMProfilesAndBox=renderUI({NULL}) 
      corvariables$portionlist_boxes=NULL
      toplot$profileAndBoxes$x=NULL
      toplot$profileAndBoxes$y=NULL
      toplot$profileAndBoxes$px=NULL
      toplot$profileAndBoxes$py=NULL
      sendSweetAlert(
        session = session,
        title = "Too many bins",
        text = "The selected number of bins exceeds the maximum allowed, because some ranges have length < nbins. Try with a lower number",
        type = "error"
      )
    }
  }else{
    output$profileProfilesAndBox<-renderPlot({NULL})
    output$boxByROIProfilesAndBox<-renderPlot({NULL})
    output$saveboxdataProfANDbox=renderUI({NULL})
    output$boxByBAMProfilesAndBox<-renderPlot({NULL})
    output$saveprofileProfilesAndBox=renderUI({NULL})
    output$saveboxByROIProfilesAndBox=renderUI({NULL})
    output$saveboxByBAMProfilesAndBox=renderUI({NULL})
    output$corProfilesAndBox <-renderPlot({NULL})
    output$savecorProfilesAndBox=renderUI({NULL})    
    output$pcorProfilesAndBox <-renderPlot({NULL})
    output$savepcorProfilesAndBox=renderUI({NULL}) 
    corvariables$portionlist_boxes=NULL
    toplot$profileAndBoxes$x=NULL
    toplot$profileAndBoxes$y=NULL
    toplot$profileAndBoxes$px=NULL
    toplot$profileAndBoxes$py=NULL
    sendSweetAlert(
      session = session,
      title = "ROI or enrichment not selected",
      text = "You must provide at least one ROI and one enrichment associated",
      type = "error"
    )
  }
},ignoreInit=TRUE)





#Respond to PDF download button for profiles
output$saveprofileProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('Profiles.pdf', sep='')
  },
  content=function(file) {
        pdf(file)

        if(input$Log2BoxProfilesAndBox){
          portionlist_profile=lapply(1:length(toplot$profileAndBoxes$portionlist_profile),function(i) {log2(toplot$profileAndBoxes$portionlist_profile[[i]])})
          yl=paste("Log2",toplot$profileAndBoxes$yl)
        }else{
          yl=toplot$profileAndBoxes$yl
          portionlist_profile=toplot$profileAndBoxes$portionlist_profile
        }
        names(portionlist_profile)=names(toplot$profileAndBoxes$portionlist_profile)
        maxval=max(unlist(lapply(portionlist_profile,max)))
        if(min(unlist(lapply(portionlist_profile,min))) <0){
          minval=min(unlist(lapply(portionlist_profile,min)))
        }else{
          minval=0
        }
        par(mar=c(10,4,1,1),xpd=TRUE)
        plot(1, type="n", xlab="",xaxt="n", ylab=yl, xlim=c(1, length(portionlist_profile[[1]])), ylim=c(minval, maxval))
        axis(1,at=c(1,length(portionlist_profile[[1]])/2 +0.5,length(portionlist_profile[[1]])),labels=c("start","center","end"))
        for(i in 1:length(portionlist_profile)){
          lines(portionlist_profile[[i]],lwd=2,col=toplot$profileAndBoxes$color_distinct[i])
        }
        legend("bottom",inset=c(0,-0.25),legend=names(portionlist_profile),col=toplot$profileAndBoxes$color_distinct,lty=rep(1,length(toplot$profileAndBoxes$color_distinct)),cex=0.6,bg="transparent")   

        dev.off()
  } 
)





#PDF download button boxplot by ROI
output$saveboxByROIProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('BoxplotByROI.pdf', sep='')
  },
  content=function(file) {
        pdf(file)        

        factor_add=rep(0:(toplot$profileAndBoxes$roinumber-1),each=toplot$profileAndBoxes$bamnumber)
        addingfactor=1:length(toplot$profileAndBoxes$portionlist_boxes)+factor_add
        if(input$Log2BoxProfilesAndBox){
          portionlist_boxes=lapply(1:length(toplot$profileAndBoxes$portionlist_boxes),function(i) {log2(toplot$profileAndBoxes$portionlist_boxes[[i]])})
          yl=paste("Log2",toplot$profileAndBoxes$yl)
        }else{
          yl=toplot$profileAndBoxes$yl
          portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
        }

        isgrouped=input$GroupColorsProfilesAndBox
        #if grouping colors, adjust legend and colors by groups
        if(!isgrouped){
          newcols=toplot$profileAndBoxes$color_distinct
          newnames=names(toplot$profileAndBoxes$portionlist_boxes)
        }else{
          #in this case, take n° ROI colors, repeated for n° BAMs
          #in the legend, put n° ROI colors and tell which ROI with number
          #in xlab, put BAMs
          #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(portionlist_boxes)))
          newnames=unique(sapply(strsplit(names(toplot$profileAndBoxes$portionlist_boxes),split=";"),"[[",1))
          #remove "ROI:" from newnames
          newnames=sapply(strsplit(newnames,split="ROI:"),"[[",2)
          newcols=rep(toplot$profileAndBoxes$color_distinct[1:toplot$profileAndBoxes$bamnumber],toplot$profileAndBoxes$roinumber)
        }

        par(mar=c(15,4,1,1),xpd=TRUE)
        suppressWarnings(boxplot(portionlist_boxes,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE,las=2))
        ats=c()
        for(i in 1:toplot$profileAndBoxes$roinumber){
          window=addingfactor[(((i-1)*toplot$profileAndBoxes$bamnumber)+1):(i*toplot$profileAndBoxes$bamnumber)]
          currentvalue=(window[length(window)]-window[1])/2
          currentvalue=window[1]+currentvalue
          ats=c(ats,currentvalue)
        }
        if(!isgrouped){
          axis(1,at=ats,label=toplot$profileAndBoxes$roiname,las=2)
          legend("bottom",inset=c(0,-0.25),legend=newnames,col=newcols,cex=0.6,bg="transparent",pch=rep(19,length(portionlist_boxes)))          
        }else{
          axis(1,at=ats,label=newnames,las=2)
          legend("bottom",inset=c(0,-0.25),legend=toplot$profileAndBoxes$bamname,col=newcols,cex=0.6,bg="transparent",pch=rep(19,length(portionlist_boxes)))
        }
        dev.off()
  } 
)



#observer for download of xls table of the boxplot data (by ROI)
output$saveenrichmentBoxProfANDboxdata<- downloadHandler(
  filename=function() {
      paste('boxplot_data.xls', sep='')
  },
  content=function(file) {
      factor_add=rep(0:(toplot$profileAndBoxes$roinumber-1),each=toplot$profileAndBoxes$bamnumber)
      addingfactor=1:length(toplot$profileAndBoxes$portionlist_boxes)+factor_add
      if(input$Log2BoxProfilesAndBox){
        portionlist_boxes=lapply(1:length(toplot$profileAndBoxes$portionlist_boxes),function(i) {log2(toplot$profileAndBoxes$portionlist_boxes[[i]])})
      }else{
        portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
      }
      names(portionlist_boxes)=names(toplot$profileAndBoxes$portionlist_boxes)
      maxvalues=max(unlist(lapply(portionlist_boxes,length)))
      arr=matrix(rep("",maxvalues*length(portionlist_boxes)),ncol=length(portionlist_boxes))
      for(i in 1:length(portionlist_boxes)){
        arr[1:length(portionlist_boxes[[i]]),i]=portionlist_boxes[[i]]
      }
      colnames(arr)=names(portionlist_boxes)
      colnames(arr)=gsub(" ","_",colnames(arr))
      write.table(arr,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  } 
  
)





#PDF download boxplot by BAM
output$saveboxByBAMProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('BoxplotByBAM.pdf', sep='')
  },
  content=function(file) {
        pdf(file)
        #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
        #invert the order
        newlist=list()
        newcols=c()
        newnames=c()

        #if grouping colors, adjust legend and colors by groups
        if(!input$GroupColorsProfilesAndBox){
          for(i in 1:toplot$profileAndBoxes$bamnumber){
            pos_inverted=seq(i,length(toplot$profileAndBoxes$portionlist_boxes),toplot$profileAndBoxes$bamnumber)
            newlist=c(newlist,toplot$profileAndBoxes$portionlist_boxes[pos_inverted])
            newcols=c(newcols,toplot$profileAndBoxes$color_distinct[pos_inverted])
            newnames=c(newnames,names(toplot$profileAndBoxes$portionlist_boxes)[pos_inverted])
          }          
        }else{
          #in this case, take n° ROI colors, repeated for n° BAMs
          #in the legend, put n° ROI colors and tell which ROI with number
          #in xlab, put BAMs
          for(i in 1:toplot$profileAndBoxes$bamnumber){
            pos_inverted=seq(i,length(toplot$profileAndBoxes$portionlist_boxes),toplot$profileAndBoxes$bamnumber)
            newlist=c(newlist,toplot$profileAndBoxes$portionlist_boxes[pos_inverted])
          }
          #numbers_rois=unique(gsub(".*\\(|\\).*", "", names(toplot$profileAndBoxes$portionlist_boxes)))
          newnames=unique(sapply(strsplit(names(toplot$profileAndBoxes$portionlist_boxes),split=";"),"[[",1))
          newcols=rep(toplot$profileAndBoxes$color_distinct[1:toplot$profileAndBoxes$roinumber],toplot$profileAndBoxes$bamnumber)
        }

        factor_add=rep(0:(toplot$profileAndBoxes$bamnumber-1),each=toplot$profileAndBoxes$roinumber)
        addingfactor=1:length(newlist)+factor_add
        par(mar=c(15,4,1,1),xpd=TRUE)
        if(input$Log2BoxProfilesAndBox){
          newlist=lapply(1:length(newlist),function(i) {log2(newlist[[i]])})
          yl=paste("Log2",toplot$profileAndBoxes$yl)
        }else{
          yl=toplot$profileAndBoxes$yl
        }
        suppressWarnings(boxplot(newlist,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE,las=2))
        ats=c()
        for(i in 1:toplot$profileAndBoxes$bamnumber){
          window=addingfactor[(((i-1)*toplot$profileAndBoxes$roinumber)+1):(i*toplot$profileAndBoxes$roinumber)]
          currentvalue=(window[length(window)]-window[1])/2
          currentvalue=window[1]+currentvalue
          ats=c(ats,currentvalue)
        }
        axis(1,at=ats,label=toplot$profileAndBoxes$bamname,las=2)
        legend("bottom",inset=c(0,-0.25),legend=newnames,col=newcols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$profileAndBoxes$portionlist_boxes)))
        dev.off()
  } 
)





#PDF download button for correlation heatmap
output$savecorProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('Correlation_heatmap.pdf', sep='')
  },
  content=function(file) {
          pdf(file)

          if(input$Log2BoxProfilesAndBox){
            portionlist_boxes=lapply(toplot$profileAndBoxes$portionlist_boxes,log2)
          }else{
            portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
          }
          mat=do.call(cbind,portionlist_boxes)
          #if log2, 0 will be -Inf. Correct them
          mat[is.infinite(mat) &mat<0 ]=0
          #colnames(mat)=bamselected[xleft:xright]
          colnames(mat)=toplot$profileAndBoxes$bamname
          #accroding to the input, use pearson or spearman
          correlation_total=cor(mat,method=toplot$profileAndBoxes$choosecorMethodProfilesAndBox)    
          trasp_cor=t(correlation_total)
          brk=c( seq( -1 , 1,0.01))
          my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

          par(mar=c(12,12,1,1),xpd=TRUE)
          image(0:nrow(trasp_cor), 0:ncol(trasp_cor),trasp_cor[,ncol(trasp_cor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
          axis( 2, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= rev(colnames( trasp_cor )), las= 2 )
          axis( 1, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= colnames( trasp_cor ), las= 2 )
          for (x in (nrow(correlation_total)-1+0.5):0.5  )
            for (y in 0.5: ((ncol(correlation_total)-1+0.5)   ))
              text(y,x, round(correlation_total[ncol(correlation_total)-x+0.5,y+0.5],2),col="blue")

          dev.off()
  } 
)






#PDF download of partial correlation heatmap
output$savepcorProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('PartialCorrelation_heatmap.pdf', sep='')
  },
  content=function(file) {
          
          if(input$Log2BoxProfilesAndBox){
            portionlist_boxes=lapply(toplot$profileAndBoxes$portionlist_boxes,log2)
          }else{
            portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
          }
          mat=do.call(cbind,portionlist_boxes)
          #colnames(mat)=bamselected[xleft:xright]
          colnames(mat)=toplot$profileAndBoxes$bamname
          mat[is.infinite(mat) &mat<0 ]=0
          #accroding to the input, use pearson or spearman
          correlation_partial=pcor(mat,method=toplot$profileAndBoxes$choosecorMethodProfilesAndBox)$estimate
          trasp_pcor=t(correlation_partial)
          brk=c( seq( -1 , 1,0.01))
          my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

          pdf(file)
          par(mar=c(12,12,1,1),xpd=TRUE)
          image(0:nrow(trasp_pcor), 0:ncol(trasp_pcor),trasp_pcor[,ncol(trasp_pcor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
          axis( 2, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= rev(colnames( trasp_pcor )), las= 2 )
          axis( 1, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= colnames( trasp_pcor ), las= 2 )
          for (x in (nrow(correlation_partial)-1+0.5):0.5  )
            for (y in 0.5: ((ncol(correlation_partial)-1+0.5)   ))
              text(y,x, round(correlation_partial[ncol(correlation_partial)-x+0.5,y+0.5],2),col="blue")
        
          dev.off()
  } 
)





#PDF download button (reactive) for scatterplot
output$savescatterProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('Scatterplot_enrichments.pdf', sep='')
  },
  content=function(file) {
    pdf(file)

    my.cols <- rev(brewer.pal(5, "RdYlBu"))
    z <- kde2d(toplot$profileAndBoxes$array_x, toplot$profileAndBoxes$array_y, n=50)
    plot(toplot$profileAndBoxes$array_x,toplot$profileAndBoxes$array_y,xlab=toplot$profileAndBoxes$name_x,ylab=toplot$profileAndBoxes$name_y,main="Correlation (random 2000)")
    contour(z, drawlabels=FALSE, nlevels=5, col=my.cols, add=TRUE, lwd=2)
    legend("bottomright",legend=toplot$profileAndBoxes$corlegend)
    dev.off()
  } 
)







##################################################################
##################################################################
##################################################################
##################################################################
# respond to dynamics
##################################################################
##################################################################
##################################################################
##################################################################


###react to help buttons:
#parameters
observeEvent(input$msg_dynamicsOnGenes_parameters, {
  boxHelpServer(msg_dynamicsOnGenes_parameters)
})
#profiles
observeEvent(input$msg_dynamicsOnGenes_profiles, {
  boxHelpServer(msg_dynamicsOnGenes_profiles)
})
#TSS enrichment
observeEvent(input$msg_dynamicsOnGenes_TSSenrichment, {
  boxHelpServer(msg_dynamicsOnGenes_TSSenrichment)
})
#genebodies enrichment
observeEvent(input$msg_dynamicsOnGenes_genebodiesEnrichment, {
  boxHelpServer(msg_dynamicsOnGenes_genebodiesEnrichment)
})
#TES enrichment
observeEvent(input$msg_dynamicsOnGenes_TESenrichment, {
  boxHelpServer(msg_dynamicsOnGenes_TESenrichment)
})
#TSS ranked
observeEvent(input$msg_dynamicsOnGenes_TSSranked, {
  boxHelpServer(msg_dynamicsOnGenes_TSSranked)
})
#Genebodies ranked enrichments
observeEvent(input$msg_dynamicsOnGenes_genebodiesRanked, {
  boxHelpServer(msg_dynamicsOnGenes_genebodiesRanked)
})
#stalling Index
observeEvent(input$msg_dynamicsOnGenes_stallingIndexRanked, {
  boxHelpServer(msg_dynamicsOnGenes_stallingIndexRanked)
})


observeEvent(input$plotDynamics,{
  #check length ROIvariables$listROI (>0),
  #then read input$genelistsforDynamics
  #then read input$BAMforDynamics
  #at this point, genelist slected has al the triplet p,t,t and if BAM selected, both
  #have the BAM. 
  if(length(ROIvariables$listROI)>0 & length(input$genelistsforDynamics)>0 & 
      length(input$BAMforDynamics)>0 & isvalid(input$binsforDynamics)){
    set.seed(123)
    toplot$dynamics$roinames=input$genelistsforDynamics
    toplot$dynamics$bamnames=input$BAMforDynamics
    toplot$dynamics$nbin=input$binsforDynamics
    #at this point, at least one valid genelist selected and at least one BAM file
    #select promoters, transcripts and TES of genelists selected

    if(input$chooseNormalizationforDynamics=="totread"){
      Snormalization=FALSE
      ylabel="Normalized reads (rpm)"
    }else{
      Snormalization=TRUE
      ylabel="Read density (rpm/bp)"
    }
    nomi=unlist(lapply(ROIvariables$listROI,getName))

    totallist_boxplot=rep(list(NA),length(input$BAMforDynamics))
    names(totallist_boxplot)=input$BAMforDynamics
    totallist_boxplot=rep(list(totallist_boxplot),length(input$genelistsforDynamics))
    names(totallist_boxplot)=input$genelistsforDynamics 
    totallist_boxplot=rep(list(totallist_boxplot),3)
    names(totallist_boxplot)=c("TSS","GB","TES")

    totallist_SI=totallist_boxplot
    names(totallist_SI)=c("TSS","GB","SI")

    totallist_profile=rep(list(NA),length(input$BAMforDynamics))
    names(totallist_profile)=input$BAMforDynamics
    totallist_profile=rep(list(totallist_profile),length(input$genelistsforDynamics))
    names(totallist_profile)=input$genelistsforDynamics        
    #TO BE PARALLELIZED IF NC>1
    lengths_rois=rep(NA,length(input$genelistsforDynamics))
    lengths_uniquegeneIDs=rep(NA,length(input$genelistsforDynamics))

    for(i in 1:length(input$genelistsforDynamics)){
      #change this if "genelist " is before the actual name of the genelist
      #current=strsplit(input$genelistsforDynamics[i],split=" ")[[1]][2]
      current=input$genelistsforDynamics[i]

      searchname_promoters=paste("promoters_genelist_",current,sep="")
      searchname_transcripts=paste("transcripts_genelist_",current,sep="")
      searchname_TES=paste("TES_genelist_",current,sep="")
      #find promoters, transcripts, TS for that ROI, unify the strand for the + strand
      pos_promoters= which(searchname_promoters==nomi)
      roi_promoters= unifyStrand(ROIvariables$listROI[[pos_promoters]])
      pos_transcripts= which(searchname_transcripts==nomi)
      #not unify strand. Invert matrixes of negative strands later, it's faster.
      roi_transcripts= ROIvariables$listROI[[pos_transcripts]]
      pos_TES= which(searchname_TES==nomi)
      roi_TES= unifyStrand(ROIvariables$listROI[[pos_TES]])   
      #get BAM for each roi of that genelist 
      BAM_promoters=names(getBAMlist(roi_promoters))
      BAM_transcripts=names(getBAMlist(roi_transcripts))      
      BAM_TES=names(getBAMlist(roi_TES)) 

      for(k in 1:length(input$BAMforDynamics)){
        #extract baseCoverage for promoters/transcripts/TES for that gene list   
        BAMnametosearch=input$BAMforDynamics[k]
        pos_BAM_promoters=which(BAMnametosearch==BAM_promoters)
        pos_BAM_transcripts=which(BAMnametosearch==BAM_transcripts)
        pos_BAM_TES=which(BAMnametosearch==BAM_TES)
        content_promoters=getBAMlist(roi_promoters)[[pos_BAM_promoters]]
        content_transcripts=getBAMlist(roi_transcripts)[[pos_BAM_transcripts]]
        content_TES=getBAMlist(roi_TES)[[pos_BAM_TES]]
        #for this genelist (p,t,t) and BAM, for promoter, transcripts and TES, calculate
        #profiles and box on the promoters,gb,TES
        #calculate content_transcripts_mod for gb only, where you cut 30% AND TSS down and TES up
        #to remove 30%, thirtypercentToRemove= l*15/80, where l is the length already increased by 30%

        #PROBLEM: if no elements on plus strand?!
        range_transcripts=getRange(roi_transcripts)
        smalllength=width(range_transcripts)
        biglength=unlist(lapply(content_transcripts,length))
        thirtypercent=(biglength-smalllength)/2

        #repeated for each BAM, but the important thing is that is the same i (ROI)
        #required for the number of ranges each genelist have (for the legend)
        lengths_rois[i]=length(range_transcripts)
        lengths_uniquegeneIDs[i]=length(unique(as.character(range_transcripts$gene_id)))

        #now cut TSS down
        #find on + strand, the fix-end(TSS)
        range_TSS=getRange(roi_promoters)
        plusstrand= as.logical(strand(range_TSS)=="+")
        part_TSS=unique(end(range_TSS[plusstrand])- start(getFixed(roi_promoters))[plusstrand])
        #cut TES up
        #find on + strans, the start(TES)-fix
        range_TES=getRange(roi_TES)
        plusstrand=as.logical(strand(range_TES)=="+")
        part_TES=unique( start(getFixed(roi_TES))[plusstrand]-start(range_TES[plusstrand]))

        #remove transcripts, promoters, TES whose length(GB)< smalllength-(part_TSS+part_TES)
        tokeep=smalllength-(part_TSS+part_TES) > input$binsforDynamics
        smalllength=smalllength[tokeep]
        content_transcripts=content_transcripts[tokeep]
        biglength=biglength[tokeep]
        thirtypercent=thirtypercent[tokeep]
        range_transcripts=range_transcripts[tokeep]
        pos_plus=as.character(strand(range_transcripts))!="-"
        pos_minus=as.character(strand(range_transcripts))=="-"
        thirtypercent_plus=thirtypercent[pos_plus]
        thirtypercent_minus=thirtypercent[pos_minus]
        content_transcripts_plus=content_transcripts[pos_plus]
        content_transcripts_minus=content_transcripts[pos_minus]
        smalllength_plus=smalllength[pos_plus]
        smalllength_minus=smalllength[pos_minus]
        smalllength=c(smalllength_plus,smalllength_minus)

        widthTES_reorder=width(getRange(roi_TES)[tokeep])
        widthTES_reorder_plus=widthTES_reorder[pos_plus]
        widthTES_reorder_minus=widthTES_reorder[pos_minus]
        widthTES_reorder=c(widthTES_reorder_plus,widthTES_reorder_minus)

        #remove the 30% up/downstream and part of TSS and TES
        start_pos_plus=thirtypercent_plus+1+part_TSS
        end_pos_plus=unlist(lapply(content_transcripts_plus,length))-thirtypercent_plus-part_TES

        start_pos_minus=thirtypercent_minus+1+part_TES
        end_pos_minus=unlist(lapply(content_transcripts_minus,length))-thirtypercent_minus-part_TSS
        
        #PROFILES
        
        #calc. matrixes for profiles, and invert negative strand
        matcov_plus=makeMatrixFrombaseCoverage(content_transcripts_plus,Nbins=toplot$dynamics$nbin,Snorm=Snormalization)
        matcov_minus=makeMatrixFrombaseCoverage(content_transcripts_minus,Nbins=toplot$dynamics$nbin,Snorm=Snormalization)
        matcov_minus <- matcov_minus[ , ncol(matcov_minus):1]
        matcov=rbind(matcov_plus,matcov_minus)
        #calculate median instead of mean
        totallist_profile[[i]][[k]]=apply(matcov,2,input$chooseMetricforDynamics)

        #here we sum the part of the transcripts that will be used:
        #    TSS      transcript    TES
        #.....|.....-------------....|....
        #-30%   TSS part      TES part  +30%
        #we remove the 30% up and down from original BAM files from transcripts, plus the part
        #of the TSS and TES that will not be considered as "transcript". We keep only the internal part...

        #to improve efficiency, use the custom Rcpp function

        content_transcripts_mod_plus=cutAndSumTranscripts(GRbaseCoverageOutput=content_transcripts_plus,
                                                      StartingPositions=start_pos_plus,
                                                      EndingPositions=end_pos_plus)
        content_transcripts_mod_minus=cutAndSumTranscripts(GRbaseCoverageOutput=content_transcripts_minus,
                                                      StartingPositions=start_pos_minus,
                                                      EndingPositions=end_pos_minus)

        #join mod plus and minus
        content_transcripts_mod=c(content_transcripts_mod_plus,content_transcripts_mod_minus)

        #BOXPLOTS
        #for boxplots, sum everything inside TSS,TES and the part of the transcripts cut before

        totallist_boxplot[[1]][[i]][[k]]=sapply(content_promoters[tokeep],sum)
        totallist_boxplot[[2]][[i]][[k]]=content_transcripts_mod
        totallist_boxplot[[3]][[i]][[k]]=sapply(content_TES[tokeep],sum)

        #STALLING INDEX
        #use content_transcripts_mod for gb previously calculated.
        #it depends on total reads, so the division of boxplot list must be carried out later
        if(input$chooseNormalizationforDynamics=="totread"){
          totallist_SI[[2]][[i]][[k]]=totallist_boxplot[[2]][[i]][[k]]+totallist_boxplot[[3]][[i]][[k]]#calculate with content_transcripts_mod, + totallist_boxplot[[3]][[i]][[k]] (TES part)
        #else, divide by the GB + TES interval:
        # TSS      GB           TES
        #....|-------------|-----.-----|
        }else{
          totallist_SI[[2]][[i]][[k]]=(totallist_boxplot[[2]][[i]][[k]]+totallist_boxplot[[3]][[i]][[k]])/ ( (smalllength-(part_TSS+part_TES))+ widthTES_reorder     )
        }

        #if read density, divide by length
        if(input$chooseNormalizationforDynamics!="totread"){
          totallist_boxplot[[1]][[i]][[k]]=totallist_boxplot[[1]][[i]][[k]]/width(getRange(roi_promoters)[tokeep])
          totallist_boxplot[[2]][[i]][[k]]=totallist_boxplot[[2]][[i]][[k]]/(smalllength-(part_TSS+part_TES)) # reduced block
          totallist_boxplot[[3]][[i]][[k]]=totallist_boxplot[[3]][[i]][[k]]/width(getRange(roi_TES)[tokeep])     
        }

        totallist_SI[[1]][[i]][[k]]=totallist_boxplot[[1]][[i]][[k]]
        
        #calculate good pseudocount: one tenth of the first percentile of non-zero totallist_SI[[2]] values
        pseudocount=totallist_SI[[2]][[i]][[k]]
        pseudocount=pseudocount[pseudocount!=0]
        pseudocount=quantile(pseudocount,0.0001)/10
        totallist_SI[[3]][[i]][[k]]=totallist_SI[[1]][[i]][[k]]/(totallist_SI[[2]][[i]][[k]]+pseudocount)#calculate stalling index

      }
    }

    randomcolors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    cols=sample(randomcolors,length(input$genelistsforDynamics)*length(input$BAMforDynamics))
    toplot$dynamics$cols=cols
    toplot$dynamics$ylabel=ylabel
    toplot$dynamics$totallist_profile=totallist_profile
    toplot$dynamics$totallist_boxplot=totallist_boxplot
    toplot$dynamics$totallist_SI=totallist_SI
    totalnames=c()
    for(i in 1:length(totallist_profile)){
      provv=names(totallist_profile[[i]])
      #add the number for each genelist (lengths_rois in position i)
      provv=paste(names(totallist_profile)[i],"(",lengths_rois[i],";",lengths_uniquegeneIDs[i],"unique ID); ",provv,sep="")
      #remove "genelist " from names
      # for(k in 1:length(provv)){
      #   provv[k]=substr(provv[k],10,nchar(provv[k]))
      # }
      totalnames=c(totalnames,provv)
    }
    toplot$dynamics$totalnames=totalnames
    #plot the profiles, boxplots, stalling index


    print("Plotting dynamics on genes")
    #profiles
    output$plotProfileDynamics<-renderPlot({
      #find max and min of the profiles 
      mins=c()
      maxs=c()
      for(i in 1:length(totallist_profile)){
        for(k in 1:length(totallist_profile[[i]])){
          if(input$islogforDynamics){
            maxs=c(maxs,max( log2(totallist_profile[[i]][[k]])))
            mins=c(mins,min( log2(totallist_profile[[i]][[k]])[!is.infinite(log2(totallist_profile[[i]][[k]]))]  ))
          }else{
            maxs=c(maxs,max(totallist_profile[[i]][[k]]))
            mins=c(mins,min(totallist_profile[[i]][[k]]))            
          }

        }
      }
      maxvalue=max(maxs)
      minvalue=min(mins[!is.infinite(mins)])
      if(input$islogforDynamics){
        ylabel=paste("log2",ylabel)
      }

      par(mar=c(5,5,1,4))
      plot(1, type="n", xlab="",xaxt="n", ylab=ylabel, xlim=c(1, toplot$dynamics$nbin), ylim=c(minvalue, maxvalue))
      axis(1,at=c(1,toplot$dynamics$nbin/2 +0.5,toplot$dynamics$nbin),labels=c("TSS-30%","genes","TES+30%"))
      count=1
      for(i in 1:length(totallist_profile)){
        for(k in 1:length(totallist_profile[[i]])){
          if(input$islogforDynamics){
            lines( log2(totallist_profile[[i]][[k]]),lwd=2,col= toplot$dynamics$cols[count],pch=".",cex=2,xaxt="n")
          }else{
            lines(totallist_profile[[i]][[k]],lwd=2,col= toplot$dynamics$cols[count],pch=".",cex=2,xaxt="n")
          }
          count=count+1
        }
      }

      legend("topright",legend=totalnames,col= toplot$dynamics$cols,lty=rep(1,length(totalnames)),bg="transparent")  

    })


    #boxplots
    output$plotboxTSSDynamics<-renderPlot({
      currentlist=totallist_boxplot[[1]]
      list2=list()
      #transform 2-order list in 1 order list
      for(i in 1:length(currentlist)){
        list2=c(list2,currentlist[[i]])
      }

      #HERE put the log2 of input$islogforDynamics
      #if Infinite values, those won't be drown
      if(input$islogforDynamics){
        for(i in 1:length(list2)){
          list2[[i]]=log2(list2[[i]])
        }
      }

      #boxplot list2. space at each ROI
      numroi=length(currentlist)
      numbam=length(currentlist[[1]])
      factor_add=rep(0:(numroi-1),each=numbam)
      addingfactor=1:length(list2)+factor_add
      ats=c()
      for(i in 1:numroi){
        window=addingfactor[(((i-1)*numbam)+1):(i*numbam)]
        currentvalue=(window[length(window)]-window[1])/2
        currentvalue=window[1]+currentvalue
        ats=c(ats,currentvalue)
      }

      if(input$islogforDynamics){
        ylabel=paste("log2",ylabel)
      }
      par(mar=c(10,4,1,1),xpd=TRUE)
      suppressWarnings(boxplot(list2,at=addingfactor,col=toplot$dynamics$cols,ylab=ylabel,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
      axis(1,at=ats,label=toplot$dynamics$roinames)
      legend("bottom",inset=c(0,-0.4),legend=totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(totalnames)))

    })

    #download data of boxplot of TSS
    output$saveboxdatadynamicsTSS=renderUI({downloadButton('saveenrichmentBoxdynamicsTSSdata', 'Download data')})

    output$plotboxGBDynamics<-renderPlot({
      currentlist=totallist_boxplot[[2]]
      list2=list()
      #transform 2-order list in 1 order list
      for(i in 1:length(currentlist)){
        list2=c(list2,currentlist[[i]])
      }


      #HERE put the log2 of input$islogforDynamics
      #if Infinite values, those won't be drown
      if(input$islogforDynamics){
        for(i in 1:length(list2)){
          list2[[i]]=log2(list2[[i]])
        }
      }
      #boxplot list2. space at each ROI
      numroi=length(currentlist)
      numbam=length(currentlist[[1]])
      factor_add=rep(0:(numroi-1),each=numbam)
      addingfactor=1:length(list2)+factor_add
      ats=c()
      for(i in 1:numroi){
        window=addingfactor[(((i-1)*numbam)+1):(i*numbam)]
        currentvalue=(window[length(window)]-window[1])/2
        currentvalue=window[1]+currentvalue
        ats=c(ats,currentvalue)
      }
      if(input$islogforDynamics){
        ylabel=paste("log2",ylabel)
      }
      par(mar=c(10,4,1,1),xpd=TRUE)
      suppressWarnings(boxplot(list2,at=addingfactor,col=toplot$dynamics$cols,ylab=ylabel,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
      axis(1,at=ats,label=toplot$dynamics$roinames)

      legend("bottom",inset=c(0,-0.4),legend=totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(totalnames)))
    })

    #download data of boxplot of GB
    output$saveboxdatadynamicsGB=renderUI({downloadButton('saveenrichmentBoxdynamicsGBdata', 'Download data')})

    output$plotboxTESDynamics<-renderPlot({
      currentlist=totallist_boxplot[[3]]
      list2=list()
      #transform 2-order list in 1 order list
      for(i in 1:length(currentlist)){
        list2=c(list2,currentlist[[i]])
      }

      #HERE put the log2 of input$islogforDynamics
      #if Infinite values, those won't be drown
      if(input$islogforDynamics){
        for(i in 1:length(list2)){
          list2[[i]]=log2(list2[[i]])
        }
      }

      #boxplot list2. space at each ROI
      numroi=length(currentlist)
      numbam=length(currentlist[[1]])
      factor_add=rep(0:(numroi-1),each=numbam)
      addingfactor=1:length(list2)+factor_add
      ats=c()
      for(i in 1:numroi){
        window=addingfactor[(((i-1)*numbam)+1):(i*numbam)]
        currentvalue=(window[length(window)]-window[1])/2
        currentvalue=window[1]+currentvalue
        ats=c(ats,currentvalue)
      }
      if(input$islogforDynamics){
        ylabel=paste("log2",ylabel)
      }
      par(mar=c(10,4,1,1),xpd=TRUE)
      suppressWarnings(boxplot(list2,at=addingfactor,col=toplot$dynamics$cols,ylab=ylabel,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
      axis(1,at=ats,label=toplot$dynamics$roinames)

      legend("bottom",inset=c(0,-0.4),legend=totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(totalnames)))
    })
    #download data of boxplot of TES
    output$saveboxdatadynamicsTES=renderUI({downloadButton('saveenrichmentBoxdynamicsTESdata', 'Download data')})
    
    output$plotSITSSDynamics<-renderPlot({
      currentlist=totallist_SI[[1]]
      totalnames=c()
      for(i in 1:length(currentlist)){
        provv=names(currentlist[[i]])
        provv=paste(names(currentlist)[i],provv)
        totalnames=c(totalnames,provv)
      }
      list2=list()
      #transform 2-order list in 1 order list
      minvals=c()
      maxvals=c()
      for(i in 1:length(currentlist)){
        list2=c(list2,currentlist[[i]])
      }



      outlayer_thresh=input$percentageOutlayerCumulPlots
      for(i in 1:length(list2)){
        mincurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],outlayer_thresh )
        maxcurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],1-outlayer_thresh )
        minvals=c(minvals,mincurrent)
        maxvals=c(maxvals,maxcurrent)
      }


      
      maxvals=max(maxvals)
      minvals=min(minvals)

      lengths=sapply(list2,length)
      par(mar=c(10,4,1,1))
      plot(1, type="n", ylab="Ranked Genes", xlab=paste("log2",ylabel), ylim=c(0, quantile(1:max(lengths),1-(2*outlayer_thresh))), xlim=c(minvals, maxvals))
      for(i in 1:length(list2)){
        #remove a fraction of outlayers
        #celan -Inf
        vals=sort(log2(list2[[i]]))
        vals=vals[!is.infinite(vals)]
        vals=vals[vals<quantile(vals,1-outlayer_thresh) & vals>quantile(vals,outlayer_thresh)]
        lines(vals,1:length(vals),col=toplot$dynamics$cols[i],lwd=2)
      }


      legend("bottom",inset=c(0,-0.4),legend=totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(totalnames)))

    })

    output$plotSIGBDynamics<-renderPlot({
      currentlist=totallist_SI[[2]]
      totalnames=c()
      for(i in 1:length(currentlist)){
        provv=names(currentlist[[i]])
        provv=paste(names(currentlist)[i],provv)
        totalnames=c(totalnames,provv)
      }
      list2=list()
      #transform 2-order list in 1 order list
      minvals=c()
      maxvals=c()
      for(i in 1:length(currentlist)){
        list2=c(list2,currentlist[[i]])
      }



      outlayer_thresh=input$percentageOutlayerCumulPlots
      for(i in 1:length(list2)){
      	mincurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],outlayer_thresh )
      	maxcurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],1-outlayer_thresh )
        minvals=c(minvals,mincurrent)
        maxvals=c(maxvals,maxcurrent)
      }
      
      maxvals=max(maxvals)
      minvals=min(minvals)

      lengths=sapply(list2,length)
      
      par(mar=c(10,4,1,1))
      plot(1, type="n", ylab="Ranked Genes", xlab=paste("log2",ylabel), ylim=c(0, quantile(1:max(lengths),1-(2*outlayer_thresh))), xlim=c(minvals, maxvals))
      for(i in 1:length(list2)){
      	#remove a fraction of outlayers
      	#celan -Inf
      	vals=sort(log2(list2[[i]]))
      	vals=vals[!is.infinite(vals)]
      	vals=vals[vals<quantile(vals,1-outlayer_thresh) & vals>quantile(vals,outlayer_thresh)]
        lines(vals,1:length(vals),col=toplot$dynamics$cols[i],lwd=2)
      }

      legend("bottom",inset=c(0,-0.4),legend=totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(totalnames)))

    })





    output$plotSISIDynamics<-renderPlot({
      currentlist=totallist_SI[[3]]
      totalnames=c()
      for(i in 1:length(currentlist)){
        provv=names(currentlist[[i]])
        provv=paste(names(currentlist)[i],provv)
        totalnames=c(totalnames,provv)
      }
      list2=list()
      #transform 2-order list in 1 order list
      minvals=c()
      maxvals=c()
      for(i in 1:length(currentlist)){
        list2=c(list2,currentlist[[i]])
      }


      outlayer_thresh=input$percentageOutlayerCumulPlots
      for(i in 1:length(list2)){
        mincurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],outlayer_thresh )
        maxcurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],1-outlayer_thresh )
        minvals=c(minvals,mincurrent)
        maxvals=c(maxvals,maxcurrent)
      }

 
      maxvals=max(maxvals)
      minvals=min(minvals)
      
      lengths=sapply(list2,length)

      par(mar=c(10,4,1,1))
      plot(1, type="n", ylab="Ranked Genes", xlab="stalling index", ylim=c(0, quantile(1:max(lengths),1-(2*outlayer_thresh))), xlim=c(minvals, maxvals))
      for(i in 1:length(list2)){
        #remove a fraction of outlayers
        #celan -Inf
        vals=sort(log2(list2[[i]]))
        vals=vals[!is.infinite(vals)]
        vals=vals[vals<quantile(vals,1-outlayer_thresh) & vals>quantile(vals,outlayer_thresh)]
        lines(vals,1:length(vals),col=toplot$dynamics$cols[i],lwd=2)
      }
      legend("bottom",inset=c(0,-0.4),legend=totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(totalnames)))
    })

    #buttons for downloads of PDF
    output$saveprofileDynamics=renderUI({downloadButton('saveprofileDynamicsbutton', 'Get PDF')})
    output$saveboxTSSDynamics=renderUI({downloadButton('saveboxTSSDynamicsbutton', 'Get PDF')})
    output$saveboxGBDynamics=renderUI({downloadButton('saveboxGBDynamicsbutton', 'Get PDF')})
    output$saveboxTESDynamics=renderUI({downloadButton('saveboxTESDynamicsbutton', 'Get PDF')})
    output$saveSITSSDynamics=renderUI({downloadButton('saveSITSSDynamicsbutton', 'Get PDF')})
    output$saveSIGBDynamics=renderUI({downloadButton('saveSIGBDynamicsbutton', 'Get PDF')})
    output$saveSISIDynamics=renderUI({downloadButton('saveSISIDynamicsbutton', 'Get PDF')})
  }else{
    #no gene list selected and/or no BAM available/selected
    output$plotProfileDynamics<-renderPlot({NULL})
    output$plotboxTSSDynamics<-renderPlot({NULL})
    output$plotboxGBDynamics<-renderPlot({NULL})
    output$plotboxTESDynamics<-renderPlot({NULL})
    output$plotSITSSDynamics<-renderPlot({NULL})
    output$plotSIGBDynamics<-renderPlot({NULL})
    output$plotSISIDynamics<-renderPlot({NULL})
    output$saveprofileDynamics=renderUI({NULL})
    output$saveboxTSSDynamics=renderUI({NULL})
    output$saveboxGBDynamics=renderUI({NULL})
    output$saveboxTESDynamics=renderUI({NULL})
    output$saveSITSSDynamics=renderUI({NULL})
    output$saveSIGBDynamics=renderUI({NULL})
    output$saveSISIDynamics=renderUI({NULL})
    sendSweetAlert(
      session = session,
      title = "Genelist or enrichment not selected",
      text = "You must provide at least one genelist (promoters+transcripts+TES) and one enrichment associated",
      type = "error"
    )
  }

},ignoreInit=TRUE)






#download PDF of profile dynamics
output$saveprofileDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('Profile_dynamics.pdf', sep='')
  },
  content=function(file) {
    pdf(file,width=15,height=8)

    #find max and min of the profiles 
    mins=c()
    maxs=c()
    for(i in 1:length(toplot$dynamics$totallist_profile)){
      for(k in 1:length(toplot$dynamics$totallist_profile[[i]])){
        if(input$islogforDynamics){
          maxs=c(maxs,max( log2(toplot$dynamics$totallist_profile[[i]][[k]])))
          mins=c(mins,min( log2(toplot$dynamics$totallist_profile[[i]][[k]])[!is.infinite(log2(toplot$dynamics$totallist_profile[[i]][[k]]))]  ))
        }else{
          maxs=c(maxs,max(toplot$dynamics$totallist_profile[[i]][[k]]))
          mins=c(mins,min(toplot$dynamics$totallist_profile[[i]][[k]]))            
        }

      }
    }
    maxvalue=max(maxs)
    minvalue=min(mins[!is.infinite(mins)])
    if(input$islogforDynamics){
      toplot$dynamics$ylabel=paste("log2",toplot$dynamics$ylabel)
    }

    par(mar=c(5,5,1,4))
    plot(1, type="n", xlab="",xaxt="n", ylab=toplot$dynamics$ylabel, xlim=c(1, toplot$dynamics$nbin), ylim=c(minvalue, maxvalue))
    axis(1,at=c(1,toplot$dynamics$nbin/2 +0.5,toplot$dynamics$nbin),labels=c("TSS-30%","genes","TES+30%"))
    count=1
    for(i in 1:length(toplot$dynamics$totallist_profile)){
      for(k in 1:length(toplot$dynamics$totallist_profile[[i]])){
        if(input$islogforDynamics){
          lines( log2(toplot$dynamics$totallist_profile[[i]][[k]]),lwd=2,col= toplot$dynamics$cols[count],pch=".",cex=2,xaxt="n")
        }else{
          lines(toplot$dynamics$totallist_profile[[i]][[k]],lwd=2,col= toplot$dynamics$cols[count],pch=".",cex=2,xaxt="n")
        }
        count=count+1
      }
    }
    legend("topright",legend=toplot$dynamics$totalnames,col= toplot$dynamics$cols,lty=rep(1,length(toplot$dynamics$totalnames)),bg="transparent")  

    dev.off()
  } 
)

#save PDF of TSS boxplot
output$saveboxTSSDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('TSS_boxplot.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    
    currentlist=toplot$dynamics$totallist_boxplot[[1]]
    list2=list()
    #transform 2-order list in 1 order list
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }

    #HERE put the log2 of input$islogforDynamics
    #if Infinite values, those won't be drown
    if(input$islogforDynamics){
      for(i in 1:length(list2)){
        list2[[i]]=log2(list2[[i]])
      }
    }

    #boxplot list2. space at each ROI
    numroi=length(currentlist)
    numbam=length(currentlist[[1]])
    factor_add=rep(0:(numroi-1),each=numbam)
    addingfactor=1:length(list2)+factor_add
    ats=c()
    for(i in 1:numroi){
      window=addingfactor[(((i-1)*numbam)+1):(i*numbam)]
      currentvalue=(window[length(window)]-window[1])/2
      currentvalue=window[1]+currentvalue
      ats=c(ats,currentvalue)
    }

    if(input$islogforDynamics){
      toplot$dynamics$ylabel=paste("log2",toplot$dynamics$ylabel)
    }
    par(mar=c(10,4,1,1),xpd=TRUE)
    suppressWarnings(boxplot(list2,at=addingfactor,col=toplot$dynamics$cols,ylab=toplot$dynamics$ylabel,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
    axis(1,at=ats,label=toplot$dynamics$roinames)

    legend("bottom",inset=c(0,-0.4),legend=toplot$dynamics$totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$dynamics$totalnames)))


    dev.off()
  } 
)





#save PDF of GB boxplot
output$saveboxGBDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('GB_boxplot.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    
    currentlist=toplot$dynamics$totallist_boxplot[[2]]
    list2=list()
    #transform 2-order list in 1 order list
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }
    #HERE put the log2 of input$islogforDynamics
    #if Infinite values, those won't be drown
    if(input$islogforDynamics){
      for(i in 1:length(list2)){
        list2[[i]]=log2(list2[[i]])
      }
    }
    #boxplot list2. space at each ROI
    numroi=length(currentlist)
    numbam=length(currentlist[[1]])
    factor_add=rep(0:(numroi-1),each=numbam)
    addingfactor=1:length(list2)+factor_add
    ats=c()
    for(i in 1:numroi){
      window=addingfactor[(((i-1)*numbam)+1):(i*numbam)]
      currentvalue=(window[length(window)]-window[1])/2
      currentvalue=window[1]+currentvalue
      ats=c(ats,currentvalue)
    }
    if(input$islogforDynamics){
      toplot$dynamics$ylabel=paste("log2",toplot$dynamics$ylabel)
    }
    par(mar=c(10,4,1,1),xpd=TRUE)
    suppressWarnings(boxplot(list2,at=addingfactor,col=toplot$dynamics$cols,ylab=toplot$dynamics$ylabel,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
    axis(1,at=ats,label=toplot$dynamics$roinames)

    legend("bottom",inset=c(0,-0.4),legend=toplot$dynamics$totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$dynamics$totalnames)))
    dev.off()
  } 
)



#save PDF of TES boxplot
output$saveboxTESDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('TES_boxplot.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    
    currentlist=toplot$dynamics$totallist_boxplot[[3]]
    list2=list()
    #transform 2-order list in 1 order list
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }

    #HERE put the log2 of input$islogforDynamics
    #if Infinite values, those won't be drown
    if(input$islogforDynamics){
      for(i in 1:length(list2)){
        list2[[i]]=log2(list2[[i]])
      }
    }

    #boxplot list2. space at each ROI
    numroi=length(currentlist)
    numbam=length(currentlist[[1]])
    factor_add=rep(0:(numroi-1),each=numbam)
    addingfactor=1:length(list2)+factor_add
    ats=c()
    for(i in 1:numroi){
      window=addingfactor[(((i-1)*numbam)+1):(i*numbam)]
      currentvalue=(window[length(window)]-window[1])/2
      currentvalue=window[1]+currentvalue
      ats=c(ats,currentvalue)
    }
    if(input$islogforDynamics){
      toplot$dynamics$ylabel=paste("log2",toplot$dynamics$ylabel)
    }
    par(mar=c(10,4,1,1),xpd=TRUE)
    suppressWarnings(boxplot(list2,at=addingfactor,col=toplot$dynamics$cols,ylab=toplot$dynamics$ylabel,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
    axis(1,at=ats,label=toplot$dynamics$roinames)

    legend("bottom",inset=c(0,-0.4),legend=toplot$dynamics$totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$dynamics$totalnames)))
    dev.off()
  } 
)




#observer for download button data of boxplot (TSS,GB,TES)
output$saveenrichmentBoxdynamicsTSSdata<- downloadHandler(
  filename=function() {
      paste('boxplot_data_TSS.xls', sep='')
  },
  content=function(file) {
    currentlist=toplot$dynamics$totallist_boxplot[[1]]
    list2=list()
    #transform 2-order list in 1 order list
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }

    #HERE put the log2 of input$islogforDynamics
    #if Infinite values, those won't be drown
    if(input$islogforDynamics){
      for(i in 1:length(list2)){
        list2[[i]]=log2(list2[[i]])
      }
    }

    names(list2)=toplot$dynamics$totalnames
    maxvalues=max(unlist(lapply(list2,length)))
    arr=matrix(rep("",maxvalues*length(list2)),ncol=length(list2))
    for(i in 1:length(list2)){
      arr[1:length(list2[[i]]),i]=list2[[i]]
    }
    colnames(arr)=names(list2)
    colnames(arr)=gsub(" ","_",colnames(arr))
    write.table(arr,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  } 
)


output$saveenrichmentBoxdynamicsGBdata<- downloadHandler(
  filename=function() {
      paste('boxplot_data_GB.xls', sep='')
  },
  content=function(file) {
    currentlist=toplot$dynamics$totallist_boxplot[[2]]
    list2=list()
    #transform 2-order list in 1 order list
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }

    #HERE put the log2 of input$islogforDynamics
    #if Infinite values, those won't be drown
    if(input$islogforDynamics){
      for(i in 1:length(list2)){
        list2[[i]]=log2(list2[[i]])
      }
    }

    names(list2)=toplot$dynamics$totalnames
    maxvalues=max(unlist(lapply(list2,length)))
    arr=matrix(rep("",maxvalues*length(list2)),ncol=length(list2))
    for(i in 1:length(list2)){
      arr[1:length(list2[[i]]),i]=list2[[i]]
    }
    colnames(arr)=names(list2)
    colnames(arr)=gsub(" ","_",colnames(arr))
    write.table(arr,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  } 
)

output$saveenrichmentBoxdynamicsTESdata<- downloadHandler(
  filename=function() {
      paste('boxplot_data_TES.xls', sep='')
  },
  content=function(file) {
    currentlist=toplot$dynamics$totallist_boxplot[[3]]
    list2=list()
    #transform 2-order list in 1 order list
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }

    #HERE put the log2 of input$islogforDynamics
    #if Infinite values, those won't be drown
    if(input$islogforDynamics){
      for(i in 1:length(list2)){
        list2[[i]]=log2(list2[[i]])
      }
    }

    names(list2)=toplot$dynamics$totalnames
    maxvalues=max(unlist(lapply(list2,length)))
    arr=matrix(rep("",maxvalues*length(list2)),ncol=length(list2))
    for(i in 1:length(list2)){
      arr[1:length(list2[[i]]),i]=list2[[i]]
    }
    colnames(arr)=names(list2)
    colnames(arr)=gsub(" ","_",colnames(arr))
    write.table(arr,file=file,row.names=FALSE,sep="\t",quote=FALSE   ) 
  } 
)



#save PDF of SI on TSS
output$saveSITSSDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('TSS_SI_dynamics.pdf', sep='')
  },
  content=function(file) {
    
    currentlist=toplot$dynamics$totallist_SI[[1]]

    list2=list()
    #transform 2-order list in 1 order list
    minvals=c()
    maxvals=c()
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }

	outlayer_thresh=input$percentageOutlayerCumulPlots
    for(i in 1:length(list2)){
	  mincurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],outlayer_thresh )
	  maxcurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],1-outlayer_thresh )
	  minvals=c(minvals,mincurrent)
	  maxvals=c(maxvals,maxcurrent)
	}
	  
	maxvals=max(maxvals)
	minvals=min(minvals)

    lengths=sapply(list2,length)

    pdf(file)
    par(mar=c(10,4,1,1))
    plot(1, type="n", ylab="Ranked Genes", xlab=paste("log2",toplot$dynamics$ylabel), ylim=c(0, quantile(1:max(lengths),1-(2*outlayer_thresh))), xlim=c(minvals, maxvals))
    for(i in 1:length(list2)){
      vals=sort(log2(list2[[i]]))
      vals=vals[!is.infinite(vals)]
      vals=vals[vals<quantile(vals,1-outlayer_thresh) & vals>quantile(vals,outlayer_thresh)]
      lines(vals,1:length(vals),col=toplot$dynamics$cols[i],lwd=2)
    }

    legend("bottom",inset=c(0,-0.4),legend=toplot$dynamics$totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$dynamics$totalnames)))      
    dev.off()
  } 
)






#save PDF of SI on GB
output$saveSIGBDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('GB_SI_dynamics.pdf', sep='')
  },
  content=function(file) {
    
    
    currentlist=toplot$dynamics$totallist_SI[[2]]

    list2=list()
    #transform 2-order list in 1 order list
    minvals=c()
    maxvals=c()
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }



    outlayer_thresh=input$percentageOutlayerCumulPlots
    for(i in 1:length(list2)){
      mincurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],outlayer_thresh )
      maxcurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],1-outlayer_thresh )
      minvals=c(minvals,mincurrent)
      maxvals=c(maxvals,maxcurrent)
    }
    maxvals=max(maxvals)
    minvals=min(minvals)

    lengths=sapply(list2,length)

    pdf(file)
    par(mar=c(10,4,1,1))

    plot(1, type="n", ylab="Ranked Genes", xlab=paste("log2",toplot$dynamics$ylabel), ylim=c(0, quantile(1:max(lengths),1-(2*outlayer_thresh)) ), xlim=c(minvals, maxvals))
    for(i in 1:length(list2)){
      vals=sort(log2(list2[[i]]))
      vals=vals[!is.infinite(vals)]
      vals=vals[vals<quantile(vals,1-outlayer_thresh) & vals>quantile(vals,outlayer_thresh)]
      lines(vals,1:length(vals),col=toplot$dynamics$cols[i],lwd=2)
    }

    legend("bottom",inset=c(0,-0.4),legend=toplot$dynamics$totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$dynamics$totalnames)))      

    dev.off()
  } 
)




#save PDF of SI on SI
output$saveSISIDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('SI_SI_dynamics.pdf', sep='')
  },
  content=function(file) {
    pdf(file)
    
    currentlist=toplot$dynamics$totallist_SI[[3]]

    list2=list()
    #transform 2-order list in 1 order list
    minvals=c()
    maxvals=c()
    for(i in 1:length(currentlist)){
      list2=c(list2,currentlist[[i]])
    }


    outlayer_thresh=input$percentageOutlayerCumulPlots
    for(i in 1:length(list2)){
      mincurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],outlayer_thresh )
      maxcurrent=quantile(log2(list2[[i]])[!is.infinite(log2(list2[[i]]))],1-outlayer_thresh )
      minvals=c(minvals,mincurrent)
      maxvals=c(maxvals,maxcurrent)
    }


    maxvals=max(maxvals)
    minvals=min(minvals)

    lengths=sapply(list2,length)
    par(mar=c(10,4,1,1))
    plot(1, type="n", ylab="Ranked Genes", xlab="stalling index", ylim=c(0, quantile(1:max(lengths),1-(2*outlayer_thresh))), xlim=c(minvals, maxvals))
    for(i in 1:length(list2)){
      #remove a fraction of outlayers
      #celan -Inf
      vals=sort(log2(list2[[i]]))
      vals=vals[!is.infinite(vals)]
      vals=vals[vals<quantile(vals,1-outlayer_thresh) & vals>quantile(vals,outlayer_thresh)]
      lines(vals,1:length(vals),col=toplot$dynamics$cols[i],lwd=2)
    }
    legend("bottom",inset=c(0,-0.4),legend=toplot$dynamics$totalnames,col=toplot$dynamics$cols,cex=0.6,bg="transparent",pch=rep(19,length(toplot$dynamics$totalnames)))      

    dev.off()
  } 
)


