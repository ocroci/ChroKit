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



#help buttons options
observeEvent(input$help_singleEvaluation_normalizationtotalread, {boxHelpServer(help_singleEvaluation_normalizationtotalread)})
observeEvent(input$help_singleEvaluation_normalizationreaddensity, {boxHelpServer(help_singleEvaluation_normalizationreaddensity)})





observeEvent(input$plotSingleEval,{
	set.seed(123)
	if (!is.null(input$ROIchooseSingleEval) & length(ROIvariables$listROI)>0) {
		nomi=unlist(lapply(ROIvariables$listROI,getName))
    	if ("promoters"%in% nomi & "transcripts"%in% nomi) {
      		pos=match(input$ROIchooseSingleEval,nomi)
      		roi=uniqueROI(ROIvariables$listROI[[pos]]) 
          duplicated_pos_range=!duplicated(getRange(ROIvariables$listROI[[pos]]))
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
            
            #build the interface to piechartSingleEval_options options (appear only when plot)
            output$piechartSingleEval_options<-renderUI({selectInput("chooseColorPaletteSingleEval","Choose color palette:",choices=c(
                                                  "red/blue/grey"="black_red_blue_gray",
                                                  "red/orange/green"="black_red_orange_green"

                                                ))})

            
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
              #make it reactive to input$choosecolor.... otherwise default values if option still does not exist
              if(isvalid(input$chooseColorPaletteSingleEval)){
                colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]][2:4]
              }else{
                colors=strsplit("black_red_blue_gray",split="_")[[1]][2:4]
              }
              
              plot_singleeval_pie(elements=elements,colors=colors,title=maintitle,label=label)
            })


            #barplot plot
            output$viewDistributionBarSingleEval<-renderPlot({
              #make it reactive to input$choosecolor....
              if(isvalid(input$chooseColorPaletteSingleEval)){
                colors=strsplit(input$chooseColorPaletteSingleEval,split="_")[[1]][2:4]
              }else{
                colors=strsplit("black_red_blue_gray",split="_")[[1]][2:4]
              }
              plot_singleeval_bar(elements=elements,colors=colors,title=maintitle,label=label)
            })
            #buttons for download
            output$saveviewDistributionPieSingleEval=renderUI({downloadButton('saveviewDistributionPieSingleEvalbutton', 'Get PDF')})
            output$saveviewDistributionBarSingleEval=renderUI({downloadButton('saveviewDistributionBarSingleEvalbutton', 'Get PDF')})
            

            print(paste("Single evaluation of",input$ROIchooseSingleEval))


            #width distribution
            if (length(range)>2){
              quant=0.95
              cuttedga=range[width(range)<quantile(width(range),quant)]
              if (length(range_promo)>2&length(cuttedga)>2){
                #build the interface to densitySingleEval_options options (appear only when plot density)
                output$densitySingleEval_options<-renderUI({selectInput("chooseColorPaletteSingleEval_density","Choose color palette:",choices=c(
                                                  "black/red/blue/grey"="black_red_blue_gray",
                                                  "black/red/orange/green"="black_red_orange_green"
                                                ))})
                #density plot can be done
                output$widthDistributionSingleEval<-renderPlot({
                  if(isvalid(input$chooseColorPaletteSingleEval_density)){
                    colors=strsplit(input$chooseColorPaletteSingleEval_density,split="_")[[1]]
                  }else{
                    colors=strsplit("black_red_blue_gray",split="_")[[1]]
                  }
                  
                  plot_singleeval_density(range=range,range_promo=range_promo,range_intra=range_intra,
                                          range_inter=range_inter,colors=colors)
                })
                #download button for PDF density plot
                output$savewidthDistributionSingleEval<-renderUI({downloadButton('savewidthDistributionSingleEvalbutton', 'Get PDF')})

              }else{
                #only total plot, but at least one between promoters, intra, inter must be >0 by definition
                output$widthDistributionSingleEval<-renderPlot({NULL})
                output$densitySingleEval_options<-renderUI({NULL})
                #no PDF button, because no plot
                output$savewidthDistributionSingleEval<-renderUI({NULL})
              }

            }else{
              #range has length==0
              output$widthDistributionSingleEval<-renderPlot({NULL})
              output$densitySingleEval_options<-renderUI({NULL})
              #no PDF button, because no plot
              output$savewidthDistributionSingleEval<-renderUI({NULL})
            }

            #extract enrichment from correct position of the ROI selected
            rawvals=Enrichlist$rawcoverage[[pos]]
            keyvals=Enrichlist$decryptkey[[pos]]
            normvals=Enrichlist$normfactlist[[pos]]
      		  getbam=names(rawvals)
            
        		if (!is.null(getbam)& length(getbam)>0){
        		  pos2=match(input$BAMchooseSingleEval,getbam)
      			  bam_orig=rawvals[[pos2]]
              keyvals_orig=keyvals[[pos2]]
              norm_orig=normvals[[pos2]]
              
              #here MUST remove duplicated ranges (because before we used "uniqueROI")
              bam_orig=bam_orig[duplicated_pos_range]
              keyvals_orig=keyvals_orig[duplicated_pos_range]

      			  #calculate enrichments for boxplots
      			  bam_promo_orig=bam_orig[ov_range>0]
              key_promo_orig=keyvals_orig[ov_range>0]
      			  bam_ws_orig=bam_orig[ov_range==0]
              key_ws_orig=keyvals_orig[ov_range==0]
      			  bam_intra_orig=bam_ws_orig[ov_transcripts>0]
              key_intra_orig=key_ws_orig[ov_transcripts>0]
      			  bam_inter_orig=bam_ws_orig[ov_transcripts==0]
              key_inter_orig=key_ws_orig[ov_transcripts==0]
 
              #not very efficient, because the matrix could be calculated once, but fast enough
              bam=as.integer(makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bam_orig,Nbins=1,Snorm=FALSE,key=keyvals_orig,norm_factor=norm_orig))
              bam_promo=bam[ov_range>0]
              bam_ws=bam[ov_range==0]
              bam_intra=bam_ws[ov_transcripts>0]
              bam_inter=bam_ws[ov_transcripts==0]

              toplot$viewDistributionPieSingleEval$bam=bam
              toplot$viewDistributionPieSingleEval$bam_promo=bam_promo
              toplot$viewDistributionPieSingleEval$bam_intra=bam_intra
              toplot$viewDistributionPieSingleEval$bam_inter=bam_inter



              listtoprofile=list(bam_orig,bam_promo_orig,bam_intra_orig,bam_inter_orig)
              keystoprofile=list(keyvals_orig,key_promo_orig,key_intra_orig,key_inter_orig)
              rangs=list(range,range_promo,range_intra,range_inter)
              lengths_sampled=as.list(rep(NA,4))
              for(i in 1:length(listtoprofile)){
                set.seed(123)
                smp=sample(1:length(listtoprofile[[i]]),floor(length(listtoprofile[[i]])/10),replace=FALSE)
                listtoprofile[[i]]=listtoprofile[[i]][smp]
                keystoprofile[[i]]=keystoprofile[[i]][smp]
                lengths_sampled[[i]]=width(rangs[[i]])[smp]
              }
              #find matrix in bins for promo,intra,inter
              matrixes=list()
              for(i in 1:length(listtoprofile)){
                if(length(listtoprofile[[i]])>0){
                  matrixes[[i]]=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=listtoprofile[[i]],Nbins=50,Snorm=FALSE,key=keystoprofile[[i]],norm_factor=norm_orig)
                }else{
                  matrixes[[i]]=matrix(rep(0,50),nrow=1)
                }
              }
              toplot$viewDistributionPieSingleEval$matrixes=matrixes
              toplot$viewDistributionPieSingleEval$lengths_sampled=lengths_sampled


              #plots will be plotted: add UI options for color and normalization
              output$boxSingleEval_options<-renderUI({
                                              list(
                                                selectInput("chooseColorPaletteSingleEval_box","Choose color palette:",choices=c(
                                                  "black/red/blue/grey"="black_red_blue_gray",
                                                  "black/red/orange/green"="black_red_orange_green")),
                                                radioButtons("chooseNormalizationSingleEval_box","Choose normalization:",
                                                      choiceNames=list(
                                                        htmlhelp("Total reads","help_singleEvaluation_normalizationtotalread"),
                                                        htmlhelp("Read density (reads/bp)","help_singleEvaluation_normalizationreaddensity")
                                                      ),
                                                      choiceValues=list(
                                                        "totread",
                                                        "readdensity"
                                                      )                                                  
                                                      )

                                                )


                                               })
              output$profileSingleEval_options<-renderUI({
                                                list(
                                                selectInput("chooseColorPaletteSingleEval_profile","Choose color palette:",choices=c(
                                                  "black/red/blue/grey"="black_red_blue_gray",
                                                  "black/red/orange/green"="black_red_orange_green")),
                                                radioButtons("chooseNormalizationSingleEval_profile","Choose normalization:",choices=c(
                                                  "Total reads"="totread",
                                                  "Read density (reads/bp)"="readdensity"
                                                    ))
                                                )

                                               })



              #plot enrchment box (reactive to color option and normalization)
              output$enrichmentBoxSingleEval<-renderPlot({
                if(isvalid(input$chooseColorPaletteSingleEval_box)){
                    colors=strsplit(input$chooseColorPaletteSingleEval_box,split="_")[[1]]
                }else{
                    colors=strsplit("black_red_blue_gray",split="_")[[1]]
                }
                if(isvalid(input$chooseNormalizationSingleEval_box)){
                  normalization=input$chooseNormalizationSingleEval_box
                }else{
                  normalization="totread"
                }

                plot_singleeval_boxplot(normalization=normalization,bam=bam,bam_promo=bam_promo,bam_intra=bam_intra,bam_inter=bam_inter,
                                        range=range,range_promo=range_promo,range_intra=range_intra,range_inter=range_inter,colors=colors)
              })


              #plot enrichment profile (reactive to color option and normalization)
              output$TSSprofileSingleEval<-renderPlot({
                if(isvalid(input$chooseColorPaletteSingleEval_profile)){
                    colors=strsplit(input$chooseColorPaletteSingleEval_profile,split="_")[[1]]
                }else{
                    colors=strsplit("black_red_blue_gray",split="_")[[1]]
                }
                if(isvalid(input$chooseNormalizationSingleEval_profile)){
                  normalization=input$chooseNormalizationSingleEval_profile
                  if(normalization=="totread"){
                    ylab="Total reads"
                  }else{
                    ylab="Read density (reads/bp)"
                    #divide for length of ranges
                    for (i in 1:length(matrixes)){
                      #equivalent. Dividing by an array with l=number of rows of the matrix, will divide each col by these numbers
                      if (length(lengths_sampled[[i]])>0){
                        matrixes[[i]]=matrixes[[i]]/lengths_sampled[[i]]
                      }
                      
                    }
                  }
                }else{
                  normalization="totread"
                  ylab="Total reads"
                }

                profile_to_plot=lapply(matrixes,function(i) {apply(i,2,mean)})
                names(profile_to_plot)=c("all","promoters","genebody","intergenic")
                plot_singleeval_profile(colors=colors,profile_to_plot=profile_to_plot,ylab=ylab)

              })

              #buttons for PDFs and data
              output$saveboxdataSingleEval=renderUI({downloadButton('saveenrichmentBoxSingleEvaldata', 'Save data')})
              output$saveenrichmentBoxSingleEval=renderUI({downloadButton('saveenrichmentBoxSingleEvalbutton', 'Get PDF')})
              output$saveenrichmentProfileSingleEval=renderUI({downloadButton('saveenrichmentProfileSingleEvalbutton', 'Get PDF')})
                          

      			}else{
      				#roi doesn't have bam...
              output$enrichmentBoxSingleEval<-renderPlot({plot_text(text="you need to\n associate an \nenrichment file",cex=1.4)})
              output$TSSprofileSingleEval<-renderPlot({plot_text(text="you need to\n associate an  \nenrichment file",cex=1.4)})
              toplot$viewDistributionPieSingleEval$Colors_profile=NULL
              toplot$viewDistributionPieSingleEval$Colors_box=NULL
              toplot$viewDistributionPieSingleEval$bam=NULL
              toplot$viewDistributionPieSingleEval$bam_promo=NULL
              toplot$viewDistributionPieSingleEval$bam_intra=NULL
              toplot$viewDistributionPieSingleEval$bam_inter=NULL
              toplot$viewDistributionPieSingleEval$matrixes=NULL
              output$saveboxdataSingleEval=renderUI({NULL})
              output$saveenrichmentBoxSingleEval=renderUI({NULL})
              output$saveenrichmentProfileSingleEval=renderUI({NULL})
              output$boxSingleEval_options<-renderUI({NULL})
              output$profileSingleEval_options<-renderUI({NULL})
            }
      		}else{
            #not easy to be here. if roi is null, something should go wrong upstream...
            output$viewDistributionPieSingleEval<-renderPlot({NULL})
            output$viewDistributionBarSingleEval<-renderPlot({NULL})
            output$saveviewDistributionPieSingleEval=renderUI({NULL})
            output$saveviewDistributionBarSingleEval=renderUI({NULL})
            output$savewidthDistributionSingleEval<-renderUI({NULL})
            output$widthDistributionSingleEval<-renderPlot({NULL})
            output$enrichmentBoxSingleEval<-renderPlot({plot_text(text="you need to\n associate an \nenrichment file",cex=1.4)})
            output$TSSprofileSingleEval<-renderPlot({plot_text(text="you need to\n associate an  \nenrichment file",cex=1.4)})
      			#logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">ROI not found...<br></font>'
            toplot$viewDistributionPieSingleEval$Colors=NULL
            toplot$viewDistributionPieSingleEval$Colors_density=NULL
            toplot$viewDistributionPieSingleEval$Colors_profile=NULL
            toplot$viewDistributionPieSingleEval$Colors_box=NULL
            toplot$viewDistributionPieSingleEval$elements=NULL
            toplot$viewDistributionPieSingleEval$maintitle=NULL
            toplot$viewDistributionPieSingleEval$label=NULL
            toplot$viewDistributionPieSingleEval$nomi=NULL
            toplot$viewDistributionPieSingleEval$range=NULL
            toplot$viewDistributionPieSingleEval$range_promo=NULL
            toplot$viewDistributionPieSingleEval$range_intra=NULL
            toplot$viewDistributionPieSingleEval$range_inter=NULL
            toplot$viewDistributionPieSingleEval$matrixes=NULL
            output$piechartSingleEval_options<-renderUI({NULL})
            output$densitySingleEval_options<-renderUI({NULL})
            output$boxSingleEval_options<-renderUI({NULL})
            output$profileSingleEval_options<-renderUI({NULL})
            output$saveboxdataSingleEval=renderUI({NULL})
            output$saveenrichmentBoxSingleEval=renderUI({NULL})
            output$saveenrichmentProfileSingleEval=renderUI({NULL})
      		}

    	}else{
        output$viewDistributionPieSingleEval<-renderPlot({NULL})
        output$viewDistributionBarSingleEval<-renderPlot({NULL})
        output$saveviewDistributionPieSingleEval=renderUI({NULL})
        output$saveviewDistributionBarSingleEval=renderUI({NULL})
        output$savewidthDistributionSingleEval<-renderUI({NULL})
        output$widthDistributionSingleEval<-renderPlot({NULL})
        output$enrichmentBoxSingleEval<-renderPlot({plot_text(text="you need to\n associate an \nenrichment file",cex=1.4)})
        output$TSSprofileSingleEval<-renderPlot({plot_text(text="you need to\n associate an  \nenrichment file",cex=1.4)})
    		#logvariables$msg[[length(logvariables$msg)+1]]= '<font color="red">promoters/transcripts not found... ask to database<br></font>'
        sendSweetAlert(
          session = session,
          title = "Annotated elements not found",
          text = "Promoters, transcripts and TES of a specific genome assembly not found, but you need them: go to 'Databases' section and choose a genome assembly",
          type = "error"
        )        
        toplot$viewDistributionPieSingleEval$Colors=NULL
        toplot$viewDistributionPieSingleEval$Colors_density=NULL
        toplot$viewDistributionPieSingleEval$Colors_profile=NULL
        toplot$viewDistributionPieSingleEval$Colors_box=NULL
        toplot$viewDistributionPieSingleEval$elements=NULL
        toplot$viewDistributionPieSingleEval$maintitle=NULL
        toplot$viewDistributionPieSingleEval$label=NULL
        toplot$viewDistributionPieSingleEval$nomi=NULL
        toplot$viewDistributionPieSingleEval$range=NULL
        toplot$viewDistributionPieSingleEval$range_promo=NULL
        toplot$viewDistributionPieSingleEval$range_intra=NULL
        toplot$viewDistributionPieSingleEval$range_inter=NULL
        toplot$viewDistributionPieSingleEval$matrixes=NULL
        output$piechartSingleEval_options<-renderUI({NULL})
        output$densitySingleEval_options<-renderUI({NULL})
        output$boxSingleEval_options<-renderUI({NULL})
        output$profileSingleEval_options<-renderUI({NULL})
        output$saveboxdataSingleEval=renderUI({NULL})
        output$saveenrichmentBoxSingleEval=renderUI({NULL})
        output$saveenrichmentProfileSingleEval=renderUI({NULL})
    	}

	}else{
		#roi not found 
    output$viewDistributionPieSingleEval<-renderPlot({NULL})
    output$viewDistributionBarSingleEval<-renderPlot({NULL})
    output$saveviewDistributionPieSingleEval=renderUI({NULL})
    output$saveviewDistributionBarSingleEval=renderUI({NULL})
    output$savewidthDistributionSingleEval<-renderUI({NULL})
    output$widthDistributionSingleEval<-renderPlot({NULL})
    output$enrichmentBoxSingleEval<-renderPlot({plot_text(text="you need to\n associate an \nenrichment file",cex=1.4)})
    output$TSSprofileSingleEval<-renderPlot({plot_text(text="you need to\n associate an  \nenrichment file",cex=1.4)})
    toplot$viewDistributionPieSingleEval$Colors=NULL
    toplot$viewDistributionPieSingleEval$Colors_density=NULL
    toplot$viewDistributionPieSingleEval$Colors_profile=NULL
    toplot$viewDistributionPieSingleEval$Colors_box=NULL
    toplot$viewDistributionPieSingleEval$elements=NULL
    toplot$viewDistributionPieSingleEval$maintitle=NULL
    toplot$viewDistributionPieSingleEval$label=NULL
    toplot$viewDistributionPieSingleEval$nomi=NULL
    toplot$viewDistributionPieSingleEval$range=NULL
    toplot$viewDistributionPieSingleEval$range_promo=NULL
    toplot$viewDistributionPieSingleEval$range_intra=NULL
    toplot$viewDistributionPieSingleEval$range_inter=NULL
    toplot$viewDistributionPieSingleEval$matrixes=NULL
    output$saveboxdataSingleEval=renderUI({NULL})
    output$piechartSingleEval_options<-renderUI({NULL})
    output$densitySingleEval_options<-renderUI({NULL})
    output$boxSingleEval_options<-renderUI({NULL})
    output$profileSingleEval_options<-renderUI({NULL})
    output$saveboxdataSingleEval=renderUI({NULL})
    output$saveenrichmentBoxSingleEval=renderUI({NULL})
    output$saveenrichmentProfileSingleEval=renderUI({NULL})
	}

},ignoreInit=TRUE)






output$saveviewDistributionPieSingleEvalbutton<- downloadHandler(
  filename=function() {
      paste('Pie_location_single.pdf', sep='')
  },
  content=function(file) {
      pdf(file)
      plot_singleeval_pie(elements=toplot$viewDistributionPieSingleEval$elements,colors=toplot$viewDistributionPieSingleEval$Colors[2:4],
        title=toplot$viewDistributionPieSingleEval$maintitle,label=toplot$viewDistributionPieSingleEval$label)
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
      plot_singleeval_bar(elements=toplot$viewDistributionPieSingleEval$elements,colors=toplot$viewDistributionPieSingleEval$Colors_density,
                      title=toplot$viewDistributionPieSingleEval$maintitle,label=toplot$viewDistributionPieSingleEval$label)
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
        plot_singleeval_density(range=toplot$viewDistributionPieSingleEval$range,range_promo=toplot$viewDistributionPieSingleEval$range_promo,
                                range_intra=toplot$viewDistributionPieSingleEval$range_intra,range_inter=toplot$viewDistributionPieSingleEval$range_inter)
        dev.off()
  } 
)



#button for saving boxplot of SingleEval
output$saveenrichmentBoxSingleEvalbutton<- downloadHandler(
  filename=function() {
      paste('Boxplot_single.pdf', sep='')
  },
  content=function(file) {

      pdf(file)
      plot_singleeval_boxplot(normalization=toplot$viewDistributionPieSingleEval$normalization_box,bam=toplot$viewDistributionPieSingleEval$bam,
                    bam_promo=toplot$viewDistributionPieSingleEval$bam_promo,bam_intra=toplot$viewDistributionPieSingleEval$bam_intra,
                    bam_inter=toplot$viewDistributionPieSingleEval$bam_inter,range=toplot$viewDistributionPieSingleEval$range,
                    range_promo=toplot$viewDistributionPieSingleEval$range_promo,range_intra=toplot$viewDistributionPieSingleEval$range_intra,
                    range_inter=toplot$viewDistributionPieSingleEval$range_inter,colors=toplot$viewDistributionPieSingleEval$Colors_box)
      dev.off()
  } 
)

#observer for download of xls table of the boxplot data of SingleEval
output$saveenrichmentBoxSingleEvaldata<- downloadHandler(
  filename=function() {
      paste('boxplot_data.xls', sep='')
  },
  content=function(file) {
      if(toplot$viewDistributionPieSingleEval$normalization_box=="totread"){
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
        arr[1:length(toplot$viewDistributionPieSingleEval$bam),1]=toplot$viewDistributionPieSingleEval$bam/width(toplot$viewDistributionPieSingleEval$range)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_promo),2]=toplot$viewDistributionPieSingleEval$bam_promo/width(toplot$viewDistributionPieSingleEval$range_promo)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_intra),3]=toplot$viewDistributionPieSingleEval$bam_intra/width(toplot$viewDistributionPieSingleEval$range_intra)
        arr[1:length(toplot$viewDistributionPieSingleEval$bam_inter),4]=toplot$viewDistributionPieSingleEval$bam_inter/width(toplot$viewDistributionPieSingleEval$range_inter)
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
    plot_singleeval_profile(colors=toplot$viewDistributionPieSingleEval$Colors_profile,profile_to_plot=toplot$viewDistributionPieSingleEval$profile_to_plot,
                          ylab=toplot$viewDistributionPieSingleEval$ylabprofile)

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
#overlap + enrichment box
observeEvent(input$msg_pairwiseOverlaps_box, {boxHelpServer(msg_pairwiseOverlaps_box)})
#overlap + enrichment scatter
observeEvent(input$msg_pairwiseOverlaps_scatter, {boxHelpServer(msg_pairwiseOverlaps_scatter)})
#overlap + enrichment calibration
observeEvent(input$msg_pairwiseOverlaps_calibration, {boxHelpServer(msg_pairwiseOverlaps_calibration)})



#react to parameters help button
observeEvent(input$help_pairwiseOverlaps_parameters_enrich1, {boxHelpServer(help_pairwiseOverlaps_parameters_enrich1)})
observeEvent(input$help_pairwiseOverlaps_parameters_enrich2, {boxHelpServer(help_pairwiseOverlaps_parameters_enrich2)})
observeEvent(input$help_pairwiseOverlaps_parameters_minbpoverlap, {boxHelpServer(help_pairwiseOverlaps_parameters_minbpoverlap)})

observeEvent(input$help_singleEvaluation_normalizationtotalread_cmpbox, {boxHelpServer(help_singleEvaluation_normalizationtotalread)})
observeEvent(input$help_singleEvaluation_normalizationreaddensity_cmpbox, {boxHelpServer(help_singleEvaluation_normalizationreaddensity)})
observeEvent(input$help_singleEvaluation_normalizationtotalread_cmpscatter, {boxHelpServer(help_singleEvaluation_normalizationtotalread)})
observeEvent(input$help_singleEvaluation_normalizationreaddensity_cmpscatter, {boxHelpServer(help_singleEvaluation_normalizationreaddensity)})



observeEvent(input$plotCmp,{

  #verify to have at least input$ROI2chooseCmp and input$ROI1chooseCmp 
  if (length(ROIvariables$listROI)>0&length(input$ROI2chooseCmp)>0 & length(input$ROI1chooseCmp)>0 & !is.na(input$minOverlapCmp) ){
    n1=input$ROI1chooseCmp
    n2=input$ROI2chooseCmp

    #further test if ROI have BAM selected
    #some bug in the observer
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos1=match(input$ROI1chooseCmp,nomi)
    BAM1_present=input$BAM1chooseCmp%in%names(Enrichlist$rawcoverage[[pos1]])
    pos2=match(input$ROI2chooseCmp,nomi)
    BAM2_present=input$BAM2chooseCmp%in%names(Enrichlist$rawcoverage[[pos2]])
    if(!isvalid(BAM1_present)){
      BAM1_choose=""
    }else{
      BAM1_choose=input$BAM1chooseCmp
    }
    if(!isvalid(BAM2_present)){
      BAM2_choose=""
    }else{
      BAM2_choose=input$BAM2chooseCmp
    }

    #find the corresponding ROI objects from their names, as usual. Keep only unique ranges (for example,
    #when promoters/transcripts of a genelist). Remember to remove duplicate ranges positions from enrichments
    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos1=match(input$ROI1chooseCmp,nomi)
    roi1=uniqueROI(ROIvariables$listROI[[pos1]])
    pos2=match(input$ROI2chooseCmp,nomi)
    roi2=uniqueROI(ROIvariables$listROI[[pos2]])
    not_dup_range1=!duplicated(getRange(ROIvariables$listROI[[pos1]]))
    not_dup_range2=!duplicated(getRange(ROIvariables$listROI[[pos2]]))
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


    #IMPORTANT: without "unique", we keep all duplicated ranges and consider them as different.
    #this is the case for example promoters/tranascripts of a genelist
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

      #range1_ov=unique(range1[qh])
      #range2_ov=unique(range2[sh])
      range1_ov=unique(range1[qh])
      range2_ov=unique(range2[sh])
      range1_notov=unique(range1[-qh])
      range2_notov=unique(range2[-sh])    
      perc1=  round( (length(range1_ov)/length(range1))*100 )
      perc2=  round( (length(range2_ov)/length(range2))*100 )
      jaccard=round(  length(common)  /length(totalunion),2)
    }

    #make the matrix for barplot and barplot of the overlap (the most correct representation)
    matbar=as.matrix(data.frame(n1=c(length(range1_ov),length(range1_notov),0,0),
                      n2=c(0,0,length(range2_ov),length(range2_notov))))
    colnames(matbar)=paste(c("ROI 1","ROI 2")," (",c(length(range1_ov),length(range2_ov)),"; ",c(perc1,perc2),"%)",sep="")
    
    toplot$cmp$matbar=matbar
    toplot$cmp$jaccard=jaccard
    toplot$cmp$n1=n1
    toplot$cmp$n2=n2

    output$pairwiseoverlaps_overlap_options<-renderUI({
                  selectInput("chooseColorPaletteCmp_overlap","Choose color palette:",choices=c(
                                                  "red/grey"="red_gray_red4_grey20",
                                                  "blue/green"="blue_green_blue4_green4"
                                                ))
    })

    output$viewBarplotCmp<-renderPlot({
      if(isvalid(input$chooseColorPaletteCmp_overlap)){
        Colors=strsplit(input$chooseColorPaletteCmp_overlap,split="_")[[1]]
      }else{
        Colors=strsplit("red_gray_red4_grey20",split="_")[[1]]
      }
      plot_cmp_barplot(matbar=matbar,n1=n1,n2=n2,colors=Colors,jaccard=jaccard)
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
      if(isvalid(input$chooseColorPaletteCmp_overlap)){
        Colors=strsplit(input$chooseColorPaletteCmp_overlap,split="_")[[1]]
      }else{
        Colors=strsplit("red_gray_red4_grey20",split="_")[[1]]
      }
      plot_cmp_venn(area1=area1,area2=area2,common=common,n1=n1,n2=n2,colors=Colors)
    })


    #download button for PDF Barplot of overlap
    output$saveviewBarplotCmp=renderUI({downloadButton('saveviewBarplotCmpbutton', 'Get PDF')})
    #download button for Venn
    output$saveviewVennCmp=renderUI({downloadButton('saveviewVennCmpbutton', 'Get PDF')})
    print(paste("Pairwise comparison between",input$ROI1chooseCmp,"and",input$ROI2chooseCmp))



    #now, check BAMlists of the 2 ROIs selected
    rawvals1=Enrichlist$rawcoverage[[pos1]]
    keyvals1=Enrichlist$decryptkey[[pos1]]
    normvals1=Enrichlist$normfactlist[[pos1]]
    rawvals2=Enrichlist$rawcoverage[[pos2]]
    keyvals2=Enrichlist$decryptkey[[pos2]]
    normvals2=Enrichlist$normfactlist[[pos2]]

    getbam1=names(rawvals1)  
    getbam2=names(rawvals2)   

    toplot$cmp$getbam1=getbam1
    toplot$cmp$getbam2=getbam2
    toplot$cmp$BAM1chooseCmp=BAM1_choose
    toplot$cmp$BAM2chooseCmp=BAM2_choose


    #if there is at least one overlap and at least one BAMfile has been selectd by the user
    #(meaning that is associated with at least one of the 2 ROIs selected), proceed
    if ( (nchar(BAM1_choose)>0 |nchar(BAM2_choose)>0)  & (!is.null(getbam1) | !is.null(getbam2)) & length(ov)>0){
      ov_range1=findOverlaps(common,range1,minoverlap=input$minOverlapCmp)
      ov_range2=findOverlaps(common,range2,minoverlap=input$minOverlapCmp)
      lab_axis="Total reads"
      #if BAM1 is present and selected, proceed
      if(nchar(BAM1_choose)>0 ){
        pos_bam1_1=match(BAM1_choose,getbam1)
        bam1_1=rawvals1[[pos_bam1_1]]
        keys1_1=keyvals1[[pos_bam1_1]]
        norm1_1=normvals1[[pos_bam1_1]]
        #remove dups range1
        bam1_1=bam1_1[not_dup_range1]
        keys1_1=keys1_1[not_dup_range1]

        #decrypt and sum
        bam1_1=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bam1_1,Nbins=1,Snorm=FALSE,key=keys1_1,norm_factor=norm1_1)
        if(!is.null(bam1_1)){
          range1only_bam1=bam1_1[-qh]
          #extract BAM1 of range1, in positions of range1 that overlaps with range2 
          bam1_common=bam1_1[subjectHits(ov_range1)]

          #   #here, we can have that multiple range1 overlap with a single range 2 in "common" regions.
          #   #in this case, for each overlapping AREA, find the MEAN of signals of BAM1 of ranges1
          #   #that overlap with this single range2:
          #   #range1: ------ ---------           ---
          #   #range2:  ----------------------------
          #   #the average of signals of BAM1 calculated on range1 is the signal that
          #   #will appear as "common regions"

          #if BAM1 is associated also with range2, we can have the reads of BAM1
          #but in the exclusive regions of ROI2!!
          pos_bam1_2=match(BAM1_choose,getbam2)
          if(!is.na(pos_bam1_2) & length(pos_bam1_2)>0){
            bam1_2=rawvals2[[pos_bam1_2]]
            keys1_2=keyvals2[[pos_bam1_2]]
            norm1_2=normvals2[[pos_bam1_2]]
            #remove dups range2
            bam1_2=bam1_2[not_dup_range2]
            keys1_2=keys1_2[not_dup_range2]

            #decrypt and sum
            bam1_2=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bam1_2,Nbins=1,Snorm=FALSE,key=keys1_2,norm_factor=norm1_2)
            range2only_bam1=bam1_2[-sh]       
          }else{
            #bam2 for range2, but not bam1 for range2
            bam1_2=NA
            range2only_bam1=NA
          }          
        }else{
          bam1_common=NA
          range1only_bam1=NA
          range2only_bam1=NA          
        }
        
      }else{
        bam1_common=NA
        range1only_bam1=NA
        range2only_bam1=NA
      }

      #if BAM2 is present and selected, proceed
      if(nchar(BAM2_choose)>0){
        pos_bam2_2=match(BAM2_choose,getbam2)
        bam2_2=rawvals2[[pos_bam2_2]]
        keys2_2=keyvals2[[pos_bam2_2]]
        norm2_2=normvals2[[pos_bam2_2]]
        #remove dups range2
        bam2_2=bam2_2[not_dup_range2]
        keys2_2=keys2_2[not_dup_range2]
        
        #decrypt and sum
        bam2_2=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bam2_2,Nbins=1,Snorm=FALSE,key=keys2_2,norm_factor=norm2_2)

        if(!is.null(bam2_2)){
          range2only_bam2=bam2_2[-sh]
          bam2_common=bam2_2[subjectHits(ov_range2)] 
          #if BAM2 is associated also with ROI 1,
          pos_bam2_1=match(BAM2_choose,getbam1)
          if(!is.na(pos_bam2_1) & length(pos_bam2_1)>0){
            bam2_1=rawvals1[[pos_bam2_1]]
            keys2_1=keyvals1[[pos_bam2_1]]
            norm2_1=normvals1[[pos_bam2_1]]
            #remove dups range1
            bam2_1=bam2_1[not_dup_range1]
            keys2_1=keys2_1[not_dup_range1]

            #decrypt and sum
            bam2_1=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=bam2_1,Nbins=1,Snorm=FALSE,key=keys2_1,norm_factor=norm2_1)
            range1only_bam2=bam2_1[-qh]

          } else{
            #bam1 for range1, but not bam2 for range1
            bam2_1=NA
            range1only_bam2=NA
          }          
        }else{
          range2only_bam2=NA
          range1only_bam2=NA
          bam2_common=NA          
        }
         
      }else{
        #no bam2 for range2. => no bam2 for range2 and range1
        range2only_bam2=NA
        range1only_bam2=NA
        bam2_common=NA
      }
      toplot$cmp$range1=range1
      toplot$cmp$range2=range2
      toplot$cmp$range1_notov=range1_notov
      toplot$cmp$range2_notov=range2_notov
      toplot$cmp$ov_range1=ov_range1
      toplot$cmp$ov_range2=ov_range2
      toplot$cmp$range2only_bam1=range2only_bam1
      toplot$cmp$range1only_bam1=range1only_bam1
      toplot$cmp$bam1_common=bam1_common
      toplot$cmp$range1only_bam2=range1only_bam2
      toplot$cmp$range2only_bam2=range2only_bam2
      toplot$cmp$bam2_common=bam2_common

      toplot$cmp$BAM1chooseCmp=BAM1_choose
      toplot$cmp$BAM2chooseCmp=BAM2_choose

      #make the boxplot with data available. If both bam files associated with
      #both selected ROIs, the boxplot will be complete
      #as common regions, the signal of BAM1 is the BAM1 signal of range1 inside
      #common areas (the average for each single overlapping area)
      #the same for BAM2: is the average of BAM2 signals inside ranges2 in common areas

      output$pairwiseoverlaps_box_options<-renderUI({
                list(
                  selectInput("chooseColorPaletteCmp_box","Choose color palette:",choices=c(
                                                  "red/grey"="red_gray_red4_grey20",
                                                  "blue/green"="blue_green_blue4_green4"
                                                )),
                  radioButtons("chooseNormalizationCmp_box","Choose normalization:",
                                      choiceNames=list(
                                        htmlhelp("Total reads","help_singleEvaluation_normalizationtotalread_cmpbox"),
                                        htmlhelp("Read density (reads/bp)","help_singleEvaluation_normalizationreaddensity_cmpbox")
                                      ),
                                      choiceValues=list(
                                        "totread",
                                        "readdensity"
                                      )
                              ),
                  checkboxInput("islogCmp_box", label="log2",value = FALSE, width = NULL)
                )
      })      
      #here put specific options for the boxplot and the boxplot
      output$viewBoxplotCmp<-renderPlot({

        #if not present, start with default colors
        if (isvalid(input$chooseColorPaletteCmp_box)){
          Colors=strsplit(input$chooseColorPaletteCmp_box,split="_")[[1]]
        }else{
          Colors=strsplit("red_gray_red4_grey20",split="_")[[1]]
        }
        #if not present, start with default normalization
        if (isvalid(input$chooseNormalizationCmp_box)){
          normalization=input$chooseNormalizationCmp_box
        }else{
          normalization=strsplit("totread",split="_")[[1]]
        } 
        if (isvalid(input$islogCmp_box)){
          islog=input$islogCmp_box
        }else{
          islog=FALSE
        }       

        if(normalization=="readdensity"){
          range1only_bam1=range1only_bam1/width(range1_notov)
          width_range1_common=width(range1[subjectHits(ov_range1)])
          bam1_common=bam1_common/width_range1_common
          range2only_bam1=range2only_bam1/width(range2_notov)
          
          range2only_bam2=range2only_bam2/width(range2_notov)
          width_range2_common=width(range2[subjectHits(ov_range2)])
          bam2_common=bam2_common/width_range2_common
          range1only_bam2=range1only_bam2/width(range1_notov)  
        }

        if (!all(is.na(bam1_common))){
          common_bam1=tapply(bam1_common,queryHits(ov_range1),mean)
        }else{
          common_bam1=NA
        }
        if (!all(is.na(bam2_common))){
          common_bam2=tapply(bam2_common,queryHits(ov_range2),mean)
        }else{
          common_bam2=NA
        }         

        plot_cmp_box(islog=islog,normalization=normalization,range2only_bam1=range2only_bam1,
                    range1only_bam1=range1only_bam1,common_bam1=common_bam1,range1only_bam2=range1only_bam2,range2only_bam2=range2only_bam2,
                    common_bam2=common_bam2,n1=n1,n2=n2,BAM1_choose=BAM1_choose,BAM2_choose=BAM2_choose,colors=Colors)
      })


      #pdf button for boxplot
      output$saveviewBoxplotCmp=renderUI({downloadButton('saveviewBoxplotCmpbutton', 'Get PDF')})
      #now the download data button will appear
      output$saveboxdataCmp=renderUI({downloadButton('saveenrichmentBoxCmpdata', 'Save data')})



      quantiles=seq(0,0.99,0.05) #will be x coordinates for the final plot
      #here, at least bam1_1 or bam2_2 must exist by definition
      if(exists("bam1_1") & !exists("bam2_2") ){
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

      #options for calibration plot
      output$pairwiseoverlaps_calibration_options<-renderUI({
                  selectInput("chooseColorPaletteCmp_calibration","Choose color palette:",choices=c(
                                                  "red/grey"="red_gray_red4_grey20",
                                                  "blue/green"="blue_green_blue4_green4"
                                                ))
      })
      #plot the calibration of overlap
      output$viewCalibrationCmp<-renderPlot({
        if (isvalid(input$chooseColorPaletteCmp_calibration)){
          Colors=strsplit(input$chooseColorPaletteCmp_calibration,split="_")[[1]]
        }else{
          Colors=strsplit("red_gray_red4_grey20",split="_")[[1]]
        }
        
        plot_cmp_calibration(quantiles=quantiles,only1=only1,only2=only2,combined1=combined1,combined2=combined2,
                              n1=n1,n2=n2,colors=Colors)
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
      if(all(!is.na(range2only_bam1)) & all(!is.na(range1only_bam1)) & all(!is.na(range2only_bam2)) & all(!is.na(range1only_bam2)) 
        & all(!is.na(bam1_common)) & all(!is.na(bam2_common)) ){
         

        #######################################################################
        #creation of df for scatterplot, if both bam files are present 
        common_block=unique(queryHits(ov_range1))
        df=as.data.frame(matrix( rep(NA,(length(range2only_bam2)+length(range1only_bam1)+length(common_block))*7),ncol=7 ))
        labs=c(rep("2only",length(range2only_bam2)), rep("1only",length(range1only_bam1)),rep("common",length(common_block)) )
        df[,3]=labs
        colnames(df)=c("bam1","bam2","label","width_range1","width_range2","times_factor_hits_range1","times_factor_hits_range2")

  

        #fill df with enrichment values with the proper belonging subset
        if(length(range2only_bam1)>0){
          df[1:length(range2only_bam2),1]=range2only_bam1
          df[1:length(range2only_bam2),2]=range2only_bam2 
          df[1:length(range2only_bam2),5]=width(range2_notov)
        }
        if(length(range1only_bam1)>0){
          df[(length(range2only_bam2)+1):(length(range2only_bam2)+length(range1only_bam2)),1]=range1only_bam1
          df[(length(range2only_bam2)+1):(length(range2only_bam2)+length(range1only_bam2)),2]=range1only_bam2
          df[(length(range2only_bam2)+1):(length(range2only_bam2)+length(range1only_bam2)),4]=width(range1_notov)       
        }
        if(length(common_block)>0){
          width_range1_common=width(range1[subjectHits(ov_range1)])
          width_range2_common=width(range2[subjectHits(ov_range2)])
          df[(nrow(df)-length(common_block)+1):nrow(df),1]=tapply(bam1_common,queryHits(ov_range1),sum)
          df[(nrow(df)-length(common_block)+1):nrow(df),4]=tapply(width_range1_common,queryHits(ov_range1),sum)
          df[(nrow(df)-length(common_block)+1):nrow(df),6]=tapply(bam1_common,queryHits(ov_range1),length)
          df[(nrow(df)-length(common_block)+1):nrow(df),2]=tapply(bam2_common,queryHits(ov_range2),sum) 
          df[(nrow(df)-length(common_block)+1):nrow(df),5]=tapply(width_range2_common,queryHits(ov_range2),sum)  
          df[(nrow(df)-length(common_block)+1):nrow(df),7]=tapply(bam2_common,queryHits(ov_range2),length)     
        }
        #######################################################################



        toplot$cmp$completedf=df
        #random sample df if more than 5.000
        if (nrow(df)>5000){
          set.seed(123) # set the seed, so the result will be the same every run with same parameters
          df=df[sort(sample(nrow(df), 5000,replace=FALSE)), ]
        }

        toplot$cmp$df=df
        # toplot$cmp$cor_common=cor_common
        # toplot$cmp$cor_all=cor_all
        

        #here try to put the checkbox input with the choice of what to show in the scatterplot
        listscatterchoice=list("exclusive ROI-1 ranges"="1only","exclusive ROI-2 ranges"="2only","common ranges"="common")
        

        output$pairwiseoverlaps_scatter_options<-renderUI({
                  list(
                    selectInput("chooseColorPaletteCmp_scatter","Choose color palette:",choices=c(
                                                    "red/grey"="red_gray_red4_grey20",
                                                    "blue/green"="blue_green_blue4_green4"
                                                  )),
                    radioButtons("chooseNormalizationCmp_scatter","Choose normalization:",
                                      choiceNames=list(
                                        htmlhelp("Total reads","help_singleEvaluation_normalizationtotalread_cmpscatter"),
                                        htmlhelp("Read density (reads/bp)","help_singleEvaluation_normalizationreaddensity_cmpscatter")
                                      ),
                                      choiceValues=list(
                                        "totread",
                                        "readdensity"
                                      )
                                                        ),
                    checkboxInput("islogCmp_scatter", label="log2",value = FALSE, width = NULL),
                    checkboxGroupInput("scatterplotChoice",label="Show only:",
                                        choices=listscatterchoice,selected=unlist(listscatterchoice) )
                  )
        }) 
        output$viewScatterplotCmp<-renderPlot({
          if (isvalid(input$chooseColorPaletteCmp_scatter)){
            Colors=strsplit(input$chooseColorPaletteCmp_scatter,split="_")[[1]]
          }else{
            Colors=strsplit("red_gray_red4_grey20",split="_")[[1]]
          }

          if (isvalid(input$chooseNormalizationCmp_scatter)){
            normalization=input$chooseNormalizationCmp_scatter
          }else{
            normalization="totread"
          }

          if(isvalid(input$scatterplotChoice)){
            scatterplotchoice=input$scatterplotChoice
          }else{
            scatterplotchoice=c("1only","2only","common")
          }

          if(isvalid(input$islogCmp_scatter)){
            islog=input$islogCmp_scatter
          }else{
            islog=FALSE
          }
          

          plot_cmp_scatter(df=df,normalization=normalization,subsetToShow=scatterplotchoice,islog=islog,
                    BAM2chooseCmp=BAM2_choose,BAM1chooseCmp=BAM1_choose,
                    n1=n1,n2=n2,colors=Colors,insets=c(-0.75,-1.2))
        })
        #button for PDF of the scatterplot
        output$saveviewScatterplotCmp=renderUI({downloadButton('saveviewScatterplotCmpbutton', 'Get PDF')})
        #button to download data about the scatterplot
        output$saveScatterdataCmp=renderUI({downloadButton('saveenrichmentScatterCmpdata', 'Save data')})

      }else{
        output$viewScatterplotCmp<-renderPlot({plot_text(text="you need to associate\ntwo enrichment files for both ROIs",cex=1.4)})
        output$pairwiseoverlaps_scatter_options<-renderUI({NULL})
        output$saveScatterdataCmp=renderUI({NULL})
        output$saveviewScatterplotCmp=renderUI({NULL})
      }
        
    }else{
      #plot enrichment stuff (box & scatter) are NULL
      output$viewCalibrationCmp<-renderPlot({plot_text(text="you need to associate at least\nan enrichment file to one of\nthe two ROIs",cex=1.4)})
      output$viewBoxplotCmp<-renderPlot({plot_text(text="you need to associate at least\nan enrichment file to one of\nthe two ROIs",cex=1.4)})
      output$viewScatterplotCmp<-renderPlot({plot_text(text="you need to associate\ntwo enrichment files for both ROIs",cex=1.4)})
      output$saveboxdataCmp=renderUI({NULL})
      output$saveviewBoxplotCmp=renderUI({NULL})
      output$saveScatterdataCmp=renderUI({NULL})
      output$saveviewScatterplotCmp=renderUI({NULL})
      output$saveviewCalibrationCmp=renderUI({NULL})
      output$pairwiseoverlaps_box_options<-renderUI({NULL})
      output$pairwiseoverlaps_scatter_options<-renderUI({NULL})
      output$pairwiseoverlaps_calibration_options<-renderUI({NULL})
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
    output$saveScatterdataCmp=renderUI({NULL})
    output$saveviewScatterplotCmp=renderUI({NULL})
    output$saveviewCalibrationCmp=renderUI({NULL})
    output$pairwiseoverlaps_overlap_options<-renderUI({NULL})
    output$pairwiseoverlaps_box_options<-renderUI({NULL})
    output$pairwiseoverlaps_scatter_options<-renderUI({NULL})
    output$pairwiseoverlaps_calibration_options<-renderUI({NULL})
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
    Colors=strsplit(input$chooseColorPaletteCmp_overlap,split="_")[[1]]
    #for PDF, better to split legend to next page
    plot_cmp_barplot(matbar=toplot$cmp$matbar,n1=toplot$cmp$n1,n2=toplot$cmp$n2,colors=Colors,
                  jaccard=toplot$cmp$jaccard,splitlegend=TRUE)
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
    Colors=strsplit(input$chooseColorPaletteCmp_overlap,split="_")[[1]]
    plot_cmp_venn(toplot$cmp$area1,toplot$cmp$area2,toplot$cmp$common,toplot$cmp$n1,toplot$cmp$n2,Colors)
    dev.off()
  } 
)



#observer to get PDF button for boxplot cmp
output$saveviewBoxplotCmpbutton<- downloadHandler(
  filename=function() {
      paste('Boxplot_enrichment_overlap.pdf', sep='')
  },
  content=function(file) {
        Colors=strsplit(input$chooseColorPaletteCmp_box,split="_")[[1]]
        range1only_bam1=toplot$cmp$range1only_bam1
        ov_range1=toplot$cmp$ov_range1
        range1_notov=toplot$cmp$range1_notov
        range2_notov=toplot$cmp$range2_notov
        range2only_bam1=toplot$cmp$range2only_bam1
        ov_range2=toplot$cmp$ov_range2
        range2only_bam2=toplot$cmp$range2only_bam2
        range1only_bam2=toplot$cmp$range1only_bam2
        range1=toplot$cmp$range1
        range2=toplot$cmp$range2
        bam1_common=toplot$cmp$bam1_common
        bam2_common=toplot$cmp$bam2_common
        

        if(input$chooseNormalizationCmp_box=="readdensity"){
          lab_axis="Read density (reads/bp)"
          range1only_bam1=range1only_bam1/width(range1_notov)
          width_range1_common=width(range1[subjectHits(ov_range1)])
          bam1_common=bam1_common/width_range1_common
          range2only_bam1=range2only_bam1/width(range2_notov)
          
          range2only_bam2=range2only_bam2/width(range2_notov)
          width_range2_common=width(range2[subjectHits(ov_range2)])
          bam2_common=bam2_common/width_range2_common
          range1only_bam2=range1only_bam2/width(range1_notov)  
        }

        if (!all(is.na(bam1_common))){
          common_bam1=tapply(bam1_common,queryHits(ov_range1),mean)
        }else{
          common_bam1=NA
        }
        if (!all(is.na(bam2_common))){
          common_bam2=tapply(bam2_common,queryHits(ov_range2),mean)
        }else{
          common_bam2=NA
        } 

        pdf(file)
        plot_cmp_box(islog=input$islogCmp_box,normalization=input$chooseNormalizationCmp_box,range2only_bam1=range2only_bam1,
                    range1only_bam1=range1only_bam1,common_bam1=common_bam1,range1only_bam2=range1only_bam2,range2only_bam2=range2only_bam2,
                    common_bam2=common_bam2,n1=toplot$cmp$n1,n2=toplot$cmp$n2,BAM1_choose=toplot$cmp$BAM1chooseCmp,BAM2_choose=toplot$cmp$BAM2chooseCmp,colors=Colors)

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
        Colors=strsplit(input$chooseColorPaletteCmp_calibration,split="_")[[1]]
        plot_cmp_calibration(quantiles=toplot$cmp$quantiles,only1=toplot$cmp$only1,only2=toplot$cmp$only2,
                          combined1=toplot$cmp$combined1,combined2=toplot$cmp$combined2,n1=toplot$cmp$n1,
                          n2=toplot$cmp$n2,insets=c(0,-0.5),colors=Colors)
        dev.off()
  } 
)


#observer for download data of boxplot
output$saveenrichmentBoxCmpdata<- downloadHandler(
  filename=function() {
      paste('boxplot_data.xls', sep='')
  },
  content=function(file) {
      range1only_bam1=toplot$cmp$range1only_bam1
      ov_range1=toplot$cmp$ov_range1
      range1_notov=toplot$cmp$range1_notov
      range2_notov=toplot$cmp$range2_notov
      range2only_bam1=toplot$cmp$range2only_bam1
      ov_range2=toplot$cmp$ov_range2
      range2only_bam2=toplot$cmp$range2only_bam2
      range1only_bam2=toplot$cmp$range1only_bam2
      range1=toplot$cmp$range1
      range2=toplot$cmp$range2
      bam1_common=toplot$cmp$bam1_common
      bam2_common=toplot$cmp$bam2_common

      if(input$chooseNormalizationCmp_box=="readdensity"){
        range1only_bam1=range1only_bam1/width(range1_notov)
        width_range1_common=width(range1[subjectHits(ov_range1)])
        bam1_common=bam1_common/width_range1_common
        range2only_bam1=range2only_bam1/width(range2_notov)
        
        range2only_bam2=range2only_bam2/width(range2_notov)
        width_range2_common=width(range2[subjectHits(ov_range2)])
        bam2_common=bam2_common/width_range2_common
        range1only_bam2=range1only_bam2/width(range1_notov)  
      }

      if (!all(is.na(bam1_common))){
        common_bam1=tapply(bam1_common,queryHits(ov_range1),mean)
      }else{
        common_bam1=NA
      }
      if (!all(is.na(bam2_common))){
        common_bam2=tapply(bam2_common,queryHits(ov_range2),mean)
      }else{
        common_bam2=NA
      } 

      if(input$islogCmp_box){
        arraytoplot=list(log2(range2only_bam1),log2(range1only_bam1),log2(common_bam1),
                log2(range1only_bam2),log2(range2only_bam2),log2(common_bam2))
      }else{
        arraytoplot=list(range2only_bam1,range1only_bam1,common_bam1,
                range1only_bam2,range2only_bam2,common_bam2)
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
        Colors=strsplit(input$chooseColorPaletteCmp_scatter,split="_")[[1]]
        pdf(file)
        plot_cmp_scatter(df=df,normalization=input$chooseNormalizationCmp_scatter,subsetToShow=input$scatterplotChoice,islog=input$islogCmp_scatter,
                    BAM2chooseCmp=toplot$cmp$BAM2chooseCmp,BAM1chooseCmp=toplot$cmp$BAM1chooseCmp,
                    n1=toplot$cmp$n1,n2=toplot$cmp$n2,colors=Colors)
        dev.off()
  } 
)



#respond to download data about scatterplot if present:
output$saveenrichmentScatterCmpdata<- downloadHandler(
  filename=function() {
      paste('scatter_data.xls', sep='')
  },
  content=function(file) {
      completedf=toplot$cmp$completedf
      if(input$chooseNormalizationCmp_scatter=="readdensity"){
        #determine 1 and 2 columns of df based on width. For common ranges, divide sums of signals and width
        completedf[completedf$label=="2only",1]=completedf[completedf$label=="2only",1]/completedf$width_range2[completedf$label=="2only"]
        completedf[completedf$label=="2only",2]=completedf[completedf$label=="2only",2]/completedf$width_range2[completedf$label=="2only"]
        completedf[completedf$label=="1only",1]=completedf[completedf$label=="1only",1]/completedf$width_range1[completedf$label=="1only"]
        completedf[completedf$label=="1only",2]=completedf[completedf$label=="1only",2]/completedf$width_range1[completedf$label=="1only"]
        completedf[completedf$label=="common",1]=completedf[completedf$label=="common",1]/completedf$width_range1[completedf$label=="common"]
        completedf[completedf$label=="common",2]=completedf[completedf$label=="common",2]/completedf$width_range2[completedf$label=="common"]

      }else{
        #here divide by times only, not for width
        completedf[completedf$label=="common",1]=completedf[completedf$label=="common",1]/completedf$times_factor_hits_range1[completedf$label=="common"]
        completedf[completedf$label=="common",2]=completedf[completedf$label=="common",2]/completedf$times_factor_hits_range2[completedf$label=="common"]

      }
      if(input$islogCmp_scatter){
        completedf[,1]=log2(completedf[,1])
        completedf[,2]=log2(completedf[,2])
      }

      #rename labels in a better way
      df=completedf
      df[df[,3]=="2only",][,3]=toplot$cmp$n2
      df[df[,3]=="1only",][,3]=toplot$cmp$n1

      #keep only values (remove widths etc. because already used for normalizations above)
      df=df[,1:3]
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


#observer help buttons parameters
observeEvent(input$help_digitalHeatmap_parameters_masterROI, { boxHelpServer(help_digitalHeatmap_parameters_masterROI)})
observeEvent(input$help_digitalHeatmap_parameters_ROItoview, { boxHelpServer(help_digitalHeatmap_parameters_ROItoview)})
observeEvent(input$help_digitalHeatmap_parameters_ROIordering, { boxHelpServer(help_digitalHeatmap_parameters_ROIordering)})
observeEvent(input$help_digitalHeatmap_parameters_ROIforcluster, { boxHelpServer(help_digitalHeatmap_parameters_ROIforcluster)})
observeEvent(input$help_digitalHeatmap_parameters_clusterKmeans, { boxHelpServer(help_digitalHeatmap_parameters_clusterKmeans)})
observeEvent(input$help_digitalHeatmap_parameters_clusterHierarchical, { boxHelpServer(help_digitalHeatmap_parameters_clusterHierarchical)})
observeEvent(input$help_digitalHeatmap_parameters_clusternumber, { boxHelpServer(help_digitalHeatmap_parameters_clusternumber)})
observeEvent(input$help_digitalHeatmap_parameters_nbins, { boxHelpServer(help_digitalHeatmap_parameters_nbins)})
observeEvent(input$help_digitalHeatmap_parameters_strandspecific, { boxHelpServer(help_digitalHeatmap_parameters_strandspecific)})
observeEvent(input$help_digitalHeatmap_parameters_randomsample, { boxHelpServer(help_digitalHeatmap_parameters_randomsample)})
observeEvent(input$help_digitalHeatmap_parameters_positionaloverlap, { boxHelpServer(help_digitalHeatmap_parameters_positionaloverlap)})




observeEvent(input$confirmUpdateDigitalHeat1,{
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
    rawvals=Enrichlist$rawcoverage[pos]
    keyvals=Enrichlist$decryptkey[pos]
    normvals=Enrichlist$normfactlist[pos]
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
      bigkeyvalslist=list()
      bignormlist=list()
      fixes=list()
      ranges_forclick=list()
      #loop through all master ROIs selected
      for(i in 1:length(roi)){
      	pos_notdup_range=!duplicated(getRange(roi[[i]]))
        unifROI=uniqueROI(roi[[i]])
        rangeroi=getRange(unifROI)
        bigrange[[i]]=granges(rangeroi)
        bigrange[[i]]$label=nameROI[i]

        ranges_forclick[[i]]=rangeroi


    		bigbamlist[[i]]=rawvals[[i]]
        bigkeyvalslist[[i]]=keyvals[[i]]
    		bamlist2=list()
        keylist2=list()
    		if(length(bigbamlist[[i]])>0){
    		  	for(j in 1:length(bigbamlist[[i]])){
    		    	element_i=bigbamlist[[i]][[j]]
    		    	#remove duplicated positions
    		    	bamlist2[[j]]=element_i[pos_notdup_range]
              element_j=bigkeyvalslist[[i]][[j]]
              keylist2[[j]]=element_j[pos_notdup_range]
    		    }		    
    		}
    		bigbamlist[[i]]=bamlist2
        bigkeyvalslist[[i]]=keylist2

        names(bigbamlist[[i]])=names(bigkeyvalslist[[i]])=names(rawvals[[i]])

        
        bignormlist[[i]]=normvals[[i]]
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

      bams=toplot$digital$ROIsForDigitalHeat
      lista=list()

      if(length(bams)>0){
        output$showoptioncolorsforDigitalHeat<-renderUI({
            radioButtons("optioncolorsforDigitalHeat",label="Select colors:",choiceNames=c("global color","custom colors"),
                            choiceValues=c("global","custom"),selected="global")
        }) 
 
      }else{
        #bams are NULL. show nothing
        output$showoptioncolorsforDigitalHeat<-renderUI({NULL})
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

        if (isvalid(input$optioncolorsforDigitalHeat)){
          optioncolors= input$optioncolorsforDigitalHeat
        }else{
          optioncolors="global"
        }

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
        if(isolate(optioncolors)=="global"){
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
            #if optioncolors is custom, palette should have length==nbams
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
      #set.seed(123)
      #color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      
      toplot$digital$clusternumbermatrix=NULL
      if(length(input$ROImaster) ==1 & length(input$ROIforClusteringDigitalHeat)>0 ){
        #color_distinct_cluster=sample(colors_list, toplot$digital$clustnumDigitalHeat)
        color_distinct_cluster=colors_list[1:toplot$digital$clustnumDigitalHeat]
        if(toplot$digital$clustnumDigitalHeat<=433 & toplot$digital$clustnumDigitalHeat>0){        
          if(!is.null(toplot$digital$clustering$clustobj)){

            #if custering will be drown, it means that we have 1 master ROI and clustering is ok.
            #therefore, keep info for subsequent subselection based on the click on cluster            
            #remember that this ROI must be re-annotated
            #1-define the starting "material", obtained at the beginning of the heat button (for the first
            #and the only master ROI selected):
            range_tokeep=ranges_forclick[[1]]
            bams_tokeep=bigbamlist[[1]]
            keyvals_tokeep=bigkeyvalslist[[1]]
            normfact_tokeep=bignormlist[[1]]
            fix_tokeep=fixes[[1]]
            #2-subsample them (1 single ROI)
            range_tokeep=range_tokeep[Digital_sample_pos]
            bams_tokeep=lapply(bams_tokeep,function(bamob){bamob[Digital_sample_pos]})
            keyvals_tokeep=lapply(keyvals_tokeep,function(kvalob){kvalob[Digital_sample_pos]})
            fix_tokeep=fix_tokeep[Digital_sample_pos]
            #3-re-ordering based on clustering (toplot$digital$clustering$clustobj$ord)
            range_tokeep=range_tokeep[toplot$digital$clustering$ord]
            bams_tokeep=lapply(bams_tokeep,function(bamob){bamob[toplot$digital$clustering$ord]})
            keyvals_tokeep=lapply(keyvals_tokeep,function(kvalob){kvalob[toplot$digital$clustering$ord]})
            fix_tokeep=fix_tokeep[toplot$digital$clustering$ord]
            #4-save in temporary variables to be used by clicking the clsuter of digital
            toplot$digital$range_tokeep=range_tokeep
            toplot$digital$bams_tokeep=bams_tokeep
            toplot$digital$keyvals_tokeep=keyvals_tokeep
            toplot$digital$fix_tokeep=fix_tokeep
            #add also norm factors (not subsampled nor reordered, but keep them)
            toplot$digital$normfact_tokeep=normfact_tokeep

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
            output$textNameClustDigitalHeat<-renderPlot({
              par(mar = rep(0, 4))
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i")     
              text(x=0.5,labels="clusters",y=0.5,srt=90,cex=1.4)           
            })  

          }else{
            output$clustersImageLeftDigital<-renderPlot({NULL})
            output$textNameClustDigitalHeat<-renderPlot({NULL})
            toplot$digital$color_distinct_cluster=NULL
            toplot$digital$clusternumbermatrix=NULL
            toplot$digital$range_tokeep=NULL
            toplot$digital$bams_tokeep=NULL
            toplot$digital$keyvals_tokeep=NULL
            toplot$digital$normfact_tokeep=NULL
            toplot$digital$fix_tokeep=NULL
          }
        }else{
          output$clustersImageLeftDigital<-renderPlot({NULL})
          output$textNameClustDigitalHeat<-renderPlot({NULL})
          toplot$digital$color_distinct_cluster=NULL
          toplot$digital$clusternumbermatrix=NULL
          toplot$digital$range_tokeep=NULL
          toplot$digital$bams_tokeep=NULL
          toplot$digital$keyvals_tokeep=NULL
          toplot$digital$normfact_tokeep=NULL
          toplot$digital$fix_tokeep=NULL
        }

      }else{
        output$clustersImageLeftDigital<-renderPlot({NULL})
        output$textNameClustDigitalHeat<-renderPlot({NULL})
        toplot$digital$clusternumbermatrix=NULL
        toplot$digital$color_distinct_cluster=NULL
        toplot$digital$range_tokeep=NULL
        toplot$digital$bams_tokeep=NULL
        toplot$digital$keyvals_tokeep=NULL
        toplot$digital$normfact_tokeep=NULL
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
        color_distinct_bias=colors_list[1: (length(freq_overlap)*length(freq_overlap[[1]])) ]
        #color_distinct_bias=sample(color, length(freq_overlap)*length(freq_overlap[[1]]))
        toplot$digital$color_distinct_bias=color_distinct_bias


        output$frequencyDigitalHeat_options<-renderUI({
          checkboxInput("FracToPercDigitalHeat", label=list("positional overlap %",htmlhelp("","help_digitalHeatmap_parameters_positionaloverlap")),value = TRUE, width = NULL)
        })


        output$frequencyOvDigitalHeat <- renderPlot({

          if(isvalid(input$FracToPercDigitalHeat)){
            frequencyperc=input$FracToPercDigitalHeat
          }else{
            frequencyperc=TRUE
          } 
          #interactively change freq/% if button is pressed:
          ovtoplot=freq_overlap
          if(frequencyperc){
            for(i in 1:length(freq_overlap)){
              ovtoplot[[i]]=lapply(freq_overlap[[i]],function(k){k/length_subsets[[i]]})
            }
            maxval=1
          }else{
            maxval=max(length_subsets)
          }

          plot_digital_frequency(ovtoplot=ovtoplot,nbin=nbin,fraction=frequencyperc,maxval=maxval,
                                  colors=color_distinct_bias)
        
        })
        output$savefrequencyOvDigitalHeat=renderUI({downloadButton('savefrequencyOvDigitalHeatbutton', 'Get PDF')})
      }else{
        output$frequencyOvDigitalHeat <- renderPlot({plot_text(text="The number of bins must be >1",cex=1.4)})
        output$savefrequencyOvDigitalHeat=renderUI({NULL})
        output$frequencyDigitalHeat_options<-renderUI({NULL})
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
          output$JaccardDigital<-renderPlot({plot_text(text="jaccard index plot will appear\nonly if a single master ROI is selected",cex=1.4)})
          output$saveJaccardDigital=renderUI({NULL})
          toplot$digital$matlistTotal=NULL
          #plot NULL, because only one ROI considered
        }
      }else{
        #plot NULL
        output$JaccardDigital<-renderPlot({plot_text(text="jaccard index plot will appear\nonly if a single master ROI is selected",cex=1.4)})
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
      output$showcolorsDigitalheat<-renderUI({NULL})
      output$saveheatmapDigital=renderUI({NULL})
      output$textNameDigitalHeat<-renderPlot({NULL}) 
      output$JaccardDigital<-renderPlot({NULL})  
      output$saveJaccardDigital=renderUI({NULL})
      output$showoptioncolorsforDigitalHeat<-renderUI({NULL})
      toplot$digital$matlistTotal=NULL 
      toplot$digital$clusternumbermatrix=NULL
      output$clustersImageLeftDigital<-renderPlot({NULL})
      output$textNameClustDigitalHeat<-renderPlot({NULL})
      toplot$digital$range_tokeep=NULL
      toplot$digital$bams_tokeep=NULL
      toplot$digital$keyvals_tokeep=NULL
      toplot$digital$normfact_tokeep=NULL
      toplot$digital$fix_tokeep=NULL
      output$newROIfromDigitalHeat_out<-renderUI({NULL})
      output$frequencyDigitalHeat_options<-renderUI({NULL})
      output$textfractionelementsDigitalHeat<-renderText({NULL})
      output$frequencyOvDigitalHeat <- renderPlot({NULL})
      output$savefrequencyOvDigitalHeat=renderUI({NULL})
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
    output$showcolorsDigitalheat<-renderUI({NULL})
    output$textNameDigitalHeat<-renderPlot({NULL})
    output$JaccardDigital<-renderPlot({NULL})
    output$saveJaccardDigital=renderUI({NULL})
    toplot$digital$matlistTotal=NULL
    toplot$digital$clusternumbermatrix=NULL
    output$clustersImageLeftDigital<-renderPlot({NULL})
    output$textNameClustDigitalHeat<-renderPlot({NULL})
    toplot$digital$range_tokeep=NULL
    toplot$digital$bams_tokeep=NULL
    toplot$digital$keyvals_tokeep=NULL
    toplot$digital$normfact_tokeep=NULL
    toplot$digital$fix_tokeep=NULL
    output$newROIfromDigitalHeat_out<-renderUI({NULL})
    output$frequencyDigitalHeat_options<-renderUI({NULL})
    output$textfractionelementsDigitalHeat<-renderText({NULL})
    output$showoptioncolorsforDigitalHeat<-renderUI({NULL})
    output$frequencyOvDigitalHeat <- renderPlot({NULL})
    output$savefrequencyOvDigitalHeat=renderUI({NULL})
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
          }
          else{
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


      if(!is.null(toplot$digital$clusternumbermatrix)){
        par(mar = c(20,12,0,0))
        image(0:nrow(toplot$digital$clusternumbermatrix), 0:ncol(toplot$digital$clusternumbermatrix),toplot$digital$clusternumbermatrix,col=toplot$digital$color_distinct_cluster,axes = FALSE, xlab = "", ylab = "",ylim=c(0,ncol(trasp)))
        axis( 2,at=pos_vertical,labels=paste(rev(toplot$digital$ROImaster),"\n",labels_vert),las=1,tick=FALSE)
        axis(1,0.5,labels="cluster",las= 2,tick=FALSE)
        par(mar = c(20,0,0,0))
        image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
        #REMOVE x,y lim if cordinates don't match
              ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp))   
              )            
      }else{
        par(mar = c(20,12,0,0))
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
    }else{
      maxval=max(toplot$digital$length_subsets)
    }

    pdf(file)

    plot_digital_frequency(ovtoplot=ovtoplot,nbin=toplot$digital$nbin,fraction=input$FracToPercDigitalHeat,
                              maxval=maxval,colors=toplot$digital$color_distinct_bias,splitlegend=TRUE)
  
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
          keyvals_sel=lapply(toplot$digital$keyvals_tokeep,function(bamob){bamob[toplot$gadgetdigital$ybottom:toplot$gadgetdigital$ytop]})
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
          Enrichlist$rawcoverage[[input$newROIfromDigitalHeat]]=bams_sel
          Enrichlist$decryptkey[[input$newROIfromDigitalHeat]]=keyvals_sel
          Enrichlist$normfactlist[[input$newROIfromDigitalHeat]]=toplot$digital$normfact_tokeep
          ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                name=input$newROIfromDigitalHeat,
                                range=range_sel,
                                fixed=fix_sel,
                                BAMlist=list(),
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


#observers for all help buttons of parameters
observeEvent(input$help_analogicHeatmap_parameters_ROI, {boxHelpServer(help_analogicHeatmap_parameters_ROI)})
observeEvent(input$help_analogicHeatmap_parameters_enrichments, {boxHelpServer(help_analogicHeatmap_parameters_enrichments)})
observeEvent(input$help_analogicHeatmap_parameters_enrichmentorder, {boxHelpServer(help_analogicHeatmap_parameters_enrichmentorder)})
observeEvent(input$help_analogicHeatmap_parameters_ranking, {boxHelpServer(help_analogicHeatmap_parameters_ranking)})
observeEvent(input$help_analogicHeatmap_parameters_clustering, {boxHelpServer(help_analogicHeatmap_parameters_clustering)})
observeEvent(input$help_analogicHeatmap_parameters_clusterKmeans, {boxHelpServer(help_analogicHeatmap_parameters_clusterKmeans)})
observeEvent(input$help_analogicHeatmap_parameters_clusterHierarchical, {boxHelpServer(help_analogicHeatmap_parameters_clusterHierarchical)})
observeEvent(input$help_analogicHeatmap_parameters_clusternumber, {boxHelpServer(help_analogicHeatmap_parameters_clusternumber)})
observeEvent(input$help_analogicHeatmap_parameters_nbins, {boxHelpServer(help_analogicHeatmap_parameters_nbins)})
observeEvent(input$help_analogicHeatmap_parameters_subsample, {boxHelpServer(help_analogicHeatmap_parameters_subsample)})
observeEvent(input$help_analogicHeatmap_parameters_uniform, {boxHelpServer(help_analogicHeatmap_parameters_uniform)})
observeEvent(input$help_analogicHeatmap_parameters_individual, {boxHelpServer(help_analogicHeatmap_parameters_individual)})


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



# toListenAnalogHeat <- reactive({
#     list(input$confirmUpdateAnalogHeat,input$confirmUpdateAnalogHeat2)
# })
observeEvent(input$confirmUpdateAnalogHeat,{
  
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
  ################################################################################################








  if (length(ROIvariables$listROI)>0 & length(input$ROIsForAnalogHeat)>0 & length(toplot$analogic$BAMsForAnalogHeat)>0 & input$binsAnalogHeat>0 & samplerandom>0 & isvalid(input$sampleRandomAnalogHeat) & isvalid(input$binsAnalogHeat)){
   
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

    nomi=unlist(lapply(ROIvariables$listROI,getName))
    pos=match(input$ROIsForAnalogHeat,nomi)
    #all rois selected
    roi=ROIvariables$listROI[pos]

    rawvals=Enrichlist$rawcoverage[pos]
    keyvals=Enrichlist$decryptkey[pos]
    normvals=Enrichlist$normfactlist[pos]

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
      normlisttouse=list() #for usage, same as bigbamlist
      normlisttostore=list() #for extraction, same as BAMs
      keylisttouse=list()
      keylisttostore=list()
      completerange=list()
      for(i in 1:length(roi)){
        #may cause problems: slow and when re-build ROI from heatmap selection,
        #strands "-" are inverted, while should not be...
        #unifROI=unifyStrand(roi[[i]])
        unifROI=roi[[i]]
        pos_notdup_range=!duplicated(getRange(unifROI))

        unifROI=uniqueROI(unifROI)
        fixes[[i]]=granges(getFixed(unifROI))
        completerange[[i]]=getRange(unifROI)
        #we have to collapse => all ranges must have the same columns
        bigrange[[i]]=granges(completerange[[i]])
        bigrange[[i]]$label=nameROI[i]      
        
        #when collapsing ranges (for example alternative promoters identical),
        #we must also collapse enrichments. uniqueROI function works only for ranges and fixes
        #because BAMlist attribute has been removed in the last version. for enrichments, we removed duplicated positions 
        #later
        totbams=rawvals[[i]]
        totkeys=keyvals[[i]]
    		bamlist2=list()
        keyvals2=list()
    		if(length(totbams)>0){
    		  	for(j in 1:length(totbams)){
    		    	element_i=totbams[[j]]
    		    	#remove duplicated positions
    		    	bamlist2[[j]]=element_i[pos_notdup_range]
              element_k=totkeys[[j]]
              keyvals2[[j]]=element_k[pos_notdup_range]
    		    }		    
    		}
    		totbams=bamlist2
        totkeys=keyvals2
    		BAMs[[i]]=totbams
        keylisttostore[[i]]=totkeys
        normlisttostore[[i]]=normvals[[i]]
        #order even if BAM reordered, should be fine
        bamnames[[i]]=names(rawvals[[i]])
        pos2=match(toplot$analogic$BAMsForAnalogHeat,bamnames[[i]])
        bigbamlist[[i]]=totbams[pos2]
        keylisttouse[[i]]=totkeys[pos2]
        normlisttouse[[i]]=normvals[[i]][pos2]
        
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
              finalbam=finalkey=list()
              for(i in 1:length(BAMs)){
                finalbam[[i]]=finalkey[[i]]=as.list(rep(NA,length(BAMs[[i]])))
              }
              cumulatesum=0
              #for each ROI
              for(i in 1:length(BAMs)){
                #for each BAM, subselect only for those selected in that label:
                #it's a RELATIVE position from the beginning of that label, so subtract what is previous
                for (k in 1:length(BAMs[[i]])){
                  finalbam[[i]][[k]]=BAMs[[i]][[k]][positions_splitted[[i]]-cumulatesum]
                  finalkey[[i]][[k]]=keylisttostore[[i]][[k]][positions_splitted[[i]]-cumulatesum]
                }
                names(finalbam[[i]])=names(finalkey[[i]])=bamnames[[i]]
                #the same of ROIs (complete, if have different columns):
                completerange[[i]]=completerange[[i]][positions_splitted[[i]]-cumulatesum]
                fixes[[i]]=fixes[[i]][positions_splitted[[i]]-cumulatesum]
                cumulatesum=cumulatesum+length(BAMs[[i]][[1]])
              }
              names(finalbam)=names(finalkey)=input$ROIsForAnalogHeat

              #here, the normalizations should keep the same as before: normlisttostore, because they 
              #do not undergo random sample (one single value for each enrichment!)
              #split also finalrange_sampled and fixes according to "labelstosplit"
              # finalrange_sampled_split=split( finalrange_sampled ,factor(labelstosplit,levels=unique(labelstosplit)))
              # fixes_sampled_split=split( fixes ,factor(labelstosplit,levels=unique(labelstosplit)))          
              ############################################################

              #print("prepared...")
              #combine lists from different ROIs (c)
              subselectedBAMlist=list()
              subselectkeylist=list()
              for(i in 1:length(bigbamlist[[1]])){
                provv=list()
                provvkey=c()
                #for each ROI
                for (k in 1:length(bigbamlist)){
                  provv=c(provv,bigbamlist[[k]][[i]])
                  provvkey=c(provvkey,keylisttouse[[k]][[i]])
                }   
                #subselect "provv": provv is the big union of all ROI for a specific BAM
                provv=provv[finalBAMs_sample_pos]
                provvkey=provvkey[finalBAMs_sample_pos]
                #decrypt (now light because already subsampled)
                #provv=decryptcov( list(provv,provvkey),chunk=length(provv))
                subselectkeylist[[i]]=provvkey
                subselectedBAMlist[[i]]=provv
              }
              #use only finalBAMs_sample_pos positions sample. This is because the "random" sample
              #number refers to the total number of rows plotted, and not for each range.
              #this means that we have to merge matrixes first and then subsample. This will keep 
              #the proportion between ROIs intact.


              #re-split the matlists according to the correct label (sampled with the same position)
              #ROIs must be in the higher order hierarchy of list 
              slicedbamlist=as.list(rep(NA,length(toplot$analogic$BAMsForAnalogHeat)))
              slicedbamlist=rep(list(slicedbamlist),length(input$ROIsForAnalogHeat))
              slicedkeylist=slicedbamlist
              for(i in 1:length(subselectedBAMlist)){
                slicedbamlist2=split(subselectedBAMlist[[i]],factor(labelstosplit,levels=unique(labelstosplit)))
                slicedkeylist2=split(subselectkeylist[[i]],factor(labelstosplit,levels=unique(labelstosplit)))
                #for each ROI obtained after split
                for(k in 1:length(slicedbamlist2)){
                  slicedbamlist[[k]][[i]]=slicedbamlist2[[k]]
                  slicedkeylist[[k]][[i]]=slicedkeylist2[[k]]
                }      
              }
              #here, check if normlisttouse follows the same hierarchy of slicedbamlist (ROI/enrich)
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
                                        return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE,norm_factor=normlisttouse[[i]][[k]],key=slicedkeylist[[i]][[k]]))
                                      })
                  }else{
                    decision=0
                    tryCatch({
                      matlist=mclapply(1:length(slicedbamlist[[i]]),function(k) {
                       
                                        return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE,norm_factor=normlisttouse[[i]][[k]],key=slicedkeylist[[i]][[k]]))
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
                                        return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE,norm_factor=normlisttouse[[i]][[k]],key=slicedkeylist[[i]][[k]]))
                                      })
                    }
                  }

                  ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                  ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
                  #pos_toinvert=as.character(strand_sampled[[i]])=="-"
                  #for k in matlist (k should be the bam)
                  # for(k in 1:length(matlist)){
                  #   #for each of the BAM file in this current ROI;
                  #   matprovv=matlist[[k]]
                  #   if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
                  #     matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
                  #   }
                  # }


                  return(matlist)    
                })
              }else{

                #if ROIs>BAMs, parallelize at ROI level
                if(nc==1){
                  matlists=lapply(1:length(slicedbamlist),function(i) {
                    #in case, filter for transcript flag if mem problems
                    matlist=lapply(1:length(slicedbamlist[[i]]),function(k) {
                                          return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE,norm_factor=normlisttouse[[i]][[k]],key=slicedkeylist[[i]][[k]]))
                                        }) 
                    ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                    ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
                    return(matlist)    
                  })                
                }else{
                  decision=0
                  tryCatch({
                    matlists=mclapply(1:length(slicedbamlist),function(i) {
                      #in case, filter for transcript flag if mem problems
                      matlist=lapply(1:length(slicedbamlist[[i]]),function(k) {
                        
                                            return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE,norm_factor=normlisttouse[[i]][[k]],key=slicedkeylist[[i]][[k]]))
                                          }) 
                      ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                      ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
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
                                            return(makeMatrixFrombaseCoverage(slicedbamlist[[i]][[k]],Nbins=nbintouse,Snorm=TRUE,norm_factor=normlisttouse[[i]][[k]],key=slicedkeylist[[i]][[k]]))
                                          }) 
                      ###HERE invert "-" strand. They are in strand_sampled list. each element of this list
                      ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
                      return(matlist)    
                    }) 

                  }
                }
              }


              #return matlists in any case. If fixed size is the same for all, return complete matrix

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
                  if(!inherits(input$distmethodAnalogHeat,"character")  |!inherits(input$clustmethodAnalogHeat,"character")  ){
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
                  finalkey[[i]][[k]]=finalkey[[i]][[k]][toorder_final[[i]]]
                }
              }
              
              toplot$analogic$finalrange=completerange
              toplot$analogic$finalbam=finalbam
              toplot$analogic$finalkey=finalkey
              toplot$analogic$finalnorm=normlisttostore
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
              toplot$analogic$matrixes_processed=matrixes_processed
              toplot$analogic$labelstosplit=labelstosplit
              toplot$analogic$clusterTypeAnalogHeat=input$clusterTypeAnalogHeat

              #here, put the code to generate clustering number arrays to plot
              #and set a variable "clustering" if all conditions are met, to draw the side colors
              isclusteringok= input$chooseOrderingAnalogHeat=="clustering" & length(input$ROIsForAnalogHeat) ==1 & length(input$BAMsForClusteringAnalogHeat)>0 &
                              !is.na(heatvariables$clustnumAnalogHeat) & heatvariables$clustnumAnalogHeat<=433 & heatvariables$clustnumAnalogHeat>0
              
              toplot$analogic$clusternumbermatrix=NULL
              if(isclusteringok){
                #set.seed(123)
                #color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
                #color_distinct_cluster=sample(color, heatvariables$clustnumAnalogHeat)
                color_distinct_cluster=colors_list[1:heatvariables$clustnumAnalogHeat]

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



              #plotting elments. 

              #plot cluster color heatmap on the left;
              if(isclusteringok){
                output$clustersImageLeft<-renderPlot({
                  par(mar = c(0,0,0,0))
                  image(0:nrow(matrix_like), 0:ncol(matrix_like),matrix_like,col=color_distinct_cluster,axes = FALSE, xlab = "", ylab = "",ylim=c(0,ncol(matrix_like)))
                }) 
                output$textNameClustAnalogHeat<-renderPlot({
                  par(mar = rep(0, 4))
                  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n',xaxs="i")     
                  text(x=0.5,labels="clusters",y=0.5,srt=90,cex=1.4)           
                })             
              }else{
                output$clustersImageLeft<-renderPlot({NULL}) 
                output$textNameClustAnalogHeat<-renderPlot({NULL})
              }




              #save heatmap data xls button
              output$showsaveAnalogHeatdata=renderUI({downloadButton('saveAnalogHeatdata', 'Download matrix data')})
              #clear button for new ROI: we have another analysis
              output$newROIfromAnalogHeat_out<-renderUI({NULL})



              #reset the color menu:

              bams=heatvariables$BAMsForAnalogHeat
              lista=list()

              if(length(bams)>0){
                output$showoptioncolorsforAnalogHeat<-renderUI({
                  radioButtons("optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("global color","custom colors"),
                            choiceValues=c("global","custom"),selected="global")
                })
                output$showcolorsheat<-renderUI({
                  selectInput("colorCustomAnalogHeat_global",label="Choose global color:",c("white/red"="white_red4",
                                                          "white/blue"="white_blue",
                                                          "rainbow"="rainbow",
                                                          "blue/white/red"="blue_white_red",
                                                          "exponential blue"="white_white_white_blue_blue4",
                                                          "gray scale"="white_grey90_grey80_grey70_grey50_black"))
                })

                output$showchooseQuantileMethodAnalogHeat<-renderUI({
                  list(
                    HTML("<b>Quantile threshold</b>"),
                    radioButtons("chooseQuantileMethodAnalogHeat",label=NULL,
                        choiceNames=list(
                          htmlhelp("Uniform","help_analogicHeatmap_parameters_uniform"),
                          htmlhelp("Individual","help_analogicHeatmap_parameters_individual")
                        ),
                        choiceValues=list("allBAM","eachBAM")
                    )
                  )
                })
                output$showquantileThreshAnalogHeat<-renderUI({
                  sliderInput('quantileThreshAnalogHeat',label=NULL,min = 0.1, max = 1, value = 0.9,step=0.002)
                })

              }else{
                #bams are NULL. show nothing
                output$showchooseQuantileMethodAnalogHeat<-renderUI({NULL})
                output$showquantileThreshAnalogHeat<-renderUI({NULL})  
                output$showoptioncolorsforAnalogHeat<-renderUI({NULL})
              }



              
              output$heatmapAnalog<-renderPlot({

                if (isvalid(isolate(input$optioncolorsforAnalogHeat))){
                  optioncolors= isolate(input$optioncolorsforAnalogHeat)
                }else{
                  optioncolors="global"
                }

                if (isvalid(input$chooseQuantileMethodAnalogHeat)){
                  method=input$chooseQuantileMethodAnalogHeat
                }else{
                  method="allBAM"
                }

                if (isvalid(input$quantileThreshAnalogHeat)){
                  quantThresh=input$quantileThreshAnalogHeat
                }else{
                  quantThresh=0.9
                }                

                if(numbersamples<10){
                  factormult=10/numbersamples
                }else{
                  factormult=1
                }
                matProc_analogic=toplot$analogic$matrixes_processed

                #modify max values above quantile threshold. one for each BAM (not fixed value)
                nbin=heatvariables$binsAnalogHeat

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

                if(optioncolors=="global"){
                  #here, only one palette for everything
                  if(palette_col=="rainbow"){
                    palette <- rev(colorRampPalette(brewer.pal(9, 'Spectral'))(n =brk-1))
                  }else{
                    colorsplitted=strsplit(palette_col,split="_")[[1]]
                    palette=colorRampPalette(colorsplitted)(n=brk-1)
                  }
                  
                }else{

                  if(length(palette_col)==1 ){
                    #because it can be only ONE enrichment to show => if custom color, "red", and
                    #should be treated differently from "white_red" or example
                    if(palette_col=="white_red4"){
                      colorsplitted=strsplit(palette_col,split="_")[[1]]
                      palette=colorRampPalette(colorsplitted)(n =brk-1)  
                    }
                    else{
                      palette=c()
                      for(i in 1:length(palette_col)){
                        palette=c(palette,colorRampPalette(c("white",palette_col[i]))(n =brk-1))
                      }
                    }
                              
                  }else{
                    #if optioncolors is custom, palette should have length==nbams
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
                if(!inherits(trasp,"matrix")  ){
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
                if (isvalid(isolate(input$optioncolorsforAnalogHeat))){
                  optioncolors= isolate(input$optioncolorsforAnalogHeat)
                }else{
                  optioncolors="global"
                }
                if(length(optioncolors)>0){
                  if (optioncolors=="global"){
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
                }else{
                  NULL
                }

              })

              #when pressing the button (new analysis), all profiles/boxes will be reset
              output$profileAnalogHeat<-renderPlot({NULL})
              output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
              output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
              output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
              output$saveprofileAnalogHeat=renderUI({NULL})
              output$boxplotByROIAnalogHeat<-renderPlot({NULL})
              output$saveboxplotByROIAnalogHeat=renderUI({NULL})
              output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
              output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
              output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
              output$showboxAnalogHeat_logOptions<-renderUI({NULL})
              output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
              output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})



            }else{
              heatvariables$BAMsForAnalogHeat=NULL
              output$heatmapAnalog<-renderPlot({NULL})
              #clean also profiles and boxes from heat selection
              output$profileAnalogHeat<-renderPlot({NULL})
              output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
              output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
              output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
              output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
              output$boxplotByROIAnalogHeat<-renderPlot({NULL})
              output$corAnalogHeat<-renderPlot({NULL})
              output$pcorAnalogHeat<-renderPlot({NULL})
              output$saveprofileAnalogHeat=renderUI({NULL})
              output$saveboxplotByROIAnalogHeat=renderUI({NULL})
              output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
              output$savecorAnalogHeat=renderUI({NULL})
              output$savepcorAnalogHeat=renderUI({NULL})
              output$showoptioncolorsforAnalogHeat<-renderUI({NULL})
              output$showchooseQuantileMethodAnalogHeat<-renderUI({NULL})
              output$showquantileThreshAnalogHeat<-renderUI({NULL}) 
              output$showcolorsheat<-renderUI({NULL})
              output$textfractionelementsAnalogHeat<-renderText({NULL})
              output$colorScaleAnalogHeat<-renderPlot({NULL})
              output$textNameAnalogHeat <- renderPlot({NULL})
              output$textNameClustAnalogHeat<-renderPlot({NULL})
              toplot$analogic$finalrange=NULL
              toplot$analogic$finalbam=NULL
              toplot$analogic$finalkey=NULL
              toplot$analogic$finalnorm=NULL
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
            #clean also profiles and boxes from heat selection
            output$profileAnalogHeat<-renderPlot({NULL})
            output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
            output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
            output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
            output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
            output$boxplotByROIAnalogHeat<-renderPlot({NULL})
            output$corAnalogHeat<-renderPlot({NULL})
            output$pcorAnalogHeat<-renderPlot({NULL})
            output$saveprofileAnalogHeat=renderUI({NULL})
            output$textNameClustAnalogHeat<-renderPlot({NULL})
              output$saveboxplotByROIAnalogHeat=renderUI({NULL})
              output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
              output$savecorAnalogHeat=renderUI({NULL})
              output$savepcorAnalogHeat=renderUI({NULL})
output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
            output$showoptioncolorsforAnalogHeat<-renderUI({NULL})
            output$showchooseQuantileMethodAnalogHeat<-renderUI({NULL})
            output$showquantileThreshAnalogHeat<-renderUI({NULL}) 
            output$showcolorsheat<-renderUI({NULL})
            output$textfractionelementsAnalogHeat<-renderText({NULL})
            output$colorScaleAnalogHeat<-renderPlot({NULL})
            output$textNameAnalogHeat <- renderPlot({NULL})
            toplot$analogic$finalrange=NULL
            toplot$analogic$finalbam=NULL
            toplot$analogic$finalkey=NULL
            toplot$analogic$finalnorm=NULL
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
          #clean also profiles and boxes from heat selection
          output$profileAnalogHeat<-renderPlot({NULL})
          output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
          output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
          output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
          output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
          output$boxplotByROIAnalogHeat<-renderPlot({NULL})
          output$corAnalogHeat<-renderPlot({NULL})
          output$pcorAnalogHeat<-renderPlot({NULL})
          output$saveprofileAnalogHeat=renderUI({NULL})
          output$textNameClustAnalogHeat<-renderPlot({NULL})
              output$saveboxplotByROIAnalogHeat=renderUI({NULL})
              output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
              output$savecorAnalogHeat=renderUI({NULL})
              output$savepcorAnalogHeat=renderUI({NULL})
output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
          output$showoptioncolorsforAnalogHeat<-renderUI({NULL})
          output$showchooseQuantileMethodAnalogHeat<-renderUI({NULL})
          output$showquantileThreshAnalogHeat<-renderUI({NULL}) 
          output$showcolorsheat<-renderUI({NULL})
          output$textfractionelementsAnalogHeat<-renderText({NULL})
          output$colorScaleAnalogHeat<-renderPlot({NULL})
          output$textNameAnalogHeat <- renderPlot({NULL})
          toplot$analogic$finalrange=NULL
          toplot$analogic$finalbam=NULL
          toplot$analogic$finalkey=NULL
          toplot$analogic$finalnorm=NULL
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
        #clean also profiles and boxes from heat selection
        output$profileAnalogHeat<-renderPlot({NULL})
        output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
        output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
        output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
        output$boxplotByROIAnalogHeat<-renderPlot({NULL})
        output$corAnalogHeat<-renderPlot({NULL})
        output$pcorAnalogHeat<-renderPlot({NULL})
        output$saveprofileAnalogHeat=renderUI({NULL})
        output$textNameClustAnalogHeat<-renderPlot({NULL})
              output$saveboxplotByROIAnalogHeat=renderUI({NULL})
              output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
              output$savecorAnalogHeat=renderUI({NULL})
              output$savepcorAnalogHeat=renderUI({NULL})
output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
        output$showoptioncolorsforAnalogHeat<-renderUI({NULL})
        output$showchooseQuantileMethodAnalogHeat<-renderUI({NULL})
        output$showquantileThreshAnalogHeat<-renderUI({NULL}) 
        output$showcolorsheat<-renderUI({NULL})
        output$textfractionelementsAnalogHeat<-renderText({NULL})
        output$colorScaleAnalogHeat<-renderPlot({NULL})
        output$textNameAnalogHeat <- renderPlot({NULL})
        toplot$analogic$finalrange=NULL
        toplot$analogic$finalbam=NULL
        toplot$analogic$finalkey=NULL
        toplot$analogic$finalnorm=NULL
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
      #clean also profiles and boxes from heat selection
      output$profileAnalogHeat<-renderPlot({NULL})
      output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
      output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
      output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
      output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
      output$boxplotByROIAnalogHeat<-renderPlot({NULL})
      output$corAnalogHeat<-renderPlot({NULL})
      output$pcorAnalogHeat<-renderPlot({NULL})
      output$saveprofileAnalogHeat=renderUI({NULL})
      output$textNameClustAnalogHeat<-renderPlot({NULL})
              output$saveboxplotByROIAnalogHeat=renderUI({NULL})
              output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
              output$savecorAnalogHeat=renderUI({NULL})
              output$savepcorAnalogHeat=renderUI({NULL})
      output$showoptioncolorsforAnalogHeat<-renderUI({NULL})
output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
      output$showchooseQuantileMethodAnalogHeat<-renderUI({NULL})
      output$showquantileThreshAnalogHeat<-renderUI({NULL}) 
      output$showcolorsheat<-renderUI({NULL})
      output$textfractionelementsAnalogHeat<-renderText({NULL})
      output$colorScaleAnalogHeat<-renderPlot({NULL})
      output$textNameAnalogHeat <- renderPlot({NULL})
      toplot$analogic$finalrange=NULL
      toplot$analogic$finalbam=NULL
      toplot$analogic$finalkey=NULL
      toplot$analogic$finalnorm=NULL
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
    #clean also profiles and boxes from heat selection
    output$profileAnalogHeat<-renderPlot({NULL})
    output$showprofileAnalogHeat_logOptions<-renderUI({NULL})
    output$showprofileAnalogHeat_colorschemeOptions<-renderUI({NULL})
    output$showprofileAnalogHeat_colorlistOptions<-renderUI({NULL})
    output$boxplotByBAMAnalogHeat<-renderPlot({NULL})
    output$boxplotByROIAnalogHeat<-renderPlot({NULL})
    output$corAnalogHeat<-renderPlot({NULL})
    output$pcorAnalogHeat<-renderPlot({NULL})
    output$saveprofileAnalogHeat=renderUI({NULL})
    output$textNameClustAnalogHeat<-renderPlot({NULL})
              output$saveboxplotByROIAnalogHeat=renderUI({NULL})
              output$saveboxplotByBAMAnalogHeat=renderUI({NULL})
              output$savecorAnalogHeat=renderUI({NULL})
              output$savepcorAnalogHeat=renderUI({NULL})
output$showboxAnalogHeat_colorlistOptions<-renderUI({NULL})
        output$showboxAnalogHeat_logOptions<-renderUI({NULL})
output$showboxAnalogHeat_colorschemeOptions<-renderUI({NULL})
output$showboxAnalogHeat_groupcolOptions<-renderUI({NULL})
    output$showoptioncolorsforAnalogHeat<-renderUI({NULL})
    output$showchooseQuantileMethodAnalogHeat<-renderUI({NULL})
    output$showquantileThreshAnalogHeat<-renderUI({NULL}) 
    output$showcolorsheat<-renderUI({NULL})
    output$textfractionelementsAnalogHeat<-renderText({NULL})
    output$colorScaleAnalogHeat<-renderPlot({NULL})
    output$textNameAnalogHeat <- renderPlot({NULL})
    toplot$analogic$finalrange=NULL
    toplot$analogic$finalbam=NULL
    toplot$analogic$finalkey=NULL
    toplot$analogic$finalnorm=NULL
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
            if(any(palette_col=="rainbow")){
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
              }
              else{
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
            par(mar = c(17,12,0,0))
            image(0:nrow(toplot$analogic$clusternumbermatrix), 0:ncol(toplot$analogic$clusternumbermatrix),toplot$analogic$clusternumbermatrix,col=toplot$analogic$color_distinct_cluster,axes = FALSE, xlab = "", ylab = "",ylim=c(0,ncol(trasp)))
            axis( 2,at=pos_vertical,labels=paste(rev(heatvariables$ROIsForAnalogHeat),"\n",labels_vert),las=1,tick=FALSE)
            axis(1,0.5,labels="cluster",las= 2,tick=FALSE)
            par(mar = c(17,0,0,0))
            image(0:nrow(trasp), 0:ncol(trasp),trasp[,ncol(trasp):1,drop=FALSE],axes = FALSE, xlab = "", ylab = "",col=palette
            #REMOVE x,y lim if cordinates don't match
                  ,xlim=c(0,nrow(trasp)*factormult),ylim=c(0,ncol(trasp))   
                  )            
          }else{
            par(mar = c(17,12,0,0))
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
      if(input$colorscheme_profileAnalogHeat=="custom"){
        tosearch=paste("colorCustomAnalogHeat_profile",1:length(toplot$gadgetanalogic$portionlist_profile),sep="")              
        if(length(tosearch)>0){
          listinputcols=list()
          for(i in 1:length(tosearch)){
            listinputcols[[i]]=input[[tosearch[i]]]
          }
          listinputcols=unlist(listinputcols)
          if(length(listinputcols)==length(tosearch)){
            colors=listinputcols
          }else{
            colors=toplot$gadgetanalogic$color_distinct
          }
        }else{
          colors=toplot$gadgetanalogic$color_distinct
        }
      }else{
        #random colors
        colors=toplot$gadgetanalogic$color_distinct
      }
      plot_analog_profile(islog2=input$isLog2_profileAnalogHeat,portionlist_profile=toplot$gadgetanalogic$portionlist_profile,
                      colors=colors,ispdf=TRUE)
      dev.off()
  } 
)




#observer for PDF to download boxplot by ROI
output$saveboxplotByROIAnalogHeatbutton<- downloadHandler(
  filename=function() {
      paste('BoxplotByROI_analog_heatmap.pdf', sep='')
  },
  content=function(file) {
      

      islog=input$isLog2_boxAnalogHeat
      grouped=input$GroupColorsAnalogHeat_box

      if (input$colorscheme_boxAnalogHeat=="random"){
        if(!grouped){
          newcols=toplot$gadgetanalogic$color_distinct
          newnames=names(toplot$gadgetanalogic$portionlist_boxes)
        }else{
          newnames=unique(sapply(strsplit(names(toplot$gadgetanalogic$portionlist_boxes),split=";"),"[[",1))
          newcols=rep(toplot$gadgetanalogic$color_distinct[1:toplot$gadgetanalogic$bamnumber],toplot$gadgetanalogic$roinumber)
        }

      }else{
        #forse grouped to FALSE. Can happen sometimes that even with color=custom, group is still TRUE
        grouped=FALSE
        #here is valid color scheme but not random: extract single colors
        tosearch=paste("colorCustomAnalogHeat_box",1:length(toplot$gadgetanalogic$portionlist_boxes),sep="")              
        if(length(tosearch)>0){
          listinputcols=list()
          for(i in 1:length(tosearch)){
            listinputcols[[i]]=input[[tosearch[i]]]
          }
          listinputcols=unlist(listinputcols)
          if(length(listinputcols)==length(tosearch)){
            newcols=listinputcols
            newnames=names(toplot$gadgetanalogic$portionlist_boxes)
          }else{
            newcols=toplot$gadgetanalogic$color_distinct
            newnames=names(toplot$gadgetanalogic$portionlist_boxes)
          }
        }else{
          newcols=toplot$gadgetanalogic$color_distinct
          newnames=names(toplot$gadgetanalogic$portionlist_boxes)
        }
      }

      pdf(file)
      plot_analog_boxByROI(materialtoplot=toplot$gadgetanalogic$portionlist_boxes,roiname=toplot$gadgetanalogic$roiname,
                                bamname=toplot$gadgetanalogic$bamname,newnames=newnames,islog=islog,
                                isgrouped=grouped,colors=newcols,ispdf=TRUE)

      dev.off()
  } 
)





#observer for boxplot by BAM PDF download button
output$saveboxplotByBAMAnalogHeatbutton<- downloadHandler(
  filename=function() {
      paste('BoxplotByBAM_analog_heatmap.pdf', sep='')
  },
  content=function(file) {
      
      newlist=list()
      newcols=c()
      newnames=c()

      #if random, check if group or not
      if (input$colorscheme_boxAnalogHeat=="random"){
        #if grouping colors, adjust legend and colors by groups
        if(!input$GroupColorsAnalogHeat_box){
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
      }else{
        #extract all colors,because color scheem is custom
        tosearch=paste("colorCustomAnalogHeat_box",1:length(toplot$gadgetanalogic$portionlist_boxes),sep="")              
        if(length(tosearch)>0){
          listinputcols=list()
          for(i in 1:length(tosearch)){
            listinputcols[[i]]=input[[tosearch[i]]]
          }
          listinputcols=unlist(listinputcols)
          if(length(listinputcols)==length(tosearch)){
            newcols=listinputcols
            for(i in 1:toplot$gadgetanalogic$bamnumber){
              pos_inverted=seq(i,length(toplot$gadgetanalogic$portionlist_boxes),toplot$gadgetanalogic$bamnumber)
              newlist=c(newlist,toplot$gadgetanalogic$portionlist_boxes[pos_inverted])
              newnames=c(newnames,names(toplot$gadgetanalogic$portionlist_boxes)[pos_inverted])
            }
          }else{
            for(i in 1:toplot$gadgetanalogic$bamnumber){
              pos_inverted=seq(i,length(toplot$gadgetanalogic$portionlist_boxes),toplot$gadgetanalogic$bamnumber)
              newlist=c(newlist,toplot$gadgetanalogic$portionlist_boxes[pos_inverted])
              newcols=c(newcols,toplot$gadgetanalogic$color_distinct[pos_inverted])
              newnames=c(newnames,names(toplot$gadgetanalogic$portionlist_boxes)[pos_inverted])
            }
          }
        }else{
          for(i in 1:toplot$gadgetanalogic$bamnumber){
            pos_inverted=seq(i,length(toplot$gadgetanalogic$portionlist_boxes),toplot$gadgetanalogic$bamnumber)
            newlist=c(newlist,toplot$gadgetanalogic$portionlist_boxes[pos_inverted])
            newcols=c(newcols,toplot$gadgetanalogic$color_distinct[pos_inverted])
            newnames=c(newnames,names(toplot$gadgetanalogic$portionlist_boxes)[pos_inverted])
          }
        }
      }            



      pdf(file)

      plot_analog_boxByBAM(materialtoplot=newlist,roinumber=toplot$gadgetanalogic$roinumber,newcols=newcols,newnames=newnames,
                                bamname=toplot$gadgetanalogic$bamname,islog=input$isLog2_boxAnalogHeat,colors=newcols,ispdf=TRUE)

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
            correct_keys=toplot$analogic$finalkey[[pos]]
            correct_norm=toplot$analogic$finalnorm[[pos]]
            newBAMlist=newkeylist=as.list(rep(NA,length(correct_bams)))
            for(i in 1:length(newBAMlist)){
              newBAMlist[[i]]=correct_bams[[i]][newstart:newstop]
              newkeylist[[i]]=correct_keys[[i]][newstart:newstop]
            }
            names(newBAMlist)=names(newkeylist)=names(correct_bams)
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
              Enrichlist$rawcoverage[[input$newROIfromAnalogHeat]]=newBAMlist
              Enrichlist$decryptkey[[input$newROIfromAnalogHeat]]=newkeylist
              Enrichlist$normfactlist[[input$newROIfromAnalogHeat]]=correct_norm
              ROIvariables$listROI[[length(ROIvariables$listROI)+1]]=new("RegionOfInterest",
                                              name=input$newROIfromAnalogHeat,
                                              range=correct_range,
                                              fixed=correct_fix,
                                              BAMlist=list(),
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


#respond to help buttons for parameters
observeEvent(input$help_enrichmentInRois_parameters_ROIs, {boxHelpServer(help_enrichmentInRois_parameters_ROIs)})
observeEvent(input$help_enrichmentInRois_parameters_enrichments, {boxHelpServer(help_enrichmentInRois_parameters_enrichments)})
observeEvent(input$help_enrichmentInRois_parameters_bins, {boxHelpServer(help_enrichmentInRois_parameters_bins)})

observeEvent(input$help_enrichmentInRois_profile_totalread, {boxHelpServer(help_singleEvaluation_normalizationtotalread)})
observeEvent(input$help_enrichmentInRois_profile_readdensity, {boxHelpServer(help_singleEvaluation_normalizationreaddensity)})
observeEvent(input$help_enrichmentInRois_box_totalread, {boxHelpServer(help_singleEvaluation_normalizationtotalread)})
observeEvent(input$help_enrichmentInRois_box_readdensity, {boxHelpServer(help_singleEvaluation_normalizationreaddensity)})
observeEvent(input$help_enrichmentInRois_heat_totalread, {boxHelpServer(help_singleEvaluation_normalizationtotalread)})
observeEvent(input$help_enrichmentInRois_heat_readdensity, {boxHelpServer(help_singleEvaluation_normalizationreaddensity)})











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
    rawvals=Enrichlist$rawcoverage[pos]
    keyvals=Enrichlist$decryptkey[pos]
    normvals=Enrichlist$normfactlist[pos]
    maxbinstohave=min(unlist(lapply(roi,checkMaxBins)))

    toplot$profileAndBoxes$maxbinstohave=maxbinstohave
    
    #execute everything only if number of bins acceptable
    if(input$binsProfilesAndBox<=maxbinstohave){
      lengthROIs=list()
      nameROI=unlist(lapply(roi,getName))
      bigbamlist=list()
      bigkeyvalslist=list()
      bignormlist=list()
      for(i in 1:length(roi)){
        lengthROIs[[i]]=list(width(getRange(roi[[i]])))
        allbams=names(rawvals[[i]] )
        pos2=match(input$BAMsForProfilesAndBox,allbams)
        #order even if BAM reordered, should be fine
        bigbamlist[[i]]=rawvals[[i]][pos2]
        bigkeyvalslist[[i]]=keyvals[[i]][pos2]
        bignormlist[[i]]=normvals[[i]][pos2]
        names(bigbamlist[[i]])=names(bigkeyvalslist[[i]])=names(bignormlist[[i]])=input$BAMsForProfilesAndBox
      }
      names(lengthROIs)=names(bigbamlist)=names(bigkeyvalslist)=names(bignormlist)=nameROI




      #from bigbamlist (1st order: roi, 2nd order: BAM) do the binning from baseCOverage to matrix
      matlists=list()
    
      if(nc==1){
        matlists=lapply(1:length(bigbamlist),function(i) {
          matlist=lapply(1:length(bigbamlist[[i]]),function(k) {
                      provv=bigbamlist[[i]][[k]]
                      #provv=decryptcov( list(provv,bigkeyvalslist[[i]][[k]]),chunk=length(provv))
                      return(makeMatrixFrombaseCoverage(GRbaseCoverageOutput=provv,Nbins=binstouse,Snorm=FALSE,key=bigkeyvalslist[[i]][[k]],norm_factor=bignormlist[[i]][[k]]))
                    }) 
          names(matlist)=names(bigbamlist[[i]])
          # ###HERE invert "-" strand. They are in strands list. each element of this list
          # ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
          #pos_toinvert=as.character(strands[[i]])=="-"
          # for(k in 1:length(matlist)){
          #   #for each of the BAM file in this current ROI;
          #   matprovv=matlist[[k]]
          #   if(!all(pos_toinvert)==FALSE | dim(matprovv)[2]>1){
          #     matlist[[k]][pos_toinvert,]=matprovv[pos_toinvert,][,ncol(matprovv[pos_toinvert,]):1]
          #   }
          # }
          return(matlist)
        })
      }else{

        decision=0
        tryCatch({
          #if number of ROIs is > than number of BAMs, parallelize ROIs
          if(length(bigbamlist)>length(bigbamlist[[i]])){
            matlists=mclapply(1:length(bigbamlist),function(i) {
              matlist=lapply(1:length(bigbamlist[[i]]),function(k) {
                      provv=bigbamlist[[i]][[k]]
                      #provv=decryptcov( list(provv,bigkeyvalslist[[i]][[k]]),chunk=length(provv))
                      return(makeMatrixFrombaseCoverage(GRbaseCoverageOutput=provv,Nbins=binstouse,Snorm=FALSE,key=bigkeyvalslist[[i]][[k]],norm_factor=bignormlist[[i]][[k]]))
                        }) 
              names(matlist)=names(bigbamlist[[i]])

              # ###HERE invert "-" strand. They are in strands list. each element of this list
              # ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
              return(matlist)
            },mc.cores=nc)            
          
          }else{
            matlists=lapply(1:length(bigbamlist),function(i) {
              matlist=mclapply(1:length(bigbamlist[[i]]),function(k) {
                      provv=bigbamlist[[i]][[k]]
                      #provv=decryptcov( list(provv,bigkeyvalslist[[i]][[k]]),chunk=length(provv))
                      return(makeMatrixFrombaseCoverage(GRbaseCoverageOutput=provv,Nbins=binstouse,Snorm=FALSE,key=bigkeyvalslist[[i]][[k]],norm_factor=bignormlist[[i]][[k]]))
                        },mc.cores=nc) 
              names(matlist)=names(bigbamlist[[i]])

              # ###HERE invert "-" strand. They are in strands list. each element of this list
              # ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
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
                      provv=bigbamlist[[i]][[k]]
                      #provv=decryptcov( list(provv,bigkeyvalslist[[i]][[k]]),chunk=length(provv))
                      return(makeMatrixFrombaseCoverage(GRbaseCoverageOutput=provv,Nbins=binstouse,Snorm=FALSE,key=bigkeyvalslist[[i]][[k]],norm_factor=bignormlist[[i]][[k]]))
                      }) 
            names(matlist)=names(bigbamlist[[i]])

            # ###HERE invert "-" strand. They are in strands list. each element of this list
            # ###is a ROI. With "bamsignals" package, negative strand is already reversed at the beginning
            return(matlist)
          })
        }
      }
   
      names(matlists)=nameROI
      
      #print("reordering lists")
      #obtain list form list of lists
      portionlist=list()
      listforbox=list()
      keysforbox=list()
      normforbox=list()
      lengthROIs_unrolled=list()
      for(i in 1:length(matlists)){
        provv=matlists[[i]]
        provv_box=bigbamlist[[i]]
        provv_keybox=bigkeyvalslist[[i]]
        provv_normbox=bignormlist[[i]]
        l_ROI=rep(lengthROIs[[i]],length(provv))
        names(l_ROI)=names(provv)=names(provv_box)=names(provv_normbox)=paste("ROI:",names(matlists)[i]," (",nrow(matlists[[i]][[1]]),"); ",names(matlists[[i]]),sep="" )
        portionlist=c(portionlist,provv)
        lengthROIs_unrolled=c(lengthROIs_unrolled,l_ROI)
        listforbox=c(listforbox,provv_box)
        keysforbox=c(keysforbox,provv_keybox)
        normforbox=c(normforbox,provv_normbox)
      }

      #while values for profiles were already normalized, values for boxes are not

      #do not apply mean. Wait the plot to decide is Snorm 
      portionlist_profile=portionlist


      #apply cumulative enrichment for boxplots.
      #need to use the original bamlist, because in binning we loose the remaining
      #of the division
      portionlist_boxes=lapply(1:length(listforbox),function(i){
        block=listforbox[[i]]
        key=keysforbox[[i]]
        norm=normforbox[[i]]
        len=sapply(block,length)
        #decrypt and normalize data for boxes
        mat_boxes=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=block,Nbins=1,Snorm=FALSE,key=key,norm_factor=norm) 
        return(mat_boxes) 
      })
      names(portionlist_profile)=names(portionlist_boxes)=names(portionlist)

      n=length(portionlist_profile)
      color_distinct=colors_list[1:n]


      #kep the values for the scatterplot, when you select cells in the correlation heatmap

      names(portionlist_profile)=names(portionlist_boxes)=names(portionlist)
      toplot$profileAndBoxes$portionlist_profile=portionlist_profile
      toplot$profileAndBoxes$portionlist_boxes=portionlist_boxes
      toplot$profileAndBoxes$lengthROIs_unrolled=lengthROIs_unrolled
      toplot$profileAndBoxes$color_distinct=color_distinct

      roinumber=length(matlists)
      bamnumber=length(matlists[[1]])

      toplot$profileAndBoxes$roinumber=roinumber
      toplot$profileAndBoxes$bamnumber=bamnumber
      toplot$profileAndBoxes$bamname=names(bigbamlist[[1]])
      toplot$profileAndBoxes$roiname=names(bigbamlist)



      #custom options for profile
      output$showprofileProfileAndBox_logOptions<-renderUI({
        checkboxInput("isLog2_profileProfileAndBox", label="log2",value = FALSE, width = NULL)
      })

      output$showprofileProfileAndBox_colorschemeOptions<-renderUI({
        radioButtons("colorscheme_profileProfileAndBox",label="Choose color scheme:",choices=c(
                                                  "Random colors"="random",
                                                  "Custom colors"="custom"
                                                        ),selected="random")       
      }) 
      output$showprofileProfileAndBox_isdensityOptions<-renderUI({
        radioButtons("normalization_profileProfileAndBox",label="Normalization method:",
                                      choiceNames=list(
                                        htmlhelp("Total reads","help_enrichmentInRois_profile_totalread"),
                                        htmlhelp("Read density (reads/bp)","help_enrichmentInRois_profile_readdensity")
                                      ),
                                      choiceValues=list("totread","readdensity"),selected="readdensity")
      })     


      print("Drawing profiles and boxplots")
      ##PLOT profiles and boxplots using the obtained matrix lists
      output$profileProfilesAndBox<-renderPlot({
        orig_names=names(portionlist_profile)
        if(isvalid(input$isLog2_profileProfileAndBox)){
          islog2=input$isLog2_profileProfileAndBox
        }else{
          islog2=FALSE
        }

        if (isvalid(input$normalization_profileProfileAndBox)){
          norm_method=input$normalization_profileProfileAndBox
        }else{
          norm_method="readdensity"
        }

        if (isvalid(input$colorscheme_profileProfileAndBox)){
          if(input$colorscheme_profileProfileAndBox=="custom"){
            #only the first if custom
            if (isvalid(input$colorCustomProfileAndBox_profile1)){
              #extract colors based on choice (loop through them):
              tosearch=paste("colorCustomProfileAndBox_profile",1:length(portionlist_profile),sep="")              
              if(length(tosearch)>0){
                listinputcols=list()
                for(i in 1:length(tosearch)){
                  listinputcols[[i]]=input[[tosearch[i]]]
                }
                listinputcols=unlist(listinputcols)
                if(length(listinputcols)==length(tosearch)){
                  colors=listinputcols
                }else{
                  colors=color_distinct
                }
              }else{
                colors=color_distinct
              }
            }else{
              colors=color_distinct
            }
          }else{
            colors=color_distinct
          }
        }else{
          colors=color_distinct
        }


        yl="Total reads"
        #here Snorm:
        if (norm_method=="readdensity"){
          portionlist_profile=lapply(1:length(portionlist_profile),function(i){
            mat=portionlist_profile[[i]]
            mat2=mat/lengthROIs_unrolled[[i]]
            return(mat2)
          })
          yl="Read density (reads/bp)"
        }#else, do not do anything

        portionlist_profile=lapply(1:length(portionlist_profile),function(i){
          mat=portionlist_profile[[i]]
          mat_profile=apply(mat,2,mean)
        })
        names(portionlist_profile)=orig_names
        plot_profilesAndBox_profile(islog2=islog2,portionlist_profile=portionlist_profile,
                                colors=colors,yl=yl,ispdf=FALSE) 
      })

      #PDF download button for profiles
      output$saveprofileProfilesAndBox=renderUI({downloadButton('saveprofileProfilesAndBoxbutton', 'Get PDF')})






      #custom options for boxes
      output$showBoxProfileAndBox_logOptions<-renderUI({
        checkboxInput("isLog2_boxProfileAndBox", label="log2",value = FALSE, width = NULL)
      })

      output$showBoxProfileAndBox_colorschemeOptions<-renderUI({
        radioButtons("colorscheme_boxProfileAndBox",label="Choose color scheme:",choices=c(
                                                  "Random colors"="random",
                                                  "Custom colors"="custom"
                                                        ),selected="random")       
      }) 
      output$showBoxProfileAndBox_isdensityOptions<-renderUI({
        radioButtons("normalization_boxProfileAndBox",label="Normalization method:",
                                                choiceNames=list(
                                                  htmlhelp("Total reads","help_enrichmentInRois_box_totalread"),
                                                  htmlhelp("Read density (reads/bp)","help_enrichmentInRois_box_readdensity")
                                                ),choiceValues=list("totread","readdensity"),selected="readdensity")
      })  





      #plot boxplots of regions selected. Use portionlist_boxes
      output$boxByBAMProfilesAndBox<-renderPlot({

        if (isvalid(input$isLog2_boxProfileAndBox)){
          islog=input$isLog2_boxProfileAndBox
        }else{
          islog=FALSE
        }
        if (isvalid(input$normalization_boxProfileAndBox)){
          normalization=input$normalization_boxProfileAndBox
        }else{
          normalization="readdensity"
        }
        if (isvalid(input$GroupColorsProfileAndBox_box)){
          grouped=input$GroupColorsProfileAndBox_box
        }else{
          grouped=FALSE
        }

        yl="Total reads"
        if(normalization=="readdensity"){
          for (i in 1:length(portionlist_boxes)){
            portionlist_boxes[[i]]=portionlist_boxes[[i]]/lengthROIs_unrolled[[i]]
          }
          yl="Read density (reads/bp)"
        } 

        #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
        #invert the order
        newlist=list()
        newcols=c()
        newnames=c()


        if (isvalid(input$colorscheme_boxProfileAndBox)){
          if (input$colorscheme_boxProfileAndBox=="random"){
            if (!grouped){
              #if random, not grouped colors, default simple color block
              for(i in 1:bamnumber){
                pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                newlist=c(newlist,portionlist_boxes[pos_inverted])
                newcols=c(newcols,color_distinct[pos_inverted])
                newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
              } 
            }else{
              for(i in 1:bamnumber){
                pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                newlist=c(newlist,portionlist_boxes[pos_inverted])
              }
              newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
              newcols=rep(color_distinct[1:toplot$profileAndBoxes$roinumber],toplot$profileAndBoxes$bamnumber)
            }

          }else{
            #color not random: extract all colors from input variables
            tosearch=paste("colorCustomProfileAndBox_box",1:length(portionlist_boxes),sep="")
            if (length(tosearch)>0){
              listinputcols=list()
              for(i in 1:length(tosearch)){
                listinputcols[[i]]=input[[tosearch[i]]]
              }
              listinputcols=unlist(listinputcols)
              if(length(listinputcols)==length(tosearch)){
                newcols=listinputcols
                for(i in 1:bamnumber){
                  pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                  newlist=c(newlist,portionlist_boxes[pos_inverted])
                  newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                }
              }else{
                for(i in 1:bamnumber){
                  pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                  newlist=c(newlist,portionlist_boxes[pos_inverted])
                  newcols=c(newcols,color_distinct[pos_inverted])
                  newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
                } 
              }
            }else{
              for(i in 1:bamnumber){
                pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
                newlist=c(newlist,portionlist_boxes[pos_inverted])
                newcols=c(newcols,color_distinct[pos_inverted])
                newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
              } 
            }
          }

        }else{
          #option not existing => default = random colors!
          for(i in 1:bamnumber){
            pos_inverted=seq(i,length(portionlist_boxes),bamnumber)
            newlist=c(newlist,portionlist_boxes[pos_inverted])
            newcols=c(newcols,color_distinct[pos_inverted])
            newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
          } 
        }

        plot_profilesAndBox_boxByBAM(materialtoplot=newlist,colors=newcols,yl=yl,newnames=newnames,
                                    roinumber=roinumber,bamname=names(bigbamlist[[1]]),islog=islog)

      })
 


      output$boxByROIProfilesAndBox<-renderPlot({
        if (isvalid(input$isLog2_boxProfileAndBox)){
          islog=input$isLog2_boxProfileAndBox
        }else{
          islog=FALSE
        }
        if (isvalid(input$normalization_boxProfileAndBox)){
          normalization=input$normalization_boxProfileAndBox
        }else{
          normalization="readdensity"
        }
        if (isvalid(input$GroupColorsProfileAndBox_box)){
          grouped=input$GroupColorsProfileAndBox_box
        }else{
          grouped=FALSE
        }
        yl="Total reads"
        if(normalization=="readdensity"){
          for (i in 1:length(portionlist_boxes)){
            portionlist_boxes[[i]]=portionlist_boxes[[i]]/lengthROIs_unrolled[[i]]
          }
          yl="Read density (reads/bp)"
        } 


        if (isvalid(input$colorscheme_boxProfileAndBox)){
          if (input$colorscheme_boxProfileAndBox=="random"){
            if (!grouped){
              newcols=color_distinct
              newnames=names(portionlist)
            }else{
              newnames=unique(sapply(strsplit(names(portionlist),split=";"),"[[",1))
              newnames=sapply(strsplit(newnames,split="ROI:"),"[[",2)
              newcols=rep(color_distinct[1:toplot$profileAndBoxes$bamnumber],toplot$profileAndBoxes$roinumber)
            }
          }else{
            grouped=FALSE
            #here is valid color scheme but not random: extract single colors
            tosearch=paste("colorCustomProfileAndBox_box",1:length(portionlist_boxes),sep="") 
            if(length(tosearch)>0){
              listinputcols=list()
              for(i in 1:length(tosearch)){
                listinputcols[[i]]=input[[tosearch[i]]]
              }
              listinputcols=unlist(listinputcols)
              if (length(listinputcols)==length(tosearch)){
                newcols=listinputcols
                newnames=names(portionlist)
              }else{
                newcols=color_distinct
                newnames=names(portionlist)
              }
            }else{
              newcols=color_distinct
              newnames=names(portionlist)
            }

          }
        }else{
          grouped=FALSE
          newcols=color_distinct
          newnames=names(portionlist)
        }

        plot_profilesAndBox_boxByROI(materialtoplot=portionlist_boxes,roiname=names(bigbamlist),bamname=names(bigbamlist[[1]]),
                              newnames=newnames,islog=islog,yl=yl,isgrouped=grouped,
                              colors=newcols,ispdf=FALSE)
      })  




      #PDF download buttons for boxplots by ROI and by BAM
      output$saveboxByROIProfilesAndBox=renderUI({downloadButton('saveboxByROIProfilesAndBoxbutton', 'Get PDF')})
      output$saveboxByBAMProfilesAndBox=renderUI({downloadButton('saveboxByBAMProfilesAndBoxbutton', 'Get PDF')})
      #download data of boxplot by ROI
      output$saveboxdataProfANDbox=renderUI({downloadButton('saveenrichmentBoxProfANDboxdata', 'Download data')})



      if(bamnumber>1 & roinumber==1){



        #specific options for cor heatmap (readdensity, log2, cormethod)
        #log2 present only if cormethod==pearson

        output$showCorProfileAndBox_isdensityOptions<-renderUI({
          radioButtons("normalization_CorProfileAndBox",label="Normalization method:",
                                                choiceNames=list(
                                                  htmlhelp("Total reads","help_enrichmentInRois_heat_totalread"),
                                                  htmlhelp("Read density (reads/bp)","help_enrichmentInRois_heat_readdensity")
                                                ),choiceValues=list("totread","readdensity"),selected="readdensity")
        })  
        output$showCorProfileAndBox_corMethodOptions<-renderUI({
          radioButtons("corMethodProfilesAndBox_Cor",label="Type of correlation",choices=c(
                                                  "Pearson"="pearson",
                                                  "Spearman"="spearman"
                                                        ),selected="pearson")
        })  

        output$corProfilesAndBox<-renderPlot({
          if (isvalid(input$normalization_CorProfileAndBox)){
            normalization=input$normalization_CorProfileAndBox
          }else{
            normalization="readdensity"
          }

          if (isvalid(input$isLog2_CorProfileAndBox)){
            islog=input$isLog2_CorProfileAndBox
          }else{
            islog=FALSE
          }

          if (isvalid(input$corMethodProfilesAndBox_Cor)){
            #force islog=FALSE if cor method=spearman
            cormethod=input$corMethodProfilesAndBox_Cor
            if (!cormethod=="pearson"){
              islog=FALSE
            }
          }else{
            cormethod="pearson"
          }
          

          #respond to Snorm and change accordingly
          if(normalization=="readdensity"){
            for (i in 1:length(portionlist_boxes)){
              portionlist_boxes[[i]]=portionlist_boxes[[i]]/lengthROIs_unrolled[[i]]
            }
          } 
          if(islog){
            portionlist_boxes=lapply(portionlist_boxes,log2)
          }
          plot_profilesAndBox_corheat(portionlist_boxes=portionlist_boxes,bamname=names(bigbamlist[[1]]),
                        cormethod=cormethod)
        })
        


        #PDF download button for cor matrix
        output$savecorProfilesAndBox=renderUI({downloadButton('savecorProfilesAndBoxbutton', 'Get PDF')})
        output$scatterProfilesAndBox<-renderPlot({plot_text(text="click a cell in the correlation or\npartial correlation heatmap\nto view the scatterlpot",cex=1.4)})



        #at least 3 for partial correlation, otherwise problems with lm
        if(bamnumber>2){


          output$pcorProfilesAndBox<-renderPlot({
            if (isvalid(input$normalization_CorProfileAndBox)){
              normalization=input$normalization_CorProfileAndBox
            }else{
              normalization="readdensity"
            }

            if (isvalid(input$isLog2_CorProfileAndBox)){
              islog=input$isLog2_CorProfileAndBox
            }else{
              islog=FALSE
            }

            if (isvalid(input$corMethodProfilesAndBox_Cor)){
              #force islog=FALSE if cor method=spearman
              cormethod=input$corMethodProfilesAndBox_Cor
              if (!cormethod=="pearson"){
                islog=FALSE
              }
            }else{
              cormethod="pearson"
            }
            #respond to Snorm and change accordingly
            if(normalization=="readdensity"){
              for (i in 1:length(portionlist_boxes)){
                portionlist_boxes[[i]]=portionlist_boxes[[i]]/lengthROIs_unrolled[[i]]
              }
            } 

	          if(islog){
	            portionlist_boxes=lapply(portionlist_boxes,log2)
	          }

            plot_profilesAndBox_pcorheat(portionlist_boxes=portionlist_boxes,bamname=names(bigbamlist[[1]]),
                                cormethod=cormethod)
          })
  
          #PDf download button for partial correlation
          output$savepcorProfilesAndBox=renderUI({downloadButton('savepcorProfilesAndBoxbutton', 'Get PDF')})


        }else{
        	output$pcorProfilesAndBox <-renderPlot({plot_text(text="you need to select a single ROI\nwith at least\n3 enrichment files associated",cex=1.4)})
          output$savepcorProfilesAndBox=renderUI({NULL})
          output$scatterProfilesAndBox<-renderPlot({plot_text(text="click a cell in the correlation or\npartial correlation heatmap\nto view the scatterlpot",cex=1.4)})
	        toplot$profileAndBoxes$px=NULL
	        toplot$profileAndBoxes$py=NULL
        }


      }else{
        output$corProfilesAndBox <-renderPlot({plot_text(text="you need to select a single ROI\nwith at least\n2 enrichment files associated",cex=1.4)})
        output$savecorProfilesAndBox=renderUI({NULL})
        output$scatterProfilesAndBox<-renderPlot({NULL}) 
        output$pcorProfilesAndBox <-renderPlot({plot_text(text="you need to select a single ROI\nwith at least\n3 enrichment files associated",cex=1.4)})
        output$savepcorProfilesAndBox=renderUI({NULL})
        output$showCorProfileAndBox_logOptions<-renderUI({NULL})
        output$showCorProfileAndBox_isdensityOptions<-renderUI({NULL})
        output$showCorProfileAndBox_corMethodOptions<-renderUI({NULL})
        output$savescatterProfilesAndBox=renderUI({NULL})
        toplot$profileAndBoxes$x=NULL
        toplot$profileAndBoxes$y=NULL
        toplot$profileAndBoxes$px=NULL
        toplot$profileAndBoxes$py=NULL
      }
    }else{
      output$profileProfilesAndBox<-renderPlot({NULL})
      output$showprofileProfileAndBox_logOptions<-renderUI({NULL})
      output$showprofileProfileAndBox_colorschemeOptions<-renderUI({NULL})
      output$showprofileProfileAndBox_isdensityOptions<-renderUI({NULL})
      output$showprofileProfileAndBox_colorlistOptions<-renderUI({NULL})
      output$showBoxProfileAndBox_logOptions<-renderUI({NULL})
      output$showBoxProfileAndBox_colorschemeOptions<-renderUI({NULL})
      output$showBoxProfileAndBox_isdensityOptions<-renderUI({NULL})
      output$showBoxProfileAndBox_colorlistOptions<-renderUI({NULL})  
      output$showBoxProfileAndBox_groupcolOptions<-renderUI({NULL})
      output$showCorProfileAndBox_logOptions<-renderUI({NULL})
      output$showCorProfileAndBox_isdensityOptions<-renderUI({NULL})
      output$showCorProfileAndBox_corMethodOptions<-renderUI({NULL})
      output$savescatterProfilesAndBox=renderUI({NULL})
      output$boxByROIProfilesAndBox<-renderPlot({NULL})
      output$saveboxdataProfANDbox=renderUI({NULL})
      output$boxByBAMProfilesAndBox<-renderPlot({NULL})
      output$scatterProfilesAndBox<-renderPlot({NULL})   
      output$corProfilesAndBox <-renderPlot({NULL})
      output$savecorProfilesAndBox=renderUI({NULL})    
      output$pcorProfilesAndBox <-renderPlot({NULL})
      output$savepcorProfilesAndBox=renderUI({NULL})
      output$saveprofileProfilesAndBox=renderUI({NULL})
      output$saveboxByROIProfilesAndBox=renderUI({NULL})
      output$saveboxByBAMProfilesAndBox=renderUI({NULL}) 
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
    output$scatterProfilesAndBox<-renderPlot({NULL})
    output$showprofileProfileAndBox_logOptions<-renderUI({NULL})
    output$showprofileProfileAndBox_colorschemeOptions<-renderUI({NULL})
    output$showprofileProfileAndBox_isdensityOptions<-renderUI({NULL})
    output$showprofileProfileAndBox_colorlistOptions<-renderUI({NULL})
    output$showBoxProfileAndBox_logOptions<-renderUI({NULL})
    output$showBoxProfileAndBox_colorschemeOptions<-renderUI({NULL})
    output$showBoxProfileAndBox_isdensityOptions<-renderUI({NULL})
    output$showBoxProfileAndBox_colorlistOptions<-renderUI({NULL})  
    output$showBoxProfileAndBox_groupcolOptions<-renderUI({NULL})
    output$showCorProfileAndBox_logOptions<-renderUI({NULL})
    output$showCorProfileAndBox_isdensityOptions<-renderUI({NULL})
    output$showCorProfileAndBox_corMethodOptions<-renderUI({NULL})
    output$savescatterProfilesAndBox=renderUI({NULL})
    output$corProfilesAndBox <-renderPlot({NULL})
    output$savecorProfilesAndBox=renderUI({NULL})    
    output$pcorProfilesAndBox <-renderPlot({NULL})
    output$savepcorProfilesAndBox=renderUI({NULL}) 
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
        portionlist_profile=toplot$profileAndBoxes$portionlist_profile
        
        #colors
        if(input$colorscheme_profileProfileAndBox=="custom"){
          #extract colors based on choice (loop through them):
          tosearch=paste("colorCustomProfileAndBox_profile",1:length(portionlist_profile),sep="")              
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if(length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=toplot$profileAndBoxes$color_distinct
            }
          }else{
            colors=toplot$profileAndBoxes$color_distinct
          }
        }else{
          colors=toplot$profileAndBoxes$color_distinct
        }
  


        #here Snorm:
        yl="Total reads"
        if (input$normalization_profileProfileAndBox=="readdensity"){
          portionlist_profile=lapply(1:length(portionlist_profile),function(i){
            mat=portionlist_profile[[i]]
            mat2=mat/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
            return(mat2)
          })
          yl="Read density (reads/bp)"
        }#else, do not do anything

        portionlist_profile=lapply(1:length(portionlist_profile),function(i){
          mat=portionlist_profile[[i]]
          mat_profile=apply(mat,2,mean)
        })
        names(portionlist_profile)=names(toplot$profileAndBoxes$portionlist_profile)
        


        pdf(file)
        plot_profilesAndBox_profile(islog2=input$isLog2_profileProfileAndBox,portionlist_profile=portionlist_profile,
                                colors=colors,yl=yl,ispdf=TRUE)
        dev.off()
  } 
)





#PDF download button boxplot by ROI
output$saveboxByROIProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('BoxplotByROI.pdf', sep='')
  },
  content=function(file) {
        portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
        yl="Total reads"
        if(input$normalization_boxProfileAndBox=="readdensity"){
          for (i in 1:length(portionlist_boxes)){
            portionlist_boxes[[i]]=portionlist_boxes[[i]]/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
          }
          yl="Read density (reads/bp)"
        } 
        grouped=input$GroupColorsProfileAndBox_box

        if (input$colorscheme_boxProfileAndBox=="random"){
          if (!grouped){
            newcols=toplot$profileAndBoxes$color_distinct
            newnames=names(portionlist_boxes)
          }else{
            newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
            newnames=sapply(strsplit(newnames,split="ROI:"),"[[",2)
            newcols=rep(toplot$profileAndBoxes$color_distinct[1:toplot$profileAndBoxes$bamnumber],toplot$profileAndBoxes$roinumber)
          }
        }else{
          #here is valid color scheme but not random: extract single colors
          tosearch=paste("colorCustomProfileAndBox_box",1:length(portionlist_boxes),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)

            if (length(listinputcols)==length(tosearch)){
              newcols=listinputcols
              newnames=names(portionlist_boxes)
            }else{
              newcols=toplot$profileAndBoxes$color_distinct
              newnames=names(portionlist_boxes)
            }
          }else{
            newcols=toplot$profileAndBoxes$color_distinct
            newnames=names(portionlist_boxes)
          }

        }


        pdf(file)
        plot_profilesAndBox_boxByROI(materialtoplot=portionlist_boxes,roiname=toplot$profileAndBoxes$roiname,
                              bamname=toplot$profileAndBoxes$bamname,newnames=newnames,islog=input$isLog2_boxProfileAndBox,
                              yl=yl,isgrouped=grouped,colors=newcols,ispdf=TRUE)
        dev.off()
  } 
)



#observer for download of xls table of the boxplot data (by ROI)
output$saveenrichmentBoxProfANDboxdata<- downloadHandler(
  filename=function() {
      paste('boxplot_data.xls', sep='')
  },
  content=function(file) {

      portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
      if(input$normalization_boxProfileAndBox=="readdensity"){
        for (i in 1:length(portionlist_boxes)){
          portionlist_boxes[[i]]=portionlist_boxes[[i]]/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
        }
      } 

      factor_add=rep(0:(toplot$profileAndBoxes$roinumber-1),each=toplot$profileAndBoxes$bamnumber)
      addingfactor=1:length(portionlist_boxes)+factor_add
      if(input$isLog2_boxProfileAndBox){
        portionlist_boxes=lapply(1:length(portionlist_boxes),function(i) {log2(portionlist_boxes[[i]])})
      }else{
        portionlist_boxes=portionlist_boxes
      }
      names(portionlist_boxes)=names(portionlist_boxes)
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

        portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
        yl="Total reads"
        if(input$normalization_boxProfileAndBox=="readdensity"){
          for (i in 1:length(portionlist_boxes)){
            portionlist_boxes[[i]]=portionlist_boxes[[i]]/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
          }
          yl="Read density (reads/bp)"
        } 

        #will be BAM1 roi1 roi2 roi3... BAM2 roi1 roi2 roi3
        #invert the order
        newlist=list()
        newcols=c()
        newnames=c()



        if (input$colorscheme_boxProfileAndBox=="random"){
          if (!input$GroupColorsProfileAndBox_box){
            #if random, not grouped colors, default simple color block
            for(i in 1:toplot$profileAndBoxes$bamnumber){
              pos_inverted=seq(i,length(portionlist_boxes),toplot$profileAndBoxes$bamnumber)
              newlist=c(newlist,portionlist_boxes[pos_inverted])
              newcols=c(newcols,toplot$profileAndBoxes$color_distinct[pos_inverted])
              newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
            } 
          }else{
            for(i in 1:toplot$profileAndBoxes$bamnumber){
              pos_inverted=seq(i,length(portionlist_boxes),toplot$profileAndBoxes$bamnumber)
              newlist=c(newlist,portionlist_boxes[pos_inverted])
            }
            newnames=unique(sapply(strsplit(names(portionlist_boxes),split=";"),"[[",1))
            newcols=rep(toplot$profileAndBoxes$color_distinct[1:toplot$profileAndBoxes$roinumber],toplot$profileAndBoxes$bamnumber)
          }

        }else{
          #color not random: extract all colors from input variables
          tosearch=paste("colorCustomProfileAndBox_box",1:length(portionlist_boxes),sep="")
          if (length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if(length(listinputcols)==length(tosearch)){
              newcols=listinputcols
              for(i in 1:toplot$profileAndBoxes$bamnumber){
                pos_inverted=seq(i,length(portionlist_boxes),toplot$profileAndBoxes$bamnumber)
                newlist=c(newlist,portionlist_boxes[pos_inverted])
                newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
              }
            }else{
              for(i in 1:toplot$profileAndBoxes$bamnumber){
                pos_inverted=seq(i,length(portionlist_boxes),toplot$profileAndBoxes$bamnumber)
                newlist=c(newlist,portionlist_boxes[pos_inverted])
                newcols=c(newcols,toplot$profileAndBoxes$color_distinct[pos_inverted])
                newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
              } 
            }
          }else{
            for(i in 1:toplot$profileAndBoxes$bamnumber){
              pos_inverted=seq(i,length(portionlist_boxes),toplot$profileAndBoxes$bamnumber)
              newlist=c(newlist,portionlist_boxes[pos_inverted])
              newcols=c(newcols,toplot$profileAndBoxes$color_distinct[pos_inverted])
              newnames=c(newnames,names(portionlist_boxes)[pos_inverted])
            } 
          }
        }


        pdf(file)
        plot_profilesAndBox_boxByBAM(materialtoplot=newlist,colors=newcols,yl=yl,newnames=newnames,
                                    roinumber=toplot$profileAndBoxes$roinumber,bamname=toplot$profileAndBoxes$bamname,
                                    islog=input$isLog2_boxProfileAndBox,ispdf=TRUE)
        dev.off()
  } 
)





#PDF download button for correlation heatmap
output$savecorProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('Correlation_heatmap.pdf', sep='')
  },
  content=function(file) {
          portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
          #respond to Snorm and change accordingly
          if(input$normalization_CorProfileAndBox=="readdensity"){
            for (i in 1:length(portionlist_boxes)){
              portionlist_boxes[[i]]=portionlist_boxes[[i]]/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
            }
          } 

          if(input$isLog2_CorProfileAndBox){
            portionlist_boxes=lapply(portionlist_boxes,log2)
          }else{
            portionlist_boxes=portionlist_boxes
          }

          pdf(file)
          plot_profilesAndBox_corheat(portionlist_boxes=portionlist_boxes,bamname=toplot$profileAndBoxes$bamname,
                                cormethod=input$corMethodProfilesAndBox_Cor)
          dev.off()
  } 
)






#PDF download of partial correlation heatmap
output$savepcorProfilesAndBoxbutton<- downloadHandler(
  filename=function() {
      paste('PartialCorrelation_heatmap.pdf', sep='')
  },
  content=function(file) {
          portionlist_boxes=toplot$profileAndBoxes$portionlist_boxes
          #respond to Snorm and change accordingly
          if(input$normalization_CorProfileAndBox=="readdensity"){
            for (i in 1:length(portionlist_boxes)){
              portionlist_boxes[[i]]=portionlist_boxes[[i]]/toplot$profileAndBoxes$lengthROIs_unrolled[[i]]
            }
          } 

          if(input$isLog2_CorProfileAndBox){
            portionlist_boxes=lapply(portionlist_boxes,log2)
          }

          pdf(file)
          plot_profilesAndBox_pcorheat(portionlist_boxes=portionlist_boxes,bamname=toplot$profileAndBoxes$bamname,
                                cormethod=input$corMethodProfilesAndBox_Cor)
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


#help buttons for parameters
observeEvent(input$help_dynamicsOnGenes_parameters_genelist, {boxHelpServer(help_dynamicsOnGenes_parameters_genelist)})
observeEvent(input$help_dynamicsOnGenes_parameters_ROIassociated, {boxHelpServer(help_dynamicsOnGenes_parameters_ROIassociated)})
observeEvent(input$help_dynamicsOnGenes_parameters_ernichments, {boxHelpServer(help_dynamicsOnGenes_parameters_ernichments)})
observeEvent(input$help_dynamicsOnGenes_parameters_nbins, {boxHelpServer(help_dynamicsOnGenes_parameters_nbins)})
observeEvent(input$help_dynamicsOnGenes_parameters_fractionexclude, {boxHelpServer(help_dynamicsOnGenes_parameters_fractionexclude)})

observeEvent(input$help_dynamicsOnGenes_parameters_totread, {boxHelpServer(help_singleEvaluation_normalizationtotalread)})
observeEvent(input$help_dynamicsOnGenes_parameters_readdensity, {boxHelpServer(help_singleEvaluation_normalizationreaddensity)})












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
    nomi=unlist(lapply(ROIvariables$listROI,getName))

    totallist_boxplot=rep(list(NA),length(input$BAMforDynamics))
    names(totallist_boxplot)=input$BAMforDynamics
    totallist_boxplot=rep(list(totallist_boxplot),length(input$genelistsforDynamics))
    names(totallist_boxplot)=input$genelistsforDynamics 
    totallist_boxplot=rep(list(totallist_boxplot),3)
    names(totallist_boxplot)=c("TSS","GB","TES")
    #create list for lengths (normalization)
    totallist_lengths_forNorm_profile=as.list(rep(NA,length(input$genelistsforDynamics)))
    totallist_lengths_forNorm_SI_GB=as.list(rep(NA,length(input$genelistsforDynamics)))
    totallist_lengths_forNorm_box_SI_prom=as.list(rep(NA,length(input$genelistsforDynamics)))
    totallist_lengths_forNorm_box_GB=as.list(rep(NA,length(input$genelistsforDynamics)))
    totallist_lengths_forNorm_box_TES=as.list(rep(NA,length(input$genelistsforDynamics)))
    totallist_SI=totallist_boxplot
    names(totallist_SI)=c("TSS","GB","SI")

    totallist_profile=rep(list(NA),length(input$BAMforDynamics))
    names(totallist_profile)=input$BAMforDynamics
    totallist_profile=rep(list(totallist_profile),length(input$genelistsforDynamics))
    names(totallist_profile)=names(totallist_lengths_forNorm_profile)=names(totallist_lengths_forNorm_SI_GB)=
    names(totallist_lengths_forNorm_box_SI_prom)=
    names(totallist_lengths_forNorm_box_GB)=
    names(totallist_lengths_forNorm_box_TES)= input$genelistsforDynamics        
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
      #get enrichment for each roi of that genelist 
      rawvals_promoters=Enrichlist$rawcoverage[[pos_promoters]]
      keyvals_promoters=Enrichlist$decryptkey[[pos_promoters]]
      normvals_promoters=Enrichlist$normfactlist[[pos_promoters]]
      rawvals_transcripts=Enrichlist$rawcoverage[[pos_transcripts]]
      keyvals_transcripts=Enrichlist$decryptkey[[pos_transcripts]]
      normvals_transcripts=Enrichlist$normfactlist[[pos_transcripts]]
      rawvals_TES=Enrichlist$rawcoverage[[pos_TES]]
      keyvals_TES=Enrichlist$decryptkey[[pos_TES]]
      normvals_TES=Enrichlist$normfactlist[[pos_TES]]
      BAM_promoters=names(rawvals_promoters)
      BAM_transcripts=names(rawvals_transcripts)      
      BAM_TES=names(rawvals_TES) 

      #calculate all ranges and 30% increment for transcripts of that genelit
      range_transcripts=getRange(roi_transcripts)
      smalllength=width(range_transcripts)
      #biglength=unlist(lapply(content_transcripts,length))
      thirtypercent=round((smalllength/10)*3)
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
      thirtypercent=thirtypercent[tokeep]
      range_transcripts=range_transcripts[tokeep]

      pos_plus=as.character(strand(range_transcripts))!="-"
      pos_minus=as.character(strand(range_transcripts))=="-"
      thirtypercent_plus=thirtypercent[pos_plus]
      thirtypercent_minus=thirtypercent[pos_minus]
      smalllength_plus=smalllength[pos_plus]
      smalllength_minus=smalllength[pos_minus]
      smalllength=c(smalllength_plus,smalllength_minus)
      thirtypercent=c(thirtypercent_plus,thirtypercent_minus)

      widthTES_reorder=width(getRange(roi_TES)[tokeep])
      widthTES_reorder_plus=widthTES_reorder[pos_plus]
      widthTES_reorder_minus=widthTES_reorder[pos_minus]
      widthTES_reorder=c(widthTES_reorder_plus,widthTES_reorder_minus)

      #remove the 30% up/downstream and part of TSS and TES
      #for new function, remove +1, remove -1 in cpp function, multiply by 2 if key=TRUE
      start_pos_plus=thirtypercent_plus+1+part_TSS
      end_pos_plus=smalllength_plus+thirtypercent_plus-part_TES
      start_pos_minus=thirtypercent_minus+1+part_TES
      end_pos_minus=smalllength_minus+thirtypercent_minus-part_TSS

      #the length of transcripts+ 30% for interval normalization of profiles (both before TSS and after TES)
      #all bp widths normalizers
      length_transcripts_filtered_plusFraction=smalllength+2*thirtypercent
      totallist_lengths_forNorm_profile[[i]]=length_transcripts_filtered_plusFraction
      totallist_lengths_forNorm_SI_GB[[i]]=( (smalllength-(part_TSS+part_TES))+ widthTES_reorder )
      totallist_lengths_forNorm_box_SI_prom[[i]]=width(range_TSS[tokeep])
      totallist_lengths_forNorm_box_GB[[i]]=(smalllength-(part_TSS+part_TES))
      totallist_lengths_forNorm_box_TES[[i]]=width(range_TES[tokeep])
      for(k in 1:length(input$BAMforDynamics)){
        #extract baseCoverage for promoters/transcripts/TES for that gene list   
        BAMnametosearch=input$BAMforDynamics[k]
        pos_BAM_promoters=which(BAMnametosearch==BAM_promoters)
        pos_BAM_transcripts=which(BAMnametosearch==BAM_transcripts)
        pos_BAM_TES=which(BAMnametosearch==BAM_TES)
        content_promoters=rawvals_promoters[[pos_BAM_promoters]]
        contentkeys_promoters=keyvals_promoters[[pos_BAM_promoters]]
        contentnorm_promoters=normvals_promoters[[pos_BAM_promoters]]
        content_transcripts=rawvals_transcripts[[pos_BAM_transcripts]]
        contentkeys_transcripts=keyvals_transcripts[[pos_BAM_transcripts]]
        contentnorm_transcripts=normvals_transcripts[[pos_BAM_transcripts]]
        content_TES=rawvals_TES[[pos_BAM_TES]]
        contentkeys_TES=keyvals_TES[[pos_BAM_TES]]
        contentnorm_TES=normvals_TES[[pos_BAM_TES]]
        #for this genelist (p,t,t) and BAM, for promoter, transcripts and TES, calculate
        #profiles and box on the promoters,gb,TES
        #calculate content_transcripts_mod for gb only, where you cut 30% AND TSS down and TES up
        #to remove 30%, thirtypercentToRemove= l*15/80, where l is the length already increased by 30%

        #here I have to decrypt, because length of raw data already depends on decryption key
        content_transcripts=content_transcripts[tokeep]
        contentkeys_transcripts=contentkeys_transcripts[tokeep]
        content_transcripts_plus=content_transcripts[pos_plus]
        content_transcripts_minus=content_transcripts[pos_minus]
        contentkeys_transcripts_plus=contentkeys_transcripts[pos_plus]
        contentkeys_transcripts_minus=contentkeys_transcripts[pos_minus]

        #PROFILES
        #calc. matrixes for profiles, and invert negative strand
        #Snormalization will be done in the plot for interactivity
        matcov_plus=makeMatrixFrombaseCoverage(content_transcripts_plus,Nbins=toplot$dynamics$nbin,Snorm=FALSE,key=contentkeys_transcripts_plus,norm_factor=contentnorm_transcripts)
        matcov_minus=makeMatrixFrombaseCoverage(content_transcripts_minus,Nbins=toplot$dynamics$nbin,Snorm=FALSE,key=contentkeys_transcripts_minus,norm_factor=contentnorm_transcripts)
        #With "bamsignals" package, negative strand is already reversed at the beginning
        #matcov_minus <- matcov_minus[ , ncol(matcov_minus):1]
        matcov=rbind(matcov_plus,matcov_minus)
        #for interval norm, keep list [i] for each genelist for lengths (length_transcripts_filtered_plusFraction)
        totallist_profile[[i]][[k]]=matcov



        #here we sum the part of the transcripts that will be used:
        #    TSS      transcript    TES
        #.....|.....-------------....|....
        #-30%   TSS part      TES part  +30%
        #we remove the 30% up and down from original BAM files from transcripts, plus the part
        #of the TSS and TES that will not be considered as "transcript". We keep only the internal part...
        #to improve efficiency, use the custom Rcpp function
        content_transcripts_mod_plus=cutAndSumTranscripts(GRbaseCoverageOutput=content_transcripts_plus,
                                                      StartingPositions=start_pos_plus,
                                                      EndingPositions=end_pos_plus,
                                                      norm_factor=contentnorm_transcripts,
                                                      key=contentkeys_transcripts_plus
                                                      )
        content_transcripts_mod_minus=cutAndSumTranscripts(GRbaseCoverageOutput=content_transcripts_minus,
                                                      StartingPositions=start_pos_minus,
                                                      EndingPositions=end_pos_minus,
                                                      norm_factor=contentnorm_transcripts,
                                                      key=contentkeys_transcripts_minus)
        #join mod plus and minus
        content_transcripts_mod=c(content_transcripts_mod_plus,content_transcripts_mod_minus)

        #BOXPLOTS
        #for boxplots, sum everything inside TSS,TES and the part of the transcripts cut before
        content_promoters=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=content_promoters,Nbins=1,Snorm=FALSE,key=contentkeys_promoters,norm_factor=contentnorm_promoters)
        content_TES=makeMatrixFrombaseCoverage(GRbaseCoverageOutput=content_TES,Nbins=1,Snorm=FALSE,key=contentkeys_TES,norm_factor=contentnorm_TES)

        totallist_boxplot[[1]][[i]][[k]]=content_promoters[tokeep]
        totallist_boxplot[[2]][[i]][[k]]=content_transcripts_mod
        totallist_boxplot[[3]][[i]][[k]]=content_TES[tokeep]


        #STALLING INDEX
        #use content_transcripts_mod for gb previously calculated.
        #it depends on total reads, so the division of boxplot list must be carried out later
        
        totallist_SI[[2]][[i]][[k]]=totallist_boxplot[[2]][[i]][[k]]+totallist_boxplot[[3]][[i]][[k]]#calculate with content_transcripts_mod, + totallist_boxplot[[3]][[i]][[k]] (TES part)
        #else, divide by the GB + TES interval:
        # TSS      GB           TES
        #....|-------------|-----.-----|


        totallist_SI[[1]][[i]][[k]]=totallist_boxplot[[1]][[i]][[k]]

      }
    }


    cols=colors_list[1:(length(input$genelistsforDynamics)*length(input$BAMforDynamics))]

    toplot$dynamics$cols=cols
    toplot$dynamics$totallist_profile=totallist_profile
    toplot$dynamics$totallist_boxplot=totallist_boxplot
    toplot$dynamics$totallist_SI=totallist_SI
    #widths for normalization
    toplot$dynamics$totallist_lengths_forNorm_profile=totallist_lengths_forNorm_profile
    toplot$dynamics$totallist_lengths_forNorm_SI_GB=totallist_lengths_forNorm_SI_GB
    toplot$dynamics$totallist_lengths_forNorm_box_SI_prom=totallist_lengths_forNorm_box_SI_prom
    toplot$dynamics$totallist_lengths_forNorm_box_GB=totallist_lengths_forNorm_box_GB
    toplot$dynamics$totallist_lengths_forNorm_box_TES=totallist_lengths_forNorm_box_TES
    totalnames=c()
    for(i in 1:length(totallist_profile)){
      provv=names(totallist_profile[[i]])
      #add the number for each genelist (lengths_rois in position i)
      provv=paste(names(totallist_profile)[i],"(",lengths_rois[i],";",lengths_uniquegeneIDs[i],"unique ID); ",provv,sep="")
      totalnames=c(totalnames,provv)
    }
    toplot$dynamics$totalnames=totalnames






    #plot the profiles, boxplots, stalling index
    print("Plotting metagene profile on genes")

    #general options for all plots: log2, normalization, color
    output$show_islogforDynamics<-renderUI({
      checkboxInput("islogforDynamics", label="log2",value = FALSE, width = NULL)
    })


    output$show_chooseNormalizationforDynamics<-renderUI({
          radioButtons("chooseNormalizationforDynamics","Choose normalization:",
                                                choiceNames=list(
                                                  htmlhelp("Total reads","help_dynamicsOnGenes_parameters_totread"),
                                                  htmlhelp("Read density (reads/bp)","help_dynamicsOnGenes_parameters_readdensity")
                                                  ),choiceValues=list("totread","readdensity"),selected="readdensity"
                      )
    })
    output$show_colorschemeDynamics<-renderUI({
      radioButtons("colorschemeDynamics",label="Choose color scheme:",choices=c(
                                                "Random colors"="random",
                                                "Custom colors"="custom"
                                                      ),selected="random")       
    }) 



    #custom option for profile (median or mean)
    output$show_chooseMetricforDynamics<-renderUI({
        radioButtons("chooseMetricforDynamics",label="Mean or median?",choices=c(
                                                "median"="median",
                                                "mean"="mean"
                                                      ),selected="mean")
    })
    #profiles
    output$plotProfileDynamics<-renderPlot({
      #here normalize totallist_profile elements (i) based on lengths
      #totallist_lengths_forNorm_profile. Each element of totallist_profile[i] must be divided by
      #the array: totallist_lengths_forNorm_profile[[i]]. Then mean or median
      profileList_to_plot=totallist_profile
      if (isvalid(input$chooseMetricforDynamics)){
        metric=input$chooseMetricforDynamics
      }else{
        metric="mean"
      }
      if (isvalid(input$islogforDynamics)){
        islog=input$islogforDynamics
      }else{
        islog=FALSE
      }
      if (isvalid(input$chooseNormalizationforDynamics)){
        normalization=input$chooseNormalizationforDynamics
      }else{
        normalization="readdensity"
      }

      if (isvalid(input$colorschemeDynamics)){
        if (input$colorschemeDynamics=="random"){
          colors=cols
        }else{
          tosearch=paste("colorCustomDynamics",1:length(totalnames),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if (length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=cols
            }
          }else{
            colors=cols
          }
        }
      }else{
        colors=cols
      }

      if(normalization=="totread"){
        ylabel="Total reads"
      }else{
        for (i in 1:length(profileList_to_plot)){
          for(k in 1:length(profileList_to_plot[[i]])){
            profileList_to_plot[[i]][[k]]=profileList_to_plot[[i]][[k]]/totallist_lengths_forNorm_profile[[i]]
          }
        }
        ylabel="Read density (reads/bp)"
      }

      for (i in 1:length(profileList_to_plot)){
        for(k in 1:length(profileList_to_plot[[i]])){
          profileList_to_plot[[i]][[k]]=apply(profileList_to_plot[[i]][[k]],2,metric)
        }
      }

      #find max and min of the profiles 
      mins=c()
      maxs=c()
      for(i in 1:length(profileList_to_plot)){
        for(k in 1:length(profileList_to_plot[[i]])){
          if(islog){
            maxs=c(maxs,max( log2(profileList_to_plot[[i]][[k]])))
            mins=c(mins,min( log2(profileList_to_plot[[i]][[k]])[!is.infinite(log2(profileList_to_plot[[i]][[k]]))]  ))
          }else{
            maxs=c(maxs,max(profileList_to_plot[[i]][[k]]))
            mins=c(mins,min(profileList_to_plot[[i]][[k]]))            
          }

        }
      }
      maxvalue=max(maxs)
      minvalue=min(mins[!is.infinite(mins)])
      if(islog){
        ylabel=paste("log2",ylabel)
      }

      plot_dynamics_profile(profileList_to_plot=profileList_to_plot,nbin=toplot$dynamics$nbin,islog=islog,
                        minval=minvalue,maxval=maxvalue,ylabel=ylabel,totalnames=totalnames,colors=colors) 
    
    })





    #boxplots
    output$plotboxTSSDynamics<-renderPlot({
      currentlist=totallist_boxplot[[1]]
      if (isvalid(input$islogforDynamics)){
        islog=input$islogforDynamics
      }else{
        islog=FALSE
      }
      if (isvalid(input$chooseNormalizationforDynamics)){
        normalization=input$chooseNormalizationforDynamics
      }else{
        normalization="readdensity"
      }
      if (isvalid(input$colorschemeDynamics)){
        if (input$colorschemeDynamics=="random"){
          colors=cols
        }else{
          tosearch=paste("colorCustomDynamics",1:length(totalnames),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if (length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=cols
            }
          }else{
            colors=cols
          }
        }
      }else{
        colors=cols
      }

      if(normalization=="totread"){
        ylabel="Total reads"
      }else{
        for (i in 1:length(currentlist)){
          for(k in 1:length(currentlist[[i]])){
            currentlist[[i]][[k]]=currentlist[[i]][[k]]/totallist_lengths_forNorm_box_SI_prom[[i]]
          }
        }
        ylabel="Read density (reads/bp)"
      }

      plot_dynamics_box(currentlist=currentlist,islog=islog,ylabel=ylabel,colors=colors,
                    totalnames=totalnames,main="TSS")
    })
    #download data of boxplot of TSS
    output$saveboxdatadynamicsTSS=renderUI({downloadButton('saveenrichmentBoxdynamicsTSSdata', 'Download data')})



    output$plotboxGBDynamics<-renderPlot({
      currentlist=totallist_boxplot[[2]]

      if (isvalid(input$islogforDynamics)){
        islog=input$islogforDynamics
      }else{
        islog=FALSE
      }
      if (isvalid(input$chooseNormalizationforDynamics)){
        normalization=input$chooseNormalizationforDynamics
      }else{
        normalization="readdensity"
      }
      if (isvalid(input$colorschemeDynamics)){
        if (input$colorschemeDynamics=="random"){
          colors=cols
        }else{
          tosearch=paste("colorCustomDynamics",1:length(totalnames),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if (length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=cols
            }
          }else{
            colors=cols
          }
        }
      }else{
        colors=cols
      }

      if(normalization=="totread"){
        ylabel="Total reads"
      }else{
        for (i in 1:length(currentlist)){
          for(k in 1:length(currentlist[[i]])){
            currentlist[[i]][[k]]=currentlist[[i]][[k]]/totallist_lengths_forNorm_box_GB[[i]]
          }
        }
        ylabel="Read density (reads/bp)"
      }
      plot_dynamics_box(currentlist=currentlist,islog=islog,ylabel=ylabel,colors=colors,
                    totalnames=totalnames,main="Genebody")
    })

    #download data of boxplot of GB
    output$saveboxdatadynamicsGB=renderUI({downloadButton('saveenrichmentBoxdynamicsGBdata', 'Download data')})




    output$plotboxTESDynamics<-renderPlot({
      currentlist=totallist_boxplot[[3]]
      main="TES"
      if (isvalid(input$islogforDynamics)){
        islog=input$islogforDynamics
      }else{
        islog=FALSE
      }
      if (isvalid(input$chooseNormalizationforDynamics)){
        normalization=input$chooseNormalizationforDynamics
      }else{
        normalization="readdensity"
      }
      if (isvalid(input$colorschemeDynamics)){
        if (input$colorschemeDynamics=="random"){
          colors=cols
        }else{
          tosearch=paste("colorCustomDynamics",1:length(totalnames),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if (length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=cols
            }
          }else{
            colors=cols
          }
        }
      }else{
        colors=cols
      }

      if(normalization=="totread"){
        ylabel="Total reads"
      }else{
        for (i in 1:length(currentlist)){
          for(k in 1:length(currentlist[[i]])){
            currentlist[[i]][[k]]=currentlist[[i]][[k]]/totallist_lengths_forNorm_box_TES[[i]]
          }
        }
        ylabel="Read density (reads/bp)"
      }
      plot_dynamics_box(currentlist=currentlist,islog=islog,ylabel=ylabel,colors=colors,
                    totalnames=totalnames,main="TES")
    })
    #download data of boxplot of TES
    output$saveboxdatadynamicsTES=renderUI({downloadButton('saveenrichmentBoxdynamicsTESdata', 'Download data')})
    



    #specific option for cumulative plots (quantile threshold)
    output$show_percentageOutlayerCumulPlots<-renderUI({
      list(
        list(HTML("<b>Fraction of outliers to exclude in cumulative plots:</b>"),htmlhelp("","help_dynamicsOnGenes_parameters_fractionexclude")),
        sliderInput('percentageOutlayerCumulPlots',label=NULL,min = 0, max = 0.3, value = 0.05,step=0.01)
      )
    })


    output$plotSITSSDynamics<-renderPlot({
      if(isvalid(input$percentageOutlayerCumulPlots)){
        outlayer_thresh=input$percentageOutlayerCumulPlots
      }else{
        outlayer_thresh=0.05
      }
      if (isvalid(input$islogforDynamics)){
        islog=input$islogforDynamics
      }else{
        islog=FALSE
      }
      if (isvalid(input$chooseNormalizationforDynamics)){
        normalization=input$chooseNormalizationforDynamics
      }else{
        normalization="readdensity"
      }
      if (isvalid(input$colorschemeDynamics)){
        if (input$colorschemeDynamics=="random"){
          colors=cols
        }else{
          tosearch=paste("colorCustomDynamics",1:length(totalnames),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if (length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=cols
            }
          }else{
            colors=cols
          }
        }
      }else{
        colors=cols
      }

      currentlist=totallist_SI[[1]]
      if(normalization=="totread"){
        ylabel="Total reads"
      }else{
        for (i in 1:length(currentlist)){
          for(k in 1:length(currentlist[[i]])){
            currentlist[[i]][[k]]=currentlist[[i]][[k]]/totallist_lengths_forNorm_box_SI_prom[[i]]
          }
        }
        ylabel="Read density (reads/bp)"
      }


      plot_dynamics_cumulative(currentlist=currentlist,islog=islog,ylabel=ylabel,
                              outlayer_thresh=outlayer_thresh,colors=colors,main="TSS")
    })




    output$plotSIGBDynamics<-renderPlot({
      if(isvalid(input$percentageOutlayerCumulPlots)){
        outlayer_thresh=input$percentageOutlayerCumulPlots
      }else{
        outlayer_thresh=0.05
      }
      if (isvalid(input$islogforDynamics)){
        islog=input$islogforDynamics
      }else{
        islog=FALSE
      }
      if (isvalid(input$chooseNormalizationforDynamics)){
        normalization=input$chooseNormalizationforDynamics
      }else{
        normalization="readdensity"
      }
      if (isvalid(input$colorschemeDynamics)){
        if (input$colorschemeDynamics=="random"){
          colors=cols
        }else{
          tosearch=paste("colorCustomDynamics",1:length(totalnames),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if (length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=cols
            }
          }else{
            colors=cols
          }
        }
      }else{
        colors=cols
      }

      currentlist=totallist_SI[[2]]
      
      if(normalization=="totread"){
        ylabel="Total reads"
      }else{
        for (i in 1:length(currentlist)){
          for(k in 1:length(currentlist[[i]])){
            currentlist[[i]][[k]]=currentlist[[i]][[k]]/totallist_lengths_forNorm_SI_GB[[i]]
          }
        }
        ylabel="Read density (reads/bp)"
      }
      plot_dynamics_cumulative(currentlist=currentlist,islog=islog,ylabel=ylabel,
                              outlayer_thresh=outlayer_thresh,colors=colors,main="Genebody")
    })





    output$plotSISIDynamics<-renderPlot({
      if(isvalid(input$percentageOutlayerCumulPlots)){
        outlayer_thresh=input$percentageOutlayerCumulPlots
      }else{
        outlayer_thresh=0.05
      }
      if (isvalid(input$islogforDynamics)){
        islog=input$islogforDynamics
      }else{
        islog=FALSE
      }
      if (isvalid(input$chooseNormalizationforDynamics)){
        normalization=input$chooseNormalizationforDynamics
      }else{
        normalization="readdensity"
      }
      if (isvalid(input$colorschemeDynamics)){
        if (input$colorschemeDynamics=="random"){
          colors=cols
        }else{
          tosearch=paste("colorCustomDynamics",1:length(totalnames),sep="") 
          if(length(tosearch)>0){
            listinputcols=list()
            for(i in 1:length(tosearch)){
              listinputcols[[i]]=input[[tosearch[i]]]
            }
            listinputcols=unlist(listinputcols)
            if (length(listinputcols)==length(tosearch)){
              colors=listinputcols
            }else{
              colors=cols
            }
          }else{
            colors=cols
          }
        }
      }else{
        colors=cols
      }

      currentlistTSS=totallist_SI[[1]]
      currentlistGB=totallist_SI[[2]]

      currentlist=rep(list(NA),length(currentlistTSS[[1]]))
      names(currentlist)=names(currentlistTSS[[1]])
      currentlist=rep(list(currentlist),length(currentlistTSS))
      names(currentlist)=names(currentlistTSS)

      if(normalization=="totread"){
        for (i in 1:length(currentlist)){
          for(k in 1:length(currentlist[[i]])){
            pseudocount=currentlistGB[[i]][[k]]
            pseudocount=pseudocount[pseudocount!=0]
            pseudocount=quantile(pseudocount,0.0001)/10
            currentlist[[i]][[k]]=currentlistTSS[[i]][[k]]/(currentlistGB[[i]][[k]]+pseudocount)#calculate stalling index            
          }
        }
        
      }else{
        for (i in 1:length(currentlist)){
          for(k in 1:length(currentlist[[i]])){
            normTSS=currentlistTSS[[i]][[k]]/totallist_lengths_forNorm_box_SI_prom[[i]]
            normGB=currentlistGB[[i]][[k]]/totallist_lengths_forNorm_SI_GB[[i]]
            pseudocount=normGB
            pseudocount=pseudocount[pseudocount!=0]
            pseudocount=quantile(pseudocount,0.0001)/10            
            currentlist[[i]][[k]]=normTSS/(normGB+pseudocount)#calculate stalling index            
          }
        }
      }      

      plot_dynamics_cumulative(currentlist=currentlist,islog=islog,ylabel="stalling index",
                              outlayer_thresh=outlayer_thresh,colors=colors,main="Stalling index")
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
    output$saveboxdatadynamicsTSS=renderUI({NULL})
    output$saveboxdatadynamicsGB=renderUI({NULL})
    output$saveboxdatadynamicsTES=renderUI({NULL})
    output$show_percentageOutlayerCumulPlots<-renderUI({NULL})
    output$show_chooseMetricforDynamics<-renderUI({NULL})
    output$show_islogforDynamics<-renderUI({NULL})
    output$show_chooseNormalizationforDynamics<-renderUI({NULL})
    output$show_colorschemeDynamics<-renderUI({NULL})
    output$show_colorsDynamics<-renderUI({NULL})
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
    profileList_to_plot=toplot$dynamics$totallist_profile
    #colors:
    if (input$colorschemeDynamics=="random"){
      colors=toplot$dynamics$cols
    }else{
      tosearch=paste("colorCustomDynamics",1:length(toplot$dynamics$totalnames),sep="") 
      if(length(tosearch)>0){
        listinputcols=list()
        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }
        listinputcols=unlist(listinputcols)
        if (length(listinputcols)==length(tosearch)){
          colors=listinputcols
        }else{
          colors=toplot$dynamics$cols
        }
      }else{
        colors=toplot$dynamics$cols
      }
    }
    #normalization
    if(input$chooseNormalizationforDynamics=="totread"){
      ylabel="Total reads"
    }else{
      for (i in 1:length(profileList_to_plot)){
        for(k in 1:length(profileList_to_plot[[i]])){
          profileList_to_plot[[i]][[k]]=profileList_to_plot[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_profile[[i]]
        }
      }
      ylabel="Read density (reads/bp)"
    }
    #mean vs median
    for (i in 1:length(profileList_to_plot)){
      for(k in 1:length(profileList_to_plot[[i]])){
        profileList_to_plot[[i]][[k]]=apply(profileList_to_plot[[i]][[k]],2,input$chooseMetricforDynamics)
      }
    }

    #find max and min of the profiles 
    mins=c()
    maxs=c()
    for(i in 1:length(profileList_to_plot)){
      for(k in 1:length(profileList_to_plot[[i]])){
        if(input$islogforDynamics){
          maxs=c(maxs,max( log2(profileList_to_plot[[i]][[k]])))
          mins=c(mins,min( log2(profileList_to_plot[[i]][[k]])[!is.infinite(log2(profileList_to_plot[[i]][[k]]))]  ))
        }else{
          maxs=c(maxs,max(profileList_to_plot[[i]][[k]]))
          mins=c(mins,min(profileList_to_plot[[i]][[k]]))            
        }

      }
    }
    maxvalue=max(maxs)
    minvalue=min(mins[!is.infinite(mins)])
    if(input$islogforDynamics){
      ylabel=paste("log2",ylabel)
    }


    pdf(file,width=15,height=8)
    plot_dynamics_profile(profileList_to_plot=profileList_to_plot,nbin=toplot$dynamics$nbin,islog=input$islogforDynamics,
                        minval=minvalue,maxval=maxvalue,ylabel=ylabel,totalnames=toplot$dynamics$totalnames,colors=colors)
    dev.off()
  } 
)




#save PDF of TSS boxplot
output$saveboxTSSDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('TSS_boxplot.pdf', sep='')
  },
  content=function(file) {
    currentlist=toplot$dynamics$totallist_boxplot[[1]]
    if (input$colorschemeDynamics=="random"){
      colors=toplot$dynamics$cols
    }else{
      tosearch=paste("colorCustomDynamics",1:length(toplot$dynamics$totalnames),sep="") 
      if(length(tosearch)>0){
        listinputcols=list()
        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }
        listinputcols=unlist(listinputcols)
        if (length(listinputcols)==length(tosearch)){
          colors=listinputcols
        }else{
          colors=toplot$dynamics$cols
        }
      }else{
        colors=toplot$dynamics$cols
      }
    }

    if(input$chooseNormalizationforDynamics=="totread"){
      ylabel="Total reads"
    }else{
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_SI_prom[[i]]
        }
      }
      ylabel="Read density (reads/bp)"
    }

    pdf(file)
    plot_dynamics_box(currentlist=currentlist,islog=input$islogforDynamics,ylabel=ylabel,
          colors=colors,totalnames=toplot$dynamics$totalnames,main="TSS",ispdf=TRUE)
    dev.off()
  } 
)





#save PDF of GB boxplot
output$saveboxGBDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('GB_boxplot.pdf', sep='')
  },
  content=function(file) {
    
    currentlist=toplot$dynamics$totallist_boxplot[[2]]
    if (input$colorschemeDynamics=="random"){
      colors=toplot$dynamics$cols
    }else{
      tosearch=paste("colorCustomDynamics",1:length(toplot$dynamics$totalnames),sep="") 
      if(length(tosearch)>0){
        listinputcols=list()
        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }
        listinputcols=unlist(listinputcols)
        if (length(listinputcols)==length(tosearch)){
          colors=listinputcols
        }else{
          colors=toplot$dynamics$cols
        }
      }else{
        colors=toplot$dynamics$cols
      }
    }

    if(input$chooseNormalizationforDynamics=="totread"){
      ylabel="Total reads"
    }else{
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_GB[[i]]
        }
      }
      ylabel="Read density (reads/bp)"
    }

    pdf(file)
    plot_dynamics_box(currentlist=currentlist,islog=input$islogforDynamics,ylabel=ylabel,
          colors=colors,totalnames=toplot$dynamics$totalnames,main="Genebody",ispdf=TRUE)
    dev.off()
  } 
)



#save PDF of TES boxplot
output$saveboxTESDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('TES_boxplot.pdf', sep='')
  },
  content=function(file) {
    
    
    currentlist=toplot$dynamics$totallist_boxplot[[3]]


    if (input$colorschemeDynamics=="random"){
      colors=toplot$dynamics$cols
    }else{
      tosearch=paste("colorCustomDynamics",1:length(toplot$dynamics$totalnames),sep="") 
      if(length(tosearch)>0){
        listinputcols=list()
        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }
        listinputcols=unlist(listinputcols)
        if (length(listinputcols)==length(tosearch)){
          colors=listinputcols
        }else{
          colors=toplot$dynamics$toplot$dynamics$cols
        }
      }else{
        colors=toplot$dynamics$cols
      }
    }

    if(input$chooseNormalizationforDynamics=="totread"){
      ylabel="Total reads"
    }else{
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_TES[[i]]
        }
      }
      ylabel="Read density (reads/bp)"
    }

    pdf(file)
    plot_dynamics_box(currentlist=currentlist,islog=input$islogforDynamics,ylabel=ylabel,
          colors=colors,totalnames=toplot$dynamics$totalnames,main="TES",ispdf=TRUE)
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

    if(input$chooseNormalizationforDynamics=="readdensity"){
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_SI_prom[[i]]
        }
      }
    }

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

    if(input$chooseNormalizationforDynamics=="readdensity"){
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_GB[[i]]
        }
      }
    }

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

    if(input$chooseNormalizationforDynamics=="readdensity"){
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_TES[[i]]
        }
      }
    }

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

    if (input$colorschemeDynamics=="random"){
      colors=toplot$dynamics$cols
    }else{
      tosearch=paste("colorCustomDynamics",1:length(toplot$dynamics$totalnames),sep="") 
      if(length(tosearch)>0){
        listinputcols=list()
        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }
        listinputcols=unlist(listinputcols)
        if (length(listinputcols)==length(tosearch)){
          colors=listinputcols
        }else{
          colors=toplot$dynamics$cols
        }
      }else{
        colors=toplot$dynamics$cols
      }
    }

    if(input$chooseNormalizationforDynamics=="totread"){
      ylabel="Total reads"
    }else{
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_SI_prom[[i]]
        }
      }
      ylabel="Read density (reads/bp)"
    }


    pdf(file)
    plot_dynamics_cumulative(currentlist=currentlist,islog=input$islogforDynamics,ylabel=ylabel,
                              outlayer_thresh=input$percentageOutlayerCumulPlots,colors=colors,main="TSS",ispdf=TRUE)     
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

    if (input$colorschemeDynamics=="random"){
      colors=toplot$dynamics$cols
    }else{
      tosearch=paste("colorCustomDynamics",1:length(toplot$dynamics$totalnames),sep="") 
      if(length(tosearch)>0){
        listinputcols=list()
        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }
        listinputcols=unlist(listinputcols)
        if (length(listinputcols)==length(tosearch)){
          colors=listinputcols
        }else{
          colors=toplot$dynamics$cols
        }
      }else{
        colors=toplot$dynamics$cols
      }
    }

    if(input$chooseNormalizationforDynamics=="totread"){
      ylabel="Total reads"
    }else{
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          currentlist[[i]][[k]]=currentlist[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_SI_GB[[i]]
        }
      }
      ylabel="Read density (reads/bp)"
    }


    pdf(file)
    plot_dynamics_cumulative(currentlist=currentlist,islog=input$islogforDynamics,ylabel=ylabel,
                              outlayer_thresh=input$percentageOutlayerCumulPlots,colors=colors,main="Genebody",ispdf=TRUE)     
    dev.off()
  } 
)




#save PDF of SI on SI
output$saveSISIDynamicsbutton<- downloadHandler(
  filename=function() {
      paste('SI_SI_dynamics.pdf', sep='')
  },
  content=function(file) {
    currentlistTSS=toplot$dynamics$totallist_SI[[1]]
    currentlistGB=toplot$dynamics$totallist_SI[[2]]

    if (input$colorschemeDynamics=="random"){
      colors=toplot$dynamics$cols
    }else{
      tosearch=paste("colorCustomDynamics",1:length(toplot$dynamics$totalnames),sep="") 
      if(length(tosearch)>0){
        listinputcols=list()
        for(i in 1:length(tosearch)){
          listinputcols[[i]]=input[[tosearch[i]]]
        }
        listinputcols=unlist(listinputcols)
        if (length(listinputcols)==length(tosearch)){
          colors=listinputcols
        }else{
          colors=toplot$dynamics$cols
        }
      }else{
        colors=toplot$dynamics$cols
      }
    }

    currentlist=rep(list(NA),length(currentlistTSS[[1]]))
    names(currentlist)=names(currentlistTSS[[1]])
    currentlist=rep(list(currentlist),length(currentlistTSS))
    names(currentlist)=names(currentlistTSS)
    if(input$chooseNormalizationforDynamics=="totread"){
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          pseudocount=currentlistGB[[i]][[k]]
          pseudocount=pseudocount[pseudocount!=0]
          pseudocount=quantile(pseudocount,0.0001)/10
          currentlist[[i]][[k]]=currentlistTSS[[i]][[k]]/(currentlistGB[[i]][[k]]+pseudocount)#calculate stalling index            
        }
      }
      
    }else{
      for (i in 1:length(currentlist)){
        for(k in 1:length(currentlist[[i]])){
          normTSS=currentlistTSS[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_box_SI_prom[[i]]
          normGB=currentlistGB[[i]][[k]]/toplot$dynamics$totallist_lengths_forNorm_SI_GB[[i]]
          pseudocount=normGB
          pseudocount=pseudocount[pseudocount!=0]
          pseudocount=quantile(pseudocount,0.0001)/10            
          currentlist[[i]][[k]]=normTSS/(normGB+pseudocount)#calculate stalling index            
        }
      }
    }      


    pdf(file)
    plot_dynamics_cumulative(currentlist=currentlist,islog=input$islogforDynamics,ylabel="stalling index",
                              outlayer_thresh=input$percentageOutlayerCumulPlots,colors=colors,main="Stalling index",ispdf=TRUE)    
    dev.off()
  } 
)
















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




#observer for help buttons of parameters
observeEvent(input$help_goAnalysis_parameters_fromROI, {boxHelpServer(help_goAnalysis_parameters_fromROI)})
observeEvent(input$help_goAnalysis_parameters_fromgenelist, {boxHelpServer(help_goAnalysis_parameters_fromgenelist)})
observeEvent(input$help_goAnalysis_parameters_kindofID, {boxHelpServer(help_goAnalysis_parameters_kindofID)})
observeEvent(input$help_goAnalysis_parameters_nearestgenes, {boxHelpServer(help_goAnalysis_parameters_nearestgenes)})
observeEvent(input$help_goAnalysis_parameters_genewindow, {boxHelpServer(help_goAnalysis_parameters_genewindow)})
observeEvent(input$help_goAnalysis_parameters_signatures, {boxHelpServer(help_goAnalysis_parameters_signatures)})
observeEvent(input$help_goAnalysis_parameters_orderresults, {boxHelpServer(help_goAnalysis_parameters_orderresults)})
observeEvent(input$help_goAnalysis_parameters_minsize, {boxHelpServer(help_goAnalysis_parameters_minsize)})
observeEvent(input$help_goAnalysis_parameters_maxsize, {boxHelpServer(help_goAnalysis_parameters_maxsize)})
observeEvent(input$help_goAnalysis_parameters_generatio, {boxHelpServer(help_goAnalysis_parameters_generatio)})
observeEvent(input$help_goAnalysis_parameters_padjthresh, {boxHelpServer(help_goAnalysis_parameters_padjthresh)})







# toListenGO <- reactive({
#     list(input$doTheGO)#,input$doTheGO2)
# })

observeEvent(input$doTheGO,{
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
    if(isolate(input$chooseSourceGO)=="fromROI"){

      if(!isvalid(isolate(input$selectROIGO))){
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
          text = "Be sure to select at least one ROI and that a genome assembly is active. Go to 'Assembly' section for that",
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


        #here put parameters only for heatmap (color scale and color)
        output$show_colorScaleGO<-renderUI({
          selectInput("colorScaleGO",label="Choose a color scale:",c("white/red"="white_red4",
                                                              "white/blue"="white_blue",
                                                              "white/green"= "white_green4"))
        })
        output$show_scaleQuantileGO<-renderUI({
          sliderInput('scaleQuantileGO',label="Quantile threshold for padj colorscale",min = 0.1, max = 1, value = 0.9,step=0.002)
        })
        
                        

        output$plotOntology<-renderPlot({
          if (isvalid(input$colorScaleGO)){
            colorpal=input$colorScaleGO
          }else{
            colorpal="white_red4"
          }
          if(isvalid(input$scaleQuantileGO)){
            scaleQuantile=input$scaleQuantileGO
          }else{
            scaleQuantile=0.9
          }

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
          if (isvalid(input$colorScaleGO)){
            palette_col=input$colorScaleGO
          }else{
            palette_col="white_red4"
          }
          if(isvalid(input$scaleQuantileGO)){
            scaleQuantile=input$scaleQuantileGO
          }else{
            scaleQuantile=0.9
          }   

          trasp[trasp>quantile(trasp,scaleQuantile)]=quantile(trasp,scaleQuantile)          
          brk=201
            
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
        #here only one block: plot barplot and not heatmap
        output$show_colorScaleGO<-renderUI({NULL})
        output$show_scaleQuantileGO<-renderUI({NULL})
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
      output$plotOntology<-renderPlot({plot_text(text="There are no significant results,\nor thresholds (p adjusted, gene ratio) are too stringent.\nTry to relax them.",cex=1.4)})
      output$tableOntology <- renderDataTable({NULL})
      output$textNameGO <- renderPlot({NULL})
      output$tableGOdownloadButton<-renderUI({NULL})
      output$plotMaterialLeft<-renderPlot({NULL})
      output$colorScaleGO<-renderPlot({NULL})
      output$saveheatmapGO<-renderUI({NULL})
      output$show_colorScaleGO<-renderUI({NULL})
      output$show_scaleQuantileGO<-renderUI({NULL})

    }    
  }else{
    output$plotOntology<-renderPlot({NULL})
    output$tableOntology <- renderDataTable({NULL})   
    output$textNameGO <- renderPlot({NULL}) 
    output$tableGOdownloadButton<-renderUI({NULL})
    output$plotMaterialLeft<-renderPlot({NULL})
    output$colorScaleGO<-renderPlot({NULL})
    output$saveheatmapGO<-renderUI({NULL})
    output$show_colorScaleGO<-renderUI({NULL})
    output$show_scaleQuantileGO<-renderUI({NULL})

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







