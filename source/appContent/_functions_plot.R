###############################################################################
# special: plot only text
###############################################################################

#lines can be separated with \n
plot_text<-function(text,cex,color="black") {
	par(mar = c(0, 0, 0, 0)) 
	plot(x = 0:1,                   # Create empty plot
     y = 0:1,
     ann = F,
     bty = "n",
     type = "n",
     xaxt = "n",
     yaxt = "n")

	text(x = 0.5,                   # Add text to empty plot
     y = 0.5,
     text, 
     cex = cex,
     col=color)
}


###############################################################################
# single evaluation
###############################################################################



#single evaluation pie
plot_singleeval_pie<-function(elements,colors,title,label){
	par(mar=c(7,7,7,7))
	pie(elements,col=colors,main=title,cex.main=1.2,labels=label,cex=1)		
}

#single evaluation bar
plot_singleeval_bar<-function(elements,colors,title,label){
    par(mar=c(10,6,5,5))
    names(elements)=label
    barplot(elements,col=colors,main=title,cex.main=1.2,las=2)
    title(ylab = "Interval number", cex.lab = 1,
              line = 4)
}

#single evaluation density
plot_singleeval_density<-function(range,range_promo,range_intra,range_inter,colors) {
	quant=0.95
    cuttedga=range[width(range)<quantile(width(range),quant)]
    cuttedprom=range_promo[width(range_promo)<quantile(width(range_promo),quant)]
    cuttedintra=range_intra[width(range_intra)<quantile(width(range_intra),quant)]
    cuttedinter=range_inter[width(range_inter)<quantile(width(range_inter),quant)]

	if (length(range_intra)>2){
	  if(length(range_inter)>2){
	    xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x,density(log2(width(cuttedinter)))$x), 
	            max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x,density(log2(width(cuttedinter)))$x)  )
	    ylim=c(0, max(max(density(log2(width(cuttedprom)))$y),max(density(log2(width(cuttedintra)))$y),max(density(log2(width(cuttedinter)))$y) ) )
	    plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 width",cex.lab=1,cex.axis=1,ylab="probability",ylim=ylim,xlim=xlim,
	      main=paste("Interval width (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=3,type="l",col=colors[1])
	    lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=3)
	    lines(density(log2(width(cuttedintra)))$x,density(log2(width(cuttedintra)))$y,cex=1,col=colors[3],ylim=ylim,xlim=xlim,lty=1,lwd=3)
	    lines(density(log2(width(cuttedinter)))$x,density(log2(width(cuttedinter)))$y,cex=1,col=colors[4],ylim=ylim,xlim=xlim,lty=1,lwd=3)
	  }else{
	    xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x), 
	            max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x,density(log2(width(cuttedintra)))$x  ))
	    ylim=c(0, max(max(density(log2(width(cuttedprom)))$y),max(density(log2(width(cuttedintra)))$y) ) )
	    plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 width",cex.lab=1,cex.axis=1,ylab="probability",ylim=ylim,xlim=xlim,
	      main=paste("Interval width (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=3,type="l",col=colors[1])
	    lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=3)
	    lines(density(log2(width(cuttedintra)))$x,density(log2(width(cuttedintra)))$y,cex=1,col=colors[3],ylim=ylim,xlim=xlim,lty=1,lwd=3)
	  }
	}else{
	  xlim=c(min(density(log2(width(cuttedga)))$x ,density(log2(width(cuttedprom)))$x), 
	            max(density(log2(width(cuttedga)))$x,density(log2(width(cuttedprom)))$x  ))
	  ylim=c(0, max(max(density(log2(width(cuttedprom)))$y) ) )
	  plot(density(log2(width(cuttedga)))$x,density(log2(width(cuttedga)))$y,pch=".",cex=1,xlab="log2 width",cex.lab=1,cex.axis=1,ylab="probability",ylim=ylim,xlim=xlim,
	    main=paste("Interval width (0 - ",quant*100,"%)",sep=''),cex.main=1.2,lwd=3,type="l",col=colors[1])
	  lines(density(log2(width(cuttedprom)))$x,density(log2(width(cuttedprom)))$y,cex=1,col=colors[2],ylim=ylim,xlim=xlim,lty=1,lwd=3)
	}
	legend("topright" , c("All intervals","Promoter","Genebody","Intergenic") , col=colors , lty=rep(1,4),lwd=3,cex=1)	
}

#single evaluation boxplot
plot_singleeval_boxplot<-function(normalization,bam,bam_promo,bam_intra,bam_inter,
									range,range_promo,range_intra,range_inter,colors) {
    par(mar=c(8,5,5,1))
    if (normalization=="totread"){
      suppressWarnings(boxplot(bam,bam_promo,bam_intra,bam_inter,notch=TRUE,outline=FALSE,varwidth=TRUE,col=colors,xaxt="n",ylab="Total reads"))
    }else{
      bam=round( bam/width(range),3)
      bam_promo=round( bam_promo/width(range_promo),3)
      bam_intra=round( bam_intra/width(range_intra),3)
      bam_inter=round( bam_inter/width(range_inter),3)
      suppressWarnings(boxplot(bam,bam_promo,bam_intra,bam_inter,notch=TRUE,outline=FALSE,varwidth=TRUE,col=colors,xaxt="n",ylab="Read density (reads/bp)"))
    }
    axis(1,at=1:4,labels=c("All ranges","Promoter","Genebody","Intergenic"),las=2)
}

#single evaluation profile
plot_singleeval_profile<-function(colors,profile_to_plot,ylab) {
    par(mar=c(6,6,3,1))
    mins=Reduce(min,profile_to_plot)
    maxs=Reduce(max,profile_to_plot)

    plot(0,type='n',xaxt="n",ylim=c(mins,maxs),xlim=c(0,50),ylab=ylab,xlab="Genomic Window",
                        main="Shape profile",cex.main=1.2,cex=1,xaxt="n",cex.lab=1.2,cex.axis=1)
    for(i in 1:length(profile_to_plot)){
      lines(profile_to_plot[[i]],pch=".",lwd=3,col=colors[i])
    }

    axis(1,at=seq(1,length(profile_to_plot[[1]]),length(profile_to_plot[[1]])/2-1),
                      labels=c( "start","midpoint","end" ),cex.axis=1) 
}










###############################################################################
# pairwise comparison
###############################################################################



plot_cmp_barplot<-function(matbar,n1,n2,colors,jaccard,margins=c(14,6,5,5),splitlegend=FALSE) {
    par(xpd=TRUE,mar=margins)
    barplot(matbar,col=c(colors[3],colors[1],colors[4],colors[2]),main=paste("Jaccard idx:",jaccard),cex.main=1.2,las=2)
    title(ylab = "Interval number", cex.lab = 1,line = 4)
    if(splitlegend){
    	plot.new()
    	legend("topleft",legend=c(paste(n1,"alone"),
                                                      paste(n1,"overlapping"),
                                                      paste(n2,"overlapping"),
                                                      paste(n2,"alone")),pch=19,
        			col=c(colors[1],colors[3],colors[4],colors[2]))	
    }else{
    	legend("bottom",inset=c(0,-1.5),legend=c(paste(n1,"alone"),
                                                      paste(n1,"overlapping"),
                                                      paste(n2,"overlapping"),
                                                      paste(n2,"alone")),pch=19,
        			col=c(colors[1],colors[3],colors[4],colors[2]))	
    }
}



plot_cmp_venn<-function(area1,area2,common,n1,n2,colors) {
	griglia<-draw.pairwise.venn(area1=area1+length(common), area2=area2+length(common), 
                cross.area=length(common),c(n1,n2),cex=1,ext.dist=c(.01,-0.04),
                ext.line.lwd=0,ext.length=0,ext.pos=90,alpha=c(0.1,0.1),col=c(colors[1],colors[2]),fill=c(colors[1],colors[2]),
                label.col=c(colors[1],colors[3],colors[4]),cat.pos=c(310,50),cat.dist=c(0.1,0.1),lwd=c(2,2),cat.cex=c(1.2,1.2),
                cat.col=c(colors[1],colors[4]),margin=0.3,main="Intervals overlap",main.cex=1.2,main.col="black",main.pos=c(1,1),filename=NULL)	
}



plot_cmp_box<-function(islog,normalization,range2only_bam1,range1only_bam1,common_bam1,
						range1only_bam2,range2only_bam2,common_bam2,n1,n2,BAM1_choose,BAM2_choose,colors) {


    #normalization must be "readdensity" or the other
    if(normalization=="readdensity"){
       lab_axis="Read density (reads/bp)"
    }else{
       lab_axis="Total reads"
    }
    if(islog){
      arraytoplot=list(log2(range2only_bam1),log2(range1only_bam1),log2(common_bam1),
            log2(range1only_bam2),log2(range2only_bam2),log2(common_bam2))
      laby=paste("log2",lab_axis)
    }else{
      arraytoplot=list(range2only_bam1,range1only_bam1,common_bam1,
            range1only_bam2,range2only_bam2,common_bam2)
      laby=lab_axis
    }
    par(xpd=TRUE,mar=c(16,4,2,2))
    suppressWarnings(boxplot(arraytoplot,col=c(colors[2],colors[1],colors[3],colors[1],colors[2],colors[4]),
            outline=FALSE,notch=TRUE,varwidth=TRUE,at =c(1,2,3, 5,6,7),xaxt="n",ylab=laby))

    axis(1,at=c(2,6),labels=c("enrich. 1","enrich. 2"),las=2,cex=1)
    legend("bottom",inset=c(0,-1),legend=c(paste(n1,"alone"),
                                                  paste(n1,"overlapping"),
                                                  paste(n2,"overlapping"),
                                                  paste(n2,"alone")),pch=19,
          col=c(colors[1],colors[3],colors[4],colors[2]))
    #another legend for complete name of enrichments
    legend("bottom",inset=c(0,-1.5),legend=c(paste("enrich. 1: ",BAM1_choose,sep=""),paste("enrich. 2: ",BAM2_choose,sep="")),box.col = "black",bg = "white")
}




plot_cmp_scatter<-function(df,normalization,subsetToShow,islog,BAM2chooseCmp,BAM1chooseCmp,n1,n2,insets=c(-0.5,-0.8),colors) {
	#subsetToShow is which subset (all, common, only 1, only 2) of the scatter to show
    if(normalization=="readdensity"){
      #determine 1 and 2 columns of df based on width. For common ranges, divide sums of signals and width
      df[df$label=="2only",1]=df[df$label=="2only",1]/df$width_range2[df$label=="2only"]
      df[df$label=="2only",2]=df[df$label=="2only",2]/df$width_range2[df$label=="2only"]
      df[df$label=="1only",1]=df[df$label=="1only",1]/df$width_range1[df$label=="1only"]
      df[df$label=="1only",2]=df[df$label=="1only",2]/df$width_range1[df$label=="1only"]
      df[df$label=="common",1]=df[df$label=="common",1]/df$width_range1[df$label=="common"]
      df[df$label=="common",2]=df[df$label=="common",2]/df$width_range2[df$label=="common"]
      lab_axis="Read density (reads/bp)"
    }else{
      #here divide by times only, not for width
      df[df$label=="common",1]=df[df$label=="common",1]/df$times_factor_hits_range1[df$label=="common"]
      df[df$label=="common",2]=df[df$label=="common",2]/df$times_factor_hits_range2[df$label=="common"]
      lab_axis="Total reads"
    }


    if(islog){
      toplot1=log2(df[,1])
      toplot2=log2(df[,2])
      #set -inf values to min
      toplot1[is.infinite(toplot1)]=min(toplot1[!is.infinite(toplot1)])
      toplot2[is.infinite(toplot2)]=min(toplot2[!is.infinite(toplot2)])
      xlims=c(min(toplot1),max(toplot1))
      ylims=c(min(toplot2),max(toplot2))
      laby=paste("log2",BAM2chooseCmp,lab_axis)
      labx=paste("log2",BAM1chooseCmp,lab_axis)
    }else{
      toplot1=df[,1]
      toplot2=df[,2]
      xlims=c(0,max(toplot1))
      ylims=c(0,max(toplot2))
      laby=paste(BAM2chooseCmp,lab_axis)
      labx=paste(BAM1chooseCmp,lab_axis)
    }

	cor_all=round(cor(toplot1,toplot2),2)
    cor_common=round( cor(toplot1[df[,3]=="common"],toplot2[df[,3]=="common"] ),2)
	pos=df[,3] %in% subsetToShow
    df2=df[pos,]
    toplot1=toplot1[pos]
    toplot2=toplot2[pos]

    par(xpd=TRUE,mar=c(15,5,2,2))
    plot(toplot1,toplot2,col=ifelse(df2[,3]=="1only",colors[1],ifelse(df2[,3]=="2only",colors[2],"black")), xlab=labx,ylab=laby,
                        xlim=xlims,ylim=ylims)
    legend("bottom",inset=c(0,insets[1]),legend=c(paste("cor common intervals: ",cor_common,sep=""),paste("cor all intervals: ",cor_all,sep="")),box.col = "black",bg = "white")
    legend("bottom",inset=c(0,insets[2]),legend=c(paste(n2,"intervals"),paste(n1,"intervals"),"common"),box.col = "black",bg = "white",pch=19,col=c(colors[2],colors[1],"black"))	
}


plot_cmp_calibration<-function(quantiles,only1,only2,combined1,combined2,n1,n2,insets=c(0,-0.8),colors) {
    par(xpd=TRUE,mar=c(13,5,2,2))
    plot(1, type="n", xlab="Enrichment threshold (quantile)", ylab="fraction of overlap", xlim=c(0, 1), ylim=c(0, 1),xaxt="n")
    
    if(!is.null(only1)& is.null(only2)){
      lines(quantiles,only1,type="p",col=colors[1])
      legend("bottomleft",inset=c(0,-0.8),legend=paste(n1,"overlap"),col=colors[1],bg="transparent",pch=19)
    }else if(is.null(only1)& !is.null(only2)){
      lines(quantiles,only2,type="p",col=colors[2])
      legend("bottomleft",inset=c(0,-0.8),legend=paste(n2,"overlap"),col=colors[2],bg="transparent",pch=19)
    }else if(!is.null(only1)& !is.null(only2)){
      lines(quantiles,only1,type="p",col=colors[1])
      lines(quantiles,only2,type="p",col=colors[2])
      lines(quantiles,combined1,type="p",col=colors[3])
      lines(quantiles,combined2,type="p",col=colors[4])
      legend("bottom",inset=insets,
              legend=c(paste(n1,"overlap"),paste(n2,"overlap"),paste(n1,"overlap combined"),paste(n2,"overlap combined")),
                  col=c(colors[1],colors[2],colors[3],colors[4]),pch=19,bg="transparent" )
    }
    axis(1,quantiles)	
}








###############################################################################
# digital heatmap
###############################################################################


plot_digital_frequency<-function(ovtoplot,nbin,fraction=TRUE,maxval,colors,splitlegend=FALSE) {
  
  if(fraction){
    ylabfreq="Fraction of overlaps"   
  }else{
    ylabfreq="Number of overlaps"
  }
  #now plot (either % or absolute numbers)
  par(xpd=TRUE,mar=c(13,5,2,2))
  plot(0,type='n',xaxt="n",ylim=c(0,maxval),xlim=c(1,nbin),ylab=ylabfreq,xlab="Genomic Window",
                main="Overlaps bias",cex.main=1.2,cex=1,xaxt="n",cex.lab=1.2,cex.axis=1)
  count=1
  legnames=c()
  for(i in 1:length(ovtoplot)){
    for(k in 1:length(ovtoplot[[i]])){
      #draw the lines/points
      lines(ovtoplot[[i]][[k]],col=colors[count],lwd=3)
      legnames[count]=paste(names(ovtoplot[[i]])[k],"in",names(ovtoplot)[i])
      count=count+1
    }
  }
  axis(1,at=seq(1,nbin,(nbin-1)/2),
              labels=c( "start","midpoint","end" ),cex.axis=1)
  #legend for each line
  if (splitlegend){
    plot.new()
    legend("topleft",legend=legnames,col=colors,bg="transparent",pch=19,cex=0.8) 
  }else{
    legend("bottom",inset=c(0,-0.8),legend=legnames,col=colors,bg="transparent",pch=19,cex=0.8) 
  }
     
}






###############################################################################
# analogic heatmap
###############################################################################


plot_analog_profile<-function(islog2,portionlist_profile,colors,ispdf=FALSE) {
  if(islog2){
    portionlist_profile2=lapply(portionlist_profile,log2)
    yl="Log2 Read density (reads/bp)"
  }else{
    portionlist_profile2=portionlist_profile
    yl="Read density (reads/bp)"
  }
  maxval=max(unlist(lapply(portionlist_profile2,max)))
  minval=min(unlist(lapply(portionlist_profile2,min)))


  if (length(portionlist_profile2)<12){
    bottommargin=length(portionlist_profile2)+4*1.2
  }else{
    if(ispdf){
      bottommargin=5
    }else{
      bottommargin=12+4*1.2
    }
    
  }
  
  par(mar=c(bottommargin,4,2,2),xpd=TRUE)
  plot(1, type="n", xlab="",xaxt="n", ylab=yl, xlim=c(1, length(portionlist_profile2[[1]])), ylim=c(minval, maxval))
  axis(1,at=c(1,length(portionlist_profile2[[1]])/2 +0.5,length(portionlist_profile2[[1]])),labels=c("start","center","end"))

  for(i in 1:length(portionlist_profile2)){
    lines(portionlist_profile2[[i]],lwd=3,col=colors[i])
  }
  #if plotted into PDF, something changes with margins, we have to adapt
  if (!ispdf){
    legend("bottom",inset=c(0,-((length(portionlist_profile2))+4)/20),legend=names(portionlist_profile),col=colors,lty=rep(1,length(colors)),cex=0.7,bg="transparent",lwd=3)  
  }else{
    if (length(portionlist_profile2)<12){
      legend("bottom",inset=c(0,-((length(portionlist_profile2))+2)/20),legend=names(portionlist_profile),col=colors,lty=rep(1,length(colors)),cex=0.7,bg="transparent",lwd=3) 
    }else{
      plot.new()
      legend("topleft",legend=names(portionlist_profile2),col=colors,lty=rep(1,length(colors)),cex=0.7,bg="transparent",lwd=3) 
    }
    
  }
  
}






plot_analog_boxByBAM<-function(materialtoplot,roinumber,newcols,newnames,bamname,islog,colors,ispdf=FALSE){

  factor_add=rep(0:(length(bamname)-1),each=roinumber)
  addingfactor=1:length(materialtoplot)+factor_add

  if(islog){
    newlist2=lapply(materialtoplot,log2)
    yl="Log2 Read density (reads/bp)"
  }else{
    newlist2=materialtoplot
    yl="Read density (reads/bp)"
  }


  if (!ispdf){
    if (length(newnames)<12){
      bottommargin=length(newnames)+4*1.2
    }else{
      bottommargin=12+4*1.2
    }
  }else{
    bottommargin=14
  }

  par(mar=c(bottommargin,4,1,1),xpd=TRUE)
  suppressWarnings(boxplot(newlist2,at=addingfactor,col=newcols,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
  ats=c()
  for(i in 1:length(bamname)){
    window=addingfactor[(((i-1)*roinumber)+1):(i*roinumber)]
    currentvalue=(window[length(window)]-window[1])/2
    currentvalue=window[1]+currentvalue
    ats=c(ats,currentvalue)
  }
  axis(1,at=ats,label=bamname,las=2)
  if (!ispdf){
    legend("bottom",inset=c(0,-((length(newnames))+4)/20),legend=newnames,col=newcols,cex=0.7,bg="transparent",pch=rep(19,roinumber*length(bamname)))
  }else{
    plot.new()
    legend("topleft",legend=newnames,col=newcols,cex=0.7,bg="transparent",pch=rep(19,roinumber*length(bamname)))
  }
  

}





plot_analog_boxByROI<-function(materialtoplot,roiname,bamname,newnames,islog,isgrouped,colors,ispdf=FALSE) {
  factor_add=rep(0:(length(roiname)-1),each=length(bamname))
  addingfactor=1:length(materialtoplot)+factor_add
  if(islog){
    portionlist_boxes2=lapply(materialtoplot,log2)
    yl="Log2 Read density (reads/bp)"
  }else{
    portionlist_boxes2=materialtoplot
    yl="Read density (reads/bp)"
  }


  if (!isgrouped){
    tocheck=newnames
  }else{
    tocheck=bamname
  }

  if (!ispdf){
    if (length(tocheck)<12){
      bottommargin=length(tocheck)+4*1.2
    }else{
      bottommargin=12+4*1.2
    }
  }else{
    bottommargin=14
  }


  par(mar=c(bottommargin,4,1,1),xpd=TRUE)
  suppressWarnings(boxplot(portionlist_boxes2,at=addingfactor,col=colors,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
  ats=c()
  for(i in 1:length(roiname)){
    window=addingfactor[(((i-1)*length(bamname))+1):(i*length(bamname))]
    currentvalue=(window[length(window)]-window[1])/2
    currentvalue=window[1]+currentvalue
    ats=c(ats,currentvalue)
  }


  if (!ispdf){
    if(!isgrouped){
      axis(1,at=ats,label=roiname,las=2)
      legend("bottom",inset=c(0,-((length(tocheck))+4)/20),legend=newnames,col=colors,cex=0.7,bg="transparent",pch=rep(19,length(materialtoplot)))          
    }else{
      axis(1,at=ats,label=newnames,las=2)
      legend("topright",legend=bamname,col=colors,cex=0.7,bg="transparent",pch=rep(19,length(materialtoplot)))
    }    
  }else{

      if(!isgrouped){
        axis(1,at=ats,label=roiname,las=2)
        plot.new()
        legend("topleft",legend=newnames,col=colors,cex=0.7,bg="transparent",pch=rep(19,length(materialtoplot)))
      }else{
        axis(1,at=ats,label=newnames,las=2)
        plot.new()
        legend("topleft",legend=bamname,col=colors,cex=0.7,bg="transparent",pch=rep(19,length(materialtoplot)))        
      }

  }

}







###############################################################################
# profiles and box
###############################################################################


plot_profilesAndBox_profile<-function(islog2,portionlist_profile,colors,yl,ispdf=FALSE) {
  if(islog2){
    portionlist_profile=lapply(portionlist_profile,log2)
    yl=paste("Log2",yl)
  }else{
  }


  if (length(portionlist_profile)<12){
    bottommargin=length(portionlist_profile)+4*1.2
  }else{
    if(ispdf){
      bottommargin=5
    }else{
      bottommargin=12+4*1.2
    }
    
  }


  maxval=max(unlist(lapply(portionlist_profile,max)))
  if(min(unlist(lapply(portionlist_profile,min))) <0){
    minval=min(unlist(lapply(portionlist_profile,min)))
  }else{
    minval=0
  }
  par(mar=c(bottommargin,4,1,1),xpd=TRUE)
  plot(1, type="n", xlab="",xaxt="n", ylab=yl, xlim=c(1, length(portionlist_profile[[1]])), ylim=c(minval, maxval))
  axis(1,at=c(1,length(portionlist_profile[[1]])/2 +0.5,length(portionlist_profile[[1]])),labels=c("start","center","end"))
  for(i in 1:length(portionlist_profile)){
    lines(portionlist_profile[[i]],lwd=3,col=colors[i])
  }


  #if plotted into PDF, something changes with margins, we have to adapt
  if (!ispdf){
    legend("bottom",inset=c(0,-((length(portionlist_profile))+4)/20),legend=names(portionlist_profile),col=colors,lty=rep(1,length(colors)),cex=0.6,bg="transparent",lwd=3)  
  }else{
    if (length(portionlist_profile)<12){
      legend("bottom",inset=c(0,-((length(portionlist_profile))+2)/20),legend=names(portionlist_profile),col=colors,lty=rep(1,length(colors)),cex=0.6,bg="transparent",lwd=3) 
    }else{
      plot.new()
      legend("topleft",legend=names(portionlist_profile),col=colors,lty=rep(1,length(colors)),cex=0.6,bg="transparent",lwd=3) 
    }
    
  }

}





plot_profilesAndBox_boxByBAM<-function(islog,materialtoplot,colors,yl,newnames,roinumber,bamname,ispdf=FALSE) {
  factor_add=rep(0:(length(bamname)-1),each=roinumber)
  addingfactor=1:length(materialtoplot)+factor_add

  if (!ispdf){
    if (length(newnames)<12){
      bottommargin=length(newnames)+4*1.2
    }else{
      bottommargin=12+4*1.2
    }
  }else{
    bottommargin=14
  }


  par(mar=c(bottommargin,4,1,1),xpd=TRUE)
  if(islog){
    materialtoplot=lapply(materialtoplot,log2)
    yl=paste("Log2",yl)
  }else{
  }
  suppressWarnings(boxplot(materialtoplot,at=addingfactor,col=colors,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE))
  ats=c()
  for(i in 1:length(bamname)){
    window=addingfactor[(((i-1)*roinumber)+1):(i*roinumber)]
    currentvalue=(window[length(window)]-window[1])/2
    currentvalue=window[1]+currentvalue
    ats=c(ats,currentvalue)
  }
  axis(1,at=ats,label=bamname,las=2)
  #if is pdf plot 
  if (!ispdf){
    legend("bottom",inset=c(0,-((length(newnames))+4)/20),legend=newnames,col=colors,cex=0.6,bg="transparent",pch=rep(19,roinumber*length(bamname)))  
  }else{
    plot.new()
    legend("topleft",legend=newnames,col=colors,cex=0.6,bg="transparent",pch=rep(19,roinumber*length(bamname)))
  }

}






plot_profilesAndBox_boxByROI<-function(materialtoplot,roiname,bamname,newnames,islog,yl,isgrouped,colors,ispdf=FALSE) {

  factor_add=rep(0:(length(roiname)-1),each=length(bamname))
  addingfactor=1:length(materialtoplot)+factor_add

  if(islog){
    materialtoplot=lapply(materialtoplot,log2)
    yl=paste("Log2",yl)
  }

  if (!isgrouped){
    tocheck=newnames
  }else{
    tocheck=bamname
  }

  if (!ispdf){
    if (length(tocheck)<12){
      bottommargin=length(tocheck)+4*1.2
    }else{
      bottommargin=12+4*1.2
    }
  }else{
    bottommargin=14
  }




  par(mar=c(bottommargin,4,1,1),xpd=TRUE)
  suppressWarnings(boxplot(materialtoplot,at=addingfactor,col=colors,ylab=yl,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE,las=2))
  ats=c()
  for(i in 1:length(roiname)){
    window=addingfactor[(((i-1)*length(bamname))+1):(i*length(bamname))]
    currentvalue=(window[length(window)]-window[1])/2
    currentvalue=window[1]+currentvalue
    ats=c(ats,currentvalue)
  }


  if (!ispdf){
    if(!isgrouped){
      axis(1,at=ats,label=roiname,las=2)
      legend("bottom",inset=c(0,-((length(tocheck))+4)/20),legend=newnames,col=colors,cex=0.6,bg="transparent",pch=rep(19,length(materialtoplot)))          
    }else{
      axis(1,at=ats,label=newnames,las=2)
      legend("topright",legend=bamname,col=colors,cex=0.6,bg="transparent",pch=rep(19,length(materialtoplot)))
    }    
  }else{

      if(!isgrouped){
        axis(1,at=ats,label=roiname,las=2)
        plot.new()
        legend("topleft",legend=newnames,col=colors,cex=0.6,bg="transparent",pch=rep(19,length(materialtoplot)))
      }else{
        axis(1,at=ats,label=newnames,las=2)
        plot.new()
        legend("topleft",legend=bamname,col=colors,cex=0.6,bg="transparent",pch=rep(19,length(materialtoplot)))        
      }

  }
}




plot_profilesAndBox_corheat<-function(portionlist_boxes,bamname,cormethod) {
  mat=do.call(cbind,portionlist_boxes)
  #if log2, 0 will be -Inf. Correct them
  mat[is.infinite(mat) &mat<0 ]=0
  colnames(mat)=bamname
  #according to the input, use pearson or spearman
  correlation_total=cor(mat,method=cormethod)    
  trasp_cor=t(correlation_total)
  brk=c( seq( -1 , 1,0.01))
  my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

  par(mar=c(12,12,1,1),xpd=TRUE)
  image(0:nrow(trasp_cor), 0:ncol(trasp_cor),trasp_cor[,ncol(trasp_cor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
  axis( 2, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= rev(colnames( trasp_cor )), las= 2,cex.axis=0.8 )
  axis( 1, at=seq(0.5,ncol(trasp_cor)+0.5-1,1 ), labels= colnames( trasp_cor ), las= 2,cex.axis=0.8 )
  for (x in (nrow(correlation_total)-1+0.5):0.5  )
    for (y in 0.5: ((ncol(correlation_total)-1+0.5)   ))
      text(y,x, round(correlation_total[ncol(correlation_total)-x+0.5,y+0.5],2),col="blue",cex=0.8) 
}



plot_profilesAndBox_pcorheat<-function(portionlist_boxes,bamname,cormethod) {
  mat=do.call(cbind,portionlist_boxes)
  #colnames(mat)=bamselected[xleft:xright]
  colnames(mat)=bamname
  mat[is.infinite(mat) &mat<0 ]=0
  #accroding to the input, use pearson or spearman
  correlation_partial=pcor(mat,method=cormethod)$estimate
  colnames(correlation_partial)=rownames(correlation_partial)= colnames(mat)
  #correction for likely a bug in pcor function

  trasp_pcor=t(correlation_partial)
  brk=c( seq( -1 , 1,0.01))
  my_palette <- colorRampPalette(c("darkred","red", "white", "green","darkgreen"))(n = length(brk)-1 )  

  par(mar=c(12,12,1,1),xpd=TRUE)
  image(0:nrow(trasp_pcor), 0:ncol(trasp_pcor),trasp_pcor[,ncol(trasp_pcor):1,drop=FALSE],axes=FALSE,xaxt="n",yaxt="n", xlab = "", ylab = "",col=my_palette,breaks=brk)
  axis( 2, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= rev(colnames( trasp_pcor )), las= 2,cex.axis=0.8 )
  axis( 1, at=seq(0.5,ncol(trasp_pcor)+0.5-1,1 ), labels= colnames( trasp_pcor ), las= 2,cex.axis=0.8 )
  for (x in (nrow(correlation_partial)-1+0.5):0.5  )
    for (y in 0.5: ((ncol(correlation_partial)-1+0.5)   ))
      text(y,x, round(correlation_partial[ncol(correlation_partial)-x+0.5,y+0.5],2),col="blue",cex=0.8)  
}











###############################################################################
# dynamics
###############################################################################

plot_dynamics_profile<-function(profileList_to_plot,nbin,minval,maxval,ylabel,totalnames,colors,islog) {
  par(mar=c(5,5,1,4))
  plot(1, type="n", xlab="",xaxt="n", ylab=ylabel, xlim=c(1, nbin), ylim=c(minval, maxval))
  axis(1,at=c(1,nbin/2 +0.5,nbin),labels=c("TSS-30%","genes","TES+30%"))
  count=1
  for(i in 1:length(profileList_to_plot)){
    for(k in 1:length(profileList_to_plot[[i]])){
      if(islog){
        lines( log2(profileList_to_plot[[i]][[k]]),lwd=3,col= colors[count],pch=".",cex=2,xaxt="n")
      }else{
        lines(profileList_to_plot[[i]][[k]],lwd=3,col= colors[count],pch=".",cex=2,xaxt="n")
      }
      count=count+1
    }
  }
  legend("topright",legend=totalnames,col= colors,lty=rep(1,length(totalnames)),bg="transparent",lwd=3)  
}




plot_dynamics_box<-function(currentlist,islog,ylabel,colors,totalnames,main,ispdf=FALSE) {
  list2=list()
  #transform 2-order list in 1 order list
  for(i in 1:length(currentlist)){
    list2=c(list2,currentlist[[i]])
  }

  #HERE put the log2 of input$islogforDynamics
  #if Infinite values, those won't be drown
  if(islog){
    for(i in 1:length(list2)){
      list2[[i]]=log2(list2[[i]])
    }
    ylabel=paste("log2",ylabel)
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

  par(mar=c(5,4,3,1),xpd=TRUE)
  suppressWarnings(boxplot(list2,at=addingfactor,col=colors,ylab=ylabel,xaxt="n",notch=TRUE,varwidth=TRUE,outline=FALSE,main=main))
  axis(1,at=ats,label=names(currentlist))
  if (ispdf){
    plot.new()
    legend("topleft",legend=totalnames,col=colors,cex=1.2,bg="transparent",pch=rep(19,length(totalnames))) 
  }else{
    #legend("bottom",inset=c(0,-0.4),legend=totalnames,col=colors,cex=0.6,bg="transparent",pch=rep(19,length(totalnames))) 
  }
   
}





plot_dynamics_cumulative<-function(currentlist,islog,ylabel,outlayer_thresh,colors,main,ispdf=FALSE) {
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
    tmp=currentlist[[i]]
    if (islog){
      ylabel=paste("log2",ylabel)
      tmp=lapply(tmp,log2)
    }
    list2=c(list2,tmp)
  }
  for(i in 1:length(list2)){
    mincurrent=quantile(list2[[i]][!is.infinite(list2[[i]])],outlayer_thresh )
    maxcurrent=quantile(list2[[i]][!is.infinite(list2[[i]])],1-outlayer_thresh )
    minvals=c(minvals,mincurrent)
    maxvals=c(maxvals,maxcurrent)
  }

  maxvals=max(maxvals)
  minvals=min(minvals)

  lengths=sapply(list2,length)
  par(mar=c(5,4,3,1))
  plot(1, type="n", ylab="Ranked Genes", xlab=ylabel, ylim=c(0, quantile(1:max(lengths),1-(2*outlayer_thresh))), xlim=c(minvals, maxvals),main=main)
  for(i in 1:length(list2)){
    #remove a fraction of outlayers
    #celan -Inf
    vals=sort(list2[[i]])
    vals=vals[!is.infinite(vals)]
    vals=vals[vals<quantile(vals,1-outlayer_thresh) & vals>quantile(vals,outlayer_thresh)]
    lines(vals,1:length(vals),col=colors[i],lwd=3)
  }
  if (ispdf){
    plot.new()
    legend("topleft",legend=totalnames,col=colors,cex=1.2,bg="transparent",lty=1,lwd=3)
  }else{
    legend("bottom",inset=c(0,-0.4),legend=totalnames,col=colors,cex=0.6,bg="transparent",lty=1,lwd=3)
  }
}