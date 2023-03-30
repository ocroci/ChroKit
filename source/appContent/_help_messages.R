#messages for help buttons - ChroKit

#define common schemes (tip, warning...)
Field<-function(name,body){
	return(
		tags$p(
			tags$b(name),
			tags$br(),
			HTML(body)
		)
	)
}

Comment<-function(body,plural=F){
	if (plural==F){
		comm="COMMENT"
	}else{
		comm="COMMENTS"
	}	
	return(
		tags$i(
			comm,
			tags$br(),
			HTML(body)
		)
	)
}

Warning<-function(body,color="#b22222",plural=F){
	if (plural==F){
		warn="WARNING"
	}else{
		warn="WARNINGS"
	}
	return(
		list(
			HTML(paste("<font color=\"",color,"\">",warn,"</font>",sep="")),
			tags$br(),
			HTML(paste("<font color=\"#b22222\">",body,"</font>",sep=""))
		)
	)	
}

Tip<-function(body,plural=F){
	if(plural==F){
		tp="TIP"
	}else{
		tp="TIPS"
	}
	return(
		tags$i(
			tp,
			tags$br(),
			HTML(body)
		)
	)	
}




############################################################
# function for inline help (actionbutton in HTML block)
############################################################

htmlhelp<-function(text,id) {
	return(
		list(
			HTML(text),
			actionButton(id,"?",style="height:20px;width:14px;padding:0px;background-color: #b7b7b7;font-weight: bold")
		)
	)
}

htmlwarning<-function(text,id) {
	return(
		list(
			HTML(text),
			actionButton(id,"!",style="height:20px;width:14px;padding:0px;background-color: #ff0000;font-weight: bold")
		)
	)
}


############################################################
############################################################
############################################################
# Global_overview
############################################################
############################################################
############################################################
msg_global_overview<-list(
	title="Global ChroKit overview",
	text=list(
		HTML("ChroKit session starts by defining new ROIs. These are genomic coordinates that can be derived from BED/GTF files or from lists of genes.
			Optionally, enrichment files (bam or bigWig format) can be imported for further analyses. ROIs can be manipulated or filtered, and then a variety of plots can be generated.
			<br><img src='global_scheme_1.png' alt='Res' width='100%'><br><br>"),
		HTML("The sidebar menu contains the main ChroKit funcionalities: <br><img src='global_scheme_2.png' alt='Res' width='100%'><br><br>"),
		HTML("ROIs can be imported, viewed or downloaded: <br><img src='global_scheme_3.png' alt='Res' width='100%'><br><br>"),
		HTML("Different kind of visualizations can be used for analyse and explore NGS experiments: <br><img src='global_scheme_4.png' alt='Res' width='100%'><br><br>")
	)
)



############################################################
############################################################
############################################################
# ROI coordinate files
############################################################
############################################################
############################################################

############################################################
# import
############################################################


msg_coordinateFiles_chooseCoordinates<-list(
	title="Import a new ROI",
	text=list(
		HTML("Here you can import new ROIs (Regions of Interest). A ROI can be:<br>
		<li>A set of genomic coordinates</li>
		<li>A list of genes</li><br> 
		You can import or create new ROIs from files (bed/gtf/gff format), from lists of genes and from sequence pattern occurrences in the genome"),
		tags$br(),tags$br(),
		Comment("If an assembly is already loaded, the imported ROI will be automatically annotated with the genes nearest to its genomic ranges")	
	)
)


help_BED_fromfiles<-list(
	title="Import a ROI from file",
	text=list(HTML("Here you can import a ROI from a bed, gtf or gff file. It must be a tab-delimited text file with 3 (optionally 4) columns with the format :<br><br>
					Column 1: <i>chr</i>&emsp;Column 2: <i>start</i>&emsp;Column 3: <i>end</i>&emsp;Column 4: <i>strand</i><br><br>
					The following is an example of the content of the text file:<br><img src='File_ROI_format.png' alt='Res'><br><br>
					<i>'chr'</i> column should be in format 'chrN', where 'N' is the number of the chromosome (1,2,...X,Y,M)<br>
	 				<i>'start, end'</i> columns are numbers indicating the starting and ending point of each genomic range<br>
	 				<i>'strand'</i> column (optional) is the strand for each genomic range(can be either +, -, or * if strand is not determined)"),
		tags$br(),	
		tags$br(),
		Warning("Files must be reachable from the machine on which ChroKit is running.")

	)
)

help_BED_fromgenelist<-list(
	title="Import a genelist",
	text=list("Here you can import a genelist by pasting the gene ID or symbols, or loading a text file containing IDs or symbols.",
		tags$br(),
		tags$br(),
		Comment("When a gene list is loaded, the promoters, transcripts, TES of the genes associated to that gene lists are loaded in memory as new ROIs. 
				All annotated isoforms are loaded as well"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. The file containing the gene list must be a text file in which each row is a gene symbol/gene ID.<br> 
				<b>2</b>. It must be readable from the system in which ChroKit is running.<br>
				<b>3</b>. A genome assembly must be loaded to import a gene list (To load a genome assembly, go to the 'Assembly' section).",plural=T)
	)
)


help_BED_kindofID<-list(
	title="Kind of gene IDs",
	text=list("Select which kind of identifiers are in the gene list you are importing: Symbols, ENTREZ IDs, ...")
)
help_BED_maxtranscriptlen<-list(
	title="Maximum transcript length",
	text=list("Select the maximum length of the transcript allowed in the gene list you are going to import.
					Genes with transcripts length above that threshold will not be loaded",
		tags$br(),
		tags$br(),
		Warning("Association of enrichments to long transcripts may cause memory problems")

	)

)


help_BED_frompatterngenome<-list(
	title="Import a ROI from pattern occurrences",
	text=list(
		"Import a new ROI from a sequence pattern inside the genome. Coordinates of the created ROI will be the positions of occurrences
			 of the pattern in the genome",
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. A genome assembly must be loaded to search for patterns. To load a genome assembly, use 'Assembly' section.<br>
				<b>2</b>. A BSgenome database must be also loaded: this database contains the information about the sequence 
				in a specific genomic region in the genome assembly in use. If this database is not installed in the 
				system, a new button will appear to enable the User to download it from bioConductor. 
				Make sure your internet connection is working",plural=T)			 
	)
)


help_BED_frompatterngenome_warning<-list(
	title="Warning",
	text=list(
		HTML(paste("<font color=\"#b22222\">This could require a substantial amount of computational time</font>",sep=""))
	 
	)
)


help_BED_IUPACpattern<-list(
	title="IUPAC nomenclature",
	text=list(
				HTML("Select the pattern you want to extract. Use the IUPAC nomenclature:<br>
				<b>A,C,G,T,U</b>: known bases<br>
				<b>R</b>: A or G<br>
				<b>Y</b>: C or T<br>
				<b>S</b>: G or C<br>
				<b>W</b>: A or T<br>
				<b>K</b>: G or T<br>
				<b>M</b>: A or C<br>
				<b>B</b>: C or G or T<br>
				<b>D</b>: A or G or T<br>
				<b>H</b>: A or C or T<br>
				<b>V</b>: A or C or G<br>
				<b>N</b>: any base<br>
				<b>. or -</b>: gap<br>")
	)

)

help_BED_headeroption<-list(
	title="Is the header present?",
	text=list("If the file has a header (the column names are present), check this option")
)


help_BED_linesskip<-list(
	title="Skip lines",
	text=list("Set the number of lines that should be ignored at the beginning of the file. If you have extra-lines at the beginning of your file, you can skip them by
		 			inserting their number here")
)

help_BED_filepreview<-list(
	title="View the file before importing",
	text=list(
		"Here you can interactively explore the opened file. If you are satisfied and you want to import it as new ROI, click 'Confirm and import as ROI' button"
	)
)


############################################################
# view ROI
############################################################
msg_quickviewROIs<-list(
	title="Quick look at imported ROIs",
	text="Select one or more ROIs to view the number of genomic ranges or the ranges width distribution"
)


help_BED_viewoptions<-list(
	title="Which info to display",
	text=HTML("Select which information you want to display of the selected ROIs.<br>
		<b>Distribution of genomic ranges width</b>: This shows the distribution of the log2 of the width of the ranges<br>
		<b>Number of ranges</b>: This shows a barplot of the number of genomic ranges for each of the selected ROIs")
)



############################################################
# delete ROI
############################################################

#Delete ROIs
msg_deleteRois_deleteRois<-list(
	title="Remove ROIs",
	text=list(
		"Displays all the ROIs in the current session. Click to select the ROIs to be deleted.",
		tags$br(),
		tags$br(),
		Tip("Remove unnecessary ROIs if your session becomes too heavy (too many ROIs with too many enrichments associated to them). This 
			will save memory and make the session files smaller"),
		tags$br(),
		tags$br(),
		Warning("When you remove a ROI, all the enrichments associated to it will be removed, as well")
	)
)



############################################################
# get ROI
############################################################

msg_getRois_BOX<-list(
	title="Export a ROI to a file",
	text=list(
		"Download a ROI and its associated information."
	)
)

help_BED_getroi_eachGR<-list(
	title="Information for each genomic range of the ROI",
	text=list("Displays the genomic coordinates and/or nearest gene IDs and/or associated enrichments 
				for each genomic range of the selected ROI.",
				tags$br(),
				tags$br(),
				Warning("A genome assembly must be loaded to get the annotations. To load a genome assembly, 
				use 'Assembly' section. Once an assembly has been loaded, every genomic range of a 
				ROI is annotated to the nearest gene")
		)
)

help_BED_getroi_genomicWindow<-list(
	title="Annotated genes inside genomic window",
	text=list("It retrieves the list of genes (IDs and/or symbols) 
					found at a user-defined distance from the midpoint of all genomic ranges of the ROI.",

		tags$br(),
		tags$br(),
		tags$br(),
		Warning("A genome assembly must be loaded to get the annotations. To load a genome assembly, 
				use 'Assembly' section")
		)
)




help_BED_getroi_windowvalue<-list(
	title="Number of bp for genomic window",
	text=list("The number of base pair upstream and downstream the midpoint of each genomic range of the selected ROI.
			All the genes inside this window will be displayed.")
)




############################################################
############################################################
# enrichment files
############################################################
############################################################

#import enrichment file
msg_enrichmentFiles_importEnrichment<-list(
	title="Import an enrichment file (BAM/bigWIG)",
	text=list(
		"This section is used to associate enrichment files to a ChroKit session.",
		tags$br(),
		tags$br(),
		Field("Choose one or more files from filesystem","Select the enrichment file(s) (BAM or bigWIG) by exploring the filesystem"),
		Field("Manually type the path of a file","Type here the path to the enrichment file (BAM or bigWIG)"),
		tags$br(),
		Comment("Formats supported: BAM and bigWIG. When an enrichment file is selected and imported, the program keeps in memory the link (path)
				 to that file"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. The file must be reachable from the machine in which ChroKit is running.<br>

				<b>2</b>. If a BAM file is chosen, a BAM index (.bai) file should be present in the same directory 
				and must have the same name of the bam, plus the '.bai' extension. 
				For example, if the file is 'enrichment.bam', the corresponding index should be 'enrichment.bam.bai'. 
				bigWIG files are not currently supported by Windows operating systems",plural=T),
		tags$br(),
		tags$br(),
		Tip("To create a BAM index, use samtools (outside ChroKit) with the following syntax:<br><br> 
			samtools index &lt;FileName&gt;.bam")				

	)
)
#delete enrichment file
msg_enrichmentFiles_deleteEnrichment<-list(
	title="Delete enrichment files",
	text=list(
		"Displays all the enrichment files loaded in the current session. 
		Click to select the links to the enrichment files (file paths) to be deleted.",
		tags$br(),
		tags$br(),

		Comment("Note that only the links to the enrichment files (file paths) will be deleted. 
				Original files won't be deleted from the system"),
		tags$br(),
		tags$br(),

		Warning("The enrichments already associated to ROIs won't be deleted"),
		tags$br(),
		tags$br(),

		Tip("This feature is useful when the location or the name of the enrichment files has changed in the system")				
	)
)
#rename enrichment file
msg_enrichmentFiles_renameEnrichment<-list(
	title="Rename enrichment files",
	text=list("This will change the name of an enrichment file in all the menus and plots",
		tags$br(),
		tags$br(),	
		Warning("This only changes the name of the enrichment file within ChroKit, while the original file in the system will not be renamed")
	)
)




# #ROI selection
# msg_getRois_roiSelection<-list(
# 	title="Get all the info associated to a ROI",
# 	text=list(
# 		"Select the ROI for which you want to explore or download the information. 
# 		Once a ROI is selected, the following information can be retrieved",
# 		tags$br(),
# 		tags$br(),
# 		tags$h4("Features for each genomic range:"),
# 		tags$br(),
# 		"Displays the genomic coordinates and/or nearest gene IDs and/or associated enrichments for each genomic range of the selected ROI",
# 		tags$br(),
# 		tags$br(),
# 		Field("View Ranges","Genomic coordinates of the ROI and the strand"),
# 		Field("Select annotations to show","Here you can include IDs and/or symbols 
# 											of the annotated genes"),
# 		tags$br(),
# 		Warning("A genome assembly must be loaded to get the annotations. To load a genome assembly, 
# 				use 'Assembly' section. Once an assembly has been loaded, every genomic range of a 
# 				ROI is annotated to the nearest gene"),
# 		tags$br(),
# 		tags$br(),
# 		Field("Enrichments","These are the signal enrichments associated to the ROI. It shows the sum 
# 				of the pileup of the reads for the selected enrichment within each genomic range of the ROI"),
# 		tags$br(),
# 		tags$h4("Gene list inside genomic window:"),
# 		tags$br(),
# 		"Displays the list of gene IDs and/or symbols inside a genomic window from the midpoint of each genomic range of the selected ROI",
# 		tags$br(),
# 		tags$br(),
# 		Field("Get annotated genes within a genomic window","It retrieves the list of genes (IDs and/or symbols) 
# 					found at a user-defined distance from the genomic ranges of the ROI. User defined 
# 					distance is set in '<b>Select the genomic window (bp)</b>'"),
# 		Field("Select the genomic window (bp)","Select the size of the genomic window flanking the genomic ranges 
# 					of a ROI. The width is measured from  the center of the genomic ranges"),		
# 		tags$br(),
# 		tags$h4("Edit notes of the ROI:"),
# 		tags$br(),
# 		"Displays the notes of the selected ROI. These notes can be edited and saved or downloaded as text file",
# 		tags$br()

# 	)
# )













#Rename ROIs
msg_deleteRois_renameRois<-list(
	title="Change the name of a ROI",
	text=list(
		"Select a ROI and type the new name",
		tags$br(),
		tags$br(),
		Warning("You cannot use 'promoters', 'transcripts', 'TES' as names. Also names cannot begin with 'promoters_genelist_', 'transcripts_genelist_' or 'TES_genelist_'")
	)
)
#Reorder ROIs
msg_deleteRois_reorderRois<-list(
	title="Change the order of ROIs",
	text=list(
		"For each ROI in the list, select the number corresponding to the new position 
		in the ranking, and press the 'Reorder!' button",
		tags$br(),
		tags$br(),
		Comment("This is the order used for some displayed items, as for example position-based heatmaps or menus")
	)
)



############################################################
############################################################
# Assembly
############################################################
############################################################


help_txdb_windowupdownstream<-list(
	title="base pair for defining TSS and TES regions",
	text=list(
		HTML("By indicating the number of base pairs upstream and downstream from TSS and TES, 
									you can choose the size of the genomic window encompassing the TSS and TES"),
		tags$br(),
		tags$br(),
		Tip("By choosing 2000 upstream and 1000 downstream  of  a TSS,  you  will define promotorial regions of 3000bp, 
					starting at  -2000 and ending at +1000 from the annotated TSS")
	)
)

#Extract annotated elements from database
msg_databases_extractAnnotatedElements<-list(
	title="Select an available assembly",
	text=list(

		Field("Choose assembly to use","Choose a genome assembly to extract annotated genomic elements (promoters, transcripts, TES)"),
		tags$br(),

		Comment("1. The selected assembly will be used to automatically annotate all the ROIs that will be imported later.<br>
				2. The selected assembly will be used to obtain a set of promoters, transcripts and TES from a custom list of genes ('ROIs' section). <br>
				3. ROIs already imported into the program will be automatically re-annotated according to the newly selected genome assembly 
				(except those derived from a previous assembly, for example the promoters of a list of genes obtained from another assembly).<br>
				4. Additional genome assemblies can be downloaded from bioconductor using the '... or download an assembly' box 
				on the right. ",plural=T)
		
	)	
)

#Download databases
msg_databases_downloadDatabases<-list(
	title="Download annotation packages from bioConductor",
	text=list(
		"Download the packages of the desired genome assembly (if not already present in the system) 
			directly from bioConductor",
		tags$br(),
		tags$br(),
		Comment("This will download the databases associated to the selected genome assembly (for example, 
				for mm10 genome assembly it will download TxDb.Mmusculus.UCSC.mm10.knownGene and org.Mm.eg.db databases)"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. Be sure your internet connection is working! <br>
				<b>2</b>. This may require some time",plural=T)		
	
	)
)



############################################################
############################################################
# ROI preparation
############################################################
############################################################


############################################################
# main menu
############################################################
msg_ROImanipulation<-list(
	title="Prepare ROIs for analyses and visualization",
	text=list("A set of tools for preparing ROIs for data visualization")
)

msg_prepare_ultraeasy<-list(
	title="Prepare a ROI for visualization",
	text=list("The easiest way to prepare a ROI for analyses. This step simply associates one or more enrichment files to a ROI.
			This is the minimum requirement for basic quantitative analyses.",
		tags$br(),
		tags$br(),
		Warning("You need to import enrichment files for this operation")
	)
)

msg_prepare_forheat<-list(
	title="Prepare a ROI for heatmaps",
	text=list("Create a new modified ROI from a pre-existent one, that can be used for heatmaps visualization and analyses.
			If enrichment files are not loaded, the resulting ROI can be used only for position-based heatmaps."
	)
)

msg_prepare_formetagene<-list(
	title="Prepare a genelist for metagene profiles",
	text=list("Associates enrichments to imported genelists for metagene profile analyses.",
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. You need to import enrichment files<br>
				<b>2</b>. You need to import the correct genome assembly<br>
				<b>3</b>. You need to import genelists<br>",plural=T)	
	)
)

msg_prepare_manual<-list(
	title="In-depth ROI editing",
	text=list("Manually edit ROIs. This is a 'swiss army knife' for modifying genomic ranges, associating enrichments, filtering genomic intervals in detail."
	)
)



############################################################
# prepare ultraeasy
############################################################
help_ultraeasy_enrichmentAssoc<-list(
	title="Associate enrichments",
	text=list("Associate one or more enrichments to the ROI. This is required for quantitative analyses on the reads enrichment of NGS experiments",
		tags$br(),
		tags$br(),
		Comment("This step is required to perform quantitative analyses"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>.You need to import enrichment files for this operation<br>
				<b>2</b>. bigWIG files are not supported in Windows operating systems",plural=T)
	)
)



############################################################
# prepare for heatmaps
############################################################
help_prepheat_subsample<-list(
	title="Subsample ROI",
	text=list("This value determines how many of the genomic ranges of the selected ROI will be randomly subsampled.
			For instance, for a ROI of 10000 ranges, by selecting 30% you will generate a random subset of 3000 genomic regions.")

)

help_prepheat_summit<-list(
	title="Center the genomic ranges on summits",
	text=list("If 'Yes', set the midpoint of each genomic range on the summit (i.e. position with the highest amount of reads) of a specific enrichment.
			The ranges shown in heatmaps will be centered on this summits.",
		tags$br(),
		tags$br(),
		Warning("You need to import enrichment files for this operation")
		)
)

help_prepheat_resize<-list(
	title="Resize genomic ranges",
	text=list("Resize each genomic range to a fixed window, by selecting the number of bp upstream and downstream the center of genomic ranges.")
)


help_prepheat_enrichmentAssoc<-list(
	title="Associate enrichments",
	text=list("Associate one or more enrichments to the ROI. This is required for quantitative analyses on the reads enrichment of NGS experiments",
		tags$br(),
		tags$br(),
		Comment("This step is required to generate enrichment-based heatmaps and to perform quantitative analyses"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>.You need to import enrichment files for this operation<br>
				<b>2</b>. bigWIG files are not supported in Windows operating systems",plural=T)	
	)

)

############################################################
# prepare for genelists
############################################################


help_prepgenelist_ROIassociated<-list(
	title="ROIs which constitute the selected genelists",
	text=list("The list of ROIs (promoters, transcripts and TES) that, together, form the selected genelists")
)

help_prepgenelist_enrichmentAssoc<-list(
	title="Associate enrichments",
	text=list("Associate one or more enrichments to genelists. 
				Selected enrichments will be associated to both the promoters, transcripts and TES that constitute the genelists.",
		tags$br(),
		tags$br(),
		Warning("<b>1</b>.You need to import enrichment files for this operation<br>
				<b>2</b>. bigWIG files are not supported in Windows operating systems",plural=T)	
	)
)



############################################################
# manual ROI editing
############################################################

#overlaps
help_roimanual_overlaps<-list(
	title="Overlap two or more ROIs",
	text=list(HTML("Here you can set the rules for defining overlapping or non-overlapping genomic ranges and generate new ROIs.<br><br>
		<img src='overlaps.png' alt='Res' width=60% height=60%>")
	)
)


help_roimanual_overlaps_selectref<-list(
	title="Select reference ROI",
	text=list(HTML("The ROI(s) that will be used as reference. If more than one ROI is selected, 
				a combination of them will be used as starting point for the overlaps.<br><br>
				<img src='overlaps_referenceroi.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),
		Comment("1. If multiple ROIs are selected, the program will calculate 
				the union or intersection of their genomic ranges and will produce a single reference.<br>
				2. The new ROI will keep the enrichments associated to the reference ROI only if a single reference ROI is selected",plural=T)
	)
)

help_roimanual_overlaps_intersectionref<-list(
	title="Intersection of reference ROIs",
	text=HTML("the aggregated reference ROI will be constituted by 
					the genomic ranges common to all the ROIs selected. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated reference ROI will be: 
					A&cap;B&cap;C")
)

help_roimanual_overlaps_unionref<-list(
	title="Union of reference ROIs",
	text=HTML("The aggregated reference ROI will be the union of all the genomic 
					ranges contained in the selected ROIs. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated reference ROI will be:
					A&cup;B&cup;C")
)


help_roimanual_overlaps_selectoverlapwith<-list(
	title="Select overlapping contrast ROI",
	text=list(HTML("Select ROIs that have to overlap with the reference ROI<br><br>
		<img src='overlaps_overlappingroi.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),	
		Comment("Only genomic ranges of the reference ROI that overlap with those selected here will be kept")
	)
)

help_roimanual_overlaps_intersectionoverlapwith<-list(
	title="Intersection of overlapping contrast ROI",
	text=HTML("The aggregated contrast ROI will be the genomic ranges common to all the ROIs selected. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be: A&cap;B&cap;C")
)

help_roimanual_overlaps_unionoverlapwith<-list(
	title="Union of overlapping contrast ROI",
	text=HTML("The aggregated contrast ROI will be the union of all the genomic 
					ranges contained in the selected ROIs. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be: A&cup;B&cup;C")
)

help_roimanual_overlaps_selectnotoverlapwith<-list(
	title="Select not overlapping contrast ROI",
	text=list(HTML("Select ROIs that do not have to overlap with the reference ROI<br><br>
		<img src='overlaps_notoverlappingroi.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),	
		Comment("The resulting ROI will contain only the genomic ranges of the reference ROI that 
					do not overlap with those selected here")
	)
)

help_roimanual_overlaps_intersectionnotoverlapwith<-list(
	title="Intersection of overlapping contrast ROI",
	text=HTML("The aggregated contrast ROI will be the genomic ranges common to all the ROIs selected. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be: A&cap;B&cap;C")
)

help_roimanual_overlaps_unionnotoverlapwith<-list(
	title="Union of overlapping contrast ROI",
	text=HTML("The aggregated contrast ROI will be the union of all the genomic 
					ranges contained in the selected ROIs. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be: A&cup;B&cup;C")
)


help_roimanual_overlaps_minimumbp<-list(
	title="Minimum number of bp to consider for overlaps",
	text="Two genomic ranges will be deemed overlapping if 
						the number of shared bases is equal or greater than the value indicated here"

)

help_roimanual_overlaps_strandspecific<-list(
	title="Strand-specific overlaps",
	text="If selected, the overlaps will be determined considering the strand 
						information: only genomic regions with concordant strand information (i.e. same strand) 
						will be tested. If strand information is not available, Chrokit will ignore this option"
)




#resize
help_roimanual_resize<-list(
	title="Resize ROI boundaries",
	text=list(HTML("Change the width of the genomic ranges present in a ROI.<br><br>
		<img src='resize.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),
		Warning("The new ROI won't keep the enrichments associated to the old ROI")

	)
)

help_roimanual_pattern<-list(
	title="Extract sequence pattern from a ROI",
	text=list(HTML("You can extract user-provided sequence patterns from a ROI.
			The output is a new ROI composed of genomic ranges that contain the user-defined pattern. 
			The new ROI obtained after the pattern search will be strand-specific, based on the pattern found. 
			The genomic ranges will be centred on the pattern.<br><br>
			<img src='pattern.png' alt='Res' width=60% height=60%>"),
			tags$br(),
			tags$br(),
			Warning("<b>1</b>. A genome assembly must be loaded to search for patterns. To load a genome assembly, use 'Assembly' section.<br>
				<b>2</b>. A BSgenome database must be also loaded: this database contains the information about the sequence 
				in a specific genomic region in the genome assembly in use. If this database is not installed in the 
				system, a new button will appear to enable the User to download it from bioConductor. 
				Make sure your internet connection is working.<br>
				<b>3</b>. The new ROI won't keep the enrichments associated to the old ROI.",plural=T)


	)
)

help_roimanual_pattern_IUPAC<-list(
	title="Select pattern (IUPAC nomenclature)",
	text=HTML("Select the pattern you want to extract. Use the IUPAC nomenclature:<br>
				<b>A,C,G,T,U</b>: known bases<br>
				<b>R</b>: A or G<br>
				<b>Y</b>: C or T<br>
				<b>S</b>: G or C<br>
				<b>W</b>: A or T<br>
				<b>K</b>: G or T<br>
				<b>M</b>: A or C<br>
				<b>B</b>: C or G or T<br>
				<b>D</b>: A or G or T<br>
				<b>H</b>: A or C or T<br>
				<b>V</b>: A or C or G<br>
				<b>N</b>: any base<br>
				<b>. or -</b>: gap<br>"
	)
)

help_roimanual_pattern_bothstrands<-list(
	title="Search on both strands",
	text="The sequence pattern will be searched on both strands for all genomic ranges of the ROI."
)

help_roimanual_pattern_strandspecific<-list(
	title="Strand-specific pattern search",
	text=list("The sequence pattern will be searched only in the strand defined for each genomic range of the ROI.",
		tags$br(),
		tags$br(),	
		Comment("For example, if one of the genomic ranges in the ROI is defined as '-' strand, the pattern will be searched only in the negative strand for this range")	
	)
)


#summit
help_roimanual_summit<-list(
	title="Center a ROI on the summit of an enrichment",
	text=list(HTML("For each genomic range of a ROI, it calculates the summit of a signal enrichment (i.e. the coordinate with the maximum value of pileup reads). 
		Then, each range of the ROI will be centered on that summit<br><br>
		<img src='summit_detection.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. At least one enrichment must be associated to the ROI<br>
				<b>2</b>. The new ROI won't keep the enrichments associated to the old ROI",plural=T)		
	)
)


help_roimanual_summit_enrichmenttouse<-list(
	title="Enrichment to use for summit detection",
	text="This enrichment is used to calculate the point of maximum reads pileup (summit)"
)


#association enrichments
help_roimanual_enrichmentAssoc<-list(
	title="Associate or remove enrichments to ROI(s)",
	text=list(HTML("Associate one or more enrichments (from BAM/bigWIG files imported) to one or more ROIs"),
		tags$br(),
		tags$br(),
		Comment(		
		"This operation will calculate the pileup of the reads for each base of each range of each ROI selected.
 		The information of the reads come from the enrichment files (either BAM or bigWIG) selected.
 		This is a computationally intensive task, it may require time (from several seconds to few minutes, depending on the size of the ROIs, 
 		the size of the enrichment files, the number of ROIs and the computational power). <br><br>
 		<img src='enrichment_association.png' alt='Res' width=80% height=80%>
 		<br><br>Signals will be normalized only for BAM files: for bigWIG files, the normalization coefficient will be 1."),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. You need to import enrichment files for the association<br>
				<b>2</b>. bigWIG files are not supported in Windows operating systems",plural=T)	
	)
)


help_roimanual_enrichmentAssoc_enrichments<-list(
	title="Select enrichment(s) to associate to selected ROI(s)",
	text=list("Select one or more enrichments (BAM/bigWIG) 
				to associate. This operation will calculate the pileup of the reads 
				from the enrichment files (either BAM or bigWIG) selected for each base of each range of the ROI(s)",
		tags$br(),
		tags$br(),
		Tip("To import enrichment files go to 'Enrichment files' in the 'Import data' section"),
		tags$br(),
		Comment("The normalization is available only for BAM files. bigWIG files will not be normalized.")
	)
)

help_roimanual_enrichmentAssoc_NOnorm<-list(
	title="Do not normalize",
	text="The number of reads will not be normalized (obtaining raw pileup number for each base pair of each range)."
)

help_roimanual_enrichmentAssoc_librarysize<-list(
	title="Library-size normalization",
	text="The number of reads will be divided by the size of their library and multiplied by 1 million."
)

help_roimanual_enrichmentAssoc_customnorm<-list(
	title="Use custom normalizer",
	text="The number of reads will be divided by the library size of another enrichment file present in the list and multiplied by 1 million."
)

help_roimanual_enrichmentAssoc_spikein<-list(
	title="Spike-in normalization",
	text=list("The number of reads will be normalized by the library size of another enrichment file present in the list (spike-in);
						this number will be divided by the number of reads of a control and multiplied by the number of 
						reads of the spike-in of the control. Finally, the number will be multiplied by 1 million.",
		tags$br(),
		tags$br(),
		Comment("The normalizing factor in spike-in normalization will be:<br><br> (&ltspike-in reads in control&gt*1000000)/(&ltspike-in reads&gt*&ltcontrol reads&gt)<br><br>
					This normalization is particularly useful for ChIP-Seq experiments.")
	)
)

help_roimanual_enrichmentAssoc_numbercores<-list(
	title="Set CPU cores",
	text=list("Choose the number of CPU cores to use for enrichment association",
		tags$br(),
		tags$br(),
		Comment("The more cores you use, the faster the association will be, but more memory will be required.
			Leave '1' unless you have a lot of RAM or if you are unsure.")
	)
)


#ROI sample
help_roimanual_subsample<-list(
	title="Random subset a ROI",
	text=list(HTML("Create a new ROI from a random sample of genomic rages of another ROI<br><br>
		<img src='random_sample.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),
		Tip("Subsetting is useful to reduce memory usage. This is strongly suggested if you are 
			dealing with large ROIs (with above 50/60000 genomic ranges), since it reduces processing 
			time and memory usage during some downstream analysis (such as the association of enrichments), 
			but preserve a statistically significant sample representative of the original ROI")				

	)
)

#filter for width
help_roimanual_filterwidth<-list(
	title="Subset a ROI based on the width of its genomic ranges",
	text=list(HTML("ROIs can be subsetted according to the width of their genomic ranges. Only genomic ranges above/below certain width thresholds
				will be considered.<br><br>
				<img src='filterwidth.png' alt='Res' width=60% height=60%>")

	)
)

help_roimanual_filterwidth_min<-list(
	title="Min width threshold",
	text="Genomic ranges below this width will be discarded"
)

help_roimanual_filterwidth_max<-list(
	title="Max width threshold",
	text="Genomic ranges above this width will be discarded"
)


#filter for enrichment
help_roimanual_filterenrich<-list(
	title="Subset a ROI based on a specific enrichment",
	text=list(HTML("ROIs can be subsetted by their signal enrichment.
		Only genomic ranges above/below certain enrichment thresholds
				will be considered.<br><br>
				<img src='filterenrichment.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),		
		Warning("At least one enrichment must be associated to the ROI")

	)
)

help_roimanual_filterenrich_min<-list(
	title="Min enrichment threshold",
	text="Genomic ranges below this enrichment value will be discarded"
)

help_roimanual_filterenrich_max<-list(
	title="Max enrichment threshold",
	text="Genomic ranges above this enrichment value will be discarded"
)


help_roimanual_renameroi<-list(
	title="Rename a ROI",
	text=list("Change the name of a ROI. The new name will appear in all the menus and plots."

	)
)

help_roimanual_editnotes<-list(
	title="Edit ROI notes",
	text=list("View or edit a short report for a specific ROI"

	)
)

help_roimanual_renameenrich<-list(
	title="Rename associated enrichments",
	text=list("Rename enrichments associated to a ROI",
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. The name of the enrichment will change in all the other ROIs to which it had been associated<br>
			<b>2</b>. At least one enrichment must be associated to the ROI",plural=T)

	)
)

help_roimanual_reorderroi<-list(
	title="Reorder ROi list",
	text=list("Change the order in which ROIs appear in menus and plots.",
		tags$br(),
		tags$br(),
		Comment("This is useful when you want to change the order of the ROIs in which they appear in plots",)		

	)
)

help_roimanual_reorderroi_index<-list(
	title="Change the order in which ROI appear",
	text="For each ROI, set the number corresponding to the new position in the ranking."
)

help_roimanual_reorderenrich<-list(
	title="Reorder associated enrichments",
	text=list("Change the order in which the enrichments appear in a ROI.",
		tags$br(),
		tags$br(),
		Comment("1. This is useful when you want to change the order of the enrichments in which they appear in plots. <br>
				2. The order of the enrichments will be changed only in the selected ROI",plural=T),
		tags$br(),
		tags$br(),		
		Warning("At least one enrichment must be associated to the ROI")

	)
)

help_roimanual_reorderenrich_index <-list(
	title="Change the order in which the enrichments appear in the ROI",
	text="For each of the enrichments associated to the selected ROI, 
			set the number corresponding to the new position in the ranking."
)








############################################################
############################################################
# Genomics
############################################################
############################################################


#### Single evaluation

# (parameters)
msg_singleEvaluation_parameters<-list(
	title="Parameters for single ROI evaluation",
	text=list(
		"Set parameters to view the genomic distribution of the genomic ranges of a ROI (and its enrichments)",
		tags$br(),
		tags$br(),
		Warning("A genome assembly must be loaded. To load a genome assembly, go to the 'Assembly' section")
	)
)
#Distribution
msg_singleEvaluation_distribution<-list(
	title="Distribution in annotated genomic elements",
	text="Displays the location of genomic ranges of the selected ROI in the annotated regions 
	(promoters, genebodies, intergenic regions). Results can be shown both as barplots or piecharts"
)
#Width distribution
msg_singleEvaluation_widthDistribution<-list(
	title="Width distribution",
	text=list(
		"Displays the distribution of the width of genomic ranges of the selected ROI stratified 
		by class of annotated regions (promoters, genebodies, intergenic regions)",
		tags$br(),
		tags$br(),
		Warning("To avoid outliers, only 0-95% of the distribution is shown")		
	)
)
#Enrichment boxplot
msg_singleEvaluation_enrichmentBoxplot<-list(
	title="Enrichment in annotated genomic elements",
	text=list(
		"Displays a boxplot 
		of the reads (or the read density) stratified by the annotated regions (promoters, genebodies, intergenic regions)",
		tags$br(),
		tags$br(),
		Warning("The enrichment file must be associated to the ROI. To associate enrichments to a ROI, 
			go to 'ROI preparation' and select 'Prepare ROI basic'")		
	)
)
#Peak average profile
msg_singleEvaluation_peakProfile<-list(
	title="Profile of the enrichment in annotated genomic elements",
	text=list(
		"Displays the profile of an enrichment in selected ROI stratified by the annotated regions 
		(promoters, genebodies, intergenic regions), centred at the midpoint of the genomic ranges",
		tags$br(),
		tags$br(),
		Warning("The enrichment file must be associated to the selected ROI. To associate enrichments to a ROI, 
			go to 'ROI preparation' and select 'Prepare ROI basic'")		
	)
)



help_singleEvaluation_normalizationtotalread<-list(
	title="Show total reads",
	text="Show the number of reads. (Depending on the enrichment association, these can be normalized for library size or not)"
)
help_singleEvaluation_normalizationreaddensity<-list(
	title="Show read density",
	text="Show the number of reads, 
		 			normalized by the width of the genomic ranges."
)







#### Pairwise overlaps

#Parameters
msg_pairwiseOverlaps_parameters<-list(
	title="Parameters for pairwise overlaps",
	text=list("Analyze the pairwise overlap between two ROIs (and their associated enrichments)",
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. This analysis can be performed only if two ROIs are selected<br>
				<b>2</b>. Enrichments can be analyzed only if they are associated 
				to ROI-1, ROI-2 or both. If enrichments are not associated to either ROI-1 or ROI-2, only 
				the fraction of overlapping/not overlapping genomic ranges will be displayed.",plural=T)								
	)
)


help_pairwiseOverlaps_parameters_minbpoverlap<-list(
	title="Minimum number of bp for overlap",
	text="Set the minimum number of bases 
					used to define the overlap between genomic ranges of the two ROIs."
)


help_pairwiseOverlaps_parameters_ROIuniverse<-list(
	title="ROI to be used as background",
	text="The background ROI is needed for the computation of the statistics of the overlap.
			This background will be the union of the selected ROI(s), plus ROI-1 and ROI-2.
			The background ROI will be the universe of the hypergeometric test."

)

help_pairwiseOverlaps_parameters_enrich1<-list(
	title="Choose enrichment of ROI-1",
	text="Choose one of the enrichments associated to the first ROI to 
				evaluate the signal of reads in the subset of genomic ranges of ROI-1 alone 
				or those overlapping with ROI-2"
)

help_pairwiseOverlaps_parameters_enrich2<-list(
	title="Choose enrichment of ROI-2",
	text="Choose one of the enrichments associated to the second ROI to evaluate 
				the signal of reads in the subset of genomic ranges of ROI-2 alone or those overlapping 
				with ROI-1"
)



#use them for read density
# help_singleEvaluation_normalizationtotalread
# help_singleEvaluation_normalizationreaddensity



#Overlap
msg_pairwiseOverlaps_overlap<-list(
	title="Displays the overlap of the two ROIs",
	text=list(
		"The overlap of the two ROIs can be displayed as a barplot or with a Venn diagram",
		tags$br(),
		tags$br(),
		Field("Barplot","It displays the fraction of the genomic ranges of the ROI-1 that overlap 
			with the ROI-2 and the fraction of ranges exclusively present only in ROI-1 (and vice versa)"),
		Field("Venn","It shows the overlap between the ROI-1 and ROI-2. The common area in the diagram 
				represents the number of regions that are overlapping between the two ROIs"),
		tags$br(),
		Comment("In the Venn diagram, the number of common genomic ranges is an approximation: those common 
				regions are all the contiguous genomic intervals in which ROI-1 and ROI-2 overlap")


			# tags$br(),
			# tags$br(),
			# "ROI1: -----   ----     ----     ---------     ----------",
			# tags$br(),
			# "ROI2: ---      ----    ----------------      ----    ----",
			# tags$br(),
			# "In this example, the number of overlapping regions will be 4."

	)
)





msg_pairwiseOverlaps_box<-list(
	title="Boxplot",
	text=list(HTML("It shows the values (number of reads) of the selected enrichments (i.e. enrichment-1 and enrichment-1) within:<br>
				1) regions of ROI-1 not overlapping with ROI-2 (ROI-1 alone)<br>
				2) regions of ROI-2 not overlapping with ROI-1 (ROI-2 alone)<br>
				3) common ROI-1 and ROI-2 regions (ROI-1 overlapping and ROI-2 overlapping) <br>
				<img src='scheme_enrichment_pairwise.png' alt='Res' width=100% >")
	)

)

msg_pairwiseOverlaps_scatter<-list(
	title="Scatterplot",
	text=list("If enrichment-1 and enrichment-2 are associated to both ROI-1 and ROI-2, it shows the 
				correlation between the two enrichments. Overlapping and exclusive (ROI-1 alone, ROI-2 alone) genomic 
				ranges are highlighted.",
			tags$br(),
			tags$br(),
			Warning("This works only if both enrichment-1 and enrichment-2 are associated with both ROI-1 and ROI-2")
	)

)

msg_pairwiseOverlaps_calibration<-list(
	title="Calibration",
	text=list("It shows the fraction of the overlapping regions between ROI-1 and ROI-2 at increasing values 
			of enrichments thresholds. The threshold is applied either to ROI-1 or ROI-2 or to both ROIs together (combined).",
			tags$br(),
			tags$br(),
			Comment("This plot recalculates the overlap between ROI-1 and ROI-2 if we consider only those genomic ranges 
				with an enrichment above an increasing threshold. 'combined' means that the enrichment threshold 
				is applied to both ROI-1 and ROI-2 simultaneously; otherwise, the threshold is applied only to 
				one of the two ROIs, thus considering all the genomic ranges of the other ROI for the analysis.<br><br>
				With this plot, it is possible to assess which is the % of overlap between A (ROI-1) and B (ROI-2) considering 
				only the top 20% enriched genomic ranges in A, or considering the top 20% enriched regions in A and B. 
				If, in a ChIP-seq experiment, the fraction of overlap increases when considering only 
				the top enriched ranges, it may suggest that common genomic ranges are high affinity sites for 
				both ChIP-seq A and ChIP-seq B.")					
	)

)







#### Digital heatmap

#Parameters
msg_digitalHeatmap_parameters<-list(
	title="Parameters for position-based heatmap",
	text=list(
		"Set the parameters for a heatmap showing overlaps between multiple ROIs"
	)
)


help_digitalHeatmap_parameters_masterROI<-list(
title="Select the master ROI",
text=list("The reference ROI in which to see overlapping regions. The overlaps will be calculated against the genomic ranges of the master ROI",
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. If more than one is selected, it won't be possible to extract a 
				new ROI from heatmap clusters.<br> 
				<b>2</b>. The jaccard index matrix won't 
				be shown.",plural=T)	

	)

)
help_digitalHeatmap_parameters_ROItoview<-list(
title="The ROIs to investigate for overlapping regions",
text=list("Select one or more ROIs for which you want to test the overlaps within the master ROI(s) genomic ranges.
		Each genomic range of the master ROI(s) will be overlapped with all these ROIs"

)

)
help_digitalHeatmap_parameters_ROIordering<-list(
title="ROI ordering",
text=list("You can set the order by which ROIs are displayed in the heatmap."

)

)
help_digitalHeatmap_parameters_ROIforcluster<-list(
title="ROIs for cluster",
text=list("Select which is/are the ROI(s) that will drive the clustering of the rows 
			(genomic ranges) of the heatmap. The clustering will be based on the overlap (0 or 1) 
			in all the bins of each genomic range"

)

)
help_digitalHeatmap_parameters_clusterKmeans<-list(
title="K-means clustering",
text=list("Performs a K-means clustering on the rows of the heatmap. This allows the identification of specific patterns of overlap 
					between the different ROIs. User-defined K-means parameters are: number of clusters, 
					number of starting points and number of iterations"

)

)
help_digitalHeatmap_parameters_clusterHierarchical<-list(
title="Hierarchical clustering",
text=list("Performs a hierarchical clustering on the rows of the heatmap. User-defined parameters for hierarchical clustering 
					are: number of clusters, distance method and clustering method"

)

)
help_digitalHeatmap_parameters_clusternumber<-list(
title="Set number of clusters",
text=list("The number of clusters represents how many different overlapping patterns will be found by the clustering.
		Increase this number to detect lowly-represented overlapping patterns.",
		tags$br(),
		tags$br(),
		Warning("Maximum number allowed: 433")
)

)#identifies more subsets with a common overlap pattern
help_digitalHeatmap_parameters_nbins<-list(
title="Set the number of bins",
text=list("The number of bins represents
	the number of intervals each genomic range of the master ROI is divided into for calculating overlaps.", 
	tags$br(),
	tags$br(),
	Tip("Increase this number to improve resolution."),
	tags$br(),
	tags$br(),
	Tip("When you want to extract the most common patterns of overlap between a large number of ROIs, 
			set the number of bins =1, a k-mean clustering driven by all the shown ROIs and a number of clusters 
			k=2^n, where n=number of ROIs. In this way, all the possible combinations of overlaps between ROIs can 
			be displayed at once"),
	tags$br(),
	tags$br(),
	Warning("The plot 'Positional distribution of overlaps' will be shown only if the number of bins is > 2")

	)

)
help_digitalHeatmap_parameters_strandspecific<-list(
title="Strand specific overlaps",
text=list("Only the 
			overlaps in the same strand between master ROI and all other ROIs will be considered.
			If master ROI does not have strand information, this option will be ignored."

	)

)
help_digitalHeatmap_parameters_randomsample<-list(
	title="Random sample of genomic ranges to show",
	text=list("The number of genomic ranges of the master ROI to show in the heatmap (default: 2000 ranges). The smaller the subset, the faster the heatmap will be shown.",
		tags$br(),
		tags$br(),
		Tip("1. Usually, a random sample of 2000 genomic ranges is sufficient to have a good representation of 
			overlaps in the heatmap.<br>
			2. If the aim of the analysis is to extract new ROIs from heatmap clusters for downstream analyses, 
			include all the genomic ranges of the master ROI (set a value >= the sum of ranges of the master ROIs selected)",plural=T)
	)

)


help_digitalHeatmap_parameters_positionaloverlap<-list(
	title="positional overlap %",
	text="If checked, the plot shows the fraction of overlapping events instead of the absolute number"
)


#Heatmap
msg_digitalHeatmap_heatmap<-list(
	title="Position-based heatmap",
	text=list(
		"An interactive heatmap that shows, for each genomic range of the master ROI selected (each row) 
		the overlap of the other ROI(s) selected. Note that overlaps are calculated for each bin of the genomic ranges.
		Genomic ranges of the ROIs selected are displayed as colors.
		The fraction of the overlapping genomic range for each ROI compared to the total ranges of the master 
		ROI is displayed on the names in the x axis. Clusters are shown on the left side of the heatmap, as colored bars.",
		tags$br(),
		tags$br(),		
		"If a cluster bar is clicked, the number of genomic ranges belonging to the selected clusters are displayed, 
		along, with the cluster number. The field “New ROI from selection” will also appear, giving you the possibility 
		to create a new ROI from the genomic ranges of the cluster selected.",
		tags$br(),
		tags$br(),
		Comment("The new ROI will keep the enrichments associated to the old ROI"),
		tags$br(),
		tags$br(),
		Warning("Cluster's bars appear only if you have selected a single master ROI"),
		tags$br(),
		tags$br(),
		Tip("If the aim of the clustering analysis is to identify all the regions of the genome that share a defined 
			combination of overlaps, it is important to include all genomic ranges of the master ROI and not a random subset")			
	)
)


#help for cluster, fraction ranges and extract from cluster
help_digitalHeatmap_clickinfo<-list(
	title="Click clusters to get info and extract ROIs",
	text=list(HTML("Click one of the coloured blocks in the cluster bar. On the right, you will see how many genomic regions 
		constitute the selected cluster and you can extract a new ROI having these genomic ranges.<br><br>
		<img src='digital_click_cluster.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),
		Tip("Use this to extract ROIs with interesting overlapping patterns for downstream analyses, for example the gene ontology.")	
	)	
)

help_digitalHeatmap_fractionranges<-list(
	title="The number of genomic ranges displayed",
	text=list("This is the number of genomic ranges shown in the heatmap / the total genomic ranges (the sum of 
			ranges of all the master ROI(s) selected).",
		tags$br(),
		tags$br(),			
		Comment("This is set in 'Random sample of genomic ranges' parameter. These numbers are green when all the available ranges
			in the master ROI(s) are shown.")
	)
)

help_digitalHeatmap_extractROI<-list(
	title="Extract new ROIs",
	text=list("Here you can import a new ROI from the genomic ranges in the selected cluster.",
		tags$br(),
		tags$br(),	
		Comment("1. To extract new ROI, it's strongly suggested to use all the ranges of the master ROI ('Random sample of genomic ranges' parameter).
				In this way, no ranges are randomly lost for downstream analyses.<br>
				2. The new ROI will preserve all the enrichments associated to the master ROI, if present.")

	)
)






#Jaccard idx
msg_digitalHeatmap_jaccardIdx<-list(
	title="Shows the heatmap of Jaccard index ",
	text=list(
		"This plot shows the  heatmap of Jaccard index of all the pairwise combinations of the ROIs displayed 
		in the heatmap.",
		tags$br(),
		tags$br(),
		Comment("For a given pair of ROIs shown in the position-based heatmap, JI=(number of genomic ranges of the 
			master ROI overlapping both with the first and the second ROI) / (number of genomic ranges of the 
			master ROI overlapping only with first ROI + number of genomic ranges of the master ROI overlapping 
			only with second ROI + number of genomic ranges of the master ROI overlapping both with the first and 
			the second ROI)"),
		tags$br(),
		tags$br(),
		Warning("This plot is generated only if a single master ROI has been selected and if the number of 
				ROIs shown in the heatmap is > 1")			
	)
)


#Overlap frequency bias
msg_digitalHeatmap_overlapBias<-list(
	title="Shows the positional distribution of overlaps",
	text=list(
		"This plot is a meta-representation that shows the frequency of overlap of a ROI over the 
		genomic ranges of a reference ROI. It can be used to evaluate if genomic regions overlap with 
		a preferential positional pattern.",
		tags$br(),
		tags$br(),
		Comment("For example, when analysing the overlap of a ChIP-seq transcription factor on a series of 
			transcripts, this plot can be used to evaluate if the overlaps are preferentially found at the 
			5’ or 3’ end of the genes"),
		tags$br(),
		tags$br(),
		Warning("This plot is displayed only if the number of bins is >2."),
		tags$br(),
		tags$br(),
		Tip("1. To obtain a sharper profile of the plot, increase the number of bins.<br>
			2. You can choose to show the fraction of overlapping genomic ranges instead of the absolute number 
			by checking '<b>positional overlap %</b>' option in parameters. This is useful when comparing multiple master 
			ROIs (it normalizes for the number of genomic ranges of each master ROI)",plural=T)

	)
)


#### Analogic heatmap

#Parameters
msg_analogicHeatmap_parameters<-list(
	title="Parameters for enrichment-based heatmap",
	text=list(
		"Set the parameters for a heatmap showing reads enrichments inside ROIs"
	)

)



help_analogicHeatmap_parameters_ROI<-list(
	title="Select one or more ROIs",
	text=list("The ROIs in which investigate enrichments. The corresponding genomic ranges will be the 'rows' of the heatmap.",
		tags$br(),
		tags$br(),
		Warning("To plot the heatmap, enrichments must be associated to the ROI(s) selected")
	)
)


help_analogicHeatmap_parameters_enrichments<-list(
	title="Select enrichments to show",
	text=list("Select one or more enrichments to show in the heatmap.",
		tags$br(),
		tags$br(),
		Tip("To associate enrichments to ROIs, go to 'ROI preparation' section"),
		tags$br(),
		tags$br(),
		Warning("An enrichment appears only if associated to ALL the ROIs selected")
	)
)



help_analogicHeatmap_parameters_enrichmentorder<-list(
	title="Enrichment order",
	text=list("You can set the order by which enrichments are displayed in the heatmap."
	)
)


help_analogicHeatmap_parameters_ranking<-list(
	title="Rank enrichments",
	text=list("Genomic ranges (rows of the heatmap) will be put in descending order based on a selected enrichment."
	)
)

help_analogicHeatmap_parameters_clustering<-list(
	title="Cluster enrichments",
	text=list("The order of heatmap rows (genomic ranges) will be calculated with a clustering algorithm. The enrichment
			that will drive the clustering have to be chosen.",
			Warning("Depending on your hardware, hierarchical clustering may lead to crashes if you cluster more than 20000 genomic regions")
	)
)



help_analogicHeatmap_parameters_clusterKmeans<-list(
	title="K-means clustering",
	text=list("Performs a K-means clustering on the rows of the heatmap. This allows the identification of specific patterns of enrichment. 
		User-defined K-means parameters are: number of clusters, 
					number of starting points and number of iterations"
	)
)
help_analogicHeatmap_parameters_clusterHierarchical<-list(
	title="Hierarchical clustering",
	text=list("Performs a hierarchical clustering on the rows of the heatmap. User-defined parameters for hierarchical clustering 
					are: number of clusters, distance method and clustering method"
	)
)


help_analogicHeatmap_parameters_clusternumber<-list(
	title="Set number of clusters",
	text=list("The number of clusters represents how many different patterns of enrichment will be detected by the clustering.
		Increase this number to detect lowly-represented patterns.",
		tags$br(),
		tags$br(),
		Warning("Maximum number allowed: 433")
	)
)



help_analogicHeatmap_parameters_nbins<-list(
	title="Set the number of bins",
	text=list("The number of bins represents
		the number of intervals each genomic range of the ROI(s) is divided into for calculating enrichments.", 
		tags$br(),
		tags$br(),
		Tip("Increase this number to improve resolution.")

	)
	
)

help_analogicHeatmap_parameters_subsample<-list(
	title="Random sample of genomic ranges to show",
	text=list("The number of genomic ranges of the ROI(s) to show in the heatmap (default: 2000 ranges). The smaller the subset, the faster the heatmap will be shown.",
		tags$br(),
		tags$br(),
		Tip("1. Usually, a random sample of 2000 genomic ranges is sufficient to have a good representation of 
			enrichments in the ROI(s).<br>
			2. If the aim of the analysis is to extract new ROIs from heatmap clusters for downstream analyses, 
			include all the genomic ranges of the ROI (set a value >= the sum of ranges of all the ROIs selected)",plural=T)
	)
)


help_analogicHeatmap_parameters_uniform<-list(
	title="Global color scale",
	text=list("The threshold for the color scale saturation will be a quantile the 
					enrichments considered together.",
			tags$br(),
			tags$br(),
			Comment("Use this option if you have to compare signals of different enrichments.")
	)
)
help_analogicHeatmap_parameters_individual<-list(
	title="Individual color scale",
	text=list("Set a different threshold of the color scale 
					saturation for each of the signal enrichments displayed.",
			tags$br(),
			tags$br(),
			Comment("Use this option if you include different enrichments in the same heatmap with values 
				that have different order of magnitudes and do not have to be compared.")
	)
)


#Heatmap
msg_analogicHeatmap_heatmap<-list(
	title="Enrichment-based heatmap",
	text=list(
		"Plots an interactive heatmap with each row representing genomic ranges in the ROI(s) selected 
		and each column the different enrichments. Signal intensity will be depicted in the specified color scale. 
		If clustering has been performed, different clusters will be shown at the left of the heatmap as bars 
		composed by different colors",
		tags$br(),
		tags$br(),
		Comment("1. The user can select a portion of the heatmap (by dragging the pointer) or a cluster bar 
				(by clicking the bar); these selections will update all the plots shown below (i.e. profiles, boxplots, 
				correlation heatmaps). <br>
				2. The number of the genomic ranges selected will be shown at the upper-right 
				corner of the heatmap.<br> 
				3. Upon selection the new field 'New ROI from selection' will appear: this will allow the user to create a new ROI containing the genomic ranges selected in the heatmap.<br> 
				4. The new ROI will keep the enrichments associated to the old ROI.",plural=T),
		tags$br(),
		tags$br(),
		Warning("A single ROI must be selected in order to extract a new ROI from the Heatmap 
				and to have the barplot of the clusters"),
		tags$br(),
		tags$br(),

		Tip("To extract new ROIs from clusters or selected areas of the heatmap, make sure to use all the genomic 
			ranges of the ROI, and not a random sample (Set 'Random sample of:' option with a value equal to or greater the number of genomic ranges in the ROI(s) selected)")
		

	)
)



#help for cluster, fraction ranges and extract from cluster
help_analogicHeatmap_clickinfo<-list(
	title="Click clusters/heatmap area to get info and extract ROIs",
	text=list(HTML("Select an area of interest in the heatmap by dragging and dropping the mouse,
		 or click one of the coloured blocks in the cluster bar. On the right, you will see how many genomic regions 
		constitute the selection and you can extract a new ROI having these genomic ranges. On the bottom, enrichments of the selected area will
		be shown as profiles and boxplots, and updated on the fly if you change the selected area.<br><br>
		<img src='analog_click_cluster.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),
		Tip("Use this to extract ROIs with interesting enrichment patterns for downstream analyses, for example the gene ontology.")	
	)	
)


help_analogicHeatmap_clickinfoonlyselection<-list(
	title="Select a heatmap area to get info and extract ROIs",
	text=list(HTML("Select an area of interest in the heatmap by dragging and dropping the mouse. On the right, you will see how many genomic regions 
		constitute the selection and you can extract a new ROI having these genomic ranges. On the bottom, enrichments of the selected area will
		be shown as profiles and boxplots, and updated on the fly if you change the selected area.<br><br>
		<img src='analog_click_area.png' alt='Res' width=60% height=60%>"),
		tags$br(),
		tags$br(),
		Tip("Use this to extract ROIs with interesting enrichment patterns for downstream analyses, for example the gene ontology.")	
	)	
)

help_analogicHeatmap_fractionranges<-list(
	title="The number of genomic ranges displayed",
	text=list("This is the number of genomic ranges shown in the heatmap / the total genomic ranges (the sum of 
			ranges of all the ROI(s) selected).",
		tags$br(),
		tags$br(),			
		Comment("This is set in 'Random sample of genomic ranges' parameter. These numbers are green when all the available ranges
			in the master ROI(s) are shown.")
	)
)

help_analogicHeatmap_extractROI<-list(
	title="Extract new ROIs",
	text=list("Here you can import a new ROI from the genomic ranges in the selected cluster or area.",
		tags$br(),
		tags$br(),	
		Comment("1. It's strongly suggested to use all the ranges of the ROI ('Random sample of genomic ranges' parameter).
				In this way, no ranges are randomly lost for downstream analyses.<br>
				2. The new ROI will preserve all the enrichments associated to the reference ROI, if present.")

	)
)











#Profiles
msg_analogicHeatmap_profiles<-list(
	title="Enrichment profiles",
	text=list(
		"Plots the profile of the enrichment of the selected area of the heatmap or the selected cluster",
		tags$br(),
		tags$br(),
		Tip("For a better resolution of the profiles, increase the number of bins in the heatmap")		
	)
)
#Enrichments
msg_analogicHeatmap_enrichments<-list(
	title="Enrichment plots",
	text=list(
		HTML("Plots the enrichment of the selected area of the heatmap, or the selected cluster.<br> 
			Plotting options are:"),
		tags$br(),
		tags$br(),
		Field("Boxplot by ROI/cluster","The enrichments are shown as boxplot stratified by the ROIs 
			provided. In the case of a clustered heatmap from a single ROI, the boxplot will be stratified by clustering"),
		Field("Boxplot by enrichment","The enrichments are shown in different ROIs as boxplot stratified by enrichment")	

	)
)




#### Enrichment in ROIs

#Parameters
msg_enrichmentInRois_parameters<-list(
	title="Parameters for visualizing enrichments in ROIs",
	text=list(
		"Define parameters for the analysis of the enrichments in 
		the ROIs: signal distribution (profiles), box-plots, signal correlations matrix and scatterplot"
	)
)


help_enrichmentInRois_parameters_ROIs<-list(
	title="Select ROI(s)",
	text=list("Choose one or more ROI(s) for the enrichment analysis.",
		tags$br(),
		tags$br(),
		Warning("For the analysis, enrichments must be associated to the ROI(s) selected.")

	)
)

help_enrichmentInRois_parameters_enrichments<-list(
	title="Select enrichments to show",
	text=list("Choose the enrichments to be shown associated to the selected ROI(s).",
		tags$br(),
		tags$br(),
		Tip("To associate enrichments to ROIs, go to 'ROI preparation' section."),
		tags$br(),
		tags$br(),
		Warning("An enrichment appears only if associated to ALL the ROIs selected.")
	)
)


help_enrichmentInRois_parameters_bins<-list(
	title="Set the number of bins",
	text=list("The number of bins represents
		the number of intervals each genomic range of the ROI(s) is divided into for calculating enrichments.", 
		tags$br(),
		tags$br(),
		Tip("Increase this number to improve resolution."),
		tags$br(),
		tags$br(),
		Comment("This parameter only affects the profile plot.")		

	)
)


#total reads and read density are the same messages:
# help_singleEvaluation_normalizationtotalread
# help_singleEvaluation_normalizationreaddensity



#Profiles
msg_enrichmentInRois_profiles<-list(
	title="Enrichment profiles",
	text=list(
		"Plots the profiles of the enrichment of the selected ROI(s) and enrichment(s)",
		tags$br(),
		tags$br(),
		Tip("To increase the resolution of the profiles, increase the number of bins")		
	)
)

#Boxplots
msg_enrichmentInRois_boxplots<-list(
	title="Enrichment boxplots",
	text=list(
		"Shows the boxplots of the enrichment of the selected ROI(s) and enrichment(s), 
		grouped by either ROI or enrichment",
		tags$br(),
		tags$br(),
		Field("Boxplot by ROI","The boxplot is grouped by the ROIs"),
		Field("Boxplot by enrichment","The boxplot is grouped by the enrichments"),
		tags$br(),
		Comment("The raw data of the enrichments in ROIs can be downloaded in tabular format 
				by clicking 'download data'")
	)
)
#Correlations
msg_enrichmentInRois_correlations<-list(
	title="Correlation heatmaps",
	text=list(
		HTML("Shows the correlation or the partial correlation between enrichments in the selected ROI.<br>
		Clicking a square in the matrix will show the corresponding scatterplot of the pairwise 
		(partial)correlation in the box on the right.<br><br>
		<img src='cormatrix_click.png' alt='Res' width=70% height=70%>"),
		tags$br(),
		tags$br(),
		Field("Cor-Heatmap","Correlation Heatmap of the enrichments in a given ROI"),
		Field("Pcor-Heatmap","Partial correlation heatmap of the enrichments in a given ROI"),
		tags$br(),
		Comment("The partial correlation between enrichment A and B will be computed using, 
					as covariates, all the other enrichments selected"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. The correlation and partial correlation heatmaps will appear only when a single ROI 
				has been selected.<br> 
				<b>2</b>. The correlation heatmap will appear only when at least 2 enrichment 
				were selected.<br> 
				<b>3</b>. Partial correlation heatmap only when at least 3 enrichments were selected",plural=T)
	)
)


#Scatterplot
msg_enrichmentInRois_scatterplot<-list(
	title="Scatterplot of pairwise correlation",
	text=list(
		"To display the scatterplot of the pairwise correlation/partial correlation, 
		click the corresponding cell in the heatmap of correlations",
		tags$br(),
		tags$br(),
		Warning("A scatterplot can be generated only if a correlation matrix has been generated.
			For computational reasons, a random sample of 2000 points is plotted")		
	)
)


#### Dynamics on genes

#Parameters
msg_dynamicsOnGenes_parameters<-list(
	title="Parameters for meta-gene representation",
	text=list(
		"Here you can define the parameters for the analyses of metagenes"
	)

)




help_dynamicsOnGenes_parameters_genelist<-list(
	title="Select genelist(s)",
	text=list("Choose the gene lists in which you want to perform the analysis",
		tags$br(),
		tags$br(),
		Tip("To import genelists, go to 'ROIs' -> 'Get promoters, transcripts, TES coordinates of a list of genes'"),
		tags$br(),
		tags$br(),
		Warning("Enrichment files must be associated to the selected genelist(s). To associate enrichments,
					go to 'ROI preparation' -> 'Prepare genelists for metagene profile'")
	)
)

help_dynamicsOnGenes_parameters_ROIassociated<-list(
	title="ROIs constituting the genelist(s)",
	text="This panel shows all the ROIs that constitute 
			the selected gene list(s). A genelist is composed by 3 ROIs: promoters, transcripts, TES of the genelist."
)

help_dynamicsOnGenes_parameters_ernichments<-list(
	title="Select enrichments to show",
	text=list("Choose the enrichments to be shown associated to the selected genelist(s).",
		tags$br(),
		tags$br(),
		Tip("To associate enrichments to genelists, go to 'ROI preparation' -> 'Prepare genelists for metagene profile', or 
				do it manually in 'ROI preparation' -> 'Manual ROI management' (for advanced users)"),
		tags$br(),
		tags$br(),
		Warning("<b>1</b>. An enrichment appears only if associated to ALL the genelists selected<br>
				<b>2</b>. An enrichment is considered associated to a genelist when it was associated to both promoters, transcripts and TES of that genelist",plural=T)
	)
)

help_dynamicsOnGenes_parameters_nbins<-list(
	title="Set the number of bins",
	text=list("The number of bins represents
		the number of intervals each gene of the list is divided into for calculating enrichments.", 
		tags$br(),
		tags$br(),
		Tip("Increase this number to improve resolution."),
		tags$br(),
		tags$br(),
		Comment("This parameter only affects the profile plot.")
	)	
)

# help_singleEvaluation_normalizationtotalread
# help_singleEvaluation_normalizationreaddensity

help_dynamicsOnGenes_parameters_fractionexclude<-list(
	title="Fraction of outlayers to exclude",
	text=list("Fraction of outliers to exclude in cumulative plots",
		tags$br(),
		tags$br(),
		Tip("Increase this value to emphasize the differences between the curves in the cumulative plots")
	)
)





	
#Profiles
msg_dynamicsOnGenes_profiles<-list(
	title="Metagene analysis",
	text=list(
		"This analysis generates different plots on groups of genes",
		tags$br(),
		tags$br(),
		Field("Metagene profile","Shows the metagene profile of enrichments in the selected gene list(s).
		This plot is particularly useful to evaluate RNApol2 distribution along genes"),
		Comment("The metagene plot starts at –30% of transcripts length from TSS and ends at +30% from TES"),
		tags$br(),
		tags$br(),
		HTML("A metagene representation is useful to analyse the behaviour of a signal in a set of genes. Each transcript length is increased by 30% 
				and divided in an equal number of bins (in the following example, 9 bins). Then, the signals from different genes
				belonging to the same bin are averaged together to produce the metagene signal. <br><img src='scheme_metagene.png' alt='Res' width=100% height=80%>"),
		tags$br(),
		tags$br(),
		Field("Enrichment boxplots","Shows the the enrichments as boxplots at the promoters, transcripts and TES of the selected gene list(s)."),				
		Comment("The upstream and downstream limits to define TSS and TES regions are set the 
				first time the genome assembly is loaded"),
		tags$br(),	
		tags$br(),	
		Field("Ranked enrichments","For each gene in the gene list(s) selected, it shows the ranked enrichments at promoters, at genebodies or their ranked stalling index."),
		Comment("The stalling index is defined as number of reads at TSS / number of reads at genebodies.
				This is typically used when calculating the Pol2 stalling index"),
		tags$br(),	
		tags$br(),	
		Tip("To change the size of intervals around TSS or TES , go to the 'Assembly' section and  
			change the default up/downstream values for the definition of TSS and TES regions, 
			then import the gene list again ('ROIs' section -> 'Get promoters, transcripts, TES coordinates of a list of genes')")

	)
)








#### GO analyses

#parameters
msg_goAnalysis_parameters<-list(
	title="Parameters for gene ontology",
	text=list(
		"Set the parameters used for the gene ontology analysis",
		tags$br(),
		tags$br(),
		Comment("'<b>Variables</b>' tab: basic parameters for the analysis; '<b>Filtering</b>' tab: parameters for filtering the results")
	)
)

help_goAnalysis_parameters_fromROI<-list(
	title="Analyze genes annotated to ROIs",
	text=list("Perform the gene ontology analysis on the genes annotated to select ROI(s).",
		tags$br(),
		tags$br(),
		Comment("If a single ROI is selected, the results will be displayed as a barplot, otherwise as a heatmap (rows=ROIs; columns=significant signatures in at least one ROI)"),
		tags$br(),	
		tags$br(),	
		Warning("A genome assembly must be loaded for the annotation to ROIs. To load a genome assembly, use 'Assembly' section")
	)
)



help_goAnalysis_parameters_fromgenelist<-list(
	title="Analyze genes from a pasted genelist",
	text=list("Perform the gene ontology analysis on a list of genes given on the fly, one per line.",
		tags$br(),tags$br(),
		Comment("If genes are lowercase (for example murine), they will be converted in uppercase (Human) to be consistent 
			with signature lists available in the program. The Output is a barplot showing the p adjusted of the hypergeometric test")
	)
)

help_goAnalysis_parameters_kindofID<-list(
	title="Kind of gene IDs",
	text=list("Which kind of identifiers are the input genes provided? Symbols, ENTREZ IDs, ...")
)

help_goAnalysis_parameters_nearestgenes<-list(
	title="Take the nearest gene",
	text="For each genomic range of selected ROI(s), the program takes the nearest gene from the midpoint.
		Therefore, the number of input genes of a ROI will be = number of genomic ranges"
)

help_goAnalysis_parameters_genewindow<-list(
	title="Take the genes inside genomic window",
	text="For each genomic range of selected ROI(s), the program takes all the genes within a certain genomic window 
		from the midpoint of the range. The amplitude of the window is set by the user"
)

help_goAnalysis_parameters_signatures<-list(
	title="Select signatures (database)",
	text=list("Select the gene signatures(s) to use for the gene ontology analysis. If a custom universe is not specified, the union of all the genes of 
			the signatures selected will be considered as the background for the hypergeometric test.",
			tags$br(),tags$br(),
			Comment("The gene signatures are lists of genes in gmt format present in the directory appContent/signatures/ of the program. 
				Gene signatures can be loaded from the Molecular Signature Database (MSigDB) or similar resources"),
			tags$br(),tags$br(),
			
			HTML("Preloaded gene signatures are:<br>
			<ul>
			<li><b>c1.all.v7.4_symbols.gmt</b>: positional genesets. Gene sets corresponding to human chromosome cytogenetic bands. </li>
			<li><b>c2.all.v7.4_symbols.gmt</b>: curated genesets. Gene sets in this collection are curated from various sources, including online pathway databases and the biomedical literature.
				<ul>
				<li><b>c2.cgp.v7.4_symbols.gmt</b>: Gene sets represent expression signatures of genetic and chemical perturbations. A number of these gene sets come in pairs: xxx_UP (and xxx_DN) gene set representing genes induced (and repressed) by the perturbation.</li>
				<li><b>c2.cp.v7.4_symbols.gmt</b>: Gene sets from pathway databases.</li>
				</ul>
			</li>
			<li><b>c3.all.v7.4_symbols.gmt</b>: Gene sets representing potential targets of regulation by transcription factors or microRNAs. The sets consist of genes grouped by elements they share in their non-protein coding regions.</li>
				<ul>
				<li><b>c3.mir.v7.4_symbols.gmt</b>: All miRNA target prediction gene sets.</li>
				<li><b>c3.tft.v7.4_symbols.gmt</b>: All transcription factor target prediction gene sets.</li>
				</ul>
			<li><b>c4.all.v7.4_symbols.gmt</b>: Computational gene sets defined by mining large collections of cancer-oriented microarray data.</li>
			<li><b>c5.all.v7.4_symbols.gmt</b>: Gene sets that contain genes annotated by the same ontology term.</li>
				<ul>
				<li><b>c5.go.v7.4_symbols.gmt</b>: Gene sets derived from all the GO processes</li>
				<li><b>c5.go.bp.v7.4_symbols.gmt</b>: Gene sets derived from the GO Biological Process ontology.</li>
				<li><b>c5.go.cc.v7.4_symbols.gmt</b>: Gene sets derived from the GO Cellular Component ontology.</li>
				<li><b>c5.go.mf.v7.4_symbols.gmt</b>: Gene sets derived from the GO Molecular Function ontology.</li>
				<li><b>c5.hpo.v7.4_symbols.gmt</b>: Gene sets derived from the Human Phenotype ontology.</li>
				</ul>
			<li><b>c6.all.v7.4_symbols.gmt</b>: Gene sets that represent signatures of cellular pathways which are often dis-regulated in cancer.</li>
			<li><b>c7.all.v7.4_symbols.gmt</b>: Gene sets that represent cell states and perturbations within the immune system.</li>
				<ul>
				<li><b>c7.immunesigdb.v7.4_symbols.gmt</b>: Gene sets representing chemical and genetic perturbations of the immune system generated by manual curation of published studies in human and mouse immunology.</li>
				<li><b>c7.vax.v7.4_symbols.gmt</b>: Gene sets curated by the Human Immunology Project Consortium (HIPC) describing human transcriptomic immune responses to vaccinations.</li>
				</ul>
			<li><b>c8.all.v7.4_symbols.gmt</b>: Gene sets that contain curated cluster markers for cell types identified in single-cell sequencing studies of human tissue.</li>
			<li><b>h.all.v7.4_symbols.gmt</b>: Hallmark gene sets summarize and represent specific well-defined biological states or processes and display coherent expression.</li>
			</ul>
			<br>
			For further info, go to MSigDB website: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C1
			")

	)
)

help_goAnalysis_parameters_orderresults<-list(
	title="How to order results",
	text="Choose the order of the significant hits resulting from GO analysis. If a single ROI has 
			been selected, or a custom gene list has been used as input, the results will be ranked based on 
			-log10 p adjusted by default. If multiple ROIs have been selected, the user can choose whether 
			to rank the heatmap based on the best p adjusted 
			or whether to order based on clustering (K-means or hierarchical)"
)

help_goAnalysis_parameters_minsize<-list(
	title="Minimum signature size",
	text="Set the minimum size of the gene signatures used for the ontology analysis. 
			All signatures with a size below this value will not be considered"
)

help_goAnalysis_parameters_maxsize<-list(
	title="Maximum signature size",
	text="Set the maximum size of the gene signatures used for the ontology analysis. 
			All signatures with a size above this value will not be considered"
)


help_goAnalysis_parameters_customuniverse<-list(
	title="Custom universe genes",
	text="Paste a list of gene symbols to be used as custom universe for the hypergeometric test in the 
		gene ontology analysis"
)



help_goAnalysis_parameters_generatio<-list(
	title="gene ratio threshold",
	text="All the hits below this gene ratio threshold are hidden from the plot. Gene ratio is defined as: genes associated to the 
			ROI or in the custom gene list / all the genes in the signature"
)

help_goAnalysis_parameters_padjthresh<-list(
	title="-log10 padj threshold",
	text="All the hits below this significance (-log10 p adjusted) threshold are hidden from the plot"
)







#GO plot
msg_goAnalysis_goPlot<-list(
	title="Gene ontology plot",
	text=list(
		"It displays the plot of the GO analysis (heatmap or barplot). 
		For each signature, genes annotated to a specific ROI can be displayed by clicking the corresponding cell or bar",
		tags$br(),		
		tags$br(),
		Comment("The p adjusted dispayed is the level of significance of the hypergeometric test.")	
	)
)



help_goAnalysis_clickinfoheat<-list(
	title="Click cells in the heatmap to get info",
	text=list(HTML("Click one of the cells in the heatmap. On the right, you will see the ontological term, the ROI and all the genes
			that are both annotated to the ROI and belonging to the ontological term.<br><br>
			<img src='GO_click_heat.png' alt='Res' width=80% height=80%>")
	)
)


help_goAnalysis_clickinfobar<-list(
	title="Click bars in the barplot get info",
	text=list("Click one of the bars. On the right, you will see the ontological term and all the genes given
			in the input that belong to that ontological term.")
)




#GO table
msg_goAnalysis_goTable<-list(
	title="Gene ontology table",
	text="GO results are displayed in a tabular format and can be explored interactively.
		You can download the tab-delimited text file for downstream processing"
)




############################################################
############################################################
# Save/Load
############################################################
############################################################

msg_saveLoad_save<-list(
	title="Save the session",
	text="Save the current working session. All ROIs created and the enrichments associated
	with them will be saved into an *rds file. This session file can then be shared with other users
	and opened with ChroKit, to reproduce all the data. In order to keep the sessions light and 
	computationally tractable, try not to go over 1Gb of session files."
)

msg_saveLoad_load<-list(
	title="Load the session",
	text="Load a session from an *rds file. Rds file should have been created with ChroKit."
)




