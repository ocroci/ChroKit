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
############################################################
# coordinate files
############################################################
############################################################

#choose coordinates
msg_coordinateFiles_chooseCoordinates<-list(
	title="From file",
	text=list(
		tags$h3("File to import"),
		"Choose the BED/GTF/GFF file to import",
		tags$br(),
		tags$br(),
		tags$h4("Parameters:"),
		"Choose the parameters to use to open a ROI from file",
		tags$br(),
		Field("Header","Tick if the file has a header"),
		Field("Lines to skip","If you have extra-lines at the beginning of your file, you can skip them by
			inserting their number here"),

		tags$br(),
		tags$h4("Open file:"),
		"You can use two different ways to open a BED/GTF/GFF file",
		tags$br(),
		Field("Select a file...","Select the BED/GFF/GTF file by exploring the filesystem"),
		Field("...or choose a file path","Type the complete path of the BED/GFF/GTF file"),
		tags$br(),

		Comment("If the file is opened correctly, the content of the table can be explored interactively 
						in 'File preview' section (see below)"),

		tags$br(),
		tags$br(),
		Warning("1. This file must be reachable from the machine on which ChroKit is running.<br>
				2. It must be a tab-delimited text file with 3 (optionally 4) columns with the format :<br><br>

				Column 1: chr&emsp;Column 2: start&emsp;Column 3: end&emsp;Column 4: strand<br><br>
				
				<b>'chr'</b> column should be in format 'chrN', where 'N' should be the number of the chromosome (1,2,...X,Y,M)<br>
				<b>'start, end'</b> columns are numbers indicating the starting and ending point of the each region, respectively<br>
				<b>'strand'</b> column (optional) represents the strand for each genomic range(can be either +, -, or * if strand is not determined)		
				",plural=T),
				
		tags$br(),
		tags$br(),
		tags$br(),
		tags$h3("File preview"),
		"Here you can check the file you opened"		
	)
)



#Import genelist
msg_genelists_importGenelist<-list(
	title="From genelist",
	text=list(
		tags$h3("Genes to import"),
		"You can import promoters, transcripts and TES from a custom list of genes (one gene per line) in 3 different ways",
		tags$br(),
		tags$br(),
		Field("Open a text file...","Select a text file containing the list of genes by exploring the filesystem"),
		Field("...or select a path...","Type the path of the text file containing the list of genes"),
		Field("...or put IDs/symbols here","Insert here the list of gene symbols or IDs to extract, and the name of the new gene list you are going to import"),

		tags$br(),
		Comment("When a gene list is opened, the promoters, transcripts, TES of the genes associated to that gene lists are loaded in memory as new ROI. 
				All annotated isoforms are loaded as well"),
		tags$br(),
		tags$br(),
		Warning("1. The file containing the gene list must be a text file in which each row is a gene symbol/gene ID.<br> 
				2. It must be readable from the system in which ChroKit is running.<br>
				3. A genome assembly must be loaded to import a gene list (To load a genome assembly, go to the 'Assembly' section).",plural=T),

		tags$br(),
		tags$br(),
		tags$h3("Parameters"),
		tags$br(),
		Field("What kind of identifiers are you importing?","Select which kind of identifiers are in the gene list you are importing: Symbols, ENTREZ IDs, ..."),
		Field("Max length for transcripts","Select the maximum length of the transcript allowed in the gene list you are going to import.  
				Genes with transcripts length above that threshold will not be loaded"),
		tags$br(),
		Warning("Association of enrichments to long transcripts may cause memory problems")

	)
)


#Delete ROIs
msg_deleteRois_deleteRois<-list(
	title="Remove ROIs",
	text=list(
		"Delete one or more ROIs",
		tags$br(),
		tags$br(),
		Tip("Remove unnecessary ROIs if your session becomes too heavy (too many ROIs with too many enrichments associated to them). This 
			will save memory and make the session files smaller"),
		tags$br(),
		tags$br(),
		Warning("When you remove a ROI, all the enrichments associated to it will be removed, as well")
	)
)
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
# enrichment files
############################################################
############################################################

#import enrichment file
msg_enrichmentFiles_importEnrichment<-list(
	title="Import an enrichment file (BAM/WIG)",
	text=list(
		"This section is used to associate enrichment files to a ChroKit session.",
		tags$br(),
		tags$br(),
		Field("Choose a file...","Select the enrichment file (BAM or WIG) by exploring the filesystem"),
		Field("...or select the path","Type here the path to the enrichment file (BAM or WIG)"),
		tags$br(),
		Comment("Formats supported: BAM and WIG. When an enrichment file is selected and imported, the program keeps in memory the link (path)
				 to that file"),
		tags$br(),
		tags$br(),
		Warning("1. The file must be reachable from the machine in which ChroKit is running.<br>

				2. If a BAM file is chosen, a BAM index (.bai) file should be present in the same directory 
				and must have the same name of the bam, plus the '.bai' extension. 
				For example, if the file is 'enrichment.bam', the corresponding index should be 'enrichment.bam.bai'. 
				WIG files are not currently supported by Windows operating systems",plural=T),
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

############################################################
############################################################
# Assembly
############################################################
############################################################

#Extract annotated elements from database
msg_databases_extractAnnotatedElements<-list(
	title="Extract annotated promoters, transcripts, TES coordinates of a genome assembly",
	text=list(

		Field("Choose assembly to use","Here you can choose the genome assembly to use for extracting annotated promoters, 
										transcripts and TES (transcription end sites)"),
		tags$br(),

		Comment("1. The selected assembly will be used to automatically annotate all the ROIs that will be imported later.<br>
				2. Additional genome assemblies can be downloaded from bioconductor using the 'Download databases' box 
				on the right. <br>
				3. ROIs already imported into the program will be automatically re-annotated according to the newly selected genome assembly.",plural=T),
		
		tags$br(),
		tags$br(),
		Field("TSS and TES Upstream/Downstream","By indicating the number of base pairs selected upstream and downstream from TSS and TES, 
									you can choose the size of the genomic window encompassing the TSS and TES"),
		tags$br(),
		Tip("By choosing 2000 upstream and 1000 downstream  of  a TSS,  you  will define promotorial regions of 3000bp, 
					starting at  -2000 and ending at +1000 from the annotated TSS"),
		
		tags$br(),
		tags$br(),
		Field("Transcripts for annotation must contain (click)","Tick here to restrict the annotation of your ROI(s) 
								only to those transcripts annotated by the indicated identifier"),
		tags$br(),
		Tip("To select only promoters/transcripts/TES which have the corresponding gene symbol and also 
					a RefSeq ID, click 'SYMBOL' and 'RefSeq' buttons")
		
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
		Warning("1. Be sure your internet connection is working! <br>
				2. This may require some time",plural=T)		
	
	)
)






############################################################
############################################################
# ROI management
# IMAGES NEEDED FOR OVERLAPS
############################################################
############################################################

#### Overlaps
msg_newRois_options<-list(
	title="Options for the overlaps",
	text=list(
		"Here you can set the rules for defining overlapping or non-overlapping genomic ranges",
		tags$br(),
		tags$br(),
		Field("Minimum number of bp to consider for overlaps","Two genomic ranges will be deemed overlapping if 
						the number of shared bases is equal or greater than the value indicated here"),
		Field("Strand-specific overlaps","If selected, the overlaps will be determined considering the strand 
						information: only genomic regions with concordant strand information (i.e. same strand) 
						will be tested. If strand information is not available, Chrokit will ignore this option"),
		tags$br(),
		Warning("The enrichments associated to the original ROI will be kept 
				only if all the genomic ranges of the new ROI are narrower than those of the original ROI."),
		tags$br(),
		tags$br(),
		Field("Name of the ROI","The name of the new ROI that will be created by the overlap analysis"),
		tags$br(),
		Warning("You cannot use 'promoters', 'transcripts', 'TES' as names. The name cannot begin with 
					'promoters_genelist_', 'transcripts_genelist_' or 'TES_genelist_'")
		
	)
)


#ROI combination
msg_newRois_ROIcombination<-list(
	title="Roles and rules for the overlap analysis",
	text=list(
		"You can build new ROIs according to combination of overlaps with other ROIs",
		tags$br(),
		tags$br(),
		Field("Choose the reference ROI","Here you can define the ROI that will be used as reference 
				for calculating the overlaps"),
		Field("Select ROI(s)","Select one or more ROIs. If more than one ROI is selected, 
				a combination of them will be used as reference"),
		Field("Criteria for building the aggregated reference ROI","If more than one ROI is selected, choose how 
				multiple reference ROIs are aggregated together: the resulting reference ROI can be either 
				the union or the intersection of the ROIs selected<br>
				<li>Intersection: the aggregated reference ROI will be constituted by 
					the genomic ranges common to all the ROIs selected. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated reference ROI will be: 
					A&cap;B&cap;C</li>
				<li>Union: The aggregated reference ROI will be the union of all the genomic 
					ranges contained in the selected ROIs. For example, if 3 ROIs 
					are selected (i.e. A, B and C), the aggregated reference ROI will be:
					A&cup;B&cup;C</li>"),
		tags$br(),
		Comment("1. If multiple ROIs are aggregated together to form a reference, the program will calculate 
				the union or intersection of genomic ranges of those ROIs, and will produce a single reference.<br>
				2. The new ROI will keep the enrichments associated to the reference ROI only if a single reference ROI is selected",plural=T),
		tags$br(),
		tags$br(),
		tags$br(),
		Field("...that overlaps with contrast ROI","Here you can define the rules to build the contrast ROI 
				that will be used for the overlap analysis with the reference ROI. The aggregated contrast ROI 
				is built from the ROIs selected which are then aggregated following the rule defined in 
				'<b>Criteria for building the aggregated contrast ROI</b>'"),
		Field("Select ROI(s)","Select one or more contrast ROIs to be overlapped with the reference ROI"),
		Field("Criteria for building the aggregated contrast ROI","If multiple contrast ROIs are selected here, 
				they will be aggregated according to different criteria:
				<li>Intersection of the contrast ROIs: the aggregated contrast ROI will be constituted by 
					the genomic ranges common to all the ROIs selected. For example, if 3 contrast ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be: 
					A&cap;B&cap;C</li>

				<li>Union of the contrast ROIs: The aggregated contrast ROI will be the union of all the genomic 
					ranges contained in the selected ROIs. For example, if 3 contrast ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be:
					A&cup;B&cup;C</li>"),
		tags$br(),
		Comment("Only genomic ranges of the reference ROI that overlap with those of the contrast ROI will be kept"),
		tags$br(),
		tags$br(),
		tags$br(),
		Field("...that doesn't overlap with contrast ROI","Here you can define the rules to build the contrast ROI 
				that will be used to identify genomic ranges non-overlapping with the reference ROI. 
				For more information refer to:'<b>Criteria for building the aggregated contrast ROI</b>'"),
		Field("Select ROI(s)","Select one or more contrast ROIs"),
		Field("Criteria for building the aggregated contrast ROI","If multiple contrast ROIs are selected here, 
				they will be aggregated according to different criteria:
				<li>Intersection of the contrast ROIs:  The aggregated contrast ROI will be constituted by the 
					genomic ranges common to all the ROIs selected. For example, if 3 contrast ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be: 
					A&cap;B&cap;C</li>
				<li>Union of the contrast ROIs :  The aggregated contrast ROI will be the union 
				of all the genomic ranges contained in the selected ROIs. For example, if 3 contrast ROIs 
					are selected (i.e. A, B and C), the aggregated contrast ROI will be:
					A&cup;B&cup;C</li>"),
		tags$br(),
		Comment("The resulting ROI will contain only the genomic ranges of the reference ROI that 
				do not overlap with those of contrast ROI")
		

	)
)




#### Modify ROIs

#Resize
msg_modifyRois_resize<-list(
	title="Resize a ROI",
	text=list(
		"This lets you choose the width of the genomic ranges present in a ROI. The center of each genomic range 
		is used as reference point, except when promoters are chosen. If the ROI is a set of promoters, 
		then the center of the genomic ranges will correspond to the TSS (transcription start site)",
		tags$br(),
		tags$br(),
		Field("Select ROI to resize","Select the ROI you want to resize"),
		Field("Select Upstream/Downstream intervals","Set the number of base 
						pairs upstream and downstream the center of the genomic regions"),
		Field("Name of the ROI","Choose the name of the resized ROI"),
		tags$br(),
		Warning("The new ROI will keep the enrichments associated to the old ROI only if 
				all genomic ranges of the old ROI decrease their width")
	)
)

#Center on summit
msg_modifyRois_summit<-list(
	title="Center a ROI on the summit of an enrichment",
	text=list(
		"This section calculates the summit of the signal in each genomic range of a ROI and 
		then center the genomic range on the summit of a signal. For each range of the ROI, 
		the summit of the signal enrichment (i.e. the position with the maximum value of pileup reads) 
		will be calculated. A new ROI, centred on the summit and with each range of width = 1bp, 
		will be created (this corresponds  to coordinates of the summits)",
		tags$br(),
		tags$br(),
		Field("Select ROI to center on summit","Select the ROI for which you want to extract the summit"),
		Field("Select enrichment to use for summit","Select the enrichment used to calculate the summit"),
		Field("Name of the ROI","Choose the name of the ROI centred on summit"),
		tags$br(),
		Warning("1. The enrichment used for calculating the summit must have been previously 
				associated to the selected ROI (use 'Associate enrichments' tab). <br>
				2. The new 
				ROI won't keep the enrichments associated to the old ROI",plural=T)
		
	)
)



#Random sample
msg_modifyRois_sample<-list(
	title="Random subset a ROI",
	text=list(
		"Create a new ROI from a random sample of genomic rages of another ROI",
		tags$br(),
		tags$br(),
		Field("Select ROI to sample","Select the ROI to subsample from the list. The new ROI 
							will be a random subset of the original one"),
		Field("Fraction to keep","Select the number or fraction of genomic ranges. 
								You can scroll the bar to select the fraction of ranges to keep, 
								or directly put the number in the field below"),
		Field("Name of the ROI","Choose the name of the subsetted ROI"),
		tags$br(),
		Tip("Subsetting is useful to reduce memory usage. This is strongly suggested if you are 
			dealing with large ROIs (with above 20/30000 genomic ranges), since it reduces processing 
			time and memory usage during some downstream analysis (such as the association of enrichments), 
			but preserve a statistically significant sample representative of the original ROI"),
		tags$br(),
		tags$br(),
		Comment("The new ROI will keep the enrichments associated to the old ROI")
		
	)
)
#Filter for width
msg_modifyRois_width<-list(
	title="Subset a ROI based on the width of its ranges",
	text=list(
		"ROIs can be subsetted according to the width of their genomic ranges. Thresholds for 
		selection can be absolute (by inserting min and max values) or relative (by defining quantiles 
		of the witdh distribution)",
		tags$br(),
		tags$br(),
		Field("Select ROI to filter","Select the ROI to filter from the list"),
		Field("Select quantiles","Genomic Ranges with a width quantile below the minimum 
								and above maximum will be discarded"),
		Field("Min width/Max width","Genomic Ranges with a width below the minimum and above maximum will be discarded"),
		Field("Name of the ROI","Choose the name of the subsetted ROI"),
		tags$br(),
		Comment("The new ROI will keep the enrichments associated to the old ROI")
		
	)
)
#Filter for enrichment
msg_modifyRois_enrichment<-list(
	title="Subset a ROI based on a specific enrichment",
	text=list(
		"ROIs can be subsetted by their signal enrichment.
		Thresholds for subsetting can be absolute (by inserting min and max values) or relative 
		(by defining quantiles of the signal distribution). The new ROI will have only the genomic 
		ranges that have the enrichment levels in the range defined by the user",
		tags$br(),
		tags$br(),
		Field("Select ROI to filter","Select the ROI to filter for enrichment from the list"),
		tags$br(),
		Warning("At least one enrichment must be associated to that ROI"),
		tags$br(),
		tags$br(),
		Field("Select enrichment to use for filtering","Select one of the enrichments associated to the selected ROI"),
		tags$br(),
		Tip("To associate enrichments to ROIs, use the 'Associate enrichments' tab"),
		tags$br(),
		tags$br(),
		Field("Select quantiles","Genomic ranges with an enrichment quantiles below the minimum and 
									above maximum will be discarded"),
		tags$br(),
		Field("Min enrichment/Max enrichment","Genomic Ranges with an enrichment below the 
										minimum and above maximum values will be discarded"),
		Field("Name of the ROI","Choose the name of the subsetted ROI"),
		tags$br(),
		Comment("The new ROI will keep the enrichments associated to the old ROI")

	)
)

#Extract patterns
msg_modifyRois_pattern<-list(
	title="Extract sequence patterns",
	text=list(
		"You can extract user-provided sequence patterns, either from a ROI or from the entire genome.
			The output is a new ROI composed of genomic ranges that contain the user-defined pattern. 
			The new ROI obtained after the pattern search will be strand-specific, based on the pattern found. 
			The genomic ranges will be centred on the pattern.",
		tags$br(),
		tags$br(),
		Field("Choose where to search for the pattern","This menu allows to choose if the pattern is searched 
					in a existing ROI or from the entire genome. The genomic ranges of the new ROI will be all 
					the occurrences of the pattern:
					<li>From a ROI: select this to extract a pattern from an available ROI. This will let you 
						choose the ROI from which extract the desired pattern</li>
					<li>From the entire genome (SLOW): select this to extract a pattern from the entire genome 
						(using the assembly in use)</li>"),
		tags$br(),
		Warning("1. A genome assembly must be loaded to search for patterns. To load a genome assembly, use 'Assembly' section.<br>
				2. A BSgenome database must be also loaded: this database contains the information about the sequence 
				in a specific genomic region in the genome assembly in use. If this database is not installed in the 
				system, a new button will appear to enable the User to download it from bioConductor. 
				Make sure your internet connection is working",plural=T),
		tags$br(),
		tags$br(),
		Field("Select pattern (IUPAC nomenclature)","Select the pattern you want to extract. Use the IUPAC nomenclature:
				<li><b>A,C,G,T,U</b>: known bases</li>
				<li><b>R</b>: A or G</li>
				<li><b>Y</b>: C or T</li>
				<li><b>S</b>: G or C</li>
				<li><b>W</b>: A or T</li>
				<li><b>K</b>: G or T</li>
				<li><b>M</b>: A or C</li>
				<li><b>B</b>: C or G or T</li>
				<li><b>D</b>: A or G or T</li>
				<li><b>H</b>: A or C or T</li>
				<li><b>V</b>: A or C or G</li>
				<li><b>N</b>: any base</li>
				<li><b>. or -</b>: gap</li>"),
		Field("Strand selection","If a ROI with strand information has been selected, then a menu will appear to let
				 you decide whether to perform the pattern search on both strands or on a single strand"),
		Field("Name of the ROI","Choose the name of the ROI centred on the pattern"),
		tags$br(),
		Warning("The new ROI won't keep the enrichments associated to the old ROI")

	)
)



#### View ROIs

#ROI selection
msg_viewRois_selectRoi<-list(
	title="Select the ROI",
	text="Select one or more ROIs for which you want to obtain information (width/number of ranges and/or how it was created)"
)
#Visualization
msg_viewRois_visualization<-list(
	title="Visualize features of selected ROI(s)",
	text=list(
		Field("width distribution","It shows the distribution of the widths of the genomic ranges of the ROI(s)"),
		Field("intervals number","It shows a barplot representing the number of genomic ranges for each of the ROI(s)")

	)
)

#Information
msg_viewRois_information<-list(
	title="View info of a selected ROI",
	text=list(
		Field("Select quantiles of the width distribution","Shows how many genomic ranges 
							have a width greater than the selected width quantile"),
		Field("Where does this ROI come from?","Shows how the selected ROI was created (summarizes all the steps)"),
		tags$br(),
		Warning("To view information, only one ROI must be selected from the 'ROI selection' menu")

	)
)

#### Get ROI

#ROI selection
msg_getRois_roiSelection<-list(
	title="Get all the info associated to a ROI",
	text=list(
		"Select the ROI for which you want to explore or download the information. 
		Once a ROI is selected, the following information can be retrieved",
		tags$br(),
		tags$br(),
		Field("View Ranges","Genomic coordinates of the ROI and the strand"),
		Field("Select annotations to show","Here you can include IDs and/or symbols 
											of the annotated genes"),
		tags$br(),
		Warning("A genome assembly must be loaded to get the annotations. To load a genome assembly, 
				use 'Assembly' section. Once an assembly has been loaded, every genomic range of a 
				ROI is annotated to the nearest gene"),
		tags$br(),
		tags$br(),
		Field("Enrichments","These are the signal enrichments associated to the ROI. It shows the sum 
				of the pileup of the reads for the selected enrichment within each genomic range of the ROI"),
		tags$br(),
		Field("Get annotated genes within a genomic window","It retrieves the list of genes (IDs and/or symbols) 
					found at a user-defined distance from the genomic ranges of the ROI. User defined 
					distance is set in '<b>Select the genomic window (bp)</b>'"),
		Field("Select the genomic window (bp)","Select the size of the genomic window flanking the genomic ranges 
					of a ROI. The width is measured from  the center of the genomic ranges")

	)
)
#Preview
msg_getRois_preview<-list(
	title="Preview of the ROI table",
	text="Lets you preview the matrix file of a ROI with all the information selected"
)





#### Associate enrichments

#ROI selection
msg_associateEnrichments_associateRemove<-list(
	title="Associate or remove enrichments to ROI(s)",
	text=list(
		"Associate one or more enrichments (BAM/WIG) to one or more ROIs",
		tags$br(),
		tags$br(),
		Field("Choose ROI(s)","Choose one or more ROIs to which you want to associate or remove enrichments"),

		tags$br(),
		tags$h3("Associate enrichment(s) to ROI(s)"),
		"Select one or more enrichment(s) to associate to the selected ROI(s).",
		tags$br(),
		tags$br(),
		Field("Select enrichment(s) to associate to selected ROI(s)","Select one or more enrichments (BAM/WIG) 
				to associate to the selected ROI(s). This operation will calculate the pileup of the reads 
				from the enrichment files (either BAM or WIG) selected  for each base of each range of each ROI selected"),
		tags$br(),
		Tip("To import enrichment files go to 'Enrichment files' in the 'Import data' section"),
		tags$br(),
		tags$br(),
		Warning("This is a computationally intensive task (it may require from several seconds to few minutes)"),
		tags$br(),
		tags$br(),
		Field("Number of cores","Choose the number of cores to use for enrichment association"),
		tags$br(),
		Comment("The more cores you use, the faster the operation will be, and more memory will be required"),
		tags$br(),
		tags$br(),
		Tip("Increase the number of cores to speed up this operation"),
		tags$br(),
		tags$br(),
		Field("Normalization method","Select how to normalize the reads:
				<li>library-size normalization: the reads will be normalized by the size of their library</li>
				<li>custom normalizer: the reads will be normalized by the library size of a user defined enrichment file present in the list</li>"),
						
		
		tags$p(
			tags$b("Enrichment files (BAM/WIG) to associate to selected ROI(s):"),
			tags$br(),
			"Select one or more enrichments to associate to selected ROI(s). To import those files go to 'Enrichment files'
			in the 'Import data' section"
		),		
		"This operation will calculate the pileup of the reads for each base of each range of each ROI selected.
 		The information of the reads come from the enrichment files (either BAM or WIG) selected.
 		This is a computationally intensive task, it may require time (from several seconds to few minutes, depending on the size of the ROIs, 
 		of the enrichment files and on the number of ROIs)",
		tags$br(),
		tags$br(),
		tags$p(
			tags$b("Number of cores"),
			tags$br(),
			"Choose the number of cores to use for enrichment association. The more cores you use, 
			the faster the operation, but the more memory it will require."
		),



		tags$br(),
		tags$br(),
		tags$br(),
		tags$h3("Remove enrichments from ROI(s)"),
		"Select the enrichments to be removed from the selected ROI(s). 
		The list  comprises all the signal enrichment loaded and associated to ROIs in the ChroKit session loaded",
		tags$br(),
		tags$br(),
		Tip("This operation will free some memory"),
		tags$br(),
		tags$br(),
		Warning("If an enrichment file is eliminated, its association to any of 
				the selected ROIs in the ongoing session is eliminated as well")

	)
)

#### Rename/Order enrichments

#Rename enrichments
msg_renameEnrichments_renameOrderEnrichments<-list(
	title="Change the name or reorder associated enrichments",
	text=list(
		Field("Select ROI","Select the ROI in which enrichments must be renamed or reordered"),
		tags$br(),
		
		tags$h3("Rename enrichments"),
		"Rename enrichments associated to ROIs",
		tags$br(),
		tags$br(),
		Field("Select enrichment to rename","Choose the enrichment associated to the selected ROI to be renamed"),
		Field("New enrichment name","Choose the new name for the enrichment selected"),
		tags$br(),
		Warning("The name of the enrichment will change in all the other ROIs to which it had been associated"),
		tags$br(),
		tags$br(),
		tags$br(),


		tags$h3("Reorder enrichments"),
		"Change the order in which the enrichments appear in the ROI. For each of the enrichments associated to the selected ROI, 
			set the number corresponding to the new position in the ranking, and finally press 'Reorder!' button",

		tags$br(),
		tags$br(),
		Comment("1. This is useful when you want to change the order of the enrichments in which they appear in plots. <br>
				2. The order of the enrichments will be changed only in the selected ROI",plural=T)

	)
)


#### GO analyses

#parameters
msg_goAnalysis_parameters<-list(
	title="Parameters for gene ontology",
	text=list(
		"Allows to set the parameters used for the hypergeometric test on genes provided by the user",
		tags$br(),
		tags$br(),
		Comment("The gene signatures are lists of genes in gmt format present in the directory appContent/signatures/ of the program. 
				Gene signatures can be loaded from the Molecular Signature Database (MSigDB) or similar resources"),
		tags$br(),
		tags$br(),
		Warning("Gene symbols in signatures must be in uppercase"),
		tags$br(),
		tags$br(),		

		tags$h3("Variables"),
		tags$br(),
		Field("Select the source","Allows you to choose whether genes are 
				provided as a gene list (loaded by the user) or as genes associated to ROIs.
				The menu is context-dependent:"),
		HTML("<li>From ROI: this option will let you select the ROIs from which annotated 
			genes will be extracted. If multiple ROIs are selected, the result will be displayed 
			as a heatmap: ROIs selected in rows and the significant signatures in columns (the colour 
			will represent the statistical significance). If a single ROI is selected, the result 
			will be displayed as a barplot of the p adjusted for each significant signature 
			
			</li>"),
		tags$br(),
		
		Warning("A genome assembly must be loaded. To load a genome assembly, use 'Assembly' section"),
		tags$br(),
		tags$br(),
		HTML("<li>From gene list: this option will let you input a list of gene symbols (one per line)
			. If genes are lowercase (murine), they will be converted in uppercase (Human) to be consistent 
			with signature lists available in the program. The Output is a barplot representing the p adjusted</li>"),
		tags$br(),
		
		Warning("The symbols given as input will be converted in uppercase (“humanized”) to match 	those within 
				the gmt signatures"),
		tags$br(),
		tags$br(),
		Field("Select signature(s)","Select the gene signatures(s) to use for the gene ontology analysis. The 
			background for the hypergeometric test (the 'universe')  will be the union of all the genes of 
			the signatures selected"),
		Field("How to order results","Choose the order of the results of GO analysis. If a single ROI has 
			been selected, or a custom gene list has been used as input, the results will be ranked based on 
			-log10 p adjusted by default. If multiple ROIs have been selected, the user can choose whether 
			to rank the heatmap based on the best p adjusted 
			or whether to rank based on clustering"),
		Field("Cluster type","Choose the clustering method for the heatmap. This is based on the p-adjusted of the signatures"),
		Field("Min/Max signature size","Set the minimum and maximum size of the gene signatures used for the ontology analysis. 
			All signatures with a size below the minimum or above the maximum will not be considered"),
		tags$br(),
		tags$br(),		
		tags$br(),
		tags$h3("Filtering"),
		tags$br(),
		Field("Quantile threshold for padj color scale","Quantile to use to set the maximum for the color scale of the p adjusted"),
		tags$br(),
		Tip("The lower the value, the more saturated the colors will be"),
		tags$br(),		
		tags$br(),
		Field("Gene ratio threshold","Filter the results using a threshold on the gene ratio (genes associated to the 
			ROI or in the custom gene list / all the genes in the signature)"),
		Field("-log10 padj threshold","Filter the results keeping the hits with a significance (-log10 p adjusted) above a certain threshold"),
		Field("Top significant hits","Filter the results only for the top n significant hits (n=number of hit, this value is user-defined)"),
		Field("Choose a color scale","Choose the color scale for the heatmap")		

	)
)

#GO plot
msg_goAnalysis_goPlot<-list(
	title="Gene ontology plot",
	text=list(
		"It displays the plot of the GO analysis (heatmap or barplot). 
		For each signature, genes annotated to a specific ROI can be displayed by clicking the corresponding cell",
		tags$br(),		
		tags$br(),
		Tip("To optimize the visualization of the results, change parameters such as the color scale")		
	)
)

#GO table
msg_goAnalysis_goTable<-list(
	title="Gene ontology table",
	text="GO results are displayed in a tabular format and can be explored interactively.
		You can download the tab-delimited text file for downstream processing"
)



#### Predefined ROI preparation for heatmaps
msg_PredefPipeline_parameters<-list(
	title="Prepare a ROI for heatmap",
	text=list(
		"This pipeline is dedicated to the processing of ROIs in a format optimal for heatmaps
		generated in the section 'Genomics'",
		tags$br(),
		tags$br(),
		Field("Select ROI for preparation","Select the ROI to be modified for heatmap visualization"),
		Field("% genomic ranges to keep","Select a random subset of genomic ranges"),
		tags$br(),
		Tip("Subsetting is suggested when dealing with large ROIs (i.e. >30/40 000 genomic ranges)"),
		tags$br(),
		tags$br(),
		Field("Do you want to center on summit?","Choose “Yes” if you want to center each genomic 
					range of the ROI to the summit of a specific enrichment"),
		tags$br(),
		Tip("If “Yes” is chosen, a menu will appear, allowing to select the enrichment to use for summit detection"),
		tags$br(),
		tags$br(),
		Field("Upstream/Downstream","Define the size of each genomic range by indicating the number of upstream/downstream 
				base pairs from its the center. The center is user-defined (midpoint or summit)"),
		Field("Select enrichment to associate to the new ROI","Select the enrichments to associate to the ROI"),
		Field("Number of cores","Choose the number of cores to use for enrichment association"),
		tags$br(),
		Comment("The more cores you use, the faster the operation will be, and more memory will be required"),
		tags$br(),
		tags$br(),
		Tip("If you are running ChroKit in a powerful machine, increase the number of cores to speed up this operation"),
		tags$br(),
		tags$br(),
		Field("New ROI name","Type the name of the new ROI here")										
	)
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
		"Parameters for the genomic analysis of a single ROI",
		tags$br(),
		tags$br(),
		Warning("A genome assembly must be loaded. To load a genome assembly, go to the 'Assembly' section"),
		tags$br(),
		tags$br(),
		Field("Select ROI","Choose the ROI"),
		Field("Select enrichment","If the selected ROI has enrichments 
				associated, choose here the enrichment you want to use for the analysis"),
		tags$br(),
		Comment("If present, the enrichment in each annotated genomic element (promoters, 
				genebodies and intergenic regions) will be shown for all the genomic ranges of the ROI"),
		tags$br(),
		tags$br(),
		Field("Choose normalization","Select what kind of signal to show in the plots:<br>
				<li><b>Total reads (rpm)</b>: The number of library-normalized reads</li>
				<li><b>Read density (rpm/bp)</b>: The number of library-normalized reads, 
					normalized also by the length of the genomic ranges</li>")	,
		Field("Choose color palette","Choose the color palette to show in the plots")	

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
		Warning("To avoid outlayers, only 0-95% of the distribution is shown")		
	)
)
#Enrichment boxplot
msg_singleEvaluation_enrichmentBoxplot<-list(
	title="Enrichment in annotated genomic elements",
	text=list(
		"If an enrichment associated to the selected ROI has been selected, it displays a boxplot 
		of the reads (or the read density) stratified by the annotated regions (promoters, genebodies, intergenic regions)",
		tags$br(),
		tags$br(),
		Warning("The enrichment file must be associated to the ROI. To associate enrichments to a ROI, 
			see 'Associate enrichments' tab in 'ROI management' section")		
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
		Warning("The enrichment file must be associated to the selected ROI. To associate enrichments, 
			see 'Associate enrichments' tab in 'ROI management' section")		
	)
)

#### Pairwise overlaps

#Parameters
msg_pairwiseOverlaps_parameters<-list(
	title="Parameters for pairwise overlaps",
	text=list(
		Field("Select ROI-1","Choose the first ROI for the overlap"),
		Field("Select ROI-2","Choose the second ROI for the overlap"),
		tags$br(),
		Warning("This analysis can be performed only if two ROIs are selected"),
		tags$br(),
		tags$br(),
		Field("Minimum number of bp for overlap","Set the minimum number of bases 
					used to define the overlap between genomic ranges of the two ROIs"),
		Field("Choose enrichment-1","Choose one of the enrichments associated to ROI-1 to 
				evaluate the signal of reads in the subset of genomic ranges of ROI-1 alone 
				or those overlapping with ROI-2"),
		Field("Choose enrichment-2","Choose one of the enrichments associated to ROI-2 to evaluate 
				the signal of reads in the subset of genomic ranges of ROI-2 alone or those overlapping 
				with ROI-1"),
		tags$br(),
		Warning("The possibility to select the enrichments will appear only if enrichments are associated 
				to ROI-1, ROI-2 or both. If enrichments are not associated to either ROI-1 or ROI-2, only 
				the fraction of overlapping/not overlapping genomic ranges will be displayed."),
		tags$br(),
		tags$br(),
		Field("log2","Used for log2 transformation of the signals of the enrichments"),
		tags$br(),
		Tip("Check 'log2' for a better visualization in the scatterplot"),
		tags$br(),
		tags$br(),
		Field("Choose normalization","Choose what value is shown in the plots:<br>
				<li><b>Total reads (rpm)</b>: The number of library-normalized reads</li>
				<li><b>Read density (rpm/bp)</b>: The number of library-normalized reads, 
					normalized also by the length of the genomic ranges</li>"),
		Field("Choose color palette","Choose the color palette used in the plots"),
		Field("In scatterplot, show:","This menu will appear only if the scatterplot is present. Choose which are the genomic ranges to show in the scatterplot: <br>
				<li>exclusive ROI-1 ranges: genomic ranges present only in ROI-1</li>
				<li>exclusive ROI-2 ranges: genomic ranges present only in ROI-2</li>
				<li>common ranges: genomic ranges common to ROI-1 and ROI-2</li>")			
					
	)
)

#Overlap
msg_pairwiseOverlaps_overlap<-list(
	title="Displays the overlap of the two ROIs",
	text=list(
		"he overlap of the two ROIs can be displayed as a barplot or with a Venn diagram",
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
#box/scatter/Calibration
msg_pairwiseOverlaps_overlapAndEnrichment<-list(
	title="Enrichments in overlaps",
	text=list(
		"A series of analysis tool to evaluate the signal enrichments in genomic ranges shared by two ROIs (ROI-1 and ROI-2)",
		tags$br(),
		tags$br(),
		Field("Boxplot","It shows the signal of the selected enrichment (i.e. enrichment-1) in the regions of ROI-1 which 
				either overlap or not with ROI-2.
				If the enrichment is associated also to ROI-2, it shows also enrichment-1 signals in the genomic ranges 
				of ROI-2 that are not overlapping with ROI-1. 
				The same is true for ROI-2 and enrichment-2."),
		Field("Scatterplot","If enrichment-1 and enrichment-2 are associated to both ROI-1 and ROI-2, it shows the 
				correlation between the two enrichments. Overlapping and exclusive (ROI-1 alone, ROI-2 alone) genomic 
				ranges  are highlighted"),
		tags$br(),
		Warning("This works only if both enrichment-1 and enrichment-2 are associated with both ROI-1 and ROI-2"),
		tags$br(),
		tags$br(),
		Tip("For a better representation of this plot, select 'log2' in the parameters"),
		tags$br(),
		tags$br(),				
		Field("Calibration","It shows the fraction of the overlapping regions between ROI-1 and ROI-2 at increasing values 
			of enrichments thresholds. The threshold is applied either to ROI-1 or ROI-2 or to both ROIs together (combined)"),
		tags$br(),
		Comment("This plot recalculates the overlap between ROI-1 and ROI-2 if we consider only those genomic ranges 
				with an enrichment above an increasing threshold. 'combined' means that the enrichment threshold 
				is applied to both ROI-1 and ROI-2 simultaneously; otherwise, the threshold is applied only to 
				one of the two ROIs, thus preserving all the genomic ranges of the other ROI for the analysis.<br><br>
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
		"Set the parameters for a heatmap showing overlaps between multiple ROIs",
		tags$br(),
		tags$br(),		
		tags$h3("Variables"),
		tags$br(),
		Field("Master ROI(s)","One or more ROIs can be selected to create the master ROI"),
		tags$br(),
		Warning("1. If more than one master ROI is selected, it won't be possible to extract a 
				new ROI from the single clusters of the heatmap.<br> 
				2. The jaccard index matrix won't 
				be shown, as well.",plural=T),
		tags$br(),
		tags$br(),
		Field("ROIs to view","Select one or more ROIs for which you want to test the overlaps within the master ROI(s)"),
		Field("ROIs ordering","You can set the order by which ROIs are displayed in the heatmap"),
		Field("ROIs for cluster","Select which is/are the ROI(s) that will drive the clustering of the rows 
			(genomic ranges) of the heatmap. The cluster will be based on the information of the overlap (0 or 1) 
			in all the bins of each genomic range"),
		Field("Clustering type","The clustering type menu will appear if at least one ROI has been selected for the clustering. 
				You can select two kind of clustering algorithms:<br>
				<li><b>K-means</b>: performs a K-means clustering on the rows of the heatmap, using the information of 
					the overlap of all ROIs selected. This allows the identification of specific patterns of overlap 
					of the different ROIs. User-defined K-means parameters are: number of clusters, 
					number of starting points and number of iterations</li>
				<li><b>Hierarchical</b>: performs a hierarchical clustering on the rows of the heatmap, using the 
					information of the overlap of all ROIs selected. User-defined parameters for hierarchical clustering 
					are: number of clusters, distance method and clustering method</li>"),
		tags$br(),
		Comment("The fields for advanced parameters of different clustering algorithms will appear depending on the clustering type selected"),
		tags$br(),
		tags$br(),
		Field("Number of bins","Set the number of bins each genomic range will be divided into. More bins means more resolution in determining 
			the positions of overlaps in each genomic range of the master ROIs"),
		tags$br(),
		Warning("The plot 'Overlap frequency bias' will be shown only if the number of bins set is > 2"),
		tags$br(),
		tags$br(),
		Tip("When you want to extract the most common patterns of overlap between a large number of ROIs, 
			set the number of bins =1, a k-mean clustering driven by all the shown ROIs and a number of clusters 
			k=2^n, where n=number of ROIs. In this way, all the possible combinations of overlaps between ROIs can 
			be displayed at once"),
		tags$br(),
		tags$br(),	
		Field("Strand-specific overlaps","If the master ROI has strand information and this option is set, only the 
			overlaps in the same strand of each genomic range of the master ROI will be considered"),
		tags$br(),
		tags$br(),
		tags$br(),		

		tags$h3("Advanced"),
		tags$br(),
		Field("Random sample of genomic ranges to show","Reduce the number of genomic ranges of selected 
				ROI(s) to a user-defined number (default: 2000 ranges)"),
		tags$br(),
		Tip("1. Usually, a random sample of 2000 genomic ranges is sufficient to have a good representation of 
			overlaps in the heatmap.<br>
			2. If the aim of the  analysis is to identify all the regions of the genome that share a defined combination 
			of overlaps, it is important to include all the genomic ranges of the master ROI and not a random subset",plural=T),
		tags$br(),
		tags$br(),
		Warning("Increasing this number will result in a slower calculation of the heatmap"),
		tags$br(),
		tags$br(),
		Field("Select colors","Choose the colors for showing the overlapping bins:<br>
				<li><b>global color</b>: set a unique color for all ROIs selected</li>
				<li><b>custom colors</b>: set a custom color for each ROI selected</li>"),
		Field("positional overlap %","If checked, the plot of positional distribution of 
				overlaps will show the fraction of overlapping events, instead of the absolute number")				

	)
)

#Heatmap
msg_digitalHeatmap_heatmap<-list(
	title="Position-based heatmap",
	text=list(
		"This plots an interactive heatmap that shows, for each genomic range of the master ROI selected (each row) 
		the overlap of the other ROI(s). Note that overlaps are calculated for each bin of the genomic ranges. Overlapping regions are displayed as colored. 
		The fraction of the overlapping genomic range for each ROI compared to the total ranges of the master 
		ROI is displayed on the x axis. Clusters are shown on the left side of the heatmap, as colored bars.",
		tags$br(),
		tags$br(),		
		"If a cluster bar is clicked, the number of genomic ranges belonging to the selected clusters are displayed, 
		along, with the cluster number. The field “New ROI from selection” will also appear, giving you the possibility 
		to create a new ROI.",
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
			2. You can choose to show the fraction of overlapping ranges instead of the absolute number 
			by checking '<b>positional overlap %</b>' option in parameters. This is useful when comparing multiple master 
			ROIs (it normalizes for the number of genomic ranges of each master ROI)",plural=T)

	)
)


#### Analogic heatmap

#Parameters
msg_analogicHeatmap_parameters<-list(
	title="Parameters for enrichment-based heatmap",
	text=list(
		"Set the parameters for the heatmap",
		tags$br(),
		tags$br(),
		tags$h3("Variables"),
		tags$br(),
		Field("Select ROI(s)","Select one or more ROIs: the corresponding genomic ranges will be the 'rows' of the heatmap"),
		tags$br(),
		Warning("To plot the heatmap, enrichments must be associated to the ROI(s) selected"),
		tags$br(),
		tags$br(),
		Field("Select enrichments to show","Select one or more enrichments to show in the heatmap"),
		tags$br(),
		Warning("An enrichment appears as available option if only the enrichment is associated to ALL the ROIs selected"),
		tags$br(),
		tags$br(),
		Tip("To associate enrichments to ROIs, go to 'Associate enrichments' tab in 'ROI management’ section"),
		tags$br(),	
		tags$br(),
		Field("Enrichment order","If multiple enrichments have been selected, you can set their order of appearance in the heatmap"),
		Field("Clustering/ranking","Select how to organize the heatmap rows. You can rank the rows by one of the selected enrichments, 
				or you can cluster the rows. You can select two kind of clustering algorithms:<br>
				<li><b>K-means</b>: performs a K-means clustering on the rows of the heatmap, using the values of the enrichments 
					selected for clustering. This allows the identification of specific patterns of enrichments in the different ROIs. 
					User-defined K-means parameters are: number of clusters, number of starting points and number of iterations</li>
				<li><b>Hierarchical</b>: performs a hierarchical clustering on the rows of the heatmap, using the information of 
					the enrichments selected for clustering. User-defined parameters for hierarchical clustering are: number of 
					clusters, distance method and clustering method</li>"),
		tags$br(),
		Comment("The fields for advanced parameters of different clustering algorithms will appear depending on the clustering type selected"),
		tags$br(),
		tags$br(),
		Field("Number of bins","Set the number of bins each enrichment will be divided into. Increasing the number of bins will 
				increase the resolution of the heatmap"),
		tags$br(),
		tags$br(),
		tags$br(),
		
		tags$h3("Advanced"),
		tags$br(),
		Field("Random sample of genomic ranges to show","It reduces the number of genomic ranges of selected ROI(s) to a 
					user-defined number (default: 2000 ranges)"),
		tags$br(),
		Warning("Increasing this number will result in a slower calculation of the heatmap"),
		tags$br(),
		tags$br(),		
		Tip("Usually, a random sample of 2000 genomic ranges is sufficient to have a good representation of 
			overlaps in the heatmap.<br><br>
			If the aim is to extract new ROI(s) from a selected area or cluster, it is important to include all 
			genomic ranges of the master ROI and not a random subset"),
		tags$br(),
		tags$br(),
		Field("log2","If selected, the signals shown as profiles and boxplots will be log2 transformed"),
		Field("Quantile threshold","Used to set the maximum quantile value for the color scale of the 
				signals. Two different kind of color scales can be set:<br>
				<li><b>Uniform</b>: the threshold for the color scale saturation will be a quantile of all 
					enrichments considered together</li>
				<li><b>Individual</b>: used to set a different threshold of the color scale 
					saturation for each of the signal enrichments displayed</li>"),
		tags$br(),
		Comment("Use 'Uniform' option if you have to compare signals of different enrichments. 
				Use 'Individual' if you include different enrichments in the same heatmap with values 
				that have different order of magnitudes)"),
		tags$br(),
		tags$br(),
		Tip("The lower the value of the threshold, the more saturated the colors will be"),	
		tags$br(),
		tags$br(),
		Field("Select colors","Choose the colors of the enrichments shown:<br>
				<li><b>default color</b>: set a unique color scale for all enrichment</li>
				<li><b>custom colors</b>: set a custom color scale for each enrichment selected</li>"),
		Field("Group colors (boxes)","If checked, group the colors in the boxplots of the selected area 
				of the heatmap according to the stratification (ROIs, clusters or enrichments)")	


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
		Warning("A single ROI must be selected in order to extract a new ROI from the selected area of the Heatmap 
				or to have the barplot of the clusters"),
		tags$br(),
		tags$br(),

		Tip("To extract new ROIs from clusters or selected areas of the heatmap, make sure to use all the genomic 
			ranges of the ROI, and not a random sample (this is done in '<b>Random sample of genomic ranges to show</b>' 
			(change in GUI) by inputing a value equal to or greater than the number of genomic ranges in the ROI selected)")
		

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
		"Plots the enrichment of the selected area of the heatmap, or the selected cluster.<br> 
			Plotting options are:",
		tags$br(),
		tags$br(),
		Field("Boxplot by ROI/cluster","The enrichments are shown as boxplot stratified by the ROIs 
			provided. In the case of a clustered heatmap from a single ROI, the boxplot will be stratified by clustering"),
		Field("Boxplot by enrichment","The enrichments are shown in different ROIs as boxplot stratified by enrichment"),
		tags$br(),
		Tip("To highlight the stratification made be the boxplots (ROIs/clusters or enrichments), 
				check the '<b>Group colors (boxes)</b>' field in parameters"),
		tags$br(),
		tags$br(),
		Field("cor","Heatmap showing the pairwise Pearson correlation between the enrichments in a given ROI"),
		Field("pcor","Heatmap showing the Pearson pairwise partial correlation between the enrichments in a given ROI"),
		tags$br(),
		Comment("The partial correlation between enrichment A and B will be computed using, as covariates, 
				all the other enrichments selected"),
		tags$br(),
		tags$br(),
		Warning("1. The correlation and partial correlation heatmaps are available only if a single ROI is selected. <br>
				2. The correlation heatmap is available only if at least 2 enrichments are selected. <br>
				3. The partial correlation heatmap is available only if at least 3 enrichments are selected.",plural=T)		

	)
)

#### Enrichment in ROIs

#Parameters
msg_enrichmentInRois_parameters<-list(
	title="Parameters for visualizing enrichments in ROIs",
	text=list(
		"This section is dedicated to the analysis of the enrichments in 
		the ROIs: signal distribution (profiles), box-plots, signal correlations matrix and scatterplot",
		tags$br(),
		tags$br(),
		Field("Select ROI(s)","Choose one or more ROI(s) for the analysis"),
		Field("Select enrichments to show","Choose the enrichments to be shown associated to the selected ROI(s)"),
		tags$br(),
		Warning("An enrichment will appear as available option only if the enrichment is associated to ALL the ROIs selected"),
		tags$br(),
		tags$br(),
		Tip("To associate enrichments to ROIs, go to 'Associate enrichments' tab in 'ROI management' section")	,
		tags$br(),
		tags$br(),			
		Field("Number of bins","Set the number of bins each enrichment will be divided into, 
				in order to generate the enrichment profile plot"),
		tags$br(),
		Tip("Increase this number to improve the resolution of the profiles"),
		tags$br(),
		tags$br(),		
		Field("Normalization method for the enrichments","Select how to normalize the signals of the enrichments:<br>
				<li><b>Total reads (rpm)</b>: The number of library-normalized reads</li>
				<li><b>Read density (rpm/bp)</b>: The number of library-normalized reads, 
					normalized also by the length of the genomic ranges</li>"),
		Field("log2","If checked, the signals of enrichments will be log2 transformed"),
		Field("Type of correlation","Select which correlation will be shown in 
				correlation/partial correlation heatmaps (Spearman or Pearson)"),
		Field("Group colors (boxes)","Select the colors of boxplots according to the ROI or the enrichment")

	)
)

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
				by clicking 'download data'"),
		tags$br(),
		tags$br(),		
		Tip("To color according to the grouping of the boxplots (ROIs or enrichments), tick the 
			'Group colors (boxes)' field in 'parameters' to color according to the grouping of the 
			boxplots (ROIs or enrichments), tick the 'Group colors (boxes)' field in the parameters")	
	)
)
#Correlations
msg_enrichmentInRois_correlations<-list(
	title="Correlation heatmaps",
	text=list(
		"Shows the correlation or the partial correlation between enrichments in the selected ROI.<br>
		Clicking a square in the matrix will show the corresponding scatterplot of the pairwise 
		(partial)correlation in the scatterplot on the right",
		tags$br(),
		tags$br(),
		Field("Cor-Heatmap","Correlation Heatmap of the enrichments in a given ROI"),
		Field("Pcor-Heatmap","Partial correlation heatmap of the enrichments in a given ROI"),
		tags$br(),
		Comment("The partial correlation between enrichment A and B will be computed using, 
					as covariates, all the other enrichments selected"),
		tags$br(),
		tags$br(),
		Warning("1. The correlation and partial correlation heatmaps will appear only when a single ROI 
				has been selected.<br> 
				2. The correlation heatmap will appear only when at least 2 enrichment 
				were selected.<br> 
				3. Partial correlation heatmap only when at least 3 enrichments were selected",plural=T),
		tags$br(),
		tags$br(),				
		Tip("To change the correlation type (Spearman or Pearson) use the option in the parameters")

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
		"Here you can define the parameters for the metagene plot",
		tags$br(),
		tags$br(),
		Field("Select the gene list","Choose the gene lists in which you want to perform the analysis"),
		tags$br(),
		Tip("To import a gene list, go to 'ROI' section"),
		tags$br(),
		tags$br(),
		Field("ROIs associated to selected gene lists","This panel shows all the ROIs that constitute 
			the selected gene list(s) (promoters, transcript, TES of the genelist(s) selected)"),
		Field("Select enrichments to show","Choose the enrichments you want to analyze. These 
			enrichments must be associated both with promoters, transcripts and TES ROIs of the gene list"),
		tags$br(),
		Tip("To associate enrichments to a gene list, use the 'Associate enrichments' tab in 'ROI management' section"),
		tags$br(),
		tags$br(),
		Field("Number of bins","Set the number of bins for the gene body portion of the metagene plot"),
		tags$br(),
		Tip("To improve the resolution of the metagene profile, increase the number of bins"),
		tags$br(),
		tags$br(),
		Field("Mean or median?","Choose whether to show the mean or the median of the signal in the metagene profile"),
		Field("log2","Choose whether to transform the signals of the enrichments using log2"),
		Field("Choose normalization","Choose the normalization method for the enrichments:<br>
		 		<li><b>Total reads (rpm)</b>: The number of library-normalized reads</li>
				<li><b>Read density (rpm/bp)</b>: The number of library-normalized reads, normalized also by the length of the genomic ranges</li>")				

	)

)
	
#Profiles
msg_dynamicsOnGenes_profiles<-list(
	title="Metagene profile",
	text=list(
		"Shows the metagene profile of enrichments in the selected gene list(s).
		This plot is particularly useful to evaluate RNApol2 distribution along genes",
		tags$br(),
		tags$br(),
		Comment("The metagene plot starts at –30% of transcripts length from TSS and ends at +30% from TES")		
	)
)
#TSS enrichment
msg_dynamicsOnGenes_TSSenrichment<-list(
	title="Enrichments at promoters",
	text=list(
		"Shows a boxplot of the enrichments at the promoters of the selected gene list(s).",
		tags$br(),
		tags$br(),
		Comment("The upstream and downstream limits from TSS to define promoters are set the 
				first time the genome assembly is loaded"),
		tags$br(),
		tags$br(),	
		Tip("To change the size of intervals around TSS  or TES , go to the 'Assembly' section and  
			change the default up/downstream values for the definition of TSS and TES intervals, 
			then reload the gene list in 'ROI' section")		
	)
)
#genebodies enrichment
msg_dynamicsOnGenes_genebodiesEnrichment<-list(
	title="Enrichments at genebodies",
	text=list(
		"Shows a boxplot of the enrichments at genebodies of the selected gene list(s)",
		tags$br(),
		tags$br(),
		Comment("A genebody is defined as the annotated gene, minus the TSS and the TES intervals"),
		tags$br(),
		tags$br(),
		Tip("To change the size of intervals around TSS  or TES , go to the 'Assembly' section and  
			change the default up/downstream values for the definition of TSS and TES intervals, 
			then reload the gene list in 'ROI' section")				
	)
)
#TES enrichment
msg_dynamicsOnGenes_TESenrichment<-list(
	title="Enrichments at TES",
	text=list(
		"Shows a boxplot of the enrichments at the TES of the selected gene list(s)",
		tags$br(),
		tags$br(),
		Comment("The upstream and downstream limits from TES are defined the first time the genome assembly is loaded"),
		tags$br(),
		tags$br(),
		Tip("To change the size of intervals around TSS  or TES , go to the 'Assembly' section and  
			change the default up/downstream values for the definition of TSS and TES intervals, 
			then reload the gene list in 'ROI' section")		

	)
)
#TSS ranked enrichments
msg_dynamicsOnGenes_TSSranked<-list(
	title="Promoters ranked enrichments",
	text=list(
		"For each gene in the gene list(s) selected, it shows the ranked enrichments at promoters",
		tags$br(),
		tags$br(),
		Comment("The upstream and downstream limits from TSS were defined the first time the genome assembly was loaded"),
		tags$br(),
		tags$br(),		
		Tip("To change the size of intervals around TSS  or TES , go to the 'Assembly' section and  
			change the default up/downstream values for the definition of TSS and TES intervals, 
			then reload the gene list in 'ROI' section")			
	)
)

#Genebodies ranked enrichments
msg_dynamicsOnGenes_genebodiesRanked<-list(
	title="Genebodies ranked enrichments",
	text=list(
		"For each gene in the gene list(s) selected, it shows the ranked enrichments at genebodies",
		tags$br(),
		tags$br(),
		Comment("A genebody is defined as the annotated gene, minus the TSS and the TES intervals"),
		tags$br(),
		tags$br(),
		Tip("To change the size of intervals around TSS  or TES , go to the 'Assembly' section and  
			change the default up/downstream values for the definition of TSS and TES intervals, 
			then reload the gene list in 'ROI' section")				
	)
)
#Stalling Index
msg_dynamicsOnGenes_stallingIndexRanked<-list(
	title="Stalling index (ranked)",
	text=list(
		"It calculates the stalling index (= number of reads at TSS / number of reads at genebodies) and plot each gene in ascending order",
		tags$br(),
		tags$br(),
		Comment("This is typically used to calculate the Pol2 stalling index")		
	)
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





































