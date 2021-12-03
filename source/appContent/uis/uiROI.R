tabROI<-tabItem (tabName = "ROIblock",
  tabsetPanel(type="pills",
    tabPanel("Overlaps",
      fluidRow(
        box(width=3,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_newRois_options",title="Options"),
          #minimum overlap (true for all the blocks)
          HTML("<b>Minimum number of bp to consider for overlaps:</b>"),
          numericInput(inputId = 'minOverlapNEWROI',label=NULL,min = 1, step = 5,value=1),
          checkboxInput("StrandSpecOverlapNEWROI", label="Strand-specific overlaps",value = FALSE, width = NULL),

          textInput("ROIname",label="Name of the ROI",placeholder="type new ROI name here"),
          # resize width from the center (work also on the summit, when BAM is available)
          #button for "create ROI" (check if some BED files are selected, otherwise msg in logs: not possible)
          actionButton("maketheROI", "Build the ROI!")            
        ), 
             
        box(width=9,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title=boxHelp(ID="msg_newRois_ROIcombination",title="Combination of ROIs"),
          fluidRow(
            column(width=4,
              HTML("<h4><b>Choose the reference ROI:</b></h4>"),
              HTML("Select ROI(s):"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
                checkboxGroupInput("selectedROIs",NULL,NULL)
              ),
              # radiobutton for union/intersection
              radioButtons("choiceROI",label="Criteria for building the aggregated reference ROI:" ,
                                        choices=c("Intersection"="intersection",
                                                  "Union"="union"),
                                        selected="union"
                          )            
            ),
            column(width=4,
              HTML("<h4><b>...that overlaps with contrast ROI:</b></h4>"),
              HTML("Select ROI(s):"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
                checkboxGroupInput("overlapROIs",NULL,NULL)
              ),
              radioButtons("choiceoverlapROI",label="Criteria for building the aggregated contrast ROI:" ,
                                      choices=c("Intersection of the contrast ROIs"="stringent",
                                                #"Overlap of the contrast ROIs"="permissive",
                                                "Union of the contrast ROIs"="allofthem"),
                                      selected="allofthem"
              )            
            
            ),
            column(width=4,
              HTML("<h4><b>...that doesn't overlap with contrast ROI:</b></h4>"),
              HTML("Select ROI(s):"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 400px; max-width: 300px; background-color: #ffffff;",
                checkboxGroupInput("notoverlapROIs",NULL,NULL)
              ),
              radioButtons("choicenotoverlapROI",label="Criteria for building the aggregated contrast ROI:" ,
                                      choices=c("Intersection of the contrast ROIs"="intersection",
                                                "Union of the contrast ROIs"="union"),
                                      selected="union"
              )            
            
            )
          )  

        )
 
      )
    ),
    


    # RESIZE ROIs
    tabPanel("Modify ROIs",

      fluidRow(
        box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_modifyRois_resize",title="Resize"),
          selectInput("selectROItoresize",label="Select ROI to resize:",NULL),

          #other options if needed

          fluidRow(
            column(width=4,HTML("<b>Upstream</b>")),
            column(width=4),
            column(width=4,HTML("<b>Downstream</b>"))
          ),
          fluidRow(
            column(width=4,
              numericInput(inputId = 'sliderUpstreamROI',label=NULL,min = 0, max = 20000, step = 50,value=1000) 
            ),
            column(width=4,
              HTML("<h4>Range center</h4>")
            ),
            column(width=4,
              numericInput(inputId = 'sliderDownstreamROI',label=NULL,min = 0, max = 20000, step = 50,value=1000) 
            )
          ),
          textInput("ROInameResize",label="Name of the ROI",placeholder="type new ROI name here"),
          actionButton("resizeROI","Create ROI")
        ),

        box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_modifyRois_summit",title="Center on summit"),

          selectInput("selectROItoCenterSummit",label="Select ROI to center on summit:",NULL),
          selectInput("selectBAMtoCenterSummit",label="Select enrichment to use for summit:",NULL),
          textInput("ROInameSummit",label="Name of the ROI",placeholder="type new ROI name here"),
          actionButton("SummitROI","Create ROI")
        ),

        box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_modifyRois_sample",title="Random sample"),

          selectInput("selectROItoSample",label="Select ROI to sample:",NULL),
          sliderInput('quantileThreshSample',label="Fraction to keep:",min = 0, max = 1, value = 0.3,step=0.05),
          numericInput(inputId = 'numberSample',label=NULL,min = 0, max = 0, step = 50,value=0),
          textInput("ROInameSample",label="Name of the ROI",placeholder="type new ROI name here"),
          actionButton("SampleROI","Create ROI")
        )

        

      ),




      fluidRow(

        box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_modifyRois_width",title="Filter for width"),

          selectInput("selectROItoFilterWIDTH",label="Select ROI to filter:",NULL),
          sliderInput('quantileThreshFilterWIDTH',label="Quantile intervals:",min = 0, max = 1, value = c(0,0.9),step=0.002),
          fluidRow(
            column(width=4,HTML("<b>MIN width</b>")),
            column(width=4),
            column(width=4,HTML("<b>MAX width</b>"))
          ),
          fluidRow(
            column(width=4,
              numericInput(inputId = 'absoluteFilter1WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0)
            ),
            column(width=4,
              HTML("<h4>      &lt; --- &gt; </h4>")
            ),
            column(width=4,
              numericInput(inputId = 'absoluteFilter2WIDTH',label=NULL,min = 0, max = 0, step = 50,value=0) 
            )
          ),

          textInput("ROInameFilterWIDTH",label="Name of the ROI",placeholder="type new ROI name here"),
          actionButton("FilterROIWIDTH","Create ROI")

        ),


        box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_modifyRois_enrichment",title="Filter for enrichment"),

          selectInput("selectROItoFilter",label="Select ROI to filter:",NULL),
          selectInput("selectBAMtoFilter",label="Select enrichment to use for filtering:",NULL),
          sliderInput('quantileThreshFilter',label="Quantile intervals:",min = 0, max = 1, value = c(0,0.9),step=0.002),
          fluidRow(
            column(width=4,HTML("<b>MIN enrichment</b>")),
            column(width=4),
            column(width=4,HTML("<b>MAX enrichment</b>"))
          ),
          fluidRow(
            column(width=4,
              numericInput(inputId = 'absoluteFilter1',label=NULL,min = 0, max = 0, step = 50,value=0)
            ),
            column(width=4,
              HTML("<h4>      &lt; --- &gt; </h4>")
            ),
            column(width=4,
              numericInput(inputId = 'absoluteFilter2',label=NULL,min = 0, max = 0, step = 50,value=0) 
            )
          ),

          textInput("ROInameFilter",label="Name of the ROI",placeholder="type new ROI name here"),
          actionButton("FilterROI","Create ROI")

        ),

        box(width=4,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_modifyRois_pattern",title="Extract patterns"),

          radioButtons("choiceWherePattern",label="Choose where to search for the pattern:" ,
                                   choices=c("From a ROI"="fromROI",
                                             "From the entire genome (SLOW!)"="fromGenome"),
                                   selected="fromROI"
          ),
          #selectInput("selectROItoExtractPattern",label="ROI to extract pattern from:",NULL),
          uiOutput("selectWherePattern"),
          #warning in case BSgenome DB not present
          uiOutput("showWarningBSgenome"),
          #motif text input
          textInput("PatternToSearch",label="Select pattern (IUPAC nomenclature)",placeholder="ATCNYGG"),
          #menu or text to choose if both strands or strand in range
          uiOutput("showStrandOptsPattern"),
          #name of the new motif
          textInput("ROInamePattern",label="Name of the ROI",placeholder="type new ROI name here"),
          actionButton("ExtractPatternROI","Create ROI")

        )




      )


    ),






        
    #for each ROI, show features (density plot of the width and /or barplot of ranges number)
    tabPanel("View ROI",
      fluidRow(
        box(width=3,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_viewRois_selectRoi",title="ROI selection"),

          HTML("<b>Select ROI to view:</b>"),
          wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 300px; max-width: 300px; background-color: #ffffff;",
            checkboxGroupInput("confirmviewROI", label=NULL,choices=NULL)
          ),
          
          actionButton("updatechoiceROI", "Update choice")
        ),
        box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_viewRois_visualization",title="Visualization"),

          tabBox(width=12,
            tabPanel("width distribution",
              plotOutput('viewROIwidth'),
              htmlOutput("saveviewpeakswidthROI")
            ),
            tabPanel("intervals number",
              plotOutput("viewpeaksnumberROI"),
              htmlOutput("saveviewpeaksnumberROI")
             
            ) 
          )
        ), 
        
        box(width=3,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_viewRois_information",title="Information"),

            HTML("<br>"),
            sliderInput("quantileROIwidth",label="Select quantile of the width distribution:",min=0,max=1,value=.5,step=0.01),
            htmlOutput("viewROIstat"),
            HTML("<br>"),
            HTML("<b>Where does this ROI come from?</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
              htmlOutput("viewROIsource")
            )
        )
      )
    ),


    #for each ROI, show features (density plot of the width and /or barplot of ranges number)
    tabPanel("Get ROI",
      fluidRow(
        box(width=3,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_getRois_roiSelection",title="ROI selection"),

          HTML("<b>Select ROI:</b>"),
          selectInput("listgetROI",NULL,choices=NULL),
          uiOutput("showROIoptionsToViewRANGE"),
          uiOutput("showROIoptionsToViewMETADATA"),
          uiOutput("showROIoptionsToViewENRICHMENTS"),
          #checkboxInput("putEnrichments", label="Put available enrichments",value = TRUE, width = NULL),
          actionButton("showdataframeROI", "Preview ROI"),

          HTML("<br><br><br><br>"),
          uiOutput("showWindowAnnotation")

        ),

        
        box(width=9,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_getRois_preview",title="Preview"),

          htmlOutput("previewROItodownload"),
          htmlOutput("previewROItodownloadbutton")
          
        )
      )
    ),




    #BAM association tab
    tabPanel("Associate enrichments",

      fluidRow(


        box(width=7,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_associateEnrichments_associateRemove",title="Associate/remove enrichments"),

          fluidRow(
            column(width=6,
              HTML("<b>Choose ROI(s):</b><br>"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 500px; max-width: 300px; background-color: #ffffff;",
                checkboxGroupInput("selectROItoBAMassociate",label=NULL,NULL)
              )
            
            ),
            column(width=6,
              HTML("<h3>Associate enrichment(s) to ROI(s)</h3><br>"),
              HTML("<b>Select enrichment(s) to associate to selected ROI(s):</b><br>"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 180px; max-width: 300px; background-color: #ffffff;",
                checkboxGroupInput("selectBAMtoassociate", NULL,choices=NULL)
              ),
              uiOutput("radioForNorm"),
              uiOutput("menuForNorm"),
              HTML("<b>Number of cores:</b>"),
              numericInput(inputId = 'coresCoverage',label=NULL,min = 1, max = nc, step = 1,value=1),
              HTML("<i>WARNING</i>: time consuming (some minutes)<br>"),
              actionButton("confirmBAMassociate", "Associate!"),
              HTML("<br><br>"),
              HTML('<hr size=3>'),
              HTML("<br>"),
              HTML("<h3>Enrichments associated to ROI(s)</h3><br>"),
              HTML("<b>Select enrichment(s) to remove from selected ROI(s):</b><br>"),
              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 180px; max-width: 300px; background-color: #ffffff;",
                checkboxGroupInput("selectBAMtoDeassociate", NULL,choices=NULL)
              ),
              actionButton("confirmBAMDeassociate", "Remove")            
            )
          )
        ),

        box(width=5,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_renameEnrichments_renameOrderEnrichments",title="Rename/reorder enrichments"),
          selectInput("selectROIforBAMrename",label="Select ROI for rename/reorder enrichments:",NULL),
          HTML("<br>"),
          HTML("<h3>Rename enrichments</h3><br>"),
          #HTML("<br>"),
          selectInput("selectedBAMtoRename",label="Select enrichment to rename:",NULL),
          #HTML("<br><br>"),
          textInput("newBAMname","Select new name:",placeholder="type new enrichment name here",value=""),
          actionButton("renameBAM", "Rename"),
          HTML("<br><br>"),
          HTML('<hr size=3>'),
          HTML("<br>"),
          HTML("<h3>Reorder enrichments</h3><br>"),
          wellPanel(id = "logPanelBAM",style = "overflow-y:scroll; max-height: 400px",
            fluidPage(
              uiOutput("dinamicBAM")
            )
          ),
          #button to confirm to reorder
          actionButton("reorderBAM","Reorder!")
        )
      )

    ),






    tabPanel("GO analyses",
      fluidRow (
        
        column(width=3,

          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_goAnalysis_parameters",title="Parameters"),

            tabBox(width=12,
              tabPanel("Variables",
                HTML("<br><b>Select the source</b></br>"),
                radioButtons("chooseSourceGO",label=NULL,choices=c(
                              "From ROI"="fromROI",
                              "From gene list"="fromGeneList"
                            )),
                #here, the UI (checkboxGroup if from ROI, textInput if from custom list)
                uiOutput("viewSelectGenesGO"),
                #here, put choice of kind of ID of the genes (symbols, etrez, ensembl): only symbols if database (promoters) is not present
                uiOutput("additionalparametersGO"),
                uiOutput("chooseWindowROIGO"),
                #here, we serve all possible genesets using GenesetsGMT global variable
                HTML("<b>Select signature(s):</b>"),
                wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
                  checkboxGroupInput(inputId="selectedGenesetsGO",label=NULL,choices=names(GenesetsGMT))
                ),

                #radiobutton to choose if ranking or clustering the results
                uiOutput("chooseOrderingGO_widget"),
                
                # radioButtons("chooseOrderingGO","How to order results",choices=c(
                #                             "Ranking best padj"="ranking",
                #                             "Custering"="clustering"
                #                       ),selected="ranking") ,
                #if "clsutering" is selected, show radiobutton of the type of clustering ("kmean" or "hierarchical")
                uiOutput("clustertypeGO_widget"),
                uiOutput("clusternumbershowGO_widget"),
                uiOutput("clusterHDistMethodGO_widget"),
                uiOutput("clusterHClustMethodGO_widget"),
                uiOutput("clusterKstartsGO_widget"),
                uiOutput("clusterKiterationsGO_widget"),


                HTML("<b>Min signature size:</b>"),
                numericInput(inputId = 'minSizeGO',label=NULL,min = 1, step = 5,value=15),
                HTML("<b>Max signature size:</b>"),
                numericInput(inputId = 'maxSizeGO',label=NULL,min = 1, step = 5,value=500),              

                actionButton("doTheGO","GO!")
              ),

              tabPanel("Filtering",
                #here, sliderInputs of various thresholds for the analyses
                sliderInput('scaleQuantileGO',label="Quantile threshold for padj colorscale",min = 0.1, max = 1, value = 0.9,step=0.002),
                #geneRatio:
                sliderInput('quantileGeneRatioGO',label="Gene ratio threshold:",min = 0, max = 1, value = 0,step=0.05),
                #-log10Padj:
                sliderInput('log10padjGO',label="-log10 padj threshold:",min = 1, max = 50, value = 2,step=1),
                #top n statistically significant:
                sliderInput('topNGO',label="Top significant hits:",min = 1, max = 100, value = 10,step=1),
                #color scale
                selectInput("colorScaleGO",label="Choose a color scale:",c("white/red"="white_red4",
                                                              "white/blue"="white_blue",
                                                              "white/green"= "white_green4"))
              )

              
            )

          )



        ),
        


        column(width=9,
          fluidRow(
            #put plot (barplot/heatmap)
            box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
              title=boxHelp(ID="msg_goAnalysis_goPlot",title="GO plot"),

              fluidRow(
                column(width=9,#style='padding:0px;',
                  fluidRow(
                    column(width=3,
                      plotOutput("plotMaterialLeft")
                    ),
                    column(width=9,
                      plotOutput("plotOntology",click="GO_click",brush=brushOpts(id="GO_brush",delayType="debounce",delay=300,resetOnNew=TRUE)),#,height=750,width=600),
                      plotOutput("textNameGO"),
                      htmlOutput("saveheatmapGO")
                    #width=600),
                    )
                  )
                ),
                column(width=3,#style='padding:0px;',
                  plotOutput("colorScaleGO",height=100),
                  uiOutput("showTermClicked"),
                  uiOutput("showGenesClicked")
                ) 
              )

            )
            
          ),

          fluidRow(
            #put table to download
            box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
              title=boxHelp(ID="msg_goAnalysis_goTable",title="GO table"),

              wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; background-color: #ffffff;",
                dataTableOutput("tableOntology")
              ),
              uiOutput("tableGOdownloadButton")
            )
          
          )
        )


      )
    ),






    tabPanel("Prepare ROI for heatmap",
      fluidRow(
        column(width=6,
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                  title=boxHelp(ID="msg_PredefPipeline_parameters",title="ROI preparation for heatmap"),

            #ROI menu
            selectInput("selectROIpredefPipeline",label="Select ROI for preparation:",NULL),
            HTML("<br>"),
            HTML('<hr size=3>'),
            HTML("<br>"),
            #subsample %
            sliderInput('quantileThreshPredefPipeline',label="% genomic ranges to keep:",min = 0, max = 1, value = 0.3,step=0.05),
            #this below allows absolute number input
            #numericInput(inputId = 'numberSamplePredefPipeline',label=NULL,min = 0, max = 0, step = 50,value=0),              

            #output UI with absolute number calculated (give warning if > 50k)
            uiOutput("textNumRangesPredefPipeline"),
            HTML("<br>"),
            HTML('<hr size=3>'),
            HTML("<br>"),
            #decide if center on summit (Yes/No); if Yes, menu with enrichment associated
            radioButtons("choiceSummitPredefPipeline",label="Do you want to center on summit?" ,
                                      choices=c("Yes"="Yes",
                                                "No"="No"),
                                      selected="No"
            ),
            uiOutput("menuSummitPredefPipeline"),       
            #selectInput("selectBAMtoCenterPredefPipeline",label="Enrichment to use for summit:",NULL),
            HTML("<br>"),
            HTML('<hr size=3>'),
            HTML("<br>"),
            #set intervals from midpoint/summit
            fluidRow(
              column(width=4,HTML("<b>Upstream</b>")),
              column(width=4),
              column(width=4,HTML("<b>Downstream</b>"))
            ),
            fluidRow(
              column(width=4,
                numericInput(inputId = 'sliderUpstreamPredefPipeline',label=NULL,min = 0, max = 20000, step = 50,value=2000) 
              ),
              column(width=4,
                HTML("<h4>Center</h4>")
              ),
              column(width=4,
                numericInput(inputId = 'sliderDownstreamPredefPipeline',label=NULL,min = 0, max = 20000, step = 50,value=2000) 
              )
            ),        

            HTML("<br>"),
            HTML('<hr size=3>'),
            HTML("<br>"),
            #select enrichment files for final ROI (from enrichmnt opened, "re-do" association)
            uiOutput("menuEnrichPredefPipeline"),
            # wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 180px; max-width: 300px; background-color: #ffffff;",
            #         checkboxGroupInput("selectBAMassociatePredefPipeline", NULL,choices=NULL)
            # ),
            #number of cores for this association:
            HTML("<b>Number of cores:</b>"),
            numericInput(inputId = 'coresPredefPipeline',label=NULL,min = 1, max = nc, step = 1,value=nc),

            HTML("<br>"),
            HTML('<hr size=3>'),
            HTML("<br>"), 


            HTML("<b>New ROI name:</b>"),
            textInput("ROInamePredefPipeline",label=NULL,placeholder="type new ROI name here"),
            actionButton("PrepareROIpredefPipeline","Create ROI")
          )      
        ),
        column(width=6)
        
      )
      
     

    )    



  )         
)




