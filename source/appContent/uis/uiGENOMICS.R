
tabGENOMICS<-tabItem (tabName = "GENOMICSblock",
  tabsetPanel(type="pills",


  	#std plots...
    tabPanel("Single evaluation",
      #double radio button : only 2 ranges, 
      #Venn diagrams, 
 
      fluidRow(

        column(width=3,
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_singleEvaluation_parameters",title="Parameters"),

            selectInput("ROIchooseSingleEval", "Select ROI:",NULL),
            selectInput("BAMchooseSingleEval", "Select enrichment:",NULL),
            HTML("<br><br>"),
            radioButtons("chooseNormalizationSingleEval","Choose normalization:",choices=c(
                                                "Total reads (rpm)"="totread",
                                                "Read density (rpm/bp)"="readdensity"
                                                      )),
            selectInput("chooseColorPaletteSingleEval","Choose color palette:",choices=c(
                                                  "black/red/blue/grey"="black_red_blue_gray",
                                                  "black/red/orange/green"="black_red_orange_green"

                                                )),
            HTML("<br>"),
            actionButton("plotSingleEval","Update plot")
          )

        ),

        column(width=9,
          fluidRow(
            fluidRow(

              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_distribution",title="Distribution"),

                tabBox(width=12,id="Peaks location",
                  tabPanel("Piechart",id="piechartSingleEval",
                    plotOutput('viewDistributionPieSingleEval'),
                    htmlOutput("saveviewDistributionPieSingleEval")
                  ),
                  tabPanel("Barplot",id="barplotSingleEval",
                    plotOutput("viewDistributionBarSingleEval"),
                    htmlOutput("saveviewDistributionBarSingleEval")
                  ) 
                )
              ),




              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_widthDistribution",title="Width distribution"),

                plotOutput("widthDistributionSingleEval"),
                htmlOutput("savewidthDistributionSingleEval")
              )
            ),

            fluidRow(
              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_enrichmentBoxplot",title="Enrichment boxplot"),

                plotOutput("enrichmentBoxSingleEval"),
                htmlOutput("saveenrichmentBoxSingleEval"),
                htmlOutput("saveboxdataSingleEval")
                #downloadButton('saveenrichmentBoxSingleEvaldata', 'Save data')

              ),
              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_peakProfile",title="Peak average profile"),

                plotOutput("TSSprofileSingleEval"),
                htmlOutput("saveenrichmentProfileSingleEval")
              )
            )

          )
        )

      )
      
    ),






    #std cmp...
    tabPanel("Pairwise overlaps",

      fluidRow(

        column(width=3,
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_pairwiseOverlaps_parameters",title="Parameters"),

            selectInput("ROI1chooseCmp", "Choose ROI-1:",NULL),
            selectInput("ROI2chooseCmp", "Choose ROI-2:",NULL),
            HTML("<br>"),
            HTML("<b>Minimum number of bp for overlap:</b>"),
            numericInput(inputId = 'minOverlapCmp',label=NULL,min = 1, step = 5,value=1),

            uiOutput("BAMmenuchooseCmp1"),
            uiOutput("BAMmenuchooseCmp2"),

            # HTML("<b>Choose enrichment1:</b>"),
            # selectInput(inputId="BAM1chooseCmp", label=NULL,choices=character(0)),
            # HTML("<b>Choose enrichment2:</b>"),
            # selectInput(inputId="BAM2chooseCmp", label=NULL,choices=character(0)),

            HTML("<br>"),
            checkboxInput("islogCmp", label="log2",value = FALSE, width = NULL),
            radioButtons("chooseNormalizationCmp","Choose normalization:",choices=c(
                                                "Total reads (rpm)"="totread",
                                                "Read density (rpm/bp)"="readdensity"
                                                      )),
            selectInput("chooseColorPaletteCmp","Choose color palette:",choices=c(
                                                  "red/grey"="red_gray_red4_grey20",
                                                  "blue/green"="blue_green_blue4_green4"
                                                )),
            uiOutput("showScatterChoice"),
            HTML("<br>"),
            actionButton("plotCmp","Update plot")
          )

        ),

        column(width=9,

          box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_pairwiseOverlaps_overlap",title="Overlap"),

            tabBox(width=12,id="Overlap",
              tabPanel("Barplot",id="barplotCmp",
                plotOutput("viewBarplotCmp"),
                htmlOutput("saveviewBarplotCmp")
              ),
              tabPanel("Venn",id="vennCmp",
                plotOutput('viewVennCmp'),
                htmlOutput("saveviewVennCmp")
              )

            )
          ),


          box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_pairwiseOverlaps_overlapAndEnrichment",title="Overlap and enrichment"),

            tabBox(width=12,id="Enrichments",
              tabPanel("Boxplot",id="boxplotCmp",
                plotOutput('viewBoxplotCmp'),
                htmlOutput("saveviewBoxplotCmp"),
                htmlOutput("saveboxdataCmp")
              ),
              tabPanel("Scatterplot",id="scatterplotCmp",
                plotOutput("viewScatterplotCmp"),
                htmlOutput("saveviewScatterplotCmp"),
                htmlOutput("saveScatterdataCmp")
              ),
              tabPanel("Calibration",id="CalibrationCmp",
                plotOutput("viewCalibrationCmp"),
                htmlOutput("saveviewCalibrationCmp")
              ) 
            )
          )
        )
      )
    ),


    #Digital heatmap
    tabPanel("Position-based Heatmap",
      fluidRow (

        column(width=3,


          box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
            title=boxHelp(ID="msg_digitalHeatmap_parameters",title="Parameters"),

            tabBox(width=12,
              tabPanel("Variables",
                HTML("<br>"),
                actionButton("confirmUpdateDigitalHeat1", "Update plot"),
                HTML("<br><br>"),
                #ROI to select that guide the heatmap (master ROI)
                HTML("<b>Master ROI(s):</b>"),
                wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
                  checkboxGroupInput("ROImaster",NULL,NULL)
                ),
                #ROIs available to be viewed
                HTML("<b>ROIs to view:</b>"),
                wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
                  checkboxGroupInput("ROIsForDigitalHeat",NULL,NULL)
                ),
                
                #reorder ROI in digital heatmap
                uiOutput("reorderROImenuDigitalHeat"),

                HTML("<b>ROIs for cluster:</b>"),
                wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
                  checkboxGroupInput("ROIforClusteringDigitalHeat",NULL,choices=NULL)
                ),

                uiOutput("clustertypeDigitalHeat"),
                uiOutput("clusternumbershowDigitalHeat"),
                uiOutput("clusterHDistMethodDigitalHeat"),
                uiOutput("clusterHClustMethodDigitalHeat"),
                uiOutput("clusterKstartsDigitalHeat"),
                uiOutput("clusterKiterationsDigitalHeat"),

                HTML("<b>Number of bins:</b>"),
                numericInput(inputId = 'binsDigitalHeat',label=NULL,min = 1, max = 50, step = 1,value=10),
                #check box for strand-specific overlap:
                checkboxInput("StrandSpecOverlap", label="Strand-specific overlaps",value = FALSE, width = NULL)
                

                
              ),
              tabPanel("Advanced",
                HTML("<br>"),
                actionButton("confirmUpdateDigitalHeat2", "Update plot"),
                HTML("<br><br>"),
                HTML("<b>Random sample of genomic ranges to show:</b>"),
                numericInput(inputId = 'sampleRandomDigitalHeat',label=NULL,min = 0, max = 0, step = 1000,value=0),

                radioButtons("optioncolorsforDigitalHeat",label="Select colors:",choiceNames=c("global color","custom colors"),choiceValues=c("global","custom"),selected="global"),
                uiOutput("showcolorsDigitalheat"),

                checkboxInput("FracToPercDigitalHeat", label="positional overlap %",value = TRUE, width = NULL)
                #HTML("<b>Universe for overlaps:</b>"),
                #radioButtons("chooseOrderingDigitalHeat",label=NULL,choices=c("All acessible sites"="global","Sites of this ROI"="local")),
                # selectInput("distmethodDigitalHeat",label="Distance method:",c("Euclidean"="euclidean",
                #                                                                 "Manhattan"="manhattan",
                #                                                                 #"Canberra"="canberra",
                #                                                                 "Minkowski"="minkowski")),
                # selectInput("clustmethodDigitalHeat",label="Clustering method:",c("Average"="average",
                #                                                                  "Complete"="complete",
                #                                                                  "Median"="median",
                #                                                                  "Centroid"="centroid")),
                
              )        
            )

          )


        ),
        


        column(width=9,
          fluidRow(
            # box(width=1,height=600,title=" ",
            #   plotOutput("clustersImageLeft",height=350,click="rowdendrogram_click_....."),
            #   plotOutput("emptyplot",height=250)
            # ),
            box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_digitalHeatmap_heatmap",title="Heatmap"),
              #mini-fluid page which puts together heatmap, color scale....
              fluidRow(
                column(width=9,#style='padding:0px;',
                  fluidRow(
                    column(width=1,
                      plotOutput("clustersImageLeftDigital",click="rowdendrogram_click_Digital")
                    ),
                    column(width=11,
                      plotOutput("heatmapDigital",click="heatmapDigital_click",brush=brushOpts(id="heatmapDigital_brush",delayType="debounce",delay=300,resetOnNew=TRUE)),
                      plotOutput("textNameDigitalHeat",height=200),

                      fluidRow(
                        column(width=6,
                          htmlOutput("saveheatmapDigital")
                        ),
                        column(width=6#,
                          #htmlOutput("showsaveDigitalHeatdata")
                        )
                      )
                    )
                  )
                ),
                column(width=3,#style='padding:0px;',
                  htmlOutput("textfractionelementsDigitalHeat"),
                  htmlOutput("textselectedelementsDigitalHeat"),
                  uiOutput("newROIfromDigitalHeat_out")
                  # htmlOutput("textfractionelementsDigitalHeat"),
                  # htmlOutput("textselectedelementsDigitalHeat"),
                  # plotOutput("colorScaleDigitalHeat",height=90),
                  # uiOutput("newROIfromDigitalHeat_out")
                )  
              )
            )
          ),


          fluidRow(
            box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_digitalHeatmap_jaccardIdx",title="Jaccard idx"),
              plotOutput("JaccardDigital"),
              htmlOutput("saveJaccardDigital")
            ),
            box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_digitalHeatmap_overlapBias",title="Positional distribution of overlaps"),
              plotOutput("frequencyOvDigitalHeat"),
              htmlOutput("savefrequencyOvDigitalHeat")
            )  
          )
        )


      )

    ),






    #analogic heatmaps...
    tabPanel (  "Enrichment-based Heatmap",


      fluidRow (

        column(width=3,



          box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
            title=boxHelp(ID="msg_analogicHeatmap_parameters",title="Parameters"),

            tabBox(width=12,
              tabPanel("Variables",
                HTML("<br>"),
                actionButton("confirmUpdateAnalogHeat", "Update plot"),
                HTML("<br><br>"),
                HTML("<b>Select ROI(s):</b>"),
                wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
                  checkboxGroupInput("ROIsForAnalogHeat",NULL,NULL)
                ),

                HTML("<b>Select enrichments to show:</b>"),
                wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
                  checkboxGroupInput("BAMsForAnalogHeat",NULL,NULL)
                ),
                #menu for order the selected enrichments
                uiOutput("reorderBAMmenuAnalogHeat"),
                HTML("<b>Clustering/ranking</b>"),
                radioButtons("chooseOrderingAnalogHeat",label=NULL,choices=c("Ranking"="ranking","Clustering"="clustering")),
                uiOutput("orderingAnalogHeat"),
                uiOutput("clustertypeAnalogHeat"),
                uiOutput("clusternumbershowAnalogHeat"),
                uiOutput("clusterHDistMethodAnalogHeat"),
                uiOutput("clusterHClustMethodAnalogHeat"),
                uiOutput("clusterKstartsAnalogHeat"),
                uiOutput("clusterKiterationsAnalogHeat"),

                numericInput(inputId = 'binsAnalogHeat',label="Number of bins:",min = 1, max = 200, step = 1,value=50)

            
                # HTML("<b>New ROI from heatmap:</b>"),

                #,
                #uiOutput("confirmImportROIfromAnalogHeat_out")
                # textInput("newROIfromAnalogHeat",NULL,value="",placeholder = "my_ROI"),
                # actionButton("confirmImportROIfromAnalogHeat", "Import ROI")

              ),

              tabPanel("Advanced",
                HTML("<br>"),
                actionButton("confirmUpdateAnalogHeat2", "Update plot"),
                HTML("<br><br>"),    
                numericInput(inputId = 'sampleRandomAnalogHeat',label="Random sample of genomic ranges to show:",min = 0, max = 0, step = 1000,value=0),


                checkboxInput("Log2BoxAnalogHeat", label="log2",value = FALSE, width = NULL),


                HTML("<b>Quantile threshold</b>"),
                radioButtons("chooseQuantileMethodAnalogHeat",label=NULL,choices=c("Uniform"="allBAM","Individual"="eachBAM")),
                sliderInput('quantileThreshAnalogHeat',label=NULL,min = 0.1, max = 1, value = 0.9,step=0.002),
          
                radioButtons("optioncolorsforAnalogHeat",label="Select colors:",choiceNames=c("default color","custom colors"),choiceValues=c("global","custom"),selected="global"),
                
                uiOutput("showcolorsheat"),

                HTML("<b></b>"),
                checkboxInput("GroupColorsAnalogHeat", label="Group colors (boxes)",value = FALSE, width = NULL)

              )
              
            )
          )
        ),
        
        column(width=9,
          fluidRow(
            # box(width=1,height=600,title=" ",
            #   plotOutput("clustersImageLeft",height=350,click="rowdendrogram_click_Analog"),
            #   plotOutput("emptyplot",height=250)
            # ),
            box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_analogicHeatmap_heatmap",title="Heatmap"),

              # div(
              #   dropdownButton(        
              #       .......      
              #       circle = TRUE, status = "danger", icon = icon("gear"),right=TRUE
              #     ),class = "pull-right",id = "moveme"
              # ),
              #mini-fluid page which puts together heatmap, color scale....
              fluidRow(
                column(width=9,#style='padding:0px;',
                  fluidRow(
                    column(width=1,
                      plotOutput("clustersImageLeft",click="rowdendrogram_click_Analog")
                    ),
                    column(width=11,
                      plotOutput("heatmapAnalog",click="heatmap_click",brush=brushOpts(id="heatmap_brush",delayType="debounce",delay=300,resetOnNew=TRUE)),#,height=750,width=600),
                      plotOutput("textNameAnalogHeat",height=200),
                      fluidRow(
                        column(width=6,
                          htmlOutput("saveheatmapAnalog")
                        ),
                        column(width=6,
                          htmlOutput("showsaveAnalogHeatdata")
                        )
                      )
                    )
                    #width=600),
                  )
                ),
                column(width=3,#style='padding:0px;',
                  htmlOutput("textfractionelementsAnalogHeat"),
                  htmlOutput("textselectedelementsAnalogHeat"),
                  plotOutput("colorScaleAnalogHeat",height=90),
                  uiOutput("newROIfromAnalogHeat_out")
                )  
              )
            )
          ),


          fluidRow(
            box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_analogicHeatmap_profiles",title="Profiles"),
              plotOutput("profileAnalogHeat"),
              htmlOutput("saveprofileAnalogHeat")
            ),

            box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_analogicHeatmap_enrichments",title="Enrichments"),

              tabBox(width=12,

                tabPanel("Boxplot by ROI/cluster",
                  plotOutput("boxplotByROIAnalogHeat"),
                  htmlOutput("saveboxplotByROIAnalogHeat")
                ),
                tabPanel("Boxplot by enrichment",
                  plotOutput("boxplotByBAMAnalogHeat"),
                  htmlOutput("saveboxplotByBAMAnalogHeat")
                ),
                tabPanel("cor",
                  plotOutput("corAnalogHeat"),
                  htmlOutput("savecorAnalogHeat")
                ),
                tabPanel("pcor",
                  plotOutput("pcorAnalogHeat"),
                  htmlOutput("savepcorAnalogHeat")
                )          
              ) 
            )
          )
        )
      )

    ),



    #profiles and boxplots
    tabPanel("Enrichments in ROIs",
      fluidRow (

        column(width=3,
          box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
            title=boxHelp(ID="msg_enrichmentInRois_parameters",title="Parameters"),

            HTML("<b>Select ROI(s):</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("ROIsForProfilesAndBox",NULL,NULL)
            ),

            HTML("<b>Select enrichments to show:</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("BAMsForProfilesAndBox",NULL,NULL)
            ),

            HTML("<b>Number of bins:</b>"),
            numericInput(inputId = 'binsProfilesAndBox',label=NULL,min = 1, max = 500, step = 1,value=50),
            HTML("<br>"),
            radioButtons("chooseNormalizationProfilesAndBox",label="Normalization method for enrichments:",choices=c(
                                                  "Total reads (rpm)"="totread",
                                                  "Read density (rpm/bp)"="readdensity"
                                                        ),selected="readdensity"),
            HTML("<b></b>"),
            checkboxInput("Log2BoxProfilesAndBox", label="log2",value = FALSE, width = NULL),
            HTML("<br>"),
            
            radioButtons("choosecorMethodProfilesAndBox",label="Type of correlation",choices=c(
                                                  "Pearson"="pearson",
                                                  "Spearman"="spearman"
                                                        ),selected="pearson"),
            HTML("<b></b>"),
            checkboxInput("GroupColorsProfilesAndBox", label="Group colors (boxes)",value = FALSE, width = NULL),
            HTML("<b></b>"),
            actionButton("confirmUpdateProfilesAndBox", "Update plot")
          )

        ),

        column(width=9,
          fluidRow(
            fluidRow(
              box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_enrichmentInRois_profiles",title="Profiles"),

                plotOutput("profileProfilesAndBox"),
                htmlOutput("saveprofileProfilesAndBox")
              ),

              box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_enrichmentInRois_boxplots",title="Boxplots"),

                tabBox(width=12,height=500,
                  tabPanel("Boxplot by ROI",
                    plotOutput("boxByROIProfilesAndBox"),
                    htmlOutput("saveboxByROIProfilesAndBox"),
                    htmlOutput("saveboxdataProfANDbox")
                  ),

                  tabPanel("Boxplot by enrichment",
                    plotOutput("boxByBAMProfilesAndBox"),
                    htmlOutput("saveboxByBAMProfilesAndBox")
                  )            
                )
              )
            ),

            fluidRow(


              box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_enrichmentInRois_correlations",title="Correlations"),

                tabBox(width=12,height=500,
                  tabPanel("Cor-Heatmap",
                    plotOutput("corProfilesAndBox",click="cor_click"),
                    htmlOutput("savecorProfilesAndBox")
                  ),

                  tabPanel("Pcor-Heatmap",
                    plotOutput("pcorProfilesAndBox",click="pcor_click"),
                    htmlOutput("savepcorProfilesAndBox")
                  )            
                )
              ),



              box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_enrichmentInRois_scatterplot",title="Scatterplot"),

                plotOutput("scatterProfilesAndBox"),
                htmlOutput("savescatterProfilesAndBox")
              )
            )

          )
        )
          

      )
    ),



    #dynamics on genes
    tabPanel("Metagene profiles",
      #The USER has to be careful to use down/upstream distances from the TSS while retrieving the
      #TSS/TES. If not satisfied with this values, the user shoult re-construct original promoters
      #and transcripts and TES from the database
      fluidRow(

        column(width=3,
          box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
            title=boxHelp(ID="msg_dynamicsOnGenes_parameters",title="Parameters"),

            HTML("<b>Select the gene list:</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("genelistsforDynamics",NULL,NULL)
            ),
            
            uiOutput("showROItriadGeneList"),
            HTML("<br>"),
            HTML("<b>Select enrichments to show:</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("BAMforDynamics",NULL,NULL)
            ),
            
            HTML("<br>"),
            HTML("<b>Number of bins:</b>"),
            numericInput(inputId = 'binsforDynamics',label=NULL,min = 1, max = 300, step = 1,value=100),
            HTML("<br>"),
            radioButtons("chooseMetricforDynamics",label="Mean or median?",choices=c(
                                                "median"="median",
                                                "mean"="mean"
                                                      ),selected="mean"),
            HTML("<br>"),
            checkboxInput("islogforDynamics", label="log2",value = FALSE, width = NULL),
            radioButtons("chooseNormalizationforDynamics","Choose normalization:",choices=c(
                                                "Total reads (rpm)"="totread",
                                                "Read density (rpm/bp)"="readdensity"
                                                      ),selected="readdensity"),
            HTML("<br>"),
            HTML("<b>Fraction of outliers to exclude in cumulative plots:</b>"),
            sliderInput('percentageOutlayerCumulPlots',label=NULL,min = 0, max = 0.3, value = 0.05,step=0.01),
            HTML("<br>"),
            actionButton("plotDynamics","Update plot")
          )

        ),

        column(width=9,

          #profile
          fluidRow(
            column(width=12,
              box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_profiles",title="Metagene profile"),

                plotOutput("plotProfileDynamics"),
                htmlOutput("saveprofileDynamics")
              )
                           
            )
          ),
          #boxplot
          fluidRow(
            column(width=4,
              box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_TSSenrichment",title="TSS enrichment (boxplot)"),

                plotOutput("plotboxTSSDynamics"),
                htmlOutput("saveboxTSSDynamics"),
                htmlOutput("saveboxdatadynamicsTSS")
              )
            ),
            column(width=4,
              box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_genebodiesEnrichment",title="Genebodies enrichment (boxplot)"),

                plotOutput("plotboxGBDynamics"),
                htmlOutput("saveboxGBDynamics"),
                htmlOutput("saveboxdatadynamicsGB")
              )
            ),
            column(width=4,
              box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_TESenrichment",title="TES enrichment (boxplot)"),

                plotOutput("plotboxTESDynamics"),
                htmlOutput("saveboxTESDynamics"),
                htmlOutput("saveboxdatadynamicsTES")
              )
            )
          ),
          #stallingindex
          fluidRow(
            column(width=4,
              box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_TSSranked",title="TSS ranked enrichments"),

                plotOutput("plotSITSSDynamics"),
                htmlOutput("saveSITSSDynamics")
              )
            ),
            column(width=4,
              box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_genebodiesRanked",title="Genebodies ranked enrichments"),

                plotOutput("plotSIGBDynamics"),
                htmlOutput("saveSIGBDynamics")
              )
            ),
            column(width=4,
              box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_stallingIndexRanked",title="Stalling Index"),

                plotOutput("plotSISIDynamics"),
                htmlOutput("saveSISIDynamics")
              )
            )
          )

        )




      )
    )


  )
)
