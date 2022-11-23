
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
            selectInput("ROIchooseSingleEval", "1) Select ROI:",NULL),
            #selectInput("BAMchooseSingleEval", "2) Select enrichment (optional):",NULL),
            uiOutput("BAMmenuchoose_singleeval"),
            HTML("<br>"),
            actionButton("plotSingleEval","Update plot")
          )

        ),

        column(width=9,
          fluidRow(
            fluidRow(

              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_distribution",title="Distribution"),
                column(width=8,
                  tabBox(width=12,id="Peaks location",
                    tabPanel("Piechart",id="piechartSingleEval", 
                      plotOutput('viewDistributionPieSingleEval'),
                      htmlOutput("saveviewDistributionPieSingleEval"),
                    ),
                    tabPanel("Barplot",id="barplotSingleEval",
                      plotOutput("viewDistributionBarSingleEval"),
                      htmlOutput("saveviewDistributionBarSingleEval")
                    )
                  )

                ),column(width=4,
                  uiOutput("piechartSingleEval_options")
                ) 

              ),




              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_widthDistribution",title="Width distribution"),

                plotOutput("widthDistributionSingleEval"),
                uiOutput("densitySingleEval_options"),
                htmlOutput("savewidthDistributionSingleEval")
              )
            ),

            fluidRow(
              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_enrichmentBoxplot",title="Enrichment boxplot"),
                column(width=8,
                  plotOutput("enrichmentBoxSingleEval"),
                  htmlOutput("saveenrichmentBoxSingleEval"),
                  htmlOutput("saveboxdataSingleEval")
                ),column(width=4,
                  uiOutput("boxSingleEval_options")
                )
              ),
              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_peakProfile",title="Peak average profile"),
                column(width=8,
                  plotOutput("TSSprofileSingleEval"),
                  htmlOutput("saveenrichmentProfileSingleEval")
                ),column(width=4,
                  uiOutput("profileSingleEval_options")
                )

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
            uiOutput("show_minoverlapcmp"),

            uiOutput("BAMmenuchooseCmp1"),
            uiOutput("BAMmenuchooseCmp2"),

            HTML("<br>"),
            actionButton("plotCmp","Update plot")
          )

        ),

        column(width=9,

          fluidRow(

            fluidRow(
              box(width=7,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_pairwiseOverlaps_overlap",title="Overlap"),
                column(width=8,
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
                ),column(width=4,
                  uiOutput("pairwiseoverlaps_overlap_options")
                )

              ),


              box(width=5,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                title=boxHelp(ID="msg_singleEvaluation_enrichmentBoxplot",title="Enrichment boxplot"),
                column(width=8,
                  title=boxHelp(ID="msg_pairwiseOverlaps_box",title="Overlap and enrichment"),
                  plotOutput('viewBoxplotCmp'),
                  htmlOutput("saveviewBoxplotCmp"),
                  htmlOutput("saveboxdataCmp")
                ),column(width=4,
                  uiOutput("pairwiseoverlaps_box_options")
                )
              )
            ),

            fluidRow(
              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title=boxHelp(ID="msg_pairwiseOverlaps_scatter",title="Enrichment scatterplot"),
                column(width=8,
                  plotOutput("viewScatterplotCmp"),
                  htmlOutput("saveviewScatterplotCmp"),
                  htmlOutput("saveScatterdataCmp")
                ),column(width=4,
                  uiOutput("pairwiseoverlaps_scatter_options")
                )
              ),

              box(width=6,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                      title=boxHelp(ID="msg_pairwiseOverlaps_calibration",title="Enrichment calibration"),
                column(width=8,
                  plotOutput("viewCalibrationCmp"),
                  htmlOutput("saveviewCalibrationCmp")
                ),column(width=4,
                  uiOutput("pairwiseoverlaps_calibration_options")
                )
              )
            )

          )

        )
      ),


    ),




    #Digital heatmap
    tabPanel("Position-based Heatmap",
      fluidRow (

        column(width=3,


          box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
            title=boxHelp(ID="msg_digitalHeatmap_parameters",title="Parameters"),

            #ROI to select that guide the heatmap (master ROI)
            list(HTML("<b>Master ROI(s):</b>"),htmlhelp("","help_digitalHeatmap_parameters_masterROI")),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("ROImaster",NULL,NULL)
            ),


            #ROIs available to be viewed. SHould appear after selecting the master ROI
            uiOutput("ROIsForDigitalHeat_menu"),
            #reorder ROI in digital heatmap. 
            uiOutput("reorderROImenuDigitalHeat"),
            #show ROI for clustering. SHould appear after selecting the ROIs available to be viewed
            uiOutput("ROIsForClusterDigital_menu"),
            uiOutput("clustertypeDigitalHeat"),
            uiOutput("clusternumbershowDigitalHeat"),
            uiOutput("clusterHDistMethodDigitalHeat"),
            uiOutput("clusterHClustMethodDigitalHeat"),
            uiOutput("clusterKstartsDigitalHeat"),
            uiOutput("clusterKiterationsDigitalHeat"),
            uiOutput("showbinsDigitalHeat"),
            #check box for strand-specific overlap:
            uiOutput("showStrandSpecOverlap"),
            #checkboxInput("StrandSpecOverlap", label="Strand-specific overlaps",value = FALSE, width = NULL),
            uiOutput("sampleRandomROIDigital"),
            HTML("<br>"),
            uiOutput("show_confirmUpdateDigitalHeat1")   

          )


        ),
        


        column(width=9,
          fluidRow(

            box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_digitalHeatmap_heatmap",title="Heatmap"),
              #mini-fluid page which puts together heatmap, color scale....
              fluidRow(
                column(width=12,
                  uiOutput("show_clickmessage_digital")
                )
              ),
              fluidRow(
                column(width=9,#style='padding:0px;',
                  fluidRow(
                    column(width=1,
                      plotOutput("clustersImageLeftDigital",click="rowdendrogram_click_Digital"),
                      plotOutput("textNameClustDigitalHeat",height=250)
                    ),
                    column(width=11,
                      plotOutput("heatmapDigital",click="heatmapDigital_click",brush=brushOpts(id="heatmapDigital_brush",delayType="debounce",delay=300,resetOnNew=TRUE)),
                      plotOutput("textNameDigitalHeat",height=250),

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
                  #here put specific graphical options for heat (colors(global/custom) and menu)
                  uiOutput("textfractionelementsDigitalHeat"),
                  uiOutput("textselectedelementsDigitalHeat"),
                  uiOutput("newROIfromDigitalHeat_out"),
                  HTML("<br><br>"),
                  uiOutput("showoptioncolorsforDigitalHeat"),
                  uiOutput("showcolorsDigitalheat")
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
              column(width=8,
                plotOutput("frequencyOvDigitalHeat"),
                htmlOutput("savefrequencyOvDigitalHeat")
              ),column(width=4,
                uiOutput("frequencyDigitalHeat_options")
              )
               
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

                list(HTML("<b>Select ROI(s):</b>"),htmlhelp("","help_analogicHeatmap_parameters_ROI")),
                wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 200px; max-width: 300px; background-color: #ffffff;",
                  checkboxGroupInput("ROIsForAnalogHeat",NULL,NULL)
                ),
                #dynamic menu for parameters for analog heatmap
                uiOutput("showBAMsForAnalogHeat"),
                uiOutput("reorderBAMmenuAnalogHeat"),
                uiOutput("showrankingmethod"),
                uiOutput("orderingAnalogHeat"),
                uiOutput("clustertypeAnalogHeat"),
                uiOutput("clusternumbershowAnalogHeat"),
                uiOutput("clusterHDistMethodAnalogHeat"),
                uiOutput("clusterHClustMethodAnalogHeat"),
                uiOutput("clusterKstartsAnalogHeat"),
                uiOutput("clusterKiterationsAnalogHeat"),
                uiOutput("showbinsAnalogHeat"),
                uiOutput("showsampleRandomAnalogHeat"),
                HTML("<br>"),
                uiOutput("show_confirmUpdateAnalogHeat")
          )
        ),
        
        column(width=9,

          fluidRow(
            box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_analogicHeatmap_heatmap",title="Heatmap"),
              fluidRow(
                column(width=10,
                  uiOutput("show_clickmessage_analogic")
                ),
                column(width=2,
                  plotOutput("colorScaleAnalogHeat",height=80)
                )                
              ),
              #mini-fluid page which puts together heatmap, color scale....
              fluidRow(
                column(width=9,#style='padding:0px;',
                  fluidRow(
                    column(width=1,
                      plotOutput("clustersImageLeft",click="rowdendrogram_click_Analog"),
                      plotOutput("textNameClustAnalogHeat",height=250),
                    ),
                    column(width=11,
                      plotOutput("heatmapAnalog",click="heatmap_click",brush=brushOpts(id="heatmap_brush",delayType="debounce",delay=300,resetOnNew=TRUE)),#,height=750,width=600),
                      plotOutput("textNameAnalogHeat",height=250),
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
                  
                  
                  uiOutput("textfractionelementsAnalogHeat"),
                  uiOutput("textselectedelementsAnalogHeat"),
                  uiOutput("newROIfromAnalogHeat_out"),
                  HTML("<br><br>"),
                  uiOutput("showoptioncolorsforAnalogHeat"),                
                  uiOutput("showcolorsheat"),
                  uiOutput("showchooseQuantileMethodAnalogHeat"),
                  uiOutput("showquantileThreshAnalogHeat")


                )  
              )
            )
          ),

          #profile analog heat box
          fluidRow(
            box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_analogicHeatmap_profiles",title="Profiles"),



              
              fluidRow(
                column(width=12,
                  plotOutput("profileAnalogHeat"),
                  htmlOutput("saveprofileAnalogHeat")
                )
              ),
              #now fluidRow with options
              fluidRow(
                column(width=4,
                  uiOutput("showprofileAnalogHeat_logOptions"),
                  uiOutput("showprofileAnalogHeat_colorschemeOptions")
                ),
                column(width=8,
                  uiOutput("showprofileAnalogHeat_colorlistOptions")
                )
              )
                
     
              
            ),

            box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
              title=boxHelp(ID="msg_analogicHeatmap_enrichments",title="Enrichments"),

              fluidRow(
                column(width=12,
                  tabBox(width=12,
                    tabPanel("Boxplot by ROI/cluster",
                      plotOutput("boxplotByROIAnalogHeat"),
                      htmlOutput("saveboxplotByROIAnalogHeat")
                    ),
                    tabPanel("Boxplot by enrichment",
                      plotOutput("boxplotByBAMAnalogHeat"),
                      htmlOutput("saveboxplotByBAMAnalogHeat")
                    )
                    # tabPanel("cor",
                    #   plotOutput("corAnalogHeat"),
                    #   htmlOutput("savecorAnalogHeat")
                    # ),
                    # tabPanel("pcor",
                    #   plotOutput("pcorAnalogHeat"),
                    #   htmlOutput("savepcorAnalogHeat")
                    # )          
                  )
                )
              ),
              fluidRow(
                column(width=4,
                  uiOutput("showboxAnalogHeat_logOptions"),
                  uiOutput("showboxAnalogHeat_colorschemeOptions")
                ),
                column(width=8,
                  uiOutput("showboxAnalogHeat_groupcolOptions"),
                  uiOutput("showboxAnalogHeat_colorlistOptions")
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

            list(HTML("<b>Select ROI(s):</b>"),htmlhelp("","help_enrichmentInRois_parameters_ROIs")),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 150px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("ROIsForProfilesAndBox",NULL,NULL)
            ),
            uiOutput("showBAMsforProfilesAndBox"),
            uiOutput("showbinsforProfilesAndBox"),
            HTML("<br>"),
            uiOutput("show_confirmUpdateProfilesAndBox")
          )

        ),

        column(width=9,
          fluidRow(
            fluidRow(
              box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_enrichmentInRois_profiles",title="Profiles"),

                fluidRow(
                  column(width=12,
                    plotOutput("profileProfilesAndBox"),
                    htmlOutput("saveprofileProfilesAndBox")
                  )
                ),
                #now fluidRow with options
                fluidRow(
                  column(width=4,
                    uiOutput("showprofileProfileAndBox_logOptions"),
                    uiOutput("showprofileProfileAndBox_colorschemeOptions"),
                    uiOutput("showprofileProfileAndBox_isdensityOptions")
                  ),
                  column(width=8,
                    uiOutput("showprofileProfileAndBox_colorlistOptions")
                  )
                )

                
              ),

              box(width=6,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_enrichmentInRois_boxplots",title="Boxplots"),


              fluidRow(
                column(width=12,

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
                column(width=4,
                  uiOutput("showBoxProfileAndBox_logOptions"),
                  uiOutput("showBoxProfileAndBox_colorschemeOptions"),
                  uiOutput("showBoxProfileAndBox_isdensityOptions")
                ),
                column(width=8,
                  uiOutput("showBoxProfileAndBox_groupcolOptions"),
                  uiOutput("showBoxProfileAndBox_colorlistOptions")
                )
              )              



              )
            ),

            fluidRow(


              box(width=7,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_enrichmentInRois_correlations",title="Correlations"),


                column(width=9,

                  tabBox(width=12,#height=400,
                    tabPanel("Cor-Heatmap",
                      fluidRow(
                        column(width=12,
                          uiOutput("show_clickmessage_cor")
                        )
                      ),
                      plotOutput("corProfilesAndBox",click="cor_click"),
                      htmlOutput("savecorProfilesAndBox")
                    ),
                    tabPanel("Pcor-Heatmap",
                      fluidRow(
                        column(width=12,
                          uiOutput("show_clickmessage_pcor")
                        )
                      ),
                      plotOutput("pcorProfilesAndBox",click="pcor_click"),
                      htmlOutput("savepcorProfilesAndBox")
                    )            
                  )
                ),
                column(width=3,
                  uiOutput("showCorProfileAndBox_isdensityOptions"),
                  uiOutput("showCorProfileAndBox_corMethodOptions"),
                  uiOutput("showCorProfileAndBox_logOptions")
                )


              ),



              box(width=5,solidHeader = TRUE,status = "primary",collapsible = TRUE,
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
        column(width=3,style='padding:0px;',
          box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
            title=boxHelp(ID="msg_dynamicsOnGenes_parameters",title="Parameters"),

            uiOutput("show_genelistsforDynamics"),

            
            uiOutput("showROItriadGeneList"),
            HTML("<br>"),
            uiOutput("show_BAMforDynamics"),

            
            HTML("<br>"),
            uiOutput("show_binsforDynamics"), 

            HTML("<br>"),
            uiOutput("show_plotDynamics")
          )

        ),

        column(width=9,style='padding:0px;',
          box(width=12,solidHeader = TRUE,status = "primary",collapsible = TRUE,
                title=boxHelp(ID="msg_dynamicsOnGenes_profiles",title="Metagene profile"), 
            column(width=9,
              #profile
              fluidRow(
                column(width=12,
                    plotOutput("plotProfileDynamics"),
                    fluidRow(
                      column(width=4,
                        htmlOutput("saveprofileDynamics")
                      ),
                      column(width=4,
                        uiOutput("show_chooseMetricforDynamics")
                      )
                    )
                                   
                )
              ),

              HTML('<hr size=3>'),
              #boxplot
              fluidRow(
                column(width=4,
                    plotOutput("plotboxTSSDynamics"),
                    htmlOutput("saveboxTSSDynamics"),
                    htmlOutput("saveboxdatadynamicsTSS")
                ),
                column(width=4,
                    plotOutput("plotboxGBDynamics"),
                    htmlOutput("saveboxGBDynamics"),
                    htmlOutput("saveboxdatadynamicsGB")
                ),
                column(width=4,
                    plotOutput("plotboxTESDynamics"),
                    htmlOutput("saveboxTESDynamics"),
                    htmlOutput("saveboxdatadynamicsTES")
                )
              ),
              HTML('<hr size=3>'),
              #stallingindex
              fluidRow(
                column(width=4,
                    plotOutput("plotSITSSDynamics"),
                    htmlOutput("saveSITSSDynamics")
                ),
                column(width=4,
                    plotOutput("plotSIGBDynamics"),
                    htmlOutput("saveSIGBDynamics")
                ),
                column(width=4,
                    plotOutput("plotSISIDynamics"),
                    htmlOutput("saveSISIDynamics")
                )
                
              ),
              fluidRow(
                column(width=6,
                  uiOutput("show_percentageOutlayerCumulPlots")
                ),
                column(width=6)
              )

            ),
            column(width=3,
              uiOutput("show_islogforDynamics"),
              uiOutput("show_chooseNormalizationforDynamics"),
              uiOutput("show_colorschemeDynamics"),
              uiOutput("show_colorsDynamics")
            )   




          )

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
                HTML("<b>Select the source</b></br>"),
                radioButtons("chooseSourceGO",label=NULL,
                                  choiceNames=list(
                                    htmlhelp("From ROI","help_goAnalysis_parameters_fromROI"),
                                    htmlhelp("From gene list","help_goAnalysis_parameters_fromgenelist")
                                  ),choiceValues=list("fromROI","fromGeneList")
                            ),
                #here, the UI (checkboxGroup if from ROI, textInput if from custom list)
                uiOutput("viewSelectGenesGO"),
                #here, put choice of kind of ID of the genes (symbols, etrez, ensembl): only symbols if database (promoters) is not present
                uiOutput("additionalparametersGO"),
                uiOutput("chooseWindowROIGO"),
                #here, we serve all possible genesets using GenesetsGMT global variable

                uiOutput("show_selectedGenesetsGO"),

                #radiobutton to choose if ranking or clustering the results
                uiOutput("chooseOrderingGO_widget"),
                
                uiOutput("clustertypeGO_widget"),
                uiOutput("clusternumbershowGO_widget"),
                uiOutput("clusterHDistMethodGO_widget"),
                uiOutput("clusterHClustMethodGO_widget"),
                uiOutput("clusterKstartsGO_widget"),
                uiOutput("clusterKiterationsGO_widget"),
                uiOutput("show_minmaxSizeGO"),
            
                uiOutput("show_doTheGO")
              ),

              tabPanel("Filtering",
                #geneRatio:
                sliderInput('quantileGeneRatioGO',label=list("Gene ratio threshold:",htmlhelp("","help_goAnalysis_parameters_generatio")),min = 0, max = 1, value = 0,step=0.05),
                #-log10Padj:
                sliderInput('log10padjGO',label=list("-log10 padj threshold:",htmlhelp("","help_goAnalysis_parameters_padjthresh")),min = 1, max = 50, value = 2,step=1),
                #top n statistically significant:
                sliderInput('topNGO',label="Top significant hits:",min = 1, max = 100, value = 10,step=1)
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
                column(width=12,
                  uiOutput("show_clickmessage_GO")
                )
              ),

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
                  #here, sliderInputs of various thresholds for the analyses
                  uiOutput("show_scaleQuantileGO"),
                  #color scale
                  uiOutput("show_colorScaleGO"),
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
    )




  )
)
