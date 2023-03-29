####################################################################################
# TAB BED files
####################################################################################

tabBED <- tabItem(tabName = "BEDblock",

  tabsetPanel(id="newROItabset",type="pills",

    tabPanel("Generate new ROIs",
      fluidRow(  
        column(width=12,style='padding:0px;',
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                  title=boxHelp(ID="msg_coordinateFiles_chooseCoordinates",title="New ROI"),
            #how to include an icon into radiobutton. the png must be into www/ directory
            radioButtons("importROImainchoice",label="Select how to get or create a new ROI",
                                    choiceNames=list(
                                      htmlhelp("Import genomic coordinates (bed/gtf/gff files)","help_BED_fromfiles"),#&nbsp&nbsp<img src='resizeico.png' alt='Res' width='30' height='10'>"),
                                      htmlhelp("Get promoters, transcripts, TES coordinates of a list of genes","help_BED_fromgenelist"),
                                      list(htmlhelp("Generate ROI from a sequence pattern in the genome","help_BED_frompatterngenome"),htmlwarning("","help_BED_frompatterngenome_warning"))
                                    ),
                                    choiceValues=list(
                                      "fromfile",
                                      "fromgenelist",
                                      "frompattern"
                                    ),
                                      selected=character(0)),            

            HTML("<br><br>"),
            uiOutput("importROIwindowToShow")
          )
        )#,

  
      )#,


    ),


    tabPanel("Managing ROIs",value="managingROItabPanel",

      fluidRow(
        column(width=9,style='padding:0px;',
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_quickviewROIs",title="Quick ROI preview"),
            fluidRow(
              column(width=4,
                uiOutput("show_confirmviewROI")
                # HTML("<b>Select ROI to view:</b>"),
                # wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 300px; max-width: 300px; background-color: #ffffff;",
                #   checkboxGroupInput("confirmviewROI", label=NULL,choices=NULL)
                # )
              ),

              column(width=8,
                uiOutput("show_chooseROIvisualiz"),
                plotOutput('viewROImaterial'),
                htmlOutput("saveviewpeaksROImaterial")


              )
            )
          )
        ),
        column(width=3,style='padding:0px;',
          box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
            title=boxHelp(ID="msg_deleteRois_deleteRois",title="Loaded ROIs"),

            
            HTML("<b>Available ROIs:</b>"),
            uiOutput("show_selectedCustomROItoRemove"),

            HTML("Select ROIs to delete<br><br>"),
            actionButton("deleteROI", "Delete")
          ) 
        )
      )
    ),


    tabPanel("Download ROIs",value="downloadROItabPanel",
      fluidRow(
        box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
          title=boxHelp(ID="msg_getRois_BOX",title="Download ROI"),
          fluidRow(
            column(width=4,
              uiOutput("show_choosegetROImenu"),
              uiOutput("showROIoptionsToGET")    
            ),
            column(width=8, 
              htmlOutput("previewROItodownload"),
              htmlOutput("previewROItodownloadbutton")
            )
          )
        )
      )

    )

    



  )



          
)







