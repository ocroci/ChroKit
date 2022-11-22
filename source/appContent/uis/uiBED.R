####################################################################################
# TAB BED files
####################################################################################

tabBED <- tabItem(tabName = "BEDblock",

  fluidRow(  
    column(width=8,style='padding:0px;',
      box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
              title=boxHelp(ID="msg_coordinateFiles_chooseCoordinates",title="New ROI"),
        #how to include an icon into radiobutton. the png must be into www/ directory
        radioButtons("importROImainchoice",label="Select how to get or create a new ROI",
                                choiceNames=list(
                                  htmlhelp("Import genomic coordinates (bed/gtf/gff files)","help_BED_fromfiles"),#&nbsp&nbsp<img src='resizeico.png' alt='Res' width='30' height='10'>"),
                                  htmlhelp("Get promoters, transcripts, TES coordinates of a list of genes","help_BED_fromgenelist"),
                                  htmlhelp("Generate ROI from a sequence pattern in the genome","help_BED_frompatterngenome")
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
    ),

    column(width=4,style='padding:0px;',
      box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
        title=boxHelp(ID="msg_deleteRois_deleteRois",title="Loaded ROIs"),

        HTML("<b>Available ROIs:</b>"),
        wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 300px; background-color: #ffffff;",
          checkboxGroupInput("selectedCustomROItoRemove",NULL,NULL)
        ),
        HTML("Select ROIs to delete<br><br>"),
        actionButton("deleteROI", "Delete")
      )     
    )   
  ),


  fluidRow(
    column(width=6,style='padding:0px;',
      box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
        title=boxHelp(ID="msg_quickviewROIs",title="Quick ROI preview"),
        fluidRow(
          column(width=4,
            HTML("<b>Select ROI to view:</b>"),
            wellPanel(id = "logPanel",style = "overflow-y:scroll; overflow-x:scroll; max-height: 300px; max-width: 300px; background-color: #ffffff;",
              checkboxGroupInput("confirmviewROI", label=NULL,choices=NULL)
            )
          ),

          column(width=8,
            uiOutput("show_chooseROIvisualiz"),
            plotOutput('viewROImaterial'),
            htmlOutput("saveviewpeaksROImaterial")


          )
        )
      )
    ),





    #download the ROI table viewing interactively all the infos
    column(width=6,style='padding:0px;',
      box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
        title=boxHelp(ID="msg_getRois_BOX",title="Download ROI"),
        fluidRow(
          column(width=4,
            HTML("<b>Select ROI:</b>"),
            selectInput("listgetROI",NULL,choices=NULL),
            radioButtons("choosegetROItype",label=NULL,
                              choiceNames=list(
                                htmlhelp("Features for each genomic range","help_BED_getroi_eachGR"),
                                htmlhelp("Gene list inside genomic window","help_BED_getroi_genomicWindow")
                              ),
                              choiceValues=list(
                                "eachRange",
                                "genesWindow"
                              ),selected="eachRange"),   
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
