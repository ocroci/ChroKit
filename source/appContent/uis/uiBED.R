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


    )#,



  )



          
)







