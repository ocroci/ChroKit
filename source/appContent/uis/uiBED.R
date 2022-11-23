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
      uiOutput("boxdeleteroi")     
    )   
  ),


  fluidRow(
    column(width=6,style='padding:0px;',
      uiOutput("boxviewroi")
    ),





    #download the ROI table viewing interactively all the infos
    column(width=6,style='padding:0px;',
      uiOutput("boxgetroi")
    )
  )

        


          
)
