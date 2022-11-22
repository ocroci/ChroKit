####################################################################################
# TAB MANIPULATE ROI
####################################################################################

tabMANIPULATEROI <- tabItem(tabName = "MANIPULATEROIblock",
  
  fluidRow(
    box(width=12,collapsible = TRUE,status = "primary",solidHeader = TRUE,
                  title=boxHelp(ID="msg_ROImanipulation",title="ROI manipulation"),
      #main menu: different manipulation pipelines: automatic or manual   
      fluidRow(  
        column(width=8,
          #how to include an icon into radiobutton. the png must be into www/ directory
          radioButtons("choose_ROImanipulation_pipeline",label="Select how to prepare a ROI",
                                  choiceNames=list(
                                    htmlhelp("Prepare ROI basic (ULTRA-EASY)","msg_prepare_ultraeasy"),
                                    htmlhelp("Prepare ROI for heatmaps (EASY)","msg_prepare_forheat"),#&nbsp&nbsp<img src='resizeico.png' alt='Res' width='30' height='10'>"),
                                    htmlhelp("Prepare genelists for metagene profile (EASY)","msg_prepare_formetagene"),
                                    htmlhelp("Manual ROI management (ADVANCED USERS)","msg_prepare_manual")
                                  ),
                                  choiceValues=list(
                                    "prepare_easy",
                                    "prepare_heat",
                                    "prepare_metagene",
                                    "prepare_manual"
                                  ),
                                    selected=character(0)),
        ),
        column(width=4)
      ),
    ),
    uiOutput("show_ROImanipulationUI")
 
  )



          
)