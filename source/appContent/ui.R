### global UI file for the main structure of the GUI

x=rownames(installed.packages())


####################################################################################
####################################################################################
####################################################################################
# HEADER of the dashboard (with name of the program)
####################################################################################
####################################################################################
####################################################################################

header <- dashboardHeader(title = "ChroKit",disable = FALSE,titleWidth = 180,
    #tags$li(class = "dropdown",

            tags$li(class = "dropdown", htmlOutput("showcurrentASSEMBLY"),style = "padding-top: 10px; padding-bottom: 4px; padding-right: 15px; color: #fff; font-size:100%;"),
            tags$li(class = "dropdown", actionButton("loadExampleData", "Load example data",style='padding:4px; font-size:80%'),style = "padding-top: 8px; padding-bottom: 0px; padding-right: 10px; color: #fff; font-size:100%;"),
            tags$li(class = "dropdown", actionButton("gototutorial", "Go to tutorial",style='padding:4px; font-size:80%',onclick ="window.open('https://ocroci.github.io/ChroKit/', '_blank')"),style = "padding-top: 8px; padding-bottom: 0px; padding-right: 10px; color: #fff; font-size:100%;"),
            tags$li(class = "dropdown", actionButton("show_overview", "Overview",style='padding:4px; font-size:80%'),style = "padding-top: 8px; padding-bottom: 0px; padding-right: 10px; color: #fff; font-size:100%;")
            
            )
             
          #)

####################################################################################
####################################################################################
####################################################################################
# SIDEBAR of the dashboard (contains the main menus)
####################################################################################
####################################################################################
####################################################################################

sidebar<- dashboardSidebar(
   width = 180,
  sidebarMenu(#style = "position: fixed; overflow: visible;",
    style="position: fixed; height: 90vh; overflow-y: auto;",
    HTML("&nbsp&nbsp;&nbsp;&nbsp© Ottavio Croci<br>"),
    #actionButton("loadExampleData", "Load example data"),
    #HTML("<br>"),
    HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp<b>1) Import data</b>"),
    menuItem("ROIs", tabName = "BEDblock", icon = icon("fas fa-file-text")),
    menuItem("Enrichment files", tabName = "BAMblock", icon = icon("fas fa-file")),
    menuItem("Assembly", tabName = "TXDBblock", icon = icon("database")),
    shinyFilesButton('loadenv', label='Load session file', 'Select rds session file to load', icon=icon("fas fa-file-export"),FALSE),
    #menuItem("Gene lists", tabName = "GENELISTSblock", icon = icon("list")),
    HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp<b>2) Data management</b>"),
    menuItem("ROI preparation",tabName="MANIPULATEROIblock",icon=icon("fas fa-hammer")),
    #menuItem("Associate enrichments", tabName = "ASSOCIATEblock", icon = icon("fas fa-paperclip")),
    #menuItem("ROI management", tabName = "ROIblock", icon = icon("sort-amount-down-alt")),
    HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp<b>3) Data visualization</b>"),
    menuItem("Genomics", tabName = "GENOMICSblock", icon = icon("chart-area")),
    #HTML("<br>"),
    HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp<b>Save your progress</b>"),
    shinySaveButton("saveWork",label="Save session","Save working environment in rds file...",icon=icon("fas fa-file-import"),filetype=list(rds="rds")),
    #menuItem("Save/Load", tabName = "SAVELOADblock", icon = icon("save")),
    HTML("<br>"),
    .busyIndicator(text="Loading..." , wait=1000 , image="gif.gif"),
    htmlOutput("showRAMusageGC"),
    plotOutput("showRAMbar",height="20",width="120"),
    #HTML("<br>"),
    
    menuItem("Log messages", tabName = "LOGSblock",icon = icon("fas fa-book"))
   

  )      
)

####################################################################################
####################################################################################
####################################################################################
# BODY of the dashboard
####################################################################################
####################################################################################
####################################################################################
### each of these will execute the corresponding R file


####################################################################################
# #TAb of managing BEDs
####################################################################################
source(file.path("uis","uiBED.R"),local=TRUE)$value

####################################################################################
# #TAb of managing BAMs
####################################################################################
source(file.path("uis","uiBAM.R"),local=TRUE)$value

####################################################################################
# #TAb of managing TXDB databases
####################################################################################
source(file.path("uis","uiTXDB.R"),local=TRUE)$value

####################################################################################
# #TAb of manipulating ROIs
####################################################################################
source(file.path("uis","uiMANIPULATEROI.R"),local=TRUE)$value


####################################################################################
# #TAb of enrichments association
####################################################################################
#source(file.path("uis","uiASSOCIATE.R"),local=TRUE)$value


####################################################################################
# #TAb of managing ROIs
####################################################################################
#source(file.path("uis","uiROI.R"),local=TRUE)$value

####################################################################################
# #TAb of genomics
####################################################################################
source(file.path("uis","uiGENOMICS.R"),local=TRUE)$value

####################################################################################
# #TAb of save/load
####################################################################################
source(file.path("uis","uiLOGS.R"),local=TRUE)$value



####################################################################################
# put all tab together in the body. 1:1 correspondence with sidebar!
####################################################################################

body<-dashboardBody(

  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  tabItems(tabBED,tabBAM,tabTXDB,tabMANIPULATEROI,tabGENOMICS,tabLOGS)
  #everything in common, to show in all the tabs (example: log strings) 
  #logs are now in uis/uiBED.R

)

###finally, build the dashboard with all the elements
dashboardPage(header, sidebar, body)


