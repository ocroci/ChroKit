### global UI file for the main structure of the GUI

x=rownames(installed.packages())

#function for the indicator "Loading...". This function was taken from
#xiaodaigh and dcurrier (https://github.com/AnalytixWare/ShinySky/blob/master/R/busy-indicator.r)
#and corrected. Maybe the real time is inside setInterval function
.busyIndicator <- function(text = "Processing..."
                        , image = "http://i.giphy.com/l3V0EQrPMh1nnfbFe.gif"
                        , wait=1000) {
  tagList(
    singleton(tags$head(
      tags$link(rel = "stylesheet"
        , type = "text/css" 
        ,href = file.path("panel","inst","extdata","busyIndicator.css")
      )))
    ,div(class = "mybusyindicator",p(text),img(src=image))
    ,tags$script(sprintf(
    " setInterval(function(){
       if ($('html').hasClass('shiny-busy')) {
        setTimeout(function() {
          if ($('html').hasClass('shiny-busy')) {
            $('div.mybusyindicator').show()
          }
        }, %d)          
      } else {
        $('div.mybusyindicator').hide()
      }
    },1000)
    ",wait)
    )
  ) 
}

####################################################################################
####################################################################################
####################################################################################
# HEADER of the dashboard (with name of the program)
####################################################################################
####################################################################################
####################################################################################

header <- dashboardHeader(title = "ChroKit",disable = FALSE,
    tags$li(class = "dropdown",
            tags$li(class = "dropdown", htmlOutput("showcurrentASSEMBLY"),style = "padding-top: 10px; padding-bottom: 10px; padding-right: 20px; color: #fff; font-size:130%;")
            )
             
          )

####################################################################################
####################################################################################
####################################################################################
# SIDEBAR of the dashboard (contains the main menus)
####################################################################################
####################################################################################
####################################################################################

sidebar<- dashboardSidebar(
  sidebarMenu(style = "position: fixed; overflow: visible;",
    HTML("&nbsp© Ottavio Croci<br><br>"),
    HTML("&nbsp;&nbsp;&nbsp;&nbsp<b>1) Import data</b>"),
    actionButton("loadExampleData", "Load example data"),
    menuItem("ROIs", tabName = "BEDblock", icon = icon("file-excel-o")),

    menuItem("Enrichment files", tabName = "BAMblock", icon = icon("file-o")),
    menuItem("Assembly", tabName = "TXDBblock", icon = icon("database")),
    shinyFilesButton('loadenv', label='Load session file', 'Select rds session file to load', icon=icon("fas fa-file-export"),FALSE),
    #menuItem("Gene lists", tabName = "GENELISTSblock", icon = icon("list")),
    HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp<b>2) Data management</b>"),
    menuItem("ROI preparation",tabName="MANIPULATEROIblock",icon=icon("fas fa-hammer")),
    #menuItem("Associate enrichments", tabName = "ASSOCIATEblock", icon = icon("fas fa-paperclip")),
    #menuItem("ROI management", tabName = "ROIblock", icon = icon("sort-amount-down-alt")),
    HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp<b>3) Data visualization</b>"),
    menuItem("Genomics", tabName = "GENOMICSblock", icon = icon("chart-area")),
    HTML("<br>"),
    HTML("<br>&nbsp;&nbsp;&nbsp;&nbsp<b>Save your progress</b>"),
    shinySaveButton("saveWork",label="Save session","Save working environment in rds file...",icon=icon("fas fa-file-import"),filetype=list(rds="rds")),
    #menuItem("Save/Load", tabName = "SAVELOADblock", icon = icon("save")),
    HTML("<br><br>"),
    htmlOutput("showRAMusageGC"),
     HTML("<br>"),
    .busyIndicator(text="Loading..." , wait=1000 , image="gif.gif")
   

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
#source(file.path("uis","uiSAVELOAD.R"),local=TRUE)$value



####################################################################################
# put all tab together in the body. 1:1 correspondence with sidebar!
####################################################################################

body<-dashboardBody(

  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  tabItems(tabBED,tabBAM,tabTXDB,tabMANIPULATEROI,tabGENOMICS), 
  #everything in common, to show in all the tabs (example: log strings) 
  HTML("<b><h3>Logs: </h3></b>"),
  wellPanel(id = "logPanel",style = "overflow-y:scroll; max-height: 250px",
      htmlOutput("showlogs")
  )

)

###finally, build the dashboard with all the elements
dashboardPage(header, sidebar, body)


