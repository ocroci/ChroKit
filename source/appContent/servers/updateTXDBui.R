
##update available and not available (missing, for the download)
##menus defined in ui TXDB.R
observe({
	updateSelectInput(session,"searchASSEMBLYforuse",label=NULL,choices=as.list(names(availDB)))
})

observe({
	updateSelectInput(session,"searchASSEMBLYfordownload",label=NULL,choices=as.list(DATABASEvariables$missingASSEMBLIES))
})



##update absolute values of upstream/downstream of TSS/TES if slidebar changed:
# observe({
#   input$upstreamTSS
#   updateNumericInput(session,inputId = 'absoluteFilterUpstreamTSS',label=NULL,min = 0, max = 0, step = 50,value=input$upstreamTSS)
# })
# observe({
#   input$downstreamTSS
#   updateNumericInput(session,inputId = 'absoluteFilterDownstreamTSS',label=NULL,min = 0, max = 0, step = 50,value=input$downstreamTSS)
# })
# observe({
#   input$upstreamTES
#   updateNumericInput(session,inputId = 'absoluteFilterUpstreamTES',label=NULL,min = 0, max = 0, step = 50,value=input$upstreamTES)
# })
# observe({
#   input$downstreamTES
#   updateNumericInput(session,inputId = 'absoluteFilterDownstreamTES',label=NULL,min = 0, max = 0, step = 50,value=input$downstreamTES)
# })




