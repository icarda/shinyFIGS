extractWCDataUI <- function(id){
  ns <- NS(id)
  tagList(
    selectInput(ns("var"), "Select climatic variable(s)", multiple = TRUE, c("Var" = "", c('average temperature'='tavg', 'minimum temperature'='tmin', 'maximum temperature'='tmax', 'precipitation'='prec', 'solar radiation'='srad', 'water vapor pressure'='vapr', 'wind speed'='wind', 'Bioclimatic variables'='bio'))),
    uiOutput(ns("bio_vars"))
  )
}

extractWCDataMod <- function(input, output, session, rv){
  
  WCdata <- reactive({
    req(rv$datasetInput)
    if(is.null(input$var)){
      showNotification("Please choose at least one climate variable", type = "warning", duration = 5)
      return()
    }
    else{
      df_cleaned <- rv$datasetInput %>%
        dplyr::mutate(
          !!rlang::sym(rv$lng) := as.numeric(!!rlang::sym(rv$lng)), 
          !!rlang::sym(rv$lat) := as.numeric(!!rlang::sym(rv$lat))
        ) %>%
        dplyr::filter(!is.na(!!rlang::sym(rv$lng)) & !is.na(!!rlang::sym(rv$lat)))
      df_cleaned <- df_cleaned %>%
        filter(PopulationType!="Genetic stock" & PopulationType!="Unreleased breeding material" & PopulationType!="Research material")
      withProgress(message = "Extracting World Clim Data ...", {
        WCdata <- icardaFIGSr::extractWCdata(df_cleaned, long = rv$lng, lat = rv$lat, var = input$var, res = 2.5)
      })
      return(WCdata)
    }
  })
  
  output$bio_vars <- renderUI({
    if("bio" %in% input$var){
      includeMarkdown("Rmd/bio_vars.Rmd")
    }
  })
  
  return(WCdata)
}

