#' @title Getting List of Crops Available in ICARDA's Genebank Documentation System
#' @description Return a list with codes and names of available crops.
#' @return A list containing all crops available in ICARDA's Genebank Documentation System.
#' @details The crop codes and names are fetched from ICARDA's online server.

getCrops <- function() {
  result <- read.csv("https://grs.icarda.org/web_services/getCrops.php")
  df <- data.frame(colnames(result)[1],colnames(result)[2])
  names(df) <- c("CropCode","CropName")
  names(result) <- c("CropCode","CropName")
  crops <- rbind(df,result)
  return(crops)
}

#'Map accessions using leaflet package
#'
#'
#'

mapAccessions <- function(map, df, long, lat, y){
  
  df_clean <- df %>%
    dplyr::mutate(
      lng = as.numeric(!!rlang::sym(long)), 
      lat = as.numeric(!!rlang::sym(lat))
    ) %>%
    dplyr::filter(!is.na(lng) & !is.na(lat))
  
  if(y == "None"){
    leaflet_map <- map %>% clearMarkers() %>%
      clearControls() %>% removeLayersControl() %>%
      leaflet::addCircleMarkers(data = df_clean, lng = ~lng, lat = ~lat,
                       color = "#2d7436",
                       radius = 1.5,
                       fill = TRUE,
                       fillColor = "#2d7436",
                       fillOpacity = 1, stroke = TRUE, weight = 0.1) 
  }
  else {
    ## omit NAs in y column
    df.na.omit <- df_clean[!is.na(df_clean[[y]]), ]

    if (is.numeric(df_clean[[y]])){
      pal <- leaflet::colorNumeric(
        palette = c("viridis"),
        domain = df_clean[[y]],
        na.color = "#808080"
      )
    }
    else {
      pal <- leaflet::colorFactor(
        palette = c("viridis"),
        domain = df_clean[[y]],
        na.color = "#808080"
      )
    }
    
    leaflet_map <- map %>% clearMarkers() %>%
      clearControls() %>% removeLayersControl() %>%
      leaflet::addCircleMarkers(data = df.na.omit,
                                lng = ~lng,
                                lat = ~lat,
                                color = "black",
                                radius = 1.5,
                                fill = TRUE,
                                fillColor = ~pal(df.na.omit[[y]]),
                                label = ~df.na.omit[[y]],
                                stroke = TRUE,
                                fillOpacity = 1,
                                weight = 0.1,
                                group = "withoutNAs") %>%
      leaflet::addCircleMarkers(data = df_clean, lng = ~lng, lat = ~lat,
                                color = "black",
                                radius = 1.5,
                                fill = TRUE,
                                fillColor = ~pal(df_clean[[y]]),
                                label = ~df_clean[[y]],
                                fillOpacity = 1, stroke = TRUE, weight = 0.1, group = "withNAs") %>%
      leaflet::addLegend("bottomright", pal = pal, values = df_clean[[y]], opacity = 1,  title = y) %>%
      addLayersControl(baseGroups = c("withNAs","withoutNAs"),
                       options = layersControlOptions(collapsed = FALSE))
  }
  leaflet_map
}

#' this function gets the country ISO3 code giving the country name
#'
#'
#'

countryCode <- function(countryName){
  if(is.null(countryName)) code <- NULL
  else code <- countries %>% distinct() %>% filter(ADMIN %in% countryName) %>% pull(ADM0_A3)
  return(code)
}

#' search a given pattern in an object
#'
#'
#'

search4pattern <- function(pattern, obj){
  matchItems <- grep(pattern = paste0(pattern, collapse = "|"), x = obj, ignore.case = T, value = T)
  return(matchItems)
}

#' Map two datasets using leaflet package
#'
#'
#'
#'

map_two_dfs <- function(map, df1, df2, lng, lat, type){
  df1_clean <- df1 %>%
    dplyr::mutate(
      lng = as.numeric(!!rlang::sym(lng)), 
      lat = as.numeric(!!rlang::sym(lat))
    ) %>%
    dplyr::filter(!is.na(lng) & !is.na(lat))
  
  df2_clean <- df2 %>%
    dplyr::mutate(
      lng = as.numeric(!!rlang::sym(lng)), 
      lat = as.numeric(!!rlang::sym(lat))
    ) %>%
    dplyr::filter(!is.na(lng) & !is.na(lat))
  
  df1_clean$Aggregated <- "Overall Data"
  df2_clean$Aggregated <- type
  d <- rbind(df1_clean, df2_clean)
  d_overall <- subset(d, Aggregated == "Overall Data")
  d_subset <- subset(d, Aggregated == type)
  
  pal <- leaflet::colorFactor(c("#2d7436", "#ED7506"), domain = c("Overall Data", type))
  
  map %>% leaflet::clearMarkers() %>%
    leaflet::clearControls() %>% leaflet::removeLayersControl() %>%
    leaflet::addMapPane(type, zIndex = 450) %>%
    leaflet::addMapPane("Overall Data", zIndex = 410) %>%
    leaflet::addCircleMarkers(data = d_overall, lng = d_overall[[lng]], lat = d_overall[[lat]],
                     radius = 2,
                     color = "black",
                     fill = TRUE,
                     fillColor = "#2d7436",
                     fillOpacity = 0.8, stroke = TRUE, weight = 0.1, group = "Overall Data") %>% 
      leaflet::addCircleMarkers(data = d_subset, lng = d_subset[[lng]], lat = d_subset[[lat]],
                     radius = 2,
                     color = "black",
                     fillColor = "#ED7506",
                     stroke = TRUE, fillOpacity = 0.8, weight = 0.1, group = type) %>%
      leaflet::addLayersControl(overlayGroups = c("Overall Data", type),
                     options = layersControlOptions(collapsed = FALSE)) 
  # %>%
  #   leaflet::addLegend(pal = pal, values = d[["Aggregated"]], opacity = 1,  title = '')
}

#' get the value that has maximun occurences in a vector, or choose a value randomly if occurences are equivalent
#'
#'
#'

max.frequency <- function(x){
  t <- table(x)
  elts <- names(t)[t == max(t)]
  if(length(elts) > 1) max.frq <- sample(elts, 1)
  else max.frq = elts
  return(max.frq)
}

#'
#'
#'

create_buttons <- function(vars) {
  lapply(
    vars,
    FUN = function(var) {
      button <- list(
        method = 'restyle',
        args = list('transforms[0].value', var),
        label = var
      )
    }
  )
}

#'Generate the dropdown list from variables for X axis
#'
#'
x_axis_dd_list <- function(df, vars) {
  lapply(
    vars,
    FUN = function(var) {
      button <- list(
        method = 'update',
        args = list(list(x = list(df[[var]]))),
        label = var
      )
    }
  )
}

make_sliders <- function(x, var) {
    rng <- range(x, na.rm = TRUE)
    sliderInput(var, var, min = floor(rng[1]), max = ceiling(rng[2]), 
                value = c(floor(rng[1]),ceiling(rng[2])))
}

render_hists <- function(output, x, var){
  plotname <- paste("histo", var, sep = "")
  output[[plotname]] <- plotly::renderPlotly({
    plotly::plot_ly(x, x = x[[var]], color = "#ff8103") %>% 
    plotly::add_histogram() %>%
    plotly::add_annotations(var, x = 0.5, y = 1, 
      xref = "paper", yref = "paper", showarrow = FALSE)
    })
}

render_prints <- function(output, x, var){
  printname <- paste("sumClimVar", var, sep = "")
  output[[printname]] <- renderPrint({
    print(var)
    summary(x[[var]])
  })
}

filter_var <- function(x, val) {
  !is.na(x) & x >= val[1] & x <= val[2]
}

update_slider <- function(session, x, var) {
  newMin <- floor(min(x, na.rm = TRUE))
  newMax <- ceiling(max(x, na.rm = TRUE))
  updateSliderInput(session, inputId = var, value = c(newMin,newMax))
}


#' functions for traits data exploration (by Tamara Ortiz)
#' 

csvDownloadButton <- function(tableId, label = "Download as CSV", filename = "data.csv") {
  htmltools::tags$button(
    label,
    onclick = sprintf("Reactable.downloadDataCSV('%s', '%s')", tableId, filename)
  )
}

# General ----
summaryPerYear <- function(df, accessions){
  summary_year <- df %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      n_igs = dplyr::n_distinct(AccessionNumber),
      coverage = n_igs / nrow(accessions) * 100,
      .groups = "drop"
    )
  plot <- ggplot(summary_year, aes(x = factor(YEAR), y = n_igs)) +
    geom_col(fill = "#ADD0B3") +
    labs(title = "Unique IGs per Year", x = "Year", y = "Number of Unique IGs") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  return(ggplotly(plot, tooltip = c("y")))
}

cummulativePerYear <- function(df, accessions){
  coverage_cum <- df %>%
    dplyr::arrange(YEAR, AccessionNumber) %>%
    dplyr::distinct(YEAR, AccessionNumber) %>%
    dplyr::mutate(first_seen = !duplicated(AccessionNumber)) %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(new_igs = sum(first_seen), .groups = "drop") %>%
    dplyr::mutate(
      cum_igs = cumsum(new_igs),
      coverage_cum = 100 * cum_igs / nrow(accessions)
    )
  
  coverage_cum <- coverage_cum %>%
    dplyr::mutate(YEAR = factor(YEAR))
  
  # as bar
  barPlot <- ggplot(coverage_cum, aes(x = YEAR, y = coverage_cum)) +
    geom_col(fill = "lightblue") +
    geom_text(aes(label = sprintf("%.1f%%", coverage_cum)),
              vjust = -0.3, size = 2
    ) +
    coord_cartesian(ylim = c(0, 105)) +
    labs(
      title = "Cumulative Collection Coverage by Year",
      x = "Year",
      y = "Coverage (%)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  return(ggplotly(barPlot, tooltip = c("y")))
}

# Numeric ----

traitSummary <- function(df, traitName, ig_trait = NULL) {
  #trait <- names(df)[12]
  safe_trait <- gsub("[^A-Za-z0-9]", "_", traitName) #remove special characters
  
  trait_summary <- df %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      mean = mean(.data[[traitName]], na.rm = TRUE),
      sd  = sd(.data[[traitName]], na.rm = TRUE),
      min  = min(.data[[traitName]], na.rm = TRUE),
      max = max(.data[[traitName]], na.rm = TRUE),
      .groups = "drop"
    )
  
  trait_summaryTable <- htmltools::browsable(
    tagList(
      csvDownloadButton(paste0("traitSummary", safe_trait),
                        "Download as CSV",
                        filename = paste0("traitSummary", safe_trait, ".csv")),
      reactable(trait_summary,
                filterable = TRUE,
                searchable = TRUE,
                elementId = paste0("traitSummary", safe_trait))
    )
  )
  
  # Maximum outliers
  maxOutliers <- df[df[[traitName]] > (mean(df[[traitName]]) + 3 * (sd(df[[traitName]]))), ]
  maxOutliers <- maxOutliers %>% dplyr::select(c('AccessionNumber',
                                                 'YEAR',
                                                 'EXPT',
                                                 'PLOT',
                                                 'PlantingDate',
                                                 traitName
                                                 ))
  
  maxOutliersTable <- htmltools::browsable(
    tagList(
      csvDownloadButton(paste0("maxOutliers", safe_trait),
                        "Download as CSV",
                        filename = paste0("maxOutliers", safe_trait, ".csv")),
      reactable(maxOutliers,
                filterable = TRUE,
                searchable = TRUE,
                elementId = paste0("maxOutliers", safe_trait))
    )
  )

  # Minimum outliers
  minOutliers <- df[df[[traitName]] < (mean(df[[traitName]]) - 3 * (sd(df[[traitName]]))), ]
  minOutliers <- minOutliers %>% dplyr::select(c('AccessionNumber',
                                                 'YEAR',
                                                 'EXPT',
                                                 'PLOT',
                                                 'PlantingDate',
                                                 traitName))
  
  minOutliersTable <- htmltools::browsable(
    tagList(
      csvDownloadButton(paste0("minOutliers", safe_trait),
                        "Download as CSV",
                        filename = paste0("minOutliers", safe_trait, ".csv")),
      reactable(minOutliers,
                filterable = TRUE,
                searchable = TRUE,
                elementId = paste0("minOutliers", safe_trait))
    )
  )
  
  # Histogram of trait values
  histPlot <- ggplot(df, aes(x = .data[[traitName]])) +
    geom_histogram(fill = "lightblue", bins = 30, color = "black") +
    geom_vline(aes(xintercept = mean(.data[[traitName]], na.rm = TRUE)),
               color = "red", linewidth = 1, linetype = "dashed") +
    geom_vline(aes(xintercept = median(.data[[traitName]], na.rm = TRUE)),
               color = "darkgreen", linewidth = 1, linetype = "dotdash") +
    labs(title = paste("Distribution of", traitName),
         x = traitName, y = "Count",
         caption = "Red = Mean | Green = Median") +
    theme_minimal()
  
  # Boxplot of trait values
  boxPlot <- ggplot(df, aes(x = factor(YEAR), y = .data[[traitName]])) +
    geom_boxplot(fill = "lightgreen") +
    geom_point(aes(text = paste("IG:", AccessionNumber,
                                "<br>Year:", YEAR,
                                "<br>", traitName, ":", .data[[traitName]])),
               alpha = 0, size = 2, stroke = 0) +
    labs(title = paste(traitName, "Distribution by Year"),
         x = "Year", y = traitName) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  # Convert to interactive plotly
  boxplot_plotly <- ggplotly(boxPlot, tooltip = "text")
  
  # Remove hover info from boxplot layer
  boxplot_plotly <- boxplot_plotly %>%
    plotly::style(hoverinfo = "none", traces = 1)
  
  return(list(
    traitSummary = trait_summaryTable,
    maxOutliers = maxOutliersTable,
    minOutliers = minOutliersTable,
    histogram = ggplotly(histPlot, tooltip = c("x")),
    boxplot = boxplot_plotly
  ))
}

# Factor ----

summaryPerYearF <- function(df, accessions){
  # unique IGs per year + coverage
  summary_year <- df %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      n_igs = dplyr::n_distinct(IG),
      coverage = n_igs / nrow(accessions) * 100,
      .groups = "drop"
    )
  
  # unique IGs per year
  plot <- ggplot(summary_year, aes(x = factor(YEAR), y = n_igs)) +
    geom_col(fill = "steelblue") +
    labs(title = "Unique IGs per Year", x = "Year", y = "Number of Unique IGs") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(ggplotly(plot, tooltip = c("y")))
}

cummulativePerYearF <- function(df, accessions){
  # cumulative coverage by year
  coverage_cum <- df %>%
    dplyr::arrange(YEAR, IG) %>%
    dplyr::distinct(YEAR, IG) %>%
    dplyr::mutate(first_seen = !duplicated(IG)) %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(new_igs = sum(first_seen), .groups = "drop") %>%
    dplyr::mutate(
      cum_igs = cumsum(new_igs),
      coverage_cum = 100 * cum_igs / nrow(accessions)
    )
  coverage_cum <- coverage_cum %>%
    dplyr::mutate(YEAR = factor(YEAR))
  
  # as bar
  barPlot <- ggplot(coverage_cum, aes(x = factor(YEAR))) +
    geom_col(aes(y = coverage_cum), fill = "lightblue") +
    geom_text(
      aes(y = coverage_cum, label = sprintf("%.1f%%", coverage_cum)),
      vjust = -0.3,
      size = 2,
      hjust = 0.5
    ) +
    ylim(0, 105) +
    labs(
      title = "Cumulative Collection Coverage by Year",
      x = "Year", y = "Coverage (%)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(ggplotly(barPlot, tooltip = c("y")))
}

traitSummaryF <- function(df, traitName, factor_trait_info){
  #trait <- names(df)[12]
  
  # Histogram of trait values
  histPlot <- ggplot(df, aes(x = .data[[traitName]])) +
    geom_bar(aes(text = after_stat(count)),
             fill = "#FFD580", color ="black") +
    labs(title = paste("Distribution of", traitName), x = paste(traitName, "Category"), y = "Count") +
    theme_minimal()
  
  fctrs <- unlist(factor_trait_info$valid_options[factor_trait_info$Trait == traitName])
  
  df$Value <- ifelse(map_lgl(str_split(df[[traitName]], pattern = "[;&]+"), 
                              ~ all(str_trim(.x) %in% as.character(fctrs))),
                     "Available", "Not Available")
  
  # Dynamic color palette
  max_colors <- length(levels(df[[traitName]]))
  base_pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
  extended_pal <- base_pal(max_colors)
  
  # Proportions by group
  # Create a df for proportions
  var_sym <- rlang::sym(traitName)
  
  df_props <- df %>%
    dplyr::count(YEAR, !!var_sym, Value, name = "count_segment") %>%
    dplyr::group_by(YEAR) %>%
    dplyr::mutate(
      count_by_year = sum(count_segment)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      pct = count_segment / count_by_year * 100,
      tooltip_text = paste0(round(pct, 1), "% (n=", count_segment, ")")
    )
  
  propGroup <- ggplot(df_props, aes(x = factor(YEAR), 
                                    y = count_segment, 
                                    fill = factor(.data[[traitName]]))) +
    geom_col(
      aes(
        text = tooltip_text,
        color = Value,
        linewidth = Value
      ),
      position = "fill"
    ) +
    scale_fill_manual(values = extended_pal) +
    scale_color_manual(values = c("Not Available" = "black", "Available" = NA), guide = "none") +
    scale_linewidth_manual(values = c("Not Available" = 0.4, "Available" = 0), guide = "none") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = paste(traitName, "distribution by Year"), x = "Year", y = "Proportion", fill = traitName) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  # Counts by group
  countsGroup <- ggplot(df, aes(x = factor(YEAR), fill = factor(.data[[traitName]]), 
                                color = Value, linewidth = Value)) +
    geom_bar( aes(text = after_stat(count)),
      position = position_dodge2(preserve = "single")) +
    scale_fill_manual(values = extended_pal) +
    scale_color_manual(values = c("Not Available" = "black", "Available" = NA), 
                       guide = "none") +
    scale_linewidth_manual(values = c("Not Available" = 0.4, "Available" = 0), 
                           guide = "none") +
    labs(title = paste(traitName, "counts by Year"), x = "Year", y = "Count", fill = traitName) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  return(list(histogram = ggplotly(histPlot, tooltip = c("x", "text")),
              proportionsByGroup = ggplotly(propGroup, tooltip = c("fill", "text")),
              countsByGroup = ggplotly(countsGroup, tooltip = c("fill", "text"))))
}

# Passport completeness ----
has_val <- function(x) !is.na(x) & as.character(x) != ""

pdci_pop_group <- function(x){
  dplyr::case_when(
    x %in% c("Wild","Weedy") ~ "wild",
    x %in% c("Landrace") ~ "landrace",
    x %in% c("Genetic stock","Research material","Unreleased breeding material") ~ "breeding",
    x %in% c("Cultivar") ~ "cultivar",
    TRUE ~ "other"
  )
}

pdci_from_columns <- function(df){
  get <- function(nm) if (nm %in% names(df)) df[[nm]] else rep(NA, nrow(df))
  
  IG  <- get("AccessionNumber")
  Crop <- get("CropName")
  Country <- get("CountryOfOrigin")
  Pop <- get("PopulationType")
  Taxon <- get("TaxonName")
  Lat <- get("Latitude")
  Lon <- get("Longitude")
  Alt <- get("Altitude")
  CollYear <- get("CollectingYear")
  
  # site description proxy (D14/COLLSITE, D29/COLLDESCR proxy-ish)
  Site <- has_val(get("Site")) | has_val(get("SiteCode")) | has_val(get("Province")) | has_val(get("Admin1stLevel"))
  
  popg <- pdci_pop_group(Pop)
  coord_pair <- has_val(Lat) & has_val(Lon)
  
  score <- rep(0, nrow(df))
  
  # Independent of population type
  # D05 GENUS + D06 SPECIES
  # Taxon present -> genus + species present
  score <- score + ifelse(has_val(Taxon), 120 + 80, 0)
  
  # D10 CROPNAME
  score <- score + ifelse(has_val(Crop), 45, 0)
  
  # D20 SAMPSTAT
  score <- score + ifelse(has_val(Pop), 80, 0)
  
  # Depending on population type
  # D13 ORIGCTY
  score <- score + ifelse(
    has_val(Country),
    ifelse(popg %in% c("wild","landrace"), 80, 40),
    0
  )
  
  # D14 COLLSITE (including missing coords)
  score <- score + ifelse(
    Site,
    dplyr::case_when(
      popg == "wild" & !coord_pair ~ 70,
      popg == "wild" ~ 20,
      popg == "landrace" & !coord_pair ~ 45,
      popg == "landrace" ~ 15,
      popg == "other" & !coord_pair ~ 20,
      popg == "other" ~ 10,
      TRUE ~ 0
    ),
    0 # cultivar is 0
  )
  
  # D15 + D16 LAT/LON (both must be present)
  score <- score + ifelse(
    coord_pair,
    dplyr::case_when(
      popg == "wild" ~ 60 + 60,
      popg == "landrace" ~ 40 + 40,
      popg == "other" ~ 15 + 15,
      TRUE ~ 0
    ),
    0 # cultivar is 0
  )
  
  # D17 ELEVATION
  score <- score + ifelse(
    has_val(Alt),
    dplyr::case_when(
      popg == "wild"     ~ 20,
      popg == "landrace" ~ 15,
      popg == "other"    ~ 5,
      TRUE ~ 0
    ),
    0
  )
  
  # D18 COLLDATE (CollectionYear)
  score <- score + ifelse(
    has_val(CollYear),
    dplyr::case_when(
      popg %in% c("wild","landrace") ~ 30,
      popg == "other" ~ 10,
      TRUE ~ 0
    ),
    0
  )
  
  # Final PDCI: sum / 100
  round(score / 100, 2)
}
