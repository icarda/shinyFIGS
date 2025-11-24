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

#' @title Extracting historical climatic data from WorldClim 2.1
#' @description Return a data frame based on specified climatic variables. 
#' @param sites object of class "data.frame" with coordinates of sites from which to extract data.
#' @param long character. Name of column from \code{sites} with longitude.
#' @param lat character. Name of column from \code{sites} with latitude.
#' @param res numeric. Spatial resolution. Default 2.5
#' @param var character. Climatic variable(s) to be extracted: 'tavg', 'tmin', 'tmax', 'prec', 'bio', 'srad', 'vapr', 'wind'
#' @return An object of class "data.frame" with specified climatic variables for coordinates in \code{sites}.

extractWCdata <- function(sites, long, lat, var, res = 2.5){
  #remove records having NA coordinates
  out <- list(
    is.na(sites[[long]]), 
    is.na(sites[[lat]])
  )
  
  outReduced <- !Reduce("|", out)
  sites <- sites[outReduced,]
  xy <- cbind(sites[[long]], sites[[lat]])
  sp <- sp::SpatialPoints(xy)
  
  for (ivar in var){
    rasterfile <- .getRasterData(var = ivar, res = res)

    for (i in 1:length(names(rasterfile))){
      f.name <- names(rasterfile)[i]
      var.name <- sub(paste(".*",res,"m_", sep = ''), "", f.name)
      print(var.name)
      sites[, var.name] <- raster::extract(rasterfile[[i]], sp, method = 'simple')
      sites[, var.name] <- round(sites[, var.name], 3)
    }
  }
  
  # test for parallel processing ...
  # for (ivar in var){
  #   rasterfile <- .getRasterData(var = ivar, res = res)
  #   n_cores <- detectCores()
  #   cluster <- makeCluster(8)
  #   registerDoParallel(cluster)
  #   clim_vals <- foreach(i=1:length(names(rasterfile)), 
  #                        .packages = c("terra","raster"),
  #                        .combine='cbind') %dopar% {
  #                          terra::extract(rasterfile[[i]], sp)
  #                        }
  #   var.name <- foreach(f.name = names(rasterfile), .combine = "c") %dopar% {
  #     sub(paste(".*","2.5m_", sep = ''), "", f.name)
  #   }
  #   parallel::stopCluster(cluster)
  #   colnames(clim_vals) <- var.name
  #   sites <- cbind(sites, clim_vals)
  # }
  
  return(sites)
}

#'
#'
#'
#'
#'

.getRasterData <- function(var, res){
  
  stopifnot(var %in% c('tavg', 'tmin', 'tmax', 'prec', 'bio', 'srad', 'vapr', 'wind'))
  
  path <- getwd()
  path <- paste(path, '/WorldClim_', res, '/', sep='')
  dir.create(path, showWarnings=FALSE)
  
  theurl = paste("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_", res,"m_", var,".zip", sep='')
  
  zip <- paste('wc2.1_', res, 'm_', var ,'.zip', sep='')
  zipfile <- paste(path, zip, sep='')
  
  if (var  != 'bio') {
    tiffiles <- paste('wc2.1_', res, 'm_', var, '_', sprintf("%02d",1:12), '.tif', sep='')
  } else {
    tiffiles <- paste('wc2.1_', res, 'm_', var,'_', 1:19, '.tif', sep='')	
  }
  
  files <- paste(path, tiffiles, sep='')
  fc <- sum(file.exists(files))
  
  if ( fc < length(files) ) {
    if (!file.exists(zipfile)) {
      .download(theurl, zipfile)
      if (!file.exists(zipfile))	{ 
        message("\n Could not download file -- perhaps it does not exist") 
      }
    }	
    utils::unzip(zipfile, exdir=dirname(zipfile))
  }
  
  st <- raster::stack(files)
  
  raster::projection(st) <- "+proj=longlat +datum=WGS84"
  return(st)
}

#'
#'
#'
#'
#'

.download <- function(url, filename) {
  fn <- paste(tempfile(), '.download', sep='')
  res <- utils::download.file(url=url, destfile=fn, method="auto", quiet = FALSE, mode = "wb", cacheOK = TRUE)
  if (res == 0) {
    w <- getOption('warn')
    on.exit(options('warn' = w))
    options('warn'=-1) 
    if (! file.rename(fn, filename) ) { 
      file.copy(fn, filename)
      file.remove(fn)
    }
  } else {
    stop('could not download the file' )
  }
}


#'Map accessions using leaflet package
#'
#'
#'

mapAccessions <- function(map, df, long, lat, y){
  
  if(y == "None"){
    leaflet_map <- map %>% clearMarkers() %>%
      clearControls() %>% removeLayersControl() %>%
      leaflet::addCircleMarkers(data = df, lng = df[[long]], lat = df[[lat]],
                       color = "#2d7436",
                       radius = 1.5,
                       fill = TRUE,
                       fillColor = "#2d7436",
                       fillOpacity = 1, stroke = TRUE, weight = 0.1) 
  }
  else {
    ## omit NAs in y column
    df.na.omit <- df[!is.na(df[[y]]), ]

    if (is.numeric(df[[y]])){
      pal <- leaflet::colorNumeric(
        palette = c("viridis"),
        domain = df[[y]],
        na.color = "#808080"
      )
    }
    else {
      pal <- leaflet::colorFactor(
        palette = c("viridis"),
        domain = df[[y]],
        na.color = "#808080"
      )
    }
    
    leaflet_map <- map %>% clearMarkers() %>%
      clearControls() %>% removeLayersControl() %>%
      leaflet::addCircleMarkers(data = df.na.omit,
                                lng = df.na.omit[[long]],
                                lat = df.na.omit[[lat]],
                                color = "black",
                                radius = 1.5,
                                fill = TRUE,
                                fillColor = ~pal(df.na.omit[[y]]),
                                label = ~df.na.omit[[y]],
                                stroke = TRUE,
                                fillOpacity = 1,
                                weight = 0.1,
                                group = "withoutNAs") %>%
      leaflet::addCircleMarkers(data = df, lng = df[[long]], lat = df[[lat]],
                                color = "black",
                                radius = 1.5,
                                fill = TRUE,
                                fillColor = ~pal(df[[y]]),
                                label = ~df[[y]],
                                fillOpacity = 1, stroke = TRUE, weight = 0.1, group = "withNAs") %>%
      leaflet::addLegend("bottomright", pal = pal, values = df[[y]], opacity = 1,  title = y) %>%
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
  df1$Aggregated <- "Overall Data"
  df2$Aggregated <- type
  d <- rbind(df1, df2)
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
      n_igs = dplyr::n_distinct(IG),
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

traitSummary <- function(df, ig_trait = NULL) {
  trait <- names(df)[12]
  safe_trait <- gsub("[^A-Za-z0-9]", "_", trait) #remove special characters
  
  trait_summary <- df %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      mean = mean(.data[[trait]], na.rm = TRUE),
      sd  = sd(.data[[trait]], na.rm = TRUE),
      min  = min(.data[[trait]], na.rm = TRUE),
      max = max(.data[[trait]], na.rm = TRUE),
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
  maxOutliers <- df[df[[12]] > (mean(df[[12]]) + 3 * (sd(df[[12]]))), ]
  maxOutliers <- maxOutliers %>% dplyr::select(c(1,2,3,4,5,12))
  
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
  minOutliers <- df[df[[12]] < (mean(df[[12]]) - 3 * (sd(df[[12]]))), ]
  minOutliers <- minOutliers %>% dplyr::select(c(1,2,3,4,5,12))
  
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
  histPlot <- ggplot(df, aes(x = .data[[trait]])) +
    geom_histogram(fill = "lightblue", bins = 30, color = "black") +
    geom_vline(aes(xintercept = mean(.data[[trait]], na.rm = TRUE)),
               color = "red", linewidth = 1, linetype = "dashed") +
    geom_vline(aes(xintercept = median(.data[[trait]], na.rm = TRUE)),
               color = "darkgreen", linewidth = 1, linetype = "dotdash") +
    labs(title = paste("Distribution of", trait),
         x = trait, y = "Count",
         caption = "Red = Mean | Green = Median") +
    theme_minimal()
  
  # Boxplot of trait values
  boxPlot <- ggplot(df, aes(x = factor(YEAR), y = .data[[trait]])) +
    geom_boxplot(fill = "lightgreen") +
    geom_point(aes(text = paste("IG:", IG,
                                "<br>Year:", YEAR,
                                "<br>", trait, ":", .data[[trait]])),
               alpha = 0, size = 2, stroke = 0) +
    labs(title = paste(trait, "Distribution by Year"),
         x = "Year", y = trait) +
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
  trait <- names(df)[12]
  
  # Histogram of trait values
  histPlot <- ggplot(df, aes(x = .data[[trait]])) +
    geom_bar(aes(text = after_stat(count)),
             fill = "#FFD580", color ="black") +
    labs(title = paste("Distribution of", trait), x = paste(trait, "Category"), y = "Count") +
    theme_minimal()
  
  fctrs <- unlist(factor_trait_info$valid_options[factor_trait_info$Trait == traitName])
  df$Value <- ifelse(as.character(df[[12]]) %in% as.character(fctrs),
                     "Available", "Not Available")
  
  # Proportions by group
  # Create a df for proportions
  var_sym <- rlang::sym(trait)
  
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
                                    fill = factor(.data[[trait]]))) +
    geom_col(
      aes(
        text = tooltip_text,
        color = Value,
        linewidth = Value
      ),
      position = "fill"
    ) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Set3")) +
    scale_color_manual(values = c("Not Available" = "black", "Available" = NA), guide = "none") +
    scale_linewidth_manual(values = c("Not Available" = 0.4, "Available" = 0), guide = "none") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = paste(trait, "distribution by Year"), x = "Year", y = "Proportion", fill = trait) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  # Counts by group
  countsGroup <- ggplot(df, aes(x = factor(YEAR), fill = factor(.data[[trait]]), 
                                color = Value, linewidth = Value)) +
    geom_bar( aes(text = after_stat(count)),
      position = position_dodge2(preserve = "single")) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Set3")) +
    scale_color_manual(values = c("Not Available" = "black", "Available" = NA), 
                       guide = "none") +
    scale_linewidth_manual(values = c("Not Available" = 0.4, "Available" = 0), 
                           guide = "none") +
    labs(title = paste(trait, "counts by Year"), x = "Year", y = "Count", fill = trait) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  return(list(histogram = ggplotly(histPlot, tooltip = c("x", "text")),
              proportionsByGroup = ggplotly(propGroup, tooltip = c("fill", "text")),
              countsByGroup = ggplotly(countsGroup, tooltip = c("fill", "text"))))
}

traitSummaryFFil <- function(df, traitNumber){
  trait <- names(df)[12]
  
  # Histogram of trait values
  histPlot <- ggplot(df, aes(x = .data[[trait]])) +
    geom_bar(fill = "lightblue", color = "black") +
    labs(title = paste("Distribution of", trait), x = "EGV category", y = "Count") +
    theme_minimal()
  
  # Proportions by group
  propGroup <- ggplot(df, aes(x = factor(YEAR), fill = factor(.data[[trait]]))) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Set3")) +
    scale_y_continuous(labels = percent_format()) +
    labs(title = paste(trait, " distribution by Year"), x = "Year", y = "Proportion", 
         fill = trait) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  # Counts by group
  countsGroup <- ggplot(df, aes(x = factor(YEAR), fill = factor(.data[[trait]]))) +
    geom_bar(position = position_dodge2(preserve = "single")) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Set3")) +
    labs(title = paste(trait, "counts by Year"), x = "Year", y = "Count", fill = trait) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  
  return(list(histogram = ggplotly(histPlot, tooltip = c("x")),
              proportionsByGroup = ggplotly(propGroup, tooltip = c("fill")),
              countsByGroup = ggplotly(countsGroup, tooltip = c("fill"))))
}

