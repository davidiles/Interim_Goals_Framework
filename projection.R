library(bbsBayes2)
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(readxl)
library(ggthemes)
library(ggpubr)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Interim_Goals_Framework")
rm(list=ls())

# ------------------------------------------------------------------------------
# Species names and abbreviations (avian core 2022)
# ------------------------------------------------------------------------------

ac <- read_xlsx("data/1_Avian_Core_20220422_FULLVERSION.xlsx")

# ------------------------------------------------------------------------------
#  Species that have experienced declines; from Marcel and Willow
#  Column definitions:
#      Species: the species name
#      StartDate: the desired "start date" for calculating historical trend; chosen subjectively by Willow
#      Survey: the best source of trend information for this species
# ------------------------------------------------------------------------------

from_Willow <- read_xlsx("data/DecliningSpeciesShortorLongTermTrendForPrediction.xlsx")

# ------------------------------------------------------------------------------
# Fit BBS models for each species; not necessary if annual indices are already provided by analysts
# ------------------------------------------------------------------------------

target_species <- subset(from_Willow, Survey == "BBS")

# for (species_name in target_species$Species){
# 
#   # 4 letter species abbreviation code
#   sp_code <- subset(ac, English_Name == species_name)$Species_ID
# 
#   # Directory where fitted model will be stored
#   filename <- paste0("fitted_models/",sp_code,".RDS")
# 
#   # Fit model for this species, save resulting indices
#   if (!file.exists(filename)){
# 
#     # Stratify data
#     s <- stratify(by = "bbs_usgs", species = species_name)
# 
#     p <- tryCatch(expr = {prepare_data(s)},
#                   error = function(e){NULL})
# 
#     # Fit model & generate indices
#     if (is.null(p)) next
#     md <- prepare_model(p, model = "gamye")
#     m <- run_model(md, save_model = FALSE)
# 
#     ## * THIS IS WHAT YOU WOULD NEED FROM ADAM
#     i <- generate_indices(model_output = m, regions = "country")
# 
#     # National indices
#     samps <- i$samples$country_Canada
# 
#     # Save results
#     sp_results <- list(samps = samps)
#     saveRDS(sp_results, file = filename)
#   }
# }

# ------------------------------------------------------------------------------
# This gigantic clunky function summarizes empirical indices, and performs population 
# projections. Output from this function is used for plotting trajectories and
# summarizing current status, relative to goals
# ------------------------------------------------------------------------------

projection_function <- function(
  
  # Annual indices for the species, which is either: 
  #    1) a matrix of samples from posterior (columns = years, rows = sample from posterior)
  #    2) a vector (where each element should be named, according to the year associated with it) 
  indices = NULL,
  
  # Year in which goals are set
  year_goals_are_set = 2022,
  
  # Target rate of population change, once goal is achieved
  target_trend = 3,    
  
  # Years until target growth rate is reached (leave as NULL if using annual_growth_rate_improvement)
  years_to_target_trend = NULL, 
  
  # Annual improvement in population growth rate (leave as NULL if using years_to_target_trend)
  annual_growth_rate_improvement = NULL,
  
  # Final year of projection
  end_of_projection = 2100,  
  
  # Years for calculating baseline population index (mean across these years)
  baseline_years = seq(1980,1984),    
  
  # How many years back to calculate "current" trend over?
  length_current_trend = 30
){           
  
  if ((is.null(years_to_target_trend) & is.null(annual_growth_rate_improvement)) | (!is.null(years_to_target_trend) & !is.null(annual_growth_rate_improvement))){
    print("Must specify either a time until target is achieved OR annual improvement in growth rate")
    break
  } 
  
  # --------------------------------
  # Matrix or vector of annual indices (columns are years, rows are samples from posterior)
  #    columns must be labeled as year
  # --------------------------------
  
  if (is.matrix(indices)){
    year_seq <- colnames(indices) %>% as.numeric()
  } else if (is.vector(indices)){
    year_seq <- names(indices) %>% as.numeric()
    indices <- matrix(indices, ncol = length(indices), nrow = 1)
    colnames(indices) <- year_seq
  }
  
  # --------------------------------
  # Loop through samples from posterior and perform analysis
  #   1) fit gams through indices to "smooth" the annual variation
  #   2) estimate trend based on smoothed indices
  #   3) calculate future projections starting in year of goal
  #   4) evaluate current status relative to projections
  # --------------------------------
  
  # Properties of data/projection (don't need to be modified)
  final_year_of_data = max(year_seq)
  year_seq_projection <- seq(year_seq[1],end_of_projection)
  
  # Empty objects to store various estimates/summary statistics
  indices_StatusQuo <- indices_Recovery <- indices_gam <- matrix(NA, nrow = nrow(indices), ncol = length(year_seq_projection))
  Index_baseline <- Trend_samples <- c()
  
  # Loop through samples from posterior
  for (i in 1:nrow(indices)){
    
    # ---------------------------------------------
    # Extract estimated historical dynamics from bbs fit
    # ---------------------------------------------
    
    # Annual indices for this draw from the posterior
    i_dat <- data.frame(log_y = log(indices[i,]),Year = year_seq) %>% 
      subset(Year <= final_year_of_data)
    
    # ---------------------------------------------
    # Fit a gam through indices to describe pattern of historical change, removing year-to-year noise
    # ---------------------------------------------
    
    # Choose knots for gam; more knots allows for more wiggliness
    # 5 years is probably a good balance of removing year-to-year fluctuations, while describing patterns of change on a scale that is relevant for conservation
    knot_distance = 5
    knots <- seq(min(year_seq), max(year_seq)+10, knot_distance)
    gam <- gam(log_y~s(Year, bs = 'cs', k = length(knots)),
               knots = list(Year = knots),
               data = i_dat)
    
    # Smoothed predictions from gam; these smoothed indices are used to calculate trends
    i_dat$gam_pred <- predict(gam, newdata = i_dat)
    indices_gam[i,1:length(i_dat$gam_pred)] <- exp(i_dat$gam_pred)
    
    # ---------------------------------------------
    # Estimate baseline abundance as average across relevant years
    # ---------------------------------------------
    
    Index_baseline[i] <- mean(exp(i_dat$gam_pred[i_dat$Year %in% baseline_years])) 
    
    # ---------------------------------------------
    # Trend estimate as as the geometric mean annual growth rate (mean of annual differences in log-scale indices)
    # Note that we could alternatively calculate trend as log-linear slope through smoothed indices, using lm(gam_pred ~ Year, data = i_dat_for_trend)
    #    - this has the advantage of removing the strong dependency on start/end years
    # ---------------------------------------------
    
    i_dat_for_trend <- subset(i_dat, Year %in% seq(year_goals_are_set,year_goals_are_set - length_current_trend))
    
    trend_log <- mean(diff(i_dat_for_trend$gam_pred))
    trend_percent <- 100*(exp(trend_log)-1) # Convert to percent change per year
    
    # Save estimate of trend for this sample from posterior
    Trend_samples <- c(Trend_samples,trend_percent)
    
    # ---------------------------------------------
    # Conduct future projections under two scenarios:
    #    1) status quo / business as usual
    #    2) achieving a specified target trajectory
    # ---------------------------------------------
    
    # Dataframe to store projection information
    projection_i <- data.frame(Year = year_seq_projection,
                               y_Obs = NA,
                               y_StatusQuo = NA,
                               y_Recovery = NA)
    
    # Observed indices
    projection_i$y_Obs[projection_i$Year %in% i_dat$Year] <- exp(i_dat$log_y)
    projection_i$y_gam[projection_i$Year %in% i_dat$Year] <- exp(i_dat$gam_pred)
    projection_i <- subset(projection_i, Year >= (year_goals_are_set - length_current_trend))
    
    # **********************
    # Historical trend
    # **********************
    
    # Fill in StatusQuo and Recovery columns (may not be necessary, but just seems more statisfying to have these filled instead of NA)
    projection_i$y_StatusQuo[1] <- projection_i$y_gam[1]
    projection_i$y_Recovery[1] <- projection_i$y_gam[1]
    
    for (y in seq((min(projection_i$Year)+1),year_goals_are_set)){
      projection_i$y_StatusQuo[projection_i$Year == y] <- projection_i$y_StatusQuo[projection_i$Year == (y-1)] *  exp(trend_log)
      projection_i$y_Recovery[projection_i$Year == y] <- projection_i$y_Recovery[projection_i$Year == (y-1)] *  exp(trend_log)
    }
    
    # **********************
    # Project the StatusQuo scenario into the future, at same rate of change
    # **********************
    
    for (y in seq(year_goals_are_set+1,end_of_projection)) projection_i$y_StatusQuo[projection_i$Year == y] <- projection_i$y_StatusQuo[projection_i$Year == (y-1)] *  exp(trend_log)
    
    # **********************
    # Project a recovery scenario into the future
    # **********************
    
    # If current trend is somehow already higher than the eventual target, start the projection using the target trend (we've already achieved it)
    if (trend_percent > target_trend) trend_percent <- target_trend 
    
    # Sequence of annual growth rates if recovery scenario is defined based on "years until target trend is achieved"
    if (!is.null(years_to_target_trend)) trend_seq <- c(seq(trend_percent,target_trend,length.out = years_to_target_trend+1),rep(target_trend,1000)) # adds a bunch of extra years just to ensure we project far enough; these get trimmed off later
    
    # Sequence of annual growth rates if recovery scenario is defined based on annual increments in growth rate
    if (!is.null(annual_growth_rate_improvement)) trend_seq <- c(seq(trend_percent,target_trend,annual_growth_rate_improvement),rep(target_trend,1000)) # adds a bunch of extra years just to ensure we project far enough; these get trimmed off later
    
    # What years do we actually need to store?
    proj_years <- seq(year_goals_are_set,end_of_projection)
    
    # Trim off unnecessary years from projection (we don't need to store 1000 years into the future)
    trend_seq <- trend_seq[1:length(proj_years)] 
    
    # Convert to log-scale annual change (instead of percent) for easier math
    trend_seq_logscale <- log(trend_seq/100 + 1) 
    
    # Calculate annual indices under recovery scenario, using the sequence of annual growth rates
    for (y in seq(year_goals_are_set+1,end_of_projection)) projection_i$y_Recovery[projection_i$Year == y] <- projection_i$y_Recovery[projection_i$Year == (y-1)] *  exp(trend_seq_logscale[proj_years == y])
    
    # Store annual indices under each projection, for this sample from the posterior 
    indices_StatusQuo[i,which(year_seq_projection %in% projection_i$Year)] <- projection_i$y_StatusQuo
    indices_Recovery[i,which(year_seq_projection %in% projection_i$Year)] <- projection_i$y_Recovery
    indices_gam[i,1:length(i_dat$gam_pred)] <- exp(i_dat$gam_pred)
    
  }
  
  # ******************************************************************
  # The code below summarizes the projections in various ways
  # and assesses current status relative to the projections.
  #   (note: probably very clunky and inefficiently coded!)
  # ******************************************************************
  
  # ----------------------------------------------------
  # Round estimates/projections to 3 decimal places
  # ----------------------------------------------------
  
  indices_gam <- round(indices_gam,3)
  indices_StatusQuo <- round(indices_StatusQuo,3)
  indices_Recovery <- round(indices_Recovery,3)
  
  # ----------------------------------------------------
  # Summarize indices each year
  # ----------------------------------------------------
  
  indices_summarized <- indices %>%
    reshape2::melt() %>%
    group_by(year) %>%
    summarize(Index_q_0.025 = quantile(value,0.025),
              Index = quantile(value,0.5),
              Index_q_0.975 = quantile(value,0.975)) %>%
    dplyr::rename(Year = year)
  
  # ----------------------------------------------------
  # Current status (relative to projections and historical abundance)
  # ----------------------------------------------------
  
  status_year_number <- which(year_seq_projection == final_year_of_data)
  
  Prob_Exceed_StatusQuo <- 100*mean(indices_gam[,status_year_number] > indices_StatusQuo[,status_year_number]) %>% round(3)
  Prob_Exceed_Recovery <- 100*mean(indices_gam[,status_year_number] > indices_Recovery[,status_year_number]) %>% round(3)
  Prob_Exceed_Baseline <- 100*mean(indices_gam[,status_year_number] > indices_gam[,1]) %>% round(3)
  
  # ---------------------------------------------------------------------
  # Calculate indices expressed as percent of baseline level
  # ---------------------------------------------------------------------
  
  percent_of_Baseline_StatusQuo <- indices_StatusQuo
  percent_of_Baseline_Recovery <- indices_Recovery
  
  for (j in 1:ncol(percent_of_Baseline_StatusQuo)){
    percent_of_Baseline_StatusQuo[,j] <- 100* percent_of_Baseline_StatusQuo[,j]/Index_baseline
    percent_of_Baseline_Recovery[,j] <- 100* percent_of_Baseline_Recovery[,j]/Index_baseline
  }
  
  percent_of_Baseline_Obs <- indices
  for (j in 1:ncol(percent_of_Baseline_Obs)) percent_of_Baseline_Obs[,j] <- 100* percent_of_Baseline_Obs[,j]/Index_baseline
  
  # -------------------------------
  # Relevant summaries of indices
  # -------------------------------
  
  # Credible intervals on gam fit to historical indices
  gam_summary <- reshape2::melt(indices_gam) %>%
    rename(samp = Var1, year_number = Var2, gam = value) %>%
    mutate(Year = year_seq_projection[year_number]) %>%
    group_by(Year) %>%
    summarize(gam_med = median(gam, na.rm = TRUE),
              gam_q_0.025 = quantile(gam,0.025, na.rm = TRUE),
              gam_q_0.975 = quantile(gam,0.975, na.rm = TRUE))
  
  # Credible intervals on status quo projection
  StatusQuo_summary <- reshape2::melt(indices_StatusQuo) %>%
    rename(samp = Var1, year_number = Var2, StatusQuo = value) %>%
    mutate(Year = year_seq_projection[year_number]) %>%
    group_by(Year) %>%
    summarize(StatusQuo_med = median(StatusQuo),
              StatusQuo_q_0.025 = quantile(StatusQuo,0.025, na.rm = TRUE),
              StatusQuo_q_0.975 = quantile(StatusQuo,0.975, na.rm = TRUE))
  
  # Credible intervals on recovery projection
  Recovery_summary <- reshape2::melt(indices_Recovery) %>%
    rename(samp = Var1, year_number = Var2, Recovery = value) %>%
    mutate(Year = year_seq_projection[year_number]) %>%
    group_by(Year) %>%
    summarize(Recovery_med = median(Recovery),
              Recovery_q_0.025 = quantile(Recovery,0.025, na.rm = TRUE),
              Recovery_q_0.975 = quantile(Recovery,0.975, na.rm = TRUE))
  
  # Credible intervals on estimates of percent of historical abundance
  percent_change_summary_Obs <- reshape2::melt(percent_of_Baseline_Obs) %>%
    rename(Year = year, Obs = value) %>%
    group_by(Year) %>%
    summarize(Obs_med = median(Obs),
              Obs_q_0.025 = quantile(Obs,0.025, na.rm = TRUE),
              Obs_q_0.975 = quantile(Obs,0.975, na.rm = TRUE))
  
  # Credible intervals on estimates of percent change, under status quo projection
  percent_change_summary_StatusQuo <- reshape2::melt(percent_of_Baseline_StatusQuo) %>%
    rename(samp = Var1, year_number = Var2, StatusQuo = value) %>%
    mutate(Year = year_seq_projection[year_number]) %>%
    group_by(Year) %>%
    summarize(StatusQuo_med = median(StatusQuo),
              StatusQuo_q_0.025 = quantile(StatusQuo,0.025, na.rm = TRUE),
              StatusQuo_q_0.975 = quantile(StatusQuo,0.975, na.rm = TRUE))
  
  # Credible intervals on estimates of percent change, under recovery projection
  percent_change_summary_Recovery <- reshape2::melt(percent_of_Baseline_Recovery) %>%
    rename(samp = Var1, year_number = Var2, Recovery = value) %>%
    mutate(Year = year_seq_projection[year_number]) %>%
    group_by(Year) %>%
    summarize(Recovery_med = median(Recovery),
              Recovery_q_0.025 = quantile(Recovery,0.025, na.rm = TRUE),
              Recovery_q_0.975 = quantile(Recovery,0.975, na.rm = TRUE))
  
  # -----------------------------------------
  # Description of recovery progress, assuming goals are met
  # -----------------------------------------
  
  suppressWarnings({
    Recovery_index_2050 <- subset(Recovery_summary, Year == 2050)
    Recovery_percent_2050 <- subset(percent_change_summary_Recovery, Year == 2050)
    
    # Year in which full recovery is achieved
    Year_of_full_Recovery <- apply(percent_of_Baseline_Recovery[,which(year_seq_projection > year_goals_are_set)],1,function(x) min(which(x >= 100), na.rm = TRUE))
    if (mean(Year_of_full_Recovery) == Inf | mean(Year_of_full_Recovery) == -Inf){
      Year_of_full_Recovery <- rep(NA,3)} else {
        
        Year_of_full_Recovery <- year_seq_projection[year_seq_projection > year_goals_are_set][Year_of_full_Recovery] %>%
          quantile(c(0.025,0.5,0.975))
      }
  })
  
  # ******************************************************************
  # This is a tabular summary of various statistics, which might be useful for
  # summarizing the results of these projections
  # ******************************************************************
  recovery_description <- data.frame(
    
    # Estimates of historical trend
    Trend_q0.025 = quantile(Trend_samples,0.025),
    Trend_q0.500 = quantile(Trend_samples,0.500),
    Trend_q0.975 = quantile(Trend_samples,0.975),
    
    # Current status of the population (as of last survey)
    Current_Prob_Exceed_StatusQuo = Prob_Exceed_StatusQuo,
    Current_Prob_Exceed_Recovery = Prob_Exceed_Recovery,
    Current_Prob_Exceed_Baseline = Prob_Exceed_Baseline,
    
    # The projected population index in 2050 if recovery target is met
    Recovery_Index_2050_q0.025 = Recovery_index_2050$Recovery_q_0.025,
    Recovery_Index_2050_q0.500 = Recovery_index_2050$Recovery_med,
    Recovery_Index_2050_q0.975 = Recovery_index_2050$Recovery_q_0.975,
    
    # The projected population index in 2050, expressed as percent of the historical baseline, if recovery target is met
    Recovery_percent_2050_q0.025 = Recovery_percent_2050$Recovery_q_0.025,
    Recovery_percent_2050_q0.500 = Recovery_percent_2050$Recovery_med,
    Recovery_percent_2050_q0.975 = Recovery_percent_2050$Recovery_q_0.975,
    
    # The projected year in which recovery to 100% of baseline would be achieved
    Year_of_Full_Recovery_q0.025 = Year_of_full_Recovery[1],
    Year_of_Full_Recovery_q0.500 = Year_of_full_Recovery[2],
    Year_of_Full_Recovery_q0.975 = Year_of_full_Recovery[3])
  
  # ******************************************************************
  # Function will return all the objects that were created within it (yikes - lots of objects!)
  # This is because I wasn't sure exactly which objects we'd want to plot
  # ******************************************************************
  
  return(as.list(environment()))
  
}

# ------------------------------------------------------------------------------
# Conduct projections for each species that Marcel and Willow sent
#  Notes: in these analysis, Willow chose a different "length of historical trend" for each species based on subjective decisions
#         the code reads in these 
# ------------------------------------------------------------------------------

# An empty table to store species results
sp_table <- data.frame()

for (species_name in target_species$Species){
  
  print(species_name)
  
  # 4 letter species abbreviation code
  sp_code <- subset(ac, English_Name == species_name)$Species_ID
  
  # File where annual national indices are stored
  filename <- paste0("fitted_models/",sp_code,".RDS")
  
  # Skip species if estimates are unavailable
  if (!file.exists(filename)) next
  
  # Load fitted model results
  sp_results <- readRDS(filename)
  
  # Annual indices for the species
  species_indices <- sp_results$samps
  
  # How far back to calculate trend over?
  start_year_for_trend <- subset(target_species, Species == species_name)$StartDate
  
  # Generate projections
  sp_projection <- projection_function(
    
    indices = species_indices,
    
    # Years for calculating baseline index for 'full recovery' (mean across these years)
    baseline_years = seq(1970,1974),    
    
    # How far back to calculate "current" trend over)
    length_current_trend = 2022-start_year_for_trend,
    
    # Year in which goals were set
    year_goals_are_set = 2022,
    
    # Percent change per year once population reaches its target growth rate
    target_trend = 3,    
    
    # Goal is to increase growth rate by 0.5% per year, until target it reached
    annual_growth_rate_improvement = 0.5,
    
    # Final year of projection
    end_of_projection = 2200
  )
  
  # Summary of projections
  sp_projection$recovery_description
  
  # Fix up the table slightly
  sp_projection$recovery_description <- sp_projection$recovery_description %>%
    mutate(Species = species_name) %>%
    relocate(Species)
  
  # Append to sp_table
  sp_table <- rbind(sp_table,sp_projection$recovery_description)
  
  # ************************************************
  # Plot time series of population indices
  # ************************************************
  
  sp_plot_index <- ggplot()+
    
    # Vertical line showing the year goals were set
    geom_vline(xintercept = sp_projection$year_goals_are_set, size=2, col = "black", alpha = 0.2)+
    geom_text(aes(x = sp_projection$year_goals_are_set+1, y = 0.01), 
              label = "<- Year goals were set", col = "black", alpha = 0.2,
              hjust=0, fontface = "bold", size = 2)+
    
    # Observed indices
    geom_errorbar(data = subset(sp_projection$indices_summarized, Year <= sp_projection$final_year_of_data),aes(x = Year, ymin = Index_q_0.025, ymax = Index_q_0.975), width = 0, col = "gray30")+
    geom_point(data = subset(sp_projection$indices_summarized, Year <= sp_projection$final_year_of_data),aes(x = Year, y = Index), col = "gray30")+
    
    # gam smooth
    geom_ribbon(data = sp_projection$gam_summary, aes(x = Year, ymin = gam_q_0.025, ymax = gam_q_0.975), alpha = 0.4, fill = "gray50")+
    geom_line(data = sp_projection$gam_summary, aes(x = Year, y = gam_med), col = "gray50", linewidth = 1)+
    
    # Historical "trend" line
    geom_ribbon(data = subset(sp_projection$StatusQuo_summary, Year <= sp_projection$year_goals_are_set), aes(x = Year, ymin = StatusQuo_q_0.025, ymax = StatusQuo_q_0.975), alpha = 0.2, fill = "black")+
    geom_line(data = subset(sp_projection$StatusQuo_summary, Year <= sp_projection$year_goals_are_set), aes(x = Year, y = StatusQuo_med), col = "black", linewidth = 1)+
    
    # Status quo projection
    geom_ribbon(data = subset(sp_projection$StatusQuo_summary, Year >= sp_projection$year_goals_are_set), aes(x = Year, ymin = StatusQuo_q_0.025, ymax = StatusQuo_q_0.975), alpha = 0.2, fill = "orangered")+
    geom_line(data = subset(sp_projection$StatusQuo_summary, Year >= sp_projection$year_goals_are_set), aes(x = Year, y = StatusQuo_med), col = "orangered", linewidth = 1)+
    
    # Recovery projection
    geom_ribbon(data = subset(sp_projection$Recovery_summary, Year >= sp_projection$year_goals_are_set), aes(x = Year, ymin = Recovery_q_0.025, ymax = Recovery_q_0.975), alpha = 0.2, fill = "dodgerblue")+
    geom_line(data = subset(sp_projection$Recovery_summary, Year >= sp_projection$year_goals_are_set), aes(x = Year, y = Recovery_med), col = "dodgerblue", linewidth = 1)+
    
    # Labels and theme
    ylab("Population Index")+
    xlab("Year")+
    theme_few()+
    ggtitle(species_name)+
    labs(subtitle = paste0("\nCurrent Status (",sp_projection$final_year_of_data,"):\n\n",
                           "Exceeds status quo trajectory: ",sp_projection$Prob_Exceed_StatusQuo,"% chance\n",
                           "Exceeds recovery trajectory: ",sp_projection$Prob_Exceed_Recovery,"% chance\n",
                           "Exceeds 1970 abundance: ",sp_projection$Prob_Exceed_Baseline,"% chance\n"))+
    
    # Plot up to 2050
    coord_cartesian(ylim=c(0,max(apply(sp_projection$indices,2,function(x) quantile(x, 0.975)))),
                    xlim=c(1970,2050))+
    scale_x_continuous(breaks = seq(1970,sp_projection$end_of_projection,10))
  print(sp_plot_index)
  
  # Save plot
  png(paste0("./summary_figures/Indices_",sp_code,".png"), units = "in", width = 8, height = 5, res = 600)
  print(sp_plot_index)
  dev.off()
  
  # ************************************************
  # Plot time series of percent population change, relative to a baseline
  # ************************************************
  
  sp_plot_percent_change <- ggplot()+
    
    # Year goals were set:
    geom_vline(xintercept = sp_projection$year_goals_are_set, size=2, col = "black", alpha = 0.2)+
    geom_text(aes(x = sp_projection$year_goals_are_set+1, y = 0), 
              label = "<- Year goals were set", col = "black", alpha = 0.2,
              hjust=0, fontface = "bold", size = 2)+
    
    # 100%
    geom_hline(yintercept = 100, linetype = 2)+
    
    # Observed indices
    geom_errorbar(data = subset(sp_projection$percent_change_summary_Obs,Year <= sp_projection$final_year_of_data), aes(x = Year, ymin = Obs_q_0.025, ymax = Obs_q_0.975), width = 0, col = "gray30")+
    geom_point(data = subset(sp_projection$percent_change_summary_Obs, Year <= sp_projection$final_year_of_data), aes(x = Year, y = Obs_med), col = "gray30")+
    
    # Historical trend
    geom_ribbon(data = subset(sp_projection$percent_change_summary_StatusQuo, Year <= sp_projection$year_goals_are_set), 
                aes(x = Year, 
                    ymin = StatusQuo_q_0.025, 
                    ymax = StatusQuo_q_0.975,
                    fill = "Historical"), alpha = 0.2)+
    geom_line(data = subset(sp_projection$percent_change_summary_StatusQuo, Year <= sp_projection$year_goals_are_set), 
              aes(x = Year, y = StatusQuo_med,
                  col = "Historical"), size = 1)+
    
    # StatusQuo projection
    geom_ribbon(data = subset(sp_projection$percent_change_summary_StatusQuo, Year >= sp_projection$year_goals_are_set), 
                aes(x = Year, 
                    ymin = StatusQuo_q_0.025, 
                    ymax = StatusQuo_q_0.975,
                    fill = "Status Quo"), 
                alpha = 0.2)+
    geom_line(data = subset(sp_projection$percent_change_summary_StatusQuo, Year >= sp_projection$year_goals_are_set), 
              aes(x = Year, y = StatusQuo_med, col = "Status Quo"), size = 1)+
    
    # Recovery projection
    geom_ribbon(data = subset(sp_projection$percent_change_summary_Recovery, Year >= sp_projection$year_goals_are_set), 
                aes(x = Year, 
                    ymin = Recovery_q_0.025, 
                    ymax = Recovery_q_0.975,
                    fill = "Recovery Target"), 
                alpha = 0.2)+
    geom_line(data = subset(sp_projection$percent_change_summary_Recovery, Year >= sp_projection$year_goals_are_set), 
              aes(x = Year, y = Recovery_med, col = "Recovery Target"), size = 1)+
    
    # Labels and theme
    scale_x_continuous(breaks = seq(1970,sp_projection$end_of_projection,10))+
    
    scale_fill_manual(values=c("gray50","dodgerblue","orangered"), name = "Scenario")+
    scale_color_manual(values=c("black","dodgerblue","orangered"), name = "Scenario")+
    
    ylab("Percent of Baseline abundance")+
    xlab("Year")+
    ggtitle(species_name)+
    coord_cartesian(ylim=c(0,200),
                    xlim=c(1970,2050))+
    labs(subtitle = paste0("\nCurrent Status (",sp_projection$final_year_of_data,"):\n\n",
                           "Exceeds status quo trajectory: ",sp_projection$Prob_Exceed_StatusQuo,"% chance\n",
                           "Exceeds recovery trajectory: ",sp_projection$Prob_Exceed_Recovery,"% chance\n",
                           "Exceeds Baseline abundance: ",sp_projection$Prob_Exceed_Baseline,"% chance\n"))+
    theme_few()
  print(sp_plot_percent_change)
  
  # Save plot
  png(paste0("./summary_figures/Percent_Change_",sp_code,".png"), units = "in", width = 8, height = 5, res = 600)
  print(sp_plot_percent_change)
  dev.off()
  
}

# Save summary table
write.csv(sp_table,file = "./results_summary/projection_table.csv",row.names = FALSE)
