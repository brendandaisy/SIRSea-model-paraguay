library(tidyverse)
library(lubridate)
library(hubValidations)
library(glue)

format_quantiles <- function(pred_quants, forecast_date) {
    pred_quants <- pred_quants %>%
        filter(horizon != -1)
    #add one week to date 
    #data_check <- data 
    #data$date <- as.Date(data$date)
    # Now, add 7 days to each date in the date column
    
    data1 <- pred_quants |>
        rename(target_end_date = date,
               output_type_id = quant_level) |>
        mutate(value = round(value)) |>
        filter(target_end_date >= forecast_date) |>
        mutate(reference_date = forecast_date,
               output_type = "quantile",
               target = "inc sari hosp",
               age_group = str_to_title(age_group),
               horizon = as.numeric(as.factor(target_end_date)) - 1,  
               #horizon = as.numeric(date - min(date, na.rm = TRUE)),
               output_type_id = format(as.numeric(sub("%", "", output_type_id)) / 100, trim = TRUE)) |>
        select(reference_date, target, horizon, target_end_date, age_group, output_type, output_type_id, value)
    data1 %>%
        mutate(output_type_id = as.numeric(output_type_id),  # Convert to numeric first, if not already
               output_type_id = gsub("0+$", "",              # Remove trailing zeros
                                     gsub("\\.$", "",        # Remove dot if it's the last character
                                          sprintf("%.3f", output_type_id))))
}

check_format <- function(submission_file) {
    sub_validation <- hubValidations::validate_submission(
        hub_path="/Users/brendaisy/Documents/Projects/Paraguay Forecasting/NCIRD-GIB-Paraguay-Forecasting",
        file_path=submission_file
    )
    ## Want all green checkmarks
    print(sub_validation)
    ## Want to make sure there are no missing required values
    print(glue("Number of missing values: {nrow(sub_validation$req_vals$missing)}"))
}
