library(tidyverse)

# estimate the date of each Saturday for the provide year and week labels
add_target_date <- function(adm, date_pulled) {
    date_end <- date_pulled - days(wday(date_pulled))
    num_dates <- distinct(adm, year, week) |> nrow()
    num_age_groups <- length(unique(adm$age_group))
    
    adm |> 
        arrange(year, week) |> # make extra sure things are in order
        mutate(
            date=rev(rep(date_end - weeks(1:num_dates - 1), each=num_age_groups)),
            .before=year
        )
}

admits0 <- read_csv("mechanistic-model/raw-sari-admissions.csv")

date_pulled <- ymd("2024-04-16") # most recent date when the data were updated

admits <- admits0 |> 
    filter(age_group != "Unknown") |> 
    mutate(age_group=fct_relevel(age_group, "3-19", "20-59", after=1)) |> 
    add_target_date(date_pulled) |> 
    filter(date < max(date)) |> # remove most recent week of data
    rename(epiweek=week)

admits |> 
    pivot_longer(inc_flu_hosp:inc_rsv_hosp) |> 
    ggplot(aes(date, value+1, col=age_group, group=age_group)) +
    geom_line() +
    facet_wrap(~name, nrow=3, scales="free_y") +
    scale_x_date(date_breaks="3 month", date_labels="%m-%Y", guide=guide_axis(angle=45)) +
    scale_y_continuous(trans="sqrt") +
    labs(x=NULL, y="new hospitalizations", col=NULL) +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank())