library(tidyverse)
library(ggtext)

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

admits0 <- read_csv("raw-sari-admissions.csv")

date_pulled <- ymd("2024-09-04") # most recent date when the data were updated

admits <- admits0 |> 
    filter(age_group != "Unknown") |> 
    # mutate(age_group=fct_relevel(age_group, "3-19", "20-59", after=1)) |> 
    add_target_date(date_pulled) |> 
    mutate(age_group=case_match(
        age_group,
        c("<3", "3-19") ~ "Pediatric",
        .default="Adult"
    )) |> 
    group_by(date, year, week, age_group) |> 
    summarize(across(inc_flu_hosp:inc_sari_hosp, sum), .groups="drop") |> 
    filter(date < max(date)) |> # remove most recent week of data
    mutate(inc_covid_hosp=ifelse(date < "2020-01-01", NA, inc_covid_hosp)) |> 
    rename(
        epiweek=week,
        Flu=inc_flu_hosp,
        Covid=inc_covid_hosp,
        RSV=inc_rsv_hosp,
        SARI=inc_sari_hosp
    )

hc1 <- "#046CE9"
hc2 <- "#F4A247"

admits |> 
    pivot_longer(Flu:RSV) |> 
    # filter(year >= 2023) |> 
    ggplot(aes(date, value+1, col=age_group, group=age_group)) +
    geom_line() +
    facet_wrap(~name, nrow=3, scales="free_y") +
    scale_color_manual(values=c(hc1, hc2)) +
    scale_x_date(date_breaks="1 year", date_labels="%Y") +
    scale_y_continuous(trans="sqrt") +
    labs(
        x=NULL, y="new confirmed hospitalizations", col=NULL, 
        title="**individual diseases**"
    ) +
    theme_minimal() +
    theme(panel.grid.minor.x=element_blank(), legend.position="none", plot.title=element_markdown())

ggsave("mechanistic-model/figs/individual-diseases.pdf", width=5, height=5)

ggplot(admits, aes(date, SARI+1, col=age_group, group=age_group)) +
    geom_line() +
    scale_color_manual(values=c(hc1, hc2)) +
    scale_x_date(date_breaks="1 year", date_labels="%Y") +
    scale_y_continuous(trans="sqrt") +
    labs(x=NULL, y="new hospitalized SARI", col=NULL, title="**target data**") +
    theme_minimal() +
    theme(panel.grid.minor.x=element_blank(), plot.title=element_markdown())

ggsave("mechanistic-model/figs/target-data.pdf", width=5, height=3.2)

#  visualize bed capacity limits in Pediatrics------------------------------------
ped_recent <- admits |> 
    filter(age_group == "Pediatric", year %in% c(2022, 2023))

ped_max <- admits |> 
    filter(age_group == "Pediatric", !(year %in% c(2020, 2024))) |> 
    mutate(prev=SARI*8/7) |> 
    group_by(year) |> 
    slice_max(prev)

p1 <- ggplot(ped_max, aes(as.factor(year), prev)) +
    geom_col(fill="gray70", width=0.5) +
    geom_hline(yintercept=457, linewidth=1.2, col="#AF333F", linetype="dashed") +
    labs(x="year", y="bed occupancy", title="Pediatric bed occupancy for SARI") +
    theme_minimal() +
    background_grid(minor="none") +
    theme(plot.title=element_text(face="bold"), axis.title=element_text(size=rel(1.1)))

p2 <- admits |> 
    filter(age_group == "Pediatric", year < 2024) |> 
    ggplot(aes(date, RSV)) +
    geom_line(col="#4A82BA") +
    scale_x_date(date_breaks="1 year", date_labels="%Y") +
    # scale_y_continuous(trans="sqrt") +
    labs(
        x=NULL, y="RSV", col=NULL, 
        title="Pediatric laboratory-confirmed RSV"
    ) +
    theme_minimal() +
    theme(panel.grid.minor.x=element_blank(), plot.title=element_text(face="bold"), axis.title=element_text(size=rel(1.1)))

plot_grid(p1, p2, nrow=2, rel_heights=c(1, 0.75))

ggsave("pediatric-bed-at-capacity.pdf", width=5.5, height=4)

#  Explore seasonal patterns and timing of peak-----------------------------------
admits0 <- read_csv("mechanistic-model/target-sari-admissions.csv")

admits <- admits0 |>
    filter(age_group %in% c("Pediatric", "Adult", "Overall")) |>
    filter(date < max(date)) |>
    rename(epiweek=week)

peaks <- admits |> 
    filter(date < "2019-12-20" | date > "2022-03-01", date < "2023-12-20") |> 
    group_by(year, age_group) |> 
    slice_max(inc_sari_hosp) |>
    slice_max(date) |> 
    ungroup() |> 
    mutate(
        age_group=fct_relevel(age_group, "Pediatric"),
        era=ifelse(year < 2020, "2015-2019", "2022-2023")
    )

ggplot(peaks, aes(year, epiweek)) +
    geom_point() +
    # scale_color_manual()
    # scale_color_gradient2(low="blue3", mid="gray40", high="tomato", midpoint=median(peaks$inc_sari_hosp)) +
    labs(title="week peak by year") +
    theme_minimal() +
    theme(legend.position="none")

ggplot(peaks, aes(year, epiweek)) +
    geom_hline(yintercept=27, linetype="13", col="#046CE9", linewidth=1.1) +
    geom_point() +
    facet_wrap(~age_group, nrow=3) +
    theme_minimal()

ggsave("mechanistic-model/supp-figs/season-peaks.pdf", width=2.1, height=3.2)
