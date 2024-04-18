library(cmdstanr)
library(posterior)
library(tidyverse)
library(bayesplot)
library(cowplot)

# add_target_date <- function(adm) {
#     # isoweek 1 for 2015 began the 4th, so remove one day to get to first Saturday for which we have data
#     date_start <- ymd("2015-01-03")
#     num_dates <- distinct(adm, year, week) |> nrow() # need to remove an extra date, assuming weeks we get are epiweeks   s
#     num_age_groups <- length(unique(adm$age_group))
#     
#     mutate(
#         adm, 
#         date=rep(date_start + weeks(1:num_dates - 1), each=num_age_groups),
#         .before=year
#     )
# }

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

# Choose subset of data to start working with-------------------------------------

admits_sub <- admits |> 
    filter(date > "2021-09-01") |>
    mutate(age_group=case_when(
        age_group %in% c("20-59", "60+") ~ "adult",
        TRUE ~ "pediatric"
    )) |> 
    group_by(date, year, epiweek, age_group) |>
    summarize(inc_sari_hosp=sum(inc_sari_hosp), .groups="drop") |> 
    arrange(age_group, date)

ggplot(admits_sub, aes(date, inc_sari_hosp)) +
    geom_line() +
    facet_wrap(~age_group, nrow=1) +
    scale_x_date(date_breaks="3 month", date_labels="%m-%Y", guide=guide_axis(angle=45)) +
    # scale_y_continuous(trans="sqrt") +
    labs(x=NULL, y="new hospitalizations", col=NULL) +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank())

# Fit a stan model----------------------------------------------------------------
ord_knots <- 3
weeks_ahead <- 5
Tmax <- length(unique(admits_sub$date)) + weeks_ahead
wk <- c(filter(admits_sub, age_group == "adult")$epiweek, 15:19) # TODO: fix this
knots <- rep(1:ceiling(Tmax/ord_knots), each=ord_knots)[1:Tmax]

ymat <- matrix(admits_sub$inc_sari_hosp, nrow=2, byrow=TRUE) # remember adult is first row

stan_dat <- list(
    `T`=length(unique(admits_sub$date)),
    H=weeks_ahead,
    G=2,
    N=c(1e6, 1e6), # let's just say a million people in each group at risk
    wk=wk,
    knots=knots,
    y=ymat
    # alpha=0.7 # inf. period of 10 days
)

init_params <- function() {
    list(
        i_init=runif(2, 1e-5, 1e-3),
        rho=runif(1, 0, 0.1),
        kappa=runif(1, 0, 0.1),
        alpha=runif(1, 0.7, 1),
        phi=rnorm(53, -4, 0.5),
        psi=matrix(rnorm(2*47, -4, 0.5), nrow=2)
        # v0=runif(1, -9, -7),
        # w0=runif(1, -9, -7),
        # v1=0,
        # w1=0
        # v1=runif(1, -0.5, 0.5),
        # w1=runif(1, -0.5, 0.5)
    )
}

exec <- cmdstan_model("mechanistic-model/arsirs2.stan")

fit <- exec$sample(
    data=stan_dat,
    chains=8,
    parallel_chains=8,
    init=init_params,
    iter_sampling=500,
    iter_warmup=1000,
    adapt_delta=0.95,
    max_treedepth=15
)

(fit_summ <- fit$summary())

post <- as_draws_df(fit$draws())

mcmc_intervals(fit$draws(variables=c("alpha", "mu", "rho", "kappa", "sd_phi", "sd_psi", "v1", "w1")))

last_date <- max(admits_sub$date)

# Save the predicted quantiles----------------------------------------------------
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

post_pred <- post |> 
    mutate(draw=.draw) |> 
    select(contains("yhat"), draw) |> 
    pivot_longer(-draw, values_to="pred_count") |> 
    mutate(
        age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "adult", "pediatric"), 
        t=as.integer(str_extract(name, "(?<=,)\\d+")),
        .keep="unused", .before=pred_count
    ) |> 
    arrange(draw, age_group, t) |> 
    group_by(draw, age_group) |> 
    mutate(date=c(unique(admits_sub$date), last_date + weeks(1:weeks_ahead))) |> 
    ungroup()

pp_overall <- post_pred |> 
    group_by(draw, date) |> 
    summarise(pred_count=sum(pred_count), .groups="drop") |> 
    mutate(age_group="overall")

post_pred_quants <- bind_rows(post_pred, pp_overall) |> 
    filter(date > last_date) |> 
    group_by(age_group, date) |> 
    reframe(enframe(quantile(pred_count, quantiles_needed), "quant_level")) |> 
    arrange(date, age_group)

write_csv(post_pred_quants, paste0("mechanistic-model/submission-files/quantiles-", today(), ".csv"))

# Prepare the weekly model summary------------------------------------------------
post_pred_summ <- post_pred |> 
    group_by(age_group, date) |> 
    reframe(enframe(quantile(pred_count, c(0.025, 0.25, 0.5, 0.75, 0.975)))) |> 
    pivot_wider() |> 
    left_join(admits_sub, by=c("age_group", "date"))

p1 <- ggplot(post_pred_summ, aes(date)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="gray80", col=NA, alpha=0.6) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="gray60", col=NA, alpha=0.6) +
    geom_point(aes(y=inc_sari_hosp), col="tomato", size=0.85) +
    facet_wrap(~age_group) +
    scale_x_date(date_breaks="3 month", date_labels="%m-%Y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="SARI (all data)", title="model fit") +
    theme_bw()

prev_counts_date <- last_date - years(1) - weeks(10)

prev_counts <- admits_sub |> 
    filter(date > prev_counts_date, date <= (last_date - years(1) + weeks(weeks_ahead))) |> 
    mutate(date=rep(last_date + weeks(-9:weeks_ahead), times=2))

p2 <- post_pred_summ |> 
    filter(date > (last_date - weeks(10))) |> 
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="gray80", col=NA, alpha=0.6) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="gray60", col=NA, alpha=0.6) +
    geom_point(aes(y=inc_sari_hosp), prev_counts, col="blue3", size=0.92) +
    geom_point(aes(y=inc_sari_hosp), col="tomato", size=0.95) +
    facet_wrap(~age_group) +
    scale_x_date(date_breaks="1 week", date_labels="%d-%b", guide=guide_axis(angle=45)) +
    labs(x=NULL, y="SARI (latest and previous year)", title="predictions") +
    theme_bw() +
    theme(panel.grid.minor.x=element_blank())

plot_grid(p1, p2, nrow=2)

ggsave(paste0("mechanistic-model/submission-files/predictions", today(), ".pdf"), width=7, height=7)

post_phi_summ <- post |> 
    select(contains("phi[")) |> 
    pivot_longer(everything(), values_to="phi") |> 
    mutate(t=parse_number(name), .keep="unused") |> 
    group_by(t) |> 
    reframe(enframe(quantile(phi, c(0.025, 0.25, 0.5, 0.75, 0.975)))) |> 
    pivot_wider()

p1 <- ggplot(post_phi_summ, aes(t)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="lightblue", col=NA, alpha=0.5) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="blue3", col=NA, alpha=0.4) +
    labs(x="week", title="seasonal transmission effect") +
    theme_bw()

post_psi_summ <- post |> 
    select(contains("psi[")) |> 
    pivot_longer(everything(), values_to="psi") |> 
    mutate(
        age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "adult", "pediatric"), 
        t=as.integer(str_extract(name, "(?<=,)\\d+"))
    ) |> 
    group_by(age_group, t) |> 
    reframe(enframe(quantile(psi, c(0.025, 0.25, 0.5, 0.75, 0.975)))) |> 
    pivot_wider()

p2 <- ggplot(post_psi_summ, aes(t)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="lightblue", col=NA, alpha=0.5) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="blue3", col=NA, alpha=0.4) +
    facet_wrap(~age_group, nrow=1) +
    labs(x="knot index", title="weekly transmission effect") +
    theme_bw()

# post_rnot_summ <- post |> 
#     select(contains("eff")) |> 
#     pivot_longer(everything(), values_to="eff_rnot") |> 
#     mutate(t=parse_number(name), .keep="unused") |> 
#     filter(t > 1) |> # t=1 isn't really constrained by data, so not useful to compare
#     group_by(t) |> 
#     reframe(enframe(quantile(eff_rnot, c(0.025, 0.25, 0.5, 0.75, 0.975)))) |> 
#     pivot_wider()
# 
# p2 <- ggplot(post_rnot_summ, aes(t)) +
#     geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="lightblue", col=NA, alpha=0.5) +
#     geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="blue3", col=NA, alpha=0.4) +
#     geom_hline(yintercept=1, linetype="dashed", col="tomato") +
#     labs(title="effective reproductive number") +
#     theme_bw()

library(cowplot)


###
post_sus_summ <- post |> 
    select(contains("S[")) |> 
    pivot_longer(everything(), values_to="S") |> 
    mutate(
        age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "adult", "pediatric"), 
        t=as.integer(str_extract(name, "(?<=,)\\d+"))
    ) |> 
    group_by(age_group, t) |> 
    reframe(enframe(quantile(S, c(0.025, 0.25, 0.5, 0.75, 0.975)))) |> 
    pivot_wider()

p3 <- ggplot(post_sus_summ, aes(t)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="lightblue", col=NA, alpha=0.5) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="blue3", col=NA, alpha=0.4) +
    geom_hline(yintercept=stan_dat$N, linetype="dashed", col="tomato") +
    # geom_vline(xintercept=stan_dat$`T`, linetype="dashed", col="blue4") +
    facet_wrap(~age_group, nrow=1) +
    labs(x="time", title="susceptible population") +
    theme_bw()

plot_grid(p1, p2, p3, nrow=3)

ggsave(paste0("mechanistic-model/submission-files/latent-effects", today(), ".pdf"), width=6, height=8)






