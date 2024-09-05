library(cmdstanr)
library(posterior)
library(tidyverse)
library(bayesplot)
library(cowplot)

# estimate the date of each Saturday for the provide year and week labels
add_target_date <- function(adm, date_pulled) {
    date_end <- date_pulled + days(7 - wday(date_pulled))
    num_dates <- distinct(adm, year, week) |> nrow()
    num_age_groups <- length(unique(adm$age_group))
    
    adm |> 
        arrange(year, week) |> # make extra sure things are in order
        mutate(
            date=rev(rep(date_end - weeks(1:num_dates - 1), each=num_age_groups)),
            .before=year
        )
}

post_array_quantiles <- function(post, var, age_dim=FALSE, q=c(0.025, 0.25, 0.75, 0.975)) {
    post_long <- post |> 
        select(contains(paste0(var, "["))) |> 
        pivot_longer(everything(), values_to=var)
    
    if (age_dim) {
        ret <- post_long |> 
            mutate(
                age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "adult", "pediatric"), 
                t=as.integer(str_extract(name, "(?<=,)\\d+"))
            ) |> 
            group_by(age_group, t)
    }
    else
        ret <- mutate(post_long, t=parse_number(name), .keep="unused") |> group_by(t)
    
    ret |> 
        reframe(enframe(quantile(.data[[var]], q))) |> 
        pivot_wider()
}

admits0 <- read_csv("mechanistic-model/target-sari-admissions.csv")
# admits0 <- read_csv("mechanistic-model/raw-sari-admissions.csv")

admits <- admits0 |>
    filter(age_group != "Unknown") |>
    filter(date < max(date)) |> # remove most recent week of data
    rename(epiweek=week)
# admits <- admits0 |> 
#     filter(age_group != "Unknown") |> 
#     mutate(age_group=fct_relevel(age_group, "3-19", "20-59", after=1)) |> 
#     add_target_date(date_pulled) |> 
#     filter(date < max(date)) |> # remove most recent week of data
#     rename(epiweek=week)

# Choose a subset of more recent data for the mechanistic model-------------------
admits_sub <- admits |>
    filter(date > "2021-09-01", age_group %in% c("Pediatric", "Adult")) |>
    filter(epiweek != 53) |> # for now, since only one instance of it in the training data
    arrange(age_group, date)
# admits_sub <- admits |> 
#     filter(date > "2021-09-01") |>
#     mutate(age_group=case_when(
#         age_group %in% c("20-59", "60+") ~ "adult",
#         TRUE ~ "pediatric"
#     )) |> 
#     group_by(date, year, epiweek, age_group) |>
#     summarize(inc_sari_hosp=sum(inc_sari_hosp), .groups="drop") |> 
#     arrange(age_group, date)

# Prepare data for stan-----------------------------------------------------------
ord_knots <- 3 # the "order" of the non-repeating transmission effect
weeks_ahead <- 5 # forecast ahead weeks 0 through 4
Tmax <- length(unique(admits_sub$date)) + weeks_ahead
wk_obs <- filter(admits_sub, age_group == "Adult")$epiweek
wk_pred <- last(wk_obs)+1:weeks_ahead
wk <- c(wk_obs, ifelse(wk_pred > 52, wk_pred %% 53 + 1, wk_pred))
knots <- rep(1:ceiling(Tmax/ord_knots), each=ord_knots)[1:Tmax] # indexes for non-repeating trans. effect

# convert observed SARI to a matrix:
ymat <- matrix(admits_sub$inc_sari_hosp, nrow=2, byrow=TRUE) # remember adult is first row

stan_dat <- list(
    `T`=length(unique(admits_sub$date)),
    H=weeks_ahead,
    G=2,
    N=c(1e6, 1e6), # let's just say a million people in each group at risk
    wk=wk,
    knots=knots,
    y=ymat,
    i_init=c(0.005, 0.005)
)

# Fit the stan model--------------------------------------------------------------
init_params <- function() {
    list(
        # i_init=runif(2, 1e-5, 1e-3),
        rho=runif(1, 0, 0.1),
        kappa=runif(1, 0, 0.1),
        alpha=runif(1, 0.7, 1),
        phi=rnorm(52, -4, 0.5),
        psi=matrix(rnorm(2*50, -4, 0.5), nrow=2)
    )
}

# compile the current model
exec <- cmdstan_model("mechanistic-model/arsirs2.stan")

# fit model with MCMC
fit <- exec$sample(
    data=stan_dat,
    chains=8,
    parallel_chains=8,
    init=init_params,
    iter_sampling=1250,
    iter_warmup=9000,
    adapt_delta=0.955,
    max_treedepth=12
)

(fit_summ <- fit$summary())
post <- as_draws_df(fit$draws())

# MCMC model diagnostics----------------------------------------------------------
fit_nuts <- nuts_params(fit)
fit_lp <- log_posterior(fit)

mcmc_trace(post, pars=c("lp__"), size=1.2, np=fit_nuts) +
    scale_color_viridis_d()

# mechanistic scalar parameters. none of these look problematic or very correlated, just
# not very well constrained by the data
mcmc_pairs(
    post, pars=c("alpha", "mu", "rho", "kappa"), np=fit_nuts,
    off_diag_args=list(size=0.75)
)

# transmission rate parameters appear most correlated with TODO??? whittle epi params down
mcmc_pairs(
    post, pars=c("alpha", "rho", "sd_phi", "sd_psi", "v1", "w1"), np=fit_nuts,
    off_diag_args=list(size=0.75)
)

mcmc_nuts_divergence(fit_nuts, fit_lp)

# Save the predicted quantiles----------------------------------------------------
last_date <- max(admits_sub$date)
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

post_pred <- post |> 
    mutate(draw=.draw) |> 
    select(contains("yhat"), draw) |> 
    pivot_longer(-draw, values_to="pred_count") |> 
    mutate(
        age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "Adult", "Pediatric"), 
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
    mutate(age_group="Overall")

post_pred_quants <- bind_rows(post_pred, pp_overall) |> 
    filter(date > last_date) |> 
    mutate(horizon=as.numeric(as.factor(date)) - 1) |> 
    group_by(age_group, date, horizon) |> 
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

ggsave(paste0("mechanistic-model/submission-files/predictions", today(), ".pdf"), width=7, height=6)

### page 2

post_phi_summ <- post_array_quantiles(post, "phi", age_dim=FALSE)

p1 <- ggplot(post_phi_summ, aes(t)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="lightblue", col=NA, alpha=0.5) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="blue3", col=NA, alpha=0.4) +
    labs(x="week", title="seasonal transmission effect") +
    theme_bw()

post_psi_summ <- post_array_quantiles(post, "psi", age_dim=TRUE)

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

post_sus_summ <- post_array_quantiles(post, "S", age_dim=TRUE)

p3 <- ggplot(post_sus_summ, aes(t)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="lightblue", col=NA, alpha=0.5) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="blue3", col=NA, alpha=0.4) +
    geom_hline(yintercept=stan_dat$N, linetype="dashed", col="tomato") +
    # geom_vline(xintercept=stan_dat$`T`, linetype="dashed", col="blue4") +
    facet_wrap(~age_group, nrow=1) +
    labs(x="time", title="susceptible population") +
    theme_bw()

plot_grid(plot_grid(NULL, p1, NULL, rel_widths=c(1, 3, 1), nrow=1), p2, p3, nrow=3)

ggsave(paste0("mechanistic-model/submission-files/latent-effects", today(), ".pdf"), width=6, height=8)


###Another one
tmp <- post |> 
    subset_draws(variable=c("C", "I")) |> 
    mutate(draw=.draw) |> 
    pivot_longer(-draw) |> 
    mutate(
        age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "Adult", "Pediatric"), 
        t=as.integer(str_extract(name, "(?<=,)\\d+"))
    ) |> 
    filter(t > 1) |> 
    mutate(
        name=str_extract(name, "\\w"),
        t=ifelse(name == "I", t+1, t)
    )

tmp |> 
    pivot_wider() |> 
    mutate(deltaI=C-)
    





