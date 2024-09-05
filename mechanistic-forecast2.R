library(cmdstanr)
library(posterior)
library(splines)
library(tidyverse)
library(bayesplot)
library(cowplot)

source("mechanistic-model/format-quantiles.R")

post_array_quantiles <- function(post, var, age_dim=FALSE, q=c(0.025, 0.25, 0.75, 0.975)) {
    post_long <- post |> 
        select(contains(paste0(var, "["))) |> 
        pivot_longer(everything(), values_to=var)
    
    if (age_dim) {
        ret <- post_long |> 
            mutate(
                age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "Pediatric", "Adult"), 
                age_group=fct_relevel(age_group, "Pediatric"),
                t=as.integer(str_extract(name, "(?<=,)\\d+"))
            ) |> 
            group_by(age_group, t)
    }
    else
        ret <- mutate(post_long, t=parse_number(name), .keep="unused") |> group_by(t)
    
    ret |> 
        summarize(
            mean=mean(.data[[var]]),
            qs = list(value = quantile(.data[[var]], probs=q)), 
            .groups="drop"
        ) |> 
        unnest_wider(qs)
}

admits0 <- read_csv("target-sari-admissions.csv")

# TODO 6/25: it makes sense to incorporate some age-specific epi params. Sure it doesn't matter much for
# data fitting, but if we trust the underlying and the effect they might have on forecasts, they can 
# indeed impact forecast trajectory in hopefully an accurate way

# TODO 7/18: it may be useful to try and scale the individual disease data to population counts, since
# it may be awkward to work with population N while there are like 50 or less of each disease
# We could try to scale the data based on estimated coverage from the sentinel sites

DROP_2_WEEKS <- FALSE

admits <- admits0 |>
    filter(age_group != "Unknown") |>
    filter(date < (if (DROP_2_WEEKS) max(date) - weeks(1) else max(date))) |> # remove most recent 2 weeks of data
    rename(epiweek=week)

# Choose a subset of more recent data for the mechanistic model-------------------
admits_sub <- admits |>
    filter(date > "2021-09-01", age_group %in% c("Pediatric", "Adult")) |>
    mutate(age_group=fct_relevel(age_group, "Pediatric")) |> 
    filter(epiweek != 53) |> # for now, since only one instance of it in the training data
    arrange(age_group, date)

# Prepare data for stan-----------------------------------------------------------
weeks_ahead <- if (DROP_2_WEEKS) 6 else 5 # forecast ahead weeks 0 through 4
Tmax <- length(unique(admits_sub$date)) + weeks_ahead
spline_df <- round(Tmax/7) # the number of knots to place evenly (approx every 7 weeks)
B <- t(bs(1:Tmax, spline_df, Boundary.knots=c(-2, Tmax+3)))

# convert observed SARI to a matrix:
ymat <- matrix(admits_sub$inc_sari_hosp, nrow=2, byrow=TRUE) # remember pediatric is first row

stan_dat <- list(
    `T`=length(unique(admits_sub$date)),
    H=weeks_ahead,
    G=2,
    N=c(1e6, 1e6), # let's just say a million people in each group at risk
    df=spline_df,
    B=B,
    y=ymat,
    # i_init=c(0.005, 0.005),
    alpha=c(0.9, 0.7), # pediatrics recover more quickly
    sd_phi_df=3,
    sd_phi_scale=5
)

# Fit the stan model--------------------------------------------------------------

# compile the current model
exec <- cmdstan_model("mechanistic-model/sirsea2.stan")

# fit model with MCMC
fit <- exec$sample(
    data=stan_dat,
    chains=8,
    parallel_chains=8,
    # init=init_params,
    iter_sampling=1000,
    iter_warmup=10000,
    adapt_delta=0.9
    # max_treedepth=12
)

(fit_summ <- fit$summary()) # want rhat all <1.01
post <- as_draws_df(fit$draws())

# MCMC model diagnostics----------------------------------------------------------
fit_nuts <- nuts_params(fit)
fit_lp <- log_posterior(fit)

mcmc_trace(post, pars=c("lp__"), size=1.2, np=fit_nuts) +
    scale_color_viridis_d()

# Posterior epi parameters-------------------------------------------------------
hc1 <- "#046CE9"
hc2 <- "#F4A247"

epi_params <- post |> 
    select(c("mu", "kappa", "w1", "w2", contains("rho"))) |> 
    pivot_longer(everything())

guide_lines <- tibble(name=unique(epi_params$name), vline=c(NA, NA, 0, 0, NA, NA))

ggplot(epi_params, aes(x=value)) +
    geom_density(fill="lightblue", col="gray40", alpha=0.7) +
    geom_vline(aes(xintercept=vline), guide_lines, col=hc2, linetype="13", linewidth=1.05) +
    facet_wrap(~name, ncol=2, scales="free") +
    theme_minimal() +
    background_grid(minor="none")

ggsave("mechanistic-model/supp-figs/post-epi-params.pdf", width=4.1, height=4.5)

# Prepare the weekly model summary------------------------------------------------
last_date <- max(admits_sub$date)
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

post_pred <- post |> 
    mutate(draw=.draw) |> 
    select(contains("yhat"), draw) |> 
    pivot_longer(-draw, values_to="pred_count") |> 
    mutate(
        age_group=ifelse(str_extract(name, "\\d+(?=,)") == "1", "Pediatric", "Adult"), 
        age_group=fct_relevel(age_group, "Pediatric"),
        t=as.integer(str_extract(name, "(?<=,)\\d+")),
        .keep="unused", .before=pred_count
    ) |> 
    arrange(draw, age_group, t) |> 
    group_by(draw, age_group) |> 
    mutate(date=c(unique(admits_sub$date), last_date + weeks(1:weeks_ahead))) |> 
    ungroup()

post_pred_summ <- post_pred |> 
    filter(t > 1) |> 
    group_by(age_group, date) |> 
    summarize(
        mean=mean(pred_count),
        qs=list(value=quantile(pred_count, c(0.025, 0.25, 0.5, 0.75, 0.975))), 
        .groups="drop"
    ) |> 
    unnest_wider(qs) |> 
    left_join(admits_sub, by=c("age_group", "date"))

library(ggtext)

theme_mod_output <- function() {
    theme_minimal() +
        theme(
            panel.grid.minor=element_blank(),
            plot.title=element_markdown(size=13),
            axis.title.y=element_markdown(),
            strip.text=element_text(size=12)
        )
}

p1 <- ggplot(post_pred_summ, aes(date)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="gray80", col=NA, alpha=0.6) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="gray60", col=NA, alpha=0.6) +
    geom_point(aes(y=inc_sari_hosp), col=hc1, size=0.83) +
    facet_wrap(~age_group) +
    coord_cartesian(ylim=c(NA, quantile(admits_sub$inc_sari_hosp, 0.99)+5)) +
    scale_x_date(date_breaks="4 month", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(
        x=NULL, 
        y=glue::glue("SARI (all <span style='color:{hc1}'>training</span> data)"), 
        title=paste0("***SIRSea* model fit and predictions** - ", format(today(), "%d %b %Y"))
    ) +
    theme_mod_output()

prev_counts_date <- last_date - years(1) - weeks(10)

prev_counts <- admits_sub |> 
    filter(date > prev_counts_date, date <= (last_date - years(1) + weeks(weeks_ahead))) |> 
    mutate(date=rep(last_date + weeks(-9:weeks_ahead), times=2))

# get the exact prediction dates so that the breaks line up
recent_dates <- post_pred_summ |> 
    filter(date > (last_date - weeks(10))) |> 
    pull(date) |> 
    unique()

p2 <- post_pred_summ |> 
    filter(date > (last_date - weeks(10))) |> 
    ggplot(aes(date)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill="gray80", col=NA, alpha=0.6) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill="gray60", col=NA, alpha=0.6) +
    geom_line(aes(y=mean), col="gray40", alpha=0.8) +
    geom_point(aes(y=inc_sari_hosp), prev_counts, col=hc2, size=0.95) +
    geom_point(aes(y=inc_sari_hosp), col=hc1, size=0.98) +
    facet_wrap(~age_group) +
    scale_x_date(breaks=recent_dates, date_labels="%d %b", guide=guide_axis(angle=45)) +
    labs(
        x=NULL, 
        y=glue::glue("SARI (<span style='color:{hc1}'>latest</span> and <span style='color:{hc2}'>previous year</span>)")
    ) +
    theme_mod_output()

plot_grid(p1, p2, nrow=2)

ggsave(paste0("mechanistic-model/submission-files/predictions-", today(), ".pdf"), width=6.8, height=5.6)

# Save the predicted quantiles----------------------------------------------------
forecast_date <- "2024-09-07" # the "official" date for this week. Saturday following submission

pp_overall <- post_pred |> 
    group_by(draw, date) |> 
    summarise(pred_count=sum(pred_count), .groups="drop") |> 
    mutate(age_group="Overall")

post_pred_quants <- bind_rows(post_pred, pp_overall) |> 
    filter(date > last_date) |> 
    mutate(horizon=as.numeric(as.factor(date)) - 1) |> # TODO: assumes NOT dropping 2 weeks atm 
    group_by(age_group, date, horizon) |> 
    reframe(enframe(quantile(pred_count, quantiles_needed), "quant_level")) |> 
    arrange(date, age_group) |> 
    format_quantiles(forecast_date) # following Saturday

fname <- paste0("UGA_flucast-SIRSea/", forecast_date, "-UGA_flucast-SIRSea.csv")
write_csv(post_pred_quants, paste0("/Users/brendaisy/Documents/Projects/Paraguay Forecasting/NCIRD-GIB-Paraguay-Forecasting/model-output/", fname))
check_format(fname)

### page 2
post_beta_summ <- post_array_quantiles(post, "beta", age_dim=TRUE) |> 
    mutate(across(mean:`97.5%`, ~.x*1e6)) |> 
    left_join(distinct(post_pred, t, age_group, date))

p1 <- ggplot(post_beta_summ, aes(date)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), col=NA, fill="gray80", alpha=0.6) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), col=NA, fill="gray60", alpha=0.6) +
    geom_line(aes(y=mean), col="gray40", alpha=0.8) +
    geom_vline(xintercept=last_date, col=hc1, linetype="13", linewidth=1.05) +
    facet_wrap(~age_group, nrow=1) +
    scale_y_continuous(labels=scales::scientific) +
    scale_x_date(date_breaks="4 month", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y=NULL, title="**transmission rate**") +
    theme_mod_output()

post_sus_summ <- post_array_quantiles(post, "S", age_dim=TRUE) |> 
    left_join(distinct(post_pred, t, age_group, date))

p2 <- ggplot(post_sus_summ, aes(date)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), col=NA, fill="gray80", alpha=0.6) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), col=NA, fill="gray60", alpha=0.6) +
    geom_line(aes(y=mean), col="gray40", alpha=0.8) +
    geom_vline(xintercept=last_date, col=hc1, linetype="13", linewidth=1.05) +
    facet_wrap(~age_group, nrow=1) +
    scale_x_date(date_breaks="4 month", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y=NULL, title="**susceptible population**") +
    theme_mod_output()

post_rnot_summ <- post_array_quantiles(post, "eff_rnot") |> 
    left_join(distinct(post_pred, t, date))

p3 <- ggplot(post_rnot_summ, aes(date)) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), col=NA, fill="gray80", alpha=0.6) +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), col=NA, fill="gray60", alpha=0.6) +
    geom_line(aes(y=mean), col="gray40", alpha=0.8) +
    geom_hline(yintercept=1, col=hc2, linetype="13", linewidth=1.05) +
    geom_vline(xintercept=last_date, col=hc1, linetype="13", linewidth=1.05) +
    # facet_wrap(~age_group, nrow=1) +
    scale_x_date(date_breaks="4 month", date_labels="%b '%y", guide=guide_axis(angle=45)) +
    labs(x=NULL, y=NULL, title="**effective reproductive number**") +
    theme_mod_output()

pbottom <- plot_grid(NULL, p3, NULL, nrow=1, rel_widths=c(1, 2, 1))
plot_grid(plot_grid(p1, p2, nrow=2), pbottom, nrow=2, rel_heights=c(2, 1))
# pright <- plot_grid(NULL, p3, NULL, ncol=1, rel_heights=c(1, 2, 1))
# plot_grid(plot_grid(p1, p2, nrow=2), pright, nrow=1, rel_widths=c(1.7, 1))

ggsave(paste0("mechanistic-model/submission-files/latent-effects", today(), ".pdf"), width=5.5, height=6)
