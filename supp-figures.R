library(tidyverse)
library(MASS, exclude="select")
library(splines)
library(cowplot)

ar1_corr_mat <- function(nt, w) {
    exponent <- abs(matrix(1:nt - 1, nt, nt, byrow=TRUE) - (1:nt - 1))
    w^exponent
}

sample_phi <- function(w, B, nrep=10, ar_sd=1) {
    Sigma <- ar1_corr_mat(nrow(B), 0.8)
    
    map_dfr(1:nrep, \(rep) {
        tibble(
            w=w,
            rep=rep,
            t=1:ncol(B),
            ts=c(t(accumulate(1:(nrow(B)-1), ~rnorm(1, w*.x, ar_sd), .init=0) %*% B))
            # ts=c(t(mvrnorm(1, rep(0, nrow(B)), ar_sd * Sigma)) %*% B)
        )
    })
}

spline_df <- 15
Tmax <- 150
B <- t(bs(1:Tmax, spline_df, Boundary.knots=c(-2, Tmax+3)))

phi <- bind_rows(sample_phi(0, B), sample_phi(-0.9, B), sample_phi(0.8, B))

ggplot(phi, aes(t, ts, group=rep, col=as.factor(rep))) +
    geom_line(alpha=0.7, linewidth=0.97) +
    facet_wrap(~w, nrow=1, labeller=label_both) +
    scale_color_viridis_d() +
    labs(x="time", y=NULL) +
    theme_minimal() +
    theme(legend.position="none")

ggsave("mechanistic-model/supp-figs/b-spline-ex.pdf", width=6, height=3.2)
