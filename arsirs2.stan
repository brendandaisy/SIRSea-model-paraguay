
data {
    int<lower=2> T; // number of observed time points
    int<lower=2> H; // forecast horizon - number of weeks to do prediction
    int G; // number of age groups
    array[G] int<lower=0> N; // population size
    // int K; // order of the AR process on beta
    
    // array[G, G] real M;
    
    array[T+H] int wk; // week index of each row of data
    array[T+H] int knots; // knot index of each row of data
    
    array[G, T] int y; // observed counts
    
    // real<lower=0, upper=1> alpha; // recovery rate
}

transformed data {
    // vector[G] zeros = rep_vector(0, G);
    int Tmax = T + H;
    int num_weeks = max(wk); // dimension of seasonal pattern
    int num_knots = max(knots); // dim. of "weekly" pattern
}

parameters {
    array[G] real<lower=0, upper=0.1> i_init;
    real<lower=0.7, upper=3.5> alpha; // recovery rate
    real<lower=0.1, upper=0.2> mu; // waning immunity rate
    real<lower=0, upper=0.5> rho; // probability an infectious case goes to hospital
    real<lower=0, upper=0.5> kappa; // between age mixing rate
    
    vector[num_weeks] phi; // seasonal transmission pattern
    array[G] vector[num_knots] psi; // age x week transmission pattern
    
    // shared (for now) params for each age's ts:
    real<lower=0> sd_phi;
    real<lower=0> sd_psi;
    real v0;
    real<lower=-1, upper=1> v1; // [-1, 1] indicates a stationary AR process
    real w0;
    real<lower=-1, upper=1> w1;
    // real b2;
}

transformed parameters {
    array[G] vector[Tmax] beta;
    
     matrix[G, Tmax] C; // new infections
     matrix[G, Tmax] S;
     matrix[G, Tmax] I;
     // matrix[G, T] R;
    
    // set up transmission matrix and initial conditions
    for (g in 1:G) {
        beta[g] = exp(phi[wk] + psi[g, knots]) / N[g];
        
        C[g, 1] = y[g, 1];
        S[g, 1] = N[g] * (1-i_init[g]);
        I[g, 1] = N[g]*i_init[g];
        // R[g, 1] = 0; // out of whole population, any existing, true immunity probably negligible
    }
    
    // simulate model
    for (t in 2:(Tmax)) {
        for (g in 1:G) {
            real prob_inf = 0;
            for (gg in 1:G) {
                if (g == gg)
                    prob_inf += beta[gg, t] * (1-kappa)^(G-1) * I[gg, t-1];
                else
                    prob_inf += beta[gg, t] * kappa * I[gg, t-1];
                // print(prob_inf);
            }
            C[g, t] = S[g, t-1] * prob_inf;
            S[g, t] = S[g, t-1] - C[g, t] + mu*(N[g]-S[g, t-1]-I[g, t-1]);
            I[g, t] = I[g, t-1] + C[g, t] - alpha*I[g, t-1];
            // R[g, t] = R[g, t-1] + alpha*I[g, t-1] - mu*R[g, t-1];
        }
        // print(beta[g]);
        // print(I[1:G, t]);
    }
    // print(beta);
    // print(I);
}

model {
    // hyperpriors
    rho ~ beta(1.2, 5);
    kappa ~ beta(1.2, 5);
    v0 ~ cauchy(0, 1.5);
    v1 ~ cauchy(0, 1.5);
    w0 ~ cauchy(0, 1.5);
    w1 ~ cauchy(0, 1.5);
    // b2 ~ cauchy(0, 2.5);
    sd_phi ~ gamma(1, 1);
    sd_psi ~ gamma(1, 1);
    
    phi[1] ~ normal(0, sd_phi);
    phi[2:num_weeks] ~ normal(v0 + v1*phi[1:(num_weeks-1)], sd_phi);
    
    for (g in 1:G) {
        psi[g, 1] ~ normal(0, sd_psi);
        psi[g, 2:num_knots] ~ normal(w0 + w1*psi[g, 1:(num_knots-1)], sd_psi);
        
        y[g, 2:T] ~ poisson(rho * C[g, 2:T]); // TODO: C indexing ineff.
    }
}

generated quantities {
    array[G, Tmax] int yhat;
    // array[G, Tmax] int eff_rnot;
    
    for (g in 1:G) {
        yhat[g, 1] = y[g, 1];
        yhat[g, 2:Tmax] = poisson_rng(rho * C[g, 2:Tmax]);
        
        // eff_rnot[g, 1:Tmax] = beta .* S / alpha;
    }
    // array[T+H] int yhat;
    // vector[T+H] eff_rnot;
    // array[H] real fS; // keep track of just the S compartment predictions
    // yhat[1] = y[1];
    // yhat[2:T] = poisson_rng(rho * C[2:T]);
    // eff_rnot[1:T] = beta .* S / alpha;
    // {
    //     vector[H] flogbeta;
    //     vector[H] fbeta;
    //     array[H] real fC;
    //     array[H] real fI;
    //     array[H] real fR;
    //     
    //     flogbeta[1] = normal_rng(b0 + b1*logbeta[T], sd_beta);
    //     fbeta[1] = exp(flogbeta[1]) / N;
    //     fS[1] = S[T] - fbeta[1]*S[T]*I[T] + mu*R[T];
    //     fI[1] = I[T] + fbeta[1]*S[T]*I[T] - alpha*I[T];
    //     fR[1] = R[T] + alpha*I[T] - mu*R[T];
    //     yhat[T+1] = poisson_rng(rho * fbeta[1]*S[T]*I[T]);
    //     eff_rnot[T+1] = fbeta[1] * fS[1] / alpha;
    //     
    //     for (t in 2:H) {
    //         real mu = b0;
    //         for (k in 1:K) {
    //             mu += b[k] * logbeta[t-k];
    //         }
    //         logbeta[t] ~ normal(mu, sd_beta);
    //         flogbeta[t] = normal_rng(b0 + b1*flogbeta[t-1], sd_beta);
    //         fbeta[t] = exp(flogbeta[t]) / N;
    //         
    //         fC[t] = fbeta[t]*fS[t-1]*fI[t-1];
    //         fS[t] = fS[t-1] - fC[t] + mu*fR[t-1];
    //         fI[t] = fI[t-1] + fC[t] - alpha*fI[t-1];
    //         fR[t] = fR[t-1] + alpha*fI[t-1] - mu*fR[t-1];
    //         
    //         yhat[T+t] = poisson_rng(rho * fC[t]);
    //         eff_rnot[T+t] = fbeta[t] * fS[t] / alpha;
    //     }
    // }
}
