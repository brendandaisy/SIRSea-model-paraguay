
data {
    int<lower=0> N; // population size
    int<lower=2> T; // number of time points (=sample size currently)
    int<lower=2> H; // forecast horizon - number of weeks to do prediction
    int K; // order of the AR process on beta
    
    array[T] int wk; // week of each row of data
    
    array[T] int y; // observed counts
}

transformed data {
    int num_weeks = max(wk);
}

parameters {
    real<lower=y[1]*1.0/N, upper=0.5> i_init;
    real<lower=0, upper=1> alpha;
    vector[num_weeks] phi;
    vector[T]
    real<lower=0, upper=1> mu;
    real<lower=0, upper=1> rho;
    
    real<lower=0> sigma_trans;
    real b0;
    array[K] real<lower=-1, upper=1> b; // uniform on [-1, 1] gives a stationary AR process
}

transformed parameters {
    vector[T] beta;
    beta = exp(logbeta) / N;
    
    vector[T] C; // new infections
    vector[T] S;
    array[T] real I;
    array[T] real R;
    C[1] = y[1];
    S[1] = N * (1-i_init);
    I[1] = N*i_init;
    R[1] = 0; // out of whole population, any existing, true immunity probably negligible
    
    // simulate model up until T
    for (t in 2:T) {
        C[t] = beta[t]*S[t-1]*I[t-1];
        S[t] = S[t-1] - C[t] + mu*R[t-1];
        I[t] = I[t-1] + C[t] - alpha*I[t-1];
        R[t] = R[t-1] + alpha*I[t-1] - mu*R[t-1];
    }
}

model {
    // hyperpriors
    rho ~ beta(1.2, 5);
    b0 ~ cauchy(0, 2.5);
    // log((1+b1) / (1-b1)) ~ cauchy(0, 2.5);
    sigma_trans ~ cauchy(0, 2.5); // TODO: make super informative
    
    logbeta[1] ~ normal(0, sigma_trans);
    // update the prior constraints for beta
    for (t in 2:T) {
        real mu = b0;
        for (k in 1:(min(t, K+1)-1)) {
            mu += b[k] * logbeta[t-k];
        }
        logbeta[t] ~ normal(mu, sigma_trans);
        
        // update likelihood contribution
        y[t] ~ poisson(rho * C[t]);
    }
    // for (t in 2:T) {
    //     
    //     logbeta[t] ~ normal(b0 + b1*logbeta[t-1], sigma_trans);
    //     
    //     
    // }
}

generated quantities {
    array[T+H] int yhat;
    vector[T+H] eff_rnot;
    array[H] real fS; // keep track of just the S compartment predictions
    yhat[1] = y[1];
    yhat[2:T] = poisson_rng(rho * C[2:T]);
    eff_rnot[1:T] = beta .* S / alpha;
    {
        vector[H] flogbeta;
        vector[H] fbeta;
        array[H] real fC;
        array[H] real fI;
        array[H] real fR;
        
        flogbeta[1] = normal_rng(b0 + b1*logbeta[T], sigma_trans);
        fbeta[1] = exp(flogbeta[1]) / N;
        fS[1] = S[T] - fbeta[1]*S[T]*I[T] + mu*R[T];
        fI[1] = I[T] + fbeta[1]*S[T]*I[T] - alpha*I[T];
        fR[1] = R[T] + alpha*I[T] - mu*R[T];
        yhat[T+1] = poisson_rng(rho * fbeta[1]*S[T]*I[T]);
        eff_rnot[T+1] = fbeta[1] * fS[1] / alpha;
        
        for (t in 2:H) {
            real mu = b0;
            for (k in 1:K) {
                mu += b[k] * logbeta[t-k];
            }
            logbeta[t] ~ normal(mu, sigma_trans);
            flogbeta[t] = normal_rng(b0 + b1*flogbeta[t-1], sigma_trans);
            fbeta[t] = exp(flogbeta[t]) / N;
            
            fC[t] = fbeta[t]*fS[t-1]*fI[t-1];
            fS[t] = fS[t-1] - fC[t] + mu*fR[t-1];
            fI[t] = fI[t-1] + fC[t] - alpha*fI[t-1];
            fR[t] = fR[t-1] + alpha*fI[t-1] - mu*fR[t-1];
            
            yhat[T+t] = poisson_rng(rho * fC[t]);
            eff_rnot[T+t] = fbeta[t] * fS[t] / alpha;
        }
    }
}
