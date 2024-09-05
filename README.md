# SIRSea

A deterministic, discrete-time SIRS model with time-varying transmission rate. 
The transmission rate is constrained in shape to allow room for other model dynamics such as
building and waning immunity, by using a B-spline penalized by an autoregressive process of order 2.

## Running the model

Use `mechanistic-forecast2.R` for the current model and nice visualizations.
