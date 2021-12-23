## causamed
Implementation of Tyler Vanderweele's methods for mediation & causal inference, from the 2015 book

> VanderWeele, T., 2015. Explanation in causal inference: methods for mediation and interaction. Oxford University Press.

Overview and changelog:

## Implemented
* Regression methods for continuous mediator/outcome (section 3.1.1)
* IPTW-based method for scenario with exposure-induced M->Y confounding (section 5.3.1)
* G-estimation based method for joint effect of multiple mediators (section 5.3.3)
* Calculation of random interventional analogues to direct and indirect effects (section 5.4.2)
* version 0.1.0 : Structural means models to deal with exposure-induced confounding of M->Y given a continuous exposure or mediator (section 5.3.5)
* version 0.1.1 : Handling of mice objects to add to resampler (slow but more stable than pooling from few imputed datasets)

## Installation

This package can be installed using `devtools::install_github`
