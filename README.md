## causamed
Implementation of Tyler Vanderweele's methods for mediation

####Everything pulled straight from the book "Explanation in Causal Inference; Methods for mediation and interaction" by Tyler Vanderweele (2015)

Currently implements
1. Regression methods for continuous mediator/outcome (section 3.1.1) <br>
2. IPTW-based method for scenario with exposure-induced M->Y confounding (section 5.3.1) <br>
3. G-estimation based method for joint effect of multiple mediators (section 5.3.3) <br>
4. Calculation of random interventional analogues to direct and indirect effects (section 5.4.2)

To do (not much work)
5. Structural means models to deal with exposure-induced confounding of M->Y given a continuous exposure or mediator
6. Extension of methods to categorical and continuous outcomes where not already the case
7. Handling of mice objects
8. Regression-based methods for multiple confounders
9. Fix (?) random interventional analogue for categorical exposure case

To do (lots of work)
9. Implementation of sensitivity analyses for regression-based methods
10. More sensitivity analyses in settings where NDE/NIE not identified
11. Handle non-parametric/semi-parametric approachs to propensity score estimation for weighting-based approaches
