## causamed
Implementation of Tyler Vanderweele's methods for mediation

####Everything pulled straight from the book "Explanation in Causal Inference; Methods for mediation and interaction" by Tyler Vanderweele (2015)

Currently implements <br>
1. Regression methods for continuous mediator/outcome (section 3.1.1) <br>
2. IPTW-based method for scenario with exposure-induced M->Y confounding (section 5.3.1) <br>
3. G-estimation based method for joint effect of multiple mediators (section 5.3.3) <br>
4. Calculation of random interventional analogues to direct and indirect effects (section 5.4.2) <br> <br>
5. Structural means models to deal with exposure-induced confounding of M->Y given a continuous exposure or mediator (section 5.3.5) <br>

To do (not much work) <br>
6. Extension of methods to categorical and continuous outcomes where not already the case <br>
7. Handling of mice objects <br>
8. Regression-based methods for multiple confounders <br>
9. Fix (?) random interventional analogue for categorical exposure case <br> <br>

To do (lots of work) <br> 
10. Implementation of sensitivity analyses for regression-based methods <br>
11. More sensitivity analyses in settings where NDE/NIE not identified <br>
12. Handle non-parametric/semi-parametric approachs to propensity score estimation for weighting-based approaches <br>
