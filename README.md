# Sample Scripts

The purpose of this repository is to share a couple examples of my past work on analytic projects. 

Because the scripts contained herein are analyses of health data, I cannot share the data in this public format and I have redacted information linking the scripts to a specific project.  

Below is a summary of the key techniques or methods applied in each script.


## Example A: Mediation Analysis with 4-Way Decomposition of Effects and Bootstrapped CIs

This script examined factors mediating the association between poverty or food insecurity and a suite of neurodevelopmental outcomes in children aged 3-5 years. The analyis implemented VanderWeele's (2014) 4-way decomposition of effects (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4220271/) and bootstrapping methods to estimate confidence intervals. Data recoding and manipulation steps are also included, as well as scripts for plots. 

Note that if opened in R Studio, this script is navigable through the drop-down menu at the bottom of the script window.


## Example B: Longitudinal Analaysis with Generalized Estimating Equations (GEE)

This script includes a series of analyses to investigate the association between insecticide exposure and antibody titres in children, and whether that association was modified by factors such as child sex, poverty, and food insecurity (i.e., Effect Measure Modification). To investigate the effect longitudinally, generalized estimating equations (GEE) were applied. Confidence intervals were bootstrapped. The script also includes investigation and plotting of non-linear dose response splines. Finally, data recoding and manipulation steps are also included. 

Note that if opened in R Studio, this script is navigable through the drop-down menu at the bottom of the script window.



