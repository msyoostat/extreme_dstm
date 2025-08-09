# Intro
This is the Github repo for the paper "Modeling high and low extremes with a novel dynamic spatio-temporal model" by Myungsoo Yoo, Likun Zhang, Christopher K. Wikle, and Thomas Opitz. The preprint is available on arXiv: https://arxiv.org/abs/2508.01481.

In this work, we propose a physically motivated spatio-temporal dynamic model in which the innovations switch between heavy- and light-tailed distributions according to a probability, capturing external drivers of extreme events. The framework accounts for spatial and temporal dependence without assuming temporal independence and allows flexible extremal dependence across space and time. It also provides a fully model-based approach for detecting extreme events with uncertainty quantification. We illustrate the model using hourly PM2.5 data over the Central U.S. in March 2024.

# real_data
## data
- pm2_5_hour_5.RData: This dataset contains PM2.5 measurements used in the paper, focusing on Arkansas, Oklahoma, and southwestern Missouri from March 19 at 00:30 to March 30 at 23:30. The original data are from NASA’s GEOS-CF v1.0 system (Keller et al., 2021).
- hourly3_basis_170_4_2_crs.RData: This dataset contains the local basis functions used for PM2.5 dataset. 
## code
- pm_2_5_full_stable.R: 
- pm_2_5_full_vg.R: 
- pm_2_5_missing_stable.R: 


# simulation

# Reference

C. A. Keller, K. E. Knowland, B. N. Duncan, J. Liu, D. C. Anderson, S. Das, R. A. Lucchesi, E. W. Lundgren, J. M. Nicely, E. Nielsen, L. E. Ott, E. Saunders, S. A. Strode, P. A. Wales, D. J. Jacob, and S. Pawson. Description of the nasa geos composition forecast modeling system geos-cf v1.0. *Journal of Advances in Modeling Earth Systems*, 13(4): e2020MS002413, 2021. doi: https://doi.org/10.1029/2020MS002413
