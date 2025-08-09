# Intro
This is the Github repo for the paper "Modeling high and low extremes with a novel dynamic spatio-temporal model" by Myungsoo Yoo, Likun Zhang, Christopher K. Wikle, and Thomas Opitz. The preprint is available on arXiv: https://arxiv.org/abs/2508.01481.

In this work, we propose a physically motivated spatio-temporal dynamic model in which the innovations switch between heavy- and light-tailed distributions according to a probability, capturing external drivers of extreme events. The framework accounts for spatial and temporal dependence without assuming temporal independence and allows flexible extremal dependence across space and time. It also provides a fully model-based approach for detecting extreme events with uncertainty quantification. We illustrate the model using hourly PM2.5 data over the Central U.S. in March 2024.

# real_data
## data
- pm2_5_hour_5.RData: This dataset contains PM2.5 measurements used in the paper, focusing on Arkansas, Oklahoma, and southwestern Missouri from March 19 at 00:30 to March 30 at 23:30. The original data are from NASA’s GEOS-CF v1.0 system (Keller et al., 2021).
- hourly3_basis_170_4_2_crs.RData: This dataset contains the local basis functions used for PM2.5 dataset. 
## code
- pm_2_5_full_stable.R: This is the code to fit the "stable" model described in Section 7.1.
- pm_2_5_full_vg.R: This is the code to fit the "variance-gamma" model described in Section 7.1.
- pm_2_5_missing_stable.R: This is the code to fit the "stable" model described in Section 7.2.


# simulation
## detection
### data
### code


# Reference

Keller, C. A., Knowland, K. E., Duncan, B. N., Liu, J., Anderson, D. C., Das, S., Lucchesi, R. A., Lundgren, E. W., Nicely, J. M., Nielsen, E., Ott, L. E., Saunders, E., Strode, S. A., Wales, P. A., Jacob, D. J., & Pawson, S. (2021). Description of the NASA GEOS Composition Forecast Modeling System GEOS-CF v1.0. Journal of Advances in Modeling Earth Systems, 13(4), e2020MS002413. https://doi.org/10.1029/2020MS002413
