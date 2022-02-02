# CoastalRisk
Estimates the impact of water risks on crop production in 6 coastal Bangladesh polders


                                  Coastal Risk Model

Risk Modelv20 developed by E.J. Barbour October 2021 with contributions in study design from J.W. Hall, E. Borgomeo, M.S.G. Adnan, M. Salehin, M.S.A Khan, and K. Paprocki. Results from this code are published in the manuscript: "The unequal distribution of water risks and adaptation benefits in coastal Bangladesh". The code draws on Matlab code developed by Borgomeo, E., Hall, J., and Salehin, M. 2018. "Avoiding the water-poverty trap: insights from a conceptual human-water dynamical model for coastal Bangladesh". International Journal of Water Resources Development. 34 (6) 900-922. https://doi.org/10.1080/07900627.2017.1331842.

It contains the following: 
1. RiskModelv20.f90: Main Program
2. poisson_fixed_time.f90: generates stochastic flood events
3. agri_mauza_module.f90: calculates crop production based on salinity, waterlogging, crop yields/ratios
4. crop_matrix_module.f90: generates matrix of crop yields for entire time series
5. sediment.f90: calculates new elevation for each mauza based on subsidence (and deposition if using).
