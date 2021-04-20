# dFe

This includes all the data from the GEOTRACES
IDP 2017, plus
historical compilations by Moore and Braucher,2008 and Tagliabue et al,
2012.

The format in this text file is
month, lat, lon, depth, dfe, gridx, gridy, gridlevel, inregion

gridx,gridy,gridlevel, refer to the gx3v7 grid, and just ignore
inregion,

Schlitzer, R., Anderson, R. F., Masferrer Dodas, E, et al., The
GEOTRACES Intermediate Data Product 2017, Chem. Geol. (2018),
https://doi.org/10.1016/j.chemgeo.2018.05.040.

Moore, J. K., & Braucher, O. (2008). Sedimentary and mineral dust
sources of dissolved iron to the world ocean. Biogeosciences, 5(3),
631-656.

Tagliabue, A., Mtshali, T., Aumont, O., Bowie, A. R., Klunder, M. B.,
Roychoudhury, A. N., & Swart, S. (2012). A global compilation of
dissolved iron measurements: Focus on distributions and processes in the
Southern Ocean. Biogeosciences, 9(6), 2333-2349.



# DOM
1) This (DOMdatabaseDec2012_gx3v7) is a dataset of (in order of tab delimit):  month, lon, lat, depth, [DON], [DOP], [DOC], pop_gx3v7 grid lon, pop_gx3v7 grid lat, pop grid_level, and an integer that assigns a data point to a ‘region’ eg. North Atlantic for color-coded plotting. It is a tab delimited text file that we use with Keith M’s IDL plotting routines. 

2) Yes µM (mmol/m3) for [DON], [DOP], [DOC]. These are bulk pools. Missing value is denoted by –1.0

3) To compare MARBL DOM pools to the obs dataset you should sum DOC + DOCr (e.g. semilable + refractory) to get bulk DOC, etc. for DON and DOP. Then you can compute C:P and N:P ratios of obs and output as desired. 

4) Yes happy to help. I had a quick read of the section on detrital OM cycling. There is also a small flux from POM remin to the refractory DOM pools. Semilabile DOP has an additional sink due to direct autotrophic utilization as a P source. I can add details on these as comments to the manuscript. 

P.S. I’m attaching the .csv file that was used to create the database file that Keith M. sent. Here missing values are denoted with -999. You will see some data points that have a month value of ‘0’. These correspond to data from HOT and BATS in which I include 1 time-series averaged depth profile from each station.

P.P.S. If you want to compare simulated semilable DOM to obs slDOM (or refractory to refractory), you can obtain this split from the obs data by: 
1) subtracting 1.8 µM from obs DON (i.e. DONr obs = 1.8 µM) 
2) subtracting 0.03 µM from obs DOP
3) subtracting a value between 38.5 – 47.1 µM for DOCr depending on ocean basin (see Table S2 from Letscher & Moore, 2015 GBC). 


2015. Letscherx, R.T., and J.K. Moore, Preferential remineralization of
dissolved organic phosphorus and non-Redfield DOM dynamics in the global
ocean: Impacts on marine productivity, nitrogen fixation, and carbon
export, Global Biogeochemical Cycles, 29, 325-340,
doi:10.1002/2014GB004904.