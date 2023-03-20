# Online repository for the PhenoPop method. 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7334379.svg)](https://doi.org/10.5281/zenodo.7334379)

Phenotypic deconvolution in heterogeneous cancer cell populations using drug screen data. 

This repository contains the Matlab implementation of the method described in the article [Phenotypic deconvolution in heterogeneous cancer cell populations using drug screening data]([url](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00028-0)), and also the necessary code and data to reproduce the results and figures in the publication. A Python package of the method is available at https://github.com/ocbe-uio/pyPhenoPop.

The required Matlab version is Matlab R2020b or newer. We are aware that the Matlab function **fmincon** from the Optimization Toolbox may cause runtime errors in versions newer than Matlab R2020b. To resolve this error, we suggest reinstalling the Optimization Toolbox and trying again (https://se.mathworks.com/matlabcentral/answers/473469-error-in-fmincon-undefined-getipoptions). 
The PhenoPop method is available as in the "Phenopop_inference" folder as the function "PhenoPop(experiment_name, NR, Concentration_array, Time_array, datafile)". 
See "example.m" in the "Phenopop_inference" folder or the Wiki pages for specifications on how to run the method. 

Author list: 
A. Köhn-Luque\*, E. M. Myklebust\*, D. S. Tadele, M. Giliberto, J. Noory, E. Harivel, P. Arsenteva, S.M. Mumenthaler, F. Schjesvold, K. Taskén, J. M. Enserink, K. Leder†, A. Frigessi†, and J. Foo†.

\*: These authors contributed equally.

†: These authors contributed equally.

Abstract: 
Tumor heterogeneity is an important driver of tumor recurrence, as treatments that initially elicit clinical responses can select for drug-tolerant tumor subpopulations, leading to the outgrowth of resistant clones and cancer treatment failure. Profiling the drug-response heterogeneity of tumor samples using traditional genomic deconvolution methods has yielded limited results, due in part to the imperfect mapping between genomic variation and functional characteristics. Here, by introducing an underlying population dynamic model of tumor subclonal response to therapy, we enable the phenotypic deconvolution of bulk drug response data into component subpopulations, and the estimation of their differential drug sensitivities and population frequencies. We used this method, called PhenoPop, to perform deconvolution on tumor drug screening data generated both experimentally and in silico. This study demonstrates how mechanistic population modeling can be leveraged to develop statistical frameworks for profiling phenotypic heterogeneity from bulk tumor samples and to perform individualized patient treatment predictions.
