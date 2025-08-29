# Replicate Section 5 from Reguly (2025) there are the following steps:

1. Download the data from  Pop-Eleches & Urquiola (2013) from the following url: https://www.aeaweb.org/articles?id=10.1257/aer.103.4.1289\
2. Run `create_FE_outcomes.R` code to de-mean the outcome variables and save it to you path
3. `hetero_analysis_only2.m` replicates the tree estimate for Figure 2
4. `hetero_analysis_only2_forest.m` replicates the forest estimate for Figure 3.
5. `hetero_analysis_manycov.m` analyze the survey data.

Additional files and folders:
`pop_elches_one_feature.m` does the marginalization.
In the **appendix** folder one can find replications which are reported in the online appendix.