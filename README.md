# ARC
Accumulation Rate Curve Estimator

data_sim.r
Code to generate data for the simulation scenarios.

main_run.r
Code to generate the point estimate using the ARC estimator, estimate count covariances and bootstrap to obtain confidence intervals.

point_est.r
Functions to calculate the point estimate using the ARC method.

var_est.r
Functions to calculate the variance of the y_k counts for the ARC method.

data_gen.r
Helper function to process data to use as inputs for the ARC method.

wagner_tree_microbiome_cleaned_example_1.RData
The OTU table used in Applied Example 1 from: Wagner MR, Lundberg DS, del Rio TG, Tringe SG, Dangl JL, Mitchell-Olds T (2016) Host genotype and age shape the leaf and root microbiomes of wild perennial plant. Nat Commun. https​://doi.org/10.1038/ncomm​s1215​1
Original data available at: https://datadryad.org/stash/dataset/doi:10.5061/dryad.g60r3

steam_catalogue_cleaned_example_2.RData
The catalogue data used in Applied Example 2.
