# Reliable-Moment
Code for the Reliable Moment model (Cayco-Gajic et al. 2018)

Dependencies:
- Minimum Probability Flow learning: https://github.com/Sohl-Dickstein/Minimum-Probability-Flow-Learning
See also Sohl-Dickstein et al., 2011 
- Dichotomised Gaussian fitting: https://github.com/mackelab/CorBinian
See also Macke et al., 2009
- minFunc: http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html

To generate figures:
- Figure 1: Run `toy_model.m`
- Figures 2-3,A1: Run `test_pairwise_model.m`. to get `test_pairwise_model.mat`, then run `plot_pairwise_model.m` (2-4) or `plot_triplet_fits.m` (A1).
- Figure 4: Run `test_pairwise_model.m` and `calculate_dissimilarity_new_vs_old.m` to get `probs_new_vs_old.mat`, then run relevant sections of `plot_pairwise_model.m`
- Figure 5: Run `get_DG_parameters.m` to get `DG_parameters.mat`, then run `test_DG.m` and `calculate_dissimilarity_DG.m` to get `dissimilarities_DG.mat` and `test_DG.mat`, finally run `plot_DG.mat`
- Figure A2: Run `test_triplets.m` to get `test_triplets.mat`, then plot figure with `plot_triplet_fits.m`
