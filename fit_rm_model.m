% fit RM model with threshold p_min, data X
% E_true - empirical sample moments
% words - patterns corresponding to those moments
% h - interaction parameters

function [E_true,words,h] = fit_rm_model(X,p_min)

% MPF parameters
maxlinesearch = 10000; minf_options = []; 
minf_options.maxIter = maxlinesearch;
minf_options.maxFunEvals = maxlinesearch;

[N,~] = size(X);

% determine which moments are >p_min
[E_true, index] = calc_all_moments(X,p_min); 
% convert to spiking patterns
words = ind_to_words(index,N);

% fit Reliable Moment model with MPF
h_init = zeros(size(words,2),1);
h = minFunc(@K_dK_rm_L1,h_init,minf_options,words,X);
