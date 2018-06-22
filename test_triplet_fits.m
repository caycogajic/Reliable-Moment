% This example fits the rm and ri models to 
% a pairwise model with sparse triplets

clear all

addpath('/path/to/minFunc')

% First generate pairwise max ent data with random interactions
N = 20; % number of units
nsamples = 10000; % number of training samples
N_iter = 50;

pmin_range = logspace(log10(.05),log10(.001),20);
ri_thresh_range = logspace(log10(.005),log10(10^-5),20);

%%
Xs = cell(N_iter,1);
hs_true = cell(N_iter,1);
words_true = cell(N_iter,1);

hs_rm = cell(N_iter,length(pmin_range));
words_rm = cell(N_iter,length(pmin_range));

hs_ri = cell(N_iter,length(ri_thresh_range));
Zs_ri = cell(N_iter,length(ri_thresh_range));
words_ri = cell(N_iter,length(ri_thresh_range));

for i = 1:N_iter
    i
    % choose a random coupling matrix to generate the test data
    % parameters chosen to make similar to cortical data
    % ie low firing rate and range of correlations
    J = 2*randn(N,N)/sqrt(N); J = (J+J')/2;
    J = J-diag(diag(J)); 
    bias = randn( N, 1 )*.5+3;
    J = J+diag(bias);

    % Convert to RM style
    words = zeros(N,N+nchoosek(N,2));
    h_true = zeros(N+nchoosek(N,2),1);
    words(1:N,1:N) = eye(N);
    h_true(1:N) = bias;
    ix = N+1;
    for j = 1:N
        for k = (j+1):N
            words(j,ix) = 1;
            words(k,ix) = 1;
            h_true(ix) = 2*J(j,k);
            ix = ix+1;
        end
    end
    
    % Add sparse network of triplet correlations
    p = .05;
    for j = 1:N
        for k = (j+1):N
            for l = (k+1):N
                if rand < p
                    words(j,ix) = 1;
                    words(k,ix) = 1;
                    words(l,ix) = 1;
                    h_true(ix) = randn;
                    ix = ix+1;
                end
            end
        end
    end
    
    X = sample_rm(h_true,words,N,10000,100*N,10*N);
    
    words_true{i} = words;
    hs_true{i} = h_true;

    %for different ranges of pmin, calculate RM model
    for j = 1:length(pmin_range)
        [E_true,words,h] = fit_rm_model(X,pmin_range(j));
        % Save trained info
        hs_rm{i,j} = h;
        words_rm{i,j} = words;
    end

    %for different ranges of RIthresh, calculate RI model
    for j = 1:length(ri_thresh_range)
        [words,h_ri,Z_ri] = fit_RI(X,ri_thresh_range(j));
        % Save trained info
        hs_ri{i,j} = h_ri;
        Zs_ri{i,j} = Z_ri;
        words_ri{i,j} = words;
    end
    
    % Save true data paramters
    Xs{i} = X; 
end

save('test_triplet_fits.mat','pmin_range','ri_thresh_range','Xs','hs_true','words_true','words_rm','hs_rm','hs_ri','Zs_ri','words_ri')
