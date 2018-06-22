% This example fits the rm and ri models to 
% a DG model

clear all

addpath('/path/to/minFunc')

load DG_parameters.mat

% First generate pairwise max ent data with random interactions
N = 20; % number of units
nsamples = 10000; % number of training samples
N_iter = 50;

pmin_range = logspace(log10(.05),log10(.001),20);
ri_thresh_range = logspace(log10(.005),log10(10^-5),20);

%%
Xs = cell(N_iter,1);
Js = cell(N_iter,1);

hs_rm = cell(N_iter,length(pmin_range));
words_rm = cell(N_iter,length(pmin_range));

hs_ri = cell(N_iter,length(ri_thresh_range));
Zs_ri = cell(N_iter,length(ri_thresh_range));
words_ri = cell(N_iter,length(ri_thresh_range));

time_rm = zeros(N_iter,length(pmin_range));
time_ri = zeros(N_iter,length(ri_thresh_range));

for i = 1:N_iter
    i
    % generate data to be fit from pre-loaded DG parameters
    R = Rs(:,:,i); g = gs(:,i);
    t = R' * randn(N,nsamples);
    X = double(t>-repmat(g,1,nsamples).*ones(N,nsamples));
    
    %[Ps_emp{i},words_emp,M] = get_empirical_probs(X_te);

    %for different ranges of pmin, calculate RM model
    parfor j = 1:length(pmin_range)
        tic
        [E_true,words,h] = fit_rm_model(X,pmin_range(j));
        % Save trained info
        hs_rm{i,j} = h;
        words_rm{i,j} = words;
        time_rm(i,j)=toc;
    end

    %for different ranges of RIthresh, calculate RI model
    parfor j = 1:length(ri_thresh_range)
        tic
        [words,h_ri,Z_ri] = fit_RI(X,ri_thresh_range(j));
        % Save trained info
        hs_ri{i,j} = h_ri;
        Zs_ri{i,j} = Z_ri;
        words_ri{i,j} = words;
        time_ri(i,j)=toc;
    end
    
    % Save true data paramters
    Xs{i} = X; 

end

save('test_DG.mat','pmin_range','ri_thresh_range','Xs','words_rm','hs_rm','hs_ri','Zs_ri','words_ri','time_rm','time_ri')
