clear all; clc

addpath('/path/to/mackelab-pop_spike-a161a8ed79e0/util/')
addpath('/path/to/mackelab-pop_spike-a161a8ed79e0/dich_gauss/')

N = 20;
nsamples = 10000;
N_iter = 50;

%%
% Parameters for distribution
m = 5; v= 4;

Rs = zeros(N,N,N_iter);
gs = zeros(N,N_iter);
for iter = 1:N_iter
    % Get firing rates
    f_rates = lognrnd(log(m^2/sqrt(v+m^2)),sqrt(log(v/m^2+1)),N,1);
    mean_activity = f_rates*.02;

    rho_specified = normrnd(.1,.05,N,N); % specified correlation matrix
    rho_specified = (rho_specified+rho_specified')/2;
    rho_specified = rho_specified - diag(diag(rho_specified)) + eye(N,N);

    % Get specified covariance from specified correlation
    cov_specified = zeros(N,N);
    for i = 1:N
        var_i = mean_activity(i)-mean_activity(i)^2;
        for j = 1:N
            var_j = mean_activity(j)-mean_activity(j)^2;
            cov_specified(i,j) = rho_specified(i,j)*(sqrt(var_i) * sqrt(var_j));
        end
    end

    [g,R] = sampleDichGauss01(mean_activity,cov_specified);
    
    Rs(:,:,iter) = R;
    gs(:,iter) = g;
    
end

save('DG_parameters','Rs','gs')

