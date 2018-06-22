load DG_parameters.mat

pmin_range = logspace(log10(.05),log10(.001),20);
ri_thresh_range = logspace(log10(.005),log10(10^-5),20);

load test_DG.mat
N=20;
nsamples = 10000; % number of training samples
N_iter = 50;

load dissimilarities_DG.mat

%% Figure 5a,b
% Plot distribution of firing rates and correlations

% First plot distributions of correlations and firing rates
J = zeros(N,N);
for i = 1:N
    for j = (i+1):N
        J(i,j) = 1;
    end
end

fr_all = zeros(N,N_iter);
C_all = zeros(nchoosek(N,2),N_iter);
for k = 1:N_iter
    fr_all(:,k) = mean(Xs{k},2)/.02;
    C = corrcoef(Xs{k}');
    C_all(:,k) = C(J==1); 
end

figure, histogram(fr_all(:))
set(gca,'FontSize',15)
xlabel('Firing rate (Hz)')
ylabel('Number')

figure, histogram(C_all(:))
set(gca,'FontSize',15)
xlabel('Correlation')
ylabel('Number')


%% Figure 5c,e
% plot d as a function of number parameters
num_params_ri = zeros(N_iter,length(pmin_range));
d_mean_ri = zeros(N_iter,length(pmin_range));
d_std_ri = zeros(N_iter,length(pmin_range));

num_params_rm = zeros(N_iter,length(pmin_range));
d_mean_rm = zeros(N_iter,length(pmin_range));
d_std_rm = zeros(N_iter,length(pmin_range));

for j = 1:N_iter
    
    for k = 1:length(pmin_range)    
        h = hs_rm{j,k};
        words = words_rm{j,k};
        which_hois = find(sum(words,1)>2);
        num_params_rm(j,k) = size(words,2);
        d_mean_rm(j,k) = mean(h(sum(words,1)>1));
    end
    
    for k = 1:length(ri_thresh_range)
        h = hs_ri{j,k};
        words = words_ri{j,k};
        which_hois = find(sum(words,1)>2);
        num_params_ri(j,k) = size(words,2);
        d_mean_ri(j,k) = mean(h(sum(words,1)>1));
    end
end

% Figure 5c
figure, semilogx(num_params_ri(:),d_test_ri(:),'.m')
hold on, semilogx(num_params_rm(:),d_test_rm(:),'.b')
hold on, plot(num_params_ri(12,:),d_test_ri(12,:),'.-','Color',[1,0,1],'LineWidth',2,'MarkerSize',15)
hold on, plot(num_params_rm(12,:),d_test_rm(12,:),'.-b','LineWidth',2,'MarkerSize',15)
set(gca,'FontSize',15)
xlim([10,20000])
xlabel('Number parameters')
ylabel('d(P_{true},F_{model})')

% Figure 5e
figure, loglog(num_params_ri(:),time_ri(:),'.m')
hold on, loglog(num_params_rm(:),time_rm(:),'.b')
hold on, plot(num_params_ri(12,:),time_ri(12,:),'.-','Color',[1,0,1],'LineWidth',2,'MarkerSize',15)
hold on, plot(num_params_rm(12,:),time_rm(12,:),'.-b','LineWidth',2,'MarkerSize',15)
set(gca,'FontSize',15)
xlim([10,20000])
xlabel('Number parameters')
ylabel('Time (s)')
%% Figure 5d
% Compare two specific cases with same number of parameters

i=12;
j_rm = 8; j_ri=9;
num_params_rm(i,j_rm), num_params_ri(i,j_ri)

X_train = Xs{i};
R = Rs(:,:,i); g = gs(:,i);
t = R' * randn(N,nsamples);
X_test = double(t>-repmat(g,1,nsamples).*ones(N,nsamples));
    
[P_train,words_train,M_train] = get_empirical_probs(X_train);
[P_test,words_test,~] = get_empirical_probs(X_test);

Z_rm = get_Z_GT(words_train,hs_rm{i,j_rm},words_rm{i,j_rm},M_train);
P_rm_test = get_rm_probs(words_test,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);

f_ri_test = get_ri_probs(words_test,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});

figure, loglog(P_test,f_ri_test,'.m','MarkerSize',10)
hold on,loglog(P_test,P_rm_test,'.b','MarkerSize',10)
hold on, plot([10^-6,1],[10^-6,1],'k')
set(gca,'FontSize',15)
xlabel('Probability (test)')
ylabel('Probability (model)')
