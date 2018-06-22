% Load data and parameters
load test_pairwise_model.mat

N = 20; % number of units
nsamples = 10000; % number of training samples
N_iter = 50;

pmin_range = logspace(log10(.05),log10(.001),20);
ri_thresh_range = logspace(log10(.005),log10(10^-5),20);

%% Plot Figure 2a,b

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
%% Plot figure 2c-e
% Specific example
j = 1;
k = 20;

X_train = Xs{j};
J = Js{j};

% Get model data
h = hs_rm{j,k};
words = words_rm{j,k};

disp(strcat('pmin:',num2str(pmin_range(k))))
disp(strcat('Number of fitted pairs:',num2str(sum(sum(words,1)==2))))
disp(strcat('Number of fitted triplets:',num2str(sum(sum(words,1)==3))))
disp(strcat('Number of fitted quadruplets:',num2str(sum(sum(words,1)==4))))
disp(strcat('Number of fitted quintuplets+:',num2str(sum(sum(words,1)>4))))

X_test = sample_ising( J, nsamples,100*N,10*N);
X_rm = sample_rm(h,words,N,10000,100*N,10*N);

% Figure 2c
% Empirical (training) vs. fitted pairwise moments
x=[]; y=[];
for i = 1:N
    for j2 = (i+1):N
        y = [y,mean(X_rm(i,:).*X_rm(j2,:))];
        x = [x,mean(X_train(i,:).*X_train(j2,:))];
    end
end
figure, loglog(x,y,'.r','MarkerSize',10)
% Empirical (training) vs. fitted moments of all orders
x=[]; y=[];
for i = 1:size(words,2)
    cells = find(words(:,i));
    y = [y,mean(prod(X_rm(cells,:)))];
    x = [x,mean(prod(X_train(cells,:)))];
end
hold on, loglog(x,y,'.b','MarkerSize',10)
hold on, plot([10^-4,.1],[10^-4,.1],'k')
set(gca,'FontSize',15)
xlabel('Moments (training data)')
ylabel('Moments (model)')

% Figure 2d - probabilities of spiking patterns
[~,words_train,M_train] = get_empirical_probs(X_train);
Z_rm = get_Z_GT(words_train,h,words,M_train);

[P_test,words_test,~] = get_empirical_probs(X_test);
P_rm_test = get_rm_probs(words_test,h,words,Z_rm);

figure, loglog(P_test,P_rm_test,'.b','MarkerSize',10)
hold on, plot([10^-4,1],[10^-4,1],'k')
set(gca,'FontSize',15)
xlabel('Probability (test data)')
ylabel('Probability (model)')

% Figure 2e - predicted correlations
% All correlations
C = corrcoef(X_rm'); C=C(:);
C_true = corrcoef(X_test'); C_true=C_true(:);
figure, plot(C_true(C~=1),C(C~=1),'.r','MarkerSize',10)
hold on, plot([-.1,.2],[-.1,.2],'k')
set(gca,'FontSize',15)
xlabel('Correlations (test data)')
ylabel('Correlations (model)')

% Correlations for moments fitted by the model
x = [];% Correlations from test data
y = [];% correlations from RM model data
ix=find(sum(words,1)==2);
C = corrcoef(X_rm');
C_true = corrcoef(X_test'); 
for i = 1:length(ix)
    cells = find(words(:,ix(i)));
    y = [y,C(cells(1),cells(2))];
    x = [x,C_true(cells(1),cells(2))];
end
hold on, plot(x,y,'.b','MarkerSize',10)
axis([-.1,.2,-.1,.2])

% Figure 2e, inset - predicted firing rates
% note .02 corresponds to 20 ms bins
figure, plot(mean(X_test,2)/.02,mean(X_rm,2)/.02,'.b','MarkerSize',10)
hold on, plot([0,6],[0,6],'k')
set(gca,'FontSize',15)
xlabel('Firing rate, Hz (test data)')
ylabel('Firing rate, Hz (model)')

%% Plot Figure 3 - spurious correlations

num_hocs_ri = zeros(N_iter,length(pmin_range));
h_mean_ri = zeros(N_iter,length(pmin_range));
h_std_ri = zeros(N_iter,length(pmin_range));

num_hocs_rm = zeros(N_iter,length(pmin_range));
h_mean_rm = zeros(N_iter,length(pmin_range));
h_std_rm = zeros(N_iter,length(pmin_range));

for j = 1:N_iter
    
    for k = 1:length(pmin_range)    
        h = hs_rm{j,k};
        words = words_rm{j,k};
        which_hois = find(sum(words,1)>2);
                
        num_hocs_rm(j,k) = length(which_hois);
        h_mean_rm(j,k) = mean(abs(h(which_hois)));
        h_std_rm(j,k) = std(h(which_hois));
    end
    
    for k = 1:length(ri_thresh_range)
        h = hs_ri{j,k};
        words = words_ri{j,k};
        which_hois = find(sum(words,1)>2);
        
        num_hocs_ri(j,k) = length(which_hois);
        h_mean_ri(j,k) = mean(abs(h(which_hois)));
        h_std_ri(j,k) = std(h(which_hois));
    end
end
figure, semilogx(num_hocs_ri(:),h_mean_ri(:),'.','Color',[1,0,1])
hold on, semilogx(num_hocs_rm(:),h_mean_rm(:),'.b')
hold on, plot(num_hocs_ri(1,:),h_mean_ri(1,:),'.-','Color',[1,0,1])
hold on, plot(num_hocs_rm(1,:),h_mean_rm(1,:),'.-b')
set(gca,'FontSize',15)
xlabel('Number fitted HOIs')
ylabel('Average magnitude')

figure, semilogx(num_hocs_ri(:),h_std_ri(:),'.m')
hold on, semilogx(num_hocs_rm(:),h_std_rm(:),'.b')
hold on, plot(num_hocs_ri(1,:),h_std_ri(1,:),'.-','Color',[1,0,1])
hold on, plot(num_hocs_rm(1,:),h_std_rm(1,:),'.-b')
set(gca,'FontSize',15)
set(gca,'FontSize',15)
xlabel('Number fitted HOIs')
ylabel('Standard deviation')
%% Appendix B Figure A1
% Normalize hoi's by magnitude of average pair
num_hocs_ri = zeros(N_iter,length(pmin_range));
num_params_ri = zeros(N_iter,length(pmin_range));
h_mean_pair_ri = zeros(N_iter,length(pmin_range));
h_mean_ratio_ri = zeros(N_iter,length(pmin_range));

num_hocs_rm = zeros(N_iter,length(pmin_range));
num_params_rm = zeros(N_iter,length(pmin_range));
h_mean_pair_rm = zeros(N_iter,length(pmin_range));
h_mean_ratio_rm = zeros(N_iter,length(pmin_range));

for j = 1:N_iter
    
    for k = 1:length(pmin_range)    
        h = hs_rm{j,k};
        words = words_rm{j,k};
        which_hois = find(sum(words,1)>2);
        
        which_hois_pair = find(sum(words,1)==2);
        
        num_hocs_rm(j,k) = length(which_hois);
        num_params_rm(j,k) = size(words,2);
        h_mean_pair_rm(j,k) = mean(abs(h(which_hois_pair))); % avg h for pairs
        h_mean_ratio_rm(j,k) = mean(abs(h(which_hois)))/h_mean_pair_rm(j,k);
    end
    
    for k = 1:length(ri_thresh_range)
        h = hs_ri{j,k};
        words = words_ri{j,k};
        which_hois = find(sum(words,1)>2);

        which_hois_pair = find(sum(words,1)==2);
        
        num_hocs_ri(j,k) = length(which_hois);
        num_params_ri(j,k) = size(words,2);
        h_mean_pair_ri(j,k) = mean(abs(h(which_hois_pair)));
        h_mean_ratio_ri(j,k) = mean(abs(h(which_hois)))/h_mean_pair_ri(j,k);
    end
end
figure, semilogx(num_hocs_ri(:),h_mean_pair_ri(:),'.','Color',[1,0,1])
hold on, semilogx(num_hocs_rm(:),h_mean_pair_rm(:),'.b')
hold on, plot(num_hocs_ri(1,:),h_mean_pair_ri(1,:),'.-','Color',[1,0,1])
hold on, plot(num_hocs_rm(1,:),h_mean_pair_rm(1,:),'.-b')
set(gca,'FontSize',15)
xlabel('Number fitted HOIs')
ylabel('Average magnitude (pairs)')

figure, semilogx(num_hocs_ri(:),h_std_ri(:),'.m')
hold on, semilogx(num_hocs_rm(:),h_std_rm(:),'.b')
hold on, plot(num_hocs_ri(1,:),h_std_ri(1,:),'.-','Color',[1,0,1])
hold on, plot(num_hocs_rm(1,:),h_std_rm(1,:),'.-b')
set(gca,'FontSize',15)
set(gca,'FontSize',15)
xlabel('Number fitted HOIs')
ylabel('Average magnitude of HOIs (norm. by pairs)')


%% Figure 4b,d
load probs_new_vs_old.mat
i=1:50
figure, loglog(num_hocs_ri(:),d_new_ri(:),'.m','MarkerSize',10)
hold on, loglog(num_hocs_rm(:),d_new_rm(:),'.b','MarkerSize',10)
loglog(num_hocs_ri(1,:),d_new_ri(1,:),'.-m','MarkerSize',15,'LineWidth',2)
loglog(num_hocs_rm(1,:),d_new_rm(1,:),'.-b','MarkerSize',15,'LineWidth',2)
set(gca,'FontSize',15)
xlabel('Number fitted HOIs')
ylabel('d(P_{true},P_{model})')
title('New words')

figure, loglog(num_hocs_ri(:),d_old_ri(:),'.m','MarkerSize',10)
hold on, loglog(num_hocs_rm(:),d_old_rm(:),'.b','MarkerSize',10)
loglog(num_hocs_ri(1,:),d_old_ri(1,:),'.-m','MarkerSize',15,'LineWidth',2)
loglog(num_hocs_rm(1,:),d_old_rm(1,:),'.-b','MarkerSize',15,'LineWidth',2)
set(gca,'FontSize',15)
xlabel('Number fitted HOIs')
ylabel('d(P_{true},P_{model})')
title('Old words')
%% Figure 4a,c
% Specific Example with same number params
i=19;
j_rm = 18; j_ri=10;
num_params_rm(i,j_rm), num_params_ri(i,j_ri)


X_train = Xs{i};
X_test = sample_ising( Js{i}, nsamples,100*N,10*N);

[P_train,words_train,M_train] = get_empirical_probs(X_train);
[P_test,words_test,~] = get_empirical_probs(X_test);

[words_common,ix_common,~] = intersect(words_test',words_train','rows');
[words_new,ix_new] = setdiff(words_test',words_train','rows');

words_common = words_common'; 
words_new = words_new'; 

Z_pw = get_Z_pw(Js{i});
P_test = get_pw_probs(words_test,Js{i},Z_pw);
P_old = get_pw_probs(words_common,Js{i},Z_pw);
P_new = get_pw_probs(words_new,Js{i},Z_pw);

Z_rm = get_Z_GT(words_train,hs_rm{i,j_rm},words_rm{i,j_rm},M_train);
P_rm_test = get_rm_probs(words_test,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);
P_rm_old = get_rm_probs(words_common,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);
P_rm_new = get_rm_probs(words_new,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);


f_ri_test = get_ri_probs(words_test,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});
f_ri_old = get_ri_probs(words_common,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});
f_ri_new = get_ri_probs(words_new,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});


figure, loglog(P_old,f_ri_old,'.m','MarkerSize',10)
hold on,loglog(P_old,P_rm_old,'.b','MarkerSize',10)
set(gca,'FontSize',15)
xlabel('Frequency (ground truth)')
ylabel('Frequency (model)')
title('Old patterns')

figure, loglog(P_new,f_ri_new,'.m','MarkerSize',10)
hold on,loglog(P_new,P_rm_new,'.b','MarkerSize',10)
set(gca,'FontSize',15)
xlabel('Frequency (ground truth)')
ylabel('Frequency (model)')
title('New patterns')



