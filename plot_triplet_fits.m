load test_triplet_fits.mat

N = 20; % number of units
nsamples = 10000; % number of training samples
N_iter = 50;

pmin_range = logspace(log10(.05),log10(.001),20);
ri_thresh_range = logspace(log10(.005),log10(10^-5),20);

%% First plot distributions of correlations and firing rates
% to make sure they are similar to Figure 2a,b

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

%% Plot Figure A2
% Compare ground truth triplets with their inferred parameters 
% in rm and ri models

figure, hold on


for j = 1:N_iter
    
    ix_triplets = find(sum(words_true{j},1) == 3);
    h_triplet_true = hs_true{j}(ix_triplets);
    
    for k = 1:length(ri_thresh_range)
        h = hs_ri{j,k};
        words = words_ri{j,k};
        
        h_triplet_ri = nan(size(h_triplet_true));
        for i = 1:length(ix_triplets)
            theword = words_true{j}(:,ix_triplets(i));
            ix_triplet_word = find(words'*theword==3 & sum(words',2)==3);
            if ~isempty(ix_triplet_word)
                h_triplet_ri(i) = h(ix_triplet_word);
            end
        end

        plot(h_triplet_true,-h_triplet_ri,'.','Color',[1,0,1])
        % Note difference in sign between ri and rm models
    end
end
for j = 1:N_iter
    
    ix_triplets = find(sum(words_true{j},1) == 3);
    h_triplet_true = hs_true{j}(ix_triplets);
    
    for k = 1:length(pmin_range)    
        h = hs_rm{j,k};
        words = words_rm{j,k};
        
        h_triplet_rm = nan(size(h_triplet_true));
        for i = 1:length(ix_triplets)
            theword = words_true{j}(:,ix_triplets(i));
            ix_triplet_word = find(words'*theword==3 & sum(words',2)==3);
            if ~isempty(ix_triplet_word)
                h_triplet_rm(i) = h(ix_triplet_word);
            end
        end
        
        plot(h_triplet_true,h_triplet_rm,'.b')
    end
end
set(gca,'FontSize',15)
xlabel('Triplet (true)')
ylabel('Triplet (model)')
