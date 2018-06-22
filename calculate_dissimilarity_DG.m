% Calculateds d(P_true,F_model) for RM and RI models for DG data

d = @(P,Q) nansum(P.*abs(log2(P./Q)));

d_test_rm = zeros(N_iter,length(pmin_range));
d_test_ri = zeros(N_iter,length(pmin_range));

for i = 1:N_iter
    i
    X_train = Xs{i};
    
    R = Rs(:,:,i); g = gs(:,i);
    t = R' * randn(N,nsamples);
    X_test = double(t>-repmat(g,1,nsamples).*ones(N,nsamples));
    [P_train,words_train,M_train] = get_empirical_probs(X_train);
    [P_test,words_test,~] = get_empirical_probs(X_test);

    for j_rm = 1:length(pmin_range)
        Z_rm = get_Z_GT(words_train,hs_rm{i,j_rm},words_rm{i,j_rm},M_train);
        P_rm_test = get_rm_probs(words_test,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);
        
        d_test_rm(i,j_rm) = d(P_test,P_rm_test);
    end

    for j_ri = 1:length(ri_thresh_range)
        f_ri_test = get_ri_probs(words_test,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});
        
        d_test_ri(i,j_ri) = d(P_test,f_ri_test);
    end
end

save('dissimilarities_DG','d_test_ri','d_test_rm')