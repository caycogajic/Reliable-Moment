% Calculateds d(P_true,F_model) for RM and RI models
% both for new and old spiking patterns
% rather slow because calculates exact ground truth probabilities

d = @(P,Q) nansum(P.*abs(log2(P./Q)));

d_new_rm = zeros(N_iter,length(pmin_range));
d_old_rm = zeros(N_iter,length(pmin_range));

d_new_ri = zeros(N_iter,length(pmin_range));
d_old_ri = zeros(N_iter,length(pmin_range));

for i = 1:N_iter

    i
    
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

    for j_rm = 1:length(pmin_range)
        Z_rm = get_Z_GT(words_train,hs_rm{i,j_rm},words_rm{i,j_rm},M_train);
        P_rm_test = get_rm_probs(words_test,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);
        P_rm_old = get_rm_probs(words_common,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);
        P_rm_new = get_rm_probs(words_new,hs_rm{i,j_rm},words_rm{i,j_rm},Z_rm);
        
        d_new_rm(i,j_rm) = d(P_new,P_rm_new);
        d_old_rm(i,j_rm) = d(P_old,P_rm_old);
    end

    for j_ri = 1:length(ri_thresh_range)
        f_ri_test = get_ri_probs(words_test,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});
        f_ri_old = get_ri_probs(words_common,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});
        f_ri_new = get_ri_probs(words_new,hs_ri{i,j_ri},words_ri{i,j_ri},Zs_ri{i,j_ri});
        
        d_new_ri(i,j_ri) = d(P_new,f_ri_new);
        d_old_ri(i,j_ri) = d(P_old,f_ri_old);
    end
end

save('probs_new_vs_old','d_new_ri','d_old_ri','d_new_rm','d_old_rm')