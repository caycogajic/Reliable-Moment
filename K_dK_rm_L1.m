% MPF learning for reliable moment model
% Modified from K_dK_ising_L1.m
% See Sohl-Dickstein et al. 2011
% https://github.com/Sohl-Dickstein/Minimum-Probability-Flow-Learning

function [K, dK] = K_dK_rm_L1( h, words, X )
    h=h';
    w = length(h);
    lambda = 5.9948e-04;
    [N, T] = size( X );
    
    % Single bit flip
    activewords = ( words'*X == sum(words,1)'*ones(1,T) );
    Etemp = h*activewords;
    E = ones(N,1)*Etemp;
    
    Eflip = zeros(N,T);
    dK_single = zeros(1,length(h));
    Kfull = zeros(N,T);
    
    word_orders = sum(words,1)';
    for k = 1:N
        Xflip = X;
        Xflip(k,:) = 1 - Xflip(k,:);
        activewords_flip = bsxfun(@eq,word_orders,words'*Xflip);
        Eflip(k,:) = h*activewords_flip;
        Kfull(k,:) = exp(.5*(Etemp - Eflip(k,:)));
        dK_single = dK_single + Kfull(k,:)*(activewords - activewords_flip)';
    end
    K_single = sum(Kfull(:));
    
    % Full bit flip
    notX = 1-X;
    
    E = E(1,:);
    activewords_all = ( words'*notX == sum(words,1)'*ones(1,T) );
    Eflip = h*activewords_all;
    
    Kfull = exp(.5*(E-Eflip));
    K_all = sum(Kfull(:));
    
    dK_all = Kfull*(activewords - activewords_all)';
    
    K = (K_single + K_all) /T;
    dK = (dK_single + dK_all) /T;
    
    wsum = sum(abs(h(:))); K = K + lambda*wsum;
    dsum = sign(h);     dK = dK + lambda*dsum; dK= dK';
        
