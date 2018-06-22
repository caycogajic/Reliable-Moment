% Return empirical probabilities of words 
% X - data
% P - probability of all words in the data
% words - all words in the data
% P_ones - probability of words only observed once

function [P,words,P_ones] = get_empirical_probs(X)

    [~,T] = size(X);
    words = unique(X','rows')';
    [~,w] = size(words);
    
    P = zeros(1,w);
    for k = 1:w
        theword = words(:,k);
        temp = ismember(X',theword','rows')';
        P(k) = sum(temp);
        X = X(:,~temp);
    end
    
    P_ones = sum(P==1)/T;
    
    
    P = P/T;
    
    % Sort
    [~,ix] = sort(P,'descend');
    P = P(ix); words = words(:,ix);
    