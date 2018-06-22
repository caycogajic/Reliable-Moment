% Approximates Z using Good Turing
% words_emp - words from data 
% h - interaction parameters
% words_fit_rm - words fit with rm model
% M - missing mass

function Z = get_Z_GT(words_emp,h,words_fit_rm,M)
    
    [~,w] = size(words_emp);
    
    P_rm = zeros(1,w);
    for k = 1:w
        theword = words_emp(:,k);
        activewords = ( words_fit_rm'*theword == sum(words_fit_rm,1)' );
        P_rm(k) = exp(-h'*activewords);
    end
        
    % Approximate partition function
    Z = sum(P_rm)/(1-M);

    