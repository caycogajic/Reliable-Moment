% Return reliable moment probabilities for  spiking patterns in words_emp
% words_emp - words from data 
% h - interaction parameters
% words_fit_rm - words fit with rm model
% Z - partition function

function P_rm = get_rm_probs(words_emp,h,words_fit_rm,Z)
    
    [~,w] = size(words_emp);
    
    P_rm = zeros(1,w);
    for k = 1:w
        theword = words_emp(:,k);
        activewords = ( words_fit_rm'*theword == sum(words_fit_rm,1)' );
        P_rm(k) = exp(-h'*activewords);
    end
    
    P_rm = P_rm/Z;
    