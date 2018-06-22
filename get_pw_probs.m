function P_pw = get_pw_probs(words_emp,J,Z)
    
    [~,w] = size(words_emp);
    
    P_pw = zeros(1,w);
    for k = 1:w
        theword = words_emp(:,k);
        P_pw(k) = exp(-theword'*J*theword);
    end
    
    P_pw = P_pw/Z;