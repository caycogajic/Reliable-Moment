% get frequencies of ri model for words in words_emp

function f_ri = get_ri_probs(words_emp,h,words_fit_ri,Z)

    w = size(words_emp,2);
    f_ri = zeros(1,size(words_emp,2));

    for i=1:w

        %get the data word
        theword = words_emp(:,i)';

        %"active bits" are the sub-words whose interactions we need to subtract
        active_bits = find(theword*words_fit_ri == sum(words_fit_ri));
        interact_sum = 0;
        for jj = 1:length(active_bits)
            interact_sum = interact_sum + h(active_bits(jj));
        end

        f_ri(i) = (exp(interact_sum))/Z;

    end