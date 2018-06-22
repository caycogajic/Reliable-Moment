% turns index for moments into words (/spiking patterns)

function words = ind_to_words(index,N)

num_words = 0;
for k = 1:size(index,2)
    num_words = num_words + size(index{k},1);
end

words = zeros(N,num_words); word_ix = 0;
for k = 1:size(index,2)
    for j = 1:size(index{k},1)
        word_ix = word_ix+1;
        words(index{k}(j,:),word_ix) = 1;
    end
end