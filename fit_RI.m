%Reliable interaction fitting code
%Implements Ganmor et al. 2011 PNAS Algorithm
%To find the RI model describing some binary data
%Joel Zylberberg, joelzy@uw.edu 2013
%output is partition function Z, reliable_interactions
%and reliable_words
%
%function [reliable_words reliable_interactions Z] = fit_RI(spiketrains,RI_thresh)


%%%%%%%%%%%%%%%%%%
%spiketrains are numbins X numcells in size
%%%%%%%%%%%%%%%%%%%%%

function [reliable_words reliable_interactions Z] = fit_RI(spiketrains,RI_thresh)

%%%%%%%%%%%%%%%% initialize the variables and count words
%get sizes
spiketrains = spiketrains';
[numbins d] = size(spiketrains);
nneur = d; %size of words to consider

%initialize some data
word = reshape(spiketrains(1,1:nneur),nneur,1);
words(:,1) = word;
wordcounts(1) = 1;
numwords = 1;

%count the words in the data
for i=2:numbins

    %pull a word out of the data and compare to my list
    word = reshape(spiketrains(i,1:nneur),nneur,1);
    diffs = sum(abs(words-repmat(word,1,numwords)));
    whichword = find(diffs==0);

    %if it's not on the list, add it. Otherwise, update counter
    if length(whichword) ==0;
        words(:,numwords+1) = word;
        wordcounts(numwords+1) = 1;
        numwords = numwords+1;

    else
        wordcounts(whichword) = wordcounts(whichword) + 1;
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%
listofstates = words';
numspikes = transpose(sum(listofstates'));




%now sort in order of word frequency
[numspikes word_IDs] = sort(numspikes);
counts = wordcounts(word_IDs);
wordprobs = counts/numbins;
listofstates = listofstates(word_IDs,:);


%%%%%%%%%%%%%%%%%%%%%%%fit the model
%now that we have a list of word frequencies, let's fit the model...

%how many reliable patterns to consider...
whichtoconsider = (wordprobs>RI_thresh);
numtoconsider = sum(whichtoconsider)-1;

%%%%%%

interactions = zeros(numwords,1);

reliable_words = zeros(numtoconsider,nneur);
reliable_interactions = zeros(numtoconsider,1);

%first get the partition function
Z = numbins/counts(1); %1/P(0000...). 

wordcounter = 1;
for i = 2:numwords

    if whichtoconsider(i)
        %and what is the bitstring of that word?
        theword = listofstates(i,:);
        reliable_words(wordcounter,:) = theword;

        %first, put in the "naive" interaction, which doesn't subtract out the
        %part already accounted for by lower-order interactions
        interactions(i) = log(Z) + log(wordprobs(i));

        %then subtract out the stuff due to lower order interactions     
        if numspikes(i)>1
            %"active bits" are the sub-words whose interactions we need to subtract
            active_bits = find(theword*listofstates' == numspikes');

            for jj = 1:length(active_bits)
             %skip the "full word", and the "empty word" since we already
             %accounted for those
                if (active_bits(jj) ~= i) && (active_bits(jj) ~= 1)
                     interactions(i) = interactions(i) - interactions((active_bits(jj)));        
                end

            end
        end

        reliable_interactions(wordcounter) = interactions(i);
        wordcounter = wordcounter + 1;
    end
end
     
%%%%%%%%%%%%%%%%%%%
    
reliable_words = reliable_words';
 