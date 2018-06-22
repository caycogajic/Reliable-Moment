% Gibbs sampling for Reliable Moment model
% Modified from sample_ising.m
% h - interaction parameters
% rm_words - words(ie spiking patterns) corresponding to reliable moments 

function X_out = sample_rm( h, rm_words, ndims, nsamples, burnin, independent_steps )
    nsamplingsteps = burnin + floor((nsamples-1)*(independent_steps));
	sigmoid = @(x) 1./(1+exp(-x));
    
    % choose dimensions to gibbs sample
    upd_i = floor( rand( nsamplingsteps, 1 ) * ndims ) + 1;
    % precalculate the random nubers for comparison
    uni_rand = rand( nsamplingsteps, 1 );

    X_out = zeros( ndims, nsamples );
    x = floor( rand( ndims, 1 ) * 2 );

    i_out = 1;
    next_sample = burnin;

    for si = 1:nsamplingsteps
        x(upd_i(si)) = 0;
        x(upd_i(si)) = 1; 
        whichwords = (rm_words(upd_i(si),:)==1) & (x'*rm_words == sum(rm_words,1));
        E_act = sum(h(whichwords));
        
        p_act = sigmoid(-E_act); % ****** sampling 0/1 here
        if p_act > uni_rand(si)
            x(upd_i(si)) = 1;
        else
            x(upd_i(si)) = 0;
        end	   
        
        if si == next_sample % copy to the output array if appropriate
            next_sample = si + independent_steps;
            X_out(:,i_out) = x;
            i_out = i_out + 1;
            if mod(i_out,100) == 0
                i_out
            end
        end
    end    