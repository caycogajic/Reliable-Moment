% MPF learning for pairwise models
% Borrowed from Sohl-Dickstein et al. 2011
% https://github.com/Sohl-Dickstein/Minimum-Probability-Flow-Learning

function [K, dK] = K_dK_ising_L1( J, X )
    lambda = 5.9948e-04*.8;
    [ndims, nbatch] = size( X );

    J = reshape( J, [ndims, ndims] );
    J = (J + J')/2;

    Y = J*X;
    diagJ = diag(J);
    % XnotX contains (X - [bit flipped X])
    XnotX = 2*X-1;
    % Kfull is a [ndims, nbatch] matrix containing the contribution to the objective function from flipping each bit in the rows, for each datapoint on the columns
    Kfull = exp(XnotX .* Y - (1/2)*diagJ(:,ones(1,nbatch))); 
    K = sum(Kfull(:));
    
    lt = Kfull.*XnotX;
    dJ = lt * X';
    dJ = dJ - (1/2)*diag( sum(Kfull, 2) );
             
    
    %% add all bit flip comparison case
    % calculate the energies for the data states
    EX = sum( X.*Y );
    % calculate the energies for the states where all bits are flipped relative to data states
    notX = 1-X;
    notY = J*notX;
    EnotX = sum( notX.*notY );
    % calculate the contribution to the MPF objective function from all-bit-flipped states
    K2full = exp( (EX - EnotX)/2 );
    K2 = sum(K2full);
    % calculate the gradient contribution from all-bit-flipped states
    dJ2 = bsxfun( @times, X, K2full ) * X'/2 - bsxfun( @times, notX, K2full ) * notX'/2;
    % add all-bit-flipped contributions on to full objective
    K = K + K2;
    dJ = dJ + dJ2;
    
    % symmetrize coupling matrix
    dJ = (dJ + dJ')/2;
    dK = dJ(:);

    % average over batch
    K  = K  / nbatch;
    dK = dK / nbatch;
    
    wsum = sum(abs(J(:)));
    K = K + lambda*wsum;
    
    dsum = sign(J(:));
    dK = dK + lambda*dsum;
    
    