% Brute force calculation of partition function for pw model
function Z_pw = get_Z_pw(J)
    N= size(J,1);
    
    Z_pw = 0;
    for i = 1:2^N
        x = zeros(N,1);
        x_ = de2bi(i-1);
        x(1:length(x_)) = x_;
        Z_pw = Z_pw+exp(-x'*J*x);
    end
    