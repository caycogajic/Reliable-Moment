% determine with are moments are over threshold
% X - data
% p_min - threshold probability

function [E_fit,indices_fit] = calc_all_moments(X,p_min)

[N,~] = size(X);

E_temp = mean(X,2); indices_1 = find(E_temp>p_min);
E_fit = cell(1,1); E_fit{1,1} = E_temp(indices_1); 
indices_fit = cell(1,1); indices_fit{1,1} = indices_1; 
stop = 0; k=2; 			% k gives the order of interaction

while stop == 0
    indices_tocheck = nchoosek(indices_1,k); 
    if k>2
        indices_tocheck = checkifinlist(indices_tocheck,indices_old);
    end
    E_fit_temp = [];
    indices_old = [];
    for j = 1:size(indices_tocheck,1)
        % Moment / probability of subset of neurons spiking together
        E_temp = mean(prod(X(indices_tocheck(j,:),:),1));
        if E_temp >= p_min
            E_fit_temp = [E_fit_temp;E_temp];
            indices_old = [indices_old; indices_tocheck(j,:)];
        end
    end
    if isempty(indices_old) || k == N
        stop = 1;
    elseif ~isempty(indices_old)
        E_fit{1,k} = E_fit_temp;
        indices_fit{1,k} = indices_old;        
    end
    k = k+1;
end