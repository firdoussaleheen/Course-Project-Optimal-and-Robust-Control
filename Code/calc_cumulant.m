%%% RETURNS z = ARRAY OF CUMULANTS
%%%         kap = ARRAY OF CUMULANTS / i! FOR CONTINUED FRACTIONS APPROXIMATION

function [D_vect z cumulants] = calc_cumulant(G, W, X, k, x0, n)
    
    D_vect = zeros(k,1);                % vector of values of D (*1.32)
    [m nn] = size(X);                   % m = # rows of X (number of time values computed)
    cumulants = zeros(k,1);             % vector of cumulants, ith space is ith cumulant
    z = zeros(k, 1);                    % vector of ki/L!, where ki is the ith cumulant
    
    for i = 1:k
        
        Htrace = zeros(m,1);            % initialize Htrace vector
        H = zeros(n);                   % initialize H matrix
        
        for j = 1:m
            H = reshape(X(j,(i-1)*n+1:i*n),n,n);        % reshape H to n x n matrix
            Htrace(j) = trace(H*G*W*G');  % Htrace(j) = tr[H(j)GWG']
        end
            
        Hint = trapz(Htrace);           % integration approximation of (*1.32)
        D_vect(i) = Hint;
        
        cumulants(i) = factorial(i)*2^(i-1)*x0'*H*x0 + factorial((i-1))*2^(i-1)*D_vect(i);
        z(i) = cumulants(i)/factorial(i);
    end
end