function dXdt = mRiccati(t, X, A, B, Q, R, G, W, n, gamma, k)
    X = reshape(X, k*n, n); %Convert from "k*n^2"-by-1 to "k*n"-by-"n"

        % Kbar = Hr for r = 1 in (*4.33)
    K_bar = X(1:n,1:n);

        % loop to compute sum(gamma*Hr) in (*4.33)
    for i=2:k
        K_bar = K_bar + gamma(i)*X((i-1)*n+1:i*n,1:n);
    end

        % create K in (4.33)
    K = -inv(R)*B'*K_bar;

        % initialize dXdt
    dXdt = zeros(k*n,n);

        % construct dXdt k = 1 case of coupled differential Riccati equations
    dXdt(1:n,1:n) = (A+B*K)'*X(1:n,1:n) + X(1:n,1:n)*(A+B*K) + K'*R*K + Q;

        % construct dXdt which is matrix of k > 1 coupled differential Riccati
        % equations
    for i = 2:k
        for j = 1:i-1
            final = 2*factorial(i)/(factorial(j)*factorial(i-j))*X((j-1)*n+1:j*n,1:n)*G*W*G'*X(((i-j)-1)*n+1:(i-j)*n,1:n);
        end
        dXdt((i-1)*n+1:i*n,1:n) = (A + B*K)'*X((i-1)*n+1:i*n,1:n) + X((i-1)*n+1:i*n,1:n)*(A+B*K) + final;
    end

    dXdt = dXdt(:); % Convert from "k*n"-by-"n" to "k*n^2"-by-1
end


