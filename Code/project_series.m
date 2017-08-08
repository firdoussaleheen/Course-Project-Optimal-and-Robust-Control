%%% All * equation numbers refer to Statistical Control Book
%%% All ** equation numbers refer to Density Shaping Paper

clear
clc

global n %% System matrix is n x n
global k %% k = number of cumulants
global L G

syms z s t

n = 1;
k = 2;
gamma = zeros(1,k);
gamma(1)=1;
gamma(2)=1;
gamma(3)=1;
gamma(4)=1;
gamma(5)=1;
gamma(6)=1;
gamma(7)=1;
gamma(8)=1;
gamma(9)=1;
gamma(10)=1;

%%% Initialize System matrices
A = 1;              %% n x n
B = 1;              %% n x m
Q = 1;              %% n x n
R = 1;              %% m x m
W = 0.2;
%W = normrnd(1,1);   %% n x n
G = 0.5;            %% n x n
x0 = 1;             %% n x 1 initial conditions
L = k;              %% counter for continued fractions

%%% WON'S SOLUTION
%%% couples Riccati equations for MCV case (k = 2)
%syms M V
%M1 = A'*M + M*A + Q - M*B*inv(R)*B'*M + gamma(1)^2*V*B*inv(R)*B'*V;
%M2 = 4*M*G*W*G'*M + A'*V + V*A - M*B*inv(R)*B'*V - V*B*inv(R)*B'*M - 2*gamma(1)*V*B*inv(R)*B'*V;

%%% Solution of couple MCV Riccati equations
%S=solve(M1,M2);
%M=S.M;
%V=S.V;

%%% MY SOLUTION
%%% Solve general system of k coupled differential Riccati equations

tspan = [0 5];           % integration time span
X_final = zeros(k*n,n);  % terminating conditions
[T X] = ode45(@mRiccati, tspan, X_final, [], A, B, Q, R, G, W, n, gamma, k);
                        % X = solution of differentials over tspan
                        % X is a integrated backwards in time
                        % thus first row of X is terminal values
                        % last row of X is initial values
X = flipud(X);          % flips X
                        
% returns kap - a cell of cumulants
%         z - a cell of ki/i!, for continued fractions approximation (**22)
%%%%not needed???%%%%%%%%%%%[c z kap]= kcc(L+1, x0);
[D_vect z cumulants] = calc_cumulant(G, W, X, k, x0, n);

% computes A from continued fractions approximation
contA = Amat(z, L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create first characteristic function from second (cumulant series)
syms s
%D_vect = flipud(D_vect);
firstchar = 1;
for l = 1:k
    firstchar = firstchar*exp(s^(l)*z(l));
end
%series = tf(D_new,1);
%f_cs = exp(series);
f_cs = simplify(eval(firstchar));
syms x
f_cs = subs(f_cs, s, -1i*x);

% inverse fourier transform of first char function to get PDF (cumulant series)
%costfn_cs = ifourier(f_cs, x, t);
%costfn_cs = ifft(f_cs, x, t);

% PLOT PROBABILITY DENSITY FUNCTION (cumulant series)
t = (z(1)-50):0.01:(z(1)+50);
ezplot(costfn_cs, t)
title('Density function of the cost (Cumulant Series Method), J')
xlabel('J')
ylabel('density function')
axis auto
% ylim([0 .5])
% xlim([-z(1)-20 -z(1)+20])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%y = ifft(costfn_cs, t)
syms t
pdf = eval(int(costfn_cs, t, -inf, inf));
pdf
%syms u
%f = 1/(2*pi) * int(f_cs*exp(1i*x*u),x,-inf,inf);
