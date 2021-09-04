% This function is used to extract the matrices A,B,C in the state-space
% representation.
% 
% Weighted matrix has been considered in the B and C matrices.
%
% The coordinate system is same as that in oss1.m

% The inputs are same as that in oss1.m (except the type of operator is excluded)


function [A,B,C,dUdy] = oss_operator(N, kx, kz, Re, F)
%% ----- derivative matrices ----- 
[D,y] = cheb(N);
D1 = D(2:N,2:N);
D2 = D^2;
D2 = D2(2:N,2:N);
S = diag([0; 1 ./(1-y(2:N).^2); 0]);
D4 = (diag(1-y.^2)*D^4 - 8*diag(y)*D^3 - 12*D^2)*S;
D4 = D4(2:N,2:N);

%% ----- weight matrix -----
[~,W] = clenCurt(N);
W = sqrtm(diag(W(2:N)));
W = [W, zeros(N-1), zeros(N-1);zeros(N-1), W, zeros(N-1);zeros(N-1), zeros(N-1), W];

%% ----- mean flow ----- 
if F == 1
    U = diag(1 - y(2:N) .* y(2:N));                       % Poiseuille flow
    dUdy = - diag(2 * y(2:N));
    dU2dy2 = -diag(2*ones(length(y(2:N)),1));
else
    U = diag(y(2:N));                            % Couette flow
    dUdy = eye(length(y(2:N)));
    dU2dy2 = zeros(length(y(2:N)));
end

%% ----- form operators -----
k2 = (kx * kx + kz * kz);
I = eye(N-1);

% inverting Laplacian requires BC
L = D2 - k2 * I;
L2 = D4 + k2 * k2 * I - 2 * k2 * D2;

% form OS, Squire and coupling operators
% and assemble into OSS operator

A1 = -1i * kx * U * L + 1i * kx * dU2dy2 + L2 ./ Re;     % Oss Sommerfeld operator  
Ao = L \ A1;

As = L ./ Re -1i * kx * U;                               % Squire operator
As1 = -1i * kz * dUdy;

A = [Ao zeros(N-1);As1 As];

B = [L \ (-1i * kx * D1), L \ (-k2 * I), L \ (-1i * kz * D1); 1i * kz * I, zeros(N-1), -1i * kx * I];

C = [1i * kx * D1, -1i * kz * I; k2 * I, zeros(N-1); 1i * kz * D1, 1i * kx * I] ./ k2;

%% add weighted matrix
B = B / W;
C = W * C;