% This function is used to plot the input-output resolvent norm of the 
% Plane Poiseuille/Couette flow. 

% The state-space representation is shown below:
%       dx/dt = Ax + Bu
%           y = Cx
%       state x -- [v, \eta]^T (wall normal perturbation velocity & wall normal vorticity perturbation)
%       input u -- [dx, dy, dz]^T (non-linear terms & external forces)
%       output y -- [u, v, w]^T (perturbation velocities in three directions)

% In order to fit the energy norm, the input matrix and the output
% matrix should be changed using Clenshaw-Curtis quadrature. Then the altered input 
% matrix is B = B / W, the altered output matrix is C = W * C, where W can be expressed
% as W^{H} * W = M, M is the weighted matrix attained from clenCurt.m.

% In order to achieve accurate results, the number of modes (N) should be
% around 100.

%% specify the parameters
N = 100;
kx = 1;
kz = 0.25;
Re = 1000;
F = 0; % (1 for Poiseuille flow, 0 for Couette flow)

%% calculate the resolvent norm and the resolvent limit
% The resolvent limit only considers the eigenvalues of the matrix A.
% Because of the non-normality of the A matrix, the resolvent norm is
% greater than the resolvent limit.
[A,B,C] = oss_operator(N, kx, kz, Re, F);
EA = diag(eig(A));                  % eigenvalue decomposition of the matrix 'A'

omega = [-1.5:0.01:1.5];
ff = length(omega);
sig = zeros(ff,1);
sie = zeros(ff,1);
w = omega(1);
dw = omega(2) - omega(1);

for j = 1:ff
    D1 = -i * w * eye(2*N-2) - A;
    D2 = -i * w * eye(2*N-2) - EA;      % only consider the spectrum of matrix A
    R1 = C / D1 * B;
    R2 = C / D2 * B;
    sig(j) = norm(R1);
    sie(j) = norm(R2);
    w = w + dw;
end 

%% plot the result
semilogy(omega,sig,'Color','k','LineWidth',2)
hold on
semilogy(omega,sie,'Color','r','LineWidth',1)

axis([-1.5 1.5 1 100])
grid on
xlabel('$\omega$','interpreter','latex')
ylabel('$||(C(i\omega I - A))^{-1}B||$','interpreter','latex')
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
set(gca,'Fontsize',16)
