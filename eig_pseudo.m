% This code is used to generate Fig.6 in 'analysis of fluid systems
% stability receptivity sensitivity'.
% There are three main parts in this code.
%   1. eigenvalues of the O-S operator 
%   2. pseudospectra of the O-S operator
%   3. numerical range (function 'eigenvalue_selection.m' will be used)

%% specify the parameters
N = 100;
kx = 1;
kz = 0.25;
Re = 1000;
F = 0;


%% eigenvalue calculation
% The eigenvalue problem is of the form: Av = (\lambda)v = (-i*\omega)v,
% the eigenvalues are calculated first then comes the frequency.
% BE CAREFUL! To get the accurate eigenvalues, the number of modes should
% be greater than 50.
[A,B,C] = oss_operator(N,kx,kz,Re,F);
w = eig(A);                             % compute the eigenvalues of the operator
wr = -imag(w);
wi = real(w);
%%
figure(1)
x_left = -1.2; x_right = 1.2; y_low = -2; y_high = 0.2;
scatter(wr,wi,'filled'), hold on
yline(0,'Color','k','LineWidth',2)
patch('Faces',[1 2 3 4],'Vertices',[x_left 0;x_right 0;x_right y_high;x_left y_high],'EdgeColor','none','FaceColor','black','FaceAlpha',.2)
axis([x_left x_right y_low y_high])
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
xlabel('$\omega_r$','interpreter','latex')
ylabel('$\omega_i$','interpreter','latex')
set(gca,'Fontsize',16)
%% pseudospectra calculation
% The interval should be set small if accurate pseudospectra contours are
% required.
x = [-1.1: 0.02: 1.1];                  % specify the coordinate range of the plot
y = [-1.5: 0.02: 0.4];

ss = zeros(length(y),length(x));        % store the minimum singular value of each point

for m = 1: length(x)
    for n = 1: length(y)
        z = -1i * (x(m) + 1i*y(n));
        ss(n,m) = min(svd(z * eye(2*N-2) - A));         % Be careful of the positions of 'm' and 'n'
    end
end

%% numerical range calculation
% The algorithm is used. (section 3.3 in 'analysis of fluid systems stability receptivity sensitivity')

theta = [0:0.01:2*pi];
zz = zeros(length(theta),1);

% Since the numerical range includes all the eigenvalues, and some
% eigenvalue-eigenvector pairs are unecessary plus far away from the
% mainstream, so only a few eigenvalues which are near to the
% stable-unstable delimiting line are selected. So there is a minimum
% imagnary part value (IMIN).
IMIN = -2.5;        
nA = eigenvalue_selection(A,IMIN,C); % see function 'elgenvalue_selection.m'

% exceute the algorithm
for k = 1 : length(theta)
    R = exp(1i * theta(k));
    MN = R * nA;
    MNbar = (MN + MN')/2;
    [V,D] = eig(MNbar);
   
    [dmax,mm] = max(diag(D));
    vmax = V(:,mm);
    zz(k) = 1i * (vmax' * nA * vmax) / (vmax' * vmax);
end


%% plot the result
figure(2)
hold on
axis([x(1) x(end) y(1) y(end)])         % specify the coordinate range

% plot the eigenvalues
plot(wr,wi,'b.','MarkerSize',10)

% plot the numerical range
plot(real(zz),imag(zz),'Color','r','LineWidth',2)   
% plot(real(zz(1)),imag(zz(1)),'r*','MarkerSize',10)

% plot the pseudospectra 
v = [0.000001,0.00001,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.15,0.2]; % specify the values of the contours
contour(x,y,ss,v,'LineWidth',0.5,'ShowText','off','Color',[0.4 0.4 0.4])

% plot the line 'y=0' delimiting the stable and unstable areas
line([x(1) x(end)],[0 0],'LineStyle','-','Color','Black','LineWidth',1)

patch('Faces',[1 2 3 4],'Vertices',[x(1) 0;x(end) 0;x(end) y(end);x(1) y(end)],'EdgeColor','none','FaceColor','black','FaceAlpha',.2)

ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
xlabel('$\omega_r$','interpreter','latex')
ylabel('$\omega_i$','interpreter','latex')
set(gca,'Fontsize',16)