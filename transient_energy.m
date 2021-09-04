% This code is used to generate plots of transient energy growth.

% G(t) = ||C * exp(At) * B||^2

%% specify the values
N = 100;
kx = 1;
kz = 0.25;
Re = 2000;
F = 1;

[A,B,C] = oss_operator(N,kx,kz,Re,F);
t = [0:1:80];
G = zeros(length(t),1);

[~,~,V1] = svd(C * expm(3.*A) * B);
[~,~,V2] = svd(C * expm(15.*A) * B);
[~,~,V3] = svd(C * expm(45.*A) * B);

G1 = zeros(length(t),1);
G2 = zeros(length(t),1);
G3 = zeros(length(t),1);

for k = 1:length(t)
    G(k) = norm(C * expm(t(k).*A) * B)^2;
    G1(k) = norm(C * expm(t(k).*A) * B * V1(:,1))^2;
    G2(k) = norm(C * expm(t(k).*A) * B * V2(:,1))^2;
    G3(k) = norm(C * expm(t(k).*A) * B * V3(:,1))^2;
end


%% plot the result
G(1) = 1;
semilogy(t,G,'Color','k','LineWidth',2),hold on
semilogy(t,G1,'LineWidth',1)
semilogy(t,G2,'LineWidth',1)
semilogy(t,G3,'LineWidth',1)
axis([0 50 1 20])
grid on 
xlabel('t')
ylabel('G')

title('Poiseuille flow, \alpha = 1, \beta = 0.25, Re = 2000')
