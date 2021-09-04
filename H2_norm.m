% This code is used to generate the figure 5 in 'Componentwise energy amplication in channel flows'

% There are two variables, kx (streamwise wave number) and kz (spanwise wave number)

% The coordinates are logarithmically gridded in both directions.

% The neurograph is plotted using the command 'imagesc'
% The command 'imagesc' sets the reversed vertical coordinate by default.
% To alter this, go to menu
% edit--axes properties...--RULERS--YDir, choose 'normal'

% BE CAREFUL! This program is time-consuming!

%% set the ranges of wave numbers
kxe = [-4 : 0.01 : 0];      % the exponent of wave number
kze = [-2 : 0.01 : 1];

kx = 10 .^ kxe;             % the value of wave number
kz = 10 .^ kze;

lkx = length(kxe);
lkz = length(kze);

H = zeros(lkx,lkz);         % store the H_2 norm

for a = 1 : lkx
   for b = 1 : lkz
       [A,B,C] = oss_operator(30,kx(a),kz(b),2000,1);
       P = ss(A,B,C,0);
       H(a,b) = log10(norm(P,2));
   end
end

%% plot the result
imagesc(H);
xticks([1 101 201 301])
xticklabels([10^(-2) 10^(-1) 10^(0) 10^(1)])
yticks([1 101 201 301 401])
yticklabels([10^(-4) 10^(-3) 10^(-2) 10^(-1) 1])
xlabel('k_z'),ylabel('k_x'),title('log_{10}||H||_2')