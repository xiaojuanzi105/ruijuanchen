clear all;
clc;
%  close all;
%% real data set
L=0*2*10^(-10);
N=500;

a=1;
%%
X=rand(10,10);
y=10*rand(10,1);
for i=1:10
    if y(i)<1
        y(i)=0;
    else
        y(i)=1;
    end
end
%%
W=X;
H=y; 
d=size(W,2);
%% initial conditions
z11=rand(d,1); z22=zeros(d,1); z33=1;
z0=[z11;z22;z33];
%%
SW1=1;
SW2=2;
ImF21=zeros(N+1,1);
ImF22=zeros(N+1,1);
ImF23=zeros(N+1,1);

ExF_NAG=zeros(N+1,1);
ExF_GD=zeros(N+1,1);
%% 
% h1=0.02;
% h2=0.1;
% h3=0.2;
% h_NAG=0.02;
% h_GD=0.02;
%%
h1=0.2;
h2=0.02;
h3=0.001;
h_NAG=0.01;
h_GD=0.01;
p1=2;
p2=4;
p3=5;

ImX21=gauss_bfgs_crj_diff(a, h2, N, z0, d, W, H, L, p1);
ImX22=gauss_bfgs_crj_diff(a, h2, N, z0, d, W, H, L, p2);
ImX23=gauss_bfgs_crj_diff(a, h2, N, z0, d, W, H, L, p3);

X_NAG=nesterov(h_NAG, N, z0, d, W, H, L, SW2);
X_GD=nesterov(h_GD, N, z0, d, W, H, L, SW1);
for j=1:N+1
%% loss 
ImF21(j) = Fu(ImX21(:, j), W, H, L);
ImF22(j) = Fu(ImX22(:, j), W, H, L);
ImF23(j) = Fu(ImX23(:, j), W, H, L);
% ExFF(j) = Fu(Exp_X(:, j), W, H, L);
ExF_NAG(j) = Fu(X_NAG(:, j), W, H, L);
ExF_GD(j) = Fu(X_GD(:, j), W, H, L);
end
% 
figure
semilogy(1:N+1, ImF21,'r-','LineWidth', 1.5);hold on
semilogy(1:N+1, ImF22,'b-','LineWidth', 1.5);hold on
semilogy(1:N+1, ImF23,'r--','LineWidth', 1.5);hold on
% semilogy(1:N+1, ExFF,'b--','LineWidth', 1.5);
semilogy(1:N+1, ExF_NAG,'k:','LineWidth', 1.5);hold on;
semilogy(1:N+1, ExF_GD,'k-.','LineWidth', 1.5);
xlabel('Iterations', 'FontSize',16);
ylabel('Objective','FontSize',16);
%title('Minimizing regularized quadratic function on  set','FontSize',16);
legend({'Im p=2','Im p=3','Im p=5','NAG','GD'},'FontSize',16); %'Ex-2s2o-Kutta'
set(gca,'FontSize',16);