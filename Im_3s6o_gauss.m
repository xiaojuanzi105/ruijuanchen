% clear all;
% close all;
 function Im_X=Im_3s6o_gauss(a, h, N, z0, d, W, H, L)
% L=0.1;
% function z=Gauss(a, b, y0, d, h)
% Gauss(0, 5, [2*pi;0], 2, 0.01)
%% 3-stage 6-order Gauss: implicit Runge-Kutta method
%%
% h=0.10;
% N=200;
% a=h;
% % b=5;
% % d= 10;

% % x0=ones(d,1);
% % y0=ones(d,1);
% z11=rand(d,1);
% z22=zeros(d,1);
% 
% z0=[z11;z22];
% W = randn(d);
% H = ones(d,1);

%% 3-stage gauss
s=3;
IA=[5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;...
     1/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;...
     5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
IB=[5/18, 4/9, 5/18];
IC=[1/2-sqrt(15)/10; 1/2; 1/2+sqrt(15)/10];
%% 3-order Heun for initial conditions
A1=[0 0 0; 1/3 0 0; 0 2/3 0];
B1=[1/4 0 3/4];
C1=[0;1/3;2/3];
%%
t0=a;

IZ=zeros(2*d+1,N+1);
IX=zeros(d,N+1);

IZ(:,1)=z0;
% Im_Err=zeros(1,N+1);
IF=zeros(N+1,1);

IX(:,1)=z0(1:d);
% IF(1)=Fu(d,IX(:,1), W, H);
% X_True=inv(W'*W+L*eye(d))*(W'*H);
% X_True=0;
% Im_Err(1)=norm(IX(:,1)-X_True);
% i=2;
for i=2:N+1
%     i=(n-a)/h+1;
err1=1;
% err2=1;
t=kron(ones(s,1),t0)+h*IC;
%% Calculate the first RK variable
tc11=t0; 
tc12=t0+IC(1)*h*C1(2); 
tc13=t0+IC(1)*h*C1(3);
Z11=z0;
Z12= z0+(h*IC(1))*(kron(A1(2,1),eye(2*d+1))*Obj(d, tc11, Z11, W, H, L));
Z13= z0+h*IC(1)*kron(A1(3,1:2),eye(2*d+1))...
    *[Obj(d, tc11, Z11, W, H, L); Obj(d, tc12, Z12, W, H, L)];
Z10= z0+h*IC(1)*kron(B1(1:3),eye(2*d+1))...
    *[Obj(d, tc11, Z11, W, H, L); ...
    Obj(d, tc12, Z12, W, H, L);...
    Obj(d, tc13, Z13, W, H, L)];

%% Calculate the SECOND RK variable
tc21=t0; 
tc22=t0+IC(2)*h*C1(2); 
tc23=t0+IC(2)*h*C1(3);
Z21= z0;
Z22= z0+h*IC(2)*kron(A1(2,1),eye(2*d+1))...
    *Obj(d, tc21, Z21, W, H, L);
Z23= z0+h*IC(2)*kron(A1(3,1:2),eye(2*d+1))...
    *[Obj(d, tc21, Z21, W, H, L); Obj(d, tc22, Z22, W, H, L)];
Z20= z0+h*IC(2)*kron(B1(1:3),eye(2*d+1))...
    *[Obj(d, tc21, Z21, W, H, L); Obj(d, tc22, Z22, W, H, L); Obj(d, tc23, Z23, W, H, L)];

%% 
%% Calculate the third RK variable
tc31=t0; 
tc32=t0+IC(3)*h*C1(2); 
tc33=t0+IC(3)*h*C1(3);
Z31= z0;
Z32= z0+h*IC(3)*kron(A1(2,1),eye(2*d+1))...
    *Obj(d, tc31, Z31, W, H, L);
Z33= z0+h*IC(3)*kron(A1(3,1:2),eye(2*d+1))...
    *[Obj(d, tc31, Z31, W, H, L); Obj(d, tc32, Z32, W, H, L)];
Z30= z0+h*IC(3)*kron(B1(1:3),eye(2*d+1))...
    *[Obj(d, tc31, Z31, W, H, L); Obj(d, tc32, Z32, W, H, L); Obj(d, tc33, Z33, W, H, L)];
%According to the above two calculations, 
%the initial value conditions of the whole calculation process are obtained
% for j=1:N
Z0=[Z10;Z20;Z30];%The dimension: 2*d*s
H0=eye(s*(2*d+1));
Hk=H0;
Zk=Z0;

F=[Obj(d, t(1), Z0(1:2*d), W, H, L);...
     Obj(d, t(2), Z0(2*d+1:4*d), W, H, L);...
     Obj(d, t(3), Z0(4*d+1:2*s*d), W, H, L)];
G0=Z0-kron(ones(s,1),z0)-h*kron(IA,eye(2*d+1))*F;
Gk=G0;
j=1;
while err1>0.001 
    %Start the Newton method, F here,
% is the non-linear function in the Newton method, 
% and the dimension is s*d
%%
%  F=[Obj(d, t(1), Z10, W, H);...
%      Obj(d, t(2), Z20, W, H)];
% G0=Z0-kron(ones(s,1),z0)-h*kron(A,eye(2*d))*F;
p = -Hk*Gk;
alpha = 1/(2^j);
Zk1=Zk+ alpha *p;
 Zk11 = Zk1(1:2*d+1);
 Zk12 = Zk1(2*d+2:4*d+2);
 Zk13 = Zk1(4*d+3:2*d*s+s);

 Fk1= [Obj(d, t(1), Zk11, W, H, L);...
          Obj(d, t(2), Zk12, W, H, L);...
          Obj(d, t(3), Zk13, W, H, L)];
Gk1=Zk1-kron(ones(s,1),z0)-h*kron(IA,eye(2*d+1))*Fk1;
Sk=Zk1-Zk;
Yk=Gk1-G0;
 err1=sqrt(Sk'*Sk);
j=j+1;
end
Z0=Zk1;
Z10=Z0(1:2*d+1);
Z20=Z0(2*d+1+1:4*d+2);
Z30=Z0(4*d+2+1:s*2*d+3);

z1=z0+h*kron(IB, eye(2*d+1))...
    *[Obj(d, t(1), Z10, W, H, L);...
      Obj(d, t(2), Z20, W, H, L);...
      Obj(d, t(3), Z30, W, H, L)];
IZ(:,i)=z1;

t0=t0+h;
z0=z1;
IX(:,i)=z0(1:d);
% IF(i)=Fu(d, IX(:,i), W, H, L);
% Im_Err(i)=norm(IX(:,i)-X_True);

end
Im_X=IX;
% figure;
% semilogy(1:N+1,Im_Err,'or');hold on;
% figure;
% figure;
% plot(1:N+1,IF);
% 
% E1=XX(:,end);
% E2=X_True;
