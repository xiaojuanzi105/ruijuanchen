function Ex_X = gauss_bfgs_crj_diff(a, h, N, z0, d, W, H, L, cp)
% % 
% % Gauss(0, 5, [2*pi;0], 2, 0.01)

%% 2-stage gauss
s=2;
A=[1/4 1/4-sqrt(3)/6; 1/4+sqrt(3)/6 1/4];
B=[1/2 1/2];
C=[1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
%% 3-order Heun for initial condition
A1=[0 0 0; 1/3 0 0; 0 2/3 0];
B1=[1/4 0 3/4];
C1=[0;1/3;2/3];
%%
% L=100;
t0=a;

Zz=zeros(2*d+1,N+1);
XX=zeros(d,N+1);

Zz(:,1)=z0;
% Err=zeros(1,N+1);
FF=zeros(N+1,1);
XX(:,1)=z0(1:d);
% FF(1)=Fu(XX(:,1), W, H, L);
% Fu(x, A, b, L)
% Err(1)=norm(XX(:,1)-X_True);
% i=2;
for i=2:N+1
%     i=(n-a)/h+1;
err1=1;
% err2=1;
t=kron(ones(s,1),t0)+h*C;
%% Calculate the first RK variable
tc11=t0; 
tc12=t0+C(1)*h*C1(2); 
tc13=t0+C(1)*h*C1(3);
Z11=z0;
% Obj(d, t, y, A, b)
Z12= z0+(h*C(1))*(kron(A1(2,1),eye(2*d+1))*Obj(d, tc11, Z11, W, H, L, cp));
% f = Obj(d, t, z, A, b, L)
Z13= z0+h*C(1)*kron(A1(3,1:2),eye(2*d+1))...
    *[Obj(d, tc11, Z11, W, H, L, cp); Obj(d, tc12, Z12, W, H, L, cp)];

Z10= z0+h*C(1)*kron(B1(1:3),eye(2*d+1))...
    *[Obj(d, tc11, Z11, W, H, L, cp); Obj(d, tc12, Z12, W, H, L, cp); Obj(d, tc13, Z13, W, H, L, cp)];
%% Calculate the SECOND RK variable
tc21=t0; tc22=t0+C(2)*h*C1(2); tc23=t0+C(2)*h*C1(3);
Z21= z0;
Z22= z0+h*C(2)*kron(A1(2,1),eye(2*d+1))...
    *Obj(d, tc21, Z21, W, H, L, cp);
Z23= z0+h*C(2)*kron(A1(3,1:2),eye(2*d+1))...
    *[Obj(d, tc21, Z21, W, H, L, cp); Obj(d, tc22, Z22, W, H, L, cp)];
Z20= z0+h*C(2)*kron(B1(1:3),eye(2*d+1))...
    *[Obj(d, tc21, Z21, W, H, L, cp); Obj(d, tc22, Z22, W, H, L, cp); Obj(d, tc23, Z23, W, H, L, cp)];
%% 
%According to the above two calculations, 
%the initial value conditions of the whole calculation process are obtained
% for j=1:N
Z0=[Z10;Z20];%The dimension: (2*d+1)*s
H0=eye(2*(2*d+1));
Hk=H0;
Zk=Z0;

F=[Obj(d, t(1), Zk(1:2*d+1), W, H, L, cp);...
     Obj(d, t(2), Zk(2*d+2:4*d+2), W, H, L, cp)];
 
G0=Z0-kron(ones(s,1),z0)-h*kron(A,eye(2*d+1))*F;
Gk=G0;
j=1;k=1;
while err1>0.01 
    %Start the Newton method, F here,
% is the non-linear function in the Newton method, 
% and the dimension is s*d
%%
%  F=[Obj(d, t(1), Z10, W, H);...
%      Obj(d, t(2), Z20, W, H)];
% G0=Z0-kron(ones(s,1),z0)-h*kron(A,eye(2*d))*F;
p = -Hk*Gk;
alpha = 1/(2^j);
Zk1=Zk+ h*alpha *p;
 Zk11 = Zk1(1:2*d);
 Zk12 = Zk1(2*d+1:4*d);
 Fk1=[Obj(d, t(1), Zk11, W, H, L, cp);...
     Obj(d, t(2), Zk12, W, H, L, cp)];
Gk1=Zk1-kron(ones(s,1),z0)-h*kron(A,eye(2*d+1))*Fk1;
Sk=Zk1-Zk;
Yk=Gk1-Gk;
%%
  gammak = 1/(Yk'*Sk);
    skykT = Sk * Yk';
    skskT = Sk * Sk';
    E = eye(2*(2*d+1));
    Hk = (E-gammak*skykT)*Hk*(E-gammak*skykT) + gammak*skskT;
%%
   
k = k + 1;
    Zk = Zk1;
    Gk = Gk1;
 err1=sqrt(Sk'*Sk);
j=j+1;
end
Z0=Zk1;
Z10=Z0(1:2*d+1);
Z20=Z0(2*d+2:4*d+2);
z1=z0+h*kron(B, eye(2*d+1))...
    *[Obj(d, t(1), Z10, W, H, L, cp);...
    Obj(d, t(2), Z20, W, H, L, cp)];
Zz(:,i)=z1;
t0=t0+h;
z0=z1;
XX(:,i)=z0(1:d);
% FF(i)=Fu(XX(:,i), W, H, L);
% Err(i)=norm(XX(:,i)-X_True);
end

Ex_X=XX;
end
% end
% figure;
% semilogy(1:N+1,Err,'+b');hold on;
% figure;
% plot(1:N+1,FF);

% E1=XX(:,end);
% E2=X_True;
%%
% function Hk = BFGS(d, Hk,sk,yk)
%     gammak = 1/(yk'*sk);
%     skykT = sk * yk';
%     skskT = sk * sk';
%     E = eye(2*2*d);
%     Hk = (E-gammak*skykT)*Hk*(E-gammak*skykT) + gammak*skskT;
% end

%% %   Z(:,k)=Zk;
% r=Z0-kron(ones(s,1),z0)-h*kron(A,eye(2*d))*F;
%Grad_Obj(d, t, y, A, b)
% df=[Grad_Obj(d, t(1), Z10, W, H) zeros(2*d,2*d); ...
% zeros(2*d,2*d) 
% Z1=Z0-pinv(eye(s*2*d)-h*kron(A,eye(2*d))*df)*r;
% err1=norm(Z1-Z0);
% err2=norm(r);