function Im1_X = Im_1s1o_trape(a, h, N, z0, d, W, H, L)
% % 
%%  1-stage gauss
s=1;  A=1/2; B=1; C=1/2;err1=1;
% s=1;  A=0.001; B=1; C=0.1;err1=1;
%cancer data set
% %for 1-stage implicit method, the best parameters 
%with the step size h=0.01; 
%%  3-order Heun for initial condition
A1=[0 0 0; 1/3 0 0; 0 2/3 0]; B1=[1/4 0 3/4]; C1=[0;1/3;2/3];
%%
t0=a;
Zz=zeros(2*d+1,N+1);
XX=zeros(d,N+1);
Zz(:,1)=z0;
XX(:,1)=z0(1:d);
err1=1;
err2=1;
for i=2:N+1
t=kron(ones(s,1),t0)+h*C;
%% Calculate the first RK variable
tc11=t0; 
tc12=t0+C*h*C1(2); 
tc13=t0+C*h*C1(3);
Z11=z0;
Z12= z0+(h*C)*(kron(A1(2,1),eye(2*d+1))*Obj(d, tc11, Z11, W, H, L));
Z13= z0+h*C*kron(A1(3,1:2),eye(2*d+1))...
    *[Obj(d, tc11, Z11, W, H, L); Obj(d, tc12, Z12, W, H, L)];
Z10= z0+h*C*kron(B1(1:3),eye(2*d+1))...
    *[Obj(d, tc11, Z11, W, H, L); ...
    Obj(d, tc12, Z12, W, H, L); Obj(d, tc13, Z13, W, H, L)];
%% 
%According to the above two calculations, 
%the initial value conditions of the whole calculation process are obtained
% for j=1:N
Z0=Z10;%The dimension: (2*d+1)*(s=1)
while err1>=10^(-12) & err2>=10^(-12) 

F=Obj(d, t(1), Z0, W, H, L);
G0=Z0-kron(ones(s,1),z0)-h*kron(A,eye(2*d+1))*F;
dF = Grad_Obj(d, t, Z0, W, H, L);
% f = Grad_Obj(d, t, ~, A, ~)
Z1= Z0 - inv(eye(2*d+1) - h *kron(A,eye(2*d+1))*dF)*G0;
err1=norm(Z1-Z0);
err2=norm(G0);
Z0=Z1;
end
% while err1>0.001 
%     %Start the Newton method, F here,
% % is the non-linear function in the Newton method, 
% % and the dimension is s*d
% %%
% p = -Hk*Gk;
% alpha = 1/(2^j);
% Zk1=Zk+ h*alpha *p;
%  Fk1=Obj(d, t(1), Zk1, W, H, L);
% Gk1=Zk1-kron(ones(s,1),z0)-h*kron(A,eye(2*d+1))*Fk1;
% Sk=Zk1-Zk;
% Yk=Gk1-Gk;
% %%
%   gammak = 1/(Yk'*Sk);
%     skykT = Sk * Yk';
%     skskT = Sk * Sk';
%     E = eye(2*d+1);
%     Hk = (E-gammak*skykT)*Hk*(E-gammak*skykT) + gammak*skskT;
% %%
% k = k + 1;
%     Zk = Zk1;
%     Gk = Gk1;
%  err1=sqrt(Sk'*Sk);
% j=j+1;
% end

Z1=Z0(1:2*d+1);
z1=z0+h*kron(B, eye(2*d+1))...
    *Obj(d, t(1), Z1, W, H, L);
Zz(:,i)=z1;
t0=t0+h;
z0=z1;
XX(:,i)=z0(1:d);
end
Im1_X=XX;
end


