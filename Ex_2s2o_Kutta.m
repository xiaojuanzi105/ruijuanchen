function Exp_X=Ex_2s2o_Kutta(a, h, N, z0, d, W, H, L)

% L=0;
% d=10;
% N=1000;
% h=0.1; a=h;
% 
% z11=ones(d,1);
% z22=zeros(d,1);
% % ExError = zeros(N+1,1);
% y0=[z11;z22];
% 
% W = randn(d);
% H = ones(d,1);
% d1=d/2;
% H(d1+1:d)=zeros(d1,1);

Ex_y=zeros(d,N+1);
% Err=zeros(1,N+1);
% FK=zeros(1,N);
% X_True= inv(W'*W+L*eye(d))*(W'*H);
% Err(1) = norm(z11-X_True);
%% 2级2阶的显示RK方法
s=2;
A=zeros(s);
A(2, 1)=1;
b=[1/2 1/2];
c=[0;1];
%%
t0=a;
y0=z0;
% N=(t2-t1)/h;
% % y=zeros(2*d+1, N+1);
% % y(:,1)=z0;
% FK=zeros(1,N+1);
t=zeros(s);
%   FK(1)=Fu(z0(1:d), W, H, L);

Ex_y( :, 1) = z0(1:d);
for i=2:N+1
    for j=1:s
        t(j)=t0+c(j)*h;
    end
    
    Y1=y0;
    Y2=y0+h*kron(A(2,1), eye(2*d+1))*Obj(d, t(1), Y1, W, H, L);%f(t(1),Y1)
    
    y1=y0+h*(kron(b(1:2), eye(2*d+1))...
        *[Obj(d, t(1), Y1, W, H, L); Obj(d, t(2), Y2, W, H, L)]);
    Ex_y(:,i)=y0(1:d);
%     Err(i) =norm( Ex_y(:,i) -X_True);
%     FK(i)=Fu(Ex_y(:,i), W, H, L);
    
    t0=t0+h;
    y0=y1;
end
Exp_X=Ex_y;
% semilogy

% subplot(2,1,1)
% plot(1:N+1,FK,'b-o')
% subplot(2,1,2)
% semilogy(1:N+1,Err,'b--')
end