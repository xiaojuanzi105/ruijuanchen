function Exp_X=Ex_1s1o_Kutta(a, h, N, z0, d, W, H, L)


%% 1级1阶的显示RK方法
s=1;
A=zeros(s);

b=1;
c=0;
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
        t(j)=t0+c*h;
    end
    
    Y1=y0;

    y1=y0+h*(kron(b, eye(2*d+1))...
        *Obj(d, t(1), Y1, W, H, L));
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