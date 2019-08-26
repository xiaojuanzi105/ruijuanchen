function Ex_X=Ex_3s3o_Heun(a, h, N, z0, d, W, H, L)

Ex_y=zeros(d,N+1);

%% 4级4阶的显示RK方法
s=3;
A=zeros(s);
A(2, 1)=1/3;
A(3,2)=2/3;
b=[1/4 0 3/4];
c=[0;1/3;2/3];
%%
t0=a;
% y0=z0;
%y=zeros(2*d+1,N+1);
y0=z0;
t=zeros(s);
Ex_y( :, 1) = z0(1:d);
for i=2:N+1
    for j=1:s
        t(j)=t0+c(j)*h;
    end
    Y1=y0;
    Y2=y0+h*kron(A(2,1), eye(2*d+1))...
        *Obj(d, t(1), Y1, W, H, L);
    Y3=y0+h*kron(A(3,1:2), eye(2*d+1))...
        *[Obj(d, t(1), Y1, W, H, L); ...
        Obj(d, t(2), Y2, W, H, L)];
    
    y1=y0+h*kron(b(1:3), eye(2*d+1))...
        *[Obj(d, t(1), Y1, W, H, L); ...
        Obj(d, t(2), Y2, W, H, L);...
        Obj(d, t(3), Y3, W, H, L)];
       
    t0=t0+h;
    y0=y1;
   
    Ex_y(:,i)=y0(1:d);

end
Ex_X=Ex_y;
end