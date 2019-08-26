function F=Fu(x, A, b, L)
% L=0.1;
% x=z(1:d);
F=0.5*x'*(A'*A)*x - b'*A*x + 0.5*(b'*b) + 0.5*L*(x'*x);
end