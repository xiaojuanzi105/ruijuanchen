%% ?Direct Runge-Kutta Discretization Achieves Acceleration
function f= dF(x, A, b)
L=0.1;
f = (A'*A)*x - A'*b + L*x;
end
