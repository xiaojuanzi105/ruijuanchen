function grad_F=grad_Fu(x, A, b, L)
  grad_F=(A'*A)*x - A' *b+ L*x;
end