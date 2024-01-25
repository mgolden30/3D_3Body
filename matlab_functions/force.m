function [f, df] = force(x)
  %   f - force
  %  df - Jacobian

  r1 = x(1:3);
  r2 = x(4:6);

  r12 = r1-r2;
  r13 = 2*r1 + r2;
  r23 = 2*r2 + r1;

  d12 = norm(r12);
  d13 = norm(r13);
  d23 = norm(r23);

  f  = [-r12/d12^3 - r13/d13^3;
        +r12/d12^3 - r23/d23^3;];

  % Jacobian of force
  df = zeros(6,6);
  I  = eye(3);
  dg = @(r,d) -I/d^3 + 3*r*r'/d^5; %macro for subexpression
  
  df(1:3, 1:3) =  dg(r12,d12) + 2*dg(r13,d13);
  df(1:3, 4:6) = -dg(r12,d12) +   dg(r13,d13);
  
  df(4:6, 1:3) = -dg(r12,d12) +   dg(r23,d23);
  df(4:6, 4:6) =  dg(r12,d12) + 2*dg(r23,d23);
end