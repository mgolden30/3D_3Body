

function [f, df] = objective(z, N)
  %{
  PURPOSE:
  Compute the objective function of periodic orbits
  %}

  [xf, dxf] = symplectic_steps( z(1:12), z(13)/N, N );
  
  f = xf - z(1:12); %objective function

  df  = zeros(12,13);   %Jacobian
  
  df(:,1:12) = dxf - eye(12);
  v = [xf(7:12); force(xf)];
  df(:,13) = v; %approximate the derivative with respect to period as the state space velocity
end
