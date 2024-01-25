function [xp, dxp] = symplectic_steps(x, dt, N)
  % Compute the forward time evolution xp of an initial state x
  % And its exact Jacobian dxp
  
  dxp = eye(12); %initialize Jacobian
  xp  = x;       %initialize position

  %The Jacobian for updating position is a constant matrix
  %initialize it once
  A = eye(12);

  % A = [I, dt/2*I; 0, I];
  
  %The Jacobian for momnetum updating will be upper/lower triangular, but
  %the identity part of it can be initialized a single time
  B = eye(12);

  c = zeros(3,1);
  d = zeros(3,1);

  %{
  %third order method
  c(1) = (sqrt(7/3) - 1)/2;
  c(2) = -c(1);
  c(3) = 1;

  d(1) = (1 + c(1))/(1 - c(1))/2;
  d(2) = -c(1)/2;
  d(3) = 1-d(1)-d(2);
  %}

  c = zeros(4,1);
  d = zeros(4,1);
  crt = 2.0^(1.0/3.0); %cube root 2

  c(1) = 1.0/ 2.0 /(2.0 - crt ); %c1
  c(2) = (1.0 - crt) / 2.0 / (2 - crt); %c2
  c(3) = c(2);
  c(4) = c(1);

  d(1) = 1.0/(2.0 - crt );
  d(2) = -crt/(2.0-crt);
  d(3) = d(1);
  d(4) = 0.0;



  for i = 1:N
    for s = 1:numel(c)
      %update position and Jacobian
      xp(1:6) = xp(1:6) + xp(7:12)*dt*c(s);
      A(1:6,7:12) = eye(6)*dt*c(s);
      dxp     = A*dxp;

      %update momentum and Jacobian
      [f, df]     = force(xp); % Compute force and its derivatives
      B(7:12,1:6) = dt*df*d(s);
      xp(7:12)    = xp(7:12) + dt*f*d(s);
      dxp         = B*dxp;            %update Jacobian
    end
  end
end