function f = spectral_objective(, q)
  %{
  PURPOSE:
  This defines a function f(q) such that f(q) corresponds to periodic
  orbits of the three body problem.

  INPUT:
  m - a huge vector of size 144*(N+1) corresponding to the Floquet matrix
      in time. 

  q - a huge vector of size 6*N+4. It is made up of the following:
      qs - 6*N is the timeseries of [r1; r2] where r1 and r2 are size [3,N]
      T  - 1 is the period
      alpha - 3 are the SO(3) generators to allow for RPOs

  OUTPUT:
  f - 

  %}
  
  N = (numel(q)-4)/6;

  f     = 0*q;
  qs    = reshape(q(1:6*N), [6,N]);
  T     = q(6*N+1);
  alpha = q(6*N+(2:4)); %drift rate in SO(3)

  %Construct rotation generators
  Jz = [0,1,0;-1,0,0;0,0,0];
  Jy = [0,0,-1;0,0,0;1,0,0];
  Jx = [0,0,0;0,1,0;0,-1,0];

  %rotation generator
  g = alpha(1)*Jx + alpha(2)*Jy + alpha(3)*Jz;

  r1 = qs(1:3,:);
  r2 = qs(4:6,:);

  r12 =   r1 -   r2;
  r13 = 2*r1 +   r2;
  r23 =   r1 + 2*r2;
  
  d12 = vecnorm(r12);
  d13 = vecnorm(r13);
  d23 = vecnorm(r23);  

  F = 0*qs;
  F(1:3,:) = -r12./d12.^3 - r13./d13.^3;
  F(4:6,:) =  r12./d12.^3 - r23./d23.^3;
  
  
  a = 2*pi/T; %conversion factor
  k = 0:N-1; k(k>N/2) = k(k>N/2) - N;

  f1 = @(x) fft(x,N,2);
  f2 = @(x) real(ifft(x,N,2));

  %Define an operator to undo two time derivatives 
  % (and preserve the zero mode)
  k2_inv = -1./k.^2;
  k2_inv(1) = 1; %preserve the zero mode

  %compute the time derivative of qs
  dqs = a * f2( 1i*k.*f1(qs) );

  %momentum is the time derivative with rotation subtracted off
  ps = dqs - [g*qs(1:3,:); g*qs(4:6,:) ];

  %Make it a 6x6 matrix
  g = [g, zeros(3,3); zeros(3,3), g];

  %Compute inertial momentum 
  ps_inertial = ps - g*qs;
  
  %Add rotation terms to force
  rhs = F + 2*g*ps + g*g*qs;

  obj = a*a*(qs - mean(qs, 2)) - f2( k2_inv .* f1(rhs) ) ;

  f(1:6*N) = reshape( obj, [6*N,1] );

  %Set up slicing conditions
  %First enforce mean H = -1
  %Note that even though we are in a rotated frame, the Hamiltonian is the
  %same.
  p1 = ps(1:3,:);
  p2 = ps(4:6,:);
  H = sum(p1.^2 + p2.^2 + p1.*p2) - 1./d12 - 1./d23 - 1./d13;
  f(6*N+1) = mean(H) + 1;

  p1 = ps_inertial(1:3,1);

  %Force r1 to point along x axis
  f(6*N+2) = r1(2);
  f(6*N+3) = r1(3);

  %force p1 to be in the x-y plane
  f(6*N+4) = p1(3);
end