function H = mean_H(q)
  N  = (numel(q) - 4)/6;
  qs = reshape(q(1:6*N), [6,N]);
  T  = q(6*N+1);
  alpha = q(6*N+(2:4));

  %Need Hamiltonian
  k = 0:N-1; k(k>N/2) = k(k>N/2) - N;

  f1 = @(x) fft(x,N,2);
  f2 = @(x) real(ifft(x,N,2));
  
  ps = f2( 1i*k*(2*pi)/T.*f1(qs) );

  %Construct rotation generators
  Jx = [0,1,0;-1,0,0;0,0,0];
  Jy = [0,0,-1;0,0,0;1,0,0];
  Jz = [0,0,0;0,1,0;0,-1,0];

  %rotation generator
  alpha = q(6*N + (2:4));
  g = alpha(1)*Jx + alpha(2)*Jy + alpha(3)*Jz;

  %Compute inertial momentum
  ps_inertial = ps - [g, zeros(3,3); zeros(3,3), g]*qs;

  r1 = qs(1:3, :);
  r2 = qs(4:6, :);
  
  d12 = vecnorm(   r1 -   r2 );
  d13 = vecnorm( 2*r1 +   r2 );
  d23 = vecnorm(   r1 + 2*r2 );

  p1 = ps_inertial(1:3, :);
  p2 = ps_inertial(4:6, :);

  H = sum(p1.^2 + p2.^2 + p1.*p2) - 1./d12 - 1./d13 - 1./d23;

  H = mean(H);
end