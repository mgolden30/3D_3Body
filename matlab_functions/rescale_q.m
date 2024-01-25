function q = rescale_q(q)
  %{
  PURPOSE:
  We wish to fix the Hamiltonian H=-1 for periodic orbits. We can use the
  rescaling symmetry of n-body dynamics.
  %}

  %Decompose the q vector into useful physical quantities.
  [qs, ps, T, alpha, N, g] = unpack_q( q );

  r1 = qs(1:3,:);
  r2 = qs(4:6,:);

  r12 =   r1 -   r2;
  r13 = 2*r1 +   r2;
  r23 =   r1 + 2*r2;
  
  d12 = vecnorm(r12);
  d13 = vecnorm(r13);
  d23 = vecnorm(r23);  
  
  %Make it a 6x6 matrix
  g = [g, zeros(3,3); zeros(3,3), g];

  %Set up slicing conditions
  %First enforce mean H = -1
  %Note that even though we are in a rotated frame, the Hamiltonian is the
  %same.
  p1 = ps(1:3,:);
  p2 = ps(4:6,:);
  H = sum(p1.^2 + p2.^2 + p1.*p2) - 1./d12 - 1./d23 - 1./d13;

  l = sqrt(abs(H));
  lm = mean(l); 

  %Rescale each step so that H=-1 to numerical precision.
  q(       1:6*N ) = reshape( qs.*l.*l, [6*N, 1] );
  q(6*N + (1:6*N)) = reshape( ps./l,   [6*N, 1] );
  
  q(12*N+1)     = T*lm*lm*lm;
  q(12*N+(2:4)) = alpha/lm/lm/lm;
end