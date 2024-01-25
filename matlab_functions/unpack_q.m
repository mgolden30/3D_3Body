function [qs, ps, T, alpha, N, g] = unpack_q( q )
  %{
  PURPOSE: standardize the deconstruction of a state since it is used quite
  a bit.
  %}

  N = (numel(q)-4)/12;

  %read position and momentum
  qs    = reshape(q(0*N + (1:6*N)), [6,N]);
  ps    = reshape(q(6*N + (1:6*N)), [6,N]);
  T     = q(12*N+1);
  alpha = q(12*N+(2:4)); %drift rate in SO(3)

  %Construct rotation generators
  Jz = [0,1,0;-1,0,0;0,0,0];
  Jy = [0,0,-1;0,0,0;1,0,0];
  Jx = [0,0,0;0,1,0;0,-1,0];

  %rotation generator
  g = alpha(1)*Jx + alpha(2)*Jy + alpha(3)*Jz;

end