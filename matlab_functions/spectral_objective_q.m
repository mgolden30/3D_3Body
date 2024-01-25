function f = spectral_objective_q( q )
  %{
  PURPOSE:
  This defines a function f(q) such that f(q) corresponds to periodic
  orbits of the three body problem.

  INPUT:
  q - see pack_q

  OUTPUT:
  f - equations of motion for an RPO in integral form
  %}
  

  [qs, ps, T, alpha, N, g] = unpack_q( q );

  r1 = qs(1:3,:);
  r2 = qs(4:6,:);

  r12 =   r1 -   r2;
  r13 = 2*r1 +   r2;
  r23 =   r1 + 2*r2;
  
  %vecnorm is not allowed with vpa
  %d12 = vecnorm(r12);
  %d13 = vecnorm(r13);
  %d23 = vecnorm(r23);  

  d12 = sqrt( sum(r12.^2) );
  d23 = sqrt( sum(r23.^2) );
  d13 = sqrt( sum(r13.^2) );

  F = 0*qs;
  F(1:3,:) = -r12./d12.^3 - r13./d13.^3;
  F(4:6,:) =  r12./d12.^3 - r23./d23.^3;
  
  a = 2*pi/T; %time conversion factor
  k = 0:N-1; k(k>N/2) = k(k>N/2) - N;

  f1 = @(x) fft(x,N,2);
  f2 = @(x) real(ifft(x,N,2));

  %Define an operator to undo a time derivative 
  % (and preserve the zero mode)
  dt_inv = 1./(1i*k);
  dt_inv(1) = 1; %preserve the zero mode

  %Make it a 6x6 matrix
  g = [g, zeros(3,3); zeros(3,3), g];
  
  %Add rotation terms usual right hand sides
  rhs_q = ps + g * qs;
  rhs_p = F  + g * ps;

  %Let fq and fp be the components of f corresponding to position and
  %momentum, repectively.
  fq = a*(qs - mean(qs,2)) - f2( dt_inv .* f1( rhs_q )) ;
  fp = a*(ps - mean(ps,2)) - f2( dt_inv .* f1( rhs_p ));

  %Set up slicing conditions
  %First enforce mean H = -1
  %Note that even though we are in a rotated frame, the Hamiltonian is the
  %same.
  p1 = ps(1:3,:);
  p2 = ps(4:6,:);
  H = sum(p1.^2 + p2.^2 + p1.*p2) - 1./d12 - 1./d23 - 1./d13;
  fT = mean(H) + 1;

  p1 = ps(1:3,1);

  %Rotate r1 to point along x axis
  %and p1 to be in the x-y plane
  %falpha = [r1(2); r1(3); p1(3)];

  %The above does not work well in some situations
  falpha = alpha;

  f = pack_q( fq, fp, fT, falpha, N );
end