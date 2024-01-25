function q = fix_rotation_q( q )
  %{
  PURPOSE:
  For any exact solution to the equations of motion, there is a degree of
  freedom 
  %}


  %get position and momentum out of q in an idiot-proof fashion
  [qs, ps, T, ~, N, g] = unpack_q( q );

  %Define a rotation macro to take a 3x3 matrix and apply it to a 12xN matrix
  rotate_traj = @(R, traj) reshape( R * reshape( traj, [3,2*N]) , [6,N]);
 
  %take intial r1 and rotate it to point along x axis and make p1 point in
  %x-y plane. QR decomposition makes this easy!
  r1 = qs(1:3,1);
  p1 = ps(1:3,1);
  A = [r1,p1];
  [Q,~] = qr(A);

  %Now actually rotate everything
  qs = rotate_traj( Q.', qs );
  ps = rotate_traj( Q.', ps );

  %Matrices need multiplied on both sides
  g = Q.' * g * Q;

  %See unpack_state to make sure this is right
  alpha =  [g(2,3); -g(1,3); g(1,2)];
  
  %pack everyhting back up in the rotated frame
  q = pack_q( qs, ps, T, alpha, N );
end