function xs = plot_traj( x, dt, N)
  xs = zeros(12,N);
  xs(:,1) = x;
  for i = 2:N
    [xs(:,i), ~] = symplectic_steps(xs(:,i-1), dt, 1);
  end 

  r1 = xs(1:3,:);
  r2 = xs(4:6,:);
  
  %Find a rotation to diagonalize a rank-2 tensor
  V = find_rotation( r1, r2 );
  
  r1 = V.'*r1;
  r2 = V.'*r2;
  r3 = -r1-r2; %Don't need to rotate since we already rotated r1 and r2
  
  %grey = [1 1 1] * 0.25;
  plt  = @(rr, c) plot3( rr(1,:), rr(2,:), rr(3,:), 'color', c, 'linewidth', 5 );
  %plt2 = @(r)    plot3( r(1,:), r(2,:), r(3,:), 'color', grey, 'linewidth', 3, 'LineStyle', '--' );

  q = 0.4;

  plt(r1(:,1:i), [1 q q]);
  hold on
  plt(r2(:,1:i), [q q 1]);
  plt(r3(:,1:i), [q 1 q]);
  hold off

  xlim([ -3 3 ]);
  ylim([ -3 3 ]);
  zlim([ -3 3 ]);
  
  pbaspect([1 1 1]);
end