function visualize_q(q)
  %{
  PURPOSE:
  Visualize everything you need to know about this orbit. Make it sexy.
  %}

  %
  [qs, ps, T, ~, N, g] = unpack_q( q );

  r1 = qs(1:3,:);
  r2 = qs(4:6,:);
  r3 = -r1-r2;

  plt  = @(rr, c) plot3( rr(1,:), rr(2,:), rr(3,:), 'color', c, 'linewidth', 5 );
  
  qq = 0.4;
  tiledlayout(2,2)

  nexttile
  plt(r1, [1 qq qq]);
  hold on
  plt(r2, [qq qq 1]);
  plt(r3, [qq 1 qq]);
  hold off

  xlim([ -3 3 ]);
  ylim([ -3 3 ]);
  zlim([ -3 3 ]);
  
  pbaspect([1 1 1]);
  title('corotating frame');

  %% Lab frame orbit
  nexttile
  
  r1_lab = r1;
  r2_lab = r2;

  dt = T/N;
  for i = 2:N
    rot = expm( -dt*(i-1)*g ); %matrix exponential to get the full rotation
    r1_lab(:,i) = rot*r1_lab(:,i);
    r2_lab(:,i) = rot*r2_lab(:,i);
  end
  r3_lab = -r1_lab-r2_lab;

  plt(r1_lab, [1 qq qq]);
  hold on
  plt(r2_lab, [qq qq 1]);
  plt(r3_lab, [qq 1 qq]);
  hold off

  xlim([ -3 3 ]);
  ylim([ -3 3 ]);
  zlim([ -3 3 ]);
  
  pbaspect([1 1 1]);
  title('lab frame');


  %% spectral scaling of coefficients
  nexttile

  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;
  qs2 = fft( qs, N, 2 );
  ps2 = fft( qs, N, 2 );

  semilogy( k, abs(qs2) );
  hold on
  semilogy( k, abs(ps2) );
  hold off
  xlabel('n');
  ylabel('|c_n|')

  xlim([0, N/2]);
  ylim( [ 1e-16, max(ylim)] );
  title('spectral scaling (in corotating)')
  %legend();

  %% Energy Conservation
  nexttile
  
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

  semilogy( abs(H + 1) );
  ylabel("|H+1|");
  xlabel("timestep");
  title("energy conservation");
end