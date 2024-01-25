%{
PURPOSE:
The C code finds many periodic orbits, but it does not check for
redundancy. Let's fix that here.
%}

clear;
restoredefaultpath();
addpath("matlab_functions/");


states = dir("states/*.txt");

%For each state I will compute a small set of time averaged observables to
%determine if the state has been found yet.


num_unique = 1;        %increment this every time a unique solution is found
obs = inf * ones(4,1); %eigenvalues of Q and 2-norm of l

mkdir("./unique_solutions")

for j = 1:numel( states )
  j/numel(states)
  filename = states(j).name;

  z = readmatrix( "states/" + filename );
  z = rescale_state(z); %rescale to H=-1
  N = 1024; %number of timesteps
  T = z(13);

  zs = zeros( 12, N );
  zs(:,1) = z(1:12);
  for i = 2:N
    zs(:,i) = symplectic_steps( zs(:,i-1), T/N, 1 );
  end

  r1 = zs(1:3, :);
  r2 = zs(4:6, :);
  r3 = -r1-r2;

  p1 = zs(7:9, :);
  p2 = zs(10:12,:);
  p3 = -p1-p2;

  Q = zeros(3,3);

  %Compute time averaged quantities to identify unique solutions
  for i = 1:3
    for k = 1:3
      Q(i,k) = mean( r1(i,:).*r1(k,:) + r2(i,:).*r2(k,:) + r3(i,:).*r3(k,:), 2 );
    end
  end

  l = mean( cross(r1,p1) + cross(r2,p2) + cross(r3,p3), 2 );
  
  %compute observables for this solution
  this_obs = [eigs(Q); norm(l)];

  if min(vecnorm( obs - this_obs )) < 0.05
    %we already know this solution
    continue;
  end

  %otherwise, we have a new solution
  obs(:, num_unique) = this_obs;
  system("cp states/" + filename + " unique_solutions/" + num_unique + ".txt" );
  num_unique = num_unique + 1;
end



function H = hamiltonian(x)
  r1 = x(1:3);
  r2 = x(4:6);
  p1 = x(7:9);
  p2 = x(10:12);

  r12 = r1-r2;
  r13 = 2*r1 + r2;
  r23 = 2*r2 + r1;

  d12 = norm(r12);
  d13 = norm(r13);
  d23 = norm(r23);

  H = sum(p1.^2) + sum(p2.^2) + sum(p1.*p2) - 1/d12 - 1/d13 - 1/d23;
end


function z = rescale_state(z)
  l = sqrt(abs(hamiltonian(z)));
  z(1:6)  = z(1:6) * l^2;
  z(7:12) = z(7:12) / l;
  z(13)   = z(13) * l^3;
end
