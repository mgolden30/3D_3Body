%{
Sample from a distribution

 \int d^nx \, \rho_0(x) = \int d^nx \, \delta(H+1)

%}

R = 3;
P = 3;
N = 1024;

[samples, weights] = generate_samples( R, P, N );

%Compute the mean angular momentum
r1 = samples(1:3, :);
r2 = samples(4:6, :);
p1 = samples(7:9, :);
p2 = samples(10:12, :);

r3 = -r1-r2;
p3 = -p1-p2;

L = cross(r1,p1) + cross(r2,p2) + cross(r3,p3);

L_ave = sum(L.*weights, 2)/N

function [samples, weights] = generate_samples( R, P, N )
  %{ 
  Randomly sample a natural distribution for the three body problem.
  %}

  samples = zeros( 12, N );
  weights = zeros(  1, N );

  i = 1;
  while(i<=N)
    r1 = random_vector( R, 3 );
    r2 = random_vector( R, 3 );
    r3 = -r1-r2;
    if norm(r3)>R
      continue;
    end

    p1 = random_vector( P,3 );
    p2 = zeros(3,1);
    p2(1:2) = random_vector( P,2 ); %solve for z component

    %Get the last component of momentum by fixing energy to -1
    V = -1/norm(r1-r2) - 1/norm(r1-r3) - 1/norm(r2-r3);
    
    %Hamiltonian is p1^2 + p2^2 + p1*p2 + V
    %This is a quadratic equation for p2z
    a = 1;
    b = p1(3);
    c = sum(p1.^2) + sum(p1(1:2).*p2(1:2)) + sum(p2(1:2).^2) + V;

    if b*b - 4*a*c < 0
      continue; %there is no solution for H=-1
    end

    sgn = sign( 2*rand() - 1 );

    p2(3) = (-b + sgn*sqrt(b*b - 4*a*c))/(2*a);

    p3 = -p1-p2;
    if( norm(p3) > P )
      continue;
    end

    %If you made it this far, we have a valid sample
    samples(:,i) = [r1;r2;p1;p2];
    weights(i)   = 1/2/abs(p1(3) + p2(3));
    i = i + 1;
  end
end


function v = random_vector( R, dim )
  %Generate a random vector in the sphere of radius R
  
  v = R*(2*rand(dim,1) - 1);
  if( norm(v) > R )
    %try again
    v = random_vector( R, dim );
  end
end