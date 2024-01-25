function V = find_rotation( r1, r2 )
  r3 = -r1-r2;
  
  N = size(r1,2);
  r = zeros(3,N,3);
  r(:,:,1) = r1;
  r(:,:,2) = r2;
  r(:,:,3) = r3;

  
  %compute the mean r_i r_j tensor
  Q = zeros(3,3);
  for i = 1:3
    for j = 1:3
        for k = 1:3
          Q(i,j) = Q(i,j) + mean( r(i,:,k).*r(j,:,k) );
        end
    end
  end

  %rotate state to diagonalize Q
  [V,~] = eigs(Q,3);
end