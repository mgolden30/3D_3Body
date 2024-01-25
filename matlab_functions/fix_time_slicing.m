function q = fix_time_slicing( q )
  %{
  The goal of this function is to perform the nontrivial operation of
  shifting q in time such that q(t=0) corresponds to maximal hyperradius

  R = r1^2 + r2^2 + r3^2
  %}


  N  = (numel(q)-4)/6;
  qs = reshape( q(1:6*N), [6,N] );

  tiledlayout(1,2);
  nexttile

  r1 = qs(1:3,:);
  r2 = qs(4:6,:);
  r3 = -r1-r2;

  %Compute gridded R
  R = sum(r1.^2 + r2.^2 + r3.^3);
 
  %find the maximum value
  [~, idx] = max(R);

  %Do a linear fit near R and interpolate
  RR = @(i) R(mod(i-1,N)+1); %To allow indices to be arbitrary integers
  
  fit_x = -10:10;
  Rs = RR(idx + fit_x);

  %plot(fit_x, Rs)

  p  = polyfit( fit_x, Rs, 2);
  dp = [2*p(1), p(2)];

  idx_true = idx + roots(dp);

  if( abs(roots(dp)) > 10)
    %Don't trust. No shift needed.
    return;
  end

  %Shift using Fourier Transform
  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;

  qs = fft(qs, N, 2);
  %size(qs)
  qs = exp( 1i*k*idx_true/N) .* qs;
  %size(qs)
  qs = real(ifft(qs,N,2));
  %size(qs)

  %Update q with the translated values
  q(1:6*N) = reshape( qs, [6*N,1] );
end