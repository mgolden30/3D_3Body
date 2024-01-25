function q = change_sampling(q, N_new)
  %{
  PURPOSE:
  Perhaps you need more or fewer points in time to resolve a solution.
  This function allows you to painlessly change resolution.
  %}

  [qs, ps, T, alpha, N, ~] = unpack_q( q );

  if( N == N_new )
    %Nothing to do
    return;
  end

  %Only allow even N_new
  assert( mod(N_new,2) == 0);

  %k is the relevant number of modes we will copy
  k = min(N,N_new)/2;

  qs_new = zeros(6, N_new);
  ps_new = zeros(6, N_new);

  %Take fft of function input
  qs = fft( qs, N, 2 );
  ps = fft( ps, N, 2 );

  qs_new( :,  1:k      ) = qs(:, 1:k );
  qs_new( :, end-k+1:end  ) = qs(:, end-k+1:end );
  ps_new( :,  1:k      ) = ps(:, 1:k );
  ps_new( :, end-k+1:end ) = ps(:, end-k+1:end );
    
  qs_new = real(ifft(qs_new,N_new,2));
  ps_new = real(ifft(ps_new,N_new,2));

  %Need to rescale by ratio of N's
  qs_new = qs_new/N*N_new;
  ps_new = ps_new/N*N_new;

  q = pack_q( qs_new, ps_new, T, alpha, N_new );
end