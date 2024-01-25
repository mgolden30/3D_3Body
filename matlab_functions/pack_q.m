function q = pack_q( qs, ps, T, alpha, N )
  %{
  PURPOSE: standardize the deconstruction of a state since it is used quite
  a bit.
  %}

  q = [ reshape(qs, [6*N,1]); 
        reshape(ps, [6*N,1]); 
        T; 
        alpha ];
end