function [z, converge_flag] = super_newton(z,N,maxit,tol,maxit2)
  converge_flag = true; %default

  for i = 1:maxit
    [f, df, ddf] = objective(z, N);

    %Use regular Newton output as initial guess guess
    [dz, ~] = lsqr( df, -f, tol );

    for j = 1:32
       
      %Keep the intermediate calculation of multiplying hessian by step
      df2 = reshape(reshape(ddf, [8*9,9]) * dz, [8,9] ); 
      
      %second order Taylor expansion
      f2  = f + df*dz + 1/2 * df2*dz;
      
      if( j == 1)
        err0 = norm(f2);
      end

      total_df = df + df2;
      [ddz, ~] = lsqr( total_df, -f2, tol);
      dz = dz + ddz; %Newton step
    end
    err = norm(f2);

    fprintf("step %d: |f| = %e\t err0 = %e\t err = %e\t |step| = %e\n", i, norm(f), err0, err, norm(dz));

    z = z + dz;

    
    damp = 1;
    for k = 1:maxit2
      z2 = z + damp*dz;
      f2 = objective(z2,N);

      if( norm(f2) < norm(f) )
        break;
      end
      
      damp = damp/2;
    end
    z = z2;

  end
end