function H = hamiltonian(x)
  r1 = x(1:2);
  r2 = x(3:4);
  p1 = x(5:6);
  p2 = x(7:8);

  r12 = r1-r2;
  r13 = 2*r1 + r2;
  r23 = 2*r2 + r1;

  d12 = norm(r12);
  d13 = norm(r13);
  d23 = norm(r23);

  H = sum(p1.^2) + sum(p2.^2) + sum(p1.*p2) - 1/d12 - 1/d13 - 1/d23;
end