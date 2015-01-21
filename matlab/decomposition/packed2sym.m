%lapack-like packed vector to symetric matrix
function S = packed2sym(p,UpLo)
l = length(p);
n = sqrt(1/4+2*l) - 1/2;
assert(n==floor(n), 'invalid vector size');
S = zeros(n,n);
if (UpLo=='l' || UpLo=='L')
  k = 0;
  for i=1:n
    S(i:n,i) = p(k+1:k+(n-i+1));
    S(i,i:n) = S(i:n,i)';
    k = k + n - i + 1;
  end
elseif (UpLo=='u' || UpLo=='U')
  k=0;
  for i=1:n
    S(1:i,i) = p(k+1:k+i);
    S(i,1:i) = S(1:i,i)';
    k = k+i;
  end
else
  assert(false);
end
end