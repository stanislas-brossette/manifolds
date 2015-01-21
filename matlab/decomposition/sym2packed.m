%symetric matrix to lapack-like packed vector
function p = sym2packed(S,UpLo)
[m,n] = size(S);
assert(m==n);
p = zeros(n*(n+1)/2,1);
if (UpLo=='l' || UpLo=='L')
  k = 0;
  for i=1:n
    p(k+1:k+(n-i+1)) = S(i:n,i);
    k = k + n - i + 1;
  end
elseif (UpLo=='u' || UpLo=='U')
  k=0;
  for i=1:n
    p(k+1:k+i) = S(1:i,i);
    k = k+i;
  end
else
  assert(false);
end
end