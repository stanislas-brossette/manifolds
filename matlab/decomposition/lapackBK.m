function [P,L,B,Ps,Ls] = lapackBK(A)
n = size(A,1);
p = sym2packed(A,'l');
[p,ipiv,info] = dsptrf(p,'l');

%reconstruct from packed result
k=0;
l=0;
R = packed2sym(p,'l');
Ls = {};
Ps = {};
B = zeros(n,n);
while k<n
  Ps{end+1} = eye(n);
  kpiv = ipiv(k+1);
  if kpiv<0
    s=2;
    if -kpiv~=k+2
      Ps{end}([k+2,-kpiv],[k+2,-kpiv]) = [0 1; 1 0];
    end
  else
    s=1;
    if kpiv~=k+1
      Ps{end}([k+1,kpiv],[k+1,kpiv]) = [0 1; 1 0];
    end
  end
  Ls{end+1} = [ eye(k)              zeros(k,s)        zeros(k,n-s-k);...
                zeros(s,k)            eye(s)          zeros(s,n-s-k);...
                zeros(n-s-k,k)  R(k+s+1:end,k+1:k+s)    eye(n-s-k)  ];
  B(k+1:k+s,k+1:k+s) = R(k+1:k+s,k+1:k+s);
  k = k+s;
end
[P,L] = buildPLBrute(Ps,Ls);
end


function [P,L] = buildPLBrute(Ps,Ls)
L = eye(size(Ps{1},1));
P = L;
for i=1:length(Ps)
  L = Ps{i}*L*Ps{i}'*Ls{i};
  P = Ps{i}*P;
end
end