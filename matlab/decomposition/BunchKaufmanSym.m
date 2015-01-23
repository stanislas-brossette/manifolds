function [P,L,B,t,Ls,Bi] = BunchKaufmanSym(A)
n = size(A,1);
Bi = A;
B = zeros(n);
Ls = {};
s = {};
t = 1:n;
k = 0;
while k<n;
  [ti,Bi,si] = BKStepSym(Bi,k);
  B(k+1:k+si,k+1:k+si) = Bi(k+1:k+si,k+1:k+si); if si==2, B(k+1,k+2)=B(k+2,k+1); end
  Ls{end+1} = [eye(k) zeros(k,n-k); zeros(si,k) eye(si) zeros(si,n-si-k); ...
                    zeros(n-k-si,k) Bi(k+si+1:n,k+1:k+si) eye(n-k-si)];
  t(k+1:k+si) = ti;
  k = k+si;
  s{end+1} = si;
end
[P,L] = buildPLBrute(t,Bi,s);
end

function [P,L] = buildPLBrute(t,A,s)
n = size(A,1);
L = eye(n);
P = eye(n);
k=0;
for i=1:length(s)
  Pi = eye(n);
  si = s{i};
  for j=1:si
    tmp = Pi(k+j,:); Pi(k+j,:) = Pi(t(k+j),:); Pi(t(k+j),:) = tmp;
  end
  L(k+si+1:n,k+1:k+si) = A(k+si+1:n,k+1:k+si);
  P = Pi*P;
  k = k+si;
end
end

