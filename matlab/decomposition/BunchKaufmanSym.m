function [P,L,B,t,Ls] = BunchKaufmanSym(A)
n = size(A,1);
Bi = A;
B = zeros(n);
Ls = {};
s = {};
t = 1:n;
k = 0;
while k<n;
  [ti,Li,Bi,E,si] = BKStepSym(Bi,k);
  B(k+1:k+si,k+1:k+si) = E;
  Ls{end+1} = [eye(k) zeros(k,n-k); zeros(n-k,k) Li];
  t(k+1:k+si) = ti;
  k = k+si;
  s{end+1} = si;
end
[P,L] = buildPLBrute(t,Ls,s);
end

function [P,L] = buildPLBrute(t,Ls,s)
n = size(Ls{1},1);
L = eye(n);
P = eye(n);
k=0;
for i=1:length(Ls)
  Pi = eye(n);
  si = s{i};
  for j=1:si
    tmp = Pi(k+j,:); Pi(k+j,:) = Pi(t(k+j),:); Pi(t(k+j),:) = tmp;
  end
  L = Pi*L*Pi'*Ls{i};
  P = Pi*P;
  k = k+si;
end
end

