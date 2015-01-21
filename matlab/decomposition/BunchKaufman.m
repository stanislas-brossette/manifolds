function [P,L,B,Ps,Ls] = BunchKaufman(A)
n = size(A,1);
Bi = A;
B = zeros(n);
Ls = {};
Ps = {};
k = 0;
while k<n;
  [Pi,Li,Bi,E,s] = BKStep(Bi);
  B(k+1:k+s,k+1:k+s) = E;
  Ls{end+1} = [eye(k) zeros(k,n-k); zeros(n-k,k) Li];
  Ps{end+1} = [eye(k) zeros(k,n-k); zeros(n-k,k) Pi];
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

% function R = verif(P,L,B,Bi)
% k = length(L);
% ni = size(B,1)-size(Bi,1);
% R = B;
% R(ni+1:end,ni+1:end) = Bi;
% M = eye(size(B,1));
% for i=1:k
%   M = M*P{i}'*L{i};
% end
% R = M*R*M';
% end