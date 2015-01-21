function [P,L,B,E,s] = BKStep(A)
n = size(A,1);
if n==1
  P = 1;
  L = 1;
  E = A;
  B = [];
  s = 1;
  return;
end
alpha = (1+sqrt(17))/8;
[lambda,r] = max(abs(A(2:n,1)));
r = r+1;
if lambda>0
  if abs(A(1,1))>=alpha*lambda
    s = 1; P = eye(n);
  else
    [sigma,p] = max(abs([A(1:r-1,r);A(r+1:n,r)]));
    if p>=r, p=p+1; end
    if sigma*abs(A(1,1)) >= alpha*lambda^2
      s = 1; P = eye(n);
    elseif abs(A(r,r)) >= alpha*sigma
      s = 1; P = eye(n); P([1,r],[1,r]) = [0 1;1 0];
    else
      s = 2; P = eye(n); 
      if p==1 && r~=2, P([2,r],[2,r]) = [0 1;1 0];
      elseif p==1 && r==2;
      elseif p~=1 && r==2, P([1,p],[1,p]) = [0 1;1 0];
      elseif p==2, P([1,2,r],[1,2,r]) = [0,0,1;1,0,0;0,1,0];
      else P([1,p],[1,p]) = [0 1;1 0]; P([2,r],[2,r]) = [0 1;1 0];
      end
    end
  end
  
  Ap = P*A*P';
  E = Ap(1:s,1:s);
  C = Ap(s+1:n,1:s);
  B = Ap(s+1:n,s+1:n);
  L = [eye(s), zeros(s,n-s); C/E, eye(n-s)];
  B = B - (C/E)*C';
  %B
else
  s = 1;
  P = eye(n);
  L = eye(n);
  E = A(1,1);
  B = A(2:n,2:n);
end
end