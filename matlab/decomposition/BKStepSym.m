function [t,B,s] = BKStepSym(A,k)
n = size(A,1);
if n==1
  t = 1;
  L = 1;
  B = [];
  s = 1;
  return;
end
alpha = (1+sqrt(17))/8;
[lambda,r] = max(abs(A(k+2:n,k+1)));
r = r+k+1;
if lambda>0
  if abs(A(k+1,k+1))>=alpha*lambda
    s = 1; t = k+1;
  else
    [sigma,p] = max(abs([A(r,k+1:r-1)';A(r+1:n,r)]));
    p = p+k;
    if p>=r, p=p+1; end
    if sigma*abs(A(k+1,k+1)) >= alpha*lambda^2
      s = 1; t = k+1;
    elseif abs(A(r,r)) >= alpha*sigma
      s = 1; t = r;
    else
      s = 2; 
      if p==k+1 && r~=k+2, t = [k+1,r];
      elseif p==k+1 && r==k+2, t=[k+1,k+2];
      elseif p~=k+1 && r==k+2, t=[p,k+2];
      elseif p==k+2, t=[r,r];
      else t = [p,r];
      end
    end
  end
  
  Ap = A;
  for i=1:s
    %tmp = A(:,i); A(:,i) = A(:,t(i)); A(:,t(i))=tmp;
    %tmp = A(i,:); A(i,:) = A(t(i),:); A(t(i),:)=tmp;
    tmp=Ap(t(i)+1:end,k+i); Ap(t(i)+1:end,k+i)=Ap(t(i)+1:end,t(i)); Ap(t(i)+1:end,t(i))=tmp;
    tmp=Ap(k+i,1:k+i-1); Ap(k+i,1:k+i-1)=Ap(t(i),1:k+i-1); Ap(t(i),1:k+i-1)=tmp;
    tmp=Ap(k+i+1:t(i)-1,k+i); Ap(k+i+1:t(i)-1,k+i)=Ap(t(i),k+i+1:t(i)-1)'; Ap(t(i),k+i+1:t(i)-1)=tmp';
    tmp=Ap(k+i,k+i); Ap(k+i,k+i)=Ap(t(i),t(i)); Ap(t(i),t(i))=tmp;
  end
  %Ap-A
  E = Ap(k+1:k+s,k+1:k+s); if s==2, E(1,2)=E(2,1); end
  C = Ap(k+s+1:n,k+1:k+s);
  B = Ap(k+s+1:n,k+s+1:n);
  L = C/E;
  I = repmat([1:n-k-s]',1,n-k-s);
  J = repmat([1:n-k-s],n-k-s,1);
  B1 = B';
  B2 = B;
  B(I<J) = B1(I<J);
  B = B - (C/E)*C';
  B(I<J) = B2(I<J);
  Ap(k+s+1:n,k+s+1:n) = B;
  Ap(k+s+1:n,k+1:k+s) = L;
  B = Ap;
else
  s = 1;
  t = k+1;
  L = zeros(n-k-1,1);
  B = A;
end
B
end