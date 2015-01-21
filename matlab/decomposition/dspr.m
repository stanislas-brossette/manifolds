function Ap = dspr(Ap,alpha,x,incx,UpLo)
l = length(Ap);
n = sqrt(1/4+2*l) - 1/2;
info = 0;

%check entry
upper = UpLo=='u' || UpLo=='U';
if (~upper && UpLo~='l' && UpLo~='L')
  info = -1;
elseif n~=floor(n)
  info = -2;
elseif incx==0
  info = -5;
end
if info~=0
  assert(false);
end

%quick return if possible
if n==0 || alpha==0, return; end

%set the start point in x if the increment is not unity.
if incx <= 0
  kx = 1 - (n-1)*incx;
elseif incx~=1
  kx = 1;
end

%Start the operations. In this version the elements of the array Ap
%are accessed sequentially with one pass through Ap.

kk=1;
if upper
  %Form A when upper triangle is stored in Ap
  if incx==1
    for j=1:n
      if x(j) ~= 0
        temp = alpha*x(j);
        k = kk;
        for i=1:j
          Ap(k) = Ap(k) + x(i)*temp;
          k = k+1;
        end
      end
      kk = kk + j;
    end
  else
    jx = kx;
    for j=1:n
      if (x(jx)~=0)
        temp = alpha*x(jx);
        ix = kx;
        for k = kk:kk+j-1
          Ap(k) = Ap(k) + x(ix)*temp;
          ix = ix + incx;
        end
      end
      jx = jx + incx;
      kk = kk+j;
    end
  end
else
  %Form A when lower triangle is stored in AP.
  if incx==1
    for j=1:n
      if x(j) ~= 0
        temp = alpha*x(j);
        k = kk;
        for i=j:n
          Ap(k) = Ap(k) + x(i)*temp;
          k = k+1;
        end
      end
      kk = kk + n-j+1;
    end
  else
    jx = kx;
    for j=1:n
      if x(jx) ~= 0
        temp = alpha*x(jx);
        ix = jx;
        for k = kk:kk+n-j
          Ap(k) = Ap(k) + x(ix)*temp;
          ix = ix + incx;
        end
      end
      jx = jx + incx;
      kk = kk+n-j+1;
    end
  end
end
end