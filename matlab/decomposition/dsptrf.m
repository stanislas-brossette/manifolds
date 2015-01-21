function [Ap,ipiv,info] = dsptrf(Ap, UpLo)
l = length(Ap);
n = sqrt(1/4+2*l) - 1/2;
info = 0;
%check entry
upper = UpLo=='u' || UpLo=='U';
if (~upper && UpLo~='l' && UpLo~='L')
  info = -1;
elseif n~=floor(n)
  info = -2;
end
if info~=0
  ipiv = [];
  return
end

ipiv = zeros(1,n);
alpha = (1+sqrt(17))/8;

if upper
  assert(false, 'non implemented yet')
else
  k = 1;
  kc = 1;
  npp = n*(n+1)/2;
  
  while k<=n
    knc = kc;
    kstep = 1;
    
    %Determine rows and columns to be interchanged and whether
    %a 1-by-1 or 2-by-2 pivot block will be used
    
    absakk = abs(Ap(kc));
    
    %IMAX is the row-index of the largest off-diagonal element in
    %column K, and COLMAX is its absolute value
    if k<n
      [colmax,imax] = max(abs(Ap(kc+1:kc+n-k)));
      imax = k + imax;
    else
      colmax = 0;
    end
    
    if max(absakk,colmax) == 0
      %Column K is zero: set INFO and continue
      if info == 0
        info = k;
      end
      kp = k;
    else
      if absakk >= alpha*colmax
        %no interchange, use 1-by-1 pivot block
        kp = k;
      else
        %JMAX is the column-index of the largest off-diagonal
        %element in row IMAX, and ROWMAX is its absolute value
        rowmax = 0;
        kx = kc + imax - k;
        for j=k:imax-1
          if abs(Ap(kx))>rowmax
            rowmax = abs(Ap(kx));
            jmax = j;
          end
          kx = kx + n- j;
        end
        kpc = npp - (n-imax+1)*(n-imax+2)/2 + 1;
        if imax<n
          [rowmax2,jmax] = max(abs(Ap(kpc+1:kpc+n-imax)));
          rowmax = max(rowmax,rowmax2);
          jmax = imax + jmax; %jmx is not used later, we can remove this
        end
        
        if absakk >= alpha*colmax*(colmax/rowmax)
          %no interchange, use 1-by-1 pivot block
          kp = k;
        elseif abs(Ap(kpc))>=alpha*rowmax
          %interchange rows and columns K and IMAX, use 1-by-1 pivot block
          kp = imax;
        else
          %interchange rows and columns K+1 and IMAX, use 2-by-2 pivot
          %block
          kp = imax;
          kstep = 2;
        end
      end
      
      kk = k+kstep-1;
      if kstep==2
        knc = knc + n - k + 1;
      end
      if kp~=kk
        %Interchange rows and columns KK and KP in the trailing submatrix A(k:n,k:n)
        if kp < n
          tmp = Ap(kpc+1:kpc+n-kp);                         %CALL DSWAP
          Ap(kpc+1:kpc+n-kp) = Ap(knc+kp-kk+1:knc+n-kk);
          Ap(knc+kp-kk+1:knc+n-kk) = tmp;
        end
        
        kx = knc + kp-kk;
        for j=kk+1:kp-1
          kx = kx + n-j+1;
          t = Ap(knc+j-kk);
          Ap(knc+j-kk) = Ap(kx);
          Ap(kx) = t;
        end
        t = Ap(knc);
        Ap(knc) = Ap(kpc);
        Ap(kpc) = t;
        if kstep == 2
          t = Ap(kc+1);
          Ap(kc+1) = Ap(kc+kp-k);
          Ap(kc+kp-k) = t;
        end
      end
    
      %Update the trailing submatrix
      if kstep == 1
        %1-by-1 pivot block D(k): column k now holds W(k) = L(k)*D(k)
        %where L(k) is the k-th column of L
        if k<n
          %Perform a rank-1 update of A(k+1:n,k+1:n) as
          %A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T
          r1 = 1/Ap(kc);
          Ap(kc+n-k+1:end) = dspr(Ap(kc+n-k+1:end),-r1,Ap(kc+1:kc+n-k),1,UpLo);
          Ap(kc+1:kc+n-k) = r1*Ap(kc+1:kc+n-k); %CALL DSCAL
        end
      else
        %2-by-2 pivot block D(k): columns K and K+1 now hold
        %( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
        %where L(k) and L(k+1) are the k-th and (k+1)-th columns of L
        if k<n-1
          %Perform a rank-2 update of A(k+2:n,k+2:n) as
          %A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T
          %   = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T
          %where L(k) and L(k+1) are the k-th and (k+1)-th
          %columns of L
          D21 = Ap(k+1+(k-1)*(2*n-k)/2);
          D11 = Ap(k+1+k*(2*n-k-1)/2)/D21;
          D22 = Ap(k+(k-1)*(2*n-k)/2)/D21;
          t = 1/(D11*D22-1);
          D21 = t/D21;

          for j=k+2:n
            wk = D21*(D11*Ap(j+(k-1)*(2*n-k)/2) - Ap(j+k*(2*n-k-1)/2));
            wkp1 = D21*(D22*Ap(j+k*(2*n-k-1)/2) - Ap(j+(k-1)*(2*n-k)/2));

            for i=j:n
              Ap(i+(j-1)*(2*n-j)/2) = Ap(i+(j-1)*(2*n-j)/2) ...
                                    - Ap(i+(k-1)*(2*n-k)/2)*wk ...
                                    - Ap(i+k*(2*n-k-1)/2)*wkp1;
            end
            Ap(j+(k-1)*(2*n-k)/2) = wk;
            Ap(j+k*(2*n-k-1)/2) = wkp1;
          end
        end
      end
    end
  
    %Store details of the interchanges in IPIV
    if kstep==1
      ipiv(k) = kp;
    else
      ipiv(k) = -kp;
      ipiv(k+1) = -kp;
    end
    
    %ncrease K and return to the start of the main loop
    k = k+kstep;
    kc = knc+n-k+2;
  end
end
end