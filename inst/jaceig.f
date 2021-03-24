c
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        code for the evaluation of eigenvalues and eigenvectors
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine jaceig(a,n,eps,ifvect,uout,nelems)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(2,2),u(2,2),rlam(4),uout(n,n)
c 
c        This subroutine uses Jacobi iterations to find the
c        spectrum of a real symmetric matrix. It is a primitive
c        and robust code. It is also mercifully short.
c        IMPORTANT! The dimensionality of the matrix must be
c        .ge. 3. Otherwise, a terrible bomb will happen!
c 
c                 Input parameters:
c 
c  a - the matrix whose spectrum is to be found; destroyed by the
c        subroutine
c  n - the dimensionality of a; must be .ge. 3
c  eps - the accuracy (relative) to which the calculations are to
c        be conducted
c  ifvect - the integer parameter telling the subroutine whether it
c        should compute the eigenvectors (or only the eigenvalues)
c        ifvect=1 will cause the eigenvectors to be computed
c        ifvect=0 will suppress the computation of the eigenvectors
c 
c                 Output parameters:
c 
c  a - the first n elements of a contain the (ordered) spectrum of a
c  uout - the eigenvectors of a corresponding to the eigenvalues
c        returned in the first column of a. These are only returned if
c        on input, ifvect was set to 1
c  nelems - the number of jacobi rotations performed by the subroutine;
c        provided solely for the amusement of a sophisticated reader
c 
c        . . . if the user so requested, initialized the matrix u
c 
        if(ifvect .eq. 0) goto 1150
c 
        do 1120 i=1,n
        do 1110 j=1,n
        uout(j,i)=0
 1110 continue
c 
        uout(i,i)=1
 1120 continue
c 
 1150 continue
c 
        maxiter=1000
        d=0
        do 1400 i=1,n
        do 1200 j=1,i
        if(abs(a(j,i)) .gt. d) d=abs(a(j,i))
 1200 continue
 1400 continue
        throut=eps*d
c 
c       one pair of off-diagonal elements after another, use
c       Jacobi rotations to diagonalize the matrix a
c 
        nelems=0
        do 4000 ijk=1,maxiter
c 
c       find the threshold after which the elements of the matrix
c       will be considered to be zero
c 
        d=0
        do 1500 i=2,n
        do 1450 j=1,i-1
        if(abs(a(j,i)) .gt. d) d=abs(a(j,i))
 1450 continue
 1500 continue
c 
cccc        thresh=d/10
        thresh=d/2
c 
        ifdone=1
        do 2800 i=2,n
        do 2600 j=1,i-1
c 
        if(abs(a(i,j)) .gt. throut) ifdone=0
        if(abs(a(j,i)) .lt. thresh) goto 2600
c 
        nelems=nelems+1
c 
        b(1,1)=a(j,j)
        b(2,2)=a(i,i)
        b(1,2)=a(j,i)
        b(2,1)=b(1,2)
c 
c       diagonalize this 2 \times 2 - matrix
c 
        call oneu(b,rlam,u)
c 
c       apply matrix u to the rows i,j of the matrix a
c 
        do 1600 k=1,n
c 
        dj=u(1,1)*a(j,k)+u(1,2)*a(i,k)
        di=u(2,1)*a(j,k)+u(2,2)*a(i,k)
c 
        a(j,k)=dj
        a(i,k)=di
 1600 continue
c 
c       apply matrix u^* to the columns i,j of the matrix a
c 
        do 1800 k=1,n
c 
        dj=u(1,1)*a(k,j)+u(1,2)*a(k,i)
        di=u(2,1)*a(k,j)+u(2,2)*a(k,i)
c 
        a(k,j)=dj
        a(k,i)=di
 1800 continue
c 
        if(ifvect .eq. 0) goto 2200
c 
c       if the user so requested, apply matrix u^* to the
c       columns i,j of the matrix u
c 
        do 2000 k=1,n
c 
        dj=u(1,1)*uout(k,j)+u(1,2)*uout(k,i)
        di=u(2,1)*uout(k,j)+u(2,2)*uout(k,i)
c 
        uout(k,j)=dj
        uout(k,i)=di
 2000 continue
c 
 2200 continue
c 
 2600 continue
 2800 continue
ccccccccccccc        call prinf('in jaceig, ijk=*',ijk,1)
c 
        do 3400 i=2,n
        do 3200 j=1,n-1
        a(i,j)=a(j,i)
 3200 continue
 3400 continue
c 
        if(ifdone .eq. 1) goto 4200
c 
 4000 continue
c 
 4200 continue
c 
c      copy the eigenvalues into the first column of a and sort them
c 
        do 4400 i=1,n
        a(i,1)=a(i,i)
 4400 continue
c 
c        put the first three largest eigenvalues in the first three
c        positions
c 
  
         do 5200 i=1,3
c 
         d=-1.0d30
         do 4600 j=i,n
c 
         if(a(j,1) .lt. d) goto 4600
         jmax=j
         d=a(j,1)
 4600 continue
c 
         if(ifvect .eq. 0) goto 5000
         do 4800 j=1,n
  
        d=uout(j,jmax)
        uout(j,jmax)=uout(j,i)
        uout(j,i)=d
 4800 continue
c 
 5000 continue
c 
        d=a(i,1)
        a(i,1)=a(jmax,1)
        a(jmax,1)=d
c 
 5200 continue
c 
c       sort the remaining eigenvalues
c 
        do 5600 i=4,n
        a(i,2)=i
 5600 continue
c 
ccccccccccccccccc        call prin2('before sorting, a(1)=*',a(1,1),n)
c 
        do 6000 i=1,n
c 
        do 5800 j=4,n-1
        if(a(j,1) .gt. a(j+1,1)) goto 5800
c 
        d=a(j,1)
        a(j,1)=a(j+1,1)
        a(j+1,1)=d
c 
        d=a(j,2)
        a(j,2)=a(j+1,2)
        a(j+1,2)=d
 5800 continue
 6000 continue
c 
c        if the user so requested, reorder the eigenvectors to corresponfd
c        to the reordered eigenvalues
c 
        if(ifvect .eq. 0) return
        do 6400 i=4,n
c 
        i1=i
        i2=a(i,2) +0.1
c 
        do 6200 j=1,n
        a(j,i1)=uout(j,i2)
 6200 continue
c 
 6400 continue
c 
        do 6800 i=4,n
        do 6600 j=1,n
        uout(j,i)=a(j,i)
 6600 continue
c 
 6800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine oneu(b,rlam,u)
        implicit real *8 (a-h,o-z)
        save
        dimension u(2,2),v(2,2),rlam(2),
     1      b(2,2),z(2,2),rlams(2,2)
        data eps/1.0d-10/,done/1.0d0/
c 
c       this subroutine produces the eigendecomposition
c       of a 2 * 2 real symmetric matrix b, so that on exit,
c 
c       a = u^* d u
c 
c       with u orthogonal, and d a diagonal matrix with the
c       vector rlam=(rlam(1),rlam(2)) on the diagonal.
c 
c       now, diagonalize the simmetrized matrix
c 
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        if(dabs(tg2phi) .lt. 1.0d-3) goto 1500
c 
        dd=dsqrt(4/tg2phi**2+4)
        if(tg2phi .lt. 0) tgphi=(-2/tg2phi-dd)/2
        if(tg2phi .gt. 0) tgphi=(-2/tg2phi+dd)/2
c 
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1500 continue
c 
        beta=tg2phi/2
        alpha=dsqrt(done-beta**2)
c 
        goto 1800
 1600 continue
c 
c       denominator is too small. use taylor expansion
c 
        alpha=dsqrt((1-den/rnum)/2)
        beta=dsqrt((1+den/rnum)/2)
 1800 continue
c 
c       construct the diagonalizing matrix
c       and the resulting diagonal
c 
        v(1,1)=alpha
        v(1,2)=beta
        v(2,1)=-beta
        v(2,2)=alpha
c 
        u(1,1)=alpha
        u(1,2)=-beta
        u(2,1)=beta
        u(2,2)=alpha
c 
c       finally, compute the resulting diagonal elements
c 
        call jac_prod2(v,b,z)
        call jac_prod2(z,u,rlams)
        rlam(1)=rlams(1,1)
        rlam(2)=rlams(2,2)
c 
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine jac_prod2(a,b,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
