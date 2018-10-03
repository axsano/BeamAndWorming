subroutine inv_gauss_jordan(a,ia,n)
    implicit none

    integer, intent(in) :: n
    real(8), intent(in)  :: a(n,n)
    real(8), intent(out) :: ia(n,n)
    real(8), dimension(n,n) :: ones,pe,pp,ps
    real(8), allocatable,dimension(:,:,:) :: e,p,s
    integer i, j, m, nn

    nn = (n**2-n)/2

    allocate(e(n,n,nn),p(n,n,nn),s(n,n,n))

    ones(:,:)=0.0D0
    forall(i=1:n,j=1:n) ones(i,j)=(i/j)*(j/i)

    ia(:,:)=a(:,:)

    m = 1
    do j=1,n-1
        do i=j+1,n
            e(:,:,m) = ones(:,:)
            e(i,j,m) =-ia(i,j)/ia(j,j)
            ia = matmul(e(:,:,m),ia)
            m = m + 1
        enddo
    enddo

    m = 1
    do j=n,2,-1
        do i=j-1,1,-1
            p(:,:,m) = ones(:,:)
            p(i,j,m) =-ia(i,j)/ia(j,j)
            ia = matmul(p(:,:,m),ia)
            m = m + 1
        enddo
    enddo

    do i=1,n
        s(:,:,i) = ones(:,:)
        s(i,i,i) = 1/ia(i,i)
        ia = matmul(s(:,:,i),ia)
    enddo

    pe=e(:,:,nn)
    do i=nn-1,1,-1
        pe = matmul(pe,e(:,:,i))
    enddo

    pp=p(:,:,nn)
    do i=nn-1,1,-1
        pp = matmul(pp,p(:,:,i))
    enddo

    ps=s(:,:,n)
    do i=n-1,1,-1
        ps = matmul(ps,s(:,:,i))
    enddo
   
    ia = matmul(ps,matmul(pp,pe)) 

    deallocate(e,p,s)

end subroutine
