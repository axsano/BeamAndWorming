subroutine thomas_lu(at,bt,ct,rt,qq,n,idx)
    implicit none

    integer,intent(in) :: n,idx
    real(8),dimension(n,n,idx),intent(in) :: at,bt,ct
    real(8),dimension(n,idx),  intent(in) :: rt
    real(8),intent(out),dimension(n,idx)  :: qq

    integer i
    real(8),dimension(n,n,idx) :: xx
    real(8),dimension(n,idx)   :: yy 
    real(8),dimension(n,n)     :: ax
     
    ax = bt(:,:,1)
    call inv_gauss_jordan(ax,ax,n)
    xx(:,:,1) =    matmul(ax,ct(:,:,1))
    yy(:,1)   =    matmul(ax,rt(:,1))

    do i=2,idx
        ax = bt(:,:,i) - matmul(at(:,:,i),xx(:,:,i-1))
        call inv_gauss_jordan(ax,ax,n)           
        xx(:,:,i) = matmul(ax,ct(:,:,i))
        yy(:,i)   = matmul(ax,(rt(:,i) - matmul(at(:,:,i),yy(:,i-1))))
    enddo
  
    qq(:,idx) = yy(:,idx)

    do i=idx-1,1,-1
        qq(:,i) = yy(:,i) - matmul(xx(:,:,i),qq(:,i+1))           
    enddo

end subroutine
