subroutine mesh()
    use mod_global
    implicit none
    
    select case (selcase)
        case(1)
            filexgrid = 'xgrid_ni21nj12.dat'
            fileygrid = 'ygrid_ni21nj12.dat'
            nqsi = 21
            neta = 12
        case(2)
            filexgrid = 'xgrid_ni63nj36.dat'
            fileygrid = 'ygrid_ni63nj36.dat'
            nqsi = 63
            neta = 36
        case(3)
            filexgrid = 'xgrid_ni105nj60.dat'
            fileygrid = 'ygrid_ni105nj60.dat'
            nqsi = 105
            neta = 60
        case(4)
            filexgrid = 'xgrid_ni63nj12.dat'
            fileygrid = 'ygrid_ni63nj12.dat'
            nqsi = 63
            neta = 12
        case(5)
            filexgrid = 'xgrid_ni105nj12.dat'
            fileygrid = 'ygrid_ni105nj12.dat'
            nqsi = 105
            neta = 12
    end select
    
    imax = nqsi
    jmax = neta
        
    allocate(x(imax,jmax),y(imax,jmax))
    allocate( jac(imax,jmax))
    allocate(xksi(imax,jmax),yksi(imax,jmax))
    allocate(xeta(imax,jmax),yeta(imax,jmax))
    allocate(ksix(imax,jmax),ksiy(imax,jmax))
    allocate(etax(imax,jmax),etay(imax,jmax))
    
    open(100, file=filexgrid, form='formatted')
    open(101, file=fileygrid, form='formatted')
    do i=1,imax
        read(100,*) (x(i,j), j=1,jmax)
        read(101,*) (y(i,j), j=1,jmax)
    enddo
    
    do j=1,jmax
        do i=1,imax
            if(i==1) then
                xksi(i,j) = (x(i+1,j) - x(i,j)) !ok
                yksi(i,j) = (y(i+1,j) - y(i,j)) !ok 

            elseif(i==imax) then
                xksi(i,j) = (x(i,j) - x(i-1,j)) !ok
                yksi(i,j) = (y(i,j) - y(i-1,j)) !ok
            else
                xksi(i,j) = 0.5D0*(x(i+1,j) - x(i-1,j)) !ok
                yksi(i,j) = 0.5D0*(y(i+1,j) - y(i-1,j)) !ok
            endif

            if(j==1) then
                xeta(i,j) = (x(i,j+1) - x(i,j)) !ok
                yeta(i,j) = (y(i,j+1) - y(i,j)) !ok
            elseif(j==jmax) then
                !xeta(i,j) = 0.5D0*( x(i,j-3) - x(i,j-1))
                !yeta(i,j) = 0.5D0*(-y(i,j-3) - y(i,j-1))
                
                xeta(i,j) = (x(i,j) - x(i,j-1))
                yeta(i,j) = (y(i,j) - y(i,j-1))
            else
                xeta(i,j) = 0.5D0*(x(i,j+1) - x(i,j-1)) !ok
                yeta(i,j) = 0.5D0*(y(i,j+1) - y(i,j-1)) !ok
            endif     
        enddo
    enddo

    jac(:,:) = 1/(xksi(:,:)*yeta(:,:) - xeta(:,:)*yksi(:,:))
   
    ksix(:,:) =  jac(:,:)*yeta(:,:)
    ksiy(:,:) = -jac(:,:)*xeta(:,:)
    etax(:,:) = -jac(:,:)*yksi(:,:)
    etay(:,:) =  jac(:,:)*xksi(:,:)

    allocate(gout(imax,jmax,9))

    gout(:,:,1)  = jac(:,:)
    gout(:,:,2)  = xksi(:,:)
    gout(:,:,3)  = yksi(:,:)
    gout(:,:,4)  = xeta(:,:)
    gout(:,:,5)  = yeta(:,:)
    gout(:,:,6)  = ksix(:,:)
    gout(:,:,7)  = ksiy(:,:)
    gout(:,:,8)  = etax(:,:)
    gout(:,:,9)  = etay(:,:)

    !write mesh and metrics 
    open(102, file='plot2d.xyz', form='unformatted')
    write(102)1
    write(102)imax,jmax
    write(102)((x(i,j), i=1,imax), j=1,jmax), &
              ((y(i,j), i=1,imax), j=1,jmax)
    close(102)

    open(103, file='gout.q', form ='unformatted')
    write(103)1
    write(103)imax,jmax,9
    write(103)(((gout(i,j,l), i=1,imax), j=1,jmax), l=1,9)
    close(103)

    deallocate(gout)
    
end subroutine
    
