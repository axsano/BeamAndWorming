subroutine initialize()
    use mod_global
    implicit none
    
    real(kind=dp) :: aux
    
    rtot = Ptot/(Rspecific*Ttot)
    aid = sqrt(gamma*Ptot/rtot)

    acr = sqrt(2*gamma*(gamma-1)/(gamma+1)*cv*Ttot)
    Radm = Rspecific*Ttot/acr**2

    aux = (1-(gamma-1)/(gamma+1)*(1+tan(alpha)**2)*(uin/acr)**2) 

    tin = (Ttot*aux)/Ttot
    pin = (Ptot*aux**(gamma/(gamma-1)))/(rtot*acr**2)
    !rin = pin/(Radm*tin)
    rin = (rtot*aux**(1/(gamma-1)))/rtot
    ein = (ptot/(gamma-1) + 0.5D0*rtot*(uin**2 + vin**2))/(rtot*acr**2)

    uout = 0.0D0
    vout = 0.0D0
    tout = tin
    pout = pin/3.0D0
    rout = rin !check this boundary condition
    eout = ((ptot/3)/(gamma-1) + 0.5D0*rtot*(uin**2 + vin**2))/(rtot*acr**2)

    write(*,*)
    write(*,*)"Initial conditions:"
    write(*,*)
    write(*,900)"uin = ",uin," [m/s]"
    write(*,900)"vin = ",vin," [m/s]"
    write(*,900)"alpha = ",alpha," [rad]"
    write(*,900)"Ptot = ",Ptot," [Pa]"
    write(*,900)"Ttot = ",Ttot," [K]"
    write(*,900)"rtot = ",rtot," [kg/m^3]"
    write(*,900)"acr = ",acr," [m/s]"
    write(*,900)"aid = ",aid," [m/s]"
    write(*,900)
    write(*,*)"Adimensional initial conditions:"
    write(*,*)
    write(*,*)"Flow properties:"
    write(*,900)"Radm = ",Radm,""
    write(*,*)
    write(*,*)"Inlet:"
    write(*,900)"Pin = ",pin,""
    write(*,900)"Tin = ",tin,""
    write(*,900)"rin = ",rin,""
    write(*,900)"ein = ",ein,""
    write(*,*)
    write(*,*)"Outlet:"
    write(*,900)"Pout = ",pout,""
    write(*,900)"Tout = ",tout,""
    write(*,900)"rout = ",rout,""
    write(*,900)"eout = ",eout,""
    write(*,*)
   
    do j=1,jmax
        do i=1,imax
            if(i==imax) then
                q(i,j,1) = rout
                q(i,j,2) = 0.0D0
                q(i,j,3) = 0.0D0
                q(i,j,4) = eout 
            else
                q(i,j,1) = rin
                q(i,j,2) = 0.0D0
                q(i,j,3) = 0.0D0
                q(i,j,4) = ein
            endif
        enddo
    enddo
    
    do l=1,4    
    do j=1,jmax
        do i=1,imax
            qb(i,j,l) = q(i,j,l)/jac(i,j)
        enddo
    enddo
    enddo
    
    time = 0
    
    write(filename,"(A4,I4.4,A2)")"qout",iter,".q"
    open(103,file=filename,form ='unformatted')
    write(103)imax,jmax
    write(103)mach,alpha,reyn,time
    write(103)(((qb(i,j,l)*jac(i,j), i=1,imax), j=1,jmax), l=1,4)
    close(103)
    
    900 format(A10,F12.4,A)
    
end subroutine
