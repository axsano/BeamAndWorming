program main
    use mod_global      
    implicit none

    real(kind=dp) :: t1,t2
    integer       :: nsteps

    !select the mesh
    !#1  -  nqsi = 21
    !       neta = 12
    !#2  -  nqsi = 63
    !       neta = 36 
    !#3  -  nqsi = 105
    !       neta = 60
    !#4  -  nqsi = 63
    !       neta = 12
    !#5  -  nqsi = 105
    !       neta = 12

    selcase = 1
    call mesh()
    call allocatearray()

    !initial conditions
    alpha = 0.0D0
    alpha = (alpha*pi/180.0D0)
    uin   = 0.0D0
    vin   = uin*tan(alpha)
    Ttot  = 531.2D0  !ºR
    Ptot  = 2117.0D0 !lb/ft^2

    Ttot = Ttot*5.0D0/9.0D0 !K
    Ptot = Ptot*47.880172     !Pa
    reyn = 0.0D0
    mach = 0.0D0
 
    call initialize() !qb starts here
  
    !solver
    epse           = 32.0/2
    epsi           = 3*epse   
    K2             = 0.25D0
    K4             = 1.0D0/100.0D0
    cfl            = 0.07 
    itermax        = 5000 !7500
    residualmax    = 1e-9
    nsteps         = 10
    
    allocate(residual(itermax))

    call cpu_time(t1)

    imaxr = 2
    jmaxr = 2
    lmaxr = 2

    do iter = 1,itermax
        
        call primitive()
        call fluxvectors()

        call calc_dt()

!        call calc_dissipation1()  !second order dissipation
!        call calc_dissipation2()  !second order dissipation     
        call calc_dissipation3()  !fourth order dissipation
!        call calc_dissipation4()  !pulliam 1984

!        call calc_rhs1() !residual using explicit euler
        call calc_rhs2() !residual uding block tridiagonal matrix

        call calc_residual()
!        call calc_explicit_euler()
        call calc_block_tridiagonal1()
!        call calc_block_tridiagonal2() !pulliam 1984

!        call bc_inlet1()  !first order extrapolation
!        call bc_outlet1() !first order extrapolation
!        call bc_wall1()   !first order extrapolation
        call bc_inlet2()
        call bc_outlet2()
        call bc_wall2()
        call bc_symmetry()
        
        time = iter
        if(mod(iter,nsteps)==0) then
            write(filename,"(A4,I4.4,A2)")"qout",iter,".q"
            open(103,file=filename,form ='unformatted')
            write(103)imax,jmax
            write(103)mach,alpha,reyn,time
            write(103)(((qb(i,j,l)*jac(i,j), i=1,imax), j=1,jmax), l=1,4)
            close(103)
        endif

    enddo

    call cpu_time(t2)
    
    open(104,file="residual.dat")
    do i=1,itermax
        write(104,*)i,residual(i)
    enddo

    call output()

    write(*,*)
    write(*,*)"Elapsed time:",(t2-t1)

end program

subroutine primitive()
    use mod_global
    implicit none
    
    do j=1,jmax
        do i=1,imax
            r(i,j) = qb(i,j,1)*jac(i,j)
            u(i,j) = qb(i,j,2)/qb(i,j,1)
            v(i,j) = qb(i,j,3)/qb(i,j,1)
            e(i,j) = qb(i,j,4)*jac(i,j)
            p(i,j) = (gamma-1)*(e(i,j) - 0.5D0*r(i,j)*(u(i,j)**2 + v(i,j)**2))
            a(i,j) = sqrt(gamma*p(i,j)/r(i,j))
        enddo
    enddo

    do j=1,jmax
        do i=1,imax
            uc(i,j) = ksix(i,j)*u(i,j) + ksiy(i,j)*v(i,j)
            vc(i,j) = etax(i,j)*u(i,j) + etay(i,j)*v(i,j)
        enddo
    enddo
    
end subroutine

subroutine fluxvectors()
    use mod_global
    implicit none

    do j = 1,jmax
        do i = 1,imax
            eb(i,j,1) = r(i,j)*uc(i,j)
            eb(i,j,2) = r(i,j)*u(i,j)*uc(i,j) + ksix(i,j)*p(i,j)
            eb(i,j,3) = r(i,j)*v(i,j)*uc(i,j) + ksiy(i,j)*p(i,j)
            eb(i,j,4) = (e(i,j) + p(i,j))*uc(i,j)
            
            fb(i,j,1) = r(i,j)*vc(i,j)
            fb(i,j,2) = r(i,j)*u(i,j)*vc(i,j) + etax(i,j)*p(i,j)
            fb(i,j,3) = r(i,j)*v(i,j)*vc(i,j) + etay(i,j)*p(i,j)
            fb(i,j,4) = (e(i,j) + p(i,j))*vc(i,j)
        enddo
    enddo
         
    do l = 1,4
    do j=1,jmax
        do i=1,imax
            eb(i,j,l) = eb(i,j,l)/jac(i,j)
            fb(i,j,l) = fb(i,j,l)/jac(i,j)
        end do
    end do
    enddo
            
end subroutine

subroutine calc_rhs1()
    use mod_global
    implicit none
    
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            rksi(i,j,l) = (eb(i+1,j,l) - eb(i-1,j,l))/2.0D0
            reta(i,j,l) = (fb(i,j+1,l) - fb(i,j-1,l))/2.0D0
        enddo
    enddo
    enddo
   
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            rhs(i,j,l) = rksi(i,j,l) + reta(i,j,l) &
                        -epse*(1/jac(i,j))*(dksi(i,j,l) + deta(i,j,l))
        enddo
    enddo
    enddo
 
end subroutine

subroutine calc_rhs2()
    use mod_global
    implicit none
    
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            rksi(i,j,l) = (eb(i+1,j,l) - eb(i-1,j,l))/2.0D0
            reta(i,j,l) = (fb(i,j+1,l) - fb(i,j-1,l))/2.0D0
        enddo
    enddo
    enddo
   
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            rhs(i,j,l) = -dt(i,j)*(rksi(i,j,l) + reta(i,j,l)) &
                         -dt(i,j)*epse*(1/jac(i,j))*(-dksi(i,j,l) - deta(i,j,l))
        enddo
    enddo
    enddo 
 
end subroutine

subroutine calc_dissipation1()
    use mod_global
    implicit none
    
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            dksi(i,j,l) = jac(i,j)*(+1*qb(i+1,j,l) & 
                                    -2*qb(i  ,j,l) &
                                    +1*qb(i-1,j,l))
            deta(i,j,l) = jac(i,j)*(+1*qb(i,j+1,l) &
                                    -2*qb(i,j  ,l) &
                                    +1*qb(i,j-1,l))
        enddo
    enddo
    enddo
    
end subroutine

subroutine calc_dissipation2()
    use mod_global
    implicit none
    
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            dksi(i,j,l) = (+1*jac(i+1,j)*qb(i+1,j,l) & 
                           -2*jac(i  ,j)*qb(i  ,j,l) &
                           +1*jac(i-1,j)*qb(i-1,j,l))
            deta(i,j,l) = (+1*jac(i,j+1)*qb(i,j+1,l) &
                           -2*jac(i,j  )*qb(i,j  ,l) &
                           +1*jac(i,j-1)*qb(i,j-1,l))
        enddo
    enddo
    enddo
    
end subroutine

subroutine calc_dissipation3()
    use mod_global
    implicit none
    
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            if(i==2 .or. i==imax-1) then
                dksi(i,j,l) = (+1*jac(i+1,j)*qb(i+1,j,l) & 
                               -2*jac(i  ,j)*qb(i  ,j,l) &
                               +1*jac(i-1,j)*qb(i-1,j,l))
                               
            elseif(j==2 .or. j==jmax-1) then
                deta(i,j,l) = (+1*jac(i,j+1)*qb(i,j+1,l) &
                               -2*jac(i,j  )*qb(i,j  ,l) &
                               +1*jac(i,j-1)*qb(i,j-1,l))
            
            else
                dksi(i,j,l) = -(+1*jac(i+2,j)*qb(i+2,j,l) &
                                -4*jac(i+1,j)*qb(i+1,j,l) &
                                +6*jac(i  ,j)*qb(i  ,j,l) &
                                -4*jac(i-1,j)*qb(i-1,j,l) &
                                +1*jac(i-2,j)*qb(i-2,j,l))

                deta(i,j,l) = -(+1*jac(i,j+2)*qb(i,j+2,l) &
                                -4*jac(i,j+1)*qb(i,j+1,l) &
                                +6*jac(i,j  )*qb(i,j  ,l) &
                                -4*jac(i,j-1)*qb(i,j-1,l) &
                                +1*jac(i,j-2)*qb(i,j-2,l))
            endif
        enddo
    enddo
    enddo
    
end subroutine

subroutine calc_dissipation4()
    use mod_global
    implicit none
    
    real(kind=dp) :: maux,turn
    
    do j=2,jmax-1
        do i=2,imax-1
            sig(i,j) = abs(uc(i,j)) + a(i,j)*sqrt(ksix(i,j)**2 + ksiy(i,j)**2) + &
                       abs(vc(i,j)) + a(i,j)*sqrt(etax(i,j)**2 + etay(i,j)**2)
        enddo
    enddo
 
    do j=2,jmax-1
        do i=2,imax-1
            nuksi(i,j) = abs(p(i+1,j) - 2*p(i,j) + p(i-1,j))/ &
                            (p(i+1,j) + 2*p(i,j) + p(i-1,j)) 
    
            nueta(i,j) = abs(p(i,j+1) - 2*p(i,j) + p(i,j-1))/ &
                            (p(i,j+1) + 2*p(i,j) + p(i,j-1)) 
        enddo
    enddo

    do j=2,jmax-1
        do i=2,imax-1
            eksi2(i,j) = K2*dt(i,j)*max(nuksi(i+1,j),nuksi(i,j),nuksi(i-1,j))
            eksi4(i,j) = max(0.0D0,(K4*dt(i,j) - eksi2(i,j)))

            eeta2(i,j) = K2*dt(i,j)*max(nueta(i,j+1),nueta(i,j),nueta(i,j-1))
            eeta4(i,j) = max(0.0D0,(K4*dt(i,j) - eeta2(i,j)))
        enddo
    enddo

    do j=2,jmax-1
        do i=2,imax-1
            eksii(i,j) = 3*(eksi2(i,j) + eksi4(i,j))
            eetai(i,j) = 3*(eeta2(i,j) + eeta4(i,j))
        enddo
    enddo
    
    do l=1,4
    do j=2,jmax-1
        do i=2,imax-1    
            if(i==2 .or. i==imax-1) then
            dksi(i,j,l) = (sig(i+1,j)/jac(i+1,j) - sig(i-1,j)/jac(i-1,j))* &
                           (+eksi2(i  ,j)*(+1*jac(i+1,j)*qb(i+1,j,l)  &
                                           -1*jac(i  ,j)*qb(i  ,j,l)) &
                            -eksi2(i-1,j)*(+1*jac(i  ,j)*qb(i  ,j,l)  &
                                           +1*jac(i-1,j)*qb(i-1,j,l))) 
            
            else
            dksi(i,j,l) = -(sig(i+1,j)/jac(i+1,j) - sig(i-1,j)/jac(i-1,j))* &
                           (+eksi2(i  ,j)*(+1*jac(i+1,j)*qb(i+1,j,l)  &
                                           -1*jac(i  ,j)*qb(i  ,j,l)) &
                            -eksi2(i-1,j)*(+1*jac(i  ,j)*qb(i  ,j,l)  &
                                           +1*jac(i-1,j)*qb(i-1,j,l)) &
                            -eksi4(i  ,j)*(+1*jac(i+2,j)*qb(i+2,j,l)  &
                                           -3*jac(i+1,j)*qb(i+1,j,l)  &
                                           +3*jac(i  ,j)*qb(i  ,j,l)  &
                                           -1*jac(i-1,j)*qb(i-1,j,l))*turn &
                            +eksi4(i-1,j)*(+1*jac(i+1,j)*qb(i+1,j,l)  &
                                           -3*jac(i  ,j)*qb(i  ,j,l)  &
                                           +3*jac(i-1,j)*qb(i-1,j,l)  &
                                           -1*jac(i-2,j)*qb(i-2,j,l)))
            endif
        enddo
    enddo
    enddo

    do l=1,4
    do j=2,jmax-1
        do i=2,imax-1
            if(j==2 .or. j==jmax-1) then
            deta(i,j,l) = (sig(i,j+1)/jac(i,j+1) - sig(i,j-1)/jac(i,j-1))* &
                           (+eeta2(i,j  )*(+1*jac(i,j+1)*qb(i,j+1,l)  &
                                           -1*jac(i,j  )*qb(i,j  ,l)) &
                            -eeta2(i,j-1)*(+1*jac(i,j  )*qb(i,j  ,l)  &
                                           +1*jac(i,j-1)*qb(i,j-1,l)))
            
            else
            deta(i,j,l) = -(sig(i,j+1)/jac(i,j+1) - sig(i,j-1)/jac(i,j-1))* &
                           (+eeta2(i,j  )*(+1*jac(i,j+1)*qb(i,j+1,l)  &
                                           -1*jac(i,j  )*qb(i,j  ,l)) &
                            -eeta2(i,j-1)*(+1*jac(i,j  )*qb(i,j  ,l)  &
                                           +1*jac(i,j-1)*qb(i,j-1,l)) &
                            -eeta4(i,j  )*(+1*jac(i,j+2)*qb(i,j+2,l)  &
                                           -3*jac(i,j+1)*qb(i,j+1,l)  &
                                           +3*jac(i,j  )*qb(i,j  ,l)  &
                                           -1*jac(i,j-1)*qb(i,j-1,l))*turn &
                            +eeta4(i,j-1)*(+1*jac(i,j+1)*qb(i,j+1,l)  &
                                           -3*jac(i,j  )*qb(i,j  ,l)  &
                                           +3*jac(i,j-1)*qb(i,j-1,l)  &
                                           -1*jac(i,j-2)*qb(i,j-2,l)))
            endif
        enddo
    enddo
    enddo

    do j = 2,jmax-1
        do i = 2,imax-1
            eksii(i,j) = 3*(eksi2(i,j) + eksi4(i,j))
            eetai(i,j) = 3*(eeta2(i,j) + eeta4(i,j))
        enddo
    enddo
 
end subroutine

subroutine calc_residual()
    use mod_global
    implicit none

    do l=1,4
    do j=2,jmax-1
        do i=2,imax-1
            if(abs(rhs(i,j,l)) .gt. abs(rhs(imaxr,jmaxr,lmaxr))) then
                imaxr = i
                jmaxr = j
                lmaxr = l
            endif
        enddo
    enddo
    enddo

    residual(iter) = abs(rhs(imaxr,jmaxr,lmaxr)) 

    if(mod(iter,1)==0) then
        write(*,*)
        write(*,fmt3) iter,residual(iter)
    endif

end subroutine

subroutine calc_dt()
    use mod_global
    implicit none
    
    real(kind=dp) :: aux1, aux2
    
    do j = 1,jmax
        do i = 1,imax
            aux1 = abs(uc(i,j)) + a(i,j)*sqrt(ksix(i,j)**2 + ksiy(i,j)**2)
            aux2 = abs(vc(i,j)) + a(i,j)*sqrt(etax(i,j)**2 + etay(i,j)**2)
            
            dt(i,j) = cfl/max(aux1,aux2)
        enddo
    enddo
    
end subroutine

subroutine calc_explicit_euler()
    use mod_global
    implicit none
    
    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            qb(i,j,l) = qb(i,j,l) - dt(i,j)*rhs(i,j,l)
        enddo
    enddo
    enddo
    
end subroutine

subroutine calc_block_tridiagonal1()
    use mod_global
    implicit none

    real(kind=dp), dimension(4,4,imax-2) :: AI,BI,CI
    real(kind=dp), dimension(4,imax-2)   :: DI,RI
    real(kind=dp), dimension(4,4)      :: II,Ahf,Ahb

    real(kind=dp), dimension(4,4,jmax-2) :: AJ,BJ,CJ
    real(kind=dp), dimension(4,jmax-2)   :: DJ,RJ
    real(kind=dp), dimension(4,4)      :: Bhf,Bhb

    forall(i=1:4, j=1:4) II(i,j) = (i/j)*(j/i)

    do j = 2,jmax-1
        do i = 2,imax-1
            call jacmatrix(qb(i-1,j,:),jac(i-1,j),ksix(i-1,j),ksiy(i-1,j),Ahb)
            call jacmatrix(qb(i+1,j,:),jac(i+1,j),ksix(i+1,j),ksiy(i+1,j),Ahf)

            AI(:,:,i-1) = -dt(i,j)*(0.5D0*Ahb + epsi*(1/jac(i,j))*jac(i-1,j)*II)
            BI(:,:,i-1) = +II + 2.0D0*dt(i,j)*epsi*II
            CI(:,:,i-1) = +dt(i,j)*(0.5D0*Ahf - epsi*(1/jac(i,j))*jac(i+1,j)*II)

            do l = 1,4
                DI(l,i-1) = rhs(i,j,l)
            enddo
        enddo

        call thomas_lu(AI,BI,CI,DI,RI,4,imax-2)

        do l = 1,4
        do i = 2,imax-1
            fn(i,j,l) = RI(l,i-1)
        enddo
        enddo

    enddo 

    do i = 2,imax-1
        do j = 2,jmax-1
            call jacmatrix(qb(i,j-1,:),jac(i,j-1),etax(i,j-1),etay(i,j-1),Bhb)
            call jacmatrix(qb(i,j+1,:),jac(i,j+1),etax(i,j+1),etay(i,j+1),Bhf)

            AJ(:,:,j-1) = -dt(i,j)*(0.5D0*Bhb + epsi*(1/jac(i,j))*jac(i,j-1)*II)
            BJ(:,:,j-1) = +II + 2.0D0*dt(i,j)*epsi*II
            CJ(:,:,j-1) = +dt(i,j)*(0.5D0*Bhf - epsi*(1/jac(i,j))*jac(i,j+1)*II)

            do l = 1,4
                DJ(l,j-1) = fn(i,j,l)
            enddo
        enddo

        call thomas_lu(AJ,BJ,CJ,DJ,RJ,4,jmax-2)

        do l = 1,4
        do j = 2,jmax-1
            dqb(i,j,l) = RJ(l,j-1)
        enddo
        enddo

    enddo 

    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            qb(i,j,l) = qb(i,j,l) + dqb(i,j,l)
        enddo
    enddo
    enddo
   
end subroutine

subroutine calc_block_tridiagonal2()
    use mod_global
    implicit none

    real(kind=dp), dimension(4,4,imax-2) :: AI,BI,CI
    real(kind=dp), dimension(4,imax-2)   :: DI,RI
    real(kind=dp), dimension(4,4)      :: II,Ahf,Ahb

    real(kind=dp), dimension(4,4,jmax-2) :: AJ,BJ,CJ
    real(kind=dp), dimension(4,jmax-2)   :: DJ,RJ
    real(kind=dp), dimension(4,4)      :: Bhf,Bhb

    forall(i=1:4, j=1:4) II(i,j) = (i/j)*(j/i)

    do j = 2,jmax-1
        do i = 2,imax-1
            call jacmatrix(qb(i-1,j,:),jac(i-1,j),ksix(i-1,j),ksiy(i-1,j),Ahb)
            call jacmatrix(qb(i+1,j,:),jac(i+1,j),ksix(i+1,j),ksiy(i+1,j),Ahf)

            AI(:,:,i-1) = -dt(i,j)*(0.5D0*Ahb + eksii(i,j)*(1/jac(i,j))*jac(i-1,j)*II)
            BI(:,:,i-1) = +II + 2.0D0*dt(i,j)*eksii(i,j)*II
            CI(:,:,i-1) = +dt(i,j)*(0.5D0*Ahf - eksii(i,j)*(1/jac(i,j))*jac(i+1,j)*II)

            do l = 1,4
                DI(l,i-1) = rhs(i,j,l)
            enddo
        enddo

        call thomas_lu(AI,BI,CI,DI,RI,4,imax-2)

        do l = 1,4
        do i = 2,imax-1
            fn(i,j,l) = RI(l,i-1)
        enddo
        enddo

    enddo 

    do i = 2,imax-1
        do j = 2,jmax-1
            call jacmatrix(qb(i,j-1,:),jac(i,j-1),etax(i,j-1),etay(i,j-1),Bhb)
            call jacmatrix(qb(i,j+1,:),jac(i,j+1),etax(i,j+1),etay(i,j+1),Bhf)

            AJ(:,:,j-1) = -dt(i,j)*(0.5D0*Bhb + eetai(i,j)*(1/jac(i,j))*jac(i,j-1)*II)
            BJ(:,:,j-1) = +II + 2.0D0*dt(i,j)*eetai(i,j)*II
            CJ(:,:,j-1) = +dt(i,j)*(0.5D0*Bhf - eetai(i,j)*(1/jac(i,j))*jac(i,j+1)*II)

            do l = 1,4
                DJ(l,j-1) = fn(i,j,l)
            enddo
        enddo

        call thomas_lu(AJ,BJ,CJ,DJ,RJ,4,jmax-2)

        do l = 1,4
        do j = 2,jmax-1
            dqb(i,j,l) = RJ(l,j-1)
        enddo
        enddo

    enddo 

    do l = 1,4
    do j = 2,jmax-1
        do i = 2,imax-1
            qb(i,j,l) = qb(i,j,l) + dqb(i,j,l)
        enddo
    enddo
    enddo
   
end subroutine


subroutine bc_inlet1()
    use mod_global
    implicit none
    
    real(kind=dp) :: u_i1,v_i1,p_i1,r_i1,e_i1,a_i1,t_i1,acr_i1,aux_i1
    
    do j = 2,jmax-1
        i = 1
        
        u_i1 = qb(i+1,j,2)/qb(i+1,j,1)
        v_i1 = u_i1*tan(alpha)
        
        acr_i1 = sqrt(2*gamma*(gamma - 1)/(gamma + 1)*cv*Ttot)
        
        aux_i1 = (1 - (gamma-1)/(gamma+1)*(u_i1**2 + v_i1**2)/acr_i1**2)
        
        t_i1 = (Ttot*aux_i1)/Ttot
        
        p_i1 = (Ptot*aux_i1**(gamma/(gamma-1)))/(rtot*acr_i1**2)
        
        r_i1 = p_i1/(Radm*t_i1)
        
        qb(i,j,1) = r_i1/jac(i,j)
        qb(i,j,2) = qb(i,j,1)*u_i1
        qb(i,j,3) = qb(i,j,1)*v_i1
        qb(i,j,4) = (p_i1/(gamma-1) + 0.5D0*r_i1*(u_i1**2 + v_i1**2))/jac(i,j)
    enddo
    
end subroutine

subroutine bc_inlet2()
    use mod_global
    implicit none
    
    real(kind=dp) :: u_i1,v_i1,p_i1,r_i1,e_i1,a_i1,t_i1,acr_i1,aux_i1
    real(kind=dp) :: l4,r4,dpdu_i1,du_i1,dp_i1
    
    call primitive()
    call calc_dt()

    do j = 2,jmax-1
        i = 1

        l4 = dt(i+1,j)*(u(i+1,j) - a(i+1,j))/(x(i+1,j) - x(i,j)) 
        r4 = -(l4/(1-l4))*((p(i+1,j) - p(i,j)) - (r(i+1,j)*a(i+1,j))*(u(i+1,j) - u(i,j)))
        
        dpdu_i1 = -(1+tan(alpha)**2)*(1-(gamma-1)/(gamma+1)*(1+tan(alpha)**2)*u(i,j))**(1/(gamma-1))*u(i,j)
    
        du_i1 = r4/(dpdu_i1 - r(i+1,j)*a(i+1,j))
   
        u_i1 = u(i,j) + du_i1
        v_i1 = u_i1*tan(alpha)
        
        acr_i1 = sqrt(2*gamma*(gamma - 1)/(gamma + 1)*cv*Ttot)
        
        aux_i1 = (1 - (gamma-1)/(gamma+1)*(u_i1**2 + v_i1**2)/acr_i1**2)
        
        t_i1 = (Ttot*aux_i1)/Ttot
        
        p_i1 = (Ptot*aux_i1**(gamma/(gamma-1)))/(rtot*acr_i1**2)
        
        r_i1 = p_i1/(Radm*t_i1)
        
        qb(i,j,1) = r_i1/jac(i,j)
        qb(i,j,2) = qb(i,j,1)*u_i1
        qb(i,j,3) = qb(i,j,1)*v_i1
        qb(i,j,4) = (p_i1/(gamma-1) + 0.5D0*r_i1*(u_i1**2 + v_i1**2))/jac(i,j)
    enddo
    
end subroutine

subroutine bc_wall1()
    use mod_global
    implicit none
    
    real(kind=dp) :: uw,vw,pw,rw,ew,aw,tiw,acrw,auxw,ucw,vcw
    
    do i = 1,imax
        j = 1
        
        ucw = ksix(i,j+1)*u(i,j+1) + ksiy(i,j+1)*v(i,j+1)
        vcw = etax(i,j+1)*u(i,j+1) + etay(i,j+1)*v(i,j+1)
        
        uw = +etay(i,j+1)*ucw/jac(i,j+1)
        vw = -etax(i,j+1)*ucw/jac(i,j+1)
        
        rw = qb(i,j+1,1)*jac(i,j+1)
        ew = qb(i,j+1,4)*jac(i,j+1)
        
        pw = (gamma-1)*(ew - 0.5D0*rw*(uw**2 + vw**2))

        qb(i,j,1) = rw/jac(i,j)
        qb(i,j,2) = qb(i,j,1)*uw
        qb(i,j,3) = qb(i,j,1)*vw
        qb(i,j,4) = (pw/(gamma-1) + 0.5D0*rw*(uw**2 + vw**2))/jac(i,j)
    enddo
    
end subroutine

subroutine bc_wall2()
    use mod_global
    implicit none
    
    real(kind=dp) :: uw,vw,rw,ew,aw,tw,acrw,auxw,ucw,vcw
    real(kind=dp),dimension(imax) :: aa,bb,cc,rr,pw
   
    do i = 1,imax
        j = 1

        ucw = ksix(i,j)*u(i,j) + ksiy(i,j)*v(i,j)
        vcw = etax(i,j)*u(i,j) + etay(i,j)*v(i,j)
        
        uw = +etay(i,j)*ucw/jac(i,j)
        vw = -etax(i,j)*ucw/jac(i,j)
     
        aa(i) = -(ksix(i,j)*etax(i,j) + ksiy(i,j)*etay(i,j))
        bb(i) = -3.0D0*(etax(i,j)**2 + etay(i,j)**2)
        cc(i) = +(ksix(i,j)*etax(i,j) + ksiy(i,j)*etay(i,j))

        rr(i) = -2.0D0*r(i,j)*ucw*(ksix(i,j)*etax(i,j) + ksiy(i,j)*etay(i,j)) &
                -(etax(i,j)**2 + etay(i,j)**2)*(4.0D0*p(i,j+1) - p(i,j+2))
    enddo 

    call solve_tridiag(aa,bb,cc,rr,pw,imax)
 
    do i = 1,imax
        j = 1
        
        ucw = ksix(i,j)*u(i,j) + ksiy(i,j)*v(i,j)
        vcw = etax(i,j)*u(i,j) + etay(i,j)*v(i,j)
        
        uw = +etay(i,j)*ucw/jac(i,j)
        vw = -etax(i,j)*ucw/jac(i,j)
             
        tw = (Ttot*(1 - (gamma-1)/(gamma+1)*(uw**2 + vw**2)/acr**2))/Ttot
 
        rw = pw(i)/(Radm*tw)
        
        qb(i,j,1) = rw/jac(i,j)
        qb(i,j,2) = qb(i,j,1)*uw
        qb(i,j,3) = qb(i,j,1)*vw
        qb(i,j,4) = (pw(i)/(gamma-1) + 0.5D0*rw*(uw**2 + vw**2))/jac(i,j)
    enddo
    
end subroutine
        
subroutine bc_outlet1()
    use mod_global
    implicit none
    
    real(kind=dp) :: u_imax,v_imax,p_imax,r_imax,e_imax,a_imax,m_imax,acr_imax
    
    do j = 2,jmax-1
        i = imax
        
        u_imax = qb(i-1,j,2)/qb(i-1,j,1)
        v_imax = qb(i-1,j,3)/qb(i-1,j,1)
        
        r_imax = qb(i-1,j,1)*jac(i-1,j)
        e_imax = qb(i-1,j,4)*jac(i-1,j)
        
        p_imax = (gamma-1)*(e_imax - 0.5D0*r_imax*(u_imax**2 + v_imax**2))
        a_imax = sqrt(gamma*p_imax/r_imax)
        m_imax = sqrt(u_imax**2 + v_imax**2)/a_imax
        acr_imax = sqrt(2*gamma*(gamma - 1)/(gamma + 1)*cv*Ttot)
        
        if(m_imax < 1.0D0) then
            p_imax = ((Ptot/3)*(1-(gamma-1)/(gamma+1)*(u_imax**2 + v_imax**2)/acr_imax**2)**(gamma/(gamma-1)))/(rtot*acr_imax**2)
        else
            p_imax = (gamma-1)*(e_imax - 0.5D0*r_imax*(u_imax**2 + v_imax**2))
        endif
        
        qb(i,j,1) = r_imax/jac(i,j)
        qb(i,j,2) = qb(i,j,1)*u_imax
        qb(i,j,3) = qb(i,j,1)*v_imax
        qb(i,j,4) = (p_imax/(gamma-1) + 0.5D0*r_imax*(u_imax**2 + v_imax**2))/jac(i,j)
    enddo
    
end subroutine

subroutine bc_outlet2()
    use mod_global
    implicit none
    
    real(kind=dp) :: u_imax,v_imax,p_imax,r_imax,e_imax,a_imax,m_imax,acr_imax
    real(kind=dp) :: du_imax,dv_imax,dr_imax,dp_imax
    real(kind=dp) :: l1,l3,l4 
    real(kind=dp) :: r1,r2,r3,r4   

    call primitive()
    call calc_dt()

    do j = 2,jmax-1
        i = imax

        m_imax = sqrt((qb(i-1,j,2)/qb(i-1,j,1))**2 + &
                      (qb(i-1,j,3)/qb(i-1,j,1))**2)/a(i-1,j)
        acr_imax = sqrt(2*gamma*(gamma - 1)/(gamma + 1)*cv*Ttot)

        l1 = dt(i-1,j)*u(i-1,j)/(x(i,j) - x(i-1,j))
        l3 = dt(i-1,j)*(u(i-1,j) + a(i-1,j))/(x(i,j) - x(i-1,j))
        l4 = dt(i-1,j)*(u(i-1,j) - a(i-1,j))/(x(i,j) - x(i-1,j))

        r1 = -(l1/(1+l1))*((r(i,j) - r(i-1,j)) - (1/a(i-1,j)**2)*(p(i,j) - p(i-1,j)))
        r2 = -(l1/(1+l1))*( v(i,j) - v(i-1,j))
        r3 = -(l3/(1+l3))*((p(i,j) - p(i-1,j)) + (r(i-1,j)*a(i-1,j))*(u(i,j) - u(i-1,j)))
        r4 = -(l4/(1+l4))*((p(i,j) - p(i-1,j)) - (r(i-1,j)*a(i-1,j))*(u(i,j) - u(i-1,j)))
        
        if(m_imax < 1.0D0) then
            p_imax = ((Ptot/3)*(1-(gamma-1)/(gamma+1)*(u_imax**2 + v_imax**2)/acr_imax**2)**(gamma/(gamma-1)))/(rtot*acr_imax**2)
            u_imax = u(i,j) + r3/(r(i-1,j)*a(i-1,j))
            v_imax = v(i,j) + r2
            r_imax = r(i,j) + r1
        else
            dp_imax = (r3 + r4)/2
            p_imax = p(i,j) + dp_imax
            r_imax = r(i,j) + r1 + (1/a(i-1,j)**2)*dp_imax
            u_imax = u(i,j) + (r3 - dp_imax)/(r(i-1,j)*a(i-1,j))
            v_imax = v(i,j) + r2
        endif
        
        qb(i,j,1) = r_imax/jac(i,j)
        qb(i,j,2) = qb(i,j,1)*u_imax
        qb(i,j,3) = qb(i,j,1)*v_imax
        qb(i,j,4) = (p_imax/(gamma-1) + 0.5D0*r_imax*(u_imax**2 + v_imax**2))/jac(i,j)
    enddo
    
end subroutine

subroutine bc_symmetry()
    use mod_global
    implicit none

    do i=1,imax
        j = jmax
        qb(i,j,1) =  qb(i,j-2,1) 
        qb(i,j,2) =  qb(i,j-2,2) 
        qb(i,j,3) = -qb(i,j-2,3)
        qb(i,j,4) =  qb(i,j-2,4)
    enddo

end subroutine            
            
subroutine output()
    use mod_global
    implicit none

    character(len=100), parameter :: fileplace = "/home/alex/Dropbox/CC298/proj01/beam_and_worming_v3/data/"

    call primitive()

    open(200,file='u.dat')
    open(201,file='v.dat')    
    open(202,file='r.dat')   
    open(203,file='e.dat')  
    open(204,file='p.dat')  
    open(205,file='a.dat')

    do i=1,imax
        write(200,*)(u(i,j),j=1,jmax)
        write(201,*)(v(i,j),j=1,jmax)
        write(202,*)(r(i,j),j=1,jmax)
        write(203,*)(e(i,j),j=1,jmax)
        write(204,*)(p(i,j),j=1,jmax)
        write(205,*)(a(i,j),j=1,jmax)
    enddo

    do n = 200,205
        close(n)
    enddo

end subroutine
