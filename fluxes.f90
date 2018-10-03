subroutine fluxes(qbin,jacin,kx,ky,Xh)
    use mod_global
    implicit none

    real(kind=dp),dimension(4),intent(in)  :: qbin    
    real(kind=dp),             intent(in)  :: jacin,kx,ky
    real(kind=dp),             intent(out) :: Xh(4,4)
    
    real(kind=dp) :: phi2,theta,uaux,vaux,raux,eaux
    
    raux = qbin(1)*jacin
    uaux = qbin(2)/qbin(1)
    vaux = qbin(3)/qbin(1)
    eaux = qbin(4)*jacin

    phi2  = 0.5D0*(gamma-1)*(uaux**2 + vaux**2)
    theta = kx*uaux + ky*vaux

    Xh(1,1) = 0.0D0
    Xh(1,2) = kx
    Xh(1,3) = ky
    Xh(1,4) = 0.0D0
    
    Xh(2,1) = kx*phi2 - uaux*theta
    Xh(2,2) = theta - kx*(gamma-2)*uaux
    Xh(2,3) = ky*uaux - kx*(gamma-1)*vaux
    Xh(2,4) = kx*(gamma-1)
    
    Xh(3,1) = ky*phi2 - vaux*theta
    Xh(3,2) = kx*vaux - ky*(gamma-1)*uaux
    Xh(3,3) = theta - ky*(gamma-2)*vaux
    Xh(3,4) = ky*(gamma-1)
    
    Xh(4,1) = theta*(2*phi2 - gamma*eaux/raux)
    Xh(4,2) = kx*(gamma*eaux/raux - phi2) - (gamma-1)*uaux*theta
    Xh(4,3) = ky*(gamma*eaux/raux - phi2) - (gamma-1)*vaux*theta
    Xh(4,4) = gamma*theta
    
end subroutine
