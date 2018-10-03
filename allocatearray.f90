subroutine allocatearray()
    use mod_global
    implicit none

    allocate(q(imax,jmax,4))  
    allocate(qb(imax,jmax,4))
    allocate(eb(imax,jmax,4))
    allocate(fb(imax,jmax,4))
    allocate(dt(imax,jmax))
    allocate(rksi(imax,jmax,4))
    allocate(reta(imax,jmax,4))
    allocate(rhs(imax,jmax,4))
    allocate(dqb(imax,jmax,4))
    allocate(fn(imax,jmax,4))
    allocate(dksi(imax,jmax,4))
    allocate(deta(imax,jmax,4))

    allocate(u(imax,jmax))
    allocate(v(imax,jmax))
    allocate(p(imax,jmax))
    allocate(r(imax,jmax))
    allocate(e(imax,jmax))
    allocate(a(imax,jmax))
    allocate(uc(imax,jmax))
    allocate(vc(imax,jmax))

    allocate(eksi2(imax,jmax),eksi4(imax,jmax)) 
    allocate(eeta2(imax,jmax),eeta4(imax,jmax)) 
    allocate(nuksi(imax,jmax),nueta(imax,jmax))  
    allocate(sig(imax,jmax))
    allocate(eksii(imax,jmax),eetai(imax,jmax))

end subroutine
