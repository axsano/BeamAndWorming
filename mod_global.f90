module mod_global
    implicit none

    !flow properties
    integer, parameter      :: dp    = 8
    real(kind=dp),parameter :: pi    = 4*atan(1.0D0)
    real(kind=dp),parameter :: gamma = 1.4
    real(kind=dp),parameter :: cp    = 1.0D0*1000.0D0 ![J/(kg*K)]
    real(kind=dp),parameter :: cv    = cp/gamma       ![J/(kg*K)]
    real(kind=dp),parameter :: Rspecific = cp - cv        ![J/(kg*K] 
    real(kind=dp)::Radm

    character(len= 20),parameter :: fmt1 = "(*(F12.7))"
    character(len= 20),parameter :: fmt2 = "(*(ES12.1))"
    character(len=100),parameter :: fmt3 = "(I6,20ES16.4)"
    character(len= 20),parameter :: fmt4 = "(*(F9.6))"

    !mesh properties
    integer :: i,j,k,l,n
    integer :: imax,jmax,kmax,nmax,nout(3)
    integer :: nqsi,neta
    integer :: terminalout(11)
   
    character(len=100) :: filexgrid,fileygrid,outname
 
    real(kind=dp), allocatable, dimension(:,:)   :: x,y,jac
    real(kind=dp), allocatable, dimension(:,:)   :: xksi,yksi,xeta,yeta
    real(kind=dp), allocatable, dimension(:,:)   :: ksix,ksiy,etax,etay
    real(kind=dp), allocatable, dimension(:,:,:) :: gout

    !initial conditions
    integer :: selcase
    real(kind=dp) :: uin,vin,rin,ein,pin,tin
    real(kind=dp) :: uout,vout,rout,eout,pout,tout
    real(kind=dp) :: acr,alpha,mach,reyn,mu,time
    real(kind=dp) :: aid
    real(kind=dp) :: Ptot,Ttot,rtot

    real(kind=dp) :: Ahat(4,4), Bhat(4,4) 

    !flux vectors
    real(kind=dp),allocatable,dimension(:,:,:) :: q,qb,eb,fb,dqb

    !solver parameters
    character(len=1024) :: filename
    integer :: iter,itermax,imaxr,jmaxr,lmaxr,nstep,seldissipation
    real(kind=dp) :: cfl,epse,epsi,c,K2,K4
    real(kind=dp) :: residualmax
    real(kind=dp),allocatable,dimension(:) :: residual
    real(kind=dp),allocatable,dimension(:,:) :: dt

    real(kind=dp),allocatable,dimension(:,:,:) :: rksi,reta,rhs,fn
    real(kind=dp),allocatable,dimension(:,:,:) :: deta,dksi

    real(kind=dp),allocatable,dimension(:,:) :: u,v,p,r,e,a,uc,vc

    !artificial dissipation
    real(kind=dp),allocatable,dimension(:,:) :: eksi2,eeta2,eksi4,eeta4
    real(kind=dp),allocatable,dimension(:,:) :: nuksi,nueta,sig  
    real(kind=dp),allocatable,dimension(:,:) :: eksii,eetai

end module
