!!!************************************************************
!!!*						                	                *
!!!*                     K Epsilon 1D      	         	    *
!!!*							        	                    *
!!!*                Author: Dr. Andrea Boghi  		        *
!!!*							        	                    *
!!!*								                            *
!!!************************************************************
    
Program main_K_epsilon
    implicit none
    real*8,allocatable :: y(:),U(:),kt(:),eps(:)
    real*8,allocatable :: U0(:),kt0(:),eps0(:)
    real*8,allocatable :: nut(:),f2(:),dnutdy(:),dUdy(:)
    real*8,allocatable :: tau_mu(:),tau_R(:),Pk(:),Tk(:),Dk(:)
    real*8,allocatable :: Peps(:),Teps(:),Deps(:),epseps(:)
    integer j,ny,iter,niter
    real*8  Re_tau,sigmak,sigmae,Ce1,Ce2,Cmu,f1,alphaU,alphaKt,alphaeps
    real*8 RU,RKt,Reps,dy,aU_w,aU_e,sU,aK_w,aK_e,sK,aE_w,aE_e,sE
    logical flag

    open(1,file='imp_ke.dat')
    read(1,*) flag
    read(1,*) ny
    read(1,*) niter
    read(1,*) Re_tau
    read(1,*) alphaU
    read(1,*) alphaKt
    read(1,*) alphaeps 
    close(1)

    print*, ' niter =', niter, ' ny =', ny 
    print*, ' Re_tau =', Re_tau, ' flag =', flag 

    allocate(y(1:ny),U(1:ny),kt(1:ny),eps(1:ny),nut(1:ny),dnutdy(1:ny),dUdy(1:ny),f2(1:ny))
    allocate(U0(1:ny),kt0(1:ny),eps0(1:ny))

    call initialization(flag,ny,Re_tau,y,U,Kt,eps,dy)

    call nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu,f1)

    do iter=1,niter

        U0 = U
        Kt0 = kt
        eps0 = eps

        call nagano_takawa_k_epsilon_functions(nut,f2,ny,y,kt,eps,Cmu)
        call ddy(ny,nut,dnutdy,dy)
        call ddy(ny,U,dUdy,dy)
        U(1) = 0.d0
        Kt(1) = 0.d0
        eps(1) = 2.d0*( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2)))-dsqrt(dabs(Kt(3))))/(2.d0*dy) )**2.d0
        do j =2,ny-1
            call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),dy,Re_tau)
            U(j) =  sU + aU_e*U(j+1) + aU_w*U(j-1)
            call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),dy,sigmak)
            Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
            call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),dy,sigmae,Ce1,f1,Ce2,f2(j))
            eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
        enddo
        call u_coefficients(aU_w,aU_e,sU,nut(ny),dnutdy(ny),dy,Re_tau)
        U(ny) =  sU + (aU_e + aU_w)*U(ny-1)
        call K_coefficients(aK_w,aK_e,sK,eps(ny),nut(ny),dnutdy(ny),dUdy(ny),dy,sigmak)
        Kt(ny) = sK + (aK_e + aK_w)*Kt(ny-1)
        call E_coefficients(aE_w,aE_e,sE,eps(ny),Kt(ny),nut(ny),dnutdy(ny),dUdy(ny),dy,sigmae,Ce1,f1,Ce2,f2(ny))
        eps(ny) = sE + (aE_e + aE_w)*eps(ny-1)

        U = dabs(alphaU*U +(1.d0-alphaU)*U0)
        Kt = dabs(alphaKt*Kt +(1.d0-alphaKt)*Kt0)
        eps = dabs(alphaeps*eps +(1.d0-alphaeps)*eps0)
        print*, ' completed =', 100*real(iter)/real(niter)
    enddo

    open(1,file='point_ke.dat',form='unformatted')
    write(1) y,U,Kt,eps,nut
    close(1)

    allocate(tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny))
    allocate(Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny))

    call output_fields(ny,U,Kt,eps,nut,f1,f2,dy,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps)

    open(14,file='momentumNT.csv',form='formatted')
    write(14,*) '"y","U","k","epsilon","nut","tau_mu","tau_R","Pk","Tk","Dk","Peps","Teps","Deps","epsEps"'
    do j=1,ny
       write(14,101) y(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j),',',tau_mu(j),',',tau_R(j),',',Pk(j),',',Tk(j),',',Dk(j),',', &
       Peps(j),',',Teps(j),',',Deps(j),',',epsEps(j)
    enddo
    close(14)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10, &
    A,e18.10)

    end

!!!*************************************************
!!!*						         	             *
!!!*            initialization                        *
!!!*								                 *
!!!*************************************************

subroutine initialization(flag,ny,Re_tau,y,U,Kt,eps,dy)
    implicit none
    logical, intent(in) :: flag
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau
    real*8, intent(out) :: y(1:ny),U(1:ny),kt(1:ny),eps(1:ny),dy
    integer j
    real*8 Kappa, Cmu,nut(1:ny),y_mid

    Kappa = 4.d-1
    Cmu = 9.d-2
    y_mid = 11.635

    !initial conditions
    if (flag) then
        
        print *,'continuation'
  
        open(1,file='point_ke.dat',form='unformatted')
        read(1) y(1:ny),U(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny)
        close(1)
  
    else
  
        print *,'new run'

        dy = Re_tau/(ny-1)
  
        do j=1,ny
          y(j) = dy*(j-1)
        enddo

        do j=1,ny
            if(y(j)<=y_mid) then
                U(j) = y(j)
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y_mid)**2.d0
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu )
            else
                U(j) = (1.d0/Kappa)*dlog(y(j)) +5.5d0 
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y(j))**2.d0
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu )
            endif
        enddo
  
    endif

    end

!!!*************************************************
!!!*						         	             *
!!!*                    DDY                        *
!!!*								                 *
!!!*************************************************
    
subroutine  ddy(ny,A,DA,dy)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),dy
    real*8, intent(out) :: DA(1:ny)
    integer j

    DA(1) =  -(3.d0*A(1) -4.d0*A(2) +A(3))/(2.d0*dy)

    do j=2,ny-1
        DA(j) = (A(j+1) -A(j-1))/(2.d0*dy)
    enddo

    DA(ny) = (3.d0*A(ny) -4.d0*A(ny-1) +A(ny-2))/(2.d0*dy)

    end

!!!*************************************************
!!!*						         	             *
!!!*                    D2DY2                        *
!!!*								                 *
!!!*************************************************

subroutine  d2dy2(ny,A,D2A,dy)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),dy
    real*8, intent(out) :: D2A(1:ny)
    real*8 dy2
    integer j

    dy2 = dy*dy
    
    D2A(1) =  (12.d0*a(1) -30.d0*a(2) +24.d0*a(3) -6.d0*a(4))/(6.d0*dy2)
    
    do j=2,ny-1
        D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/dy2
    enddo
    
    D2A(ny) = (12.d0*a(ny) -30.d0*a(ny-1) +24.d0*a(ny-2) -6.d0*a(ny-3))/(6.d0*dy2)
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Nagano Takawa K - Epsilon Constants 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu,f1)
    implicit none
    real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu,f1

    sigmaK= 1.0d0
    sigmae= 1.3d0
    Ce1=1.45d0
    Ce2=1.9d0
    Cmu=0.09d0
    f1=1.d0
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Nagano Takawa K - Epsilon Functions 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  nagano_takawa_k_epsilon_functions(nut,f2,ny,y,kt,eps,Cmu)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: y(1:ny),Kt(1:ny),eps(1:ny),Cmu
    real*8, intent(out) :: nut(1:ny),f2(1:ny)
    real*8 Ret(1:ny),fmu(1:ny),dkdy(1:ny),depsdy(1:ny),dy, Ret_min, eps_min
    integer j

    eps_min = 1.d-6

    do j=1,ny
        if (eps(j) <= eps_min) then
            dy = y(2)-y(1)
            call ddy(ny,Kt,dkdy,dy)
            call ddy(ny,eps,depsdy,dy)
            Ret(j)= dabs(Kt(j)*Kt(j)/eps_min)
        else
            Ret(j)= dabs(Kt(j)*Kt(j)/eps(j))
        endif
    enddo

    Ret_min = 1.d-6
    do j=1,ny
        if (Ret(j) <= Ret_min) then
            fmu(j)= (1.d0 +4.1d0/Ret_min**0.75d0)*(1.d0 -dexp(-y(j)/26.d0))**2.d0
        else
            fmu(j)= (1.d0 +4.1d0/Ret(j)**0.75d0)*(1.d0 -dexp(-y(j)/26.d0))**2.d0
        endif
    enddo

    do j=1,ny
        nuT(j)= Cmu*fmu(j)*Ret(j)
    enddo

    do j=1,ny
        f2(j)= (1.d0 -0.3d0*dexp(-(Ret(j)/6.5d0)**2.d0))*(1.d0 -dexp(-y(j)/6.d0))**2.d0
    enddo
    
    end

!!!***************************************************
!!!*						         	               *
!!!*                U coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  u_coefficients(aU_w,aU_e,sU,nut,dnutdy,dy,Re_tau)
    implicit none
    real*8, intent(in) :: nut,dnutdy,dy,Re_tau
    real*8, intent(out) :: aU_w,aU_e,sU

    aU_w = ( 1.d0 + nut - dy*dnutdy/2.d0 )/(2.d0*(1.d0+nut))
    aU_e = ( 1.d0 + nut + dy*dnutdy/2.d0 )/(2.d0*(1.d0+nut))
    sU = (dy*dy)/(2.d0*Re_tau*(1.d0+nut))

    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdy,dUdy,dy,sigmak)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,dUdy,dy,sigmak
    real*8, intent(out) :: aK_w,aK_e,sK

    aK_w = (1.d0 + nut/sigmak - (dy/sigmak)*dnutdy/2.d0)/(2.d0*(1.d0+nut/sigmak))
    aK_e = (1.d0 + nut/sigmak + (dy/sigmak)*dnutdy/2.d0)/(2.d0*(1.d0+nut/sigmak))
    sK = (nut*dUdy*dUdy - eps)*(dy*dy)/(2.d0*(1.d0+nut/sigmak))

    end

!!!***************************************************
!!!*						         	               *
!!!*                E coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdy,dUdy,dy,sigmae,Ce1,f1,Ce2,f2)
    implicit none
    real*8, intent(in) :: eps,Kt,nut,dnutdy,dUdy,dy,sigmae,Ce1,f1,Ce2,f2
    real*8, intent(out) :: aE_w,aE_e,sE
    real*8 den
    logical method1

    method1 = .false.
    if (method1) then
        aE_w = ( 1.d0 + nut/sigmae + (dy/sigmae)*dnutdy/2.d0 )/(2.d0*(1.d0 +nut/sigmae) )
        aE_e = ( 1.d0 + nut/sigmae - (dy/sigmae)*dnutdy/2.d0 )/(2.d0*(1.d0 +nut/sigmae) )
        if (Kt==0.d0) then
            sE = 0.d0
        else
            sE = (dy*dy/(2.d0*(1.d0 +nut/sigmae)))*(Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps)*eps/Kt
        endif
    else
        den = (Kt - (dy*dy/(2.d0*(1.d0 +nut/sigmae)))*(Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps))
        aE_w = Kt*( 1.d0 + nut/sigmae + (dy/sigmae)*dnutdy/2.d0 )/(2.d0*(1.d0 +nut/sigmae) )/den
        aE_e = Kt*( 1.d0 + nut/sigmae - (dy/sigmae)*dnutdy/2.d0 )/(2.d0*(1.d0 +nut/sigmae) )/den
        sE = 0.d0 
    endif
    
    end

!!!*************************************************
!!!*						         	           *
!!!*               output fields                   *
!!!*								               *
!!!*************************************************

subroutine  output_fields(ny,U,Kt,eps,nut,f1,f2,dy,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: dy,sigmaK,sigmaE,Ce1,Ce2,f1,f2(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny)
    real*8, INTENT(OUT) :: tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny)
    real*8, INTENT(OUT) :: Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny)
    real*8 dUdy(1:ny),d2Ktdy2(1:ny),d2epsdy2(1:ny),dKtdy(1:ny),depsdy(1:ny),dnutdy(1:ny)

    call ddy(ny,U,dUdy,dy)
    call ddy(ny,nut,dnutdy,dy)

    tau_mu = dUdy
    tau_R = nut*dUdy
    Pk = tau_R*dUdy

    call d2dy2(ny,Kt,D2Ktdy2,dy)
    call ddy(ny,Kt,dKtdy,dy)

    Dk = D2Ktdy2
    Tk = (nut/sigmaK)*d2Ktdy2 + (dKtdy/sigmaK)*dnutdy

    Peps(2:ny) = f1*Ce1*(eps(2:ny)/Kt(2:ny))*Pk(2:ny)
    Peps(1) = Peps(2)

    call d2dy2(ny,eps,D2epsdy2,dy)
    call ddy(ny,eps,depsdy,dy)

    Deps = D2epsdy2
    Teps = (nut/sigmaE)*d2epsdy2 + (depsdy/sigmaE)*dnutdy
    epsEps(2:ny) = -f2(2:ny)*Ce2*(eps(2:ny)/Kt(2:ny))*eps(2:ny)
    epsEps(1) = epsEps(2)
    
    end