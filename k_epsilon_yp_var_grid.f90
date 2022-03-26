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
    real*8,allocatable :: y(:),U(:),kt(:),eps(:),detady(:),d2etady2(:)
    real*8,allocatable :: U0(:),kt0(:),eps0(:)
    real*8,allocatable :: nut(:),f2(:),dnutdy(:),dUdy(:),dnutdeta(:),dUdeta(:)
    real*8,allocatable :: tau_mu(:),tau_R(:),Pk(:),Tk(:),Dk(:)
    real*8,allocatable :: Peps(:),Teps(:),Deps(:),epseps(:)
    integer j,ny,iter,niter
    real*8  Re_tau,dy_min,sigmak,sigmae,Ce1,Ce2,Cmu,f1,alphaU,alphaKt,alphaeps
    real*8 resU,resK,resE,deta,aU_w,aU_e,sU,aK_w,aK_e,sK,aE_w,aE_e,sE, conv_fac
    logical flag

    open(1,file='imp_ke_var.dat')
    read(1,*) flag
    read(1,*) ny
    read(1,*) niter
    read(1,*) Re_tau
    read(1,*) dy_min
    read(1,*) alphaU
    read(1,*) alphaKt
    read(1,*) alphaeps 
    close(1)

    print*, ' niter =', niter, ' ny =', ny 
    print*, ' Re_tau =', Re_tau, ' flag =', flag 

    allocate(y(1:ny),U(1:ny),kt(1:ny),eps(1:ny),nut(1:ny),dnutdy(1:ny),dUdy(1:ny),f2(1:ny))
    allocate(U0(1:ny),kt0(1:ny),eps0(1:ny),detady(1:ny),d2etady2(1:ny),dnutdeta(1:ny),dUdeta(1:ny))

    call initialization(flag,ny,Re_tau,dy_min,y,detady,d2etady2,U,Kt,eps,deta)

    call nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu,f1)

    conv_fac = 1.d0
    open(11,file='residuals_momentum_var.csv')
    write(11,*) '"iter","resU","resK","resE"'

    do iter=1,niter

        U0 = U
        Kt0 = kt
        eps0 = eps

        call nagano_takawa_k_epsilon_functions(nut,f2,ny,y,kt,eps,Cmu)
        call ddeta(ny,nut,dnutdeta,deta)
        dnutdy = dnutdeta*detady
        call ddeta(ny,U,dUdeta,deta)
        dUdy = dUdeta*detady
        U(1) = 0.d0
        !call u_coefficients(aU_w,aU_e,sU,nut(1),dnutdy(1),deta,Re_tau,d2etady2(1),detady(1))
        !U(1) =  sU + (aU_e + aU_w)*U(2) + 2.d0*aU_w*dy 
        Kt(1) = 0.d0
        eps(1) = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2)))-dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
        do j =2,ny-1
            call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j))
            U(j) =  sU + aU_e*U(j+1) + aU_w*U(j-1)
            call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),deta,sigmak,d2etady2(j),detady(j))
            Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
            call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae,Ce1,f1,Ce2,f2(j),d2etady2(j), &
            detady(j))
            eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
        enddo
        call u_coefficients(aU_w,aU_e,sU,nut(ny),dnutdy(ny),deta,Re_tau,d2etady2(ny),detady(ny))
        U(ny) =  sU + (aU_e + aU_w)*U(ny-1)
        call K_coefficients(aK_w,aK_e,sK,eps(ny),nut(ny),dnutdy(ny),dUdy(ny),deta,sigmak,d2etady2(ny),detady(ny))
        Kt(ny) = sK + (aK_e + aK_w)*Kt(ny-1)
        call E_coefficients(aE_w,aE_e,sE,eps(ny),Kt(ny),nut(ny),dnutdy(ny),dUdy(ny),deta,sigmae,Ce1,f1,Ce2,f2(ny),d2etady2(ny), &
        detady(ny))
        eps(ny) = sE + (aE_e + aE_w)*eps(ny-1)

        call residuals(ny,U,U0,resU)
        call residuals(ny,Kt,Kt0,resK)
        call residuals(ny,eps,eps0,resE)
        write(11,102) conv_fac*iter,',',resU,',',resK,',',resE

        U = dabs(alphaU*U +(1.d0-alphaU)*U0)
        Kt = dabs(alphaKt*Kt +(1.d0-alphaKt)*Kt0)
        eps = dabs(alphaeps*eps +(1.d0-alphaeps)*eps0)
        print*, ' completed =', 100*real(iter)/real(niter), ' resU = ', resU, ' resK = ', resK, ' resE = ', resE
        
    enddo
    close(11)

    open(1,file='point_ke_var.dat',form='unformatted')
    write(1) y,detady,d2etady2,U,Kt,eps,nut,f2
    close(1)

    allocate(tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny))
    allocate(Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny))

    call output_fields(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps,detady,d2etady2)

    open(14,file='momentumNT_var.csv',form='formatted')
    write(14,*) '"y","U","k","epsilon","nut","tau_mu","tau_R","Pk","Tk","Dk","Peps","Teps","Deps","epsEps"'
    do j=1,ny
       write(14,101) y(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j),',',tau_mu(j),',',tau_R(j),',',Pk(j),',',Tk(j),',',Dk(j),',', &
       Peps(j),',',Teps(j),',',Deps(j),',',epsEps(j)
    enddo
    close(14)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10, &
    A,e18.10)

    102 format(e18.10,A,e18.10,A,e18.10,A,e18.10)

    end

!!!*************************************************
!!!*						         	             *
!!!*            initialization                        *
!!!*								                 *
!!!*************************************************

subroutine initialization(flag,ny,Re_tau,dy_min,y,detady,d2etady2,U,Kt,eps,deta)
    implicit none
    logical, intent(in) :: flag
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),U(1:ny),kt(1:ny),eps(1:ny),deta
    integer j
    real*8 Kappa, Cmu,nut(1:ny),f2(1:ny),y_mid

    Kappa = 4.d-1
    Cmu = 9.d-2
    y_mid = 11.635

    !initial conditions
    if (flag) then
        
        print *,'continuation'
  
        open(1,file='point_ke_var.dat',form='unformatted')
        read(1) y(1:ny),detady(1:ny),d2etady2(1:ny),U(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),f2(1:ny)
        close(1)
  
    else
  
        print *,'new run'

        call grid(ny,dy_min,Re_tau,y,detady,d2etady2,deta)

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
!!!*                    ddeta                        *
!!!*								                 *
!!!*************************************************
    
subroutine  ddeta(ny,A,DA,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),deta
    real*8, intent(out) :: DA(1:ny)
    integer j

    DA(1) =  -(3.d0*A(1) -4.d0*A(2) +A(3))/(2.d0*deta)

    do j=2,ny-1
        DA(j) = (A(j+1) -A(j-1))/(2.d0*deta)
    enddo

    DA(ny) = (3.d0*A(ny) -4.d0*A(ny-1) +A(ny-2))/(2.d0*deta)

    end

!!!*************************************************
!!!*						         	             *
!!!*                    d2deta2                        *
!!!*								                 *
!!!*************************************************

subroutine  d2deta2(ny,A,D2A,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),deta
    real*8, intent(out) :: D2A(1:ny)
    real*8 deta2
    integer j

    deta2 = deta*deta
    
    D2A(1) =  (12.d0*a(1) -30.d0*a(2) +24.d0*a(3) -6.d0*a(4))/(6.d0*deta2)
    
    do j=2,ny-1
        D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/deta2
    enddo
    
    D2A(ny) = (12.d0*a(ny) -30.d0*a(ny-1) +24.d0*a(ny-2) -6.d0*a(ny-3))/(6.d0*deta2)
    
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
    real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
    integer j

    eps_min = 1.d-12

    do j=1,ny
        if (eps(j) <= eps_min) then
            Ret(j)= dabs(Kt(j)*Kt(j)/eps_min)
        else
            Ret(j)= dabs(Kt(j)*Kt(j)/eps(j))
        endif
    enddo

    Ret_min = 1.d-12
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
subroutine  u_coefficients(aU_w,aU_e,sU,nut,dnutdy,deta,Re_tau,d2etady2,detady)
    implicit none
    real*8, intent(in) :: nut,dnutdy,deta,Re_tau,d2etady2,detady
    real*8, intent(out) :: aU_w,aU_e,sU
    real*8 dev

    dev = deta*( (1.d0+nut)*d2etady2 + dnutdy*detady )/(4.d0*(1.d0+nut)*(detady)**2.d0)

    aU_w = 5.d-1 - dev
    aU_e = 5.d-1 + dev
    sU = (deta*deta)/(2.d0*Re_tau*(1.d0+nut)*(detady)**2.d0)

    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdy,dUdy,deta,sigmak,d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,dUdy,deta,sigmak,d2etady2,detady
    real*8, intent(out) :: aK_w,aK_e,sK
    real*8 dev

    dev = deta*( (sigmak+nut)*d2etady2 + dnutdy*detady )/(4.d0*(sigmak+nut)*detady**2.d0)

    aK_w = 5.d-1 - dev
    aK_e = 5.d-1 + dev
    sK = (nut*dUdy*dUdy - eps)*(deta*deta)/(2.d0*(1.d0+nut/sigmak)*detady**2.d0)

    end

!!!***************************************************
!!!*						         	               *
!!!*                E coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdy,dUdy,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,Kt,nut,dnutdy,dUdy,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady
    real*8, intent(out) :: aE_w,aE_e,sE
    real*8 K_min, Kb, dev
    logical method1

    dev = deta*( (sigmae+nut)*d2etady2 + dnutdy*detady )/(4.d0*(sigmae+nut)*detady**2.d0)

    K_min = 1.d-12

    Kb = (Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps)*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detady**2.d0))

    method1 = .true.
    if (method1) then
        aE_w = 5.d-1 - dev
        aE_e = 5.d-1 + dev
        if (Kt<=K_min) then
            sE = 0.d0*Kb*eps/K_min
        else
            sE = Kb*eps/Kt
        endif
    else
        aE_w = (5.d-1 - dev)/(1.d0 - Kb/Kt)
        aE_e = (5.d-1 + dev)/(1.d0 - Kb/Kt)
        sE = 0.d0 
    endif
    
    end

!!!*************************************************
!!!*						         	           *
!!!*               output fields                   *
!!!*								               *
!!!*************************************************

subroutine  output_fields(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps,detady, &
    d2etady2)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: deta,sigmaK,sigmaE,Ce1,Ce2,f1,f2(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny),detady(1:ny),d2etady2(1:ny)
    real*8, INTENT(OUT) :: tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny)
    real*8, INTENT(OUT) :: Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny)
    real*8 dUdy(1:ny),d2Ktdeta2(1:ny),d2epsdeta2(1:ny),dKtdeta(1:ny),depsdeta(1:ny),dnutdy(1:ny),dUdeta(1:ny),dnutdeta(1:ny)

    call ddeta(ny,U,dUdeta,deta)
    dUdy = dUdeta*detady
    call ddeta(ny,nut,dnutdeta,deta)
    dnutdy = dnutdeta*detady

    tau_mu = dUdy
    tau_R = nut*dUdy
    Pk = tau_R*dUdy

    call d2deta2(ny,Kt,d2Ktdeta2,deta)
    call ddeta(ny,Kt,dKtdeta,deta)

    Dk = d2Ktdeta2*detady**2.d0 + dKtdeta*d2etady2
    Tk = (nut/sigmaK)*(d2Ktdeta2*detady**2.d0 + dKtdeta*d2etady2) + (dKtdeta*detady/sigmaK)*dnutdy

    Peps(2:ny) = f1*Ce1*(eps(2:ny)/Kt(2:ny))*Pk(2:ny)
    Peps(1) = Peps(2)

    call d2deta2(ny,eps,D2epsdeta2,deta)
    call ddeta(ny,eps,depsdeta,deta)

    Deps = d2epsdeta2*detady**2.d0 + depsdeta*d2etady2
    Teps = (nut/sigmaE)*(d2epsdeta2*detady**2.d0 + depsdeta*d2etady2) + (depsdeta*detady/sigmaE)*dnutdy
    epsEps(2:ny) = -f2(2:ny)*Ce2*(eps(2:ny)/Kt(2:ny))*eps(2:ny)
    epsEps(1) = epsEps(2)
    
    end

!!!*************************************************
!!!*						         	           *
!!!*               residuals                   *
!!!*								               *
!!!*************************************************

subroutine  residuals(ny,A,A0,resA)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),A0(1:ny)
    real*8, INTENT(OUT) :: resA
    real*8 sumN, sumD
    integer j

    sumN = 0
    do j=1,ny
        sumN = sumN +dabs(A(j)- A0(j))
    enddo

    sumD = 0
    do j=1,ny
        sumD = sumD + dabs(A0(j))
    enddo

    resA = sumN/sumD

    end

!!!*************************************************
!!!*						         	             *
!!!*            grid                       *
!!!*								                 *
!!!*************************************************

subroutine grid(ny,dy_min,Re_tau,y,detady,d2etady2,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),deta
    integer j
    real*8 a, b, eta

    deta = Re_tau/(ny-1)

    a = (deta - dy_min)/(Re_tau*deta - deta*deta)
    b = (dy_min*Re_tau - deta*deta)/(Re_tau*deta - deta*deta)

    do j=1,ny
        eta = deta*(j-1)
        y(j) = a*eta**2.d0 + b*eta
        detady(j) = 1.d0/(2.d0*a*eta + b)
        d2etady2(j) = -2.d0*a/(2.d0*a*eta + b)**3.d0
    enddo

    print*, ' dy_max =', y(ny)-y(ny-1), ' dy_min =', y(2)-y(1), ' ratio =', (y(ny)-y(ny-1))/(y(2)-y(1))

    end
