!!!************************************************************
!!!*						                	                *
!!!*                     Thermal 1D      	         	    *
!!!*							        	                    *
!!!*                Author: Dr. Andrea Boghi  		        *
!!!*							        	                    *
!!!*								                            *
!!!************************************************************
    
Program main_theraml
    implicit none
    real*8,allocatable :: y(:),U(:),kt(:),eps(:)
    real*8,allocatable :: T(:),Th2(:),T0(:),Th20(:)
    real*8,allocatable :: nut(:),f2(:),dnutdy(:),dUdy(:)
    integer j,ny,iter,niter
    real*8  Re_tau,sigmaT,sigmaTh2,Ce1,Ce2,alphaT,alphaTh2
    real*8 resT,resTh2,dy,aT_w,aT_e,sT,aTh2_w,aTh2_e,sTh2, conv_fac
    logical flag

    open(1,file='imp_ther.dat')
    read(1,*) flag
    read(1,*) ny
    read(1,*) niter
    read(1,*) Re_tau
    read(1,*) alphaT
    read(1,*) alphaTh2
    close(1)

    print*, ' niter =', niter, ' ny =', ny 
    print*, ' Re_tau =', Re_tau, ' flag =', flag 

    allocate(y(1:ny),U(1:ny),kt(1:ny),eps(1:ny),nut(1:ny),dnutdy(1:ny),dUdy(1:ny),f2(1:ny))
    allocate(T(1:ny),Th2(1:ny),T0(1:ny),Th20(1:ny))

    open(1,file='point_ke.dat',form='unformatted')
    read(1) y(1:ny),U(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),f2(1:ny)
    close(1)

    open(14,file='momentum_test.csv',form='formatted')
    write(14,*) '"y","U","k","epsilon","nut"'
    do j=1,ny
       write(14,101) y(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j)
    enddo
    close(14)

    goto 1001

    conv_fac = 1.d0
    open(11,file='residuals_thermal.csv')
    write(11,*) '"iter","resT","resTh2"'

    do iter=1,niter

        T0 = T
        Th20 = Th2

        call ddy(ny,nut,dnutdy,dy)
        call ddy(ny,U,dUdy,dy)
        T(1) = 0.d0
        Th2(1) = 0.d0
        do j =2,ny-1
            call T_coefficients(aT_w,aT_e,sT,nut(j),dnutdy(j),dy,Re_tau,sigmaT)
            T(j) =  sT + aT_e*T(j+1) + aT_w*T(j-1)
            call Th2_coefficients(aTh2_w,aTh2_e,sTh2,eps(j),nut(j),dnutdy(j),dUdy(j),dy,sigmaTh2)
            Th2(j) = sTh2 + aTh2_e*Th2(j+1) + aTh2_w*Th2(j-1)
        enddo
        call T_coefficients(aT_w,aT_e,sT,nut(ny),dnutdy(ny),dy,Re_tau,sigmaT)
        T(ny) =  sT + (aT_e + aT_w)*T(ny-1)
        call Th2_coefficients(aTh2_w,aTh2_e,sTh2,eps(ny),nut(ny),dnutdy(ny),dUdy(ny),dy,sigmaTh2)
        Th2(ny) = sTh2 + (aTh2_e + aTh2_w)*Th2(ny-1)

        call residuals(ny,T,T0,resT)
        call residuals(ny,Th2,Th20,resTh2)

        write(11,102) conv_fac*iter,',',resT,',',resTh2

        T = dabs(alphaT*T +(1.d0-alphaT)*T0)
        Th2 = dabs(alphaTh2*Th2 +(1.d0-alphaTh2)*Th20)
        print*, ' completed =', 100*real(iter)/real(niter), ' resT = ', resT, ' resTh2 = ', resTh2
        
    enddo
    close(11)

    open(1,file='point_thermal.dat',form='unformatted')
    write(1) y,T,Th2
    close(1)

    open(14,file='momentum_test.csv',form='formatted')
    write(14,*) '"y","U","k","epsilon","nut"'
    do j=1,ny
       write(14,101) y(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j)
    enddo
    close(14)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    102 format(e18.10,A,e18.10,A,e18.10)

    1001 continue

    end
