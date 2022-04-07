Program TISE_Euler

    ! This code numercially solves the Time Independent 1-D Schrodinger Equation Numerically.Euler method.
    ! looking at radial part
    
    implicit none

    integer,parameter :: numsteps = 1E5 !***wht is maxsteps 1E7 the largest value before i get an error? 
    double precision, dimension(0:numsteps) :: Veff, potl, centrifugal, r, U , UP, UPP
    double precision :: Sum=0.0, Rydberg=0.5, Energy, tau , hbar=1.0 , mass=1.0 , e=1.0 , k=1.0 ,Zr,lr,nr ! tau = size of increments , h = planks constant
    double precision :: percent, r0 , expect , difference , trunc
    integer :: a, i, j , n , l , nmax, lmax, zmax, iunit , iunitmax , numfiles , Z
    integer :: rr0 ! principle qn, l = anglar momentum 
    character :: namevar2*40 , namevar*40, nvar*2, lvar*2, zvar*2  !16 
    
    double precision :: Pi=3.14159265358979 ! 3238462643383279502884197169399373510
  
    !b=1000.0
    
    write(*,*)"n="
    read(*,*) n
    write(*,*)"L="
    read(*,*)  l
    write(*,*)"Z="
    read(*,*)  Z
    
    Zr=real(Z)
    lr=real(l)
    nr=real(n)
    rr0=real(r0)
    
    r0 = 10.0 * real(n)**2.0*real(z)**(-1.0) * (1.0 + (0.5 - real(l)*(real(l)+1.0)*0.5*real(n)**(-2.0)))
    !r0 = 350.00
    
    rr0=real(r0)
    
    tau = r0/numsteps !b/numsteps
    r(numsteps) = r0 !b
    U(numsteps) = 1.0D-150  !1.0! !
    U(numsteps-1) = 1.0D-150 !1.0 !
    UP(numsteps) = 0.0
    Energy=-(Rydberg*(Zr)**(2.0))/n**2

    !tau = (b - a) * numsteps**(-1.0) 
    !write(*,*) 'numsteps=' , numsteps
    !write(*,*) 'Numsteps * tau *r0 = ', tau*numsteps*r0
    write(*,*) 'Energy= ', Energy
    
    !do i = numsteps, 1 , -1 !Euler Method 
    
    do i = numsteps, 2 , -1   ! Verlet Method
    
      r(i-1)=r0-tau*(numsteps - i +1.0) ! r0-tau*(numsteps-i)
    
      !centrifugal(i-1) = (hbar**2.0*lr*(lr+1.0))/(2.0*mass*r(i-1)**2.0)
      !potl(i-1) = -(1.0*Zr*k*e**2.0)/r(i-1)  ! from int coulombs force
      !Veff(i-1) = potl(i-1) + centrifugal(i-1)
    
      if(upp(i).ne.upp(i)) then
        write(*,*) 'error, upp(i) does not equal upp(i)'
    
        write(*,*) 'veff(i)',veff(i) 
        write(*,*) 'r(i)',r(i) 
        write(*,*) 'u(i)', u(i)
        write(*,*) 'up(i)', up(i)
        write(*,*) 'hbar', hbar
        write(*,*) 'mass', mass
        write(*,*) 'Energy' ,Energy
    
        stop
        
      endif
    
     ! UP(i-1) = UP(i) - tau * UPP(i)!Eulers
     ! U(i-1) =  U(i) - tau * UP(i) !Eulers
     ! UPP(i) = (-2.0*mass/hbar**2.0)*(Energy-Veff(i))*U(i) !Eulers
    
     ! UPP(i-1) = (-2.0*mass/hbar**2.0)*(Energy-Veff(i-1))*U(i-1) ! Verlet
     ! U(i-2)=2.0*U(i-1)-U(i)+tau**(2.0)*UPP(i-1)! Verlet
    
      UPP(i-1) = ((-2.0*mass/hbar**2.0)*(Energy-((-(1.0*Zr*k*e**2.0)/r(i-1)) + &
      &((hbar**2.0*lr*(lr+1.0))/(2.0*mass*r(i-1)**2.0)))))*U(i-1) ! Verlet
    
      U(i-2)=2.0*U(i-1)-U(i)+UPP(i-1)*tau**(2.0)
    
      !r(i-1) = r(i) - tau
     
    enddo
    r(0)=0.0
    
    write(*,*) 'Pre-normalization : u(0)',u(0)
    write(*,*) 'Pre-normalization : u(r0/20)',u(rr0/20)
    write(*,*) 'Pre-normalization : u(r0/10)',u(rr0/10)
    write(*,*) 'Pre-normalization : u(r0)',u(rr0)
        
    !trunc=0.5
    !if(l>7)then
      trunc=r0/15.0
    !endif
  
    !normailization
    sum = 0.0
    do i=numsteps,0 , -2 !Integrate from 0 -> numsteps in steps of 2
      if ((l>0).and.(abs(u(i)).gt.abs(u(i+1))).and.(r(i).lt.trunc)) then
    
        write (6,*) "Truncating at r(",i,")=",r(i)
       u(i)=0.0
        u(i+1)=0.0
      else
        sum = sum + U(i)**(2.0)+4.0*U(i+1)**(2.0)+U(i+2)**(2.0)
      end if
    enddo
    sum = (tau/3.0) * sum
    
    write(*,*) 'Post-normalization : u(0)',u(0)
    write(*,*) 'Post-normalization : u(r0/20)',u(rr0/20)
    write(*,*) 'Post-normalization : u(r0/10)',u(rr0/10)
    write(*,*) 'Post-normalization : u(r0)',u(rr0)
    
    if(u(0)<0.0)then
    do i=0, numsteps, 1
     U(i)=Sum**(-0.5)*U(i)*(-1.0) ! MULTIPLIED BY -1 to flip sign! 
     !!!!U(i)=Sum**(-0.5)*U (i)/sign(one,U(1))
    enddo
    else
    u=u/sqrt(sum)
    endif 
    
    write(*,*) 'Post-normalization final ? : u(0)',u(0)
    write(*,*) 'Post-normalization final ? : u(r0/20)',u(rr0/20)
    write(*,*) 'Post-normalization final ? : u(r0/10)',u(rr0/10)
    write(*,*) 'Post-normalization final ? : u(r0)',u(rr0)
    
    !write(6,*) 'sum' , sum
    !expectation value
    Sum=0.0
    do i=0 , numsteps , 2 
     Sum = Sum + r(i)*U(i)**2.0 + 4.0*r(i)*U(i+1)**2.0 + r(i)*U(i+2)**2.0
    enddo
    Sum = (tau/3.0)*Sum
    
    expect=r0/10
    
    write(*,*) 'r0', r0
    write(*,*) '<r>_actual=', expect
    write(*,*) '<r>_code=' ,sum
    write(*,*) 'Difference' , 100.0*(sum-expect)/expect
    
    if (n >= 10) then
        write (nvar, "(I2)"),  n
      end if
    
      if (n < 10) then
        write (nvar, "(A1,I1)"), "0", n
      end if
    
      if (l >= 10) then
        write (lvar, "(I2)"),  l
      end if
          
      if (l < 10) then 
        write (lvar, "(A1,I1)"), "0", l 
      end if
    
      if (z >= 10) then
        write (zvar, "(I2)"),  z
      end if
    
      if (z < 10) then
        write (zvar, "(A1,I1)"), "0", z 
      end if
    
      write (namevar, "(A9, A2, A1, A2, A1, A2, A4)") "U_verlet_", nvar, "_", lvar, "_", zvar, ".dat" 
    write(*,*) "writing file:",namevar
      open(unit=100, file=namevar)
    
    !Writing Out U's
    !write(100,*)  'r(i)     ',  'U(i)     ' ,'centrifugal(i)     ' ,'potl(i)    ', 'UPP(i)    '
    do i = 0, numsteps
      write(100,*)  r(i), (4.0)*Pi*U(i)**(2.0)  !u(i)! ,centrifugal(i) ,potl(i), UPP(i)
    enddo
    close(100)

    end program TISE_Euler
