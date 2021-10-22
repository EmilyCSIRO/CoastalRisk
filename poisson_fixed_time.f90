SUBROUTINE poisson_fixed_time (time, MC, fldP, timing)
   !Generates stochastic flood events
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: time, MC
   REAL, INTENT(OUT) :: fldP (6,time)
   INTEGER, INTENT(OUT) :: timing (time)
   REAL :: event_time(time), randomnum, t, dt, floodsize, rantemp
   REAL :: lambda, aParam, bParam, aParamP(6), bParamP(6), flood(time)
   INTEGER :: n, seed, i, ecount,pcount
   
   INTEGER :: usechangeRand=0  !change to 1 if using changing random number
   
   !variables for random number function
   REAL :: rand                 ! Declare the type of the rand() function
   INTEGER :: k,j, id, S                 ! Counts random numbers
   INTEGER, DIMENSION(:), ALLOCATABLE :: S_ARRAY
   integer*4 timeArray(3)    ! Holds the hour, minute, and second 
   REAL :: random_beta
   LOGICAL :: FIRST
	  
    !Assign parameter values	  
     lambda=1.0/27.0   !Flood frequency assumed to be the same across the coastal zone, but with varying magnitude 
     aParamP=[0.9,0.8,1.0,0.8,0.5,0.7]
     bParamP=[2.9,8.4,1.8,8.3,7.8,8.0]
      
     n = 0
     t = 0.0
     dt=0.0
	 
	 DO i=1, time
	   event_time(i)=0.0
	 END DO
	 
	 DO i=1, time
	   flood(i) = 0.0
	   timing(i)=0
	 END DO

	 IF (usechangeRand==0) THEN          !use a constant random number generator which stays the same each time it is run   
	   !Poisson random number generator using constant value to determine flood timing
	   CALL random_seed(seed)
	   DO WHILE (t < time)  
         !call function for random number generator below which changes each time
         CALL random_number(randomnum)
	     dt=-log(randomnum)/lambda
	     n=n+1
	     event_time(n)=t   
	     t=t+dt
	   END DO	 
	   event_time=NINT(event_time)
	 
	  !Random number generator using beta distribution to simulate flood magnitude
       !assign parameters based on polder
       DO pcount=1,6   !polder 1 to 6	 
	     aParam=aParamP(pcount)
	     bParam=bParamP(pcount)	  

	     floodsize=0.0
	     ecount=0
	
	     DO i=1, time
	       IF(event_time(i)>0) THEN
	         n=NINT(event_time(i))
		     ecount=ecount+1   !need to initialise the beta distribution routine at the start of each round
	         FIRST = ecount .eq. 1	 
		     floodsize=random_beta(aParam,bParam,FIRST)
		     flood(n)=floodsize                  
		     timing(n)=n                             
	       END IF
	     END DO
	     fldP(pcount,:)=flood(:)

	   END DO  !polder 1 to 6 	 
	 END IF
	 
	 IF (usechangeRand==1) THEN   !Use a changing random number generation
	   !Poisson distribution for flood timing
	   IF (MC+9999<0.001) THEN   !not being run using Monte Carlo, use time stamp for random number generator
	     call itime(timeArray)     ! Get the current time
         k =   timeArray(1)+timeArray(2)+timeArray(3)
	   ELSE
	     k=MC
	   END IF
       randomnum = rand(k)   !this initialises the seed
	   DO WHILE (t < time)  
	     randomnum=rand()             !call function for random number generator for each iteration using the seed generated earlier
		 dt=-log(randomnum)/lambda
	     n=n+1
	     event_time(n)=t   
	     t=t+dt
	   END DO	 
	   event_time=NINT(event_time)	 
	 
      !Beta distribution for flood magnitude   NOTE this generates the same random number each round. Tried to implement the time-based seed but wasn't successful
	
       DO pcount=1,6   !polder 1 to 6	 
	     aParam=aParamP(pcount)
	     bParam=bParamP(pcount)	  

   	     floodsize=0.0
	     ecount=0
	
	     DO i=1, time
	       IF(event_time(i)>0) THEN
	         n=NINT(event_time(i))
		     ecount=ecount+1   !need to initialise the beta distribution routine at the start of each round
	         FIRST = ecount .eq. 1
		 		 
	         floodsize=random_beta(aParam,bParam,FIRST)
		 
		     flood(n)=floodsize                  
		     timing(n)=n                             
	       END IF
	     END DO
	 	 fldP(pcount,:)=flood(:)
	   END DO  ! End do polder 1 to 6 
	 	 
	 END IF
	 
 END SUBROUTINE poisson_fixed_time
 


 !*******************BETA DISTRIBUTION RANDOM NUMBER GENERATOR************
 !From: http://users.cla.umn.edu/~erm/data/qr09/codes/random.f90
 !Ellen McGrattan, Department of Economics, University of Minnesota

real FUNCTION random_beta(aa, bb, first)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,1]
! FROM A BETA DISTRIBUTION WITH DENSITY
! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
! USING CHENG'S LOG LOGISTIC METHOD.

!     AA = SHAPE PARAMETER FROM DISTRIBUTION (0 < real)
!     BB = SHAPE PARAMETER FROM DISTRIBUTION (0 < real)

IMPLICIT NONE
real, INTENT(IN)    :: aa, bb
LOGICAL, INTENT(IN) :: first


!     Local variables
real, PARAMETER     :: aln4 = 1.3862944
real                :: a, b, d, f, h, t, c, g, r, s, x, y, z
LOGICAL             :: swap

real :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                             vsmall = TINY(1.0), vlarge = HUGE(1.0)
							 
							 
							

SAVE d, f, h, t, c, swap

IF (aa <= zero .OR. bb <= zero) THEN
  PRINT *, 'IMPERMISSIBLE SHAPE PARAMETER VALUE(S)'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = aa
  b = bb
  swap = b > a
  IF (swap) THEN
    g = b
    b = a
    a = g
  END IF
  d = a/b
  f = a+b
  IF (b > one) THEN
    h = SQRT((two*a*b - f)/(f - two))
    t = one
  ELSE
    h = b
    t = one/(one + (a/(vlarge*b))**b)
  END IF
  c = a+h
END IF


DO
     

  CALL RANDOM_NUMBER(r)   
  CALL RANDOM_NUMBER(x)   
  
 
  s = r*r*x
  
  IF (r < vsmall .OR. s <= zero) CYCLE
  IF (r < t) THEN
    x = LOG(r/(one - r))/h
    y = d*EXP(x)
    z = c*x + f*LOG((one + d)/(one + y)) - aln4
    IF (s - one > z) THEN
      IF (s - s*z > one) CYCLE
      IF (LOG(s) > z) CYCLE
    END IF
    random_beta = y/(one + y)
  ELSE
    IF (4.0*s > (one + one/d)**f) CYCLE
    random_beta = one
  END IF
  EXIT
END DO

!WRITE(*,*) random_beta

  IF (swap) random_beta = one - random_beta
  RETURN
END FUNCTION random_beta

 