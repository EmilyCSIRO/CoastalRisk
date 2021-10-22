SUBROUTINE crop_matrix_module(yld, nyear, numcrop, optim, cropYield)
  !Calculates matrix of crop yields for entire simulation based on cropping patterns.
  IMPLICIT NONE

  REAL, INTENT(IN):: yld(numcrop)
  INTEGER, INTENT(IN) :: nyear, numcrop, optim
  REAL, INTENT(OUT) :: cropYield(numcrop,12*nyear)
  
  INTEGER ::  k,i,t,n
  REAL, ALLOCATABLE :: crop_matrix(:,:), yield_matrix(:,:) 
  INTEGER, ALLOCATABLE :: cpmonth_st(:), cpmonth_end(:), crop_calendar(:,:) 
  
  ALLOCATE (crop_matrix(numcrop,12*nyear),yield_matrix(numcrop,12*nyear))
  ALLOCATE (cpmonth_st(numcrop), cpmonth_end(numcrop),crop_calendar(numcrop,12))

  !crops: local aus, hyv aus, local aman; hyv aman, local boro, hybrid boro, hyv boro, wheat, maize, jute, water gourd, winter pumpkin, summer pumpkin 
  	cpmonth_st = (/4,4,8,8,11,1,1,11,11,5,7,11,3/)
    cpmonth_end = (/8,8,12,12,5,5,5,4,5,9,3, 3, 11/)

  ! Generate crop calendar showing 1 when crops are growing and -9999 when they aren't 
  DO i=1,numcrop
    DO t=1,12
	  IF (cpmonth_st(i) < cpmonth_end(i)) THEN
	    IF (cpmonth_st(i)<=t .AND. t<=cpmonth_end(i)) THEN
          crop_calendar(i,t)=1
		ELSE
		  crop_calendar(i,t)=-9999
        END IF	
	  ELSE IF (t>=cpmonth_st(i) .OR. t<=cpmonth_end(i)) THEN  
	    crop_calendar(i,t)=1
      ELSE
        crop_calendar(i,t)=-9999	  
	  END IF  
	END DO
  END DO
  IF (optim==2) THEN 
    OPEN(unit=55, file='Outputs/CropCalendar.txt', status='unknown')   !*turn off for optimisation 
    DO i=1, numcrop
      WRITE(55,*) (crop_calendar(i,t),t=1,12)
    END DO
    CLOSE(unit=55)
  ENDIF
    
  yield_matrix = SPREAD(yld,2,12*nyear)  !copy yields for all years
    
  n=0
  DO k=1,nyear
    DO i=1,numcrop
	  DO t=1,12
	  	    crop_matrix(i,t+n)=crop_calendar(i,t)  !Repeat crop calendar for all years
	  END DO 
	END DO
  	n=n+12
  END DO
  DO i=1,numcrop
    DO t=1,12*nyear
	  cropYield(i,t)=yield_matrix(i,t)*crop_matrix(i,t)
	  IF (cropYield(i,t) <0) THEN
	    cropYield(i,t)=-9999
	  END IF
	END DO
  END DO  
  IF (optim==2) THEN    !***turn off for optimisation
    OPEN(unit=56, file='Outputs/MaxCropYield.txt', status='unknown')
    DO i=1,numcrop                         
      WRITE(56,*) (cropYield(i,t), t=1,12*nyear)
    END DO
    CLOSE(unit=56)
   END IF
  

END SUBROUTINE crop_matrix_module