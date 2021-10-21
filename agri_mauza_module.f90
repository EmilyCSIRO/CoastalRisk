SUBROUTINE agri_mauza_module(saltM, wtlg, nummz, nummzM, nyear,t, numcrop, crop_yield, area_hh,  crop_r, Yhhmax_m)
  IMPLICIT NONE
  
  REAL, INTENT(IN) :: saltM(nummz), wtlg(nummz), crop_yield(numcrop),  area_hh(4),  crop_r(numcrop)
  INTEGER, INTENT(IN) :: nyear, t, numcrop, nummz,nummzM
  REAL, INTENT(OUT) :: Yhhmax_m(nummzM,numcrop,4)  !Production in ha per household type per mauza  
  REAL :: soil,waterlog
  REAL, ALLOCATABLE :: Y_a(:,:)
  INTEGER ::  k,j,i,c
  
  ALLOCATE (Y_a(numcrop,4))

  DO c=1, nummz
    soil=MAX(saltM(c),0.0)
	waterlog=MAX(wtlg(c),0.0)
	!soil=0.0  !temp for testing impact of salinity
	!waterlog=0.0  !temp for testing impact of waterlogging
	
	DO i=1, SIZE(area_hh)	
	  DO j=1,numcrop !Crop 1 to 10
	    IF (crop_yield(j)>0) THEN  ! Adjust crop yield based on salinity and waterlogging
		  Y_a(j,i)= crop_yield(j) -(crop_yield(j)*MIN(0.5*soil+0.5*waterlog,1.0))  !t/ha per crop
	    ELSE
	      Y_a(j,i)=-9999
	    END IF
        !Calculate crop production based on the ratio and area for each crop per HH type (assuming no flooding)
		Yhhmax_m(c,j,i)=Y_a(j,i)*area_hh(i)*crop_r(j)  
		   
	    IF (Yhhmax_m(c,j,i)<0) THEN
	      Yhhmax_m(c,j,i) = -9999
	    END IF
	    IF (i==1) THEN
	      IF (j <=10) THEN
		    Yhhmax_m(c,j,i)=-9999    !Subsistence farmer, take only the last three crops	    
		  END IF
	    ELSE IF (i>1) THEN
	      IF (j >10) THEN
		    Yhhmax_m(c,j,i)=-9999    !Farm holdings, remove the last three crops
		  END IF
	    END IF	
	  END DO	  
	END DO 
  END DO
  
END SUBROUTINE agri_mauza_module