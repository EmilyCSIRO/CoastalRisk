 SUBROUTINE sediment (flooding,m_frfl, m_elev, p_frfl, elvband, nummz,nbands, elev_new, eband_new,t,p)
   !Calculates new elevation for each mauza based on subsidence and deposition (if included)
   IMPLICIT NONE
   REAL, INTENT(IN):: flooding, m_frfl(nummz), m_elev(nummz), p_frfl
   INTEGER, INTENT(IN) :: elvband(nbands), nummz,nbands,t,p
   REAL, INTENT(OUT) :: elev_new (nummz)
   INTEGER, INTENT(OUT) :: eband_new(nummz)
   
   REAL :: subs, deps
   INTEGER :: c, b
   
   subs=0.02/12.0   !Subsidence of 2cm/yr (Auerbach et al)
   !deps=0.37    !Deposition reported by Auerbach with an average of 37 cm and as high as 60-70cm in the two years after the storm. Turn on if using for deposition scenario
   deps=0.0
   
   DO c=1,nummz
     IF (m_frfl (c)>0.5 .AND. flooding > 0.001) THEN       !Mauza is flooded (frfl equal to 1. the 0.5 is arbitrary)
       elev_new(c)= m_elev(c)-subs+ deps*p_frfl
     ELSE
	   elev_new(c)= m_elev(c)-subs     !subsidence only
	 END IF
	 DO b=1,nbands
	   IF (elev_new(c) < elvband(b)) eband_new(c) = nbands+1-b !nbands+1-b to account for the order of the categories with category 10 read in first (ie b=1)
	   IF (elev_new(c) > elvband(b) .AND. b==1) WRITE(*,*) 'New elevation greater than existing band'
	 END DO
   END DO
    
 END SUBROUTINE sediment