!**************************************************************************************************************
!
!                          Coastal Risk Model
!
!  Risk Modelv20 developed by E.J. Barbour October 2021 with contributions in study design from
!  J.W. Hall, E. Borgomeo, M.S.G. Adnan, M. Salehin, M.S.A Khan, and K. Paprocki.
!  Results from this code are published in the manuscript: 
!  "The unequal distribution of water risks and adaptation benefits in coastal Bangladesh"
!
!*************************************************************************************************************

INCLUDE 'poisson_fixed_time.f90'   !Generates stochastic flood events
INCLUDE 'agri_mauza_module.f90'    !Calculates crop production based on salinity, waterlogging, and crop yields/ratios
INCLUDE 'crop_matrix_module.f90'   !Calculates matrix of crop yields for entire simulation based on cropping patterns
INCLUDE 'sediment.f90'             !Calculates new elevation for each mauza based on subsidence and deposition (if included)

MODULE RiskModel

! Reconcile totincomeHH (doesn't consider mauzas), TotIncomeBy_HH, and SumMauza_APhhAllCrop (not by polder)

IMPLICIT NONE

  !***********PARAMETERS**************************************************************
  !These are used to generate results during the Monte Carlo Simulation
  !*** Flooding ******
  !floodp=floodevent by month for each polder; floodm=floodevent per mauza (max 62) for each polder
  !pold_maxCropFlood = maximum crop area flooded for each polder and each year  
   REAL :: floodp(720,6), floodm(6,62),pold_maxCropFlood(6,60)
  
  !***Embankment and Salinity****
  !pold_eR=embankment reliability by polder and month
  !pold_sL=salinity for each polder and month; TInvCost=total invetment cost
  !confail is an index to identify if eR goes below zero. If it does, a value of >0 is assigned.
  !costfail is an index to identify if the total cost of investment goes above a user specified budget.
  REAL :: pold_eR(6,720), pold_sL(6,720)   
  INTEGER :: TInvCost
  INTEGER ::confail, costfail
  
  !*** Agriculture ***
  !SDamage=total damange for each HH type for all years and all mauzas; Damage=total damage by polder and HH
  !CropInc=crop income per polder, HH and year; CropIncBy_CropHH= crop income by polder, crop and HH
  !pold_avYieldBy_crop=average crop yield per polder and crop
  !AvIncomePerYrBy_HH= average crop income per polder and HH; AvIncomePerYrBy_HH_numhold=Divides by number of holdings per HH type by polder
  !poldermauzaAllYears_income=total crop income per polder and per mauza (max 62 mauzas)
  REAL :: SDamage(4), Damage(6,4), CropInc(6,4,60), CropIncBy_CropHH(6,13,4), pold_avYieldBy_crop(6,13)
  REAL :: AvIncomePerYrBy_HH(6,4),AvIncomePerYrBy_HH_numhold(6,4), poldermauzaAllYears_income(6,62)   

  !***FLAGS/SWITCHES*********************************************************
  
  !These can be set by the user to change the mode of operation
  ! inv=flag for which investment strategy to use (1 if time-based, otherwise based on embankment reliability), invYr=investment decision variable
  !Optimisation flag -set to 1 if on and 2 if off. If on, outputs to file will not be produced. 
  !Sed is flag for sediment subroutine. Set to 1 to activate, any other value to deactivate.
  !Historical simulation flag set to 1 if turned on  
  !ck= error check mode, set to 1 to check results against base values and ensure MCrun=1
  
  INTEGER ::  inv=0, opt=2, sed=1, his=0, ck=1  
  
  !***********************************************************************************************************
  
CONTAINS

SUBROUTINE MainProcess 
IMPLICIT NONE
  
  !Flooding module 
       !  ****PARAMETERS****
	   !   A_F=Initial fraction flooded at time of flood;
       !   Beta=Infrastructure damage parameter when flood occurs 
	   !   poldfrflood = function that calculates the fraction flooded for each polder for each timestep
       !   floodPold = time series of flood events external to each polder
	   !   mauzaflood = binary, indicates if a mauza is flooded (1) or unflooded (0)
       !   fr_fl (PolderType)=fraction flooded after flood event accounting for recession of floodwater 
	   !   timing_month=a time series of months of when a flood occurred external to the polder	   
	   !  ******************
   REAL :: A_F, beta=20.0, poldfrflood
   REAL, ALLOCATABLE :: floodPold(:,:)
   INTEGER, ALLOCATABLE :: mauzaflood(:),timing_month(:)
  
  !Embankment/drainage reliability
       !  ****PARAMETERS****
       !   i_e is the size of investment in embankment condition
	   !   alpha_0 and alpha_1=drainage infrastructure deterioration (0 is initial; 1 is during investment period)
       !   w_i_start and w_i = drainage condition (start = initial condition) 
       !   eR_start = initial embankment condition for all polders
	   !   InvestThresh = embankment condition used to trigger an investment in reliability-based investment (specified by user)
       !   MaxInvest = difference between maximum and current embankment reliability
	   !   TotInvest = total improvement in embankments across all polders
	   !   Budget & InvCost = maximum budget for investment in $US million; actual investment cost
	   !   Investment is the user defined investment for reliability-based, and is the fraction of needed investment that is used (1.0 = back to maximum reliability)
	   !   inv_start & inv_end = month for start and end of investment period (1-40 years excl warmup)
	   !   InvCost = !US million based on embankment length (estimated from Bangladesh Water Development Board, 2013, 
	   !      Bangladesh Water Development Board. Development Project Proforma/Proposal (DPP) for Blue Gold Program (BWDB Component). 80-83 .
	   !  ******************
   REAL :: i_e, alpha_0=-0.005, alpha_1=-0.005, w_i_start=0.8, eR_start=0.6  
   REAL :: InvestThresh, MaxInvest, TotInvest, Budget=100
   REAL :: Investment
   INTEGER :: inv_start=12, inv_end=492,InvCost(6)= [31,26,33,34,42,38] 
   REAL, ALLOCATABLE :: w_i(:)
   
  !Agriculture
       !  ****PARAMETERS****
       !   yield = maximum crop yield based on BBS values
	   !   crop_ratio = ratio of different types of crops for a polder
	   !   crop_yield = array of maximum crop yield based on BBS values assuming areas unflooded and repeated over entire simulation. 
	   !   price = taka per tonne per crop
       !   Prodmax_m_temp/Prodmax_m_tempABS = maximum crop production accounting for salinity and waterlogging by mauza, crop, farm type, month
       !   mauza_ProdFlood = crop production accounting for flooding, salinity, waterlogging, by mauza, crop, farm type, and month
       !   mauza_ProdYearMin = minimum production for an entire year by mauza, crop, farm type and year
	   !   mauza_ProdYearPoldTotHH = Production by farm type summed for an entire polder by crop, farm type, and year
       !   pold_avProdBy_crop = Average crop production over simulation by crop: summed over all years (shown for individual crops) and divided by number of years
	   !   pold_ProdBy_CropYr = Sum of crop production across all households by polder (for individual crops and each year) - same as polder(p)%ProdBy_CropYr but different format for calcs
       !   pold_avYieldBy_CropYr = Average crop yield for each crop and each polder per year
	   !   mauza_CropIncome (mauza, crop, HH, month) = production x crop price per mauza, crop, farm type and year
       !   SumMauza_CropIncome (crop, HH, year) = sum of crop income for all mauzas within a polder by crop, farm type, and year
	   !   mauza_income_HHYr (mauza, HH, year) = income for each mauza by farm type and year for all crops 
       !   mauza_totincomeHH (mauza, HH, year) = Total income per mauza for each HH and each year summed over number of each holding type per mauza
       !   mauza_totincome (mauza,year) = Total income per mauza for each year summed over number of each holding type per mauza
	   !   mauzaAllYears_income (mauza) = Total crop income per mauza
	   !   maxinc = the estimated maximum crop income that would be obtained for each HH type if there was no flooding or salinity
	   !  ******************  
   REAL, ALLOCATABLE :: yield(:),crop_ratio(:), crop_yield(:,:), price(:)
   REAL, ALLOCATABLE :: MaxCropProd(:,:),SumMaxCropProd(:),PropMaxProd(:,:)
   REAL, ALLOCATABLE :: Prodmax_m_temp(:,:,:,:), Prodmax_m_tempABS(:,:,:,:), mauza_ProdFlood(:,:,:,:),mauza_ProdYearMin(:,:,:,:)
   REAL, ALLOCATABLE :: mauza_ProdYearPoldTotHH(:,:,:,:), pold_avProdBy_crop(:,:), pold_ProdBy_CropYr(:,:,:)
   REAL, ALLOCATABLE :: pold_avYieldBy_CropYr(:,:,:)
   REAL, ALLOCATABLE :: mauza_CropIncome(:,:,:,:), SumMauza_CropIncome(:,:,:), mauza_income_HHYr(:,:,:)
   REAL, ALLOCATABLE :: mauza_totincomeHH(:,:,:), mauza_totincome(:,:),mauzaAllYears_income(:)
   INTEGER, ALLOCATABLE :: maxinc(:,:)
 
  !Geospatial
       !  ****PARAMETERS****
       !   PoldCropArea = total area (ha) being cropped for a polder across all mauzas and HH types
	   !   PoldCropAreaSubs = total area being cropped for a polder across all mauzas for subsistence farms
	   !   PoldCropAreaHold = total area being cropped for a polder across all mauzas for farm holdings
	   !   Elv_pold_sum = sum of elevations within a polder accounting for area (number of holding types in each mauza)
	   !   Elv_pold = average elevation for a polder based on that of mauzas within the polder
	   !   elevnew = new elevation for an individual mauza
   	   !   CropAreaBy_Crop = total area under each crop 
	   !   SumHoldingsHH = Total number of holdings by holding type
	   !   elevband = elevation that defines each elevation category
	   !   ecat_new = new elevation category
   REAL :: PoldCropArea, PoldCropAreaSubs, PoldCropAreaHold, Elv_pold_sum, Elv_pold  
   INTEGER :: numbands = 9 !Allocate number of elevation bands based on SRTM DEM used
   REAL, ALLOCATABLE ::  elevnew(:), CropAreaBy_Crop(:,:)     
   INTEGER, ALLOCATABLE :: SumHoldingsHH(:), elevband(:), ecat_new(:)
    
  !Other environmental factors
       !  ****PARAMETERS****
       !   salinity_profile = monthly pattern of salinity,    
       !   saltFlood = salinity resulting from a flood event
	   !   saltDecl = assumed rate of decline in salinity following a flood event
	   !   sumsalin = used to calculate total salinity in a polder to estimate average polder salinity based on mauza salinity
	   !   saltInt = gradual increase in baseline salinity over time due to saline intrusion
	   !   riverlevel_profile = monthly river levels
	   !   wlg_pold_sum = sum of waterlogging in a polder accounting for area to determine the average polder waterlogging
	   !   wlg_pold = average waterlogging for a polder
	   !   sp = time series of baseline salinity for entire simulation based on the salinity profile
	   !   eW = time series of baseline river levels for entire simulation based on the river level profile
	   !   saltCon & saltCon_start = initial salinity concentration increased over time based on saline intrusion
	   !   saltRiv = salinity concentration in the river for each polder
   REAL :: salinity_profile(12), saltFlood, saltDecl, sumsalin,saltInt
   REAL :: riverlevel_profile(12), wlg_pold_sum
   REAL, ALLOCATABLE :: sp(:), saltCon(:), eW(:),wlg_pold(:,:)
   INTEGER, ALLOCATABLE :: saltCon_start(:), saltRiv(:)
   
  !Household variables
   REAL :: area_hh(4)   ! Area of each HH type
  
  !General temporary variables
   !i= polder 1 to 6; t = time (months); tt= time excluding warm up; numpolders = total polders in available data
   !nummauzaP = number of mauzas in each polder; spold = number of case study polders
   REAL :: invDV !Investment Decision Variable - timing of investment
   INTEGER :: i, t, tt, numpolders, nummauzaM, nummauzaP, spold 
   INTEGER :: nyear, k, n, y, m,j, ncrop, c, b
   INTEGER :: t_it, t_it2, t_it3, t1,  p, numupazilas,d
   INTEGER :: poldnum(6), filenum, nummauzaAll(6), fileElevation(6), fileHoldings(6)
   INTEGER :: MCnum, numevents
   CHARACTER(20) :: dummy2
   CHARACTER(40) :: dummy1,filename,filename2, filename3, filenameMauzInc
   
  !Error Checks
  !poldfrfl is the base polder(p)%fr_fl, wlg_pold_Check is the check for average polder waterlogging, 
  !saltCheck is for polder salinity polder(p)%sL(t), Check is a generic check between expected and actual
  !CropIncCheck is for yearly total crop income per polder and HH type, recorded in CropIncomeAllMauzaAllCrops.txt
   REAL:: poldfrfl_Check(6,732),wlg_pold_Check(6,732), saltCheck(6,732), Check, embCheck(6,732)

  !Variables grouped as types  
  TYPE PolderType
    CHARACTER(10) :: polderID
	CHARACTER(40) :: polderName
	CHARACTER(30) :: district, upazila
	REAL, ALLOCATABLE :: fr_fl(:), sL(:), eR(:)
	REAL, ALLOCATABLE :: invest(:)
    REAL, ALLOCATABLE :: mauza_frfl(:)
    REAL, ALLOCATABLE :: mauzaArea(:)
	REAL, ALLOCATABLE :: TotProdBy_Yr(:)
	REAL, ALLOCATABLE :: ProdBy_CropHHYr(:,:,:), ProdBy_HHYr(:,:), ProdBy_CropYr(:,:)
	REAL, ALLOCATABLE :: Elev(:), Elev_start(:)
	INTEGER, ALLOCATABLE ::mauzaHoldingID(:), mauzaHoldings(:,:), mauzaID(:)
	INTEGER, ALLOCATABLE :: ElevCat_frfl(:), ElevCat(:)
	REAL :: d_r
	REAL :: floodband(10), CropFlood
	INTEGER :: numholdings(4), upazID, invYr(3),Cost
  END TYPE PolderType
  TYPE (PolderType), ALLOCATABLE :: polder(:)
  
  TYPE MauzaType
	REAL :: saltM, saltM_scale, waterLgM
	REAL, ALLOCATABLE :: Yhhmax(:,:,:), Yhhmax_temp2(:,:,:)
  END TYPE MauzaType
  TYPE (MauzaType), ALLOCATABLE :: mauza(:)
  
  TYPE CropRatioType
    CHARACTER(30) :: district, upazila
	REAL :: ratio(13)    
	INTEGER :: upazilaID, districtID
  END TYPE CropRatioType
  TYPE (CropRatioType), ALLOCATABLE :: cropR(:)
  
  !Intialise parameter values
  nyear=61   
  numpolders=98  !total number of polders for reading in data
  spold=6  !number of test case polders
  numupazilas=60
  nummauzaM=62  !set to largest number of nummauza for allocation purposes
  ncrop=13
  nummauzaAll = [62, 37, 6, 9, 37, 21]
  poldnum = [31,32,35,36,65,86] !polders 29(31), 30 (32), 32(35), 33(36), 43/1(65), 54(86)  
  numevents=nyear-1 !Excludes warm up year
  
  !*********Confirm parameters used with user*************************************  
  IF (opt==2) THEN  !Write to screen parameters to check using correct values
    WRITE(*,*) "Number years is    ", numevents
	WRITE(*,*) "DEM used: SRTM"
	IF (inv==1) THEN
	  WRITE(*,*) "Investment strategy: Time-based; investment at   "
	ELSE
	   WRITE(*,*) "Investment strategy: Reliability-based; investment at   ", InvestThresh, "   of  ", Investment
	END IF
	IF (his==1) THEN
	  WRITE(*,*) "WARNING - Running in historical mode"
	END IF
  END IF
  !*******************************************************************************


 !***********ALLOCATE ARRAYS*************************   
  !Flooding module
  ALLOCATE (floodPold(spold,12*nyear)) 
  ALLOCATE (timing_month(12*nyear))
  ALLOCATE (mauzaflood(nummauzaM))  
  
  !Embankment/drainage reliability 
   ALLOCATE (w_i(12*nyear)) 
  
  !for allocation purposes set arrays to largest number of nummuaza for all sample polders 
  
  !Agriculture
  ALLOCATE (maxinc(spold,ncrop),price(ncrop), crop_ratio(ncrop), cropR(62))
  ALLOCATE (yield(ncrop),crop_yield(ncrop,12*nyear),MaxCropProd(spold,ncrop),PropMaxProd(spold,numevents))
  ALLOCATE (SumMaxCropProd(spold))  
  ALLOCATE (mauza_CropIncome(nummauzaM,ncrop,4,12*nyear), SumMauza_CropIncome(ncrop,4,numevents))
  ALLOCATE (mauza_income_HHYr(nummauzaM,4,nyear))
  ALLOCATE (mauza_totincome(nummauzaM,nyear),mauzaAllYears_income(nummauzaM))
  ALLOCATE (mauza_totincomeHH(nummauzaM,4,nyear))
  ALLOCATE (mauza_ProdYearPoldTotHH(nummauzaM,ncrop,4,nyear))
  ALLOCATE (SumHoldingsHH(4),elevband(numbands),ecat_new(nummauzaM),elevnew(nummauzaM))
  ALLOCATE (Prodmax_m_temp(nummauzaM,ncrop,4,12*nyear),Prodmax_m_tempABS(nummauzaM,ncrop,4,12*nyear))
  ALLOCATE (mauza_ProdFlood(nummauzaM,ncrop,4,12*nyear), mauza_ProdYearMin(nummauzaM,ncrop,4,nyear))
  ALLOCATE (pold_avProdBy_crop(numpolders,ncrop),pold_ProdBy_CropYr(numpolders,ncrop,numevents))  
  ALLOCATE (pold_avYieldBy_CropYr(numpolders,ncrop,numevents)) 
  
  !Geospatial
  ALLOCATE (polder(numpolders),mauza(nummauzaM), CropAreaBy_Crop(numpolders,ncrop))  
 
  !Other environmental factors
  ALLOCATE (sp(12*nyear))
  ALLOCATE (saltCon_start(spold), saltCon(12*nyear),saltRiv(spold))
  ALLOCATE (eW(12*nyear), wlg_pold(spold,12*nyear))  
  
 DO i=1, spold
   p = poldnum(i)
   nummauzaP = nummauzaAll(i)
   ALLOCATE (polder(p)%fr_fl(12*nyear), polder(p)%eR(12*nyear), polder(p)%sL(12*nyear))
   ALLOCATE (polder(p)%invest(12*nyear))
   ALLOCATE (polder(p)%ElevCat_frfl(10))
   ALLOCATE (polder(p)%ElevCat(nummauzaP),polder(p)%mauzaID(nummauzaP),polder(p)%mauzaHoldingID(nummauzaP))
   ALLOCATE (polder(p)%Elev(nummauzaP),polder(p)%Elev_start(nummauzaP))
   ALLOCATE (polder(p)%mauzaHoldings(nummauzaP,4),  polder(p)%mauza_frfl(nummauzaP),  polder(p)%mauzaArea(nummauzaP))
   ALLOCATE (polder(p)%TotProdBy_Yr(nyear),polder(p)%ProdBy_CropHHYr(ncrop,4,nyear), polder(p)%ProdBy_HHYr(4,numevents))   
   ALLOCATE (polder(p)%ProdBy_CropYr(ncrop,nyear))
   END DO
 
 DO i=1, nummauzaM
   ALLOCATE (mauza(i)%Yhhmax(ncrop,4,12*nyear), mauza(i)%Yhhmax_temp2(ncrop,4,12*nyear))
 END DO 


 !******************READ PARAMETER FILES*************************
         !********Investment Parameters******************** 
   DO i=1, spold  !Assign investment cost to each polder
     p = poldnum(i)
	 polder(p)%Cost=InvCost(i)
   END DO
   IF (inv==1) THEN   ! Time-based investment
     OPEN(unit=60, file='InvParams001.txt', status='old') !Specified when investment should occur for all polders
	 DO j=1,3  !Three investments can be selected for each polder
	   DO i=1,spold
	     p = poldnum(i)
	     READ(60,*) invDV  !read as index from 0 to 1 from optimisation, where 0-0.11=5yrs, 0.12-0.22=10 years etc >0.8 is no investment
	     IF (invDV<=0.1111) THEN
  	       polder(p)%invYr(j) = 5*12+12   !convert 5 years to months plus warm up period of 12 months
	     ELSE IF (invDV<=0.2222) THEN
	       polder(p)%invYr(j) = 10*12+12 
	     ELSE IF (invDV<=0.3333) THEN
	       polder(p)%invYr(j) = 15*12+12 
	     ELSE IF (invDV<=0.4444) THEN
	       polder(p)%invYr(j) = 20*12+12 
	     ELSE IF (invDV<=0.5556) THEN
	       polder(p)%invYr(j) = 25*12+12 
	     ELSE IF (invDV<=0.6667) THEN
	       polder(p)%invYr(j) = 30*12+12 
	     ELSE IF (invDV<=0.7778) THEN
	       polder(p)%invYr(j) = 35*12+12
	     ELSE IF (invDV<=0.8889) THEN
	       polder(p)%invYr(j) = 40*12+12
	     ELSE 
	       polder(p)%invYr(j) = 999999   ! No polder investment, set to time which the code will not reach 
	     END IF
	   END DO
	 END DO

   ELSE   !Reliability-based investment (investment based on embankment reliability)
     OPEN(unit=60, file='RiskParams001.txt', status='old')
	 READ(60,'(E30.16)') InvestThresh
     READ (60,*) Investment
   END IF
   CLOSE(unit=60)
         !******************************************* 	 
		 
         !********Geospatial and Erosion********************   
   !Read in crop ratios to polders based on Upazila-scale data, and assign high or low erosion rates to each polder 
   OPEN(unit=71, file='Inputs/MajorCropRatioUpazila.txt', status='old')
   READ(71,*)
   DO i=1, numupazilas
     READ(71,*) cropR(i)%upazilaID, cropR(i)%ratio(:), cropR(i)%district, cropR(i)%districtID, cropR(i)%upazila 
   END DO   
   CLOSE(unit=71)
 
   OPEN (unit=63, file='Inputs/PolderErosion.txt', status='old')
   OPEN (unit=64, file='Inputs/PolderUpazila.txt', status='old')
   READ(63,*)
   READ(64,*)  
   DO i=1, numpolders  !read in data for all 98 polders
	 READ(63,*) dummy1, dummy2
	 IF (dummy2=='N') THEN
	   polder(i)%d_r= 0.0025                 !half the rate of deterioration
	 ELSE IF (dummy2=='Y') THEN
	   polder(i)%d_r= 0.005                  !standard polder deterioriation rate set by Edoardo - anything higher sets eR to zero
	 ELSE 
	   WRITE(*,*) 'Error: no erosion label entered for polder ', i 
	 END IF
	 READ(64,*) dummy1, polder(i)%upazID, polder(i)%district, polder(i)%upazila
   END DO  
   CLOSE(unit=63)
   CLOSE(unit=64)
   !*******************************************
   
   !********Elevation*********************************************
   !This is the starting elevation for each mauza, read in as a category (elevations divided into 9 categories)
   !Cat1=-8 to 2m; Cat2=2 to 3m; Cat3=3 to 3m; Cat4=4-5m etc; Cat 9=9 to 23m
   IF (his==1) THEN   !Running in historical mode, add 1m to elevation (Auerbach)
     OPEN (unit=72, file='Inputs/ElevationCat_Mauza29_SRTM_hist.txt', status='old')
     OPEN (unit=76, file='Inputs/ElevationCat_Mauza30_SRTM_hist.txt', status='old')
     OPEN (unit=77, file='Inputs/ElevationCat_Mauza32_SRTM_hist.txt', status='old')
     OPEN (unit=78, file='Inputs/ElevationCat_Mauza33_SRTM_hist.txt', status='old')
     OPEN (unit=79, file='Inputs/ElevationCat_Mauza43_SRTM_hist.txt', status='old')
     OPEN (unit=80, file='Inputs/ElevationCat_Mauza54_SRTM_hist.txt', status='old')
   ELSE
     OPEN (unit=72, file='Inputs/ElevationCat_Mauza29_SRTM.txt', status='old')
     OPEN (unit=76, file='Inputs/ElevationCat_Mauza30_SRTM.txt', status='old')
     OPEN (unit=77, file='Inputs/ElevationCat_Mauza32_SRTM.txt', status='old')
     OPEN (unit=78, file='Inputs/ElevationCat_Mauza33_SRTM.txt', status='old')
     OPEN (unit=79, file='Inputs/ElevationCat_Mauza43_SRTM.txt', status='old')
     OPEN (unit=80, file='Inputs/ElevationCat_Mauza54_SRTM.txt', status='old') 
   END IF	 

   READ(72,*)
   READ(76,*)
   READ(77,*)
   READ(78,*)
   READ(79,*)
   READ(80,*)

   fileElevation = [72,76,77,78,79,80]
   
   DO i=1, spold
     nummauzaP = nummauzaAll(i)
	 p = poldnum(i)
	 filenum = fileElevation(i)
	 DO c = 1, nummauzaP   
       READ(filenum,*) polder(p)%mauzaID(c), dummy1, polder(p)%Elev(c), polder(p)%ElevCat(c)	  
	   polder(p)%Elev_start(c)=polder(p)%Elev(c)  !Assign starting elevation category to each mauza in each polder
	 END DO
	 CLOSE (unit=filenum)
   END DO

   OPEN (unit=67, file='Inputs/ElevationArea_SRTM.txt', status='old') !relates fraction flooded to elevation bands
   OPEN (unit=68, file='Inputs/ElevationBands_SRTM.csv', status='old') !the 9 elevation bands
   READ(68,*)
   DO i=1,2
     READ(67,*)
   END DO
   DO i=1,numbands
     READ(67,*) polder(31)%floodband(i), polder(32)%floodband(i), polder(35)%floodband(i), polder(36)%floodband(i), &
	       polder(65)%floodband(i), polder(86)%floodband(i)
     READ(68,*) elevband(i)
   END DO  
   CLOSE(unit=67)
   CLOSE(unit=68)  
   !***********************************************************************************************************
   
   !*********Farm Holdings ************************************************************************************
   !This reads in the number of small, medium, large and non-farm holdings for each mauza for each polder
    OPEN(unit=73, file='Inputs/MauzaFarmHoldings29.txt', status='old')
	OPEN(unit=81, file='Inputs/MauzaFarmHoldings30.txt', status='old')
	OPEN(unit=82, file='Inputs/MauzaFarmHoldings32.txt', status='old')
	OPEN(unit=83, file='Inputs/MauzaFarmHoldings33.txt', status='old')
	OPEN(unit=84, file='Inputs/MauzaFarmHoldings43.txt', status='old')
	OPEN(unit=85, file='Inputs/MauzaFarmHoldings54.txt', status='old')
	READ(73,*)
	READ(81,*)
	READ(82,*)
	READ(83,*)
	READ(84,*)
	READ(85,*)
		
	fileHoldings=[73, 81,82,83,84,85]
		
	DO i=1,spold
	  nummauzaP = nummauzaAll(i)
	  p = poldnum(i)
	  filenum = fileHoldings(i)
	  DO c=1, nummauzaP
	    READ(filenum,*)polder(p)%mauzaHoldingID(c), dummy1, polder(p)%mauzaArea(c), polder(p)%mauzaHoldings(c,:)   			
   	    !Error check for inputs
	   	IF(polder(p)%mauzaID(c)-polder(p)%mauzaHoldingID(c)>0.001) THEN
	      WRITE(*,*) "Mauza IDs do not match - check ElevationCat and MauzaFarmHolding text files"
	        PAUSE
	    END IF	  	  
	  END DO
	  CLOSE (unit=filenum)
	END DO
	!*********************************************************************************************************


  !***************ASSIGN PARAMETER VALUES**********************************************************************
  !
  !**************Matrix of crop yields*********************************************************
    !crops: local aus, hyv aus, local aman; hyv aman, local boro, hybrid boro, hyv boro, wheat, maize, jute, water gourd, winter pumpkin, summer pumpkin 
    !price is a rough estimate, insufficient data to justify a time series at present.
    !price is in taka per tonne, no price found for hybrid boro, assumed to be the same as hyv boro> Jute is in taka per bale (180kg)
    !area_hh is the size of land owned in ha for different types of households (polder 29), subsistence farmer; small farmer; medium; large	
	!sizes are calculated based on 2008 BBS agricultural statistics using total operated area and total number of holdings for Dumuria and Batiaghata upazillas
	!maxinc is the estimated maximum crop income that would be obtained for each HH type if there was no flooding or salinity
    !saved under C:\Emily\University_Oxford\REACH\Khulna Observatory\WP4_CoastalRiskModelling\RiskModelCode\V4\Outputs\Run5_Calibrated\CropYieldCalcs_P29
	yield =(/2.5,5.5,2.5,5.5,2.5,10.0,7.0,5.0,8.0,15.0,28.0, 35.0, 35.0/) !t/ha, except jute which is in bales
    price=[13360,13150,13680,15660,14330,16220,16220,20450, 16560,8557,13300, 13240, 13240]       
    area_hh = (/0.01, 0.51, 2.0, 6.6/)    
    maxinc(1,:)=[3901000,333957000,339716000,123859000] 
    maxinc(2,:)=[1522000,109667000,143508000,53761000]
    maxinc(3,:)=[4116000,83658000,86132000,68443000]
    maxinc(4,:)=[4759000,148566000,161293000,84563000]
    maxinc(5,:)=[7194000,161674000,204000000,113461000]
    maxinc(6,:)=[4437000,95685000,172664000,109330000]  
	
    CALL crop_matrix_module (yield(:),nyear, ncrop, opt, crop_yield(:,:))  !Generate crop calendar and maxtrix of maximum yields for entire simulation
  !*************************************************************************************************** 
  
  !***************Salinity time series**************************************************************** 
    salinity_profile = [ 0.15,0.17,0.8,0.85,0.9,0.3,0.15,0.05,0.05,0.05,0.05,0.05 ] ! seasonal salinity profile from Lazar et al [2015]
    riverlevel_profile=[1.9,2.0,2.1,2.3,2.4,2.4,2.6,2.6,2.6,2.5,2.3,2.1]  !m from three observed gauges    
	!Generate a time series of salinity levels and river levels DO t_it=0,nyear-1
    DO t_it=0,nyear-1
      DO t_it2=1,12
        sp(t_it*12 + t_it2)=salinity_profile(t_it2)
	    eW(t_it*12 + t_it2)=riverlevel_profile(t_it2)
      END DO
    END DO  
  !***************************************************************************************************************
  
  !*************Embankments and Drainage**************************************************************************
   !i_e is the investment made to improve embankment reliability
   !TotInvest is the total embankment investment 
   !TInvCost is the total embankment investment cost
   i_e = 0  !initialise to zero so investments aren't made automatically
   TotInvest=0.0
   TInvCost=0
   confail = 0   !confail is an index to identify if eR goes below zero. If it does, a value of >0 is assigned.
   costfail = 0 !costfail is an index to identify if the total cost of investment goes above a user specified budget.
   w_i = 1d0   !drainage infrastructure
   w_i(:) = w_i(:)*w_i_start
   !**************************************************************************************************************

   !*************Flooding and agriculture ************************************************************************
   floodp=0d0
   floodm=0d0 !EJB 28/10/2020
   MaxCropProd = 0d0
   pold_maxCropFlood=0d0  !EJB 1/10/2021
   !**************************************************************************************************************

   !*************Miscellaneous - investment, flood, timestep******************************************************
   DO i = 1, spold
     p = poldnum(i)
	 polder(p)%invest=0d0
!	 polder(p)%poldermauza_totincome=0d0
	 polder(p)%CropFlood=0d0
	 polder(p)%sL=0d0
	 polder(p)%eR(:)=1d0
	 polder(p)%eR(:)=polder(p)%eR(:)*eR_start
   END DO

   t1=13  ! First year (first 12 months) used as warm up
   !*************************************************************************************************************


   !-------------------------------------------------------------------------------------------------------------
   !                         ***START MAIN PROCESS***
   !-------------------------------------------------------------------------------------------------------------

  
   !******************FLOODING MODULE*************************************************************************** 
    !****Calculate external flood events*******
    !Generate different flood event sequence for each polder based on varying magnitude
	CALL poisson_fixed_time (12*nyear, MCnum, floodPold(:,:), timing_month(:))  
 
    DO t=1, 12*numevents      !save flood series excluding warm-up period	
      floodp(t,:)=floodPold(:,t+t1-1)
    END DO 
   !************************************************************************************************************
 


   !*****************FLOOD IMPACT ON CROP PRODUCTION************************************************************
     !loop through each polder and each time step to calculate the crop area flooded for each mauza, the impact of salinity,
	 !    and the subsequent impact on crop production and crop income.

  DO i = 1, spold !numpolders 
    !Initialise parameter values
	p = poldnum(i)
	nummauzaP = nummauzaAll(i)
	mauzaflood=0d0
	y=1
    m=1
	mauza_ProdYearMin=1d0
	mauza_ProdYearMin=9999*mauza_ProdYearMin
	mauza_ProdYearPoldTotHH=0d0  
	SumMauza_CropIncome=0d0
	mauza_income_HHYr=0d0
	mauza_CropIncome=0d0
	mauza_totincome=0d0
	mauza_totincomeHH=0d0
	mauzaAllYears_income=0d0   
	ecat_new=0d0
	elevnew=0d0
	Prodmax_m_temp=0d0
	Prodmax_m_tempABS=0d0


	!Assign crop area ratios to each polder by upazila	
    crop_ratio=0d0
	DO d=1, numupazilas
	 IF (polder(p)%upazID==cropR(d)%upazilaID) THEN
		crop_ratio(:)= cropR(d)%ratio(:)
	  END IF
	END DO
	
	!Calculate polder crop areas including by HH type 
	PoldCropArea=0.0
	PoldCropAreaSubs=0.0  
	PoldCropAreaHold=0.0  
	DO c=1,nummauzaP
	  DO k=1,4  !HH1 to 4
        PoldCropArea = PoldCropArea + polder(p)%mauzaHoldings(c,k)*area_hh(k)			  
        IF (k == 1) PoldCropAreaSubs = PoldCropAreaSubs +  polder(p)%mauzaHoldings(c,k)*area_hh(k)    
        IF (k /= 1) PoldCropAreaHold = PoldCropAreaHold + polder(p)%mauzaHoldings(c,k)*area_hh(k)		  
	  END DO
	END DO
    
	!Calculate the total maximum potential crop production based on area and yield to compare with actual production
	DO j=1, ncrop  !EJB whole loop
      IF (j <11) THEN
	    CropAreaBy_Crop(i,j)=PoldCropAreaHold*crop_ratio(j)    
      ELSE
	    CropAreaBy_Crop(i,j)=PoldCropAreaSubs*crop_ratio(j)    
      END IF
        MaxCropProd(i,j) = yield(j)*CropAreaBy_Crop(i,j)   !tonnes
    END DO
    SumMaxCropProd(:)=SUM(MaxCropProd(:,:),DIM=2) !Total production for all crops
	 
    !Open files as loop through polders to write out mauza results
    IF (opt==2) THEN 
	  WRITE(filenameMauzInc,'(a,i4.4,a)') "Outputs/TotMauzaIncome",p,".txt"
      OPEN(unit=840+p, file=filenameMauzInc, status='unknown')
	  WRITE(filename2,'(a,i4.4,a)') "Outputs/elevation_mauza",p,".csv"
	  OPEN(unit=590+p, file=filename2, status='unknown')
	  WRITE(filename3,'(a,i4.4,a)') "Outputs/elevation_mauzaCat",p,".csv"  !Elevation category for each mauza
	  OPEN(unit=401+p, file=filename3, status='unknown')
    END IF	
	
 
    !************Loop through each month from time t=13 to total number of years
    DO t=t1, 12*nyear   
      tt=t-12  ! to start at tt=1, t=13

	  !*******EMBANKMENT AND DRAINAGE RELIABILITY****************************************************************
	  !*******Calculate embankment condition (reliability) based on flood events and any investments************
	  
	  !Reliability-based investment (inv/=1) is when investment is triggered by embankment condition falling below a user-specified value
	  !The user also specifies whether the embankment is fully restored (Investment=1) or only maintained (Investment<1)
	  IF (t>= inv_start .AND. t<= inv_end .AND. inv/=1) THEN  !only calculate investment when within investment period
	    IF (polder(p)%eR(t-1)<investThresh) THEN
	      MaxInvest=1.0-polder(p)%eR(t-1)
	      i_e=MaxInvest*Investment  !where investment is a fraction between 0 and 1 
	      polder(p)%invest(t)=i_e
	      TotInvest=TotInvest + polder(p)%invest(t)
	    END IF
	  END IF
	  
	  !Time-based investment where the user specifies when an investment occurs (or can be generated by optimisation)
	  !It is assumed that the polder is completely rehabilitated, hence i_e is 1 minus the current condition
	  !Three investments are allowed for year polder, hence a max of invYr(3)
	  IF (inv==1) THEN
	    IF (t== polder(p)%invYr(1) .OR. t== polder(p)%invYr(2) .OR. t== polder(p)%invYr(3)) THEN
	      i_e = 1.0 - polder(p)%eR(t-1)
	      TInvCost = TInvCost + polder(p)%Cost   !US million
		  IF (TInvCost>Budget) costfail=costfail+1
	    ELSE
	      i_e = 0.0  !allow invYr > run time time to represent no investment
	    END IF	
	  END IF


      !Calculate embankment reliability and drainage condition 
	  !Considers flood events (floodPold), ongoing deterioration (d_r & alpha), and any investment (for embankments) (i_e)
	  IF (t==1) THEN
	    polder(p)%eR(t) = eR_start - polder(p)%d_r*eR_start - floodPold(i,t)/beta
	    w_i(t) = w_i_start +alpha_0*w_i_start
	  ELSE IF (t>= inv_start .AND. t<= inv_end) THEN
	    polder(p)%eR(t)=min(polder(p)%eR(t-1) - polder(p)%d_r*polder(p)%eR(t-1) - floodPold(i,t)/beta + i_e,1.0)   !calculation of eR seems to be working correctly, difficult to directly compare due to random number generation in flood
	    w_i(t) = min(w_i(t-1) + alpha_1*w_i(t-1),1.0)  ! correct, exactly the same as Edoardo's for all 40 years
      ELSE 
	    polder(p)%eR(t) = polder(p)%eR(t-1) - polder(p)%d_r*polder(p)%eR(t-1) - floodPold(i,t)/beta
        w_i(t) = w_i(t-1) +alpha_0*w_i(t-1)
	  END IF
	  IF (polder(p)%eR(t)<0) THEN ! interrupt simulation if embankment reliability goes below zero
        WRITE(*,*)'flood defence reliability is less than zero for polder ',p, 'run the script again'
	    polder(p)%eR(t)=0.0
	    confail = confail + 1
      END IF 
	 	
	  pold_eR(i,tt)=polder(p)%eR(t)
	
	  i_e=0.0 !reset to zero for next investment to be made, otherwise will keep investing every year.
      !********************************************************************************************

      !*********INTERNAL POLDER FLOODING**********************************************************
      !Calculate the fraction flooded for each polder using the function poldfrflood
	   polder(p)%fr_fl(t)=poldfrflood(timing_month(:),polder(p)%eR(:), floodPold(i,:), t, nyear)
	 
	  !Estimate which mauzas are flooded based on total polder fraction flooded
       DO c=1, SIZE(polder(p)%ElevCat_frfl)  !Set all mauzas to 1 (flooded) to initialise
         polder(p)%ElevCat_frfl(c) = 1
       END DO
       DO b=1,numbands  !Identify which polder elevation band is flooded (1 if flooded, 0 if not flooded)
         IF (polder(p)%fr_fl(t)<polder(p)%floodband(b)) polder(p)%ElevCat_frfl(numbands+1-b)=0 !11-b to account for the order of the categories with category 10 read in first (ie b=1)
       END DO 
	 
	   polder(p)%CropFlood=0d0
	   DO c=1, nummauzaP  	 !Assign flooded/not flooded category to each mauza based on elevation band
	     DO d=1, numbands
	       IF (ABS(polder(p)%ElevCat(c) - d)<0.001) THEN !e.g. if mauza cat3 elevation it will be assigned a 1 if that band is flooded or 0 otherwise
		     polder(p)%mauza_frfl(c) = polder(p)%ElevCat_frfl(d)
		     floodm(i,c)=floodm(i,c) + polder(p)%mauza_frfl(c) !EJB 28/10/2020
		   END IF
	     END DO
     		       		 
	     DO k = 1, size(area_hh) !calculate the total crop area flooded based on number mauzas flooded, number of each holding, and the size of holding type
	       polder(p)%CropFlood = polder(p)%CropFlood + polder(p)%mauza_frfl(c)*REAL(polder(p)%mauzaHoldings(c,k)) &
		       *area_hh(k)
	     END DO		 	   
	   END DO	   
  	   pold_maxCropFlood(i,y)= max(pold_maxCropFlood(i,y), polder(p)%CropFlood) !maximum crop area flooded for each polder for each year
      !********************************************************************************************************
	 
	  !************ELEVATION AND WATERLOGGING******************************************************************
       !Update elevation based on sediment deposition from flood	 
	  IF (sed==1) THEN  !Only activate when sed set to 1 by user, otherwise function with no subsidence
	    mauza(:)%waterLgM=0d0  !reset to zero for each time step 
	  
	    !Calculate the new elevation for each mauza considering subsidence
	    CALL sediment(floodPold(i,t), polder(p)%mauza_frfl(:),polder(p)%Elev(:), polder(p)%fr_fl(t), elevband(:), & 
	        nummauzaP, numbands, elevnew(:), ecat_new(:),t,p)
	
	    !Update elevation data
	    Elv_pold_sum=0.0
	    wlg_pold_sum=0.0
	   
	    DO c=1,nummauzaP
	      polder(p)%Elev(c)=elevnew(c)	  
	    ! Update elevation category for mauza based on sediment deposition/subsidence
	      IF (ABS(polder(p)%ElevCat(c)-REAL(ecat_new(c))) > 0.5) THEN
	        polder(p)%ElevCat(c) = REAL(ecat_new(c))  !elevation band has changed	
          END IF	  
	  	  !Calculate waterlogging based on river level (eW), polder elevation, and drainage condition (w_i) 
	      IF (eW(t)<polder(p)%Elev(c)) THEN
		    mauza(c)%waterLgM=0.5*(eW(t)/polder(p)%Elev(c)) + 0.5*(1.0-w_i(t))
		  ELSE
		    mauza(c)%waterLgM = 1.0  !River is higher than the polder 
		  END IF				
	      DO k=1,4
 	        Elv_pold_sum=Elv_pold_sum + polder(p)%Elev(c)*polder(p)%mauzaHoldings(c,k)*area_hh(k)
		    wlg_pold_sum=wlg_pold_sum + mauza(c)%waterLgM * polder(p)%mauzaHoldings(c,k)*area_hh(k)
	      END DO
	    END DO
	    Elv_pold=Elv_pold_sum/(PoldCropArea)  !Average elevation of the polder based on the revised elevation of each mauza (crop area only)
	    wlg_pold(i,t)=wlg_pold_sum/PoldCropArea  !Average waterlogging for the polder
	  	  
	    !Write new mauza elevation for each polder plus the elevation category
	    IF (opt==2) THEN 		 
	      WRITE(590+p,'(62(3x,f14.4))') polder(p)%Elev
		  WRITE(401+p,'(62(3x,i2.2))') polder(p)%ElevCat
	    END IF   
	  END IF !End if sed=1
	  !********************************************************************************************************
	
	  !************ SALINITY **********************************************************************************
	  !w_i is the state of the drainage infrastructure 
	  !sp is a monthly salinty profile which is the same for every year	
	
	  IF(his==1) THEN
	    saltCon_start=[4,3,8,8,3,4]   !estimated historical concentration 60 years before observed data
	  ELSE
	    saltCon_start=[9,9,14,14,8,9]   !average dS/m at start of simulation
	  END IF
	  saltRiv=[17,28,28,28,2,28]
	  saltInt=0.09/12.0   !Historical increase in soil salinity roughly 0.09 dS/m per year due to saline intrusion. 
	  saltDecl= 0.03  !0.03   !Decline in salt concentration following a flood. Assumes takes roughly one year to halve
	  mauza(:)%saltM = 0d0
	  mauza(:)%saltM_scale = 0d0
	  
	  ! Calculate salinity based on seasonal salinity profile and gradual increase based on saline intrusion 
	  IF (t==13) THEN 
	    saltCon(t)=REAL(saltCon_start(i))
	  ELSE
        saltCon(t)= saltCon(t-1) + saltInt	   
      END IF
      polder(p)%sL(t) = saltCon(t)*sp(t)
      mauza(:)%saltM = polder(p)%sL(t)
 
	  ! Calculate salinity due to flood inundation for the entire polder and take the maximum of this and the salinity profile
	  IF (timing_month(t)>0 .AND. polder(p)%fr_fl(t)>0) THEN 
	    DO c=1,nummauzaP
	      IF (polder(p)%mauza_frfl(c)>0.5) THEN!Mauza flooded	
	        mauzaflood(c) = 1	        			
	        saltFlood=REAL(saltRiv(i))*exp(-saltDecl*REAL(mauzaflood(c)))
	        mauza(c)%saltM = MAX(polder(p)%sL(t), saltFlood)
	      ELSE IF (mauzaflood(c)>0) THEN   !flood event but mauza not flooded, previous flood at mauza
            mauzaflood(c) = mauzaflood(c) + 1
		    saltFlood=REAL(saltRiv(i))*exp(-saltDecl*REAL(mauzaflood(c)))
		    mauza(c)%saltM = MAX(polder(p)%sL(t), saltFlood)
	      END IF			
	    END DO
	  ELSE
	    DO c=1,nummauzaP
	      IF(mauzaflood(c)>0) THEN   !Previous flood event at mauza
		    mauzaflood(c) = mauzaflood(c) + 1
	        saltFlood=REAL(saltRiv(i))*exp(-saltDecl*REAL(mauzaflood(c)))
		    mauza(c)%saltM = MAX(polder(p)%sL(t), saltFlood)
		  END IF
	    END DO
	  END IF
	  
!	  CALL salt(t, i, spold, timing_month(:), nyear, polder(p)%fr_fl(t), nummauzaP, polder(p)%mauza_frfl(:), sp(t), & 
!	       saltCon(t), mauza(:)%saltM, polder(p)%sL(t))
	  	 
	  !Calculate polder average salinity
	  sumsalin=0.0
	  DO c=1, nummauzaP
	    DO k=1, 4
         sumsalin=sumsalin + mauza(c)%saltM*polder(p)%mauzaHoldings(c,k)*area_hh(k)
	    END DO
	  END DO
	  polder(p)%sL(t) = sumsalin/(PoldCropArea)  
	  
	  !Scale salinity to a value between 0 and 1 where >= 16 dS/m =1
	  DO c=1,nummauzaP
	    mauza(c)%saltM_scale=MIN((mauza(c)%saltM/16.0),1.0)
	  END DO
	  polder(p)%sL(t)=MIN((polder(p)%sL(t)/16.0),1.0)
	  pold_sL(i,tt)= polder(p)%sL(t)
	  !********************************************************************************************************

	  !****CROP PRODUCTION*************************************************************************************
	  !Calculate the impact of flooding, waterlogging, and salinity on crop production
      !
      !Calculate crop production based on salinity and waterlogging for each crop type and each farm holding type at a mauza scale
	  CALL agri_mauza_module(mauza(:)%saltM_scale, mauza(:)%waterLgM,nummauzaP,nummauzaM, nyear,t, ncrop, crop_yield(:,t),&
	      area_hh(:), crop_ratio(:),  Prodmax_m_temp(:,:,:,t))
		 
	  Prodmax_m_tempABS(:,:,:,t) = ABS(Prodmax_m_temp(:,:,:,t))  !Convert and -9999 flags to 99999 where no crop production for a particular crop
  	  !Calculate crop production based on any flooding  	
	  DO c= 1, nummauzaP    
        DO k=1, size(area_hh)
	      DO j=1, ncrop	
	        IF (polder(p)%mauza_frfl(c)<0.01) THEN  !not flooded (i.e. if =0)
		      mauza_ProdFlood(c,j,k,t)=Prodmax_m_tempABS(c,j,k,t)
		    ELSE IF (polder(p)%mauza_frfl(c)>0.01) THEN !flooded	(i.e. if =1) 		  
		      IF (Prodmax_m_tempABS(c,j,k,t)>9000) THEN !outside season for that crop
		        mauza_ProdFlood(c,j,k,t)=9999
		      ELSE
		        mauza_ProdFlood(c,j,k,t)=0.0    !Crop production for that mauza becomes zero due to flooding
		      END IF			  	  
		    END IF
		    mauza_ProdYearMin(c,j,k,y)=min(mauza_ProdYearMin(c,j,k,y), mauza_ProdFlood(c,j,k,t))	  !This takes the minimum monthly yield to account for any flooding which may have occurred
	      END DO
	    END DO
	  END DO	  
      !Calculate the crop income at the end of each year (m=12)
	  IF (m==12) THEN
	    DO k=1, size(area_hh)
	      DO j=1,ncrop	
		    DO c=1, nummauzaP
		      IF(mauza_ProdYearMin(c,j,k,y)<9999) THEN
		        mauza_CropIncome(c,j,k,y)=mauza_ProdYearMin(c,j,k,y)*price(j)
		        mauza_income_HHYr(c,k,y)=mauza_income_HHYr(c,k,y) + mauza_CropIncome(c,j,k,y)	!Sum income across all crops for a HH and muaza			
			    SumMauza_CropIncome(j,k,y)=SumMauza_CropIncome(j,k,y)+REAL(polder(p)%mauzaHoldings(c,k))*mauza_CropIncome(c,j,k,y)  !Sum income across all mauzas
			    mauza_ProdYearPoldTotHH(c,j,k,y)=REAL(polder(p)%mauzaHoldings(c,k))*mauza_ProdYearMin(c,j,k,y) !Total production for all holdings by holding type for each crop, year, and mauza			    	  
			  END IF				  
		    END DO !End c=1, nummauza	
		  END DO
	    END DO
	    polder(p)%ProdBy_CropHHYr(:,:,y)= SUM(mauza_ProdYearPoldTotHH(:,:,:,y), DIM=1) !Polder production across all mauzas
	    polder(p)%ProdBy_HHYr(:,y)=SUM(polder(p)%ProdBy_CropHHYr(:,:,y),DIM=1) !Polder production across all mauzas and all crops
	    polder(p)%ProdBy_CropYr(:,y)=SUM(polder(p)%ProdBy_CropHHYr(:,:,y),DIM=2)    !Sum of crop production across all households by polder (for individual crops and each year).  Need in this format for later calcs 	 	
	    pold_ProdBy_CropYr(i,:,y)=polder(p)%ProdBy_CropYr(:,y)  !Set as array for calc of yields, same as ProdBy_CropYr minus warm up. 
	    polder(p)%TotProdBy_Yr(y)=SUM(polder(p)%ProdBy_HHYr(:,y),DIM=1) !Total production for a year for each polder
	    PropMaxProd(i,y) = polder(p)%TotProdBy_Yr(y)/SumMaxCropProd(i)  !Proportion of actual vs maximum crop production each year
	    !Calculate totals based on number of holdings within a polder
		DO k=1, size(area_hh)    
		  DO c=1, nummauzaP			    
		    mauza_totincomeHH(c,k,y)= REAL(polder(p)%mauzaHoldings(c,k))*mauza_income_HHYr(c,k,y) !Total income per mauza for each HH type and each year
		    mauza_totincome(c,y)= mauza_totincome(c,y) + mauza_totincomeHH(c,k,y)  !Total income per mauza for each year
		  END DO		
		END DO
        !Sum crop income for each mauza
		 DO c=1, nummauzaP
		   mauzaAllYears_income(c)=mauzaAllYears_income(c)+ mauza_totincome(c,y) !sum of crop income per mauza
		   poldermauzaAllYears_income(i,c) = mauzaAllYears_income(c)  !sum of crop income per polder and mauza
         END DO		

		    IF (opt==2) THEN	!Test
		      WRITE(840+p,*)p, mauza_totincome(:,y)  !turn off for optimisation EJB
		    END IF   !if opt==2
		
		
		y=y+1
		m=0
	  END IF  !End if m==12
	  
	  m=m+1  !m is month 1 to 12
	  	   	    
    END DO   !End DO t=t1, 12*nyear
  
    !Calculate crop income, production and yield statistics (totals, averages)
    CropIncBy_CropHH(i,:,:)=SUM(SumMauza_CropIncome(:,:,:), DIM=3)  !total crop income for each crop and household type
    CropInc(i,:,:)=SUM(SumMauza_CropIncome(:,:,:), DIM=1) !sums all crops for income by HH and year
	pold_avProdBy_crop(i,:)=SUM(polder(p)%ProdBy_CropYr(:,:),DIM=2)/(REAL(numevents)) !Average crop production over simulation by crop: summed over all years (shown for individual crops) and divided by number of years
	DO j=1,ncrop
	  pold_avYieldBy_crop(i,j)=pold_avProdBy_crop(i,j)/(CropAreaBy_Crop(i,j))  !Average crop yield for each crop and each polder over entire simulation
	  pold_avYieldBy_CropYr(i,j,:)= pold_ProdBy_CropYr(i,j,:)/(CropAreaBy_Crop(i,j)) !Average crop yield for each crop and each polder per year
	END DO

	AvIncomePerYrBy_HH(i,:)= SUM(CropInc(i,:,:), DIM=2)/REAL(numevents)  !Average income for each polder and each HH type
	SumHoldingsHH(:)=SUM(polder(p)%mauzaHoldings(:,:),DIM=1) !Total number of holdings by holding type
	AvIncomePerYrBy_HH_numhold(i,:)= AvIncomePerYrBy_HH(i,:)/REAL(SumHoldingsHH(:)) !Divides by number of holdings per HH type by polder

	!Calculate damage as a percentage of the maximum possible yield given no floods
	DO k=1,4
	  Damage(i,k)= (REAL(maxinc(i,k))-AvIncomePerYrBy_HH(i,k))/(REAL(maxinc(i,k)))
    END DO
			
  
    IF (opt==2) THEN !Write the total income per mauza for each polder 
	  WRITE(filename,'(a,i4.4,a)') "Outputs/MauzaIncomeAllYearsAllHH",p,".txt"
	  OPEN(unit=750+p, file=filename, status='unknown')	
      DO c=1, nummauzaP	  
	    WRITE(750+p,*) polder(p)%mauzaID(c), mauzaAllYears_income(c) 	
	  END DO

	  CLOSE (unit=750+p)
	  CLOSE (unit=840+p)
	  CLOSE (unit=590+p)
	  CLOSE (unit=401+p)
    END IF
  
   
  END DO   !End DO p=1, numpolder
  	     
  SDamage(:)=SUM(Damage(:,:), DIM=1) !Sum of damage across polders for each HH type (SDamage(4))
  !*******************************END MAIN PROCESS***********************************************************
 
  
!************WRITE OUTPUTS TO FILE***************************************************************************  
  IF(opt==2) THEN !Write to file only if not using optimisation	 
    !*****FLOODING & EMBANKMENT OUTPUTS*****
    OPEN(unit=54, file='Outputs/reliability.txt', status='unknown')  !Embankment reliability by month
	OPEN (unit=40, file='Outputs/floodeventAllP.txt', status='unknown')  !Flood events external to polders by month
	OPEN (unit=41, file='Outputs/fractionfloodedAllP.txt', status='unknown')  !Internal fraction flooded by month
    OPEN (unit = 115, file = 'Outputs/MaxMauzaCropFloodedYearAllP.txt', status='unknown') !Crop area flooded by year
	WRITE(54,'(6(3x,a))') '29', '30', '32', '33', '43/1', '54'
	   
	DO m=1,12*nyear  
	  WRITE(54,'(6(3x,f8.4))')polder(31)%eR(m),polder(32)%eR(m),polder(35)%eR(m),polder(36)%eR(m),polder(65)%eR(m),polder(86)%eR(m)
	  WRITE(41,'(6(3x,f8.4))')polder(31)%fr_fl(m),polder(32)%fr_fl(m),polder(35)%fr_fl(m),polder(36)%fr_fl(m),polder(65)%fr_fl(m), & 
           polder(86)%fr_fl(m)
	END DO 
	DO m=1,12*numevents
	  WRITE(40,'(6(3x,f8.4))') (floodp(m,p), p=1,spold)
	END DO   	   
	DO y=1, numevents
      WRITE(115,'(i4,6(3x,4(3x,f14.4)))')y,pold_maxCropFlood(:,y)		
	END DO  !End numevents	   	   
	CLOSE (unit=54)     	   
	CLOSE (unit=40)
	CLOSE (unit=41)  
	CLOSE (unit=115)
	 
	!****SALINITY AND WATERLOGGING*****
  	OPEN(unit=53, file='Outputs/salinitylevel.txt', status='unknown')     !Salinity by month
	OPEN (unit=61, file='Outputs/PolderWaterlog.txt', status='unknown')   !Waterlogging by month	   
	DO m=1,12*nyear
	  WRITE(53,'(6(3x,f8.4))') polder(31)%sL(m),polder(32)%sL(m),polder(35)%sL(m),polder(36)%sL(m),polder(65)%sL(m), &
		     polder(86)%sL(m)    
	  WRITE(61,'(6(3x,f8.4))') (wlg_pold(p,m), p=1,spold)
	END DO      
    CLOSE (unit=53)
	CLOSE (unit=61)
	   
	!****AGRICULTURAL PRODUCTION******
    OPEN (unit=113, file='Outputs/TotProdBy_Yr.txt', status='unknown') !Total production for a year for each polder
	OPEN (unit=119, file='Outputs/ProdBy_HHYr.txt', status='unknown')  !Polder production across all mauzas and all crops
	OPEN (unit=120, file='Outputs/avYieldBy_CropYr.txt', status='unknown')  !Average crop yield for each crop and each polder per year
	OPEN (unit=121, file='Outputs/PropProdBy_PoldYr.txt', status='unknown')  !Proportion of actual vs maximum crop production each year
	WRITE(119,'(a,3x, 24(3x,a))') 'Year','P29HH1','P29HH2', 'P29HH3','P29HH4','P30HH1','P30HH2', 'P30HH3','P30HH4',&
		    'P32HH1','P32HH2', 'P32HH3','P32HH4','P33HH1','P33HH2', 'P33HH3','P33HH4', &
		    'P43HH1','P43HH2', 'P43HH3','P43HH4','P54HH1','P54HH2', 'P54HH3','P54HH4'						  
	WRITE(120,'(a,3x, 3(3x,a))') 'Year','Polder','Crop','Yield'	!EJB
	DO y=1, numevents
      WRITE(113,'(i4,3x,6(3x,f12.4))') y,polder(31)%TotProdBy_Yr(y),polder(32)%TotProdBy_Yr(y),polder(35)%TotProdBy_Yr(y), & 
		     polder(36)%TotProdBy_Yr(y),polder(65)%TotProdBy_Yr(y),polder(86)%TotProdBy_Yr(y)
	  WRITE(119,'(i4,6(3x,4(3x,f14.4)))') y,(polder(31)%ProdBy_HHYr(k,y), k=1,4),(polder(32)%ProdBy_HHYr(k,y),k=1,4), & 
		    (polder(35)%ProdBy_HHYr(k,y), k=1,4), (polder(36)%ProdBy_HHYr(k,y),k=1,4),(polder(65)%ProdBy_HHYr(k,y), k=1,4),&
            (polder(86)%ProdBy_HHYr(k,y),k=1,4)
	  WRITE(121, '(i4,6(3x, f8.5))') y,(PropMaxProd(i,y),i=1,6)	    
	  DO i=1,spold  !EJB - loop
	    DO j=1,ncrop
		  WRITE(120,'(3(i4,3x),f18.7)') y,i,j,pold_avYieldBy_CropYr(i,j,y)
		END DO !End ncrop
	  END DO	!End spold				
	END DO  !End numevents	   	   
	CLOSE (unit=113)
    CLOSE (unit=119)
    CLOSE (unit=120) !EJB
    CLOSE (unit=121)
	
	!****CROP INCOME *********
	OPEN (unit=114, file='Outputs/CropIncomeAllMauzaAllCrops.txt', status='unknown') !crop income for each holding time for each year, summed across crops and mauzas
	OPEN (unit=112, file='Outputs/CropIncomeBy_CropHH.txt', status='unknown') !total crop income for each crop and household type 
	write(114,'(a,3x, 24(3x,a))') 'Year','P29HH1','P29HH2', 'P29HH3','P29HH4','P30HH1','P30HH2', 'P30HH3','P30HH4',&
		  'P32HH1','P32HH2', 'P32HH3','P32HH4','P33HH1','P33HH2', 'P33HH3','P33HH4', &
		  'P43HH1','P43HH2', 'P43HH3','P43HH4','P54HH1','P54HH2', 'P54HH3','P54HH4'	
	write(112,'(a,3x, 24(3x,a))') 'Crop','P29HH1','P29HH2', 'P29HH3','P29HH4','P30HH1','P30HH2', 'P30HH3','P30HH4',&
		  'P32HH1','P32HH2', 'P32HH3','P32HH4','P33HH1','P33HH2', 'P33HH3','P33HH4', &
		  'P43HH1','P43HH2', 'P43HH3','P43HH4','P54HH1','P54HH2', 'P54HH3','P54HH4'						  
	DO y=1, numevents
	  WRITE(114,'(i4,6(3x,4(3x,f14.4)))') y, (CropInc(1,k,y),k=1,4),(CropInc(2,k,y),k=1,4),(CropInc(3,k,y),k=1,4), &
		    (CropInc(4,k,y),k=1,4),(CropInc(5,k,y),k=1,4),(CropInc(6,k,y),k=1,4)		
	END DO  !End numevents
	DO j=1, ncrop
      WRITE(112,'(i4,6(3x,4(3x,f16.2)))') j, (CropIncBy_CropHH(1,j,k), k=1,4),(CropIncBy_CropHH(2,j,k),k=1,4), &
		    (CropIncBy_CropHH(3,j,k),k=1,4),(CropIncBy_CropHH(4,j,k),k=1,4),(CropIncBy_CropHH(5,j,k),k=1,4), &
	        (CropIncBy_CropHH(6,j,k),k=1,4)
	  WRITE(filename,'(a,i4.4,a)') "Outputs/cropyield",j,".txt"
	  OPEN(unit=1000+j, file=filename, status='unknown') !Polder production across all mauzas
	  write(1000+j,'(a,3x, 24(3x,a))') 'Year','P29HH1','P29HH2', 'P29HH3','P29HH4','P30HH1','P30HH2', 'P30HH3','P30HH4',&
		     'P32HH1','P32HH2', 'P32HH3','P32HH4','P33HH1','P33HH2', 'P33HH3','P33HH4', &
		     'P43HH1','P43HH2', 'P43HH3','P43HH4','P54HH1','P54HH2', 'P54HH3','P54HH4'
	  DO y=1,numevents
	    WRITE(1000+j,'(i4,6(3x,4(3x,f14.4)))') y,(polder(31)%ProdBy_CropHHYr(j,k,y), k=1,4), &
		      (polder(32)%ProdBy_CropHHYr(j,k,y), k=1,4),&
			  (polder(35)%ProdBy_CropHHYr(j,k,y), k=1,4),(polder(36)%ProdBy_CropHHYr(j,k,y), k=1,4), &
			  (polder(65)%ProdBy_CropHHYr(j,k,y), k=1,4),(polder(86)%ProdBy_CropHHYr(j,k,y), k=1,4)
	  END DO
	  CLOSE (unit=1000+j)
	END DO	
	CLOSE (unit=112)
    CLOSE (unit=114)
  END IF  !End IF (opt==2)
     

	 DO i=1, spold
	   p=poldnum(i)
	   nummauzaP=nummauzaAll(i)
	 END DO
  !**********************************************************************************************************
 
  !********ERROR CHECKS**************************************************************************************
	!Ensure MCruns set to 1
  IF (ck==1 .AND. opt==2 .AND. inv==0) THEN
    OPEN (unit= 411, file='Outputs/ErrorCheck/fractionfloodedAllP.txt', status='old') !Fraction Flooded
    DO i=1,12*nyear
	  READ(411,*) poldfrfl_Check (i,1:6)
	  DO j=1,spold
	    p=poldnum(j)
		Check=ABS(polder(p)%fr_fl(i) - poldfrfl_Check(i,j))
		IF (Check>0.0001) THEN
     	  WRITE(*,*) 'polder(p)%fr_fl does not match base for polder', j, 'at time', i
		  WRITE(*,*) 'Expected', poldfrfl_Check (i,j), 'Actual', polder(p)%fr_fl(i)
		  PAUSE
		END IF
	  END DO
	END DO
	CLOSE (unit=411)
	OPEN (unit=611, file='Outputs/ErrorCheck/PolderWaterlog.txt', status='old')   !Waterlogging by month
	DO i=1,12*nyear
	  READ(611,*) wlg_pold_Check(1:spold,i)
	  DO j=1, spold
		Check=ABS(wlg_pold(j,i)-wlg_pold_Check(j,i))
		IF (Check>0.0001) THEN
     	  WRITE(*,*) 'wlg_pold does not match base for polder', j, 'at time', i
		  WRITE(*,*) 'Expected', wlg_pold_Check(j,i), 'Actual', wlg_pold(j,i)
		  PAUSE
		END IF 
	  END DO
	END DO
	CLOSE (unit=611)
	OPEN (unit=531, file='Outputs/ErrorCheck/salinitylevel.txt', status='old') !Salinity level by month
	DO i=1,12*nyear
	  READ(531,*) saltCheck (1:spold,i)
	  DO j=1,spold
		p=poldnum(j)
		Check=ABS(polder(p)%sL(i) - saltCheck (j,i))
		IF (Check>0.0001) THEN
     	  WRITE(*,*) 'polder salinity (sL) does not match base for polder', j, 'at time', i
		  WRITE(*,*) 'Expected', saltCheck(j,i), 'Actual', polder(p)%sL(i)
		  PAUSE
		END IF 
      END DO
	END DO
	CLOSE (unit=631)
	OPEN(unit=541, file='Outputs/ErrorCheck/reliability.txt', status='old') !Embankment reliability by month
	READ(541,*)
	DO i=1, 12*nyear
	  READ(541,*)embCheck(1:spold,i)
	  DO j=1, spold
		p=poldnum(j)
		Check=ABS(polder(p)%eR(i)-embCheck(j,i))
		IF (Check>0.0001) THEN
     	  WRITE(*,*) 'embankment reliability (eR) does not match base for polder', j, 'at time', i
		  WRITE(*,*) 'Expected', embCheck(j,i), 'Actual', polder(p)%eR(i)
		  PAUSE
		END IF 		   
	  END DO
	END DO
	CLOSE (unit=541)
  END IF
  !********END ERROR CHECKS**********************************************************************************
  
  !*******OPTIMISATION OBJECTIVES FOR RELIABILITY-BASED******************************************************
  !Not currently in use, need revising objectives if used. May need moving to Program section.
  IF (inv/= 1) THEN  !Reliability-based
    OPEN(unit=57, file='Obj001.txt', status='unknown')
	WRITE (57,*) TotInvest
	WRITE (57,*) confail
  END IF
  !**********************************************************************************************************
  
END SUBROUTINE MainProcess
END MODULE RiskModel 
!************************************************************************************************************
  
 !***********************************************************************************************************
 !
 !                                                MAIN PROGRAM
 !
 !*********************************************************************************************************** 
 PROGRAM CoastalRisk
 USE RiskModel
 IMPLICIT NONE
 
   !***********PARAMETERS***********************************************
   !m=MC run number; k=HH number; t=timestep;  j=number crops; y= year; i=polder number
   !MCruns = total number of stochastic runs defined by the user
   !cMax = the maximum number of mauzas in a particular polder
   !MCnummauzaAll = array containing max number of mauzas for all six polders
   !MCcostfail = total investment cost failure for each MC run
   !TotMCcostfail = investment cost failure summed over all MC runs
   INTEGER :: m,k,t,i,j,y,c 
   INTEGER :: MCruns,cMax, MCnummauzaAll(6), TotMCcostfail 
   INTEGER, ALLOCATABLE :: MCcostfail(:)
   
   !*****Flood, salinity, waterlogging**************
   !MCpold_sL = monthly time series of salinity for each polder and each MC run
   !MCfloodp=floodevent by month for each polder; MCfloodm=floodevent per mauza (max 62) for each polder
   !MCpold_maxCropFlood = maximum crop area flooded for each polder and each year
   REAL, ALLOCATABLE :: MCpold_sL(:,:,:), MCfloodp(:,:,:), MCfloodm(:,:,:) 
   REAL, ALLOCATABLE :: MCpold_maxCropFlood(:,:,:)
   
   !*****Embankments and Investment*****************
   !MCpold_eR = monthly time series of embankment reliability for each polder and each MC run
   REAL, ALLOCATABLE :: MCpold_eR(:,:,:)
   
   !*****Agriculture********************************
   !TotMCDam = Running sum of damage for each HH type across all MC runs
   !AvMCDamage = average damage for each HH type across all MC runs
   !TotAvMCDamage = sum of average damage across HH types for all MC runs: ****OBJECTIVE FUNCTION****
   !HH1AvMCDamage = average damage for HH1 only (subsistence farmers) across all MC runs: ****OBJECTIVE FUNCTION****
   !MCpoldermauzaAllYears_income=total crop income per polder and per mauza (max 62 mauzas) for each MC run
   !TotMCpoldmauza_income = sum of crop income per polder and mauza, summed over MC runs
   !AvMCpoldmauza_income = Average income per mauza and polder across MC runs
   !MCDamage = crop damage by polder and HH type for all MC runs
   !MCSDamage = sum of crop damage by HH type across all polders for each MC run
   !MCAvIncomePerYrBy_HH= average crop income per polder and HH for each MC run
   !MCAvIncomePerYrBy_HH_numhold=Divides by number of holdings per HH type by polder, shown for each  MC run
   !MCAvYield = average crop yield for each polder and each crop for each MC run  
   REAL :: TotMCDam(4), AvMCDamage(4), TotAvMCDamage, HH1AvMCDamage
   REAL :: TotMCpoldmauza_income(6,62),AvMCpoldmauza_income(6,62)
   REAL, ALLOCATABLE :: MCpoldermauzaAllYears_income(:,:,:)
   REAL, ALLOCATABLE :: MCDamage(:,:,:), MCSDamage(:,:)
   REAL, ALLOCATABLE :: MCAvIncomePerYrBy_HH(:,:,:), MCAvIncomePerYrBy_HH_numhold(:,:,:) 
   REAL, ALLOCATABLE :: MCAvYield(:,:,:) 
   
   !*****Checks************************************
   REAL :: MCAvIncomePerYrBy_HH_Check(4,6), Check, MCDamagePerYrBy_HH_Check(4,6)   

   !*****Initialise Parameters*********************   
   MCruns=1   !*****USER INPUT********* Select the number of MC runs**********  
   TotMCDam=0d0
   
   !*****Allocate arrays***************************
   ALLOCATE (MCcostfail(MCruns), MCpold_sL(MCruns,6,720), MCfloodp(MCruns,720,6))
   ALLOCATE (MCfloodm(MCruns,6,62), MCpold_maxCropFlood(MCruns,6,60), MCpold_eR(MCruns,6,720))  
   ALLOCATE (MCpoldermauzaAllYears_income(MCruns,6,62))
   ALLOCATE (MCDamage(MCruns,6,4),MCSDamage(MCruns,4))
   ALLOCATE (MCAvIncomePerYrBy_HH(MCruns,6,4), MCAvIncomePerYrBy_HH_numhold(MCruns,6,4))
   ALLOCATE (MCAvYield(MCruns,6,13)) 

   !*****Initialise Parameters*********************
   MCpold_eR=0d0
   MCpold_sL=0d0
   MCpold_maxCropFlood=0d0
   TotMCpoldmauza_income=0d0  
   AvMCpoldmauza_income=0d0    
   MCnummauzaAll = [62, 37, 6, 9, 37, 21] !Number of mauzas in each polder
   !***********************************************************************************

   !***********************************************************************************
   !                              MAIN CALCULATIONS
   !***********************************************************************************
   DO m=1,MCruns
     IF (opt==2) WRITE(*,*) "Run number:  ", m
     CALL MainProcess   !Call main module with all subroutines
	 MCDamage(m,:,:)=Damage(:,:)
	 MCSDamage(m,:)= SDamage(:)
	 TotMCDam(:)=TotMCDam(:)+ MCSDamage(m,:)
	 MCpold_eR(m,:,:)=pold_eR(:,:)
	 MCpold_sL(m,:,:)=pold_sL(:,:)
	 MCpold_maxCropFlood(m,:,:)=pold_maxCropFlood(:,:)
	 MCAvYield(m,:,:)=pold_avYieldBy_crop(:,:)   
	 MCAvIncomePerYrBy_HH(m,:,:)= AvIncomePerYrBy_HH(:,:)
	 MCAvIncomePerYrBy_HH_numhold(m,:,:)= AvIncomePerYrBy_HH_numhold(:,:)
	 MCpoldermauzaAllYears_income(m,:,:) = poldermauzaAllYears_income(:,:) 
	 MCfloodp(m,:,:)=floodp(:,:)
	 MCfloodm(m,:,:)=floodm(:,:) 
	 MCcostfail(m)=costfail
   END DO
	 
   !Calculate averages and totals
   AvMCDamage(:)=TotMCDam(:)/(REAL(MCruns)) !average damage for each HH type across all MC runs
   TotAVMCDamage=SUM(AvMCDamage(:),DIM=1)   !sum of average damage across HH types for all MC runs - objective function
   HH1AvMCDamage=AVMCDamage(1)              !average damage for HH1 only (subsistence farmers) across all MC runs - objective function
   TotMCcostfail=SUM(MCcostfail(:),DIM=1)   !investment cost failure summed over all MC runs
   TotMCpoldmauza_income(:,:)=SUM(MCpoldermauzaAllYears_income(:,:,:),DIM=1) !Sum income per mauza and polder over all MC runs   
   AvMCpoldmauza_income(:,:)=TotMCpoldmauza_income(:,:)/(REAL(MCruns)) !Average income per mauza and polder across MC runs
  
   !Error Check - make sure set to reliability-based with no investment
	IF (ck==1 .AND. MCruns==1 .AND. opt==2 .AND. inv==0) THEN
	   OPEN (unit= 91, file='Outputs/MC/ErrorCheck/MCAvIncomePerYrBy_HH.txt', status='old')
	   READ(91,*)
	   DO k=1,4
		 READ(91,*) m,j,MCAvIncomePerYrBy_HH_Check(k,1:6)
		 DO i=1,6
		   Check=ABS(MCAvIncomePerYrBy_HH (1,i,k) - MCAvIncomePerYrBy_HH_Check(k,i))
		   IF (Check>0.000001) THEN
     		 WRITE(*,*) 'MCAvIncomePerYrBy_HH does not match base'
			 WRITE(*,*) 'Expected', MCAvIncomePerYrBy_HH_Check(k,i), 'Actual', MCAvIncomePerYrBy_HH (1,i,k)
		     PAUSE
		   END IF
		 END DO
	   END DO
	   CLOSE (unit=91)
	   OPEN(unit= 511, file='Outputs/MC/ErrorCheck/MCDamagePerYrBy_HH.txt', status='old')
	   READ(511,*)
	   DO k=1,4
	     READ(511,*) m,j, MCDamagePerYrBy_HH_Check(k,1:6)
		 DO i=1,6
		   Check=ABS(MCDamage(1,i,k) - MCDamagePerYrBy_HH_Check(k,i))
		   IF (Check>0.0001) THEN
     		 WRITE(*,*) 'MCDamage does not match base'
			 WRITE(*,*) 'Expected', MCDamagePerYrBy_HH_Check(k,i), 'Actual', MCDamage(1,i,k)
		     PAUSE
		   END IF
		 END DO
	   END DO
	   CLOSE(unit=511)
    END IF
	   
   !********WRITE OBJECTIVE FUNCTIONS TO OUTPUT***************************************
   IF (inv==1) THEN   !Time-based investment
     OPEN(unit=57, file='Obj001.txt', status='unknown')
	 !Select EITHER TotAvMCDamage OR HH1AvMCDamage depending on if optimising for all HH types or HH1 only (comment other one out)
	 !If Single objective, also comment out TInvCost and just optimise for minimum damage
	 !If single objective with a cost constraint, use TotMCcostfail and either TotAvMCDamage or HH1AvMCDamage (set 'Budget' parameter)
	 WRITE(57,*) TotAVMCDamage
	 WRITE(57,*) HH1AvMCDamage
	 WRITE(57,*) TInvCost  
	! WRITE (57,*) TotMCcostfail
	 CLOSE (unit=57)
   END IF
   !**********************************************************************************
   
   !********WRITE OUTPUTS TO FILE*****************************************************
   IF (opt==2) THEN !Test
     OPEN (unit=49, file='Outputs/MC/floodeventMC.txt', status='unknown')
	 OPEN (unit=58, file='Outputs/MC/MCAvIncomePerYrBy_HH.txt', status='unknown')
	 OPEN (unit=582, file='Outputs/MC/MCAvIncomePerYrBy_HH_numhold.txt', status='unknown')
	 OPEN (unit=51, file='Outputs/MC/MCDamagePerYrBy_HH.txt', status='unknown')
	 OPEN (unit=583, file='Outputs/MC/MCReliabilityBy_Polder.txt', status='unknown')
	 OPEN (unit=584, file='Outputs/MC/MCSalinityBy_Polder.txt', status='unknown')
	 OPEN (unit=585, file='Outputs/MC/MCMaxCropFloodBy_Polder.txt', status='unknown')
	 OPEN (unit=586, file='Outputs/MC/MCAvYieldBy_PolderCrop.txt', status='unknown')  !EJB
	 OPEN (unit=587, file='Outputs/MC/MCAvpoldmauza_income.txt', status='unknown')  !EJB
	 OPEN (unit=588, file='Outputs/MC/MCMauzaFloodTotMonths.txt', status='unknown')  !EJB 28/10/2020

	 WRITE(49, '(a,3x, 7(3x,a))') 'MCrun', 'Time', 'P29','P30','P32','P33','P43','P54'
	 WRITE(58, '(a,3x, 7(3x,a))') 'MCrun', 'HH', 'P29','P30','P32','P33','P43','P54'
	 WRITE(582, '(a,3x, 7(3x,a))') 'MCrun', 'HH', 'P29','P30','P32','P33','P43','P54'
	 WRITE(51, '(a,3x, 7(3x,a))') 'MCrun', 'HH', 'P29','P30','P32','P33','P43','P54'
	 WRITE(583, '(a,3x, 7(3x,a))') 'MCrun', 'Month', 'P29','P30','P32','P33','P43','P54'
	 WRITE(584, '(a,3x, 7(3x,a))') 'MCrun', 'Month', 'P29','P30','P32','P33','P43','P54'
	 WRITE(585, '(a,3x, 7(3x,a))') 'MCrun', 'Month', 'P29','P30','P32','P33','P43','P54'
	 WRITE(586, '(a,3x, 3(3x,a))') 'MCrun', 'Polder', 'Crop','Production'  !EJB
	 WRITE(587, '(a,3x, 2(3x,a))') 'Polder', 'Mauza','MCAvIncome'  !EJB
	 WRITE(588, '(a,3x, 3(3x,a))') 'MCrun', 'Polder', 'Mauza', 'TotalMonthsFlooded' !EJB 28/10/2020
	 DO m=1, MCruns
	   DO t=1,720
		    WRITE(49,'(i4,3x,i4,6(3x,4(3x,f18.4)))') m, t, MCfloodp(m,t,1),MCfloodp(m,t,2),MCfloodp(m,t,3),&
		    MCfloodp(m,t,4),MCfloodp(m,t,5),MCfloodp(m,t,6)	
            WRITE(583,'(i4,3x,i4,6(3x,4(3x,f18.4)))') m, t, MCpold_eR(m,1,t),MCpold_eR(m,2,t),MCpold_eR(m,3,t),&
		    MCpold_eR(m,4,t),MCpold_eR(m,5,t),MCpold_eR(m,6,t)	
            WRITE(584,'(i4,3x,i4,6(3x,4(3x,f18.4)))') m, t, MCpold_sL(m,1,t),MCpold_sL(m,2,t),MCpold_sL(m,3,t),&
		    MCpold_sL(m,4,t),MCpold_sL(m,5,t),MCpold_sL(m,6,t)			
	   END DO
	   
	   DO y=1,60
	     WRITE(585,'(i4,3x,i4,6(3x,4(3x,f18.4)))') m, y, MCpold_maxCropFlood(m,1,y),MCpold_maxCropFlood(m,2,y),&
	     MCpold_maxCropFlood(m,3,y), MCpold_maxCropFlood(m,4,y),MCpold_maxCropFlood(m,5,y),MCpold_maxCropFlood(m,6,y)
	   END DO
   
	   DO k=1,4
			WRITE(58,'(i4,3x,i4,6(3x,4(3x,f18.4)))') m, k, MCAvIncomePerYrBy_HH(m,1,k),MCAvIncomePerYrBy_HH(m,2,k),&
		      MCAvIncomePerYrBy_HH(m,3,k),MCAvIncomePerYrBy_HH(m,4,k),MCAvIncomePerYrBy_HH(m,5,k),&
			  MCAvIncomePerYrBy_HH(m,6,k)
			WRITE(582,'(i4,3x,i4,6(3x,4(3x,f18.4)))') m, k, MCAvIncomePerYrBy_HH_numhold(m,1,k),&
		      MCAvIncomePerYrBy_HH_numhold(m,2,k),MCAvIncomePerYrBy_HH_numhold(m,3,k),&
			  MCAvIncomePerYrBy_HH_numhold(m,4,k),MCAvIncomePerYrBy_HH_numhold(m,5,k), &
			  MCAvIncomePerYrBy_HH_numhold(m,6,k)
		 WRITE(51,'(i4,3x,i4,6(3x,4(3x,f18.4)))') m, k, MCDamage(m,1,k),MCDamage(m,2,k),MCDamage(m,3,k),&
		    MCDamage(m,4,k),MCDamage(m,5,k),MCDamage(m,6,k)			
       END DO
	   
	   DO i=1,6  
	     DO j=1,13
		   WRITE(586,'(3(i4,3x),f18.7)') m,i,j,MCAvYield(m,i,j)
		 END DO
	   END DO  
	   DO i=1,6 
         cMax=MCnummauzaAll(i)
         DO c=1,cMax
		   WRITE(588,'(3(i4,3x),f18.7)')m,i,c,MCfloodm(m,i,c) 
	     END DO
	   END DO	   	   
	 END DO  !m=1 to MC runs
	 
	 DO i=1,6
     cMax=MCnummauzaAll(i)
       DO c=1,cMax
	     WRITE(587,'(2(i4,3x),f18.7)') i,c,AvMCpoldmauza_income(i,c)
	   END DO
	 END DO
	 
	 CLOSE (unit=49)
	 CLOSE (unit=58)
	 CLOSE (unit=582)
	 CLOSE (unit=51)
	 CLOSE (unit=583)
	 CLOSE (unit=584)
	 CLOSE (unit=585)
	 CLOSE (unit=586)  
	 CLOSE (unit=587)  
	 CLOSE (unit=588) 
   END IF  !IF opt==2
   !**********************************************************************
  
 END PROGRAM CoastalRisk
 !*************************************************************************************************

    !*************************************************************************
	!                          FUNCTION POLDFRFLOOD
	!*************************************************************************
	!Calculates the fraction of each polder flooded based on the magnitude of the flood and the embankment reliability
	!Flood recession is also calculated
	REAL FUNCTION poldfrflood (timeflood, poldER, poldflood, tmonth, numyears)
	  IMPLICIT NONE
	  REAL, INTENT(IN):: poldER (12*numyears), poldflood (12*numyears)  !embankment reliability; magnitude of flood event for a polder
	  INTEGER, INTENT(IN):: timeflood(12*numyears), tmonth, numyears  !current monthly timestep, polder number, total number of years
	  REAL :: A_F
		   !Flooding of polder area and retreat of waters 	 
	  IF (timeflood(tmonth)>0) THEN ! a flood happens
        A_F = MIN(1.0, 1.0/((MAX(poldER(tmonth),0.0001)))*(poldflood(tmonth))**3  ) ! F is the flood severity
        poldfrflood = A_F
      else if (timeflood(tmonth-1)>0) then 
        A_F = min(1.0, 1.0/((MAX(poldER(tmonth-1),0.0001)))*(poldflood(tmonth-1))**3  )! F is the flood severity
        poldfrflood = A_F*exp(-1.0/5.0)
      else if (timeflood(tmonth-2)>0) then 
        A_F = min( 1.0, 1.0/((MAX(poldER(tmonth-2),0.0001)))*(poldflood(tmonth-2))**3  )! F is the flood severity
        poldfrflood = A_F*exp(-2.0/5.0)
      else if (timeflood(tmonth-3) >0) then 
        A_F = min( 1.0, 1.0/((MAX(poldER(tmonth-3),0.0001)))*(poldflood(tmonth-3))**3  )! F is the flood severity
        poldfrflood = A_F*exp(-3.0/5.0)
      else if (timeflood(tmonth-4) >0) then 
        A_F = min( 1.0, 1.0/((MAX(poldER(tmonth-4),0.0001)))*(poldflood(tmonth-4))**3  )! F is the flood severity
        poldfrflood = A_F*exp(-4.0/5.0)
      else if (timeflood(tmonth-5) >0) then 
         A_F = min(1.0, 1.0/((MAX(poldER(tmonth-5),0.0001)))*(poldflood(tmonth-5))**3  )! F is the flood severity
         poldfrflood = A_F*exp(-5.0/5.0)
      else if (timeflood(tmonth-6)>0) then 
	 	A_F = min( 1.0, 1.0/((MAX(poldER(tmonth-6),0.0001)))*(poldflood(tmonth-6))**3  )! F is the flood severity
         poldfrflood = A_F*exp(-6.0/5.0)
      else
         poldfrflood =0
      end if	
    END FUNCTION poldfrflood
    !***************END FUNCTION*******************************************************************