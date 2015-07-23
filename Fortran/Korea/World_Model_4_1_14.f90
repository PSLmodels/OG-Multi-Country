
! ********************************************************************************************************************************** 
! *         Module With Global Declarations                                                                                        *
! **********************************************************************************************************************************
! Countries are always signed by the same number, i.e. 1=US, 2=EU, 3=Japan, 4=China, 5=India, 6=Russia
! FIRST_COUNTRY is the first country chosen for the simulation 
! LAST_COUNTRY is the last country chosen for the simulation
! If FIRST_COUNTRY is equal to LAST_COUNTRY, then the simulation is for the closed economy of the respective country
! NOTE: Not all combinations of countries is possible until now, only the following are available:
! US, EU, Japan, China, US-EU, EU-Japan, Japan-China, US-EU-Japan, EU-Japan-China, US-EU-Japan-China

      MODULE GLOBAL_DATA
           IMPLICIT NONE
           SAVE
           CHARACTER (LEN = 100), DIMENSION(3) :: INFILE_NAME
           CHARACTER (LEN = 100), DIMENSION(19) :: OUTFILE_NAME
           INTEGER, PARAMETER :: COUNTRIES = 6, FIRST_COUNTRY = 1, LAST_COUNTRY = 6, FIRST_FERTILITY_YEAR = 23, &
                & LAST_FERTILITY_YEAR = 45, FIRST_WORK_YEAR = 21, GENS = 90, LAST_EDUCATION_YEAR = 20, Y_CLASSES = 2, & 
                & YEARS = 300, FIRST_YEAR = 2008
           INTEGER :: FIRST_SOLUTION_YEAR, IRUN, TRANS_YEAR, ITER, NITER           
           INTEGER, DIMENSION(-1:YEARS, COUNTRIES) :: RETIREMENT_AGE
           REAL*8 :: ALP, ALPHA, DAMP, DAMPK, DAMPR, DEL, EXOG_NPOP, GAMMA, HOURS, RHO, SIGFIG, TECH, THETA
           REAL*8, DIMENSION(COUNTRIES) ::  GOVS
           REAL*8, DIMENSION(FIRST_COUNTRY:LAST_COUNTRY) :: POP_SCALE, TOTAL_FERTILITY_SS
           REAL*8, DIMENSION(-1:YEARS) :: AGG_ASSETS_WORLD, AP, DD_WORLD, FOREIGN_ASSETS_WORLD, RG, TRADE_BALANCE_WORLD, YY_WORLD, &
                & TOT_GOV_ENDOW_SHARE, TOT_GOV_ENDOW_REVENUE !Add extra endowment variables, national income variables
           REAL*8, DIMENSION(-1:YEARS) :: NET_ENDOW_INVEST !%$%$%$ Define net endow invest globally
           REAL*8, DIMENSION(-GENS:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: AVG_BIRTH_AGE, TOTAL_FERTILITY
           REAL*8, DIMENSION(Y_CLASSES, COUNTRIES) :: C_SHARE           
           REAL*8, DIMENSION(-1:YEARS, COUNTRIES) :: DEBT_LEVEL, DISABILITY_BENEFITS_IND, INHERITANCE_TAX_RATE, MU_1, MU_2_GOV, &
                & MU_2_TAX, MU_3, MU_4, CONTRIBUTION_CEILING, GOVERNMENT_DISCRETIONARY_SPENDING
           REAL*8, DIMENSION(-1:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: AGG_ASSETS, AGG_DISABILITY_TAX_RATE, AGG_HEALTH_TAX_RATE, &
                & AGG_PENSION_TAX_RATE, AVG_LABOR_EARNINGS, CAPITAL, CAPITAL_TAX_RATE, CC, CONSUMP_PRICE, CORP_TAX, DD, DEBT, &
                & DISABILITY_BENEFITS, EDUCATION_EXPENDITURES, FOREIGN_ASSETS, GOVERNMENT_EXPENDITURES, HEALTH_BENEFITS, INVEST, &
                & NPOP, PENSION_BENEFITS, POP_EFFICIENT, R, TRADE_BALANCE, YY, YY_0, DEFICIT, TRF, GOV_ENDOW_SHARE, TOTAL_EXPENDITURES !%%%Add Gov share, other national income stuff
           REAL*8, DIMENSION(0:LAST_EDUCATION_YEAR, FIRST_COUNTRY:LAST_COUNTRY) :: EDUCATION_EXPENDITURES_IND
           REAL*8, DIMENSION(0:GENS, 1:6) :: HEALTH_EXPENDITURES_PROFILE
           REAL*8, DIMENSION(0:GENS, COUNTRIES) :: MIGRANTS
           REAL*8, DIMENSION(Y_CLASSES) :: BETA                
           REAL*8, DIMENSION(-1:YEARS, Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) :: AGG_ASSETS_FOR_BEQUESTS, AVG_INDEXED_EARNINGS, &
                & C_SHARE_MIGRANTS, WAGE_INDEX_CLASS
           REAL*8, DIMENSION(-1:YEARS, Y_CLASSES, COUNTRIES) :: PRODUCTIVITY_C
           REAL*8, DIMENSION(-1:YEARS, Y_CLASSES + 1, FIRST_COUNTRY:LAST_COUNTRY) :: AGG_ASSETS_MIGRANTS, &
                & AGG_AVG_DISABILITY_TAX_RATE, AGG_AVG_HEALTH_TAX_RATE, AGG_AVG_PENSION_TAX_RATE, AGG_AVG_WAGE_TAX_RATE, &
                & AGG_MARG_DISABILITY_TAX_RATE, AGG_MARG_HEALTH_TAX_RATE, AGG_MARG_PENSION_TAX_RATE, AGG_MARG_WAGE_TAX_RATE, LABOR
           REAL*8, DIMENSION(0:GENS, -1:YEARS, COUNTRIES) :: MORTALITY
           REAL*8, DIMENSION(0:GENS, -1:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: HEALTH_BENEFITS_IND, SURVIVAL_PROBABILITY
           REAL*8, DIMENSION(1:GENS, -1:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: MIGRATION_SCALE
           REAL*8, DIMENSION(3, -1:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: ALEPH, BETH
           REAL*8, DIMENSION(6, -1:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: TAX_REVENUES
           REAL*8, DIMENSION(FIRST_FERTILITY_YEAR:LAST_FERTILITY_YEAR, -GENS:YEARS, COUNTRIES) :: FERTILITY
           REAL*8, DIMENSION(FIRST_FERTILITY_YEAR:LAST_FERTILITY_YEAR, COUNTRIES) :: FERTILITY_SS
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, -1:YEARS, COUNTRIES) :: DELTA
           REAL*8, DIMENSION(0:GENS, -1:YEARS, Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) :: PENSION_REPLACEMENT_RATE
           REAL*8, DIMENSION(0:GENS + 1, -1:YEARS, Y_CLASSES + 1, FIRST_COUNTRY:LAST_COUNTRY) :: POP
           REAL*8, DIMENSION(3, -1:YEARS, COUNTRIES) :: ENDOGENOUS_TAX_RATIO             
           REAL*8, DIMENSION(4, 2, -1:YEARS, COUNTRIES) :: TAX_RATE
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, -1:YEARS, Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) :: AGE_EFFICIENCY, ASSETS, &
                & AVG_DISABILITY_TAX_RATE, AVG_HEALTH_TAX_RATE, AVG_PENSION_TAX_RATE, AVG_WAGE_TAX_RATE, CONSUMP, CONSUMP_KIDS, &
                & LEISURE, MARG_DISABILITY_TAX_RATE, MARG_HEALTH_TAX_RATE, MARG_PENSION_TAX_RATE, MARG_WAGE_TAX_RATE, &
                & PENSION_BENEFITS_IND, SHADOW_WAGE, UTILITY, TRANSFER
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, -1:YEARS, Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY, FIRST_COUNTRY:LAST_COUNTRY) :: &
                & H_TRANSFER     
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, 0:YEARS, Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) :: UTILITY_0
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, Y_CLASSES, COUNTRIES) :: ASSETS_INITIAL_YEAR, BASE_ASSETS
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) :: ZXI     
           			!!!!!Create a flow of endowments in , PVENDOWMENT will be constructed from this $%$%$%
		   REAL*8, DIMENSION(-1:YEARS):: ENDOWMENT
           REAL*8, DIMENSION(-1:YEARS+1):: PVENDOWMENT
           REAL*8, DIMENSION(-1:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: GDP, NATIONAL_INCOME, COUNTRY_ENDOW_SHARE !THESE ARE NEW VARIABLES FOR CACULATING NATIONAL STUFF
      END MODULE GLOBAL_DATA

! ********************************************************************************************************************************** 
! *                                                       Main program.                                                            *
! **********************************************************************************************************************************

! ********************************************************************************************************************************** 
! *   NAME: PROGRAM DEMO                                                                                                           *
! *   PURPOSE: The main program unit.                                                                                              *
! **********************************************************************************************************************************

      PROGRAM WORLD_MODEL 
           USE GLOBAL_DATA
           IMPLICIT NONE
           IRUN = 0
!          NITER = 0 -> Baseline path only
!          NITER = 1 -> Baseline path and policy reform
           NITER = 1
!  
           CALL INITIALIZE
           CALL GET_INITIAL_STEADY_STATE
           CALL GET_TRANSITION           
           CALL RESULTS
           IF(NITER == 0) STOP
!             
           IRUN = 1
           CALL INITIALIZE_SCENARIO
           CALL GET_TRANSITION
           !CALL RESULTS
           CALL PRINT_HICKSIAN_EV
      END

! ********************************************************************************************************************************** 
! *    NAME: SUBROUTINE INITIALIZE                                                                                                 *
! *    PURPOSE: Set parameters for the initial base case                                                                           *
! **********************************************************************************************************************************

      SUBROUTINE INITIALIZE
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, IK, IM21, J, Y_CLASS, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8, DIMENSION(13) :: GBA
           REAL*8, DIMENSION(COUNTRIES) :: EDUCATION_EXPENDITURES_SCALE, HEALTH_EXPENDITURE_SCALE
           REAL*8, DIMENSION(0:LAST_EDUCATION_YEAR) :: EDUCATION_EXPENDITURES_PROFILE
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, -1:YEARS, COUNTRIES) :: OAP
           DATA EDUCATION_EXPENDITURES_PROFILE / 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 24.3513120D0, 47.256560D0, &
                &  49.4795920D0, 53.2609330D0, 53.5932940D0, 54.5947520D0, 53.8163270D0, 52.4825070D0, 51.147230D0, 50.9256560D0, &
                &  50.8148690D0, 49.5379010D0, 45.6005830D0, 33.5160350D0, 20.9854230D0 /
           INFILE_NAME(1) = 'POP_RUS_01_20_14.DAT' !(fertility rates for Russia adjusted (multyplied) for those not in 23-45 age group) 


		   INFILE_NAME(2) = 'ENDOWMENT_1_11_14.DAT' !$%$%$%$ Choose your Endowment File
           
           OUTFILE_NAME(1) = 'POPMOD1.OUT'
           OUTFILE_NAME(2) = 'POPMOD2.OUT'
           OUTFILE_NAME(3) = 'ERGEB1.OUT'
           OUTFILE_NAME(4) = 'ERGEB2.OUT'
           OUTFILE_NAME(5) = 'INITIAL_EQUILIBRIUM.OUT'
           OUTFILE_NAME(12) = 'BASELINE_DATA.DAT'
           OUTFILE_NAME(13) = 'RESULTS_USA.OUT'
           OUTFILE_NAME(14) = 'RESULTS_EU.OUT'
           OUTFILE_NAME(15) = 'RESULTS_JAPAN.OUT'
           OUTFILE_NAME(16) = 'RESULTS_CHINA.OUT'
           OUTFILE_NAME(17) = 'RESULTS_INDIA.OUT'
           OUTFILE_NAME(18) = 'RESULTS_RUSSIA.OUT'
           OUTFILE_NAME(19) = 'ENDOWMENT_DATA.OUT'

           
           OPEN (UNIT = 5, FILE = OUTFILE_NAME(3), STATUS = 'REPLACE')
!          Initialize Parameters
           SIGFIG = 0.0000001       
    
           DAMP = 0.7
           DAMPK = 0.05
           DAMPR = 0.2

!          POPULATION PARAMETERS
!          Gewichtung der Bevölkerung (z.B.China mit kleinerer Bevölkerung laufen lassen)
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                POP_SCALE(COUNTRY) = 1.0D0
                IF(COUNTRY == 1) POP_SCALE(COUNTRY) = 1.01D0
                IF(COUNTRY == 2) POP_SCALE(COUNTRY) = 0.81D0
                IF(COUNTRY == 3) POP_SCALE(COUNTRY) = 1.59D0
                IF(COUNTRY == 4) POP_SCALE(COUNTRY) = 0.98D0
                IF(COUNTRY == 5) POP_SCALE(COUNTRY) = 0.99D0  
                IF(COUNTRY == 6) POP_SCALE(COUNTRY) = 1.00D0 
           ENDDO
!          Anteil der Einkommensklassen an der Bevölkerung (nur für Jahr 0)
           C_SHARE(1, 1) = 0.70D0
           C_SHARE(2, 1) = 0.30D0
           C_SHARE(1, 2) = 0.70D0
           C_SHARE(2, 2) = 0.30D0
           C_SHARE(1, 3) = 0.70D0
           C_SHARE(2, 3) = 0.30D0
           C_SHARE(1, 4) = 0.75D0
           C_SHARE(2, 4) = 0.25D0
           C_SHARE(1, 5) = 0.75D0
           C_SHARE(2, 5) = 0.25D0

           C_SHARE(1, 6) = 0.70D0
           C_SHARE(2, 6) = 0.30D0
!
!          Transition Years
           TRANS_YEAR = 50
!          Effizienzniveau der Einkommensklassen
           DO YEAR = -1, YEARS
!               Assume constant productivity in developed countries             
                DO COUNTRY = 1, 3
                     PRODUCTIVITY_C(YEAR, 1, COUNTRY) = 1.0D0
                     PRODUCTIVITY_C(YEAR, 2, COUNTRY) = 1.0D0
                ENDDO
     
                     PRODUCTIVITY_C(YEAR, 1, 6) = 1.0D0
                     PRODUCTIVITY_C(YEAR, 2, 6) = 1.0D0

           ENDDO
!           DO YEAR = -1, 0
!               Assume increasing productivity in developing countries                
!				The product
                PRODUCTIVITY_C(-1, 1, 2) = 0.6D0
                PRODUCTIVITY_C(-1, 2, 2) = 0.6D0
                PRODUCTIVITY_C(-1, 1, 3) = 0.47D0
                PRODUCTIVITY_C(-1, 2, 3) = 0.47D0
                PRODUCTIVITY_C(-1, 1, 4) = 0.06D0
                PRODUCTIVITY_C(-1, 2, 4) = 0.06D0
                PRODUCTIVITY_C(-1, 1, 5) = 0.035D0
                PRODUCTIVITY_C(-1, 2, 5) = 0.035D0
                PRODUCTIVITY_C(0, 1, 2) = 0.6D0
                PRODUCTIVITY_C(0, 2, 2) = 0.6D0
                PRODUCTIVITY_C(0, 1, 3) = 0.47D0
                PRODUCTIVITY_C(0, 2, 3) = 0.47D0
                PRODUCTIVITY_C(0, 1, 4) = 0.06D0
                PRODUCTIVITY_C(0, 2, 4) = 0.06D0
                PRODUCTIVITY_C(0, 1, 5) = 0.035D0
                PRODUCTIVITY_C(0, 2, 5) = 0.035D0

!$$$$$$                 PRODUCTIVITY_C(-1, 1, 6) = 0.6D0
!$$$$$$                 PRODUCTIVITY_C(-1, 2, 6) = 0.6D0
!$$$$$$                 PRODUCTIVITY_C(0, 1, 6) = 0.6D0
!$$$$$$                 PRODUCTIVITY_C(0, 2, 6) = 0.6D0
                
				PRODUCTIVITY_C(-1, 1, 6) = 0.36D0
                PRODUCTIVITY_C(-1, 2, 6) = 0.36D0
                PRODUCTIVITY_C(0, 1, 6) = 0.36D0
                PRODUCTIVITY_C(0, 2, 6) = 0.36D0
				
!           ENDDO                     
           PRODUCTIVITY_C(10, 1, 2) = 1.0D0
           PRODUCTIVITY_C(10, 2, 2) = 1.0D0
           PRODUCTIVITY_C(10, 1, 3) = 1.0D0
           PRODUCTIVITY_C(10, 2, 3) = 1.0D0
           PRODUCTIVITY_C(40, 1, 6) = 1.0D0 !PRODUCTIVITY_C(10, 1, 6) = 1.0D0
           PRODUCTIVITY_C(40, 2, 6) = 1.0D0 !PRODUCTIVITY_C(10, 2, 6) = 1.0D0
           
           PRODUCTIVITY_C(15, 1, 4) = 1.0D0 !PRODUCTIVITY_C(15, 1, 4) = 1.0D0
           PRODUCTIVITY_C(15, 2, 4) = 1.0D0 !PRODUCTIVITY_C(15, 2, 4) = 1.0D0
           
           PRODUCTIVITY_C(75, 1, 5) = 1.0D0
           PRODUCTIVITY_C(75, 2, 5) = 1.0D0          
!           
           DO YEAR = 1, YEARS
                DO Y_CLASS = 1, Y_CLASSES
                     IF (YEAR <= 9) &               
                          & PRODUCTIVITY_C(YEAR, Y_CLASS, 2) = PRODUCTIVITY_C(YEAR - 1, Y_CLASS, 2) + &
                          & (PRODUCTIVITY_C(10, Y_CLASS, 2) - PRODUCTIVITY_C(0, Y_CLASS, 2)) / 11.D0  
                     IF (YEAR >= 11) PRODUCTIVITY_C(YEAR, Y_CLASS, 2) = PRODUCTIVITY_C(10, Y_CLASS, 2)
                     IF (YEAR <= 9) &               
                          & PRODUCTIVITY_C(YEAR, Y_CLASS, 3) = PRODUCTIVITY_C(YEAR - 1, Y_CLASS, 3) + &
                          & (PRODUCTIVITY_C(10, Y_CLASS, 3) - PRODUCTIVITY_C(0, Y_CLASS, 3)) / 12.D0  
                     IF (YEAR >= 11) PRODUCTIVITY_C(YEAR, Y_CLASS, 3) = PRODUCTIVITY_C(10, Y_CLASS, 3)                                  
                     IF (YEAR <= 14) &
                          & PRODUCTIVITY_C(YEAR, Y_CLASS, 4) = PRODUCTIVITY_C(YEAR - 1, Y_CLASS, 4) + &
                          & (PRODUCTIVITY_C(15, Y_CLASS, 4) - PRODUCTIVITY_C(0, Y_CLASS, 4)) / 16.D0
                     IF (YEAR >= 16) PRODUCTIVITY_C(YEAR, Y_CLASS, 4) = PRODUCTIVITY_C(15, Y_CLASS, 4)
                     IF (YEAR <= 74) &                            
                          & PRODUCTIVITY_C(YEAR, Y_CLASS, 5) = PRODUCTIVITY_C(YEAR - 1, Y_CLASS, 5) + &
                          & (PRODUCTIVITY_C(75, Y_CLASS, 5) - PRODUCTIVITY_C(0, Y_CLASS, 5)) / 76.D0
                     IF (YEAR >= 76) PRODUCTIVITY_C(YEAR, Y_CLASS, 5) = PRODUCTIVITY_C(75, Y_CLASS, 5)  

                     IF (YEAR <= 39) &
                          & PRODUCTIVITY_C(YEAR, Y_CLASS, 6) = PRODUCTIVITY_C(YEAR - 1, Y_CLASS, 6) + &
                          & (PRODUCTIVITY_C(40, Y_CLASS, 6) - PRODUCTIVITY_C(0, Y_CLASS, 6)) / 41.D0
                     IF (YEAR >= 41) PRODUCTIVITY_C(YEAR, Y_CLASS, 6) = PRODUCTIVITY_C(40, Y_CLASS, 6)  
                ENDDO    
           ENDDO  
!          Population growth rate after TRANS_YEAR
           EXOG_NPOP = 0.0D0
!          Veränderung des Immigrationsparameters MIGRATION_SCALE (muss ab TRANS_YEAR gleichbleiben!)
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = -1, YEARS
                     DO GEN = 1, GENS
                          MIGRATION_SCALE(GEN, YEAR, COUNTRY) = 1.D0
                     ENDDO
                ENDDO
           ENDDO
!          Veränderung des Anteils der Migranten in den Einkommensklassen (muss ab TRANS_YEAR gleichbleiben!)
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = -1, YEARS
                     DO Y_CLASS = 1, Y_CLASSES
                          C_SHARE_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = C_SHARE(Y_CLASS, COUNTRY)
                     ENDDO
                ENDDO
           ENDDO
!          Parameter for Technical Change / Progress
           TECH = 0.01D0           
           CALL GET_POPULATION_DEVELOPMENT
!          HOUSHOLD PARAMETERS
!          Time Preference Rate           
           DO YEAR = -1, YEARS
                DO GEN = FIRST_WORK_YEAR, GENS
                     DELTA(GEN, YEAR, 1) = 0.01D0
                     DELTA(GEN, YEAR, 2) = -0.02D0
                     DELTA(GEN, YEAR, 3) = -0.02D0
                     DELTA(GEN, YEAR, 4) = -0.07D0
                     DELTA(GEN, YEAR, 5) = -0.026D0
                     
                     DELTA(GEN, YEAR, 6) = -0.01D0
                ENDDO
           ENDDO
           DO COUNTRY = 4, 4
                IF(COUNTRY == 4)DELTA(FIRST_WORK_YEAR, 0, COUNTRY) = -0.07D0
                IF(COUNTRY == 4)DELTA(FIRST_WORK_YEAR, 25, COUNTRY) = 0.01D0
                DO YEAR = 1, YEARS
                   IF(COUNTRY == 4) THEN
                     IF(YEAR <= 24) DELTA(FIRST_WORK_YEAR, YEAR, COUNTRY) = DELTA(FIRST_WORK_YEAR, YEAR -1, COUNTRY) + &
                          & (DELTA(FIRST_WORK_YEAR, 25, COUNTRY) - DELTA(FIRST_WORK_YEAR, 0, COUNTRY)) / 26.0D0 * 1.0D0
                     IF(YEAR >= 26) DELTA(FIRST_WORK_YEAR, YEAR, COUNTRY) = DELTA(FIRST_WORK_YEAR, 25, COUNTRY) * 1.0D0 
                   ENDIF
                ENDDO
                DO GEN = 22, GENS
                     DELTA(GEN, 0, COUNTRY) = DELTA(FIRST_WORK_YEAR, 0, COUNTRY) * 1.0D0
                ENDDO
!               Um sicherzugehen, auch die anderen Alter besetzen
                DO GEN = 22, GENS
                     DO J = GEN, GENS
                          IK = GET_YEAR_BECOMING_J(0, GEN, J)
                          DELTA(J, IK, COUNTRY) = DELTA(GEN, 0, COUNTRY) * 1.0D0
                     ENDDO
                ENDDO
                DO YEAR = 0, YEARS
                     DO J = FIRST_WORK_YEAR, GENS
                          IK = GET_YEAR_BECOMING_J(YEAR, FIRST_WORK_YEAR, J)
                          DELTA(J, IK, COUNTRY) = DELTA(FIRST_WORK_YEAR, YEAR, COUNTRY) * 1.0D0
                     ENDDO
                ENDDO
           ENDDO
!          Intertemporal Elasticity of Substitution
           GAMMA = 0.250D0
!          Intratemporal Elasticity of Substitution
           RHO = 0.40D0
!          Parameter zur Gewichtung des Kindernutzens Parameter for wieghing the utility of children
           THETA = 1.0D0
!          Time Endowment
           HOURS = 1.0D0
!          Leisure Preference Parameter
           ALP = 1.50D0


!           
!          Efficiency Profile; Formeln (17), (38)
           DO COUNTRY = 1, COUNTRIES
                DO YEAR = -1, YEARS
                     DO GEN = FIRST_WORK_YEAR, GENS
                          IF(COUNTRY == 1 .AND. GEN <= 62) THEN
                               OAP(GEN, YEAR, COUNTRY) = 1.0D0
                          ELSEIF(COUNTRY == 1 .AND. GEN >= 63) THEN
                               OAP(GEN, YEAR, COUNTRY) = 0.0D0
                          ELSEIF(COUNTRY == 2 .AND. GEN <= 59) THEN
                               OAP(GEN, YEAR, COUNTRY) = 1.0D0
                          ELSEIF(COUNTRY == 2 .AND. GEN >= 60) THEN
                              OAP(GEN, YEAR, COUNTRY) = 0.0D0
!     
                          ELSEIF(COUNTRY == 3 .AND. GEN <= 59) THEN
                                OAP(GEN, YEAR, COUNTRY) = 1.0D0
                          ELSEIF(COUNTRY == 3 .AND. GEN >= 60) THEN
                                OAP(GEN, YEAR, COUNTRY) = 0.0D0
!
                          ELSEIF(COUNTRY == 4 .AND. GEN <= 59) THEN
                                OAP(GEN, YEAR, COUNTRY) = 1.0D0
                          ELSEIF(COUNTRY == 4 .AND. GEN >= 60) THEN
                                OAP(GEN, YEAR, COUNTRY) = 0.0D0
                                
                          ELSEIF(COUNTRY == 5 .AND. GEN <= 59) THEN
                                OAP(GEN, YEAR, COUNTRY) = 1.0D0
                          ELSEIF(COUNTRY == 5 .AND. GEN >= 60) THEN
                                OAP(GEN, YEAR, COUNTRY) = 0.0D0


                          ELSEIF(COUNTRY == 6 .AND. GEN <= 59) THEN
                                OAP(GEN, YEAR, COUNTRY) = 1.0D0
                          ELSEIF(COUNTRY == 6 .AND. GEN >= 60) THEN
                                OAP(GEN, YEAR, COUNTRY) = 0.0D0      
                          ENDIF
                     ENDDO
                ENDDO
           ENDDO
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = -1, YEARS
                     DO Y_CLASS = 1, Y_CLASSES                          
                          DO GEN = FIRST_WORK_YEAR, GENS
                               IF (YEAR <= 0 .OR. (YEAR > 0 .AND. GEN == FIRST_WORK_YEAR)) THEN
                                    AGE_EFFICIENCY(GEN, YEAR, Y_CLASS, COUNTRY) = PRODUCTIVITY_C(YEAR, Y_CLASS, COUNTRY) * &
                                         & EXP(4.47 + 0.033 * (GEN - LAST_EDUCATION_YEAR) - 0.00067 * &
                                         & (GEN - LAST_EDUCATION_YEAR) ** 2) / EXP(4.50233) * &
                                         & (1. + TECH) ** (GEN - FIRST_WORK_YEAR) * OAP(GEN, YEAR, COUNTRY)
                               ELSEIF (YEAR > 0 .AND. GEN > FIRST_WORK_YEAR) THEN
                                    IM21 = GET_YEAR_BECOMING_J(YEAR, GEN, FIRST_WORK_YEAR)
                                    AGE_EFFICIENCY(GEN, YEAR, Y_CLASS, COUNTRY) = PRODUCTIVITY_C(IM21, Y_CLASS, COUNTRY) * &
                                         & EXP(4.47 + 0.033 * (GEN - LAST_EDUCATION_YEAR) - 0.00067 * &
                                         & (GEN - LAST_EDUCATION_YEAR) ** 2) / EXP(4.50233) * &
                                         & (1. + TECH) ** (GEN - FIRST_WORK_YEAR) * OAP(GEN, YEAR, COUNTRY)
                               ENDIF
                          ENDDO 
                     ENDDO
                ENDDO
           ENDDO
!          PRODUCTION PARAMETERS
!          Capital Share
           ALPHA = 0.35D0
!          Depreciation rate
           DEL = 0.075D0
!		   Weight of income-class specific labor input in production function
           BETA(1) = 0.4D0
           BETA(2) = 0.25D0
!           
!          GOVERNMENT PARAMETERS
				

!Government Endowment Shares: These must sum to less than 1 in every year, the remainder is owned by households $%$%$%$%
           GOV_ENDOW_SHARE(:,:) = 0.0D0

           !Russia gets this percentage of world endowment in every year
           !GOV_ENDOW_SHARE(:,6) = 0.245D0 !GOV_ENDOW_SHARE(:,6) = 0.20D0 !We think this is closesish to data, but problems
		   GOV_ENDOW_SHARE(:,6) = 0.234D0
            
           GOV_ENDOW_SHARE(:,1) = 0.0195D0
           GOV_ENDOW_SHARE(:,2) = 0.0315D0
           GOV_ENDOW_SHARE(:,4) = 0.0335D0
           GOV_ENDOW_SHARE(:,5) = 0.0125D0


           COUNTRY_ENDOW_SHARE(:, 1) = 0.1938D0
           COUNTRY_ENDOW_SHARE(:, 2) = 0.0793D0
           COUNTRY_ENDOW_SHARE(:, 3) = 0.00D0
           COUNTRY_ENDOW_SHARE(:, 4) = 0.0838D0
           COUNTRY_ENDOW_SHARE(:, 5) = 0.0310D0
           COUNTRY_ENDOW_SHARE(:, 6) = 0.6120D0
            !UPDATE THIS WITH BETTER DATA!


           DO YEAR = -1, YEARS
! HIER ANPASSEN!             

              !Government Resource Revenues:

                !Calculate Total Government Share $%$%$%$%

       		DO COUNTRY = 1, COUNTRIES     
				TOT_GOV_ENDOW_SHARE(YEAR) = TOT_GOV_ENDOW_SHARE(YEAR) + GOV_ENDOW_SHARE(YEAR, COUNTRY)
          	ENDDO
        		Write (*,*) TOT_GOV_ENDOW_SHARE(YEAR)
                
!			   Consumption Tax
                TAX_RATE(1, 1, YEAR, 1) = -1.D0
                TAX_RATE(1, 2, YEAR, 1) = 0.D0
                TAX_RATE(1, 1, YEAR, 2) = -1.D0
                TAX_RATE(1, 2, YEAR, 2) = 0.D0
                TAX_RATE(1, 1, YEAR, 3) = -1.D0                
                TAX_RATE(1, 2, YEAR, 3) = 0.D0
                TAX_RATE(1, 1, YEAR, 4) = -1.D0                
                TAX_RATE(1, 2, YEAR, 4) = 0.D0
                TAX_RATE(1, 1, YEAR, 5) = -1.D0                
                TAX_RATE(1, 2, YEAR, 5) = 0.D0

                TAX_RATE(1, 1, YEAR, 6) = 0.090D0   !TAX_RATE(1, 1, YEAR, 6) = 0.057D0              
                TAX_RATE(1, 2, YEAR, 6) = 0.D0
                
                ENDOGENOUS_TAX_RATIO(1, YEAR, 1) = 0.5D0
                ENDOGENOUS_TAX_RATIO(1, YEAR, 2) = 0.65D0
                ENDOGENOUS_TAX_RATIO(1, YEAR, 3) = 0.58D0
                ENDOGENOUS_TAX_RATIO(1, YEAR, 4) = 0.83D0
                ENDOGENOUS_TAX_RATIO(1, YEAR, 5) = 0.83D0

                ENDOGENOUS_TAX_RATIO(1, YEAR, 6) = 0.0D0
!               Wage Tax
                TAX_RATE(2, 1, YEAR, 1) = -1.0D0
                TAX_RATE(2, 2, YEAR, 1) = 0.009D0                
                TAX_RATE(2, 1, YEAR, 2) = -1.0D0
                TAX_RATE(2, 2, YEAR, 2) = 0.009D0                
                TAX_RATE(2, 1, YEAR, 3) = -1.0D0
                TAX_RATE(2, 2, YEAR, 3) = 0.009D0                
                TAX_RATE(2, 1, YEAR, 4) = -1.0D0
                TAX_RATE(2, 2, YEAR, 4) = 0.009D0                
                TAX_RATE(2, 1, YEAR, 5) = -1.0D0
                TAX_RATE(2, 2, YEAR, 5) = 0.009D0

                TAX_RATE(2, 1, YEAR, 6) = -1.0D0
                TAX_RATE(2, 2, YEAR, 6) = 0.009D0
                ENDOGENOUS_TAX_RATIO(2, YEAR, 1) = 0.5D0
                ENDOGENOUS_TAX_RATIO(2, YEAR, 2) = 0.35D0
                ENDOGENOUS_TAX_RATIO(2, YEAR, 3) = 0.42D0
                ENDOGENOUS_TAX_RATIO(2, YEAR, 4) = 0.17D0
                ENDOGENOUS_TAX_RATIO(2, YEAR, 5) = 0.17D0

                ENDOGENOUS_TAX_RATIO(2, YEAR, 6) = 1.0D0
!               Capital Tax
                TAX_RATE(3, 1, YEAR, 1) = 0.110D0  
                TAX_RATE(3, 2, YEAR, 1) = 0.0D0
                TAX_RATE(3, 1, YEAR, 2) = 0.140D0
                TAX_RATE(3, 2, YEAR, 2) = 0.0D0
                TAX_RATE(3, 1, YEAR, 3) = 0.080D0                
                TAX_RATE(3, 2, YEAR, 3) = 0.0D0
                TAX_RATE(3, 1, YEAR, 4) = 0.03D0
                TAX_RATE(3, 2, YEAR, 4) = 0.0D0
                TAX_RATE(3, 1, YEAR, 5) = 0.03D0
                TAX_RATE(3, 2, YEAR, 5) = 0.0D0

                TAX_RATE(3, 1, YEAR, 6) = 0.018D0
                TAX_RATE(3, 2, YEAR, 6) = 0.0D0
                ENDOGENOUS_TAX_RATIO(3, YEAR, 1) = 0.D0
                ENDOGENOUS_TAX_RATIO(3, YEAR, 2) = 0.D0
                ENDOGENOUS_TAX_RATIO(3, YEAR, 3) = 0.D0
                ENDOGENOUS_TAX_RATIO(3, YEAR, 4) = 0.D0
                ENDOGENOUS_TAX_RATIO(3, YEAR, 5) = 0.D0

                ENDOGENOUS_TAX_RATIO(3, YEAR, 6) = 0.D0
!               Inheritance Tax
                INHERITANCE_TAX_RATE(YEAR, 1) = 0.D0
                INHERITANCE_TAX_RATE(YEAR, 2) = 0.D0
                INHERITANCE_TAX_RATE(YEAR, 3) = 0.D0
                INHERITANCE_TAX_RATE(YEAR, 4) = 0.D0
                INHERITANCE_TAX_RATE(YEAR, 5) = 0.D0

                INHERITANCE_TAX_RATE(YEAR, 6) = 0.D0
!               Corporate tax
				CORP_TAX(YEAR, 1) = 0.32D0                
                CORP_TAX(YEAR, 2) = 0.3D0
                CORP_TAX(YEAR, 3) = 0.3D0
                CORP_TAX(YEAR, 4) = 0.3D0
                CORP_TAX(YEAR, 5) = 0.3D0

                CORP_TAX(YEAR, 6) = 0.190D0 !CORP_TAX(YEAR, 6) = 0.175D0 

  
                
!               GOVPY is the part of the Pension Benefits contributed by the government (the other part is financed by taxes)
                MU_1(YEAR, 1) = 0.0D0
                MU_1(YEAR, 2) = 0.20D0
                MU_1(YEAR, 3) = 0.0D0
                MU_1(YEAR, 4) = 0.0D0
                MU_1(YEAR, 5) = 0.0D0

                MU_1(YEAR, 6) = 0.24D0
!               MU_2_GOV is the part of the Health Benefits treated as government consumption
                MU_2_GOV(YEAR, 1) = 0.73D0
                MU_2_GOV(YEAR, 2) = 0.73D0
                MU_2_GOV(YEAR, 3) = 0.73D0
                MU_2_GOV(YEAR, 4) = 0.73D0
                MU_2_GOV(YEAR, 5) = 0.73D0

                MU_2_GOV(YEAR, 6) = 1.0D0
!               MU_2_TAX is the part of Health Benefits financed by general taxes
                MU_2_TAX(YEAR, 1) = 0.77D0
                MU_2_TAX(YEAR, 2) = 0.20D0
                MU_2_TAX(YEAR, 3) = 0.25D0
                MU_2_TAX(YEAR, 4) = 0.0D0
                MU_2_TAX(YEAR, 5) = 0.0D0 
                
                MU_2_TAX(YEAR, 6) = 1.0D0            

                !What is MU_3? This is a parameter that relates to how disability benefits (double check,
                !seems to be the part of disability benefits paid for through general taxes)
                MU_3(YEAR, 1) = 0.0D0
                MU_3(YEAR, 2) = 0.20D0
                MU_3(YEAR, 3) = 0.0D0
                MU_3(YEAR, 4) = 0.0D0
                MU_3(YEAR, 5) = 0.0D0
                MU_3(YEAR, 6) = 0.0D0
!               MU_4 is the share of corporate tax revenues which is lump-sum transferred to households                
                MU_4(YEAR, 1) = 0.6D0!2.D0/3.D0
                MU_4(YEAR, 2) = 0.4D0
                MU_4(YEAR, 3) = 0.35D0
                MU_4(YEAR, 4) = 0.0D0
                MU_4(YEAR, 5) = 0.2D0   

                MU_4(YEAR, 6) = 0.0D0 
           ENDDO
!          Schuldenstand
           DO YEAR = -1, YEARS
                DEBT_LEVEL(YEAR, 1) = 0.70D0
                DEBT_LEVEL(YEAR, 2) = 0.76D0
                DEBT_LEVEL(YEAR, 3) = 1.46D0
                DEBT_LEVEL(YEAR, 4) = 0.21D0
                DEBT_LEVEL(YEAR, 5) = 0.72D0

                DEBT_LEVEL(YEAR, 6) = 0.1055D0!DEBT_LEVEL(YEAR, 6) = 0.104D0
           ENDDO
!          GOVS is the exogenous government spending, stands for other expenditures
           GOVS(1) = 0.034D0 * PRODUCTIVITY_C(0, 1, 1)    
           GOVS(2) = 0.075D0 * PRODUCTIVITY_C(0, 1, 2)
           GOVS(3) = 0.075D0 * PRODUCTIVITY_C(0, 1, 3)
           GOVS(4) = 0.34D0 * PRODUCTIVITY_C(0, 1, 4)
           GOVS(5) = 0.2D0 * PRODUCTIVITY_C(0, 1, 5)

           GOVS(6) = 0.142D0 * PRODUCTIVITY_C(0, 1, 6) !GOVS(6) = 0.160D0 * PRODUCTIVITY_C(0, 1, 6) 
!          Pension System
           DO YEAR = -1, YEARS
!               Renteneintrittsalter
                RETIREMENT_AGE(YEAR, 1) = 63
                RETIREMENT_AGE(YEAR, 2) = 60
                RETIREMENT_AGE(YEAR, 3) = 60
                RETIREMENT_AGE(YEAR, 4) = 60
                RETIREMENT_AGE(YEAR, 5) = 60

                RETIREMENT_AGE(YEAR, 6) = 60
!               Parameter replacement rate !!!It is not a tax rate but just a parameter for replacement rate calculating
                TAX_RATE(4, 1, YEAR, 1) = 0.65D0
                TAX_RATE(4, 2, YEAR, 1) = 0.10D0
                TAX_RATE(4, 1, YEAR, 2) = 0.6D0
                TAX_RATE(4, 2, YEAR, 2) = 0.15D0
                TAX_RATE(4, 1, YEAR, 3) = 0.4D0                
                TAX_RATE(4, 2, YEAR, 3) = 0.15D0
                TAX_RATE(4, 1, YEAR, 4) = 0.3D0
                TAX_RATE(4, 2, YEAR, 4) = 0.0D0
                TAX_RATE(4, 1, YEAR, 5) = 0.3D0
                TAX_RATE(4, 2, YEAR, 5) = 0.0D0

                !TAX_RATE(4, 1, YEAR, 6) = 0.62D0
                !TAX_RATE(4, 2, YEAR, 6) = 0.05D0
                TAX_RATE(4, 1, YEAR, 6) = 0.62D0
                TAX_RATE(4, 2, YEAR, 6) = 0.02D0
!               Faktor der Beitragsbemessungsgrenze
                CONTRIBUTION_CEILING(YEAR, 1) = 2.9D0
                CONTRIBUTION_CEILING(YEAR, 2) = 2.D0
                CONTRIBUTION_CEILING(YEAR, 3) = 1.5D0
                CONTRIBUTION_CEILING(YEAR, 4) = 3.0D0
                CONTRIBUTION_CEILING(YEAR, 5) = 3.0D0

                CONTRIBUTION_CEILING(YEAR, 6) = 1000.0D0
           ENDDO
!          Health System
!          HEALTH_EXPENDITURES_PROFILE = altersabhängiges Ausgabenprofil;
!          HEALTH_EXPENDITURE_SCALE = gibt das Niveau der Gesundheitsausgaben an
           HEALTH_EXPENDITURE_SCALE(1) = 0.032D0 * PRODUCTIVITY_C(0, 1, 1)
           HEALTH_EXPENDITURE_SCALE(2) = 0.145D0 * PRODUCTIVITY_C(0, 1, 2)
           HEALTH_EXPENDITURE_SCALE(3) = 0.027D0 * PRODUCTIVITY_C(0, 1, 3)
           HEALTH_EXPENDITURE_SCALE(4) = 0.017D0 * PRODUCTIVITY_C(0, 1, 4)
           HEALTH_EXPENDITURE_SCALE(5) = 0.010D0 * PRODUCTIVITY_C(0, 1, 5)
           
		   HEALTH_EXPENDITURE_SCALE(6) = 0.0385D0 * PRODUCTIVITY_C(0, 1, 6)!HEALTH_EXPENDITURE_SCALE(6) = 0.045D0 * PRODUCTIVITY_C(0, 1, 6)
  DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = -1, YEARS
                     DO GEN = 0, GENS
                          IF(COUNTRY <= 3) THEN
                               HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_EXPENDITURES_PROFILE(GEN, COUNTRY) * &
                                    & HEALTH_EXPENDITURE_SCALE(COUNTRY)
                          ELSEIF(COUNTRY == 6) THEN
                               HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_EXPENDITURES_PROFILE(GEN, COUNTRY) * &
                                    & HEALTH_EXPENDITURE_SCALE(COUNTRY)
                          ELSEIF(COUNTRY == 4) THEN
                               HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_EXPENDITURES_PROFILE(GEN, 3) * &
                                    & HEALTH_EXPENDITURE_SCALE(COUNTRY)
                          ELSEIF(COUNTRY == 5) THEN
                               HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_EXPENDITURES_PROFILE(GEN, 3) * &
                                    & HEALTH_EXPENDITURE_SCALE(COUNTRY)
                          ENDIF          
                     ENDDO     
                ENDDO
           ENDDO
!          Disability Insurance (USA, EU, JAPAN)
           DO YEAR = -1, YEARS
                DISABILITY_BENEFITS_IND(YEAR, 1) = 0.015D0 * PRODUCTIVITY_C(0, 1, 1)
                DISABILITY_BENEFITS_IND(YEAR, 2) = 0.031D0 * PRODUCTIVITY_C(0, 1, 2)
                DISABILITY_BENEFITS_IND(YEAR, 3) = 0.011D0 * PRODUCTIVITY_C(0, 1, 3)
                DISABILITY_BENEFITS_IND(YEAR, 4) = 0.0D0
                DISABILITY_BENEFITS_IND(YEAR, 5) = 0.0D0

                DISABILITY_BENEFITS_IND(YEAR, 6) = 0.0D0
           ENDDO
!          Education System
           EDUCATION_EXPENDITURES_SCALE(1) = 0.0033D0 * PRODUCTIVITY_C(0, 1, 1)
           EDUCATION_EXPENDITURES_SCALE(2) = 0.0047D0 * PRODUCTIVITY_C(0, 1, 2)
           EDUCATION_EXPENDITURES_SCALE(3) = 0.0038D0 * PRODUCTIVITY_C(0, 1, 3)
           EDUCATION_EXPENDITURES_SCALE(4) = 0.0042D0 * PRODUCTIVITY_C(0, 1, 4)
           EDUCATION_EXPENDITURES_SCALE(5) = 0.0023D0 * PRODUCTIVITY_C(0, 1, 5)

           EDUCATION_EXPENDITURES_SCALE(6) = 0.00405D0 * PRODUCTIVITY_C(0, 1, 6) !EDUCATION_EXPENDITURES_SCALE(6) = 0.0035D0 * PRODUCTIVITY_C(0, 1, 6)
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO GEN = 0, LAST_EDUCATION_YEAR
                     EDUCATION_EXPENDITURES_IND(GEN, COUNTRY) = EDUCATION_EXPENDITURES_PROFILE(GEN) * &
                          & EDUCATION_EXPENDITURES_SCALE(COUNTRY)
                ENDDO
           ENDDO
!          Initialize Variables
           CALL INITIALIZE_VARIABLES
!
!          Set initial assets
!          USA
           DO Y_CLASS = 1, Y_CLASSES
                ASSETS_INITIAL_YEAR(21, Y_CLASS, 1) = 0.D0
                ASSETS_INITIAL_YEAR(28, Y_CLASS, 1) = 65.9D0 * 0.01025D0 * PRODUCTIVITY_C(0, Y_CLASS, 1) * &
                     & 0.8260D0 * 0.40D0
                ASSETS_INITIAL_YEAR(35, Y_CLASS, 1) = 196.2D0 * 0.01025D0 * PRODUCTIVITY_C(0, Y_CLASS, 1) * &
                     & 0.8260D0 * 0.40D0
                ASSETS_INITIAL_YEAR(47, Y_CLASS, 1) = 362.7D0 * 0.01025D0 * PRODUCTIVITY_C(0, Y_CLASS, 1) * &
                     & 0.8260D0 * 0.40D0
                ASSETS_INITIAL_YEAR(64, Y_CLASS, 1) = 500.2D0 * 0.01025D0 * PRODUCTIVITY_C(0, Y_CLASS, 1) * &
                     & 0.8260D0 * 0.40D0
                ASSETS_INITIAL_YEAR(74, Y_CLASS, 1) = 465.5D0 * 0.01025D0 * PRODUCTIVITY_C(0, Y_CLASS, 1) * &
                     & 0.8260D0 * 0.40D0
                ASSETS_INITIAL_YEAR(82, Y_CLASS, 1) = 310.2D0 * 0.01025D0 * PRODUCTIVITY_C(0, Y_CLASS, 1) * &
                     & 0.8260D0 * 0.40D0
           ENDDO
           DO Y_CLASS = 1, Y_CLASSES
                GBA(1) = (ASSETS_INITIAL_YEAR(28, Y_CLASS, 1) - ASSETS_INITIAL_YEAR(21, Y_CLASS, 1)) / 7.0
                GBA(2) = (ASSETS_INITIAL_YEAR(35, Y_CLASS, 1) - ASSETS_INITIAL_YEAR(28, Y_CLASS, 1)) / 7.0
                GBA(3) = (ASSETS_INITIAL_YEAR(47, Y_CLASS, 1) - ASSETS_INITIAL_YEAR(35, Y_CLASS, 1)) / 12.0
                GBA(4) = (ASSETS_INITIAL_YEAR(64, Y_CLASS, 1) - ASSETS_INITIAL_YEAR(47, Y_CLASS, 1)) / 17.0
                GBA(5) = (ASSETS_INITIAL_YEAR(74, Y_CLASS, 1) - ASSETS_INITIAL_YEAR(64, Y_CLASS, 1)) / 10.0
                GBA(6) = (ASSETS_INITIAL_YEAR(82, Y_CLASS, 1) - ASSETS_INITIAL_YEAR(74, Y_CLASS, 1)) / 8.0
                GBA(7) = (0.0 - ASSETS_INITIAL_YEAR(82, Y_CLASS, 1)) / 9.0
                DO GEN = 22, 27
                    ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 1) = ASSETS_INITIAL_YEAR(GEN-1, Y_CLASS, 1) + GBA(1)
                ENDDO
                DO GEN = 29, 34
                    ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 1) = ASSETS_INITIAL_YEAR(GEN-1, Y_CLASS, 1) + GBA(2)
                ENDDO
                DO GEN = 36, 46
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 1) = ASSETS_INITIAL_YEAR(GEN-1, Y_CLASS, 1) + GBA(3)
                ENDDO
                DO GEN = 48, 63
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 1) = ASSETS_INITIAL_YEAR(GEN-1, Y_CLASS, 1) + GBA(4)
                ENDDO
                DO GEN = 65, 73
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 1) = ASSETS_INITIAL_YEAR(GEN-1, Y_CLASS, 1) + GBA(5)
                ENDDO
                DO GEN = 75, 81
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 1) = ASSETS_INITIAL_YEAR(GEN-1, Y_CLASS, 1) + GBA(6)
                ENDDO
                DO GEN = 83, GENS
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 1) = ASSETS_INITIAL_YEAR(GEN-1, Y_CLASS, 1) + GBA(7)
                ENDDO
           ENDDO
!          EU
           DO Y_CLASS = 1, Y_CLASSES
                ASSETS_INITIAL_YEAR(21, Y_CLASS, 2) = 0.0
                ASSETS_INITIAL_YEAR(22, Y_CLASS, 2) = 51.2730 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(25, Y_CLASS, 2) = 65.3330 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(30, Y_CLASS, 2) = 89.0130 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(35, Y_CLASS, 2) = 136.1440 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(40, Y_CLASS, 2) = 187.2660 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(45, Y_CLASS, 2) = 233.0530 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(50, Y_CLASS, 2) = 279.5810 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(55, Y_CLASS, 2) = 318.0230 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(60, Y_CLASS, 2) = 289.0890 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(65, Y_CLASS, 2) = 251.6210 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(70, Y_CLASS, 2) = 200.5960 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(75, Y_CLASS, 2) = 186.5430 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(80, Y_CLASS, 2) = 150.860 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS, 2) * 0.84090 * &
                     & 0.40D0
                ASSETS_INITIAL_YEAR(90, Y_CLASS, 2) = 50.290 * 0.0150 * PRODUCTIVITY_C(0, Y_CLASS,2 ) * 0.84090 * &
                     & 0.40D0
           ENDDO
           DO Y_CLASS = 1, Y_CLASSES
                GBA(1) = (ASSETS_INITIAL_YEAR(25, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(22, Y_CLASS, 2)) / 3.0
                GBA(2) = (ASSETS_INITIAL_YEAR(30, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(25, Y_CLASS, 2)) / 5.0
                GBA(3) = (ASSETS_INITIAL_YEAR(35, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(30, Y_CLASS, 2)) / 5.0
                GBA(4) = (ASSETS_INITIAL_YEAR(40, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(35, Y_CLASS, 2)) / 5.0
                GBA(5) = (ASSETS_INITIAL_YEAR(45, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(40, Y_CLASS, 2)) / 5.0
                GBA(6) = (ASSETS_INITIAL_YEAR(50, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(45, Y_CLASS, 2)) / 5.0
                GBA(7) = (ASSETS_INITIAL_YEAR(55, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(50, Y_CLASS, 2)) / 5.0
                GBA(8) = (ASSETS_INITIAL_YEAR(60, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(55, Y_CLASS, 2)) / 5.0
                GBA(9) = (ASSETS_INITIAL_YEAR(65, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(60, Y_CLASS, 2)) / 5.0
                GBA(10) = (ASSETS_INITIAL_YEAR(70, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(65, Y_CLASS, 2)) / 5.0
                GBA(11) = (ASSETS_INITIAL_YEAR(75, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(70, Y_CLASS, 2)) / 5.0
                GBA(12) = (ASSETS_INITIAL_YEAR(80, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(75, Y_CLASS, 2)) / 5.0
                GBA(13) = (ASSETS_INITIAL_YEAR(90, Y_CLASS, 2) - ASSETS_INITIAL_YEAR(80, Y_CLASS, 2)) / 10.0
                DO GEN = 23, 24
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(1)
                ENDDO
                DO GEN = 26, 29
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(2)
                ENDDO
                DO GEN = 31, 34
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(3)
                ENDDO
                DO GEN = 36, 39
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(4)
                ENDDO
                DO GEN = 41, 44
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(5)
                ENDDO
                DO GEN = 46, 49
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(6)
                ENDDO
                DO GEN = 51, 54
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(7)
                ENDDO
                DO GEN = 56, 59
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(8)
                ENDDO
                DO GEN = 61, 64
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(9)
                ENDDO
                DO GEN = 66, 69
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(10)
                ENDDO
                DO GEN = 71, 74
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(11)
                ENDDO
                DO GEN = 76, 79
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(12)
                ENDDO
                DO GEN = 81, 89
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 2) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 2) + GBA(13)
                ENDDO
           ENDDO
!          JAPAN
           DO Y_CLASS = 1, Y_CLASSES
                ASSETS_INITIAL_YEAR(21, Y_CLASS, 3) = 0.0
                ASSETS_INITIAL_YEAR(22, Y_CLASS, 3) = 29.020 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(26, Y_CLASS, 3) = 58.030 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(32, Y_CLASS, 3) = 107.170 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(37, Y_CLASS, 3) = 156.320 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(42, Y_CLASS, 3) = 218.280 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(47, Y_CLASS, 3) = 280.240 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(52, Y_CLASS, 3) = 337.760 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(57, Y_CLASS, 3) = 395.270 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(62, Y_CLASS, 3) = 448.930 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(67, Y_CLASS, 3) = 502.590 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(70, Y_CLASS, 3) = 527.010 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
                ASSETS_INITIAL_YEAR(90, Y_CLASS, 3) = 105.40 * 0.00970 * PRODUCTIVITY_C(0, Y_CLASS, 3) * 0.92010 * 0.40D0
           ENDDO
           DO Y_CLASS = 1, Y_CLASSES
                GBA(1) = (ASSETS_INITIAL_YEAR(26, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(22, Y_CLASS, 3)) / 4.0
                GBA(2) = (ASSETS_INITIAL_YEAR(32, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(26, Y_CLASS, 3)) / 6.0
                GBA(3) = (ASSETS_INITIAL_YEAR(37, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(32, Y_CLASS, 3)) / 5.0
                GBA(4) = (ASSETS_INITIAL_YEAR(42, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(37, Y_CLASS, 3)) / 5.0
                GBA(5) = (ASSETS_INITIAL_YEAR(47, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(42, Y_CLASS, 3)) / 5.0
                GBA(6) = (ASSETS_INITIAL_YEAR(52, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(47, Y_CLASS, 3)) / 5.0
                GBA(7) = (ASSETS_INITIAL_YEAR(57, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(52, Y_CLASS, 3)) / 5.0
                GBA(8) = (ASSETS_INITIAL_YEAR(62, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(57, Y_CLASS, 3)) / 5.0
                GBA(9) = (ASSETS_INITIAL_YEAR(67, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(62, Y_CLASS, 3)) / 5.0
                GBA(10) = (ASSETS_INITIAL_YEAR(70, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(67, Y_CLASS, 3)) / 3.0
                GBA(11) = (ASSETS_INITIAL_YEAR(90, Y_CLASS, 3) - ASSETS_INITIAL_YEAR(70, Y_CLASS, 3)) / 20.0
                DO GEN = 23, 25
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(1)
                ENDDO
                DO GEN = 27, 31
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(2)
                ENDDO
                DO GEN = 33, 36
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(3)
                ENDDO
                DO GEN = 38, 41
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(4)
                ENDDO
                DO GEN = 43, 46
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(5)
                ENDDO
                DO GEN = 48, 51
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(6)
                ENDDO
                DO GEN = 53, 56
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(7)
                ENDDO
                DO GEN = 58, 61
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(8)
                ENDDO
                DO GEN = 63, 66
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(9)
                ENDDO
                DO GEN = 68, 69
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(10)
                ENDDO
                DO GEN = 71, 89
                     ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) = ASSETS_INITIAL_YEAR(GEN -1, Y_CLASS, 3) + GBA(11)
                ENDDO
           ENDDO
!          CHINA and INDIA and RUSSIA
           DO Y_CLASS = 1, Y_CLASSES
                DO GEN = FIRST_WORK_YEAR, GENS
                          ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 4) = ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) * 0.10
                          ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 5) = ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) * 0.10
                          ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 6) = ASSETS_INITIAL_YEAR(GEN, Y_CLASS, 3) * 0.10
!                              & * 0.129150
                ENDDO
           ENDDO           


           RETURN
      END SUBROUTINE INITIALIZE    
      
! ********************************************************************************************************************************** 
! *    NAME: SUBROUTINE INITIALIZE_SCENARIO                                                                                        *
! *    PURPOSE: Set parameters for the reform run.                                                                                 *
! **********************************************************************************************************************************

      SUBROUTINE INITIALIZE_SCENARIO
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, Y_CLASS, YEAR
           REAL*8, EXTERNAL :: GET_LIFETIME_UTILITY
!           
!          Speichern des Nutzens zur Bestimmung der HICKS_EQUIVALENT_VARIATIONs
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = 0, YEARS
                     YY_0(YEAR, COUNTRY) = YY(YEAR, COUNTRY)
                ENDDO  
                DO Y_CLASS = 1, Y_CLASSES
                     DO GEN = FIRST_WORK_YEAR, GENS
                          UTILITY_0(GEN, 0, Y_CLASS, COUNTRY) = GET_LIFETIME_UTILITY(GEN, 0, Y_CLASS, COUNTRY)
                     ENDDO
                     DO YEAR = 1, YEARS
                          UTILITY_0(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY) = &
                               & GET_LIFETIME_UTILITY(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY)
                     ENDDO
                ENDDO
           ENDDO
!          Reformwerte eingeben
!          POPULATION PARAMETERS
!          Population growth rate after TRANS_YEAR
!$$$$$$            EXOG_NPOP = 0.00           
!          Veränderung des Immigrationsparameters MIGRATION_SCALE
!          und Migrantenanteils an den Einkommensklassen C_SHARE_MIGRANTS
!$$$$$$            DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
!$$$$$$                 DO YEAR = 0, YEARS
!$$$$$$                      DO GEN = 1, GENS
!$$$$$$                           MIGRATION_SCALE(GEN, YEAR, COUNTRY) = 0.0D0
!$$$$$$                      ENDDO
!$$$$$$                 ENDDO
!$$$$$$            ENDDO
!          Veränderung des Anteils der Migranten in den Einkommensklassen (muss ab TRANS_YEAR gleichbleiben!)
!$$$$$$            DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
!$$$$$$                 DO YEAR = 0, YEARS
!$$$$$$                      DO Y_CLASS = 1, Y_CLASSES
!$$$$$$                           C_SHARE_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = C_SHARE(Y_CLASS, COUNTRY)
!$$$$$$                      ENDDO
!$$$$$$                 ENDDO
!$$$$$$            ENDDO
!$$$$$$            CALL GET_POPULATION_DEVELOPMENT
!
           DO YEAR = 2010 - FIRST_YEAR, YEARS
!               
				!Set Corporate Taxes to Zero to simplify the Endowment Case $%$%$%$
                CORP_TAX(YEAR, 1) = 0.D0                
                CORP_TAX(YEAR, 2) = 0.D0
                CORP_TAX(YEAR, 3) = 0.D0
                CORP_TAX(YEAR, 4) = 0.D0
                CORP_TAX(YEAR, 5) = 0.D0
!                
                MU_4(YEAR, 1) = 0.D0
                MU_4(YEAR, 2) = 0.D0
                MU_4(YEAR, 3) = 0.D0
                MU_4(YEAR, 4) = 0.D0
                MU_4(YEAR, 5) = 0.D0                                
!                
!               Consumption Tax
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     TAX_RATE(1, 1, YEAR, COUNTRY) = ALEPH(1, YEAR, COUNTRY)
                     TAX_RATE(1, 2, YEAR, COUNTRY) = 0.D0
                     ENDOGENOUS_TAX_RATIO(1, YEAR, COUNTRY) = 0.D0
!               Wage Tax
!$$$$$$                      TAX_RATE(2, 1, YEAR, COUNTRY) = ALEPH(2, YEAR, COUNTRY)
!$$$$$$                      TAX_RATE(2, 2, YEAR, COUNTRY) = BETH(2, YEAR, COUNTRY)                      
                     ENDOGENOUS_TAX_RATIO(2, YEAR, COUNTRY) = 1.D0
!               Corporate tax                
!$$$$$$                 CORP_TAX(YEAR, 1) = 0.D0
!$$$$$$                 MU_4(YEAR, 1) = 0.D0
                ENDDO
           ENDDO  
!
!
           RETURN
      END SUBROUTINE INITIALIZE_SCENARIO
      
! ********************************************************************************************************************************** 
! *   NAME:  SUBROUTINE GET_POPULATION_DEVELOPMENT                                                                                 *
! *   PURPOSE:                                                                                                                     *
! **********************************************************************************************************************************

      SUBROUTINE GET_POPULATION_DEVELOPMENT
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, Y_CLASS, YEAR, PAR_AGE
           
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = 0, YEARS
                     IF(YEAR >= 1) POP(91, YEAR, Y_CLASSES+1, COUNTRY) = 0.
                     DO GEN = 0, GENS
                          IF(YEAR >= 1) POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) = 0.
                          DO Y_CLASS = 1, Y_CLASSES
                               IF(YEAR >= 1) POP(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          ENDDO
                          SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY) = 1.
                     ENDDO
                ENDDO
           ENDDO
           IF(IRUN == 0) CALL READ_POPULATION_DATA
!          Ueberlebenswahrscheinlichkeiten; Formel (2)
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = 0, YEARS
                     DO GEN = 0, GENS
                          IF(YEAR == 0) THEN
                               SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY) = 1.0
                          ELSE
                               IF(GEN == 0) THEN
                                    SURVIVAL_PROBABILITY(0, YEAR, COUNTRY) = (1.-MORTALITY(0, YEAR, COUNTRY))
                               ELSE
                                    SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY) = SURVIVAL_PROBABILITY(GEN-1, YEAR-1, COUNTRY) * &
                                         & (1.-MORTALITY(GEN, YEAR, COUNTRY))
                               ENDIF
                          ENDIF
                     ENDDO
                ENDDO
           ENDDO
!
! Now assume child mortality in INDIA
           DO YEAR = 0, YEARS
             DO GEN = 0, 20
               IF(YEAR == 0)MORTALITY(GEN, YEAR, 5) = 0.01D0
               IF(YEAR >=1 .AND. YEAR <=50)THEN
                      MORTALITY(GEN, YEAR, 5) = MORTALITY(GEN, YEAR-1, 5) - 0.007D0/50.D0
               ENDIF
               IF(YEAR > 50)MORTALITY(GEN, YEAR, 5) = 0.003D0
             ENDDO
           ENDDO      

           ! Now assume child mortality in Russia
           DO YEAR = 0, YEARS
             DO GEN = 0, 18
               IF(YEAR == 0)MORTALITY(GEN, YEAR, 6) = 0.0022D0
               IF(YEAR >=1 .AND. YEAR <=42)THEN
                      MORTALITY(GEN, YEAR, 6) = MORTALITY(GEN, YEAR-1, 6) - 0.00095D0/42.D0
               ENDIF
               IF(YEAR > 42)MORTALITY(GEN, YEAR, 6) = 0.0011D0
             ENDDO
           ENDDO

		              
!           
!$$$$$$           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
!$$$$$$                DO YEAR = 0, YEARS
!$$$$$$                     DO GEN = 0, GENS
!$$$$$$                          IF(GEN == 0) THEN
!$$$$$$                               SURVIVAL_PROBABILITY(0, YEAR, COUNTRY) = (1. - MORTALITY(0, YEAR, COUNTRY))
!$$$$$$                          ELSE
!$$$$$$                               SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY) = SURVIVAL_PROBABILITY(GEN -1, YEAR, COUNTRY) * &
!$$$$$$                                    & (1. - MORTALITY(GEN, YEAR, COUNTRY))
!$$$$$$                          ENDIF
!$$$$$$                     ENDDO
!$$$$$$                ENDDO
!$$$$$$           ENDDO               
!          Berechnen der Bevölkerung ab dem Jahr 1
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = 1, YEARS
                     DO GEN = GENS, 0, -1
                          IF(GEN > 0) THEN
                               DO Y_CLASS = 1, Y_CLASSES
                                    IF(YEAR <= TRANS_YEAR) THEN
                                         POP(GEN, YEAR, Y_CLASS, COUNTRY) = (POP(GEN-1, YEAR-1, Y_CLASS, COUNTRY) + &
                                              & MIGRATION_SCALE(GEN, YEAR, COUNTRY) * C_SHARE_MIGRANTS(YEAR, Y_CLASS, COUNTRY) * &
                                              & MIGRANTS(GEN, COUNTRY)) * (1.-MORTALITY(GEN, YEAR, COUNTRY))
                                    ELSE
                                         POP(GEN, YEAR, Y_CLASS, COUNTRY) = (POP(GEN-1, YEAR-1, Y_CLASS, COUNTRY) + &
                                              & MIGRATION_SCALE(GEN, YEAR, COUNTRY) * C_SHARE_MIGRANTS(YEAR, Y_CLASS, COUNTRY) * &
                                              & MIGRANTS(GEN, COUNTRY) * (1.+EXOG_NPOP) ** (YEAR-TRANS_YEAR)) * &
                                              & (1.-MORTALITY(GEN, YEAR, COUNTRY))
                                    ENDIF                               
                                    POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) = POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) + &
                                         & POP(GEN, YEAR, Y_CLASS, COUNTRY)
                               ENDDO
                          ELSE
                               DO Y_CLASS = 1, Y_CLASSES
                                    DO PAR_AGE = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                                         IF(YEAR <= TRANS_YEAR) THEN
                                              POP(0, YEAR, Y_CLASS, COUNTRY) = POP(0, YEAR, Y_CLASS, COUNTRY) + &
                                                   & POP(PAR_AGE, YEAR, Y_CLASS, COUNTRY) * &
                                                   & FERTILITY(PAR_AGE, YEAR, COUNTRY)* (1.-MORTALITY(0, YEAR, COUNTRY))
                                         ELSE
                                              POP(0, YEAR, Y_CLASS, COUNTRY) = (1.+EXOG_NPOP) * POP(0, YEAR-1, Y_CLASS, COUNTRY)
                                         ENDIF
                                    ENDDO
                                    POP(0, YEAR, Y_CLASSES+1, COUNTRY) = POP(0, YEAR, Y_CLASSES+1, COUNTRY) + &
                                         & POP(0, YEAR, Y_CLASS, COUNTRY)
                               ENDDO
                          ENDIF
                          POP(91, YEAR, Y_CLASSES+1, COUNTRY) = POP(91, YEAR, Y_CLASSES+1, COUNTRY) + &
                               & POP(GEN, YEAR, Y_CLASSES+1, COUNTRY)
                     ENDDO
                ENDDO
           ENDDO
!
           !CALL PRINT_POPULATION_SUMMARY            
!           
!          Calculate population growth rates           
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                NPOP(-1, COUNTRY) = EXOG_NPOP
                NPOP(0, COUNTRY) = POP(22, 0, Y_CLASSES+1, COUNTRY) / POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, COUNTRY) - 1.
                DO YEAR = 1, YEARS
                     NPOP(YEAR, COUNTRY) = POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, COUNTRY) / &
                          & POP(FIRST_WORK_YEAR, YEAR-1, Y_CLASSES+1, COUNTRY) - 1.         
                ENDDO
           ENDDO     
!           
!          Calculate fertility rates endogenously from year trans_year+1 onwards
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY                
                DO YEAR = TRANS_YEAR+1, YEARS
                     TOTAL_FERTILITY(YEAR, COUNTRY) = 0.
                     DO GEN = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                          FERTILITY(GEN, YEAR, COUNTRY) = (1.+EXOG_NPOP) * FERTILITY(GEN, YEAR-1, COUNTRY) * &
                               & POP(GEN, YEAR-1, Y_CLASSES+1, COUNTRY) / POP(GEN, YEAR, Y_CLASSES+1, COUNTRY)
                          TOTAL_FERTILITY(YEAR, COUNTRY) = TOTAL_FERTILITY(YEAR, COUNTRY) + FERTILITY(GEN, YEAR, COUNTRY)
                     ENDDO  
                ENDDO  
           ENDDO    
!           
!          Now shift population counts according to the base year in the transition
           IF (FIRST_YEAR /= 2000) THEN
                DO COUNTRY = FIRST_COUNTRY, 5
                     DO YEAR = 0, YEARS - FIRST_YEAR + 2000
                          POP(91, YEAR, Y_CLASSES+1, COUNTRY) = POP(91, YEAR + FIRST_YEAR - 2000, Y_CLASSES+1, COUNTRY)
                          DO GEN = 0, GENS
                               POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) = POP(GEN, YEAR + FIRST_YEAR - 2000, Y_CLASSES+1, COUNTRY)
                               MORTALITY(GEN, YEAR, COUNTRY) = MORTALITY(GEN, YEAR + FIRST_YEAR - 2000, COUNTRY)
                               SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY) = &
                                    & SURVIVAL_PROBABILITY(GEN, YEAR + FIRST_YEAR - 2000, COUNTRY)
                               DO Y_CLASS = 1, Y_CLASSES
                                    POP(GEN, YEAR, Y_CLASS, COUNTRY) = POP(GEN, YEAR + FIRST_YEAR - 2000, Y_CLASS, COUNTRY)
                               ENDDO
                          ENDDO                  
                     ENDDO  
                     DO YEAR = -GENS, YEARS - FIRST_YEAR + 2000
                          TOTAL_FERTILITY(YEAR, COUNTRY) = TOTAL_FERTILITY(YEAR + FIRST_YEAR - 2000, COUNTRY)
                          DO GEN = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                               FERTILITY(GEN, YEAR, COUNTRY) = FERTILITY(GEN, YEAR + FIRST_YEAR - 2000, COUNTRY)
                          ENDDO  
                     ENDDO  
!                    Fill the remaining years with the original year-300 values             
                     DO YEAR = YEARS - FIRST_YEAR + 2001, YEARS
                          POP(91, YEAR, Y_CLASSES+1, COUNTRY) = POP(91, YEARS, Y_CLASSES+1, COUNTRY) * (1. + EXOG_NPOP) ** &
                               & (YEAR - (YEARS - FIRST_YEAR + 2000))
                          DO GEN = 0, GENS
                               POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) = POP(GEN, YEARS, Y_CLASSES+1, COUNTRY) * (1. + EXOG_NPOP) ** &
                                    & (YEAR - (YEARS - FIRST_YEAR + 2000))
                               MORTALITY(GEN, YEAR, COUNTRY) = MORTALITY(GEN, YEARS, COUNTRY)
                               SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY) = SURVIVAL_PROBABILITY(GEN, YEARS, COUNTRY)
                               DO Y_CLASS = 1, Y_CLASSES
                                    POP(GEN, YEAR, Y_CLASS, COUNTRY) = POP(GEN, YEARS, Y_CLASS, COUNTRY) * &
                                         & (1. + EXOG_NPOP) ** (YEAR - (YEARS - FIRST_YEAR + 2000))
                               ENDDO
                          ENDDO 
                          TOTAL_FERTILITY(YEAR, COUNTRY) = TOTAL_FERTILITY(YEARS, COUNTRY)
                          DO GEN = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                               FERTILITY(GEN, YEAR, COUNTRY) = FERTILITY(GEN, YEARS, COUNTRY)
                          ENDDO                 
                     ENDDO 
                ENDDO                                    
           ENDIF     
           CALL PRINT_POPULATION_SUMMARY 
!             
!               Adopt population in -1 from YEARS
           IF(IRUN == 0) THEN
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     POP(91, -1, Y_CLASSES+1, COUNTRY) = POP(91, YEARS, Y_CLASSES+1, COUNTRY)
                     DO GEN = 0, GENS
                          POP(GEN, -1, Y_CLASSES+1, COUNTRY) = POP(GEN, YEARS, Y_CLASSES+1, COUNTRY)
                          MORTALITY(GEN, -1, COUNTRY) = MORTALITY(GEN, YEARS, COUNTRY)
                          SURVIVAL_PROBABILITY(GEN, -1, COUNTRY) = SURVIVAL_PROBABILITY(GEN, YEARS, COUNTRY)
                          DO Y_CLASS = 1, Y_CLASSES
                               POP(GEN, -1, Y_CLASS, COUNTRY) = POP(GEN, YEARS, Y_CLASS, COUNTRY)
                          ENDDO
                     ENDDO
                     TOTAL_FERTILITY_SS(COUNTRY) = TOTAL_FERTILITY(YEARS, COUNTRY)
                     DO GEN = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                          FERTILITY_SS(GEN, COUNTRY) = FERTILITY(GEN, YEARS, COUNTRY)
                     ENDDO
                ENDDO
!
           ENDIF     
!           
!               Calculate population growth rates           
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     NPOP(-1, COUNTRY) = EXOG_NPOP
                     NPOP(0, COUNTRY) = POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, COUNTRY) / POP(22, 0, Y_CLASSES+1, COUNTRY) - 1.
                     DO YEAR = 1, YEARS
                          NPOP(YEAR, COUNTRY) = POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, YEAR -1, Y_CLASSES+1, COUNTRY) - 1.         
                     ENDDO
                ENDDO               
!
!
!          Calculate sum of efficient overall population of age .ge. first_work_year
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = -1, YEARS
                     POP_EFFICIENT(YEAR, COUNTRY) = 0.
                     DO Y_CLASS = 1, Y_CLASSES
                          DO GEN = FIRST_WORK_YEAR, GENS
                               POP_EFFICIENT(YEAR, COUNTRY) = POP_EFFICIENT(YEAR, COUNTRY) + &
                                    & POP(GEN, YEAR, Y_CLASS, COUNTRY) * (1.+TECH) ** (FIRST_WORK_YEAR - GEN)                                    
                          ENDDO
                     ENDDO
                ENDDO
           ENDDO 
!           
!           
           RETURN
      END SUBROUTINE GET_POPULATION_DEVELOPMENT
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_INITIAL_STEADY_STATE                                                                                    *
! *   PURPOSE: Calculate the Initial Steady State                                                                                  *
! **********************************************************************************************************************************

      SUBROUTINE GET_INITIAL_STEADY_STATE
           USE GLOBAL_DATA
           IMPLICIT NONE
           CHARACTER*6 :: COUNTRY_NAME
           INTEGER :: COUNTRY, Y_CLASS
           INTEGER, PARAMETER :: MAX_ITER = 500
           REAL*8 :: H_LABOR
           REAL*8, EXTERNAL :: GET_FRACTION_PVE_CONSUMED, GET_PV_EARNINGS, GET_GOODS_MARKET
!
!          Set initial guesses for interest rate, wages and producer prices
           RG(-1) = 0.2
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                WAGE_INDEX_CLASS(-1, 1, COUNTRY) = 1.
                WAGE_INDEX_CLASS(-1, 2, COUNTRY) = 2.
           ENDDO 
!
           
!
           DO ITER = 0, MAX_ITER
!
!               Calculate AP such that the wage rate of the first income class in the US is one
                AP(-1) = 1.D0 / (BETA(1) * CAPITAL(-1, FIRST_COUNTRY)**ALPHA * LABOR(-1, 1, FIRST_COUNTRY)**(BETA(1)-1.) * &
                     & LABOR(-1, 2, FIRST_COUNTRY)**BETA(2))               
!
!               World interest rate
                RG(-1) = (1.-CORP_TAX(-1, FIRST_COUNTRY)) * (R(-1, FIRST_COUNTRY) - DEL)
!
!               Capital demand 
                IF(LAST_COUNTRY /= FIRST_COUNTRY) THEN
                     DO COUNTRY = FIRST_COUNTRY+1, LAST_COUNTRY
                          R(-1, COUNTRY) = RG(-1) / (1.-CORP_TAX(-1, COUNTRY)) + DEL
!
                          H_LABOR = 1.                          
                          DO Y_CLASS = 1, Y_CLASSES
                               H_LABOR = H_LABOR * LABOR(-1, Y_CLASS, COUNTRY)**BETA(Y_CLASS)
                               
                          ENDDO  
!
                          CAPITAL(-1, COUNTRY) = (R(-1, COUNTRY) / (ALPHA * AP(-1) * H_LABOR))**(1./(ALPHA-1.))
                     ENDDO  
                ENDIF  
                                                     
!
!               Household side
!
!               Calculate individual demands                  
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     DO Y_CLASS = 1, Y_CLASSES
                          CONSUMP(FIRST_WORK_YEAR, -1, Y_CLASS, COUNTRY) = &
                               & GET_PV_EARNINGS(FIRST_WORK_YEAR, -1, Y_CLASS, COUNTRY) * &
                               & GET_FRACTION_PVE_CONSUMED(FIRST_WORK_YEAR, -1, Y_CLASS, COUNTRY)                                   
                          CALL GET_INDIVIDUAL_DEMANDS(FIRST_WORK_YEAR, -1, Y_CLASS, COUNTRY)
                          CALL GET_SHADOW_WAGES(FIRST_WORK_YEAR, -1, Y_CLASS, COUNTRY)
                     ENDDO
                ENDDO
!                    
                CALL GET_AGGREGATE_VARIABLES(-1)                
!
!               Production side
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     CALL GET_MARGINAL_PRODUCTS(-1, COUNTRY)                
                ENDDO
!
!               Government sector           
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     CALL GET_TAXES(-1, COUNTRY)                     
                ENDDO 
!
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     CALL HEALTH_SYSTEM(-1, COUNTRY)
                     IF(COUNTRY <= 3) CALL DISABILITY_INSURANCE(-1, COUNTRY)
                     CALL PENSION_SYSTEM(-1, COUNTRY)
                ENDDO
!
                WRITE(*, 102) ITER, (CAPITAL(-1, COUNTRY), COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
102             FORMAT(I4, 1X, 5(F13.6, 1X))                            
!           
                IF (GET_GOODS_MARKET(-1) == 1) EXIT                                  
!
           ENDDO  
           IF(ITER == MAX_ITER) THEN
                     WRITE(*,110) 
110                  FORMAT('Caution! Convergence not reached!!', I3)
           ENDIF   
!
           WRITE( * , 200) ITER, (YY_WORLD(-1) - DD_WORLD(-1))
200        FORMAT(2X, 'ITER, Excess Supply = ', I6, F14.10)        
!                     
!          Print out results.
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                IF(COUNTRY == 1) COUNTRY_NAME = 'USA'
                IF(COUNTRY == 2) COUNTRY_NAME = 'EUROPA'
                IF(COUNTRY == 3) COUNTRY_NAME = 'JAPAN'
                IF(COUNTRY == 4) COUNTRY_NAME = 'CHINA'
                IF(COUNTRY == 5) COUNTRY_NAME = 'INDIA'
                IF(COUNTRY == 6) COUNTRY_NAME = 'RUSSIA'
                WRITE(5, 300) COUNTRY_NAME
300             FORMAT(24X, A6 / 15X, 'INITIAL STEADY STATE' / )
                CALL OUTPUT(-1, 1, Y_CLASSES, COUNTRY)
           ENDDO
!           
           RETURN
      END SUBROUTINE GET_INITIAL_STEADY_STATE
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_TRANSITION                                                                                              *
! *   PURPOSE: Calculate the transition                                                                                            *
! **********************************************************************************************************************************

      SUBROUTINE GET_TRANSITION
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, NFLAG, Y_CLASS, YEAR, J, Q
           INTEGER, PARAMETER :: MAX_ITER = 1500
           REAL*8 :: H_LABOR, OLD_RG, HELP_AGCO(COUNTRIES), HELP_AGC(COUNTRIES), HELP_AGW, OLD_PVE
           REAL*8, EXTERNAL :: GET_FRACTION_PVE_CONSUMED, GET_GOODS_MARKET, GET_PV_EARNINGS, GET_EFFICIENT_POPULATION
!           
		   REAL*8, DIMENSION(COUNTRIES) :: INITIAL_ASSET_SHARE   
           REAL*8 :: TEMP, ASSET_PERCENT
           IF (IRUN == 1) OPEN (UNIT = 5, FILE = OUTFILE_NAME(4), STATUS = 'REPLACE')
           IF (IRUN == 0) FIRST_SOLUTION_YEAR = 0
           IF (IRUN == 1) FIRST_SOLUTION_YEAR = 2010 - FIRST_YEAR
         
		   SIGFIG = 0.0001
           DAMP = 0.1
           DAMPK = 0.05
!$$$$$$            IF (IRUN == 1)DAMPK = 0.01
           DAMPR = 0.2
!
!
           IF (IRUN == 0) THEN
                DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                     RG(YEAR) = RG(-1)
                     DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                          WAGE_INDEX_CLASS(YEAR, 1, COUNTRY) = WAGE_INDEX_CLASS(-1, 1, 1)
                          WAGE_INDEX_CLASS(YEAR, 2, COUNTRY) = WAGE_INDEX_CLASS(-1, 2, 1)
                          CONSUMP_PRICE(YEAR, COUNTRY) = CONSUMP_PRICE(-1, COUNTRY)
                     ENDDO
                ENDDO
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                          DO Y_CLASS = 1, Y_CLASSES
                               DO J = 21, 90
                                    SHADOW_WAGE(J, YEAR, Y_CLASS, COUNTRY) = SHADOW_WAGE(J, -1, Y_CLASS, COUNTRY)
                                    PENSION_BENEFITS_IND(J, YEAR, Y_CLASS, COUNTRY) = PENSION_BENEFITS_IND(J, -1, Y_CLASS, COUNTRY)
                                    AVG_WAGE_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = AVG_WAGE_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                                    AVG_PENSION_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = AVG_PENSION_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                                    AVG_HEALTH_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = AVG_HEALTH_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                                    AVG_DISABILITY_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = &
                                         & AVG_DISABILITY_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                                    MARG_WAGE_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = MARG_WAGE_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                                    MARG_PENSION_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = &
                                         & MARG_PENSION_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                                    MARG_HEALTH_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = MARG_HEALTH_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                                    MARG_DISABILITY_TAX_RATE(J, YEAR, Y_CLASS, COUNTRY) = &
                                         & MARG_DISABILITY_TAX_RATE(J, -1, Y_CLASS, COUNTRY)
                               ENDDO  
                          ENDDO  
                     ENDDO   
                ENDDO
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                          DO Y_CLASS = 1, Y_CLASSES
                               AGG_ASSETS_FOR_BEQUESTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_FOR_BEQUESTS(-1, Y_CLASS, COUNTRY)
                          ENDDO  
                     ENDDO
                ENDDO
!
! Sum assets in the initial year
           HELP_AGW = 0.D0
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                HELP_AGCO(COUNTRY) = 0.D0
                DO Y_CLASS = 1, Y_CLASSES
                     DO GEN = FIRST_WORK_YEAR, GENS
                          HELP_AGCO(COUNTRY) = HELP_AGCO(COUNTRY) + ASSETS_INITIAL_YEAR(GEN, Y_CLASS, COUNTRY) * &
                               & GET_EFFICIENT_POPULATION(GEN, 0, Y_CLASS, COUNTRY) / (1. - MORTALITY(GEN, 0, COUNTRY))
                     ENDDO                                     
                ENDDO                        
                HELP_AGW = HELP_AGW + HELP_AGCO(COUNTRY)
           ENDDO        
!
           HELP_AGW = HELP_AGW * 4.6D0           
!           
! Now distribute world assets
!$$$$$$ 
!$$$$$$             HELP_AGC(1) = HELP_AGW * 0.3D0
!$$$$$$            HELP_AGC(2) = HELP_AGW * 0.3D0
!$$$$$$            HELP_AGC(3) = HELP_AGW * 0.1D0
!$$$$$$            HELP_AGC(4) = HELP_AGW * 0.1D0
!$$$$$$            HELP_AGC(5) = HELP_AGW * 0.1D0
!$$$$$$            HELP_AGC(6) = HELP_AGW * 0.1D0



           HELP_AGC(1) = HELP_AGW * 0.15D0
           HELP_AGC(2) = HELP_AGW * 0.2D0
           HELP_AGC(3) = HELP_AGW * 0.15D0
           HELP_AGC(4) = HELP_AGW * 0.15D0
           HELP_AGC(5) = HELP_AGW * 0.09D0
           HELP_AGC(6) = HELP_AGW * 0.16D0
!           
! Now scale individual assets accordingly           
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO Y_CLASS = 1, Y_CLASSES
                     DO GEN = FIRST_WORK_YEAR, GENS
                          ASSETS_INITIAL_YEAR(GEN, Y_CLASS, COUNTRY) = ASSETS_INITIAL_YEAR(GEN, Y_CLASS, COUNTRY) * &
                               & HELP_AGC(COUNTRY)/HELP_AGCO(COUNTRY)
                     ENDDO                                        
                ENDDO                        
           ENDDO
!
           ENDIF
!

           !!!!!!!!!!To start outside the steady state, make adjustments here:
           
           CALL GET_AGGREGATE_VARIABLES(-1)
           
           !!!!!!!!Change asset level
           
           ASSETS(:,-1,:,:) = ASSETS(:,-1,:,:)*1.0D0
           
           !!!!!!!!!Scale so that countries get the right percentage of world assets
           !!!!!!!!!The Parameters below should sum to 1
           
           CALL GET_AGGREGATE_VARIABLES(-1)
           
!$$$$$$            INITIAL_ASSET_SHARE(1) = .395D0
!$$$$$$            INITIAL_ASSET_SHARE(2) = .33D0
!$$$$$$            INITIAL_ASSET_SHARE(3) = .13D0
!$$$$$$            INITIAL_ASSET_SHARE(4) = .12D0
!$$$$$$            INITIAL_ASSET_SHARE(5) = .0175D0
!$$$$$$            INITIAL_ASSET_SHARE(6) = .0075D0
                      
           INITIAL_ASSET_SHARE(1) = .339D0
           INITIAL_ASSET_SHARE(2) = .301D0
           INITIAL_ASSET_SHARE(3) = .168D0
           INITIAL_ASSET_SHARE(4) = .117D0
           INITIAL_ASSET_SHARE(5) = .055D0
           INITIAL_ASSET_SHARE(6) = .020D0
           
            DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
             ASSET_PERCENT = AGG_ASSETS(-1,COUNTRY)/AGG_ASSETS_WORLD(-1)    
             TEMP = INITIAL_ASSET_SHARE(COUNTRY)/ASSET_PERCENT
             ASSETS_INITIAL_YEAR(:,:,COUNTRY) = ASSETS_INITIAL_YEAR(:,:,COUNTRY)*TEMP
             ASSETS(:,-1,:,COUNTRY) = ASSETS(:,-1,:,COUNTRY)*TEMP
            ENDDO
            
            CALL GET_AGGREGATE_VARIABLES(-1)
           
            DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                WRITE(*,*) AGG_ASSETS(-1,COUNTRY), AGG_ASSETS_WORLD(-1)
            ENDDO
!!!!!!!!!!To start outside the steady state, make adjustments here:
!!!!!!!!!!In this example I start russia with extra initial assets to
!!!!!!!!!!represent the oil wealth of its citizens
! BASE_ASSETS(:,:,6) = BASE_ASSETS(:,:,6)*1.2D0
           
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO Y_CLASS = 1, Y_CLASSES
                     DO GEN = FIRST_WORK_YEAR, GENS
                          IF(FIRST_SOLUTION_YEAR == 0) BASE_ASSETS(GEN, Y_CLASS, COUNTRY) = &
                               & ASSETS_INITIAL_YEAR(GEN, Y_CLASS, COUNTRY) 
                          IF(FIRST_SOLUTION_YEAR /= 0) BASE_ASSETS(GEN, Y_CLASS, COUNTRY) = &
                               & ASSETS(GEN, FIRST_SOLUTION_YEAR, Y_CLASS, COUNTRY)
                     ENDDO                                        
                ENDDO                        
           ENDDO


!
!           
           DO ITER = 0, MAX_ITER
!
!               Calculate AP such that the wage rate of the first income class in the US is one in year 0
                AP(0) = 1.D0 / (BETA(1) * CAPITAL(0, FIRST_COUNTRY)**ALPHA * LABOR(0, 1, FIRST_COUNTRY)**(BETA(1)-1.) * &
                     & LABOR(0, 2, FIRST_COUNTRY)**BETA(2))               
                DO YEAR = 1, YEARS
                     AP(YEAR) = AP(0)
                ENDDO                
!                       
!               World interest rate
                DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                     OLD_RG = RG(YEAR)
                     RG(YEAR) = (1.-CORP_TAX(YEAR, FIRST_COUNTRY)) * (R(YEAR, FIRST_COUNTRY) - DEL)
                     RG(YEAR) = OLD_RG + DAMPR * (RG(YEAR) - OLD_RG)
                ENDDO     


!
!               Capital demand 
                IF(LAST_COUNTRY /= FIRST_COUNTRY) THEN
                     DO COUNTRY = FIRST_COUNTRY+1, LAST_COUNTRY
                          DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                               R(YEAR, COUNTRY) = RG(YEAR) / (1.-CORP_TAX(YEAR, COUNTRY)) + DEL
!
                               H_LABOR = 1.                          
                               DO Y_CLASS = 1, Y_CLASSES
                                    H_LABOR = H_LABOR * LABOR(YEAR, Y_CLASS, COUNTRY)**BETA(Y_CLASS)                               
                               ENDDO  
!
                               CAPITAL(YEAR, COUNTRY) = (R(YEAR, COUNTRY) / (ALPHA * AP(YEAR) * H_LABOR))**(1./(ALPHA-1.))
                          ENDDO     
                     ENDDO  
                ENDIF         
!                      
!               Household side
!               Calculate taxes and consumer prices
!               Calculate individual demands                  
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
!                    Compute present values of earnings of FIRST_WORK_YEAR - year - old natives for every year.
                     DO Y_CLASS = 1, Y_CLASSES
                          DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                               CONSUMP(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY) = &
                                    & GET_PV_EARNINGS(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY) * &
                                    & GET_FRACTION_PVE_CONSUMED(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY)                                                   
                               CALL GET_INDIVIDUAL_DEMANDS(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY)                                                                        
                               CALL GET_SHADOW_WAGES(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY)
                          ENDDO                          
!                         Compute present values of earnings of those natives already alive in year 0
!                         (older than FIRST_WORK_YEAR years)
                          DO GEN = 22, GENS
                               CONSUMP(GEN, FIRST_SOLUTION_YEAR, Y_CLASS, COUNTRY) = &
                                    & GET_PV_EARNINGS(GEN, FIRST_SOLUTION_YEAR, Y_CLASS, COUNTRY) * &
                                    & GET_FRACTION_PVE_CONSUMED(GEN, FIRST_SOLUTION_YEAR, Y_CLASS, COUNTRY)
                               CALL GET_INDIVIDUAL_DEMANDS(GEN, FIRST_SOLUTION_YEAR, Y_CLASS, COUNTRY)
                               CALL GET_SHADOW_WAGES(GEN, FIRST_SOLUTION_YEAR, Y_CLASS, COUNTRY)
                          ENDDO                          
                     ENDDO                     
                ENDDO
!                                    
                DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                     CALL GET_AGGREGATE_VARIABLES(YEAR)
                ENDDO                
!
!               Production side
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                          CALL GET_MARGINAL_PRODUCTS(YEAR, COUNTRY)
                     ENDDO     
                ENDDO    
                
!
!
!               Government sector
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                           CALL GET_TAXES(YEAR, COUNTRY)                     
                           CALL HEALTH_SYSTEM(YEAR, COUNTRY)
                           IF(COUNTRY <= 3) CALL DISABILITY_INSURANCE(YEAR, COUNTRY)
                           CALL PENSION_SYSTEM(YEAR, COUNTRY)
                     ENDDO      
                ENDDO 

                
                NFLAG = 0
                
                DO YEAR = FIRST_SOLUTION_YEAR, YEARS
                     NFLAG = NFLAG + GET_GOODS_MARKET(YEAR)

                ENDDO

                
!
                WRITE(*,100) ITER, CAPITAL(1, FIRST_COUNTRY), CAPITAL(30, FIRST_COUNTRY), CAPITAL(50, FIRST_COUNTRY), &
                     & CAPITAL(YEARS-1, FIRST_COUNTRY), NFLAG
100             FORMAT(I4, 1X, 4(F13.6, 1X), I3)
!
                IF (NFLAG >= YEARS - (FIRST_SOLUTION_YEAR + 1) .AND. ITER >= 1) EXIT                                  
!
                IF(ITER == MAX_ITER) THEN
                     WRITE(*,110) NFLAG
110                  FORMAT('Caution! Convergence not reached!! NFLAG = ', I3)
                ENDIF                  

           ENDDO
!
           IF(NFLAG >= YEARS-(FIRST_SOLUTION_YEAR + 1))THEN
                WRITE(*, 120) NFLAG, ITER
120             FORMAT('Convergence reached in ', I3, ' markets after ', I5, ' iterations!')
           ENDIF           
!
!           
!          Print out results.
           DO YEAR = 0, YEARS
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     IF(YEAR <= 60 .OR. YEAR >= 297) CALL OUTPUT(YEAR, 1, Y_CLASSES, COUNTRY)
                ENDDO
           ENDDO
!          
!          Print initial equilibrium
           CALL PRINT_INITIAL_EQUILIBRIUM           
!
          ! write(*,*)AP(0)           
!           
           CLOSE(5)
           
            CALL GET_AGGREGATE_VARIABLES(-1)
           
            DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                WRITE(*,*) AGG_ASSETS(-1,COUNTRY), AGG_ASSETS_WORLD(-1)
                WRITE(*,*) AGG_ASSETS(1,COUNTRY), AGG_ASSETS_WORLD(1)
            ENDDO
           
           RETURN
      END SUBROUTINE GET_TRANSITION
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_PV_EARNINGS(GEN, YEAR, Y_CLASS, COUNTRY)                                                                  *
! *   PURPOSE: Compute present value of earnings of individuals age GEN in YEAR                                                    *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_PV_EARNINGS(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: IK, IKP1, J
           INTEGER, INTENT(IN) :: COUNTRY, GEN, YEAR, Y_CLASS
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8 :: S
           REAL*8, EXTERNAL :: GET_INHERITANCES, KIDS_HEALTH_BENEFITS, SUM_WAGE
!          GET_PV_EARNINGS =  present value of lifetime earnings
           GET_PV_EARNINGS = 0.
!          Set initial assets
           IF(GEN > FIRST_WORK_YEAR) GET_PV_EARNINGS = BASE_ASSETS(GEN, Y_CLASS, COUNTRY) * &
                & (1.+RG(YEAR) * CAPITAL_TAX_RATE(YEAR, COUNTRY))

!          Discount Rate
           S = 1.
!          Compute index for year corresponding to age J; Formel (22)
           DO J = GEN, GENS
                IK = GET_YEAR_BECOMING_J(YEAR, GEN, J)
                IKP1 = GET_YEAR_BECOMING_J(YEAR, GEN, J + 1)
                
                GET_PV_EARNINGS = GET_PV_EARNINGS + (HOURS * SUM_WAGE(J, IK, Y_CLASS, COUNTRY) * &
                     & (AVG_WAGE_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - AVG_PENSION_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & AVG_HEALTH_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - AVG_DISABILITY_TAX_RATE(J, IK, Y_CLASS, COUNTRY))) / S
                     
                GET_PV_EARNINGS = GET_PV_EARNINGS + ((1. - INHERITANCE_TAX_RATE(IK, COUNTRY) + RG(IK) * &
                     & CAPITAL_TAX_RATE(IK, COUNTRY)) * GET_INHERITANCES(J, IK, Y_CLASS, COUNTRY)) / S
                     
                GET_PV_EARNINGS = GET_PV_EARNINGS + (HEALTH_BENEFITS_IND(J,IK, COUNTRY) + &
                     & KIDS_HEALTH_BENEFITS(J, IK, Y_CLASS, COUNTRY)) * (1. - MU_2_GOV(IK, COUNTRY)) / S
                     
                GET_PV_EARNINGS = GET_PV_EARNINGS + TRANSFER(J, IK, Y_CLASS, COUNTRY) / S
                     
                IF(J >= FIRST_WORK_YEAR .AND. J <= 64) GET_PV_EARNINGS = GET_PV_EARNINGS + DISABILITY_BENEFITS_IND(IK, COUNTRY) / S
                  
                IF(J >= RETIREMENT_AGE(IK, COUNTRY)) GET_PV_EARNINGS = GET_PV_EARNINGS + &
                     & PENSION_BENEFITS_IND(J, IK, Y_CLASS, COUNTRY) / S
                     
                S = S * (1 + RG(IKP1) * CAPITAL_TAX_RATE(IKP1, COUNTRY))
           ENDDO
           IF(GET_PV_EARNINGS < 0) WRITE(*,*) COUNTRY, Y_CLASS, GET_PV_EARNINGS
           RETURN
      END FUNCTION GET_PV_EARNINGS
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_FRACTION_PVE_CONSUMED(GEN, YEAR, Y_CLASS, COUNTRY)                                                        *
! *   PURPOSE: Compute GET_FRACTION_PVE_CONSUMED = fraction of PVE consumed at age FIRST_WORK_YEAR.                                *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_FRACTION_PVE_CONSUMED(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: IK, IKP1, J
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8 :: H, H1, PC1, S
           REAL*8, EXTERNAL :: KIDS, SUM_WAGE
           GET_FRACTION_PVE_CONSUMED = 0.
           S = 1.
           DO J = GEN, GENS             
                IK = GET_YEAR_BECOMING_J(YEAR, GEN, J)
                IKP1 = GET_YEAR_BECOMING_J(YEAR, GEN, J + 1)
                IF(J == GEN) PC1 = CONSUMP_PRICE(IK, COUNTRY)
!                  
                H = (1. + ALP ** RHO * (SUM_WAGE(J, IK, Y_CLASS, COUNTRY) * (MARG_WAGE_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_PENSION_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - MARG_HEALTH_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_DISABILITY_TAX_RATE(J, IK, Y_CLASS, COUNTRY)) / CONSUMP_PRICE(IK, COUNTRY)) ** (1. - RHO)) ** &
                     & ((RHO - GAMMA) / (1. - RHO))
                     
                IF(J == GEN) H1 = H
                  
!               ZXI = fraction of c_ais / c_jts; Formel (30)
                ZXI(J, Y_CLASS, COUNTRY) = (S / (1 + DELTA(GEN, YEAR, COUNTRY)) ** (J - GEN) * PC1 / CONSUMP_PRICE(IK, COUNTRY) * &
                     & SURVIVAL_PROBABILITY(J, IK, COUNTRY) / SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY)) ** GAMMA * H / H1
                     
!               Marginale Konsumneigung(großes Gamma) = fraction of c_jts / W_jts and W_jts = GET_PV_EARNINGS; Formel (33)
                GET_FRACTION_PVE_CONSUMED = GET_FRACTION_PVE_CONSUMED + ZXI(J, Y_CLASS, COUNTRY) * &
                     & CONSUMP_PRICE(IK, COUNTRY) * (1. + THETA * KIDS(J, IK, Y_CLASS, COUNTRY) / H + &
                     & SUM_WAGE(J, IK, Y_CLASS, COUNTRY) * (AVG_WAGE_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & AVG_PENSION_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - AVG_HEALTH_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & AVG_DISABILITY_TAX_RATE(J, IK, Y_CLASS, COUNTRY)) / CONSUMP_PRICE(IK, COUNTRY) * &
                     & (SUM_WAGE(J, IK, Y_CLASS, COUNTRY) * (MARG_WAGE_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_PENSION_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - MARG_HEALTH_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_DISABILITY_TAX_RATE(J, IK, Y_CLASS, COUNTRY)) / &
                     & (ALP * CONSUMP_PRICE(IK, COUNTRY))) ** ( - RHO)) / S
                     
                S = S * (1. + RG(IKP1) * CAPITAL_TAX_RATE(IKP1, COUNTRY))
           ENDDO
           
           GET_FRACTION_PVE_CONSUMED = 1. / GET_FRACTION_PVE_CONSUMED
           
           RETURN
      END FUNCTION GET_FRACTION_PVE_CONSUMED
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_INDIVIDUAL_DEMANDS(GEN, YEAR, Y_CLASS, COUNTRY)                                                         *
! *   PURPOSE: Compute lifetime consumption, leisure, and asset path for person age GEN in YEAR with parents age GG at the         *
! *   person's birth for rest of his life.                                                                                         *
! **********************************************************************************************************************************

      SUBROUTINE GET_INDIVIDUAL_DEMANDS(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: IK, IKM1, IM1, J
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8 :: C1, H
           REAL*8, EXTERNAL :: GET_INHERITANCES, KIDS, KIDS_HEALTH_BENEFITS, SUM_WAGE
!          Set initial assets at beginning of year IYR.
           IF(GEN > FIRST_WORK_YEAR) ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) = BASE_ASSETS(GEN, Y_CLASS, COUNTRY)
           IF(GEN == FIRST_WORK_YEAR .AND. YEAR >= 0) ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) = 0.0
             
!          C1 is consumption in YEAR for this cohort.  The nature of the utility fct. is such that consumption in all future years
!          can be derived from C1 by multiplying it by a factor, ZXIN(GEN) (GEN represents age), a function of tax rates,
!          interest rate and wages. C1 = c_jts
           C1 = CONSUMP(GEN, YEAR, Y_CLASS, COUNTRY)
           
           DO J = GEN, GENS
                IK = GET_YEAR_BECOMING_J(YEAR, GEN, J)
!               
                CONSUMP(J, IK, Y_CLASS, COUNTRY) = ZXI(J, Y_CLASS, COUNTRY) * C1
!
                LEISURE(J, IK, Y_CLASS, COUNTRY) = CONSUMP(J, IK, Y_CLASS, COUNTRY) * ALP ** RHO * &
                     & (SUM_WAGE(J, IK, Y_CLASS, COUNTRY) * (MARG_WAGE_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_PENSION_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_HEALTH_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_DISABILITY_TAX_RATE(J, IK, Y_CLASS, COUNTRY)) / &
                     & CONSUMP_PRICE(IK, COUNTRY)) ** (-RHO)    
!
                H = (1. + ALP ** RHO * (SUM_WAGE(J, IK, Y_CLASS, COUNTRY) * (MARG_WAGE_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_PENSION_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - MARG_HEALTH_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                     & MARG_DISABILITY_TAX_RATE(J, IK, Y_CLASS, COUNTRY)) / CONSUMP_PRICE(IK, COUNTRY)) ** &
                     & (1. - RHO)) ** ((RHO - GAMMA) / (1. - RHO))
!               
                IF(J >= FIRST_FERTILITY_YEAR .AND. J <= 65) CONSUMP_KIDS(J, IK, Y_CLASS, COUNTRY) = THETA * &
                     & CONSUMP(J, IK, Y_CLASS, COUNTRY) / H
!
                IF(J > GEN) THEN
                     IKM1 = GET_YEAR_BECOMING_J(YEAR, GEN, J-1)
                     IM1 = J-1
!                    New assets = old assets + interest + net earnings - consumption; Formel (21)
                     ASSETS(J, IK, Y_CLASS, COUNTRY) = (1. + RG(IKM1) * CAPITAL_TAX_RATE(IKM1, COUNTRY)) * &
                          & (ASSETS(IM1, IKM1, Y_CLASS, COUNTRY) + GET_INHERITANCES(IM1, IKM1, Y_CLASS, COUNTRY)) - &
                          & INHERITANCE_TAX_RATE(IKM1, COUNTRY) * GET_INHERITANCES(IM1, IKM1, Y_CLASS, COUNTRY) + &
                          & (HOURS - LEISURE(IM1, IKM1, Y_CLASS, COUNTRY)) * SUM_WAGE(IM1, IKM1, Y_CLASS, COUNTRY) * &
                          & (AVG_WAGE_TAX_RATE(IM1, IKM1, Y_CLASS, COUNTRY) - AVG_PENSION_TAX_RATE(IM1, IKM1, Y_CLASS, COUNTRY) - &
                          & AVG_HEALTH_TAX_RATE(IM1, IKM1, Y_CLASS, COUNTRY) - &
                          & AVG_DISABILITY_TAX_RATE(IM1, IKM1, Y_CLASS, COUNTRY)) + &
                          & PENSION_BENEFITS_IND(IM1, IKM1, Y_CLASS, COUNTRY) - CONSUMP(IM1, IKM1, Y_CLASS, COUNTRY) * &
                          & CONSUMP_PRICE(IKM1, COUNTRY) - CONSUMP_PRICE(IKM1, COUNTRY) * KIDS(IM1, IKM1, Y_CLASS, COUNTRY) * &
                          & CONSUMP_KIDS(IM1, IKM1, Y_CLASS, COUNTRY)
                     ASSETS(J, IK, Y_CLASS, COUNTRY) = ASSETS(J, IK, Y_CLASS, COUNTRY) + &
                          & (HEALTH_BENEFITS_IND(IM1,IKM1, COUNTRY) + KIDS_HEALTH_BENEFITS(IM1, IKM1, Y_CLASS, COUNTRY))  * &
                          & (1. - MU_2_GOV(IKM1, COUNTRY)) + TRANSFER(IM1, IKM1, Y_CLASS, COUNTRY)
                     IF(J >= 22 .AND. J <= 65) ASSETS(J, IK, Y_CLASS, COUNTRY) = &
                          & ASSETS(J, IK, Y_CLASS, COUNTRY) + DISABILITY_BENEFITS_IND(IKM1, COUNTRY)
                ENDIF
                
           ENDDO
           
           RETURN
      END SUBROUTINE GET_INDIVIDUAL_DEMANDS
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_SHADOW_WAGES(GEN, YEAR, Y_CLASS, COUNTRY)                                                               *
! *   PURPOSE: Compute shadow wages for cohort age GEN in YEAR.                                                                    *
! **********************************************************************************************************************************

      SUBROUTINE GET_SHADOW_WAGES(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: IK, J, KK
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8 :: DAMP_SHW, OLDSH
           REAL*8, EXTERNAL :: GET_WAGE         
             
           DAMP_SHW = 0.6D0
!           IF(ITER >= 10) DAMP_SHW = 0.9D0
             
           DO J = GEN, GENS
                IK = GET_YEAR_BECOMING_J(YEAR, GEN, J)
                IF(J >= RETIREMENT_AGE(IK, COUNTRY)) EXIT
                IF(LEISURE(J, IK, Y_CLASS, COUNTRY) > HOURS) EXIT
                SHADOW_WAGE(J, IK, Y_CLASS, COUNTRY) = DAMP_SHW * SHADOW_WAGE(J, IK, Y_CLASS, COUNTRY)
           ENDDO
           
           IF (J >= RETIREMENT_AGE(IK, COUNTRY) .OR. LEISURE(J, IK, Y_CLASS, COUNTRY) > HOURS) THEN
                KK = J
!               Compute shadow wage for all ages
                DO J = KK, GENS
                     IK = GET_YEAR_BECOMING_J(YEAR, GEN, J)
                     OLDSH = SHADOW_WAGE(J, IK, Y_CLASS, COUNTRY)
                     SHADOW_WAGE(J, IK, Y_CLASS, COUNTRY) = CONSUMP_PRICE(IK, COUNTRY) * ALP / &
                          & (MARG_WAGE_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - MARG_PENSION_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - &
                          & MARG_HEALTH_TAX_RATE(J, IK, Y_CLASS, COUNTRY) - MARG_DISABILITY_TAX_RATE(J, IK, Y_CLASS, COUNTRY)) * &
                          & (HOURS / CONSUMP(J, IK, Y_CLASS, COUNTRY)) ** (-1. / RHO) - GET_WAGE(J, IK, Y_CLASS, COUNTRY)


!                    Damp shadow wage
                     SHADOW_WAGE(J, IK, Y_CLASS, COUNTRY) = SHADOW_WAGE(J, IK, Y_CLASS, COUNTRY) * DAMP_SHW + (1.-DAMP_SHW) * OLDSH
                ENDDO
           END IF
           
           RETURN
      END SUBROUTINE GET_SHADOW_WAGES      
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_AGGREGATE_VARIABLES(YEAR)                                                                               *
! *   PURPOSE: Damp capital and labor for iterating in transition.                                                                 *
! **********************************************************************************************************************************

      SUBROUTINE GET_AGGREGATE_VARIABLES(YEAR)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, H_COUNTRY, Y_CLASS, Q, IPL1
           INTEGER, INTENT(IN) :: YEAR
           REAL*8 :: HHB, HHK, OLD_CAP
           REAL*8, DIMENSION(Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) :: OLD_LABOR
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, SUM_INDIVIDUAL_VARIABLES
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           IPL1 = GET_YEAR_BECOMING_J(YEAR, 0, 1)
!           
           AGG_ASSETS_WORLD(YEAR) = 0.
           TRADE_BALANCE_WORLD(YEAR) = 0.
           FOREIGN_ASSETS_WORLD(YEAR) = 0.
           HHB = 0.
           HHK = 0.
!           
           OLD_CAP = CAPITAL(YEAR, FIRST_COUNTRY)           
!
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY                
                CC(YEAR, COUNTRY) = SUM_INDIVIDUAL_VARIABLES(YEAR, 1, COUNTRY)                
!                
                DO Y_CLASS = 1, Y_CLASSES                  
                     OLD_LABOR(Y_CLASS, COUNTRY) = LABOR(YEAR, Y_CLASS, COUNTRY)
                ENDDO     
                LABOR(YEAR, Y_CLASSES+1, COUNTRY) = SUM_INDIVIDUAL_VARIABLES(YEAR, 2, COUNTRY)
                DO Y_CLASS = 1, Y_CLASSES
                     LABOR(YEAR, Y_CLASS, COUNTRY) = OLD_LABOR(Y_CLASS, COUNTRY) + DAMP * (LABOR(YEAR, Y_CLASS, COUNTRY) - &                     
                          & OLD_LABOR(Y_CLASS, COUNTRY))                          
                ENDDO                          
!                
                AGG_ASSETS(YEAR, COUNTRY) = SUM_INDIVIDUAL_VARIABLES(YEAR, 3, COUNTRY)
                AGG_ASSETS_WORLD(YEAR) = AGG_ASSETS_WORLD(YEAR) + AGG_ASSETS(YEAR, COUNTRY)
           ENDDO

			!Calculate the Present Value of the endowment flow to households $%$%$%$%$
			
			!In the SS calc there is no endowment
			PVENDOWMENT(-1) = 0

!$$$$$$             IF(YEAR == -1) THEN
!$$$$$$             PVENDOWMENT(-1)=0
!$$$$$$             ENDOWMENT(-1) = ENDOWMENT(0)
!$$$$$$             PVENDOWMENT(-1) = (1.0D0 - TOT_GOV_ENDOW_SHARE(-1))*ENDOWMENT(-1)*((1+RG(-1))/RG(-1))
!$$$$$$             ENDIF
           
			IF(YEAR > -1) THEN
            PVENDOWMENT(YEAR)=0
            DO Q = 0, YEARS-YEAR
           		PVENDOWMENT(YEAR) = ((1/(1+RG(YEARS-Q)))*(PVENDOWMENT(YEAR))) + &
                & (1.0D0 - TOT_GOV_ENDOW_SHARE((YEARS-Q)))*ENDOWMENT(YEARS-Q)
        	ENDDO
            ENDIF



			!Write some additional stuff about this
        !   write(*,*) PVENDOWMENT(YEAR), RG(YEAR), YEAR
            
!
!          HHK = sum of cap.demand of SMOPECS, HHB = Gov.debt = part of K
           IF(FIRST_COUNTRY /= LAST_COUNTRY) THEN
                DO COUNTRY = FIRST_COUNTRY + 1, LAST_COUNTRY
                     HHK = HHK + CAPITAL(YEAR, COUNTRY)
                     HHB = HHB + DEBT(YEAR, COUNTRY)
                ENDDO
           ENDIF
           !Subtract PVENDOWMENT from AGG_ASSETS to get Physical Capital %$%$%$%$%$
           CAPITAL(YEAR, FIRST_COUNTRY) = AGG_ASSETS_WORLD(YEAR) - HHK - HHB - DEBT(YEAR, FIRST_COUNTRY)  - PVENDOWMENT(YEAR)

		!	write(*,*) YEAR, 'This is the agg assets, pvendowment, and usa capital', AGG_ASSETS_WORLD(YEAR), PVENDOWMENT(YEAR)
        !  	write(*,*) 'Look Here ', CAPITAL(YEAR, FIRST_COUNTRY)


           CAPITAL(YEAR, FIRST_COUNTRY) = OLD_CAP + DAMPK * (CAPITAL(YEAR, FIRST_COUNTRY) - OLD_CAP)
!
!          Sum transfers to foreign regions           
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY

                TRF(YEAR, COUNTRY) = 0.
                DO H_COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     IF (H_COUNTRY /= COUNTRY) THEN
                          DO Y_CLASS = 1, Y_CLASSES
                               DO GEN = FIRST_WORK_YEAR, GENS
                                    TRF(YEAR, COUNTRY) = TRF(YEAR, COUNTRY) + H_TRANSFER(GEN, YEAR, Y_CLASS, H_COUNTRY, COUNTRY) * &
                                         & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, H_COUNTRY) - &
                                         & H_TRANSFER(GEN, YEAR, Y_CLASS, COUNTRY, H_COUNTRY)* &
                                         & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                               ENDDO
                          ENDDO               
                     ENDIF     
                ENDDO  
           ENDDO  
!           
!          Foreign Assets and Trade Balance
           IF(FIRST_COUNTRY /= LAST_COUNTRY) THEN
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     FOREIGN_ASSETS(YEAR, COUNTRY) = AGG_ASSETS(YEAR, COUNTRY) - DEBT(YEAR, COUNTRY) - &
                          & CAPITAL(YEAR, COUNTRY)
                     IF(YEAR == -1 .OR. YEAR == YEARS) THEN
                            TRADE_BALANCE(YEAR, COUNTRY) = ((1.+TECH) * &
                               & (1.+NPOP(YEAR, FIRST_COUNTRY)) - (1.+RG(YEAR))) * FOREIGN_ASSETS(YEAR, COUNTRY) + &
                               & TRF(YEAR, COUNTRY)
                     ENDIF                             
                     IF(YEAR >= 1 .AND. YEAR > FIRST_SOLUTION_YEAR) THEN  
                            TRADE_BALANCE(YEAR-1, COUNTRY) = (1.+TECH) * &
                               & (1.+NPOP(YEAR, FIRST_COUNTRY)) * FOREIGN_ASSETS(YEAR, COUNTRY) - (1.+RG(YEAR -1)) * &
                               & FOREIGN_ASSETS(YEAR -1, COUNTRY) + TRF(YEAR-1, COUNTRY)
                     ENDIF     
                     TRADE_BALANCE_WORLD(YEAR) = TRADE_BALANCE_WORLD(YEAR) + TRADE_BALANCE(YEAR, COUNTRY)                
                     FOREIGN_ASSETS_WORLD(YEAR) = FOREIGN_ASSETS_WORLD(YEAR) + FOREIGN_ASSETS(YEAR, COUNTRY) 
                ENDDO
           ENDIF           
					FOREIGN_ASSETS_WORLD(YEAR) = FOREIGN_ASSETS_WORLD(YEAR) - PVENDOWMENT(YEAR)

                    TRADE_BALANCE_WORLD(YEAR) = TRADE_BALANCE_WORLD(YEAR) - &
                    &((1.+TECH) * (1.+NPOP(IPL1, FIRST_COUNTRY))*PVENDOWMENT(IPL1) - PVENDOWMENT(YEAR)) +&
                    &   PVENDOWMENT(YEAR)*RG(YEAR) 
               
           RETURN
      END SUBROUTINE GET_AGGREGATE_VARIABLES           
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION SUM_INDIVIDUAL_VARIABLES(YEAR, JD, COUNTRY)                                                                   *
! *   PURPOSE: General routine to compute aggregate consumption (JD = 1), labor (JD = 2), assets (JD = 3), and wages (JD = 4)      *
! **********************************************************************************************************************************

      REAL*8 FUNCTION SUM_INDIVIDUAL_VARIABLES(YEAR, JD, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: GEN, IK, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, JD, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, GET_WAGE, KIDS
           
           SUM_INDIVIDUAL_VARIABLES = 0.
!
           IF(JD == 3) AGG_ASSETS_MIGRANTS(YEAR, Y_CLASSES+1, COUNTRY) = 0.                      
!
           DO Y_CLASS = 1, Y_CLASSES           
!
                IF(JD == 2)LABOR(YEAR, Y_CLASS, COUNTRY) = 0.
                IF(JD == 3) THEN
                     AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = 0.           
                     AGG_ASSETS_FOR_BEQUESTS(YEAR, Y_CLASS, COUNTRY) = 0.
                ENDIF                  
!        
                DO GEN = FIRST_WORK_YEAR, GENS
             
                     IK = GET_YEAR_BECOMING_J(YEAR, GEN, GEN + 1)
!                                  
                     IF(JD == 1) THEN
                          SUM_INDIVIDUAL_VARIABLES = SUM_INDIVIDUAL_VARIABLES + (CONSUMP(GEN, YEAR, Y_CLASS, COUNTRY) + &
                               & KIDS(GEN, YEAR, Y_CLASS, COUNTRY) * CONSUMP_KIDS(GEN, YEAR, Y_CLASS, COUNTRY)) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)                         
                     ENDIF
!                     
                     IF(JD == 2) THEN
                          IF(LEISURE(GEN, YEAR, Y_CLASS, COUNTRY) <= HOURS .OR. GEN >= RETIREMENT_AGE(YEAR, COUNTRY)) &
                               & LABOR(YEAR, Y_CLASS, COUNTRY) = LABOR(YEAR, Y_CLASS, COUNTRY) + &
                               & AGE_EFFICIENCY(GEN, YEAR, Y_CLASS, COUNTRY) * (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY)) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ENDIF
!                
                     IF(JD == 3) THEN
                          AGG_ASSETS_FOR_BEQUESTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_FOR_BEQUESTS(YEAR, Y_CLASS, COUNTRY) + &
                               & ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) * (1. + TECH) ** (FIRST_WORK_YEAR - GEN) * &
                               & POP(GEN, YEAR, Y_CLASS, COUNTRY) * MORTALITY(GEN, YEAR, COUNTRY) / &
                               & (1. - MORTALITY(GEN, YEAR, COUNTRY))
                          SUM_INDIVIDUAL_VARIABLES = SUM_INDIVIDUAL_VARIABLES + ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY) / (1. - MORTALITY(GEN, YEAR, COUNTRY))
                          
!$$$$$$                           IF(YEAR == -1) THEN
!$$$$$$                                AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) + &
!$$$$$$                                     & ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) * MIGRANTS(GEN, COUNTRY) * &
!$$$$$$                                     & C_SHARE_MIGRANTS(YEAR, Y_CLASS, COUNTRY) * MIGRATION_SCALE(GEN, YEAR, COUNTRY) * &
!$$$$$$                                     & (1. - MORTALITY(GEN, YEAR, COUNTRY)) * &
!$$$$$$                                     & (1. + EXOG_NPOP) ** (YEARS - TRANS_YEAR+1) * (1. + TECH) ** (22 - GEN) / &
!$$$$$$                                     & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
!$$$$$$                           ELSEIF(YEAR >= 0 .AND. YEAR < TRANS_YEAR) THEN
!$$$$$$                                AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) + &
!$$$$$$                                     & ASSETS(GEN, IK, Y_CLASS, COUNTRY) * C_SHARE_MIGRANTS(IK, Y_CLASS, COUNTRY) * &
!$$$$$$                                     & MIGRATION_SCALE(GEN, IK, COUNTRY) * MIGRANTS(GEN, COUNTRY) * &
!$$$$$$                                     & (1. - MORTALITY(GEN, IK, COUNTRY)) * (1. + TECH) ** (22 - GEN) / &
!$$$$$$                                     & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
!$$$$$$                           ELSEIF(YEAR >= TRANS_YEAR .AND. YEAR <= YEARS) THEN
!$$$$$$                                AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) + &
!$$$$$$                                     & ASSETS(GEN, IK, Y_CLASS, COUNTRY) * C_SHARE_MIGRANTS(IK, Y_CLASS, COUNTRY) * &
!$$$$$$                                     & MIGRATION_SCALE(GEN, IK, COUNTRY) * MIGRANTS(GEN, COUNTRY) * &
!$$$$$$                                     & (1. - MORTALITY(GEN, IK, COUNTRY)) * (1. + EXOG_NPOP) ** (IK - TRANS_YEAR) * &
!$$$$$$                                     & (1. + TECH) ** (22 - GEN) / POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
!$$$$$$                           ENDIF
!                          
                               IF(YEAR == -1) THEN
                                    AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) + &
                                         & ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) * MIGRANTS(GEN, COUNTRY) * &
                                         & C_SHARE_MIGRANTS(YEAR, Y_CLASS, COUNTRY) * MIGRATION_SCALE(GEN, YEAR, COUNTRY) * &
                                         & (1. - MORTALITY(GEN, YEAR, COUNTRY)) * (1. + EXOG_NPOP) ** &
                                         & (YEARS - (TRANS_YEAR - FIRST_YEAR + 2000) + 1) * (1. + TECH) ** (22 - GEN) / &
                                         & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
                               ELSEIF(YEAR >= 0 .AND. YEAR < (TRANS_YEAR - FIRST_YEAR + 2000)) THEN
                                    AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) + &
                                         & ASSETS(GEN, IK, Y_CLASS, COUNTRY) * C_SHARE_MIGRANTS(IK, Y_CLASS, COUNTRY) * &
                                         & MIGRATION_SCALE(GEN, IK, COUNTRY) * MIGRANTS(GEN, COUNTRY) * &
                                         & (1. - MORTALITY(GEN, IK, COUNTRY)) * (1. + TECH) ** (22 - GEN) / &
                                         & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
                               ELSEIF(YEAR >= (TRANS_YEAR - FIRST_YEAR + 2000) .AND. YEAR <= YEARS) THEN
                                    AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) = AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY) + &
                                         & ASSETS(GEN, IK, Y_CLASS, COUNTRY) * C_SHARE_MIGRANTS(IK, Y_CLASS, COUNTRY) * &
                                         & MIGRATION_SCALE(GEN, IK, COUNTRY) * MIGRANTS(GEN, COUNTRY) * &
                                         & (1. - MORTALITY(GEN, IK, COUNTRY)) * (1. + EXOG_NPOP) ** &
                                         & (IK - (TRANS_YEAR - FIRST_YEAR + 2000)) * (1. + TECH) ** (22 - GEN) / &
                                         & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
                               ENDIF                          
                     ENDIF
!                
                     IF(JD == 4) THEN
                          IF(LEISURE(GEN, YEAR, Y_CLASS, COUNTRY) <= HOURS .OR. GEN >= RETIREMENT_AGE(YEAR, COUNTRY)) &
                               & SUM_INDIVIDUAL_VARIABLES = SUM_INDIVIDUAL_VARIABLES + &
                               & ((HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY)) * GET_WAGE(GEN, YEAR, Y_CLASS, COUNTRY)) ** 2 * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ENDIF
!
                ENDDO
!
                IF(JD == 2) SUM_INDIVIDUAL_VARIABLES = SUM_INDIVIDUAL_VARIABLES + LABOR(YEAR, Y_CLASS, COUNTRY)                
                IF(JD == 3) AGG_ASSETS_MIGRANTS(YEAR, Y_CLASSES+1, COUNTRY) = AGG_ASSETS_MIGRANTS(YEAR, Y_CLASSES+1, COUNTRY) + &
                     & AGG_ASSETS_MIGRANTS(YEAR, Y_CLASS, COUNTRY)                                        
!
           ENDDO           
!           
           RETURN
      END FUNCTION SUM_INDIVIDUAL_VARIABLES
      
! ********************************************************************************************************************************** 
! *   NAMES: SUBROUTINE GET_MARGINAL_PRODUCTS(YEAR, COUNTRY)                                                                       *
! *   PURPOSE: Compute wages for all countries, capital productivity only for BOPEC from production function.                      *
! **********************************************************************************************************************************

      SUBROUTINE GET_MARGINAL_PRODUCTS(YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: H_Y_CLASS, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, YEAR
           REAL*8 :: H_LABOR
!
           DO Y_CLASS = 1, Y_CLASSES
                H_LABOR = 1.
                DO H_Y_CLASS = 1, Y_CLASSES
                     IF(H_Y_CLASS /= Y_CLASS) H_LABOR = H_LABOR * LABOR(YEAR, H_Y_CLASS, COUNTRY)**BETA(H_Y_CLASS)
                ENDDO       
                WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY) = BETA(Y_CLASS) * AP(YEAR) * &
                          & CAPITAL(YEAR, COUNTRY)**ALPHA * LABOR(YEAR, Y_CLASS, COUNTRY)**(BETA(Y_CLASS)-1.) * H_LABOR
           ENDDO  
!
           IF(COUNTRY == FIRST_COUNTRY) THEN
                H_LABOR = 1.
                DO Y_CLASS = 1, Y_CLASSES
                     H_LABOR = H_LABOR * LABOR(YEAR, Y_CLASS, COUNTRY)**BETA(Y_CLASS) 
                ENDDO  
!
                R(YEAR, COUNTRY) = ALPHA * AP(YEAR) * CAPITAL(YEAR, COUNTRY)**(ALPHA-1.) * H_LABOR
           ENDIF     
!                 
           RETURN
      END SUBROUTINE GET_MARGINAL_PRODUCTS
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_TAXES(YEAR, COUNTRY)                                                                                    *
! *   PURPOSE: Calculate average labor earnings of working generations in period YEAR in country LA.                               *
! **********************************************************************************************************************************

      SUBROUTINE GET_TAXES(YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: GEN, I, IPASS, IPL1, KK, Y_CLASS, H_COUNTRY
           INTEGER, INTENT(IN) :: COUNTRY, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8 :: REVENUES_NEEDED, REVPRO, REVX, TDAMP, Z
           REAL*8, DIMENSION(FIRST_COUNTRY:LAST_COUNTRY) :: TAX_BASE_SQUARED
           REAL*8, DIMENSION(3, FIRST_COUNTRY:LAST_COUNTRY) :: ENDOGENOUS_TAX_BASE, TAX_BASE
           REAL*8, DIMENSION(Y_CLASSES + 1, COUNTRIES) :: BASE
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, GET_WAGE, SUM_INDIVIDUAL_VARIABLES, GET_INHERITANCES


           TOTAL_EXPENDITURES(YEAR, COUNTRY) = 0.D0  
              
           TDAMP = 0.1D0
           IF(ITER >= 20)TDAMP = 0.1D0  
           IPL1 = GET_YEAR_BECOMING_J(YEAR, 1, 2)
           EDUCATION_EXPENDITURES(YEAR, COUNTRY) = 0.
           
           DO GEN = 0, LAST_EDUCATION_YEAR
                IF(YEAR > 0) THEN
                     EDUCATION_EXPENDITURES(YEAR, COUNTRY) = EDUCATION_EXPENDITURES(YEAR, COUNTRY) + &
                          & EDUCATION_EXPENDITURES_IND(GEN, COUNTRY) * (YY_0(YEAR, COUNTRY) / YY(0, COUNTRY) * &
                          & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                          & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(0, COUNTRY) / &
                          & POP_EFFICIENT(YEAR, COUNTRY)) * (1.+TECH) ** (FIRST_WORK_YEAR - GEN) * &
                          & POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) / POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
                ELSEIF(YEAR <= 0) THEN
                     EDUCATION_EXPENDITURES(YEAR, COUNTRY) = EDUCATION_EXPENDITURES(YEAR, COUNTRY) + &
                          & EDUCATION_EXPENDITURES_IND(GEN, COUNTRY)  * (1.+TECH) ** (FIRST_WORK_YEAR - GEN) * &
                          & POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) / POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)                          
                ENDIF                
           ENDDO
!
           IF(YEAR <= 0) THEN
                GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) = GOVS(COUNTRY) * POP(91, YEAR, Y_CLASSES+1, COUNTRY) / &
                     & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) + EDUCATION_EXPENDITURES(YEAR, COUNTRY) + &
                     & MU_2_GOV(YEAR, COUNTRY) * HEALTH_BENEFITS(YEAR, COUNTRY)

                GOVERNMENT_DISCRETIONARY_SPENDING(YEAR, COUNTRY) = GOVS(COUNTRY) * POP(91, YEAR, Y_CLASSES+1, COUNTRY) / &
                     & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)      
           ELSEIF(YEAR > 0) THEN
!                               
               GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) = 0.D0
               GOVERNMENT_DISCRETIONARY_SPENDING(YEAR, COUNTRY) = 0.D0
               DO GEN = 0, GENS
                   GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) = GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) + & 
                               & GOVS(COUNTRY) * (YY_0(YEAR, COUNTRY) / YY(0, COUNTRY) * &
                               & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(0, COUNTRY) / &
                               & POP_EFFICIENT(YEAR, COUNTRY)) * &
                               & POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) / POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)                               
               ENDDO
			   GOVERNMENT_DISCRETIONARY_SPENDING(YEAR, COUNTRY) = GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) 
               
!               
               GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) = GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) + &
                   & EDUCATION_EXPENDITURES(YEAR, COUNTRY) + MU_2_GOV(YEAR, COUNTRY) * HEALTH_BENEFITS(YEAR, COUNTRY)
!
           ENDIF
!                
           DEBT(YEAR, COUNTRY) = DEBT_LEVEL(YEAR, COUNTRY) * YY(YEAR, COUNTRY)                
!
           REVENUES_NEEDED = GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) + RG(YEAR) * DEBT(YEAR, COUNTRY) + MU_1(YEAR, COUNTRY) * &
                & PENSION_BENEFITS(YEAR, COUNTRY) + (MU_2_TAX(YEAR, COUNTRY) - MU_2_GOV(YEAR, COUNTRY)) * &
                & HEALTH_BENEFITS(YEAR, COUNTRY) + MU_3(YEAR, COUNTRY) * DISABILITY_BENEFITS(YEAR, COUNTRY) - &
                & ((1. + NPOP(IPL1, FIRST_COUNTRY)) * (1.+TECH) * DEBT(IPL1, COUNTRY) - DEBT(YEAR, COUNTRY))
                
                
!
!          Inheritance tax revenues                
!
           TAX_REVENUES(4, YEAR, COUNTRY) = 0.           
           DO Y_CLASS = 1, Y_CLASSES
                DO GEN = FIRST_WORK_YEAR, GENS
                     TAX_REVENUES(4, YEAR, COUNTRY) = TAX_REVENUES(4, YEAR, COUNTRY) + INHERITANCE_TAX_RATE(YEAR, COUNTRY) * &
                          & GET_INHERITANCES(GEN, YEAR, Y_CLASS, COUNTRY)*GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                ENDDO     
           ENDDO
!           
           TAX_REVENUES(6, YEAR, COUNTRY) = TAX_REVENUES(4, YEAR, COUNTRY)
!
!          Corporate tax revenues                
!
           TAX_REVENUES(5, YEAR, COUNTRY) = YY(YEAR, COUNTRY) - AGG_ASSETS_MIGRANTS(YEAR, Y_CLASSES+1, COUNTRY)
           DO Y_CLASS = 1, Y_CLASSES
                TAX_REVENUES(5, YEAR, COUNTRY) = TAX_REVENUES(5, YEAR, COUNTRY) - WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY) * &
                     & LABOR(YEAR, Y_CLASS, COUNTRY)
           ENDDO
           TAX_REVENUES(5, YEAR, COUNTRY) = (TAX_REVENUES(5, YEAR, COUNTRY) - DEL * CAPITAL(YEAR, COUNTRY)) * &
                & CORP_TAX(YEAR, COUNTRY)
!
!          Calculate corporate tax transfers
           IF (COUNTRY == FIRST_COUNTRY) THEN
                DO H_COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                     DO Y_CLASS = 1, Y_CLASSES
                          DO GEN = FIRST_WORK_YEAR, GENS
                               TRANSFER(GEN, YEAR, Y_CLASS, H_COUNTRY) = 0.
                          ENDDO
                     ENDDO
                ENDDO            
           ENDIF  
!
           DO H_COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO Y_CLASS = 1, Y_CLASSES
                     DO GEN = FIRST_WORK_YEAR, GENS
                          H_TRANSFER(GEN, YEAR, Y_CLASS, H_COUNTRY, COUNTRY) = ASSETS(GEN, YEAR, Y_CLASS, H_COUNTRY) / &
                               & AGG_ASSETS_WORLD(YEAR) / (1. - MORTALITY(GEN, YEAR, H_COUNTRY)) * &
                               & TAX_REVENUES(5, YEAR, COUNTRY) * MU_4(YEAR, COUNTRY)
                          TRANSFER(GEN, YEAR, Y_CLASS, H_COUNTRY) = TRANSFER(GEN, YEAR, Y_CLASS, H_COUNTRY) + &
                               & H_TRANSFER(GEN, YEAR, Y_CLASS, H_COUNTRY, COUNTRY)
                     ENDDO            
                ENDDO
           ENDDO
!          
           TAX_REVENUES(5, YEAR, COUNTRY) = TAX_REVENUES(5, YEAR, COUNTRY) * (1.D0 - MU_4(YEAR, COUNTRY))
!           
           TAX_REVENUES(6, YEAR, COUNTRY) = TAX_REVENUES(6, YEAR, COUNTRY) + TAX_REVENUES(5, YEAR, COUNTRY)   

           !Add Resource Endowment Revenues
           
           TAX_REVENUES(6, YEAR, COUNTRY) = TAX_REVENUES(6, YEAR, COUNTRY) +&
           & GOV_ENDOW_SHARE(YEAR, COUNTRY)*ENDOWMENT(YEAR) 
			      
!
!          Calculate tax bases
!       
           TAX_BASE(1, COUNTRY) = CC(YEAR, COUNTRY)
!                
		   TAX_BASE(2, COUNTRY) = 0.
           DO Y_CLASS = 1, Y_CLASSES
                TAX_BASE(2, COUNTRY) = TAX_BASE(2, COUNTRY) + WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY) * &
                     & LABOR(YEAR, Y_CLASS, COUNTRY)
           ENDDO  
!                     
           TAX_BASE_SQUARED(COUNTRY) = SUM_INDIVIDUAL_VARIABLES(YEAR, 4, COUNTRY)
!                
           TAX_BASE(3, COUNTRY) = AGG_ASSETS(YEAR, COUNTRY) * RG(YEAR)
!           
           DO Y_CLASS = 1, Y_CLASSES + 1
                AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                BASE(Y_CLASS, COUNTRY) = 0.
           ENDDO
!
           DO I = 1, 3 
                ENDOGENOUS_TAX_BASE(I, COUNTRY) = 0.
           ENDDO     
!           
           DO IPASS = 1, 2
!               KK = Steuerart
                DO KK = 1, 3
!                    ALEPH is the proportional term, BETH the progressive
!                    REVPRO is the proportion of needed revenues covered by the endogenous tax
                     ALEPH(KK, YEAR, COUNTRY) = TAX_RATE(KK, 1, YEAR, COUNTRY)
                     BETH(KK, YEAR, COUNTRY) = TAX_RATE(KK, 2, YEAR, COUNTRY)  
                     REVPRO = ENDOGENOUS_TAX_RATIO(KK, YEAR, COUNTRY)
                     IF(IPASS == 1 .AND. (ALEPH(KK, YEAR, COUNTRY) == -1 .OR. BETH(KK, YEAR, COUNTRY) == -1)) THEN
                          ENDOGENOUS_TAX_BASE(KK, COUNTRY) = TAX_BASE(KK, COUNTRY)                            
                     ENDIF     
                     IF((IPASS == 1 .AND. (ALEPH(KK, YEAR, COUNTRY) /= -1 .AND. BETH(KK, YEAR, COUNTRY) /= -1.)) .OR. &
                          & (IPASS == 2 .AND. (ALEPH(KK, YEAR, COUNTRY) == -1 .OR. BETH(KK, YEAR, COUNTRY) == -1.))) THEN
                          IF(IPASS == 2 .AND. ALEPH(KK, YEAR, COUNTRY) == -1.) then
                            ALEPH(KK, YEAR, COUNTRY) = (REVX * REVPRO - &
                               & BETH(KK, YEAR, COUNTRY) * TAX_BASE_SQUARED(COUNTRY) / 2.) / &
                               & ENDOGENOUS_TAX_BASE(KK, COUNTRY)
        !                       write(*,*)kk,ALEPH(KK, YEAR, COUNTRY)
                          endif     
                          IF(IPASS == 2 .AND. BETH(KK, YEAR, COUNTRY) == -1.) BETH(KK, YEAR, COUNTRY) = 2. * (REVX * REVPRO - &
                               & ALEPH(KK, YEAR, COUNTRY) * ENDOGENOUS_TAX_BASE(KK, COUNTRY)) / TAX_BASE_SQUARED(COUNTRY)
!                               
                          IF(KK == 1) THEN
                               CONSUMP_PRICE(YEAR, COUNTRY) = CONSUMP_PRICE(YEAR, COUNTRY) * (1. - TDAMP) + &
                                    & (1. + ALEPH(KK, YEAR, COUNTRY)) * TDAMP
                               TAX_REVENUES(KK, YEAR, COUNTRY) = (CONSUMP_PRICE(YEAR, COUNTRY) -1.) * CC(YEAR, COUNTRY)
                               TAX_REVENUES(6, YEAR, COUNTRY) = TAX_REVENUES(6, YEAR, COUNTRY) + TAX_REVENUES(KK, YEAR, COUNTRY)
!
                          ELSEIF(KK == 2) THEN
                               TAX_REVENUES(KK, YEAR, COUNTRY) = 0.
                               DO Y_CLASS = 1, Y_CLASSES
                                    DO GEN = FIRST_WORK_YEAR, GENS
                                         Z = (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY)) * GET_WAGE(GEN, YEAR, Y_CLASS, COUNTRY)
                                         MARG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = TDAMP * &
                                              & (1. - ALEPH(KK, YEAR, COUNTRY) - BETH(KK, YEAR, COUNTRY) * Z) + (1. - TDAMP) * &
                                              & MARG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY)
                                         AVG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = TDAMP * &
                                              & (1. - ALEPH(KK, YEAR, COUNTRY) - BETH(KK, YEAR, COUNTRY) * Z / 2.) + (1. - TDAMP) * &
                                              & AVG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY)
                                         AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                                              & AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + (ALEPH(KK, YEAR, COUNTRY) + &
                                              & BETH(KK, YEAR, COUNTRY) * Z) * Z * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)                                                   
                                         AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = &
                                              & AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                                              & AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY)                                                   
                                         AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                                              & AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + (ALEPH(KK, YEAR, COUNTRY) + &
                                              & BETH(KK, YEAR, COUNTRY) * Z / 2.) * Z * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                                         AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = &
                                              & AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                                              & AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY)
                                         BASE(Y_CLASS, COUNTRY) = BASE(Y_CLASS, COUNTRY) + Z * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                                         BASE(Y_CLASSES+1, COUNTRY) = BASE(Y_CLASSES+1, COUNTRY) + BASE(Y_CLASS, COUNTRY)
                                         TAX_REVENUES(KK, YEAR, COUNTRY) = TAX_REVENUES(KK, YEAR, COUNTRY) + &
                                              & (1. - AVG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY)) * Z * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                                    ENDDO
                               ENDDO
                               
                               TAX_REVENUES(6, YEAR, COUNTRY) = TAX_REVENUES(6, YEAR, COUNTRY) + TAX_REVENUES(KK, YEAR, COUNTRY)
!                               
                          ELSEIF(KK == 3) THEN
                               CAPITAL_TAX_RATE(YEAR, COUNTRY) = CAPITAL_TAX_RATE(YEAR, COUNTRY) * (1. - TDAMP) + &
                                    & (1. - ALEPH(KK, YEAR, COUNTRY)) * TDAMP
                               TAX_REVENUES(KK, YEAR, COUNTRY) = (1. - CAPITAL_TAX_RATE(YEAR, COUNTRY)) * &
                                    & AGG_ASSETS(YEAR, COUNTRY) * RG(YEAR)
                               TAX_REVENUES(6, YEAR, COUNTRY) = TAX_REVENUES(6, YEAR, COUNTRY) + TAX_REVENUES(KK, YEAR, COUNTRY)
                          ENDIF
                     ENDIF
                ENDDO
                
                IF(IPASS == 1) REVX = REVENUES_NEEDED - TAX_REVENUES(6, YEAR, COUNTRY)
                  
           ENDDO
           
           DO Y_CLASS = 1, Y_CLASSES+1
                AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & BASE(Y_CLASS, COUNTRY)
                AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & BASE(Y_CLASS, COUNTRY)
           ENDDO           
!
           RETURN
      END SUBROUTINE GET_TAXES      
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE HEALTH_SYSTEM(YEAR, COUNTRY)                                                                                *
! *   PURPOSE: Calculate components of the country's health system.                                                                *
! **********************************************************************************************************************************

      SUBROUTINE HEALTH_SYSTEM(YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: GEN, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, YEAR
           REAL*8 :: MAX_TAXABLE_EARNINGS, TDAM, TDAMP
           REAL*8, DIMENSION(Y_CLASSES + 1) :: EARNINGS_BASE_PAYROLL
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, Y_CLASSES) :: EARNINGS_BASE_PAYROLL_IND
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, KIDS_HEALTH_BENEFITS, SUM_WAGE
           
           TDAMP = 1.D0
           TDAM = 1.D0
!           
           HEALTH_BENEFITS(YEAR, COUNTRY) = 0.
           EARNINGS_BASE_PAYROLL(Y_CLASSES+1) = 0.
!           
           CALL GET_AVG_LABOR_EARNINGS(YEAR, COUNTRY)
           
           MAX_TAXABLE_EARNINGS = CONTRIBUTION_CEILING(YEAR, COUNTRY) * AVG_LABOR_EARNINGS(YEAR, COUNTRY)
           
           DO GEN = 0, GENS
             
                IF(COUNTRY <= 3 .OR. COUNTRY == 6)THEN
                     IF(YEAR >= 0 .AND. YEAR <= 25) THEN
                          HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_BENEFITS_IND(GEN, 0, COUNTRY) * &
                               & (YY_0(YEAR, COUNTRY) / YY(0, COUNTRY) * POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(0, COUNTRY) / &
                               & POP_EFFICIENT(YEAR, COUNTRY)) * 1.02D0 ** YEAR
                     ELSEIF(YEAR >= 26 .AND. YEAR <= 35) THEN          
                          HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_BENEFITS_IND(GEN, 25, COUNTRY) * &
                               & (YY_0(YEAR, COUNTRY) / YY_0(25, COUNTRY) * &
                               & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, 25, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(25, COUNTRY) / &
                               & POP_EFFICIENT(YEAR, COUNTRY)) * 1.01D0 ** (YEAR - 25)
                     ELSEIF(YEAR >= 36) THEN          
                          HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_BENEFITS_IND(GEN, 35, COUNTRY) * &
                               & (YY_0(YEAR, COUNTRY) / YY_0(35, COUNTRY) * &
                               & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, 35, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(35, COUNTRY) / &
                               & POP_EFFICIENT(YEAR, COUNTRY))          
                     ENDIF
                ELSEIF(COUNTRY >= 4 .AND. COUNTRY<=5) THEN
                     IF(YEAR >= 0 .AND. YEAR <= 40) THEN
                          HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_BENEFITS_IND(GEN, 0, COUNTRY) * &
                               & (YY_0(YEAR, COUNTRY) / YY(0, COUNTRY) * &
                               & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(0, COUNTRY) / &
                               & POP_EFFICIENT(YEAR, COUNTRY)) * 1.04D0 ** YEAR
                     ELSEIF(YEAR >= 41)THEN          
                          HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_BENEFITS_IND(GEN, 41, COUNTRY) * &
                               & (YY_0(YEAR, COUNTRY) / YY_0(41, COUNTRY) * &
                               & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, 41, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(41, COUNTRY) / &
                               & POP_EFFICIENT(YEAR, COUNTRY))          
                     ENDIF          
                ENDIF                                               
!
                IF(YEAR == -1) HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) = HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY)           
!                                  
           ENDDO
           
           DO Y_CLASS = 1, Y_CLASSES
                DO GEN = FIRST_WORK_YEAR, GENS
                     HEALTH_BENEFITS(YEAR, COUNTRY) = HEALTH_BENEFITS(YEAR, COUNTRY) + (HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) + &
                          & KIDS_HEALTH_BENEFITS(GEN, YEAR, Y_CLASS, COUNTRY)) * &
                          & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                ENDDO          
           ENDDO            
!            
           DO Y_CLASS = 1, Y_CLASSES
                EARNINGS_BASE_PAYROLL(Y_CLASS) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY)-1
                     EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) = SUM_WAGE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY))
                     IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) < 0.) EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) = 0.
                     IF(COUNTRY == 1) THEN
                          EARNINGS_BASE_PAYROLL(Y_CLASS) = EARNINGS_BASE_PAYROLL(Y_CLASS) + &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ELSE
                          IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) <= MAX_TAXABLE_EARNINGS) THEN
                               EARNINGS_BASE_PAYROLL(Y_CLASS) = EARNINGS_BASE_PAYROLL(Y_CLASS) + &
                                    & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                                    & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          ELSE
                               EARNINGS_BASE_PAYROLL(Y_CLASS) = EARNINGS_BASE_PAYROLL(Y_CLASS) + MAX_TAXABLE_EARNINGS * &
                                    & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          ENDIF
                     ENDIF
                ENDDO
                EARNINGS_BASE_PAYROLL(Y_CLASSES+1) = EARNINGS_BASE_PAYROLL(Y_CLASSES+1) + EARNINGS_BASE_PAYROLL(Y_CLASS)                
           ENDDO
           
           
           IF(EARNINGS_BASE_PAYROLL(Y_CLASSES+1) /= 0) THEN
                AGG_HEALTH_TAX_RATE(YEAR, COUNTRY) = (1. - MU_2_TAX(YEAR, COUNTRY)) * (1. - MU_2_GOV(YEAR, COUNTRY) + &
                     & MU_2_GOV(YEAR, COUNTRY)) * HEALTH_BENEFITS(YEAR, COUNTRY) / EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           ELSE
                WRITE( * , 100) COUNTRY, YEAR
100             FORMAT(' PY = 0, COUNTRY = ', I1, ', YEAR = ', I3)
           ENDIF
           
           AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = 0.
           AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = 0.
           DO Y_CLASS = 1, Y_CLASSES
                AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) -1
                     IF(COUNTRY == 1) THEN
                          MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                               & MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAMP) + &
                               & TDAMP * AGG_HEALTH_TAX_RATE(YEAR, COUNTRY)
                          AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                               & AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAMP) + &
                               & TDAMP * AGG_HEALTH_TAX_RATE(YEAR, COUNTRY)
                          AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                               & AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ELSE
                          IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) <= MAX_TAXABLE_EARNINGS) THEN
                               IF(Y_CLASS < Y_CLASSES) THEN
                                    MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                                         & MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAMP) + &
                                         & TDAMP * AGG_HEALTH_TAX_RATE(YEAR, COUNTRY)
                               ELSE
                                    MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               ENDIF
                               AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                                    & AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAMP) + &
                                    & TDAMP * AGG_HEALTH_TAX_RATE(YEAR, COUNTRY)
                               AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                                    & AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                                    & MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                                    & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                                    & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                               AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                                    & AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                                    & AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                                    & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                                    & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          ELSE
                               MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                                    & AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAM) + &
                                    & TDAM * AGG_HEALTH_TAX_RATE(YEAR, COUNTRY) * MAX_TAXABLE_EARNINGS / &
                                    & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS)
                               AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                                    & AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                                    & MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * MAX_TAXABLE_EARNINGS * &
                                    & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                               AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                                    & AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                                    & AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * MAX_TAXABLE_EARNINGS * &
                                    & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          ENDIF
                     ENDIF
                          
                     IF(AGG_HEALTH_TAX_RATE(YEAR, COUNTRY) == 0.) THEN
                          MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                     ENDIF
                          
                ENDDO
                
                AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                     & AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY)
                AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                     & AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY)
                AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASS)
                AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASS)
                     
           ENDDO
           
           AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) / &
                & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) / &
                & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           
           RETURN
      END SUBROUTINE HEALTH_SYSTEM
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE DISABILITY_INSURANCE(YEAR, COUNTRY)                                                                         *
! *   PURPOSE: Calculate components of the country's disability insurance system.                                                  *
! **********************************************************************************************************************************

      SUBROUTINE DISABILITY_INSURANCE(YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: GEN, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, YEAR
           REAL*8 :: MAX_TAXABLE_EARNINGS, TDAM, TDAMP
           REAL*8, DIMENSION(Y_CLASSES + 1) :: EARNINGS_BASE_PAYROLL
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, Y_CLASSES) :: EARNINGS_BASE_PAYROLL_IND
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, SUM_WAGE
           
           TDAMP = 1.D0
           TDAM = 1.D0
           
           DISABILITY_BENEFITS(YEAR, COUNTRY) = 0.
           EARNINGS_BASE_PAYROLL(Y_CLASSES+1) = 0.
           
!          Nur einschalten, wenn HEALTH_SYSTEM fehlt
!$$$$$$            CALL GET_AVG_LABOR_EARNINGS(YEAR, COUNTRY)

           MAX_TAXABLE_EARNINGS = CONTRIBUTION_CEILING(YEAR, COUNTRY) * AVG_LABOR_EARNINGS(YEAR, COUNTRY)
!
           IF(YEAR > 0) DISABILITY_BENEFITS_IND(YEAR, COUNTRY) = DISABILITY_BENEFITS_IND(0, COUNTRY) * (YY_0(YEAR, COUNTRY) / &
                & YY(0, COUNTRY) * POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(0, COUNTRY) / &
                & POP_EFFICIENT(YEAR, COUNTRY))
                            
           DO GEN = FIRST_WORK_YEAR, 64
                DISABILITY_BENEFITS(YEAR, COUNTRY) = DISABILITY_BENEFITS(YEAR, COUNTRY) + DISABILITY_BENEFITS_IND(YEAR, COUNTRY) * &
                     & (1.+TECH) ** (FIRST_WORK_YEAR - GEN) * POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) / &
                     & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
           ENDDO           
!           
           DO Y_CLASS = 1, Y_CLASSES
                EARNINGS_BASE_PAYROLL(Y_CLASS) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) -1
                     EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) = SUM_WAGE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY))
                     IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) < 0.) EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) = 0.
                     IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) <= MAX_TAXABLE_EARNINGS) THEN
                          EARNINGS_BASE_PAYROLL(Y_CLASS) = EARNINGS_BASE_PAYROLL(Y_CLASS) + &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ELSE
                          EARNINGS_BASE_PAYROLL(Y_CLASS) = EARNINGS_BASE_PAYROLL(Y_CLASS) + &
                               & MAX_TAXABLE_EARNINGS * GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ENDIF
                ENDDO
                EARNINGS_BASE_PAYROLL(Y_CLASSES+1) = EARNINGS_BASE_PAYROLL(Y_CLASSES+1) + EARNINGS_BASE_PAYROLL(Y_CLASS)
           ENDDO
           
           IF(EARNINGS_BASE_PAYROLL(Y_CLASSES+1) /= 0) THEN
                AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY) = (1. - MU_3(YEAR, COUNTRY)) * DISABILITY_BENEFITS(YEAR, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           ELSE
                WRITE(*, 100) COUNTRY, YEAR
100             FORMAT('DISABILITY, PY = 0, COUNTRY = ', I1, ', YEAR = ', I3)
           ENDIF
           
           AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = 0.
           AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = 0.
           DO Y_CLASS = 1, Y_CLASSES
                AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) -1
                     IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) <= MAX_TAXABLE_EARNINGS) THEN
                          IF(Y_CLASS < Y_CLASSES) THEN
                               MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                                    & MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAMP) + &
                                    & TDAMP * AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY)
                          ELSE
                               MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          ENDIF
                          AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                               & AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1. - TDAMP) + &
                               & TDAMP * AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY)
                          AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                               & AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                               & AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ELSE
                          MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                               & AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                               & (1.-TDAM) + TDAM * AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY) * MAX_TAXABLE_EARNINGS / &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS)
                          AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                               & AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * MAX_TAXABLE_EARNINGS * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = &
                               & AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * MAX_TAXABLE_EARNINGS * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ENDIF
                     
                     IF(AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY) == 0.) THEN
                          MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                     ENDIF
                          
                ENDDO
                
                AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = &
                     & AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                     & AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY)
                AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = &
                     & AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                     & AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY)
                AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASS)
                AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASS)
                     
           ENDDO
           
           AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) / &
                & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) / &
                & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           
           RETURN
      END SUBROUTINE DISABILITY_INSURANCE
      
! **********************************************************************************************************************************
! *   NAME: SUBROUTINE PENSION_SYSTEM(YEAR, COUNTRY)                                                                               *
! *   PURPOSE: Calculate components of the country's pension system.                                                               *
! **********************************************************************************************************************************

      SUBROUTINE PENSION_SYSTEM(YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER GEN, IK, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_OF_RETIREMENT
           REAL*8 :: MAX_TAXABLE_EARNINGS, TDAMP, TDAM
           REAL*8, DIMENSION(Y_CLASSES + 1) :: EARNINGS_BASE_PAYROLL
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, Y_CLASSES) :: EARNINGS_BASE_PAYROLL_IND
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, SUM_WAGE
           
           TDAMP = 0.2D0
           TDAM = 0.2D0                          
           
!          Compute average indexed labor earnings
           CALL GET_AVG_INDEXED_EARNINGS(YEAR, COUNTRY)
           
!          Nur einschalten, wenn HEALTH_SYSTEM fehlt
!           CALL GET_AVG_LABOR_EARNINGS(YEAR, COUNTRY)

           DO Y_CLASS = 1, Y_CLASSES
                DO GEN = FIRST_WORK_YEAR, GENS                    
                     PENSION_REPLACEMENT_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                     PENSION_BENEFITS_IND(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                ENDDO
           ENDDO           
           
           PENSION_BENEFITS(YEAR, COUNTRY) = 0.
           EARNINGS_BASE_PAYROLL(Y_CLASSES+1) = 0.
           MAX_TAXABLE_EARNINGS = CONTRIBUTION_CEILING(YEAR, COUNTRY) * AVG_LABOR_EARNINGS(YEAR, COUNTRY)
           
           DO Y_CLASS = 1, Y_CLASSES
                EARNINGS_BASE_PAYROLL(Y_CLASS) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY)-1
                     EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) = SUM_WAGE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY))
                     IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) < 0.) EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) = 0.
                     IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) <= MAX_TAXABLE_EARNINGS) THEN
                          EARNINGS_BASE_PAYROLL(Y_CLASS) = EARNINGS_BASE_PAYROLL(Y_CLASS) + &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ELSE
                          EARNINGS_BASE_PAYROLL(Y_CLASS) = EARNINGS_BASE_PAYROLL(Y_CLASS) + MAX_TAXABLE_EARNINGS * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ENDIF
                ENDDO
                EARNINGS_BASE_PAYROLL(Y_CLASSES+1) = EARNINGS_BASE_PAYROLL(Y_CLASSES+1) + EARNINGS_BASE_PAYROLL(Y_CLASS)
           ENDDO          
!            
           DO Y_CLASS = 1, Y_CLASSES
                DO GEN = RETIREMENT_AGE(YEAR, COUNTRY), GENS
                     IK = GET_YEAR_OF_RETIREMENT(GEN, YEAR, COUNTRY)
                     IF(YEAR >= 0 .AND. IK < 0) IK = 0
                     PENSION_REPLACEMENT_RATE(GEN, IK, Y_CLASS, COUNTRY) = TAX_RATE(4, 1, IK, COUNTRY) - &
                          & TAX_RATE(4, 2, IK, COUNTRY) * AVG_INDEXED_EARNINGS(IK, Y_CLASS, COUNTRY)
                     IF(COUNTRY == 1 .AND. GEN < 63) PENSION_REPLACEMENT_RATE(GEN, IK, Y_CLASS, COUNTRY) = 0.0
                     IF(PENSION_REPLACEMENT_RATE(GEN, IK, Y_CLASS, COUNTRY) < 0) &
                          & PENSION_REPLACEMENT_RATE(GEN, IK, Y_CLASS, COUNTRY) = 0.
                     PENSION_BENEFITS_IND(GEN, YEAR, Y_CLASS, COUNTRY) = PENSION_REPLACEMENT_RATE(GEN, IK, Y_CLASS, COUNTRY) * &
                          & AVG_INDEXED_EARNINGS(IK, Y_CLASS, COUNTRY) * YY_0(IK, COUNTRY) / YY(IK, COUNTRY)
                     PENSION_BENEFITS(YEAR, COUNTRY) = PENSION_BENEFITS(YEAR, COUNTRY) + &
                          & PENSION_BENEFITS_IND(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                ENDDO
           ENDDO
!                
           IF(EARNINGS_BASE_PAYROLL(Y_CLASSES+1) /= 0.) THEN
                AGG_PENSION_TAX_RATE(YEAR, COUNTRY) = (1. - MU_1(YEAR, COUNTRY)) * PENSION_BENEFITS(YEAR, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           ELSE
                WRITE( * , 100) COUNTRY, YEAR
100             FORMAT(' PY = 0, COUNTRY = ', I1, ', YEAR = ', I3)
           ENDIF                
!           
           AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = 0.
           AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = 0.
           DO Y_CLASS = 1, Y_CLASSES
                AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) -1
                     IF(EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) <= MAX_TAXABLE_EARNINGS) THEN
                          IF(Y_CLASS < Y_CLASSES) THEN
                               MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                                    & MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1 -TDAMP) + &
                                    & TDAMP * AGG_PENSION_TAX_RATE(YEAR, COUNTRY)
                          ELSE
                               MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          ENDIF
                          AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                               & AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAMP) + &
                               & TDAMP * AGG_PENSION_TAX_RATE(YEAR, COUNTRY)
                          AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS) * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ELSE
                          MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = &
                               & AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * (1.-TDAM) + &
                               & TDAM * AGG_PENSION_TAX_RATE(YEAR, COUNTRY) * MAX_TAXABLE_EARNINGS / &
                               & EARNINGS_BASE_PAYROLL_IND(GEN, Y_CLASS)
                          AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * MAX_TAXABLE_EARNINGS * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                          AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) + &
                               & AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) * MAX_TAXABLE_EARNINGS * &
                               & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                     ENDIF
                     IF(AGG_PENSION_TAX_RATE(YEAR, COUNTRY) == 0.) THEN
                          MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                     ENDIF
                ENDDO
                AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                     & AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY)
                AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) + &
                     & AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY)
                AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASS)
                AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) / &
                     & EARNINGS_BASE_PAYROLL(Y_CLASS)
           ENDDO
           AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) / &
                & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)
           AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) = AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) / &
                & EARNINGS_BASE_PAYROLL(Y_CLASSES+1)          
           
           RETURN
      END SUBROUTINE PENSION_SYSTEM
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_AVG_INDEXED_EARNINGS(YEAR, COUNTRY)                                                                     *
! *   PURPOSE: Calculate average indexed earnings of agents in period YEAR in country LA.                                          *
! **********************************************************************************************************************************

      SUBROUTINE GET_AVG_INDEXED_EARNINGS(YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: GEN, IK, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J 
           REAL*8, EXTERNAL :: SUM_WAGE
           
           DO Y_CLASS = 1, Y_CLASSES            
                AVG_INDEXED_EARNINGS(YEAR, Y_CLASS, COUNTRY) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) -1                                     
                     IK = GET_YEAR_BECOMING_J(YEAR, RETIREMENT_AGE(YEAR, COUNTRY), GEN)
                     IF(YEAR >= 0 .AND. IK < 0) IK = 0                           
                     AVG_INDEXED_EARNINGS(YEAR, Y_CLASS, COUNTRY) = AVG_INDEXED_EARNINGS(YEAR, Y_CLASS, COUNTRY) + &
                          & (HOURS - LEISURE(GEN, IK, Y_CLASS, COUNTRY)) * SUM_WAGE(GEN, IK, Y_CLASS, COUNTRY)                                    
                ENDDO
                AVG_INDEXED_EARNINGS(YEAR, Y_CLASS, COUNTRY) = AVG_INDEXED_EARNINGS(YEAR, Y_CLASS, COUNTRY) / &
                     & (RETIREMENT_AGE(YEAR, COUNTRY) - FIRST_WORK_YEAR)
           ENDDO
           
           RETURN
      END SUBROUTINE GET_AVG_INDEXED_EARNINGS

! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE GET_AVG_LABOR_EARNINGS(YEAR, COUNTRY)                                                                       *
! *   PURPOSE: Calculate average labor earnings of working generations in period YEAR in country LA.                               *
! **********************************************************************************************************************************

      SUBROUTINE GET_AVG_LABOR_EARNINGS(YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: GEN, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, YEAR
           REAL*8 :: WORKPOP
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, SUM_WAGE
           
           AVG_LABOR_EARNINGS(YEAR, COUNTRY) = 0.
           WORKPOP = 0.
           
           DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) -1
                DO Y_CLASS = 1, Y_CLASSES
                     AVG_LABOR_EARNINGS(YEAR, COUNTRY) = AVG_LABOR_EARNINGS(YEAR, COUNTRY) + &
                          & (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY)) * &
                          & SUM_WAGE(GEN, YEAR, Y_CLASS, COUNTRY) * GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                ENDDO
                WORKPOP = WORKPOP + POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) / POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
           ENDDO
           
           AVG_LABOR_EARNINGS(YEAR, COUNTRY) = AVG_LABOR_EARNINGS(YEAR, COUNTRY) / WORKPOP
           
           RETURN
      END SUBROUTINE GET_AVG_LABOR_EARNINGS

! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_GOODS_MARKET(YEAR)                                                                                        *
! *   PURPOSE: Calculate aggregate supply, demand, and investment.                                                                 *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_GOODS_MARKET(YEAR)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, IPL1, Y_CLASS

           INTEGER, INTENT(IN) :: YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           IPL1 = GET_YEAR_BECOMING_J(YEAR, 0, 1)
           GET_GOODS_MARKET = 0.
           DD_WORLD(YEAR) = 0.
           YY_WORLD(YEAR) = 0.
           NET_ENDOW_INVEST(YEAR) = 0


		
		!What follows is our Market clearing condition%$%$%$%$%

		   DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
!               Aggregate Supply
                YY(YEAR, COUNTRY) = 0.
                NATIONAL_INCOME(YEAR, COUNTRY) = 0.
                GDP(YEAR, COUNTRY) = 0.


                
                DO Y_CLASS = 1, Y_CLASSES                  
                     YY(YEAR, COUNTRY) = YY(YEAR, COUNTRY) + LABOR(YEAR, Y_CLASS, COUNTRY) * &
                          & WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY)
                          
					 NATIONAL_INCOME(YEAR, COUNTRY) = NATIONAL_INCOME(YEAR, COUNTRY) +&
                     &LABOR(YEAR, Y_CLASS, COUNTRY) * WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY)
                     
!$$$$$$                      GDP(YEAR, COUNTRY) = GDP(YEAR, COUNTRY) +&
!$$$$$$                      &LABOR(YEAR, Y_CLASS, COUNTRY) * WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY)
!$$$$$$                           
                ENDDO     
                YY(YEAR, COUNTRY) = YY(YEAR, COUNTRY) + CAPITAL(YEAR, COUNTRY) * R(YEAR, COUNTRY) + &
                     & AGG_ASSETS_MIGRANTS(YEAR, Y_CLASSES+1, COUNTRY)

				NATIONAL_INCOME(YEAR, COUNTRY) = NATIONAL_INCOME(YEAR, COUNTRY) + &
                &AGG_ASSETS(YEAR, COUNTRY) * RG(YEAR) + &
                &AGG_ASSETS_MIGRANTS(YEAR, Y_CLASSES+1, COUNTRY) +&
                &GOV_ENDOW_SHARE(YEAR, COUNTRY)*ENDOWMENT(YEAR) -&
                &RG(YEAR)*DEBT(YEAR, COUNTRY)

                
                     
                     
!$$$$$$                 GDP(YEAR, COUNTRY) = GDP(YEAR, COUNTRY) +&
!$$$$$$                      &CAPITAL(YEAR, COUNTRY) * R(YEAR, COUNTRY) + &
!$$$$$$                      &COUNTRY_ENDOW_SHARE(YEAR, COUNTRY)*ENDOWMENT(YEAR) 


                     
                IF(IRUN == 0)YY_0(YEAR, COUNTRY) = YY(YEAR, COUNTRY)
!               Bruttoinvestitionen, Formel (86):
                INVEST(YEAR, COUNTRY) = (1.+TECH) * (1.+NPOP(IPL1, FIRST_COUNTRY)) * CAPITAL(IPL1, COUNTRY) - &
                     & (1.-DEL) * CAPITAL(YEAR, COUNTRY)
!               Aggregate Demand, Formel (84):
                DD(YEAR, COUNTRY) = CC(YEAR, COUNTRY) + INVEST(YEAR, COUNTRY) + GOVERNMENT_EXPENDITURES(YEAR, COUNTRY)
                DD_WORLD(YEAR) = DD_WORLD(YEAR) + DD(YEAR, COUNTRY)
                YY_WORLD(YEAR) = YY_WORLD(YEAR) + YY(YEAR, COUNTRY)
           ENDDO
           
               NET_ENDOW_INVEST(YEAR) = (1.+TECH) * (1.+NPOP(IPL1, FIRST_COUNTRY)) * PVENDOWMENT(IPL1) - PVENDOWMENT(YEAR)
               TOT_GOV_ENDOW_REVENUE(YEAR) = TOT_GOV_ENDOW_SHARE(YEAR) * ENDOWMENT(YEAR)


           !Get Markets to Balance:
           
		   YY_WORLD(YEAR) = YY_WORLD(YEAR) + PVENDOWMENT(YEAR)*RG(YEAR) + TOT_GOV_ENDOW_REVENUE(YEAR)

           DD_WORLD(YEAR) = DD_WORLD(YEAR) + NET_ENDOW_INVEST(YEAR)
           
			IF(DABS(YY_WORLD(YEAR) - DD_WORLD(YEAR)) < SIGFIG) THEN
             GET_GOODS_MARKET = 1
             write(*,*) 'Good Year', YEAR      
           ELSE
           !  write(*,*) YEAR, NET_ENDOW_INVEST(YEAR), TECH
             write(*,*) 'Bad Year, This Much Excess Supply:', (YY_WORLD(YEAR) - DD_WORLD(YEAR))
           !  write(*,*) 'This Much PVENDOWMENT and ENDOWMENT:', PVENDOWMENT(YEAR), ENDOWMENT(YEAR)
           !  write(*,*) 'This Much Interest and Change In Previous Endowment', RG(YEAR), (PVENDOWMENT(IPL1) - PVENDOWMENT(YEAR))
            ENDIF


        

      END FUNCTION GET_GOODS_MARKET
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)                                                         *
! *   PURPOSE:                                                                                                                     *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, Y_CLASS, YEAR
           
           GET_EFFICIENT_POPULATION = (1.+TECH) ** (FIRST_WORK_YEAR - GEN) * POP(GEN, YEAR, Y_CLASS, COUNTRY) / &
                & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY)
                
           RETURN
      END FUNCTION GET_EFFICIENT_POPULATION
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION KIDS(GEN, YEAR, Y_CLASS, COUNTRY)                                                                             *
! *   PURPOSE: Compute number of kids of parents age GEN in YEAR; Formel (12)                                                      *
! **********************************************************************************************************************************

      REAL*8 FUNCTION KIDS(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: J, UU, MM
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
!           
           KIDS = 0.
!           
           IF (GEN <= 46) THEN
                UU = 1
           ELSE
                UU = GEN - 45
           ENDIF
           IF (GEN - FIRST_FERTILITY_YEAR < 20) THEN
                MM = GEN - FIRST_FERTILITY_YEAR
           ELSE
                MM = 20
           ENDIF           
           
           IF(GEN <= 22 .OR. GEN >= 66) KIDS = 0.
             
           IF(GEN >= FIRST_FERTILITY_YEAR .AND. GEN <= 65) THEN
                IF (GEN <= 45 .AND. YEAR > -1) KIDS = FERTILITY(GEN, YEAR, COUNTRY)                  
                IF (GEN <= 45 .AND. YEAR == -1) KIDS = FERTILITY_SS(GEN, COUNTRY)                      
                DO J = UU, MM
                     IF (YEAR > -1) KIDS = KIDS + POP(J, YEAR, Y_CLASS, COUNTRY) * FERTILITY(GEN-J, YEAR-J, COUNTRY) / &
                          & TOTAL_FERTILITY(YEAR-J, COUNTRY) / POP(GEN, YEAR, Y_CLASS, COUNTRY)
                     IF (YEAR == -1) KIDS = KIDS + POP(J, YEAR, Y_CLASS, COUNTRY) * FERTILITY_SS(GEN-J, COUNTRY) / &
                          & TOTAL_FERTILITY_SS(COUNTRY) / POP(GEN, YEAR, Y_CLASS, COUNTRY)     
                ENDDO
           ENDIF
           
           RETURN
      END FUNCTION KIDS
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION SUM_WAGE(J, IK, Y_CLASS, COUNTRY)                                                                             *
! *   PURPOSE: SUM_WAGE adds the wage and the shadow wage.                                                                         *
! **********************************************************************************************************************************

      REAL*8 FUNCTION SUM_WAGE(J, IK, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: COUNTRY, IK, J, Y_CLASS
           REAL*8, EXTERNAL :: GET_WAGE
           
           SUM_WAGE = GET_WAGE(J, IK, Y_CLASS, COUNTRY) + SHADOW_WAGE(J, IK, Y_CLASS, COUNTRY)
           
           RETURN
      END FUNCTION SUM_WAGE

! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_WAGE(GEN, YEAR, Y_CLASS, COUNTRY)                                                                         *
! *   PURPOSE: Compute wage profile.                                                                                               *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_WAGE(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
           
           GET_WAGE = WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY) * AGE_EFFICIENCY(GEN, YEAR, Y_CLASS, COUNTRY)
           
           RETURN
      END FUNCTION GET_WAGE
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_INHERITANCES(GEN, YEAR, Y_CLASS, COUNTRY)                                                                 *
! *   PURPOSE: Calculate inheritances                                                                                              *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_INHERITANCES(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
           REAL*8 SHR(21:49)
           DATA SHR/0.137D0,0.139D0,0.256D0,0.451D0,0.757D0,1.207D0,1.834D0,2.652D0,3.650D0,4.784D0,&
                5.969D0,7.091D0,8.018D0,8.632D0,8.847D0,8.632D0,8.018D0,7.091D0,5.969D0,4.784D0,&
                3.650D0,2.652D0,1.834D0,1.207D0,0.757D0,0.451D0,0.256D0,0.139D0,0.136D0/
           
           IF(GEN >= 50) THEN
                GET_INHERITANCES = 0.
           ELSE
                GET_INHERITANCES = SHR(GEN)/100 * AGG_ASSETS_FOR_BEQUESTS(YEAR, Y_CLASS, COUNTRY) / &
                     & POP(GEN, YEAR, Y_CLASS, COUNTRY) * (1.+TECH) ** (GEN - FIRST_WORK_YEAR)
           ENDIF
           
           RETURN
      END FUNCTION GET_INHERITANCES
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION KIDS_HEALTH_BENEFITS(GEN, YEAR, Y_CLASS, COUNTRY)                                                             *
! *   PURPOSE: Compute health benefits for kids of parents age GEN in YEAR                                                         *
! **********************************************************************************************************************************

      REAL*8 FUNCTION KIDS_HEALTH_BENEFITS(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: J, UU, MM
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
           
           KIDS_HEALTH_BENEFITS = 0.
!
           IF (GEN <= 46) THEN
                UU = 1
           ELSE
                UU = GEN - 45
           ENDIF
           IF (GEN - FIRST_FERTILITY_YEAR < 20) THEN
                MM = GEN - FIRST_FERTILITY_YEAR
           ELSE
                MM = 20
           ENDIF
           
           IF(GEN <= 22 .OR. GEN >= 66) KIDS_HEALTH_BENEFITS = 0.
             
           IF(GEN >= FIRST_FERTILITY_YEAR .AND. GEN <= 65) THEN
                IF (GEN <= 45 .AND. YEAR > -1) KIDS_HEALTH_BENEFITS = FERTILITY(GEN, YEAR, COUNTRY) * &
                     & HEALTH_BENEFITS_IND(0, YEAR, COUNTRY)                  
                IF (GEN <= 45 .AND. YEAR == -1) KIDS_HEALTH_BENEFITS = FERTILITY_SS(GEN, COUNTRY) * &
                     & HEALTH_BENEFITS_IND(0, YEAR, COUNTRY)                      
                DO J = UU, MM
                     IF (YEAR > -1) KIDS_HEALTH_BENEFITS = KIDS_HEALTH_BENEFITS + POP(J, YEAR, Y_CLASS, COUNTRY) * &
                          & FERTILITY(GEN-J, YEAR-J, COUNTRY) / TOTAL_FERTILITY(YEAR-J, COUNTRY) / &
                          & POP(GEN, YEAR, Y_CLASS, COUNTRY) * HEALTH_BENEFITS_IND(J, YEAR, COUNTRY)
                     IF (YEAR == -1) KIDS_HEALTH_BENEFITS = KIDS_HEALTH_BENEFITS + POP(J, YEAR, Y_CLASS, COUNTRY) * &
                          & FERTILITY_SS(GEN-J, COUNTRY) / TOTAL_FERTILITY_SS(COUNTRY) / POP(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & HEALTH_BENEFITS_IND(J, YEAR, COUNTRY)
                ENDDO
           ENDIF
           
           RETURN
      END FUNCTION KIDS_HEALTH_BENEFITS
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE READ_POPULATION_DATA                                                                                        *
! *   PURPOSE:                                                                                                                     *
! **********************************************************************************************************************************

      SUBROUTINE READ_POPULATION_DATA
           USE GLOBAL_DATA
           IMPLICIT NONE
           REAL*8, DIMENSION(0:GENS + 1,COUNTRIES) :: HPOP
           INTEGER COUNTRY, GEN, Y_CLASS, YEAR
           
           OPEN (UNIT = 3, FILE = INFILE_NAME(1), STATUS = 'OLD')
           
           DO COUNTRY = 1, COUNTRIES
                DO GEN = 0, 91
                     IF(GEN <= GENS) MIGRANTS(GEN, COUNTRY) = 0.
                     DO YEAR = -1, YEARS
                          IF(GEN <= GENS) MORTALITY(GEN, YEAR, COUNTRY) = 0.
                     ENDDO
                ENDDO
                DO YEAR = - GENS, YEARS
                     DO GEN = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                          FERTILITY(GEN, YEAR, COUNTRY) = 0.
                     ENDDO
                ENDDO
           ENDDO
           
!          Einlesen der Gesamtbevölkerung im Jahr 0
           READ(3, 100)
100        FORMAT( / )
           DO GEN = 0, GENS
                READ(3, 110) HPOP(GEN, 1), HPOP(GEN, 2), HPOP(GEN, 3), HPOP(GEN, 4), HPOP(GEN, 5), HPOP(GEN,6)
110             FORMAT(7X, 3(F11.6, 5X), 3(F12.6, 5X))
           ENDDO
!           
!          Einlesen der Immigranten
           READ(3, 200)
200        FORMAT( // )
1000	   FORMAT (7X, 3(F11.6, 5X), 3(F12.6, 5X))
           DO GEN = 1, 65
                READ(3, 1000) MIGRANTS(GEN, 1), MIGRANTS(GEN, 2), MIGRANTS(GEN, 3), MIGRANTS(GEN, 4)&
                &, MIGRANTS(GEN, 5), MIGRANTS(GEN,6)
           ENDDO
!
!          Einlesen der Fertilitätsraten
           READ(3, 300)
300        FORMAT(4( / ))
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 1), GEN = FIRST_FERTILITY_YEAR, 30)
310             FORMAT(5X, 8(F7.6, 1X))
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 1), GEN = 31, 38)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 1), GEN = 39, LAST_FERTILITY_YEAR)
           ENDDO
           READ(3, 320)
320        FORMAT(3( / ))
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 2), GEN = FIRST_FERTILITY_YEAR, 30)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 2), GEN = 31, 38)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 2), GEN = 39, LAST_FERTILITY_YEAR)
           ENDDO
!                        
           READ(3, 320)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 3), GEN = FIRST_FERTILITY_YEAR, 30)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 3), GEN = 31, 38)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 3), GEN = 39, LAST_FERTILITY_YEAR)
           ENDDO
           READ(3, 320)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 4), GEN = FIRST_FERTILITY_YEAR, 30)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 4), GEN = 31, 38)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 4), GEN = 39, LAST_FERTILITY_YEAR)
           ENDDO
           READ(3, 320)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 5), GEN = FIRST_FERTILITY_YEAR, 30)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 5), GEN = 31, 38)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 5), GEN = 39, LAST_FERTILITY_YEAR)
           ENDDO

           READ(3, 320)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 6), GEN = FIRST_FERTILITY_YEAR, 30)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 6), GEN = 31, 38)
           ENDDO
           READ(3, 100)
           DO YEAR = -48, TRANS_YEAR
                READ(3, 310) (FERTILITY(GEN, YEAR, 6), GEN = 39, LAST_FERTILITY_YEAR)
           ENDDO

          FERTILITY(:,:,6) = .90D0*FERTILITY(:,:,6) ! Russian fertility devided by 2 because original data is woman fertility rates
!           
!          Einlesen der Sterbewahrscheinlichkeiten
           READ(3, 300)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 1), GEN = 68, 75)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 1), GEN = 76, 83)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 1), GEN = 84, GENS)
           ENDDO
           READ(3, 320)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 2), GEN = 68, 75)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 2), GEN = 76, 83)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 2), GEN = 84, GENS)
           ENDDO          
           READ(3, 320)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 3), GEN = 68, 75)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 3), GEN = 76, 83)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 3), GEN = 84, GENS)
           ENDDO
           READ(3, 320)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 4), GEN = 68, 75)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 4), GEN = 76, 83)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 4), GEN = 84, GENS)
           ENDDO
           READ(3, 320)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 5), GEN = 68, 75)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 5), GEN = 76, 83)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 5), GEN = 84, GENS)
           ENDDO


	       READ(3, 320)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 6), GEN = 68, 75)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 6), GEN = 76, 83)
           ENDDO
           READ(3, 100)
           DO YEAR = 0, TRANS_YEAR
                READ(3, 310) (MORTALITY(GEN, YEAR, 6), GEN = 84, GENS)
           ENDDO
!           
!          Einlesen des Gesundheitsausgabenprofils
           READ(3, 200)
           DO GEN = 0, GENS
                READ(3, 400) HEALTH_EXPENDITURES_PROFILE(GEN, 1), HEALTH_EXPENDITURES_PROFILE(GEN, 2), &
                     & HEALTH_EXPENDITURES_PROFILE(GEN, 3), HEALTH_EXPENDITURES_PROFILE(GEN, 6)
400             FORMAT(7X, F5.2, 6x, F9.7, 6X, F9.7, 7X, F9.7)
           ENDDO
           
!          ANPASSUNGEN
           DO COUNTRY = 1, COUNTRIES
                DO GEN = 0, GENS
                     HPOP(GEN, COUNTRY) = HPOP(GEN, COUNTRY)*1000.
                     IF(GEN <= 65)MIGRANTS(GEN, COUNTRY) = MIGRANTS(GEN, COUNTRY)*100.
                     DO YEAR = TRANS_YEAR + 1, YEARS
                          MORTALITY(GEN, YEAR, COUNTRY) = MORTALITY(GEN, TRANS_YEAR, COUNTRY)
                     ENDDO
                ENDDO
                DO GEN = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                     DO YEAR = - GENS, - 49
                          FERTILITY(GEN, YEAR, COUNTRY) = FERTILITY(GEN, - 48, COUNTRY)
                     ENDDO
                     DO YEAR = TRANS_YEAR + 1, YEARS
                          FERTILITY(GEN, YEAR, COUNTRY) = FERTILITY(GEN, TRANS_YEAR, COUNTRY)
                     ENDDO
                ENDDO
           ENDDO    
           
!          Berechnung der Fertilitätsanteile und des durchschnittlichen Gebäralters
!          TOTAL_FERTILITY = Zus.gefaßte Geb.ziffer = Summe der altersspezifischen Geburtsziffern in YEAR( = Zahl der Kinder, die
!          1000 MENSCHEN zwischen im Laufe ihres Lebens gebären, u.d.A.daß sich die altersspezifischen Geb.ziffern nicht ändern)
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = -GENS, YEARS
                     AVG_BIRTH_AGE(YEAR, COUNTRY) = 0.
                     TOTAL_FERTILITY(YEAR, COUNTRY) = 0.
                     DO GEN = FIRST_FERTILITY_YEAR, LAST_FERTILITY_YEAR
                          TOTAL_FERTILITY(YEAR, COUNTRY) = TOTAL_FERTILITY(YEAR, COUNTRY) + FERTILITY(GEN, YEAR, COUNTRY)
                          AVG_BIRTH_AGE(YEAR, COUNTRY) = AVG_BIRTH_AGE(YEAR, COUNTRY) + FERTILITY(GEN, YEAR, COUNTRY) * GEN
                     ENDDO
                     AVG_BIRTH_AGE(YEAR, COUNTRY) = AVG_BIRTH_AGE(YEAR, COUNTRY) / TOTAL_FERTILITY(YEAR, COUNTRY)
                ENDDO
           ENDDO
           
!          Aufteilung der Bevölkerung im Jahr 0 auf die Einkommensklassen
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO GEN = 0, GENS
                     POP(GEN, 0, Y_CLASSES+1, COUNTRY) = HPOP(GEN, COUNTRY) * POP_SCALE(COUNTRY)
                     DO Y_CLASS = 1, Y_CLASSES
                          POP(GEN, 0, Y_CLASS, COUNTRY) = C_SHARE(Y_CLASS, COUNTRY) * HPOP(GEN, COUNTRY) * POP_SCALE(COUNTRY)
                     ENDDO
                ENDDO
           ENDDO
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                POP(91, 0, Y_CLASSES+1, COUNTRY) = 0.
                DO GEN = 0, GENS
                     POP(91, 0, Y_CLASSES+1, COUNTRY) = POP(91, 0, Y_CLASSES+1, COUNTRY) + POP(GEN, 0, Y_CLASSES+1, COUNTRY)
                ENDDO
           ENDDO 
!
         

!!!!!!!!!!READ IN ENDOWMENTS %$%$%$%$%$%$
           OPEN (UNIT = 7, FILE = INFILE_NAME(2), STATUS = 'OLD')
  
			ENDOWMENT(:) = 0

		   READ(7, 2000)
2000       FORMAT( / )
           
           DO YEAR = 0, 300
			READ(7, 1100) ENDOWMENT(YEAR)

1100       FORMAT(7X, F11.6)
 		   ENDDO
  
           !!!Scale endowments: Change this value to add more oil to the economy!!!
           ENDOWMENT(:) = ENDOWMENT(:)*4.30D0 !ENDOWMENT(:) = ENDOWMENT(:)*4.70D0

           !Reduce endow in IC calc
		  ! ENDOWMENT(0) = 1.0D0
           
          
			

CLOSE(7)

  CLOSE(3)       
           
           RETURN
      END SUBROUTINE READ_POPULATION_DATA

! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_YEAR_IN_POP_CALC(YEAR, GEN, JJ)                                                                           *
! *   PURPOSE: Calculate the year in which someone age GEN in year YEAR will turn JJ                                               *
! **********************************************************************************************************************************

      INTEGER FUNCTION GET_YEAR_IN_POP_CALC(YEAR, GEN, JJ)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: GEN, JJ, YEAR
           
           GET_YEAR_IN_POP_CALC = YEAR + JJ - GEN
           IF(GET_YEAR_IN_POP_CALC > YEARS) GET_YEAR_IN_POP_CALC = YEARS
             
           RETURN
      END FUNCTION GET_YEAR_IN_POP_CALC

! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE PRINT_POPULATION_SUMMARY                                                                                    *
! *   PURPOSE:                                                                                                                     *
! **********************************************************************************************************************************

      SUBROUTINE PRINT_POPULATION_SUMMARY
           USE GLOBAL_DATA
           IMPLICIT NONE
           CHARACTER (LEN = 13) :: COUNTRY_NAME
           INTEGER :: COUNTRY, GEN, YEAR
           REAL*8, DIMENSION(0:YEARS, COUNTRIES) :: LE, TDM
           REAL*8, EXTERNAL :: GET_POPULATION_SHARE
           
           IF(IRUN == 0) THEN
                OPEN ( UNIT = 4, FILE = OUTFILE_NAME(1), STATUS = 'REPLACE')
           ELSE
                OPEN ( UNIT = 4, FILE = OUTFILE_NAME(2), STATUS = 'REPLACE')
           ENDIF
           
!          Berechnung der Lebenserwartung
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = 0, YEARS
                     LE(YEAR, COUNTRY) = 0.
                     DO GEN = 0, GENS
                          LE(YEAR, COUNTRY) = LE(YEAR, COUNTRY) + SURVIVAL_PROBABILITY(GEN, YEAR, COUNTRY)
                     ENDDO
                ENDDO
           ENDDO
           
!          Berechnung der Gesamtzahl der Immigranten
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                DO YEAR = 0, YEARS
                     TDM(YEAR, COUNTRY) = 0.
                     DO GEN = 1, 65
                          IF(YEAR <= TRANS_YEAR) THEN
                               TDM(YEAR, COUNTRY) = TDM(YEAR, COUNTRY) + MIGRATION_SCALE(GEN, YEAR, COUNTRY) * &
                                    & MIGRANTS(GEN, COUNTRY)
                          ELSE
                               TDM(YEAR, COUNTRY) = TDM(YEAR, COUNTRY) + (1.+EXOG_NPOP) ** (YEAR - TRANS_YEAR) * &
                                    & MIGRATION_SCALE(GEN, YEAR, COUNTRY) * MIGRANTS(GEN, COUNTRY)
                          ENDIF
                     ENDDO
               ENDDO
           ENDDO
           
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                IF(COUNTRY == 1) COUNTRY_NAME = 'USA'
                IF(COUNTRY == 2) COUNTRY_NAME = 'EUROPA'
                IF(COUNTRY == 3) COUNTRY_NAME = 'JAPAN'
                IF(COUNTRY == 4) COUNTRY_NAME = 'CHINA'
                IF(COUNTRY == 5) COUNTRY_NAME = 'INDIA'  
                IF(COUNTRY == 6) COUNTRY_NAME = 'RUSSIA'
                WRITE(4, 112) COUNTRY_NAME
112             FORMAT('                              BEVÖLKERUNGSENTWICKLUNG ', A8 / '                               ' &
                     & '=  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = ' / )
!               OUTPUTTABELLE
                WRITE(4, 1221)
1221            FORMAT('YEAR          2008      2010      2020      2030      2040      2050      2060      2070      2100')
!               Life EXpectancies in every year
                WRITE(4, 1219)
1219            FORMAT( / 'LIFE EXPECTANCY AT BIRTH')
!$$$$$$                 WRITE(4, 1222) LE(FIRST_YEAR-2000, COUNTRY), (LE(YEAR, COUNTRY), YEAR = 10, 70, 10), LE(100, COUNTRY)
				WRITE(4, 1222) LE(0, COUNTRY), (LE(YEAR, COUNTRY), YEAR = 2, 62, 10), LE(92, COUNTRY)
1222            FORMAT('MODELL       ', 9(F5.2, 5X))
!               Periodenfertilität pro Frau
                WRITE(4, 1320)
1320            FORMAT( / 'TOTAL FERTILITY RATE (Periodenfertilität)')
!$$$$$$                 WRITE(4, 1424) 2 * TOTAL_FERTILITY(FIRST_YEAR-2000, COUNTRY), (2 * TOTAL_FERTILITY(YEAR, COUNTRY), &
!$$$$$$                      & YEAR = 10, 70, 10), 2 * TOTAL_FERTILITY(100, COUNTRY)
!$$$$$$ 1424            FORMAT('MODELL  ', 9(F10.2))
				WRITE(4, 1424) 2 * TOTAL_FERTILITY(0, COUNTRY), (2 * TOTAL_FERTILITY(YEAR, COUNTRY), &
                     & YEAR = 2, 62, 10), 2 * TOTAL_FERTILITY(92, COUNTRY)
1424            FORMAT('MODELL  ', 9(F10.2))           
!               Durchschnittliches Gebäralter
                WRITE(4, 1110)
1110            FORMAT( / 'DURCHSCHNITTLICHES GEBÄRALTER')
!$$$$$$                 WRITE(4, 1427) AVG_BIRTH_AGE(FIRST_YEAR-2000, COUNTRY), (AVG_BIRTH_AGE(YEAR, COUNTRY), YEAR = 10, 70, 10), &
!$$$$$$                      & AVG_BIRTH_AGE(100, COUNTRY)
                WRITE(4, 1427) AVG_BIRTH_AGE(0, COUNTRY), (AVG_BIRTH_AGE(YEAR, COUNTRY), YEAR = 2, 62, 10), &
                     & AVG_BIRTH_AGE(92, COUNTRY)     
1427            FORMAT('MODELL    ', 9(F8.2, 2X))
!               Bevölkerungsanteile und Dependency Ratios
                WRITE(4, 1930)
1930            FORMAT( / '                ANTEILE BESTIMMTER ALTERSGRUPPEN AN DER GESAMTBEVOELKERUNG' / )
                WRITE(4, 1830) GET_POPULATION_SHARE(FIRST_YEAR, 0, 14, COUNTRY), &
                     & (GET_POPULATION_SHARE(2000+YEAR, 0, 14, COUNTRY), YEAR = 10, 70, 10), &
                     & GET_POPULATION_SHARE(2100, 0, 14, COUNTRY)
                     IF (COUNTRY<=5) THEN
               WRITE(4, 1831)GET_POPULATION_SHARE(FIRST_YEAR, 23, 45, COUNTRY), &
                    & (GET_POPULATION_SHARE(2000+YEAR, 23, 45, COUNTRY), YEAR = 10, 70, 10), &
                    & GET_POPULATION_SHARE(2100, 23, 45, COUNTRY)
                    ELSEIF (COUNTRY ==6) THEN
               WRITE(4, 1831)GET_POPULATION_SHARE(FIRST_YEAR, 23, 45, COUNTRY), &
                    & (GET_POPULATION_SHARE(2000+YEAR, 23, 45, COUNTRY), YEAR = 10, 70, 10), &
                    & GET_POPULATION_SHARE(2100, 23, 45, COUNTRY)
                    ENDIF
               WRITE(4, 1832)GET_POPULATION_SHARE(FIRST_YEAR, 65, 90, COUNTRY), &
                    & (GET_POPULATION_SHARE(2000+YEAR, 65, 90, COUNTRY), YEAR = 10, 70, 10), &
                    & GET_POPULATION_SHARE(2100, 65, 90, COUNTRY)
!              Old age dependency ratio
               WRITE(4, 1231)GET_POPULATION_SHARE(FIRST_YEAR, 60, 90, COUNTRY) / &
                    & GET_POPULATION_SHARE(FIRST_YEAR, 15, 59, COUNTRY)*100, (GET_POPULATION_SHARE(2000+YEAR, 60, 90, COUNTRY) / &
                    & GET_POPULATION_SHARE(2000+YEAR, 15, 59, COUNTRY)*100, YEAR = 10, 70, 10), &
                    & GET_POPULATION_SHARE(2100, 60, 90, COUNTRY) / GET_POPULATION_SHARE(2100, 15, 59, COUNTRY)*100

!$$$$$$                     IF (COUNTRY<=5) THEN
!$$$$$$                WRITE(4, 1232) POP(91, FIRST_YEAR-2000, Y_CLASSES+1, COUNTRY) / 1000000, &
!$$$$$$                     & (POP(91, YEAR, Y_CLASSES+1, COUNTRY) / 1000000, YEAR = 10, 70, 10), &
!$$$$$$                     & POP(91, 100, Y_CLASSES+1, COUNTRY) / 1000000
!$$$$$$                     ELSEIF (COUNTRY ==6) THEN
               WRITE(4, 1232) POP(91, 0, Y_CLASSES+1, COUNTRY) / 1000000, &
                    & (POP(91, YEAR, Y_CLASSES+1, COUNTRY) / 1000000, YEAR = 2, 62, 10), &
                    & POP(91, 92, Y_CLASSES+1, COUNTRY) / 1000000
                    !ENDIF
1830           FORMAT('< 15    ', 9F10.2)
1831           FORMAT('23-45   ', 9F10.2)
1832           FORMAT('> 64    ', 9F10.2)
1232           FORMAT('POP(MIO)', 9F10.2)
1231           FORMAT('DEP.RAT ', 9F10.2 / )
               WRITE(4, 1220)
1220           FORMAT('                                      MIGRATION    ' / )
               WRITE(4, 5221)
5221           FORMAT('YEAR        2000       2010       2020       2030       2040       2050       2060       2070       2100')
               WRITE(4, 122) TDM(FIRST_YEAR-2000, COUNTRY) / 1000, ((TDM(YEAR, COUNTRY) / 1000), YEAR = 10, 70, 10), &
                    & TDM(100, COUNTRY) / 1000
122            FORMAT('TDM(TSD)', 1X, F7.2, 8(F11.2) / )
           ENDDO
           
           CLOSE(4)
           
           RETURN
      END SUBROUTINE PRINT_POPULATION_SUMMARY

! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_POPULATION_SHARE(IYR, IAU, IAO, COUNTRY)                                                                  *
! *   PURPOSE:                                                                                                                     *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_POPULATION_SHARE(IYR, IAU, IAO, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: IAYR, I
           INTEGER, INTENT(IN) :: IYR, IAU, IAO, COUNTRY
           
!          Model is normalized to 2000,  2000 = 0
           IAYR = IYR - 2000
           GET_POPULATION_SHARE = 0.0
           DO I = IAU, IAO
                GET_POPULATION_SHARE = GET_POPULATION_SHARE + POP(I, IAYR, Y_CLASSES+1, COUNTRY)
           ENDDO
           GET_POPULATION_SHARE = GET_POPULATION_SHARE / POP(91, IAYR, Y_CLASSES+1, COUNTRY)*100
           
           RETURN
      END FUNCTION GET_POPULATION_SHARE  

! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_YEAR_BECOMING_J(YEAR, GEN, JJ)                                                                            *
! *   PURPOSE: GET_YEAR_BECOMING_J is the year in which someone age GEN in year YEAR will turn JJ.                                 *
! **********************************************************************************************************************************

      INTEGER FUNCTION GET_YEAR_BECOMING_J(YEAR, GEN, JJ)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: GEN, JJ, YEAR
           
           GET_YEAR_BECOMING_J = YEAR + JJ - GEN
           IF(GET_YEAR_BECOMING_J > YEARS) GET_YEAR_BECOMING_J = YEARS
           IF(GET_YEAR_BECOMING_J < -1) GET_YEAR_BECOMING_J = -1
           IF(YEAR == -1) GET_YEAR_BECOMING_J = -1
           IF(YEAR == YEARS) GET_YEAR_BECOMING_J = YEARS
             
           RETURN
      END FUNCTION GET_YEAR_BECOMING_J
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_YEAR_OF_RETIREMENT(GEN, YEAR, COUNTRY)                                                                    *
! *   PURPOSE: GET_YEAR_OF_RETIREMENT is the year in which someone age GEN in year YEAR has retired.                               *
! **********************************************************************************************************************************

      INTEGER FUNCTION GET_YEAR_OF_RETIREMENT(GEN, YEAR, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: IM, IK1, IPL, KK
           INTEGER, INTENT(IN) :: COUNTRY, GEN, YEAR
           
           IF(GEN >= RETIREMENT_AGE(YEAR, COUNTRY)) THEN
                IM = GEN - 1
                IK1 = YEAR - 1                                    
                DO KK = IM, FIRST_WORK_YEAR, -1
                     IF(IK1 < -1)IK1 = -1
                     IF(RETIREMENT_AGE(IK1, COUNTRY) <= KK) IK1 = IK1 - 1
                ENDDO
                GET_YEAR_OF_RETIREMENT = IK1 + 1
           ELSE
                IPL = GEN + 1
                IK1 = YEAR + 1
                DO KK = IPL, GENS
                     IF(RETIREMENT_AGE(IK1, COUNTRY) >= KK) IK1 = IK1 + 1
                ENDDO
                GET_YEAR_OF_RETIREMENT = IK1 - 1
           ENDIF
           
           IF(GET_YEAR_OF_RETIREMENT > YEARS) GET_YEAR_OF_RETIREMENT = YEARS
           IF(GET_YEAR_OF_RETIREMENT < -1) GET_YEAR_OF_RETIREMENT = -1
           IF(YEAR == -1) GET_YEAR_OF_RETIREMENT = -1
           IF(YEAR == YEARS) GET_YEAR_OF_RETIREMENT = YEARS
             
           RETURN
      END FUNCTION GET_YEAR_OF_RETIREMENT      
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE OUTPUT(YEAR, I1, I2, COUNTRY)                                                                               *
! *   PURPOSE: Write output.                                                                                                       *
! **********************************************************************************************************************************

      SUBROUTINE OUTPUT(YEAR, I1, I2, COUNTRY)  
           USE GLOBAL_DATA
           IMPLICIT NONE
           CHARACTER (LEN = 10) :: NNAME
           CHARACTER (LEN = 10) :: COUNTRY_NAME
           INTEGER :: GEN, IPL1, JJ, Y_CLASS
           INTEGER, INTENT(IN) :: COUNTRY, I1, I2, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8 :: BUGO, DIFF, POPR, REPLA
           REAL*8, DIMENSION(Y_CLASSES+1) :: INC, POPP
           REAL*8, DIMENSION (FIRST_COUNTRY:LAST_COUNTRY) :: HPOP
           REAL*8, EXTERNAL :: GET_INHERITANCES, GET_LIFETIME_UTILITY, GET_WAGE, KIDS, KIDS_HEALTH_BENEFITS, SUM_WAGE
           
           NNAME = 'YEAR'
           IF(COUNTRY == 1) COUNTRY_NAME = 'USA'
           IF(COUNTRY == 2) COUNTRY_NAME = 'EUROPA'
           IF(COUNTRY == 3) COUNTRY_NAME = 'JAPAN'
           IF(COUNTRY == 4) COUNTRY_NAME = 'CHINA'
           IF(COUNTRY == 5) COUNTRY_NAME = 'INDIA'
           IF(COUNTRY == 6) COUNTRY_NAME = 'RUSSIA'
           IPL1 = GET_YEAR_BECOMING_J(YEAR, 0, 1)
           
!          BUGO = Budget Government; Formel (61)
           BUGO = TAX_REVENUES(6, YEAR, COUNTRY) + (1.+NPOP(IPL1, FIRST_COUNTRY)) * (1.+TECH) * DEBT(IPL1, COUNTRY) - &
                & (1.+RG(YEAR)) * DEBT(YEAR, COUNTRY) - GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) - MU_1(YEAR, COUNTRY) * &
                & PENSION_BENEFITS(YEAR, COUNTRY) - MU_3(YEAR, COUNTRY) * DISABILITY_BENEFITS(YEAR, COUNTRY) - &
                & (MU_2_TAX(YEAR, COUNTRY) - MU_2_GOV(YEAR, COUNTRY)) * HEALTH_BENEFITS(YEAR, COUNTRY)

           IF(YEAR > -1) WRITE(5, 10) COUNTRY_NAME
10         FORMAT(28X, A6)

           HPOP(COUNTRY) = 0.
           DO GEN = FIRST_WORK_YEAR, GENS
                HPOP(COUNTRY) = HPOP(COUNTRY) + POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) * (1.+TECH) ** (FIRST_WORK_YEAR - GEN)
           ENDDO
           
           WRITE(5, 20) NNAME, YEAR, TRADE_BALANCE_WORLD(YEAR), FOREIGN_ASSETS_WORLD(YEAR), AGG_ASSETS_WORLD(YEAR) / &
                & YY_WORLD(YEAR), CAPITAL(YEAR, COUNTRY) / YY(YEAR, COUNTRY)          
20         FORMAT( / A4, '          TBW       QW   AGW/Y     K/Y' / I4, 5X, F8.2, 1x, 3F8.2 / )
!              
            WRITE(5, 30) YY(YEAR, COUNTRY),DD(YEAR, COUNTRY), CC(YEAR, COUNTRY), GOVERNMENT_EXPENDITURES(YEAR, COUNTRY), &
                & INVEST(YEAR, COUNTRY), GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) - MU_2_GOV(YEAR, COUNTRY) * &
                & HEALTH_BENEFITS(YEAR, COUNTRY), TRADE_BALANCE(YEAR, COUNTRY), YY(YEAR, COUNTRY) - DD(YEAR, COUNTRY) - &
                & TRADE_BALANCE(YEAR, COUNTRY) + GOV_ENDOW_SHARE(YEAR, COUNTRY)*ENDOWMENT(YEAR) 
30         FORMAT('GMA:           YY       DD    CONS       G       I    G-HB      TB YY-DD-TB' / 9X, F8.2, 1X, 6F8.2, 1x, F8.4 / )
!                
           WRITE(5, 40)AGG_ASSETS(YEAR, COUNTRY), DEBT(YEAR, COUNTRY), DEBT(YEAR, COUNTRY) / YY(YEAR, COUNTRY), &
                & AGG_ASSETS(YEAR, COUNTRY) - DEBT(YEAR, COUNTRY), CAPITAL(YEAR, COUNTRY), &
                & AGG_ASSETS(YEAR, COUNTRY) - DEBT(YEAR, COUNTRY) - CAPITAL(YEAR, COUNTRY), FOREIGN_ASSETS(YEAR, COUNTRY)
40         FORMAT('CAMA:     SAVINGS     DEBT  DEBT/Y   K-SUP   K-DEM SUP-DEM       Q' / 9X, F8.2, 1X, 12F8.2 / ) 
!                
           WRITE(5, 41)
41         FORMAT('LABMA:        LAB     LAB1   LAB2')           
!
           WRITE(5, 42)LABOR(YEAR, Y_CLASSES+1, COUNTRY), (LABOR(YEAR, Y_CLASS, COUNTRY), Y_CLASS = 1, Y_CLASSES)
42         FORMAT(9x, 3F8.2/)
!
           WRITE(5, 80) (WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY), Y_CLASS = 1, Y_CLASSES), R(YEAR, COUNTRY), RG(YEAR)
80         FORMAT('PRICES:      MPL1     MPL2      MPK     ZINS' / 9x, 4(F8.4,1x) / )
!
           WRITE(5, 50)(TAX_REVENUES(JJ, YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, JJ = 1, 6), BUGO, &
                & CONSUMP_PRICE(YEAR, COUNTRY) - 1., AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), &
                & 1.-CAPITAL_TAX_RATE(YEAR, COUNTRY), INHERITANCE_TAX_RATE(YEAR, COUNTRY), CORP_TAX(YEAR, COUNTRY), &
                & AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY)
50         FORMAT('TAXES:       TAXC     TAXY    CAPT    TAXB   CORP    TAXG   BUDGET' / 'REV:', &
                & 5X, F8.2, 1X, 5F8.2, 1X, F8.5 / 'RATE:', 1X, 'AV.', 3X, F5.3, 4X, F5.3, 3X, F5.3, 3X, F5.3, 3X, &
                & F5.3 / 6X, 'MARG.', 10X, F5.3 / )               
!                
           INC(Y_CLASSES+1) = 0.
           POPP(Y_CLASSES+1) = 0.
           DO Y_CLASS = 1, Y_CLASSES
                INC(Y_CLASS) = 0.
                POPP(Y_CLASS) = 0.
                DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY)-1
                     INC(Y_CLASS) = INC(Y_CLASS) + GET_WAGE(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY)) * POP(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & (1.+TECH) ** (FIRST_WORK_YEAR - GEN)                
                     POPP(Y_CLASS) = POPP(Y_CLASS) + POP(GEN, YEAR, Y_CLASS, COUNTRY)
                ENDDO
                INC(Y_CLASSES+1) = INC(Y_CLASSES+1) + INC(Y_CLASS)
                POPP(Y_CLASSES+1) = POPP(Y_CLASSES+1) + POPP(Y_CLASS)
                INC(Y_CLASS) = INC(Y_CLASS) / POPP(Y_CLASS)
           ENDDO
!           
           INC(Y_CLASSES+1) = INC(Y_CLASSES+1) / POPP(Y_CLASSES+1)
!
           POPR = 0.
           REPLA = 0.
           DO GEN = RETIREMENT_AGE(YEAR, COUNTRY), GENS
                POPR = POPR + POP(GEN, YEAR, Y_CLASSES+1, COUNTRY)
           ENDDO
!           
           REPLA = PENSION_BENEFITS(YEAR, COUNTRY) * POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, COUNTRY) / POPR / INC(Y_CLASSES+1)   
!           
           WRITE(5, 60) AVG_LABOR_EARNINGS(YEAR, COUNTRY), PENSION_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                & MU_1(YEAR, COUNTRY) * PENSION_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, REPLA
60         FORMAT('PENSYST:     AVLE     PB/Y  GOVANT    REPL' / 9X, F8.4, 1X, F8.2, F8.2, 3X, F5.3 / )

           WRITE(5, 70) HEALTH_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, MU_2_GOV(YEAR, COUNTRY) * &
                & HEALTH_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, DISABILITY_BENEFITS(YEAR, COUNTRY) / &
                & YY(YEAR, COUNTRY) * 100, MU_3(YEAR, COUNTRY) * DISABILITY_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                & EDUCATION_EXPENDITURES(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100                
70         FORMAT('H/D/E-SYST:  HB/Y GOVANT_H     DB/Y GOVANT_D    TRE/Y' / 9X, 5(F8.2, 1X) / )
!                
           WRITE(5, 90) (INC(Y_CLASS), AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY), &
                & AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY), AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY), &
                & AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY), AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY), &
                & AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY), AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY), &
                & AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY), Y_CLASS = 1, Y_CLASSES), INC(Y_CLASSES+1), &
                & AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), &
                & AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), &
                & AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), &
                & AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY), AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY)
90         FORMAT(12X, 'AVINC', 6X, 'TAXY', 9X, 'SST', 10X, 'HT', 11X, 'DIT' / 21X, 4('AV.', 1X, 'MARG.', 4X) / 'IC1:', 5X, &
                & F8.4, 2X, 4(F5.3, 1X, F5.3, 2X) / 'IC2:', 5X, F8.4, 2X, 4(F5.3, 1X, F5.3, 2X) / 'ALL:', 5X, F8.4, 2X, &
                & 4(F5.3, 1X, F5.3, 2X) / )                     

           DO Y_CLASS = I1, I2
                WRITE(5, 120) Y_CLASS
120             FORMAT( / 'INCOME CLASS:', 2X, I1)
                UTILITY(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY) = &
                     & GET_LIFETIME_UTILITY(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY)
                WRITE (5, 130) UTILITY(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY)
130             FORMAT( / 'UTIL:', 2X, F18.2 / 'AGE KIDS       WAGE     SHAWAG       LEIS' &
                     & '       KONS     KIKONS       ASSET        BEQ   TRANSFER      PENSP         DIFF')
!               Formel (17)
                DO GEN = FIRST_WORK_YEAR, GENS
                     IPL1 = GET_YEAR_BECOMING_J(YEAR, GEN, GEN + 1)
!                        
                     DIFF = (1 + RG(YEAR) * CAPITAL_TAX_RATE(YEAR, COUNTRY)) * ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) + &
                          & SUM_WAGE(GEN, YEAR, Y_CLASS, COUNTRY) * (AVG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) - &
                          & AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) - AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) - &
                          & AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY)) * &
                          & (HOURS - LEISURE(GEN, YEAR, Y_CLASS, COUNTRY)) + (1.-INHERITANCE_TAX_RATE(YEAR, COUNTRY) + &
                          & RG(YEAR) * CAPITAL_TAX_RATE(YEAR, COUNTRY)) * GET_INHERITANCES(GEN, YEAR, Y_CLASS, COUNTRY) + &
                          & PENSION_BENEFITS_IND(GEN, YEAR, Y_CLASS, COUNTRY) - CONSUMP(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & CONSUMP_PRICE(YEAR, COUNTRY) - KIDS(GEN, YEAR, Y_CLASS, COUNTRY) * &
                          & CONSUMP_KIDS(GEN, YEAR, Y_CLASS, COUNTRY) * CONSUMP_PRICE(YEAR, COUNTRY)
                     IF (GEN < GENS) DIFF = DIFF - ASSETS(GEN + 1, IPL1, Y_CLASS, COUNTRY)
                          
                     DIFF = DIFF + (HEALTH_BENEFITS_IND(GEN, YEAR, COUNTRY) + KIDS_HEALTH_BENEFITS(GEN, YEAR, Y_CLASS, COUNTRY)) * &
                          & (1. - MU_2_GOV(YEAR, COUNTRY)) + TRANSFER(GEN, YEAR, Y_CLASS, COUNTRY)
                          
                     IF(GEN >= FIRST_WORK_YEAR .AND. GEN <= 64) DIFF = DIFF + DISABILITY_BENEFITS_IND(YEAR, COUNTRY)
                       
                     WRITE(5, 140) GEN, KIDS(GEN, YEAR, Y_CLASS, COUNTRY), GET_WAGE(GEN, YEAR, Y_CLASS, COUNTRY), &
                          & SHADOW_WAGE(GEN, YEAR, Y_CLASS, COUNTRY), LEISURE(GEN, YEAR, Y_CLASS, COUNTRY), &
                          & CONSUMP(GEN, YEAR, Y_CLASS, COUNTRY), CONSUMP_KIDS(GEN, YEAR, Y_CLASS, COUNTRY), &
                          & ASSETS(GEN, YEAR, Y_CLASS, COUNTRY), GET_INHERITANCES(GEN, YEAR, Y_CLASS, COUNTRY), &
                          & TRANSFER(GEN, YEAR, Y_CLASS, COUNTRY), PENSION_BENEFITS_IND(GEN, YEAR, Y_CLASS, COUNTRY), DIFF
140                  FORMAT(I3, 1X, F4.2, 1X, 5(F10.6, 1X), F11.6, 1X, 3(F10.6, 1X), F12.10)
                          
                ENDDO

           ENDDO                          
!
           WRITE(5, 150)
150        FORMAT(' ---------------------------------------------------------------------------------------------- ')

           RETURN
      END SUBROUTINE OUTPUT
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE PRINT_INITIAL_EQUILIBRIUM                                                                                   *
! *   PURPOSE: Write data of the initial equilibrium                                                                               *
! **********************************************************************************************************************************

      SUBROUTINE PRINT_INITIAL_EQUILIBRIUM 
           USE GLOBAL_DATA
           IMPLICIT NONE
           CHARACTER (LEN = 6), DIMENSION(FIRST_COUNTRY:LAST_COUNTRY) :: COUNTRY_NAME
           INTEGER :: COUNTRY
!
           IF(IRUN == 0)OPEN (UNIT = 9, FILE = OUTFILE_NAME(5), STATUS = 'REPLACE')
           IF(IRUN /= 0)OPEN (UNIT = 9, FILE = OUTFILE_NAME(5), STATUS = 'OLD', POSITION='APPEND')  
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                IF(COUNTRY == 1) COUNTRY_NAME(COUNTRY) = '   USA'
                IF(COUNTRY == 2) COUNTRY_NAME(COUNTRY) = 'EUROPE'
                IF(COUNTRY == 3) COUNTRY_NAME(COUNTRY) = ' JAPAN'
                IF(COUNTRY == 4) COUNTRY_NAME(COUNTRY) = ' CHINA'
                IF(COUNTRY == 5) COUNTRY_NAME(COUNTRY) = ' INDIA'
                IF(COUNTRY == 6) COUNTRY_NAME(COUNTRY) = ' RUSSIA'
           ENDDO             
           WRITE(9, 10)(COUNTRY_NAME(COUNTRY), COUNTRY = FIRST_COUNTRY, LAST_COUNTRY), &
                & (POP(91, 0, Y_CLASSES+1, COUNTRY) / 1000000, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
10         FORMAT(37x, 6(A6, 2x) / 'Population (in mio.)', 19x, 6(F5.0, 3x) // 'National Income')           
           WRITE(9, 20)(CC(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
20         FORMAT('Private consumption', 20x, 6(F4.1, 4x))
           WRITE(9, 30)(GOVERNMENT_EXPENDITURES(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
30         FORMAT('Government purchases', 19x, 6(F4.1, 4x))           
           WRITE(9, 40)((GOVERNMENT_EXPENDITURES(0, COUNTRY) - EDUCATION_EXPENDITURES(0, COUNTRY) - MU_2_GOV(0, COUNTRY) * &
                & HEALTH_BENEFITS(0, COUNTRY)) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
40         FORMAT(5x, 'General public expenditures', 7x, 6(F4.1, 4x))           
           WRITE(9, 50)(EDUCATION_EXPENDITURES(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
50         FORMAT(5x, 'Aggregate education outlays', 7x, 6(F4.1, 4x))
           WRITE(9, 60)(MU_2_GOV(0, COUNTRY) * HEALTH_BENEFITS(0, COUNTRY) / YY(0, COUNTRY) * 100, &
                & COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
60         FORMAT(5x, 'Aggregate health benefits', 9x, 6(F4.1, 4x))
           WRITE(9, 70)(INVEST(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
70         FORMAT('Investment', 29x, 6(F4.1, 4x))           
           WRITE(9, 80)(TRADE_BALANCE(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
80         FORMAT('Trade balance', 25x, 6(F5.1, 3x))
           WRITE(9, 120)((TRADE_BALANCE(0, COUNTRY) + RG(0) * FOREIGN_ASSETS(0, COUNTRY)) / YY(0, COUNTRY) * 100, &
                & COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
120        FORMAT('Current account', 23x, 6(F5.1, 3x))
           WRITE(9, 130)
130        FORMAT(/ 'Government indicators')
           WRITE(9, 140)(((1. - MU_1(0, COUNTRY)) * PENSION_BENEFITS(0, COUNTRY) + (1. - MU_2_TAX(0, COUNTRY)) * &
                & HEALTH_BENEFITS(0, COUNTRY) + (1. - MU_3(0, COUNTRY)) * DISABILITY_BENEFITS(0, COUNTRY)) / &
                & YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
140        FORMAT('Social contributions received', 10x, 6(F4.1, 4x))
           WRITE(9, 141)((1. - MU_1(0, COUNTRY)) * PENSION_BENEFITS(0, COUNTRY) / YY(0, COUNTRY) * &
                & 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
141        FORMAT(5x, 'Pension revenues', 18x, 6(F4.1, 4x))
           WRITE(9, 142)((1. - MU_2_TAX(0, COUNTRY)) * HEALTH_BENEFITS(0, COUNTRY) / YY(0, COUNTRY) * &
                & 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
142        FORMAT(5x, 'Health revenues', 19x, 6(F4.1, 4x))
           WRITE(9, 143)((1. - MU_3(0, COUNTRY)) * DISABILITY_BENEFITS(0, COUNTRY) / YY(0, COUNTRY) * &
                & 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
143        FORMAT(5x, 'Disability revenues', 15x, 6(F4.1, 4x))
           WRITE(9, 149)((PENSION_BENEFITS(0, COUNTRY) + HEALTH_BENEFITS(0, COUNTRY) + &
                & DISABILITY_BENEFITS(0, COUNTRY)) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
149        FORMAT('Social benefits', 24x, 6(F4.1, 4x))
           WRITE(9, 150)(PENSION_BENEFITS(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
150        FORMAT(5x, 'Aggregate pension benefits', 8x, 6(F4.1, 4x))
           WRITE(9, 160)(HEALTH_BENEFITS(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
160        FORMAT(5x, 'Aggregate health benefits', 9x, 6(F4.1, 4x))
           WRITE(9, 170)(DISABILITY_BENEFITS(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
170        FORMAT(5x, 'Aggregate disability benefits', 5x, 6(F4.1, 4x))
           WRITE(9, 180)(AGG_PENSION_TAX_RATE(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
180        FORMAT('Pension contribution rate', 14x, 6(F4.1, 4x))
           WRITE(9, 190)(AGG_HEALTH_TAX_RATE(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
190        FORMAT('Health-care contribution rate', 10x, 6(F4.1, 4x))
           WRITE(9, 200)(AGG_DISABILITY_TAX_RATE(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
200        FORMAT('Disability-insurance contribution rate', 1x, 6(F4.1, 4x))
           WRITE(9, 210)(RG(0) * DEBT(0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
210        FORMAT('Interest payment on public debt', 8x, 6(F4.1, 4x))
           WRITE(9, 220)(TAX_REVENUES(6, 0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
220        FORMAT(/ 'Tax revenues', 27x, 6(F4.1, 4x))
           WRITE(9, 230)((TAX_REVENUES(2, 0, COUNTRY) + TAX_REVENUES(3, 0, COUNTRY) + &
                & TAX_REVENUES(5, 0, COUNTRY)) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
230        FORMAT(5x, 'Direct taxes', 22x, 6(F4.1, 4x))
           WRITE(9, 240)(TAX_REVENUES(2, 0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
240        FORMAT(10x, 'Wage taxes', 19x, 6(F4.1, 4x))
           WRITE(9, 250)(TAX_REVENUES(3, 0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
250        FORMAT(10x, 'Capital taxes', 16x, 6(F4.1, 4x))
           WRITE(9, 251)(TAX_REVENUES(5, 0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
251        FORMAT(10x, 'Corporate taxes', 14x, 6(F4.1, 4x))
           WRITE(9, 260)((TAX_REVENUES(1, 0, COUNTRY) + TAX_REVENUES(4, 0, COUNTRY)) / YY(0, COUNTRY) * 100, &
                & COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
260        FORMAT(5x, 'Indirect taxes', 20x, 6(F4.1, 4x))
           WRITE(9, 261)(TAX_REVENUES(1, 0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
261        FORMAT(10x, 'Consumption taxes', 12x, 6(F4.1, 4x))
           WRITE(9, 262)(TAX_REVENUES(4, 0, COUNTRY) / YY(0, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
262        FORMAT(10x, 'Inheritance taxes', 12x, 6(F4.1, 4x))
           WRITE(9, 270)
270        FORMAT(/ 'Wage tax rates', 25x, 6(F4.1, 4x))
           WRITE(9, 280)(AGG_AVG_WAGE_TAX_RATE(0, Y_CLASSES+1, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
280        FORMAT(5x, 'Average', 27x, 6(F4.1, 4x))
           WRITE(9, 290)(AGG_MARG_WAGE_TAX_RATE(0, Y_CLASSES+1, COUNTRY) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
290        FORMAT(5x, 'Marginal', 26x, 6(F4.1, 4x))
           WRITE(9, 300)(CAPITAL(0, COUNTRY) / YY(0, COUNTRY), COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
300        FORMAT('Capital-output ratio', 19x, 6(F4.1, 4x))
           WRITE(9, 310)(CAPITAL(0, COUNTRY) / LABOR(0, Y_CLASSES+1, COUNTRY), COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
310        FORMAT('Capital-labor ratio', 20x, 6(F4.1, 4x))
           WRITE(9, 320)(RG(0) * 100, COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
320        FORMAT('Interest rate', 26x, 6(F4.1, 4x))
           WRITE(9, 330)(WAGE_INDEX_CLASS(0, 1, COUNTRY), COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
330        FORMAT('Wage rate' / 5x, 'Income class 1', 20x, 6(F4.1, 4x))
           WRITE(9, 340)(WAGE_INDEX_CLASS(0, 2, COUNTRY), COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
340        FORMAT(5x, 'Income class 2', 20x, 6(F4.1, 4x))
           WRITE(9, 360)(YY(0, COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY) / &
                & (YY(0, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), &
                & COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
360        FORMAT('GDP relative', 27x, 6(F4.2, 4x))
           WRITE(9, 380)(((YY(0+1, COUNTRY) * POP(FIRST_WORK_YEAR, 0+1, Y_CLASSES+1, FIRST_COUNTRY) / &
                & (YY(0, COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY))) - 1)*100, &
                & COUNTRY = FIRST_COUNTRY, LAST_COUNTRY)
380        FORMAT('GDP growth rate', 24x, 6(F4.1, 4x))
!
           CLOSE(9)
!           
           RETURN
      END SUBROUTINE PRINT_INITIAL_EQUILIBRIUM
      
! ********************************************************************************************************************************** 
! *    NAME: SUBROUTINE RESULTS                                                                                                    *
! *    PURPOSE:                                                                                                                    *
! **********************************************************************************************************************************

      SUBROUTINE RESULTS
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, JJ, Y_CLASS, YEAR, YEAR1, GEN
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           CHARACTER (LEN = 10), DIMENSION(FIRST_COUNTRY:LAST_COUNTRY) :: COUNTRY_NAME
           REAL*8, DIMENSION(0:YEARS, FIRST_COUNTRY:LAST_COUNTRY) :: BEQ
           REAL*8, DIMENSION(0:YEARS, 1:Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) :: AVERAGE_WAGE_RATE_CLASS, EFFICIENT_WORK_POP
           REAL*8, EXTERNAL :: GET_EFFICIENT_POPULATION, GET_INHERITANCES
!
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
 	            IF(IRUN == 0) THEN
                     IF(COUNTRY == 1) THEN
                          COUNTRY_NAME(COUNTRY) = 'USA'
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(13), STATUS = 'REPLACE')
                     ELSEIF(COUNTRY == 2) THEN
                          COUNTRY_NAME(COUNTRY) = 'EUROPE'
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(14), STATUS = 'REPLACE')
                     ELSEIF(COUNTRY == 3) THEN 
                          COUNTRY_NAME(COUNTRY) = 'JAPAN'
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(15), STATUS = 'REPLACE')
                     ELSEIF(COUNTRY == 4) THEN
                          COUNTRY_NAME(COUNTRY) = 'CHINA'
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(16), STATUS = 'REPLACE')
                     ELSEIF(COUNTRY == 5) THEN
                          COUNTRY_NAME(COUNTRY) = 'INDIA'
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(17), STATUS = 'REPLACE')
                     ELSEIF(COUNTRY == 6) THEN
                          COUNTRY_NAME(COUNTRY) = 'RUSSIA'
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(18), STATUS = 'REPLACE')
                     ENDIF     
!
                     WRITE(9, 10) COUNTRY_NAME(COUNTRY)
10                   FORMAT(A6 / 291('-') / 'BASELINE PATH' /)
                     WRITE(9, 20)
20                   FORMAT('Macroeconomic variables:', 56x, 'Prices:', 22x, 'Tax rates:', 54x, 'Tax revenues:(/Y)', 33x, &
                          & 'Government indicators:' / 84('-'), 4x, 25('-'), 4x, 60('-'), 4x, 46('-'), 4x, 60('-') / &
                          & 'YEAR',  6x, 'K', 7x, 'L', 6x, 'L1', 7x, 'L2', 7x, 'GDP', 6x, 'GNI', 6x, 'GDP',&
                          & 5x, 'C/NI', 4x, &
                          & 'CA/GDP', 3x, 'TB/GDP', 4x, 'TRF/GDP', 1x, 'I/GDP', 3x, 'BEQ/GDP', 4x, 'YEAR', 2x,&
                          & 'RG', 6x, 'W1', 6x, 'W2', &
                          & 8x, 'YEAR', 4x, 'OAS', 5x, 'HI', 5x, 'DI', 5x, 'CONS', 4x, 'WAGE', 4x, 'CAP', 5x, 'INH', 4x, 'CORP', &
                          & 7x, 'YEAR', 1x, 'CONS/GDP', 1x, 'WAGE/GDP', 2x, 'CAP/GDP', 2x, 'INH/GDP', 1x, 'CORP/GDP',&
                          & 1x, 'OIL/GDP', 2X, 'SUM/GDP', 1X, 'PENS/GDP', 1X,  'SUM2/GDP', 5x, 'YEAR', 2x, &
                          & 'PB/GDP', 3x, 'HB/GDP', 3x, 'DB/GDP', 3x, 'EDU/GDP',  1x, 'OTH/GDP', 3x, 'INT/GDP',&
                          & 1x, 'B/GDP', 5x, 'DEF/GDP',&
                          & 1x , 'TOTAL/GDP')
!
                     DO YEAR = 0, 100
!  
!                              Aggregate bequests
                               BEQ(YEAR, COUNTRY) = 0.
                               DO Y_CLASS = 1, Y_CLASSES
                                    DO GEN = FIRST_WORK_YEAR, GENS
                                         BEQ(YEAR, COUNTRY) = BEQ(YEAR, COUNTRY) + &
                                              & GET_INHERITANCES(GEN, YEAR, Y_CLASS, COUNTRY) * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                                    ENDDO          
                               ENDDO           
!                              Deficit   
                               YEAR1 = GET_YEAR_BECOMING_J(YEAR, 1, 2)                       
                               DEFICIT(YEAR, COUNTRY) = (1.+TECH) * (1. + NPOP(YEAR1, FIRST_COUNTRY)) * &
                                    & DEBT(YEAR1, COUNTRY) - DEBT(YEAR, COUNTRY)
!                              Efficient working population of income class
                               DO Y_CLASS = 1, Y_CLASSES
                                    EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY) = 0.
                                    AVERAGE_WAGE_RATE_CLASS(YEAR, Y_CLASS, COUNTRY) = 0.
                                    DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) - 1
                                         EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY) = EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY) + &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)                                              
                                    ENDDO  
!                              Average wage rate of income class
                                    DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) - 1
                                         AVERAGE_WAGE_RATE_CLASS(YEAR, Y_CLASS, COUNTRY) = &
                                              & AVERAGE_WAGE_RATE_CLASS(YEAR, Y_CLASS, COUNTRY) + &
                                              & WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY) * &
                                              & AGE_EFFICIENCY(GEN, YEAR, Y_CLASS, COUNTRY) * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY) / &
                                              & EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY)
                                    ENDDO  
                               ENDDO                                      
!
                               WRITE(9, 3001) FIRST_YEAR + YEAR, CAPITAL(YEAR, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (CAPITAL(0, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)),&
                                    & LABOR(YEAR, Y_CLASSES+1, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (LABOR(0, Y_CLASSES+1, FIRST_COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), (LABOR(YEAR, Y_CLASS, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (LABOR(0, Y_CLASS, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), &
                                    & Y_CLASS = 1, Y_CLASSES), YY(YEAR, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (YY(0, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), &
                                    
                                    &(NATIONAL_INCOME(YEAR, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (NATIONAL_INCOME(0, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY))),&

                                    &(YY(YEAR, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (YY(0, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY))),&
                                    
                                    & CC(YEAR, COUNTRY) / NATIONAL_INCOME(YEAR, COUNTRY) * 100,&
                                    & (TRADE_BALANCE(YEAR, COUNTRY) + RG(YEAR) * &
                                    & FOREIGN_ASSETS(YEAR, COUNTRY)) / YY(YEAR, COUNTRY) * 100,&
                                    
                                    &TRADE_BALANCE(YEAR, COUNTRY)/ YY(YEAR, COUNTRY) * 100,&
                                    & TRF(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & INVEST(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & BEQ(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100,&




                                    
                                    & FIRST_YEAR + YEAR, RG(YEAR) * 100, (WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY), &
                                    & Y_CLASS = 1, Y_CLASSES), &




                                    
                                    & FIRST_YEAR + YEAR, AGG_PENSION_TAX_RATE(YEAR, COUNTRY) * 100, &
                                    & AGG_HEALTH_TAX_RATE(YEAR, COUNTRY) * 100, AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY) * 100, &
                                    & (CONSUMP_PRICE(YEAR, COUNTRY) - 1.) * 100, &
                                    & AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) * 100, &
                                    & (1. - CAPITAL_TAX_RATE(YEAR, COUNTRY)) * 100, INHERITANCE_TAX_RATE(YEAR, COUNTRY) * 100, &
                                    & CORP_TAX(YEAR, COUNTRY) * 100,&




                                    
                                    & FIRST_YEAR + YEAR, (TAX_REVENUES(JJ, YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, JJ = 1, 5), &
                                    & (GOV_ENDOW_SHARE(YEAR, COUNTRY)*ENDOWMENT(YEAR)*100)&
                                    & /YY(YEAR, COUNTRY),& 
                                    & (TAX_REVENUES(6, YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100), &
                                    
                                    &((1-MU_1(YEAR, COUNTRY))*PENSION_BENEFITS(YEAR, COUNTRY))&
                                    &/ YY(YEAR, COUNTRY) * 100,&

                                    &((TAX_REVENUES(6, YEAR, COUNTRY)+&
                                    & (1-MU_1(YEAR, COUNTRY))*PENSION_BENEFITS(YEAR, COUNTRY))&
                                    &/ YY(YEAR, COUNTRY) * 100), &




                                    
                                    & FIRST_YEAR + YEAR, PENSION_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & HEALTH_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & DISABILITY_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & EDUCATION_EXPENDITURES(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & GOVERNMENT_DISCRETIONARY_SPENDING(YEAR, COUNTRY)/NATIONAL_INCOME(YEAR, COUNTRY) * 100,&
                                    & RG(YEAR) * DEBT(YEAR, COUNTRY)/YY(YEAR, COUNTRY) * 100,&
                                    & DEBT(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, DEFICIT(YEAR, COUNTRY) / &
                                    & YY(YEAR, COUNTRY) * 100, & 
                               &(PENSION_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100) &
                               &+(HEALTH_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100)&
                               &+(DISABILITY_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100)&
                               &+(EDUCATION_EXPENDITURES(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100)&
                               &+(GOVERNMENT_DISCRETIONARY_SPENDING(YEAR, COUNTRY)/YY(YEAR, COUNTRY) * 100)&
                               &+(RG(YEAR) * DEBT(YEAR, COUNTRY)/YY(YEAR, COUNTRY) * 100)
                                
                               
                               
                               
                                
                               
! 
3001                             FORMAT(I4, 1X, 8(F7.3,2x), 3(F7.3,2x), 2(F6.3,2x), 3x, I4, 1x, 3(F6.3, 2x), 3x, I4, 1x, &
                                    & 8(F6.3, 2x), 3x, I4, 1x, 9(F6.3, 3x), 3x, I4, 1x, 8(F6.3, 3x))

!$$$$$$ 3001      FORMAT(I4, 1X, 8(F7.3,2x), 3(F7.3,2x), 2(F6.3,2x), 3x, I4, 1x, 3(F6.3, 2x), 3x, I4, 1x, &
!$$$$$$                                      & 8(F6.3, 2x))      
                                    
                     ENDDO     
                     WRITE(9,40)
40                   FORMAT(84('-'), 4x, 25('-'), 4x, 60('-'), 4x, 46('-'), 4x, 60('-'))    
					 CLOSE(9)
!                                                             
!                     
                ELSEIF(IRUN == 1) THEN
                     IF(COUNTRY == 1) THEN
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(13), STATUS = 'OLD', POSITION='APPEND')
                     ELSEIF(COUNTRY == 2) THEN
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(14), STATUS = 'OLD', POSITION='APPEND')
                     ELSEIF(COUNTRY == 3) THEN 
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(15), STATUS = 'OLD', POSITION='APPEND')
                     ELSEIF(COUNTRY == 4) THEN
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(16), STATUS = 'OLD', POSITION='APPEND')
                     ELSEIF(COUNTRY == 5) THEN
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(17), STATUS = 'OLD', POSITION='APPEND')
                     ELSEIF(COUNTRY == 6) THEN
                          OPEN (UNIT = 9, FILE = OUTFILE_NAME(18), STATUS = 'OLD', POSITION='APPEND')
                     ENDIF
!
!                    Write variables of policy reform as indexed values
                     WRITE(9, 50)
50                   FORMAT(/ 'POLICY REFORM (indexed values)' /)
					 WRITE(9,20)
                     DO YEAR = 0, 100
!  
!                              Aggregate bequests
                               BEQ(YEAR, COUNTRY) = 0.
                               DO Y_CLASS = 1, Y_CLASSES
                                    DO GEN = FIRST_WORK_YEAR, GENS
                                         BEQ(YEAR, COUNTRY) = BEQ(YEAR, COUNTRY) + &
                                              & GET_INHERITANCES(GEN, YEAR, Y_CLASS, COUNTRY) * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)
                                    ENDDO          
                               ENDDO           
!                              Deficit   
                               YEAR1 = GET_YEAR_BECOMING_J(YEAR, 1, 2)                       
                               DEFICIT(YEAR, COUNTRY) = (1.+TECH) * (1. + NPOP(YEAR1, FIRST_COUNTRY)) * &
                                    & DEBT(YEAR1, COUNTRY) - DEBT(YEAR, COUNTRY)
!                              Efficient working population of income class
                               DO Y_CLASS = 1, Y_CLASSES
                                    EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY) = 0.
                                    AVERAGE_WAGE_RATE_CLASS(YEAR, Y_CLASS, COUNTRY) = 0.
                                    DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) - 1
                                         EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY) = EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY) + &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY)                                              
                                    ENDDO  
!                              Average wage rate of income class
                                    DO GEN = FIRST_WORK_YEAR, RETIREMENT_AGE(YEAR, COUNTRY) - 1
                                         AVERAGE_WAGE_RATE_CLASS(YEAR, Y_CLASS, COUNTRY) = &
                                              & AVERAGE_WAGE_RATE_CLASS(YEAR, Y_CLASS, COUNTRY) + &
                                              & WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY) * &
                                              & AGE_EFFICIENCY(GEN, YEAR, Y_CLASS, COUNTRY) * &
                                              & GET_EFFICIENT_POPULATION(GEN, YEAR, Y_CLASS, COUNTRY) / &
                                              & EFFICIENT_WORK_POP(YEAR, Y_CLASS, COUNTRY)
                                    ENDDO  
                               ENDDO                                      
!
                               WRITE(9, 3001) FIRST_YEAR + YEAR, CAPITAL(YEAR, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (CAPITAL(0, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), &
                                    & LABOR(YEAR, Y_CLASSES+1, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (LABOR(0, Y_CLASSES+1, FIRST_COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), (LABOR(YEAR, Y_CLASS, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (LABOR(0, Y_CLASS, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), &
                                    & Y_CLASS = 1, Y_CLASSES), YY(YEAR, COUNTRY) * &
                                    & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) * (1.+TECH) ** YEAR / &
                                    & (YY(0, FIRST_COUNTRY) * POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY)), &
                                    & CC(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, (TRADE_BALANCE(YEAR, COUNTRY) + RG(YEAR) * &
                                    & FOREIGN_ASSETS(YEAR, COUNTRY)) / YY(YEAR, COUNTRY) * 100, TRADE_BALANCE(YEAR, COUNTRY) / &
                                    & YY(YEAR, COUNTRY) * 100, TRF(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & INVEST(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & BEQ(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & FIRST_YEAR + YEAR, RG(YEAR) * 100, (WAGE_INDEX_CLASS(YEAR, Y_CLASS, COUNTRY), &
                                    & Y_CLASS = 1, Y_CLASSES), &
                                    & FIRST_YEAR + YEAR, AGG_PENSION_TAX_RATE(YEAR, COUNTRY) * 100, &
                                    & AGG_HEALTH_TAX_RATE(YEAR, COUNTRY) * 100, AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY) * 100, &
                                    & (CONSUMP_PRICE(YEAR, COUNTRY) - 1.) * 100, &
                                    & AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASSES+1, COUNTRY) * 100, &
                                    & (1. - CAPITAL_TAX_RATE(YEAR, COUNTRY)) * 100, INHERITANCE_TAX_RATE(YEAR, COUNTRY) * 100, &
                                    & CORP_TAX(YEAR, COUNTRY) * 100, &
                                    & FIRST_YEAR + YEAR, (TAX_REVENUES(JJ, YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, JJ = 1, 6), &
                                    & FIRST_YEAR + YEAR, PENSION_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & HEALTH_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & DISABILITY_BENEFITS(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & EDUCATION_EXPENDITURES(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & (GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) - EDUCATION_EXPENDITURES(YEAR, COUNTRY) &
                                    &- DISABILITY_BENEFITS(YEAR, COUNTRY) - HEALTH_BENEFITS(YEAR, COUNTRY) &
                                    & - PENSION_BENEFITS(YEAR, COUNTRY) - RG(YEAR) * DEBT(YEAR, COUNTRY)) / &
                                    & YY(YEAR, COUNTRY) * 100, GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, &
                                    & DEBT(YEAR, COUNTRY) / YY(YEAR, COUNTRY) * 100, DEFICIT(YEAR, COUNTRY) / &
                                    & YY(YEAR, COUNTRY) * 100, GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) + & 
                               & GOVS(COUNTRY) * (YY_0(YEAR, COUNTRY) / YY(0, COUNTRY) * &
                               & POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) / &
                               & POP(FIRST_WORK_YEAR, 0, Y_CLASSES+1, FIRST_COUNTRY) * POP_EFFICIENT(0, COUNTRY) / &
                               & POP_EFFICIENT(YEAR, COUNTRY)) * &
                               & POP(GEN, YEAR, Y_CLASSES+1, COUNTRY) / POP(FIRST_WORK_YEAR, YEAR, Y_CLASSES+1, FIRST_COUNTRY) 
                     ENDDO 
                     WRITE(9,40)
                     WRITE(9,80)     
80                   FORMAT(291('-'))                     
                     CLOSE(9)              
                ENDIF     
           ENDDO

			   		OPEN (UNIT = 50, FILE = OUTFILE_NAME(19), STATUS = 'REPLACE')

                    WRITE(50, 200)
200                   FORMAT('YEAR',  6x, 'Total Endowment Value/World Assets',&
						& 6x, 'Total Government Oil Revenue/World GDP',2x,&
                          & 'Total Household Oil Revenue/World GDP')
					
					DO YEAR = 0, 100
                    WRITE(50, 300) FIRST_YEAR + YEAR, PVENDOWMENT(YEAR)/AGG_ASSETS_WORLD(YEAR)&
                    &, (TOT_GOV_ENDOW_SHARE(YEAR)*ENDOWMENT(YEAR))/YY_WORLD(YEAR),&
                    & ((1-TOT_GOV_ENDOW_SHARE(YEAR))*ENDOWMENT(YEAR))/YY_WORLD(YEAR)
300                   FORMAT(I4, 6X, 3(F7.6,35x))    
                    ENDDO
					 WRITE(50,400)
400                   FORMAT(84('-'), 4x, 25('-'), 4x, 60('-'), 4x, 46('-'), 4x, 60('-'))    
					 CLOSE(50)

      
           RETURN
      END SUBROUTINE RESULTS
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE PRINT_HICKSIAN_EV                                                                                           *
! *   PURPOSE: Bestimmung der Hick'schen äquivalenten Variation (HICKS_EQUIVALENT_VARIATION)                                       *
! *        HICKS_EQUIVALENT_VARIATION fragt nach der Einkommensänderung, die man einem Individuum bei den ursprünglichen Preisen   *
! *        geben müßte, damit es genauso gestellt wird wie nach einer ökonomischen Änderung, d.h. bei den neuen Preisen            *
! **********************************************************************************************************************************

      SUBROUTINE PRINT_HICKSIAN_EV
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, Y_CLASS, YEAR
           REAL*8, DIMENSION(FIRST_WORK_YEAR:GENS, 0:YEARS, Y_CLASSES, FIRST_COUNTRY:LAST_COUNTRY) ::  PHI
           REAL*8, EXTERNAL :: GET_LIFETIME_UTILITY
!
           DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                IF(COUNTRY == 1) THEN
                     OPEN (UNIT = 9, FILE = OUTFILE_NAME(13), STATUS = 'OLD', POSITION='APPEND')
                ELSEIF(COUNTRY == 2) THEN
                     OPEN (UNIT = 9, FILE = OUTFILE_NAME(14), STATUS = 'OLD', POSITION='APPEND')
                ELSEIF(COUNTRY == 3) THEN 
                     OPEN (UNIT = 9, FILE = OUTFILE_NAME(15), STATUS = 'OLD', POSITION='APPEND')
                ELSEIF(COUNTRY == 4) THEN
                     OPEN (UNIT = 9, FILE = OUTFILE_NAME(16), STATUS = 'OLD', POSITION='APPEND')
                ELSEIF(COUNTRY == 5) THEN
                     OPEN (UNIT = 9, FILE = OUTFILE_NAME(17), STATUS = 'OLD', POSITION='APPEND')
                ELSEIF(COUNTRY == 6) THEN
                     OPEN (UNIT = 9, FILE = OUTFILE_NAME(18), STATUS = 'OLD', POSITION='APPEND')
                ENDIF
!                
!               Welfare measure as in Auerbach/Kotlikoff (1987)
                DO Y_CLASS = 1, Y_CLASSES
                     DO GEN = FIRST_WORK_YEAR, GENS
                          PHI(GEN, 0, Y_CLASS, COUNTRY) = ((GET_LIFETIME_UTILITY(GEN, 0, Y_CLASS, COUNTRY) / &
                               & UTILITY_0(GEN, 0, Y_CLASS, COUNTRY)) ** (GAMMA / (GAMMA -1.)) -1.) 
                     ENDDO
                     DO YEAR = 1, YEARS
                          PHI(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY) = &
                               & ((GET_LIFETIME_UTILITY(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY) / &
                               & UTILITY_0(FIRST_WORK_YEAR, YEAR, Y_CLASS, COUNTRY)) ** (GAMMA / (GAMMA-1.)) -1.)
                     ENDDO
                ENDDO
!
                WRITE(9, 50)
                WRITE(9, 20)
20              FORMAT(/ 'WELFARE EFFECTS - POLICY REFORM', 12x / 'BIRTH YEAR', 6X, '1', 7X, '2')
                DO GEN = FIRST_YEAR - 1920, FIRST_YEAR - 1980, - 5
                     WRITE(9, 40) FIRST_YEAR - GEN, (PHI(GEN, 0, Y_CLASS, COUNTRY)*100, &
                          & Y_CLASS = 1, Y_CLASSES)
40                   FORMAT(I4, 6X, 2(F7.3, 1X))
                ENDDO
                DO GEN = 1990 - FIRST_YEAR, 2100 - FIRST_YEAR, 5
                     WRITE(9, 40) GEN + FIRST_YEAR, (PHI(FIRST_WORK_YEAR, FIRST_WORK_YEAR + GEN, Y_CLASS, COUNTRY)*100, &
                          & Y_CLASS = 1, Y_CLASSES)
                ENDDO
                WRITE(9, 50)
50              FORMAT(31('-'), 10x)
                CLOSE(9)
           ENDDO             
           
           RETURN
      END SUBROUTINE PRINT_HICKSIAN_EV
      
! ********************************************************************************************************************************** 
! *   NAME: FUNCTION GET_LIFETIME_UTILITY(GEN, YEAR, Y_CLASS, COUNTRY)                                                             *
! *   PURPOSE: Compute lifetime utility                                                                                            *
! **********************************************************************************************************************************

      REAL*8 FUNCTION GET_LIFETIME_UTILITY(GEN, YEAR, Y_CLASS, COUNTRY)
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: IK, J
           INTEGER, INTENT(IN) :: COUNTRY, GEN, Y_CLASS, YEAR
           INTEGER, EXTERNAL :: GET_YEAR_BECOMING_J
           REAL*8, EXTERNAL :: KIDS
           
           GET_LIFETIME_UTILITY = 0.
           
!          Formel (13)
           DO J = GEN, GENS
                IK = GET_YEAR_BECOMING_J(YEAR, GEN, J)
                GET_LIFETIME_UTILITY = GET_LIFETIME_UTILITY + (1./(1.+DELTA(GEN, YEAR, COUNTRY))) ** (J - GEN) * &
                     & SURVIVAL_PROBABILITY(J, IK, COUNTRY) * ((CONSUMP(J, IK, Y_CLASS, COUNTRY) ** (1.-(1./RHO)) + &
                     & ALP * LEISURE(J, IK, Y_CLASS, COUNTRY) ** (1.- (1./RHO))) ** ((1.-(1./GAMMA)) / (1.-(1./RHO))))
                     
                IF(J >= 23 .AND. J<= 65) GET_LIFETIME_UTILITY = GET_LIFETIME_UTILITY + (1./(1.+DELTA(GEN, YEAR, COUNTRY))) ** &
                     & (J - GEN) * SURVIVAL_PROBABILITY(J, IK, COUNTRY) * (THETA * KIDS(J, IK, Y_CLASS, COUNTRY) * &
                     & CONSUMP_KIDS(J, IK, Y_CLASS, COUNTRY) ** (1.- (1./GAMMA)))         
           ENDDO
           
           GET_LIFETIME_UTILITY = GET_LIFETIME_UTILITY / (1.-(1./GAMMA))
           
           RETURN
      END FUNCTION GET_LIFETIME_UTILITY
      
! ********************************************************************************************************************************** 
! *   NAME: SUBROUTINE INITIALIZE_VARIABLES                                                                                        *
! *   PURPOSE:                                                                                                                     *
! **********************************************************************************************************************************

      SUBROUTINE INITIALIZE_VARIABLES
           USE GLOBAL_DATA
           IMPLICIT NONE
           INTEGER :: COUNTRY, GEN, J, Y_CLASS, YEAR, H_COUNTRY
           
           DO YEAR = -1, YEARS
                AGG_ASSETS_WORLD(YEAR) = 0.
                TRADE_BALANCE_WORLD(YEAR) = 0.
                RG(YEAR) = 0.1
                AP(YEAR) = 1.
                DO COUNTRY = FIRST_COUNTRY, LAST_COUNTRY                     
                     R(YEAR, COUNTRY) = 0.1
                     AGG_ASSETS(YEAR, COUNTRY) = 1.
                     INVEST(YEAR, COUNTRY) = 0.
                     YY(YEAR, COUNTRY) = 1.
                     YY_0(YEAR, COUNTRY) = 1.
                     DD(YEAR, COUNTRY) = 1.
                     CAPITAL_TAX_RATE(YEAR, COUNTRY) = 1.
                     DEBT(YEAR, COUNTRY) = 0.
                     GOVERNMENT_EXPENDITURES(YEAR, COUNTRY) = 0.
                     EDUCATION_EXPENDITURES(YEAR, COUNTRY) = 0.
                     DISABILITY_BENEFITS(YEAR, COUNTRY) = 0.
                     HEALTH_BENEFITS(YEAR, COUNTRY) = 0.
                     PENSION_BENEFITS(YEAR, COUNTRY) = 0.
                     TRADE_BALANCE(YEAR, COUNTRY) = 0.
                     FOREIGN_ASSETS(YEAR, COUNTRY) = 0.
                     AVG_LABOR_EARNINGS(YEAR, COUNTRY) = 0.
                     CAPITAL(YEAR, COUNTRY) = 100.
                     CONSUMP_PRICE(YEAR, COUNTRY) = 1.
                     AGG_DISABILITY_TAX_RATE(YEAR, COUNTRY) = 0.
                     TRF(YEAR, COUNTRY) = 0.
                     DO Y_CLASS = 1, Y_CLASSES
                          LABOR(YEAR, Y_CLASS, COUNTRY) = 10.
                          AGG_ASSETS_FOR_BEQUESTS(YEAR, Y_CLASS, COUNTRY) = 0.
                     ENDDO     
                     DO Y_CLASS = 1, Y_CLASSES+1
                          AGG_AVG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_MARG_WAGE_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_AVG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_MARG_PENSION_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_AVG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_MARG_HEALTH_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_AVG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                          AGG_MARG_DISABILITY_TAX_RATE(YEAR, Y_CLASS, COUNTRY) = 0.
                     ENDDO     
                     CC(YEAR, COUNTRY) = 1.                
                     DO J = 1, 6
                          TAX_REVENUES(J, YEAR, COUNTRY) = 0.
                     ENDDO
                     DO GEN = FIRST_WORK_YEAR, GENS
                          DO Y_CLASS = 1, Y_CLASSES
                               CONSUMP(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               CONSUMP_KIDS(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               LEISURE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               ASSETS(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               SHADOW_WAGE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.0000001
                               PENSION_BENEFITS_IND(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               PENSION_REPLACEMENT_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               AVG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 1.
                               MARG_WAGE_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 1.
                               AVG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               MARG_PENSION_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               AVG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               MARG_HEALTH_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               AVG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               MARG_DISABILITY_TAX_RATE(GEN, YEAR, Y_CLASS, COUNTRY) = 0.
                               TRANSFER(GEN, YEAR, Y_CLASS, COUNTRY) = 0.05
                               DO H_COUNTRY = FIRST_COUNTRY, LAST_COUNTRY
                                    H_TRANSFER(GEN, YEAR, Y_CLASS, COUNTRY, H_COUNTRY) = 0.05
                               ENDDO     
                          ENDDO
                     ENDDO
                ENDDO
           ENDDO
!           
           RETURN
      END SUBROUTINE INITIALIZE_VARIABLES

! Added block that helps to change initial assets easily
! new paarmeters for replacement rates calculating
! bigger health benefits as a share of productivity and smaller consumption tax parameter
! Russian initial productivity lowered to 33% of US's and time to converge with US's set to 40 years
! fertility rates for Russia adjusted (multyplied) for those not in 23-45 age group