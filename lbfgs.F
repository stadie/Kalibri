

      SUBROUTINE LVMINI(NPAR,MVEC,NFCN, AUX)
*     Initialization of Limited-memory Variable-metrik Minimization
      IMPLICIT NONE

*     argument variables -----------------------------------------------  
      INTEGER NPAR,MVEC,NFCN,JRET
      DOUBLE PRECISION X(*),F,AUX(*)
      REAL EPS,WLF1,WLF2

*     common variables -------------------------------------------------
      INTEGER       N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER

*     call variables ---------------------------------------------------
      INTEGER IRET,NRANK
      REAL WOLFC1,WOLFC2
      DOUBLE PRECISION ALPHA

*     functions --------------------------------------------------------
      INTEGER LVMDIM
      DOUBLE PRECISION DBDOT

*     internal variables -----------------------------------------------
      LOGICAL PRINT
      INTEGER I,J,IJ,NF,LITER,MF,NCALS
      REAL DFEPS,DFWLF1,DFWLF2,CONVEG
      DOUBLE PRECISION FPITER,FCLAST,FCRELT,ANGLGS,AGRADL,DGINIT,
     +   FCDACT,ALPHA0,FOPT,FPLUS,FMINUS,HIJ,HSTEP,DIAG,PROD,
     +   SRCH,ARGA,ACUT,PI  

      PARAMETER (PI=3.141 592 653 589 793 D0)

*     ==================================================================

*     ...

      PARAMETER (DFEPS=1.0E-4,DFWLF1=1.0E-4,DFWLF2=0.9) ! default const.

*     recognition of bad function values -------------------------------
      DOUBLE PRECISION TINY,FU
      PARAMETER (TINY=1.0D-10)
      LOGICAL NANINF
      DATA CONVEG/DFEPS/,WOLFC1/DFWLF1/,WOLFC2/DFWLF2/
      NANINF(FU)=(.NOT.FU+TINY.GT.0.0D0.AND..NOT.FU-TINY.LT.0.0D0).OR.
     +             ABS(FU).GT.1.0D+30
      SAVE
*     ...
      N=ABS(NPAR)             ! number of parameters 
      PRINT=.FALSE.
      IF(NPAR.LT.0) PRINT=.TRUE.

      LAST=LVMDIM(NPAR,MVEC)  ! necessary dimension

      M =MAX(1,ABS(MVEC))     ! number of history pairs, at least =1
      NF=NFCN                 ! maximum number of function evaluations
      IF(NF.LE.0) THEN        ! default value
         NF=91+(909*N*N)/(100+N*N)  ! 100 ... 1000 function evaluations 
         IF(MVEC.LT.0) NF=NF+N+N    ! + 2 * n for cov. matrix
      END IF 
      NCAL=0                  ! count number of function evaluations      
      ITER=0                  ! iteration counter 
      ICPT=0                  ! current pair-index (and index k)
      NCPT=0                  ! number of stored sy pairs
      IDER=0                  ! reset derivative flag
 
      ANGLGS=0.0D0

      DO I=1,LAST  
       AUX(I)=0.0D0 ! reset aux complete array
      END DO 
      RETURN
      
*     algorithmic constants --------------------------------------------

      ENTRY LVMEPS(EPS,WLF1,WLF2) ! change fit constants     
      IF(EPS.LE.0.0) THEN
         CONVEG=DFEPS   ! default value
      ELSE
         CONVEG=EPS     ! user value
      END IF

      IF(WLF1.LE.0.0) THEN
         WOLFC1=DFWLF1  ! default value
      ELSE
         WOLFC1=WLF1    ! user value
      END IF

      IF(WLF2.LE.0.0) THEN
         WOLFC2=DFWLF2  ! default value
      ELSE
         WOLFC2=WLF2    ! user value
      END IF
      RETURN
     
*     ------------------------------------------------------------------
      ENTRY LVMFUN(X,F,JRET,AUX)        ! called with new F and gradient
      NCAL=NCAL+1 ! increase counter for function evaluations
      IF(NCAL.GE.NF) GOTO 91            ! maximum number of functions
      AGRADL=SQRT(DBDOT(N,AUX,AUX))     ! actual |gradient| length
      IF(IDER.NE.0) GOTO 60     ! function for covariance calculation 
c      WRITE(*,*) 'new function is ',F
      IF(NCAL.NE.1) GOTO 30     ! function for minimization
*
*     first function value ---------------------------------------------
*      
      NB=0
      DO I=1,N
       IF(AUX(N+I).NE.0) NB=1 ! diagonal element not zero
       AUX(ISCA+I)=1.0D0      ! default scale factor
c       AUX(ISPT+I)=AUX(I)     ! for min derivative
c       AUX(IYPT+I)=AUX(I)     ! for max derivative 
       AUX(ISPT+I)=ABS(AUX(I))/(1.0D0+ABS(F))     ! for min derivative         
       AUX(IYPT+I)=ABS(AUX(I))/(1.0D0+ABS(F))     ! for max derivative
       AUX(IPAR+I)  =X(I)     ! initial parameter values
      END DO
      MF=1
      IF(NB.EQ.0) MF=3 

      IF(PRINT) THEN
         WRITE(*,*) ' '
         WRITE(*,*) 'Limited-memory BFGS routine LVMINI initialized',
     +              '========================'     
         WRITE(*,*)   '   with N =',N, '  dimension parameter'
         WRITE(*,*)   '        M =',M, '  number of history pairs'
         WRITE(*,102) '       NF =',NF,
     +                '   max number of function calls'   
         WRITE(*,101) '      eps =',CONVEG,'   convergence parameter'
         WRITE(*,101) '      c_1 =',WOLFC1,'   first Wolfe parameter'
         WRITE(*,101) '      c_2 =',WOLFC2,'   second Wolfe parameter'
         IF(NB.EQ.1) THEN
            WRITE(*,*)   'Diagonal of Hessian provided' 
         END IF
         WRITE(*,*) 'Requires AUX array with at least',LAST,' words'
         CALL LVMPRT(6,AUX,-1) ! print start values
         WRITE(*,*) ' '
         WRITE(*,103) 
      END IF

      IF(NANINF(F)) THEN        ! true for bad function value
         IF(PRINT) WRITE(*,*) 'Bad first function value: ',F
         GOTO 97
      END IF

*     first search vector is negative normalized gradient --------------

      IF(AGRADL.EQ.0.0D0) GOTO 92 ! zero gradient

      DO I=1,N
       AUX(IPGR+I)=AUX(I)             ! store gradient 
       AUX(ISEA+I)=-AUX(I)/AGRADL     ! normalized gradient
      END DO 

      FPITER=-0.5D0*DBDOT(N,AUX,AUX(ISEA+1)) ! gradient * step vector/2
      ALPHA=1.0D0 

      IF(PRINT) THEN
         WRITE(*,104) ITER,NCAL,NCPT,F
      END IF
      ITER=1
      FCLAST=F
      
*     ******************************************************************
*     *                        Minimization                            *
*     ******************************************************************
*
*     (re)initialize "Linesearch Algorithm with Sufficient Decrease"
 10   IF(NCAL.GE.NF) GOTO 91           ! maximum number of functions
      CALL MTLDEF(1.0E-3,WOLFC1,WOLFC2,1.0E-17,1.0E-20,1.0E+20, MF,10)
      MF=1
      DGINIT=DBDOT(N,AUX,AUX(ISEA+1))  ! initial slope
      IF(DGINIT.GE.0.0D0) THEN         ! ... is not negative
         WRITE(*,*) 'Huch ',DGINIT
         WRITE(*,*) 'AUX ',AUX(1),AUX(2)
         WRITE(*,*) 'D   ',AUX(ISEA+1),AUX(ISEA+2)
         GOTO 93 
      END IF
      GOTO 40

*     prepare next function-value/gradient calculation
 20   DO I=1,NB*N+N
       AUX(I)=0.0D0   ! reset AUX
      END DO      
      JRET=-1
      RETURN ! for next function-value/gradient evaluation
*     ... and at next entry ... 
*     new function-value/gradient calculated
 30   FCDACT=FCLAST-F                      ! actual function difference
      FCRELT=FCDACT/(FPITER+1.0D-6*CONVEG) ! ratio actual/expected
      IF(FCRELT.GT.(+99999.99D0)) FCRELT=+99999.99D0 ! for ...
      IF(FCRELT.LT.(-99999.99D0)) FCRELT=-99999.99D0 ! ... printing only

*     update scale factors ---------------------------------------------

      IF(NCPT.EQ.0) THEN     ! (s,y) pairs not yet stored 
         IF(NB.EQ.0) THEN    ! no diagonal elements present
            DO I=1,N
              AUX(ISPT+I)=AUX(ISPT+I)+ABS(AUX(I))/(1.0D0+ABS(F))
              AUX(ISCA+I)=DFLOAT(NCAL)/AUX(ISPT+I)
            END DO
         END IF 
      END IF ! (s,y) pairs not yet stored

      IF(PRINT.AND.ITER.NE.LITER) THEN
         IF(ITER.EQ.1) THEN
            WRITE(*,104) ITER,NCAL,NCPT,F,FPITER,FCRELT,ALPHA
         ELSE
            WRITE(*,104) ITER,NCAL,NCPT,F,FPITER,FCRELT,ALPHA,ANGLGS
         END IF
      ELSE IF(PRINT) THEN 
         WRITE(*,105)      NCAL,NCPT,F,   FCRELT,ALPHA
      END IF
      LITER=ITER   ! remember last iteration number

 40   CALL MTLINE(N,X,F,AUX(1),AUX(ISEA+1),ALPHA, AUX(IPAR+1),IRET)
*     value of IRET=-1   return for another function evaluation
*                  = 0   N or STP bad
*                  = 1   conditions hold
*                  = 2   interval at limit
*                  = 3   max number of calls
*                  = 4   lower bound
*                  = 5   upper bound
*                  = 6   rounding error limitation
*                  = 7   all below 10^-5

*     recognize fast convergence (without pair update) 
      IF(NCAL.EQ.NCALS+1.AND.ITER.GE.M.AND.MAX(FPITER,FCDACT).LE.CONVEG)
     +   THEN
         IF(FCDACT.LT.0.0D0) THEN ! step was increasing
            F=FCLAST
            DO I=1,N
             X(I)=AUX(IPAR+I)     ! back to optimal X
            END DO
         END IF
         GOTO 50
      END IF

      IF(IRET.LT.0) GOTO 20 ! ... for another function evaluation
      IF(PRINT.AND.IRET.NE.1) CALL MTINFO(IRET) ! print error info
      IF(IRET.EQ.0) GOTO 95 ! bad N or bad step in MTLINE 
      IF(IRET.EQ.1) THEN
*        sufficient decrease conditions hold - store next pair ---------
         CALL LVMUPD(ALPHA,AUX(ISEA+1),AUX,AUX(IPGR+1),AUX)
      END IF
*
*     L-BFGS two-loop recursion ----------------------------------------
*
      DO I=1,N
       AUX(ISEA+I)=AUX(I)             ! copy gradient
      END DO
      CALL LVMSOL(AUX(ISEA+1),AUX)    ! next search vector
      DO I=1,N
       AUX(ISEA+I)=-AUX(ISEA+I)       ! reverse sign of search vector
       AUX(IPGR+I)=AUX(I)             ! store last gradient
      END DO

      PROD=DBDOT(N,AUX,AUX(ISEA+1))
      SRCH=SQRT(DBDOT(N,AUX(ISEA+1),AUX(ISEA+1))) 
      ARGA=-PROD/(AGRADL*SRCH)
      IF(.NOT.ARGA.LT.1.0D0) ARGA=0.0D0 ! change illegal value 
      ANGLGS=180.0D0*ASIN(ARGA)/PI ! 90 degrees - ANGLGS is 

*     angle between the gradient direction and the search direction
*        is 90 degrees minus ANGLGS
*        minimization is inefficient for extremely small ANGLGS
*        If the angle is extremely small, the search direction is
*        slighly changed to the gradient direction                 

      ACUT=1.75D-6 ! corresponds to 0.0001 degree, ie. the angle
                   ! between gradient and search is 89.9999 degrees            
      IF(ARGA.LT.ACUT) THEN
         DO I=1,N 
          AUX(ISEA+I)=AUX(ISEA+I)-AUX(I)/AGRADL*SRCH*ACUT
         END DO
      END IF


*     start new line search iteration ----------------------------------

      FPITER=-0.5D0*PROD                     ! gradient * step vector/2

      ALPHA=1.0D0                     ! default start value
      IF(ITER.GT.0) THEN
         ALPHA0=(FCLAST-F)/FPITER     ! evtl. new smaller start value
         IF(ALPHA0.GT.0.0D0) ALPHA=MIN(1.0D0,1.01D0*ALPHA0)
c         IF(ALPHA0.GT.0.0D0) WRITE(*,*) 'ALPHA=',ALPHA
      END IF
      ITER=ITER+1
      FCLAST=F        ! previous function value

*
*     test for convergence ---------------------------------------------
*
      NCALS=NCAL ! number of function evaluations at new line-search
      IF(ITER.LE.M.OR.MAX(FPITER,FCDACT).GT.CONVEG) GOTO 10
*
*     convergence of minimization --------------------------------------
*

 50   AUX(IFUN)=F              ! save optimal function value
      DO J=1,N  ! use pair-info to estimate the inverse Hessian
       DO I=1,N
        AUX(ISEA+I)=0.0D0
       END DO                              ! unit vector e_j
       AUX(ISEA+J)=1.0D0
       CALL LVMSOL(AUX(ISEA+1),AUX)        ! solve -> matrix column
       AUX(IAPP+J)=SQRT(ABS(AUX(ISEA+J)))  ! store sqrt(diag)
      END DO
      IF(NCAL+N+N.GT.NF) INGR=0            ! skip covariance matrix
      IF(PRINT.AND.INGR.EQ.0.) THEN
         CALL LVMPRT(6,AUX,2) ! and parameters with approx. errors
         CALL LVMPRT(6,AUX,0) ! final function value ...
      END IF
      JRET=0
      IF(INGR.EQ.0) GOTO 90   ! finished (no covariance matrix)

*     ******************************************************************
*     *                    Covariance matrix                           *
*     ******************************************************************
*

*     derivatives ------------------------------------------------------
*
*     using: IPAR ... optimal parameter values, with FOPT
*            IPGR ... previous gradient (optimal gradient)
*            ISEA ... search vector (gradient at +)  
*            INGR ... gradient, numerical
*            INHS ... Hessian, numerical

 60   IF(IDER.EQ.0) THEN          

*        prepare error calculation ------------------------------------- 

         FOPT=F              ! save optimal function value
         DO I=1,N
          AUX(IPAR+I)=X(I)   ! save optimal values
          AUX(IPGR+I)=AUX(I) ! save gradient at optimal values
          AUX(INGR+I)=0.0D0  ! reset gradient, numerical  
         END DO
         DO I=1,(N*N+N)/2
          AUX(INHS+I)=0.0D0  ! reset Hessian, numerical
         END DO
         IF(PRINT) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'derivative calculation:'
            WRITE(*,*) '  par_i   step(+-)     secder   dF_minus',
     +      '    dF_plus'
         END IF
      ELSE IF(IDER.GT.0) THEN
         FPLUS=F                           ! F_+
         DO I=1,N
          AUX(ISEA+I)=AUX(I)               ! g_+                   
         END DO
         X(IDER)=AUX(IPAR+IDER)-HSTEP      ! X - h_j
         IDER=-IDER
         GOTO 20
      ELSE IF(IDER.LT.0) THEN
         FMINUS=F                          ! F_-
         IDER=-IDER
         X(IDER)=AUX(IPAR+IDER)            ! restore X 

         DO I=1,N
          HIJ=0.5D0*(AUX(ISEA+I)-AUX(I))/HSTEP
          IF(I.LT.IDER) THEN
             IJ=INHS+I+(IDER*IDER-IDER)/2
             AUX(IJ)=0.5D0*(AUX(IJ)+HIJ)     ! average
          ELSE IF(I.EQ.IDER) THEN
             AUX(INHS+(I*I+I)/2)=HIJ         ! diagonal element
          ELSE IF(I.GT.IDER) THEN
             AUX(INHS+IDER+(I*I-I)/2)=HIJ    ! store
          END IF
         END DO
         AUX(INDI+IDER)=(FPLUS-2.0D0*FOPT+FMINUS)/HSTEP**2 

         IF(PRINT) THEN 
            WRITE(*,106)IDER,HSTEP,AUX(INDI+IDER),FMINUS-FOPT,FPLUS-FOPT
         END IF
 
      END IF 

 70   IF(IDER.EQ.N) THEN      ! matrix is complete
*        restore
         F=FOPT               ! restore optimal function value
         DO I=1,N
          AUX(I)=AUX(IPGR+I)  ! ... optimal gradient
         END DO
*        matrix inversion:               diagonal    pointer
         CALL LVMINV(AUX(INHS+1),N,NRANK,AUX(ISEA+1),AUX(IPGR+1))
         DO I=1,N
          DIAG=AUX(INHS+(I*I+I)/2)             ! diagonal element
          AUX(INGR+I)=0.0D0
          IF(DIAG.GT.0.0D0) AUX(INGR+I)=SQRT(DIAG)  ! error
          DIAG=DIAG*AUX(ISEA+I)          
          AUX(INDI+I)=0.0D0
          IF(DIAG.GE.1.0D0) AUX(INDI+I)=SQRT(1.0D0-1.0D0/DIAG)!glcorr
         END DO
         IF(PRINT) CALL LVMPRT(6,AUX,6)
         JRET=0     
         GOTO 90
      END IF

      IDER=IDER+1                     ! next derivative
      HSTEP=AUX(IAPP+IDER) 
c      WRITE(*,*) 'HSTEP =',HSTEP
 
      AUX(INDI+IDER)=0.0D0            ! second derivative
      IF(HSTEP.EQ.0.0D0) THEN     ! fixed parameter
         IF(PRINT) THEN
            WRITE(*,106) IDER,HSTEP
         END IF
         GOTO 70
      END IF  
      X(IDER)=AUX(IPAR+IDER)+HSTEP    ! X + h_j
      GOTO 20

*
*     Error exits ------------------------------------------------------
*
 97   JRET=JRET+1  ! bad first function value  
 96   JRET=JRET+1  
 95   JRET=JRET+1  ! bad N or bad step in MTLINE  
 94   JRET=JRET+1  
 93   JRET=JRET+1  ! DGINIT is not negative 
 92   JRET=JRET+1  ! gradient zero 
 91   JRET=JRET+1  ! maximum number of function evaluations reached  
*
*     Final exit -------------------------------------------------------
*
 90   IF(PRINT) THEN
         WRITE(*,*) 'Limited-memory BFGS routine ending ',
     +              '==================================='
         WRITE(*,*) ' '
      END IF 

*     Format statements ------------------------------------------------ 

 101  FORMAT(1X,A,G10.3,A)
 102  FORMAT(1X,A,I10  ,A)
 103  FORMAT('Iter Nfcn Npair',6X,
     +       'Fcn_value   dFex      dFac/dFex   fstep',3X,' degr')
 104  FORMAT(I4,I5,I6,F15.4,    E11.2,F11.3,F8.3,F8.4) 
 105  FORMAT(4X,I5,I6,F15.4,11X,F11.3,F8.3)
 106  FORMAT(I8,G11.2,E11.2,2E11.2)
      END 


      INTEGER FUNCTION LVMDIM(NPAR,MVEC) ! pointer and total memory
      IMPLICIT NONE
      INTEGER NPAR,MVEC,NP,MV
*     common variables -------------------------------------------------
      INTEGER       N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER

      SAVE
*     ... 
      NP=ABS(NPAR)      
      MV=ABS(MVEC)

*     offset                  nr_of_words   description
*     ------------------      -----------   -------------------------- 
*          0                !     N         gradient 
*          NP               !     N         diagonal elements
      IPAR=NP+NP            !     N         parameter vector
      IPGR=IPAR+NP          !     N         previous gradient
      IAPP=IPGR+NP          !     N         approximate errors
      ISCA=IAPP+NP          !     N         scale factors
      ISEA=ISCA+NP          !     N         search vector
      IFUN=ISEA+NP+1        !     1         function value  
      IRHO=IFUN             !     M         rho(k)
!          IRHO+MV          !     M         alpha(k)
      ISPT=IRHO+MV+MV       !     M*N       s vectors
      IYPT=ISPT+MV*NP       !     M*N       y vectors

      INGR=0                !     unused
      INDI=0                !     unused
      INHS=0                !     unused
      LAST=IYPT+MV*NP       !     total space

      IF(MVEC.LT.0) THEN
         INGR=IYPT+MV*NP    !     N         gradient, numerical
         INDI=INGR+NP       !     N         diagonal, numerical 
         INHS=INDI+NP       !  (N*N+N)/2    Hessian, numerical              
         LAST=INHS+(NP*NP+NP)/2
      END IF
      LVMDIM=LAST
      END 


      INTEGER FUNCTION LVMIND(IARG) ! return index for ...
*     IARG = 0   function value
*          = 1   parameter vector
*          = 2   approximate errors
*          = 3   accurate errors
*          = 4   global correlations
*          = 5   covariance matrix
      IMPLICIT NONE
      INTEGER IARG
*     common variables -------------------------------------------------
      INTEGER       N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER
      SAVE
*     ...
      LVMIND=0         
      IF(IARG.EQ.(-1)) THEN
         LVMIND=N
      ELSE IF(IARG.EQ.(-2)) THEN
         LVMIND=NCAL
      ELSE IF(IARG.EQ.0) THEN
         LVMIND=IFUN
      ELSE IF(IARG.EQ.1) THEN 
         LVMIND=IPAR
      ELSE IF(IARG.EQ.2) THEN 
         LVMIND=IAPP
      ELSE IF(IARG.EQ.3) THEN 
         LVMIND=INGR
      ELSE IF(IARG.EQ.4) THEN 
         LVMIND=INDI
      ELSE IF(IARG.EQ.5) THEN 
         LVMIND=INHS
      END IF
      END

      SUBROUTINE LVMPRT(LUN,AUX,JARG) ! print on unit LUN ... 
      IMPLICIT NONE
      DOUBLE PRECISION AUX(*),RHO,DEI,DEJ
      REAL    CM(10)
      INTEGER LVMIND,LUN,JARG,LUP,NFUN,NPAR,IFUN,IPAR,IAPP,INGR,INDI,
     +        INHS,IARG1,IARG2,IARG,I,J,L,M
      SAVE

*          =-1  initial parameter values
*     IARG = 0  function value
*          = 1  parameter values
*          = 2  parameter values with approximate errors
*          = 3  parameter values with errors
*          = 4  parameter values with global correlations and errors
*          = 5  correlation matrix
*          = 6  means 0, 2, 4, 5  
*     ...
      LUP=LUN
      IF(LUP.EQ.0) LUP=6 ! standard output
      NFUN=LVMIND(-2)   
      NPAR=LVMIND(-1)
      IFUN=LVMIND(0)   ! get indices
      IPAR=LVMIND(1) 
      IAPP=LVMIND(2)
      INGR=LVMIND(3) 
      INDI=LVMIND(4)
      INHS=LVMIND(5) 

      IF(JARG.NE.6) THEN      
         IARG1=JARG
         IARG2=JARG
      ELSE
         IARG1=0
         IARG2=6  
      END IF

      DO IARG=IARG1,IARG2,2
      WRITE(LUP,*) ' '
      IF(IARG.EQ.0) THEN
         WRITE(LUP,*) NPAR,'-parameter function value =',AUX(IFUN),
     +                     ' after',NFUN,' function evaluations' 
      ELSE IF(IARG.EQ.(-1)) THEN  ! parameter values
         WRITE(LUP,*) '  Initial parameter values:'
         WRITE(LUP,101)    
         WRITE(LUP,102) (I,AUX(IPAR+I),I=1,NPAR)
      ELSE IF(IARG.EQ.1) THEN  ! parameter values 
         WRITE(LUP,101) 
         WRITE(LUP,102) (I,AUX(IPAR+I),I=1,NPAR)
      ELSE IF(IARG.EQ.2) THEN  ! parameter values with appr. errors
         WRITE(LUP,*) '  Optimization results ',
     +                 'with approximate errors:'
         WRITE(LUP,103)
         WRITE(LUP,104) (I,AUX(IPAR+I),AUX(IAPP+I),I=1,NPAR)
      ELSE IF(IARG.EQ.3.OR.IARG.EQ.4) THEN
         WRITE(LUP,*) '    Optimization results:'  
         IF(IARG.EQ.3) THEN
            WRITE(LUP,*)'    par_i         value             error'
         ELSE IF(IARG.EQ.4) THEN
            WRITE(LUP,*)'    par_i         value             error',
     +      '        glcorr'
         END IF
         DO I=1,NPAR 
          IF(IARG.EQ.3) THEN
             WRITE(LUP,105) I,AUX(IPAR+I),AUX(INGR+I)
          ELSE IF(IARG.EQ.4) THEN
             WRITE(LUP,105) I,AUX(IPAR+I),AUX(INGR+I),AUX(INDI+I) 
          END IF
         END DO 
      ELSE IF(IARG.EQ.5.OR.IARG.EQ.6) THEN         
         WRITE(LUP,*) '    Correlation coefficients [per mille]',
     +              ' between parameters i and j:'
         WRITE(LUP,106) (I,I=1,MIN(10,NPAR))
         DO I=1,NPAR
          DEI=AUX(INGR+I)
          IF(DEI.GT.0.0D0) THEN
             M=0
             DO J=1,I
              RHO=0.0D0
              DEJ=AUX(INGR+J)
              IF(DEI*DEJ.GT.0.0D0) THEN
                 RHO=AUX(INHS+J+(I*I-I)/2)/(DEI*DEJ)
              END IF
              M=M+1
              CM(M)=RHO
              IF(M.EQ.10) THEN
                 IF(J.LE.10) THEN
                    WRITE(LUP,107) I,CM
                 ELSE
                    WRITE(LUP,108)   CM
                 END IF
                 M=0
              END IF
             END DO
             IF(M.NE.0) THEN
                IF(J.LE.10) THEN
                   WRITE(LUP,107) I,(CM(L),L=1,M)
                ELSE
                   WRITE(LUP,108)   (CM(L),L=1,M)
                END IF
                M=0
             END IF
          ELSE
             WRITE(LUP,107) I
          END IF
         END DO
      END IF

      END DO 
      WRITE(LUP,*) ' '

*     Formats 
 101  FORMAT(3('   par_i       value    '))
 102  FORMAT(3(I8,G16.6))
 103  FORMAT(2('   par_i      value         error    '))
 104  FORMAT(2(I8,G15.6,' +-',G11.3))
 105  FORMAT(I10,G18.6,'+-',G16.6,F10.3)
 106  FORMAT('     par_i  \ j=',10I6)
 107  FORMAT(I10,6X,10F6.3)
 108  FORMAT(16X,   10F6.3)
      END 


      SUBROUTINE LVMUPD(FAC,SVEC,YVECP,YVECM,AUX)          ! pair update
      IMPLICIT NONE
      DOUBLE PRECISION FAC,SVEC(*),YVECP(*),YVECM(*),AUX(*),DBDOT,D(*)
      DOUBLE PRECISION DOTYS,DOTSBY,BETA

      INTEGER I,J

*     common variables -------------------------------------------------
      INTEGER       N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER
      COMMON/LMCOMM/N,M,NB,IPAR,IPGR,IAPP,ISCA,ISEA,IFUN,IRHO,ISPT,IYPT,
     +              INGR,INDI,INHS,LAST,NCPT,ICPT,ICP,NPT,NCAL,ITER,IDER

      SAVE
*     ... 

*     update set of vector pairs ---------------------------------------

      NPT =ICPT*N           ! offset before increasing ICPT
      NCPT=MIN(NCPT+1,M)    ! number of stored vectors
      ICPT=ICPT+1           ! increase pointer   ICPT=POINT
      IF(ICPT.EQ.M) ICPT=0       
      ICP =ICPT             ! pointer
      IF(ICP.EQ.0) ICP=M

      IF(NB.EQ.0) THEN      ! no diagonal 
         DO I=1,N
          AUX(ISPT+NPT+I)=FAC*SVEC(I)/AUX(ISCA+I)
          AUX(IYPT+NPT+I)=(YVECP(I)-YVECM(I))*AUX(ISCA+I)
         END DO
         DOTYS=DBDOT(N,AUX(ISPT+NPT+1),AUX(IYPT+NPT+1))
         DOTSBY=DOTYS/DBDOT(N,AUX(IYPT+NPT+1),AUX(IYPT+NPT+1))
      ELSE 
         DO I=1,N    ! new vectors with scaling factors
          AUX(ISPT+NPT+I)=FAC*SVEC(I)         ! param diff.
          AUX(IYPT+NPT+I)=(YVECP(I)-YVECM(I)) ! gradient diff.
         END DO
         DOTYS=DBDOT(N,AUX(ISPT+NPT+1),AUX(IYPT+NPT+1))  !!! should not be -
      END IF 

      AUX(IRHO+ICP)=0.0D0 ! rho(k)=1/(s(k)*y(k))
      IF(DOTYS.NE.0.0D0) AUX(IRHO+ICP)=1.0D0/DOTYS
      RETURN

*     ------------------------------------------------------------------

      ENTRY LVMSOL(D,AUX)       ! calculate B (vector d)
      IF(NB.EQ.0) THEN
         DO I=1,N
          D(I)=D(I)*AUX(ISCA+I)    ! scaling
         END DO
      END IF 

*     backward: alpha(i) := rho(i) * s(i)T*y(i); q := q-alpha(i)*y(i)

      ICP=ICPT 
      DO J=1,NCPT               ! backward loop
       ICP=ICP-1                ! ... with step -1 
       IF(ICP.LT.0) ICP=M-1     ! calculate alpha(j)    
       AUX(IRHO+M+ICP+1)=AUX(IRHO+ICP+1)*DBDOT(N,AUX(ISPT+ICP*N+1),D)
       CALL DBAXPY(N,-AUX(IRHO+M+ICP+1),AUX(IYPT+ICP*N+1),D)
      END DO

*     solve with approx Hessian ---------------------------------------- 

      IF(NB.EQ.0) THEN          ! no second derivative
         DO I=1,N
          D(I)=DOTSBY*D(I)      ! scale gradient
         END DO 
      ELSE IF(NB.EQ.1) THEN     ! diagonal matrix
         DO I=1,N               ! solve for approx. search direction
          IF(AUX(N+I).EQ.0.0D0) THEN
             D(I)=0.0D0         ! diagonal element zero 
          ELSE                  ! using aprox. Hessian ...
             D(I)=1.333D0*D(I)/ABS(AUX(N+I))
          END IF
         END DO 
      ELSE IF(NB.GT.1) THEN     ! band matrix
         STOP 'not implemented'
      END IF

*     forward loop -----------------------------------------------------

      DO J=1,NCPT               ! forward loop
       BETA=AUX(IRHO+ICP+1)*DBDOT(N,AUX(IYPT+ICP*N+1),D) 
       CALL DBAXPY(N,AUX(IRHO+M+ICP+1)-BETA,AUX(ISPT+ICP*N+1),D)
       ICP=ICP+1                ! ... with step +1 
       IF(ICP.EQ.M) ICP=0    
      END DO
      IF(NB.NE.0) RETURN

      DO I=1,N
       D(I)=D(I)*AUX(ISCA+I)    ! scaling back
      END DO

      END

*     matrix/vector products -------------------------------------------

      DOUBLE PRECISION FUNCTION DBDOT(N,DX,DY)   ! dot product
*     Dot product of two vectors: DX * DY
      IMPLICIT NONE
      INTEGER N,I 
      DOUBLE PRECISION DX(*),DY(*),DTEMP
*     ...
      DTEMP=0.0D0
      DO I = 1,MOD(N,5)
       DTEMP=DTEMP+DX(I)*DY(I)
      END DO
      DO I =MOD(N,5)+1,N,5
       DTEMP=DTEMP+DX(I)*DY(I)+DX(I+1)*DY(I+1)+DX(I+2)*DY(I+2)
     +                        +DX(I+3)*DY(I+3)+DX(I+4)*DY(I+4)
      END DO
      DBDOT=DTEMP
      END


      SUBROUTINE DBAXPY(N,DA,DX,DY)              ! multiply, addition
*     Constant times vector added to a vector: DY:=DY+DA*DX
      IMPLICIT NONE
      DOUBLE PRECISION DX(*),DY(*),DA
      INTEGER I,N
*     ...
      DO I=1,MOD(N,4)
        DY(I)=DY(I)+DA*DX(I)
      END DO
      DO I=MOD(N,4)+1,N,4
        DY(I  )=DY(I  )+DA*DX(I  )
        DY(I+1)=DY(I+1)+DA*DX(I+1)
        DY(I+2)=DY(I+2)+DA*DX(I+2)
        DY(I+3)=DY(I+3)+DA*DX(I+3)
      END DO
      END

*     matrix inversion ------------------------------------------------- 

      SUBROUTINE LVMINV(V,N,NRANK,DIAG,NEXT)   ! matrix inversion  
*           V = symmetric N-by-N matrix in symmetric storage mode
*               V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
*               replaced by inverse matrix
*
*     DIAG(N) =  double precision scratch array
*     NEXT(N) =  aux array (double precision instead of integer)
*
      IMPLICIT NONE
      INTEGER I,J,JJ,IJ,K,KK,JK,L,JL,LK,N,NRANK,NEXT0,LAST
      DOUBLE PRECISION V(*),DIAG(N),VKK,VJK,EPS,NEXT(N)
      PARAMETER (EPS=1.0D-10)
*     ...
      NEXT0=1                     ! start index
      DO I=1,N
       NEXT(I)=I+1                ! set "next" pointer
       DIAG(I)=ABS(V((I*I+I)/2))  ! save abs of diagonal elements
      END DO
      NEXT(N)=-1                  ! end flag

      NRANK=0
      DO I=1,N                    ! start of loop
       K  =0
       VKK=0.0D0

       J=NEXT0
       LAST=0 
 05    IF(J.GT.0) THEN
          JJ=(J*J+J)/2
          IF(ABS(V(JJ)).GT.MAX(ABS(VKK),EPS*DIAG(J))) THEN
             VKK=V(JJ)              
             K=J
             L=LAST
          END IF
          LAST=J
          J=NEXT(LAST)
          GOTO 05   
       END IF 

       IF(K.NE.0) THEN            ! pivot found
          KK=(K*K+K)/2
          IF(L.EQ.0) THEN
             NEXT0=NEXT(K)
          ELSE 
             NEXT(L)=NEXT(K)
          END IF
          NEXT(K)=0.0D0           ! index is used, reset
          NRANK=NRANK+1           ! increase rank and ...
          VKK    =1.0D0/VKK
          V(KK)  =-VKK
          JK     =KK-K
          JL     =0
          DO J=1,N                ! elimination
           IF(J.EQ.K) THEN
              JK=KK
              JL=JL+J
           ELSE
              IF(J.LT.K) THEN
                 JK=JK+1
              ELSE
                 JK=JK+J-1
              END IF
              VJK  =V(JK)
              V(JK)=VKK*VJK
              LK   =KK-K
              DO L=1,J
               JL=JL+1
               IF(L.EQ.K) THEN
                  LK=KK
               ELSE
                  IF(L.LT.K) THEN
                     LK=LK+1
                  ELSE
                     LK=LK+L-1
                  END IF
                  V(JL)=V(JL)-V(LK)*VJK
               END IF
              END DO
           END IF
          END DO
       ELSE
          DO K=1,N
           IF(NEXT(K).NE.0.0D0) THEN
              DO J=1,K
               IF(NEXT(J).NE.0.0D0) V((K*K-K)/2+J)=0.0D0  ! clear row/col
              END DO
           END IF
          END DO
          GOTO 10
       END IF
      END DO             ! end of loop
 10   DO IJ=1,(N*N+N)/2
       V(IJ)=-V(IJ)      ! finally reverse sign of all matrix elements
      END DO
      END


      SUBROUTINE LVMIDI(NPAR,PAR,FUNCT, AUX) ! minimization of [Minuit] function
      IMPLICIT NONE
      INTEGER NPAR,MVEC,IRET
      DOUBLE PRECISION PAR(*),FSUM,AUX(*)
      EXTERNAL FUNCT
      SAVE
*     ...      
      MVEC=MIN(ABS(NPAR),12)          ! number of delta pairs
      CALL LVMINI(NPAR,-MVEC,0, AUX) ! initialization

*     function and gradient evaluation
 10   CALL FUNCT(ABS(NPAR),AUX,FSUM,PAR,2)            ! calculate F, gradient
      CALL LVMFUN(PAR,FSUM, IRET,AUX)
      IF(IRET.LT.0)          GOTO 10
      END


*     line search subprograms ------------------------------------------

      SUBROUTINE MTLINE(N,X,F,G,S,STP, WA,INFO)
*     Linesearch Algorithms with "Guaranteed Sufficient Decrease"
*     by JJ More and D Thuente, ACM Transactions on Mathenatical
*     Software 20 (1994), pp 286-307 - modified by V.B.
*
*        ===============================================================
*        The algorithm is designed to find a step which satisfies 
*        the sufficient decrease condition
*              F(X+STEP*S) <= F(X) + FTOL*STP*[gradF(X).S]
*        and the curvature condition
*              ABS[gradF(X+STP*S).S] <= GTOL*ABS[gradF(X).S]
*
*        Typical values for Constants: 
*                               FTOL = 1.0E-4
*                               GTOL = 0.9  for Newton, Quasi-Newton
*                                    = 0.1  for CG
*
*        Conditions: 0 < FTOL < GTOL <  1   for Newton, Quasi-Newton
*                    0 < FTOL < GTOL < 1/2  for CG
*        ===============================================================
*      
*     CALL MTLDEF(...)      ! initialization, see end of routine
*     NC=0                  ! reset
*     define STP =  (in units of S(.), usually 1.0) 
* 01  calculate     (note: first value before first MTLINE call) 
*             F    = function value,
*             G(.) = gradient,
*             S(.) = search vector, at X(.) 
*     NC=NC+1               ! count function calls  
*     CALL MTLINE(N,X,F,G,S,STP, WA, INFO)
*        WA(.) =work array, of dimension N
*        INFO  =returned flag (not to be modified) 
*     IF(INFO.LT.0) GOTO 01 ! recalculate function ...  
*     IF(INFO.NE.1)         ! error condition
*     IF(INFO.EQ.1)         ! sufficient decrease and deriv. condition
*                           ! i.e. Wolfe conditions hold 
*     CALL MTINFO(INFO)     ! to print end condition
*     ------------------------------------------------------------------ 

      IMPLICIT NONE
      INTEGER J,N,MINF,MAXF,INFO,INFOC,NFEV,MINFE,MAXFE,NSF,NSS
      REAL PS,FTOLE,GTOLE,XTOLE,STMINE,STMAXE
      DOUBLE PRECISION F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
      DOUBLE PRECISION X(N),G(N),S(N),WA(N),DBDOT
      DOUBLE PRECISION SF(4),SFM,SEPS,TINY,FU,SS(10) 
      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
     *       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
     *       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO 
      LOGICAL NANINF,BRACKT,STAGE1
      PARAMETER (TINY=1.0D-10)
      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
      NANINF(FU)=(.NOT.FU+TINY.GT.0.0D0.AND..NOT.FU-TINY.LT.0.0D0).OR.
     +             ABS(FU).GT.1.0D+30 ! true if NaN or too big
      SAVE
*     ...
      SS(10)=0.0D0
*     -------------------------  initialization  ----------------------
      IF(NFEV.LT.0) THEN
         NFEV =0     ! reset  
         INFOC=1
         INFO =0 
         IF(N.LE.0.OR.STP.LE.ZERO) THEN 
            WRITE(*,*) 'MTLINE call with N=',N,', STP=',STP
            RETURN ! check N and step value
         END IF 
         DGINIT=DBDOT(N,G,S)

         IF(DGINIT.GE.ZERO) THEN
            WRITE(*,*) 'Initial DG is ge zero:',DGINIT,N
            RETURN ! check initial search direction (INFO=0)
         END IF 

         BRACKT= .FALSE.         ! local variables
         STAGE1= .TRUE.

         FINIT = F               ! initial function value
         DGTEST= FTOL*DGINIT     ! initial direct. derivative
         WIDTH = STPMAX - STPMIN ! width 
         WIDTH1= WIDTH/P5        
         DO J = 1, N
          WA(J) = X(J)           ! save initial parameter values
         END DO

         STX=ZERO    ! step is zero             at best step
         FX =FINIT   ! function value
         DGX=DGINIT  ! directional derivative

         STY=ZERO    ! step is zero             at other endpoint
         FY =FINIT   ! function value
         DGY=DGINIT  ! directional derivative

         SF(1)=0.0D0     ! sum of (F-FINIT) i.e. fct rel.to start
         SF(2)=0.0D0     ! sum of (F-FINIT)**2
         NSF=1

         SS(1)=0.0D0      ! sum of alpha     i.e. slope fit
         SS(2)=0.0D0      ! sum of alpha**2
         SS(3)=DGINIT     ! sum of DG
         SS(4)=0.0D0      ! sum of DG*alpha
         SS(5)=DGINIT**2  ! sum of DG**2
         NSS=1

      ELSE

*        ---------------- new function and gradient --------------------
         INFO=0
         NFEV=NFEV+1

         DG=DBDOT(N,G,S)

         IF(.NOT.NANINF(F)) THEN
            SS(1)=SS(1)+STP
            SS(2)=SS(2)+STP*STP
            SS(3)=SS(3)+DG
            SS(4)=SS(4)+DG*STP
            SS(5)=SS(5)+DG*DG
            NSS=NSS+1
            IF(NFEV.NE.1) THEN
               SF(1)=SF(1)+(F-FINIT)
               SF(2)=SF(2)+(F-FINIT)**2
               NSF=NSF+1
            END IF
         END IF

*        STP=        ! current step             at current step
*        F=          ! function value
*        DG=         ! directional derivative       
         FTEST1 = FINIT + STP*DGTEST

         IF((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     +       .OR. INFOC .EQ. 0)
     +                            INFO = 6 ! rounding error limitation 
         IF(STP .EQ. STPMAX .AND.
     +        F .LE. FTEST1 .AND. DG .LE. DGTEST)
     +                            INFO = 5 ! step at the upper bound
         IF(STP .EQ. STPMIN .AND.
     +       (F .GT. FTEST1 .OR. DG .GE. DGTEST))
     +                            INFO = 4 ! step at the lower bound
         IF(NFEV.GE.MAXF)         INFO = 3 ! max number of calls reached

         IF(BRACKT.AND.STMAX-STMIN.LE.XTOL*STMAX)
     +                            INFO = 2 ! relative interval width limit

         IF(F.LE.FTEST1.AND.ABS(DG).LE.GTOL*(-DGINIT).
     *       AND.NFEV.GE.MINF)    INFO = 1 ! min conditions hold ! 

         IF(NSF.GT.1) THEN
            SF(3)=SF(1)/DFLOAT(NSF)
            SF(4)=SQRT
     +            ((SF(2)-SF(1)*SF(1)/DFLOAT(NSF))/DFLOAT(NSF-1))
            SFM=MAX(ABS(SF(3)),SF(4))
         END IF

         SS(10)=0.0D0 
         IF(NSS.GT.2) THEN
            SS(6)=DFLOAT(NSS)*SS(2)-SS(1)*SS(1)         ! determinant
            SS(7)=(SS(2)*SS(3)-SS(1)*SS(4))/SS(6)       ! intercept
            SS(8)=(DFLOAT(NSS)*SS(4)-SS(1)*SS(3))/SS(6) ! slope
            SS(9)=SQRT(MAX(0.0D0,SS(5)-SS(3)*SS(7)
     +           -SS(4)*SS(8))/DFLOAT(NSS-2))           ! sigma(slope)
            SS(10)=-SS(7)/SS(8)                         ! minimum?
*           These are results from the slopes only, no function values
*           involved. Near the minimum the function values may be 
*           distorted by round-off errors (non-smooth curve). However
*           the slopes are smooth and approximately a straight line.
*           The minimum SS(10) is the value of zero slope from the
*           least-squares straight-line fit. Unused so far!
*           It can be probably used, if sigma small    
c            WRITE(*,*) 'Line ',SS(7),SS(8),' sigma=',SS(9),
c     +                 ' ->',SS(10)
         END IF

         IF(INFO.NE.1.AND.NSF.GE.MAX(MINF,MIN(3,MAXF/2))) THEN
            IF(SFM.LT.SEPS) THEN     ! round-off tess
               INFO=1
               IF(STP.NE.STX) THEN
                  STP=STX
                  DO J=1,N
                   X(J)=WA(J)+STP*S(J)
                  END DO
                  F=FX
               END IF  
            END IF
         END IF 

         IF(INFO.EQ.1) RETURN        ! convergence - min conditions hold

         IF(INFO.NE.0) THEN ! termination 
            STP=STX
            DO J=1,N
             X(J)=WA(J)+STP*S(J)
            END DO 
            RETURN
         END IF

*        ------------------------ define STP ---------------------------

*        first stage" seek a step for which the modified function has
*        a non-positive value and non-negative derivative

         IF(STAGE1.AND.F.LE.FTEST1.AND.
     +   DG.GE.MIN(FTOL,GTOL)*DGINIT) STAGE1=.FALSE.

         IF(STAGE1.AND.F.LE.FX.AND.F.GT.FTEST1) THEN

*           use modified function to predict a step, if no step was
*           obtained for which the modified function has a non-positive
*           function value and non-negative derivative, and if   
*           a lower function value has been obtained, but with
*           insufficient decrease

            FM   = F - STP*DGTEST    ! modify function value ...
            FXM  = FX - STX*DGTEST
            FYM  = FY - STY*DGTEST
            DGM  = DG - DGTEST      ! ... and gradient values
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST

*           update the interval and compute new step

            CALL CSTEPM(STX,FXM,DGXM,
     +                  STY,FYM,DGYM,
     +                  STP,FM ,DGM,
     +                  BRACKT,STMIN,STMAX,INFOC)

            FX  = FXM + STX*DGTEST  ! reset function value ...
            FY  = FYM + STY*DGTEST
            DGX = DGXM + DGTEST    ! ... and gradient values
            DGY = DGYM + DGTEST

         ELSE

*           update the interval and compute new step

            CALL CSTEPM(STX,FX,DGX,
     +                  STY,FY,DGY,
     +                  STP,F ,DG,                       ! update interval
     +                  BRACKT,STMIN,STMAX,INFOC)        ! and get new step
         END IF

         IF(BRACKT) THEN ! force sufficient decrease of size of interval
            IF(ABS(STY-STX).GE.P66*WIDTH1) STP=STX+P5*(STY-STX)
            WIDTH1=WIDTH
            WIDTH =ABS(STY-STX)
         END IF

      END IF

*     ---------- prepare calculation of next function value ------------

      IF(BRACKT) THEN               ! interval of uncertainty defined      
         STMIN=MIN(STX,STY)             ! minimum step
         STMAX=MAX(STX,STY)             ! maximum step
      ELSE                          ! interval of uncertainty undefined 
         STMIN=STX                      ! minimum step
         STMAX=STP+XTRAPF*(STP-STX)     ! maximum step
      END IF

      STP=MAX(STPMIN,MIN(STP,STPMAX)) ! STPMIN < STP < STPMAX

*     test unusual termination - STP = lowest point so far

      IF((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     +   .OR. NFEV .GE. MAXF-1 .OR. INFOC .EQ. 0
     +   .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) THEN
*        unusual termination
         STP = STX
      END IF
      INFO=-1 ! return to compute function value
      IF(SS(10).NE.0.0D0) THEN
c         WRITE(*,*) INFO,' sigma and steps: ',SFM, SS(9),SS(10),STP
         IF(SFM.LT.1.0D0) STP=SS(10)
      END IF
      DO J = 1, N          ! calculate the parameter values
       X(J)=WA(J)+STP*S(J) ! for step STP
      END DO
      RETURN

*     MTINFO print error information -----------------------------------

      ENTRY MTINFO(INFO)  ! print information about errors
c      IF(INFO.GE.0) WRITE(*,*) ' '
      IF(INFO.EQ.0) THEN
         WRITE(*,101) INFO,'N or STP lt 0 or step not descending'
      ELSE IF(INFO.EQ.1) THEN
         WRITE(*,101) INFO,'Linesearch convergence conditions reached'
      ELSE IF(INFO.EQ.2) THEN
         WRITE(*,101) INFO,'interval of uncertainty at lower limit'
      ELSE IF(INFO.EQ.3) THEN
         WRITE(*,101) INFO,'max nr of line search calls reached'
      ELSE IF(INFO.EQ.4) THEN
         WRITE(*,101) INFO,'step at the lower bound' 
      ELSE IF(INFO.EQ.5) THEN
         WRITE(*,101) INFO,'step at the upper bound'
      ELSE IF(INFO.EQ.6) THEN
         WRITE(*,101) INFO,'rounding error limitation'  
      ELSE IF(INFO.EQ.7) THEN
         WRITE(*,101) INFO,'function noise limit'
      END IF 
c      IF(INFO.GE.0) WRITE(*,*) ' '
 101  FORMAT('     Info =',I3,3X,A)
      RETURN

      ENTRY MTLDEF(PS,FTOLE, GTOLE, XTOLE,  STMINE, STMAXE, MINFE,MAXFE)
*        example:  2* 1.0E-4 1.0E-1 1.0E-17 1.0E-20 1.0E+20 1     40  

      SEPS=PS       ! round-off error tolerance 
      FTOL=FTOLE    ! tolerance C1  for sufficient decrease of function value
      GTOL=GTOLE    ! tolerance C2  for directional derivative condition
      XTOL=XTOLE    ! tolerance for width of intercal
      STPMIN=STMINE ! lower bound for the step
      STPMAX=STMAXE ! upper bound for the step
      MINF=MINFE    ! min number of linesearch function calls
      MAXF=MAXFE    ! max number of linesearch function calls (>1)
      IF(FTOL.LT.ZERO.OR.GTOL.LT.ZERO.OR.XTOL.LT.ZERO.OR.
     +   STMIN.LT.ZERO.OR.STPMAX.LT.STPMIN.OR.
     +   MINF.GT.MAXF.OR.MAXF.LE.1) THEN
         WRITE(*,*) 'Improper MTLDEF arguments - stop'
         STOP
      END IF 
      NFEV=-1               ! reset counter
      END


      SUBROUTINE CSTEPM(STX,FX,DX,  ! step, fvalue, der at best step
     +                  STY,FY,DY,  !                   at other endpt
     +                  STP,FP,DP,  !                   at current step 
     +                  BRACKT,        ! bracketed true/false
     +                  STPMIN,STPMAX, ! bounds for ste step
     +                  INFO)          ! returned: INFO=0 means error 
*                       INFO=1 ... 4   ! used method of interpolation
*     compute a safeguarded step for a linesearch
*     update an interval of uncertainty for the function minimum
      IMPLICIT NONE

      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      DOUBLE PRECISION TINY,FU,REDUCE
      PARAMETER (TINY=1.0D-10,REDUCE=0.1D0)
      LOGICAL NANINF,BRACKT,BOUND 
      NANINF(FU)=(.NOT.FU+TINY.GT.0.0D0.AND..NOT.FU-TINY.LT.0.0D0).OR.
     +             ABS(FU).GT.1.0D+30               ! NaN or infinity
*     ...
      INFO=0
      IF((BRACKT.AND.(STP.LE.MIN(STX,STY).OR.      ! check errors
     *    STP.GE.MAX(STX,STY))).OR.
     *    DX*(STP-STX).GE.0.0.OR.STPMAX.LT.STPMIN) RETURN

      SGND=DP*SIGN(1.0D0,DX)

      IF(NANINF(FP)) THEN
         INFO=5 
         BOUND=.TRUE.
         BRACKT=.TRUE. 
         STPF=REDUCE*STP
         FP=1.0D10
      ELSE IF(FP.GT.FX) THEN
         INFO = 1
*        higher function value: minimum is bracketed.
*        If the cubic step is closer to STX than the quadratic step,
*        the cubic step is taken, else the average of the two steps. 
         BOUND=.TRUE.
         THETA= 3.0D0*(FX - FP)/(STP - STX) + DX + DP
         S=MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA=S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF(STP.LT.STX) GAMMA = -GAMMA
         P=(GAMMA-DX)+THETA
         Q=((GAMMA-DX)+GAMMA)+DP
         R=P/Q
         STPC=STX+R*(STP-STX)
         STPQ=STX+((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP-STX)
         IF(ABS(STPC-STX).LT.ABS(STPQ-STX)) THEN
            STPF=STPC
         ELSE
            STPF=STPC+(STPQ-STPC)/2
         END IF
         BRACKT = .TRUE.
      ELSE IF(SGND.LT.0.0) THEN
         INFO = 2
*        lower function value and opp. derivative signs: minimum bracketed.
*        If the cubic step is closer to STX than the quadratic step,
*        the cubic step is taken, else the quadratic step.          
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF(ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
         END IF
         BRACKT = .TRUE.
      ELSEIF(ABS(DP).LT.ABS(DX)) THEN
         INFO = 3
*        lower function value and same derivative signs, der. decreases:
*        The cubic step is only used, if the cubic tends to infinity in
*        the driection of the step, or if the minimum of the cuic is beyond
*        STP. Otherwise the cubic step is defined to be either STPMIN or
*        STPMAX. The quadratic step also computed and if the minimum is
*        bracketed, then the step closest to STX is taken, else the step
*        farthest away is taken.
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF(R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
         END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF(BRACKT) THEN
            IF(ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         ELSE
            IF(ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         END IF
      ELSE
         INFO = 4
*        lower function value and same derivative signs, no der. decrease:
*        If the minimum is not bracketed, the step is either STPMIN or
*        STPMAX, else the cubic step is taken. 
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
         END IF
      END IF

*     update the interval of uncertainty (does not depend on new step) 
      IF(FP.GT.FX) THEN
         STY=STP
         FY =FP
         DY =DP
      ELSE
         IF(SGND.LT.0.0) THEN
            STY=STX
            FY =FX
            DY =DX
         END IF
         STX=STP
         FX =FP
         DX =DP
      END IF

      STPF=MAX(STPMIN,MIN(STPF,STPMAX))      
      STP =STPF                       ! new step
      IF(BRACKT.AND.BOUND) THEN       ! safeguard the step
         IF(STY.GT.STX) THEN
            STP= MIN(STX+0.66D0*(STY-STX),STP)
         ELSE
            STP= MAX(STX+0.66D0*(STY-STX),STP)
         END IF
      END IF
      END



 
