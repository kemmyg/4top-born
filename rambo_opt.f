c  Downloaded 17.10.99 from
c  http://www.cpc.cs.qub.ac.uk/cpc/cgi-bin/list_summary.pl?CatNumber=AAFU
c  Modifications:
c    Comment out sample program.   G. Cowan 17 October, 1999
c    Replace RN by RANMAR (CERNLIB V113).  G. Cowan, 17 October, 1999

c AAFURAMBO.  A NEW MONTE CARLO TREATMENT OF MULTIPARTICLE PHASE SPACE AT
c HIGH ENERGIES.  R. KLEISS, W.J. STIRLING, S.D. ELLIS.
c REF. IN COMP. PHYS. COMMUN. 40 (1986) 359
      SUBROUTINE RAMBO(N,ET,XM,P,WT,LW)                                  
C------------------------------------------------------                  
C                                                                        
C                       RAMBO                                            
C                                                                        
C             RA(NDOM)  M(OMENTA)  BO(OSTER)                             
C                                                                        
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR                   
C    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING                    
C                                                                       
C    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)                
C    ET = TOTAL CENTRE-OF-MASS ENERGY                                   
C    XM = PARTICLE MASSES ( DIM=N )                                     
C    P  = PARTICLE MOMENTA ( DIM=(4,N) )                                
C    WT = WEIGHT OF THE EVENT                                           
C    LW = FLAG FOR EVENT WEIGHTING:                                     
C         LW = 0 WEIGHTED EVENTS                                        
C         LW = 1 UNWEIGHTED EVENTS ( FLAT PHASE SPACE )                 
C------------------------------------------------------                 
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION XM(100),P(4,100),Q(4,100),Z(100),R(4),                      
     .   B(3),P2(100),XM2(100),E(100),V(100),IWARN(5)                   
      DATA ACC/1.D-14/,ITMAX/6/,IBEGIN/0/,IWARN/5*0/
clk      DATA ACC/1.D-8/,ITMAX/6/,IBEGIN/0/,IWARN/5*0/
C     LK
      IBEGIN = 0
c      print*, 'iwarn',iwarn
C      print*, 'wt begin', wt
c      iwarn(1)=5
C INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT            
      IF(IBEGIN.NE.0) GOTO 103
c      print*, 'entered ibegin if'
c      IBEGIN=1                                                          
      TWOPI=8.*DATAN(1.D0)                                              
      PO2LOG=DLOG(TWOPI/4.)                                             
      Z(2)=PO2LOG                                                       
      DO 101 K=3,100                                                    
  101 Z(K)=Z(K-1)+PO2LOG-2.*DLOG(DFLOAT(K-2))                           
      DO 102 K=3,100                                                    
  102 Z(K)=(Z(K)-DLOG(DFLOAT(K-1)))                                     
C                                                                       
C CHECK ON THE NUMBER OF PARTICLES                                      
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104                                  
      PRINT 1001,N                                                      
      STOP                                                              
C                                                                       
C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES        
  104 XMT=0.                                                            
      NM=0                                                              
      DO 105 I=1,N                                                      
      IF(XM(I).NE.0.D0) NM=NM+1                                         
  105 XMT=XMT+DABS(XM(I))                                               
      IF(XMT.LE.ET) GOTO 106                                            
      PRINT 1002,XMT,ET                                                 
      STOP                                                              
C                                                                       
C CHECK ON THE WEIGHTING OPTION                                         
  106 IF(LW.EQ.1.OR.LW.EQ.0) GOTO 201                                   
      PRINT 1003,LW                                                     
      STOP                                                              
C                                                                       
C THE PARAMETER VALUES ARE NOW ACCEPTED                                 
C     
C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE 
  201 DO 202 I=1,N                                                      
      C=2.*RN(1)-1.                                                     
      S=DSQRT(1.-C*C)                                                   
      F=TWOPI*RN(2)
      Q(4,I)=-DLOG(RN(3)*RN(4)) 
c      print*, -DLOG(RN(3)*RN(4)), RN(3), RN(4),RN(3)*RN(4)
      Q(3,I)=Q(4,I)*C                                                   
      Q(2,I)=Q(4,I)*S*DCOS(F)                                           
  202 Q(1,I)=Q(4,I)*S*DSIN(F)                                           
C               
C      print*, 'wt before conformal', wt
C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION              
      DO 203 I=1,4                                                      
  203 R(I)=0.                                                           
      DO 204 I=1,N                                                      
      DO 204 K=1,4                                                      
  204 R(K)=R(K)+Q(K,I)
c      print*,'Q',(Q(4,i),i=1,4)
c      print*,'R', R
      RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)   
c      print*, 'RMAS', RMAS
      DO 205 K=1,3                                                      
  205 B(K)=-R(K)/RMAS                                                   
      G=R(4)/RMAS                                                       
      A=1./(1.+G)
      X=ET/RMAS                                                         
c      print*, 'X',X
c      print*, 'RMAS', RMAS
C                                                                       
C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S  
c      print*, 'wt q to p', wt
      DO 207 I=1,N                                                      
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)                            
      DO 206 K=1,3                                                      
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))                              
  207 P(4,I)=X*(G*Q(4,I)+BQ)                                            
C                                                                       
C RETURN FOR UNWEIGHTED MASSLESS MOMENTA                                
      WT=1.D0                                                           
      IF(NM.EQ.0.AND.LW.EQ.1) RETURN                                    
C                                                                       
C CALCULATE WEIGHT AND POSSIBLE WARNINGS                                
      WT=PO2LOG
c      print*, 'wt weight possible warnings', wt
      IF(N.NE.2) WT=(2.*N-4.)*DLOG(ET)+Z(N)                             
      IF(WT.GE.-180.D0) GOTO 208  
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
cl    k-------------------------
 208  IF(WT.LE. 174.D0) GOTO 209                                        
      IF(IWARN(2).LE.5) PRINT 1005,WT                                   
      IWARN(2)=IWARN(2)+1                                               
C                                                                       
C RETURN FOR WEIGHTED MASSLESS MOMENTA                                  
  209 IF(NM.NE.0) GOTO 210                                              
      WT=DEXP(WT)                                                       
      RETURN                                                            
C                                                                       
C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X                  
  210 XMAX=DSQRT(1.-(XMT/ET)**2)  
      DO 301 I=1,N                                                      
      XM2(I)=XM(I)**2                                                   
  301 P2(I)=P(4,I)**2                                                  
      ITER=0                                                           
      X=XMAX                                                           
      ACCU=ET*ACC                                                      
  302 F0=-ET
c      print*, 'F0,ACC and ACCU beginnig',F0,ACC,ACCU
      G0=0.                                                            
      X2=X*X                                                           
      DO 303 I=1,N   
c     lk here is the knackpunkt. E(I) goes NaN sometimes due to X2 and P2(I) 
c      print*,'X2,P2(I) ',X2,P2(I)
c,XM2,P2(I)   
      E(I)=DSQRT(XM2(I)+X2*P2(I))
c     print*, 'E(I)', E(I)
c      print*, 'F0 before',F0                                      
      F0=F0+E(I)          
c      print*, 'F0 after', F0
  303 G0=G0+P2(I)/E(I) 
c      print*,'dabs(F0), accu 303',(F0),ACCU
      IF(DABS(F0).LE.ACCU) GOTO 305                                    
      ITER=ITER+1                                                      
      IF(ITER.LE.ITMAX) GOTO 304                                       
      PRINT 1006,ITMAX                                                 
      GOTO 305                                                         
  304 X=X-F0/(X*G0)                                                    
      GOTO 302                                                         
  305 DO 307 I=1,N                                                     
      V(I)=X*P(4,I)                                                    
      DO 306 K=1,3                                                     
  306 P(K,I)=X*P(K,I)                                                  
  307 P(4,I)=E(I)                                                      
C                                                                      
C CALCULATE THE MASS-EFFECT WEIGHT FACTOR                              
      WT2=1.                                                           
      WT3=0.                                                           
      DO 308 I=1,N                                                     
      WT2=WT2*V(I)/E(I)                                                
  308 WT3=WT3+V(I)**2/E(I)                                             
      WTM=(2.*N-3.)*DLOG(X)+DLOG(WT2/WT3*ET)                           
      IF(LW.EQ.1) GOTO 401                                             
C                                                                      
      WT=WT+WTM                                                        
      IF(WT.GE.-180.D0) GOTO 309                                       
      IF(IWARN(3).LE.5) PRINT 1004,WT
      IWARN(3)=IWARN(3)+1                                              
  309 IF(WT.LE. 174.D0) GOTO 310                                       
      IF(IWARN(4).LE.5) PRINT 1005,WT                                  
      IWARN(4)=IWARN(4)+1                                              
  310 WT=DEXP(WT)                                                      
      RETURN                                                           
C                                                                      
C UNWEIGHTED MASSIVE MOMENTA REQUIRED: ESTIMATE MAXIMUM WEIGHT         
  401 WT=DEXP(WTM)                                                     
      IF(NM.GT.1) GOTO 402                                             
C                                                                      
C ONE MASSIVE PARTICLE                                                 
      WTMAX=XMAX**(4*N-6)                                              
      GOTO 405                                                         
  402 IF(NM.GT.2) GOTO 404                                             
C                                                                      
C TWO MASSIVE PARTICLES                                                
      SM2=0.                                                           
      PM2=0.                                                           
      DO 403 I=1,N                                                     
      IF(XM(I).EQ.0.D0) GOTO 403                                       
      SM2=SM2+XM2(I)                                                   
      PM2=PM2*XM2(I)                                                   
  403 CONTINUE                                                         
      WTMAX=((1.-SM2/(ET**2))**2-4.*PM2/ET**4)**(N-1.5)                
      GOTO 405                                                         
C                                                                      
C MORE THAN TWO MASSIVE PARTICLES: AN ESTIMATE ONLY                    
  404 WTMAX=XMAX**(2*N-5+NM)                                           
C                                                                      
C DETERMINE WHETHER OR NOT TO ACCEPT THIS EVENT                        
  405 W=WT/WTMAX                                                       
      IF(W.LE.1.D0) GOTO 406                                           
      IF(IWARN(5).LE.5) PRINT 1007,WTMAX,W                             
      IWARN(5)=IWARN(5)+1                                              
  406 CONTINUE                                                         
      IF(W.LT.RN(5)) GOTO 201                                          
      WT=1.D0                                                          
      RETURN                                                           
 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')    
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',             
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)                           
 1003 FORMAT(' RAMBO FAILS: LW=',I3,' IS NOT AN ALLOWED OPTION')       
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')    
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')    
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',        
     . ' DESIRED ACCURACY =',D15.6)                                    
 1007 FORMAT(' RAMBO WARNS: ESTIMATE FOR MAXIMUM WEIGHT =',D15.6,      
     . '     EXCEEDED BY A FACTOR ',D15.6)                             
      END                                                              
C  REPLACE RN BY RANMAR (CERNLIB V113).  G. COWAN, 17 OCTOBER, 1999
      REAL*8 FUNCTION RN(IDMY)
      IMPLICIT NONE
      INTEGER         IDMY
      REAL            RVEC(1)
c 222  CALL RANMAR (RVEC, 1)
c     only use random_number, ranmar not good
 222  call random_number(rvec)
      RN = RVEC(1)
      IF (RN .eq. 0d0) THEN 
c         PRINT*, 'RN = 0D0',RN
         GOTO 222 
      ENDIF
c      PRINT*, 'RN(',IDMY,')= ',RN
      RETURN
      END

      SUBROUTINE RANMAR(RVEC,LENV)
C ----------------------------------------------------------------------
C<<<<<FUNCTION RANMAR(IDUMM)
C CERNLIB V113, VERSION WITH AUTOMATIC DEFAULT INITIALIZATION
C     TRANSFORMED TO SUBROUTINE TO BE AS IN CERNLIB
C     AM.LUTZ   NOVEMBER 1988, FEB. 1989
C
C!UNIVERSAL RANDOM NUMBER GENERATOR PROPOSED BY MARSAGLIA AND ZAMAN
C IN REPORT FSU-SCRI-87-50
C        MODIFIED BY F. JAMES, 1988 AND 1989, TO GENERATE A VECTOR
C        OF PSEUDORANDOM NUMBERS RVEC OF LENGTH LENV, AND TO PUT IN
C        THE COMMON BLOCK EVERYTHING NEEDED TO SPECIFY CURRRENT STATE,
C        AND TO ADD INPUT AND OUTPUT ENTRY POINTS RMARIN, RMARUT.
C
C     UNIQUE RANDOM NUMBER USED IN THE PROGRAM
C ----------------------------------------------------------------------
      COMMON / INOUT / INUT,IOUT
      DIMENSION RVEC(*)
      COMMON/RASET1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        DEFAULT INITIALIZATION. USER HAS CALLED RANMAR WITHOUT RMARIN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RMARIN(IJKLIN, NTOTIN,NTOT2N)
C         INITIALIZING ROUTINE FOR RANMAR, MAY BE CALLED BEFORE
C         GENERATING PSEUDORANDOM NUMBERS WITH RANMAR. THE INPUT
C         VALUES SHOULD BE IN THE RANGES:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C TO GET THE STANDARD VALUES IN MARSAGLIA-S PAPER, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          ALWAYS COME HERE TO INITIALIZE
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
c      WRITE(IOUT,201) IJKL,NTOT,NTOT2
c 201  FORMAT(1X,' RANMAR INITIALIZED: ',I10,2X,2I10)
      DO 2 II= 1, 97
      S = 0.
      T = .5
      DO 3 JJ= 1, 24
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = 0.5*T
    2 U(II) = S
      TWOM24 = 1.0
      DO 4 I24= 1, 24
    4 TWOM24 = 0.5*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       COMPLETE INITIALIZATION BY SKIPPING
C            (NTOT2*MODCNS + NTOT) RANDOM NUMBERS
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
       WRITE (IOUT,'(A,I15)') ' RMARIN SKIPPING OVER ',NOW
       DO 40 IDUM = 1, NTOT
       UNI = U(I97)-U(J97)
       IF (UNI .LT. 0.)  UNI=UNI+1.
       U(I97) = UNI
       I97 = I97-1
       IF (I97 .EQ. 0)  I97=97
       J97 = J97-1
       IF (J97 .EQ. 0)  J97=97
       C = C - CD
       IF (C .LT. 0.)  C=C+CM
   40  CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          NORMAL ENTRY TO GENERATE LENV RANDOM NUMBERS
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. 0.)  UNI=UNI+1.
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. 0.)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. 0.) UNI=UNI+1.
C        REPLACE EXACT ZEROES BY UNIFORM DISTR. *2**-24
         IF (UNI .EQ. 0.)  THEN
         UNI = TWOM24*U(2)
C             AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LETS BE SAFE.
         IF (UNI .EQ. 0.) UNI= TWOM24*TWOM24
         ENDIF
      RVEC(IVEC) = UNI
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           ENTRY TO OUTPUT CURRENT STATUS
      ENTRY RMARUT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END

c      real*8 FUNCTION RN(IDMY)                                                 AAFU0198
c      DATA I/65539/                                                     AAFU0199
cc      DATA C/Z39200000/                                                 AAFU0200
c      DATA C/1/                                                          AAFU0200
cc      RN = 1.
c    1 IDMY=IDMY                                                         AAFU0201
c      I=I*69069                                                         AAFU0202
c      IF(I.LE.0) I=I + 2147483647 + 1                                   AAFU0203
c      J=I/256                                                           AAFU0204
c      RN=C*FLOAT(256*J)                                                 AAFU0205
c      IF(RN.NE.0.) RETURN                                               AAFU0206
c      PRINT 2,I,C,J,RN                                                  AAFU0207
c    2 FORMAT(' RN WARNING: WRONG VALUE OCCURRED'/,                      AAFU0208
c     .       ' I,C,J,RN  :',I15,D15.6,I15,D15.6)                        AAFU0209
c      GOTO 1                                                            AAFU0210
c      END                                                               AAFU0211

c      IMPLICIT REAL*8(A-H,O-Z)                                          AAFU0212
c      DIMENSION XM(30),P(4,30)                                          AAFU0213
c      DO 1 I=1,30                                                       AAFU0214
c    1 XM(I)=0.D0                                                        AAFU0215
c      ET=100.D0                                                         AAFU0216
c      CALL RAMBO(30,ET,XM,P,W,0)                                        AAFU0217
c      PRINT 1001,ET,(P(K,1),K=1,4)                                      AAFU0218
c 1001 FORMAT('1',50('*')/,                                              AAFU0219
c     .       ' THIS IS A RAMBO TESTRUN  '/,                             AAFU0220
c     .       ' AT A TOTAL INVARIANT MASS OF'/,D30.8/,                   AAFU0221
c     .       '0THE FIRST MOMENTUM:'/,(D30.8))                           AAFU0222
c      PRINT 1002,W                                                      AAFU0223
c 1002 FORMAT('0THE PHASE SPACE VOLUME FOR 30 MASSLESS PARTICLES IS'/,   AAFU0224
c     .       D30.8)                                                     AAFU0225
c      DO 2 I=1,30                                                       AAFU0226
c    2 XM(I)=1.D0                                                        AAFU0227
c      VOL=0.D0                                                          AAFU0228
c      DO 3 I=1,1000                                                     AAFU0229
c      CALL RAMBO(30,ET,XM,P,W,0)                                        AAFU0230
c    3 VOL=VOL+W                                                         AAFU0231
c      VOL=VOL/1000.D0                                                   AAFU0232
c      PRINT 1003,VOL                                                    AAFU0233
c 1003 FORMAT('0THE PHASESPACE VOLUME FOR 30 ',                          AAFU0234
c     . 'PARTICLES WITH MASS=1 IS'/,D30.8,'  (ESTIMATED)'/,'0',50('*'))  AAFU0235
c      STOP                                                              AAFU0236
c      END                                                               AAFU0237
                                                                        AAFU****
