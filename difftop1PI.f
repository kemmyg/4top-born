C
C
C ABSTRACT:
C        DiffTop (vers. 1.0.0) program computes total cross section or
C        differential distributions in pT and y for hadronic
C        ttbar reactions at LO, approx NLO, approx NNLO,
C        by using techniques of QCD threshold resummation
C        at LL, NLL and NNLL accurracy, in pico barn.
C
C AUTHOR(S):
C           Sven-Olaf Moch and Marco Guzzi, DESY HAMBURG
C
C
C
C PAPER(S):
C            'Top-quark pair production at hadron colliders:
C             differential cross section and phenomenological applications with DiffTop',
C             by Marco Guzzi, Katerina Lipka and Sven-Olaf Moch,
C             e-Print: arXiv:1406.0386 [hep-ph].


      program difftop
      implicit none

      integer ptychoice, naccur
      integer nqqborn, nggborn, nggnlo, nqqnlo, nggnnlo, nqqnnlo
      integer iloop, jloop
      integer kount, max
      integer ipdf, iset
      integer PDFflag
      integer ptsi, itsi, dimi, ptse, itse, dime, nel
      integer ptsb, itsb, dimb
      integer nalfas, Oalphas
      integer sclset             ! Flag for the scale choice (0) scale = m; (1) scale = sqrt(m^2+pt^2)                            
      integer Npt, Nyt           ! Num. points in Pt and rapidity                                                                 
      integer Nptp1, Nytp1
      integer N,cin,calls                   ! Number of particles for RAMBO
      integer i,j,iter
      integer choice, N_Ehat, N_beta,N_E
      real*8 mq,muR,kR,muF,mu2      ! Ren./Fact. scales                                                                              
      real*8 m2, scms, dscms, scmsroot, scale, alphas
      real*8 nlf
      real*8 epsi, epse, epsb, delta
      real*8 shad, shadinel , shadel, betacap
      real*8 pi, ca, cf, co, ck, cqed, tf, nf
      real*8 mfac,pt2,yyrap,pt2rt
      real*8 pt2rtmin, pt2rtmax      ! Min and Max values of pT                                                                   
      real*8 yyrapmin, yyrapmax      ! MIn and Max valueas of y  
      real*8 shaderr, shdinrr , shdelrr
      real*8 qcdl1, qcdl2, qcdl3     ! Flag for using pdf database          
      real*8 s1,s2,s3,s4
      real*8 xm,p,w,wt,ET,Vol,rndm
      real*8 p2,E
      real*8 wgt,E_i
      real*8 Ehat,beta,wdth_Ehat,wdth_beta,wdth_E
c      dimension xm(4),p(4,4)

      


      character*10 q

      parameter (pi = 3.14159265359d0)
      parameter (N = 100)
      common/calls/ptsb,cin,calls
      common/result/s1,s2,s3,s4
      common/group/ca, cf, co, ck, cqed, tf, nf
      common/m2only/m2
      common/sonly/scms
      common/include/ptychoice, naccur, nqqborn, nggborn,
     #     nggnlo, nqqnlo, nggnnlo, nqqnnlo

      common/coupling/scale, alphas, nlf
      common/sclchoice/mq, muR, kR, muF, sclset
      common/quark/q
      common/pdfdb/qcdl1,qcdl2,qcdl3,nalfas
      common/plus/delta
      common/channel/nel
      common/ptyblock/pt2,yyrap
      common/PDFname/PDFflag
      common/pdfiset/iset
      common/pdf/ipdf
      common/Ehatbeta/Ehat,beta,E_i
      common/choice/choice
      

     
c The function ff gives the total cross section
c to be integrated by VEGAS

      external ff



c      open (unit = 11, file ='difftop.in', status = 'old')
      open (unit = 11, file ='difftop.in')  

c     open (unit = 12, file ='difftop.out')

      open (unit = 12,
     &     file ='histo_4top_Ehat_13_173_gg_mu2mq_abmp16_lumi_
     &20_bins.dat')
      open (unit = 13,
     &     file ='histo_4top_beta_13_173_gg_mu2mq_abmp16_lumi_
     &20_bins.dat')
      open (unit = 14, file =
c     & 'total_4top_qqb-anll_a25_1-14_173_mu2mq_abmp16.dat')
     & 'TEST.dat')
      
      rewind 11
      rewind 12
      rewind 13
      rewind 14

      read(11,*) ptychoice      ! it decides which distrbution or Total Xsec                                                      

      read(11,*) PDFflag        ! it decides whether or not initializing LHAPDFs                                                  
      read(11,*) iset           ! it specifies the PDF set xxx.LHgrid from the input file                                         
      read(11,*) nalfas
      read(11,*) qcdl1          ! LambdaQCD 1-loop                                                                                
      read(11,*) qcdl2          ! LambdaQCD 2-loop                                                                                
      read(11,*) qcdl3          ! LambdaQCD 3-loop                                                                                

      read(11,*) ipdf           ! identifies standalone PDFs when PDFflag=0 or LHAPDF when PDFflag=1                              
      read(11,*) naccur
      naccur=0
c      nel = 0
      read(11,*) nggborn
      read(11,*) nqqborn
      read(11,*) nggnlo
      read(11,*) nqqnlo
      read(11,*) nggnnlo
      read(11,*) nqqnnlo
      read(11,*) sclset           ! scale choice (0) scale = m; (1) scale = sqrt(m^2+pt^2)                                        

      read(11,*) mfac             ! it varies the factorization scale                                                             
      read(11,*) kR               ! it varies the renormalization scale                                                           
      read(11,*) mq
      read(11,*) scmsroot
      read(11,*) nlf
      read(11,*) q

      read(11,*) pt2rtmin        ! Min val in pT grid                                                                             
      read(11,*) pt2rtmax        ! Max val in pT grid                                                                             
      read(11,*) Npt             ! Num points for (Npt+1) pT-grid                                                                 

      read(11,*) yyrapmin        ! Min val in y grid                                                                              
      read(11,*) yyrapmax        ! Max val in y grid                                                                              
      read(11,*) Nyt             ! Num points for (Ny+1) y-grid                                                                   

      read(11,*) delta           ! cut-off parameter                                                                              
      read(11,*) choice

 87   format(a90)

      print*,'*      WELCOME to DIFF4TOP       *'

      
      if(PDFflag.eq.1) then
         print*,'Initializing LHA PDFs...'
         call InitialPDF
         call GetOrderAs(Oalphas)
         print*,'=============================='
         print*,'Alphas order from LHAPDF grid:', Oalphas
         print*,'=============================='
      endif



      print*,'*****************************************'
      print*,'*                                       *'
      print*,'*   PARAMETERS READING difftop.in file  *'
      print*,'*                                       *'
      print*,'*****************************************'


      print*,' choice of distrib. /pt/y/tot-pT/tot-y/ -> /0/1/2/3/:',
     >     ptychoice

      print*,' PDFs initialization: PDFflag = ',PDFflag
      print*,' parton PDF Set: .LHgrid number = ',iset
      print*,' nalfas flag for alphas evolution = ',nalfas !LHAPDF evol. nalfas=(99); internal alphas evol. nalfas=(2,3)          
      print*,' ipdf flag for parton evolution = ',ipdf
      print*,' accurracy:  LL (0), NLL (1), NNLL (2): ', naccur
      print*,' include gg born (1) not (0): nggborn = ',nggborn
      print*,' include qq born (1) not (0): nqqborn = ',nqqborn
      print*,' include one loop gg (1) not (0): nggnlo = ',nggnlo
      print*,' include one loop qq (1) not (0): nqqnlo = ',nqqnlo
      print*,' include two loop gg (1) not (0): nggnnlo = ',nggnnlo
      print*,' include two loop qq (1) not (0): nqqnnlo = ',nqqnnlo
      print*,' scale choice: (0) m; (1) sqrt(m^2+pt^2) = ',sclset
      print*,' Factorization muF/scale = kF = ',mfac
      print*,' Renormalization muR/scale = kR = ',kR
      print*,' Mass of the produced heavy quark: mQ = ',mq
      print*,' center of mass energy: sqrt(S) = ', scmsroot
      print*,' the number of light flavors: nlf = ', nlf
      print*,' the heavy-quark produced is ',q
      print*,' Npt, pTmin, pTmax',Npt,pt2rtmin,pt2rtmax
      print*,' Nyt, ymin, ymax',Nyt,yyrapmin,yyrapmax
      print*,' the cutoff parameter delta =',delta
            print*, 'choice : ', choice

      print*,'======================================='
      print*,'======================================='
      print*,'======================================='

c group structure constants ca = N of SU(N) = 3                                                                                   
      CA = 3.d0
      CF = 4.d0/3.d0
      CO = 24.d0
      CK = 8.d0/3.d0
      CQED = 80.d0/9.d0
      TF = 0.5d0
c put nf to number of light flavours                                                                                              
      nf = nlf

      m2 = mq*mq
      scms = scmsroot**2
      if (scms .le. 4.d0*m2) then
         print*,'variables out of range'
c         goto 39
      endif

c     VEGAS parameters
      ptsi = 30000
      itsi = 10
      epsi = 1d-20
      dimi = 2

      ptse = 10000
      itse = 10
      epse = 1d-20
      dime = 2

      ptsb = 10000
      itsb = 10
      epsb = 1d-20
      dimb = 2

      scale = 2d0*dsqrt(m2)

!     Scale setting
      muR = scale*kR
      mu2 =  muR*muR

C     We set the factorization scale           
      muF=mfac*scale

c     implement the various loops
      shad = 0d0
      shaderr = 0d0
c     choice = 1: partonic cross section in terms of sqrt(s)
c     choice = 2: partonic cross section in terms of beta
c     choice = 3: hadronic total  cross section
c     choice = 4: Luminosity i terms of sqrt(s)
c     choice = 5: Luminosity i terms of beta
c     choice = 6: saturation distribution smax
c     choice = 7: several scms values at once
c     choice = 5

      
      
      if (choice .eq. 1 .or. choice .eq. 4 .or. choice .eq. 6) then
         N_Ehat = 100       
         wdth_Ehat = (dsqrt(scms)-4d0*mq)/dble(N_Ehat)
         Do i = 0,N_Ehat -1
            Ehat = 4d0*mq + wdth_Ehat*(dble(i) + 0.5)
            print*, 'Ehat = ', Ehat
            
            if (Ehat < 4d0*mq) then
               goto 11
            endif

!     Scale setting
            muR = scale*kR
            mu2 =  muR*muR

C     We set the factorization scale           
            muF=mfac*scale

c     implement the various loops

            shadel   = 0.d0
            shdelrr= 0.d0


            call vegas(ff,epsb,dimb,ptsb,itsb,1,0)

            shadel = shadel + s1
            shdelrr = shdelrr + s2
 11         shad =  shadel
            shaderr =  shdelrr
            
            write(12,100) 4d0*mq + dble(i)*wdth_Ehat,Ehat, shad, shaderr

            print*,'Bin, Ehat [GeV], Result [pb], Err', 
     &           4d0*mq + i*wdth_Ehat,Ehat, shad, shaderr
         enddo
         
      endif

      if (choice .eq. 2 .or. choice .eq. 5) then
         N_beta = 20
         wdth_beta = 1d0/dble(N_beta)

         Do i = 0, N_beta-1
            beta = 0d0 + wdth_beta*(dble(i) + 0.5)
            print*, beta
         enddo
         Do i = 0, N_beta-1
            shadel = 0d0
            shdelrr = 0d0
            beta = 0d0
            beta = 0d0 + wdth_beta*(dble(i) + 0.5)
            print*, 'beta = ',beta
         
            shadel   = 0.d0
            shdelrr= 0.d0

            call vegas(ff,epsb,dimb,ptsb,itsb,1,0)

            shadel = shadel + s1
            shdelrr = shdelrr + s2
            shad =  shadel
            shaderr =  shdelrr
            write(13,100) dble(i)*wdth_beta, beta, shad, shaderr

            print*,'Bin, Beta, Result [pb], Err',i*wdth_beta,beta,
     &           shad, shaderr
         enddo
         
      endif
      
      if (choice .eq. 3) then
        
         call vegas(ff,epsb,dimb,ptsb,itsb,1,0)

         shadel   = 0.d0
         shdelrr= 0.d0
         
         shadel = shadel + s1
         shdelrr = shdelrr + s2
         shad =  shadel
         shaderr =  shdelrr
         
         write(14,100) shad, shaderr


         print*,'Energy [GeV], Result [pb], Err [pb]',
     &        dsqrt(scms), shad, shaderr
      endif

      if (choice .eq. 7) then
         N_E = 27
c     wdth_s = (sqrt(scms) - 1000d0)/(dble(N))
         wdth_E = 500d0
         Do i=0,N_E-1
            E_i = 1000d0 + wdth_E*(i)
            print*, i,E_i

            call vegas(ff,epsb,dimb,ptsb,itsb,1,0)

            shadel   = 0.d0
            shdelrr= 0.d0
            
            shadel = shadel + s1
            shdelrr = shdelrr + s2
            shad =  shadel
            shaderr =  shdelrr
            
            write(14,100) E_i,shad, shaderr


            print*,'Energy [GeV], Result [pb], Err [pb]',
     &           E_i, shad, shaderr
         enddo
      endif
      
      
 100  format (4x,4(e12.5,4x))
 999  stop
      end
