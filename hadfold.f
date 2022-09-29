c  This routine gives the value for alpha_s and the parton densities to be
c  convoluted with the partonic cross sections to give the structure
c  functions or cross sections.
*             (msbar), Lambda = 0.326 GeV(nf=4)
*             (msbar), Lambda = 0.226 GeV(nf=5) 3 loop
*             (msbar), Lambda = 0.186 GeV(nf=5) 2 loop
*             (msbar), Lambda = 0.0878 GeV(nf=5) 1 loop
      subroutine hadfold(xa,xb)
      implicit none

      integer ncount
      integer ipdf, iset
      integer PDFflag
      integer nalfas
      integer sclset         ! Flag for the scale choice (0) scale = m; (1) scale = sqrt(m^2+pt^2)

      real*8 xa, xb
      real*8 scale, alphas, mu2, nlf
      real*8 pi
      real*8 gluexa, gluexb
      real*8 ggflux, sumqqbflux, sumqbqflux, sumqgflux, sumgqflux,
     #     sumqbgflux, sumgqbflux
      real*8 uxa, ubxa, dxa, dbxa,
     #     cxa, cbxa, sxa, sbxa, bxa, bbxa
      real*8 uxb, ubxb, dxb, dbxb,
     #     cxb, cbxb, sxb, sbxb, bxb, bbxb
      real*8 Ctq5Pdf
      real*8 alphasLH, parton
      real*8 qcdl1,qcdl2, qcdl3
      real*8 mq,muR,kR,muF


      character*10 q
      parameter( pi = 3.14159265359d0 )

      common/coupling/scale, alphas, nlf
      common/sclchoice/mq, muR, kR, muF, sclset
      common/fluxes/ggflux, sumqqbflux, sumqbqflux,
     #       sumqgflux, sumgqflux, sumqbgflux, sumgqbflux
      common/quark/q
      common/pdfdb/qcdl1,qcdl2,qcdl3,nalfas
      common/PDFname/PDFflag
      common/pdfiset/iset
      common/pdf/ipdf

      external alphasLH
      external parton


      data ncount/0/
      save ncount


!     Scale setting
      muR = scale*kR
      mu2 =  muR*muR
c      print*, 'entered hadfold, shit'

      if (mu2.le.1d0/16d0*mq**2.or.mu2.ge.16d0*mq**2) then
        print*,'ERROR: Extreme choice of muR variation!'
        STOP
      endif

!  calculation of alpha_s by using Lambda QCD,
!  only used for making tests.
!  Checked vs RunDec by K.G. Chetyrkin, Johann H. Kuhn, M. Steinhauser, Comput.Phys.Commun. 133 (2000) 43-65

      if (nalfas.eq.1) then
        alphas = 12.d0*pi/( (33.d0-2.d0*nlf)*dlog(mu2/qcdl1**2) )
      elseif (nalfas.eq.2) then
        alphas = 12.d0*pi/( (33.d0-2.d0*nlf)*dlog(mu2/qcdl2**2) )
     >    *( 1.d0 - 6.d0*(153.d0-19.d0*nlf)/(33.d0-2.d0*nlf)**2
     >    *dlog(dlog(mu2/qcdl2**2))
     >    /dlog(mu2/qcdl2**2)               )
      elseif (nalfas.eq.3) then
        alphas =
     >    4*0.3141592653589793D1*(1/(11-2.D0/3.D0*nlf)/dlog(mu2/qcdl3**
     >    2)-(102-38.D0/3.D0*nlf)/(11-2.D0/3.D0*nlf)**3*dlog(dlog(mu2
     >    /qcdl3**2))/dlog(mu2/qcdl3**2)**2+1/((11-2.D0/3.D0*nlf)**5)
     >    /dlog(mu2/qcdl3**2)**3*((102-38.D0/3.D0*nlf)**2*dlog(dlog(mu2
     >    /qcdl3**2))**2-(102-38.D0/3.D0*nlf)**2*dlog(dlog(mu2/qcdl3**2)
     >    )+(2857.D0/2.D0-5033.D0/18.D0*nlf+325.D0/54.D0*nlf**2)*(11-2
     >    .D0/3.D0*nlf)-(102-38.D0/3.D0*nlf)**2))

! calculation of alpha_s by LHAPDF
      elseif (nalfas.eq.99) then

! scale setting
        muR = kR*scale

        if (muR.le.1d0/4d0*mq.or.muR.ge.4d0*mq) then
          print*,'ERROR: Extreme choice of muR variation!'
          STOP
        endif

        if(PDFflag .eq. 1) then

          alphas=alphasLH(muR)

        elseif(PDFflag .eq. 1 .and. nalfas .ne. 99) then
          print*,'ERROR: Wrong assignment of PDFflag, ipdf and nalfas'
          STOP
        endif

      endif




      if (ncount .eq. 0 .and. PDFflag .eq. 0) then
         print*,'alphas = ',alphas
         if (nalfas .eq. 2) then
            print*,'qcdl2 = ',qcdl2
         else
            print*,'qcdl3 = ',qcdl3
         endif
         print*,'mu2 = ',mu2
         print*,'scale = ',scale


         if (ipdf .eq. 1) then
            print*,'ERROR: Stand alon module not active'
            stop
         else
            print*,'ERROR: improper value of ipdf'
            stop
         endif
      elseif(ncount .eq. 0 .and. PDFflag .eq. 1) then
         print*,'alphas value is read from the LHAPDF inteface:'
         print*,'Q[GeV], alphas(Q) = ', muR, alphas

      endif


!  Parton densities from LHAPDF
      if (PDFflag.eq.1.and.ipdf.eq.0) then

         gluexa = parton(0,xa,muF)/xa
         gluexb = parton(0,xb,muF)/xb


         uxa  = parton(1,xa,muF)/xa
         ubxa = parton(-1,xa,muF)/xa
         dxa  = parton(2,xa,muF)/xa
         dbxa = parton(-2,xa,muF)/xa
         sxa  = parton(3,xa,muF)/xa
         sbxa = parton(-3,xa,muF)/xa
         cxa  = parton(4,xa,muF)/xa
         cbxa = parton(-4,xa,muF)/xa
         bxa  = parton(5,xa,muF)/xa
         bbxa = parton(-5,xa,muF)/xa

         uxb  = parton(1,xb,muF)/xb
         ubxb = parton(-1,xb,muF)/xb
         dxb  = parton(2,xb,muF)/xb
         dbxb = parton(-2,xb,muF)/xb
         sxb  = parton(3,xb,muF)/xb
         sbxb = parton(-3,xb,muF)/xb
         cxb  = parton(4,xb,muF)/xb
         cbxb = parton(-4,xb,muF)/xb
         bxb  = parton(5,xb,muF)/xb
         bbxb = parton(-5,xb,muF)/xb

      elseif(PDFflag.eq.1.and.nalfas.ne.99.or.
     ~        PDFflag.eq.1.and.ipdf.ne.0) then
         print*,'ERROR: wrong assignment of PDFflag, ipdf and nalfas'
         STOP 'LHAPDF interface requires PDFflag=1, ipdf=0, nalfas=99'
       endif


!  Stand alone user parton density functions
! can be implemented here

      if (PDFflag.eq.0.and.(ipdf .eq. 1)) then
        print*,'ERROR:Stand Alone PDFs Module is not active!'
        STOP
C         gluexa = stdaPARTON(0,xa,muF)/xa
C         gluexb = stdaPARTON(0,xb,muF)/xb

C         uxa  = stdaPARTON(1,xa,muF)/xa
C         ubxa = stdaPARTON(-1,xa,muF)/xa
C         dxa  = stdaPARTON(2,xa,muF)/xa
C         dbxa = stdaPARTON(-2,xa,muF)/xa
C         sxa  = stdaPARTON(3,xa,muF)/xa
C         sbxa = stdaPARTON(-3,xa,muF)/xa
C         cxa  = stdaPARTON(4,xa,muF)/xa
C         cbxa = stdaPARTON(-4,xa,muF)/xa
C         bxa  = stdaPARTON(5,xa,muF)/xa
C         bbxa = stdaPARTON(-5,xa,muF)/xa

C         uxb  = stdaPARTON(1,xb,muF)/xb
C         ubxb = stdaPARTON(-1,xb,muF)/xb
C         dxb  = stdaPARTON(2,xb,muF)/xb
C         dbxb = stdaPARTON(-2,xb,muF)/xb
C         sxb  = stdaPARTON(3,xb,muF)/xb
C         sbxb = stdaPARTON(-3,xb,muF)/xb
C         cxb  = stdaPARTON(4,xb,muF)/xb
C         cbxb = stdaPARTON(-4,xb,muF)/xb
C         bxb  = stdaPARTON(5,xb,muF)/xb
C         bbxb = stdaPARTON(-5,xb,muF)/xb

      elseif(PDFflag .eq. 1 .and. (ipdf .eq. 1)
     ~        .and. nalfas .ne. 99) then
         print*,'ERROR: wrong assignment of PDFflag, ipdf and nalfas'
         STOP
      endif



C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C PP and PPbar gg, qq and qg fluxes
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ggflux = gluexa*gluexb


c qq channel; pp collider parton flux
c depending upon the light flavors we choose the light
c quark densities for bottom, top.
      if (q .eq. 'top@LHC') then
         sumqqbflux =  uxa*ubxb + dxa*dbxb + sxa*sbxb + cxa*cbxb +
     #        bxa*bbxb
c         sumqqbflux =  uxa*ubxb
         sumqbqflux =  ubxa*uxb + dbxa*dxb + sbxa*sxb + cbxa*cxb +
     #        bbxa*bxb
      endif
c      if (q .eq. 'bottom@LHC') then
c         sumqqflux = uxa*ubxb + dxa*dbxb + sxa*sbxb + cxa*cbxb +
c     #        ubxa*uxb + dbxa*dxb + sbxa*sxb + cbxa*cxb
c      endif
c      if (q .eq. 'charm@LHC') then
c         sumqqflux = uxa*ubxb + dxa*dbxb + sxa*sbxb +
c     #        ubxa*uxb + dbxa*dxb + sbxa*sxb
c      endif

c qg channel; pp collider parton flux; note f_q/p = f_qb/pb
c all combinations added up qa.gb + ga.qbb + qba.gb + ga.qb
      if (q .eq. 'top@LHC') then
c         sumqgflux = (uxa + dxa + sxa + cxa + bxa)*gluexb +
c     #        gluexa*(ubxb + dbxb + sbxb + cbxb + bbxb) +
c     #        (ubxa + dbxa + sbxa + cbxa + bbxa)*gluexb +
c     #        gluexa*(uxb + dxb + sxb + cxb + bxb)
         sumqgflux = (uxa + dxa + sxa + cxa + bxa)*gluexb
         sumgqflux = gluexa*(uxb + dxb + sxb + cxb + bxb)
         sumqbgflux = (ubxa + dbxa + sbxa + cbxa + bbxa)*gluexb
         sumgqbflux = gluexa*(ubxb + dbxb + sbxb + cbxb + bbxb)

      endif
      if (q .eq. 'bottom@LHC') then
         sumqgflux = (uxa + dxa + sxa + cxa)*gluexb +
     #        gluexa*(ubxb + dbxb + sbxb + cbxb) +
     #        (ubxa + dbxa + sbxa + cbxa)*gluexb +
     #        gluexa*(uxb + dxb + sxb + cxb)
      endif
      if (q .eq. 'charm@LHC') then

         sumqgflux = (uxa + dxa + sxa)*gluexb +
     #        gluexa*(ubxb + dbxb + sbxb) +
     #        (ubxa + dbxa + sbxa)*gluexb +
     #        gluexa*(uxb + dxb + sxb)
      endif


c qq channel; ppbar collider parton flux; note f_q/p = f_qb/pb
c depending upon the light flavors we choose the light
c quark densities for bottom, top.
c      if (q .eq. 'bottom@TEV') then
c         sumqqflux = uxa*uxb + dxa*dxb + sxa*sxb + cxa*cxb +
c     #        ubxa*ubxb + dbxa*dbxb + sbxa*sbxb + cbxa*cbxb
c      endif
c      if (q .eq. 'top@TEV') then
c         sumqqflux =  uxa*uxb + dxa*dxb + sxa*sxb + cxa*cxb + bxa*bxb +
c     #        ubxa*ubxb + dbxa*dbxb + sbxa*sbxb + cbxa*cbxb + bbxa*bbxb
c      endif

c qg channel; ppbar collider parton flux; note f_q/p = f_qb/pb
c all combinations added up qa.gb + ga.qbb + qba.gb + ga.qb
      if (q .eq. 'bottom@TEV') then
         sumqgflux = (uxa + dxa + sxa + cxa)*gluexb +
     #        gluexa*(ubxb + dbxb + sbxb + cbxb) +
     #        (ubxa + dbxa + sbxa + cbxa)*gluexb +
     #        gluexa*(uxb + dxb + sxb + cxb)
      endif
      if (q .eq. 'top@TEV') then
         sumqgflux = (uxa + dxa + sxa + cxa + bxa)*gluexb +
     #        gluexa*(ubxb + dbxb + sbxb + cbxb + bbxb) +
     #        (ubxa + dbxa + sbxa + cbxa + bbxa)*gluexb +
     #        gluexa*(uxb + dxb + sxb + cxb + bxb)
      endif









      ncount = 1
      return
      end
