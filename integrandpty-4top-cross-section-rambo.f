      double precision function ff(yy)
c-----------------------------------------------------
      implicit double precision(a-z)
c     Integrand."yy" is 10-dimensional array
      dimension yy(10)

      real*8 yy
      real*8 m2, scms, beta, sigmaqq
      real*8 mu2, scale, pi, alphas, nlf 
      real*8 ggflux, sumqqbflux, sumqbqflux, sumqgflux, sumgqflux,
     #sumqbgflux, sumgqbflux
      real*8 mq,muR,kR,muF          
      real*8 xm,p,wt,tau, z
      real*8 gs,gs2,gs4,k,l,phi,xpar, me_qq, me_qQ2tTtT, Vol
      
      integer i,m,j,mpart,iter,N,ii,itsb,choice
      
      parameter(mpart = 4)
      parameter (pi = 3.14159265359d0)

      dimension xm(mpart),p(4,mpart),k(4,2), xpar(2)
      character*10 q
      
      common/group/ca, cf, co, ck, cqed, tf, nf
      common/m2only/m2
      common/sonly/scms
      common/coupling/scale, alphas, nlf
      common/sclchoice/mq, muR, kR, muF, sclset 
      common/fluxes/ggflux, sumqqbflux, sumqbqflux,
     #sumqgflux, sumgqflux, sumqbgflux, sumgqbflux
      common/Ehatbeta/Ehat,beta,E_i
      common/choice/choice
      external me_qQ2tTtT, me_gg2tTtT
      


c     Ehat distribution
      if (choice .eq. 1 .or. choice .eq. 4) then
         shat = Ehat**2d0
      endif
      if (choice .eq. 2 .or. choice .eq. 5) then
         shat = dble(mpart)**2d0*m2/(1d0 - beta**2d0)
      endif
      if (choice .eq. 4 .or. choice .eq. 5) then
         stilde = shat + (scms - shat)*yy(1)
      endif
      if (choice .eq. 3) then
         rho = dble(mpart)**2d0*m2/scms
      elseif (choice .eq. 7) then
         rho = dble(mpart)**2d0*m2/(E_i**2d0)
      endif
      if (choice .eq. 3 .or. choice .eq. 7) then
c     csm momentum fractions x1 and x2
         
         xpar(1)=rho+(1d0-rho)*yy(1)
         xpar(2)=-dlog(xpar(1))*yy(2)

         tau = xpar(1)
         z = xpar(2)
         
      endif
      if (choice .eq. 3) then
         shat = tau*scms
      elseif(choice .eq. 7) then
         shat = tau*E_i**2d0

      endif
      
      if (choice .eq. 6) then
         smax = Ehat**2d0
         shat = dble(mpart)**2d0*m2 + (smax - dble(mpart)**2d0*m2)*yy(1)
         stilde = shat + (scms - shat)*yy(2)
      endif

      ff = 0d0
      
      if (choice .eq. 1 .or. choice .eq. 2 .or. choice .eq. 3
     &     .or. choice .eq. 6 .or. choice .eq. 7) then
         
         gs2 = 4*pi*alphas
         gs = dsqrt(gs2)
         gs4 = gs2**2d0
         
         k(1,1)=dsqrt(shat)/2
         k(2,1)=0d0
         k(3,1)=0d0
         k(4,1)=k(1,1)

         k(1,2)=-dsqrt(shat)/2
         k(2,2)=0d0
         k(3,2)=0d0
         k(4,2)=-k(1,2)

         N=mpart
         Do i=1,N
            xm(i)=mq
         Enddo

c     call RAMBO
         
         call rambo(N,dsqrt(shat),xm,p,wt,0)
         
         wt = wt*1d0/((2d0*pi)**(dble(N)-2d0))
     &        *1d0/((2d0*pi)**(2d0*dble(N)-2d0))


c     compute invariants
         t11 = m2 - 2d0*dotprod(k,p,1,1)
         t21 = m2 - 2d0*dotprod(k,p,2,1)
         t12 = m2 - 2d0*dotprod(k,p,1,2)
         t22 = m2 - 2d0*dotprod(k,p,2,2)
         t13 = m2 - 2d0*dotprod(k,p,1,3)
         t23 = m2 - 2d0*dotprod(k,p,2,3)
         s12 = 2d0*m2 + 2d0*dotprod(p,p,1,2)
         s13 = 2d0*m2 + 2d0*dotprod(p,p,1,3)
         s23 = 2d0*m2 + 2d0*dotprod(p,p,2,3)
c     consistency check:
         cnsstncy = t11 + t21 + t12 + t22 + t13 + t23 + s13 + s23 + shat
     &        +s12 -10d0*m2
         
         if (abs(cnsstncy) .le. 1e-5) then
            goto 11
         else
            print*, 'WARNING INVARIANTS MUST BE WRONG.'
            stop
         endif
      endif
      
c     Call hadfold
 11   if (choice .eq. 1 .or. choice .eq. 2) then
         call hadfold(0d0,0d0) 
      endif

      if (choice .eq. 3 .or. choice .eq. 7) then
         call hadfold(dexp(-z),tau*dexp(z))
      endif
      
      if (choice .eq. 4 .or. choice .eq. 5 .or. choice .eq.6) then
         call hadfold(stilde/scms,shat/stilde)
      endif
c     compute matrix elements
      if (choice .eq. 1 .or. choice .eq. 2
     &     .or. choice .eq. 3 .or. choice .eq.6 .or. choice .eq. 7) then
         
         me_qqb=me_qQ2tTtT(shat,t11,t21,t12,t22,s12,t13,t23,s13,s23,mq
     &        ,gs,3d0)/144d0

         me_qbq=me_qQ2tTtT(shat,t21,t11,t22,t12,s12,t23,t13,s13,s23,mq
     &        ,gs,3d0)/144d0

         me_gg=me_gg2tTtT(shat,t11,t21,t12,t22,s12,t13,t23,s13,s23,mq                                                              
     &        ,gs,3d0)/1024d0   
      endif
c     compute partonic cross sections
      if (choice .eq. 2 .or. choice .eq. 5) then
         goto 12
      else
         beta = dsqrt(1d0 - (dble(mpart)**2d0*m2)/shat)
         goto 12
      endif

 12   if (choice .eq. 1 .or. choice .eq. 2 .or. choice .eq. 3 .or.
     &     choice .eq. 6 .or. choice .eq. 7)then
         
c     partonic cross sections
         sigmaqqb = 1d0/(2d0*shat)*wt*me_qqb 
         sigmaqbq = 1d0/(2d0*shat)*wt*me_qbq
         sigmagg = 1d0/(2d0*shat)*wt*me_gg

      endif
      if (choice .eq. 1 .or. choice .eq. 2) then
         
         ff = ff + (
c     Quark channels
     &        sigmaqqb
     &        + sigmaqbq
     &        )
c     LL Correction qqb,qbq
c     &     *(1d0 + alphas/(pi)*8d0*cf*dlog(beta)**2d0)
c     gluon channel
c     &        + sigmagg
c     LL corrections
c     &        *(1d0 + alphas/(pi)*8d0*ca*dlog(beta)**2d0)
         jacobian = 1d0
         
      endif
      
      if (choice .eq. 3 .or. choice .eq. 7) then
         
c     only Born cross section
      ff = ff + (sumqqbflux*sigmaqqb
     &        + sumqbqflux*sigmaqbq
     &        + ggflux*sigmagg
     &        )
         
c     born + leading log
c        a=1.0d0

c         ff = ff + (sumqqbflux*sigmaqqb + sumqbqflux*sigmaqbq)
c     &        *(1d0 + alphas/(pi)*8d0*cf*dlog(beta)**2d0 - a*dlog(beta))
c     &     + ggflux*sigmagg*(1d0 + alphas/(pi)*8d0*ca*dlog(beta)**2d0)

         jacobian = -(1.d0-rho)*dlog(xpar(1))
      endif
      
      if (choice .eq. 4 .or. choice .eq. 5 .or.
     &     choice .eq. 6 ) then
         
c     qqb Luminosity
         Lqqb = (1d0/scms)*(1d0/stilde)*sumqqbflux
c     qbq Luminosity
         Lqbq = (1d0/scms)*(1d0/stilde)*sumqbqflux
c     gluon luminosity
         Lgg = (1d0/scms)*(1d0/stilde)*ggflux
      endif
      if (choice .eq.  4 .or. choice .eq. 5) then

         ff = ff
     &        + Lqqb
     &        + Lqbq
c     &        + Lgg
         
         jacobian = (scms - shat)

      endif

      if (choice .eq. 6) then

         ff = ff + (Lqqb*sigmaqqb
     &        + Lqbq*sigmaqbq
     &        + Lgg*sigmagg
     &        )

         jacobian = (smax - dble(mpart)**2d0*m2)*(scms - shat)
         
      endif
      
      
      if (choice .eq. 4 .or. choice .eq. 5)then
         ff = jacobian*ff
      else
c     conversion to express the cross section in pico barn for physical cross sections
         ff =  jacobian*19733.d0*19733.d0*ff
      endif




      return 
      end


cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     Minowski-skalar product with gamma_mu_nu=diag(-1,1,-1,+1)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      double precision function dotprod(k1,k2,l,m)
      Implicit none
      real*8 ans,k1,k2
      integer i,l,m
      dimension k1(4,l), k2(4,m)

      dotprod = -k1(1,l)*k2(1,m) - k1(2,l)*k2(2,m) - k1(3,l)*k2(3,m) +
     &     k1(4,l)*k2(4,m)
      
      return 
      end
