
      Subroutine InitialPdf

      implicit none

      integer iset
      double precision value
      character*20 parm(20)     
      dimension value(20)       

      common /pdfiset/iset
     
      if (iset.ge.10000 .and. iset.le.999999) then
c     call setlhaparm('SILENT') !uncomment to printout LHAPDF statements
      call setlhaparm('NOSTAT')
      
      parm(1)='DEFAULT'
      value(1)=iset
      call pdfset(parm,value)

      
      else 
         print *,
     >        'Illegal argument in pdfset(iset) call; iset=', iset
         Stop
      endif                    
      print*, 'initialpdf'
      return 
      end





      FUNCTION Parton (LPRTN, XD, QD)
      implicit none

      integer iset, iprtn, lprtn
      double precision xold, qold, f, tmp, parton, x, q, xd, qd
      Dimension f(-6:6)     ! LHAPDF

      Common /PdfIset/iset
      Data Xold, Qold /-1d0, -1d0 /
     
      Iprtn = LPRTN
      x = XD
      Q = QD
c      print*, iset, 'pac'
      If (Iset.ge.10000 .and. Iset.le.999999) then !LHAPDF
         If (X.ne.Xold .or. Q.ne.Qold) then
            Call evolvePDF(X,Q,f)
            Xold=X
            Qold=Q
         Endif
         If(Iprtn==1 .or. Iprtn==2) then
           Tmp=f(3-Iprtn)
         Elseif(Iprtn==-1 .or. Iprtn==-2) then
           Tmp=f(-3-Iprtn)
         Else
           Tmp=f(Iprtn)
         Endif
      Else
         print *,'Iset chosen is inactive. Parton set to 0'
         Tmp = 0.D0
       EndIf                    !iset

      Parton = Tmp
c      print*, 'PARTON',parton
      RETURN
C                        ****************************
      END
 







      DOUBLE PRECISION FUNCTION ALPHASLH(Q)
C 
      implicit none
      DOUBLE PRECISION Q, alphasPDF, tmp

      tmp = alphasPDF(Q)
      alphaslh=tmp

C
      RETURN
      END
C
