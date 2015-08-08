      program bbHW
c-----------------------------------------------------------
c  This program calculates the NLO and NNLO soft 
c  corrections to FCNC top quark + gamma production
c-----------------------------------------------------------
      implicit double precision(a-z)
      integer pts, its, norder,
     & iloop,nscale
      dimension mql(1:11)
      dimension q2l(1:20)
      dimension out(2,0:20)
      common/result/s1,s2,s3,s4
      common/cross/scap
      common/mass2/m2
      common/qcdpar/q2,qcdl4
      common/parone/norder
      common/partwo/nf
      common/pdfalfa/alfas
      common/pdfmstw/isetmstw

c
      parameter(pi=3.14159265359d0)
c                                                                     
      external f
c
      open(unit=7, file='output.dat',status='old')
      open(unit=8, file='outputmass.dat',position='append',status='old')

c
c      open(unit=10,file='topgugam2a1960muzetanew.dat',status='new')

      open(unit=11, file='bbhw.dat',status='old')
      rewind 11
c
      read(11,*)pts
      read(11,*)its
      read(11,*)norder
      read(11,*)nf
      read(11,*)srcaps
      read(11,*)mq
      read(11,*)qcdl4
      read(11,*)nscale
c
      print '(''number of points :'',i7)',pts
      print '(''number of iterations :'',i7)',its
      print '(''born (0), 1st ord.(1), 2nd ord(2):'',i2)',norder
      print '(''number of flavors :'',d20.7)',nf
      print '(''square root of s for p-pbar:'',d20.7)',srcaps
      print '(''charged higgs mass:'',d20.7)',mq
      print '(''qcdl4:'',d20.7)',qcdl4
      print '(''scale (1)m (0)m/2 (2)2m:'',i2)',nscale
c
      scap = srcaps*srcaps
c
c computation of the total cross section:
c   lines 56 [m2] to 63 [endif] were originally in the code!
c      do 510 iloop= 1,20
c      mql(iloop) = 50d0*dble(iloop-1.d0) + 500.d0
c      mq=mql(iloop)        
      m2 = mq*mq
      if (nscale.eq.0) then
         q2 = m2/4.d0
         elseif (nscale.eq.1) then 
            q2 = m2
            elseif (nscale.eq.2) then
               q2 = 4.d0*m2
       endif
c-------mu dependence----------
      m2 = mq*mq
      do 510 iloop= 1,20
       if (iloop .le. 10) then
        q2l(iloop) =dble(iloop)/10.d0 
         else
        q2l(iloop)=dble(iloop)-10.d0
       endif
      q2=m2*(q2l(iloop))**2
c------end of mu dep----------------
c The above, lines 64 - 74 [ten lines ish] were originally commented out!!
         mu=dsqrt(q2)
         r = mu/mq
c
         alphaSorder=2
         alphaSMZ=0.11707d0     
         mCharm=1.4d0
         mBottom=4.75d0
         CALL INITALPHAS(alphaSorder,1.D0,91.1876D0,alphaSMZ,
     &     mCharm,mBottom,1.D10)
          alfas=ALPHAS(mu)
         print'(''****** alfas = '',d20.8)',alfas
         isetmstw=0
c
         call vegas(f,1d-4,4,pts,its,1,0)
c change to pb gev^(-2) unit  from    mub gev^(-2)
        s1  = s1*1.d6
       sigppb = s1
c        out(1,iloop)=mql(iloop)
c        out(1,iloop)=q2l(iloop)
c        out(1,iloop)=E1
c        out(2,iloop)=s1
      print'(''the ppbar correction:'',d20.8)',sigppb
      print'(''****** mq = '',d20.8)',mq
       write(7,*) r, sigppb 
       write(8,*) mq, sigppb
 510   continue
c       write(10,1500) ((out(i,j),i=1,2),j=1,15)
c       write(10,1500) ((out(i,j),i=1,2),j=1,11)
c       write(10,1500) ((out(i,j),i=1,2),j=1,20)
c 1500  format(e15.4,5x,e15.4)

       goto 550
c
c
 550     continue
c
      stop
      end
c
c*********************************************************************
c
      double precision function f(x)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit double precision (a-z)
      integer norder
      dimension x(10)
      common/result/s1,s2,s3,s4
      common/cross/scap
      common/mass2/m2
      common/qcdpar/q2,qcdl4
      common/parone/norder
      common/partwo/nf
      common/jacob/jacob1,jacob2,jacob3,jacob4
      common/jacobhe/hjacob,ejacob
      common/t1u1/t1el,u1el,shel
      common/t1u1sh/t1,u1,sh
c
      parameter(pi=3.14159265359d0)
c
       SS = scap	 
       mu=dsqrt(q2)
c
c we do the integrals in the order T,U,xb,s2 
c
c Definition of the integration variables
c
      MW2 = 80.385d0*80.385d0
c This is a placeholder mass, clearly
c
      Tmax = m2-1.d0/2.d0*(SS+m2-MW2)
     & +1.d0/2.d0*dsqrt((SS+m2-MW2)**2-4.d0*m2*SS)
      Tmin =m2-1.d0/2.d0*(SS+m2-MW2)
     & -1.d0/2.d0*dsqrt((SS+m2-MW2)**2-4.d0*m2*SS) 
      TT=Tmin+x(1)*(Tmax-Tmin)
      jacob1=Tmax-Tmin
c
      Umax = m2+m2*SS/(TT-m2)
      Umin = -SS-TT+m2+MW2
      UU=Umin+x(2)*(Umax-Umin)
      jacob2=Umax-Umin
c
      xbmax = 1.d0
      xbmin = (MW2-TT)/(SS+UU-m2)
      xb=xbmin+x(3)*(xbmax-xbmin)
      jacob3=xbmax-xbmin
c
      s4max=xb*(SS+UU-m2)+TT-MW2
      s4min=0.d0 
      ss4=s4min+x(4)*(s4max-s4min)
      jacob4=s4max-s4min
c
      xa=(ss4-m2+MW2-xb*(UU-m2))/(xb*SS+TT-m2)
      xael=(-m2+MW2-xb*(UU-m2))/(xb*SS+TT-m2)
      xbel=xb

      t1 = xa*(TT-m2)
      t1el=xael*(TT-m2) 
      u1 = xb*(UU-m2)
      u1el=u1 
      sh = xa*xb*SS
      shel=xael*xbel*SS
c
      hjacob=xa*xb/(xb*SS+TT-m2)
      ejacob=xael*xbel/(xbel*SS+TT-m2)
c---------------------------------------------------------------------
c We call the elastic and the inelastic subroutines separately.
c However, because the soft plus virtual pieces have been
c put back into the inelastic piece, the inelastic
c routines (labelled by a '2' suffix), also
c use the elastic variable. 
c In particular, the parton elastic s 'shel', the
c parton flux 'sfpartel' and the jacobian 'ejacob' are
c transported through common blocks.
c
      if(norder.ge.1)then
         f=jacob1*jacob2*jacob3*ejacob
     & *aptel(s4max,m2,xael,xbel)
     &    +jacob1*jacob2*jacob3*jacob4
     & *aptinel(ss4,m2,xa,xb)
      else
         f=jacob1*jacob2*jacob3*ejacob
     & *aptel(s4max,m2,xael,xbel)
      endif
c      print '(''f='',d20.7)',f
c
      return
      end
c
c----------------------------------------------------------------------
c
      double precision function aptel(s4max,m2,xa,xb)
      implicit double precision(a-z)
      integer norder, iset, npp
      common/parone/norder      
      common/qcdpar/q2,qcdl4
      common/partwo/nf
      common/sfel/sfpartel
      common/pdfalfa/alfas
      common/pdfmstw/isetmstw

      parameter(pi=3.14159265359d0)
c
c Compute the parton densities
c
        scale=dsqrt(q2)
        scale=dsqrt(q2)
        iset=isetmstw
c
      G1   = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XA,scale,0)
      UT1  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XA,scale,5)
      UV1  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XA,scale,11)
c
      G1XA = G1/XA
      UTXA = UT1/XA
      USXA = (UT1-UV1)/XA
      UVXA = UV1/XA
c
      G2   = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XB,scale,0)
      UT2  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XB,scale,5)
      UV2  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XB,scale,11)
c Idk if this is correct since I need b-bar...
      G2XB = G2/XB
      UTXB = UT2/XB
      USXB = (UT2-UV2)/XB
      

c      for gu -> t gamma
      SFPARTel  = G1XA*USXB+(UVXA+USXA)*G2XB 
c      for gc -> t gamma
c       SFPARTel  = G1XA*CSXB+CSXA*G2XB 
c      for all: gu, gc, -> t gamma 
c      SFPARTel  = G1XA*(UVXB+2.D0*USXB+2.D0*CSXB)
c     & +(UVXA+2.D0*USXA+2.D0*CSXA)*G2XB

c
c      print '(''SFPARTel='',d20.7)',SFPARTel
c
c the computation of the parton cross section for uu:
c  the flag norder selects the order:
c  born (0), 1st order (1), 2nd order (2)
      if(norder.eq.0)then
          call siqqborn(s5,m2)
          sigtau = s5
       elseif(norder.eq.1)then
          call siqqone1(s6,s4max,m2)
          sigtau = s6
       elseif(norder.eq.2)then
          call siqqtwo1(s7,s4max,m2)
          sigtau = s7
       endif
c
c ----------------------------------------------------
          aptel=sfpartel*sigtau
c      print '(''aptel='',d20.7)',aptel
c
        return
        end
c**********************************************************
c
      double precision function aptinel(ss4,m2,xa,xb)
      implicit double precision(a-z)
      integer norder
      common/parone/norder      
      common/qcdpar/q2,qcdl4
      common/partwo/nf
      common/sfinel/sfpartinel
      common/pdfalfa/alfas
      common/pdfmstw/isetmstw

      parameter(pi=3.14159265359d0)
c
c Compute the parton densities
        scale=dsqrt(q2)
c Commented out the below!
c      CALL MRST2(XA,q2,UV1,DV1,USEA1,DSEA1,STR1,CHM1,BTM1,G1)
c
c      G1XA = G1/XA
c      UVXA = UV1/XA
c      USXA = USEA1/XA
c      CSXA = CHM1/XA 
c
c      CALL MRST2(XB,q2,UV2,DV2,USEA2,DSEA2,STR2,CHM2,BTM2,G2)
c
c      G2XB = G2/XB
c      UVXB = UV2/XB
c      USXB = USEA2/XB
c      CSXB = CHM2/XB

c commenting the above

c      CHARACTER prefix*50,prefix1*55,flavours(-5:5)*10,
c     &     xfilename*50,qfilename*50
c 
c      prefix = "Grids/mstw2008nnlo.90cl" ! prefix for the grid files
c
c Compute the parton densities
        scale=dsqrt(q2)
        iset=isetmstw
c
      G1   = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XA,scale,0)
      UT1  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XA,scale,5)
      UV1  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XA,scale,11)
c
      G1XA = G1/XA
      UTXA = UT1/XA
      USXA = (UT1-UV1)/XA
      UVXA = UV1/XA
c
      G2   = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XB,scale,0)
      UT2  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XB,scale,5)
      UV2  = GetOnePDF("../../Grids/mstw2008nnlo.90cl",iset,XB,scale,11)
c
      G2XB = G2/XB
      UTXB = UT2/XB
      USXB = (UT2-UV2)/XB
      
c stuff above just copied.


c
c      for gu -> t gamma
      SFPARTinel  = G1XA*USXB+(UVXA+USXA)*G2XB 
c      for gc -> t gamma
c       SFPARTinel  = G1XA*CSXB+CSXA*G2XB
c      for all: gu, gc, -> t gamma
c      SFPARTinel  = G1XA*(UVXB+2.D0*USXB+2.d0*CSXB)
c     & +(UVXA+2.D0*USXA+2.D0*CSXA)*G2XB
c
c
c the computation of the parton cross section for uu:
c  the flag norder selects the order:
c  born (0), 1st order (1), 2nd order (2)
      if(norder.eq.0)then
          call siqqborn(s5,m2)
          sigtau = s5
       elseif(norder.eq.1)then
          call siqqone2(s6,ss4,m2)
          sigtau = s6
       elseif(norder.eq.2)then
          call siqqtwo2(s7,ss4,m2)
          sigtau = s7
       endif
c
c ----------------------------------------------------
c In contrast to the elastic aptel, the parton flux
c sfpartinel is used in the subroutines that make
c the double differential distributions.
          aptinel=sigtau
c      print '(''aptinel='',d20.7)',aptinel
c
        return
        end
c-----------------------------------------------------------------------
c
       subroutine siqqborn(s5,m2)
c_________________________________________________________________
cIn this subroutine the Born cross section is computed
c____________________________________________________________
      implicit double precision(a-z)
      common/t1u1/t1el,u1el,shel
c
      parameter(pi=3.14159265359d0)       
c
      sigb0=qqdtdu(shel,m2,t1el,u1el)
c
c      print '(''sigb0='',d20.7)',sigb0
c
      s5=sigb0
c
      return
      end
c

c________________________________________________________________________
c
      subroutine siqqone1(s6,s4max,m2)
c__________________________________________________________________
c In this subroutine the first order correction to the
c cross section is computed in the MSbar scheme, elastic piece
c___________________________________________________________________
      implicit double precision(a-z)
      common/qcdpar/q2,qcdl4
      common/parsqq/coef,cqq
      common/jacob/jacob1,jacob2,jacob3,jacob4
      common/pdfalfa/alfas
      common/t1u1/t1el,u1el,shel
      common/partwo/nf
c
      parameter(pi=3.14159265359d0)
c
      srs=dsqrt(shel)
      srscap=dsqrt(scap)
      mu = dsqrt(q2)
c
      lns4max=dlog(s4max/m2)
c Born piece
      sigb0=qqdtdu(shel,m2,t1el,u1el)
c 
      ca=3.d0
      cf=4.d0/3.d0
      betaz=(11.d0*ca-2.d0*nf)/3.d0
      MW2=80.385d0**2
c
      zeta2=pi*pi/6.d0
c
       tel = t1el + m2
c temporary
c       tel =ABS(tel)
       uel = u1el + m2
c
       twel=tel-MW2
       uwel=uel-MW2
c
       mq=dsqrt(m2)
c       ReGSel=cf*dlog(-u1el/mq/dsqrt(shel))
c     & +ca/2.d0*dlog(t1el/u1el)+ca/2.d0 
c c3 is good
       c3m=4.d0*cf
       c2mel=-2.d0*cf*dlog(mu**2/m2)-2.d0*cf*dlog(m2/shel)
     & -2.d0*cf*dlog(-uwel/m2)-2.d0*dlog(-twel/m2)
c
c The minus sign should be checked
c
c The NNLL MSbar \delta(s_4) term
c
      c1mel=(dlog(-twel/m2)+dlog(-uwel/m2)-3.d0/2.d0)
     & *dlog(mu**2/m2)*cf
c
c toying with the code
c      c1mel=0
c      c2mel=0
c MSbar NLL
      f1=(c3m*lns4max**2/2.d0+c2mel*lns4max+c1mel)
     & *alfas/pi*sigb0         
      print '(''f1='',d20.7)',f1

c
c Computation of  surface term:
      s6=f1
c
      RETURN 
      end
c
c________________________________________________________________________
c
      subroutine siqqone2(s6,ss4,m2)
c__________________________________________________________________
c In this subroutine the first order correction to the 
c cross section is computed in the MSbar scheme, inelastic piece
c___________________________________________________________________
      implicit double precision(a-z)
      common/qcdpar/q2,qcdl4
      common/parsqq/coef,cqq
      common/jacob/jacob1,jacob2,jacob3,jacob4
      common/jacobhe/hjacob,ejacob
      common/sfel/sfpartel
      common/sfinel/sfpartinel
      common/pdfalfa/alfas
      common/t1u1/t1el,u1el,shel
      common/t1u1sh/t1,u1,sh
      common/partwo/nf
c
      parameter(pi=3.14159265359d0)
c
      srs=dsqrt(sh)
      srscap=dsqrt(scap)
      mu = dsqrt(q2)
c
      lns4=dlog(ss4/m2)

      ca=3.d0
      cf=4.d0/3.d0
      betaz=(11.d0*ca-2.d0*nf)/3.d0
c Born piece
      sigbs4=qqdtdu(sh,m2,t1,u1)
      sigb0=qqdtdu(shel,m2,t1el,u1el)
c
       MW2=80.385d0**2
c
       tel = t1el + m2
       uel = u1el + m2       
       tinel = t1 + m2
       uinel = u1 + m2        
c
       twel=tel-MW2
       uwel=uel-MW2
       twinel=tinel-MW2
       uwinel=uinel-MW2
c
       mq=dsqrt(m2)
c       ReGSel=cf*dlog(-u1el/mq/dsqrt(shel))
c     & +ca/2.d0*dlog(t1el/u1el)+ca/2.d0 
c       ReGSin=cf*dlog(-u1/mq/dsqrt(sh))
c     & +ca/2.d0*dlog(t1/u1)+ca/2.d0
c
       c3m=4.d0*cf
       c2mel=-2.d0*cf*dlog(mu**2/m2)-2.d0*cf*dlog(m2/shel)
     & -2.d0*cf*dlog(-uwel/m2)-2.d0*dlog(-twel/m2)
       c2minel=-2.d0*cf*dlog(mu**2/m2)-2.d0*cf*dlog(m2/sh)
     & -2.d0*cf*dlog(-uwinel/m2)-2.d0*dlog(-twinel/m2)
cc
c
c       c3m=2.d0*(cf+ca)
c       c2mel=2.d0*ReGSel-cf-ca-2.d0*cf*dlog(-tel/m2)
c     & -2.d0*ca*dlog(-uel/m2)-(cf+ca)*dlog(mu**2/shel)
c       c2minel=2.d0*ReGSin-cf-ca-2.d0*cf*dlog(-tinel/m2)
c     & -2.d0*ca*dlog(-uinel/m2)-(cf+ca)*dlog(mu**2/sh)
c
c testing impact of c3m vs c2 etc
c       c2minel=0
c       c2mel=0
c
c MSbar NLL
      df1ds4inel=(c3m*lns4/ss4+c2minel/ss4)
     & *alfas/pi*sigbs4
c
      df1ds4el=(c3m*lns4/ss4+c2mel/ss4)
     & *alfas/pi*sigb0
c
cnk      print '(''df1ds4='',d20.7)',df1ds4
c Computation of integral and surface term:
        int1=df1ds4inel*sfpartinel*hjacob
     &   -df1ds4el*sfpartel*ejacob
cnk      print '(''sfpartinel='',d20.7)',sfpartinel
cnk      print '(''sfpartel='',d20.7)',sfpartel
c      print '(''int1='',d20.7)',int1
c
cnk      print '(''coef='',d20.7)',coef
      s6=int1
c
      RETURN 
      end
c
c_______________________________________________________________________
c
      subroutine siqqtwo1(s7,s4max,m2)
c----------------------------------------------------------------------
c In this subroutine the second order correction to the cross
c is computed in the MSbar scheme, elastic piece
c___________________________________________________________________
      implicit double precision (a-z)
      common/qcdpar/q2,qcdl4
      common/partwo/nf
      common/parsqq/coef,cqq
      common/jacob/jacob1,jacob2,jacob3,jacob4
      common/pdfalfa/alfas
      common/t1u1/t1el,u1el,shel
c
      parameter(pi=3.14159265359d0)
c
      mu = dsqrt(q2)
c
      lns4max=dlog(s4max/m2)
c Born piece
      sigb0=qqdtdu(shel,m2,t1el,u1el)
c
      cf=4.d0/3.d0 
      ca=3.d0
      betaz=(11.d0*ca-2.d0*nf)/3.d0
c
      betaon=34.d0/3.d0*ca**2-2.d0*(5.d0/3.d0*ca+cf)*nf
c
      zeta2=pi*pi/6.d0
      zeta3=1.2020569031d0
      zeta4=1.0823232337d0
      MW2=80.385d0**2
c
       tel = t1el + m2
       uel = u1el + m2
c
       twel=tel-MW2
       uwel=uel-MW2
c
       mq=dsqrt(m2)
c       ReGSel=cf*dlog(-u1el/mq/dsqrt(shel))
c     & +ca/2.d0*dlog(t1el/u1el)+ca/2.d0 
c
       c3m=4.d0*cf
       c2mel=-2.d0*cf*dlog(mu**2/m2)
       T2mel=-2.d0*cf*(dlog(-uwel/m2) + dlog(-twel/m2) + dlog(m2/shel))
       c2el=c2mel+T2mel
c
c The NNLL MSbar \delta(s_4) term
c
      c1mel=(dlog(-twel/m2)+dlog(-uwel/m2)-3.d0/2.d0)
     & *dlog(mu**2/m2)*cf
c      c1mel=(cf*dlog(-tel/m2)+ca*dlog(-uel/m2)-3.d0/4.d0*cf)
c     & *dlog(mu**2/m2)
c
c MSbar NLL
      f2=c3m**2/2.d0*lns4max**4/4.d0
     & +(3.d0/2.d0*c3m*c2el-betaz/4.d0*c3m)
     & *lns4max**3/3.d0
c scale+zeta NNLL terms
     & +lns4max**2/2.d0*(c3m*c1mel + c2mel**2 + betaz*c3m/4.0d0
     & *dlog(mu**2/m2) + 2.d0*c2mel*T2mel
     & -zeta2*c3m**2)
c log^2(scale)+zeta 1/s4 terms
     & + lns4max*(c2mel*c1mel - zeta2*c3m*c2el
     & + betaz/4.d0*c2mel*dlog(mu**2/m2)
     & + cf*betaz/4.d0*(dlog(mu**2/m2))**2)
c Computation of integral and surface term:
c
c      print '(''qqdtdu='',d20.7)',f2
      s7 =coef*coef*f2*sigb0
C
      RETURN 
      end
c
c_______________________________________________________________________
c
      subroutine siqqtwo2(s7,ss4,m2)
c----------------------------------------------------------------------
c In this subroutine the second order correction to the cross
c is computed in the MSbar scheme, inelastic piece
c___________________________________________________________________
      implicit double precision (a-z)
      common/qcdpar/q2,qcdl4
      common/partwo/nf
      common/parsqq/coef,cqq
      common/jacob/jacob1,jacob2,jacob3,jacob4
      common/jacobhe/hjacob,ejacob
      common/sfel/sfpartel
      common/sfinel/sfpartinel
      common/pdfalfa/alfas
      common/t1u1/t1el,u1el,shel
      common/t1u1sh/t1,u1,sh
c
      parameter(pi=3.14159265359d0)
c
      mu = dsqrt(q2)
c
      lns4=dlog(ss4/m2)
c
c Born piece
      sigbs4=qqdtdu(sh,m2,t1,u1)
      sigb0=qqdtdu(shel,m2,t1el,u1el)
c
      cf=4.d0/3.d0
      ca=3.d0
      betaz=(11.d0*ca-2.d0*nf)/3.d0
C Derivative of f2 w.r.t. s4 (DIS)
c
      zeta2=pi*pi/6.d0 
      zeta3=1.2020569031d0
c
c       tinel = t1 + m2
c       uinel = u1 + m2
c       tel = t1el + m2
c       uel = u1el + m2
c
       MW2=80.385d0**2
c
       tel = t1el + m2
       uel = u1el + m2       
       tinel = t1 + m2
       uinel = u1 + m2        
c
       twel=tel-MW2
       uwel=uel-MW2
       twinel=tinel-MW2
       uwinel=uinel-MW2
c
       mq=dsqrt(m2)
c       ReGSel=cf*dlog(-u1el/mq/dsqrt(shel))
c     & +ca/2.d0*dlog(t1el/u1el)+ca/2.d0 
c       ReGSin=cf*dlog(-u1/mq/dsqrt(sh))
c     & +ca/2.d0*dlog(t1/u1)+ca/2.d0
c
       c3m=4.d0*cf
c       c2mel=-2.d0*cf*dlog(mu**2/m2)-2.d0*cf*dlog(m2/shel)
c     & -2.d0*cf*dlog(-uwel/m2)-2.d0*dlog(-twel/m2)
c       c2minel=-2.d0*cf*dlog(mu**2/m2)-2.d0*cf*dlog(m2/sh)
c     & -2.d0*cf*dlog(-uwinel/m2)-2.d0*dlog(-twinel/m2)
c     & -2.d0*ca*dlog(-uel/m2)-(cf+ca)*dlog(mu**2/shel)
cc       T2mel=c2mel+(cf+ca)*dlog(mu**2/m2)
c       c2minel=2.d0*ReGSin-cf-ca-2.d0*cf*dlog(-tinel/m2)
c     & -2.d0*ca*dlog(-uinel/m2)-(cf+ca)*dlog(mu**2/sh)
c       T2minel=c2minel+(cf+ca)*dlog(mu**2/m2)     
c
c The NNLL MSbar inelastic \delta(s_4) term
       c2mel=-2.d0*cf*dlog(mu**2/m2)
       T2mel=-2.d0*cf*(dlog(-uwel/m2) + dlog(-twel/m2) + dlog(m2/shel))
       c2minel=-2.d0*cf*dlog(mu**2/m2)
       T2minel=-2.d0*cf*(dlog(-uwinel/m2) +
     & dlog(-twinel/m2) + dlog(m2/sh))
       c2el=c2mel+T2mel
       c2inel=c2minel+T2minel
c
c The NNLL MSbar elastic \delta(s_4) term
      c1minel=(dlog(-twinel/m2)+dlog(-uwinel/m2)-3.d0/2.d0)
     & *dlog(mu**2/m2)*cf
c
c      c1mel=(cf*dlog(-tel/m2)+ca*dlog(-uel/m2)-3.d0/4.d0*cf)
c     & *dlog(mu**2/m2)
c
      c1mel=(dlog(-twel/m2)+dlog(-uwel/m2)-3.d0/2.d0)
     & *dlog(mu**2/m2)*cf
c MSbar NLL
      df2ds4inel=c3m**2/2.d0*lns4**3/ss4
     & + (3.d0/2.d0*c3m*c2inel - betaz/4.d0*c3m)
     & *lns4**2/ss4
c scale+zeta NNLL terms above is correct.
     & +lns4/ss4*(c3m*c1minel+c2minel**2
     & +2.d0*c2minel*T2minel
     & +betaz/4.d0*c3m*dlog(mu**2/m2)
     & -zeta2*c3m**2)
c log^2(scale)+zeta 1/s4 terms, above is fixed
     & +1.d0/ss4*(c2minel*c1minel
     & +c2minel*betaz/4.d0*dlog(mu**2/m2)
     & +cf*betaz/4.d0*(dlog(mu**2/m2))**2
     & -zeta2*c2inel*c3m+zeta3*c3m**2)
c inelastic part is fixed
      df2ds4el=c3m**2/2.d0*lns4**3/ss4
     & +(3.d0/2.d0*c3m*c2el-betaz/4.d0*c3m)
     & *lns4**2/ss4
c scale+zeta NNLL terms
     & +lns4/ss4*(c3m*c1mel+c2mel**2
     & +2.d0*c2mel*T2mel
     & +betaz/4.d0*c3m*dlog(mu**2/m2)
     & -zeta2*c3m**2)
c log^2(scale)+zeta 1/s4 terms
     & +1.d0/ss4*(c2mel*c1mel
     & +c2mel*betaz/4.d0*dlog(mu**2/m2)
     & +cf*betaz/4.d0*(dlog(mu**2/m2))**2
     & -zeta2*c2mel*c3m+zeta3*c3m**2)
c fixed elastic bits, I think`
c Computation of integral and surface term:
       int1=df2ds4inel*sfpartinel*hjacob*sigbs4
     &  -df2ds4el*sfpartel*ejacob*sigb0
c
c      print '(''test='',d20.7)',df2ds4inel
      s7 =coef*coef*(int1)
C
      RETURN 
      end
c
c_______________________________________________________________________
c***********************************************************
      double precision function qqdtdu(sh,m2,t1,u1)
      implicit double precision(a-z)
c
      PARAMETER( PI = 3.14159265359D0)
      common/parsqq/coef,cqq
      common/qcdpar/q2,qcdl4
      common/pdfalfa/alfas

c  Conversion factor from gev^-2 to units of 10^-30 cm^2
c
      norm = 19.733*19.733
c  the coefficients
      ncap = 3.0d0
      ncap2 = ncap*ncap
      cf=(ncap2-1.d0)/2.d0/ncap
      cqq=cf
      coef=alfas/pi
c
      sw2=0.23120
      cw2=1-sw2
      sw4=sw2*sw2
      cw4=cw2*cw2
      c2w=cw2-sw2
      c2w2=c2w*c2w
c the c2w is cos (2theta2)
c as if this doesn't evolve lolol
      kg=0.1d0
      ee=0.31333d0
      gg2=4.d0*pi*alfas
      alfa=1.d0/128.d0
      alfa2=alfa*alfa 
      v=dsqrt(1-4.d0*m2/sh)
      v3=v*v*v
      eq=-1.d0/3.d0
      Tl=-0.5d0
      cv=Tl-2.d0*q*sw2
      ca=Tl
      tb=30d0
      cb=1.d0/tb
c uhhh about that alpha
c    
      tm = t1 + m2
      um = u1 + m2
      mz2 = 91.1876d0*91.1876d0
      mt=173.d0
      mt2=mt**2
      mt4=mt**4
      MW2=80.385d0**2
c
c      coef=pi/9.d0*alfa2*v3
c      brkt=eq**2+0.5d0*eq*cv*c2w/sw2/cw2*sh/(sh-mz2) + 
c     & 1.d0/16.d0*(cv2+ca2)*c2w2/sw4/cw4*sh**2/(sh-mz2)**2
c
      mu=dsqrt(q2)
c Matrix element squared
c
c       Matrel2g=(4*ee**2*gg2*kg**2*(m2-sh-tm)
c     & *(m2**3-m2**2*sh-2.d0*sh*tm**2+m2*tm*(3.d0*sh+tm)))
c     & /(3.d0*m2*sh*(m2-tm)**2)
c      Matrel2g=coef*brkt
c uhhh so that's bb to hh... I messed up 
c fixed sw4       
      coef=pi*alfa2*mt4*cb**2/96.d0/sw4/MW2/(um-mt2)**2
      brk1=sh+(um-MW2)*(tm-MW2)/MW2
      brk2=(1.d0-(MW2+m2)/sh)**2-4.d0*MW2*m2/sh**2
      Matrel2g=coef*brk1*dsqrt(brk2)
c  
      c1=pi*alfa2*mt4*cb**2/96.d0/sw4/MW2
c
      qqdtdu = norm*Matrel2g/16.d0/pi/sh**2
c       qqdtdu = norm*Matrel2g/sh**2
c      dq= (2*c1*(2*m2**4*mw2**2*(3*(mt2 + tm) - 2*um) + 
c     -(tm + um)**3*(2*mw2**4 + tm**3*(tm + um) 
c     -   + mt2*(mw2 - tm)*(mw2 - um)*(2*mw2 + tm + um) - 
c     -     2*mw2**3*(3*tm + 2*um) - mw2*tm*(tm**2 + 4*tm*um + um**2) + 
c     -     mw2**2*(4*tm**2 + 7*tm*um + 2*um**2)) - m2**3*mw2*(-32*mw2**3
c     - + 8*mw2**2*tm + 24*mw2*tm**2 + tm**3 + 3*(-4*mw2 + tm)**2*um + 
c     -         (-16*mw2 + 3*tm)*um**2 - um**3 + 
c     -2*mt2*(4*mw2**2 + 3*tm**2 + 7*tm*um + 3*um**2 - 4*mw2*(tm + um))) 
c     - + m2**2*(32*mw2**5 + tm*um*(tm + um)**3 - 14*mw2**4*(7*tm + 6*um)
c     -  + 2*mw2**2*tm*(7*tm**2 - 28*tm*um - 13*um**2) + 
c     -    4*mw2**3*(11*tm**2 + 33*tm*um + 14*um**2) + 
c     -    4*mw2*(2*tm**4 + 4*tm**3*um - 3*tm*um**3 - um**4) + 
c     -      mt2*(30*mw2**4 - 20*mw2**3*(tm + um) + (tm + um)**4 + 
c     -    4*mw2*(tm + um)*(tm**2 + tm*um + um**2)  -
c     -    2*mw2**2*(7*tm**2 - 8*tm*um + 7*um**2))) - 
c     -      m2*(4*mw2**5*(tm + 2*um) - 4*mw2**4*(tm + um)*(tm + 2*um) + 
c     -      tm*(tm + um)**3*(tm**2 + 4*tm*um + um**2) + 
c     -mw2*(tm + um)**2*(4*tm**3 - 27*tm**2*um - 16*tm*um**2 - 3*um**3) +
c     -2*mw2**2*(tm + um)*(9*tm**3 + 28*tm**2*um + 29*tm*um**2 + 8*um**3)
c     - - mw2**3*(23*tm**3 + 57*tm**2*um + 45*tm*um**2 + 13*um**3) + 
c     - mt2*(4*mw2**5 - 4*mw2**4*(tm + um) + (tm + um)**3*(tm**2 + um**2)
c     -    - mw2*(tm + um)**2*(tm**2 - 16*tm*um + um**2) - 
c     -    2*mw2**2*(tm + um)*(9*tm**2 + 16*tm*um + 9*um**2) + 
c     -  2*mw2**3*(9*tm**2 + 19*tm*um + 9*um**2)))))/
c     -  (mw2*sh**5*(-mt2 + tm)**3*(4*m2*mw2 
c     - - (tm + um)**2)*Sqrt((-4*m2*mw2 + (tm + um)**2)/sh**2))
      dq=(2*c1*Sqrt(-4*mh2*mw2 + (tm + um)**2)*(2*mw2**2 
     &- 2*mw2*(tm + um) + tm*(sh + 2*um)))/
     - (mw2*sh**4*(mt2 - tm)**2*Sqrt((-4*mh2*mw2 + (tm + um)**2)/sh**2))
      qqdtdu=norm*dq
c
c      print '(''qqdtdu='',d20.7)',dq
c
      return
      end
