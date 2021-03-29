program oblfcn
!      version 1.27 March 2021
!
    use param
!
!  Originally developed over 15 years ago by arnie lee van buren and updated and
!  improved many times since then. For more information see the github repository:
!  GitHub.com/MathieuandSpheroidalWaveFunctions/oblate_swf. Especially see
!  the readme file and example input and output files. A description of the methods
!  used in oblfcn is provided in the manuscript 'Accurate calculation of oblate
!  spheroidal wave functions,' available at arXiv.org, identifier 1708.07929,
!  August 2017; revised September 2019.
!
!  purpose:     To calculate the first and second kind oblate
!               radial functions and their first derivatives
!               with respect to x for a range of values of the
!               order m and a specified number of values for the
!               degree l = m, m+1,... for specified values of the
!               size parameter c and shape parameter x.
!               To calculate the first kind oblate angular
!               functions and their first derivatives with
!               respect to eta for a range of values of m, l,
!               and angle coordinate eta for a specified value of c.
!
!  Oblfcn can be run in either double precision or quadruple precision
!  arithmetic. The choice is set in the module param provided in the github
!  repository. If this is not available, then create param as follows:
!    module param
!    integer, parameter :: knd = selected_real_kind(8)
!    logical, parameter :: debug = .true.
!    logical, parameter :: warn = .true.
!    logical, parameter :: output = .true.
!    logical, parameter :: suffix = .true.
!    end module param
!  Set the value of knd in the parenthesis to either 8 for double
!  precision or 16 for quadruple precision arithmetic. The logicals
!  in param are described in the readme file and below in the discussion
!  of the output files.
!
!  Oblfcn now provides good results using double precision (real*8)
!  arithmetic for values of c up to at least 10000, m up to at least
!  1000 and essentially all values of x. It is expected that oblfcn
!  will provide good results for much higher values of m and c. A
!  detailed discussion of the expected accuracy is given in the
!  readme file.
!
!     Input Options
!
!  Input parameters are read from unit 1 in the file oblfcn.dat
!  assumed to be located in the directory of the fortran file.
!  Oblfcn.dat contains the following lines of data.
!
!       line 1:
!          mmin   : minimum value for m. (integer)
!          minc   : increment for m. (integer)
!          mnum   : number of values of m. (integer)
!          lnum   : number of values of l [l=m, l=m+1,
!                   ..., l=m+lnum-1]. (integer)
!                   if lnum is less than 2*c/pi it should
!                   be an even integer. If lnum is chosen to
!                   be odd, oblfcn will increase lnum by one.
!
!       line 2:
!          ioprad : (integer)
!                 : =0 if radial functions are not computed
!                 : =1 if radial functions of only the first
!                      kind and their first derivatives are
!                      computed
!                 : =2 if radial functions of both kinds and
!                      their first derivatives are computed
!
!          iopang : (integer)
!                 : =0 if angular functions are not computed
!                 : =1 if angular functions of the first kind
!                      are computed
!                 : =2 if angular functions of the first kind and
!                      their first derivatives are computed
!
!          iopnorm: (integer)
!                 : =0 if not scaled. The angular functions have
!                      the same norm as the corresponding associated
!                      Legendre function [i.e., we use the Meixner
!                      -Schafke normalization scheme.] This norm
!                      becomes very large as m becomes large. The
!                      angular functions are computed below as
!                      a characteristic and an exponent to avoid
!                      overflow.
!                 : =1 if angular functions of the first kind
!                      (and their first derivatives if computed)
!                      are scaled by the square root of the
!                      normalization of the corresponding
!                      associated Legendre function. The resulting
!                      scaled angular functions have unity norm.
!                      This is very useful since it removes the
!                      need to calculate a normalization factor
!                      when using the angular function values given
!                      here. It also eliminates any chance for
!                      overflow when the characteristics and exponents
!                      are combined to form the angular functions.
!
!       line 3:
!          c      : value of the size parameter (= kd/2, where k =
!                   wavenumber and d = interfocal length) [real(knd)]
!          x      : value of the radial coordinate x [real(knd)]
!                   (any value can be entered if ioprad = 0)
!
!       line 4:
!          ioparg : (integer)
!                 : =0 if both arg1 and darg are angles in degrees
!                 : =1 if arg1 and darg are dimensionless values of eta
!
!          arg1   : first value for the angle coordinate (in degrees
!                   or dimensionless if equal to eta) for which angular
!                   functions are to be computed. [real(knd)]
!
!          darg   : increment used to calculate additional desired
!                   arguments for angular functions. [real(knd)]
!
!          narg   : number of desired angle arguments. (integer)
!                   (line 4 is not read when iopang = 0)
!
!     Output files
!
!  Oblfcn offers several several output files: Fort.20 and fort.30
!  list the calculated radial and angular functions. Fort.40 and
!  fort.50 are diagnostic files. Fort.60 provides warning whenever the
!  estimated accuracy falls below a specified minimum, currently set
!  equal to 8. Writing to these files is controlled by logicals specified
!  in the module param. False suppresses the file; true enables it.
!  Debug controls fort.30 and fort.40, warn controls fort.60 and output
!  controls fort.20 and fort.30. The logical suffix controls whether the
!  accuracy estimates given in fort.20 are followed by a letter designating
!  how the accuracy was determined. 'w' indicates it is based on the
!  Wronskian and 'e' designates it is based on subtraction errors involved
!  in the calculations. Setting suffix = false suppresses the letter.
!  Information about these files as well as a discussion about accuracy,
!  expansion d coefficients and eigenvalues is given in the readme file.
!
        real(knd) arg1,c,darg,step1,step2,step3,x,xneu
!
!  open input files
        open(1, file='oblfcn.dat')
!
!  read input data
        read(1,*) mmin,minc,mnum,lnum
        read(1,*) ioprad,iopang,iopnorm
        read(1,*) c,x
        if(iopang.ne.0) read(1,*) ioparg,arg1,darg,narg
!
!  open output files
if (output) then
        open(20, file='fort.20')
        open(30, file='fort.30')
end if
if (debug) then
        open(40, file='fort.40')
        open(50, file='fort.50')
end if
if (warn) then
        open(60, file='fort.60')
end if
!
!  Here is where the user sets kindd, the value for kind that
!  corresponds to double precision data on the user's computer.
!  Similarly, this is where kindq, the value of kind for quadruple
!  precision data, is set. These values are set below to 8 and 16,
!  respectively. They should be changed to the kind values for double
!  precision and quadruple precision if those values are different than
!  these.

5       kindd=8
        kindq=16
!
!  set the minimum desired accuray minacc to 10 for real*8
!  arithmetic and to 15 for real*16 arithmetic. These can be
!  changed if desired. See the readme file.
!
        if(knd.eq.kindd) minacc=10
        if(knd.eq.kindq) minacc=15
!
!  ndec: the maximum number of decimal digits available in real(knd).
!        arithmetic.
!  nex:  the maximum exponent available in real(knd).
!
        ndec=precision(c)
        nex=range(c)-1
!
!  set array dimensions
        nbp=int(2*c/3.1416)
        maxe=max(50,nbp+30)+10
        maxe2=maxe+maxe
        if(maxe2.lt.lnum) maxe2=lnum
        maxm=mmin+minc*(mnum-1)
        maxlp=lnum+maxm+1
        maxint=2*(nbp+33)+6
        maxj=lnum+3*ndec+int(c)+5+maxm
        maxp=max(lnum+3*ndec+int(c)+5,maxlp+5)
        maxn=maxj
        maxpdr=2*int(c)+4*ndec+int(100*x)+8
        neta=30
        if(ioprad.ne.2) go to 10
          step1=0.1e0_knd
          step2=0.1e0_knd
          step3=0.8e0_knd
          nstep1=1
          nstep2=1
          nstep3=3
          ngau=100*(2+int(c/500.0e0_knd))
!
          if(x.lt.0.2e0_knd) then
          step1=x/4.0e0_knd
          step2=0.075e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=4
          ngau=100*(2+int(c/500.0e0_knd))
          end if
!
          if(x.lt.0.1e0_knd) then
          step1=x/4.0e0_knd
          step2=0.025e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=4
          ngau=200
          if(c.gt.500.0e0_knd) ngau=200*(2+int(c/1000.0e0_knd))
          end if
!
          if(x.lt.0.05e0_knd) then
          step1=x/15.0e0_knd
          step2=0.002e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=2
          ngau=300
          if(c.gt.500.0e0_knd) ngau=500
          if(c.gt.1000.0e0_knd) ngau=800
          if(c.gt.1500.0e0_knd) ngau=1000
          if(c.gt.2000.0e0_knd) ngau=1200
          if(c.gt.2500.0e0_knd) ngau=1500
          if(c.gt.3000.0e0_knd) ngau=1700
          if(c.gt.3500.0e0_knd) ngau=1900
          if(c.gt.4000.0e0_knd) ngau=2200
          if(c.gt.4500.0e0_knd) ngau=2500
          if(c.gt.5000.0e0_knd) ngau=2500+300*int((c-4500.0e0_knd)/ &
             500.0e0_knd)
          end if
!
          if(x.lt.0.01e0_knd) then
          step1=x/15.0e0_knd
          step2=0.002e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=2
          ngau=300
          if(c.gt.500.0e0_knd) ngau=600
          if(c.gt.1000.0e0_knd) ngau=800
          if(c.gt.1500.0e0_knd) ngau=1000
          if(c.gt.2000.0e0_knd) ngau=300*int(c/500.0e0_knd)
          end if
!
          if(x.lt.0.001e0_knd) then
          step1=x/15.0e0_knd
          step2=0.002e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=3
          nstep2=1
          nstep3=2
          ngau=400
          if(c.gt.500.0e0_knd) ngau=600
          if(c.gt.1000.0e0_knd) ngau=800
          if(c.gt.1500.0e0_knd) ngau=1000
          if(c.gt.2000.0e0_knd) ngau=300*int(c/500.0e0_knd)
          end if
!
        xneu=0.3e0_knd
        if(c.gt.100.0e0_knd) xneu=0.04e0_knd
        if(c.gt.600.0e0_knd) xneu=0.03e0_knd
        if(c.gt.800.0e0_knd) xneu=0.01e0_knd
!
        if(x.lt.0.01e0_knd.or.ioprad.eq.0) go to 10
          if(knd.eq.kindd) then
          if(x.ge.0.01e0_knd) maxn=2*int(25/(x*x)+300/x+3*c+ &
                                   1250*knd)+4
          if(x.ge.0.1e0_knd) maxn=2*int((lnum+c/5+0.5e0_knd*maxm+200)* &
                                  1.4e0_knd/x)+4
          if(x.ge.0.5e0_knd) maxn=2*int((lnum+c/5+0.5e0_knd*maxm+300)/x) &
                                  +4
          if(x.ge.1.0e0_knd) maxn=2*int(lnum+c/5+0.5e0_knd*maxm+300)+4
          end if
          if(knd.eq.kindq) then
          if(x.ge.0.01e0_knd) maxn=2*int(25/(x*x)+400/x+3*c+ &
                                   1250*knd)+4
          if(x.ge.0.1e0_knd) maxn=2*int((lnum+c/5+0.5e0_knd*maxm+350)* &
                                  1.4e0_knd/x)+4
          if(x.ge.0.5e0_knd) maxn=2*int((lnum+c/5+0.5e0_knd*maxm+400)/x) &
                                  +4
          if(x.ge.1.0e0_knd) maxn=2*int(lnum+c/5+0.5e0_knd*maxm+400)+4
          end if
        maxn=maxn+maxm
10      maxp=max(maxn,maxp,maxpdr)
        maxq=lnum+3*ndec+int(c)+maxm+maxm+4
        if(c.le.60.0e0_knd.and.c.ge.10.0e0_knd.and.mmin.le.40 &
            .and.x.le.0.99e0_knd) maxq=max(maxq,250-int(50*x)+maxm+ &
                                       maxm+4)
        maxdr=maxpdr/2+1
        maxd=maxn/2+1
        maxmp=maxm+maxm+5
        maxt=1
        jnebmax=30
        jnenmax=3
        if(x.le.0.05e0_knd) jnenmax=1
        if(iopang.ne.0) maxt=narg
!
        call main (mmin,minc,mnum,lnum,c,ioprad,iopang,iopnorm,minacc, &
                   x,ngau,step1,nstep1,step2,nstep2,step3,nstep3,ioparg, &
                   arg1,darg,narg,maxd,maxdr,maxe,maxe2,maxint,maxj, &
                   maxlp,maxm,maxmp,maxn,maxp,maxpdr,maxq,maxt,neta, &
                   jnenmax,jnebmax,xneu,ndec,nex,kindd,kindq)
!
        end program
!
!
        subroutine main (mmin,minc,mnum,lnum,c,ioprad,iopang,iopnorm, &
                         minacc,x,ngau,step1,nstep1,step2,nstep2,step3, &
                         nstep3,ioparg,arg1,darg,narg,maxd,maxdr,maxe, &
                         maxe2,maxint,maxj,maxlp,maxm,maxmp,maxn,maxp, &
                         maxpdr,maxq,maxt,neta,jnenmax,jnebmax,xneu, &
                         ndec,nex,kindd,kindq)
!
!  purpose:     To coordinate the calculation of both the oblate
!               spheroidal radial and angular functions and their
!               first derivatives using various algorithms.
!
!     input:    mmin   : minimum desired value of m
!               minc   : increment in m used to compute other values
!               mnum   : number of values of m that are desired
!               lnum   : desired number of values of l = m, m + 1, ...,
!                        m + lnum - 1
!               c      : size parameter
!               ioprad : equal to 0 if no radial functions are desired;
!                        equal to 1 if only radial functions of the
!                          first kind and their first derivatives are
!                          desired;
!                        equal to 2 if radial functions of both kinds
!                          and their first derivatives are desired
!               iopang : equal to 0 if no angular functions are desired;
!                        equal to 1 if only angular functions of the
!                          first kind are desired;
!                        equal to 2 if angular functions of the first
!                          kind and their first derivatives are desired
!               iopnorm: equal to 0 when the angular functions have
!                        the same norm as the corresponding associated
!                        Legendre functions;
!                        equal to 1 when the angular functions are
!                        scaled by the square root of the normalization
!                        of the corresponding Legendre function, giving
!                        them unity norm
!               minacc : desired minimum accuracy for the radial
!                        functions
!               x      : radial coordinate
!               ngau   : order of the Gaussian quadrature to be used in
!                        computing integrals in subroutine pint for use
!                        in subroutine r2int where the integal method
!                        is used to calculate r2 and r2d
!               step1  : first step in eta for the Gaussian quadrature
!                        integration
!               nstep1 : number of equal substeps that step1 is divided
!                        into
!               step2  : second step in eta for the Gaussian quadrature
!                        integration
!               nstep2 : number of equal substeps that step2 is divided
!                        into
!               step3  : third and final step in eta for the Gaussian
!                        quadrature integration
!               nstep3 : number of equal substeps that step3 is divided
!                        into
!               ioparg : =0 if both arg1 and darg are angles in degrees
!                        =1 if arg1 and darg are dimensionless values
!                           of eta
!               arg1   : first value for the angle coordinate (in
!                        degrees or dimensionless if equal to eta)
!                        for which angular functions are to be computed
!               darg   : increment used to calculate additional desired
!                        arguments for angular functions.
!               narg   : number of desired angle arguments.
!               maxd   : dimension of enr array containing ratios of
!                        the expansion d coefficients
!               maxdr  : dimension of drhor array containing special d
!                        coefficient ratios used in subroutine r2leg
!                        when computing the sum of Legendre functions of
!                        the first kind that appear in the Legendre
!                        function expansion for r2 and r2d
!               maxe   : dimension of both the vector of diagonal
!                        elements and the vector of subdiagonal elements
!                        used to obtain accurate estimates of the
!                        eigenvalues for both even and odd l - m
!               maxe2  : dimension of the vector containing all of the
!                        eigenvalue estimates, including both even and
!                        odd l - m (=2*maxe)
!               maxint : maximum number of terms available for computing
!                        r2 and r2d in the subroutine r2int; dimension
!                        of the arrays of integrals computed in
!                        subroutine pint
!               maxj   : equal to the dimension of the array of ratios
!                        of spherical Bessel functions of the first kind
!                        and the array of ratios of the first derivative
!                        of this Bessel function to the corresponding
!                        Bessel function
!               maxlp  : maximum value desired for l
!               maxm   : maximum value desired for m
!               maxmp  : maxm + 5; dimension of the integer array norme
!                        used in scaling of the Neumann functions in
!                        the integrands in subroutine pint
!               maxn   : dimension of the arrays of Neumann function
!                        ratios used in computing r2 and r2d
!               maxp   : dimension of arrays of Legendre functions of
!                        the first kind used in computing angular
!                        functions, in computing integrands in
!                        subroutine pint and in computing r2 and r2d in
!                        subroutine r2eta
!               maxpdr : dimension of the arrays of ratios of both
!                        Legendre functions of the first kind and their
!                        first derivatives used in the sum of these
!                        functions that contribute to r2 and r2d in
!                        subroutine r2leg
!               maxq   : dimension of arrays of ratios of Legendre
!                        functions of the second kind and ratios of
!                        their first derivatives used in their sum in
!                        subroutine r2leg
!               maxt   : equal to narg if angular functions are
!                        computed where it is the maximum value of the
!                        first index in the arrays of Legendre functions
!                        used in subroutine s1leg;
!                        otherwise equal to 1 to specify the
!                        first index for the Legendre functions used
!                        in the variable eta method for computing r2
!                        and r2d in subroutine r2eta
!               neta   : number of values available for eta in the
!                        variable eta method for calculating r2 and r2d
!                        (subroutine r2eta) and the variable eta method
!                        for calculating r1 and r1d (subroutine r1eta);
!                        set equal to 30
!               jnenmax: number of arrays of ratios of Legendre and
!                        Neumann functions stored as eta is varied in
!                        subroutine r2eta; set equal to 3 so that the
!                        previous 3 sets of ratios are available
!                        to use without recalculating them when one of
!                        thes previous values for eta is used again for
!                        a later value of l. If r2eta will not be
!                        called, jnenmax is set equal to unity.
!               jnebmax: number of arrays of ratios of Legendre and
!                        Bessel functions stored as eta is varied in
!                        subroutine r1eta; set equal to the number of
!                        values of eta that are available. Thus the
!                        ratios never need to be recalculated when one
!                        of the previous values used for eta is used
!                        again for a later value of l
!               xneu   : minimum value of x for which either the
!                        Neumann expansion method r2neu0 or the variable
!                        eta method r2eta will be called
!               ndec   : number of decimal digits for real(knd)
!               nex    : maximum exponent for real(knd)
!               kindd  : kind value for double precision real data
!               kindq  : kind value for quadruple precision real data

        use param
!
!  real(knd) scalars
        real(knd) aj1,aj2,apcoef,apcoefn,api,arg1,c,c2,c4,coefn,coefme, &
                  coefmo,darg,dec,dfnorm,dmlf,dmfnorm,dmlmf,dmsnorm, &
                  dmlms,dmlms1,dneg,dc01,eigval,em,etaval,fac1,fac1d, &
                  fac2,fac2d,factor,factor1,pcoefn,pi,qdm0,qdm1,qm0, &
                  qm1,rm,rm2,r1cin,r1cm,r1dcin,r11c,r1dcm,r1d1c,r1ec, &
                  r1dec,r2ic,r2dic,r2lc,r2dlc,r2l1c,r2dl1c,r2nc,r2dnc, &
                  r2ec,r2dec,step1,step2,step3,ten,termpq,t1,t2,t3,t4, &
                  t5,t6,t7,wm,x,xb,xbninp,xl,xneu,wronc,wronca,wroncb, &
                  wront
        real(knd) ang,apcoef1,etaval1,pcoefe,pcoefet,pcoefo,pdcoefe, &
                  pdcoefo,pdcoefet,pcoefe1,pcoefet1,pcoefo1,pdcoefe1, &
                  pdcoefo1,pdcoefet1
!
!  integer and real(knd) arrays with dimension lnum
        integer   iqdl(lnum),iql(lnum)
        real(knd) eig(lnum),r1c(lnum),r1dc(lnum),r2c(lnum),r2dc(lnum), &
                  qdl(lnum),ql(lnum)
        integer   ir1e(lnum),ir1de(lnum),ir2e(lnum),ir2de(lnum), &
                  match(lnum)
!
!  integer and real(knd) arrays with dimension lnum+1
        integer ifajo(lnum+1)
        real(knd) fajo(lnum+1)
!
!  real(knd) arrays with dimension maxd
        real(knd) enr(maxd),bliste(maxd),gliste(maxd), &
                  blisto(maxd),glisto(maxd)
!
!  real(knd) arrays with dimension maxdr
        real(knd) drhor(maxdr)
!
!  real(knd) arrays with dimension maxint
        real(knd) pint1(maxint),pint2(maxint),pint3(maxint), &
                  pint4(maxint),rpint1(maxint),rpint2(maxint)
!
!  real(knd) array with dimension maxj
        real(knd) sbesf(maxj),sbesdf(maxj),sbesfe(maxj),sbesdfe(maxj), &
                  sbesfsv(jnebmax,maxj),sbesdfsv(jnebmax,maxj)
!
!  integer and real(knd) arrays with dimension maxlp
        integer ibese(maxlp),ineue2(maxlp),ipnormint(maxlp), &
                ineuee(maxlp),ineuesv(jnenmax,maxlp)
        integer ibesee(maxlp),ibesesv(jnebmax,maxlp)
        real(knd) pnormint(maxlp),sbesdr(maxlp),sbesn(maxlp), &
                  sneun2(maxlp),sneudr2(maxlp),sneune(maxlp), &
                  sneudre(maxlp),sneunsv(jnenmax,maxlp), &
                  sneudrsv(jnenmax,maxlp),sbesne(maxlp),sbesdre(maxlp), &
                  sbesnsv(jnebmax,maxlp),sbesdrsv(jnebmax,maxlp)
!
!  real(knd) arrays with dimensio maxmp
        real(knd) enrneg(maxmp),qdqr(maxmp)
!
!  real(knd) arrays with dimension maxn
        real(knd) sneuf2(maxn),sneudf2(maxn),prat1(maxn),sneufe(maxn), &
                  sneudfe(maxn),sneufsv(jnenmax,maxn), &
                  sneudfsv(jnenmax,maxn)
!
!  real(knd) arrays with dimension given by maxp
        real(knd) alpha(maxp),beta(maxp),coefa(maxp),coefb(maxp), &
                  coefc(maxp),coefd(maxp),coefe(maxp),gamma(maxp), &
                  pdr(maxt,maxp),pdrat(maxt,maxp),pdratt(maxp), &
                  pr(maxt,maxp),prat(maxt,maxp),pratb(maxp), &
                  pratbsv(jnenmax,maxp),prattsv(jnenmax,maxp), &
                  pdrattsv(jnenmax,maxp),pratt(maxp)
!
!  real(knd) arrays with dimension given by maxp
        real(knd) pratb1(maxp),pratt1(maxp),pdratt1(maxp), &
                  pratbsv1(jnebmax,maxp),prattsv1(jnebmax,maxp), &
                  pdrattsv1(jnebmax,maxp)
!
!  real(knd) arrays with dimension maxpdr
        real(knd) prx(maxpdr),pdrx(maxpdr)
!
!  real(knd) arrays with dimension maxq
        real(knd) qrat(maxq),qdrat(maxq)
        real(knd) qrat1(maxq),qdrat1(maxq)
!
!  real(knd) and integer arrays with dimension maxt
        real(knd) arg(maxt),barg(maxt),etainp(maxt),pdnorm(maxt), &
                  pdnorma(maxt),pnorm(maxt),pnorma(maxt),pdtempe(maxt), &
                  pdtempo(maxt),ptempe(maxt),ptempo(maxt),s1c(maxt), &
                  s1dc(maxt),xin(maxt),xlninp(maxt)
        integer ipdnorm(maxt),ipdnorma(maxt),ipnorm(maxt), &
                ipnorma(maxt),ipdtempe(maxt),ipdtempo(maxt), &
                iptempe(maxt),iptempo(maxt),is1e(maxt),is1de(maxt), &
                naccs(maxt),naccds(maxt)
!
!  real(knd) arrays with dimension ngau
        real(knd) wr(ngau),xr(ngau)
!
!  real(knd) arrays with dimension maxe2
        real(knd) eigt(maxe2),f(maxe2)
!
!  real(knd) arrays with dimension maxe
        real(knd) d(maxe),e(maxe)
!
!  real(knd) arrays with dimension neta
        real(knd) eta(neta),xbn(neta),xln(neta),wmeta2(neta)
!
!  miscellaneous integer arrays
        integer neeb(jnenmax),neeb1(jnebmax),limpsv(jnenmax), &
                limp1sv(jnebmax),limnsv(jnenmax),jelimsv(jnenmax), &
                jelim1sv(jnebmax),limjsv(jnebmax)
        
        character chr_e, chr_w
        if (suffix) then
            chr_e = 'e'
            chr_w = 'w'
        else
            chr_e = ' '
            chr_w = ' '
        end if

        ten=10.0e0_knd
        dec=ten**(-ndec-1)
        factor=ten**(nex-20)
        factor1=1.0e0_knd/factor
        ifactor=nex-20
        jtest=ndec-minacc-2
        pi=acos(-1.0_knd)
        api=pi/180.0e0_knd
        c2=c*c
        c4=c2*c2
!
!  begin loops
          if(iopang.eq.0) go to 10
            do jarg=1,narg
            arg(jarg)=arg1+(jarg-1)*darg
            if(ioparg.eq.0) barg(jarg)=cos(arg(jarg)*api)
            if(ioparg.eq.1) barg(jarg)=arg(jarg)
            end do
10        continue
          igau=0
          ibflag1=1
          ibflag2=1
if (debug) then
          if(knd.eq.kindd.and.ioprad.ne.0) write(40,20) x,c
20        format(1x,'x = ',e23.14,/,1x,'c = ',e23.14)
          if(knd.eq.kindq.and.ioprad.ne.0) write(40,25) x,c
25        format(1x,'x = ',e38.30,/,1x,'c = ',e38.30)
end if
          wront=1.0e0_knd/(c*(x*x+1.0e0_knd))
            do 1540 mi=1,mnum
            m=mmin+minc*(mi-1)
            em=m
            m2=m+m
if (debug) then
            if(knd.eq.kindd.and.iopang.ne.0) write(50,30) c,m
30          format(1x,'c = ',e23.14,'; m = ',i5)
            if(knd.eq.kindq.and.iopang.ne.0) write(50,35) c,m
35          format(1x,'c = ',e38.30,'; m = ',i5)
            if(ioprad.ne.0) write(40,40) m
40          format(1x,'m = ',i5)
end if
if (output) then
            if(knd.eq.kindd.and.iopang.ne.0) write(30,50) c,m
50          format(1x,e23.14,i5)
            if(knd.eq.kindq.and.iopang.ne.0) write(30,55) c,m
55          format(1x,e38.30,i5)
end if
            rm=m
            rm2=m+m
            icounter=0
            limcsav=0
            jbes=3*ndec+int(c)
            iopbes=1
            iopeta1=0
            jflageta1=0
            iopint=0
            iopleg=0
            iopleg1=0
            iopneu0=0
            iopeta=0
            if(ioprad.ne.2) go to 60
              if(knd.eq.kindd) then
              if(x.le.0.99e0_knd.and.c.le.20.0e0_knd) iopleg=1
              if(x.gt.0.99e0_knd.and.c.le.20.0e0_knd) iopneu0=1
              end if
              if(knd.eq.kindq) then
              if(x.le.0.99e0_knd.and.c.le.60.0e0_knd) iopleg=1
              if(x.gt.0.99e0_knd.and.c.le.40.0e0_knd) iopneu0=1
              end if
60          continue
            jneu1max=0
            jnen=0
            incnee=1
            if(x.gt.0.4e0_knd) incnee=2
            iplflag=0
            nee=28
            if(x.gt.0.1e0_knd) nee=26
            if(x.gt.0.2e0_knd) nee=24
            if(x.gt.0.3e0_knd) nee=22
            if(x.gt.0.4e0_knd) nee=20
            if(x.gt.0.5e0_knd) nee=18
            if(x.gt.0.6e0_knd) nee=14
            if(x.gt.0.7e0_knd) nee=12
            if(x.gt.0.8e0_knd) nee=8
            if(x.gt.0.9e0_knd) nee=2
            jnen1=0
            incnee1=1
            nee1=1
            idir=0
            iflageta1=0
            nsub1p=ndec
            nsubd1p=ndec
            if(iopang.eq.0) go to 70
            limps1=lnum+3*ndec+int(c)
            if((limps1+3).gt.maxp) limps1=maxp-3
            iopd=0
            if(iopang.eq.2) iopd=1
            call pleg(m,limps1,maxp,limcsav,iopd,ndec,nex,barg,narg, &
                      maxt,pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm,alpha, &
                      beta,gamma,coefa,coefb,coefc,coefd,coefe)
            limcsav=limps1
70          limj=lnum+3*ndec+int(c)+maxm
            if(ioprad.eq.0.or.x.ne.0.0e0_knd) go to 170
!
!  calculation of factors for radial functions when x = 0
            fac1=1.0e0_knd
            ifac1=0
            if(m.eq.0) go to 100
              do 90 k=1,m
              fac1=fac1*c/(k+k+1)
              if(fac1.lt.factor) go to 80
              fac1=fac1/factor
              ifac1=ifac1+ifactor
              go to 90
80            if(fac1.gt.factor1) go to 90
              fac1=fac1*factor
              ifac1=ifac1-ifactor
90            continue
100         iterm=int(log10(fac1))
            fac1=fac1*(ten**(-iterm))
            ifac1=ifac1+iterm
            fac1d=fac1*c/(m+m+3)
            ifac1d=ifac1
            if(ioprad.eq.1) go to 170
            fac2=(m+m-1)*pi/(2.0e0_knd*c)
            ifac2=0
            if(m.eq.0) go to 130
              do 120 k=1,m
              fac2=fac2*c/(k*2.0e0_knd)
              if(fac2.lt.factor) go to 110
              fac2=fac2/factor
              ifac2=ifac2+ifactor
              go to 120
110           if(fac2.gt.factor1) go to 120
              fac2=fac2*factor
              ifac2=ifac2-ifactor
120           continue
130         iterm=int(log10(abs(fac2)))
            fac2=fac2*(ten**(-iterm))
            ifac2=ifac2+iterm
            fac2d=(m+m-1)*(m+m-3)*(m+m+1)*pi/(c*c*2.0e0_knd)
            ifac2d=0
            if(m.eq.0) go to 160
              do 150 k=1,m
              fac2d=fac2d*c/(k*2.0e0_knd)
              if(fac2d.lt.factor) go to 140
              fac2d=fac2d/factor
              ifac2d=ifac2d+ifactor
              go to 150
140           if(fac2d.gt.factor1) go to 150
              fac2d=fac2d*factor
              ifac2d=ifac2d-ifactor
150         continue
160         iterm=int(log10(abs(fac2d)))
            fac2d=fac2d*(ten**(-iterm))
            ifac2d=ifac2d+iterm
170         xb=sqrt(x*x+1.0e0_knd)
            if(ioprad.eq.0.or.x.eq.0.0e0_knd) go to 180
            prat1(1)=1.0e0_knd
            prat1(2)=rm2+1.0e0_knd
              do jp=3,limj
              aj1=jp-1
              aj2=jp-2
              prat1(jp)=(rm2+aj1)*(rm2+aj2)/(aj1*aj2)
              end do
            pcoefn=(x*x+1.0e0_knd)/(x*x)
            apcoefn=(rm/2.0e0_knd)*log10(pcoefn)
            ipcoefn=int(apcoefn)
            pcoefn=ten**(apcoefn-ipcoefn)
            if(mi.ne.1) go to 180
            call sphbes(c,x,limj,maxj,maxlp,sbesf,sbesdf,sbesn,ibese, &
                        sbesdr)
180         continue
!  obtain starting eigenvalues eigt(i) for the Bouwkamp procedure
!
            nbp=int(2*c/3.1416)
            limeig=max(67,(4*nbp)/3)
            lime=max(50,nbp+30)
            lime2=lime+lime
              do 190 i=1,lime
              xl=m+i+i-2
              d(i)=xl*(xl+1.0e0_knd)/c2-(2.0e0_knd*xl*(xl+1.0e0_knd)- &
                   2.0e0_knd*em*em-1.0e0_knd)/((2.0e0_knd*xl-1.0e0_knd)* &
                   (2.0e0_knd*xl+3.0e0_knd))
190           continue
            nm1=lime-1
              do 200 i=1,nm1
              xl=m+i+i-2
              e(i)=(-1.0e0_knd/(2.0e0_knd*xl+3.0e0_knd))* &
                   sqrt(((xl+2.0e0_knd+em)*(xl+1.0e0_knd+em)* &
                   (xl+2.0e0_knd-em)*(xl+(1.0e0_knd)-em))/ &
                   ((2.0e0_knd*xl+5.0e0_knd)*(2.0e0_knd*xl+1.0e0_knd)))
200           continue
            call tridiag(d,e,lime,maxe)
              do 210 i=1,lime
              f(i+i-1)=c2*d(i)
210           continue
              do 220 i=1,lime
              xl=m+i+i-1
              d(i)=xl*(xl+1.0e0_knd)/c2-(2.0e0_knd*xl*(xl+1.0e0_knd)- &
                   2.0e0_knd*em*em-1.0e0_knd)/((2.0e0_knd*xl-1.0e0_knd)* &
                   (2.0e0_knd*xl+3.0e0_knd))
220           continue
            nm1=lime-1
              do 230 i=1,nm1
              xl=m+i+i-1
              e(i)=(-1.0e0_knd/(2.0e0_knd*xl+3.0e0_knd))* &
                   sqrt(((xl+2.0e0_knd+em)*(xl+1.0e0_knd+em)* &
                   (xl+2.0e0_knd-em)*(xl+(1.0e0_knd)-em))/ &
                   ((2.0e0_knd*xl+5.0e0_knd)*(2.0e0_knd*xl+1.0e0_knd)))
230           continue
            call tridiag(d,e,lime,maxe)
              do 240 i=1,lime
              f(i+i)=c2*d(i)
240           continue
              do 250 i=1,lime2
              eigt(i)=f(i)
250           continue
!  determine the number of leading decimal digits of agreement
!  for lower order paired eigenvalues
            listart=1
            matlim=min(lnum,nbp+1)
            matlilim=matlim+1
              do i=2,matlim,2
              match(i)=-int(log10(abs((eigt(i)-eigt(i-1))/(eigt(i)) &
                        +ten*dec)))
              if(match(i).lt.0) match(i)=0
              if(match(i).gt.ndec) match(i)=ndec
              if(match(i).gt.minacc+1) listart=i+1
                if(match(i).lt.5) then
                matlilim=i-1
                go to 255
                end if
              end do
255         continue
              if(listart.gt.1) then
              iopleg=0
              iopneu0=0
              end if
!
!  prepare for do loop over index l
            iflag=0
            iflagint=0
            iflagneu=0
            iqlegflag=0
            itestmp=0
            ipint=0
            intlim=maxint-3
            istartint=0
            if(x.lt.0.0005e0_knd) istartint=3
            if(istartint.eq.3) iopint=0
            if(iopint.eq.0.and.x.le.0.99e0_knd) iopleg=1
            nflag=0
            iflageta=0
            ieta=0
            naccr=minacc+1
            naccneu0p=0
            naccneu0p2=0
            naccleg=0
            nacclegp=minacc
            legflag=0
            jflagleg=0
            legstart=m
            ioppsum=1
            iopqnsum=1
            jneu1=0
            jneu1s=0
            jneu0=0
            jneu0s=0
            icflag1=0
            nacintsa=0
            jeta=0
            jeta1=0
            jbes=0
            jint=0
            jlegp=0
            jleg=0
            jleg1=0
            jtest0=0
            jtest1=0
            naccetas=0
            jetaflag=0
            naccleg1=0
            ir1ep=0
            istartneu0=1
            naccflag=0
            jneu0x=0
            limax=m
if (output) then
            if(knd.eq.kindd.and.ioprad.ne.0) write(20,260) x,c,m
260         format(1x,e23.14,e23.14,i5)
            if(knd.eq.kindq.and.ioprad.ne.0) write(20,265) x,c,m
265         format(1x,e38.30,e38.30,i5)
end if
              do 1510 li=1,lnum
              l=m+(li-1)
if (output) then
              if(iopang.ne.0) write(30,270) l
270           format(1x,i6)
end if
if (debug) then
              if(iopang.ne.0) write(50,280) l
280           format(1x,'l = ',i6)
end if
              ix=l-m-2*((l-m)/2)
              if(ix.eq.0) naccflag=0
              if(li.eq.1) naccrsav=minacc
              naccr=-1
              naccint=0
              naccleg=0
              naccneu0=0
              naccneu1=0
              naccleg1=0
              nacceta=0
              naccrt=0
              jflagl=0
              if(ioprad.eq.2.and.li.eq.listart.and.iopneu0.eq.0.and. &
                  x.ge.xneu) iopneu0=1
              limdrad=3*ndec+int(c)+10
              if(ioprad.ne.0.and.li.ne.1) limdrad=jbes+jbes+20+ &
                                          int(sqrt(c))
              if(iopint.ne.0.and.li.gt.listart.and.jint.gt.jbes) &
                  limdrad=jint+jint+20+int(sqrt(c))
              limdang=3*ndec+int(c)
              if(iopang.ne.0.and.li.ne.1) limdang=jang+jang+20+ &
                                                  int(sqrt(c))
              if(iopang.eq.0) limd=limdrad
              if(ioprad.eq.0) limd=limdang
              if(iopang.ne.0.and.ioprad.ne.0) limd=max(limdang,limdrad)
              if(li.eq.1) limnorm=limdrad
              if(li.gt.1) limnorm=jnorm+jnorm+20+int(sqrt(c))
              limd=max(limd,limnorm)
              if(ioprad.ne.2) go to 290
              if(iopleg.eq.1) limdleg=l-m+3*ndec+int(c)
              if(iopleg.eq.2) limdleg=jleg+jleg+20+int(sqrt(c))
              if(iopleg.ne.0.and.li.ge.listart) limd=max(limd,limdleg)
                if(knd.eq.kindd) then
                if(x.ge.0.01e0_knd) lneu0=2*int(25/(x*x)+300/x+3*c+ &
                                          1250*knd)
                if(x.ge.0.1e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                         200)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                         200)/x)
                if(x.ge.1.0e0_knd) lneu0=2*int(l-m+c/5+0.5e0_knd*m+200)
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.01e0_knd) lneu0=2*int(25/(x*x)+400/x+3*c+ &
                                          1250*knd)
                if(x.ge.0.1e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                         350)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                         300)/x)
                if(x.ge.1.0e0_knd) lneu0=2*int(l-m+c/5+0.5e0_knd*m+300)
                end if
              if(iopneu0.eq.3.and.jtest0.gt.minacc.and.li.ne.listart+1) &
                     lneu0=jneu0max+jneu0max+40+int(sqrt(c)* &
                           (2/min(1.0e0_knd,x))+100/x)
              if(iopneu0.ne.0) limd=max(limd,lneu0)
                if(knd.eq.kindd) then
                if(x.gt.0.05e0_knd) leta=2*int(25/(x*x)+300/x+3*c+ &
                                         1250*knd)
                if(x.ge.0.1e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+200)* &
                                        1.4e0_knd/x)
                if(x.ge.0.5e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+ &
                                        300)/x)
                if(x.ge.1.0e0_knd) leta=2*int(l-m+c/5+0.4e0_knd*m+300)
                end if
                if(knd.eq.kindq) then
                if(x.gt.0.05e0_knd) leta=2*int(25/(x*x)+400/x+3*c+ &
                                         1250*knd)
                if(x.ge.0.1e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+350)* &
                                        1.4e0_knd/x)
                if(x.ge.0.5e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+ &
                                        400)/x)
                if(x.ge.1.0e0_knd) leta=2*int(l-m+c/5+0.4e0_knd*m+400)
                end if
              if(iopeta.eq.3.and.nacceta.gt.minacc) &
                  leta=jeta+jeta+40+ &
                       int(sqrt(c)*(2/min(1.0e0_knd,x))+5/x)
              if(iopeta.ne.0) limd=max(limd,leta)
290           continue
              if(limd.gt.maxn) limd=maxn
!
              if(li.le.limeig) eigval=eigt(li)
              if(li.gt.limeig) eigval=4.0e0_knd*eig(li-1)-6.0e0_knd* &
                            eig(li-2)+4.0e0_knd*eig(li-3)-eig(li-4)
300           continue
!
!  use Bouwkamp procedure to obtain accurate eigenvalue
              if(l.eq.m) ienre=(3*ndec+int(c)+2*m)/2
              if(l.eq.m) jlowe=1
              if(l.eq.m) limdle=2
              if(l.eq.m+1) ienro=(3*ndec+int(c)+2*m)/2
              if(l.eq.m+1) jlowo=1
              if(l.eq.m+1) limdlo=3
!
!  compute the coefficients in the Bouwkamp method
              if(ix.eq.1) go to 330
!
!  beta coefficients (bliste) for l-m even
              if(limdle.gt.limd) go to 360
              j=jlowe
                do 310 i=limdle,limd,2
                i2=i+i
                t1=i
                t2=i-1
                t3=m2+i
                t4=m2+i-1
                t5=m2+i2-1
                t6=m2+i2-3
                t7=m2+i2+1
                bliste(j)=c4*t1*t2*t3*t4/(t5*t5*t6*t7)
                j=j+1
310             continue
!
!  gamma coefficients (gliste) for l-m even
              j=jlowe
                do 320 i=limdle-1,limd+1,2
                i2=i+i
                t1=m+i-1
                t2=m+i
                t3=m2*m2-1
                t4=m2+i2-3
                t5=m2+i2+1
                gliste(J)=t1*t2-(0.5e0_knd)*c2*((1.0e0_knd)-t3/(t4*t5))
                j=j+1
320             continue
              go to 360
330           continue
!
!  beta coefficients (blist0) for l-m odd
              if(limdlo.gt.limd) go to 360
              j=jlowo
                do 340 i=limdlo,limd,2
                i2=i+i
                t1=i
                t2=i-1
                t3=m2+i
                t4=m2+i-1
                t5=m2+i2-1
                t6=m2+i2-3
                t7=m2+i2+1
                blisto(j)=c4*t1*t2*t3*t4/(t5*t5*t6*t7)
                j=j+1
340             continue
!
!  gamma coefficient (glist0) for l-m odd
              j=jlowo
                do 350 i=limdlo-1,limd+1,2
                i2=i+i
                t1=m+i-1
                t2=m+i
                t3=m2*m2-1
                t4=m2+i2-3
                t5=m2+i2+1
                glisto(J)=t1*t2-(0.5e0_knd)*c2*((1.0e0_knd)-t3/(t4*t5))
              j=j+1
350           continue
360           continue
              itestm=0
              if(ix.eq.0) call conver (l,m,c,limd,bliste,gliste,ndec, &
                                       maxd,ioprad,kindd,kindq,eigval, &
                                       enr,ienre,itestm)
              if(ix.eq.1) call conver (l,m,c,limd,blisto,glisto,ndec, &
                                       maxd,ioprad,kindd,kindq,eigval, &
                                       enr,ienro,itestm)
              eig(li)=eigval
              if(ix.eq.1) go to 370
              limdle=limd+2
              if(2*(limd/2).ne.limd) limdle=limd+1
              jlowe=limd/2+1
              go to 380
370           limdlo=limd+1
              if(2*(limd/2).ne.limd) limdlo=limd+2
              jlowo=(limd-1)/2+1
380           call dnorm (l,m,c,ndec,nex,limd,maxd,enr,ioprad,iopang, &
                          dc01,idc01,dfnorm,idfe,dmlf,idmlfe,dmfnorm, &
                          idmfe,dmlmf,idmlmfe,dmsnorm,idmse,dmlms, &
                          idmlmse,jnorm,jsub,ksub)
              if(l.gt.m+1) then
                if(eig(li).lt.(eig(li-2)- &
                   (1.0e-10_knd)*abs(eig(li-2)))) go to 1540
               end if
              if(ioprad.eq.0) go to 1410
              if(li.lt.listart) go to 385
              if(ioprad.eq.1) go to 385
              if(li.eq.listart.and.jsub.le.ndec-minacc.and.x.le. &
                  0.99e0_knd.and.iopleg.eq.0.and.iopleg1.eq.0.and. &
                  iopint.eq.1) iopleg=1
              if(jsub.le.ndec-minacc.and.x.le.0.99e0_knd.and.iopleg.eq.0 &
                 .and.iopleg1.eq.0.and.l.ge.legstart) iopleg=1
              if(iopleg.eq.1.and.iopint.eq.1) iopint=0
                if(li.eq.listart.and.x.ge.xneu) then
                  if(jsub.le.ndec-minacc) then
                  iopneu0=1
                  else
                  iopneu0=0
                  end if
                end if
              jsubtest=ndec-min(naccrsav,minacc)
              if(li.ne.listart.and.jsub.le.jsubtest.and.x.ge.xneu &
                .and.iopneu0.eq.0.and.(iopleg.eq.0.or.naccleg.lt.minacc) &
                .and.(iopint.eq.0.or.naccint.lt.minacc).and. &
                (iopleg1.eq.0.or.naccleg1.lt.minacc)) iopneu0=4
              if(istartneu0.eq.0) iopneu0=0
              if((li.eq.listart.or.iopeta.ne.0).and.iopneu0.eq.4) &
                  iopneu0=1
              if(li.eq.listart.and.iopneu0.eq.1.and.iopint.eq.1) &
                  iopint=0
385           continue
              if(ioprad.eq.0.or.x.ne.0.0e0_knd) go to 550
!
!  determine oblate radial functions of both kinds when x = 0
              ioppsum=0
              limdr=int(c)+2*ndec
              if(ioprad.eq.2.and.m.ne.0) call dalt(l,m,c,nex,limdr, &
                         maxdr,maxmp,ioppsum,eigval,enrneg,drhor,dneg, &
                         idneg,nsdneg,nsdrhor1)
              if(m.eq.0) dneg=1.0e0_knd
              if(m.eq.0) idneg=0
              nsdneg=0
              naccr=min(ndec-max(jsub,ksub,nsdneg)-1,itestm-1,ndec-2)
              if(naccr.lt.0) naccr=0
              naccr1=min(itestm-1,ndec-ksub-1,ndec-2)
              if(ix.eq.1) go to 420
              if(li.ne.1) fac1=-real(l-m,knd)*(l-m-1)*fac1/ &
                                (real(l+m,knd)*(l+m-1))
              iterm=int(log10(abs(fac1)))
              fac1=fac1*(ten**(-iterm))
              ifac1=ifac1+iterm
              r1c(li)=fac1*dc01/dmfnorm
              ir1e(li)=int(log10(abs(r1c(li))))
              r1c(li)=r1c(li)*(ten**(-ir1e(li)))
              ir1e(li)=ir1e(li)+ifac1+idc01-idmfe
              if(abs(r1c(li)).ge.1.0e0_knd) go to 390
              r1c(li)=r1c(li)*ten
              ir1e(li)=ir1e(li)-1
390           r1dc(li)=0.0e0_knd
              ir1de(li)=0
              if(ioprad.eq.1) go to 450
              r2dc(li)=1.0e0_knd/(c*r1c(li))
              ir2de(li)=int(log10(abs(r2dc(li))))
              r2dc(li)=r2dc(li)*(ten**(-ir2de(li)))
              ir2de(li)=ir2de(li)-ir1e(li)
              if(abs(r2dc(li)).ge.1.0e0_knd) go to 400
              r2dc(li)=r2dc(li)*ten
              ir2de(li)=ir2de(li)-1
400           if(naccr.eq.0) r2c(li)=0.0e0_knd
              if(naccr.eq.0) ir2e(li)=0
              if(li.ne.1) fac2=-real(l+m-1,knd)*(l-m-1)*fac2/ &
                                (real(l-m,knd)*(l+m))
              if(naccr.eq.0) go to 450
              r2c(li)=fac2*dfnorm*dfnorm/(dneg*dc01*dmfnorm)
              ir2e(li)=int(log10(abs(r2c(li))))
              r2c(li)=r2c(li)*(ten**(-ir2e(li)))
              ir2e(li)=ir2e(li)+ifac2-idneg-idc01+idfe+idfe-idmfe
              if(abs(r2c(li)).ge.1.0e0_knd) go to 450
              r2c(li)=r2c(li)*ten
              ir2e(li)=ir2e(li)-1
410           go to 450
420           r1c(li)=0.0e0_knd
              ir1e(li)=0
              if(li.ne.2) fac1d=-real(l-m,knd)*(l-m-1)*fac1d/ &
                                 (real(l+m,knd)*(l+m-1))
              iterm=int(log10(abs(fac1d)))
              fac1d=fac1d*(ten**(-iterm))
              ifac1d=ifac1d+iterm
              r1dc(li)=fac1d*dc01/dmfnorm
              ir1de(li)=int(log10(abs(r1dc(li))))
              r1dc(li)=r1dc(li)*(ten**(-ir1de(li)))
              ir1de(li)=ir1de(li)+ifac1d+idc01-idmfe
              if(abs(r1dc(li)).ge.1.0e0_knd) go to 430
              r1dc(li)=r1dc(li)*ten
              ir1de(li)=ir1de(li)-1
430           if(ioprad.eq.1) go to 450
              r2c(li)=-1.0e0_knd/(c*r1dc(li))
              ir2e(li)=int(log10(abs(r2c(li))))
              r2c(li)=r2c(li)*(ten**(-ir2e(li)))
              ir2e(li)=ir2e(li)-ir1de(li)
              if(abs(r2c(li)).ge.1.0e0_knd) go to 440
              r2c(li)=r2c(li)*ten
              ir2e(li)=ir2e(li)-1
440           if(naccr.eq.0) r2dc(li)=0.0e0_knd
              if(naccr.eq.0) ir2de(li)=0
              if(li.ne.2) fac2d=-real(l-m,knd)*(l+m)*fac2d/ &
                                (real(l+m-1,knd)*(l-m-1))
              if(naccr.eq.0) go to 450
              r2dc(li)=fac2d*dfnorm*dfnorm/(dneg*dc01*dmfnorm)
              ir2de(li)=int(log10(abs(r2dc(li))))
              r2dc(li)=r2dc(li)*(ten**(-ir2de(li)))
              ir2de(li)=ir2de(li)+ifac2d-idneg-idc01+idfe+idfe-idmfe
              if(abs(r2dc(li)).ge.1.0e0_knd) go to 450
              r2dc(li)=r2dc(li)*ten
              ir2de(li)=ir2de(li)-1
450           continue
if (debug) then
              if(ix.eq.0) write(40,460)
460           format(5x,'calculated accurate value for r1 ', &
                     'using nonzero term in traditional Bessel ', &
                     'function expansion')
end if
              if(ix.eq.1) write(40,470)
if (debug) then
470                  format(5x,'calculated accurate value for r1d ', &
                     'using nonzero term in traditional Bessel ', &
                     'function expansion')
              if(knd.eq.kindd) write(40,570) r1c(li),ir1e(li),r1dc(li),ir1de(li)
              if(knd.eq.kindq) write(40,575) r1c(li),ir1e(li),r1dc(li),ir1de(li)
end if
              if(ioprad.eq.1) go to 1330
if (debug) then
              if(ix.eq.0) write(40,500)
500           format(5x,'calculated r2 using nonzero term in Legendre', &
                    ' function expansion and r2d from Wronskian and r1')
              if(ix.eq.1) write(40,510)
510           format(5x,'calculated r2d using nonzero term in ', &
                     'Legendre function expansion and r2 from ', &
                     'Wronskian and r1d')
end if
if (debug) then
              if(knd.eq.kindd) write(40,520) r2c(li),ir2e(li),r2dc(li),ir2de(li)
              if(knd.eq.kindq) write(40,525) r2c(li),ir2e(li),r2dc(li),ir2de(li)
end if
520           format(10x,'r2 = ',f19.15,i5,5x,'r2d = ',f19.15,i5)
525           format(10x,'r2 = ',f34.30,i5,5x,'r2d = ',f34.30,i5)
if (debug) then
              if(ix.eq.0) write(40,530) naccr
530           format(12x,'r2 is accurate to ',I2,' decimal digits. r1,' &
                     ' r1d, and r2d are highly accurate.')
              if(ix.eq.1) write(40,540) naccr
end if
540           format(12x,'r2d is accurate to ',I2,' decimal digits. r1,' &
                     ' r1d, and r2 are highly accurate.')
              go to 1330
550           continue
!
!  determine oblate radial functions of the first kind
!    r1 calculation using traditional Bessel functions series (eta=1)
if (debug) then
              write(40,560)
560           format(4x,'r1 and r1d calculation')
end if
                if(iopbes.eq.0) then
                jbes=jnorm
                go to 580
                end if
              if(li.eq.1) limr1=3*ndec+int(c)
              if(li.ne.1) limr1=jbes+jbes+20+int(sqrt(c))
              limr1=min(limr1,limj-m-2)
              call r1bes(l,m,c,x,limr1,ndec,maxd,enr,maxj,maxn,maxlp, &
                         nex,iflag,sbesf,sbesdf,sbesn,ibese,sbesdr, &
                         prat1,pcoefn,ipcoefn,dmfnorm,idmfe,ir1ep,r11c, &
                         ir11e,r1d1c,ir1d1e,jbes1,nsub,nsubd)
              jbes=jbes1
              iopbes=2
if (debug) then
              if(knd.eq.kindd) write(40,570) r11c,ir11e,r1d1c,ir1d1e
              if(knd.eq.kindq) write(40,575) r11c,ir11e,r1d1c,ir1d1e
570           format(10x,'r1 = ', f19.15,i5,5x,'r1d = ',f19.15,i5)
575           format(10x,'r1 = ', f34.30,i5,5x,'r1d = ',f34.30,i5)
end if
              r1c(li)=r11c
              ir1e(li)=ir11e
              r1dc(li)=r1d1c
              ir1de(li)=ir1d1e
                if(iopeta1.ne.0) then
                  if(nsub.lt.2.or.nsubd.lt.2) then
                  iopeta1=0
                  else
                  iopbes=0
                  iopeta1=2
                  idir=0
                  nee1=1
                  iflageta1=0
                  nsub1p=nsub
                  nsubd1p=nsubd
                  end if
                end if
                iflagbes=0
                if(l.eq.m.and.(nsub.gt.1.or.nsubd.gt.1)) then
                iopeta1=1
                iopbes=0
                nsub1p=nsub
                nsubd1p=nsubd
                end if
              naccr1=ndec-2-max(nsub,nsubd)
580           continue

!
!    r1 calculation using the variable eta expansion
              if(iopeta1.eq.0) go to 760
              lflageta1=1
              if(nee1.eq.1) lflageta1=0
                if(ieta.eq.0) then
                  do jnet=1,neta
                  ang=jnet*pi*0.5e0_knd/(neta+1)
                  eta(jnet)=cos(ang)
                  wmeta2(jnet)=2.0e0_knd*(1.0e0_knd+eta(jnet))* &
                               (sin(0.5e0_knd*ang)**2)
                  xbn(jnet)=sqrt(x*x+wmeta2(jnet))
                  xln(jnet)=eta(jnet)*x/xbn(jnet)
                  end do
                ieta=1
                end if
              if(iopeta1.eq.1) iopeta1=2
590           if(iopeta1.eq.3) go to 690
              etaval1=eta(nee1)
600           xbninp=xbn(nee1)
              netainp=1
              etainp(1)=eta(nee1)
              xlninp(1)=xln(nee1)
              limj=lnum+3*ndec+int(c)+m
              if(limj.lt.maxlp+1) limj=maxlp+1
              if(limj.gt.maxj) limj=maxj-2
              limp=limj-m
              if(jnen1.eq.0) go to 660
              jnenlim1=jnen1
              if(jnen1.gt.jnebmax) jnenlim1=jnebmax
              limplim=limp
              limjlim=limj
                do 650 jn=1,jnenlim1
                if(nee1.ne.neeb1(jn)) go to 650
                if(limplim.gt.limp1sv(jn)) limplim=limp1sv(jn)
                if(limjlim.gt.limjsv(jn)) limjlim=limjsv(jn)
                  do 610 je=1,limplim
                  pratb1(je)=pratbsv1(jn,je)
                  pratt1(je)=prattsv1(jn,je)
                  pdratt1(je)=pdrattsv1(jn,je)
610               continue
                  do 620 je=1,limjlim
                  sbesfe(je)=sbesfsv(jn,je)
                  sbesdfe(je)=sbesdfsv(jn,je)
620               continue
                  jelim=maxlp
                  if(maxlp.gt.limj+1) jelim=limj+1
                  if(jelim1.gt.jelim1sv(jn)) jelim1=jelim1sv(jn)
                  do 630 je=1,jelim1
                  sbesne(je)=sbesnsv(jn,je)
                  sbesdre(je)=sbesdrsv(jn,je)
                  ibesee(je)=ibesesv(jn,je)
630               continue
                go to 680
650             continue
660           continue
              jnen1=jnen1+1
              jnencur1=jnen1-(jnebmax*int((jnen1-1)/jnebmax))
              neeb1(jnencur1)=nee1
              call sphbes(c,xbninp,limj,maxj,maxlp,sbesfe,sbesdfe, &
                          sbesne,ibesee,sbesdre)
                do je=1,limj
                sbesfsv(jnencur1,je)=sbesfe(je)
                sbesdfsv(jnencur1,je)=sbesdfe(je)
                limjsv(jnencur1)=limj
                end do
                jelim1=maxlp
                if(maxlp.gt.limj+1) jelim1=limj+1
                  do 670 je=1,jelim1
                  sbesnsv(jnencur1,je)=sbesne(je)
                  sbesdrsv(jnencur1,je)=sbesdre(je)
                  ibesesv(jnencur1,je)=ibesee(je)
670               continue
                  jelim1sv(jnencur1)=jelim1
              iopd=3
              call pleg(m,limp,maxp,limcsav,iopd,ndec,nex,xlninp, &
                        netainp,maxt,prat,pdrat,pdnorma,ipdnorma,pnorma, &
                        ipnorma,alpha,beta,gamma,coefa,coefb,coefc, &
                        coefd,coefe)
              limcsav=max(limcsav,limp)
                do je=1,limp
                pratt1(je)=prat(1,je)
                pdratt1(je)=pdrat(1,je)
                prattsv1(jnencur1,je)=pratt1(je)
                pdrattsv1(jnencur1,je)=pdratt1(je)
                limp1sv(jnencur1)=limp
                end do
              limpd=2*(lnum+int(c)+ndec)
              if(limpd.gt.limp) limpd=limp
              iopd=2
              call pleg(m,limpd,maxp,limcsav,iopd,ndec,nex,etainp, &
                        netainp,maxt,prat,pdrat,pdnorma,ipdnorma,pnorma, &
                        ipnorma,alpha,beta,gamma,coefa,coefb,coefc, &
                        coefd,coefe)
                do je=1,limpd
                pratb1(je)=prat(1,je)
                pratbsv1(jnencur1,je)=pratb1(je)
                end do
              pratb1(limpd+1)=0.0e0_knd
              pratb1(limpd+2)=0.0e0_knd
              pratbsv1(jnencur1,limpd+1)=0.0e0_knd
              pratbsv1(jnencur1,limpd+2)=0.0e0_knd
680           continue
              pcoefe1=(x*x+1.0e0_knd)/(xbn(nee1)*xbn(nee1))
              apcoef1=(rm/2.0e0_knd)*log10(pcoefe1)
              ipcoefe1=int(apcoef1)
              pcoefe1=ten**(apcoef1-ipcoefe1)
              pcoefo1=pcoefe1*pratt1(2)/pratb1(2)
              ipcoefo1=ipcoefe1
              pdcoefe1=pcoefe1
              if(m.ne.0) pdcoefe1=-pcoefe1*rm*xln(nee1)*xbn(nee1)* &
                          xbn(nee1)/((x*x+1.0e0_knd)*wmeta2(nee1))
              ipdcoefe1=ipcoefe1
              pdcoefo1=pdcoefe1*pdratt1(2)/pratb1(2)
              ipdcoefo1=ipdcoefe1
              if(li.eq.1) go to 690
                do jl=3,li+ix,2
                pcoefe1=pcoefe1*pratt1(jl)/pratb1(jl)
                iterm=log10(abs(pcoefe1))
                pcoefe1=pcoefe1*ten**(-iterm)
                ipcoefe1=ipcoefe1+iterm
                pdcoefe1=pdcoefe1*pdratt1(jl)/pratb1(jl)
                iterm=log10(abs(pdcoefe1))
                pdcoefe1=pdcoefe1*ten**(-iterm)
                ipdcoefe1=ipdcoefe1+iterm
                end do
              continue
              if(li.lt.3) go to 690
                do jl=4,li+1-ix,2
                pcoefo1=pcoefo1*pratt1(jl)/pratb1(jl)
                iterm=log10(abs(pcoefo1))
                pcoefo1=pcoefo1*ten**(-iterm)
                ipcoefo1=ipcoefo1+iterm
                pdcoefo1=pdcoefo1*pdratt1(jl)/pratb1(jl)
                iterm=log10(abs(pdcoefo1))
                pdcoefo1=pdcoefo1*ten**(-iterm)
                ipdcoefo1=ipdcoefo1+iterm
                end do
690           continue
              if(ix.eq.0) go to 700
              pcoefet1=pcoefo1
              ipcoefet1=ipcoefo1
              pcoefo1=pcoefo1*pratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pcoefo1)))
              pcoefo1=pcoefo1*ten**(-iterm)
              ipcoefo1=ipcoefo1+iterm
              pdcoefet1=pdcoefo1
              ipdcoefet1=ipdcoefo1
              pdcoefo1=pdcoefo1*pdratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pdcoefo1)))
              pdcoefo1=pdcoefo1*ten**(-iterm)
              ipdcoefo1=ipdcoefo1+iterm
              go to 710
700           pcoefet1=pcoefe1
              ipcoefet1=ipcoefe1
              pcoefe1=pcoefe1*pratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pcoefe1)))
              pcoefe1=pcoefe1*ten**(-iterm)
              ipcoefe1=ipcoefe1+iterm
              pdcoefet1=pdcoefe1
              ipdcoefet1=ipdcoefe1
              pdcoefe1=pdcoefe1*pdratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pdcoefe1)))
              pdcoefe1=pdcoefe1*ten**(-iterm)
              ipdcoefe1=ipdcoefe1+iterm
710           continue
              wm=wmeta2(nee1)
              limeta=l+3*ndec+int(c)
              if(iopeta1.eq.3) limeta=jeta1+jeta1+20
              if(limeta.gt.limp-2) limeta=limp-2
              if(limeta.gt.limd) limeta=limd
              call r1eta(l,m,c,x,etaval1,nee1,limeta,ndec,nex,maxd, &
                         maxlp,maxj,maxp,minacc,wm,enr,sbesfe,sbesne, &
                         ibesee,sbesdfe,sbesdre,pdratt1,pratb1,pratt1, &
                         pcoefet1,ipcoefet1,pdcoefet1,ipdcoefet1,ir1ep, &
                         r1ec,ir1ee,r1dec,ir1dee,nsub1,nsubd1,jeta1)
              if(iopbes.eq.0) jbes=jeta1
              if(iopbes.ne.0) jbes=max(jbes,jeta1)
if (debug) then
720           if(knd.eq.kindd) write(40,730) etaval1,nee1,r1ec,ir1ee,r1dec,ir1dee
              if(knd.eq.kindq) write(40,735) etaval1,nee1,r1ec,ir1ee,r1dec,ir1dee
730           format(15x,'eta = ',f12.9,'; nee1 = ',i4,/,10x,'r1 = ', &
                      f19.15,i5,5x,'r1d = ',f19.15,i5)
735                   format(15x,'eta = ',f12.9,'; nee1 = ',i4,/,10x,'r1 = ', &
                      f34.30,i5,5x,'r1d = ',f34.30,i5)
end if
                if(nsub1.le.1.or.nsubd1.le.1) then
                if(idir.eq.0) idir=-1
                iflageta1=0
                nsub1p=nsub1
                nsubd1p=nsubd1
                r1c(li)=r1ec
                ir1e(li)=ir1ee
                r1dc(li)=r1dec
                ir1de(li)=ir1dee
                if(jflageta1.ne.0) jflageta1=jflageta1-1
                if(jflageta1.eq.0) iopeta1=3
                go to 750
                end if
                if(idir.eq.0) then
                  if(nee1.eq.1) then
                    if((nsub1.ge.nsub1p.and.nsubd1.ge.nsubd1p).and. &
                        (nsub1p.lt.4.or.nsubd1p.lt.4)) then
                    iopbes=1
                    iopeta1=2
                    go to 760
                    end if
                  r1cm=r1ec
                  ir1em=ir1ee
                  r1dcm=r1dec
                  ir1dem=ir1dee
                  nsub1m=nsub1
                  nsubd1m=nsubd1
                  nee1m=1
                  end if
                  if(nee1.ne.1.and.((nsub1.lt.nsub1m.and.nsubd1.le. &
                       nsubd1m).or.(nsubd1.lt.nsubd1m.and.nsub1.le. &
                       nsub1m))) then
                  r1cm=r1ec
                  ir1em=ir1ee
                  r1dcm=r1dec
                  ir1dem=ir1dee
                  nsub1m=nsub1
                  nsubd1m=nsubd1
                  nee1m=nee1
                  end if
                  if(nee1.eq.neta.or.((nsub1.gt.nsub1m+1.or. &
                     nsubd1.gt.nsubd1m+1).and.(nsub1m.lt.4.or. &
                     nsubd1m.lt.4))) then
                  r1c(li)=r1cm
                  ir1e(li)=ir1em
                  r1dc(li)=r1dcm
                  ir1de(li)=ir1dem
                  nsub1=nsub1m
                  nsubd1=nsubd1m
                  nsub1p=nsub1
                  nsubd1p=nsubd1
                  nee1=nee1m
                  iopeta1=2
                  idir=-1
                  go to 750
                  end if
                nee1=nee1+incnee1
                iopeta1=2
                nsub1p=nsub1
                nsubd1p=nsubd1
                go to 590
                end if
740             if((nsub1.gt.nsub1p.or.nsubd1.gt.nsubd1p.or.(nsub1.eq. &
                    nsub1p.and.nsubd1.eq.nsubd1p)).and.iflageta1.eq.1) &
                    then
                nee1=nee1+incnee1
                iflageta1=0
                nsub1=nsub1p
                nsubd1=nsubd1p
                iopeta1=2
                jflageta1=2
                go to 750
                end if
              r1c(li)=r1ec
              ir1e(li)=ir1ee
              r1dc(li)=r1dec
              ir1de(li)=ir1dee
              nsub1p=nsub1
              nsubd1p=nsubd1
              iflageta1=1
              nee1=nee1-incnee1
              iopeta1=2
                if(nee1.eq.0.and.l.le.nbp.and.nsub1.gt.2.and. &
                       nsubd1.gt.2.and.lflageta1.eq.0) then
                nee1=2
                lflageta1=1
                iflageta1=0
                go to 590
                end if
                if(nee1.eq.0.and.l.le.nbp.and.(nsub1.le.2.or. &
                       nsubd1.le.2.or.lflageta1.eq.1)) then
                nee1=1
                iopeta1=3
                iflageta1=0
                go to 750
                end if
                if(nee1.eq.0.and.l.gt.nbp) then
                iopbes=1
                iopeta1=2
                go to 750
                end if
              go to 590
750           continue
              naccr1=ndec-2-max(nsub1,nsubd1)
                if(nsub1.gt.2.and.nsubd1.gt.2.and.l.gt.nbp) then
                iopbes=1
                iopeta1=2
                end if
760           continue
              if(ioprad.eq.1) go to 1330
!
!  determine oblate radial functions of the second kind
!
if (debug) then
              write(40,770)
770 format(4x,'r2 and r2d calculation')
end if
!
!  decide whether to use values of r1 and r1d for r2 and r2d
              naccrt=0
                if(li.lt.matlilim) then
                  if(li.ge.listart) then
                  if(ix.eq.0) naccrt=min(match(li+1)-2,naccr1, &
                                     itestm-2)
                    if(ix.eq.1.and.x.ge.0.01e0_knd) then
                    wronc=-r1c(li)*r1dc(li-1)*ten**(ir1e(li)+ &
                           ir1de(li-1))+r1c(li-1)*r1dc(li)* &
                           ten**(ir1e(li-1)+ir1de(li))
                    naccrt=-int(log10(abs((wronc-wront)/wront)+dec))
                    if(naccrt.gt.ndec) naccrt=ndec
                    if(naccrt.lt.0) naccrt=0
                    end if
                    if(ix.eq.1.and.x.lt.0.01e0_knd) then
                    nmatch=-int(log10(abs((eig(li)-eig(li-1))/(eig(li)) &
                        +ten*dec)))
                    if(nmatch.gt.ndec) nmatch=ndec
                    naccrt=min(nmatch-2,naccr1,itestm-2)
                    end if
                  end if
                naccr=naccrt
                  if(li.lt.listart) then
                    if(ix.eq.0) then
if (debug) then
                    write(40,1360) l,l+1
end if
                    naccrt=min(match(li+1)-2,naccr1, &
                                     itestm-2)
                    naccr=naccrt
                    end if
                    if(ix.eq.1) then
if (debug) then
                    write(40,1370) l,l-1
end if
                    end if
                  go to 1390
                  end if
                if(x.le.0.99e0_knd.and.iopleg.eq.0.and.l.ge.legstart) &
                    iopleg=1
                if(x.ge.xneu.and.x.le.0.9e0_knd.and.x.gt.0.05e0_knd.and. &
                    jsub.gt.ndec-minacc.and.iopleg.eq.0.and. &
                    iopeta.eq.0) iopeta=4
                if(li.eq.listart.and.iopeta.eq.4) iopeta=1
                if(x.gt.0.99e0_knd.and.iopneu0.eq.0) iopneu0=4
                if(li.eq.listart.and.iopneu0.eq.4) iopneu0=1
                end if
!
              r1cin=r1c(li)
              ir1ein=ir1e(li)
              r1dcin=r1dc(li)
              ir1dein=ir1de(li)
780           continue
!
!     r2 calculation using integration technique
                if(li.gt.intlim-5) then
                istartint=2
                iopint=0
                end if
              if(iopint.eq.0) go to 830
              if(iopint.gt.1) go to 800
              if(istartint.eq.1) istartint=2
              if(istartint.eq.0) istartint=3
              limint=2*(nbp+33)+2
              if(igau.eq.1) go to 790
              call gauss(ndec,ngau,xr,wr)
              igau=1
790           if(iopint.eq.1.and.ipint.eq.1) iopint=2
              if(ipint.eq.1) go to 800
              ngqs=nstep1+nstep2+nstep3
              call pint(c,m,lnum,x,limint,maxint,maxlp,maxmp,ndec,nex, &
                        ngqs,ngau,wr,xr,step1,nstep1,step2,nstep2,step3, &
                        nstep3,intlim,rpint1,rpint2,pint1,pint2,pint3, &
                        pint4,norme,pnormint,ipnormint,coefme,coefmo)
              iopint=2
800           continue
              if(iopint.eq.2) limint=2*(nbp+33)
              if(iopint.eq.3) limint=jint+jint+20+int(sqrt(c))
              if(limint.gt.maxint) limint=maxint-5
              if(2*(limint/2).ne.limint) limint=limint-1
              call r2int(l,m,c,x,limint,ndec,nex,maxd,enr,dc01,idc01, &
                         maxint,maxmp,maxlp,intlim,rpint1,rpint2, &
                         pint1,pint2,pint3,pint4,norme,pnormint, &
                         ipnormint,coefme,coefmo,ipint,r2ic,ir2ie, &
                         r2dic,ir2die,jint,coefn,icoefn,isub,isubd)
              if(ipint.eq.0) ipint=1
              naccrs=naccr
              wronca=r1c(li)*r2dic*ten**(ir1e(li)+ir2die)
              wroncb=r2ic*r1dc(li)*ten**(ir2ie+ir1de(li))
              wronc=wronca-wroncb
              naccint=-int(log10(abs((wronc-wront)/wront)+dec))
              if(naccint.lt.0) naccint=0
              if(naccint.gt.ndec) naccint=ndec
              jacc=int(log10(abs(wronca/wroncb)))
                if(jacc.gt.0) then
                isubc=isub-jacc
                naccint=naccint+max(isubc,isubd)-max(isub,isubd)
                end if
                if(jacc.lt.0) then
                isubdc=isubd+jacc
                naccint=naccint+max(isub,isubdc)-max(isub,isubd)
                end if
              if(naccint.lt.0) naccint=0
              iopint=3
              if(naccint.lt.naccrs) go to 810
              naccr=naccint
              r2c(li)=r2ic
              ir2e(li)=ir2ie
              r2dc(li)=r2dic
              ir2de(li)=ir2die
810           continue
                if(naccint.gt.minacc) then
                iflagint=1
                iflagneu=0
                iopeta=0
                iopleg1=0
                end if
                if(naccint.le.minacc) then
                if(x.le.0.99e0_knd.and.ndec-jsub.ge.max(naccint, &
                   naccrs).and.l.ge.legstart.and.iopleg.eq.0) iopleg=1
                if(x.ge.xneu) iflagneu=1
                if(iflagneu.eq.1.and.(ndec-jsub).ge.minacc.and. &
                   iopneu0.eq.0) iopneu0=4
                if(iopneu0.eq.4.and.(l.eq.listart.or.iopeta.ne.0)) &
                   iopneu0=1
                if(iflagneu.eq.1.and.iopneu0.eq.0.and.iopeta.eq.0 &
                   .and.iflageta.eq.0.and.x.gt.0.05e0_knd) iopeta=4
                if(iopeta.eq.4.and.l.eq.listart) iopeta=1
                end if
if (debug) then
              if(knd.eq.kindd) write(40,820) naccint,r2ic,ir2ie,r2dic,ir2die
              if(knd.eq.kindq) write(40,825) naccint,r2ic,ir2ie,r2dic,ir2die
820           format(15x,'accuracy in decimal digits = ',i2,/,10x, &
                     'r2 = ',f19.15,i5,5x,'r2d = ',f19.15,i5)
825                  format(15x,'accuracy in decimal digits = ',i2,/,10x, &
                     'r2 = ',f34.30,i5,5x,'r2d = ',f34.30,i5)
end if
              if(istartint.eq.2) go to 1320
830           continue
!
!     r2 calculation using a Legendre expansion and joining factor
              if(iopleg.eq.0.or.ndec-jsub.lt.naccr-2) go to 940
              if(jflagleg.eq.1) go to 900
              jflagleg=1
              limdr=c+2*ndec+50*x+2
              if(limdr.gt.maxdr) limdr=maxdr
              if(ioppsum.eq.0) go to 840
              xin(1)=x
              limpleg=limdr+limdr
              iopd=4
              call pleg(m,limpleg,maxp,limcsav,iopd,ndec,nex,xin,1,maxt, &
                        prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma, &
                        alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
              limcsav=max(limcsav,limpleg)
                do jj=1,limpleg
                prx(jj)=prat(1,jj)
                pdrx(jj)=pdrat(1,jj)
                end do
840             if(iqlegflag.eq.0) then
                limq=lnum+3*ndec+int(c)
                if(c.le.60.0e0_knd.and.c.ge.10.0e0_knd.and.m.le.40.and. &
                   x.le.0.99e0_knd) limq=max(limq,250-int(50*x))
                call qleg(m,lnum,limq,maxmp,maxq,x,ndec,nex,qdrat,qdqr, &
                     qdm1,iqdm1,qdl,iqdl,qrat,qm1,iqm1,ql,iql,termpq, &
                     itermpq,qrat1,qdrat1,qm0,qdm0)
                iqlegflag=1
                end if
              fajo(1)=c/(rm2-1.0e0_knd)
              ifajo(1)=0
              if(m.eq.0) go to 870
                do im=1,m
                fajo(1)=fajo(1)*(im+im)/c
                if(abs(fajo(1)).lt.1.0e+10_knd) go to 850
                fajo(1)=fajo(1)*(1.0e-10_knd)
                ifajo(1)=ifajo(1)+10
850             continue
                if(abs(fajo(1)).gt.1.0e-10_knd) go to 860
                fajo(1)=fajo(1)*(1.0e+10_knd)
                ifajo(1)=ifajo(1)-10
860             continue
                end do
870           continue
              if(lnum.eq.1) go to 900
              fajo(2)=-c*fajo(1)/(rm2-3.0e0_knd)
              ifajo(2)=ifajo(1)
              if(lnum.eq.2) go to 900
                do jl=3,lnum-1,2
                fajo(jl)=fajo(jl-2)*(jl+m+m-1)/(jl-2)
                ifajo(jl)=ifajo(jl-2)
                  if(abs(fajo(jl)).gt.1.0e10_knd) then
                  fajo(jl)=fajo(jl)*1.0e-10_knd
                  ifajo(jl)=ifajo(jl)+10
                  end if
                fajo(jl+1)=fajo(jl-1)*(jl+m+m-1)/(jl)
                ifajo(jl+1)=ifajo(jl-1)
                  if(abs(fajo(jl+1)).gt.1.0e10_knd) then
                  fajo(jl+1)=fajo(jl+1)*1.0e-10_knd
                  ifajo(jl+1)=ifajo(jl+1)+10
                  end if
                end do
                if(2*(lnum/2).ne.lnum.and.lnum.ge.4) then
                fajo(lnum)=fajo(lnum-2)*(lnum+m+m-1)/(lnum-2)
                ifajo(lnum)=ifajo(lnum-2)
                end if
900           continue
              limleg=l-m+3*ndec+int(c)
              limdr=int(c)+2*ndec+int(50*x)
              if(iopleg.eq.2) limleg=jleg+jleg+20+int(sqrt(c))
              if(iopleg.eq.2) limdr=jlegp+10+int(sqrt(c)/2)
              if(limdr.gt.maxdr-4) limdr=maxdr-4
              if(limleg.gt.limq-4) limleg=limq-4
              nsdrhor1=0
              nsdneg=0
              call dalt(l,m,c,nex,limdr,maxdr,maxmp,ioppsum,eigval, &
                        enrneg,drhor,dneg,idneg,nsdneg,nsdrhor1)
              kflagl=0
              ifsub=max(jsub,ksub,nsdneg)
              call r2leg(l,m,c,x,lnum,minacc,limleg,limdr,iflagp,ndec, &
                         nex,maxd,maxmp,maxpdr,maxdr,maxq,enr,enrneg, &
                         drhor,nsdrhor1,dc01,idc01,dneg,idneg,nsdneg, &
                         dfnorm,idfe,dmfnorm,idmfe,prx,pdrx,qdrat,qdqr, &
                         qdm1,iqdm1,qdl,iqdl,qrat,qm1,iqm1,ql,iql,fajo, &
                         ifajo,ifsub,jsub,termpq,itermpq,ioppsum,iopqnsum, &
                         r1cin,ir1ein,r1dcin,ir1dein,naccr1,itestm,r2lc, &
                         ir2le,r2dlc,ir2dle,jleg,jlegp,jflagl,naccleg, &
                         kflagl,isub,isubd)
              iopleg=2
              nacclegs=naccleg
              if(kflagl.eq.1) nacclegs=-nex+10
              if(naccleg.lt.0) naccleg=0
              if(naccleg.gt.ndec) naccleg=ndec
              if(naccleg.lt.naccr.or.(naccleg.eq.naccr.and.naccint.eq. &
                    naccr.and.jflagl.eq.1)) go to 910
              naccr=naccleg
              r2c(li)=r2lc
              ir2e(li)=ir2le
              r2dc(li)=r2dlc
              ir2de(li)=ir2dle
910           continue
              if(naccleg.gt.minacc+1.and.nacclegp.gt.minacc+1) &
                   istartint=3
                if(naccleg.lt.minacc) then
                if(m.le.40.and.c.le.60.0e0_knd.and.c.ge.10.0e0_knd.and. &
                   li.lt.nbp.and.iopleg1.eq.0.and.naccr.lt.minacc.and. &
                   x.gt.0.1e0_knd) iopleg1=1
                if(iopleg1.eq.0.and.x.ge.xneu.and.iopneu0.eq.0.and. &
                   istartneu0.eq.1) iopneu0=4
                  if(iopneu0.eq.4.and.(li.eq.listart.or.iopeta.ne.0)) &
                     then
                  iopneu0=1
                  iopeta=0
                  end if
                end if
                if(naccleg.lt.naccrsav.and.naccleg.lt.minacc) then
                if((iopint.ne.0.and.naccint.gt.naccleg+2.and. &
                    naccint.gt.nacclegp+2).or.iopneu0.ne.0) iopleg=0
                if(iopint.eq.0) itest=naccrsav
                if(iopint.ne.0) itest=min(naccint,naccrsav)
                  if(iopleg1.eq.0) then
                  legstart=l+itest-nacclegs
                  if(legstart.lt.l+1) legstart=l+1
                  if(legstart.eq.l+1.and.iopleg.eq.0) iopleg=1
                  end if
                if(iopleg1.ne.0) legstart=l+minacc-nacclegs
                end if
                if(naccleg.gt.minacc.and.nacclegp.gt.minacc) then
                iopleg1=0
                iopneu0=0
                istartneu0=0
                iopeta=0
                iopint=0
                end if
                if(naccleg.eq.minacc) then
                  if(iopneu0.ne.0.and.istartneu0.eq.1) then
                  iopneu0=4
                  iopeta=0
                  end if
                end if
if (debug) then
920           if(knd.eq.kindd) write(40,820) naccleg,r2lc,ir2le,r2dlc,ir2dle
              if(knd.eq.kindq) write(40,825) naccleg,r2lc,ir2le,r2dlc,ir2dle
end if
              nacclegp=naccleg
940           continue
!
!     r2 calculation using the Baber and Hasse Legendre expansion
              if(iopleg1.eq.0) go to 960
                if(iqlegflag.eq.0) then
                limq=lnum+3*ndec+int(c)
                if(c.le.60.0e0_knd.and.c.ge.10.0e0_knd.and.m.le.40.and. &
                   x.le.0.99e0_knd) limq=max(limq,250-int(50*x))
                call qleg(m,lnum,limq,maxmp,maxq,x,ndec,nex,qdrat,qdqr, &
                          qdm1,iqdm1,qdl,iqdl,qrat,qm1,iqm1,ql,iql, &
                          termpq,itermpq,qrat1,qdrat1,qm0,qdm0)
                iqlegflag=1
                end if
              if(iopleg1.eq.1) limleg1=250-int(50*x)-4
              if(iopleg1.ne.1) limleg1=jleg1+20
              if(limleg1.gt.250-int(50*x)-4) limleg1=250- &
                                         int(50*x)-4
              call r2leg1(l,m,c,x,limleg1,maxq,ndec,eigval,qrat1, &
                          qdrat1,qm0,qdm0,r1cin,ir1ein,r2l1c,ir2l1e, &
                          r2dl1c,ir2dl1e,jleg1)
              wronc=r1c(li)*r2dl1c*ten**(ir1e(li)+ir2dl1e)- &
                    r2l1c*r1dc(li)*ten**(ir2l1e+ir1de(li))
              naccleg1=-int(log10(abs((wronc-wront)/wront)+dec))
              if(naccleg1.lt.0) naccleg1=0
              if(naccleg1.gt.ndec) naccleg1=ndec
              if(iopleg.eq.2.and.naccleg.ge.naccleg1) iopleg1=0
                if(naccleg1.lt.minacc+1.and.iopleg.eq.0) then
                iopleg=1
                end if
                if(naccleg1.ge.minacc.and.iopneu0.ne.0) iopneu0=4
                if(naccleg1.gt.naccr) then
                r2c(li)=r2l1c
                ir2e(li)=ir2l1e
                r2dc(li)=r2dl1c
                ir2de(li)=ir2dl1e
                naccr=naccleg1
                iopleg1=2
                end if
if (debug) then
              if(knd.eq.kindd) write(40,820) naccleg1,r2l1c,ir2l1e,r2dl1c,ir2dl1e
              if(knd.eq.kindq) write(40,825) naccleg1,r2l1c,ir2l1e,r2dl1c,ir2dl1e
end if
960           continue
!
!     r2 calculation using Neumann function expansion with eta set
!     equal to 0
              if(iopneu0.eq.0.or.ndec-jsub.lt.naccr) go to 1080
              if(iopneu0.gt.1) go to 1050
              if(ibflag2.eq.0) go to 1040
                if(knd.eq.kindd) then
                if(x.ge.0.01e0_knd) limn=2*int(25/(x*x)+300/x+ &
                                         3*c+1250*knd)+2
                if(x.ge.0.1e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*maxm+ &
                                        200)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*maxm+ &
                                        200)/x)+2
                if(x.ge.1.0e0_knd) limn=2*int(lnum+c/5+0.5e0_knd*maxm+ &
                                        200)+2
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.01e0_knd) limn=2*int(25/(x*x)+400/x+ &
                                         3*c+1250*knd)+2
                if(x.ge.0.1e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*maxm+ &
                                        350)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*maxm+ &
                                        300)/x)+2
                if(x.ge.1.0e0_knd) limn=2*int(lnum+c/5+0.5e0_knd*maxm+ &
                                        300)+2
                end if
              limn=limn+maxm
              if(limn.gt.maxn) limn=maxn-2
              call sphneu(c,xb,limn,maxn,maxlp,sneuf2,sneun2,ineue2, &
                          sneudf2,sneudr2)
              ibflag2=0
1040          continue
              iopneu0=iopneu0+1
1050          continue
              if(iopneu0.eq.4) go to 1080
                if(knd.eq.kindd) then
                if(x.ge.0.01e0_knd) limneu0=2*int(25/(x*x)+300/x+3*c+ &
                                            1250*knd)
                if(x.ge.0.1e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                           200)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                           200)/x)
                if(x.ge.1.0e0_knd) limneu0=2*int(l-m+c/5+0.5e0_knd*m+ &
                                           200)
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.01e0_knd) limneu0=2*int(25/(x*x)+400/x+3*c+ &
                                            1250*knd)
                if(x.ge.0.1e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                           350)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+ &
                                           300)/x)
                if(x.ge.1.0e0_knd) limneu0=2*int(l-m+c/5+.4e0_knd*m+300)
                end if
              if(iopneu0.eq.3.and.jtest0.ge.minacc) &
                  limneu0=jneu0max+jneu0max+40+ &
                          int(sqrt(c)*(2/min(1.0e0_knd,x))+100/x)
              limneu0=min(limneu0,limd-2)
              iopneu0=3
              call r2neu0(l,m,c,x,limneu0,ndec,nex,maxd,maxlp,maxn, &
                          minacc,enr,sneuf2,sneun2,ineue2,sneudf2, &
                          sneudr2,dfnorm,idfe,r1dcin,ir1dein,r2nc,ir2ne, &
                          r2dnc,ir2dne,jneu0,jtest0)
              jneu0max=max(jneu0s,jneu0)
              jneu0s=jneu0
              naccrs=naccr
              wronc=r1c(li)*r2dnc*ten**(ir1e(li)+ir2dne)- &
                    r2nc*r1dc(li)*ten**(ir2ne+ir1de(li))
              naccneu0=-int(log10(abs((wronc-wront)/wront)+dec))
              if(naccneu0.lt.0) naccneu0=0
              if(naccneu0.gt.ndec) naccneu0=ndec
                if(naccneu0.ge.minacc.and.naccneu0p.ge.minacc.and. &
                   naccneu0p2.ge.minacc) then
                iopint=0
                istartint=3
                iopleg=0
                iopleg1=0
                iopeta=0
                end if
              if(naccneu0.lt.minacc.and.iopeta.eq.0.and.jetaflag.eq.1 &
                   .and.ndec-jsub.lt.minacc+1.and.x.gt.0.05e0_knd) &
                   iopeta=3
              if(naccneu0.lt.naccrs) go to 1060
              naccr=naccneu0
              r2c(li)=r2nc
              ir2e(li)=ir2ne
              r2dc(li)=r2dnc
              ir2de(li)=ir2dne
1060          continue
              naccneu0p2=naccneu0p
              naccneu0p=naccneu0
if (debug) then
              if(knd.eq.kindd) write(40,820) naccneu0,r2nc,ir2ne,r2dnc,ir2dne
              if(knd.eq.kindq) write(40,825) naccneu0,r2nc,ir2ne,r2dnc,ir2dne
end if
1080          continue
              if(iopneu0.eq.4) iopneu0=1
!
!     r2 calculation using the variable eta expansion
              if(iopeta.eq.0) go to 1310
              iopnee=0
              naccetamax=0
              neemax=nee
              naccnmax=0
              netatry=1
              naccd=0
              jetam=0
                if(ieta.eq.0) then
                  do jnet=1,neta
                  ang=jnet*pi*0.5e0_knd/(neta+1)
                  eta(jnet)=cos(ang)
                  wmeta2(jnet)=2.0e0_knd*(1.0e0_knd+eta(jnet))* &
                               (sin(0.5e0_knd*ang)**2)
                  xbn(jnet)=sqrt(x*x+wmeta2(jnet))
                  xln(jnet)=eta(jnet)*x/xbn(jnet)
                  end do
                  ieta=1
                end if
              if(iopeta.eq.1) iopeta=2
1090          if(iopeta.eq.4) go to 1320
              if(iopeta.eq.3) go to 1180
              etaval=eta(nee)
1100          xbninp=xbn(nee)
              netainp=1
              etainp(1)=eta(nee)
              xlninp(1)=xln(nee)
                if(knd.eq.kindd) then
                if(x.gt.0.05e0_knd) limn=2*int(25/(x*x)+300/x+ &
                                         3*c+1250*knd)+2
                if(x.ge.0.1e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*m+ &
                                        200)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*m+ &
                                        300)/x)+2
                if(x.ge.1.0e0_knd) limn=2*int(lnum+c/5+0.5e0_knd*m+ &
                                        300)+2
                end if
                if(knd.eq.kindq) then
                if(x.gt.0.05e0_knd) limn=2*int(25/(x*x)+400/x+ &
                                         3*c+1250*knd)+2
                if(x.ge.0.1e0_knd) limn=2*((lnum+c/5+0.5e0_knd*m+ &
                                        350)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*((lnum+c/5+0.5e0_knd*m+ &
                                        400)/x)+2
                if(x.ge.1.0e0_knd) limn=2*(lnum+c/5+0.5e0_knd*m+400) &
                                        +2
                end if
              limn=limn+m
              if(limn.gt.maxn) limn=maxn-2
              limp=limn-m
              if(jnen.eq.0) go to 1160
              jnenlim=jnen
              if(jnen.gt.jnenmax) jnenlim=jnenmax
              limplim=limp
              limnlim=limn
                do 1150 jn=1,jnenlim
                if(nee.ne.neeb(jn)) go to 1150
                if(limplim.gt.limpsv(jn)) limplim=limpsv(jn)
                if(limnlim.gt.limnsv(jn)) limnlim=limnsv(jn)
                  do 1110 je=1,limplim
                  pratb(je)=pratbsv(jn,je)
                  pratt(je)=prattsv(jn,je)
                  pdratt(je)=pdrattsv(jn,je)
1110              continue
                  do 1120 je=1,limnlim
                  sneufe(je)=sneufsv(jn,je)
                  sneudfe(je)=sneudfsv(jn,je)
1120              continue
                  jelim=maxlp
                  if(maxlp.gt.limn+1) jelim=limn+1
                  if(jelim.gt.jelimsv(jn)) jelim=jelimsv(jn)
                  do 1130 je=1,jelim
                  sneune(je)=sneunsv(jn,je)
                  sneudre(je)=sneudrsv(jn,je)
                  ineuee(je)=ineuesv(jn,je)
1130              continue
if (debug) then
                write(40,1140) etaval
1140            format(8x,'r2eta: reused expansion functions for eta =' &
                       ,f13.9,'.')
end if
                go to 1180
1150           continue
1160          continue
              jnen=jnen+1
              jnencur=jnen-(jnenmax*int((jnen-1)/jnenmax))
              neeb(jnencur)=nee
              call sphneu(c,xbninp,limn,maxn,maxlp,sneufe,sneune, &
                          ineuee,sneudfe,sneudre)
                do je=1,limn
                sneufsv(jnencur,je)=sneufe(je)
                sneudfsv(jnencur,je)=sneudfe(je)
                limnsv(jnencur)=limn
                end do
                jelim=maxlp
                if(maxlp.gt.limn+1) jelim=limn+1
                  do 1170 je=1,jelim
                  sneunsv(jnencur,je)=sneune(je)
                  sneudrsv(jnencur,je)=sneudre(je)
                  ineuesv(jnencur,je)=ineuee(je)
1170              continue
                  jelimsv(jnencur)=jelim
              iopd=3
              call pleg(m,limp,maxp,limcsav,iopd,ndec,nex,xlninp, &
                        netainp,maxt,prat,pdrat,pdnorma,ipdnorma,pnorma, &
                        ipnorma,alpha,beta,gamma,coefa,coefb,coefc, &
                        coefd,coefe)
              limcsav=max(limcsav,limp)
                do je=1,limp
                pratt(je)=prat(1,je)
                pdratt(je)=pdrat(1,je)
                prattsv(jnencur,je)=pratt(je)
                pdrattsv(jnencur,je)=pdratt(je)
                limpsv(jnencur)=limp
                end do
              limpd=2*(lnum+int(c)+ndec)
              if(limpd.gt.limp) limpd=limp
              iopd=2
              call pleg(m,limpd,maxp,limcsav,iopd,ndec,nex,etainp, &
                        netainp,maxt,prat,pdrat,pdnorma,ipdnorma,pnorma, &
                        ipnorma,alpha,beta,gamma,coefa,coefb,coefc, &
                        coefd,coefe)
                do je=1,limpd
                pratb(je)=prat(1,je)
                pratbsv(jnencur,je)=pratb(je)
                end do
              pratb(limpd+1)=0.0e0_knd
              pratb(limpd+2)=0.0e0_knd
              pratbsv(jnencur,limpd+1)=0.0e0_knd
              pratbsv(jnencur,limpd+2)=0.0e0_knd
1180          continue
              pcoefe=(x*x+1.0e0_knd)/(xbn(nee)*xbn(nee))
              apcoef=(rm/2.0e0_knd)*log10(pcoefe)
              ipcoefe=int(apcoef)
              pcoefe=ten**(apcoef-ipcoefe)
              pcoefo=pcoefe*pratt(2)/pratb(2)
              ipcoefo=ipcoefe
              pdcoefe=pcoefe
              if(m.ne.0) pdcoefe=-pcoefe*rm*xln(nee)*xbn(nee)* &
                          xbn(nee)/((x*x+1.0e0_knd)*wmeta2(nee))
              ipdcoefe=ipcoefe
              pdcoefo=pdcoefe*pdratt(2)/pratb(2)
              ipdcoefo=ipdcoefe
              if(li.lt.3) go to 1190
                do jl=3,li+ix,2
                pcoefe=pcoefe*pratt(jl)/pratb(jl)
                iterm=log10(abs(pcoefe))
                pcoefe=pcoefe*ten**(-iterm)
                ipcoefe=ipcoefe+iterm
                pdcoefe=pdcoefe*pdratt(jl)/pratb(jl)
                iterm=log10(abs(pdcoefe))
                pdcoefe=pdcoefe*ten**(-iterm)
                ipdcoefe=ipdcoefe+iterm
                end do
              continue
              if(li.lt.4) go to 1190
                do jl=4,li+1-ix,2
                pcoefo=pcoefo*pratt(jl)/pratb(jl)
                iterm=log10(abs(pcoefo))
                pcoefo=pcoefo*ten**(-iterm)
                ipcoefo=ipcoefo+iterm
                pdcoefo=pdcoefo*pdratt(jl)/pratb(jl)
                iterm=log10(abs(pdcoefo))
                pdcoefo=pdcoefo*ten**(-iterm)
                ipdcoefo=ipdcoefo+iterm
                end do
1190          if(ix.eq.0) go to 1200
              pcoefet=pcoefo
              ipcoefet=ipcoefo
              pcoefo=pcoefo*pratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pcoefo)))
              pcoefo=pcoefo*ten**(-iterm)
              ipcoefo=ipcoefo+iterm
              pdcoefet=pdcoefo
              ipdcoefet=ipdcoefo
              pdcoefo=pdcoefo*pdratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pdcoefo)))
              pdcoefo=pdcoefo*ten**(-iterm)
              ipdcoefo=ipdcoefo+iterm
              go to 1210
1200          pcoefet=pcoefe
              ipcoefet=ipcoefe
              pcoefe=pcoefe*pratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pcoefe)))
              pcoefe=pcoefe*ten**(-iterm)
              ipcoefe=ipcoefe+iterm
              pdcoefet=pdcoefe
              ipdcoefet=ipdcoefe
              pdcoefe=pdcoefe*pdratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pdcoefe)))
              pdcoefe=pdcoefe*ten**(-iterm)
              ipdcoefe=ipdcoefe+iterm
1210          continue
                if(knd.eq.kindd) then
                if(x.gt.0.05e0_knd) limeta=2*int(25/(x*x)+300/x+3*c+ &
                                           1250*knd)
                if(x.ge.0.1e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+ &
                                          200)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+ &
                                          300)/x)
                if(x.ge.1.0e0_knd) limeta=2*int(l-m+c/5+0.5e0_knd*m+300)
                end if
                if(knd.eq.kindq) then
                if(x.gt.0.05e0_knd) limeta=2*int(25/(x*x)+400/x+3*c+ &
                                           1250*knd)
                if(x.ge.0.1e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+ &
                                          350)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+ &
                                          400)/x)
                if(x.ge.1.0e0_knd) limeta=2*int(l-m+c/5+0.5e0_knd*m+400)
                end if
              if(iopeta.eq.3.and.naccrsav.gt.minacc) &
                   limeta=jeta+jeta+40+ &
                          int(sqrt(c)*(2/min(1.0e0_knd,x))+5/x)
              if(iopeta.eq.3.and.naccrsav.le.minacc) &
                              limeta=jeta+jeta+500+c
              if(iopeta.eq.2) limeta=max(limeta,jeta+jeta+500+int(c))
              if(limeta.gt.limp-2) limeta=limp-2
              if(limeta.gt.limd) limeta=limd
              wm=wmeta2(nee)
              call r2eta(l,m,c,x,etaval,nee,limeta,ndec,maxd, &
                         maxlp,maxn,maxp,minacc,wm,enr,sneufe, &
                         sneune,ineuee,sneudfe,sneudre,pdratt,pratb, &
                         pratt,pcoefet,ipcoefet,pdcoefet,ipdcoefet, &
                         r1cin,ir1ein,r1dcin,ir1dein,naccnmax,r2ec, &
                         ir2ee,r2dec,ir2dee,nacceta,jeta,naccd)
              netatry=netatry+1
              naccetas=nacceta
                if(naccetas.ge.naccetamax) then
                neemax=nee
                naccetamax=naccetas
                jetam=jeta
                end if
              if(naccetamax.ge.naccd) iopnee=1
              naccrs=naccr
              if(naccetas.lt.naccrs) go to 1260
              naccr=naccetas
              r2c(li)=r2ec
              ir2e(li)=ir2ee
              r2dc(li)=r2dec
              ir2de(li)=ir2dee
1260          continue
if (debug) then
1270          if(knd.eq.kindd) write(40,1280) naccetas,etaval,nee,r2ec,ir2ee,r2dec,ir2dee
              if(knd.eq.kindq) write(40,1285) naccetas,etaval,nee,r2ec,ir2ee,r2dec,ir2dee
1280          format(15x,'r2eta accuracy = ',i2,' decimal digits; eta', &
                     ' = ',f12.9,'; nee = ',i4,/,10x,'r2 = ', f19.15,i5, &
                     5x,'r2d = ',f19.15,i5)
1285                  format(15x,'r2eta accuracy = ',i2,' decimal digits; eta', &
                     ' = ',f12.9,'; nee = ',i4,/,10x,'r2 = ', f34.30,i5, &
                     5x,'r2d = ',f34.30,i5)
end if
              iopeta=3
                if(naccetas.lt.naccetamax-2.or.nee.eq.neta) then
                nee=neemax-incnee
                iopeta=2
                iplflag=0
                go to 1310
                end if
                jetaflag=0
                if(naccetas.ge.minacc) then
                jetaflag=1
                iopint=0
                ietacount=ietacount+1
                if(ietacount.ge.5) incnflag=1
                  if(iplflag.eq.0.and.nee.gt.incnee+incnee) then
                  nee=nee-incnee
                  iopeta=2
                  end if
                iplflag=1
                go to 1330
                end if
              iopeta=2
              if(iplflag.eq.1.and.incnflag.eq.1.and.netatry.eq.2) &
                     iopnee=0
              ietacount=0
              if(iopnee.eq.0) go to 1290
              nee=neemax-incnee
              if(nee.eq.0) nee=incnee
                if(iopint.ne.0.and.naccint.ge.naccetamax) then
                iopeta=0
                iflageta=1
                end if
              go to 1310
1290          if(nee.eq.neta) go to 1300
              nee=nee+incnee
              go to 1090
1300          continue
              iopeta=3
              if(naccetas.lt.minacc.and.nee.eq.neta) &
                     nee=nee-incnee
              if(naccetas.lt.minacc) iopeta=2
              if(naccetas.lt.minacc) iplflag=0
1310          continue
              if(naccr.ge.minacc.or.(iopeta.ne.0.and.naccetas.gt. &
                   minacc)) go to 1320
              if(istartint.eq.0) istartint=1
              if(istartint.eq.1) iopint=1
              if(istartint.eq.1) go to 780
1320          continue
              if(istartint.eq.2) istartint=3
              if(iopeta.eq.4) iopeta=1
              if(iopint.ne.0.and.naccint.gt.nacintsa) nacintsa=naccint
                if(l.eq.m.and.iopint.ne.0.and.naccint.eq.naccr) then
                if(iopeta.ne.0) iflageta=1
                iopeta=0
                end if
              if(naccr.gt.0) go to 1330
              naccr=0
              r2c(li)=1.0e0_knd
              ir2e(li)=-ndec
              r2dc(li)=1.0e0_knd
              ir2de(li)=-ndec
1330          continue
                if(ioprad.eq.1) then
if (output) then
                write(20,1340) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),naccr1
1340                 format(1x,i5,2x,2(f17.14,i6,2x),i2)
end if
                go to 1400
                end if
                if(x.eq.0.0e0_knd) then
if (output) then
                write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li), &
                         r2c(li),ir2e(li),r2dc(li),ir2de(li),naccr,chr_e
end if
                go to 1400
                end if
                if(naccleg.lt.naccrt.and.naccleg1.lt.naccrt.and.naccint &
                   .lt.naccrt.and.naccneu0.lt.naccrt.and.naccneu1.lt. &
                   naccrt.and.nacceta.lt.naccrt) then
                if(ix.eq.0) naccflag=1
                naccrsav=naccr
if (debug) then
                if(ix.eq.0) write(40,1360) l,l+1
                if(ix.eq.1) write(40,1370) l,l-1
end if
                if(ix.eq.0) go to 1400
                  if(naccflag.eq.1) then
                  r2c(li-1)=r1c(li)
                  ir2e(li-1)=ir1e(li)
                  r2dc(li-1)=r1dc(li)
                  ir2de(li-1)=ir1de(li)
                  wronc=r1c(li-1)*r2dc(li-1)*ten**(ir1e(li-1)+ &
                      ir2de(li-1))-r2c(li-1)*r1dc(li-1)* &
                      ten**(ir2e(li-1)+ir1de(li-1))
                  naccrp=-int(log10(abs((wronc-wront)/wront)+dec))
                  if(naccrp.lt.0) naccrp=0
                  naccrp=min(naccrp,naccr1p,itestmp-1)
if (output) then
                  write(20,1380) l-1,r1c(li-1),ir1e(li-1),r1dc(li-1), &
                         ir1de(li-1),r2c(li-1),ir2e(li-1),r2dc(li-1), &
                         ir2de(li-1),naccrp,chr_w
end if
if (debug) then
                  write(40,1395) naccrp,l-1
end if
                  end if
                r2c(li)=-r1c(li-1)
                ir2e(li)=ir1e(li-1)
                r2dc(li)=-r1dc(li-1)
                ir2de(li)=ir1de(li-1)
                wronc=r1c(li)*r2dc(li)*ten**(ir1e(li)+ &
                      ir2de(li))-r2c(li)*r1dc(li)* &
                      ten**(ir2e(li)+ir1de(li))
                naccr=-int(log10(abs((wronc-wront)/wront)+dec))
                if(naccr.lt.0) naccr=0
                naccr=min(naccr,naccr1,itestm-1)
if (output) then
                write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li), &
                               r2c(li),ir2e(li),r2dc(li),ir2de(li), &
                               naccr,chr_w
end if
if (debug) then
                write(40,1395) naccr,l
end if
                naccflag=0
                go to 1400
                end if
                if(ix.eq.1.and.naccflag.eq.1) then
                naccflag=0
                r2c(li-1)=r1c(li)
                ir2e(li-1)=ir1e(li)
                r2dc(li-1)=r1dc(li)
                ir2de(li-1)=ir1de(li)
                wronc=r1c(li-1)*r2dc(li-1)*ten**(ir1e(li-1)+ &
                      ir2de(li-1))-r2c(li-1)*r1dc(li-1)* &
                      ten**(ir2e(li-1)+ir1de(li-1))
                naccrp=-int(log10(abs((wronc-wront)/wront)+dec))
                if(naccrp.lt.0) naccrp=0
                naccrp=min(naccrp,naccr1p,itestmp-1)
if (output) then
                write(20,1380) l-1,r1c(li-1),ir1e(li-1),r1dc(li-1), &
                         ir1de(li-1),r2c(li-1),ir2e(li-1),r2dc(li-1), &
                         ir2de(li-1),naccrp,chr_w
end if
                end if
1360          format(8x,'Values for r2 and r2d for ','l = ',i5, &
                     ' are given by r1 and r1d for l = ',i5)
1370          format(8x,'Values for r2 and r2d for ','l = ',i5, &
                     ' are given by -r1 and -r1d for l = ',i5)
              naccr=min(naccr,naccr1,itestm-1)
                if(jflagl.eq.1.and.naccr.eq.naccleg.and. &
                    naccr.ne.naccint) then
if (output) then
                write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li), &
                         r2c(li),ir2e(li),r2dc(li),ir2de(li),naccr,chr_e
end if
                else
if (output) then
                write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li), &
                         r2c(li),ir2e(li),r2dc(li),ir2de(li),naccr,chr_w
end if
                end if
1380          format(1x,i5,2x,4(f17.14,i6,2x),i2,a)
              go to 1400
1390          continue
                if(ix.eq.1) then
                r2c(li-1)=r1c(li)
                ir2e(li-1)=ir1e(li)
                r2dc(li-1)=r1dc(li)
                ir2de(li-1)=ir1de(li)
                wronc=r1c(li-1)*r2dc(li-1)*ten**(ir1e(li-1)+ &
                      ir2de(li-1))-r2c(li-1)*r1dc(li-1)* &
                      ten**(ir2e(li-1)+ir1de(li-1))
                naccrp=-int(log10(abs((wronc-wront)/wront)+dec))
                if(naccrp.lt.0) naccrp=0
                naccrp=min(naccrp,naccr1p,itestmp-1)
if (output) then
                write(20,1380) l-1,r1c(li-1),ir1e(li-1),r1dc(li-1), &
                       ir1de(li-1),r2c(li-1),ir2e(li-1),r2dc(li-1), &
                       ir2de(li-1),naccrp,chr_w
end if
                r2c(li)=-r1c(li-1)
                ir2e(li)=ir1e(li-1)
                r2dc(li)=-r1dc(li-1)
                ir2de(li)=ir1de(li-1)
                wronc=r1c(li)*r2dc(li)*ten**(ir1e(li)+ &
                      ir2de(li))-r2c(li)*r1dc(li)* &
                      ten**(ir2e(li)+ir1de(li))
                naccr=-int(log10(abs((wronc-wront)/wront)+dec))
                if(naccr.lt.0) naccr=0
                naccr=min(naccr,naccr1,itestm-1)
if (output) then
                write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li), &
                               r2c(li),ir2e(li),r2dc(li),ir2de(li), &
                               naccr,chr_w
end if
if (debug) then
                write(40,1395) naccrp,l-1
                write(40,1395) naccr,l
1395            format(10x,'Accuracy using eigenvalue match is ',i3, &
                       ' digits for l = ',i5)
end if
                end if
1400          continue
              ir1ep=ir1e(li)
              naccrsav=naccr
              naccr1p=naccr1
              itestmp=itestm
              naccrtp=naccrt
1405          continue
if (warn) then
                if(ioprad.eq.2.and.li.lt.matlilim) then
                if(ix.eq.1.and.naccrp.lt.8) write(60,*) &
                   ' est. acc. = ',naccrp, ' digits for x = ',x,' c = ', &
                   c,' m = ',m,' l = ',l-1
                if(ix.eq.1.and.naccr.lt.8) write(60,*) &
                   ' est. acc. = ',naccr, ' digits for x = ',x,' c = ', &
                   c,' m = ',m,' l = ',l
                end if
              if(ioprad.eq.2.and.li.ge.matlilim.and.naccr.lt.8) &
                  write(60,*) ' est. acc. = ',naccr, ' digits for x = ', &
                  x,' c = ', c,' m = ',m,' l = ',l
              if(ioprad.eq.1.and.naccr1.lt.8) write(60,*) &
                 'est. r1 acc. = ',naccr1,' digits for x = ',x,' c = ', &
                 c,' m = ',m,' l = ',l
end if
              naccrp=naccr
              if(iopang.eq.0) go to 1510
!
!  determine first kind oblate angular function
1410          if(l.eq.m) lims1=3*ndec+int(c)
              if(l.ne.m) lims1=jang+jang+20+c/25
              if(lims1.gt.maxp) lims1=maxp
              call s1leg(l,m,c,iopang,iopnorm,barg,narg,lims1,ndec,nex, &
                         maxt,maxd,maxp,enr,dmlms,idmlmse,pr,pdr, &
                         pdnorm,ipdnorm,pnorm,ipnorm,pdtempe,ipdtempe, &
                         pdtempo,ipdtempo,ptempe,iptempe,ptempo,iptempo, &
                         itestm,kindd,kindq,s1c,is1e,s1dc,is1de,naccs, &
                         naccds,jang,dmlms1,idmlms1e)
                do 1500 jarg=1,narg

if (debug) then
                if(ioparg.eq.0.and.iopang.eq.1) write(50,1420) arg(jarg),naccs(jarg)
                if(ioparg.eq.0.and.iopang.eq.2) write(50,1425) arg(jarg),naccs(jarg),naccds(jarg)
1420            format(1x,'theta = ',e24.15,'   accuracy = ',i2, ' digits.')
1425            format(1x,'theta = ',e24.15,'   s1 and s1d accuracy = ', i2,' and ',i2,' digits.')

                if(ioparg.eq.1.and.iopang.eq.1) write(50,1430) barg(jarg),naccs(jarg)
                if(ioparg.eq.1.and.iopang.eq.2) write(50,1435) barg(jarg),naccs(jarg),naccds(jarg)
1430            format(1x,'eta = ',e24.15,'   accuracy = ',i2, ' digits.')
1435            format(1x,'eta = ',e24.15,'   s1 and s1d accuracy = ', i2,' and ',i2,' digits.')
end if
if (output) then
                if(ioparg.eq.0.and.iopang.eq.1) write(30,1440) arg(jarg),s1c(jarg),is1e(jarg),naccs(jarg)
                if(ioparg.eq.0.and.iopang.eq.2) write(30,1450) arg(jarg),s1c(jarg),is1e(jarg),s1dc(jarg), &
                                                               is1de(jarg),naccs(jarg),naccds(jarg)
1440            format(1x,f17.14,2x,f17.14,2x,i5,2x,', ',i2)
1450            format(1x,f17.14,2x,f17.14,2x,i5,2x,f17.14,2x,i5,2x,i2,', ',i2)
                if(ioparg.eq.1.and.iopang.eq.1) write(30,1460) barg(jarg),s1c(jarg),is1e(jarg),naccs(jarg)
                if(ioparg.eq.1.and.iopang.eq.2) write(30,1470) barg(jarg),s1c(jarg),is1e(jarg),s1dc(jarg), &
                                                               is1de(jarg),naccs(jarg),naccds(jarg)
1460            format(1x,f17.14,2x,f17.14,2x,i5,2x,', ',i2)
1470            format(1x,f17.14,2x,f17.14,2x,i5,2x,f17.14,2x,i5,2x,i2,', ',i2)
end if
if (debug) then
                if(knd.eq.kindd.and.iopang.eq.1) write(50,1480) s1c(jarg),is1e(jarg)
                if(knd.eq.kindd.and.iopang.eq.2) write(50,1490) s1c(jarg),is1e(jarg),s1dc(jarg),is1de(jarg)
                if(knd.eq.kindq.and.iopang.eq.1) write(50,1485) s1c(jarg),is1e(jarg)
                if(knd.eq.kindq.and.iopang.eq.2) write(50,1495) s1c(jarg),is1e(jarg),s1dc(jarg),is1de(jarg)
end if
1480            format(12x,'s1 = ',f17.14,2x,i5)
1485            format(12x,'s1 = ',f34.30,2x,i5)
1490            format(12x,'s1 = ',f17.14,2x,i5,5x,'s1d = ',f17.14, &
                       2x,i5)
1495            format(12x,'s1 = ',f34.30,2x,i5,5x,'s1d = ',f34.30, &
                       2x,i5)
1500            continue
1510          continue
1540        continue
        return
        end subroutine
!
!
        subroutine s1leg (l,m,c,iopang,iopnorm,barg,narg,lims1,ndec,nex, &
                          maxt,maxd,maxp,enr,dmlms,idmlmse,pr,pdr, &
                          pdnorm,ipdnorm,pnorm,ipnorm,pdtempe,ipdtempe, &
                          pdtempo,ipdtempo,ptempe,iptempe,ptempo, &
                          iptempo,itestm,kindd,kindq,s1c,is1e,s1dc, &
                          is1de,naccs,naccds,jang,dmlms1,idmlms1e)
!
!  purpose:     To calculate the oblate angular functions of the first
!               kind and their first derivatives with respect to eta.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               iopang  : index = 1 when angular functions of the
!                         first kind are calculated; = 2 when the
!                         first derivatives with respect to eta are
!                         also calculated
!               iopnorm : = 1 when the angular functions (and
!                         first derivatives) are scaled by the
!                         square root of the normalization of the
!                         corresponding Legendre function, giving them
!                         unity norm; iopnorm = 0 otherwise
!               barg    : array of eta values for which angular
!                         functions are desired
!               narg    : number of eta values
!               lims1   : approximately twice the maximum number
!                         of terms available to be taken in the
!                         sums
!               ndec    : number of decimal digits available in
!                         real arithmetic
!               nex     : maximum exponent in real arithmetic
!               maxt    : dimension of barg, pdnorm, ipdnorm, pnorm,
!                         ipnorm, pdtempe, ipdtempe, pdtempo, ipdtempo,
!                         ptempe, iptempe, ptempo, iptempo, s1c, is1e,
!                         s1dc, is1de, and naccs arrays. first dimension
!                         of the doubly dimensioned arrays pr and pdr
!               maxd    : dimension of enr array
!               maxp    : second dimension of pr and pdr arrays
!               enr     : array of d coefficient ratios
!               dmlms   : characteristic of the value of the d
!                         coefficient with index l - m using the
!                         Meixner-Schafke normalization
!               idmlmse : exponent associated with dmlms
!               pr      : array of ratios of successive first kind
!                         associated Legendre functions of the same
!                         parity
!               pdr     : array of ratios of successive derivatives of
!                         first kind associated Legendre functions of
!                         the same parity
!               pdnorm  : array of characteristics of the first
!                         derivatives of associated Legendre functions
!                         of the first kind of order m and degree m
!               ipdnorm : array of exponents corresponding to pdnorm
!               pnorm   : array of characteristics of the associated
!                         Legendre functions of the first kind of order
!                         m and degree m
!               ipnorm  : array of exponents corresponding to pnorm
!               pdtempe : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is even
!                         or odd, to the first derivative of the
!                         function of order m and degree m
!               ipdtempe: array of exponents corresponding to pdtempe
!               pdtempo : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is odd
!                         or even, to the first derivtive of the
!                         function of order m and degree m
!               ipdtempo: array of exponents corresponding to pdtempo
!               ptempe  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is even or odd, to the function
!                         of order m and degree m
!               iptempe : array of exponents corresponding to ptempe
!               ptempo  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is odd or even, to the
!                         function of order m and degree m
!               iptempo : array of exponents corresponding to ptempo
!               itestm  : number of matching digits for the forward
!                         and backward recursion for d coefficient
!                         ratios
!               kindd   : number of bytes for real data in double
!                         precision
!               kindq   : number of bytes for real data in quadruple
!                         precision
!
!
!     output:   s1c     : array of characteristics of oblate
!                         angular functions of the first kind
!               is1e    : array of exponents of oblate angular
!                         functions of the first kind
!               s1dc    : array of characteristics of derivative with
!                         respect to eta of oblate angular functions
!                         of the first kind
!               is1de   : array of exponents of derivatives with respect
!                         to eta of oblate angular functions of first
!                         kind
!               naccs   : array of integer estimates of the number of
!                         accurate decimal digits in the values obtained
!                         for s1
!               naccds  : array of integer estimates of the number of
!                         accurate decimal digits in the values obtained
!                         for s1d
!               jang    : maximum value of the index j in the forward
!                         sum for r1 and r1d, i.e., the highest enr(j)
!                         used
!               dmlms1  : characteristic of the d coefficient with index
!                         l - m when the angular functions have unity
!                         norm
!               idmlms1e: exponent associated with dmlms1
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) adec,aj,c,dcon,dec,dmlms,dmlms1,dnew,dnewd,dold, &
                  doldd,factor,fterm,rm2,rm2m1,rm2m3,rm2p1,s1,s1d,ten, &
                  teste,testeo
        real(knd) barg(maxt),enr(maxd),pdr(maxt,maxp),pdnorm(maxt), &
                  pnorm(maxt),pr(maxt,maxp),pdtemp(maxt),ptemp(maxt), &
                  pdtempe(maxt),ptempe(maxt),pdtempo(maxt),ptempo(maxt), &
                  s1c(maxt),s1dc(maxt)
!
!  integer arrays
        integer ipdnorm(maxt),ipnorm(maxt),ipdtemp(maxt),iptemp(maxt), &
                ipdtempe(maxt),iptempe(maxt),ipdtempo(maxt), &
                iptempo(maxt),is1de(maxt),is1e(maxt),naccs(maxt), &
                naccds(maxt)
!
        ten=10.0e0_knd
        dec=ten**(-ndec-1)
        dcon=dec
        adec=1000.0e0_knd*dec
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        jflag=0
        kflag=0
        rm2=m+m
        rm2m1=m+m-1
        rm2p1=m+m+1
        rm2m3=m+m-3
        if(l.gt.(m+1)) go to 30
          do 20 k=1,narg
          if(pnorm(k).eq.0.0e0_knd) go to 20
          if(l.eq.(m+1)) go to 10
          ptempe(k)=pr(k,1)
          iptempe(k)=0
          pdtempe(k)=pdr(k,1)
          ipdtempe(k)=0
          ptemp(k)=ptempe(k)
          pdtemp(k)=pdtempe(k)
          iptemp(k)=0
          ipdtemp(k)=0
          go to 20
10        ptempo(k)=pr(k,2)
          iptempo(k)=0
          pdtempo(k)=pdr(k,2)
          ipdtempo(k)=0
          ptemp(k)=ptempo(k)
          pdtemp(k)=pdtempo(k)
          iptemp(k)=0
          ipdtemp(k)=0
20        continue
30      continue
        lm2=(l-m)/2
        ix=l-m-2*lm2
        ixx=ix-1
        ixx2=ixx+2
        if(l.lt.(m+2)) go to 110
          do 100 k=1,narg
          if(pnorm(k).eq.0.0e0_knd) go to 100
          if(ix.ne.0) go to 60
          ptempe(k)=ptempe(k)*pr(k,l-m+1)
            if(abs(ptempe(k)).gt.1.0e+10_knd) then
            ptempe(k)=ptempe(k)*(1.0e-10_knd)
            iptempe(k)=iptempe(k)+10
            end if
          ptemp(k)=ptempe(k)
          iptemp(k)=iptempe(k)
          if(abs(barg(k)).lt.adec) go to 100
          pdtempe(k)=pdtempe(k)*pdr(k,l-m+1)
            if(abs(pdtempe(k)).gt.1.0e+10_knd) then
            pdtempe(k)=pdtempe(k)*(1.0e-10_knd)
            ipdtempe(k)=ipdtempe(k)+10
            end if
          pdtemp(k)=pdtempe(k)
          ipdtemp(k)=ipdtempe(k)
          go to 100
60        if(abs(barg(k)).lt.adec) go to 80
          ptempo(k)=ptempo(k)*pr(k,l-m+1)
            if(abs(ptempo(k)).gt.1.0e+10_knd) then
            ptempo(k)=ptempo(k)*(1.0e-10_knd)
            iptempo(k)=iptempo(k)+10
            end if
          ptemp(k)=ptempo(k)
          iptemp(k)=iptempo(k)
80        pdtempo(k)=pdtempo(k)*pdr(k,l-m+1)
            if(abs(pdtempo(k)).gt.1.0e+10_knd) then
            pdtempo(k)=pdtempo(k)*(1.0e-10_knd)
            ipdtempo(k)=ipdtempo(k)+10
            end if
          pdtemp(k)=pdtempo(k)
          ipdtemp(k)=ipdtempo(k)
100        continue
110     continue
        lim=lims1/2-ix
        jlow=lm2+1
        jang=0
!
!  compute the associated Legendre function normalization factor
        factor=1.0e0_knd
        ifactor=0
        if(iopnorm.eq.0) go to 210
        if(m.eq.0) go to 170
          do 160 j=1,m
          aj=j
          factor=factor*(aj+aj)*(aj+aj-1.0e0_knd)
            if(factor.gt.teste) then
            factor=factor*testeo
            ifactor=ifactor+nfac
            end if
160       continue
170     if(l.eq.m) go to 190
          do 180 j=1,l-m
          aj=j
          factor=factor*(rm2+aj)/(aj)
            if(factor.gt.teste) then
            factor=factor*testeo
            ifactor=ifactor+nfac
            end if
180       continue
190     continue
        factor=factor*2.0e0_knd/(l+l+1.0e0_knd)
          if(2*(ifactor/2).ne.ifactor) then
          factor=factor*ten
          ifactor=ifactor-1
          end if
        factor=sqrt(factor)
        ifactor=ifactor/2
        iterm=int(log10(factor))
        factor=factor*(ten**(-iterm))
        ifactor=ifactor+iterm
        dmlms1=dmlms/factor
        idmlms1e=idmlmse-ifactor
if (debug) then
        if(knd.eq.kindd) write(50,200) factor,ifactor
        if(knd.eq.kindq) write(50,205) factor,ifactor
200     format(1x,'square root of Legendre norm = ',f19.15,2x,i5)
205     format(1x,'square root of Legendre norm = ',f34.30,2x,i5)
end if
210     continue
!
!  compute the angular function s1
          do 380 k=1,narg
          if(pnorm(k).eq.0.0e0_knd) go to 220
          if((ix.eq.1).and.(abs(barg(k)).lt.adec)) go to 220
          if(((abs(abs(barg(k))-1.0e0_knd)).lt.adec) &
               .and.(m.ne.0)) go to 220
          go to 230
220       s1c(k)=0.0e0_knd
          is1e(k)=0
          naccs(k)=ndec
          go to 300
230       dold=1.0e0_knd
          s1=dold
          is1=0
          fterm=s1
          lflag=0
            do 240 j=jlow,lim
            dnew=dold*enr(j)*pr(k,j+j+ixx2)
            s1=s1+dnew
            if(abs(s1).gt.fterm) fterm=abs(s1)
            if(abs(dnew/s1).lt.dcon) go to 250
              if(abs(s1).gt.teste) then
              s1=s1*testeo
              dnew=dnew*testeo
              fterm=fterm*testeo
              is1=is1+nfac
              kflag=1
              end if
            dold=dnew
240         continue
250       if(j.gt.jang) jang=j
if (debug) then
          write(50,260) barg(k),j
260       format(8x,'s1 calculation for eta = ',f13.8,' converged in ', &
                 i6,' terms.')
end if
          if(lm2.lt.1.or.kflag.eq.1) go to 280
          dold=1.0e0_knd
          j=lm2
            do 270 jj=1,lm2
            dnew=dold/(pr(k,j+j+ixx2)*enr(j))
            s1=s1+dnew
            if(abs(s1).gt.fterm) fterm=abs(s1)
            if(abs(dnew/s1).lt.dcon) go to 280
            dold=dnew
            j=j-1
270         continue
280       s1c(k)=s1*dmlms*ptemp(k)*pnorm(k)/factor
          if(s1c(k).ne.0.0e0_knd) iterm=int(log10(abs(s1c(k))))
          if(s1c(k).eq.0.0e0_knd) iterm=0
          s1c(k)=s1c(k)*(ten**(-iterm))
          is1e(k)=is1+iptemp(k)+ipnorm(k)+iterm-ifactor+idmlmse
          if(abs(s1c(k)).ge.1.0e0_knd) go to 290
          s1c(k)=s1c(k)*ten
          is1e(k)=is1e(k)-1
290       continue
          if(s1.eq.0.0e0_knd) naccs(k)=0
            if(s1.ne.0.0e0_knd) then
            iacc=int(log10(abs(fterm/s1)+dec))
            if(iacc.lt.0) iacc=0
            if(iacc.gt.ndec) iacc=ndec
            naccs(k)=min(ndec-2,itestm-1)-iacc
            if(naccs(k).lt.0) naccs(k)=0
            end if
          if(naccs(k).gt.0) go to 300
          naccs(k)=0
          naccds(k)=0
          s1c(k)=0.0e0_knd
          is1e(k)=0
          s1dc(k)=0.0e0_knd
          is1de(k)=0
          go to 380
!
!       compute the first derivative of the anguar function when
!       iopang equals 2
300       if(iopang.ne.2) go to 380
          if(pnorm(k).eq.0.0e0_knd) go to 310
          if((ix.eq.0).and.(abs(barg(k)).lt.adec)) go to 310
          if(((abs(abs(barg(k))-1.0e0_knd)).lt.adec).and.(m.ne.0) &
              .and.(m.ne.2)) go to 310
          go to 320
310       s1dc(k)=0.0e0_knd
          is1de(k)=0
          naccds(k)=ndec
          go to 370
320       doldd=1.0e0_knd
          s1d=doldd
          is1d=0
          if(l.eq.0) s1d=0.0e0_knd
          fterm=s1d
            do 330 j=jlow,lim
            dnewd=doldd*enr(j)*pdr(k,j+j+ixx2)
            s1d=s1d+dnewd
            if(abs(s1d).gt.fterm) fterm=abs(s1d)
            if(abs(dnewd/s1d).lt.dcon) go to 340
              if(abs(s1d).gt.teste) then
              s1d=s1d*testeo
              dnewd=dnewd*testeo
              fterm=fterm*testeo
              fterm=fterm*testeo
              is1d=is1d+nfac
              end if
            doldd=dnewd
330         continue
340       if(lm2.lt.1.or.kflag.eq.1) go to 360
          doldd=1.0e0_knd
          j=lm2
          ja=lm2
          if(m.eq.0.and.ix.eq.0) ja=lm2-1
          if(ja.eq.0) go to 360
            do 350 jj=1,ja
            dnewd=doldd/(pdr(k,j+j+ixx2)*enr(j))
            s1d=s1d+dnewd
            if(abs(s1d).gt.fterm) fterm=abs(s1d)
            if(abs(dnewd/s1d).lt.dcon) go to 360
            doldd=dnewd
            j=j-1
350         continue
360       s1dc(k)=s1d*dmlms*pdtemp(k)*pdnorm(k)/factor
          if(s1dc(k).ne.0.0e0_knd) iterm=int(log10(abs(s1dc(k))))
          if(s1dc(k).eq.0.0e0_knd) iterm=0
          s1dc(k)=s1dc(k)*ten**(-iterm)
          is1de(k)=is1d+ipdtemp(k)+ipdnorm(k)+iterm-ifactor+idmlmse
          if(s1d.eq.0.0e0_knd) naccds(k)=0
            if(s1d.ne.0.0e0_knd) then
            iacc=int(log10(abs(fterm/s1d)+dec))
            if(iacc.lt.0) iacc=0
            if(iacc.gt.ndec) iacc=ndec
            naccds(k)=min(ndec-2,itestm-1)-iacc
            if(naccds(k).lt.0) naccds(k)=0
            end if
          if(abs(s1dc(k)).ge.1.0e0_knd) go to 370
          s1dc(k)=s1dc(k)*ten
          is1de(k)=is1de(k)-1
370       continue
          if(naccds(k).eq.0) s1dc(k)=0.0e0_knd
          if(naccds(k).eq.0) is1de(k)=0
380       continue
        return
        end subroutine
!
!
        subroutine r1bes(l,m,c,x,limr1,ndec,maxd,enr,maxj,maxn,maxlp, &
                         nex,iflag,sbesf,sbesdf,sbesn,ibese,sbesdr, &
                         prat1,pcoefn,ipcoefn,dmfnorm,idmfe,ir1ep,r1c, &
                         ir1e,r1dc,ir1de,jbes,nsub,ndsub)
!
!  purpose:     To calculate the oblate radial function of the
!               first kind and its first derivative with respect
!               to x, using the traditional expansion of spherical
!               Bessel functions of the first kind with argument
!               c*x, i.e., with eta = 1.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x      : x
!               limr1  : approximately twice the maximum number of
!                        terms available to be taken in the series
!               ndec   : number of decimal digits available in
!                        real arithmetic
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               maxj   : dimension of sbesf and sbesdf arrays
!               maxn   : dimension of prat1 array
!               maxlp  : dimension of the sbesdr, sbesn, and ibese
!                        arrays
!               nex    : maximum exponent available in real(knd)
!                        arithmetic
!               iflag  : integer = 1 if forward series not needed;
!                        =0 if the forward series is computed
!               sbesf  : array of ratios of successive first kind
!                        spherical Bessel functions of the same parity
!               sbesdf : array of ratios of successive derivatives of
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : array of characteristics for Bessel functions
!               ibese  : array of exponents corresponding to sbesn
!               sbesdr : value of ratio of first derivative of
!                        spherical Bessel function to the corresponding
!                        Bessel function
!               prat1  : array of ratios of successive coefficients in
!                        r1 and r1d sum
!               pcoefn : characteristic of coefficient for term in both
!                        r1 and r1d sums that contains Bessel function
!                        of order l
!               ipcoefn: exponent (to the base 10) corresponding to
!                        pcoefn
!               dmfnorm : characteristic of Morse-Feshbach normalization
!                         sum of the d coefficients. equal to the
!                         reciprocal of the value of the d coefficient
!                         d(n = l - m) using this normalization for the
!                         angular functions
!               idmfe   : exponent associated with dmfnorm
!               ir1ep  : exponent for the value of r1 for l-1
!
!     output  : r1c    : characteristic of oblate radial function
!                        of the first kind
!               ir1e   : exponent of oblate radial function of the
!                        first kind
!               r1dc   : characteristic of derivative with respect
!                        to x of oblate radial function of the first
!                        kind
!               ir1de  : exponent of derivative with respect to x of
!                        oblate radial function of the first kind
!               jbes   : maximum value of the index j in the forward
!                        sum for r1 and r1d, i.e., the highest enr(j)
!                        used
!               nsub   : subtraction error in calculating r1
!               ndsub  : subtraction error in calculating r1d
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,coef,dec,dmfnorm,dnew,dnewd,dold,doldd,em,pcoefn, &
                  r1c,r1ca,r1dc,r1dcoef,r1temp,r1tempa,r1dtemp, &
                  r1dtempa,r1top,r1topa,r1dtop,spos,sposa,sdpos,ten, &
                  term,teste,testeo,x,x2
        real(knd) enr(maxd),prat1(maxn),sbesdf(maxj),sbesdr(maxlp), &
                  sbesf(maxj),sbesn(maxlp)
!
!  integer array
        integer ibese(maxlp)
!
!  convergence ratio dec is set according to the requested accuracy
        ten=10.0e0_knd
        dec=ten**(-ndec-1)
        lm2=(l-m)/2
!
!  ix=0 for l-m even, ix=1 for l-m odd
        ix=l-m-2*lm2
        lim=limr1/2-ix
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        ir1tope=0
        mml=m+m-1+ix
        em=m
        jflag=0
        kflag=0
        nsuba=0
        term=0.0e0_knd
        x2=x*x
        r1dcoef=-em/(x*(x*x+1.0e0_knd))
        coef=-(sbesn(m+2)/(sbesn(m+1)*sbesdr(m+1)))* &
             ten**(ibese(m+2)-ibese(m+1))
        r1topa=0.0e0_knd
        sposa=0.0e0_knd
!
!  forward summation of numerator series for both r1 and r1d
        dold=1.0e0_knd
        doldd=1.0e0_knd
        spos=1.0e0_knd
        sdpos=1.0e0_knd
        r1top=dold
        r1dtop=doldd
          if(lm2.eq.0.and.ix.eq.0) then
          jflag=1
          r1topa=-x2
          sposa=0.0e0_knd
          r1dtop=coef
          sdpos=0.0e0_knd
          if(coef.gt.0.0e0_knd) sdpos=coef
          end if
          if(iflag.eq.1) then
if (debug) then
          write(40,10)
10        format(8x,'r1bes: forward series not used.')
end if
          jtop=lm2
          go to 50
          end if
          do 20 j=lm2+1,lim
          jj=j+j+ix
          dnew=-dold*enr(j)*sbesf(jj+m)*prat1(jj+1)
          dnewd=-doldd*enr(j)*sbesdf(jj+m)*prat1(jj+1)
          r1top=r1top+dnew
          r1dtop=r1dtop+dnewd
            if(jflag.eq.1) then
            r1topa=r1topa+dnew
            if(dnew.gt.0.0e0_knd) sposa=sposa+dnew
            end if
          if(dnew.gt.0.0e0_knd) spos=spos+dnew
          if(dnewd.gt.0.0e0_knd) sdpos=sdpos+dnewd
          if((abs(dnew/r1top)+abs(dnewd/r1dtop)).lt.dec) go to 30
            if(abs(r1top).gt.teste) then
            r1top=r1top*testeo
            dnew=dnew*testeo
            spos=spos*testeo
            r1dtop=r1dtop*testeo
            dnewd=dnewd*testeo
            sdpos=sdpos*testeo
            ir1tope=ir1tope+nfac
            kflag=1
              if(jflag.eq.1) then
              r1topa=r1topa*testeo
              sposa=sposa*testeo
              end if
            end if
          dold=dnew
          doldd=dnewd
20        continue
30      continue
        jtop=min(j,lim)
        jterms=jtop-lm2
if (debug) then
        write(40,40) lim,jterms
40      format(8x,'r1bes: ',i6,' total terms ', &
               'available; forward series converged in ',i6,' terms.')
end if
50    continue
!
!  backward summation of numerator series for r1 and r1d
        if(lm2.lt.1.or.kflag.eq.1) go to 80
        dold=1.0e0_knd
        doldd=1.0e0_knd
          do 70 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sbesf(jj+m)*prat1(jj+1)*enr(j))
          dnewd=-doldd/(sbesdf(jj+m)*prat1(jj+1)*enr(j))
            if(j.eq.1.and.ix.eq.0) then
            jflag=1
            term=-x2*dnew
            r1topa=r1top+term
            dnewd=coef*dnewd
            if(term.gt.0.0e0_knd) sposa=sposa+term
            end if
          r1top=r1top+dnew
          r1dtop=r1dtop+dnewd
          if(dnew.gt.0.0e0_knd) spos=spos+dnew
          if(dnewd.gt.0.0e0_knd) sdpos=sdpos+dnewd
          if(j.eq.1) go to 80
          if((abs(dnew/r1top)+abs(dnewd/r1dtop)).lt.dec) go to 80
            if(abs(r1top).gt.teste) then
            r1top=r1top*testeo
            dnew=dnew*testeo
            ir1tope=ir1tope+nfac
            r1dtop=r1dtop*testeo
            dnewd=dnewd*testeo
            spos=spos*testeo
            sdpos=sdpos*testeo
            if(spos.lt.abs(r1top)*dec) spos=abs(r1top)*dec
            if(sdpos.lt.abs(r1dtop)*dec) sdpos=abs(r1dtop)*dec
            iflag=1
            end if
60        dold=dnew
          doldd=dnewd
70        continue
80        continue
        nsub=0
        if(spos/r1top.ne.0.0e0_knd) nsub=int(log10(abs(spos/r1top)))
        if(nsub.lt.0) nsub=0
        if(nsub.gt.ndec) nsub=ndec
        nsuba=0
          if(jflag.eq.1) then
          if(term.gt.0.0e0_knd) sposa=sposa+term
          if(sposa/r1topa.ne.0.0e0_knd) nsuba= &
                                        int(log10(abs(sposa/r1topa)))
          if(nsuba.lt.0) nsuba=0
          if(nsuba.gt.ndec) nsuba=ndec
          end if
        ndsub=0
        if(sdpos/r1dtop.ne.0.0e0_knd) ndsub= &
                                      int(log10(abs(sdpos/r1dtop)))
        if(ndsub.lt.0) ndsub=0
        if(ndsub.gt.ndec) ndsub=ndec
!
!  compute r1 and r1d
        r1temp=r1top*sbesn(l+1)*pcoefn/dmfnorm
        iterm=0
        if(r1temp.ne.0.0e0_knd) iterm=int(log10(abs(r1temp)))
        ir1e=ir1tope+ibese(l+1)+ipcoefn-idmfe+iterm
        r1c=r1temp*(ten**(-iterm))
        if(abs(r1c).ge.1.0e0_knd) go to 90
        r1c=r1c*ten
        ir1e=ir1e-1
90      continue
          if(jflag.eq.0) then
          r1ca=r1c
          ir1ea=ir1e
          else
          r1tempa=r1topa*sbesn(l+1)*pcoefn/dmfnorm
          iterm=0
          if(r1tempa.ne.0.0e0_knd) iterm=int(log10(abs(r1tempa)))
          ir1ea=ir1tope+ibese(l+1)+ipcoefn-idmfe+iterm
          r1ca=r1tempa*(ten**(-iterm))
          end if
        continue
        r1dtemp=r1dcoef*r1ca
        r1dtempa=(c*r1dtop*sbesn(l+1)*sbesdr(l+1)*pcoefn/ &
                  dmfnorm)*ten**(ibese(l+1)+ipcoefn+ir1tope-idmfe- &
                  ir1ea)
        r1dc=r1dtemp+r1dtempa
        ndsub1=0
        if(r1dtemp.ne.0.0e0_knd.and.r1dc.ne.0.0e0_knd) ndsub1= &
                                     int(log10(abs(r1dtemp/r1dc)))
        if(ndsub1.lt.0) ndsub1=0
        n1=0
        if(r1dtemp.ne.0.0e0_knd.and.r1dtempa.ne.0.0e0_knd) n1= &
                                  int(log10(abs(r1dtemp/r1dtempa)))
        if(n1.gt.0.and.jflag.eq.1) ndsub=max(nsuba,ndsub-n1)+ndsub1
        if(n1.le.0.and.jflag.eq.1) ndsub=max(ndsub,nsuba+n1)+ndsub1
        if(n1.gt.0.and.jflag.eq.0) ndsub=max(nsub,ndsub-n1)+ndsub1
        if(n1.le.0.and.jflag.eq.0) ndsub=max(ndsub,nsub+n1)+ndsub1
        if(ndsub.gt.ndec) ndsub=ndec
100     iterm=0
        if(r1dc.ne.0.0e0_knd) iterm=int(log10(abs(r1dc)))
        ir1de=ir1ea+iterm
        r1dc=r1dc*ten**(-iterm)
        if(abs(r1dc).ge.1.0e0_knd) go to 110
        r1dc=r1dc*ten
        ir1de=ir1de-1
110    continue
if (debug) then
        if(nsub+ndsub.gt.0) write(40,120) nsub,ndsub
120     format(15x,'subtraction errors in r1 and r1d are ',i2,' and ', &
               i2,' digits.')
end if
        jbes=jtop
130     return
        end subroutine
!
!
        subroutine r1eta (l,m,c,x,eta,nee,limeta,ndec,nex,maxd,maxlp, &
                          maxj,maxp,minacc,wm,enr,sbesf,sbesn,ibese, &
                          sbesdf,sbesdr,pdratt,pratb,pratt,pcoefn, &
                          ipcoefn,pdcoefn,ipdcoefn,ir1ep,r1c,ir1e,r1dc, &
                          ir1de,naccs1,naccs2,jeta)
!
!  purpose:     To calculate the oblate radial function of the
!               first kind and its first derivative with respect
!               to x, using an expansion of spherical Bessel
!               functions.
!
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : c
!               x       : x
!               eta     : value for eta used in calculation
!               nee     : index in the array of eta values in the main
!                         program that corresponds to the value of eta
!                         used in r1eta calculations
!               limeta  : maximum number of terms available in the sums
!                         for r1 and r1d
!               ndec    : number of decimal digits available in
!                         real arithmetic
!               nex     : maximum exponent available in real(knd)
!                         arithmetic
!               maxd    : dimension of enr array
!               maxlp   : maximum  l value desired; dimension
!                         of the sbesn, sbesdr, and ibese arrays
!               maxj    : dimension of sbesf and sbesdf arrays
!               maxp    : dimension of pdratt, pratb, and pratt arrays
!               minacc  : minimum number of accurate decimal digits
!                         that are requested
!               wm      : value of 1 - eta*eta computed in a way that
!                         avoids the subtraction error that would occur
!                         if it were computed directly when eta is near
!                         unity
!                         subtraction error that would occur
!               enr     : array of ratios of successive d coefficients
!               sbesf   : array of ratios of successive spherical
!                         Bessel functions of the same parity
!               sbesn   : array of characteristics for Bessel functions
!               ibese   : array of exponents corresponding to sbesn
!               sbesdf  : array of ratios of successive first
!                         derivatives of spherical Bessel functions of
!                         the same parity
!               sbesdr  : array of ratios of first derivatives of the
!                         spherical Bessel functions to the
!                         corresponding functions
!               pdratt  : array of ratios of successive first
!                         derivatives of the associated Legendre
!                         functions of the first kind of the same parity
!                         (used in numerator series)
!               pratb   : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in denominator series)
!               pratt   : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in numerator series)
!               pcoefn  : characteristic of the ratio of the numerator
!                         and denominator associated Legendre functions
!                         of the first kind of order m and degree l
!               ipcoefn : exponent corresponding to pcoefn
!               pdcoefn : characteristic of the ratio of the first
!                         derivative of the associated Legendre function
!                         of the first kind in the numerator and the
!                         associated Legendre function of the first kind
!                         in the denominator, both of order m and
!                         degree l
!               ipdcoefn: exponent corresponding to pdcoefn
!               ir1ep   : exponent for the value of r1 for l-1
!
!     output:   r1c     : characteristic of the radial function of the
!                         first kind
!               irie    : exponent corresponding to r1c
!               r1dc    : characteristic of the first derivative with
!                         respect to x of the radial function of the
!                         first kind
!               ir1de   : exponent corresponding to r1dc
!               naccs1  : larger of the subtraction error for the
!                         numerator series for r1 and the subtraction
!                         error for the denominator series
!               naccs2  : larger of the subtraction error for the
!                         numerator series for r1d and the subtraction
!                         error for the denominator
!                         series
!               jeta    : maximum number of terms taken in the numerator
!                         and denominator sums for r1 and r1d
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,dec,denom,dnew,dnewd,dnewd1,dnewd2,dold,doldd1, &
                  doldd2,eta,etas,pcoefn,pdcoefn,reld12,rm,rm2,r1c, &
                  r1dc,r1dcoef1,r1dcoef2,r1dtemp,r1temp,sumdnp, &
                  sumdp,sumnp,ten,test,testd,teste,testeo,wm,xet,xets,x
        real(knd) enr(maxd),sbesdr(maxlp),sbesn(maxlp),pratb(maxp), &
                  pratt(maxp),pdratt(maxp),sbesf(maxj),sbesdf(maxj)
!
!  integer arrays
        integer ibese(maxlp)
!
        ten=10.0e0_knd
        dec=ten**(-ndec-2)
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        ir1tempe=0
        iflag=1
        rm=m
        etas=eta*eta
        xet=sqrt(x*x+wm)
        xets=xet*xet
        r1dcoef1=eta*wm/(xets*xet)
        r1dcoef2=c*x/xet
        reld12=(r1dcoef2/r1dcoef1)*sbesdr(l+1)*(pcoefn/ &
               pdcoefn)*ten**(ipcoefn-ipdcoefn)
        rm2=rm*2.0e0_knd
        lm2=(l-m)/2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        lim=limeta/2-ix
!
!  compute radial function of the first kind and its first derivative
!
!  backward series for denominator
        idenom=0
        denom=1.0e0_knd
        sumdp=1.0e0_knd
        if (lm2.eq.0) go to 20
        dold=1.0e0_knd
          do 10 j=lm2,1,-1
          jj=j+j+ix
          dnew=dold/(pratb(jj+1)*enr(j))
          denom=denom+dnew
          if(dnew.gt.0.0e0_knd) sumdp=sumdp+dnew
          if(abs(dnew/denom).lt.dec) go to 20
          dold=dnew
10        continue
20      continue
!
!  forward series for denominator
        dold=1.0e0_knd
          do 30 j=lm2+1,lim
          jj=j+j+ix
          dnew=dold*enr(j)*pratb(jj+1)
          denom=denom+dnew
          if(dnew.gt.0.0e0_knd) sumdp=sumdp+dnew
          if(abs(dnew/denom).lt.dec) go to 40
            if(abs(denom).gt.teste) then
            denom=denom*testeo
            dnew=dnew*testeo
            sumdp=sumdp*testeo
            idenom=idenom+nfac
            end if
25        dold=dnew
30        continue
40      continue
        jden=j
        ndens=0
        if(sumdp/denom.ne.0.0e0_knd) ndens=int(log10(abs(sumdp/denom)))
        if(ndens.lt.0) ndens=0
        if(ndens.gt.ndec) ndens=ndec
        iterm=int(log10(abs(denom)))
        idenom=idenom+iterm
        denom=denom*ten**(-iterm)
!
!  backward series for numerator
        dold=1.0e0_knd
        doldd1=1.0e0_knd
        doldd2=reld12
        r1temp=1.0e0_knd
        sumnp=1.0e0_knd
        r1dtemp=doldd2
        sumdnp=0.0e0_knd
        if(doldd2.gt.0.0e0_knd) sumdnp=doldd2
        if(l.ne.0) r1dtemp=r1dtemp+1.0e0_knd
        if(l.ne.0) sumdnp=sumdnp+1.0e0_knd
        if(lm2.lt.1) go to 60
          do 50 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sbesf(jj+m)*pratt(jj+1)*enr(j))
          dnewd1=-doldd1/(sbesf(jj+m)*pdratt(jj+1)*enr(j))
          dnewd2=-doldd2/(sbesdf(jj+m)*pratt(jj+1)*enr(j))
          r1temp=r1temp+dnew
          if(dnew.gt.0.0e0_knd) sumnp=sumnp+dnew
          dnewd=dnewd1+dnewd2
          r1dtemp=r1dtemp+dnewd
          if(dnewd.gt.0.0e0_knd) sumdnp=sumdnp+dnewd
          if(abs(dnew/r1temp)+abs(dnewd/r1dtemp).lt.dec) go to 60
          if(abs(r1temp).lt.teste) go to 45
          r1temp=r1temp*testeo
          dnew=dnew*testeo
          ir1tempe=ir1tempe+nfac
          r1dtemp=r1dtemp*testeo
          dnewd1=dnewd1*testeo
          dnewd2=dnewd2*testeo
          sumnp=sumnp*testeo
          iflag=0
          if(sumnp.lt.abs(r1temp*dec)) sumnp=r1temp*dec
          sumdnp=sumdnp*testeo
          if(sumdnp.lt.abs(r1dtemp*dec)) sumdnp=r1dtemp*dec
45        dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
50      continue
60      continue
        if(m.eq.0.and.jj.eq.2) r1dtemp=r1dtemp-dnewd1
!
!  forward series for numerator
          if(iflag.eq.0) then
          j=lm2
          go to 130
          end if
        dold=1.0e0_knd
        doldd1=1.0e0_knd
        doldd2=reld12
        dnewsum=0.0e0_knd
        dnewdsum=0.0e0_knd
        doldd=reld12
        if(l.ne.0) doldd=reld12+1.0e0_knd
          do 110 j=lm2+1,lim-1
          jj=j+j+ix
          dnew=-dold*enr(j)*sbesf(jj+m)*pratt(jj+1)
          dnewd1=-doldd1*enr(j)*sbesf(jj+m)*pdratt(jj+1)
          dnewd2=-doldd2*enr(j)*sbesdf(jj+m)*pratt(jj+1)
          r1temp=r1temp+dnew
          if(dnew.gt.0.0e0_knd) sumnp=sumnp+dnew
          if(dnew.ne.0.0e0_knd) test=abs(dnew/r1temp)
          dnewd=dnewd1+dnewd2
          r1dtemp=r1dtemp+dnewd
          if(dnewd.gt.0.0e0_knd) sumdnp=sumdnp+dnewd
          if(dnewd.ne.0.0e0_knd) testd=abs(dnewd/r1dtemp)
          if (test.lt.dec.and.testd.lt.dec) go to 130
            if(abs(r1temp).gt.teste) then
            r1temp=r1temp*testeo
            dnew=dnew*testeo
            sumnp=sumnp*testeo
            ir1tempe=ir1tempe+nfac
            r1dtemp=r1dtemp*testeo
            sumdnp=sumdnp*testeo
            dnewd1=dnewd1*testeo
            dnewd2=dnewd2*testeo
            end if
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
110       continue
130     jnum=j
        naccns1=0
        if(sumnp/r1temp.ne.0.0e0_knd) naccns1= &
                                      int(log10(abs(sumnp/r1temp)))
        if(naccns1.lt.0) naccns1=0
        if(naccns1.gt.ndec) naccns1=ndec
        naccns2=0
        if(sumdnp/r1dtemp.ne.0.0e0_knd) naccns2= &
                                        int(log10(abs(sumdnp/r1dtemp)))
        if(naccns2.lt.0) naccns2=0
        if(naccns2.gt.ndec) naccns2=ndec
        naccs1=max(ndens,naccns1)
        naccs2=max(ndens,naccns2)
!
!
!  combining results to form the radial function characteristics
!  r1c and r1dc and corresponding exponents ir1e and ir1de
        r1c=r1temp*sbesn(l+1)*pcoefn/denom
        iterm=0
        if(r1c.ne.0.0e0_knd) iterm=int(log10(abs(r1c)))
        ir1e=ir1tempe+ibese(l+1)+ipcoefn-idenom+iterm
        r1c=r1c*ten**(-iterm)
        r1dc=(r1dcoef1*r1dtemp*sbesn(l+1)*pdcoefn/denom)* &
             ten**(ibese(l+1)+ipdcoefn-idenom-ir1e+ir1tempe)
        iterm=0
        if(r1dc.ne.0.0e0_knd) iterm=int(log10(abs(r1dc)))
        ir1de=ir1e+iterm
     r1dc=r1dc*ten**(-iterm)
        ir1e=ir1e
        ir1de=ir1de
if (debug) then
        write(40,140) jnum,jden,lim
140     format(8x,'r1eta: numerator, denominator converged in ', &
               i6,' ,',i6,' terms; ',i6,' terms available.')
        if(naccns1+naccns2.gt.0) write(40,150) naccns1,naccns2
150     format(15x,'subtraction errors in numerator for r1 and r1d are', &
               1x,i2,' and ',i2,' digits.')
        if(ndens.gt.0) write(40,160) ndens
160     format(15x,'subtaction error in denominator is ',i2,' digits')
        if(abs(r1c).ge.1.0e0_knd) go to 170
end if
        r1c=r1c*ten
        ir1e=ir1e-1
170     continue
        if(abs(r1dc).ge.1.0e0_knd) go to 180
        r1dc=r1dc*ten
        ir1de=ir1de-1
180     continue
190     jeta=max(jden,jnum)
        return
        end subroutine
!
!
!
        subroutine r2int (l,m,c,x,limint,ndec,nex,maxd,enr,dc01,idc01, &
                          maxint,maxmp,maxlp,intlim,rpint1,rpint2, &
                          pint1,pint2,pint3,pint4,norme,pnorm,ipnorm, &
                          coefme,coefmo,ipint,r2c,ir2e,r2dc,ir2de,jint, &
                          coefn,icoefn,isub,isubd)
!
!
!  purpose:     To calculate values of the radial function of the
!               second kind and its first derivative using an integral
!               representation of the radial functions in terms of the
!               angular function of the first kind together with a
!               Neumann function kernal. The angular function is
!               expanded in a series of associated Legendre functions.
!               Gaussian quadrature is used (in subroutine pint) to
!               evaluate the resulting integrals involving associated
!               Legendre functions times the Neumann function kernel.
!               This subroutine performs the summation of the
!               integrals times d coefficients to obtain r2 and r2d.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x      : x
!               limint : approximately twice the maximum number of
!                        terms available to be taken in the series
!               ndec   : number of decimal digits available in
!                        real arithmetic
!               nex    : maximum exponent in real(knd) arithmetic
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               dc01   : characteristic of the first d coefficient,
!                        either d0 or d1, depending on whether l-m
!                        is even or odd
!               idc01  : exponent (base 10) of the first d coefficient
!               maxint : dimension of pint and rpint arrays
!               maxmp  : dimension of norme array
!               maxlp  : dimension of the pnorm and ipnorm arrays
!               intlim : highest order of Legendre function for which
!                        integrals are not totally inaccurate due to
!                        subtraction errors in their calculation.
!                        series for calculating r2 and r2d will not
!                        include contributions from terms involving
!                        integrals above order intlim
!               rpint1 : arrays of ratios of successive integrals of
!                        either the first or the third kind, depending
!                        on whether l-m is even or odd
!               rpint2 : array of ratios of successive integrals of
!                        either the second or the fourth kind,
!                        depending on whether l-m is even or odd
!               pint1  : array of scaled values for the integrals of
!                        the first kind
!               pint2  : array of scaled values for the integrals of
!                        the second kind
!               pint3  : array of scaled values for the integrals of
!                        the third kind
!               pint4  : array of scaled values for the integrals of
!                        the fourth kind
!               norme  : exponent used to scale the Neumann function
!                        of order m involved in the integrals
!               pnorm  : array of characteristics of the scaling factors
!                        used for the associated Legendre functions in
!                        the integrals to avoid overflow
!               ipnorm : array of exponents (base 10) corresponding to
!                        pnorm
!               coefme : coefficient used to multiply r2 to get one of
!                        the two contributions to r2d when l-m is even
!               coefmo : coefficient used to multiply r2 to get one of
!                        the two contributions to r2d when l-m is odd
!               ipint  : equal to zero the first time r2int is called;
!                        equal to unity otherwise
!
!     output:   r2c    : characteristic of oblate radial function
!                        of the second kind
!               ir2e   : exponent of oblate radial function of the
!                        second kind
!               r2dc   : characteristic of derivative with respect
!                        to x of oblate radial function of the second
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        oblate radial function of the second kind
!               jint   : maximum value of the index j in the forward
!                        sum for r2 and r2d, i.e., the highest enr(j)
!                        used
!               coefn  : characteristic of coefficient that is only
!                        calculated once (for l = m) and is then
!                        used for all values of l
!               icoefn : exponent for coefn
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,coefa,coefl,coefme,coefmo,coefn,dec,dcon,dnew,dnewd, &
                  dold,doldd,dc01,ri,rm,rm2,r2c,r2dc,r2dpos,rs,r2dtemp, &
                  r2pos,r2temp,ten,teste,testeo,x
        real(knd) enr(maxd),pnorm(maxlp),pint1(maxint),pint2(maxint), &
                  pint3(maxint),pint4(maxint),rpint1(maxint), &
                  rpint2(maxint)
!
!  integer arrays
        integer ipnorm(maxlp)
!
        ten=10.0e0_knd
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        rm=m
        rm2=rm+rm
        lm2=(l-m)/2
        ix=l-m-2*lm2
        ixx=ix-1
        ixx2=ixx+2
        lim=limint/2-ix-1
        if(limint.gt.intlim-2) lim=intlim/2-ix-1
!
!  compute the leading coefficient
        if(ipint.ne.0) go to 20
        icoefn=norme
        coefn=0.5e0_knd
        if(m.eq.0) go to 20
          do 10 i=1,m
          ri=i
      coefn=coefn/(ri+ri)
            if(coefn.lt.testeo) then
            coefn=coefn*teste
            icoefn=icoefn-nfac
            end if
10      continue
          iterm=int(log10(abs(coefn)))
          coefn=coefn*ten**(-iterm)
          icoefn=icoefn+iterm
20      continue
        if(ix.eq.0) coefa=(rm2+1.0e0_knd)*coefn
        if(ix.eq.1) coefa=(rm2+3.0e0_knd)*coefn
        if((ix.eq.0).and.(2*(lm2/2).ne.lm2)) coefa=-coefa
        if((ix.eq.1).and.(2*((l-m-1)/4).ne.(l-m-1)/2)) coefa=-coefa
    coefl=coefa/dc01
        icoefl=-idc01+icoefn
        dec=ten**(-ndec-1)
        dcon=dec
        jlow=lm2+1
!
!  compute the integrals involving the angular functions by summing
!  d coefficients times corresponding integrals of Legendre
!  functions
!
!  forward summation of series for r2 and r2d
        iflag=0
        jint=lim
        r2temp=0.0e0_knd
        r2dtemp=0.0e0_knd
        r2pos=0.0e0_knd
        r2dpos=0.0e0_knd
        if(jlow.gt.lim) go to 40
        dold=1.0e0_knd
        doldd=1.0e0_knd
        r2dtemp=doldd
        r2temp=dold
        r2pos=dold
        r2dpos=doldd
          do 30 j=jlow,lim
          dnew=dold*enr(j)*rpint1(j+j+ixx2)
          dnewd=doldd*enr(j)*rpint2(j+j+ixx2)
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(dnew.gt.0.0e0_knd) r2pos=r2pos+dnew
          if(dnewd.gt.0.0e0_knd) r2dpos=r2dpos+dnewd
          if((abs(dnew/r2temp)+abs(dnewd/r2dtemp)).lt.dcon) go to 40
          dold=dnew
          doldd=dnewd
30        continue
!
!  backward summation of series for r2 and r2d
40      jint=min(j,lim)
        if(j.eq.0) jint=lim
        if(lm2.lt.1.or.iflag.eq.1) go to 70
        dold=1.0e0_knd
        doldd=1.0e0_knd
        j=lm2
          do 60 jj=1,lm2
          dnew=dold/(rpint1(j+j+ixx2)*enr(j))
          dnewd=doldd/(rpint2(j+j+ixx2)*enr(j))
          if(j.gt.lim) go to 50
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(dnew.gt.0.0e0_knd) r2pos=r2pos+dnew
          if(dnewd.gt.0.0e0_knd) r2dpos=r2dpos+dnewd
          if((abs(dnew/r2temp)+abs(dnewd/r2dtemp)).lt.dcon) go to 70
50        dold=dnew
          doldd=dnewd
          j=j-1
60        continue
70      continue
        isub=int(log10(abs(r2pos/r2temp)+dec))
        if(isub.lt.0) isub=0
        isubd=int(log10(abs(r2dpos/r2dtemp)+dec))
        if(isubd.lt.0) isubd=0
    r2temp=r2temp*coefl*pnorm(l-m+1)
        if(ix.eq.0) r2temp=r2temp*pint1(l-m+1)
        if(ix.eq.1) r2temp=x*r2temp*pint3(l-m+1)
        iterm=int(log10(abs(r2temp)))
        ir2e=iterm+ipnorm(l-m+1)+icoefl
        r2c=r2temp*ten**(-iterm)
        if(abs(r2c).ge.1.0e0_knd) go to 80
        r2c=r2c*ten
        ir2e=ir2e-1
80      r2dtemp=-r2dtemp*coefl*pnorm(l-m+1)*c*x
        rs=r2dtemp
        if(ix.eq.0) r2dtemp=r2dtemp*pint2(l-m+1)+r2temp*coefme
        if(ix.eq.1) r2dtemp=r2dtemp*pint4(l-m+1)*x+r2temp*coefmo
        jsuba=0
        if(ix.eq.0.and.m.ne.0) jsuba=int(log10(abs(r2temp*coefme/ &
                                     r2dtemp)+dec))
        if(ix.eq.1) jsuba=int(log10(abs(r2temp*coefmo/r2dtemp)+dec))
        jsubb=0
        if(ix.eq.0.and.m.ne.0) jsubb=int(log10(abs(rs*pint2(l-m+1)/ &
                                     r2dtemp)+dec))
        if(ix.eq.1) jsubb=int(log10(abs(rs*x*pint4(l-m+1)/r2dtemp)+ &
                           dec))
        if(m.ne.0.or.ix.eq.1) isubd=max(isub+jsuba,isubd+jsubb,0)
if (debug) then
        write(40,90) jint,lim,isub,isubd
90      format(8x,'r2int: converged in ',i6,' terms; 'i6, &
               ' available; ',i2,' and ',i2,' digits of sub. error.')
end if
        jterm=int(log10(abs(r2dtemp)))
        ir2de=jterm+ipnorm(l-m+1)+icoefl
        r2dc=r2dtemp*ten**(-jterm)
        if(abs(r2dc).ge.1.0e0_knd) go to 100
        r2dc=r2dc*ten
        ir2de=ir2de-1
100     continue
        return
        end subroutine
!
!
        subroutine r2leg (l,m,c,x,lnum,minacc,limleg,limdr,iflagp,ndec, &
                          nex,maxd,maxmp,maxpdr,maxdr,maxq,enr,enrneg, &
                          drhor,nsdrhor1,dc01,idc01,dneg,idneg,nsdneg, &
                          dfnorm,idfe,dmfnorm,idmfe,prx,pdrx,qdr,qdqr, &
                          qdm1,iqdm1,qdl,iqdl,qr,qm1,iqm1,ql,iql,fajo, &
                          ifajo,ifsub,jsub,termpq,itermpq,ioppsum,iopqnsum, &
                          r1c,ir1e,r1dc,ir1de,naccr1,itestm,r2c,ir2e, &
                          r2dc,ir2de,jleg,jlegp,jflagl,naccleg,kflagl, &
                          nsubleg,nsubdleg)
!
!  purpose:     To evaluate the oblate radial function of the
!               second kind and its first derivative with respect
!               to x using the traditional expansion in associated
!               Legendre functions.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               x       : radial coordinate x
!               lnum    : number of l values desired
!               minacc  : desired accuracy in decimal digits
!               limleg  : approximately twice the maximum number
!                         of terms available to be taken in qsum,
!                         (sum involving q's time d coefficients)
!               limdr   : maximum number of terms available to be
!                         taken in psum (sum involving p's time
!                         d rho coefficients)
!               iflagp  : integer flag set = 1 if psum series converges
!                         fully; set = 0 otherwise
!               ndec    : number of decimal digits available in
!                         real arithmetic
!               nex     : maximum exponent is real(knd) arithmetic
!               maxd    : dimension of enr array
!               maxmp   : dimension of enrneg array
!               maxpdr  : dimension of prx and pdrx arrays
!               maxdr   : dimension of drhor array
!               maxq    : dimension of qr and qdr arrays
!               enr     : array of d coefficient ratios
!               enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               nsdrhor1: subtraction error in calculating drhor(1)
!                         from drhor(2)
!               dc01    : characteristic of the ratio of the first d
!                         coefficient with nonnegative subscript, either
!                         d0 or d1 depending on whether l-m is even or
!                         odd, to the d coefficient with subscript l - m
!               idc01   : exponent (base 10) corresponding to dc01
!               dneg    : characteristic of the ratio of the d
!                         coefficient with subscript -2m+ix to the
!                         d coefficient with subscript ix, where
!                         ix = 0 or 1 depending on whether l - m
!                         is even or odd
!               idneg   : exponent corresponding to dneg
!               nsdneg  : subtraction error in calculating dneg
!               dfnorm  : characteristic of Flammer normalization sum of
!                         d coefficients. equal to the reciprocal of
!                         the value of the d coefficient d(n = l - m) using
!                         this normalization for the angular functions
!               idfe    : exponent associated with dfnorm
!               dmfnorm : characteristic of Morse-Feshbach normalization
!                         sum of the d coefficients. equal to the
!                         reciprocal of the value of the d coefficient
!                         d(n = l - m) using this normalization for the
!                         angular functions
!               idmfe   : exponent associated with dmfnorm
!               prx     : ratios of successive Legendre functions of
!                         the first kind of the same parity
!               pdrx    : ratios of successive first derivatives of
!                         Legendre functions of the first kind of the
!                         same parity
!               qdr     : ratios of first derivatives of successive
!                         Legendre functions of the second kind
!               qdqr    : array of ratios of derivatives of associated
!                         Legendre functions of the second kind to the
!                         corresponding Legendre function for degrees
!                         from -m to m-1c
!               qdm1    : characteristic of the first derivative of
!                         the associated Legendre function of the second
!                         kind with order m and degree m-1
!               iqdm1   : exponent corresponding to qdm1
!               qdl     : array of characteristics of the first
!                         derivatives of the associated Legendre
!                         functions of the second kind with order m
!                         and degrees from m to m+lnum-1, scaled by
!                                        -m/2
!                         (2m-1)!!(x*x+1)
!               iqdl    : array of exponents corresponding to qdl
!               qr      : array of ratios of successive associated
!                         Legendre functions of the second kind
!               qm1     : characteristic of the associated Legendre
!                         function of the second kind with order m
!                         and degree m-1
!               iqm1    : exponent corresponding to qm1
!               ql      : array of characteristics of the associated
!                         Legendre function of the second kind with
!                         order m and degrees from m to m+lnum-1
!                                                  -m/2
!                         scaled by (2m-1)!!(x*x+1)
!               iql     : array of exponents corresponding to ql
!               fajo    : characteristic of the joining factor of the
!                         second kind
!               ifajo   : exponent corresponding to fajo
!               ifsub   : subtraction error in forming fajo coming from
!                         dfnorm, dmfnorm, and dneg
!               jsub    : larger of the subtraction errors in the
!                         Flammer and the Morse and Feshbach
!                         normalizations
!               termpq  : characteristic of the relative size of the
!                         maximum terms in the positive degree q series
!                         and the p series used to calculate r2 and r2d
!               itermpq : exponent corresponding to termpq
!               ioppsum : integer flag = 0 if psum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if psum is computed
!               iopqnsum: integer flag = 0 if qnsum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if qnsum is computed
!               r1c     : characteristic of the radial function of the
!                         first kind
!               ir1e    : exponent corresponding to r1c
!               r1dc    : characteristic of the first derivative of the
!                         radial function of the first kind
!               ir1de   : exponent corresponding to r1dc
!               naccr1  : estimated accuracy of r1c and r1dc
!               itestm  : number of digits of match between the forward
!                         and backward recursions in calculating enr
!
!     output:   r2c     : characteristic of oblate
!                         radial function of the second kind
!               ir2e    : exponent of oblate radial function of the
!                         second kind
!               r2dc    : characteristic of derivative with
!                         respect to x of oblate radial function
!                         of the second kind
!               ir2de   : exponent of derivative with respect to x of
!                         oblate radial function of second kind
!               jleg    : maximum number of terms taken in qsum
!               jlegp   : maximum number of terms taken in psum
!               jflagl  : equal to 1 if a more accurate value
!                         for the leading coefficient drhor(1)
!                         in psum is obtained using the Wronskian;
!                         equal to 0 otherwise.
!               naccleg : Wronskian estimate if jflagl = 0; estimate
!                         of accuracy when jflag1 = 1
!               kflagl  : equal to one if either qsum or psum
!                         becomes so large that the summation
!                         is exited and r2 and r2d are set equal
!                         to 10.0e0_knd**nex; equal to 0 otherwise
!               nsubleg : subtraction error in decimal digits
!                         in the calculation of r2c
!               nsubdleg: subtraction error in decimal digits
!                         in the calculation of r2dc
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,dconp,dconq,dconqn,dec,dfnorm,dmfnorm,dneg,dnegjf, &
                  dnew,dnewd,dold,doldd,dc01,psum,psump,pdsum,pdsump, &
                  qdm1,qndsum,qndsump,qdsum,qdsump,qm1,qnsum,qnsump, &
                  qsum,qsump,r1c,r1dc,r2c,r2dc,rm,spsum,spsump,spdsum, &
                  spdsump,ten,termpq,test,testd,testm,testdm,testp,tm, &
                  wronc,wronca,wroncb,wront,x,xden,xcoef,xrhs
        real(knd) anumt1,anumt2,anumt3,anumt4,dent1,dent2          
        real(knd) drhor(maxdr),enr(maxd),enrneg(maxmp),fajo(lnum+1), &
                  prx(maxpdr),pdrx(maxpdr),qdl(lnum),qdr(maxq), &
                  qdqr(maxmp),ql(lnum),qr(maxq)
!
!  integer arrays
        integer ifajo(lnum+1),iqdl(lnum),iql(lnum)
!       
        ten=10.0e0_knd
        dec=ten**(-ndec-1)
        dconp=dec
        testp=ten**(nex-10)
        lm2=(l-m)/2
        ix=l-m-2*lm2
        imxp=m+m+ix
        ixx=1-ix
        lim1=limleg/2-ix
        lim2=limdr-1
        wront=1.0e0_knd/(c*(x*x+1.0e0_knd))
        if(ioppsum.eq.0) lim2=0
        rm=m
        tm=rm+rm
        dconq=dec
        dconqn=dec
        dnegjf=dneg*dc01
        if(m.eq.0) dnegjf=dc01
        iterm=int(log10(abs(dnegjf)))
        dnegjf=dnegjf*ten**(-iterm)
        idnegjf=idneg+idc01+iterm
        if(m.eq.0) idnegjf=idc01+iterm
        fajo(l-m+1)=fajo(l-m+1)*dmfnorm*dnegjf/dfnorm
        iterm=int(log10(abs(fajo(l-m+1))))
        fajo(l-m+1)=fajo(l-m+1)*ten**(-iterm)
        ifajo(l-m+1)=ifajo(l-m+1)+idnegjf+iterm+idmfe-idfe
!
!  begin calculation of series for r2
!
!  calculate d*q sum over positive n using pyramid summation
!
!  backward summation
        qsum=1.0e0_knd
        qdsum=1.0e0_knd
        qsump=1.0e0_knd
        qdsump=1.0e0_knd
        if(lm2.eq.0) go to 20
        dold=1.0e0_knd
        doldd=1.0e0_knd
        j=lm2
          do 10 jj=1,lm2
          dnew=-dold/(qr(j+j+imxp)*qr(j+j+imxp-1)*enr(j))
          qsum=qsum+dnew
          if(dnew.gt.0.0e0_knd) qsump=qsump+dnew
          dnewd=-doldd/(qdr(j+j+imxp)*qdr(j+j+imxp-1)*enr(j))
          qdsum=qdsum+dnewd
          if(dnewd.gt.0.0e0_knd) qdsump=qdsump+dnewd
            if(abs(qsum).gt.testp.or.abs(qdsum).gt.testp) then
            r2c=1.0e0_knd
            r2dc=1.0e0_knd
            ir2e=nex
            ir2de=nex
            nsubleg=ndec
            nsubdleg=ndec
            jleg=j
            jlegp=0
            jnleg=0
            kflagl=1
            go to 180
            end if
          if((abs(dnew/qsum)+abs(dnewd/qdsum)).lt.dconq) go to 20
          dold=dnew
          doldd=dnewd
          j=j-1
10        continue
20      continue
!
!  forward summation
        jlow=lm2+1
        dold=1.0e0_knd
        doldd=1.0e0_knd
          do 30 j=jlow,lim1
          dnew=-dold*enr(j)*qr(j+j+imxp)*qr(j+j+imxp-1)
          qsum=qsum+dnew
          if(dnew.gt.0.0e0_knd) qsump=qsump+dnew
          dnewd=-doldd*enr(j)*qdr(j+j+imxp)*qdr(j+j+imxp-1)
          qdsum=qdsum+dnewd
          if(dnewd.gt.0.0e0_knd) qdsump=qdsump+dnewd
          if((abs(dnew/qsum)+abs(dnewd/qdsum)).lt.dconq) go to 40
          dold=dnew
          doldd=dnewd
30        continue
40      continue
        jleg=j
        if(jleg.gt.lim1) jleg=lim1
        nsqsum=0
        if(qsum*qsump.ne.0.0e0_knd) nsqsum=int(log10(abs(qsump/qsum)))
        if(nsqsum.gt.ndec) nsqsum=ndec
        if(nsqsum.lt.0) nsqsum=0
        nsqdsum=0
        if(qdsum*qdsump.ne.0.0e0_knd) nsqdsum= &
                                  int(log10(abs(qdsump/qdsum)))
        if(nsqdsum.gt.ndec) nsqdsum=ndec
        if(nsqdsum.lt.0) nsqdsum=0
        qsum=qsum*ql(l-m+1)/(fajo(l-m+1)*termpq)
        iterm=int(log10(abs(qsum)))
        qsum=qsum*(ten**(-iterm))
        iqsum=iql(l-m+1)-ifajo(l-m+1)-itermpq+iterm
        qdsum=qdsum*qdl(l-m+1)/(fajo(l-m+1)*termpq)
        iterm=int(log10(abs(qdsum)))
        qdsum=qdsum*(ten**(-iterm))
        iqdsum=iqdl(l-m+1)-ifajo(l-m+1)-itermpq+iterm
          if(2*(m/2).eq.m) then
          qsum=-qsum
          qdsum=-qdsum
          end if
          if(2*((l-m+1)/4).ne.(l-m+1)/2) then
          qsum=-qsum
          qdsum=-qdsum
          end if
45      continue
!
!  calculate d*q sum over negative n
        qnsum=0.0e0_knd
        qndsum=0.0e0_knd
        qnsump=0.0e0_knd
        qndsump=0.0e0_knd
        iqnsum=0
        iqndsum=0
        j2=0
        nmterm=0
        nsqnsum=0
        nsqndsum=0
        if(iopqnsum.eq.0.or.m.eq.0) go to 90
        nmterm=m
        qnsum=enrneg(m)
        qndsum=qnsum*qdqr(m+m)
        j2=1
        if(ix.eq.1) go to 50
        qnsum=qnsum*qr(m+m-1)
        qndsum=qnsum*qdqr(m+m-1)
50      if(qnsum.gt.0.0e0_knd) qnsump=qnsum
        if(qndsum.gt.0.0e0_knd) qndsump=qndsum
          if(m.eq.1) then
          jnleg=1
          go to 80
          end if
        dold=qnsum
          do 60 j=2,m
          dnew=-dold*enrneg(m-j+1)*qr(imxp-j-j+1)*qr(imxp-j-j+2)
          qnsum=qnsum+dnew
          if(dnew.gt.0.0e0_knd) qnsump=qnsump+dnew
          dnewd=dnew*qdqr(imxp-j-j+1)
          qndsum=qndsum+dnewd
          if(dnewd.gt.0.0e0_knd) qndsump=qndsump+dnewd
          if(abs(dnew/qnsum)+abs(dnewd/qndsum).lt.dec) go to 70
          dold=dnew
60        continue
70      jnleg=j
        if(jnleg.gt.m) jnleg=m
        nsqnsum=0
80      if(qnsum*qnsump.ne.0.0e0_knd) nsqnsum= &
                          int(log10(abs(qnsump/qnsum)))
        if(nsqnsum.gt.ndec) nsqnsum=ndec
        if(nsqnsum.lt.0) nsqnsum=0
        nsqnsum=max(nsqnsum,nsdneg)
        nsqndsum=0
        if(qndsum*qndsump.ne.0.0e0_knd) nsqndsum= &
                          int(log10(abs(qndsump/qndsum)))
        if(nsqndsum.gt.ndec) nsqndsum=ndec
        if(nsqndsum.lt.0) nsqndsum=0
        nsqndsum=max(nsqndsum,nsdneg)
        qnsum=qnsum*qm1*dc01/(fajo(l-m+1)*termpq)
        iterm=int(log10(abs(qnsum)))
        qnsum=qnsum*(ten**(-iterm))
        iqnsum=iqm1+idc01-ifajo(l-m+1)-itermpq+iterm
        qnsum=qnsum*(ten**(iqnsum-iqsum))
        qndsum=qndsum*qm1*dc01/(fajo(l-m+1)*termpq)
        iterm=int(log10(abs(qndsum)))
        qndsum=qndsum*(ten**(-iterm))
        iqndsum=iqm1+idc01-ifajo(l-m+1)-itermpq+iterm
        qndsum=qndsum*(ten**(iqndsum-iqdsum))
          if(2*(m/2).ne.m) then
          qnsum=-qnsum
          qndsum=-qndsum
          end if
90      continue
!
!       calculate d(rho|n)*p summation
        psum=0.0e0_knd
        pdsum=0.0e0_knd
        ipsum=0
        ipdsum=0
        jlegp=0
        nspsum=0
        nspdsum=0
        if(ioppsum.eq.0) go to 160
        psum=prx(ixx+1)*drhor(1)
        pdsum=pdrx(ixx+1)*drhor(1)
        dold=psum
        doldd=pdsum
        if(m.ne.0.or.ix.ne.1) go to 100
        pdsum=0.0e0_knd
        doldd=drhor(1)
100     continue
        spsum=psum
        spdsum=pdsum
        psump=0.0e0_knd
        if(psum.gt.0.0e0_knd) psump=psum
        pdsump=0.0e0_knd
        if(pdsum.gt.0.0e0_knd) pdsump=pdsum
        spsump=psump
        spdsump=pdsump
        testm=1.0e0_knd
        testdm=1.0e0_knd
        iflagp=0
        jlegpf=1
        jlegpd=1
          do 130 j=2,lim2
          dnew=dold*drhor(j)*prx(j+j-ix)
          psum=psum+dnew
          if(dnew.gt.0.0e0_knd) psump=psump+dnew
          dnewd=doldd*drhor(j)*pdrx(j+j-ix)
          pdsum=pdsum+dnewd
          if(dnewd.gt.0.0e0_knd) pdsump=pdsump+dnewd
            if(abs(psum).gt.testp.or.abs(pdsum).gt.testp) then
            r2c=1.0e0_knd
            r2dc=1.0e0_knd
            ir2e=nex
            ir2de=nex
            nsubleg=ndec
            nsubdleg=ndec
            jlegp=j
            jnleg=0
            kflagl=1
            go to 180
            end if
          test=abs(dnew/psum)
          testd=abs(dnewd/pdsum)
          if(test.gt.testm.or.test.eq.0.0e0_knd) go to 110
          testm=test
          spsum=psum
          spsump=psump
          jlegpf=j
110       if(testd.gt.testdm.or.testd.eq.0.0e0_knd) go to 120
          testdm=testd
          spdsum=pdsum
          spdsump=pdsump
          jlegpd=j
120       if(test+testd.lt.dconp) go to 140
          dold=dnew
          doldd=dnewd
130       continue
        go to 150
140     continue
        iflagp=1
150     jlegp=j
        if(jlegp.gt.lim2) jlegp=lim2
        jlegp=max(jlegpf,jlegpd)
        psum=spsum
        pdsum=spdsum
        psump=spsump
        pdsump=spdsump
        nspsum=0
        if(psum*psump.ne.0.0e0_knd) nspsum=int(log10(abs(psump/psum)))
        if(nspsum.gt.ndec) nspsum=ndec
        if(nspsum.lt.0) nspsum=0
        nspdsum=0
        if(pdsum*pdsump.ne.0.0e0_knd) nspdsum= &
                          int(log10(abs(pdsump/pdsum)))
        if(nspdsum.gt.ndec) nspdsum=ndec
        if(nspdsum.lt.0) nspdsum=0
        psum=-psum*dnegjf*termpq/fajo(l-m+1)
        iterm=0
        if(psum.ne.0.0e0_knd) iterm=int(log10(abs(psum)))
        psum=psum*(ten**(-iterm))
        ipsum=idnegjf+itermpq-ifajo(l-m+1)+iterm
        psum=psum*(ten**(ipsum-iqsum))
        pdsum=pdsum*dnegjf*termpq/fajo(l-m+1)
        if(m.ne.0) pdsum=-pdsum*rm*x/(x*x+1.0e0_knd)
        iterm=0
        if(pdsum.ne.0.0e0_knd) iterm=int(log10(abs(pdsum)))
        pdsum=pdsum*(ten**(-iterm))
        ipdsum=idnegjf+itermpq-ifajo(l-m+1)+iterm
        pdsum=pdsum*(ten**(ipdsum-iqdsum))
        if(2*((l-m)/2).eq.(l-m)) pdsum=-pdsum
160     continue
        r2c=qsum+qnsum+psum
        r2dc=qdsum+qndsum+pdsum
        wronca=r1c*r2dc*(ten**(ir1e+iqdsum))
        wroncb=r1dc*r2c*(ten**(ir1de+iqsum))
        wronc=wronca-wroncb
        naccleg=-int(log10(abs((wronc-wront)/wront)+dec))
        nstest=max(nspsum,nspdsum)
        if(naccleg.gt.ndec) naccleg=ndec
          if(nsdrhor1.gt.0.and.naccleg.lt.minacc.and.x.le.0.01e0_knd &
              .and.(naccleg.gt.1.or.nstest.lt.10)) then
          anumt1=qdsum*r1c*ten**(ir1e+iqdsum)
          anumt2=qndsum*r1c*ten**(ir1e+iqdsum)
          anumt3=qsum*r1dc*ten**(ir1de+iqsum)
          anumt4=qnsum*r1dc*ten**(ir1de+iqsum)
          numc1=-ndec
          if(abs(anumt1).ne.0.0e0_knd) numc1=int(log10(abs(anumt1/wront)))
          numc2=-ndec
          if(abs(anumt2).ne.0.0e0_knd) numc2=int(log10(abs(anumt2/wront)))
          numc3=-ndec
          if(abs(anumt3).ne.0.0e0_knd) numc3=int(log10(abs(anumt3/wront)))
          numc4=-ndec
          if(abs(anumt4).ne.0.0e0_knd) numc4=int(log10(abs(anumt4/wront)))
          nacct1=ndec-(max(ifsub,nsqdsum)+numc1)
          if(nacct1.gt.ndec) nacct1=ndec
          nacct2=ndec-(max(ifsub,nsqndsum)+numc2)
          if(nacct2.gt.ndec) nacct2=ndec
          nacct3=ndec-(max(ifsub,nsqsum)+numc3)
          if(nacct3.gt.ndec) nacct3=ndec
          nacct4=ndec-(max(ifsub,nsqnsum)+numc4)
          if(nacct4.gt.ndec) nacct4=ndec 
          naccnum=min(nacct1,nacct2,nacct3,nacct4)
          if(naccnum.lt.0) naccnum=0
          dent1=r1c*pdsum*ten**(ir1e+iqdsum)
          dent2=r1dc*psum*ten**(ir1de+iqsum)
          nratio=0
          if(abs(dent1*dent2).ne.0.0e0_knd) &
              nratio=int(log10(abs(dent1/dent2)))
            if(nratio.gt.0) then
            naccd1=ndec-nspdsum
            naccd2=ndec-nspsum+nratio
            else
            naccd2=ndec-nspsum
            naccd1=ndec-nspdsum-nratio
            end if
          nacclest=min(naccnum,naccd1,naccd2)  
          if(nacclest.gt.naccr1) nacclest=naccr1
          if(nacclest.gt.itestm-2) nacclest=itestm-2 
          if(nacclest.lt.0) nacclest=0
            if(nacclest.gt.naccleg) then
            xrhs=wront-(qdsum+qndsum)*r1c*ten**(ir1e+iqdsum)+ &
                    (qsum+qnsum)*r1dc*ten**(ir1de+iqsum)
            xden=(r1c*pdsum*ten**(ir1e+iqdsum)-r1dc*psum* &
                    ten**(ir1de+iqsum))
            xcoef=xrhs/xden
            psum=psum*xcoef
            pdsum=pdsum*xcoef
            r2c=qsum+qnsum+psum
            r2dc=qdsum+qndsum+pdsum
            jflagl=1
            naccleg=nacclest
            end if
          end if
          if(jflagl.eq.0) then
          nspsum=max(nspsum,nsdrhor1)
          nspdsum=max(nspdsum,nsdrhor1)
          end if
        nqs=0
        if(qsum/r2c.eq.0.0e0_knd) nqs=-ndec
        if(qsum/r2c.ne.0.0e0_knd) nqs=int(log10(abs(qsum/r2c)))
        nqns=0
        if(qnsum/r2c.eq.0.0e0_knd) nqns=-ndec
        if(m.ne.0.and.iopqnsum.ne.0.and.qnsum/r2c.ne.0.0e0_knd) &
                      nqns=int(log10(abs(qnsum/r2c)))
        nps=0
        if(psum/r2c.eq.0.0e0_knd) nps=-ndec
        if(ioppsum.ne.0.and.psum/r2c.ne.0.0e0_knd) &
                      nps=int(log10(abs(psum/r2c)))
        nsqsum=max(nsqsum,nsdneg)+nqs
        if(nsqsum.lt.0) nsqsum=0
        if(nsqsum.gt.ndec) nsqsum=ndec
        nspsum=nspsum+nps
        if(nspsum.lt.0) nspsum=0
        if(nspsum.gt.ndec) nspsum=ndec
        nsqnsum=max(nsqnsum,nsdneg)+nqns
        if(nsqnsum.lt.0) nsqnsum=0
        if(nsqnsum.gt.ndec) nsqnsum=ndec
        nsubleg=max(nsqsum,nsqnsum,nspsum,jsub)
        r2dc=qdsum+qndsum+pdsum
        nqds=0
        if(qdsum/r2dc.eq.0.0e0_knd) nqds=-ndec
        if(qdsum/r2dc.ne.0.0e0_knd) nqds=int(log10(abs(qdsum/r2dc)))
        nqnds=0
        if(qndsum/r2dc.eq.0.0e0_knd) nqnds=-ndec
        if(m.ne.0.and.iopqnsum.ne.0.and.qndsum/r2dc.ne.0.0e0_knd) &
                      nqnds=int(log10(abs(qndsum/r2dc)))             
        if(nqns.lt.(-ndec-1).and.nqnds.lt.(-ndec-1)) iopqnsum=0
        if(qnsum.eq.0.0e0_knd.and.qndsum.eq.0.0e0_knd) iopqnsum=0
        npds=0
        if(pdsum/r2dc.eq.0.0e0_knd) npds=-ndec
        if(ioppsum.ne.0.and.pdsum/r2dc.ne.0.0e0_knd) &
                      npds=int(log10(abs(pdsum/r2dc)))
        if(nps.lt.(-ndec-1).and.npds.lt.(-ndec-1)) ioppsum=0
        if(psum.eq.0.0e0_knd.and.pdsum.eq.0.0e0_knd) ioppsum=0
        nsqdsum=max(nsqdsum,nsdneg)+nqds
        if(nsqdsum.lt.0) nsqdsum=0
        if(nsqdsum.gt.ndec) nsqdsum=ndec
        nspdsum=nspdsum+npds
        if(nspdsum.lt.0) nspdsum=0
        if(nspdsum.gt.ndec) nspdsum=ndec
        nsqndsum=max(nsqndsum,nsdneg)+nqnds
        if(nsqndsum.lt.0) nsqndsum=0
        if(nsqndsum.gt.ndec) nsqndsum=ndec
        nsubdleg=max(nsqdsum,nsqndsum,nspdsum,jsub)
        if(jflag.eq.1) naccleg=min(naccleg,ndec-max(nsqsum,nsqnsum, &
                         nsqdsum,nsqndsum))
        iterm=int(log10(abs(r2c)))
        r2c=r2c*(ten**(-iterm))
        ir2e=iqsum+iterm
        if(abs(r2c).ge.1.0e0_knd) go to 170
        r2c=r2c*ten
        ir2e=ir2e-1
170     continue
        iterm=int(log10(abs(r2dc)))
        r2dc=r2dc*(ten**(-iterm))
        ir2de=iqdsum+iterm
        if(abs(r2dc).ge.1.0e0_knd) go to 180
        r2dc=r2dc*ten
        ir2de=ir2de-1
180     continue
if (debug) then
        if(ioppsum.eq.1.and.iopqnsum.eq.1) write(40,190) jleg,jlegp, &
                                  jnleg,lim1,lim2,m,nsubleg,nsubdleg
190     format(8x,'r2leg: qsum, psum and qnsum series converged in ',i6, &
              ',' i6,' and ',i4,' terms; ',i6,',' i6,' and ' i4, &
              ' terms avail.',/,14x,i2,' and ',i2,' digits of sub.', &
              ' error in r2 and r2d.')
        if(ioppsum.eq.1.and.iopqnsum.eq.0) write(40,200) jleg,jlegp, &
                                           lim1,lim2,nsubleg,nsubdleg
200     format(8x,'r2leg: qsum and psum series converged in ',i6, &
              ' and ',i6,' terms; ',i6,' and ',i6,' terms avail.',/, &
              14x,i2,' and ',i2,' digits of sub. error in r2 and r2d;', &
              ' qnsum is negligible.')
        if(ioppsum.eq.0.and.iopqnsum.eq.1) write(40,210) jleg,jnleg, &
                                           lim1,m,nsubleg,nsubdleg
210     format(8x,'r2leg: qsum and qnsum series converged in ',i6, &
              ' and ',i4,' terms; ',i6,' and ',i4,' terms avail.',/, &
               14x,i2,' and ',i2,' digits of sub. error in r2 and r2d;' &
               ' psum is negligible.')
        if(ioppsum.eq.0.and.iopqnsum.eq.0) write(40,220) jleg,lim1, &
                                                  nsubleg,nsubdleg
220     format(8x,'r2leg: qsum series converged in ',i6,' terms with ', &
               i6,' terms avail.; 'i2,' and ',i2,' digits of',/,15x, &
               'sub. error in r2 and r2d; psum and qnsum are ', &
               'negligible.')
        if(jflagl.eq.1) write(40,230)
230     format(15x,'Wronskian used to improve accuracy of the', &
                ' the leading psum coefficient drhor(1).')
end if
        return
        end subroutine
!
!
        subroutine r2leg1(l,m,c,x,limq,maxq,ndec,eigval,qr,qdr,qm0, &
                          qdm0,r1c,ir1e,r2c,ir2e,r2dc,ir2de,jleg1)
!
!  purpose:     To evaluate the oblate radial function of the
!               second kind and its first derivative with respect
!               to x using an expansion in associated Legendre
!               functions of the second kind. This expansion is
!               due to Baber and Hasse.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               x       : radial coordinate x
!               limq    : the maximum number of terms available
!               maxq    : dimension of qr,qdr,aratio,coef1,
!                         coef2,coef3
!               ndec    : number of decimal digits available in real
!                         arithemetic
!               eigval  : eigenvalue
!               qr      : array of ratios of successive associated
!                         Legendre functions of the second kind
!                                                 m    m
!                         beginning with qr(1) = Q  / Q . Note that
!                                                 1    0
!                         this differs from the choice used in r2leg
!                         where functions with negative degree are
!                         also used and where ratios involving degrees
!                         less than m are inverted.
!               qdr     : ratios of first derivatives of successive
!                         Legendre functions of the second kind
!                                                   m     m
!                         beginning with qdr(1) = Q'  / Q' . Note that
!                                                   1     0
!                         this differs from the choice used in r2leg
!                         where functions with negative degree are
!                         also used and where ratios involving a degree
!                         less than m are inverted.
!               qm0     : associated Legendre function of the second
!                         kind with order m and degree 0.
!               qdm0    : first derivative of the associated Legendre
!                         function of the second kind with order m and
!                         degree 0.
!               r1c     : charcteristic of the corresponding radial
!                         function of the first kind
!               ir1e    : exponent of the corresponding radial function
!                         of the first kind
!
!     output:   r2c     : characteristic of oblate
!                         radial function of the second kind
!               ir2e    : exponent of oblate radial function of the
!                         second kind
!               r2dc    : characteristic of derivative with
!                         respect to x of the oblate radial function
!                         of the second kind
!               ir2de   : exponent of derivative with respect to x of
!                         the oblate radial function of second kind
!               jleg1   : number of terms taken in the series
!
        use param
!
!  real(knd) scalars
        real(knd) arg,c,dec,eigval,em,ea,e3,qm0,qdm0,ra,r1c,r1,r2,r2c, &
                  r2d,r2dc,r3,sumi,sumdi,sumr,sumdr,ten,term,termi, &
                  termdi,termr,termdr,x,xx
!  real(knd) vectors
        real(knd) qr(maxq),qdr(maxq),aratio(maxq),coef1(maxq), &
                coef2(maxq),coef3(maxq)
!
        ten=10.0e0_knd
        em=m
        m2=m+m
        dec=ten**(-ndec-1)
        arg=c*x
        xx=x*x+1.0e0_knd
        r1=r1c*(ten**(ir1e))
        ea=sin(arg)/c
        ra=cos(arg)/c
          if(m.ne.0) then
          ea=ea/em
          ra=ra/em
          end if
        la=l-(l/4)*4+1
          if(la.eq.1) then
          r3=ra
          e3=-ea
          end if
          if(la.eq.2) then
          e3=ra
          r3=ea
          end if
          if(la.eq.3) then
          r3=-ra
          e3=ea
          end if
          if(la.eq.4) then
          e3=-ra
          r3=-ea
          end if
50      m2=m+m
          do 90 j=1,limq
          n=j-m-1
          coef1(j)=c*2*(n+m+1)*(n+m2+1)/(n+n+m2+3)
          term=(n+m)*(n+m+1)
          coef2(j)=term-eigval-c*c
          coef3(j)=c*(n+n)*(n+m)/(n+n+m2-1)
90        continue
        aratio(1)=coef2(1)/coef1(1)
        aratio(limq)=0.0e0_knd
          do 100 j=limq-1,m+2,-1
          aratio(j-1)=coef3(j)/(coef1(j)*aratio(j)-coef2(j))
          if(abs(aratio(j-1)).gt.1.0e0_knd) go to 110
100       continue
110     continue
          do 120 n=2,j-1
          aratio(n)=(coef2(n)+coef3(n)/aratio(n-1))/coef1(n)
120       continue
        sumi=1.0e0_knd
        sumdi=1.0e0_knd
        termi=1.0e0_knd
        termdi=1.0e0_knd
        sumr=aratio(1)*qr(1)
        sumdr=aratio(1)*qdr(1)
        termr=sumr
        termdr=sumdr
          do 140 j=2,limq-2,2
          termi=-termi*aratio(j-1)*aratio(j)*qr(j-1)*qr(j)
          termdi=-termdi*aratio(j-1)*aratio(j)*qdr(j-1)*qdr(j)
          termr=-termr*aratio(j)*aratio(j+1)*qr(j)*qr(j+1)
          termdr=-termdr*aratio(j)*aratio(j+1)*qdr(j)*qdr(j+1)
          sumi=sumi+termi
          sumdi=sumdi+termdi
          sumr=sumr+termr
          sumdr=sumdr+termdr
          if(abs(termi/sumi)+abs(termdi/sumdi).gt.dec) go to 140
          if(abs(termr/sumr)+abs(termdr/sumdr).gt.dec) go to 140
          go to 150
140       continue
150     continue
        jleg1=min(j+1,limq-1)
        sumi=sumi*qm0
        sumdi=sumdi*qdm0
        sumr=sumr*qm0
        sumdr=sumdr*qdm0
        r2=r3*sumi+e3*sumr
        r2d=r3*sumdi+e3*sumdr
          if(4*((m+3)/4).eq.m+3) then
          r2d=-r2d
          end if
          if(2*(m/2).ne.m) then
          r2=-r2
          end if
        ir2e=int(log10(abs(r2)))
        r2c=r2*(ten**(-ir2e))
        if(abs(r2c).ge.1.0e0_knd) go to 160
        r2c=r2c*ten
        ir2e=ir2e-1
160     r2d=r2d+c*r1
        ir2de=int(log10(abs(r2d)))
        r2dc=r2d*(ten**(-ir2de))
        if(abs(r2dc).ge.1.0e0_knd) go to 170
        r2dc=r2dc*ten
        ir2de=ir2de-1
170     continue
if (debug) then
        write(40,180) jleg1,limq
180     format(8x,'r2leg1: series converged in ',i6, &
               ' terms; ',i6,' terms avail.')
end if
        return
        end subroutine
!
!
        subroutine r2neu0 (l,m,c,x,limneu,ndec,nex,maxd,maxlp,maxn, &
                           minacc,enr,sneuf,sneun,ineue,sneudf, &
                           sneudr,dfnorm,idfe,r1dc,ir1de,r2c,ir2e,r2dc, &
                           ir2de,jneu,jtest)
!
!  purpose:     To calculate the oblate radial function of the
!               second kind and its first derivative with respect
!               to x, using a series expansion of spherical
!               Neumann functions with argument c*sqrt(x*x+1);
!               (eta = 0)
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x      : x
!               limneu : maximum number of terms to be taken in the
!                        series summations for r2 and r2d
!               ndec   : number of decimal digits available in
!                        real arithmetic
!               nex    : maximum exponent if real arithmetic
!               maxd   : dimension of enr array
!               maxlp  : dimension of the sneun, sneudn, ineue, and
!                        ineude arrays
!               maxn   : dimension of sneuf and sneudf arrays
!               minacc : number of decimal digits of desired accuracy
!                        of the resulting radial functions
!               enr    : array of ratios of successive d coefficients
!               sneuf  : array of ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : array of characteristics for Neumann functions
!               ineue  : array of exponents for Neumann functions
!               sneudf : array of ratios of successive first derivatives
!                        of spherical Neumann functions of same parity
!               sneudr : array of ratios of first derivatives of Neumann
!                        functions to the corresponding functions
!               dfnorm : characteristic of Flammer normalization sum of
!                        d coefficients. equal to the reciprocal of
!                        the value of the d coefficient d(n = l - m) using
!                        this normalization for the angular functions
!               idfe   : exponent associated with dfnorm
!               r1dc   : charcteristic of corresponding first derivative
!                        of the radial function of the first kind
!               ir1de  : exponent of corresponding first derivative of
!                        the radial function of the first kind
!
!     output:   r2c    : characteristic of oblate radial function
!                        of the second kind
!               ir2e   : exponent of oblate radial function of the
!                        second kind
!               r2dc   : characteristic of derivative with respect
!                        to x of oblate radial function of the second
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        oblate radial function of the second kind
!               jneu   : index of term where best convergence is
!                        achieved for r2 or for r2d, whichever term is
!                        larger
!               jtest  : smaller of the number of digits of convergence
!                        of the forward sums for r2 and r2d
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,con,dconb,dconf,dconi,dec,dfnorm,dnew,dnewd,dold, &
                  doldd,rj1,rj2,rm,r1dc,r2,r2c,r2d,r2dc,r2dstore, &
                  r2dtemp,r2est,r2temp,sr2temp,sr2dtemp,sumcoef,sump, &
                  sumdp,ten,test,testd,testdm,teste,testeo,testm,tx,txd, &
                  x
        real(knd) enr(maxd),sneudr(maxlp),sneun(maxlp), &
                  sneuf(maxn),sneudf(maxn)
!
!  integer arrays
        integer ineue(maxlp)
!
        ten=10.0e0_knd
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        iscale=0
        rm=m
        dec=ten**(-ndec-1)
        dconf=ten**(-ndec)
        dconi=ten**(ndec+5)
        iexp=-ir1de-ineue(l+1)
          if(iexp.lt.nfac) then
          sumcoef=(ten**(iexp))/ &
                  (c*(x*x+1.0e0_knd)*r1dc*sneun(l+1))
          r2est=abs(sumcoef*dfnorm)*ten**(idfe)
          else
          r2est=ten**nfac
          end if
        dconb=r2est/dconi
        con=x/sqrt(x*x+1.0e0_knd)
        lm2=(l-m)/2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        mml=m+m-1+ix
        lim=limneu/2-ix
!
!  compute radial function of the second kind
!
!  backward series
        r2temp=1.0e0_knd
        sump=1.0e0_knd
        r2dtemp=1.0e0_knd
        sumdp=1.0e0_knd
        if(r2est.gt.dconi) go to 20
        if (lm2.lt.1) go to 20
        dold=1.0e0_knd
        doldd=1.0e0_knd
          do 10 j=lm2,1,-1
          jj=j+j+ix
          rj1=jj-ix
          rj2=jj+mml
          dnew=dold*rj1/(rj2*sneuf(jj+m)*enr(j))
          dnewd=doldd*rj1/(rj2*sneudf(jj+m)*enr(j))
          r2temp=r2temp+dnew
          if(dnew.gt.0.0e0_knd) sump=sump+dnew
          r2dtemp=r2dtemp+dnewd
          if(dnewd.gt.0.0e0_knd) sumdp=sumdp+dnewd
          if(abs(dnew/r2temp)+abs(dnewd/r2dtemp).lt.dconb) go to 20
          dold=dnew
          doldd=dnewd
10      continue
20      continue
!
!  forward series
        dold=1.0e0_knd
        doldd=1.0e0_knd
        testm=1.0e0_knd
        testdm=1.0e0_knd
        sr2temp=r2temp
        sr2dtemp=r2dtemp
        js=lim
        jds=lim
          do 70 j=lm2+1,lim
          jj=j+j+ix
          rj1=jj-ix
          rj2=jj+mml
          dnew=dold*enr(j)*sneuf(jj+m)*rj2/rj1
          dnewd=doldd*enr(j)*sneudf(jj+m)*rj2/rj1
          r2temp=r2temp+dnew
          if(dnew.gt.0.0e0_knd) sump=sump+dnew
          r2dtemp=r2dtemp+dnewd
          if(dnewd.gt.0.0e0_knd) sumdp=sumdp+dnewd
          test=abs(dnew/r2temp)
          testd=abs(dnewd/r2dtemp)
          if(test.lt.testm) go to 30
          go to 40
30        if(test.ne.0.0e0_knd) testm=test
          sr2temp=r2temp
          tx=sump
          js=j
40        continue
          if(testd.lt.testdm) go to 50
          go to 60
50        if(testd.ne.0.0e0_knd) testdm=testd
          sr2dtemp=r2dtemp
          txd=sumdp
          jds=j
60        continue
          if(test+testd.lt.dconf) go to 90
            if(abs(r2temp).gt.teste) then
            r2temp=r2temp*testeo
            r2dtemp=r2dtemp*testeo
            sr2temp=sr2temp*testeo
            sr2dtemp=sr2dtemp*testeo
            sump=sump*testeo
            sumdp=sumdp*testeo
            dnew=dnew*testeo
            dnewd=dnewd*testeo
            tx=tx*testeo
            txd=txd*testeo
            iscale=iscale+nfac
            end if
          dold=dnew
          doldd=dnewd
70        continue
        r2temp=sr2temp
        r2dtemp=sr2dtemp
90      continue
        jneu=max(js,jds)
        nterms=min(j,lim)
        jtestm=-int(log10(testm+dec))
        if(jtestm.lt.0) jtestm=0
        if(jtestm.gt.ndec) jtestm=ndec
        jtestdm=-int(log10(testdm+dec))
        if(jtestdm.lt.0) jtestdm=0
        if(jtestdm.gt.ndec) jtestdm=ndec
        jtest=min(jtestm,jtestdm)
        naccs1=0
        if(tx.ne.0.0e0_knd) naccs1=int(log10(abs(tx/r2temp)))
        if(naccs1.lt.0) naccs1=0
        if(naccs1.gt.ndec) naccs1=ndec
        naccn2=0
        if(txd.ne.0.0e0_knd) naccs2=int(log10(abs(txd/r2dtemp)))
        if(naccs2.lt.0) naccs2=0
        if(naccs2.gt.ndec) naccs2=ndec
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2=r2temp*sneun(l+1)/dfnorm
        if(ix.eq.1) r2=r2*con
        iterm=int(log10(abs(r2)))
        ir2e=ineue(l+1)-idfe+iterm+iscale
        r2c=r2*ten**(-iterm)
        if(abs(r2c).ge.1.0e0_knd) go to 100
        r2c=r2c*ten
        ir2e=ir2e-1
100    continue
        r2d=r2dtemp*sneun(l+1)*c*con*sneudr(l+1)/dfnorm
        r2dstore=r2d*con
        ndsub=0
        if(ix.eq.1) r2d=r2dstore+r2/(x*(x*x+1.0e0_knd))
        if(ix.eq.1) ndsub=-log10(abs(r2d/r2dstore))
        if(ndsub.lt.0) ndsub=0
        naccs2=naccs2+ndsub
if (debug) then
        write(40,110) nterms,lim,js,jtestm,naccs1,jds,jtestdm,naccs2
110     format(8x,'r2neu0 (eta=0) : numerator converged in ',i6, &
               ' terms; ',i6,' terms available.',/,15x,'best r2 at ',i6, &
               ' terms with convergence to ',i2,' digits;',i3,' digits', &
               ' sub. error.',/,15x,'best r2d at ',i6,' terms with', &
               ' convergence to ',i2,' digits; ',i2,' digits sub. ' &
               'error')

        if(ix.eq.1.and.ndsub.gt.0) write(40,120) ndsub
120     format(15x,'subtraction error in forming r2d = ',i2,' digits.')
end if
        iterm=int(log10(abs(r2d)))
        ir2de=ineue(l+1)-idfe+iterm+iscale
     r2dc=r2d*ten**(-iterm)
        if(abs(r2dc).ge.1.0e0_knd) go to 130
        r2dc=r2dc*ten
        ir2de=ir2de-1
130    continue
        return
        end subroutine
!
!
        subroutine r2eta (l,m,c,x,eta,nee,limeta,ndec,maxd,maxlp,maxn, &
                          maxp,minacc,wm,enr,sneuf,sneun,ineue, &
                          sneudf,sneudr,pdratt,pratb,pratt,pcoefn, &
                          ipcoefn,pdcoefn,ipdcoefn,r1c,ir1e,r1dc,ir1de, &
                          naccnmax,r2c,ir2e,r2dc,ir2de,nacceta,jeta, &
                          naccd)
!
!  purpose:     To calculate the oblate radial function of the
!               second kind and its first derivative with respect
!               to x, using an expansion of spherical Neumann
!               functions.
!
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : c
!               x       : x
!               eta     : value for eta used in calculation
!               nee     : index in the array of eta values in the main
!                         program that corresponds to the value of eta
!                         used in r2eta calculations
!               limeta  : maximum number of terms available in the sums
!                         for r2 and r2d
!               ndec    : number of decimal digits available in
!                         real arithmetic
!               maxd    : dimension of enr array
!               maxlp   : maximum  l value desired; dimension
!                         of the sneun, sneudr, and ineue arrays
!               maxn    : dimension of sneuf and sneudf arrays
!               maxp    : dimension of pdratt, pratb, and pratt arrays
!               minacc  : minimum number of accurate decimal digits
!                         that are requested
!               wm      : value of 1 - eta*eta computed in a way that
!                         avoids the subtraction error that would occur
!                         if it were computed directly when eta is near
!                         unity
!               enr     : array of ratios of successive d coefficients
!               sneuf   : array of ratios of successive spherical
!                         Neumann functions of the same parity
!               sneun   : array of characteristics for Neumann functions
!               ineue   : array of exponents corresponding to sneun
!               sneudf  : array of ratios of successive first
!                         derivatives of spherical Neumann functions of
!                         the same parity
!               sneudr  : array of ratios of first derivatives of the
!                         spherical Neumann functions to the
!                         corresponding functions
!               pdratt  : array of ratios of successive first
!                         derivatives of the associated Legendre
!                         functions of the first kind of the same parity
!                         (used in numerator series)
!               pratb   : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in denominator series)
!               pratt   : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in numerator series)
!               pcoefn  : characteristic of the ratio of the numerator
!                         and denominator associated Legendre functions
!                         of the first kind of order m and degree l
!               ipcoefn : exponent corresponding to pcoefn
!               pdcoefn : characteristic of the ratio of the first
!                         derivative of the associated Legendre function
!                         of the first kind in the numerator and the
!                         associated Legendre function of the first kind
!                         in the denominator, both of order m and
!                         degree l
!               ipdcoefn: exponent corresponding to pdcoefn
!               r1c     : characteristic of the radial function of the
!                         first kind (calculated in r1bes)
!               irie    : exponent corresponding to r1c
!               r1dc    : characteristic of the first derivative with
!                         respect to x of the radial function of the
!                         first kind (calculated in r1bes)
!               ir1de   : exponent corresponding to r1dc
!               naccnmax: maximum accuracy (in decimal digits) obtained
!                         for the current value of l from previous
!                         r2eta calculations
!
!     output:   r2c     : characteristic of the oblate radial function
!                         of the second kind
!               ir2e    : exponent of the oblate radial function of the
!                         second kind
!               r2dc    : characteristic of the first derivative with
!                         respect to x of the oblate radial function
!                         of the second kind
!               ir2de   : exponent corresponding to r2dc
!               nacceta : estimated number of accurate decimal digits in
!                         r2 and r2d, computed from the Wronskian
!               jeta    : maximum number of terms taken in the numerator
!                         sums for r2 and r2d
!               naccd:    estimated accuracy of the denominator series
!
!     input/
!     output:   naccnmax: maximum accuracy (in decimal digits) obtained
!                         for the current value of l from all previous
!                         r2eta calculations (input) and including the
!                         curent r2eta calculation (output)
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,dcon,dconi,dec,denom,dnew,dnewd,dnewd1,dnewd2, &
                  dold,doldd1,doldd2,eta,etas,pcoefn,pdcoefn,reld12, &
                  rm,rm2,r1c,r1dc,r2c,r2dc,r2dcoef1,r2dcoef2,r2dtemp, &
                  r2est,r2temp,sr2temp,sr2dtemp,r2test,sumcoef,sumdnp, &
                  sumdp,sumnp,ten,test,testd,testdm,testm,tx,txd,wm, &
                  wronc,wronca,wroncb,wroncm,wront,xet,xets,x
        real(knd) enr(maxd),sneudr(maxlp),sneun(maxlp),pratb(maxp), &
                  pratt(maxp),pdratt(maxp),sneuf(maxn),sneudf(maxn)
!
!  integer arrays
        integer ineue(maxlp)
!
        ten=10.0e0_knd
        dec=ten**(-ndec-2)
        rm=m
        dcon=ten**(-ndec)
        dconi=ten**(ndec+5)
        etas=eta*eta
        xet=sqrt(x*x+wm)
        xets=xet*xet
        sumcoef=(ten**(-ir1de-ineue(l+1)-ipcoefn))/ &
                (c*(x*x+1.0e0_knd)*r1dc*sneun(l+1)*pcoefn)
        r2dcoef1=eta*wm/(xets*xet)
        r2dcoef2=c*x/xet
        reld12=(r2dcoef2/r2dcoef1)*sneudr(l+1)*(pcoefn/ &
               pdcoefn)*ten**(ipcoefn-ipdcoefn)
        rm2=2*m
        lm2=(l-m)/2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        lim=limeta/2-ix
!
!  compute radial function of the second kind and its first derivative
!
!  backward series for denominator
        denom=1.0e0_knd
        sumdp=1.0e0_knd
        if (lm2.lt.1) go to 20
        dold=1.0e0_knd
          do 10 j=lm2,1,-1
          jj=j+j+ix
          dnew=dold/(pratb(jj+1)*enr(j))
          denom=denom+dnew
          if(dnew.gt.0.0e0_knd) sumdp=sumdp+dnew
          if(abs(dnew/denom).lt.dec) go to 20
          dold=dnew
10        continue
20      continue
!
!  forward series for denominator
        dold=1.0e0_knd
          do 30 j=lm2+1,lim
          jj=j+j+ix
          dnew=dold*enr(j)*pratb(jj+1)
          denom=denom+dnew
          if(dnew.gt.0.0e0_knd) sumdp=sumdp+dnew
          if(abs(dnew/denom).lt.dec) go to 40
          dold=dnew
30        continue
40      continue
        jden=j
        numsub=0
        if(sumdp.ne.0.0e0_knd) numsub=int(log10(abs(sumdp/denom)))
        if(numsub.lt.0) numsub=0
        if(numsub.gt.ndec) numsub=ndec
        naccd=ndec-max(2,int(log10(c)))-numsub
        if(naccd.lt.0) naccd=0
        r2est=abs(sumcoef*denom)
        r2test=r2est*dconi
!
!  backward series for numerator
        dold=1.0e0_knd
        doldd1=1.0e0_knd
        doldd2=reld12
        r2temp=1.0e0_knd
        sumnp=1.0e0_knd
        r2dtemp=doldd2
        sumdnp=0.0e0_knd
        if(doldd2.gt.0.0e0_knd) sumdnp=doldd2
        if(l.ne.0) r2dtemp=r2dtemp+1.0e0_knd
        if(l.ne.0) sumdnp=sumdnp+1.0e0_knd
        if(lm2.lt.1) go to 60
          do 50 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sneuf(jj+m)*pratt(jj+1)*enr(j))
          dnewd1=-doldd1/(sneuf(jj+m)*pdratt(jj+1)*enr(j))
          dnewd2=-doldd2/(sneudf(jj+m)*pratt(jj+1)*enr(j))
          r2temp=r2temp+dnew
          if(dnew.gt.0.0e0_knd) sumnp=sumnp+dnew
          dnewd=dnewd1+dnewd2
          r2dtemp=r2dtemp+dnewd
          if(dnewd.gt.0.0e0_knd) sumdnp=sumdnp+dnewd
          if(abs(dnew/r2temp)+abs(dnewd/r2dtemp).lt.dec) go to 60
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
50      continue
60      continue
        if(m.eq.0.and.jj.eq.2) r2dtemp=r2dtemp-dnewd1
!
!  forward series for numerator
        dold=1.0e0_knd
        doldd1=1.0e0_knd
        doldd2=reld12
        test=1.0e0_knd
        testd=1.0e0_knd
        testm=1.0e0_knd
        testdm=1.0e0_knd
        js=lm2
        jds=lm2
        sr2temp=r2temp
        sr2dtemp=r2dtemp
        tx=sumnp
        txd=sumdnp
        dnewsum=0.0e0_knd
        dnewdsum=0.0e0_knd
        doldd=reld12
        if(l.ne.0) doldd=reld12+1.0e0_knd
          do 110 j=lm2+1,lim-1
          jj=j+j+ix
          dnew=-dold*enr(j)*sneuf(jj+m)*pratt(jj+1)
          dnewd1=-doldd1*enr(j)*sneuf(jj+m)*pdratt(jj+1)
          dnewd2=-doldd2*enr(j)*sneudf(jj+m)*pratt(jj+1)
          r2temp=r2temp+dnew
          if(dnew.gt.0.0e0_knd) sumnp=sumnp+dnew
          if(dnew.ne.0.0e0_knd) test=abs(dnew/r2temp)
          if(test.ge.testm) go to 80
          testm=test
          sr2temp=r2temp
          js=j
          tx=sumnp
80        dnewd=dnewd1+dnewd2
          r2dtemp=r2dtemp+dnewd
          if(dnewd.gt.0.0e0_knd) sumdnp=sumdnp+dnewd
          if(dnewd.ne.0.0e0_knd) testd=abs(dnewd/r2dtemp)
          if(testd.ge.testdm) go to 100
          testdm=testd
          sr2dtemp=r2dtemp
          jds=j
          txd=sumdnp
100       if ((test+testd).lt.dcon) go to 130
          if(abs(r2temp)+abs(r2dtemp).gt.r2test) go to 120
          if((dnew.eq.0.0e0_knd).and.(dnewd.eq.0.0e0_knd)) go to 120
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
110       continue
120     r2temp=sr2temp
        r2dtemp=sr2dtemp
130     jmax=j
        jeta=js
        if(jds.gt.js) jeta=jds
        jtestm=-int(log10(abs(testm)))
        if(jtestm.lt.0) jtestm=0
        if(jtestm.gt.ndec) jtestm=ndec
        jtestdm=-int(log10(abs(testdm)))
        if(jtestdm.lt.0) jtestdm=0
        if(jtestdm.gt.ndec) jtestdm=ndec
        naccns1=0
        if(tx.ne.0.0e0_knd) naccns1=int(log10(abs(tx/r2temp)))
        if(naccns1.lt.0) naccns1=0
        if(naccns1.gt.ndec) naccns1=ndec
        naccns2=0
        if(txd.ne.0.0e0_knd) naccns2=int(log10(abs(txd/r2dtemp)))
        if(naccns2.lt.0) naccns2=0
        if(naccns2.gt.ndec) naccns2=ndec
        naccn1=min(jtestm-2,ndec-2-naccns1)
        naccn2=min(jtestdm-2,ndec-2-naccns2)
        naccn=min(naccn1,naccn2)
        if(naccn.lt.0) naccn=0
        naccnmaxp=naccnmax
        if(naccn.gt.naccnmax) naccnmax=naccn
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2c=r2temp*sneun(l+1)*pcoefn/denom
        iterm=0
        if(r2c.ne.0.0e0_knd) iterm=int(log10(abs(r2c)))
        ir2e=ineue(l+1)+ipcoefn+iterm
        r2c=r2c*ten**(-iterm)
        r2dc=(r2dcoef1*r2dtemp*sneun(l+1)*pdcoefn/denom)* &
             ten**(ineue(l+1)+ipdcoefn-ir2e)
        iterm=0
        if(r2dc.ne.0.0e0_knd) iterm=int(log10(abs(r2dc)))
        ir2de=ir2e+iterm
     r2dc=r2dc*ten**(-iterm)
if (debug) then
        write(40,140) jmax,jden,lim,js,jtestm,naccns1,jds,jtestdm, &
                     naccns2,naccn,naccd
140     format(8x,'r2eta: numerator, denominator converged in ', &
               i6,' ,',i6,' terms; ',i6,' terms available.',/, &
               15x,'best r2 at ',i6,' terms with convergence to ',i2, &
               ' digits; ',i2,' digits subtr. error.',/,15x, &
               'best r2d at ',i6,' terms with convergence to ',i2, &
               ' digits; ',i2,' digits subtr. error.',/,15x, &
               'estimated numerator and denominator accuracy is ',i2, &
               ' and ',i2,' digits.')
end if
        wronca=r1c*r2dc*(ten**(ir1e+ir2de))
        wroncb=r2c*r1dc*(ten**(ir2e+ir1de))
        wronc=wronca-wroncb
        wront=1.0e0_knd/(c*(x*x+1.0e0_knd))
        wroncm=max(abs(wronca),abs(wroncb))
        nacceta=-int(log10(abs((wronc-wront)/wront)+dec))
        if(nacceta.lt.0) nacceta=0
160     if(abs(r2c).ge.1.0e0_knd) go to 170
        r2c=r2c*ten
        ir2e=ir2e-1
170     continue
        if(abs(r2dc).ge.1.0e0_knd) go to 180
        r2dc=r2dc*ten
        ir2de=ir2de-1
180     continue
        if(naccn.gt.nacceta) naccnmax=max(naccnmaxp,nacceta)
        return
        end subroutine
!
!
        subroutine tridiag(d,e,n,np)
!
!  Purpose:     To calculate the eigenvalues of a symmetric tridiagonal
!               matrix by the implicit ql method and order them in
!               ascending value.
!
! Subroutine tridiag is a translation by Brian Gladman of the algol
! procedure imtql2, num. math. 12, 377-383(1968) by martin and
! wilkinson,as modified  in num. math. 15, 450(1970) by dubrulle,
! handbook for auto. comp.,vol.ii-linear algebra, 241-248(1971).
! The option and code for finding the eigenvectors has been removed.
!
!  parameters:
!
!     input:    d  : diagonal elements of the tridiagonal matrix
!               e  : the subdiagonal elements of the matrix
!               n  : order of the matrix
!               np : dimension of the vectors d and e
!
!     output:   d  : eigenvalues in ascending value.
!
        use param
!
        real(knd) b,c,f,g,p,r,s
        real(knd) d(np),e(np)
        integer ev_no
!
        do ev_no = 1, n
!  iterate to find the next eigenvalue
          do it = 0, 30
!  look for a small sub-diagonal element
            do i = ev_no, n - 1
            p = abs(d(i)) + abs(d(i + 1))
            if(p + abs(e(i)).eq.p) exit
            end do
!  end the iteration if we have an eigenvalue
          if(i.eq.ev_no) exit
!  form an implicit shift
          g = (d(ev_no + 1) - d(ev_no)) / (2.0e0_knd * e(ev_no))
            if(abs(g).ge.1.0e0_knd) then
            r = g * sqrt(1.0e0_knd + (1.e0_knd / g) ** 2)
            else
            r = 1.0e0_knd * sqrt(1.0e0_knd + g ** 2)
            end if
          g = d(i) - d(ev_no) + e(ev_no) / (g + sign(r, g))
          s = 1.0e0_knd
          c = 1.0e0_knd
          p = 0.0e0_knd
!
!  perform a plane rotation followed by Givens
!  rotations to restore tridiagonal form
            do j = i, ev_no + 1, -1
            f = s * e(j - 1)
            b = c * e(j - 1)
            if(abs(f).ge.abs(g)) then
            r = f * sqrt(1.0e0_knd + (g / f) ** 2)
            else
            r = g * sqrt(1.0e0_knd + (f / g) ** 2)
            end if
            e(j) = r
!  recover from underflow
            if (r.eq.0.0e0_knd) exit
            s = f / r
            c = g / r
            g = d(j) - p
            r = (d(j - 1) - g) * s + 2.0e0_knd * c * b
            p = s * r
            d(j) = g + p
            g = c * r - b
!
            end do
          d(j) = d(j) - p
          if(j.eq.ev_no) e(j) = g
          e(i) = 0.0e0_knd
          end do
        end do
!  order eigenvalues
        do i = 1,n-1
        k = i
        p = d(i)
          do j = i + 1, n
            if(d(j).lt.p) then
            k = j
            p = d(j)
            end if
          end do
          if(k.ne.i) then
          d(k) = d(i)
          d(i) = p
          end if
        end do
        return
        end subroutine
!
!
        subroutine conver (l,m,c,limd,blist,glist,ndec,maxd,ioprad, &
                           kindd,kindq,eigval,enr,ienr,itestm)
!
!  purpose:     To determine a converged eigenvalue using the
!               boukwamp method.
!  parameters:
!
!     input:    l     : l
!               m     : m
!               c     : c
!               limd  : number of enr values computed
!               blist : array of coefficients used in recursion relation
!               glist : array of coefficients used in recursion relation
!               ndec  : number of decimal digits available in real
!                       arithmetic
!               maxd  : dimension of enr,blist,glist arrays
!               ioprad: integer input equal to 0 if no radial functions
!                       are desired, equal to 1 if only radial functions
!                       of the first kind are desired, or equal to 2
!                       if radial functions of both the first and second
!                       kinds are desired
!               kindd : number of bytes for real data in double
!                       precision
!               kindq : number of bytes for real data in quadruple
!                       precision
!               eigval: estimated value of the eigenvalue
!
!     output:   eigval: converged eigenvalue
!               enr   : array of scaled ratios of successive d coefficients
!               ienr  : index n of last d coefficient ratio used
!                       in computing first term in the denominator
!                       of the eigenvalue correction
!               itestm: number of matching digits for the forward and
!                       backward recursion for d coefficient ratios
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,cora,corb,de,dec,dl,eigdec,eigstart,eigval,enrc,ten
        real(knd) blist(maxd),enr(maxd),enrf(maxd),glist(maxd)
!
        ten=10.0e0_knd
        nbp=int(2*c/3.1416)
        dec=ten**(-ndec-1)
        eigdec=dec*100.0e0_knd
        lmtest=max(67,(4*nbp)/3)
        lm2=(l-m)/2
        limdb=2*ienr+20
        if(l.eq.m.or.l.eq.m+1) limdb=2*ienr
        if(limdb.gt.limd) limdb=limd
!
!  begin Bouwkamp procedure
        icon=0
        eigstart=eigval
        jnde=1
        ix=l-m-2*lm2
        ifc=1
        lim2=limdb/2-ix
        iglim=lim2+1
        irio=lm2+1
        iw1=lm2+2
40      enr(1)=eigval-glist(1)
        if(lm2.lt.1) go to 60
!
!  evaluate the continued fraction
          do 50 i=1,lm2
          enr(i+1)=-blist(i)/enr(i)-glist(i+1)+eigval
50        continue
60      enr(lim2)=-blist(lim2)/(glist(iglim)-eigval)
        iw15=lim2-1
        ip=iw1+iw15
        if(iw15.lt.iw1) go to 80
!
!  evaluate the continued fraction
          do 70 i=iw1,iw15
          ipi=ip-i
          enr(ipi)=-blist(ipi)/(glist(ipi+1)-eigval+enr(ipi+1))
70        continue
80      enrc=-blist(irio)/(glist(irio+1)-eigval+enr(irio+1))
        de=enrc*enrc/blist(irio)
        corb=de
        if(lim2.lt.iw1) go to 100
!
!  compute first sum in the denominator of the correction
          do 90 i=iw1,lim2
          de=enr(i)*enr(i)/blist(i)*de
          corb=corb+de
          if(abs(de/corb).lt.dec) go to 100
90        continue
100     ienr=i
        if(ienr.lt.lim2-10) ienr=lim2-12
        de=1.0e0_knd
        cora=de
        if(lm2.eq.0) go to 120
!
!  compute second term in the denominator of the correction
          do 110 i=1,lm2
          de=blist(irio-i)/(enr(irio-i)*enr(irio-i))*de
          cora=cora+de
          if(abs(de/cora).lt.dec) go to 120
110       continue
!
!  compute the correction to the eigenvalue
120     dl=(enrc-enr(irio))/(cora+corb)
        eigval=dl+eigval
!
!  eigenvalue accurate enough?
        if(abs(dl/eigval).lt.eigdec) go to 130
        ifc=ifc+1
        if(ifc.le.20) go to 40
130     continue
        iflag=0
        int1=-int(log10(abs(dl/eigval)+dec))
        if(int1.lt.ndec-2.and.l-m+1.le.lmtest) &
                    iflag=1
        int2=-int(log10(abs((eigval-eigstart)/eigval)+dec))
        if(int2.lt.ndec-6.and.l-m+1.le.lmtest) iflag=1
        if(iflag.eq.1) eigval=eigstart
140     continue
if (debug) then
        if(knd.eq.kindd.and.ioprad.ne.0.and.iflag.eq.0) write(40,150) l, eigval,eigstart
        if(knd.eq.kindd.and.ioprad.eq.0.and.iflag.eq.0) write(50,150) l, eigval,eigstart
        if(knd.eq.kindq.and.ioprad.ne.0.and.iflag.eq.0) write(40,155) l, eigval,eigstart
        if(knd.eq.kindq.and.ioprad.eq.0.and.iflag.eq.0) write(50,155) l, eigval,eigstart
150     format(1x,'l =',i5,6x,'eigenvalue =',e24.15, '; estimate =',e24.15)
155     format(1x,'l =',i5,6x,'eigenvalue =',e39.31, '; estimate =',e39.31)
        if(knd.eq.kindd.and.ioprad.ne.0.and.iflag.eq.1) write(40,160) l, eigstart
        if(knd.eq.kindd.and.ioprad.eq.0.and.iflag.eq.1) write(50,160) l, eigstart
        if(knd.eq.kindq.and.ioprad.ne.0.and.iflag.eq.1) write(40,165) l, eigstart
        if(knd.eq.kindq.and.ioprad.eq.0.and.iflag.eq.1) write(50,165) l, eigstart
160     format(1x,'l =',i5,6x,'eigenvalue =',e24.15,' obtained' &
                     ' from tridiagonal matrix')
165     format(1x,'l =',i5,6x,'eigenvalue =',e39.31,' obtained' &
                     ' from tridiagonal matrix')
end if
!
!  calculate the d coefficient ratios (enr)
        lim2=limd/2-ix
        iesub=0
          if(lm2.eq.0.and.c.gt.50.0e0_knd) then
          iesub=-int(log10(abs((eigval-glist(1))/eigval)+dec))
          end if
        if(c.gt.50.0e0_knd.and.iesub.lt.4) go to 200
        enr(lim2)=-blist(lim2)/(glist(lim2+1)-eigval)
          do i=lim2-1,lm2+1,-1
          enr(i)=-blist(i)/(glist(i+1)-eigval+enr(i+1))
          end do
        if(lm2.eq.0) go to 190
        enr(1)=eigval-glist(1)
        if(lm2.eq.1) go to 190
          do n=2,lm2
          enr(n)=-blist(n-1)/enr(n-1)-glist(n)+eigval
          end do
190     continue
        itestm=ndec-1
        go to 260
200     enr(lim2)=-blist(lim2)/(glist(lim2+1)-eigval)
          do 210 i=lim2-1,1,-1
          enr(i)=-blist(i)/(glist(i+1)-eigval+enr(i+1))
210       continue
        enrf(1)=eigval-glist(1)
        nlim=1
        itestm=-int(log10(abs((enrf(1)-enr(1))/enrf(1))+ &
                 dec))
        if(itestm.lt.0) itestm=0
        if(itestm.ge.ndec) go to 230
          do 220 n=2,lim2
          enrf(n)=-blist(n-1)/enrf(n-1)-glist(n)+eigval
          itest=-int(log10(abs((enrf(n)-enr(n))/enrf(n))+ &
                 dec))
          if(itest.lt.0) itest=0
          if(itest.lt.itestm-4) go to 230
          if(itest.le.itestm) go to 220
          itestm=itest
          nlim=n
          if(itestm.ge.ndec) go to 230
220       continue
230     nlimp=2*(nlim-1)+ix
if (debug) then
        if(ioprad.ne.0) write(40,240) lim2,itestm,nlimp
        if(ioprad.eq.0) write(50,240) lim2,itestm,nlimp
240     format(15x,i6,' d coefs. Forward and backward recursion ', &
               'values of d(n+2)/d(n) match to ',i2,' digits at n = ', &
               i6)
end if
          do 250 n=1,nlim
          enr(n)=enrf(n)
250       continue
260     continue
        return
        end subroutine
!
!
        subroutine dnorm (l,m,c,ndec,nex,limd,maxd,enr,ioprad,iopang, &
                          dc01,idc01,dfnorm,idfe,dmlf,idmlfe,dmfnorm, &
                          idmfe,dmlmf,idmlmfe,dmsnorm,idmse,dmlms, &
                          idmlmse,jnorm,jsub,ksub)
!
!  purpose:     To compute d coefficient ratios from n values and to
!               calculate the normalization of the d coefficients.
!
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : c
!               ndec    : number of decimal digits available in
!                         real arithmetic
!               nex     : maximum exponent in real arithmetic
!               limd    : approximately twice the maximum number
!                         of terms available to be taken in the sum
!               maxd    : dimension of enr array
!               enr     : array of ratios of scaled d coefficients
!               ioprad  : set equal to zero if no radial functions are
!                         desired, set equal to two otherwise
!               iopang  : set equal to zero if no angular functions are
!                         desired, set equal to two otherwise
!
!     output:   enr     : array of ratios of d coefficients.
!                         enr(i) = ratio of the d coefficient with
!                         subscript 2*i+ix to the d coefficient with
!                         subscript 2*(i-1)+ix. Here ix =0 when l-m is
!                         even and ix=1 when l-m is odd.
!                         If the user needs the d coefficent ratios,
!                         they are available below right before
!                         statement 20.
!               dc01    : characteristic of ratio of first d
!                         coefficient (either d0 or d1, depending on
!                         whether l-m is even or odd) to the d
!                         coefficient of order equal to l-m
!               idc01   : exponent associated with dc01
!               dfnorm  : characteristic of Flammer normalization sum of
!                         d coefficients. equal to the reciprocal of
!                         the value of the d coefficient d(n = l - m)
!                         using this normalization for the angular
!                         functions
!               idfe    : exponent associated with dfnorm
!               dmlf    : characteristic of the d coefficient with index
!                         l-m in the Flammer normalization
!               idmlfe  : exponent associated with dmlf
!               dmfnorm : characteristic of Morse-Feshbach normalization
!                         sum of the d coefficients. equal to the
!                         reciprocal of the value of the d coefficient
!                         d(n = l - m) using this normalization for the
!                         angular functions
!               idmfe   : exponent associated with dmfnorm
!               dmlmf   : characteristic of the d coefficient with index
!                         l-m in the Morse-Feshbach normalization
!               idmlmfe : exponent associated with dmlmf
!               dmsnorm : characteristic of Meixner-Schafke normalization
!                         sum of the d coefficients. Equal to the
!                         reciprocal of the square of the d coefficient
!                         d(n = l - m) using this normalization for the
!                         angular functions
!               idmse   : exponent associated with dmsnorm
!               dmlms   : characteristic of the d coefficient with index
!                         l-m in the Meixner-Schafke normalization
!               idmlmse : exponent associated with dmlms
!               jnorm   : maximum index of enr required for convergence
!                         of dfnorm and dmfnorm
!               jsub    : number of decimal digits of subtraction error
!                         incurred in calculating dfnorm
!               ksub    : number of decimal digits of subtraction error
!                         incurred in calculating dmfnorm
!
        use param
!
!  real(knd) scalars and array
        real(knd) aj,aj2,arr,c,coef,csq,dec,dfnorm,dmlf,dmfnorm,dmlmf, &
                  dmsnorm,dmlms,dc01,ea,rm2,rm2m1,rm2m3,rm2p1,sgn,sump, &
                  ten,term,teste,testeo
        real(knd) enr(maxd)
!
        ten=10.0e0_knd
        rm2=m+m
        rm2m1=m+m-1
        rm2p1=m+m+1
        rm2m3=m+m-3
        dec=ten**(-ndec-1)
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        ir1tempe=0
        nbp=int(2*c/3.1416)
        csq=c*c
        lm2=(l-m)/2
        ix=l-m-2*lm2
        mml=m+m-1+ix
        lim2=limd/2-ix
        sgn=1.0e0_knd
          do 20 i=1,lim2
          arr=ix+i+i
          ea=arr+arr+rm2
            if(i.le.lm2) then
            if(enr(i).lt.(0.0e0_knd))  sgn=sgn*(-1.0e0_knd)
            end if
10        enr(i)=-(ea-1.0e0_knd)*(ea+1.0e0_knd)*enr(i)/((arr+rm2)* &
                   (arr+rm2-1.0e0_knd)*csq)
20        continue
!
!  compute the Morse-Feshbach normalizing factor
        dmfnorm=1.0e0_knd
        sump=1.0e0_knd
        term=1.0e0_knd
        jlow=l-m+2
        jterm=lm2
        iflag=0
        idmfe=0
          do 30 j=jlow,limd,2
          aj=j
          jterm=jterm+1
          term=term*(aj+rm2)*enr(jterm)*(aj+rm2-1.0e0_knd)/ &
               (aj*(aj-1.0e0_knd))
          if(term.gt.0.0e0_knd) sump=sump+term
          dmfnorm=dmfnorm+term
          if(abs(term/dmfnorm).lt.dec) go to 40
          if(abs(dmfnorm).lt.teste) go to 30
          dmfnorm=dmfnorm*testeo
          term=term*testeo
          sump=sump*testeo
          idmfe=idmfe+nfac
          iflag=1
30        continue
40      jlow=l-m
        jmf=jterm
        if(jlow.lt.2.or.iflag.eq.1) go to 60
        term=1.0e0_knd
        jterm=lm2
          do 50 j=jlow,2,-2
          aj=j
          term=term*aj*(aj-1.0e0_knd)/((aj+rm2 &
               -1.0e0_knd)*(aj+rm2)*enr(jterm))
          if(term.gt.0.0e0_knd) sump=sump+term
          jterm=jterm-1
          dmfnorm=dmfnorm+term
          if(abs(term/dmfnorm).lt.dec) go to 60
50        continue
60      continue
        if(dmfnorm.eq.0.0e0_knd) dmfnorm=dec
        ksub=int(log10(abs(sump/dmfnorm)+dec))
        if(ksub.lt.0) ksub=0
        if(ksub.gt.ndec) ksub=ndec
        iterm=0
        if(dmfnorm.ne.0.0e0_knd) iterm=int(log10(abs(dmfnorm)))
        idmfe=idmfe+iterm
        dmfnorm=dmfnorm*ten**(-iterm)
        dmlmf=1.0e0_knd/dmfnorm
        idmlmfe=-idmfe
if (debug) then
        if(ioprad.ne.0.and.c.gt.50.0e0_knd) write(40,70) jmf,ksub
        if(ioprad.eq.0.and.c.gt.50.0e0_knd) write(50,70) jmf,ksub
70      format(17x,'Morse-Feshbach norm. converged' &
               ' in ',i6,' terms with ',i2,' digit subtr. error.')
        if(ioprad.ne.0.and.c.le.50.0e0_knd) write(40,75) lim2,jmf,ksub
        if(ioprad.eq.0.and.c.le.50.0e0_knd) write(50,75) lim2,jmf,ksub
75      format(11x,i7,' d coefs. Morse-Feshbach norm. converged' &
               ' in ',i6,' terms with ',i2,' digit subtr. error.')
end if
!
!  compute the Flammer normalizing factor
        term=1.0e0_knd
        sump=term
        dfnorm=term
        idfe=0
        iflag=0
          do 80 j=lm2+1,lim2
          jj=j+j+ix
          term=-term*enr(j)*(jj+mml)/(jj-ix)
          dfnorm=dfnorm+term
          if(term.gt.0.0e0_knd) sump=sump+term
          if(abs(term/dfnorm).lt.dec) go to 90
          if(abs(dfnorm).lt.teste) go to 80
          dfnorm=dfnorm*testeo
          term=term*testeo
          sump=sump*testeo
          idfe=idfe+nfac
          iflag=1
80        continue
90      continue
        jf=min(j,lim2)
        if(lm2.lt.1.or.iflag.eq.1) go to 110
        term=1.0e0_knd
          do 100 j=lm2,1,-1
          jj=j+j+ix
          term=-term*(jj-ix)/((jj+mml)*enr(j))
          dfnorm=dfnorm+term
          if(term.gt.0.0e0_knd) sump=sump+term
          if(abs(term/dfnorm).lt.dec) go to 110
100       continue
110     continue
        jnorm=max(jf,jmf)
          if(dfnorm.eq.0.0e0_knd) then
          dfnorm=dec
          jsub=ndec
          else
          jsub=int(log10(abs(sump/dfnorm)+dec))
          if(jsub.lt.0) jsub=0
          if(jsub.gt.ndec) jsub=ndec
          end if
        iterm=0
        if(dfnorm.ne.0.0e0_knd) iterm=int(log10(abs(dfnorm)))
        idfe=idfe+iterm
        dfnorm=dfnorm*ten**(-iterm)
        dmlf=1.0e0_knd/dfnorm
        idmlfe=-idfe
if (debug) then
        if(ioprad.ne.0) write(40,120) jf,jsub
        if(ioprad.eq.0) write(50,120) jf,jsub
120     format(17x,'Flammer norm. converged in ', &
               i6,' terms with ',i2,' digit subtr. error.')
end if
!
!  compute the d0(c|ml) or d1(c|ml)
        idc01=0
       dc01=1.0e0_knd
        if(lm2.eq.0) go to 140
          do 130 kjl=1,lm2
          kkjl=lm2-kjl+1
          dc01=dc01/enr(kkjl)
            if(abs(dc01).gt.teste) then
            dc01=dc01*testeo
            idc01=idc01+nfac
            end if
            if(abs(dc01).lt.testeo) then
            dc01=dc01*teste
            idc01=idc01-nfac
            end if
130       continue
        iterm=int(log10(abs(dc01)))
        dc01=dc01*(ten**(-iterm))
        idc01=idc01+iterm
140     continue
!
!  compute the Meixner-Schafke normalizing factor
        if(iopang.eq.0) go to 200
        jflag=0
        idmse=0
        dmsnorm=1.0e0_knd
        coef=1.0e0_knd
        jlow=l-m+2
        jterm=lm2
          do 150 j=jlow,limd,2
          aj=j
          aj2=aj+aj
          jterm=jterm+1
          coef=coef*(aj+rm2)*enr(jterm)*(aj+rm2m1)*enr(jterm) &
               *(aj2+rm2m3)/(aj*(aj-1.0e0_knd)*(aj2+rm2p1))
          dmsnorm=dmsnorm+coef
          if(abs(coef/dmsnorm).lt.dec) go to 160
            if(abs(dmsnorm).gt.teste) then
            dmsnorm=dmsnorm*testeo
            coef=coef*testeo
            idmse=idmse+nfac
            jflag=1
            end if
150       continue
160     jlow=l-m
        jn=jterm
        jnorm=max(jnorm,jn)
        if(jlow.lt.2.or.jflag.eq.1) go to 180
        coef=1.0e0_knd
        jterm=lm2
        j=jlow
          do 170 jj=2,jlow,2
          aj=j
          aj2=aj+aj
          coef=coef*aj*(aj-1.0e0_knd)*(aj2+rm2p1)/((aj2+rm2m3)* &
                  enr(jterm)*enr(jterm)*(aj+rm2)*(aj+rm2m1))
          jterm=jterm-1
          j=j-2
          dmsnorm=dmsnorm+coef
          if(abs(coef/dmsnorm).lt.dec) go to 180
170       continue
180     iterm=int(log10(dmsnorm))
        dmsnorm=dmsnorm*ten**(-iterm)
        idmse=idmse+iterm
          if(2*(idmse/2).ne.idmse) then
          idmse=idmse-1
          dmsnorm=ten*dmsnorm
          end if
        dmlms=sgn/sqrt(dmsnorm)
        idmlmse=-idmse/2
if (debug) then
        write(50,190) jn,lim2
190     format(5x,' Meixner-Schafke normalization converged in ', &
               i6,' terms; ',i6,' terms available.')
end if
200     continue
        return
        end subroutine
!
!
        subroutine dalt (l,m,c,nex,limdr,maxdr,maxmp,ioppsum,eigval, &
                         enrneg,drhor,dneg,idneg,nsdneg,nsdrhor1)
!
!  purpose:     To calculate d ratios with negative subscripts
!               and d-rho ratios.
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : c
!               nex     : maximum exponent in real(knd) arithmetic
!               limdr   : number of ratios of successive d-rho
!                         coefficients calculated
!               maxdr   : dimension of drhor array
!               maxmp   : dimension of enrneg array
!               ioppsum : integer index = 0 if no d rho coefficients
!                         are calculated (psum not needed for r2leg)
!               eigval  : eigenvalue
!
!     output:   enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               dneg    : characteristic of the ratio of the d
!                         coefficient with index -2m+ix to the
!                         d coefficient with index ix, where
!                         ix = 0 if l-m is even and ix = 1 if
!                         l-m is odd
!               idneg   : exponent (base 10) of dneg
!               nsdneg  : subtraction error in calculating dneg
!               nsdrhor1: subtraction error in calculating drhor(1)
!                         from drhor(2)
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,dneg,eigval,r,rm,rn,t,ten,teste,testeo,uterm,vterm, &
                  wterm
        real(knd) ansdneg,asub,bsub
        real(knd) enrneg(maxmp),drhor(maxdr)
!
        ten=10.0e0_knd
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
!  if l-m is even, ix=0; if l-m is odd, ix=1
        ix=(l-m)-2*((l-m)/2)
        t=m+m-ix-ix
!
!  calculate ratios of d coefficients with negative subscripts
        rm=m
          if(m.eq.0) then
          dneg=0.0e0_knd
          idneg=0
          nsdneg=0
          go to 30
          end if
          do 10 i=1,m+1
          enrneg(i)=0.0e0_knd
10        continue
        n=2-2*m+ix
        rn=n
        r=n+m+m
        uterm=r*(r-1.0e0_knd)/((rn+r+1.0e0_knd)*(rn+r-1.0e0_knd))
        r=n+m-2
        vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-1.0e0_knd)/ &
              ((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-1.0e0_knd))- &
              (r*(r+1.0e0_knd)-eigval)/(c*c)
!
!       enrneg(k) = { d(-2m+2k-2)/d(-2m+2k), l-m even    }
!                   { d(-2m+1+2k-2)/d(-2m+1+2k), l-m odd }
!
!       calculations continue up to and including
!       enrneg(k=m) = { d(-2)/d(0), l-m even }
!                     { d(-1)/d(1), l-m odd  }
!
        enrneg(1)=-uterm/vterm
        dneg=enrneg(1)
        idneg=0
        ansdneg=0.0e0_knd
        nsdneg=0
          do 20 i=2,2*m,2
          ii=i-2*m+ix
          if (ii.ge.0) go to 25
          n=ii+2
          rn=n
          r=n+m+m
          uterm=r*(r-1.0e0_knd)/((rn+r+1.0e0_knd)*(rn+r-1.0e0_knd))
          r=n+m-2
          vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-1.0e0_knd)/ &
                ((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-1.0e0_knd))- &
                (r*(r+1.0e0_knd)-eigval)/(c*c)
          r=n-4
          wterm=(r+2.0e0_knd)*(r+1.0e0_knd)/((r+r+rm+rm+3.0e0_knd)* &
                (r+r+rm+rm+1.0e0_knd))
          j=i/2
          enrneg(j+1)=-uterm/(wterm*enrneg(j)+vterm)
          dneg=dneg*enrneg(j+1)
            if(vterm*enrneg(j)*wterm.ne.0.0e0_knd) &
               then
            asub=log10(abs(vterm/(vterm+ &
                             wterm*enrneg(j))))
            if(asub.gt.0.0e0_knd) ansdneg=ansdneg+asub
            bsub=log10(abs(vterm/(wterm*enrneg(j))))
            if(bsub.gt.0.0e0_knd) ansdneg=max(0.0e0_knd,ansdneg-bsub)
            if(int(ansdneg+0.5e0_knd).gt.nsdneg) nsdneg=ansdneg+ &
                                                       0.5e0_knd
            end if
            if(abs(dneg).gt.teste) then
            dneg=dneg*testeo
            idneg=idneg+nfac
            end if
            if(abs(dneg).lt.testeo) then
            dneg=dneg*teste
            idneg=idneg-nfac
            end if
          if(nsdneg.gt.5) nsdneg=nsdneg-1
20        continue
25      iterm=int(log10(abs(dneg)))
        dneg=dneg*(ten**(-iterm))
        idneg=idneg+iterm
!
!  calculate ratios of d rho coefficients
!
!       drhor(k-m) = { d(rho|2k)/d(rh0|2k-2), l-m even  }
!                    { d(rho|2k-1)/d(rho|2k-3), l-m odd }
!
30      if(ioppsum.eq.0) go to 60
          do 40 i=1,limdr
          drhor(i)=0.0e0_knd
40        continue
!
          do 50 i=2*limdr,6,-2
          n=4-i+ix-m-m
          ii=(i-2)/2
          rn=n
          r=n+m+m
          uterm=r*(r-1.0e0_knd)/((rn+r+1.0e0_knd)*(rn+r-1.0e0_knd))
          r=n+m-2
          vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-1.0e0_knd)/ &
                ((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-1.0e0_knd))- &
                (r*(r+1.0e0_knd)-eigval)/(c*c)
          r=n-4
          wterm=(r+2.0e0_knd)*(r+1.0e0_knd)/((r+r+rm+rm+3.0e0_knd)* &
                (r+r+rm+rm+1.0e0_knd))
          drhor(ii)=-uterm/(wterm*drhor(ii+1)+vterm)
50        continue
        n=-2*m+ix
        r=n+m-2
        vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-1.0e0_knd)/ &
              ((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-1.0e0_knd))- &
              (r*(r+1.0e0_knd)-eigval)/(c*c)
        r=n-4
        wterm=(r+2.0e0_knd)*(r+1.0e0_knd)/((r+r+rm+rm+3.0e0_knd)* &
              (r+r+rm+rm+1.0e0_knd))
!
!       the final value of ii is 1;
!       drhor(1) has a special value:
!       drhor(1) = { d(rho|2m+2)/d(-2m), l-m even  }
!                  { d(rho|2m+1)/d(-2m+1), l-m odd }
!
        drhor(1)=1.0e0_knd/((t-1.0e0_knd)*(t+1.0e0_knd)*(wterm*drhor(2)+ &
                    vterm))
        nsdrhor1=0
        if(vterm*wterm*drhor(2).ne.0.0e0_knd) nsdrhor1= &
                   int(log10(abs(vterm/(vterm+wterm*drhor(2)))))
        if(nsdrhor1.lt.0) nsdrhor1=0
        if(ix.eq.1) drhor(1)=-drhor(1)
60      continue
        return
        end subroutine
!
!
    subroutine gauss (ndec,n,x,w)
!
!  purpose:     To evaluate the coordinates and weighting factors
!               for an nth order Gaussian quadrature
!
!  parameters:
!
!     input:    n  : order of quadrature
!
!     output:   x  : coordinate values for quadrature
!               w  : weighting factors
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) delta,der,pi,ri,s,t,ten,test,u,v,z
        real(knd) x(n),w(n)
!
        ten=10.0e0_knd
        test=ten**(-ndec)
        imax=(n+1)/2
        pi=acos(-1.0_knd)
          do 40 i=1,imax
          ri=i
      z=cos(pi*(ri-0.25e0_knd)/(n+0.5e0_knd))
            do 20 j=1,30
            u=0.0e0_knd
        v=1.0e0_knd
          do 10 k=1,n
          t=u
              u=v
          v=((k+k-1)*z*u-(k-1)*t)/k
10             continue
            s=z*z-1.0e0_knd
        der=n*(z*v-u)/s
        delta=-v/der-0.5e0_knd*v*v*((n*n*s-n*z*z-n)*v+ &
                  2.0e0_knd*n*z*u)/(der*der*der*s*s)
            z=z+delta
        if(abs(delta/z).lt.test) go to 30
20          continue
30        continue
      x(i)=-z
      x(n+1-i)=z
      w(i)=2.0e0_knd/((1.0e0_knd-z*z)*der*der)
      w(n+1-i)=w(i)
40      continue
    return
    end subroutine
!
!
        subroutine pleg (m,lim,maxp,limcsav,iopd,ndec,nex,barg,narg, &
                         maxt,pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm,alpha, &
                         beta,gamma,coefa,coefb,coefc,coefd,coefe)
!
!  purpose:     To calculate ratios of successive associated Legendre
!               functions of the first kind for given arguments barg.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Legendre functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lim    : two less than the number of associated Legendre
!                        function ratios calculated for given arguments
!               maxp   : dimension of alpha, beta, gamma, coefa, coefb,
!                        coefc, coefd, and coefe arrays and second
!                        dimension of pr and pdr arrays
!               limcsav: integer equal to the number of coefficients in
!                        each of the arrays alpha, beta, gamma, coefa,
!                        coefb, coefc, coefd, and coefe arrays that
!                        have already been calculated in earlier calls
!                        to pleg for this value of m and will not be
!                        calculated again. [Note that the minimum
!                        array index for the coefficients is 3 and
!                        the maximum array index is limcsav+2]
!               iopd   : integer that is set = 0 if derivatives of
!                        Legendre functions (i.e., their ratios)
!                        are not required when iopang = 1 and the
!                        first derivatives of the angular functions
!                        are not requested.
!                        iopd is set = 1 when iopang = 2 and pleg is
!                        also being used to obtain ratios of first
!                        derivatives of Legendre functions for use in
!                        computing the first derivatives of the angular
!                        functions.
!                        iopd is set = 2 when pleg is being used to
!                        compute ratios of Legendre functions for use in
!                        the calculation of the denominator term used in
!                        calculating the radial functions of the second
!                        kind and their first derivatives in r2eta. Also
!                        used in r1eta in the same way when calculating
!                        radial functions of the first kind and their
!                        first derivatives.
!                        iopd is set = 3 when pleg is being used to
!                        compute ratios of both the Legendre functions
!                        and their first derivatives for use in the
!                        calculation of the numerator terms used
!                        in r2eta to calculate the radial functions of
!                        the second kind and their first deriatives. Also
!                        used in r1eta in the same way when calculating
!                        radial functions of the first kind and their
!                        first derivatives.
!                        iopd is set = 4 when pleg is being used to
!                        compute ratios of both the Legendre functions
!                        and their first derivatives for use is the
!                        calculation of r2 and r2d in r2leg.
!               ndec   : number of decimal digits in real(knd)
!                        arithmetic
!               nex    : maximum exponent in real(knd) arithmetic
!               barg   : array of narg values of eta for which Legendre
!                        functions are to be calculated
!               narg   : number of specified values of eta in barg array
!               maxt   : dimension of barg array
!
!     output:   pr     : array of ratios of successive first kind
!                        associated Legendre functions of the same
!                        parity
!               pdr    : array of ratios of successive derivatives of
!                        first kind associated Legendre functions of
!                        the same parity
!               pdnorm : array of characteristics of the first
!                        derivatives of associated Legendre functions
!                        of the first kind of order m and degree m
!               ipdnorm: array of exponents of the first derivatives
!                        of associated Legendre functions of the first
!                        kind of order m and degree m
!               pnorm  : array of characteristics of the associated
!                        Legendre functions of the first kind of order m
!                        and degree m
!               ipnorm : array of exponents of the associated Legendre
!                        functions of the first kind of order m and
!                        degree m
!
!     input/output:
!               alpha  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               beta   : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               gamma  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               coefa  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefb  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefc  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefd  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefe  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) adec,ajterm,am2p1,anden1,anden2,an2tnp1,bargs,coef, &
                  den,rm,rm2,temp1,temp2,temp3,ten,term,teste,testeo, &
                  ta,tb,tc,t1,t2,t3
        real(knd) alpha(maxp),barg(maxt),beta(maxp),coefa(maxp), &
                  coefb(maxp),coefc(maxp),coefd(maxp),coefe(maxp), &
                  gamma(maxp),pdnorm(maxt),pdr(maxt,maxp),pdr1(maxp), &
                  pr(maxt,maxp),pnorm(maxt)
!
!  integer array
        integer ipdnorm(maxt),ipnorm(maxt)
!
        ten=10.0e0_knd
        adec=ten**(-ndec+2)
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        rm=m
        rm2=m*m
        am2p1=m+m+1
        m2=m+m
        m2m1=m2-1
        mm1=m-1
        mm2=m-2
        msqp1=2*m*m+1
        coef=1.0e0_knd
        if(iopd.eq.4) coef=-1.0e0_knd
!
!  calculate the coefficients alpha(j), beta(j), and gamma(j) for
!  the three term recursion relating the Legendre function ratios
!
!              m                m
!   pr(k,j) = p    (barg(k)) / p  (barg(k))
!              m+j-1            m+j-3
!
!  and calculate the coefficients coefa(j), coefb(j), coefc(j),
!  coefd(j), and coefe(j) in the expression used to calculate
!  ratios of Legendre function derivatives
!
!               m                 m
!   pdr(k,j) = p'    (barg(k)) / p'  (barg(k))
!               m+j-1             m+j-3
!
        if(limcsav.ge.lim) go to 30
          do 10 j=limcsav+3,lim+2
            n=m+j-3
            n2=n+n
            n2p3=n2+3
            n2p1=n2+1
            n2m1=n2-1
            nmmm2=n-mm2
            nmmm1=n-mm1
            npmm1=n+mm1
            npm=n+m
            npmm1=n+mm1
          npmp1=n+m+1
          npmp2=npmp1+1
            an2tnp1=2*n*real(n+1,knd)
            anden1=nmmm2*real(nmmm1,knd)
            anden2=n2m1*anden1
            alpha(j)=real(n2p3,knd)*n2p1/anden1
            beta(j)=real(n2p1,knd)*(real(msqp1,knd)-an2tnp1)/anden2
            gamma(j)=-real(n2p3,knd)*real(npm,knd)*real(npmm1,knd)/anden2
          coefa(j)=-real(npmp2,knd)/real(nmmm1,knd)
          coefb(j)=real(n2p3,knd)*(n+2)/anden1
          coefc(j)=-real(npmp1,knd)*npmp2/anden1
          coefd(j)=real(npmp1,knd)/real(nmmm2,knd)
          coefe(j)=-real(n+1,knd)*n2p3/anden1
10        continue
        gamma(3)=0.0e0_knd
        gamma(4)=0.0e0_knd
        term=1.0e0_knd
        iterm=0
        if(m.lt.2) go to 30
          do jm=2,m
            term=(jm+jm-1)*term
              if(term.gt.teste) then
              term=term*testeo
              iterm=iterm+nfac
            end if
          end do
        jterm=int(log10(term))
        term=term*(ten**(-jterm))
        iterm=iterm+jterm
30      continue
!
!   calculate the ratios of successive Legendre functions of the same
!   parity using the three term recursion relationship
!
!   pr(k,j) = coef*alpha(j)*barg(k)*barg(k) + beta(j) +
!             gamma(j)/pr(k,j-2)
!
!   where coef = -1 when computing functions for r2leg, = +1 otherwise
!
          do 140 k=1,narg
            pnorm(k)=term
          ipnorm(k)=iterm
          pdnorm(k)=term
          ipdnorm(k)=iterm
!
!   define the first two ratios equal to unity and (2m+1)*barg(k)
          pr(k,1)=1.0e0_knd
          pr(k,2)=am2p1*barg(k)
          jdelta=1
          if(abs(barg(k)).lt.adec) jdelta=2
          bargs=barg(k)*barg(k)
            do 40 j=3,lim+2,jdelta
            pr(k,j)=coef*alpha(j)*bargs+beta(j)+(gamma(j)/pr(k,j-2))
40          continue
!
!   for abs(eta) > 0.1, calculate the corresponding ratios of first
!   derviatives using the relationship (except for eta = unity, where
!   a special expression is used)
!
!              (coefa(j)+coef*coefb(j)*barg(k)*barg(k))*pr(k,j)+coefc(j)
!   pdr(k,j) = ----------------------------------------------------
!                  pr(k,j)+coef*coefd(j)+coef*coefe(j)*barg(k)*barg(k)
!
!   where coef = -1 when computing functions for r2leg, = +1 otherwise
!
!   for abs(eta) <= 0.1, calculate the ratios using a recursion relation
!   involving successive ratios (except for eta = 0, where a special
!   expression is used)
!
          if(iopd.eq.0.or.iopd.eq.2) go to 120
          pdr(k,1)=1.0e0_knd
          pdr(k,2)=1.0e0_knd
          if(abs(barg(k)).ge.adec) go to 50
          pdr(k,2)=am2p1
            do j=4,lim+2,2
            pdr(k,j)=-real(m2m1+j,knd)/(j-2)
            end do
          go to 140
50        if(abs(abs(barg(k))-coef).ge.adec) go to 70
          if(m.eq.0) go to 60
          if(m.ne.2) go to 130
          pdr(k,1)=-2.0e0_knd*barg(k)
          go to 80
60        temp1=1.0e0_knd
          temp2=3.0e0_knd
          pdr(k,2)=1.0e0_knd
          pdr(k,3)=3.0e0_knd*barg(k)
            do j=4,lim+2
            temp3=temp2+real(j-1,knd)
            pdr(k,j)=temp3/temp1
            temp1=temp2
            temp2=temp3
            end do
          go to 140
70        if(m.ne.0) go to 80
          pdr(k,1)=1.0e0_knd
          pdr(k,2)=1.0e0_knd
          pdr(k,3)=3.0e0_knd*barg(k)
          jlow=4
          go to 90
80        pdr(k,2)=am2p1*(+coef*(rm+1.0e0_knd)*bargs-1.0e0_knd)/ &
                   (rm*barg(k))
          jlow=3
90        continue
           if(abs(barg(k)).le.0.1e0_knd) go to 110
            do 100 j=jlow,lim+2
            den=(pr(k,j)+coefd(j)+coef*coefe(j)*bargs)
            if(den.eq.0.0e0_knd) den=1.0e-50_knd
            pdr(k,j)=((coefa(j)+coef*coefb(j)*bargs)*pr(k,j)+ &
                     coefc(j))/den
           if(iopd.eq.3.and.abs(pdr(k,j)).lt.1.0e-5_knd) go to 110
100        continue
         go to 120
110      continue
         if(m.ne.0) pdr1(1)=pdr(k,2)
         if(m.eq.0) pdr1(2)=pdr(k,3)
           do j=jlow-1,lim+1
           n=j+m-1
           t3=coef
           if(coef.eq.-1.0e0_knd.and.2*(j/2).eq.j) t3=1.0e0_knd
           t1=coef*bargs-1.0e0_knd
           t2=(n*n*t1+rm2)
           ta=j*t2
           tb=(n+n+1)*t3*barg(k)*(t2+n*t1)
           tc=-(n+m)*(t2+(n+n+1)*t1)
           pdr1(j)=(tb+tc/pdr1(j-1))/ta
           end do
           do j=jlow,lim+2
           pdr(k,j)=pdr1(j-2)*pdr1(j-1)
           end do
120       continue
          if(m.eq.0.or.iopd.eq.2.or.iopd.eq.3.or.iopd.eq.4) go to 140
          if(abs(abs(barg(k))-1.0e0_knd).lt.adec) go to 130
          ajterm=rm*log10(1.0e0_knd-bargs)/2.0e0_knd
          jterm=int(ajterm)
          ipnorm(k)=ipnorm(k)+jterm
          pnorm(k)=pnorm(k)*(ten**(ajterm-real(jterm,knd)))
          if(iopd.eq.0) go to 140
          ajterm=log10(rm*abs(barg(k)))+(rm-2.0e0_knd)* &
                 log10(1.0e0_knd-bargs)/2.0e0_knd
          jterm=int(ajterm)
          ipdnorm(k)=ipdnorm(k)+jterm
          pdnorm(k)=-pdnorm(k)*(ten**(ajterm-real(jterm,knd)))
          if(barg(k).lt.0.0e0_knd) pdnorm(k)=-pdnorm(k)
          go to 140
130       pnorm(k)=0.0e0_knd
          ipnorm(k)=0
          if(m.ne.2) pdnorm(k)=0.0e0_knd
          if(m.ne.2) ipdnorm(k)=0
140       continue
        return
        end subroutine
!
!
        subroutine qleg (m,lnum,limq,maxmp,maxq,x,ndec,nex,qdr,qdqr, &
                         qdm1,iqdm1,qdl,iqdl,qr,qm1,iqm1,ql,iql,termpq, &
                         itermpq,qr1,qdr1,qm0,qdm0)
!
!  purpose:     To calculate ratios of successive associated Legendre
!               functions of the second kind for given c,x, and m.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Legendre functions of the second
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lnum   : number of l values desired (=lmax+1);
!                        also equal to the dimension of the arrays
!                        ql, iql, qdl, and iqdl
!               limq   : the number of associated Legendre function
!                        ratios minus 2m calculated for given m,lnum,
!                        c,ndec,and x
!               maxmp  : dimension of qdqr array
!               maxq   : dimension of qr and qdr arrays
!               x      : radial coordinate x
!               ndec   : number of decimal digits in real(knd)
!                        arithmetic
!               nex    : maximum exponend in real(knd) arithmetic
!
!     output:   qdr    : array of ratios of derivatives of successive
!                        associated Legendre functions of the second
!                        kind
!               qdqr   : array of ratios of derivatives of associated
!                        Legendre functions of the second kind to the
!                        corresponding Legendre function for degrees
!                        from -m to m-1
!               qdm1   : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree m-1, scaled by
!                                       -m/2
!                        (2m-1)!!(x*x+1)
!               iqdm1  : exponent corresponding to qdm1
!               qdl    : array of characteristics of the first
!                        derivatives of the associated Legendre
!                        functions of the second kind with order m
!                        and degrees from m to m+lnum-1, scaled by
!                                       -m/2
!                        (2m-1)!!(x*x+1)
!               iqdl   : array of exponents corresponding to qdl
!               qr     : array of ratios of successive associated
!                        Legendre functions of the second kind
!               qm1    : characteristic of the associated Legendre
!                        function of the second kind with order m
!                                                                 -m/2
!                        and degree m-1, scaled by (2m-1)!!(x*x+1)
!               iqm1   : exponent corresponding to qm1
!               ql     : array of characteristics of the associated
!                        Legendre function of the second kind with
!                        order m and degrees from m to m+lnum-1
!                                                 -m/2
!                        scaled by (2m-1)!!(x*x+1)
!               iql    : array of exponents corresponding to ql
!               termpq : characteristic of the relative size of the
!                        maximum terms in the positive degree q series
!                        and the p series used to calculate r2 and r2d
!                        in subroutine r2leg
!               itermpq: exponent corresponding to termpq
!               qr1    : array of ratios of successive associated
!                        Legendre functions of the second kind
!                        reindexed for use in the Baber and Hasse
!                        expansion
!               qdr1   : array of ratios of derivatives of successive
!                        associated Legendre functions of the second
!                        kind reindexed for use in the Baber and Hasse
!                        expansion
!               qm0    : characteristic of the associated Legendre
!                        function of the second kind with order m
!                                                                 -m/2
!                        and degree 0, scaled by (2m-1)!!(x*x+1)
!               qdm0   : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree 0, scaled by
!                                       -m/2
!                        (2m-1)!!(x*x+1)
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) ajm,dec,qdm1,qlow,qlow0,qlow1,qmid,qmid0,qmid1,qm1, &
                  qupp,qupp0,qupp1,q00,q0m,q1m,q11,qmmm1,rin,rin1, &
                  rin2,rin3,rin4,rm,rmsq,ten,term,termpq,teste,testeo, &
                  tjm,tm,tmr,x,xang,xfac,xc,x1d,xsqr
        real(knd) qdl(lnum),qdr(maxq),ql(lnum),qr(maxq),qdqr(maxmp)
        real(knd) qr1(maxq),qdr1(maxq),qdm0,qm0
!
!  integer arrays
        integer iqdl(lnum),iql(lnum)
!
!  Note that the factor i involved in these functions is suppressed
!  and taken into account in calculations in the subroutine r2leg
!
        ten=10.0e0_knd
        dec=ten**(-ndec)
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        rm=m
        rmsq=rm*rm
        tm=rm+rm
        tmr=tm/(tm+1.0e0_knd)
        x1d=x*x+1.0e0_knd
        xsqr=sqrt(x1d)
        xang=asin(x/xsqr)
        q00=-atan(1.0e0_knd/x)
        if(m.eq.0) qm0=q00
        iflag=0
        if(x.lt.0.001e0_knd) iflag=1
!
!                m
!  For x < 0.1: Q   is calculated by forward recursion in m starting
!                m
!                   0       1
!  with values for Q   and Q  obtained from closed form expressions,
!                   0       1
!                           -m/2
!  scaled by (2m-1)!!(x*x+1).    For x >= 0.1: It is calculated using
!                                                n     n-1      n-2
!  forward recursion on the expression relating Q  to Q    and Q  .
!                                                m     m        m
!                          0
!  The starting value for Q  is obtained by reverse recursion on
!                          m
!                           0     0        0
!  the expression relating Q  to Q    and Q    and normalizing the
!                           n     n+1      n+2
!                                    0                  1
!  result using the known value for Q  . Similarly for Q  .
!                                    0                  m
!
        iqlterm=0
        if(m.ne.0) go to 10
        ql(1)=q00
        qr(1)=x+1.0e0_knd/q00
        go to 30
10      q11=-x1d*q00-x
        if(m.ne.1) go to 20
        ql(1)=q11
        qr(3)=3.0e0_knd*x-2.0e0_knd/q11
        go to 30
20      continue
        if(x.ge.0.1e0_knd) go to 25
        qlow=q00
        qupp=q11
        if(iflag.eq.1) qmmm1=3.0e0_knd*x*q11-2.0e0_knd
          do jm=1,m-1
          ajm=jm
          tjm=real(jm+jm,knd)/(jm+jm+1)
          ql(1)=(-x1d-tjm)*qupp-tjm*x1d*qlow
            if(iflag.eq.1) then
            qmmm1=-((ajm+ajm+2)/(ajm+ajm+1))*x1d*qmmm1+x*ql(1)
            end if
          qlow=qupp
          qupp=ql(1)
          end do
        if(iflag.eq.1) qr(m+m+1)=qmmm1/ql(1)
        go to 30
25      limm=-int(ndec/(log10(xsqr-x)))
        qupp0=0.0e0_knd
        qmid0=1.0e0_knd
        qupp1=0.0e0_knd
        qmid1=1.0e0_knd
          do jn=limm+m,m,-1
          rin1=jn+1
          rin2=jn+2
          rin3=jn+jn+3
          if(2*(jn/2).ne.jn) rin3=-rin3
          qlow0=(rin3*x*qmid0-rin2*qupp0)/rin1
          qlow1=(rin3*x*qmid1-rin1*qupp1)/rin2
          qupp0=qmid0
          qupp1=qmid1
          qmid0=qlow0
          qmid1=qlow1
          end do
        q0m=qlow0
        q1m=qlow1
        iqlterm=0
          do jn=m-1,1,-1
          rin1=jn+1
          rin2=jn+2
          rin3=jn+jn+3
          if(2*(jn/2).ne.jn) rin3=-rin3
          qlow0=(rin3*x*qmid0-rin2*qupp0)/rin1
          qlow1=(rin3*x*qmid1-rin1*qupp1)/rin2
            if(abs(qlow0).gt.teste) then
            qlow0=qlow0*testeo
            qlow1=qlow1*testeo
            qmid0=qmid0*testeo
            qmid1=qmid1*testeo
            iqlterm=iqlterm+nfac
            end if
          qupp0=qmid0
          qupp1=qmid1
          qmid0=qlow0
          qmid1=qlow1
          end do
        rin3=qlow0
        qlow0=3.0e0_knd*x*qmid0-2.0e0_knd*qupp0
        q1m=(q11/qlow1)*q1m
        q0m=(q00/qlow0)*q0m
        qlow=q0m
        qmid=q1m
          do j=1,m-1
          rin1=j+j
          rin2=j+j+1
          rin3=(m+j)*(m-j+1)
          rin4=(j+j+1)*(j+j-1)
          qupp=-(rin1*x/rin2)*qmid+(rin3*x1d/rin4)*qlow
          qlow=qmid
          qmid=qupp
          end do
          if(2*(m/4).ne.(m/2)) qupp=-qupp
          ql(1)=qupp
30      iql(1)=int(log10(abs(ql(1))))
        ql(1)=ql(1)*(ten**(-iql(1)))
        iql(1)=iql(1)-iqlterm
!
!                        m    m
!  the ratios qr(k+m) = Q  / Q    , k=m+limq to k=m+1 are calculated
!                        k    k-1
!
!  for x > 0.001 using backward recursion from qr(mxqr+2m) =
!
!   m        m
!  Q      / Q       = x - sqrt(x*x+1), where the last quantity
!   mxqr+m   mxqr+m-1
!
!  is the asymptotic limit of the ratio as mxqr approaches infinity.
!  Otherwise the ratios are calculated using forward recursion from
!  the ratio qr(m+m+1)
!
        if(iflag.eq.1) go to 40
        mxqr=limq-int(ndec/(log10(xsqr-x)))
        if(mxqr.lt.2*limq) mxqr=2*limq
        qupp=x-xsqr
          do jn=mxqr+m,limq+m+1,-1
          rin=jn
          qlow=-(rin+rm-1.0e0_knd)/(x*(rin+rin-1.0e0_knd) &
               -(rin-rm)*qupp)
          qupp=qlow
          end do
        qr(limq+m+m)=qupp
          do jn=limq+m,m+2,-1
          rin=jn
          qr(jn+m-1)=-(rin+rm-1.0e0_knd)/(x*(rin+rin-1.0e0_knd) &
                     -(rin-rm)*qr(jn+m))
          end do
        go to 50
40      continue
          do jn=m+2,limq+m
          rin=jn
          qr(jn+m)=((rin+rin-1.0e0_knd)*x+(rin+rm-1.0e0_knd)/ &
                   qr(jn+m-1))/(rin-rm)
          end do
50      continue
!
!                              m     m
!  calculate ratios qr(k+m) = Q   / Q   ,k=m-1 to k=-m+1 using
!                              k-1   k
!
!                                       m      m
!  backward recursion from qr(m+m-1) = Q    / Q     = x
!                                       m-2    m-1
!
        if(m.eq.0) go to 120
        qr(m+m-1)=x
        if(m.eq.1) go to 70
          do 60 jn=m-1,2-m,-1
          rin=jn
          qr(jn+m-1)=(x*(rin+rin-1.0e0_knd) &
                     +((rin-rm)/qr(jn+m)))/(rin+rm-1.0e0_knd)
          if(qr(jn+m-1).eq.0.0e0_knd) qr(jn+m-1)=dec
60        continue
70      continue
!
!                  m
!  calculation of Q    , m > 0 by forward division of qr ratios
!                  m-1
!
!                 m
!  starting with Q  calculated from its closed form expression.
!                 0
!
          if(2*(m/2).eq.m) then
          qm1=sin(rm*xang)
          qm0=qm1
          qdm0=rm*cos(rm*xang)/x1d
          else
          qm1=cos(rm*xang)
          qm0=qm1
          qdm0=-rm*sin(rm*xang)/x1d
          end if
        xfac=0.5e0_knd*rm*log10(x1d)
        ixfac=int(xfac)
        qm1=qm1*(ten**(xfac-ixfac))
        term=1.0e0_knd
        iterm=0
        if(m.lt.2) go to 90
          do jm=2,m
          ajm=jm
          term=term*(ajm-1.0e0_knd)/(ajm+ajm-1.0e0_knd)
            if(term.lt.testeo) then
            term=term*teste
            iterm=iterm-nfac
            end if
          end do
90      qm1=qm1*term
        jterm=int(log10(abs(qm1)))
        qm1=qm1*ten**(-jterm)
        iqm1=ixfac+iterm+jterm
          if(2*(m/2).eq.m) then
          if(2*(m/4).ne.m/2) qm1=-qm1
          else
          if(2*((m-1)/4).ne.(m-1)/2) qm1=-qm1
          end if
        if(m.lt.2) go to 110
          do jm=1,m-1
          qm1=qm1/qr(jm+m)
            if(abs(qm1).gt.teste) then
            qm1=qm1*testeo
            iqm1=iqm1+nfac
            end if
          end do
110     continue
        iterm=int(log10(abs(qm1)))
        qm1=qm1*ten**(-iterm)
        iqm1=iqm1+iterm
120     continue
!
!  calculation of ratios of the first derivatives of q with respect
!  to x for degrees >= m using the relationship:
!
!                  m    m      [kx]qr(k+m)+(k+m)
!     qdr(k+m) = Q'  / Q'   =  ----------------- , k=m+1 to k=m+lim
!                  k    k-1    [(k-m)]qr(k+m)-kx
!
!
          do jm=m+1,m+limq
          ajm=jm
          qdr(jm+m)=(ajm*x*qr(jm+m)+(ajm+rm))/((ajm-rm)*qr(jm+m)- &
                     ajm*x)
          end do
!
!                   m                       m
!  calculation of q'    from the value for q
!                   m-1                     m-1
!                       m
!  and calculation of q'   , k=0 to lnum-1, from the value
!                       m+k
!        m
!  for q'
!        m
!
        if(m.gt.0) go to 140
        qdl(1)=1.0e0_knd/x1d
        iqdl(1)=0
        qdm0=qdl(1)
        go to 150
140     qdm1=-rm*x*qm1/x1d
        iterm=int(log10(abs(qdm1)))
        qdm1=qdm1*(ten**(-iterm))
        iqdm1=iqm1+iterm
        qdl(1)=rm*(x*ql(1)-2.0e0_knd*qm1*(ten**(iqm1-iql(1))))/x1d
        iqdl(1)=iql(1)
150     continue
        m2m1=m+m-1
          do jl=2,lnum
          ql(jl)=ql(jl-1)*qr(m2m1+jl)
          iql(jl)=iql(jl-1)
            if(abs(ql(jl)).gt.teste) then
            ql(jl)=ql(jl)*testeo
            iql(jl)=iql(jl)+nfac
            end if
            if(abs(ql(jl)).lt.testeo) then
            ql(jl)=ql(jl)*teste
            iql(jl)=iql(jl)-nfac
            end if
          qdl(jl)=qdl(jl-1)*qdr(m2m1+jl)
          iqdl(jl)=iqdl(jl-1)
            if(abs(qdl(jl)).gt.teste) then
            qdl(jl)=qdl(jl)*testeo
            iqdl(jl)=iqdl(jl)+nfac
            end if
            if(abs(qdl(jl)).lt.testeo) then
            qdl(jl)=qdl(jl)*teste
            iqdl(jl)=iqdl(jl)-nfac
            end if
          end do
          do jl=1,lnum
          iterm=int(log10(abs(ql(jl))))
          ql(jl)=ql(jl)*ten**(-iterm)
          iql(jl)=iql(jl)+iterm
          iterm=int(log10(abs(qdl(jl))))
          qdl(jl)=qdl(jl)*ten**(-iterm)
          iqdl(jl)=iqdl(jl)+iterm
          end do
        termpq=rm*log10(xsqr)
        itermpq=int(termpq)
        termpq=ten**(termpq-itermpq)
!
!  Calculation of ratios of the first derivatives of q with respect
!  to x for degrees from 0 to m-1 via backward recursion using a
!  relation developed from the traditional recursion relations.
!  Here
!                 m      m
!     qdr(k+m) = Q'   / Q'  , k=m-1 to k=1
!                 k-1    k
!
!  The backward recursion is started with the value for qdr(m+m-1) =
!  (x*x*(m-1)-1)/(m*x).
!
        if(m.eq.0) go to 180
        qdr(m+m-1)=(x*x*(rm-1.0e0_knd)-1.0e0_knd)/(rm*x)
        if(qdr(m+m-1).eq.0.0e0_knd) qdr(m+m-1)=dec
        if(m.lt.3) go to 180
          do jn=m-1,2,-1
          rin=jn
          term=rin*rin*x1d-rmsq
            if(term.eq.0.0e0_knd) then
            xc=x*(1.0e0_knd+dec)
            term=rin*rin*(xc*xc+1.0e0_knd)-rmsq
            end if
          qdr(jn+m-1)=(x*(jn+jn-1)*(jn*(jn-1)*x1d-rmsq)+ &
                      ((jn-m)*((jn-1)*(jn-1)*x1d-rmsq))/ &
                      qdr(jn+m))/((jn+m-1)*term)
          if(qdr(jn+m-1).eq.0.0e0_knd) qdr(jn+m-1)=dec
          end do
180     continue
!
!  Calculation of ratios of the first derivative of q with respect
!  to x to the corresponding value for q for degrees from -m to m-1
!  Here
!                  m      m
!     qdqr(k+m) = Q'   / Q   , k=m to k=-m+1
!                  k-1    k-1
!
!
        if(m.eq.0) go to 190
        qdqr(1)=((rm-1.0e0_knd)*x+(rm+rm-1.0e0_knd)/qr(1))/x1d
          do jn=2,m+m
          rin=jn
          qdqr(jn)=(x*(-rm+rin-1.0e0_knd)-(rin-1.0e0_knd)*qr(jn-1))/x1d
          end do
190     continue
!
!  Modification of The ratios qr and qdr to obtain those used in the
!  Baber and Hasse expansion.
!                                     m    m
!  The resulting ratios are qr1(k) = Q  / Q    , k=1 to limq
!                                     k    k-1
!                 m    m
!  and qdr1(k) = Q  / Q    , k=1 to limq
!                 k    k-1
!
          do j=1,limq-m
          qr1(m+j)=qr(m+m+j)
          qdr1(m+j)=qdr(m+m+j)
          end do
          if(m.gt.1) then
            do j=1,m-1
            qr1(j)=-1.0e0_knd/qr(m+j)
            qdr1(j)=-1.0e0_knd/qdr(m+j)
            end do
          end if
          if(m.ne.0) then
          qr1(m)=-(ql(1)/qm1)*(ten**(iql(1)-iqm1))
          qdr1(m)=-(qdl(1)/qdm1)*(ten**(iqdl(1)-iqdm1))
          end if
          if(2*(m/2).eq.m) then
          if(2*(m/4).ne.m/2) qm0=-qm0
          if(2*(m/4).ne.m/2) qdm0=-qdm0
          else
          if(2*((m-1)/4).ne.(m-1)/2) qm0=-qm0
          if(2*((m-1)/4).ne.(m-1)/2) qdm0=-qdm0
          end if
        return
        end subroutine
!
!
       subroutine pint(c,m,lnum,x,limint,maxint,maxlp,maxmp,ndec,nex, &
                       ngqs,ngau,wg,xg,step1,nstep1,step2,nstep2,step3, &
                       nstep3,limi,rpint1,rpint2,pint1,pint2,pint3, &
                       pint4,norme,pnorm,ipnorm,coefme,coefmo)
!
!  purpose:     To calculate integrals of the product of associated
!               Legendre functions and kernels containing spherical
!               Neumann functions and a window function. Four
!               different kernel functions are involved leading to
!               integrals of four different types. The integrals are
!               calculated using Gaussian quadrature.
!
!  parameters:
!
!     input:    c      : c
!               m      : m
!               lnum   : number of l values desired
!               x      : x
!               limint : number of integrals of each of the four types
!                        required
!               maxint : dimension of the integral arrays
!               maxlp  : dimension of characteristic and exponent
!                        arrays of integrals
!               maxmp  : dimension of the spherical Neumann function
!                        array
!               ndec   : number of decimal digits in real(knd)
!                        arithmetic
!               nex    : maximum exponent in real(knd) arithmetic
!               ngqs   : total number of steps
!               ngau   : order of Gaussian quadrature for substeps
!               wg     : ngau Gaussian quadrature weighting factors
!               xg     : corresponding Gaussian quadrature arguments
!               step1  : size of step 1
!               nstep1 : number of substeps in step1
!               step2  : size of step 2
!               nstep2 : number of substeps in step2
!               step3  : size of step 3
!               nstep3 : number of substeps in step3
!
!     output:   limi   : highest order of Legendre function for which
!                        integrals are calculated, since higher order
!                        integrals could overflow. series in subroutine
!                        r2int will be limited to order limi
!               rpint1 : array of ratios of successive integrals of
!                        the same parity of the first type (l-m even)
!                        or of the third type (l-m odd)
!               rpint2 : array of ratios of successive integrals of
!                        the same parity of the second type (l-m even)
!                        or of the fourth type (l-m odd)
!               pint1  : array of scaled values for the integrals of
!                        the first type
!               pint2  : array of scaled values for the integrals of
!                        the second type
!               pint3  : array of scaled values for the integrals of
!                        the third type
!               pint4  : array of scaled values for the integrals of
!                        the fourth type
!               norme  : scaling exponent for the spherical Neumann
!                        functions of order m
!               pnorm  : array of characteristics for the scaling
!                        factors used for the associated Legendre
!                        functions
!               ipnorm : array of exponents for the scaling factors
!                        used for the associated Legendre functions
!               coefme : coefficient used to multiply r2 to get one of
!                        the two contributions to r2d when l-m is even
!               coefmo : coefficient used to multiply r2 to get one of
!                        the two contributions to r2d when l-m is odd
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) amo2,an,arg,argb,arn,bn,c,coef,coefme,coefmo,coefo, &
                  coef1,coef2,coef4,darg,dec,emo2,etai,etaism1,etcoef1, &
                  etcoef2,factor,rm,rn,sargb,sneu1,sneu2,sneu3,step1, &
                  step2,step3,substep1,substep2,substep3,ten,term, &
                  term1,term2,test,teste,test1,testeo,twom,twomi,x, &
                  x2,xsp1,ynorm1,ynorm2,ynorm3,zi,zl,zu
        real(knd) alpha(maxint),beta(maxint),p(maxint),pint1(maxint), &
                  pint2(maxint),pint3(maxint),pint4(maxint), &
                  pnorm(maxlp),rpint1(maxint),rpint2(maxint),step(ngqs), &
                  wg(ngau),xg(ngau)
!
!  integer array
        integer ipnorm(maxlp)
!
        ten=10.0e0_knd
        emo2=(0.5e0_knd)*m
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        ifac=nex-30
        test=ten**ifac
        factor=1.0e0_knd/test
        test1=ten**(nex-10)
        x2=x*x
        xsp1=x2+1.0e0_knd
        dec=ten**(-ndec-1)
        rm=m
        amo2=0.5e0_knd*rm
        lim=limint
        limi=lim
        coefme=rm*x/xsp1
        coefmo=(rm*x2+xsp1)/(x*xsp1)
          substep1=step1/nstep1
            do j=1,nstep1
            step(j)=substep1
            end do
          substep2=step2/nstep2
            do j=1,nstep2
            step(nstep1+j)=substep2
            end do
          substep3=step3/nstep3
            do j=1,nstep3
            step(nstep1+nstep2+j)=substep3
            end do
!
!  calculation of scaling factors for the associated Legendre functions
    pnorm(1)=1.0e0_knd
    pnorm(2)=1.0e0_knd
    ipnorm(1)=0
    ipnorm(2)=0
        if(m.eq.0) go to 50
          do 40 n=1,m
          an=n+n
          bn=n+n+1
          pnorm(1)=pnorm(1)*an
          pnorm(2)=pnorm(2)*bn
          iterm1=int(log10(pnorm(1)))
          iterm2=int(log10(pnorm(2)))
          pnorm(1)=pnorm(1)*ten**(-iterm1)
          pnorm(2)=pnorm(2)*ten**(-iterm2)
          ipnorm(1)=ipnorm(1)+iterm1
          ipnorm(2)=ipnorm(2)+iterm2
40      continue
50    twom=m+m
        pnorm(3)=pnorm(1)*(twom+2)/2
        iterm3=int(log10(pnorm(3)))
        pnorm(3)=pnorm(3)*ten**(-iterm3)
        ipnorm(3)=iterm3+ipnorm(1)
        if(lnum.lt.4) go to 70
          do 60 il=4,lnum,2
          pnorm(il)=pnorm(il-2)*(twom+il-1)/(il-1)
          pnorm(il+1)=pnorm(il-1)*(twom+il)/(il)
          iterm1=log10(pnorm(il))
          iterm2=log10(pnorm(il+1))
          ipnorm(il)=ipnorm(il-2)+iterm1
          ipnorm(il+1)=ipnorm(il-1)+iterm2
          pnorm(il)=pnorm(il)*ten**(-iterm1)
          pnorm(il+1)=pnorm(il+1)*ten**(-iterm2)
60      continue
70    continue
!
!  calculation of the coefficients in the recursion relation used
!  for the scaled associated Legendre functions
    alpha(1)=(twom+1.0e0_knd)*pnorm(1)/pnorm(2)
        alpha(1)=alpha(1)*ten**(ipnorm(1)-ipnorm(2))
        beta(1)=0.0e0_knd
        alpha(2)=(twom+3.0e0_knd)*pnorm(2)/(pnorm(3)*2.0e0_knd)
        alpha(2)=alpha(2)*ten**(ipnorm(2)-ipnorm(3))
        beta(2)=-(twom+1.0e0_knd)/(twom+2.0e0_knd)
          do 80 il=3,lim+2
          alpha(il)=alpha(il-2)*(twom+il-1)*(twom+il+il-1)* &
          (il-2)/((il-1)*(twom+il)*(twom+il+il-5))
      beta(il)=-(twom+il-1)/(twom+il)
80      continue
!
          do 90 il=1,lim+3,2
          pint1(il)=0.0e0_knd
          pint2(il)=0.0e0_knd
          pint3(il+1)=0.0e0_knd
          pint4(il+1)=0.0e0_knd
90      continue
!
!  calculation of the scaling exponents for the spherical Neumann
!  functions required for the four types of integrals
        twomi=1.0e0_knd
        if(m.eq.0) go to 110
          do 100 n=1,m
      twomi=twomi*(n+n-1)/(n+n)
100       continue
110    continue
        arg=c*sqrt(xsp1)
        darg=1.0e0_knd/arg
        ynorm1=-cos(arg)*darg
        ynorm2=ynorm1*darg-sin(arg)*darg
        normy=0
          do 120 n=3,m+3
          rn=n+n-3
          arn=n+n-1
          ynorm3=-ynorm1+darg*rn*ynorm2
            if(abs(ynorm3).gt.teste) then
            ynorm3=ynorm3*testeo
            ynorm2=ynorm2*testeo
            normy=normy+nfac
            end if
          ynorm1=ynorm2
          ynorm2=ynorm3
120       continue
          iterm=int(log10(abs(ynorm3)))
          normy=normy+iterm
          norme=normy
!
!  Gaussian quadrature integration loops. first dividing integrand
!  into ngqs steps
            do 180 k=1,ngqs
              if(k.eq.1) then
              zl=0.0e0_knd
              zu=step(1)
              end if
              if(k.gt.1) then
              zl=zu
              zu=zl+step(k)
              end if
          etcoef1=(zl+zu)/2.0e0_knd
          etcoef2=step(k)/2.0e0_knd
          coef=step(k)
!
!  Gaussian quadrature integration over each step
        do 170 i=1,ngau
         zi=etcoef1+xg(i)*etcoef2
            etai=1.0e0_knd-zi
            etaism1=zi*(2.0e0_knd-zi)
            argb=x2+etaism1
            term2=xsp1*etaism1*etaism1/argb
            term1=sqrt(term2)
            sargb=sqrt(argb)
            coefo=1.0e0_knd/sargb
            arg=c*sargb
            darg=1.0e0_knd/arg
            sneu1=-cos(arg)*darg
            sneu2=term1*(sneu1*darg-sin(arg)/arg)
            norm=ifac
              do 150 n=3,m+3
              rn=n+n-3
              arn=n+n-1
              sneu3=-term2*sneu1+term1*darg*rn*sneu2
              if(n.eq.m+3) go to 155
                if(abs(sneu3).gt.teste) then
                sneu3=sneu3*testeo
                sneu2=sneu2*testeo
                norm=norm+nfac
                end if
                if(abs(sneu3).lt.testeo) then
                sneu3=sneu3*teste
                sneu2=sneu2*teste
                norm=norm-nfac
                end if
              sneu1=sneu2
              sneu2=sneu3
150           continue
155         continue
            iterm=int(log10(abs(sneu3)))
            term=ten**(-iterm)
            sneu3=sneu3*term
            sneu2=sneu2*term
            sneu1=sneu1*term
            norm=norm+iterm
            iexp=-norme+norm
            if(iexp.le.-nex+10-ifac) go to 160
              if(iexp.gt.-nex+10) then
              term=ten**iexp
              sneu3=sneu3*term
              sneu2=sneu2*term
              sneu1=sneu1*term
              coef1=coef*sneu1*wg(i)
              coef2=(coef/term1)*coefo*sneu2*wg(i)
              coef4=(coef/term2)*coefo*coefo*etai*sneu3*wg(i)
              iflag=0
              else
              term=ten**(iexp+ifac)
              sneu3=sneu3*term
              sneu2=sneu2*term
              sneu1=sneu1*term
              iflag=1
              end if
            p(1)=twomi*factor
            p(2)=alpha(1)*etai*p(1)
              if(iflag.eq.0) then
              pint1(1)=pint1(1)+coef1*p(1)
              pint2(1)=pint2(1)+coef2*p(1)
              pint4(2)=pint4(2)+coef4*p(2)
              end if
              do il=2,limi,2
              p(il+1)=alpha(il)*etai*p(il)+beta(il)*p(il-1)
              p(il+2)=alpha(il+1)*etai*p(il+1)+beta(il+1)*p(il)
                if(iflag.eq.0) then
                pint1(il+1)=pint1(il+1)+coef1*p(il+1)
                pint2(il+1)=pint2(il+1)+coef2*p(il+1)
                pint4(il+2)=pint4(il+2)+coef4*p(il+2)
                  if(abs(pint1(il+1)).gt.test1.or.pint2(il+1).gt.test1 &
                     .or.pint4(il+2).gt.test1) then
                  limi=il
                  go to 160
                  end if
                end if
                if(abs(p(il+2)).gt.test) then
                p(il+1)=p(il+1)*factor
                p(il+2)=p(il+2)*factor
                  if(iflag.eq.0) then
                  coef1=coef1*test
                  coef2=coef2*test
                  coef4=coef4*test
                  else
                  coef1=coef*sneu1*wg(i)
                  coef2=(coef/term1)*coefo*sneu2*wg(i)
                  coef4=(coef/term2)*coefo*coefo*etai*sneu3*wg(i)
                  iflag=0
                  end if
                end if
              end do
160         continue
170         continue
180      continue
190    continue
          do 200 il=1,limi,2
          pint3(il+1)=(pint2(il+2)-beta(il+1)*pint2(il)) &
                      /alpha(il+1)
200      continue
210     continue
if (debug) then
        write(40,220) ngau,substep1,substep2,substep3,limi
220     format(15x,'Gauss quad. order =',i5,'; step sizes = ',f12.10, &
               ', ',f12.10,', 'f12.10,'.',/,15x,'integrals for ',i5, &
               ' lowest order Legendre functions will be used for r2.')
end if
!
!  calculation of ratios of integrals for ease in compution of r2 and
!  r2d in subroutine r2int
        rpint1(1)=0.0e0_knd
        rpint1(2)=0.0e0_knd
        rpint2(1)=0.0e0_knd
        rpint2(2)=0.0e0_knd
          do 240 il=3,limi,2
          rpint1(il)=pint1(il)*(twom+il-1)/(pint1(il-2)*(il-1))
          rpint2(il)=pint2(il)*(twom+il-1)/(pint2(il-2)*(il-1))
          rpint1(il+1)=pint3(il+1)*(twom+il)/(pint3(il-1)*(il))
          rpint2(il+1)=pint4(il+1)*(twom+il)/(pint4(il-1)*(il))
240      continue
        return
        end subroutine
!
!
        subroutine sphbes (c,x,limj,maxj,maxlp,sbesf,sbesdf,sbesn, &
                           ibese,sbesdr)
!
!  purpose:     To calculate ratios of successive first kind spherical
!               Bessel functions of the same parity for given c and x.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Bessel functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    c      : c
!               x      : x
!               limj   : the number of spherical Bessel function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxj   : dimension of sbesf vector
!               maxlp  : the number of scale factors that are
!                        calculated
!
!     output:   sbesf  : ratios of successive first kind spherical
!                        Bessel functions of the same parity
!               sbesdf : ratios of first derivatives of successive
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : characteristics for the spherical
!                        Bessel functions
!               ibese  : exponents for the spherical
!                        Bessel functions
!               sbesdr : ratios of first derivatives of spherical Bessel
!                        functions to the corresponding spherical
!                        spherical functions
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,cx,rn,stemp0,stemp1,ten,x
        real(knd) sbesdf(maxj),sbesdr(maxlp),sbesf(maxj),sbesn(maxlp)
!
!  integer array
        integer ibese(maxlp)
!
        ten=10.0e0_knd
        cx=c*x
        lim1=cx+cx+20
!
!  compute first kind Bessel function ratios
!        sbesf(k)= j(n=k,c*x)/j(n=k-1,c*x)
!        sbesn(k)= (j(n=k-1),c*x))*10.0e0_knd**(-ibese(k))
!
        if (cx.lt.limj) go to 20
!
!  for c*x >= lim, use forward recursion to
!  get fcn. ratios:
!       j(n+1,c*x)/j(n,c*x)=(2*n+1)/(c*x)-1/(j(n,c*x)/j(n-1,c*x))
!
        stemp0=sin(cx)/cx
        sbesf(1)=(stemp0/cx-cos(cx)/cx)/stemp0
          do 10 n=1,limj-1
          rn=n+n+1
          sbesf(n+1)=(rn/cx)-(1.0e0_knd/sbesf(n))
10        continue
        go to 60
20      continue
!
!  for c*x < lim, use backward recursion to
!  get fcn. ratios:
!       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
!
        stemp0=0.0e0_knd
        if(lim1.le.limj) go to 40
          do 30 n=lim1,limj,-1
          rn=n+n+1
          stemp1=1.0e0_knd/(rn/cx-stemp0)
          stemp0=stemp1
30        continue
40      sbesf(limj)=stemp0
          do 50 n=limj-1,1,-1
          rn=n+n+1
          sbesf(n)=1.0e0_knd/(rn/cx-sbesf(n+1))
50        continue
60      continue
!
!  for all c*x, calculate the amplitude and sign scale
!  factors by forward operation on the Bessel function
!  ratios.
        stemp0=sin(cx)/cx
        stemp1=stemp0/cx-cos(cx)/cx
        ibese(1)=log10(abs(stemp0))
        sbesn(1)=stemp0*ten**(-ibese(1))
        if(abs(sin(cx)).lt.0.5e0_knd.and.cx.gt.1.0e0_knd) go to 70
        sbesn(2)=sbesn(1)*sbesf(1)
        ibese(2)=log10(abs(sbesn(2)))
        sbesn(2)=sbesn(2)*ten**(-ibese(2))
        ibese(2)=ibese(2)+ibese(1)
        go to 80
70      ibese(2)=log10(abs(stemp1))
        sbesn(2)=stemp1*ten**(-ibese(2))
        sbesf(1)=stemp1/stemp0
80      continue
           do 90 n=3,maxlp
          sbesn(n)=sbesn(n-1)*sbesf(n-1)
          ibese(n)=log10(abs(sbesn(n)))
          sbesn(n)=sbesn(n)*ten**(-ibese(n))
          ibese(n)=ibese(n)+ibese(n-1)
90        continue
!
!  calculate the ratios of the first derivatives of successive
!  Bessel functions using corresponding function ratios
          do 100 n=1,limj
          rn=n-1
          sbesdf(n)=(cx-(rn+2.0e0_knd)*sbesf(n))/(rn-cx*sbesf(n))
100       continue
!
!  calculate the ratios of the first derivative to the corresponding
!  spherical Bessel function
          do 110 n=1,maxlp
          rn=n-1
          sbesdr(n)=(rn/cx)-sbesf(n)
110       continue
!
!  calculate the ratios of successive functions and derivatives
!  of the same parity
          do 120 n=limj,2,-1
          sbesf(n)=sbesf(n-1)*sbesf(n)
          sbesdf(n)=sbesdf(n-1)*sbesdf(n)
120       continue
        return
        end subroutine
!
!
        subroutine sphneu (c,x,limn,maxn,maxlp,sneuf,sneun,ineue,sneudf, &
                           sneudr)
!
!  purpose:     To calculate ratios of spherical Neumann functions
!               and ratios of their first derivatives for given c and x.
!               to calculate the Neumann function characteristics
!               and exponents. To calculate ratios of the first
!               derivatives to corresponding functions.
!
!  parameters:
!
!     input:    c      : c
!               x      : x
!               limn   : the number of spherical Neumann function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxn   : dimension of sneuf and sneudf arrays
!               maxlp  : the number of values of scale factors
!                        that are calculated
!
!     output:   sneuf  : ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : characteristic for the spherical
!                        Neumann functions
!               ineue  : exponent for the spherical
!                        Neumann functions
!               sneudf : ratios of first derivatives of successive
!                        spherical Neumann functions of the same parity
!               sneudr : ratios of first derivatives of spherical
!                        Neumann functions to the corresponding
!                        function
!
        use param
!
!  real(knd) scalars and arrays
        real(knd) c,cx,rn,rnn,stemp0,stemp1,ten,x
        real(knd) sneudf(maxn),sneudr(maxlp),sneuf(maxn),sneun(maxlp)
!
!  integer arrays
        integer ineue(maxlp)
!
!  compute first kind ratios of Neumann functions and ratios
!  of their first derivatives
!
!        sneuf(k)=y(n=k,c*x)/y(n=k-2,c*x)
!        sneun(k)=(y(n=k-1),c*x)*10.0e0_knd**(-ineue(k))
!        sneudf(k)=y'(n=k,c*x)/y'(n=k-2,c*x)
!        sneudr(k)=(y'(n=k-1),c*x)/y(n=k-1),c*x))
!
!  use forward recursion to compute function ratios
!
!       y(n+1,c*x)/y(n,c*x)=(2*n+1)/(c*x)-1/(y(n,c*x)/y(n-1,c*x))
!
!  compute derivative ratios at same time using function ratios.
        ten=10.0e0_knd
        cx=c*x
        stemp0=-cos(cx)/cx
        stemp1=stemp0/cx-sin(cx)/cx
        sneuf(1)=stemp1/stemp0
        sneudf(1)=-(cx-2.0e0_knd*sneuf(1))/(cx*sneuf(1))
          do 10 n=1,limn-1
          rn=n
          rnn=n+n+1
          sneuf(n+1)=rnn/cx-1.0e0_knd/sneuf(n)
          sneudf(n+1)=(cx-(rn+2.0e0_knd)*sneuf(n+1))/(rn-cx*sneuf(n+1))
10        continue
        sneuf(limn+1)=0.0e0_knd
        sneuf(limn+2)=0.0e0_knd
        sneudf(limn+1)=0.0e0_knd
        sneudf(limn+2)=0.0e0_knd
!
!  calculate the characteristic and exponent
!  by forward operation on the Neumann function ratios:
        ineue(1)=int(log10(abs(stemp0)))
        sneun(1)=stemp0*ten**(-ineue(1))
        ineue(2)=int(log10(abs(stemp1)))
        sneun(2)=stemp1*ten**(-ineue(2))
        nlimit=maxlp
        if(maxlp.gt.limn+1) nlimit=limn+1
          do 20 n=3,nlimit
          sneun(n)=sneun(n-1)*sneuf(n-1)
          ineue(n)=int(log10(abs(sneun(n))))
          sneun(n)=sneun(n)*ten**(-ineue(n))
          ineue(n)=ineue(n)+ineue(n-1)
20        continue
!
!  calculate the ratios of the first derivatives to the corresponding
!  spherical Neumann functions
          do 30 n=1,nlimit
          rn=n-1
          sneudr(n)=(rn/cx)-sneuf(n)
30        continue
!
!  calculate the ratios of successive functions and derivatives
!  of the same parity
          do 40 n=limn,2,-1
          sneuf(n)=sneuf(n-1)*sneuf(n)
          sneudf(n)=sneudf(n-1)*sneudf(n)
40        continue
        return
        end subroutine

