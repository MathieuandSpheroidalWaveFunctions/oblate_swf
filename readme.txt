                               Oblate_swf

  Oblfcn is available as both a subroutine version provided as the
  module oblate_swf and a stand alone version oblfcn. It was originally
  developed by arnie lee van buren over 15 years years ago and has been
  updated and improved many times since.

  Table of Contents
  1. Purpose
  2. Introduction
  3. Input and Output
  4. Accuracy of results
  5. Obtaining the expansion d coefficients
  6. Obtaining the eigenvalues


 1. Purpose

  To calculate the first and second kind oblate radial functions r1
  and r2 and their first derivatives r1d and r2d for a given order m,
  a range of degrees l beginning at l = m and for a specific size parameter
  c and shape parameter x. To calculate the first kind oblate angular
  functions and their first derivatives with respect to eta for the same
  values of m, l and c and for a set of values of the angular coordinate
  eta. The subroutine version oblate_swf calculates values for a single
  input value of m . The stand alone version calculates values for a
  range of values of m.

 2. Introduction

  Oblfcn is written in free format fortran. It is designed around the
  maximum number of decimal digits ndec and the maximum exponent nex
  available in real arithmetic. Procedures used in oblfcn allow for
  exponents much larger than nex since the resulting floating point
  radial function values are given as a characteristic and an integer
  exponent.

  Oblfcn can be run in either double precision or quadruple precision
  arithmetic. The choice is set in the module param provided in the github
  repository. If this is not available, then create param as follows:

    module param
    integer, parameter :: knd = selected_real_kind(8)
    logical, parameter :: debug = .true.
    logical, parameter :: warn = .true.
    logical, parameter :: output = .true.
    logical, parameter :: suffix = .true.
    end module param

  Set the value of knd in the parenthesis to either 8 for double precision
  or 16 for quadruple precision arithmetic. The logicals in param are
  described below in the discussion of the output files.

  Some computers may have more than 8 bytes for double precision
  data and more than 16 bytes for quadruple precision data or may use
  kind values that do not correspond to the number of bytes. In this
  case just use the appropriate integers for the kind parameters in
  module param. Also change the values of kindd and kindq set in
  statement 5 below below the comments section to the kind values for
  double precision data and quadruple precision data, respectively.

  A description of the methods used in oblfcn is provided in the manuscript
  'Accurate calculation of oblate spheroidal wave functions,' available at
  arXiv.org, identifier 1708.07929, August 2017; revised September 2019.

  3. Input and Output

  Following is a description of the input and output parameters in the
  call statement for the subroutine version. After that will be a
  description of the the input and output files associated with the
  stand alone version. Note that these output files, if desired, can
  also be obtained when running the subroutine version. See comments about
  this below.

  A sample input and resulting output from oblfcn is provided by the
  files oblfcndat (text version of the input file oblfcn.dat for the
  stand alone version), oblfort20 (text version of the output file
  fort.20 giving the resulting radial functions) and oblfort30 (text
  version of the output file fort.30 giving the resulting angular
  functions).

  Subroutine Version of oblfcn

    subroutine oblfcn(c,m,lnum,ioprad,x,iopang,iopnorm,narg,arg, &
                      r1c,ir1e,r1dc,ir1de,r2c,ir2e,r2dc,ir2de,naccr, &
                      s1c,is1e,s1dc,is1de,naccs)

        real(knd), intent (in)  ::  c, x, arg(narg)
        integer, intent (in)    ::  m, lnum, ioprad, iopang, iopnorm, narg
        real(knd), intent (out) ::  r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), &
                                    s1c(lnum, narg), s1dc(lnum, narg)
        integer, intent (out)   ::  ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), &
                                    is1e(lnum, narg), is1de(lnum, narg), & 
                                    naccr(lnum), naccs(lnum, narg)

      Input and output parameters appearing in the subroutine call
      statement are defined below:

          c      : desired value of the size parameter (= kd/2, where
                   k = wavenumber and d = interfocal length) [real(knd)]
          m      : desired value for the order m (integer)
          lnum   : number of values desired for the degree l equal
                   to m, m + 1, m + 2, ..., m + lnum - 1 (integer)
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed
          x      : value of the radial coordinate x [real(knd)]
          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed
          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      Legendre function [i.e., we use the Meixner and
                      Schafke normalization scheme.] This norm
                      becomes very large as m becomes large. The
                      angular functions are computed below as
                      a characteristic and an exponent to avoid
                      overflow.
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.
                      This is very useful since it removes the
                      need to calculate a normalization factor
                      when using the angular function values given
                      here. It also eliminates any chance for
                      overflow when the characteristics and exponents
                      are combined to form the angular functions.
          narg   : number of values of the angular coordinate eta for
                   which angular functions are calculated (integer)
          arg    : vector containing the narg values of eta for which
                   angular functions are desired [real(knd)]
          r1c/   : vectors of length lnum containing the characteristics
          r1dc     for the radial functions of the first kind r1 and their
                   first derivatives [real(knd)]
          ir1e/  : integer vectors of length lnum containing the
          ir1de    exponents corresponding to r1c and r1dc
          r2c/   : vectors of length lnum containing the characteristics
          r2dc     for the radial functions of the second kind r2 and their
                   first derivatives [real(knd)]
          ir2e/  : integer vectors of length lnum containing the
          ir2de    exponents corresponding to r2c and r2dc
          naccr  : integer vector of length lnum containing the estimated
                   accuracy of the radial functions
          s1c,   : two-dimensional arrays s1c(lnum,narg) and
          s1dc     s1dc(lnum,narg) that contain narg calculated
                   characteristics for the angular functions and
                   their first derivatives for each of the lnum
                   values of l [real(knd)]
                   For example, s1(10,1) is the characteristic
                   of the angular function for l = m +10 -1 and
                   for the first value of eta given by arg(1)
          is1e,    integer arrays is1e(lnum,narg) and is1de(lnum,narg)
          is1de    containing the exponents corresponding to s1c and
                   s1dc
          naccs  : two-dimensional array naccs(lnum,narg) that contains
                   narg estimated accuracy values for the angular functions
                   for each of the lnum values of l

  Stand alone Version of oblfcn

     Input Data

     Input parameters are read from unit 1 in the file oblfcn.dat
     assumed to be in the directory of oblfcn.f90. oblfcn.dat
     contains the following lines of data:

       line 1:
          mmin   : minimum value for m. (integer)
          minc   : increment for m. (integer)
          mnum   : number of values of m. (integer)
          lnum   : number of values of l [l=m, l=m+1,
                   ..., l=m+lnum-1]. (integer)
                   if lnum is less than 2*c/pi it should
                   be an even integer. If lnum is chosen to
                   be odd, oblfcn will increase lnum by one.

       line 2:
          ioprad : (integer)
                 : =0 if radial functions are not computed
                 : =1 if radial functions of only the first kind
                      and their first derivatives are computed
                 : =2 if radial functions of both kinds and
                      their first derivatives are computed

          iopang : (integer)
                 : =0 if angular functions are not computed
                 : =1 if angular functions of the first kind
                      are computed
                 : =2 if angular functions of the first kind and
                      their first derivatives are computed

          iopnorm: (integer)
                 : =0 if not scaled. The angular functions have
                      the same norm as the corresponding associated
                      Legendre function [i.e., we use the Meixner-
                      Schafke normalization scheme.]
                 : =1 if angular functions of the first kind
                      (and their first derivatives if computed)
                      are scaled by the square root of the
                      normalization of the corresponding
                      associated Legendre function. The resulting
                      scaled angular functions have unity norm.

       line 3:
          c      : value of the size parameter (= kd/2, where k =
                   wavenumber and d = interfocal length) [real(knd)]
          x      : value of the radial coordinate x [real(knd)]
                   (any value can be entered if ioprad = 0)

       line 4:
          ioparg : (integer)
                 : =0 if both arg1 and darg are angles in degrees
                 : =1 if arg1 and darg are dimensionless values of eta

          arg1   : first value for the angle coordinate (in degrees
                   or dimensionless if eta) for which angular
                   functions are to be computed. [real(knd)]

          darg   : increment used to calculate additional desired
                   arguments for angular functions. [real(knd)]

          narg   : number of desired angle arguments. (integer)
                   (line 4 is not read when iopang = 0)

     Output files

     These output files are also available using the subroutine version
     oblate_swf. Generation of each of the files is controlled by a logical
     specified in the module param located before the program. False
     suppresses the output file and true enables it. The logical debug
     controls fort.30 and fort.40, the logical output controls fort.20
     and fort.30 and warn controls fort.60. The logical suffix controls
     whether the accuracy estimates given in fort.20 are followed by a letter
     designating how the accuracy was determined. 'w' indicates it is based
     on the Wronskian and 'e' indicates it is based on subtraction errors
     involved in the calculations. Setting suffix = false suppresses the letter. 

   fort.20

     This file contains values for all radial functions that have
     been calculated.
     The first line in the file contains the values for x, c, and
     m, formatted as follows (see statements 260 and 265 in
     subroutine main):

                x      : e23.14 in real*8; e38.30 in real*16
                c      : e23.14 in real*8; e38.30 in real*16
                m      : i5

     Each subsequent line in fort.20 contains radial functions
     for given values of l. The first line contains values for l = m,
     the next for l=m+1 and continuing to l=m+lnum-1. The radial
     functions are preceeded by the value of l and followed by the
     accuracy, equal to the estimated number of accurate decimal digits
     in the radial functions. (see comments below regarding naccr).

       The output and corresponding format for each line is as follows
       (see statements 1340, 1350, 1355, 1380, and 1390 in main).

         l            : value for l (i5)
         r1c(l-m+1)   : characteristic of the oblate radial
                        function of the first kind r1 for
                        the given value of l (f17.14)
         ir1e(l-m+1)  : exponent of r1 (i6)
         r1dc(l-m+1)  : characteristic of the first derivative
                        of the oblate radial function of the
                        first kind r1d (f17.14)
         ir1de(l-m+1) : exponent of r1d (i6)
         r2c(l-m+1)   : characteristic of the oblate radial
                        function of the second kind r2 for
                        the given value of l (f17.14)
         ir2e(l-m+1)  : exponent of the oblate radial function of
                        second kind (i6). If the exponent for any
                        of the functions is greater than 99999,
                        the format can be increased to i7 or
                        higher. Note that the procedures used in
                        this program allow for exponents much
                        larger than those allowed in real*8
                        or real*16 arithmetic on the users computer
                        since the floating point function values
                        provided are given as a characteristic
                        and an integer exponent. Use of ratios in
                        calculating the functions eliminates
                        overflow and underflow during the
                        calculations.
         r2dc(l-m+1)  : characteristic of the first derivative of
                        the oblate radial function of second kind
                        r2d (f17.14)
         ir2de(l-m+1) : exponent of r2d (i6). [See comment above
                        for ir2e.]
         naccr(l-m+1) : estimated accuracy: often equal to the
                        number of decimal digits of agreement
                        between the theoretical Wronskian and the
                        calculated Wronskian (i2). This is
                        indicated by the letter w following the
                        integer naccr. Since r1 and r1d are
                        almost always more accurate than r2 and r2d,
                        unless very near a root, naccr is usually
                        an estimate of the accuracy of r2 and r2d.
                        When the Wronskian is not used for the
                        accuracy estimate but instead a more
                        conservative estimate is used, it is
                        indicated by the letter e following the
                        integer naccr.

                        When leading coefficients used in the
                        traditional Legendre function method
                        for computing r2 and r2d are obtained using
                        the Wronskian, a following letter e is used
                        to indicate this. See the comments above
                        describing these situations.

                        When x is equal to zero, naccr is approximated
                        by ndec-jsub-3 where jsub is the subtraction
                        error involved in calculating the Flammer
                        normalization factor dfnorm (see subroutine
                        dnorm for a definition of dfnorm). The reduced
                        accuracy applies only to r2 for l-m even and
                        r2d for l-m odd. When the accuracy estimate
                        is zero digits, the relevant function is set
                        equal to zero.

   fort.30

     This file fort.30 contains values for all angular functions
     that have been calculated. Its first line contains the values
     for c and m, formatted as follows (see statements 50 and 55 in
     subroutine main).

                c      : e23.14 in real*8; e38.30 in real*16
                m      : i5

     The second line in fort.30 contains the value for the first l
     (=m), formatted as follows (see statement 270 in subroutine
     main):

                l      : i6

     This is followed by a series of narg lines. Each line contains
     a desired value of angle (ioparg = 0) or angular coordinate eta
     (ioparg =1) followed by the corresponding angular functions and
     accuracy. Specific output and format for each line is as follows.

        for iopang = 1:

               arg    : for ioparg = 0, angle in degrees (f17.14; see
                        statement 1440 in subroutine main)
            or barg   ; for ioparg = 1, angular coordinate eta
                        (f17.14; see statement 1460 in subroutine main)
               s1c    : characteristic of the oblate angular function
                        of first kind (f17.14; see statements 1440 and
                        1460 in subroutine main)
               is1e   : exponent of the oblate angular function of
                        first kind (i5; see statements 1440 and 1460 in
                        subroutine main)

        for iopang = 2, each line also includes:
               s1dc   : characteristic of the first derivative of the
                        oblate angular function of first kind (f17.14;
                        see statements 1450 and 1470 in subroutine
                        main)
               is1de  : exponent of the first derivative of the
                        oblate angular function of first kind (i5;
                        see statements 1450 and 1470 in subroutine
                        main)

        for iopang = 1:
               naccs  : accuracy: estimate of the number of decimal
                        digits of accuracy in the angular function
                        (i2). It is a conservative estimate based on
                        the calculated subtraction error in the series
                        calculation of the angular function. When the
                        accuracy estimate is equal to 0, the
                        corresponding angular functions are set equal
                        to zero. (i2; see statements 1440, 1450, 1460,
                        and 1470 in subroutine main). The calculated
                        angular functions tend to be less accurate the
                        larger the value of c, the smaller the value of
                        l - m (for values less than approximately 2c
                        divided by pi), and the closer eta is to zero
                        (i.e., the closer theta is to 90 degrees). For
                        c very large and l - m small, even the angular
                        functions near eta = 1 can be very inaccurate.
                        However, the loss of accuracy (in decimal
                        digits) is due to subtraction error and is
                        accompanied by a proportional decrease in the
                        magnitude of the angular function relative to
                        its corresponding associated Legendre function.
                        The decrease in magnitude almost always results
                        in a corresponding reduction in the magnitude
                        of their contribution to the solution of
                        physical problems involving these spheroidal
                        functions. Thus the lower accuracy in some of
                        the angular functions almost always has
                        insignificant impact on calculated solutions.

        for iopang = 2:
              naccds  : accuracy: includes an estimate of the number
                        of decimal digits of accuracy in the first
                        derivative of the angular function (i2)

   fort.40 and fort.50

     These files are diagnostic files that contain information
     about specific techniques used and numbers of terms required
     for the radial function and angular function calculations,
     respectively. They are annotated and should be self
     explanatory.

   fort.60

     This file may be of interest to the user, especially when using
     this program for values of c greater than 10000 or m greater
     than 1000. Whenever the estimated accuracy falls below a
     designated integer value during the running of this program,
     the associated values of x, c, m, and l are written to fort.60.
     The integer is currently set equal to 8 in the write statements
     for this file found after the line numbered 1405 in subroutine main.

  4. Accuracy of Results

  Oblfcn provides good results using double precision arithmetic
  for values of c up to at least 10000, m up to at least 1000 and
  essentially all values of x. It is expected that oblfcn will provide
  good results for much higher values of m and c.

  The radial functions of the first kind and their first derivatives
  are almost always very accurate for all choices of real arithmetic.

  Using real*8 arithmetic, the accuracy of the radial functions of
  the second kind r2 and their first derivatives r2d is nearly always
  at least 8 decimal digits for x as small as 0.0001. Accuracies
  greater than 10 digits are usually obtained. An accuracy of fewer
  than 8 digits can sometimes occur, especially when l is near the
  breakpoint and c is very large. I define the breakpoint as that value
  of l above which the radial functions of the first kind start to
  decrease in magnitude as the same time as the radial functions of the
  second kind start to increase in magnitude. An estimate for this when
  m is small is 2*c/pi. The accuracy is usually estimated by comparing
  the theoretical value 1/[c(x1*x1+1)] for the Wronskian to its value
  calculated from r1*r2d - r2*r1d. An alternative accuracy estimate is
  provided for the case where the Wronskian is used to obtain leading
  coefficients in the traditional Legendre function expansion. An
  alternative estimate is also given for the special case x = 0 (see below).

  Oblfcn was tested using real*8 arithmetic with a precision of 15
  decimal digits. Testing included values for c equal to 0.0001, 0.001,
  0.01, 0.1, 1, 10, 50, 100, 200, 500, 1000, 2000, 5000 and 10000, for
  x equal to 0.0001, 0.001, 0.01, 0.1 and 1.0, for orders m from 0 to
  200 in unit steps and from 210 to 1000 in steps of 10, and for
  degrees l from m to m - 1 + lnum, where lnum is sufficiently large
  so that the magnitudes of r1 and r1d are less than 10**(-300). It is
  unlikely that one wwould need function values for degrees higher than
  these. However, oblfcn is expected to provide accurate results for
  much higher values of lnum.

  For c <= 1000, the accuracy estimates were at least 8 decimal digits
  except for ten 7 digit results and one 6 digit result, all of which
  occurred when one of the radial functions was near a root.

  For c = 2000, there were 43 7 digit results and five 6 digit results.
  Nearly all of the 6 and 7 digit results occurred near a root.

  For c = 5000, there were 125 7 digit results and 12 6 digit results.
  There was even one 5 digit result at m = 810. Over half of the 6 and
  7 digit results occurred at values of m well above 200. Many of them
  occurred near a root.

  For c = 10000, there were 724 7 digit results and 61 6 digit results.
  Some of the 6 and 7 digit results occurred near a root. Over half
  of the 6 and 7 digit results occurred at values of m well above 200.      

  The number of individual accuracy estimates obtained during testing for
  c = 2000 was about 3 million. Even more estimates were obtained for
  c = 5000 and 10000.

  It is expected that occasional accuracies of 5 and 6 digits is probably
  adequate for most applications of these functions. Also, when one of the
  radial functions is less accurate due to being near a root, its magnitude
  is somewhat smaller than nearby radial functions and thus makes a reduced
  contribution to the solution of problems using these functions.

     Accuracy Estimates For Small x

  Testing for values of x smaller than 0.0001 produced accuracy
  estimates similar for those discussed above. However, for x less
  than about 0.1 the accuracy of r2 for l - m even and of r2d for
  l - m odd can be less than the estimated accuracy when c is large
  and l is very near but somewhat below the breakpoint. Here as x
  becomes smaller, r2 becomes progressively smaller in magnitude than
  r1 for l - m even and r2d becomes progressively smaller than r1d for
  l - m is odd. The lower accuracy is most pronounced when the near
  equality of low-order neighboring eigenvalues is used to obtain r2
  and r2d. (See below for a brief discussion of this method). Here, a
  good accuracy estimate for r2 for l - m odd and r2d for l - m even
  is given by the number of leading digits of agreement, called match,
  between neighboring paired eigenvalues. Match is a good estimate of
  the accuracy of r2d for l - m even and r2 for l - m odd. Oblfcn
  does not use match for this estimate but instead uses the Wronskian.
  Here computed r1 and r1d values for a given l - m are combined with
  r2 and r2d values obtained using accurate r1 values from l - m + 1
  if l - m is even or from l - m - 1 is l - m is odd. A conservative
  estimate of the corresponding accuracy of r2 for l - m even and r2d
  for l - m odd is obtained by first taking the absolute value of the
  logarithm of x to the base 10, truncating it to an integer, and then
  subtracting the result from match. For example, for x = 0.000001 and
  match = 12 digits, oblfcn reports an accuracy of 12 decimal digits
  using the Wronskian. An estimate of the accuracy of r2 for l - m even
  or r2d for l - m odd is 6 digits. For example, comparison of real*8
  results with real*16 results for c = 2000 show that for a given value
  of m there were at most one or two values of l where the accuracy of
  r2 for l - m even was as low as 6 or 7 digits. And this reduced
  accuracy occurs only when r2 is smaller in magnitude than r1 by about
  3 orders of magnitude. This gives an effective accuracy of r1 of 9 or
  10 digits when forming r3 or r4 for applications involving these
  functions. The same thing occurs for r2d when l - m is odd. As
  x becomes smaller than 0.000001, the reduced accuracy for a few
  values of l is balanced by a similar reduction in relative magnitude
  of the less accurate function values. Thus their contribution to the
  solutions of physical problems involving these functions is reduced
  proportionately. Their absolute accuracy decreases but their
  effective accuracy remains sufficient for most applications.

  Using real*16 arithmetic can provide greater accuracy, typically
  but not always 25 or more decimal digits. However, using real*16
  arithmetic increases the computation time considerably.

  x = 0 is a special case of the behavior described above. Here when
  c is large, l - m is even and l is well below the breakpoint value,
  r2 can be much smaller in magnitude than r2d and less accurate
  according to its reduced magnitude. As l - m increases toward the
  breakpoint value, r2 for l - m even increases in both magnitude and
  accuracy until it is similar to that of r2d. The same behavior is
  seen for r2d when l - m is odd. It is expected that the lower
  accuracy of these relatively small function values will not have a
  significant effect on the accuracy of numerical results calculated
  for physical problems involving these functions. The accuracy
  estimate given for the radial functions when x = 0 reflects the
  reduced accuracy of r2 for l - m even and r2d for l - m odd.
  Corresponding values for r2 when l - m is odd and r2d when l - m is
  even are very accurate for all choices of real arithmetic.

  For all choices of real arithmetic, the angular functions and their
  first derivatives have reduced accuracy only for low l - m, high c,
  and eta near zero or when near a root. However, their magnitude in
  this case is small relative to angular functions for higher values
  of l and/or eta not near zero. The reduced accuracy should not
  adversely effect numerical results for physical problems using
  these functions.

  Oblfcn takes advantage of the near equality of r2 and r2d for l - m
  even to r1 and r1d for l - m + 1 when c is large and l is well below
  the breakpoint. This behavior results from the fact that neighboring
  eigenvalues are nearly equal here. Similarly r2 and r2d for l - m + 1
  are nearly equal to the negative of r1 and r1d for l - m. Since the
  radial functions of the first kind are very accurate, using r1 and
  r1d values when the eigenvalues nearly match can provide values for
  r2 and r2d of sufficient accuracy. This is done in oblfcn for all
  values of l up until the number of decimal digits that neighboring
  eigenvalues match falls below the desired number of decimal digits
  plus 2. When x is very small and l is below but very near the break-
  point, this can lead to values for r2 for l - m even and r2d for
  l - m that are less accurate than the Wronskian estimate. See the
  discussion above about the accuracy for small values of x.

  Oblfcn now has improved capability in computing radial functions
  of the second kind r2 and their first derivatives r2d using the
  traditional Legendre function method for x very small (less than or
  equal to 0.01). When c is large and l is very near but below the
  so-called breakpoint, direct calculation of both the leading
  coefficient in the sum of P Legendre functions and the d coefficient
  with negative index appearing in the joining factor can suffer
  significant subtraction error. When these errors are so large that
  r2 and r2d would be less accurate than desired, the Wronskian is
  used to obtain a more accurate value for these coefficients. A
  good estimate of the resulting accuracy for r2 and r2d is obtained
  by considering the various subtraction errors involved in the
  computations.

  An integer called minacc is used in oblfcn to designate the
  minimum number of accurate digits desired for the radial spheroidal
  functions of the second kind. The value of minacc controls which
  methods are used to calculate the radial functions. Minacc is set
  equal to 10 for real*8 arithmetic. It is recommended that this not
  be changed. Minacc is set equal to 15 for real*16 arithmetic. This
  typically corresponds to the number of digits available for real*8
  arithmetic. If more accuracy is desired, minacc can be increased.
  This will usually result in more computation time. The value for
  minacc is set in oblfcn and oblate_swf following the introductory
  comments. 

  5. Obtaining the d expansion coefficients

  The user may desire values for the d coefficients that appear in
  the expression for the angular functions as well as in many of the
  expressions used to calculate the radial functions. Ratios of
  successive d coefficients are stored in the vector enr where enr(k)
  = d(subscript 2k+ix) divided by d(subscript 2k-2+ix). The vector enr
  is calculated in the subroutine dnorm in statement 20 and passed to
  subroutine main. The number lim2 of d coefficients calculated for a
  given l is chosen to be sufficient to compute radial and angular
  functions for that l. The size of Lim2 necessary to compute r1, r1d
  and s1, s1d and the normaliation factors ranges for low l from less
  than 100 for low c and somewhat less than int(c) for large c. Lim2
  increases with increasing l, eventually becoming less than l in size.
  The size of lim2 needed to compute r2 and r2d can be comparable to
  this unless they are computed using Neumann function expansions. Then
  lim2 can be much larger than 200000 or so when x is near its lower
  limit of 0.005 where Neumann expansions may be used. Note that the
  vector enr returned by the subroutine conver contains scaled ratios
  where the scaling has been chosen to produce a symmetric matrix for
  computing  eigenvalues. The scaling factors are removed in subroutine
  dnorm to obtain the desired d coefficient ratios.

  The d coefficients themselves can be obtained starting with the value
  for d with the subscript l - m. If iopnorm is set = 0, Oblfcn uses
  the Meixner-Schafke scheme for normalizing the angular functions.
  Here they have the same norm as the corresponding associated Legendre
  functions. Computation of this normalization is very accurate since
  the series involved has only positive terms. The subroutine dnorm
  computes d(subscript l-m) for this normalization and returns it as
  a characteristic dmlms and an exponent idmlmse to subroutine main.
  Use of an exponent avoids possible overflow of d(subscript l-m) for
  extremely large c and m. When the user sets iopnorm = 1 so that the
  angular functions have unit norm, the corresponding characteristic
  and exponent for d(subscript l-m) are calculated in subroutine s1
  and returned to subroutine main as dmlms1 and idmlms1e. Values for
  the characteristics and exponents of d(subscript l-m) for the  Morse-
  Feshbach and Flammer normalizations are computed in dnorm and
  returned to main as dmlmf, idmlmfe and dmlf, idmlfe. Calculation of
  the Flammer normalization suffers subtraction errors for lower values
  of l-m and large c that increase as c increases. The value for dmlf
  will have reduced accuracy in this case.  

  6. Obtaining the eigenvalues

  The eigenvalues for the oblate functions are computed in subroutine
  conver and returned to main where they are stored in the vector eig(l+1).
  There is such a vector created for each value of m.
