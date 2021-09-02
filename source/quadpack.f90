subroutine qk61 ( f,a,b,result,abserr,resabs,resasc)!,soc,redmass,grad,gradmean,kbt) 

!*****************************************************************************80
!
!! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) A, B, the limits of integration.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!                    result is computed by applying the 61-point
!                    Kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point Gauss rule (resg).
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 4 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 4 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 30-point Gauss rule
!           resk   - result of the 61-point Kronrod rule
!           reskh  - approximation to the mean value of f
!                    over (a,b), i.e. to i/(b-a)
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) absc
  real ( kind = 4 ) abserr
  real ( kind = 4 ) b
  !double precision soc
  !double precision redmass
  !double precision grad
  !double precision gradmean
  !double precision kbt
  real ( kind = 4 ) centr
  real ( kind = 4 ) dhlgth
  real ( kind = 4 ), external :: f
  real ( kind = 4 ) fc
  real ( kind = 4 ) fsum
  real ( kind = 4 ) fval1
  real ( kind = 4 ) fval2
  real ( kind = 4 ) fv1(30)
  real ( kind = 4 ) fv2(30)
  real ( kind = 4 ) hlgth
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtw
  integer ( kind = 4 ) jtwm1
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) resg
  real ( kind = 4 ) resk
  real ( kind = 4 ) reskh
  real ( kind = 4 ) result
  real ( kind = 4 ) wg(15)
  real ( kind = 4 ) wgk(31)
  real ( kind = 4 ) xgk(31)
!
!           the abscissae and weights are given for the
!           interval (-1,1). because of symmetry only the positive
!           abscissae and their corresponding weights are given.
!
!           xgk   - abscissae of the 61-point Kronrod rule
!                   xgk(2), xgk(4)  ... abscissae of the 30-point
!                   Gauss rule
!                   xgk(1), xgk(3)  ... optimally added abscissae
!                   to the 30-point Gauss rule
!
!           wgk   - weights of the 61-point Kronrod rule
!
!           wg    - weigths of the 30-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
     xgk(9),xgk(10)/ &
       9.994844100504906E-01,     9.968934840746495E-01, &
       9.916309968704046E-01,     9.836681232797472E-01, &
       9.731163225011263E-01,     9.600218649683075E-01, &
       9.443744447485600E-01,     9.262000474292743E-01, &
       9.055733076999078E-01,     8.825605357920527E-01/
  data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),xgk(17), &
    xgk(18),xgk(19),xgk(20)/ &
       8.572052335460611E-01,     8.295657623827684E-01, &
       7.997278358218391E-01,     7.677774321048262E-01, &
       7.337900624532268E-01,     6.978504947933158E-01, &
       6.600610641266270E-01,     6.205261829892429E-01, &
       5.793452358263617E-01,     5.366241481420199E-01/
  data xgk(21),xgk(22),xgk(23),xgk(24),xgk(25),xgk(26),xgk(27), &
    xgk(28),xgk(29),xgk(30),xgk(31)/ &
       4.924804678617786E-01,     4.470337695380892E-01, &
       4.004012548303944E-01,     3.527047255308781E-01, &
       3.040732022736251E-01,     2.546369261678898E-01, &
       2.045251166823099E-01,     1.538699136085835E-01, &
       1.028069379667370E-01,     5.147184255531770E-02, &
       0.0E+00                   /
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10)/ &
       1.389013698677008E-03,     3.890461127099884E-03, &
       6.630703915931292E-03,     9.273279659517763E-03, &
       1.182301525349634E-02,     1.436972950704580E-02, &
       1.692088918905327E-02,     1.941414119394238E-02, &
       2.182803582160919E-02,     2.419116207808060E-02/
  data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),wgk(17), &
    wgk(18),wgk(19),wgk(20)/ &
       2.650995488233310E-02,     2.875404876504129E-02, &
       3.090725756238776E-02,     3.298144705748373E-02, &
       3.497933802806002E-02,     3.688236465182123E-02, &
       3.867894562472759E-02,     4.037453895153596E-02, &
       4.196981021516425E-02,     4.345253970135607E-02/
  data wgk(21),wgk(22),wgk(23),wgk(24),wgk(25),wgk(26),wgk(27), &
    wgk(28),wgk(29),wgk(30),wgk(31)/ &
       4.481480013316266E-02,     4.605923827100699E-02, &
       4.718554656929915E-02,     4.818586175708713E-02, &
       4.905543455502978E-02,     4.979568342707421E-02, &
       5.040592140278235E-02,     5.088179589874961E-02, &
       5.122154784925877E-02,     5.142612853745903E-02, &
       5.149472942945157E-02/
  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
       7.968192496166606E-03,     1.846646831109096E-02, &
       2.878470788332337E-02,     3.879919256962705E-02, &
       4.840267283059405E-02,     5.749315621761907E-02, &
       6.597422988218050E-02,     7.375597473770521E-02/
  data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/ &
       8.075589522942022E-02,     8.689978720108298E-02, &
       9.212252223778613E-02,     9.636873717464426E-02, &
       9.959342058679527E-02,     1.017623897484055E-01, &
       1.028526528935588E-01/

  centr = 5.0E-01*(b+a)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 61-point Kronrod approximation to the integral,
!  and estimate the absolute error.
!
  resg = 0.0E+00
  fc = f(centr)!,soc,redmass,grad,gradmean,kbt)
  resk = wgk(31)*fc
  resabs = abs(resk)

  do j = 1, 15
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)!,soc,redmass,grad,gradmean,kbt)
    fval2 = f(centr+absc)!,soc,redmass,grad,gradmean,kbt)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 15
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)!,soc,redmass,grad,gradmean,kbt)
    fval2 = f(centr+absc)!,soc,redmass,grad,gradmean,kbt)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0E-01
  resasc = wgk(31)*abs(fc-reskh)

  do j = 1, 30
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00 .and. abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) / (5.0E+01* epsilon ( resabs ) )) then
    abserr = max ( ( epsilon ( resabs ) *5.0E+01)*resabs, abserr )
  end if

  return
end subroutine
!----------------------------------------------------------------------
subroutine qng ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QNG estimates an integral, using non-adaptive integration.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!    The routine is a simple non-adaptive automatic integrator, based on
!    a sequence of rules with increasing degree of algebraic
!    precision (Patterson, 1968).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) A, B, the limits of integration.
!
!    Input, real ( kind = 4 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!    RESULT is obtained by applying the 21-point Gauss-Kronrod rule (RES21)
!    obtained  by optimal addition of abscissae to the 10-point Gauss rule
!    (RES10), or by applying the 43-point rule (RES43) obtained by optimal
!    addition of abscissae to the 21-point Gauss-Kronrod rule, or by 
!    applying the 87-point rule (RES87) obtained by optimal addition of
!    abscissae to the 43-point rule.
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 4 ) NEVAL, the number of times the integral was evaluated.
!
!           ier    - ier = 0 normal and reliable termination of the
!                            routine. it is assumed that the requested
!                            accuracy has been achieved.
!                    ier > 0 abnormal termination of the routine. it is
!                            assumed that the requested accuracy has
!                            not been achieved.
!                    ier = 1 the maximum number of steps has been
!                            executed. the integral is probably too
!                            difficult to be calculated by qng.
!                        = 6 the input is invalid, because
!                            epsabs < 0 and epsrel < 0,
!                            result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           centr  - mid point of the integration interval
!           hlgth  - half-length of the integration interval
!           fcentr - function value at mid point
!           absc   - abscissa
!           fval   - function value
!           savfun - array of function values which have already
!                    been computed
!           res10  - 10-point Gauss result
!           res21  - 21-point Kronrod result
!           res43  - 43-point result
!           res87  - 87-point result
!           resabs - approximation to the integral of abs(f)
!           resasc - approximation to the integral of abs(f-i/(b-a))
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) absc
  real ( kind = 4 ) abserr
  real ( kind = 4 ) b
  real ( kind = 4 ) centr
  real ( kind = 4 ) dhlgth
  real ( kind = 4 ) epsabs
  real ( kind = 4 ) epsrel
  real ( kind = 4 ), external :: f
  real ( kind = 4 ) fcentr
  real ( kind = 4 ) fval
  real ( kind = 4 ) fval1
  real ( kind = 4 ) fval2
  real ( kind = 4 ) fv1(5)
  real ( kind = 4 ) fv2(5)
  real ( kind = 4 ) fv3(5)
  real ( kind = 4 ) fv4(5)
  real ( kind = 4 ) hlgth
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ipx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) neval
  real ( kind = 4 ) result
  real ( kind = 4 ) res10
  real ( kind = 4 ) res21
  real ( kind = 4 ) res43
  real ( kind = 4 ) res87
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) reskh
  real ( kind = 4 ) savfun(21)
  real ( kind = 4 ) w10(5)
  real ( kind = 4 ) w21a(5)
  real ( kind = 4 ) w21b(6)
  real ( kind = 4 ) w43a(10)
  real ( kind = 4 ) w43b(12)
  real ( kind = 4 ) w87a(21)
  real ( kind = 4 ) w87b(23)
  real ( kind = 4 ) x1(5)
  real ( kind = 4 ) x2(5)
  real ( kind = 4 ) x3(11)
  real ( kind = 4 ) x4(22)
!
!           the following data statements contain the abscissae
!           and weights of the integration rules used.
!
!           x1      abscissae common to the 10-, 21-, 43- and 87-point
!                   rule
!           x2      abscissae common to the 21-, 43- and 87-point rule
!           x3      abscissae common to the 43- and 87-point rule
!           x4      abscissae of the 87-point rule
!           w10     weights of the 10-point formula
!           w21a    weights of the 21-point formula for abscissae x1
!           w21b    weights of the 21-point formula for abscissae x2
!           w43a    weights of the 43-point formula for absissae x1, x3
!           w43b    weights of the 43-point formula for abscissae x3
!           w87a    weights of the 87-point formula for abscissae x1,
!                   x2 and x3
!           w87b    weights of the 87-point formula for abscissae x4
!
  data x1(1),x1(2),x1(3),x1(4),x1(5)/ &
       9.739065285171717E-01,     8.650633666889845E-01, &
       6.794095682990244E-01,     4.333953941292472E-01, &
       1.488743389816312E-01/
  data x2(1),x2(2),x2(3),x2(4),x2(5)/ &
       9.956571630258081E-01,     9.301574913557082E-01, &
       7.808177265864169E-01,     5.627571346686047E-01, &
       2.943928627014602E-01/
  data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),x3(9),x3(10), &
    x3(11)/ &
       9.993333609019321E-01,     9.874334029080889E-01, &
       9.548079348142663E-01,     9.001486957483283E-01, &
       8.251983149831142E-01,     7.321483889893050E-01, &
       6.228479705377252E-01,     4.994795740710565E-01, &
       3.649016613465808E-01,     2.222549197766013E-01, &
       7.465061746138332E-02/
  data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),x4(10), &
    x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),x4(19), &
    x4(20),x4(21),x4(22)/         9.999029772627292E-01, &
       9.979898959866787E-01,     9.921754978606872E-01, &
       9.813581635727128E-01,     9.650576238583846E-01, &
       9.431676131336706E-01,     9.158064146855072E-01, &
       8.832216577713165E-01,     8.457107484624157E-01, &
       8.035576580352310E-01,     7.570057306854956E-01, &
       7.062732097873218E-01,     6.515894665011779E-01, &
       5.932233740579611E-01,     5.314936059708319E-01, &
       4.667636230420228E-01,     3.994248478592188E-01, &
       3.298748771061883E-01,     2.585035592021616E-01, &
       1.856953965683467E-01,     1.118422131799075E-01, &
       3.735212339461987E-02/
  data w10(1),w10(2),w10(3),w10(4),w10(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
  data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/ &
       3.255816230796473E-02,     7.503967481091995E-02, &
       1.093871588022976E-01,     1.347092173114733E-01, &
       1.477391049013385E-01/
  data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/ &
       1.169463886737187E-02,     5.475589657435200E-02, &
       9.312545458369761E-02,     1.234919762620659E-01, &
       1.427759385770601E-01,     1.494455540029169E-01/
  data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7), &
    w43a(8),w43a(9),w43a(10)/     1.629673428966656E-02, &
       3.752287612086950E-02,     5.469490205825544E-02, &
       6.735541460947809E-02,     7.387019963239395E-02, &
       5.768556059769796E-03,     2.737189059324884E-02, &
       4.656082691042883E-02,     6.174499520144256E-02, &
       7.138726726869340E-02/
  data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),w43b(7), &
    w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/ &
       1.844477640212414E-03,     1.079868958589165E-02, &
       2.189536386779543E-02,     3.259746397534569E-02, &
       4.216313793519181E-02,     5.074193960018458E-02, &
       5.837939554261925E-02,     6.474640495144589E-02, &
       6.956619791235648E-02,     7.282444147183321E-02, &
       7.450775101417512E-02,     7.472214751740301E-02/
  data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),w87a(7), &
    w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),w87a(13),w87a(14), &
    w87a(15),w87a(16),w87a(17),w87a(18),w87a(19),w87a(20),w87a(21)/ &
       8.148377384149173E-03,     1.876143820156282E-02, &
       2.734745105005229E-02,     3.367770731163793E-02, &
       3.693509982042791E-02,     2.884872430211531E-03, &
       1.368594602271270E-02,     2.328041350288831E-02, &
       3.087249761171336E-02,     3.569363363941877E-02, &
       9.152833452022414E-04,     5.399280219300471E-03, &
       1.094767960111893E-02,     1.629873169678734E-02, &
       2.108156888920384E-02,     2.537096976925383E-02, &
       2.918969775647575E-02,     3.237320246720279E-02, &
       3.478309895036514E-02,     3.641222073135179E-02, &
       3.725387550304771E-02/
  data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7), &
    w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14), &
    w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),w87b(21), &
    w87b(22),w87b(23)/            2.741455637620724E-04, &
       1.807124155057943E-03,     4.096869282759165E-03, &
       6.758290051847379E-03,     9.549957672201647E-03, &
       1.232944765224485E-02,     1.501044734638895E-02, &
       1.754896798624319E-02,     1.993803778644089E-02, &
       2.219493596101229E-02,     2.433914712600081E-02, &
       2.637450541483921E-02,     2.828691078877120E-02, &
       3.005258112809270E-02,     3.164675137143993E-02, &
       3.305041341997850E-02,     3.425509970422606E-02, &
       3.526241266015668E-02,     3.607698962288870E-02, &
       3.669860449845609E-02,     3.712054926983258E-02, &
       3.733422875193504E-02,     3.736107376267902E-02/
!
!  Test on validity of parameters.
!
  result = 0.0E+00
  abserr = 0.0E+00
  neval = 0

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if

  hlgth = 5.0E-01 * ( b - a )
  dhlgth = abs ( hlgth )
  centr = 5.0E-01 * ( b + a )
  fcentr = f(centr)
  neval = 21
  ier = 1
!
!  Compute the integral using the 10- and 21-point formula.
!
  do l = 1, 3

    if ( l == 1 ) then

      res10 = 0.0E+00
      res21 = w21b(6) * fcentr
      resabs = w21b(6) * abs(fcentr)

      do k = 1, 5
        absc = hlgth * x1(k)
        fval1 = f(centr+absc)
        fval2 = f(centr-absc)
        fval = fval1 + fval2
        res10 = res10 + w10(k)*fval
        res21 = res21 + w21a(k)*fval
        resabs = resabs + w21a(k)*(abs(fval1)+abs(fval2))
        savfun(k) = fval
        fv1(k) = fval1
        fv2(k) = fval2
      end do

      ipx = 5

      do k = 1, 5
        ipx = ipx + 1
        absc = hlgth * x2(k)
        fval1 = f(centr+absc)
        fval2 = f(centr-absc)
        fval = fval1 + fval2
        res21 = res21 + w21b(k) * fval
        resabs = resabs + w21b(k) * ( abs ( fval1 ) + abs ( fval2 ) )
        savfun(ipx) = fval
        fv3(k) = fval1
        fv4(k) = fval2
      end do
!
!  Test for convergence.
!
      result = res21 * hlgth
      resabs = resabs * dhlgth
      reskh = 5.0E-01 * res21
      resasc = w21b(6) * abs ( fcentr - reskh )

      do k = 1, 5
        resasc = resasc+w21a(k)*(abs(fv1(k)-reskh)+abs(fv2(k)-reskh)) &
                     +w21b(k)*(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
      end do

      abserr = abs ( ( res21 - res10 ) * hlgth )
      resasc = resasc * dhlgth
!
!  Compute the integral using the 43-point formula.
!
    else if ( l == 2 ) then

      res43 = w43b(12)*fcentr
      neval = 43

      do k = 1, 10
        res43 = res43 + savfun(k) * w43a(k)
      end do

      do k = 1, 11
        ipx = ipx + 1
        absc = hlgth * x3(k)
        fval = f(absc+centr) + f(centr-absc)
        res43 = res43 + fval * w43b(k)
        savfun(ipx) = fval
      end do
!
!  Test for convergence.
!
      result = res43 * hlgth
      abserr = abs((res43-res21)*hlgth)
!
!  Compute the integral using the 87-point formula.
!
    else if ( l == 3 ) then

      res87 = w87b(23) * fcentr
      neval = 87

      do k = 1, 21
        res87 = res87 + savfun(k) * w87a(k)
      end do

      do k = 1, 22
        absc = hlgth * x4(k)
        res87 = res87 + w87b(k) * ( f(absc+centr) + f(centr-absc) )
      end do

      result = res87 * hlgth
      abserr = abs ( ( res87 - res43) * hlgth )

    end if

    if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00 ) then
      abserr = resasc * min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
    end if

    if ( resabs > tiny ( resabs ) / ( 5.0E+01 * epsilon ( resabs ) ) ) then
      abserr = max (( epsilon ( resabs ) *5.0E+01) * resabs, abserr )
    end if

    if ( abserr <= max ( epsabs, epsrel*abs(result))) then
      ier = 0
    end if

    if ( ier == 0 ) then
      exit
    end if

  end do

  return
end subroutine
!------------------------------------------------------------------------------
subroutine qags ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGS estimates the integral of a function.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) A, B, the limits of integration.
!
!    Input, real ( kind = 4 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 4 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500

  real ( kind = 4 ) a
  real ( kind = 4 ) abseps
  real ( kind = 4 ) abserr
  real ( kind = 4 ) alist(limit)
  real ( kind = 4 ) area
  real ( kind = 4 ) area1
  real ( kind = 4 ) area12
  real ( kind = 4 ) area2
  real ( kind = 4 ) a1
  real ( kind = 4 ) a2
  real ( kind = 4 ) b
  real ( kind = 4 ) blist(limit)
  real ( kind = 4 ) b1
  real ( kind = 4 ) b2
  real ( kind = 4 ) correc
  real ( kind = 4 ) defabs
  real ( kind = 4 ) defab1
  real ( kind = 4 ) defab2
  real ( kind = 4 ) dres
  real ( kind = 4 ) elist(limit)
  real ( kind = 4 ) epsabs
  real ( kind = 4 ) epsrel
  real ( kind = 4 ) erlarg
  real ( kind = 4 ) erlast
  real ( kind = 4 ) errbnd
  real ( kind = 4 ) errmax
  real ( kind = 4 ) error1
  real ( kind = 4 ) error2
  real ( kind = 4 ) erro12
  real ( kind = 4 ) errsum
  real ( kind = 4 ) ertest
  logical extrap
  real ( kind = 4 ), external :: f
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierro
  integer ( kind = 4 ) iord(limit)
  integer ( kind = 4 ) iroff1
  integer ( kind = 4 ) iroff2
  integer ( kind = 4 ) iroff3
  integer ( kind = 4 ) jupbnd
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ksgn
  integer ( kind = 4 ) ktmin
  integer ( kind = 4 ) last
  logical noext
  integer ( kind = 4 ) maxerr
  integer ( kind = 4 ) neval
  integer ( kind = 4 ) nres
  integer ( kind = 4 ) nrmax
  integer ( kind = 4 ) numrl2
  real ( kind = 4 ) resabs
  real ( kind = 4 ) reseps
  real ( kind = 4 ) result
  real ( kind = 4 ) res3la(3)
  real ( kind = 4 ) rlist(limit)
  real ( kind = 4 ) rlist2(52)
  real ( kind = 4 ) small
!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  ierro = 0
  call qk21 ( f, a, b, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1

  if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 140
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( dres >= (1.0E+00 - 5.0E+01* epsilon ( defabs ) ) * defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
    call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) go to 15

    if ( abs ( rlist(maxerr) - area12) > 1.0E-05 * abs(area12) &
      .or. erro12 < 9.9E-01 * errmax ) go to 10

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( last > 10 .and. erro12 > errmax ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel*abs(area) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
      ier = 2
    end if

    if ( iroff2 >= 5 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
      (abs(a2)+1.0E+03* tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) go to 80
    if ( noext ) go to 90

    erlarg = erlarg-erlast

    if ( abs(b1-a1) > small ) then
      erlarg = erlarg+erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
      extrap = .true.
      nrmax = 2
    end if

!40  continue
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last

      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
          go to 90
        end if
        nrmax = nrmax+1
      end do

    end if
!
!  Perform extrapolation.
!
!60  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5 .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs,epsrel*abs(reseps))

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum
    go to 90

80  continue

    small = abs ( b - a ) * 3.75E-01
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 105
  end if

  if ( abserr > errsum ) go to 115
  if ( area == 0.0E+00 ) go to 130
  go to 110

105 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
!
!  Test on divergence.
!
110 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 130

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum > abs(area) ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

130 continue
 
  if ( 2 < ier ) then
    ier = ier - 1
  end if

140 continue

  neval = 42*last-21

  return
end subroutine
!------------------------------------------------------------------------------
subroutine qagi ( f, bound, inf, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGI estimates an integral over a semi-infinite or infinite interval.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral   
!      I = integral of F over (A, +Infinity), 
!    or 
!      I = integral of F over (-Infinity,A)
!    or 
!      I = integral of F over (-Infinity,+Infinity),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) BOUND, the value of the finite endpoint of the integration
!    range, if any, that is, if INF is 1 or -1.
!
!    Input, integer ( kind = 4 ) INF, indicates the type of integration range.
!    1:  (  BOUND,    +Infinity),
!    -1: ( -Infinity,  BOUND),
!    2:  ( -Infinity, +Infinity).
!
!    Input, real ( kind = 4 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 4 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 4 ) IER, error indicator.
!    0, normal and reliable termination of the routine.  It is assumed that 
!      the requested accuracy has been achieved.
!    > 0,  abnormal termination of the routine.  The estimates for result
!      and error are less reliable.  It is assumed that the requested
!      accuracy has not been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can 
!      allow more subdivisions by increasing the data value of LIMIT in QAGI
!      (and taking the according dimension adjustments into account).
!      However, if this yields no improvement it is advised to analyze the
!      integrand in order to determine the integration difficulties.  If the
!      position of a local difficulty can be determined (e.g. singularity,
!      discontinuity within the interval) one will probably gain from
!      splitting up the interval at this point and calling the integrator 
!      on the subranges.  If possible, an appropriate special-purpose 
!      integrator should be used, which is designed for handling the type
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.  The error may be
!      under-estimated.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    4, the algorithm does not converge.  Roundoff error is detected in the
!      extrapolation table.  It is assumed that the requested tolerance
!      cannot be achieved, and that the returned result is the best which 
!      can be obtained.
!    5, the integral is probably divergent, or slowly convergent.  It must 
!      be noted that divergence can occur with any other value of IER.
!    6, the input is invalid, because INF /= 1 and INF /= -1 and INF /= 2, or
!      epsabs < 0 and epsrel < 0.  result, abserr, neval are set to zero.
!
!  Local parameters:
!
!            the dimension of rlist2 is determined by the value of
!            limexp in QEXTR.
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least (limexp+2),
!                       containing the part of the epsilon table
!                       which is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained, it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered up
!                       to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation. i.e.
!                       before subdividing the smallest interval we
!                       try to decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true-value)
!
  implicit none

  integer ( kind = 4 ), parameter :: limit = 500

  real ( kind = 4 ) abseps
  real ( kind = 4 ) abserr
  real ( kind = 4 ) alist(limit)
  real ( kind = 4 ) area
  real ( kind = 4 ) area1
  real ( kind = 4 ) area12
  real ( kind = 4 ) area2
  real ( kind = 4 ) a1
  real ( kind = 4 ) a2
  real ( kind = 4 ) blist(limit)
  real ( kind = 4 ) boun
  real ( kind = 4 ) bound
  real ( kind = 4 ) b1
  real ( kind = 4 ) b2
  real ( kind = 4 ) correc
  real ( kind = 4 ) defabs
  real ( kind = 4 ) defab1
  real ( kind = 4 ) defab2
  real ( kind = 4 ) dres
  real ( kind = 4 ) elist(limit)
  real ( kind = 4 ) epsabs
  real ( kind = 4 ) epsrel
  real ( kind = 4 ) erlarg
  real ( kind = 4 ) erlast
  real ( kind = 4 ) errbnd
  real ( kind = 4 ) errmax
  real ( kind = 4 ) error1
  real ( kind = 4 ) error2
  real ( kind = 4 ) erro12
  real ( kind = 4 ) errsum
  real ( kind = 4 ) ertest
  logical extrap
  real ( kind = 4 ), external :: f
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierro
  integer ( kind = 4 ) inf
  integer ( kind = 4 ) iord(limit)
  integer ( kind = 4 ) iroff1
  integer ( kind = 4 ) iroff2
  integer ( kind = 4 ) iroff3
  integer ( kind = 4 ) jupbnd
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ksgn
  integer ( kind = 4 ) ktmin
  integer ( kind = 4 ) last
  integer ( kind = 4 ) maxerr
  integer ( kind = 4 ) neval
  logical noext
  integer ( kind = 4 ) nres
  integer ( kind = 4 ) nrmax
  integer ( kind = 4 ) numrl2
  real ( kind = 4 ) resabs
  real ( kind = 4 ) reseps
  real ( kind = 4 ) result
  real ( kind = 4 ) res3la(3)
  real ( kind = 4 ) rlist(limit)
  real ( kind = 4 ) rlist2(52)
  real ( kind = 4 ) small
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = 0.0E+00
  blist(1) = 1.0E+00
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
!  Determine the interval to be mapped onto (0,1).
!  If INF = 2 the integral is computed as i = i1+i2, where
!  i1 = integral of f over (-infinity,0),
!  i2 = integral of f over (0,+infinity).
!
  if ( inf == 2 ) then
    boun = 0.0E+00
  else
    boun = bound
  end if

  call qk15i ( f, boun, inf, 0.0E+00, 1.0E+00, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )

  if ( abserr <= 100.0E+00 * epsilon ( defabs ) * defabs .and. &
    errbnd < abserr ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 130
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  ktmin = 0
  numrl2 = 2
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( ( 1.0E+00 - 5.0E+01 * epsilon ( defabs ) ) * defabs <= dres ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk15i ( f, boun, inf, a1, b1, area1, error1, resabs, defab1 )
    call qk15i ( f, boun, inf, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 /= error1 .and. defab2 /= error2 ) then

      if ( abs ( rlist(maxerr) - area12 ) <= 1.0E-05 * abs ( area12 ) &
        .and. 9.9E-01 * errmax <= erro12 ) then

        if ( extrap ) then
          iroff2 = iroff2 + 1
        end if

        if ( .not. extrap ) then
          iroff1 = iroff1 + 1
        end if

      end if

      if ( 10 < last .and. errmax < erro12 ) then
        iroff3 = iroff3 + 1
      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
      ier = 2
    end if

    if ( 5 <= iroff2 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals equals LIMIT.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at some points of the integration range.
!
    if ( max ( abs(a1), abs(b2) ) <= (1.0E+00 + 1.0E+03 * epsilon ( a1 ) ) * &
    ( abs(a2) + 1.0E+03 * tiny ( a2 ) )) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with NRMAX-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) then
      small = 3.75E-01
      erlarg = errsum
      ertest = errbnd
      rlist2(2) = area
      cycle
    end if

    if ( noext ) then
      cycle
    end if

    erlarg = erlarg - erlast

    if ( small < abs ( b1 - a1 ) ) then
      erlarg = erlarg + erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then

      if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
        cycle
      end if

      extrap = .true.
      nrmax = 2

    end if

    if ( ierro == 3 .or. erlarg <= ertest ) then
      go to 60
    end if
!
!  The smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    id = nrmax
    jupbnd = last

    if ( (2+limit/2) < last ) then
      jupbnd = limit + 3 - last
    end if

    do k = id, jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
        go to 90
      end if
      nrmax = nrmax + 1
    end do
!
!  Extrapolate.
!
60  continue

    numrl2 = numrl2 + 1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres ) 
    ktmin = ktmin+1

    if ( 5 < ktmin .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs, epsrel * abs(reseps) )

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ( ier + ierro ) == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00) then
    go to 105
  end if

  if ( errsum < abserr ) then
    go to 115
  end if

  if ( area == 0.0E+00 ) then
    go to 130
  end if

  go to 110

105 continue

  if ( errsum / abs ( area ) < abserr / abs ( result )  ) then
    go to 115
  end if
!
!  Test on divergence
!
110 continue

  if ( ksgn == (-1) .and. &
  max ( abs(result), abs(area) ) <=  defabs * 1.0E-02) go to 130

  if ( 1.0E-02 > (result/area) .or. &
    (result/area) > 1.0E+02 .or. &
    errsum > abs(area)) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
  115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum
  130 continue

  neval = 30 * last - 15
  if ( inf == 2 ) then
    neval = 2 * neval
  end if

  if ( 2 < ier ) then
    ier = ier - 1
  end if

  return
end subroutine
!---------------------------------------------
subroutine qk15i ( f, boun, inf, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
!
!  Discussion:
!
!    The original infinite integration range is mapped onto the interval 
!    (0,1) and (a,b) is a part of (0,1).  The routine then computes:
!
!    i = integral of transformed integrand over (a,b),
!    j = integral of abs(transformed integrand) over (a,b).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) BOUN, the finite bound of the original integration range,
!    or zero if INF is 2.
!
!    Input, integer ( kind = 4 ) INF, indicates the type of the interval.
!    -1: the original interval is (-infinity,BOUN),
!    +1, the original interval is (BOUN,+infinity),
!    +2, the original interval is (-infinity,+infinity) and
!    the integral is computed as the sum of two integrals, one 
!    over (-infinity,0) and one over (0,+infinity).
!
!    Input, real ( kind = 4 ) A, B, the limits of integration, over a subrange of [0,1].
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained 
!    by optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 4 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 4 ) RESASC, approximation to the integral of the
!    transformated integrand | F-I/(B-A) | over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           tabsc* - transformed abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of the transformed
!                    integrand over (a,b), i.e. to i/(b-a)
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) absc
  real ( kind = 4 ) absc1
  real ( kind = 4 ) absc2
  real ( kind = 4 ) abserr
  real ( kind = 4 ) b
  real ( kind = 4 ) boun
  real ( kind = 4 ) centr
  real ( kind = 4 ) dinf
  real ( kind = 4 ), external :: f
  real ( kind = 4 ) fc
  real ( kind = 4 ) fsum
  real ( kind = 4 ) fval1
  real ( kind = 4 ) fval2
  real ( kind = 4 ) fv1(7)
  real ( kind = 4 ) fv2(7)
  real ( kind = 4 ) hlgth
  integer ( kind = 4 ) inf
  integer ( kind = 4 ) j
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) resg
  real ( kind = 4 ) resk
  real ( kind = 4 ) reskh
  real ( kind = 4 ) result
  real ( kind = 4 ) tabsc1
  real ( kind = 4 ) tabsc2
  real ( kind = 4 ) wg(8)
  real ( kind = 4 ) wgk(8)
  real ( kind = 4 ) xgk(8)
!
!  the abscissae and weights are supplied for the interval
!  (-1,1).  because of symmetry only the positive abscissae and
!  their corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                    rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule, corresponding
!                    to the abscissae xgk(2), xgk(4), ...
!                    wg(1), wg(3), ... are set to zero.
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
       9.914553711208126E-01,     9.491079123427585E-01, &
       8.648644233597691E-01,     7.415311855993944E-01, &
       5.860872354676911E-01,     4.058451513773972E-01, &
       2.077849550078985E-01,     0.0000000000000000E+00/

  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
       2.293532201052922E-02,     6.309209262997855E-02, &
       1.047900103222502E-01,     1.406532597155259E-01, &
       1.690047266392679E-01,     1.903505780647854E-01, &
       2.044329400752989E-01,     2.094821410847278E-01/

  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
       0.0000000000000000E+00,     1.294849661688697E-01, &
       0.0000000000000000E+00,     2.797053914892767E-01, &
       0.0000000000000000E+00,     3.818300505051189E-01, &
       0.0000000000000000E+00,     4.179591836734694E-01/

  dinf = min ( 1, inf )

  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  tabsc1 = boun+dinf*(1.0E+00-centr)/centr
  fval1 = f(tabsc1)
  if ( inf == 2 ) fval1 = fval1+f(-tabsc1)
  fc = (fval1/centr)/centr
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the error.
!
  resg = wg(8)*fc
  resk = wgk(8)*fc
  resabs = abs(resk)

  do j = 1, 7

    absc = hlgth*xgk(j)
    absc1 = centr-absc
    absc2 = centr+absc
    tabsc1 = boun+dinf*(1.0E+00-absc1)/absc1
    tabsc2 = boun+dinf*(1.0E+00-absc2)/absc2
    fval1 = f(tabsc1)
    fval2 = f(tabsc2)

    if ( inf == 2 ) then
      fval1 = fval1+f(-tabsc1)
      fval2 = fval2+f(-tabsc2)
    end if

    fval1 = (fval1/absc1)/absc1
    fval2 = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(j)*fsum
    resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0E-01
  resasc = wgk(8) * abs(fc-reskh)

  do j = 1, 7
    resasc = resasc + wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk * hlgth
  resasc = resasc * hlgth
  resabs = resabs * hlgth
  abserr = abs ( ( resk - resg ) * hlgth )

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc* min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) / ( 5.0E+01 * epsilon ( resabs ) ) ) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
!--------------------------------------------------------------------
subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! QSORT maintains the order of a list of local error estimates.
!
!  Discussion:
!
!    This routine maintains the descending ordering in the list of the 
!    local error estimates resulting from the interval subdivision process. 
!    At each call two error estimates are inserted using the sequential 
!    search top-down for the largest error estimate and bottom-up for the
!    smallest error estimate.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LIMIT, the maximum number of error estimates the list can
!    contain.
!
!    Input, integer ( kind = 4 ) LAST, the current number of error estimates.
!
!    Input/output, integer ( kind = 4 ) MAXERR, the index in the list of the NRMAX-th 
!    largest error.
!
!    Output, real ( kind = 4 ) ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
!
!    Input, real ( kind = 4 ) ELIST(LIMIT), contains the error estimates.
!
!    Input/output, integer ( kind = 4 ) IORD(LAST).  The first K elements contain 
!    pointers to the error estimates such that ELIST(IORD(1)) through
!    ELIST(IORD(K)) form a decreasing sequence, with
!      K = LAST 
!    if 
!      LAST <= (LIMIT/2+2), 
!    and otherwise
!      K = LIMIT+1-LAST.
!
!    Input/output, integer ( kind = 4 ) NRMAX.
!
  implicit none

  integer ( kind = 4 ) last

  real ( kind = 4 ) elist(last)
  real ( kind = 4 ) ermax
  real ( kind = 4 ) errmax
  real ( kind = 4 ) errmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) iord(last)
  integer ( kind = 4 ) isucc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jbnd
  integer ( kind = 4 ) jupbn
  integer ( kind = 4 ) k
  integer ( kind = 4 ) limit
  integer ( kind = 4 ) maxerr
  integer ( kind = 4 ) nrmax
!
!  Check whether the list contains more than two error estimates.
!
  if ( last <= 2 ) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
!
!  This part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
  errmax = elist(maxerr)

  do i = 1, nrmax-1

    isucc = iord(nrmax-1)

    if ( errmax <= elist(isucc) ) then
      exit
    end if

    iord(nrmax) = isucc
    nrmax = nrmax-1

  end do
!
!  Compute the number of elements in the list to be maintained
!  in descending order.  This number depends on the number of
!  subdivisions still allowed.
!
  jupbn = last

  if ( (limit/2+2) < last ) then
    jupbn = limit+3-last
  end if

  errmin = elist(last)
!
!  Insert errmax by traversing the list top-down, starting
!  comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i = ibeg, jbnd
    isucc = iord(i)
    if ( elist(isucc) <= errmax ) then
      go to 60
    end if
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  Insert errmin by traversing the list bottom-up.
!
60 continue

  iord(i-1) = maxerr
  k = jbnd

  do j = i, jbnd
    isucc = iord(k)
    if ( errmin < elist(isucc) ) then
      go to 80
    end if
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last
!
!  Set maxerr and ermax.
!
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end subroutine
!-----------------------------------------------------
subroutine qextr ( n, epstab, result, abserr, res3la, nres )

!*****************************************************************************80
!
!! QEXTR carries out the Epsilon extrapolation algorithm.
!
!  Discussion:
!
!    The routine determines the limit of a given sequence of approximations, 
!    by means of the epsilon algorithm of P. Wynn.  An estimate of the 
!    absolute error is also given.  The condensed epsilon table is computed.
!    Only those elements needed for the computation of the next diagonal
!    are preserved.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, indicates the entry of EPSTAB which contains
!    the new element in the first column of the epsilon table.
!
!    Input/output, real ( kind = 4 ) EPSTAB(52), the two lower diagonals of the triangular
!    epsilon table.  The elements are numbered starting at the right-hand 
!    corner of the triangle.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 4 ) ABSERR, estimate of the absolute error computed from
!    RESULT and the 3 previous results.
!
!    ?, real ( kind = 4 ) RES3LA(3), the last 3 results.
!
!    Input/output, integer ( kind = 4 ) NRES, the number of calls to the routine.  This
!    should be zero on the first call, and is automatically updated
!    before return.
!
!  Local Parameters:
!
!           e0     - the 4 elements on which the
!           e1       computation of a new element in
!           e2       the epsilon table is based
!           e3                 e0
!                        e3    e1    new
!                              e2
!           newelm - number of elements to be computed in the new
!                    diagonal
!           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!           result - the element in the new diagonal with least value
!                    of error
!           limexp is the maximum number of elements the epsilon table
!           can contain. if this number is reached, the upper diagonal
!           of the epsilon table is deleted.
!
  implicit none

  real ( kind = 4 ) abserr
  real ( kind = 4 ) delta1
  real ( kind = 4 ) delta2
  real ( kind = 4 ) delta3
  real ( kind = 4 ) epsinf
  real ( kind = 4 ) epstab(52)
  real ( kind = 4 ) error
  real ( kind = 4 ) err1
  real ( kind = 4 ) err2
  real ( kind = 4 ) err3
  real ( kind = 4 ) e0
  real ( kind = 4 ) e1
  real ( kind = 4 ) e1abs
  real ( kind = 4 ) e2
  real ( kind = 4 ) e3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ib2
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) limexp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) newelm
  integer ( kind = 4 ) nres
  integer ( kind = 4 ) num
  real ( kind = 4 ) res
  real ( kind = 4 ) result
  real ( kind = 4 ) res3la(3)
  real ( kind = 4 ) ss
  real ( kind = 4 ) tol1
  real ( kind = 4 ) tol2
  real ( kind = 4 ) tol3

  nres = nres+1
  abserr = huge ( abserr )
  result = epstab(n)

  if ( n < 3 ) then
    abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))
    return
  end if

  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = huge ( epstab(n) )
  num = n
  k1 = n

  do i = 1, newelm

    k2 = k1-1
    k3 = k1-2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    e1abs = abs(e1)
    delta2 = e2-e1
    err2 = abs(delta2)
    tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
    delta3 = e1-e0
    err3 = abs(delta3)
    tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )
!
!  If e0, e1 and e2 are equal to within machine accuracy, convergence 
!  is assumed.
!
    if ( err2 <= tol2 .and. err3 <= tol3 ) then
      result = res
      abserr = err2+err3
      abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))
      return
    end if

    e3 = epstab(k1)
    epstab(k1) = e1
    delta1 = e1-e3
    err1 = abs(delta1)
    tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )
!
!  If two elements are very close to each other, omit a part
!  of the table by adjusting the value of N.
!
    if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

    ss = 1.0E+00/delta1+1.0E+00/delta2-1.0E+00/delta3
    epsinf = abs ( ss*e1 )
!
!  Test to detect irregular behavior in the table, and
!  eventually omit a part of the table adjusting the value of N.
!
    if ( epsinf > 1.0E-04 ) go to 30

20  continue

    n = i+i-1
    exit
!
!  Compute a new element and eventually adjust the value of RESULT.
!
30  continue

    res = e1+1.0E+00/ss
    epstab(k1) = res
    k1 = k1-2
    error = err2+abs(res-e2)+err3

    if ( error <= abserr ) then
      abserr = error
      result = res
    end if

  end do
!
!  Shift the table.
!
  if ( n == limexp ) then
    n = 2*(limexp/2)-1
  end if

  if ( (num/2)*2 == num ) then
    ib = 2
  else
    ib = 1
  end if

  ie = newelm+1

  do i = 1, ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do

  if ( num /= n ) then

    indx = num-n+1

    do i = 1, n
      epstab(i)= epstab(indx)
      indx = indx+1
    end do

  end if

  if ( nres < 4 ) then
    res3la(nres) = result
    abserr = huge ( abserr )
  else
    abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
      +abs(result-res3la(1))
    res3la(1) = res3la(2)
    res3la(2) = res3la(3)
    res3la(3) = result
  end if

  abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))

  return
end subroutine
!----------------------------------------------------------------
subroutine qk21 ( f, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger, 
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 4 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 4 ) f
!      real ( kind = 4 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 4 ) A, B, the limits of integration.
!
!    Output, real ( kind = 4 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk) 
!    obtained by optimal addition of abscissae to the 10-point Gauss 
!    rule (resg).
!
!    Output, real ( kind = 4 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 4 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 4 ) RESASC, approximation to the integral | F-I/(B-A) | 
!    over [A,B].
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) absc
  real ( kind = 4 ) abserr
  real ( kind = 4 ) b
  real ( kind = 4 ) centr
  real ( kind = 4 ) dhlgth
  real ( kind = 4 ), external :: f
  real ( kind = 4 ) fc
  real ( kind = 4 ) fsum
  real ( kind = 4 ) fval1
  real ( kind = 4 ) fval2
  real ( kind = 4 ) fv1(10)
  real ( kind = 4 ) fv2(10)
  real ( kind = 4 ) hlgth
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtw
  integer ( kind = 4 ) jtwm1
  real ( kind = 4 ) resabs
  real ( kind = 4 ) resasc
  real ( kind = 4 ) resg
  real ( kind = 4 ) resk
  real ( kind = 4 ) reskh
  real ( kind = 4 ) result
  real ( kind = 4 ) wg(5)
  real ( kind = 4 ) wgk(11)
  real ( kind = 4 ) xgk(11)
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11)/ &
       9.956571630258081E-01,     9.739065285171717E-01, &
       9.301574913557082E-01,     8.650633666889845E-01, &
       7.808177265864169E-01,     6.794095682990244E-01, &
       5.627571346686047E-01,     4.333953941292472E-01, &
       2.943928627014602E-01,     1.488743389816312E-01, &
       0.000000000000000E+00/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11)/ &
       1.169463886737187E-02,     3.255816230796473E-02, &
       5.475589657435200E-02,     7.503967481091995E-02, &
       9.312545458369761E-02,     1.093871588022976E-01, &
       1.234919762620659E-01,     1.347092173114733E-01, &
       1.427759385770601E-01,     1.477391049013385E-01, &
       1.494455540029169E-01/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0E+00
  fc = f(centr)
  resk = wgk(11)*fc
  resabs = abs(resk)

  do j = 1, 5
    jtw = 2*j
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1, 5
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*5.0E-01
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end subroutine
