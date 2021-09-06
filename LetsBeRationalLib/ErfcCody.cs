using System;

namespace LetsBeRationalLib {
    class ErfcCody {
        internal static double erf_cody(double x) {
            /* -------------------------------------------------------------------- */
            /* This subprogram computes approximate values for erf(x). */
            /*   Author/date: W. J. Cody, January 8, 1985 */
            /* -------------------------------------------------------------------- */
            return calerf(x, 0);
        }

        internal static double erfc_cody(double x) {
            /* -------------------------------------------------------------------- */
            /* This subprogram computes approximate values for erfc(x). */
            /*   (see comments heading CALERF). */
            /*   Author/date: W. J. Cody, January 8, 1985 */
            /* -------------------------------------------------------------------- */
            return calerf(x, 1);
        }

        /* S    REAL FUNCTION ERFCX(X) */
        /*<       DOUBLE PRECISION FUNCTION DERFCX(X) >*/
        internal static double erfcx_cody(double x) {
            /* ------------------------------------------------------------------ */
            /* This subprogram computes approximate values for exp(x*x) * erfc(x). */
            /*   (see comments heading CALERF). */
            /*   Author/date: W. J. Cody, March 30, 1987 */
            /* ------------------------------------------------------------------ */
            return calerf(x, 2);
        }

        private static double d_int(double x) { return (x > 0) ? Math.Floor(x) : -Math.Floor(-x); }

        private static double calerf(double x, int jint) {
            double[] a = { 3.1611237438705656, 113.864154151050156, 377.485237685302021, 3209.37758913846947, 0.185777706184603153 };
            double[] b = { 23.6012909523441209, 244.024637934444173, 1282.61652607737228, 2844.23683343917062 };
            double[] c__ = { 0.564188496988670089, 8.88314979438837594, 66.1191906371416295, 298.635138197400131, 881.95222124176909, 1712.04761263407058, 2051.07837782607147, 1230.33935479799725, 2.15311535474403846e-8 };
            double[] d__ = { 15.7449261107098347, 117.693950891312499, 537.181101862009858, 1621.38957456669019, 3290.79923573345963, 4362.61909014324716, 3439.36767414372164, 1230.33935480374942 };
            double[] p = { 0.305326634961232344, 0.360344899949804439, 0.125781726111229246, 0.0160837851487422766, 6.58749161529837803e-4, 0.0163153871373020978 };
            double[] q = { 2.56852019228982242, 1.87295284992346047, 0.527905102951428412, 0.0605183413124413191, 0.00233520497626869185 };

            const double zero = 0.0;
            const double half = 0.5;
            const double one = 1.0;
            const double two = 2.0;
            const double four = 4.0;
            const double sqrpi = 0.56418958354775628695;
            const double thresh = 0.46875;
            const double sixten = 16.0;

            double y, del, ysq, xden, xnum, result;

            /* ------------------------------------------------------------------ */
            /* This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x) */
            /*   for a real argument  x.  It contains three FUNCTION type */
            /*   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX), */
            /*   and one SUBROUTINE type subprogram, CALERF.  The calling */
            /*   statements for the primary entries are: */
            /*                   Y=ERF(X)     (or   Y=DERF(X)), */
            /*                   Y=ERFC(X)    (or   Y=DERFC(X)), */
            /*   and */
            /*                   Y=ERFCX(X)   (or   Y=DERFCX(X)). */
            /*   The routine  CALERF  is intended for internal packet use only, */
            /*   all computations within the packet being concentrated in this */
            /*   routine.  The function subprograms invoke  CALERF  with the */
            /*   statement */
            /*          CALL CALERF(ARG,RESULT,JINT) */
            /*   where the parameter usage is as follows */
            /*      Function                     Parameters for CALERF */
            /*       call              ARG                  Result          JINT */
            /*     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0 */
            /*     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1 */
            /*     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2 */
            /*   The main computation evaluates near-minimax approximations */
            /*   from "Rational Chebyshev approximations for the error function" */
            /*   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This */
            /*   transportable program uses rational functions that theoretically */
            /*   approximate  erf(x)  and  erfc(x)  to at least 18 significant */
            /*   decimal digits.  The accuracy achieved depends on the arithmetic */
            /*   system, the compiler, the intrinsic functions, and proper */
            /*   selection of the machine-dependent constants. */
            /* ******************************************************************* */
            /* ******************************************************************* */
            /* Explanation of machine-dependent constants */
            /*   XMIN   = the smallest positive floating-point number. */
            /*   XINF   = the largest positive finite floating-point number. */
            /*   XNEG   = the largest negative argument acceptable to ERFCX; */
            /*            the negative of the solution to the equation */
            /*            2*exp(x*x) = XINF. */
            /*   XSMALL = argument below which erf(x) may be represented by */
            /*            2*x/sqrt(pi)  and above which  x*x  will not underflow. */
            /*            A conservative value is the largest machine number X */
            /*            such that   1.0 + X = 1.0   to machine precision. */
            /*   XBIG   = largest argument acceptable to ERFC;  solution to */
            /*            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where */
            /*            W(x) = exp(-x*x)/[x*sqrt(pi)]. */
            /*   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to */
            /*            machine precision.  A conservative value is */
            /*            1/[2*sqrt(XSMALL)] */
            /*   XMAX   = largest acceptable argument to ERFCX; the minimum */
            /*            of XINF and 1/[sqrt(pi)*XMIN]. */
            // The numbers below were preselected for IEEE .
            const double xinf = 1.79e308;
            const double xneg = -26.628;
            const double xsmall = 1.11e-16;
            const double xbig = 26.543;
            const double xhuge = 6.71e7;
            const double xmax = 2.53e307;
            /*   Approximate values for some important machines are: */
            /*                          XMIN       XINF        XNEG     XSMALL */
            /*  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15 */
            /*  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15 */
            /*  IEEE (IBM/XT, */
            /*    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8 */
            /*  IEEE (IBM/XT, */
            /*    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16 */
            /*  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17 */
            /*  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18 */
            /*  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17 */
            /*  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16 */
            /*                          XBIG       XHUGE       XMAX */
            /*  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293 */
            /*  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465 */
            /*  IEEE (IBM/XT, */
            /*    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37 */
            /*  IEEE (IBM/XT, */
            /*    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307 */
            /*  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75 */
            /*  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307 */
            /*  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38 */
            /*  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307 */
            /* ******************************************************************* */
            /* ******************************************************************* */
            /* Error returns */
            /*  The program returns  ERFC = 0      for  ARG .GE. XBIG; */
            /*                       ERFCX = XINF  for  ARG .LT. XNEG; */
            /*      and */
            /*                       ERFCX = 0     for  ARG .GE. XMAX. */
            /* Intrinsic functions required are: */
            /*     ABS, AINT, EXP */
            /*  Author: W. J. Cody */
            /*          Mathematics and Computer Science Division */
            /*          Argonne National Laboratory */
            /*          Argonne, IL 60439 */
            /*  Latest modification: March 19, 1990 */
            /* ------------------------------------------------------------------ */
            /*<       INTEGER I,JINT >*/
            /* S    REAL */
            /*<    >*/
            /*<       DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5) >*/
            /* ------------------------------------------------------------------ */
            /*  Mathematical constants */
            /* ------------------------------------------------------------------ */
            /* S    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/, */
            /* S   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/, */
            /* S   2     SIXTEN/16.0E0/ */
            /*<    >*/
            /* ------------------------------------------------------------------ */
            /*  Machine-dependent constants */
            /* ------------------------------------------------------------------ */
            /* S    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/, */
            /* S   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/ */
            /*<    >*/
            /* ------------------------------------------------------------------ */
            /*  Coefficients for approximation to  erf  in first interval */
            /* ------------------------------------------------------------------ */
            /* S    DATA A/3.16112374387056560E00,1.13864154151050156E02, */
            /* S   1       3.77485237685302021E02,3.20937758913846947E03, */
            /* S   2       1.85777706184603153E-1/ */
            /* S    DATA B/2.36012909523441209E01,2.44024637934444173E02, */
            /* S   1       1.28261652607737228E03,2.84423683343917062E03/ */
            /*<    >*/
            /*<    >*/
            /* ------------------------------------------------------------------ */
            /*  Coefficients for approximation to  erfc  in second interval */
            /* ------------------------------------------------------------------ */
            /* S    DATA C/5.64188496988670089E-1,8.88314979438837594E0, */
            /* S   1       6.61191906371416295E01,2.98635138197400131E02, */
            /* S   2       8.81952221241769090E02,1.71204761263407058E03, */
            /* S   3       2.05107837782607147E03,1.23033935479799725E03, */
            /* S   4       2.15311535474403846E-8/ */
            /* S    DATA D/1.57449261107098347E01,1.17693950891312499E02, */
            /* S   1       5.37181101862009858E02,1.62138957456669019E03, */
            /* S   2       3.29079923573345963E03,4.36261909014324716E03, */
            /* S   3       3.43936767414372164E03,1.23033935480374942E03/ */
            /*<    >*/
            /*<    >*/
            /* ------------------------------------------------------------------ */
            /*  Coefficients for approximation to  erfc  in third interval */
            /* ------------------------------------------------------------------ */
            /* S    DATA P/3.05326634961232344E-1,3.60344899949804439E-1, */
            /* S   1       1.25781726111229246E-1,1.60837851487422766E-2, */
            /* S   2       6.58749161529837803E-4,1.63153871373020978E-2/ */
            /* S    DATA Q/2.56852019228982242E00,1.87295284992346047E00, */
            /* S   1       5.27905102951428412E-1,6.05183413124413191E-2, */
            /* S   2       2.33520497626869185E-3/ */
            /*<    >*/
            /*<    >*/
            /* ------------------------------------------------------------------ */

            y = Math.Abs(x);
            if (y <= thresh) {
                /* ------------------------------------------------------------------ */
                /*  Evaluate  erf  for  |X| <= 0.46875 */
                /* ------------------------------------------------------------------ */
                ysq = zero;
                if (y > xsmall)
                    ysq = y* y;

                xnum = a[4] * ysq;
                xden = ysq;
                for (int i__ = 1; i__ <= 3; ++i__) {
                    xnum = (xnum + a[i__ - 1]) * ysq;
                    xden = (xden + b[i__ - 1]) * ysq;
                }

                result = x* (xnum + a[3]) / (xden + b[3]);
                if (jint != 0) 
                    result = one - result;

                if (jint == 2)
                    result = Math.Exp(ysq) * result;

                return result;
            }
            else if (y <= four) {
                /* ------------------------------------------------------------------ */
                /*  Evaluate  erfc  for 0.46875 <= |X| <= 4.0 */
                /* ------------------------------------------------------------------ */

                xnum = c__[8] * y;
                xden = y;

                for (int i__ = 1; i__ <= 7; ++i__) {
                    xnum = (xnum + c__[i__ - 1]) * y;
                    xden = (xden + d__[i__ - 1]) * y;
                }

                result = (xnum + c__[7]) / (xden + d__[7]);

                if (jint != 2) {
                    double d__1 = y * sixten;
                    ysq = d_int(d__1) / sixten;
                    del = (y - ysq) * (y + ysq);
                    d__1 = Math.Exp(-ysq* ysq) * Math.Exp(-del);
                    result = d__1* result;
                }
            }
            else {
                /* ------------------------------------------------------------------ */
                /*  Evaluate  erfc  for |X| > 4.0 */
                /* ------------------------------------------------------------------ */
                result = zero;
                if (y >= xbig) {
                    if (jint != 2 || y >= xmax)
                        goto L300;

                    if (y >= xhuge) {
                        result = sqrpi / y;
                        goto L300;
                    }
                }

                ysq = one / (y* y);
                xnum = p[5] * ysq;
                xden = ysq;
                for (int i__ = 1; i__ <= 4; ++i__) {
                    xnum = (xnum + p[i__ - 1]) * ysq;
                    xden = (xden + q[i__ - 1]) * ysq;
                }

                result = ysq* (xnum + p[4]) / (xden + q[4]);
                result = (sqrpi - result) / y;
                if (jint != 2) {
                    double d__1 = y * sixten;
                    ysq = d_int(d__1) / sixten;
                    del = (y - ysq) * (y + ysq);
                    d__1 = Math.Exp(-ysq* ysq) * Math.Exp(-del);
                    result = d__1* result;
                }
            }

L300:
            /* ------------------------------------------------------------------ */
            /*  Fix up for negative argument, erf, etc. */
            /* ------------------------------------------------------------------ */
            if (jint == 0) {
                result = (half - result) + half;
                if (x < zero)
                    result = -(result);
            }
            else if (jint == 1) {
                if (x < zero)
                    result = two - result;
            }
            else {
                if (x < zero) {
                    if (x < xneg) 
                        result = xinf;
                    else {
                        double d__1 = x * sixten;
                        ysq = d_int(d__1) / sixten;
                        del = (x - ysq) * (x + ysq);
                        y = Math.Exp(ysq* ysq) * Math.Exp(del);
                        result = y + y - result;
                    }
                }
            }

            return result;
        }
    }
}
