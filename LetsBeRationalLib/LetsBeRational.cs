﻿using System;

namespace LetsBeRationalLib {
    public static class LetsBeRational {
        public enum OptionType { Put, Call };
        public const double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC = -Double.MaxValue;
        public const double VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM = Double.MaxValue;
        private const int implied_volatility_maximum_iterations = 2;

        private const double TWO_PI = 6.283_185_307_179_586_476;
        private const double SQRT_PI_OVER_TWO = 1.253_314_137_315_500_251;  // sqrt(pi/2) to avoid misinterpretation.
        private const double SQRT_THREE = 1.732_050_807_568_877_294;
        private const double SQRT_ONE_OVER_THREE = 0.577_350_269_189_625_765;
        private const double TWO_PI_OVER_SQRT_TWENTY_SEVEN = 1.209_199_576_156_145_234; // 2*pi/sqrt(27)
        private const double PI_OVER_SIX = 0.523_598_775_598_298_873;
        private const double ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887;
        private const double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;
        private const double SQRT_TWO_PI = 2.506628274631000502415765284811045253006986740610;
        private const double asymptotic_expansion_accuracy_threshold = -10;

        private static readonly double SQRT_DBL_EPSILON = Math.Sqrt(Double.Epsilon);
        private static readonly double FOURTH_ROOT_DBL_EPSILON = Math.Sqrt(SQRT_DBL_EPSILON);
        private static readonly double EIGHTH_ROOT_DBL_EPSILON = Math.Sqrt(FOURTH_ROOT_DBL_EPSILON);
        private static readonly double SIXTEENTH_ROOT_DBL_EPSILON = Math.Sqrt(EIGHTH_ROOT_DBL_EPSILON);
        private static readonly double small_t_expansion_of_normalised_black_threshold = 2 * SIXTEENTH_ROOT_DBL_EPSILON;
        private static readonly double SQRT_DBL_MIN = Math.Sqrt(Double.MinValue);
        private static readonly double SQRT_DBL_MAX = Math.Sqrt(Double.MaxValue);

        // Set this to 0 if you want positive results for (positive) denormalised inputs, else to Double.MinValue.
        // Note that you cannot achieve full machine accuracy from denormalised inputs!
        private const double DENORMALISATION_CUTOFF = 0;

        // s = underlyiong price
        // K = strike
        // t = dte/365
        // r = annualized risk free rate
        // q = annualized sp500 dividen yield (trailing 12 months)
        public static double ImpliedVolatility(double option_price, double s, double K, double t, double r, double q, OptionType optionType) {
            double deflater = Math.Exp(-r * t);
            double undiscounted_option_price = option_price / deflater;
            double F = s * Math.Exp((r - q) * t);
            int intOptionType = (optionType == OptionType.Put ? -1 : 1);
            double iv = ImpliedVolatilityFromATransformedRationalGuess(undiscounted_option_price, F, K, t, intOptionType);
            return iv;
        }

        // s: underlying asset price
        // K: strike price
        // sigma: annualized standard deviation, or volatility
        // t: time to expiration in years
        // r: risk-free interest rate
        // q: annualized continuous dividend rate
        public static double Delta(double s, double K, double t, double r, double sigma, double q, OptionType optionType) {
            double d1 = D1(s, K, t, r, sigma, q);

            if (optionType == OptionType.Put)
                return -Math.Exp(-q * t) * NormalDistribution.norm_cdf(-d1);
            else
                return Math.Exp(-q * t) * NormalDistribution.norm_cdf(d1);
        }

        /// <summary>
        /// Computes European option theta.
        /// <param name="optionType">call or put</param>
        /// <param name="s">Underlying price</param>
        /// <param name="K">Strike price</param>
        /// <param name="t">Time to expiration in fraction of year</param>
        /// <param name="sigma">Volatility - usually computed implied volatility</param>
        /// <param name="r">continuously compounded risk-free annual interest rate</param>
        /// <param name="q">continuously compounded annual dividend yield</param>
        /// <returns></returns>
        public static double Theta(double s, double K, double t, double r, double sigma, double q, OptionType optionType)
        {
            if (t <= 0.0)
                return 0.0;
            double sqrt_t = Math.Sqrt(t);
            double e = Math.Exp(-q * t);
            double d1 = D1(s, K, t, r, sigma, q);
            double d2 = d1 - sigma * sqrt_t;
            double f1 = -e * s * NormalDistribution.norm_pdf(d1) * sigma / (2 * sqrt_t);
            double f2;
            if (optionType == OptionType.Call)
                f2 = -r * K * Math.Exp(-r * t) * NormalDistribution.norm_cdf(d2) + q * s * e * NormalDistribution.norm_cdf(d1);
            else
                f2 = r * K * Math.Exp(-r * t) * NormalDistribution.norm_cdf(-d2) - q * s * e * NormalDistribution.norm_cdf(-d1);
            return (f1 + f2) / 365.0;
        }

        // s: underlying asset price
        // K: strike price
        // sigma: annualized standard deviation, or volatility
        // t: time to expiration in years
        // r: risk-free interest rate
        // d: annualized continuous dividend rate
        public static double Gamma(double s, double K, double t, double r, double sigma, double d, OptionType optionType)
        {
            if (t <= 0.0)
                return 0.0;
            double d1 = D1(s, K, t, r, sigma, d);
            double t1 = Math.Exp(-d * t) / (s * sigma * Math.Sqrt(t));
            double t2 = Math.Exp(-d1 * d1 / 2.0) / SQRT_TWO_PI;
            return t1 * t2;
        }

        // s: underlying asset price
        // K: strike price
        // sigma: annualized standard deviation, or volatility
        // t: time to expiration in years
        // r: risk-free interest rate
        // q: annualized continuous dividend rate
        public static double Vega(double s, double K, double t, double r, double sigma, double q, OptionType optionType)
        {
            if (t <= 0.0)
                return 0.0;
            double d1 = D1(s, K, t, r, sigma, q);
            return s * Math.Exp(-q * t) * Math.Sqrt(t) * NormalDistribution.norm_pdf(d1) / 100.0;
        }

        // s: underlying asset price
        // K: strike price
        // sigma: annualized standard deviation, or volatility
        // t: time to expiration in years
        // r: risk-free interest rate
        // q: annualized continuous dividend rate
        public static double Rho(double s, double K, double t, double r, double sigma, double q, OptionType optionType)
        {
            if (t <= 0.0)
                return 0.0;
            double d1 = D1(s, K, t, r, sigma, q);
            double d2 = d1 - sigma * Math.Sqrt(t);
            double t1 = K * t * Math.Exp(-r * t);
            double t2;
            if (optionType == OptionType.Put)
                t2 = -NormalDistribution.norm_cdf(-d2);
            else
                t2 = NormalDistribution.norm_cdf(d2);
            return t1 * t2;
        }

        // s: underlying asset price
        // K: strike price
        // sigma: annualized standard deviation, or volatility
        // t: time to expiration in years
        // r: risk-free interest rate
        // q: annualized continuous dividend rate
        // From Espen Haug, The Complete Guide To Option Pricing Formulas
        private static double D1(double s, double K, double t, double r, double sigma, double q) {
            double numerator = Math.Log(s / K) + ((r - q) + 0.5 * sigma * sigma) * t;
            double denominator = sigma * Math.Sqrt(t);
            return numerator / denominator;
        }

        const int implied_volatility_householder_method_order = 4;
        private static double householder_factor(double newton, double halley, double hh3) {
            return implied_volatility_householder_method_order > 3 ? (1.0 + 0.5 * halley * newton) / (1.0 + newton * (halley + hh3 * newton / 6)) : (implied_volatility_householder_method_order > 2 ? 1 / (1 + 0.5 * halley * newton) : 1);
        }

        private static bool is_below_horizon(double x) { return Math.Abs(x) < DENORMALISATION_CUTOFF; } // This weeds out denormalised (a.k.a. 'subnormal') numbers.

        private static double ImpliedVolatilityFromATransformedRationalGuess(double price, double F, double K, double T, int optionType) {
            return ImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations(price, F, K, T, optionType, implied_volatility_maximum_iterations);
        }

        private static double ImpliedVolatilityFromATransformedRationalGuessWithLimitedIterations(double price, double F, double K, double T, int optionType, int N) {
            double intrinsic = Math.Abs(Math.Max((optionType < 0 ? K - F : F - K), 0.0));
#if true
            if (price < intrinsic)
            {
                price = intrinsic + 0.0001;
                //return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_BELOW_INTRINSIC;
            }
#endif
            double max_price = optionType < 0 ? K : F;
            if (price >= max_price)
                return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;

            // Map in-the-money to out-of-the-money
            double x = Math.Log(F / K);
            if (optionType * x > 0) {
                price = Math.Abs(Math.Max(price - intrinsic, 0.0));
                optionType = -optionType;
            }

            return unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(price / (Math.Sqrt(F) * Math.Sqrt(K)), x, optionType, N) / Math.Sqrt(T);
        }

        // See http://en.wikipedia.org/wiki/Householder%27s_method for a detailed explanation of the third order Householder iteration.
        //
        // Given the objective function g(s) whose root x such that 0 = g(s) we seek, iterate
        //
        //     s_n+1  =  s_n  -  (g/g') · [ 1 - (g''/g')·(g/g') ] / [ 1 - (g/g')·( (g''/g') - (g'''/g')·(g/g')/6 ) ]
        //
        // Denoting  newton:=-(g/g'), halley:=(g''/g'), and hh3:=(g'''/g'), this reads
        //
        //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ]
        //
        //
        // NOTE that this function returns 0 when beta<intrinsic without any safety checks.
        //
        private static double unchecked_normalised_implied_volatility_from_a_transformed_rational_guess_with_limited_iterations(double beta, double x, int q, int N) {
            // Subtract intrinsic.
            if (q * x > 0) {
                beta = Math.Abs(Math.Max(beta - normalised_intrinsic(x, q), 0.0));
                q = -q;
            }

            // Map puts to calls
            if (q < 0) {
                x = -x;
                q = -q;
            }

            if (beta <= 0) // For negative or zero prices we return 0.
                return 0;
            if (beta < DENORMALISATION_CUTOFF) // For positive but denormalised (a.k.a. 'subnormal') prices, we return 0 since it would be impossible to converge to full machine accuracy anyway.
                return 0;
            double b_max = Math.Exp(0.5 * x);
            if (beta >= b_max)
                return VOLATILITY_VALUE_TO_SIGNAL_PRICE_IS_ABOVE_MAXIMUM;

            int iterations = 0, direction_reversal_count = 0;
            double f = -Double.MaxValue, s = -Double.MaxValue, ds = s, ds_previous = 0, s_left = Double.MinValue, s_right = Double.MaxValue;
            // The temptation is great to use the optimised form b_c = exp(x/2)/2-exp(-x/2)·Phi(sqrt(-2·x)) but that would require implementing all of the above types of round-off and over/underflow handling for this expression, too.
            double s_c = Math.Sqrt(Math.Abs(2 * x));
            double b_c = normalised_black_call(x, s_c);
            double v_c = normalised_vega(x, s_c);

            // Four branches.
            if (beta < b_c) {
                double s_l = s_c - b_c / v_c;
                double b_l = normalised_black_call(x, s_l);
                if (beta < b_l) {
                    double f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2;
                    (f_lower_map_l, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2) = compute_f_lower_map_and_first_two_derivatives(x, s_l);
                    double r_ll = RationalCubic.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, d2_f_lower_map_l_d_beta2, true);
                    f = RationalCubic.rational_cubic_interpolation(beta, 0.0, b_l, 0.0, f_lower_map_l, 1.0, d_f_lower_map_l_d_beta, r_ll);
                    if (!(f > 0)) { // This can happen due to roundoff truncation for extreme values such as |x|>500.
                                    // We switch to quadratic interpolation using f(0)≡0, f(b_l), and f'(0)≡1 to specify the quadratic.
                        double t = beta / b_l;
                        f = (f_lower_map_l * t + b_l * (1 - t)) * t;
                    }
                    s = inverse_f_lower_map(x, f);
                    s_right = s_l;
                    //
                    // In this branch, which comprises the lowest segment, the objective function is
                    //     g(s) = 1/ln(b(x,s)) - 1/ln(beta)
                    //          ≡ 1/ln(b(s)) - 1/ln(beta)
                    // This makes
                    //              g'               =   -b'/(b·ln(b)²)
                    //              newton = -g/g'   =   (ln(beta)-ln(b))·ln(b)/ln(beta)·b/b'
                    //              halley = g''/g'  =   b''/b'  -  b'/b·(1+2/ln(b))
                    //              hh3    = g'''/g' =   b'''/b' +  2(b'/b)²·(1+3/ln(b)·(1+1/ln(b)))  -  3(b''/b)·(1+2/ln(b))
                    //
                    // The Householder(3) iteration is
                    //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ]
                    //
                    for (; iterations < N && Math.Abs(ds) > Double.Epsilon * s; ++iterations) {
                        if (ds * ds_previous < 0)
                            ++direction_reversal_count;

                        if (iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right))) {
                            // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
                            // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
                            s = 0.5 * (s_left + s_right);
                            if (s_right - s_left <= Double.Epsilon * s) break;
                            direction_reversal_count = 0;
                            ds = 0;
                        }

                        ds_previous = ds;
                        double b = normalised_black_call(x, s), bp = normalised_vega(x, s);
                        if (b > beta && s < s_right) s_right = s;
                        else if (b < beta && s > s_left) s_left = s; // Tighten the bracket if applicable.

                        if (b <= 0 || bp <= 0) // Numerical underflow. Switch to binary nesting for this iteration.
                            ds = 0.5 * (s_left + s_right) - s;
                        else {
                            double ln_b = Math.Log(b), ln_beta = Math.Log(beta), bpob = bp / b, h = x / s, b_halley = h * h / s - s / 4, newton = (ln_beta - ln_b) * ln_b / ln_beta / bpob, halley = b_halley - bpob * (1 + 2 / ln_b);
                            double b_hh3 = b_halley * b_halley - 3 * Square(h / s) - 0.25, hh3 = b_hh3 + 2 * Square(bpob) * (1 + 3 / ln_b * (1 + 1 / ln_b)) - 3 * b_halley * bpob * (1 + 2 / ln_b);
                            ds = newton * householder_factor(newton, halley, hh3);
                        }

                        s += ds = Math.Max(-0.5 * s, ds);
                    }
                    return s;
                }
                else {
                    double v_l = normalised_vega(x, s_l), r_lm = RationalCubic.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(b_l, b_c, s_l, s_c, 1 / v_l, 1 / v_c, 0.0, false);
                    s = RationalCubic.rational_cubic_interpolation(beta, b_l, b_c, s_l, s_c, 1 / v_l, 1 / v_c, r_lm);
                    s_left = s_l;
                    s_right = s_c;
                }
            }
            else {
                double s_h = (v_c > Double.MinValue) ? s_c + (b_max - b_c) / v_c : s_c;
                double b_h = normalised_black_call(x, s_h);
                if (beta <= b_h) {
                    double v_h = normalised_vega(x, s_h), r_hm = RationalCubic.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_c, b_h, s_c, s_h, 1 / v_c, 1 / v_h, 0.0, false);
                    s = RationalCubic.rational_cubic_interpolation(beta, b_c, b_h, s_c, s_h, 1 / v_c, 1 / v_h, r_hm);
                    s_left = s_c;
                    s_right = s_h;
                }
                else {
                    double f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2;
                    (f_upper_map_h, d_f_upper_map_h_d_beta, d2_f_upper_map_h_d_beta2) = compute_f_upper_map_and_first_two_derivatives(x, s_h);

                    if (d2_f_upper_map_h_d_beta2 > -SQRT_DBL_MAX && d2_f_upper_map_h_d_beta2 < SQRT_DBL_MAX) {
                        double r_hh = RationalCubic.convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(b_h, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, d2_f_upper_map_h_d_beta2, true);
                        f = RationalCubic.rational_cubic_interpolation(beta, b_h, b_max, f_upper_map_h, 0.0, d_f_upper_map_h_d_beta, -0.5, r_hh);
                    }

                    if (f <= 0) {
                        double h = b_max - b_h, t = (beta - b_h) / h;
                        f = (f_upper_map_h * (1 - t) + 0.5 * h * t) * (1 - t); // We switch to quadratic interpolation using f(b_h), f(b_max)≡0, and f'(b_max)≡-1/2 to specify the quadratic.
                    }

                    s = inverse_f_upper_map(f);
                    s_left = s_h;
                    if (beta > 0.5 * b_max) { // Else we better drop through and let the objective function be g(s) = b(x,s)-beta. 
                                              //
                                              // In this branch, which comprises the upper segment, the objective function is
                                              //     g(s) = ln(b_max-beta)-ln(b_max-b(x,s))
                                              //          ≡ ln((b_max-beta)/(b_max-b(s)))
                                              // This makes
                                              //              g'               =   b'/(b_max-b)
                                              //              newton = -g/g'   =   ln((b_max-b)/(b_max-beta))·(b_max-b)/b'
                                              //              halley = g''/g'  =   b''/b'  +  b'/(b_max-b)
                                              //              hh3    = g'''/g' =   b'''/b' +  g'·(2g'+3b''/b')
                                              // and the iteration is
                                              //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ].
                                              //
                        for (; iterations < N && Math.Abs(ds) > Double.Epsilon * s; ++iterations) {
                            if (ds * ds_previous < 0)
                                ++direction_reversal_count;

                            if (iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right))) {
                                // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
                                // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
                                s = 0.5 * (s_left + s_right);
                                if (s_right - s_left <= Double.Epsilon * s) break;
                                direction_reversal_count = 0;
                                ds = 0;
                            }

                            ds_previous = ds;
                            double b = normalised_black_call(x, s), bp = normalised_vega(x, s);
                            if (b > beta && s < s_right) s_right = s; else if (b < beta && s > s_left) s_left = s; // Tighten the bracket if applicable.
                            if (b >= b_max || bp <= Double.MinValue) // Numerical underflow. Switch to binary nesting for this iteration.
                                ds = 0.5 * (s_left + s_right) - s;
                            else {
                                double b_max_minus_b = b_max - b, g = Math.Log((b_max - beta) / b_max_minus_b), gp = bp / b_max_minus_b;
                                double b_halley = Square(x / s) / s - s / 4, b_hh3 = b_halley * b_halley - 3 * Square(x / (s * s)) - 0.25;
                                double newton = -g / gp, halley = b_halley + gp, hh3 = b_hh3 + gp * (2 * gp + 3 * b_halley);
                                ds = newton * householder_factor(newton, halley, hh3);
                            }

                            s += ds = Math.Max(-0.5 * s, ds);
                        }
                        return s;
                    }
                }
            }

            // In this branch, which comprises the two middle segments, the objective function is g(s) = b(x,s)-beta, or g(s) = b(s) - beta, for short.
            // This makes
            //              newton = -g/g'   =  -(b-beta)/b'
            //              halley = g''/g'  =    b''/b'    =  x²/s³-s/4
            //              hh3    = g'''/g' =    b'''/b'   =  halley² - 3·(x/s²)² - 1/4
            // and the iteration is
            //     s_n+1  =  s_n  +  newton · [ 1 + halley·newton/2 ] / [ 1 + newton·( halley + hh3·newton/6 ) ].
            //
            for (; iterations < N && Math.Abs(ds) > Double.Epsilon * s; ++iterations) {
                if (ds * ds_previous < 0)
                    ++direction_reversal_count;
                if (iterations > 0 && (3 == direction_reversal_count || !(s > s_left && s < s_right))) {
                    // If looping inefficently, or the forecast step takes us outside the bracket, or onto its edges, switch to binary nesting.
                    // NOTE that this can only really happen for very extreme values of |x|, such as |x| = |ln(F/K)| > 500.
                    s = 0.5 * (s_left + s_right);
                    if (s_right - s_left <= Double.Epsilon * s) break;
                    direction_reversal_count = 0;
                    ds = 0;
                }
                ds_previous = ds;
                double b = normalised_black_call(x, s), bp = normalised_vega(x, s);
                if (b > beta && s < s_right) s_right = s; else if (b < beta && s > s_left) s_left = s; // Tighten the bracket if applicable.
                double newton = (beta - b) / bp, halley = Square(x / s) / s - s / 4, hh3 = halley * halley - 3 * Square(x / (s * s)) - 0.25;
                s += ds = Math.Max(-0.5 * s, newton * householder_factor(newton, halley, hh3));
            }

            //return implied_volatility_output(iterations, s);
            return s;
        }

        private static double normalised_black_call(double x, double s) {
            if (x > 0)
                return normalised_intrinsic_call(x) + normalised_black_call(-x, s); // In the money.

            if (s <= Math.Abs(x) * DENORMALISATION_CUTOFF)
                return normalised_intrinsic_call(x); // sigma=0 -> intrinsic value.
                                                     // Denote h := x/s and t := s/2.
                                                     // We evaluate the condition |h|>|η|, i.e., h<η  &&  t < τ+|h|-|η|  avoiding any divisions by s , where η = asymptotic_expansion_accuracy_threshold  and τ = small_t_expansion_of_normalised_black_threshold .

            if (x < s * asymptotic_expansion_accuracy_threshold && 0.5 * s * s + x < s * (small_t_expansion_of_normalised_black_threshold + asymptotic_expansion_accuracy_threshold))
                return asymptotic_expansion_of_normalised_black_call(x / s, 0.5 * s);

            if (0.5 * s < small_t_expansion_of_normalised_black_threshold)
                return small_t_expansion_of_normalised_black_call(x / s, 0.5 * s);

            return normalised_black_call_with_optimal_use_of_codys_functions(x, s);
        }
        private static (double, double, double) compute_f_upper_map_and_first_two_derivatives(double x, double s) {
            double f = NormalDistribution.norm_cdf(-0.5 * s);
            double fp, fpp;

            if (is_below_horizon(x))
                return (f, -0.5, 0.0);

            double w = Square(x / s);
            fp = -0.5 * Math.Exp(0.5 * w);
            fpp = SQRT_PI_OVER_TWO * Math.Exp(w + 0.125 * s * s) * w / s;
            return (f, fp, fpp);
        }

        private static (double, double, double) compute_f_lower_map_and_first_two_derivatives(double x, double s) {
            double ax = Math.Abs(x);
            double z = SQRT_ONE_OVER_THREE * ax / s;
            double y = z * z;
            double s2 = s * s;
            double Phi = NormalDistribution.norm_cdf(-z);
            double phi = NormalDistribution.norm_pdf(z);
            double f = 0.0, fp = 1.0, fpp;

            fpp = PI_OVER_SIX * y / (s2 * s) * Phi * (8.0 * SQRT_THREE * s * ax + (3 * s2 * (s2 - 8.0) - 8.0 * x * x) * Phi / phi) * Math.Exp(2 * y + 0.25 * s2);

            if (is_below_horizon(s))
                return (0.0, 1.0, fpp);

            double Phi2 = Phi * Phi;
            fp = TWO_PI * y * Phi2 * Math.Exp(y + 0.125 * s * s);

            if (!is_below_horizon(x))
                f = TWO_PI_OVER_SQRT_TWENTY_SEVEN * ax * Phi2 * Phi;

            return (f, fp, fpp);
        }

        //
        // Introduced on 2017-02-18
        //
        //     b(x,s)  =  Φ(x/s+s/2)·exp(x/2)  -   Φ(x/s-s/2)·exp(-x/2)
        //             =  Φ(h+t)·exp(x/2)      -   Φ(h-t)·exp(-x/2)
        //             =  ½ · exp(-u²-v²) · [ erfcx(u-v) -  erfcx(u+v) ]
        //             =  ½ · [ exp(x/2)·erfc(u-v)     -  exp(-x/2)·erfc(u+v)    ]
        //             =  ½ · [ exp(x/2)·erfc(u-v)     -  exp(-u²-v²)·erfcx(u+v) ]
        //             =  ½ · [ exp(-u²-v²)·erfcx(u-v) -  exp(-x/2)·erfc(u+v)    ]
        // with
        //              h  =  x/s ,       t  =  s/2 ,
        // and
        //              u  = -h/√2  and   v  =  t/√2 .
        //
        // Cody's erfc() and erfcx() functions each, for some values of their argument, involve the evaluation
        // of the exponential function exp(). The normalised Black function requires additional evaluation(s)
        // of the exponential function irrespective of which of the above formulations is used. However, the total
        // number of exponential function evaluations can be minimised by a judicious choice of one of the above
        // formulations depending on the input values and the branch logic in Cody's erfc() and erfcx().
        //
        static double normalised_black_call_with_optimal_use_of_codys_functions(double x, double s) {
            double codys_threshold = 0.46875;
            double h = x / s;
            double t = 0.5 * s;
            double q1 = -ONE_OVER_SQRT_TWO * (h + t);
            double q2 = -ONE_OVER_SQRT_TWO * (h - t);
            double two_b;

            if (q1 < codys_threshold)
                if (q2 < codys_threshold)
                    two_b = Math.Exp(0.5 * x) * ErfcCody.erfc_cody(q1) - Math.Exp(-0.5 * x) * ErfcCody.erfc_cody(q2);
                else
                    two_b = Math.Exp(0.5 * x) * ErfcCody.erfc_cody(q1) - Math.Exp(-0.5 * (h * h + t * t)) * ErfcCody.erfcx_cody(q2);
            else
                if (q2 < codys_threshold)
                two_b = Math.Exp(-0.5 * (h * h + t * t)) * ErfcCody.erfcx_cody(q1) - Math.Exp(-0.5 * x) * ErfcCody.erfc_cody(q2);
            else
                two_b = Math.Exp(-0.5 * (h * h + t * t)) * (ErfcCody.erfcx_cody(q1) - ErfcCody.erfcx_cody(q2));

            return Math.Abs(Math.Max(0.5 * two_b, 0.0));
        }

        static double normalised_black_call_using_erfcx(double h, double t) {
            // Given h = x/s and t = s/2, the normalised Black function can be written as
            //
            //     b(x,s)  =  Φ(x/s+s/2)·exp(x/2)  -   Φ(x/s-s/2)·exp(-x/2)
            //             =  Φ(h+t)·exp(h·t)      -   Φ(h-t)·exp(-h·t) .                     (*)
            //
            // It is mentioned in section 4 (and discussion of figures 2 and 3) of George Marsaglia's article "Evaluating the
            // Normal Distribution" (available at http://www.jstatsoft.org/v11/a05/paper) that the error of any cumulative normal
            // function Φ(z) is dominated by the hardware (or compiler implementation) accuracy of exp(-z²/2) which is not
            // reliably more than 14 digits when z is large. The accuracy of Φ(z) typically starts coming down to 14 digits when
            // z is around -8. For the (normalised) Black function, as above in (*), this means that we are subtracting two terms
            // that are each products of terms with about 14 digits of accuracy. The net result, in each of the products, is even
            // less accuracy, and then we are taking the difference of these terms, resulting in even less accuracy. When we are
            // using the asymptotic expansion asymptotic_expansion_of_normalised_black_call() invoked in the second branch at the
            // beginning of this function, we are using only *one* exponential instead of 4, and this improves accuracy. It
            // actually improves it a bit more than you would expect from the above logic, namely, almost the full two missing
            // digits (in 64 bit IEEE floating point).  Unfortunately, going higher order in the asymptotic expansion will not
            // enable us to gain more accuracy (by extending the range in which we could use the expansion) since the asymptotic
            // expansion, being a divergent series, can never gain 16 digits of accuracy for z=-8 or just below. The best you can
            // get is about 15 digits (just), for about 35 terms in the series (26.2.12), which would result in an prohibitively
            // long expression in function asymptotic expansion asymptotic_expansion_of_normalised_black_call(). In this last branch,
            // here, we therefore take a different tack as follows.
            //     The "scaled complementary error function" is defined as erfcx(z) = exp(z²)·erfc(z). Cody's implementation of this
            // function as published in "Rational Chebyshev approximations for the error function", W. J. Cody, Math. Comp., 1969, pp.
            // 631-638, uses rational functions that theoretically approximates erfcx(x) to at least 18 significant decimal digits,
            // *without* the use of the exponential function when x>4, which translates to about z<-5.66 in Φ(z). To make use of it,
            // we write
            //             Φ(z) = exp(-z²/2)·erfcx(-z/√2)/2
            //
            // to transform the normalised black function to
            //
            //   b   =  ½ · exp(-½(h²+t²)) · [ erfcx(-(h+t)/√2) -  erfcx(-(h-t)/√2) ]
            //
            // which now involves only one exponential, instead of three, when |h|+|t| > 5.66 , and the difference inside the
            // square bracket is between the evaluation of two rational functions, which, typically, according to Marsaglia,
            // retains the full 16 digits of accuracy (or just a little less than that).
            //
            double b = 0.5 * Math.Exp(-0.5 * (h * h + t * t)) * (ErfcCody.erfcx_cody(-ONE_OVER_SQRT_TWO * (h + t)) - ErfcCody.erfcx_cody(-ONE_OVER_SQRT_TWO * (h - t)));
            return Math.Abs(Math.Max(b, 0.0));
        }

        //     b(x,s)  =  Φ(x/s+s/2)·exp(x/2)  -   Φ(x/s-s/2)·exp(-x/2)
        //             =  Φ(h+t)·exp(x/2)      -   Φ(h-t)·exp(-x/2)
        // with
        //              h  =  x/s   and   t  =  s/2
        static double normalised_black_call_using_norm_cdf(double x, double s) {
            double h = x / s, t = 0.5 * s, b_max = Math.Exp(0.5 * x), b = NormalDistribution.norm_cdf(h + t) * b_max - NormalDistribution.norm_cdf(h - t) / b_max;
            return Math.Abs(Math.Max(b, 0.0));
        }

        // Asymptotic expansion of
        //
        //              b  =  Φ(h+t)·exp(x/2) - Φ(h-t)·exp(-x/2)
        // with
        //              h  =  x/s   and   t  =  s/2
        // which makes
        //              b  =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
        //
        //                    exp(-(h²+t²)/2)
        //                 =  ---------------  ·  [ Y(h+t) - Y(h-t) ]
        //                        √(2π)
        // with
        //           Y(z) := Φ(z)/φ(z)
        //
        // for large negative (t-|h|) by the aid of Abramowitz & Stegun (26.2.12) where Φ(z) = φ(z)/|z|·[1-1/z^2+...].
        // We define
        //                     r
        //         A(h,t) :=  --- · [ Y(h+t) - Y(h-t) ]
        //                     t
        //
        // with r := (h+t)·(h-t) and give an expansion for A(h,t) in q:=(h/r)² expressed in terms of e:=(t/h)² .
        static double asymptotic_expansion_of_normalised_black_call(double h, double t) {
            double e = (t / h) * (t / h), r = ((h + t) * (h - t)), q = (h / r) * (h / r);

            // 17th order asymptotic expansion of A(h,t) in q, sufficient for Φ(h) [and thus y(h)] to have relative accuracy of 1.64E-16 for h <= η  with  η:=-10.
            double asymptotic_expansion_sum = (2.0 + q * (-6.0E0 - 2.0 * e + 3.0 * q * (1.0E1 + e * (2.0E1 + 2.0 * e) + 5.0 * q * (-1.4E1 + e * (-7.0E1 + e * (-4.2E1 - 2.0 * e)) +
                7.0 * q * (1.8E1 + e * (1.68E2 + e * (2.52E2 + e * (7.2E1 + 2.0 * e))) + 9.0 * q * (-2.2E1 + e * (-3.3E2 + e * (-9.24E2 + e * (-6.6E2 + e * (-1.1E2 - 2.0 * e)))) +
                1.1E1 * q * (2.6E1 + e * (5.72E2 + e * (2.574E3 + e * (3.432E3 + e * (1.43E3 + e * (1.56E2 + 2.0 * e))))) + 1.3E1 * q * (-3.0E1 + e * (-9.1E2 + e * (-6.006E3 +
                e * (-1.287E4 + e * (-1.001E4 + e * (-2.73E3 + e * (-2.1E2 - 2.0 * e)))))) + 1.5E1 * q * (3.4E1 + e * (1.36E3 + e * (1.2376E4 + e * (3.8896E4 + e * (4.862E4 +
                e * (2.4752E4 + e * (4.76E3 + e * (2.72E2 + 2.0 * e))))))) + 1.7E1 * q * (-3.8E1 + e * (-1.938E3 + e * (-2.3256E4 + e * (-1.00776E5 + e * (-1.84756E5 +
                e * (-1.51164E5 + e * (-5.4264E4 + e * (-7.752E3 + e * (-3.42E2 - 2.0 * e)))))))) + 1.9E1 * q * (4.2E1 + e * (2.66E3 + e * (4.0698E4 + e * (2.3256E5 +
                e * (5.8786E5 + e * (7.05432E5 + e * (4.0698E5 + e * (1.08528E5 + e * (1.197E4 + e * (4.2E2 + 2.0 * e))))))))) + 2.1E1 * q * (-4.6E1 + e * (-3.542E3 +
                e * (-6.7298E4 + e * (-4.90314E5 + e * (-1.63438E6 + e * (-2.704156E6 + e * (-2.288132E6 + e * (-9.80628E5 + e * (-2.01894E5 + e * (-1.771E4 + e * (-5.06E2 - 2.0 * e)))))))))) +
                2.3E1 * q * (5.0E1 + e * (4.6E3 + e * (1.0626E5 + e * (9.614E5 + e * (4.08595E6 + e * (8.9148E6 + e * (1.04006E7 + e * (6.53752E6 + e * (2.16315E6 + e * (3.542E5 + e * (2.53E4 +
                e * (6.0E2 + 2.0 * e))))))))))) + 2.5E1 * q * (-5.4E1 + e * (-5.85E3 + e * (-1.6146E5 + e * (-1.77606E6 + e * (-9.37365E6 + e * (-2.607579E7 + e * (-4.01166E7 + e * (-3.476772E7 +
                e * (-1.687257E7 + e * (-4.44015E6 + e * (-5.9202E5 + e * (-3.51E4 + e * (-7.02E2 - 2.0 * e)))))))))))) + 2.7E1 * q * (5.8E1 + e * (7.308E3 + e * (2.3751E5 + e * (3.12156E6 +
                e * (2.003001E7 + e * (6.919458E7 + e * (1.3572783E8 + e * (1.5511752E8 + e * (1.0379187E8 + e * (4.006002E7 + e * (8.58429E6 + e * (9.5004E5 + e * (4.7502E4 + e * (8.12E2 +
                2.0 * e))))))))))))) + 2.9E1 * q * (-6.2E1 + e * (-8.99E3 + e * (-3.39822E5 + e * (-5.25915E6 + e * (-4.032015E7 + e * (-1.6934463E8 + e * (-4.1250615E8 + e * (-6.0108039E8 +
                e * (-5.3036505E8 + e * (-2.8224105E8 + e * (-8.870433E7 + e * (-1.577745E7 + e * (-1.472562E6 + e * (-6.293E4 + e * (-9.3E2 - 2.0 * e)))))))))))))) + 3.1E1 * q * (6.6E1 +
                e * (1.0912E4 + e * (4.74672E5 + e * (8.544096E6 + e * (7.71342E7 + e * (3.8707344E8 + e * (1.14633288E9 + e * (2.07431664E9 + e * (2.33360622E9 + e * (1.6376184E9 +
                e * (7.0963464E8 + e * (1.8512208E8 + e * (2.7768312E7 + e * (2.215136E6 + e * (8.184E4 + e * (1.056E3 + 2.0 * e))))))))))))))) + 3.3E1 * (-7.0E1 + e * (-1.309E4 +
                e * (-6.49264E5 + e * (-1.344904E7 + e * (-1.4121492E8 + e * (-8.344518E8 + e * (-2.9526756E9 + e * (-6.49588632E9 + e * (-9.0751353E9 + e * (-8.1198579E9 + e * (-4.6399188E9 +
                e * (-1.6689036E9 + e * (-3.67158792E8 + e * (-4.707164E7 + e * (-3.24632E6 + e * (-1.0472E5 + e * (-1.19E3 - 2.0 * e))))))))))))))))) * q)))))))))))))))));
            double b = ONE_OVER_SQRT_TWO_PI * Math.Exp((-0.5 * (h * h + t * t))) * (t / r) * asymptotic_expansion_sum;
            return Math.Abs(Math.Max(b, 0.0));
        }

        // Calculation of
        //
        //              b  =  Φ(h+t)·exp(h·t) - Φ(h-t)·exp(-h·t)
        //
        //                    exp(-(h²+t²)/2)
        //                 =  --------------- ·  [ Y(h+t) - Y(h-t) ]
        //                        √(2π)
        // with
        //           Y(z) := Φ(z)/φ(z)
        //
        // using an expansion of Y(h+t)-Y(h-t) for small t to twelvth order in t.
        // Theoretically accurate to (better than) precision  ε = 2.23E-16  when  h<=0  and  t < τ  with  τ := 2·ε^(1/16) ≈ 0.21.
        // The main bottleneck for precision is the coefficient a:=1+h·Y(h) when |h|>1 .
        static double small_t_expansion_of_normalised_black_call(double h, double t) {
            // Y(h) := Φ(h)/φ(h) = √(π/2)·erfcx(-h/√2)
            // a := 1+h·Y(h)  --- Note that due to h<0, and h·Y(h) -> -1 (from above) as h -> -∞, we also have that a>0 and a -> 0 as h -> -∞
            // w := t² , h2 := h²
            double a = 1 + h * (0.5 * SQRT_TWO_PI) * ErfcCody.erfcx_cody(-ONE_OVER_SQRT_TWO * h), w = t * t, h2 = h * h;
            double expansion = 2 * t * (a + w * ((-1 + 3 * a + a * h2) / 6 + w * ((-7 + 15 * a + h2 * (-1 + 10 * a + a * h2)) / 120 + w * ((-57 + 105 * a + h2 * (-18 + 105 * a + h2 * (-1 + 21 * a + a * h2))) / 5040 +
                w * ((-561 + 945 * a + h2 * (-285 + 1260 * a + h2 * (-33 + 378 * a + h2 * (-1 + 36 * a + a * h2)))) / 362880 + w * ((-6555 + 10395 * a + h2 * (-4680 + 17325 * a + h2 * (-840 + 6930 * a + h2 * (-52 + 990 * a +
                h2 * (-1 + 55 * a + a * h2))))) / 39916800 + ((-89055 + 135135 * a + h2 * (-82845 + 270270 * a + h2 * (-20370 + 135135 * a + h2 * (-1926 + 25740 * a + h2 * (-75 + 2145 * a + h2 * (-1 + 78 * a + a * h2)))))) * w) / 6227020800.0))))));
            double b = ONE_OVER_SQRT_TWO_PI * Math.Exp((-0.5 * (h * h + t * t))) * expansion;
            return Math.Abs(Math.Max(b, 0.0));
        }

        private static double inverse_f_upper_map(double f) {
            return -2.0 * NormalDistribution.inverse_norm_cdf(f);
        }

        private static double inverse_f_lower_map(double x, double f) {
            return is_below_horizon(f) ? 0 : Math.Abs(x / (SQRT_THREE * NormalDistribution.inverse_norm_cdf(Math.Pow(f / (TWO_PI_OVER_SQRT_TWENTY_SEVEN * Math.Abs(x)), 1.0 / 3.0))));
        }

        private static double normalised_intrinsic_call(double x) {
            return normalised_intrinsic(x, 1);
        }

        private static double normalised_intrinsic(double x, int optionType) {
            if (optionType * x <= 0)
                return 0;

            double x2 = x * x;
            if (x2 < 98 * FOURTH_ROOT_DBL_EPSILON) // The factor 98 is computed from last coefficient: √√92897280 = 98.1749
                return Math.Abs(Math.Max((optionType < 0 ? -1 : 1) * x * (1 + x2 * ((1.0 / 24.0) + x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) + (1.0 / 92897280.0) * x2)))), 0.0));
            double b_max = Math.Exp(0.5 * x), one_over_b_max = 1 / b_max;
            return Math.Abs(Math.Max(optionType * (b_max - one_over_b_max), 0.0));
        }

        private static double normalised_vega(double x, double s) {
            double ax = Math.Abs(x);
            return (ax <= 0) ? ONE_OVER_SQRT_TWO_PI * Math.Exp(-0.125 * s * s) : ((s <= 0 || s <= ax * SQRT_DBL_MIN) ? 0 : ONE_OVER_SQRT_TWO_PI * Math.Exp(-0.5 * (Square(x / s) + Square(0.5 * s))));
        }

        private static double Square(double x) {
            return x * x;
        }

    }

}

