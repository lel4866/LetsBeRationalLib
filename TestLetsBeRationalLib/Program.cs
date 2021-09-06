using System;
using LetsBeRationalLib;


namespace TestLetsBeRationalLib
{
    class Program
    {
        static void Main(string[] args)
        {
            double price = 6.8;
            double r = 0.0225; // risk free rate(1 year treasury yield)
            double d = 0.0192; // trailing 12 - month sp500 dividend yield
            double t = 34.0 / 365.0; // days to expiration / days in year
            double s = 2736.27; // underlying SPX price
            double K = 2475.0;

            // EXPORT_EXTERN_C double implied_volatility_from_a_transformed_rational_guess(double price, double F, double K, double T, double q /* q=±1 */)
            double iv = LetsBeRational.ImpliedVolatility(price, s, K, t, r, d, LetsBeRational.OptionType.Put);
            int a = 1;

            Console.WriteLine($"IV = {iv}");
        }
    }
}
