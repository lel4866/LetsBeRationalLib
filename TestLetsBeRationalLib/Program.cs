using System;
using LetsBeRationalLib;


namespace TestLetsBeRationalLib
{
    class Program
    {
        static void Main(string[] args)
        {
            double price = 13.2;
            double r = 0.0016324; // risk free rate(1 year treasury yield)
            double d = 0.0194; // trailing 12 - month sp500 dividend yield
            DateTime quote_dt = new DateTime(2014, 1, 2);
            DateTime exp_dt = new DateTime(2014, 1, 31);
            int dte = (exp_dt - quote_dt).Days;
            double t = dte / 365.0; // days to expiration / days in year
            double s = 1837.73; // underlying SPX price
            double K = 1800.0;

            // EXPORT_EXTERN_C double implied_volatility_from_a_transformed_rational_guess(double price, double F, double K, double T, double q /* q=±1 */)
            double iv = LetsBeRational.ImpliedVolatility(price, s, K, t, r, d, LetsBeRational.OptionType.Put);

            Console.WriteLine($"IV = {iv}");
        }
    }
}
