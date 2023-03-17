using MultiPrecision;
using System.Numerics;

namespace HermitePolynomialRootFinding {
    internal static class PolynomialGenerator {
        static readonly List<BigInteger[]> coef_table = new();

        static PolynomialGenerator() {
            coef_table.Add(new BigInteger[] { 1 });
            coef_table.Add(new BigInteger[] { 0, 2 });
        }

        public static Polynomial<N> Table<N>(int n) where N : struct, IConstant {
            if (n < coef_table.Count) {
                return new Polynomial<N>(coef_table[n].Select(c => (MultiPrecision<N>)c).ToArray());
            }

            for (int i = coef_table.Count; i <= n; i++) {
                BigInteger[] coef = new BigInteger[i + 1];

                coef[0] = -2 * (coef_table[i - 2][0] * (i - 1));

                for (int k = 1; k <= i - 2; k++) {
                    coef[k] = 2 * (coef_table[i - 1][k - 1] - coef_table[i - 2][k] * (i - 1));
                }

                coef[i - 1] = 0;
                coef[i] = 2 * coef_table[i - 1][i - 1];

                coef_table.Add(coef);
            }

            return new Polynomial<N>(coef_table[n].Select(c => (MultiPrecision<N>)c).ToArray());
        }
    }
}
