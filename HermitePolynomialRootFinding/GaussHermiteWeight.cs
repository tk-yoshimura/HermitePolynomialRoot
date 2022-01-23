using MultiPrecision;

namespace HermitePolynomialRootFinding {
    internal static class GaussHermiteWeight<N> where N : struct, IConstant {
        public static MultiPrecision<N> Eval(MultiPrecision<N> x_root, int n, Polynomial<N> poly_nm1) {
            MultiPrecision<N> w = (MultiPrecision<N>.Gamma(n + 1) * MultiPrecision<N>.Sqrt(MultiPrecision<N>.PI))
                                   / MultiPrecision<N>.Square(n * poly_nm1.Value(x_root));

            w = MultiPrecision<N>.Ldexp(w, n - 1);

            return w;
        }
    }
}
