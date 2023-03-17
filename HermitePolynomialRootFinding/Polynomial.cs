using MultiPrecision;

namespace HermitePolynomialRootFinding {
    internal class Polynomial<N> where N : struct, IConstant {
        readonly MultiPrecision<N>[] cs, ds;

        public Polynomial(IReadOnlyList<MultiPrecision<N>> cs) {
            MultiPrecision<N>[] ds = new MultiPrecision<N>[cs.Count - 1];

            for (int i = 1; i < cs.Count; i++) {
                ds[i - 1] = cs[i] * i;
            }

            this.cs = cs.ToArray();
            this.ds = ds;
        }

        public MultiPrecision<N> Value(MultiPrecision<N> x) {
            MultiPrecision<N> s = cs[^1];

            for (int i = cs.Length - 2; i >= 0; i--) {
                s = s * x + cs[i];
            }

            return s;
        }

        public MultiPrecision<N> Diff(MultiPrecision<N> x) {
            MultiPrecision<N> s = ds[^1];

            for (int i = ds.Length - 2; i >= 0; i--) {
                s = s * x + ds[i];
            }

            return s;
        }
    }
}
