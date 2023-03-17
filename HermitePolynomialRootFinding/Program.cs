using MultiPrecision;

namespace HermitePolynomialRootFinding {
    internal class Program {
        static void Main() {
            using StreamWriter sw = new("../../../../results/roots_n64_2.csv");
            using BinaryWriter sbw = new(File.Open("../../../../results_disused/roots_n64.bin", FileMode.Create));

            sw.WriteLine("x,w,wexp(x^2),Hn(x+eps),Hn(x-eps)");

            MultiPrecision<Plus32<Pow2.N64>> dx = 0.125;

            for (int n = 4; n <= 128; n++) {
                Console.WriteLine($"H{n} root finding...");

                Polynomial<Plus32<Pow2.N64>> poly = PolynomialGenerator.Table<Plus32<Pow2.N64>>(n);
                Polynomial<Plus32<Pow2.N64>> poly_nm1 = PolynomialGenerator.Table<Plus32<Pow2.N64>>(n - 1);

                List<MultiPrecision<Plus32<Pow2.N64>>> roots = new();

                while (dx > 0) {
                    roots.Clear();

                    if ((n & 1) == 1) {
                        roots.Add(0);
                    }

                    MultiPrecision<Plus32<Pow2.N64>> prev_y = poly.Value(dx / 2);

                    for (MultiPrecision<Plus32<Pow2.N64>> x = dx * 3 / 2; x <= n / 2.0; x += dx) {
                        MultiPrecision<Plus32<Pow2.N64>> y = poly.Value(x);

                        if (prev_y.Sign != y.Sign) {
                            MultiPrecision<Plus32<Pow2.N64>> px = (x - dx) - dx * prev_y / (y - prev_y);

                            MultiPrecision<Plus32<Pow2.N64>> root = MultiPrecisionUtil.NewtonRaphsonRootFinding(
                                px,
                                poly.Value, poly.Diff,
                                x_min: roots.Count > 0 ? (roots[^1] + dx / 256) : 0,
                                null, break_overshoot: false, max_iterations: 256
                            );

                            roots.Add(root);

                            Console.WriteLine($"  root found x={root:e10}");
                        }

                        prev_y = y;
                    }

                    if (roots.Where(root => root.IsFinite).Count() >= (n + 1) / 2) {
                        break;
                    }

                    dx /= 2;
                    Console.WriteLine($"  reduce search width. dx={dx}");
                }

                sw.WriteLine($"H{n}");

                sbw.Write(n);

                foreach (MultiPrecision<Plus32<Pow2.N64>> x in roots) {
                    MultiPrecision<Pow2.N64> w = GaussHermiteWeight<Plus32<Pow2.N64>>.Eval(x, n, poly_nm1).Convert<Pow2.N64>();

                    MultiPrecision<Pow2.N64> wexp = (GaussHermiteWeight<Plus32<Pow2.N64>>.Eval(x, n, poly_nm1) * MultiPrecision<Plus32<Pow2.N64>>.Exp(x * x)).Convert<Pow2.N64>();

                    MultiPrecision<Pow2.N64> my = poly.Value(
                        MultiPrecision<Pow2.N64>.BitDecrement(x.Convert<Pow2.N64>()).Convert<Plus32<Pow2.N64>>()
                    ).Convert<Pow2.N64>();

                    MultiPrecision<Pow2.N64> py = poly.Value(
                        MultiPrecision<Pow2.N64>.BitIncrement(x.Convert<Pow2.N64>()).Convert<Plus32<Pow2.N64>>()
                    ).Convert<Pow2.N64>();

                    sw.WriteLine($"{x.Convert<Pow2.N64>()},{w},{wexp},{my:e4},{py:e4}");

                    sbw.Write(x.Convert<Pow2.N64>());
                    sbw.Write(w);
                    sbw.Write(wexp);
                }

                sw.Flush();
                sbw.Flush();
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}
