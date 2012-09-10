using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DotNumerics.LinearAlgebra;
using Metabolism.Patterns;
using Metabolism.Analysis;

namespace Metabolism.Network
{
    class Compressor
    {
        public static Matrix CalculateKernelSVD(double[,] matrix)
        {
            SingularValueDecomposition svd = new SingularValueDecomposition();

            Matrix A = new Matrix(matrix), U, S, VT;
            Vector s;
            svd.ComputeSVD(A, out S, out U, out VT);
            svd.ComputeSVD(A, out s);

            int kernel_dim = S.ColumnCount - s.Length;

            if (kernel_dim == 0)
                return null;

            Matrix res = new Matrix(S.ColumnCount, kernel_dim);
            for (int i = 0; i < S.ColumnCount; ++i)
                for (int j = 0; j < kernel_dim; ++j)
                    res[i, j] = VT[S.ColumnCount - kernel_dim + j, i];

            return res;
        }

        public static LinkedListCoupling Uncompress(IIntCoupling coupled, List<List<int>> classes, double[] multiplicators)
        {

            int n = multiplicators.Length;

            LinkedListCoupling res = new LinkedListCoupling(n);

            ICollection<int> coupling_i;
            for (int k = 0; k < coupled.Length; ++k)
            {
                coupling_i = coupled[k];
                foreach (int l in coupling_i)
                {
                    foreach (int i in classes[k])
                        foreach (int j in classes[l])
                            if (!coupled[l].Contains(k))
                                res.Add(new Couple(i, j));
                            else if (k < l)
                                res.Add(new PartialCouple(i, j));
                }
            }

            Console.WriteLine("({0})\tCompressed couples uncompressed.\n", System.DateTime.Now);

            foreach (List<int> cl in classes)
                for (int i = 0; i < cl.Count; ++i)
                    for (int j = i + 1; j < cl.Count; ++j)
                        res.Add(new FullCouple(cl[i], cl[j], multiplicators[cl[j]], multiplicators[cl[i]]));

            Console.WriteLine("({0})\tFully coupled pairs added.\n", System.DateTime.Now);

            return res;
        }

        public static LinkedList<FluxPattern> Compress(LinkedList<FluxPattern> patterns, List<List<int>> classes)
        {
            LinkedList<FluxPattern> res = new LinkedList<FluxPattern>();
            if (patterns == null || patterns.Count < 1)
                return res;

            bool[] pattern;
            foreach (FluxPattern a in patterns)
            {
                pattern = new bool[classes.Count];
                for (int i = 0; i < pattern.Length; ++i)
                    pattern[i] = a[classes[i][0]];
                res.AddLast(new FluxPattern(pattern, a.Reversible));
            }

            return res;
        }

        public static FluxPattern Uncompress(FluxPattern a, int n, List<List<int>> classes)
        {
            bool[] pattern = new bool[n];
            ICollection<int> reactions = a.Values;
            foreach (int r in reactions)
                foreach (int i in classes[r])
                    pattern[i] = true;
            return new FluxPattern(pattern, a.Reversible);
        }

        public static int Compress(double[,] matrix, bool[] rev, bool[] blocked, Matrix kernel, double tolerance, out double[,] matrix_small, out bool[] rev_small, out List<List<int>> classes, out double[] multiplicators)
        {
            int n = rev.Length, m = matrix.Length / n;
            Vector[] K = kernel.GetRowVectors();

            classes = new List<List<int>>();
            multiplicators = new double[n];

            bool found;
            for (int i = 0; i < n; ++i)
                if (!blocked[i] && !(blocked[i] = K[i].NormInf() < tolerance))
                {
                    found = false;
                    for (int j = 0; j < classes.Count && !found; ++j)
                    {
                        multiplicators[i] = Vector.DotProduct(K[classes[j][0]], K[i].Transpose()) / (Vector.DotProduct(K[classes[j][0]], K[classes[j][0]].Transpose()));
                        if (found = (K[i] - K[classes[j][0]].Multiply(multiplicators[i])).NormInf() < tolerance)
                        {
                            classes[j].Add(i);
                            if (Math.Abs(multiplicators[i] - Math.Round(multiplicators[i])) < tolerance)
                                multiplicators[i] = Math.Round(multiplicators[i]);
                        }
                    }
                    if (!found)
                    {
                        classes.Add(new List<int>());
                        classes[classes.Count - 1].Add(i);
                        multiplicators[i] = 1;
                    }

                    if (Math.Abs(multiplicators[i]) < tolerance)
                        Console.WriteLine("BÄH!!!");
                }

            bool switched;
            foreach (List<int> cl in classes)
            {
                switched = false;
                foreach (int j in cl)
                {
                    if (switched = (!switched && !rev[j] && multiplicators[j] < 0))
                        foreach (int k in cl)
                            multiplicators[k] = -multiplicators[k];
                }
            }

            matrix_small = new double[m, classes.Count];
            rev_small = new bool[classes.Count];
            bool contains_irr;
            for (int i = 0; i < classes.Count; ++i)
            {
                contains_irr = false;
                foreach (int j in classes[i])
                {
                    contains_irr |= !rev[j];
                    for (int k = 0; k < m; ++k)
                        matrix_small[k, i] += multiplicators[j] * matrix[k, j];
                }
                for (int k = 0; k < m; ++k)
                    if (Math.Abs(matrix_small[k, i] - Math.Round(matrix_small[k, i])) < tolerance)
                        matrix_small[k, i] = Math.Round(matrix_small[k, i]);
                rev_small[i] = !contains_irr;
            }

            return n - classes.Count;
        }
    }
}
