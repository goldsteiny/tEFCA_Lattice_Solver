using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Tools.IO;
using Metabolism.Analysis;
using DotNumerics.LinearAlgebra;

namespace Metabolism.Patterns
{

    public interface IPattern<P> : IComparable<P> where P : IPattern<P>
    {
        int Count { get; }
        List<int> Values { get; }
        bool[] Incidence { get; }
        bool this[int index] { get; }
    }

    public interface IFluxPattern<P> : IPattern<P> where P : IFluxPattern<P>
    {
        int Length { get; }
    }

    public interface IFluxVector<P> : IFluxPattern<P> where P : IFluxVector<P>
    {
        Vector Fluxes { get; }
    }

    public abstract class Pattern<P> : IPattern<P> where P : Pattern<P>
    {
        private readonly bool[] pattern;
        private readonly SortedSet<int> reactions;
        private string as_string;

        public Pattern(bool[] pattern)
        {
            this.pattern = (bool[])pattern.Clone();
            reactions = new SortedSet<int>();
        }

        private void CreateSet()
        {
            StringBuilder str = new StringBuilder();
            for (int i = 0; i < pattern.Length; ++i)
                if (pattern[i])
                {
                    reactions.Add(i);
                    str.Append(i + 1).Append(' ');
                }
            as_string = str.Length > 1 ? '(' + str.ToString(0, str.Length - 1) + ')' : "()";
        }

        public override string ToString()
        {
            if (as_string == null)
                CreateSet();
            return as_string;
        }

        public List<int> Values
        {
            get
            {
                if (as_string == null)
                    CreateSet();
                return reactions.ToList();
            }
        }

        public bool this[int index]
        {
            get
            {
                if (index < pattern.Length && index > -1)
                    return pattern[index];
                else
                    return false;
            }
        }

        public virtual int Count
        {
            get
            {
                if (as_string == null)
                    CreateSet();
                return reactions.Count;
            }
        }


        public static bool operator <(Pattern<P> a1, Pattern<P> a2)
        {

            if (a1.Count < a2.Count)
            {
                List<int> reactions = a1.Values;
                bool smaller = true;
                for (int i = 0; i < reactions.Count && smaller; ++i)
                    smaller = a2[(reactions[i])];
                return smaller;
            }

            return false;
        }

        public static bool operator >(Pattern<P> a1, Pattern<P> a2)
        {
            return a2 < a1;
        }

        public abstract int CompareTo(P other);


        public bool[] Incidence
        {
            get { return (bool[])this.pattern.Clone(); }
        }

    }


    public class FluxPattern : Pattern<FluxPattern>, IFluxPattern<FluxPattern>
    {
        private readonly bool rev;
        private readonly int n;

        public bool Reversible { get { return rev; } }


        public FluxPattern(bool[] pattern, bool reversible)
            : base(pattern)
        {
            this.rev = reversible;
            this.n = pattern.Length;
        }

        public static FluxPattern operator +(FluxPattern a1, FluxPattern a2)
        {

            bool rev = a1.Reversible && a2.Reversible;
            bool[] sum = new bool[a1.Length];
            for (int i = 0; i < a1.Length; ++i)
                sum[i] = a1[i] || a2[i];

            return new FluxPattern(sum, rev);
        }


        public int Length
        {
            get { return n; }
        }


        public override int CompareTo(FluxPattern other)
        {
            if (this.Count < other.Count)
                return -1;
            if (this.Count > other.Count)
                return 1;

            List<int> reac_a = this.Values, reac_b = other.Values;
            for (int i = 0; i < this.Count; ++i)
                if (reac_a.ElementAt(i) < reac_b.ElementAt(i))
                    return -1;
                else if (reac_a.ElementAt(i) > reac_b.ElementAt(i))
                    return 1;

            return 0;
        }

    }

    public class FluxVector : IFluxVector<FluxVector>
    {

        private readonly FluxPattern pattern;
        private readonly Vector vector;

        public FluxVector(double[] values, bool reversible, double tolerance)
        {
            this.vector = new Vector(VectorType.Column, values);

            bool[] pattern = new bool[values.Length];
            for (int i = 0; i < pattern.Length; ++i)
                pattern[i] = Math.Abs(values[i]) < tolerance;

            this.pattern = new FluxPattern(pattern, reversible);
        }



        public int Count
        {
            get { return pattern.Count; }
        }

        public List<int> Values
        {
            get { return pattern.Values; }
        }

        public int CompareTo(FluxVector other)
        {
            return pattern.CompareTo(other.pattern);
        }




        public Vector Fluxes
        {
            get { return new Vector(VectorType.Column, vector.ToArray()); }
        }



        public bool[] Incidence
        {
            get { return this.pattern.Incidence; }
        }

        public int Length
        {
            get { return this.pattern.Length; }
        }


        public bool this[int index]
        {
            get { return this.pattern[index]; }
        }
    }
}