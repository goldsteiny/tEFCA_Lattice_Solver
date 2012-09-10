using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gurobi;
using Tools.IO;
using Metabolism.Patterns;
using Tools.Converter.Biology;

namespace Metabolism.Analysis
{

    /// <summary>
    /// This interface you to add two coordinates I and J to a class. Afterwards you can sort objects of this type e.g. lexicographically.
    /// </summary>
    public interface IPair<T> : IComparable<T> where T : IPair<T>
    {

        int I { get; }
        int J { get; }
    }

    /// <summary>
    /// This interface allows to create a class that tests whether a certain property holds for a flux pattern "a". E.g. $5 \notin a$.
    /// </summary>
    public interface ILatticeProperty
    {
        /// <summary>
        /// Tests whether a flux pattern "a" fulfills the classes property.
        /// </summary>
        /// <param name="a">The pattern you want to test.</param>
        /// <returns>true iff "a" fulfills the property.</returns>
        bool HoldsFor<FP>(FP a) where FP : IFluxPattern<FP>;
    }

    /// <summary>
    /// A property that has to hold for a flux pattern "a" if $I \notin a$
    /// </summary>
    public interface INoIImpliesProperty : ILatticeProperty
    {

        /// <summary>
        /// The coordinate relevant for the property.
        /// </summary>
        int I { get; }

        string ToString(String[] reaction_names);
    }


    /// <summary>
    /// A pair of two coordinates. You can order a set of Pair objects lexicographically.
    /// </summary>
    public class Pair<T> : IPair<T> where T : Pair<T>
    {

        public readonly int i, j;

        /// <summary>
        /// Creates a new Pair object for the fixed coordinates "i" and "j".
        /// </summary>
        public Pair(int i, int j)
        {
            this.i = i;
            this.j = j;
        }

        public int I
        {
            get { return i; }
        }

        public int J
        {
            get { return j; }
        }

        /// <summary>
        /// Lexicographical ordering of Pair objects.
        /// </summary>
        public int CompareTo(T other)
        {
            if (this.i != other.i)
                return this.i < other.i ? -1 : 1;
            if (this.j != other.j)
                return this.j < other.j ? -1 : 1;
            return 0;
        }

        public override string ToString()
        {
            return String.Format("({0}, {1})", i + 1, j + 1);
        }

        public string ToString(String[] reaction_names)
        {
            return String.Format("({0}, {1})", reaction_names[i], reaction_names[j]);
        }

        public override int GetHashCode()
        {
            return this.ToString().GetHashCode();
        }
    }


    /// <summary>
    /// An abstract class that tests whether a flux pattern "a" fulfills a property regarding a pair (i,j) of coordinates.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public abstract class LatticePropertyPair<T> : Pair<T>, ILatticeProperty where T : LatticePropertyPair<T>
    {

        /// <summary>
        /// Creates an instance of the class for the fixed coordinates (i,j)
        /// </summary>
        public LatticePropertyPair(int i, int j) : base(i, j) { }
        public abstract bool HoldsFor<FP>(FP a) where FP : IFluxPattern<FP>;
    }


    /// <summary>
    /// A pair (i,j) of coordinates that represents flux coupling $i \rightarrow j$ for the flux lattice.
    /// </summary>
    public class Couple : LatticePropertyPair<Couple>, INoIImpliesProperty
    {

        public Couple(int i, int j) : base(i, j) { }

        public override string ToString()
        {
            return String.Format("{0} -> {1}", i + 1, j + 1);
        }

        public string ToString(String[] reaction_names)
        {
            return String.Format("{0} -> {1}", reaction_names[i], reaction_names[j]);
        }

        public override int GetHashCode()
        {
            return this.ToString().GetHashCode();
        }

        /// <summary>
        /// Tests whether $i \rightarrow j$ holds for the flux pattern "a".
        /// </summary>
        public override bool HoldsFor<FP>(FP a)
        {
            return !a[j] || a[i] || a.Count == 0;
        }
    }

    /// <summary>
    /// A pair (i,j) of coordinates that represents partial flux coupling $i \leftrightarrow j$ for the flux lattice.
    /// </summary>
    public class PartialCouple : Couple
    {
        public PartialCouple(int i, int j) : base(Math.Min(i, j), Math.Max(i, j)) { }

        public override string ToString()
        {
            return String.Format("{0} <-> {1}", i + 1, j + 1);
        }

        public string ToString(String[] reaction_names)
        {
            return String.Format("{0} -> {1}\n{1} -> {0}", reaction_names[i], reaction_names[j]);
        }


        /// <summary>
        /// Tests whether $i \leftrightarrow j$ holds for the flux pattern "a".
        /// </summary>
        public override bool HoldsFor<FP>(FP a)
        {
            return a[j] == a[i];
        }
    }

    /// <summary>
    /// A pair (i,j) of coordinates that represents fully flux coupling $i \leftrightarrow_\lambda j$ for the flux lattice.
    /// </summary>
    public class FullCouple : PartialCouple
    {
        public readonly double lambda_i, lambda_j;

        /// <summary>
        /// Creates a new instance of the class.
        /// </summary>
        /// <remarks>
        /// $lambda_i f_i = lambda_j f_j \forall f \in F$.
        /// The "HoldsFor" method only test for partial coupling!
        /// </remarks>
        public FullCouple(int i, int j, double lambda_i, double lambda_j)
            : base(i, j)
        {
            this.lambda_i = i < j ? lambda_i : lambda_j;
            this.lambda_j = i < j ? lambda_j : lambda_i;
        }

        public override string ToString()
        {
            return String.Format("{0} <-> {1} ({2:0.0###} : {3:0.0###})", i + 1, j + 1, lambda_i, lambda_j);
        }

        public string ToString(String[] reaction_names)
        {
            return String.Format("{0} <-> {1} ({2:0.0###} : {3:0.0###})", reaction_names[i], reaction_names[j], lambda_i, lambda_j);
        }

    }


}


