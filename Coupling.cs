using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Metabolism.Patterns;

namespace Metabolism.Analysis
{

    /// <summary>
    /// An interface that allows you to implement classes that represents a set of properties a certain lattice satisfies. These properties implement the INoIImpliesProperty interface, e.g. flux coupling rules for the lattice.
    /// </summary>
    public interface ILatticeNoIImpliessPropertySet<P> where P : INoIImpliesProperty
    {
        LinkedList<string> ToStrings();
        LinkedList<string> ToStrings(String[] reaction_names);
        void Add(P p);
        int Count { get; }

        /// <summary>
        /// Returns the set of all properties.
        /// </summary>
        /// <returns>A list of all properties.</returns>
        List<P> ToList();

        /// <summary>
        /// An implementation of this interfaces is a set of properties all lattice members have to fulfill if they don't contain the number "i". This returns all properties for a certain "i".
        /// </summary>
        /// <returns>$\{p Property | index \notin a \Rightarrow a fulfills p \forall a in L\}$</returns>
        List<P> this[int index] { get; }

        /// <summary>
        /// The number of properties implied by the absence of "i".
        /// </summary>
        /// <returns>The number of properties implied by the absence of "i".</returns>
        int GetCount(int i);
    }

    /// <summary>
    /// An interface for flux couplings (sets of all coupled pairs of reactions in a flux lattice).
    /// </summary>
    public interface ICoupling : ILatticeNoIImpliessPropertySet<Couple> { }


    /// <summary>
    /// A flux coupling class whose sets of properties linked to the index "i" are saved as LinkedLists. This allows fast adding, removing and filtering of elements.
    /// </summary>
    public class LinkedListCoupling : LatticeNoIImpliesPropertySet<Couple>, ICoupling
    {

        /// <summary>
        /// Creates a new LinkedListCoupling object for flux lattices $L \subseteq 2^{[n]}$.
        /// </summary>
        /// <param name="n">The number of reactions.</param>
        public LinkedListCoupling(int n) : base(n) { }
    }


    /// <summary>
    /// An implementation of the ILatticeNoIImpliessPropertySet interface where all the sets of properties linked to an index "i" are saved as LinkedLists.
    /// </summary>
    /// <typeparam name="P"></typeparam>
    public class LatticeNoIImpliesPropertySet<P> : ILatticeNoIImpliessPropertySet<P> where P : INoIImpliesProperty
    {
        private readonly LinkedList<P>[] couples;
        private readonly int n;
        private int counter;

        /// <summary>
        /// Create a new instance of the LatticeNoIImpliesPropertySet class.
        /// </summary>
        /// <param name="n">The number of indices important for the lattice properties.</param>
        public LatticeNoIImpliesPropertySet(int n)
        {
            this.n = n;

            couples = new LinkedList<P>[n];
            for (int i = 0; i < n; ++i)
                couples[i] = new LinkedList<P>();

            counter = 0;
        }


        /// <summary>
        /// Removes the first of all properties saved in the class (lexicographical orderd by the index).
        /// </summary>
        /// <returns>The property removed from the set.</returns>
        public P RemoveFirst()
        {

            foreach (LinkedList<P> line in couples)
                if (line.Count > 0)
                {
                    P res = line.First.Value;
                    line.RemoveFirst();
                    --counter;
                    return res;
                }

            return default(P);

        }

        public LinkedList<string> ToStrings()
        {
            LinkedList<string> res = new LinkedList<string>();

            SortedSet<P> couples = Pairs;
            foreach (P c in couples)
                res.AddLast(c.ToString());

            return res;
        }

        public SortedSet<P> Pairs
        {
            get
            {
                SortedSet<P> res = new SortedSet<P>();
                foreach (ICollection<P> row in couples)
                    foreach (P c in row)
                        res.Add(c);
                return res;
            }
        }

        public void Add(P item)
        {
            couples[item.I].AddLast(item);
            ++counter;
        }


        public int Count
        {
            get { return counter; }
        }


        /// <summary>
        /// Removes all properties which are not fulfilled by the lattice because they don't hold for "witness". 
        /// </summary>
        /// <param name="witness">Our witness for testing properties.</param>
        /// <returns>true iff properties have been removed.</returns>
        public bool RemoveAll<FP>(FP witness) where FP : IFluxPattern<FP>
        {
            if (witness.Count == 0)
                return false;

            int counter_before = counter;
            LinkedList<P> line;
            for (int i = 0; i < n; ++i)
                if (!witness[i])
                {
                    line = new LinkedList<P>();
                    foreach (P p in couples[i])
                        if (!p.HoldsFor(witness))
                            --counter;
                        else
                            line.AddLast(p);
                    couples[i] = line;
                }
            return counter_before > counter;
        }



        public List<P> ToList()
        {
            return Pairs.ToList();
        }


        public List<P> this[int index]
        {
            get
            {
                return this.couples[index].ToList();
            }
        }

        public int GetCount(int i)
        {
            return this.couples[i].Count;
        }


        public LinkedList<string> ToStrings(String[] reaction_names)
        {
            LinkedList<string> res = new LinkedList<string>();

            SortedSet<P> couples = Pairs;
            foreach (P c in couples){
                res.AddLast(c.ToString(reaction_names));
            }
            return res;
        }
    }
}