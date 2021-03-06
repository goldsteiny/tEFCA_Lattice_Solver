﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Metabolism.Analysis.Calculators;
using Metabolism.Patterns;

namespace Metabolism.Analysis
{

    /// <summary>
    /// An interface for flux couplings based on lists of the indices of coupled reactions.
    /// </summary>
    public interface IIntCoupling
    {
        /// <summary>
        /// $Coupling_i = \{j \in [n] | i \rightarrow j\}$
        /// </summary>
        /// <returns>A list of all reactions j coupled to i $i \rightarrow j$(especially $i \rightarrow i$)</returns>
        ICollection<int> this[int index] { get; }

        /// <summary>
        /// n (the number of reactions)
        /// </summary>
        int Length { get; }

        /// <summary>
        /// A subset of the flux lattice L we used to decide which reactions are coupled.
        /// </summary>
        LinkedList<FluxPattern> Witnesses { get; }

        /// <summary>
        /// Each finite join semi-lattice has a unique maximum. So does $L_i$, the flux lattice with reaction i blocked. The maximum of $L_i$ contains all the unblocked reactions not coupled to i.
        /// </summary>
        /// <returns>$1_{L_i}$</returns>
        FluxPattern Max(int index);


        /// <summary>
        /// Converts the index-based flux coupling to a coupling based on the coupled reaction pairs.
        /// </summary>
        /// <returns>An object implementing the reaction-pair based ICoupling interface.</returns>
        ICoupling ToCoupling();


        /// <summary>
        /// Creates a Graphviz-Dot-File of the coupling relations that can be procedeeded using Graphviz
        /// </summary>
        /// <returns>List of lines of the .dot-file</returns>
        LinkedList<String> ToDot(String[] names);
    }


    /// <summary>
    /// An index based flux coupling class.
    /// </summary>
    public class IntCoupling : IIntCoupling
    {

        private readonly int n;
        private readonly LinkedList<int>[] coupled;
        private readonly LinkedList<FluxPattern> witnesses;
        private readonly FluxPattern[] maxima;
        private readonly FluxPattern max;

        /// <summary>
        /// Does the FCA for the given network and saves the result in a new IntCoupling object.
        /// </summary>
        /// <param name="n">#reactions in the network.</param>
        /// <param name="calculator">An object that calculates the maximal set of reactions in the metabolic network constrained by disabled reactions.</param>
        /// <param name="witnesses">A subset of the flux lattice L.</param>
        /// <param name="max">The maximum of L / the set of all unblocked reactions in the metabolic network.</param>
        /// <param name="frev">The set of all fully reversible reactions in the metabolic network.</param>
        /// <param name="R">A set $R \subseteq [n]$ of disabled reactions.</param>
        /// <param name="coupling">Flux coupling we already know about.</param>
        /// <param name="start">This constructor calculated the flux coupling sets $C_i$ for $i \geq start$.</param>
        public IntCoupling(IIntFCACalculator<FluxPattern> calculator, bool show_output, ICollection<FluxPattern> witnesses, out DateTime[] time_reaction, out int[] lps_reaction, out int[] wits_reaction, out int wits_used, FluxPattern max, FluxPattern frev = null, List<int> R = null, IIntCoupling coupling = null, int start = 0)
        {
            this.max = max;
            this.n = max.Length;
            this.witnesses = new LinkedList<FluxPattern>();

            time_reaction = new DateTime[n];
            lps_reaction = new int[n];
            wits_reaction = new int[n];
            int lps_count = calculator.SolverCalls;


            this.maxima = new FluxPattern[n];
            coupled = new LinkedList<int>[n];
            for (int i = 0; i < n; ++i)
            {
                coupled[i] = new LinkedList<int>();
                coupled[i].AddLast(i);
            }

            if (R == null)
                R = new List<int>();


            LinkedList<FluxPattern> new_witnesses = new LinkedList<FluxPattern>(witnesses);
            witnesses = new LinkedList<FluxPattern>();
            bool allowed;
            foreach (FluxPattern w in new_witnesses)
            {
                allowed = true;
                for (int i = 0; allowed && i < R.Count; ++i)
                    allowed = !w[R[i]];
                if (allowed)
                    witnesses.Add(w);
            }
            wits_used = witnesses.Count;


            if (coupling == null)
                coupling = this;

            FluxPattern opt, calcs;
            R.Add(-1);
            ICollection<int> c;
            for (int i = start; i < n; ++i)
                if (maxima[i] == null)
                {
                    if (max[i])
                    {


                        R[R.Count - 1] = i;

                        opt = calculator.CalculateMax(R, n, witnesses, out new_witnesses, out calcs, max, frev, coupling);
                        lps_reaction[i] = calculator.SolverCalls - lps_count;
                        lps_count = calculator.SolverCalls;

                        c = coupling[i];
                        foreach (int r in c)
                            if (r >= i && coupling.Max(r) != null && !coupling.Max(r)[i])
                            {
                                maxima[r] = opt;
                                for (int j = 0; j < n; ++j)
                                    if (max[j] && !maxima[r][j] && r != j)
                                        coupled[r].AddLast(j);
                            }

                        if (maxima[i] == null)
                        {
                            maxima[i] = opt;
                            if (maxima[i] == null)
                                Console.WriteLine("Fehler!");
                            else
                                for (int j = 0; j < n; ++j)
                                    if (!maxima[i][j] && max[j] && i != j)
                                        coupled[i].AddLast(j);
                        }

                        foreach (FluxPattern a in new_witnesses)
                            if (a.Count > 0)
                            {
                                this.witnesses.AddLast(a);
                                witnesses.Add(a);
                                wits_reaction[i]++;
                            }

                        time_reaction[i] = DateTime.Now;
                        if (show_output)
                            Console.WriteLine("({0})\tCoupling {1} from {2} calculated. ({3} LPs solved, {4} kept.)\n", time_reaction[i], i + 1, n, lps_reaction[i], witnesses.Count);

                    }
                    else
                        maxima[i] = max;
                }
            R.RemoveAt(R.Count - 1);
        }

        public ICollection<int> this[int index]
        {
            get { return coupled[index]; }
        }

        public int Length
        {
            get { return n; }
        }

        public LinkedList<FluxPattern> Witnesses
        {
            get { return this.witnesses; }
        }


        public FluxPattern Max(int index)
        {
            return maxima[index];
        }


        public ICoupling ToCoupling()
        {
            LinkedListCoupling res = new LinkedListCoupling(n);

            ICollection<int> coupled;
            for (int i = 0; i < n; ++i)
                if (this.max[i])
                {
                    coupled = this.coupled[i];
                    foreach (int j in coupled)
                        if (this.max[j])
                        {
                            if (!maxima[j][i])
                            {
                                if (i < j)
                                    res.Add(new PartialCouple(i, j));
                            }
                            else
                                res.Add(new Couple(i, j));
                        }
                }

            return res;
        }



        public LinkedList<string> ToDot(String[] reaction_names)
        {
            LinkedList<String> res = new LinkedList<string>();

            res.AddFirst("digraph FCA {");
            res.AddLast("\tnodesep=0.7;");
            res.AddLast("\trankdir=LR;");

            if (reaction_names != null)
                for (int i = 0; i < n; ++i)
                    res.AddLast(String.Format("\t{0}[label=\"{1}\"];", i + 1, reaction_names[i]));

            // compress
            int[] component = new int[n];
            List<int>[] components = new List<int>[n];
            for (int i = 0; i < n; ++i)
            {
                component[i] = -1;
                components[i] = new List<int>();
            }

            for (int i = 0; i < n; ++i)
                if (component[i] < 0)
                {
                    component[i] = i;
                    components[i].Add(i);
                    for (int j = i + 1; j < n; ++j)
                        if (!Max(i)[j] && !Max(j)[i])
                        {
                            component[j] = i;
                            components[i].Add(j);
                        }
                }

            // add components
            for (int i = 0; i < n; ++i)
                if (component[i] == i)
                    switch (components[i].Count)
                    {
                        case 0:
                        case 1:
                            break;
                        case 2:
                            res.AddLast(String.Format("\t{0}->{1}[dir=both];", components[i][0] + 1, components[i][1] + 1));
                            break;
                        default:
                            for (int j = 0; j < components[i].Count; ++j)
                                res.AddLast(String.Format("\t{0}->{1}[dir=both];", components[i][j] + 1, components[i][(j + 1) % components[i].Count] + 1));
                            break;
                    }

            // connect components
            for (int i = 0; i < n; ++i)
                if (component[i] == i)
                    for (int j = 0; j < n; ++j)
                        if (!Max(i)[j] && component[j] == j && i != j)
                            res.AddLast(String.Format("\t{0}->{1};", i + 1, j + 1));

            res.AddLast("}");

            return res;
        }
    }


    /// <summary>
    /// An interface for classes that save results of an EFCA (double reaction knock-out simulation).
    /// </summary>
    public interface IEFCACount
    {

        /// <summary>
        /// An overview how many reactions are unblocked if reactions are pairwise diabled.
        /// </summary>
        /// <returns>#unblocked reactions if i and j are disabled.</returns>
        int this[int i, int j] { get; }

        /// <summary>
        /// #reactions in the metabolic network.
        /// </summary>
        int Length { get; }

        /// <summary>
        /// Number of non-trivial lattice elements we found doing the EFCA.
        /// </summary>
        int WitnessCount { get; }

        /// <summary>
        /// Number of (MI)LPs we had to solve to do the EFCA.
        /// </summary>
        int LPCount { get; }

        /// <summary>
        /// Converts the results of the EFCA into the rows of a CSV file.
        /// </summary>
        /// <returns>A list of the rows of the CSV file.</returns>
        LinkedList<string> ToCSV(string separator, string[] names = null);
    }


    /// <summary>
    /// An enhanced flux coupling analysis simulates all $\choose n 2$ double reaction knock-outs for a metabolic network. For each pair of disabled reactions we calculate the set of still unblocked reactions. An EFCACount object saves the number of those unblocked reactions foreach pair (i,j).
    /// </summary>
    public class EFCACount : IEFCACount
    {

        private readonly int[,] max_size, targeted_size;
        private readonly FluxPattern max;
        private readonly List<int> T; 
        public readonly int n;
        private int lp_counter = 0, witness_counter = 0;
        private readonly LinkedList<FluxPattern> witnesses_kept;
        private readonly Random rand = new Random();



        /// <summary>
        /// Simulates all $\choose n 2$ double reaction knock-outs and creates a new EFCACount object for the results.
        /// </summary>
        /// <param name="calculator">An object that calculates the maximal set of reactions in the metabolic network constrained by disabled reactions.</param>
        /// <param name="witnesses">A subset of the flux lattice L.</param>
        /// <param name="max">The maximum of L / the set of all unblocked reactions in the metabolic network.</param>
        /// <param name="frev">The set of all fully reversible reactions in the metabolic network.</param>
        /// <param name="coupling">A flux coupling for the network</param>
        public EFCACount(IIntFCACalculator<FluxPattern> calculator, bool keep_witnesses, ICollection<FluxPattern> witnesses, FluxPattern max, FluxPattern frev, out DateTime[] time_reaction, out int[] lps_reaction, out int[] wits_reaction, out int[] wits_usable, IIntCoupling coupling = null, ICollection<int> T = null, bool show_output = true)
        {
            this.max = max;
            this.n = max.Length;

            time_reaction = new DateTime[n];
            lps_reaction = new int[n];
            for (int i = 0; i < n; ++i)
                lps_reaction[i] = -1;
            wits_reaction = new int[n];
            wits_usable = new int[n];

            DateTime[] time;
            int[] lps, wits;
            int temp_counter = calculator.SolverCalls, temp_wits;
            if (coupling == null)
            {
                coupling = new IntCoupling(calculator, false, witnesses, out time, out lps, out wits, out temp_wits, max, frev);
                this.witness_counter = coupling.Witnesses.Count;
                this.lp_counter = calculator.SolverCalls - temp_counter;
                temp_counter = calculator.SolverCalls;
            }

            if (T == null)
            {
                this.T = new List<int>(n);
                for (int i = 0; i < n; ++i)
                    this.T.Add(i);
            }
            else
                this.T = new List<int>(T);

            LinkedList<FluxPattern> new_witnesses = new LinkedList<FluxPattern>();
            witnesses_kept = new LinkedList<FluxPattern>(witnesses);
            IIntCoupling row;

            max_size = new int[n, n];
            targeted_size = new int[n, n];
            for (int i = 0; i < n; ++i)
            {
                max_size[i, i] = -1;
                targeted_size[i, i] = -1;
            }

            List<int> R = new List<int>();
            ICollection<int> c;
            int wits_used;
            R.Add(-1);
            for (int i = 0; i < n; ++i)
                if (max[i] && max_size[i, i] < 0)
                {
                    R[0] = i;


                    if (keep_witnesses)
                    {
                        foreach (FluxPattern a in new_witnesses)
                            if (rand.NextDouble() < ((double)a.Count) / witnesses_kept.Count)
                                witnesses_kept.AddLast(a);
                    }

                    row = new IntCoupling(calculator, false, witnesses_kept, out time, out lps, out wits, out wits_used, max, frev, R, coupling, i);
                    wits_usable[i] = wits_used;
                    if (show_output)
                        Console.WriteLine("({0})\tRow {1} from {2} calculated of EFCA. ({3} LPs solved, {4} kept.)\n", System.DateTime.Now, i + 1, n, calculator.SolverCalls - temp_counter, witnesses_kept.Count);

                    this.lp_counter += calculator.SolverCalls - temp_counter;
                    temp_counter = calculator.SolverCalls;
                    new_witnesses = row.Witnesses;
                    this.witness_counter += row.Witnesses.Count;


                    c = coupling[i];
                    foreach (int r in c)
                        if (r >= i && !coupling.Max(r)[i])
                            for (int j = r; j < n; ++j)
                                if (max[j])
                                {
                                    max_size[r, j] = row.Max(j).Count;
                                    max_size[j, r] = max_size[r, j];

                                    targeted_size[r, j] = row.Max(j).Values.Intersect(this.T).Count();
                                    targeted_size[j, r] = targeted_size[r, j];
                                }

                    lps_reaction[i] = lp_counter;
                    wits_reaction[i] = witness_counter;
                    time_reaction[i] = DateTime.Now;
                }

            if (show_output)
                Console.WriteLine("({0})\tEFCA uncompressed.\n", System.DateTime.Now);
        }

        public int this[int i, int j]
        {
            get { return max_size[i, j]; }
        }

        public int Length
        {
            get { return n; }
        }

        public int WitnessCount
        {
            get { return witness_counter; }
        }

        public int LPCount
        {
            get { return lp_counter; }
        }


        public LinkedList<string> ToCSV(string separator, string[] names = null)
        {

            if (names == null)
            {
                names = new string[n];
                for (int i = 0; i < n; ++i)
                    names[i] = (i + 1).ToString();
            }

            LinkedList<string> res = new LinkedList<string>();
            res.AddFirst(String.Format("i{0}j{0}Number.of.maximum{0}Number.of.targeted", separator));

            for (int i = 0; i < n; ++i)
                if (this.max[i])
                    for (int j = 0; j < n; ++j)
                        if (this.max[j])
                            res.AddLast(String.Format("{1}{0}{2}{0}{3}{0}{4}", separator, names[i], names[j], this[i, j], this.targeted_size[i, j]));

            return res;
        }
    }
}