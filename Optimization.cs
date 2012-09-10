using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Metabolism.Patterns;
using Metabolism.Analysis.Calculators;
using Metabolism.Analysis;

namespace Metabolism.Optimization
{


    public class Population
    {

        private SortedList<double, Individual> population;

        public readonly int population_size;

        public readonly double recombination_rate;


        public Population(FitnessEstimator fe, int population_size,
            double recombination_rate, double mut_rate_act, double mut_rate_deact, double start_density)
        {
            this.population_size = population_size;
            this.recombination_rate = recombination_rate;

            EvolutionParams prms = new EvolutionParams(mut_rate_act, mut_rate_deact);

            population = new SortedList<double, Individual>();
            Individual temp;
            for (int i = 0; i < population_size; ++i)
            {
                temp = Individual.RandomIndividual(fe, start_density, prms);
                population.Add(temp.Fitness, temp);
            }
        }


        public static Random rand = new Random();

        public void evolve()
        {
            Individual[] individuals = Individuals;

            Individual new_individual, temp_individual;
            for (int i = 0; i < population_size; ++i)
            {

                do
                {
                    new_individual = (Individual)individuals[i].Mutate();

                    if (rand.NextDouble() < recombination_rate)
                    {
                        temp_individual = (Individual)individuals[rand.Next(population_size)].Mutate();
                        new_individual.Recombine(temp_individual);
                    }
                } while (population.ContainsKey(new_individual.Fitness));

                population.Add(new_individual.Fitness, new_individual);
            }

            Console.WriteLine("optimal: {0} (Fitness: {1}).", population.Values.ElementAt(0), population.Keys.ElementAt(0));
        }

        public Individual[] Individuals
        {
            get
            {
                Individual[] individuals = new Individual[population_size];
                int i = 0;
                foreach (Individual individual in population.Values)
                    if (i < population_size)
                        individuals[i++] = individual;
                    else
                        break;
                return individuals;
            }
        }

        public Individual FittestIndividual
        {
            get
            {
                return population.Values.First();
            }
        }

        public int Size
        {
            get { return population_size; }
        }
    }

    public class Individual
    {

        public Individual Mutate()
        {
            Individual res = null;
            do
                try
                {
                    bool[] gene_array = new bool[fe.n];
                    foreach (int i in genes)
                        gene_array[i] = true;

                    List<int> new_genes = new List<int>();
                    for (int i = 0; i < fe.n; ++i)
                    {
                        if (gene_array[i])
                        {
                            if (prms.rand.NextDouble() >= prms.mut_rate_deact)
                                new_genes.Add(i);
                        }
                        else
                        {
                            if (!too_many_genes && prms.rand.NextDouble() < prms.mut_rate_act)
                                new_genes.Add(i);
                        }
                    }
                    res = new Individual(new_genes, fe, prms);
                }
                catch (CalculationException e) { }
            while (res == null);
            return res;
        }

        public Individual Recombine(Individual b)
        {

            Individual res = null;

            do
                try
                {
                    bool[] a_genes = new bool[fe.n];
                    foreach (int i in genes)
                        a_genes[i] = true;
                    bool[] b_genes = new bool[b.fe.n];
                    foreach (int i in b.genes)
                        b_genes[i] = true;

                    int s = prms.rand.Next(fe.n), t = prms.rand.Next(fe.n);
                    int gen;
                    for (int i = 0; i < t; ++i)
                    {
                        gen = (s + i) % fe.n;
                        a_genes[gen] ^= b_genes[gen];
                        b_genes[gen] ^= a_genes[gen];
                        a_genes[gen] ^= b_genes[gen];
                    }

                    List<int> new_genes = new List<int>();
                    for (int i = 0; i < fe.n; ++i)
                        if (a_genes[i])
                            new_genes.Add(i);

                    res = new Individual(new_genes, fe, prms);
                }
                catch (CalculationException e) { }
            while (res == null);
            return res;
        }


        private readonly List<int> genes;
        private readonly FluxPattern max;
        private readonly double fitness;
        private bool too_many_genes;

        public List<int> Genes { get { return new List<int>(genes); } }
        public double Fitness { get { return fitness; } }
        public FluxPattern Max { get { return max; } }

        private readonly FitnessEstimator fe;
        private readonly EvolutionParams prms;

        public Individual(List<int> genes, FitnessEstimator fe, EvolutionParams prms)
        {
            this.fe = fe;
            this.prms = prms;
            this.fitness = fe.CalculateFitness(genes, out this.max, out this.genes, out this.too_many_genes);
        }

        public static Individual RandomIndividual(FitnessEstimator fe, double density, EvolutionParams prms)
        {
            Individual res = null;
            do
                try
                {
                    List<int> r = new List<int>();
                    for (int i = 0; i < fe.n; ++i)
                        if (prms.rand.NextDouble() < density)
                            r.Add(i);
                    res = new Individual(r, fe, prms);
                }
                catch (CalculationException e) { }
            while (res == null);
            return res;
        }

        public override string ToString()
        {
            bool[] chosen = new bool[fe.n];
            foreach (int gen in genes)
                chosen[gen] = true;
            return String.Format("{0} are possible w/o {1}.", max, new FluxPattern(chosen, false));
        }
    }

    public struct EvolutionParams
    {
        public readonly double mut_rate_act, mut_rate_deact;
        public readonly Random rand;

        public EvolutionParams(double act, double deact)
        {
            this.mut_rate_act = act;
            this.mut_rate_deact = deact;
            this.rand = new Random();
        }
    }

    public class OptimizationTree
    {

        private int opt_max_size;
        private int opt_r_size;

        private LinkedList<List<int>> opt_genes;

        private FitnessEstimator fe;
        private List<int> T, D;
        private FluxPattern max_allowed;

        private bool output;
        private int n_updates;

        private List<int> order;

        private LinkedList<String> dot_tree;

        private static String NodeName(List<int> R)
        {
            if (R.Count < 1)
                return "";

            StringBuilder res = new StringBuilder();
            foreach (int r in R)
                res.Append(r+1).Append(',');
            return res.ToString(0, res.Length - 1);
        }

        public LinkedList<String> DotTree
        {
            get {
                LinkedList<String> res = new LinkedList<string>(this.dot_tree);

                res.AddFirst("digraph {");
                res.AddLast("}");

                return res;
            }
        }

        public OptimizationTree(FitnessEstimator fe, List<int> R0, List<int> Rb, List<int> approximation = null, FluxPattern max_approx = null, bool output = true)
        {
            this.dot_tree = new LinkedList<string>();

            this.opt_genes = new LinkedList<List<int>>();
            this.opt_max_size = fe.DoIncrease ? int.MaxValue : int.MinValue;
            this.opt_r_size = fe.DoIncrease ? int.MinValue : int.MaxValue;

            this.fe = fe;
            this.T = new List<int>(R0);
            this.D = new List<int>(Rb);

            bool[] allowed = new bool[fe.n];
            for (int i = 0; i < fe.n; ++i)
                allowed[i] = true;
            foreach (int r in R0)
                allowed[r] = false;
            this.max_allowed = new FluxPattern(allowed, false);



            if (approximation != null)
                Update(approximation, max_approx);

            this.output = output;
            this.n_updates = 0;



            List<int> reactions = new List<int>();
            for (int i = 0; i < Rb.Count; ++i)
                reactions.Add(i);
            this.order = reactions.OrderBy(x => System.Guid.NewGuid()).ToList();


            Optimize(fe.DoIncrease ? this.D : new List<int>(), 0);
        }


        private void Update(List<int> R, FluxPattern max = null)
        {
            if ((fe.DoIncrease && max.Count < this.opt_max_size) || (!fe.DoIncrease && max.Count > this.opt_max_size))
            {
                this.opt_genes.Clear();
                this.opt_genes.AddFirst(R);
                this.opt_max_size = max.Count;
                this.opt_r_size = R.Count;
            }
            else if (max.Count == this.opt_max_size)
            {
                this.opt_genes.AddLast(R);
                this.opt_r_size = Math.Min(this.opt_r_size, R.Count);
            }
        }


        private void Optimize(List<int> R, int i, FluxPattern upper_bound = null)
        {
            ++n_updates;
            String name = "";

            FluxPattern max, ub, lb;
            bool tmg;
            fe.CalculateFitness(R, out max, out R, out tmg, upper_bound);



            if (this.T.Intersect(max.Values).Count() == (fe.DoIncrease ? this.T.Count : 0))
            {
                Update(R, max);
                name = (String.Format("\"{0}\" [label = \"\\\\tiny${0}$\\n\\\\tiny\\\\newline \\n\\\\tiny${1}$\", shape = rectangle, style = \"filled\", fillcolor=\"0.4166667 0.6 1\"]", NodeName(R), (max.Count > 0) ? NodeName(max.Values) : "\\\\emptyset"));
            }
            else
            {
                name = (String.Format("\"{0}\" [label = \"\\\\tiny${0}$\\n\\\\tiny\\\\newline \\n\\\\tiny${1}$\", shape = rectangle, style = \"dotted,rounded\"]", NodeName(R), (max.Count > 0) ? NodeName(max.Values) : "\\\\emptyset"));

                if (fe.DoIncrease)
                {
                    if (max.Count < this.opt_max_size)
                    {

                        // filter reactions to check
                        bool[] queue = new bool[fe.n], coupling;
                        foreach (int l in R)
                            queue[l] = true;
                        for (int l = 0; l < i; ++l)
                            queue[l] = false;
                        foreach (int l in R) if (l < i)
                            {
                                coupling = fe.Coupling.Max(l).Incidence;
                                for (int k = 0; k < fe.n; ++k)
                                    queue[k] &= queue[k];
                            }
                        // branch & bound
                        List<int> fix = new List<int>(), next, temp;
                        for (int j = 0; j < fe.n; ++j)
                            if (queue[j])
                            {
                                fix.Clear();
                                foreach (int l in R) if (l < j - 1)
                                        fix.Add(l);
                                fe.CalculateFitness(fix, out ub, out temp, out tmg);
                                ++n_updates;

                                if (this.T.Intersect(ub.Values).Count() < this.T.Count)
                                {
                                    String break_witness = (String.Format("\"{2}\" [label = \"\\\\tiny${0}$\\n\\\\tiny\\\\newline \\n\\\\tiny${1}$\", shape = rectangle, style = \"dotted,filled\", fillcolor=\"0.08333334 0.6 1\"]", NodeName(fix), (ub.Count > 0) ? NodeName(ub.Values) : "\\\\emptyset", NodeName(R) + "-b"));
                                    this.dot_tree.AddLast(break_witness);
                                    this.dot_tree.AddLast(String.Format("\"{0}\" -> \"{1}\" [style=bold]", NodeName(R) + "-b", NodeName(R)));
                                    break;
                                }
                                else
                                {
                                    next = new List<int>();
                                    foreach (int l in R) if (l != j)
                                            next.Add(l);
                                    this.dot_tree.AddLast(String.Format("\"{0}\" -> \"{1}\"", NodeName(R), NodeName(next)));
                                    Optimize(next, j, ub);
                                }
                            }
                    }
                }
                else
                {
                    if (max.Count > this.opt_max_size)
                    {

                        // filter reactions to check
                        bool[] queue = new bool[fe.n], relaxation, coupling;
                        IEnumerable<int> inter = D.Intersect(max.Values);
                        foreach (int l in inter)
                            queue[l] = true;

                        // branch & bound
                        List<int> possible = new List<int>(), next, temp;
                        for (int j = i; j < fe.n; ++j)
                            if (queue[j])
                            {
                                relaxation = max.Incidence;
                                coupling = fe.Coupling.Max(j).Incidence;
                                for (int l = 0; l < fe.n; ++l)
                                    relaxation[l] &= coupling[l];
                                ub = new FluxPattern(relaxation, false);


                                possible.Clear();
                                foreach (int l in R)
                                    possible.Add(l);
                                foreach (int l in this.D)
                                    if (l >= j)
                                        possible.Add(l);

                                fe.CalculateFitness(possible, out lb, out temp, out tmg, ub);
                                ++n_updates;

                                if (this.T.Intersect(lb.Values).Count() > 0)
                                {
                                    String break_witness = (String.Format("\"{2}\" [label = \"\\\\tiny${0}$\\n\\\\tiny\\\\newline \\n\\\\tiny${1}$\", shape = rectangle, style = \"dotted,filled\", fillcolor=\"0.08333334 0.6 1\"]", NodeName(possible), (lb.Count > 0) ? NodeName(lb.Values) : "\\\\emptyset", NodeName(R) + "-b"));
                                    this.dot_tree.AddLast(break_witness);
                                    this.dot_tree.AddLast(String.Format("\"{0}\" -> \"{1}\" [style=bold]", NodeName(R) + "-b", NodeName(R)));
                                    break;
                                }
                                else
                                {
                                    next = new List<int>(R);
                                    next.Add(j);

                                    this.dot_tree.AddLast(String.Format("\"{0}\" -> \"{1}\"", NodeName(R), NodeName(next)));
                                    Optimize(next, j, ub);
                                }
                            }
                    }
                }
            }
            this.dot_tree.AddFirst(name);
        }

        public int OptMaxSize
        {
            get { return opt_max_size; }
        }

        public int OptRSize
        {
            get { return opt_r_size; }
        }

        public LinkedList<List<int>> OptGenes
        {
            get { return new LinkedList<List<int>>(opt_genes); }
        }

        public int NumberNodes { get { return n_updates; } }

        private int Blockable(int index)
        {
            if (index < 0 || index >= D.Count)
                return -1;
            return D[order[index]];
        }

    }

    public class FitnessEstimator
    {

        private readonly List<int> R0, Rb;
        private readonly bool increase;

        private readonly Random rand;

        private readonly IIntFCACalculator<FluxPattern> fca_solver;
        public readonly int n;
        private readonly LinkedList<FluxPattern> witnesses;
        private LinkedList<FluxPattern> new_witnesses;
        private readonly FluxPattern MAX, FREV;
        private readonly IIntCoupling coupling;

        public FitnessEstimator(List<int> R0, List<int> Rb, bool increase, int n, IIntFCACalculator<FluxPattern> fca_solver, LinkedList<FluxPattern> witnesses, FluxPattern max, FluxPattern frev, IIntCoupling coupling)
        {
            this.R0 = R0;
            this.Rb = Rb;
            this.increase = increase;
            this.rand = new Random();

            this.n = n;
            this.fca_solver = fca_solver;
            this.witnesses = witnesses;
            this.MAX = max;
            this.FREV = frev;
            this.coupling = coupling;
        }

        public IIntCoupling Coupling
        {
            get { return this.coupling; }
        }

        public bool DoIncrease
        {
            get { return this.increase; }
        }

        public double CalculateFitness(List<int> R, out FluxPattern max, out List<int> R_trimmed, out bool too_many_genes, FluxPattern ub = null)
        {

            R_trimmed = R = new List<int>(R.Intersect(Rb));

            FluxPattern calcs;
            max = fca_solver.CalculateMax(R, n, witnesses, out new_witnesses, out calcs, ub == null ? MAX : ub, FREV, coupling);
            too_many_genes = max.Count == 0;

            foreach (FluxPattern a in new_witnesses)
                if (rand.NextDouble() < ((double)a.Count) / witnesses.Count)
                    witnesses.AddLast(a);

            return R.Count + (increase ? (n + 1) : -(n + 1)) * (max.Count - (n + 1) * R0.Intersect(max.Values).Count()) + 0.1 * rand.NextDouble();
        }
    }
}