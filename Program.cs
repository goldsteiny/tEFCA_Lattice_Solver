using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Metabolism.Analysis;
using Metabolism.Patterns;
using Metabolism.Analysis.Calculators;
using Tools.IO;
using Tools.Converter.Biology;
using Metabolism.Optimization;

namespace FCA
{
    public class FCA
    {

        public static void Analyze(double[,] matrix, bool[] rev, bool[] trans, LinkedList<FluxPattern> forbidden_combinations, double max_value, double tolerance, bool lp, bool fc, bool do_efca, bool keep_witnesses, bool only_internal, bool show_output, out FluxPattern max, out IIntCoupling coupling, out EFCACount efca, out int reactions_unblocked, out int reactions_frev, out int couples, out ICollection<FluxPattern> circuits, out DateTime time_start, out DateTime time_unblocked, out DateTime time_frev, out DateTime time_couples, out DateTime time_efca, out DateTime[] time_reaction, out int lps_unblocked, out int lps_frev, out int lps_couples, out int lps_efca, out int[] lps_reaction, out int wits_unblocked, out int wits_frev, out int wits_couples, out int wits_efca, out int[] wits_reaction, out int[] wits_distribution, out int[] wits_usable, String name = null)
        {


            //initialisation
            time_start = DateTime.Now;
            efca = null;
            wits_usable = null;
            int n = rev.Length;

            List<int> irr = new List<int>();
            for (int i = 0; i < n; ++i)
                if (!rev[i])
                    irr.Add(i);

            LinkedList<FluxPattern> witnesses = new LinkedList<FluxPattern>();
            LinkedList<FluxPattern> blocking_witnesses, frev_witnesses;
            IIntFCACalculator<FluxPattern> fcacalculator = !lp ? (IIntFCACalculator<FluxPattern>)new MILPIntFCACalculator(matrix, rev, max_value, tolerance) :
                !fc ? (IIntFCACalculator<FluxPattern>)new LPIntFCACalculator(matrix, rev, max_value, tolerance) :
                forbidden_combinations == null || forbidden_combinations.Count == 0 ? (IIntFCACalculator<FluxPattern>)new FCIntFCACalculator(matrix, rev, trans, max_value, tolerance) :
                (IIntFCACalculator<FluxPattern>)new FastFCIntFCACalculator(matrix, rev, forbidden_combinations, max_value, tolerance);

            circuits = (fc && forbidden_combinations == null) ? ((FCIntFCACalculator) fcacalculator).Circles : null;

            // calculate blocked reactions
            FluxPattern calcs;
            int wits_used;
            
            bool[] internal_reactions = new bool[n];
            for(int i = 0; i<n; ++i)
                internal_reactions[i] = !trans[i];
            FluxPattern pattern_internal = new FluxPattern(internal_reactions, false);
            max = fcacalculator.CalculateMax(new List<int>(), n, new LinkedList<FluxPattern>(), out blocking_witnesses, out calcs, only_internal ? pattern_internal : null, null);
            reactions_unblocked = max.Count;



            foreach (FluxPattern a in blocking_witnesses)
                if (a.Count > 0)
                    witnesses.AddLast(a);
            if (name != null)
                SaveStatistics(witnesses, n, name + "_max.stats");
            time_unblocked = DateTime.Now;
            lps_unblocked = fcacalculator.SolverCalls;
            wits_unblocked = witnesses.Count;
            if (show_output)
                Console.WriteLine("({0})\tBlocked reactions calculated: {1} of {2} blocked.\n", time_unblocked, n-reactions_unblocked, n);


            // calculate fully reversible reactions
            FluxPattern frev = fcacalculator.CalculateMax(irr, n, blocking_witnesses, out frev_witnesses, out calcs, max, null);
            reactions_frev = frev.Count;

            foreach (FluxPattern a in frev_witnesses)
                if (a.Count > 0)
                    witnesses.AddLast(a);
            if (name != null)
                SaveStatistics(witnesses, n, name + "_frev.stats");
            time_frev = DateTime.Now;
            lps_frev = fcacalculator.SolverCalls - lps_unblocked;
            wits_frev = witnesses.Count - wits_unblocked;
            if (show_output)
                Console.WriteLine("({0})\tFRev calculated.\n", time_frev);


            // FCA
            Console.WriteLine("Witnesses so far:");
            foreach (FluxPattern a in witnesses)
                Console.WriteLine("\t{0}", a);
            Console.WriteLine("\nStarting FCA.\n");


            coupling = new IntCoupling(fcacalculator, show_output, witnesses, out time_reaction, out lps_reaction, out wits_reaction, out wits_used, max, frev);

            if (name != null)
                SaveStatistics(coupling.Witnesses, n, name + "_fca.stats");

            time_couples = DateTime.Now;
            wits_couples = coupling.Witnesses.Count - wits_frev - wits_unblocked;
            lps_couples = 0; wits_couples = 0;
            for (int i = 0; i < n; ++i)
            {
                lps_couples += lps_reaction[i];
                wits_couples += wits_reaction[i];
            }
            if (show_output)
                Console.WriteLine("({0})\tCouples calculated.\n", time_couples);

            wits_distribution = new int[n + 1];
            foreach (FluxPattern a in coupling.Witnesses)
                wits_distribution[a.Count]++;


            // EFCA
            lps_efca = 0;
            wits_efca = 0;
            if (do_efca)
            {
                efca = new EFCACount(fcacalculator, keep_witnesses, coupling.Witnesses, max, frev, out time_reaction, out lps_reaction, out wits_reaction, out wits_usable, coupling, null, show_output);
                lps_efca = efca.LPCount;
                wits_efca = efca.WitnessCount;
            }
            time_efca = DateTime.Now;

            couples = coupling.ToCoupling().Count;
        }

        public static void SaveStatistics(ICollection<FluxPattern> witnesses, int n, String path)
        {
            int[,] coupling = new int[n, n];
            int[] reactions = new int[n];

            foreach (FluxPattern a in witnesses)
                for (int i = 0; i < n; ++i)
                    if (a[i])
                    {
                        ++reactions[i];
                        for (int j = 0; j < n; ++j)
                            if (!a[j])
                                coupling[i, j]++;
                    }

            int k = witnesses.Count;

            LinkedList<String> file = new LinkedList<string>();
            file.AddLast(String.Format("{0} {1} {2}", -1, -1, k));

            for (int i = 0; i < n; ++i)
            {
                file.AddLast(String.Format("{0} {1} {2}", i + 1, -1, k - reactions[i]));
                file.AddLast(String.Format("{0} {1} {2}", -1, i + 1, reactions[i]));

                for (int j = 0; j < n; ++j)
                {
                    file.AddLast(String.Format("{0} {1} {2}", j + 1, i + 1, coupling[i, j]));
                }
            }

            Writer.WriteFile(file, path);
        }

        public static IIntCoupling DoFCA(double[,] matrix, bool[] rev, bool[] trans, LinkedList<FluxPattern> forbidden_combinations, double max_value, double tolerance, bool lp, bool fc, bool do_efca, bool keep_witnesses, bool only_internal, bool show_output, bool save_stats, String name, String model, String method, LinkedList<string> res, LinkedList<string> details, out FluxPattern max, out EFCACount efca, out ICollection<FluxPattern> circuits)
        {
            LinkedList<FluxPattern> witnesses;
            int n_unblocked;
            FluxPattern frev = null;
            IIntCoupling coupling = null;
            DateTime time_start, time_stop;
            DateTime[] time_reaction;

            int reactions_unblocked, reactions_frev, couples, lps, lps_unblocked, lps_frev, lps_couples, lps_efca, wits, wits_unblocked, wits_frev, wits_couples, wits_efca;
            int lps_total, wits_total;
            int[] lps_reaction, wits_reaction, wits_distribution, wits_usable;
            DateTime time_unblocked, time_frev, time_couples, time_efca;

            Analyze(matrix, rev, trans, forbidden_combinations, max_value, tolerance, lp, fc, do_efca, keep_witnesses, only_internal, show_output, out max, out coupling, out efca, out reactions_unblocked, out reactions_frev, out couples, out circuits, out time_start, out time_unblocked, out time_frev, out time_couples, out time_efca, out time_reaction, out lps_unblocked, out lps_frev, out lps_couples, out lps_efca, out lps_reaction, out wits_unblocked, out wits_frev, out wits_couples, out wits_efca, out wits_reaction, out wits_distribution, out wits_usable, save_stats ? name : null);

            lps = lps_unblocked + lps_frev + lps_couples + lps_efca;
            wits = wits_unblocked + wits_frev + wits_couples + wits_efca;

            /*String cellcolor = "\\cellcolor{gray!10}";
                        res.AddLast("\\midrule\n\\multirowbt{4}{*}{\\Var{" + model.Replace("_", "\\_") + "}}" + String.Format("& {4}Total &{4}{0}&{4}{2}&{4}{3}&{4}{1:0.0}\\\\", rev.Length, (time_efca - time_start).TotalSeconds, lps, wits, ""));
                        res.AddLast(String.Format("\\cmidrule{4} & {5}$1_L$ &{5}{0}&{5}{2}&{5}{3}&{5}{1:0.0}\\\\", reactions_unblocked, (time_unblocked - time_start).TotalSeconds, lps_unblocked, wits_unblocked, "{2-6}", cellcolor));
                        res.AddLast(String.Format("\\cmidrule{4} & $\\Frev$ &{0}&{2}&{3}&{1:0.0}\\\\", reactions_frev, (time_frev - time_unblocked).TotalSeconds, lps_frev, wits_frev, "{2-6}"));
                        res.AddLast(String.Format("\\cmidrule{4} & {5}$\\Coupling$ &{5}{0}&{5}{2}&{5}{3}&{5}{1:0.0}\\\\", couples, (time_couples - time_frev).TotalSeconds, lps_couples, wits_couples, "{2-6}", cellcolor));
            if(do_efca)
                res.AddLast(String.Format("\\cmidrule{4} & EFCA &{1}&{2}&{3}&{1:0.0}\\\\", efca.Length, (time_efca - time_couples).TotalSeconds, lps_efca, wits_efca, "{2-6}"));
              */  
            res.AddLast(String.Format("{0} {1} {2} {3} {4} {5} {6} {7:0.000} {8:0.000} {9:0.000} {10:0.000} {11:0.000} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21}", model, method, matrix.Length / rev.Length, rev.Length, reactions_unblocked, reactions_frev, couples, (time_efca - time_start).TotalSeconds, (time_unblocked - time_start).TotalSeconds, (time_frev - time_unblocked).TotalSeconds, (time_couples - time_frev).TotalSeconds, (time_efca - time_couples).TotalSeconds, lps, lps_unblocked, lps_frev, lps_couples, lps_efca, wits, wits_unblocked, wits_frev, wits_couples, wits_efca));

            //unblocked
            details.AddLast(String.Format("{0} {1} {2} {3} {4} {5} {6} {22:0.000} {23} {24} {25} {26} {27}", model, method, matrix.Length / rev.Length, rev.Length, reactions_unblocked, reactions_frev, couples, (time_efca - time_start).TotalSeconds, (time_unblocked - time_start).TotalSeconds, (time_frev - time_unblocked).TotalSeconds, (time_couples - time_frev).TotalSeconds, (time_efca - time_couples).TotalSeconds, lps, lps_unblocked, lps_frev, lps_couples, lps_efca, wits, wits_unblocked, wits_frev, wits_couples, wits_efca, time_start, -1, lps_unblocked, lps_unblocked, wits_unblocked, wits_unblocked));
            //frev
            details.AddLast(String.Format("{0} {1} {2} {3} {4} {5} {6} {22:0.000} {23} {24} {25} {26} {27}", model, method, matrix.Length / rev.Length, rev.Length, reactions_unblocked, reactions_frev, couples, (time_efca - time_start).TotalSeconds, (time_unblocked - time_start).TotalSeconds, (time_frev - time_unblocked).TotalSeconds, (time_couples - time_frev).TotalSeconds, (time_efca - time_couples).TotalSeconds, lps, lps_unblocked, lps_frev, lps_couples, lps_efca, wits, wits_unblocked, wits_frev, wits_couples, wits_efca, time_start, 0, lps_frev, lps_unblocked + lps_frev, wits_frev, wits_unblocked + wits_frev));
            //couples
            lps_total = lps_unblocked + lps_frev;
            wits_total = wits_unblocked + wits_frev;
            int temp_i = 0;
            for (int i = 0; i < rev.Length; ++i)
                if (max[i] && lps_reaction[i] > -1)
                {
                    lps_total += lps_reaction[i];
                    wits_total += wits_reaction[i];
                    temp_i = i;
                    details.AddLast(String.Format("{0} {1} {2} {3} {4} {5} {6} {22:0.000} {23} {24} {25} {26} {27} {28} {29}", model, method, matrix.Length / rev.Length, rev.Length, reactions_unblocked, reactions_frev, couples, do_efca ? (time_efca - time_start).TotalSeconds : 0, (time_unblocked - time_start).TotalSeconds, (time_frev - time_unblocked).TotalSeconds, (time_couples - time_frev).TotalSeconds, do_efca ? (time_efca - time_couples).TotalSeconds : 0, lps, lps_unblocked, lps_frev, lps_couples, lps_efca, wits, wits_unblocked, wits_frev, wits_couples, wits_efca, i + 1, (time_reaction[i] - time_start).TotalSeconds, lps_reaction[i], lps_total, wits_reaction[i], wits_total, wits_distribution[i + 1], do_efca ? wits_usable[i] : 0));
                }
                else
                    details.AddLast(String.Format("{0} {1} {2} {3} {4} {5} {6} {22:0.000} {23} {24} {25} {26} {27} {28} {29}", model, method, matrix.Length / rev.Length, rev.Length, reactions_unblocked, reactions_frev, couples, do_efca ? (time_efca - time_start).TotalSeconds : 0, (time_unblocked - time_start).TotalSeconds, (time_frev - time_unblocked).TotalSeconds, (time_couples - time_frev).TotalSeconds, do_efca ? (time_efca - time_couples).TotalSeconds : 0, lps, lps_unblocked, lps_frev, lps_couples, lps_efca, wits, wits_unblocked, wits_frev, wits_couples, wits_efca, i + 1, (time_reaction[temp_i] - time_start).TotalSeconds, 0, lps_total, 0, wits_total, wits_distribution[i + 1], do_efca ? wits_usable[i] : 0));
            return coupling;
        }

        public static bool[] FindTrans(double[,] matrix, bool[] rev)
        {
            bool[] trans = new bool[rev.Length];
            int pos, neg;
            for (int i = 0; i < rev.Length; ++i)
            {
                pos = 0; neg = 0;
                for (int j = 0; j < matrix.Length / rev.Length; ++j)
                {
                    if (matrix[j, i] > 0)
                        ++pos;
                    if (matrix[j, i] < 0)
                        ++neg;
                }
                trans[i] = (pos * neg == 0) && (pos + neg > 0);
                //trans[i] = false;
            }
            Console.WriteLine("The exchange reactions are:\n\t{0}", new FluxPattern(trans, false));
            return trans;
        }

        public static void Main(string[] args)
        {

            //init
            string path2models = args[0], name, separator = args[1];

            bool fc = true, calculate_circuits = false, do_efca = args.Length > 5 && int.Parse(args[5]) == 1, no_output = args.Length > 4 && int.Parse(args[4]) == 0, save_stats = args.Length > 7 && int.Parse(args[7]) > 0;
            bool only_internal = false, keep_witnesses = false;

            int n_simulations = args.Length > 6 ? int.Parse(args[6]) : 1;
            double max_value = double.Parse(args[2]), tolerance = double.Parse(args[3]);
            LinkedList<FluxPattern> forbidden = new LinkedList<FluxPattern>();


            double[,] matrix; bool[] rev;

            LinkedList<string> paths = Reader.ReadFile(path2models), res = new LinkedList<string>(), details = new LinkedList<string>();
//            res.AddFirst("\\begin{table}\n\\centering\n\\begin{tabular}{llrrrr}\n\\toprule\n\\rowcolor{gray!10}\\textbf{Model}&\\textbf{Value}&\\textbf{Number}&\\textbf{LPs}&\\textbf{Witnesses}&\\textbf{Time}\\\\");
            res.AddLast("Model Method Metabolites Reactions Reactions.unblocked Reactions.frev Couples Time Time.unblocked Time.frev Time.couples Time.efca LPs LPs.unblocked LPs.frev LPs.couples LPs.efca Witnesses Witnesses.unblocked Witnesses.frev Witnesses.couples Witnesses.efca");
            details.AddLast("Model Method Metabolites Reactions Reactions.unblocked Reactions.frev Couples Reaction Time LPs LPs.total Witnesses Witnesses.total Witnesses.distribution Witnesses.usable");
            String model, method;
            foreach (string path in paths)
            {

                // read network
                forbidden.Clear();

                String[] metabolite_names = null, reaction_names = null;
                if (path.Contains("matrev"))
                {
                    List<double>[] network = Reader.ReadDoubleFileToList(path, separator);
                    NetworkReader.ExtractNetwork(network, out matrix, out rev);
                    name = path.Substring(0, path.IndexOf(".matrev"));
                }
                else if (path.Contains("meta"))
                {
                    NetworkReader.ReadMetatool(Reader.ReadFile(path), out matrix, out rev, out metabolite_names, out reaction_names);
                    name = path.Substring(0, path.IndexOf(".meta"));
                }
                else
                {
                    name = path;
                    NetworkReader.ReadDir(name, separator, out matrix, out rev, out metabolite_names, out reaction_names, out forbidden);
                }
                Writer.WriteFile(NetworkReader.NetworkToDot(matrix, rev, metabolite_names, reaction_names), name + "_network.dot");
                Writer.WriteFile(NetworkReader.NetworkToConnections(matrix, rev, metabolite_names, reaction_names), name + "_connections.dot");
                Writer.WriteFile(NetworkReader.NetworkToMatlab(matrix, rev), name + ".m");

                model = name.Substring(name.LastIndexOf('\\') + 1);

                bool[] trans = FindTrans(matrix, rev);

                // Calculations.

                FluxPattern max = null;
                IIntCoupling coupling = null;
                EFCACount efca = null;
                ICollection<FluxPattern> circuits = null;
                LinkedList<String> lines = new LinkedList<string>();
                for (int milp = 0; milp < 1; ++milp)
                {
                    method = milp < 1 ? "LP" : "MILP";
                    for (int t = 0; t < n_simulations; ++t)
                    {


                        try
                        {
                            coupling = DoFCA(matrix, rev, trans, calculate_circuits ? null : forbidden, max_value, tolerance, method == "LP", fc, do_efca, keep_witnesses, only_internal, save_stats, !no_output, name, model, method, res, details, out max, out efca, out circuits);
                            Writer.WriteFile(details, String.Format("{0}_details_{1}_{2}-{3}-{4}-{5}-{6}.log", path2models.Substring(0, path2models.Length - 4), fc ? "fc" : "all", DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, DateTime.Now.Hour, DateTime.Now.Minute));
                            Writer.WriteFile(res, String.Format("{0}_{1}_{2}-{3}-{4}-{5}-{6}.log", path2models.Substring(0, path2models.Length - 4), fc ? "fc" : "all", DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, DateTime.Now.Hour, DateTime.Now.Minute));
                        }
                        catch (CalculationException e)
                        {
                            Console.WriteLine("There's happened an error during the calculation. Repeating the simulation.");
                            --t;
                        }
                    }

                    // save results of last simulation

                    String suffix = (only_internal ? "_internal" : "") + (milp == 1 ? "_milp" : "_lp") + (fc ? (calculate_circuits  ? "_fc" : "_comb" ) : "_all");

                        lines.Clear();
                        String max_string = max.ToString();
                        lines.AddLast(max_string.Substring(1, max_string.Length - 2));
                    Writer.WriteFile(lines, name + "\\unblocked"+ suffix + ".txt");
                        lines.Clear();
                        for (int r = 0; r < rev.Length; ++r)
                            if (!max[r])
                                lines.AddLast(reaction_names[r]);
                        Writer.WriteFile(lines, name + "\\blocked" + suffix + ".txt");

                        if (!only_internal)
                        {
                            Writer.WriteFile(coupling.ToCoupling().ToStrings(), name + suffix + ".coupling");
                            Writer.WriteFile(coupling.ToCoupling().ToStrings(reaction_names), name + "\\coupled" + suffix + ".txt");
                            Writer.WriteFile(coupling.ToDot(reaction_names), name + "\\coupled" + suffix + ".dot");
                            if (calculate_circuits)
                            {
                                lines.Clear();
                                foreach (FluxPattern c in circuits)
                                    lines.AddLast(c.ToString());
                                Writer.WriteFile(lines, name+"\\circuits" + suffix + ".txt");
                            }
                        if (do_efca)
                            Writer.WriteFile(efca.ToCSV(separator, reaction_names), name + "\\efca" + suffix + ".txt");

                        Writer.WriteFile(res, String.Format("{0}_{1}_{2}-{3}-{4}-{5}-{6}.log", path2models.Substring(0, path2models.Length - 4), fc ? "fc" : "all", DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, DateTime.Now.Hour, DateTime.Now.Minute));
                    }
                }
            }

            // Finish stats.
            res.AddLast("\\bottomrule\n\\end{tabular}\n\\end{table}");
            Console.WriteLine("Done. [press Enter]");
            Console.ReadLine();
        }
    }

    public class FCAOptimization
    {


        public static void Optimize(double[,] S, bool[] rev, bool increase, LinkedList<FluxPattern> witnesses, FluxPattern max, FluxPattern frev, IIntCoupling coupling, List<int> R0, List<int> Rb, bool lp, double max_value, double tolerance, bool output, int population_size, int generations, double recomb_rate, double mut_rate_act, double mut_rate_deact, double start_density, out LinkedList<String> dot_tree)
        {

            IIntFCACalculator<FluxPattern> fcacalculator = lp ? (IIntFCACalculator<FluxPattern>)new LPIntFCACalculator(S, rev, max_value, tolerance) : (IIntFCACalculator<FluxPattern>)new MILPIntFCACalculator(S, rev, max_value, tolerance);
            FitnessEstimator fe = new FitnessEstimator(R0, Rb, increase, rev.Length, fcacalculator, witnesses, max, frev, coupling);


            List<int> approximation = null;
            FluxPattern max_approx = null;

            if (population_size > 0)
            {
                Population population = new Population(fe, population_size, recomb_rate, mut_rate_act, mut_rate_deact, start_density);

                for (int i = 1; i < generations; ++i)
                    population.evolve();

                Individual best = population.FittestIndividual;
                if (best.Fitness < 0)
                {
                    approximation = best.Genes;
                    max_approx = best.Max;
                }
            }


            OptimizationTree tree = new OptimizationTree(fe, R0, Rb, approximation, max_approx, output);

            dot_tree = tree.DotTree;
            LinkedList<List<int>> opt_solutions = tree.OptGenes;

            Console.WriteLine("\n\n\nThere are {0} optimal solutions.\n\t{1} reactions are unblocked.\n\tYou have to use {2} drugs.", opt_solutions.Count, tree.OptMaxSize, tree.OptRSize);

            if (opt_solutions.Count > 0)
            {
                Individual opt = new Individual(opt_solutions.First.Value, fe, new EvolutionParams());
                Console.WriteLine("\n\tOne optimal solution is: {0}", opt);
            }

            Console.WriteLine("\n\tI had to calculate {0} out of {1} maxima to find the optima.", tree.NumberNodes, ((long)1) << Rb.Count());

        }


        public static void main(string[] args)
        {


            //init
            string path = args[0], name, separator = args[1];

            bool no_output = args.Length > 4 && int.Parse(args[4]) == 0, evolution = args.Length > 5;

            double max_value = double.Parse(args[2]), tolerance = double.Parse(args[3]);

            int population_size = evolution ? int.Parse(args[5]) : 0, generations = evolution ? int.Parse(args[6]) : 0;
            double recomb_rate = args.Length > 7 ? double.Parse(args[7]) : 0.2, mut_rate_act = args.Length > 8 ? double.Parse(args[8]) : 0.1, mut_rate_deact = args.Length > 9 ? double.Parse(args[9]) : mut_rate_act, start_density = args.Length > 10 ? double.Parse(args[10]) : 0.1;

            bool fc = false, lp = true, save_stats = false;
            bool increase = true;
            bool only_internal = false;



            double[,] matrix; bool[] rev;


            int n_unblocked;
            LinkedList<FluxPattern> witnesses = new LinkedList<FluxPattern>(), forbidden_combinations = new LinkedList<FluxPattern>();
            ICollection<FluxPattern> circuits;
            FluxPattern max = null, frev = null;
            IIntCoupling coupling = null;
            EFCACount efca = null;
            DateTime time_start, time_stop;

            List<int>[] reaction_sets;
            List<int> R0, Rb;


            // read network
            String[] metabolite_names = null, reaction_names = null;
            if (path.Contains("matrev"))
            {
                List<double>[] network = Reader.ReadDoubleFileToList(path, separator);
                NetworkReader.ExtractNetwork(network, out matrix, out rev);
                name = path.Substring(0, path.IndexOf(".matrev"));
            }
            else
            {
                NetworkReader.ReadMetatool(Reader.ReadFile(path), out matrix, out rev, out metabolite_names, out reaction_names);
                name = path.Substring(0, path.IndexOf(".meta"));
            }
            Writer.WriteFile(NetworkReader.NetworkToMatlab(matrix, rev), name + ".m");

            bool[] trans = FCA.FindTrans(matrix, rev);

            reaction_sets = Reader.ReadIntFile(name + ".rs", separator);
            R0 = reaction_sets[1]; Rb = reaction_sets[0];

            String model = name.Substring(name.LastIndexOf('\\') + 1), method = lp ? "LP" : "MILP";
            LinkedList<string> res = new LinkedList<string>(), details = new LinkedList<string>();
            res.AddLast("Model Method Metabolites Reactions Reactions.unblocked Reactions.frev Couples Time Time.unblocked Time.frev Time.couples Time.efca LPs LPs.unblocked LPs.frev LPs.couples LPs.efca Witnesses Witnesses.unblocked Witnesses.frev Witnesses.couples Witnesses.efca");
            details.AddLast("Model Method Time.start  Metabolites Reactions Reactions.unblocked Reactions.frev Couples Reaction LPs LPs.total Witnesses Witnesses.total");

            for (byte t = 0; t < 1; ++t)
            {
                try
                {
                    coupling = FCA.DoFCA(matrix, rev, trans, forbidden_combinations, max_value, tolerance, lp, fc, false, false, only_internal, !no_output, save_stats, name, model, method, res, details, out max, out efca, out circuits);
                }
                catch (CalculationException e)
                {
                    Console.WriteLine("There's happened an error during the calculation. Repeating the simulation.");
                    --t;
                }
            }

            // save results of last simulation
            Writer.WriteFile(coupling.ToCoupling().ToStrings(), name + "_lp" + ".coupling");

            LinkedList<String> dot_tree;
            Optimize(matrix, rev, increase, witnesses, max, frev, coupling, R0, Rb, true, max_value, tolerance, !no_output, population_size, generations, recomb_rate, mut_rate_act, mut_rate_deact, start_density, out dot_tree);
            Writer.WriteFile(dot_tree, name + "_tree.gv");

            Console.ReadLine();
        }
    }
}