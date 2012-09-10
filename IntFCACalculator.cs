using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Gurobi;
using Metabolism.Patterns;


namespace Metabolism.Analysis.Calculators
{


    public interface IIntFCACalculator<P> where P : IFluxPattern<P>
    {
        P CalculateMax(List<int> R, int n, ICollection<P> witnesses, out LinkedList<P> new_witnesses, out FluxPattern calcs, P MAX, P FREV, IIntCoupling C = null);

        int SolverCalls { get; }
    }

    [Serializable()]
    public class CalculationException : System.Exception
    {
        public CalculationException() : base() { }
        public CalculationException(string message) : base(message) { }
        public CalculationException(string message, System.Exception inner) : base(message, inner) { }

        // A constructor is needed for serialization when an
        // exception propagates from a remoting server to the client. 
        protected CalculationException(System.Runtime.Serialization.SerializationInfo info,
            System.Runtime.Serialization.StreamingContext context) { }
    }

    public class FCIntFCACalculator : IIntFCACalculator<FluxPattern>
    {

        public readonly double TOLERANCE, MAXVALUE;

        public readonly int m, n;
        private readonly bool[] rev, trans;

        protected GRBModel model;
        protected GRBVar[] supp, supp_pos, supp_neg;

        private int solver_calls;

        private LPIntFCACalculator relaxation_solver;
        private LinkedList<FluxPattern> circles;

        public LinkedList<FluxPattern> Circles{ get { return new LinkedList<FluxPattern>(circles);}}

        public bool ContainsCircle(FluxPattern a)
        {
            bool res;
            List<int> reactions;
            foreach (FluxPattern c in circles)
            {
                reactions = c.Values;
                res = true;
                foreach (int r in reactions)
                    res &= a[r];
                if (res)
                    return true;
            }

            return false;
        }

        
        GRBVar[] x, x_pos, x_neg;

        public FCIntFCACalculator(double[,] matrix, bool[] rev, bool[] trans, double max_value, double tolerance)
        {

            this.relaxation_solver = new LPIntFCACalculator(matrix, rev, max_value, tolerance);
            this.solver_calls = 0;
            this.circles = new LinkedList<FluxPattern>();

            this.rev = rev;
            this.trans = trans;

            this.n = rev.Length;
            this.m = matrix.Length / n;

            this.TOLERANCE = tolerance;
            this.MAXVALUE = max_value;




            // Create variables
            GRBEnv env = new GRBEnv();
            model = new GRBModel(env);
            model.GetEnv().Set(GRB.IntParam.OutputFlag, 0);
            model.GetEnv().Set(GRB.DoubleParam.IntFeasTol, 1e-9);
            

            double[] lb, ub, obj;
            char[] type;

            lb = new double[n];
            ub = new double[n];
            obj = new double[n];
            type = new char[n];

            for (int i = 0; i < n; ++i)
            {
                lb[i] = -GRB.INFINITY;
                ub[i] = GRB.INFINITY;
                obj[i] = 0;
                type[i] = GRB.CONTINUOUS;
            }
        
            
            x = model.AddVars(lb, ub, obj, type, null);
            x_pos = model.AddVars(n, GRB.CONTINUOUS);
            x_neg = model.AddVars(n, GRB.CONTINUOUS);
            
            //GRBVar[] x = model.AddVars(lb, ub, obj, type, null);
            //GRBVar[] x_pos = model.AddVars(n, GRB.CONTINUOUS);
            //GRBVar[] x_neg = model.AddVars(n, GRB.CONTINUOUS);
        


            for (int i = 0; i < n; ++i)
            {
                lb[i] = 0;
                ub[i] = 1;
                obj[i] = 0;
                type[i] = GRB.BINARY;
            }
            supp = model.AddVars(lb, ub, obj, type, null);
            supp_pos = model.AddVars(n, GRB.BINARY);
            for (int i = 0; i < n; ++i)
            {
                ub[i] = rev[i] ? 1 : 0;
                obj[i] = 0;
            }
            supp_neg = model.AddVars(lb, ub, obj, type, null);


            // Integrate new variables

            model.Update();

            GRBLinExpr expr;


            // flux conservation
            for (int i = 0; i < m; ++i)
            {
                expr = new GRBLinExpr();


                for (int j = 0; j < n; ++j)
                    if (Math.Abs(matrix[i, j]) > tolerance)
                        expr.AddTerm(matrix[i, j], x[j]);

                model.AddConstr(expr, GRB.EQUAL, 0, null);
            }

            // supports
            for (int i = 0; i < n; ++i)
            {
                model.AddConstr(supp_pos[i], GRB.LESS_EQUAL, x_pos[i], null);
                expr = new GRBLinExpr();
                expr.AddTerm(max_value, supp_pos[i]);
                model.AddConstr(expr, GRB.GREATER_EQUAL, x_pos[i], null);

                model.AddConstr(supp_neg[i], GRB.LESS_EQUAL, x_neg[i], null);
                expr = new GRBLinExpr();
                expr.AddTerm(max_value, supp_neg[i]);
                model.AddConstr(expr, GRB.GREATER_EQUAL, x_neg[i], null);

                expr = new GRBLinExpr();
                expr.AddTerm(1, supp_pos[i]);
                expr.AddTerm(1, supp_neg[i]);
                model.AddConstr(expr, GRB.EQUAL, supp[i], null);

                expr = new GRBLinExpr();
                expr.AddTerm(1, x_pos[i]);
                expr.AddTerm(-1, x_neg[i]);
                model.AddConstr(expr, GRB.EQUAL, x[i], null);
            }

            // non-trivial
            expr = new GRBLinExpr();
            for (int i = 0; i < n; ++i)
            {
                expr.AddTerm(1, supp[i]);
            }
            model.AddConstr(expr, GRB.GREATER_EQUAL, 1, null);

            model.Update();

            ForbidCycles();
        }

        private void ForbidCycles()
        {
            for (int i = 0; i < n; ++i)
            {
                supp[i].Set(GRB.DoubleAttr.Obj, 1);
                if (trans[i])
                    supp[i].Set(GRB.DoubleAttr.UB, 0);
            }

            model.Update();

            int cycles = 0;
            FluxPattern c;
            while (ForbidNextCycle(out c)) { 
                circles.AddLast(c); 
                Console.WriteLine("{0}th cycle forbidden:\n\t{1}", ++cycles, c.ToString());
            //    model.Write(String.Format("S:\\model{0:0000}.lp", cycles));
            //    Console.ReadLine();
            }
            Console.WriteLine("All cycles calculated.");
//            Console.ReadLine();

            for (int i = 0; i < n; ++i)
            {
                supp[i].Set(GRB.DoubleAttr.Obj, 0);
                if (trans[i])
                    supp[i].Set(GRB.DoubleAttr.UB, 1);
            }
            model.Update();

        }

        private bool ForbidNextCycle(out FluxPattern cycle)
        {
            bool feasible = true;

            model.Optimize();
            feasible = model.Get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL;

            cycle = null;

            if (feasible)
            {
                bool[] efp = new bool[n];

                for (int i = 0; i < n; ++i)
                    efp[i] = Math.Abs(x[i].Get(GRB.DoubleAttr.X)) > this.TOLERANCE;

                cycle = new FluxPattern(efp, false);

                /*
                if (cycle.Count < 2)
                {
                for (int i = 0; i < n; ++i)
                {
                    if (efp[i])
                    {
                        Console.WriteLine("{0}: {1} vs. {2} --> {3} + {4} = {5}", 
                            i + 1, 
                            x_pos[i].Get(GRB.DoubleAttr.X), x_neg[i].Get(GRB.DoubleAttr.X), 
                            supp_pos[i].Get(GRB.DoubleAttr.X), supp_neg[i].Get(GRB.DoubleAttr.X), 
                            supp[i].Get(GRB.DoubleAttr.X));
                    }
                }*/
                Console.WriteLine(model.Get(GRB.DoubleAttr.ObjVal));
                //}

                int k = 0;
                GRBLinExpr expr = new GRBLinExpr();
                for (int i = 0; i < n; ++i)
                {
                    if (efp[i])
                    {
                        expr.AddTerm(1, supp[i]);
                        ++k;
                    }
                }
                model.AddConstr(expr, GRB.LESS_EQUAL, k - 1, null);
                model.Update();
            }

            return feasible;
        }


        public FluxPattern CalculateMax(FluxPattern lower_bound, FluxPattern upper_bound, IIntCoupling coupling, out LinkedList<FluxPattern> witnesses, out FluxPattern calcs)
        {

            LinkedList<FluxPattern> wits_relaxed;
            upper_bound = relaxation_solver.CalculateMax(lower_bound, upper_bound, coupling, out wits_relaxed, out calcs);

            witnesses = new LinkedList<FluxPattern>();
            foreach (FluxPattern w in wits_relaxed)
                if (!this.ContainsCircle(w))
                {
                    lower_bound += w;
                    witnesses.AddLast(w);
                }

            int n = lower_bound.Length;
            bool allrev = true;
            for (int i = 0; allrev && i < n; ++i)
                allrev = rev[i] || !upper_bound[i];


            FluxPattern a;

            bool[] ub = upper_bound.Incidence, calced = new bool[n];
            for (int i = 0; i < n; ++i)
            {
                if (!ub[i])
                    supp[i].Set(GRB.DoubleAttr.UB, 0);
            }
            model.Update();

            ICollection<int> c;
            bool[] sol = new bool[n];
            for (int i = 0; i < n; ++i)
                if (calced[i] = ub[i] && !lower_bound[i])
                {

                    supp[i].Set(GRB.DoubleAttr.LB, 1);
                    model.Update();
                    model.Optimize();
                    this.solver_calls++;

                    if (model.Get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE)
                    {
                        if (coupling == null)
                        {
                            ub[i] = false;
                            supp[i].Set(GRB.DoubleAttr.UB, 0);
                        }
                        else
                        {
                            c = coupling[i];
                            foreach (int r in c)
                            {
                                ub[r] = false;
                                supp[r].Set(GRB.DoubleAttr.UB, 0);
                            }
                        }
                    }
                    else
                    {
                        for (int j = 0; j < n; ++j)
                            sol[j] = supp[j].Get(GRB.DoubleAttr.X) > 0.5;

                        a = new FluxPattern(sol, false);
                        lower_bound += a;
                        witnesses.AddLast(a);
                    }

                    supp[i].Set(GRB.DoubleAttr.LB, 0);
                    model.Update();

                }

            for (int i = 0; i < n; ++i)
            {
                supp[i].Set(GRB.DoubleAttr.UB, 1);
            }
            model.Update();

            //Console.WriteLine(lower_bound);
            //Console.ReadLine();
            calcs = new FluxPattern(calced, false);
            return lower_bound;
        }


        public FluxPattern CalculateMax(List<int> R, int n, ICollection<FluxPattern> witnesses, out LinkedList<FluxPattern> new_witnesses, out FluxPattern calcs, FluxPattern MAX = null, FluxPattern FREV = null, IIntCoupling C = null)
        {

            FluxPattern lb = new FluxPattern(new bool[n], true), ub;
            bool keep;
            foreach (FluxPattern w in witnesses)
            {
                keep = true;
                for (int i = 0; keep && i < R.Count; ++i)
                    keep = !w[R[i]];
                if (keep)
                    lb += w;
            }

            bool[] l, u;

            if (MAX != null && FREV != null)
            {
                HashSet<int> frev = new HashSet<int>(FREV.Values), nfrev = new HashSet<int>(MAX.Values);
                foreach (int r in frev)
                    nfrev.Remove(r);
                if (R.Count == 1 && frev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in nfrev)
                        l[i] = true;
                    lb += new FluxPattern(l, false);
                }
                if (nfrev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in frev)
                        l[i] = true;
                    lb += new FluxPattern(l, true);
                }

            }

            if (MAX == null)
            {
                u = new bool[n];
                for (int i = 0; i < n; ++i)
                    u[i] = true;
                ub = new FluxPattern(u, false);
            }
            else
                u = MAX.Incidence;

            ICollection<int> c;
            if (C == null)
                foreach (int r in R)
                    u[r] = false;
            else
                foreach (int r in R)
                {
                    c = C[r];
                    foreach (int i in c)
                        u[i] = false;
                }
            ub = new FluxPattern(u, false);


            return this.CalculateMax(lb, ub, C, out new_witnesses, out calcs);
        }

        public int SolverCalls
        {
            get { return this.solver_calls; }
        }


    }

    public class FastFCIntFCACalculator : IIntFCACalculator<FluxPattern>
    {

        public readonly double TOLERANCE, MAXVALUE;

        public readonly int m, n;
        private readonly bool[] rev, trans;

        protected GRBModel model;
        protected GRBVar[] supp;

        private int solver_calls;

        private LPIntFCACalculator relaxation_solver;
        private LinkedList<FluxPattern> circles;

        public bool ContainsCircle(FluxPattern a)
        {
            bool res;
            List<int> reactions;
            foreach (FluxPattern c in circles)
            {
                reactions = c.Values;
                res = true;
                foreach (int r in reactions)
                    res &= a[r];
                if (res)
                    return true;
            }

            return false;
        }


        GRBVar[] x;

        public FastFCIntFCACalculator(double[,] matrix, bool[] rev, ICollection<FluxPattern> circles, double max_value, double tolerance)
        {

            this.relaxation_solver = new LPIntFCACalculator(matrix, rev, max_value, tolerance);
            this.solver_calls = 0;
            this.circles = new LinkedList<FluxPattern>(circles);

            this.rev = rev;

            this.n = rev.Length;
            this.m = matrix.Length / n;

            this.TOLERANCE = tolerance;
            this.MAXVALUE = max_value;




            // Create variables
            GRBEnv env = new GRBEnv();
            model = new GRBModel(env);
            model.GetEnv().Set(GRB.IntParam.OutputFlag, 0);
            model.GetEnv().Set(GRB.DoubleParam.IntFeasTol, 1e-9);


            double[] lb, ub, obj;
            char[] type;

            lb = new double[n];
            ub = new double[n];
            obj = new double[n];
            type = new char[n];

            for (int i = 0; i < n; ++i)
            {
                lb[i] = rev[i] ? -GRB.INFINITY : 0;
                ub[i] = GRB.INFINITY;
                obj[i] = 0;
                type[i] = GRB.CONTINUOUS;
            }


            x = model.AddVars(lb, ub, obj, type, null);

        
            for (int i = 0; i < n; ++i)
            {
                lb[i] = 0;
                ub[i] = 1;
                obj[i] = 0;
                type[i] = GRB.BINARY;
            }
            supp = model.AddVars(lb, ub, obj, type, null);


            // Integrate new variables

            model.Update();

            GRBLinExpr expr;


            // flux conservation
            for (int i = 0; i < m; ++i)
            {
                expr = new GRBLinExpr();


                for (int j = 0; j < n; ++j)
                    if (Math.Abs(matrix[i, j]) > tolerance)
                        expr.AddTerm(matrix[i, j], x[j]);

                model.AddConstr(expr, GRB.EQUAL, 0, null);
            }

            // supports fuer a_i = 0 --> v_i = 0
            for (int i = 0; i < n; ++i)
            {
                expr = new GRBLinExpr();
                expr.AddTerm(max_value, supp[i]);
                model.AddConstr(expr, GRB.GREATER_EQUAL, x[i], null);

                if (rev[i])
                {
                    expr = new GRBLinExpr();
                    expr.AddTerm(-max_value, supp[i]);
                    model.AddConstr(expr, GRB.LESS_EQUAL, x[i], null);
                }
            }


            model.Update();

            ForbidCycles();
        }

        private void ForbidCycles()
        {
            int k;
            GRBLinExpr expr;
            foreach (FluxPattern cycle in this.circles)
            {

                k = 0;
                expr = new GRBLinExpr();
                for (int i = 0; i < n; ++i)
                {
                    if (cycle[i])
                    {
                        expr.AddTerm(1, supp[i]);
                        ++k;
                    }
                }
                model.AddConstr(expr, GRB.LESS_EQUAL, k - 1, null);
            }
            model.Update();
        }

        private bool ForbidNextCycle(out FluxPattern cycle)
        {
            bool feasible = true;

            model.Optimize();
            feasible = model.Get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL;

            cycle = null;

            if (feasible)
            {
                bool[] efp = new bool[n];

                for (int i = 0; i < n; ++i)
                    efp[i] = Math.Abs(x[i].Get(GRB.DoubleAttr.X)) > this.TOLERANCE;

                cycle = new FluxPattern(efp, false);

                /*
                if (cycle.Count < 2)
                {
                for (int i = 0; i < n; ++i)
                {
                    if (efp[i])
                    {
                        Console.WriteLine("{0}: {1} vs. {2} --> {3} + {4} = {5}", 
                            i + 1, 
                            x_pos[i].Get(GRB.DoubleAttr.X), x_neg[i].Get(GRB.DoubleAttr.X), 
                            supp_pos[i].Get(GRB.DoubleAttr.X), supp_neg[i].Get(GRB.DoubleAttr.X), 
                            supp[i].Get(GRB.DoubleAttr.X));
                    }
                }*/
                Console.WriteLine(model.Get(GRB.DoubleAttr.ObjVal));
                //}

                int k = 0;
                GRBLinExpr expr = new GRBLinExpr();
                for (int i = 0; i < n; ++i)
                {
                    if (efp[i])
                    {
                        expr.AddTerm(1, supp[i]);
                        ++k;
                    }
                }
                model.AddConstr(expr, GRB.LESS_EQUAL, k - 1, null);
                model.Update();
            }

            return feasible;
        }


        public FluxPattern CalculateMax(FluxPattern lower_bound, FluxPattern upper_bound, IIntCoupling coupling, out LinkedList<FluxPattern> witnesses, out FluxPattern calcs)
        {

            LinkedList<FluxPattern> wits_relaxed;
            upper_bound = relaxation_solver.CalculateMax(lower_bound, upper_bound, coupling, out wits_relaxed, out calcs);


            witnesses = new LinkedList<FluxPattern>();
            foreach (FluxPattern w in wits_relaxed)
                if (!this.ContainsCircle(w))
                {
                    lower_bound += w;
                    witnesses.AddLast(w);
                }

            int n = lower_bound.Length;
            bool allrev = true;
            for (int i = 0; allrev && i < n; ++i)
                allrev = rev[i] || !upper_bound[i];


            FluxPattern a;

            bool[] ub = upper_bound.Incidence, calced = new bool[n];
            for (int i = 0; i < n; ++i)
            {
                if (!ub[i])
                    supp[i].Set(GRB.DoubleAttr.UB, 0);
            }
            model.Update();

            ICollection<int> c;
            bool[] sol = new bool[n];
            for (int i = 0; i < n; ++i)
                if (calced[i] = ub[i] && !lower_bound[i])
                {

                    x[i].Set(GRB.DoubleAttr.LB, 1);
                    model.Update();
                    model.Optimize();
                    this.solver_calls++;

                    if (rev[i] && !allrev && model.Get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE)
                    {
                        x[i].Set(GRB.DoubleAttr.LB, -GRB.INFINITY);
                        x[i].Set(GRB.DoubleAttr.UB, -1);
                        model.Update();
                        model.Optimize();
                        this.solver_calls++;
                    }

                    if (model.Get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE)
                    {
                        if (coupling == null)
                        {
                            ub[i] = false;
                            supp[i].Set(GRB.DoubleAttr.UB, 0);
                        }
                        else
                        {
                            c = coupling[i];
                            foreach (int r in c)
                            {
                                ub[r] = false;
                                supp[r].Set(GRB.DoubleAttr.UB, 0);
                            }
                        }
                    }
                    else
                    {
                        for (int j = 0; j < n; ++j)
                            sol[j] = Math.Abs(x[j].Get(GRB.DoubleAttr.X)) > TOLERANCE;

                        a = new FluxPattern(sol, false);
                        lower_bound += a;
                        witnesses.AddLast(a);
                    }

                    x[i].Set(GRB.DoubleAttr.UB, GRB.INFINITY);
                    x[i].Set(GRB.DoubleAttr.LB, rev[i] ? -GRB.INFINITY : 0);
                    model.Update();

                }

            for (int i = 0; i < n; ++i)
            {
                supp[i].Set(GRB.DoubleAttr.UB, 1);
            }
            model.Update();

            //Console.WriteLine(lower_bound);
            //Console.ReadLine();
            calcs = new FluxPattern(calced, false);
            return lower_bound;
        }


        public FluxPattern CalculateMax(List<int> R, int n, ICollection<FluxPattern> witnesses, out LinkedList<FluxPattern> new_witnesses, out FluxPattern calcs, FluxPattern MAX = null, FluxPattern FREV = null, IIntCoupling C = null)
        {

            FluxPattern lb = new FluxPattern(new bool[n], true), ub;
            bool keep;
            foreach (FluxPattern w in witnesses)
            {
                keep = true;
                for (int i = 0; keep && i < R.Count; ++i)
                    keep = !w[R[i]];
                if (keep)
                    lb += w;
            }

            bool[] l, u;

            if (MAX != null && FREV != null)
            {
                HashSet<int> frev = new HashSet<int>(FREV.Values), nfrev = new HashSet<int>(MAX.Values);
                foreach (int r in frev)
                    nfrev.Remove(r);
                if (R.Count == 1 && frev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in nfrev)
                        l[i] = true;
                    lb += new FluxPattern(l, false);
                }
                if (nfrev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in frev)
                        l[i] = true;
                    lb += new FluxPattern(l, true);
                }

            }

            if (MAX == null)
            {
                u = new bool[n];
                for (int i = 0; i < n; ++i)
                    u[i] = true;
                ub = new FluxPattern(u, false);
            }
            else
                u = MAX.Incidence;

            ICollection<int> c;
            if (C == null)
                foreach (int r in R)
                    u[r] = false;
            else
                foreach (int r in R)
                {
                    c = C[r];
                    foreach (int i in c)
                        u[i] = false;
                }
            ub = new FluxPattern(u, false);


            return this.CalculateMax(lb, ub, C, out new_witnesses, out calcs);
        }

        public int SolverCalls
        {
            get { return this.solver_calls; }
        }


    }


    public class MILPIntFCACalculator : IIntFCACalculator<FluxPattern>
    {

        public readonly double TOLERANCE, MAXVALUE;

        public readonly int m, n;
        private readonly bool[] rev;

        protected GRBModel model;
        protected GRBVar[] supp, supp_pos, supp_neg, supp_irr, x, x_pos, x_neg,
          x_irr;

        public MILPIntFCACalculator(double[,] matrix, bool[] rev, double max_value, double tolerance)
        {

            this.rev = rev;

            this.n = rev.Length;
            this.m = matrix.Length / n;

            this.TOLERANCE = tolerance;
            this.MAXVALUE = max_value;

            int n_rev, n_irr;
            List<int> rev_indices, irr_indices;
            int[] new_indices;

            rev_indices = new List<int>(n);
            irr_indices = new List<int>(n);
            new_indices = new int[n];
            for (int i = 0; i < n; ++i)
                if (rev[i])
                {
                    new_indices[i] = rev_indices.Count;
                    rev_indices.Add(i);
                }
                else
                {
                    new_indices[i] = irr_indices.Count;
                    irr_indices.Add(i);
                }


            n_rev = rev_indices.Count;
            n_irr = n - n_rev;


            // Create variables
            GRBEnv env = new GRBEnv();
            model = new GRBModel(env);
            model.GetEnv().Set(GRB.IntParam.OutputFlag, 0);

            x = model.AddVars(n, GRB.CONTINUOUS);
            x_pos = model.AddVars(n_rev, GRB.CONTINUOUS);
            x_neg = model.AddVars(n_rev, GRB.CONTINUOUS);
            x_irr = model.AddVars(n_irr, GRB.CONTINUOUS);

            double[] lb, ub, obj;
            char[] type;



            lb = new double[n];
            ub = new double[n];
            obj = new double[n];
            type = new char[n];
            for (int i = 0; i < n; ++i)
            {
                ub[i] = 1;
                obj[i] = -1;
                type[i] = GRB.BINARY;
            }
            supp = model.AddVars(lb, ub, obj, type, null);
            supp_pos = model.AddVars(n_rev, GRB.BINARY);
            supp_neg = model.AddVars(n_rev, GRB.BINARY);
            supp_irr = model.AddVars(n_irr, GRB.BINARY);


            // Integrate new variables

            model.Update();

            GRBLinExpr expr;


            // flux conservation
            for (int i = 0; i < m; ++i)
            {
                expr = new GRBLinExpr();


                for (int j = 0; j < n; ++j)
                    if (Math.Abs(matrix[i, j]) > tolerance)
                        expr.AddTerm(matrix[i, j], x[j]);

                model.AddConstr(expr, GRB.EQUAL, 0, null);
            }

            // thermodynamical constraints
            for (int j = 0; j < n; ++j)
            {
                x[j].Set(GRB.DoubleAttr.UB, MAXVALUE);
                x[j].Set(GRB.DoubleAttr.LB, rev[j] ? -MAXVALUE : 0);
            }

            // supports
            for (int i = 0; i < n_rev; ++i)
            {
                model.AddConstr(supp_pos[i], GRB.LESS_EQUAL, x_pos[i], null);
                expr = new GRBLinExpr();
                expr.AddTerm(max_value, supp_pos[i]);
                model.AddConstr(expr, GRB.GREATER_EQUAL, x_pos[i], null);

                model.AddConstr(supp_neg[i], GRB.LESS_EQUAL, x_neg[i], null);
                expr = new GRBLinExpr();
                expr.AddTerm(max_value, supp_neg[i]);
                model.AddConstr(expr, GRB.GREATER_EQUAL, x_neg[i], null);

                expr = new GRBLinExpr();
                expr.AddTerm(1, supp_pos[i]);
                expr.AddTerm(1, supp_neg[i]);
                model.AddConstr(expr, GRB.EQUAL, supp[rev_indices[i]], null);

                expr = new GRBLinExpr();
                expr.AddTerm(1, x_pos[i]);
                expr.AddTerm(-1, x_neg[i]);
                model.AddConstr(expr, GRB.EQUAL, x[rev_indices[i]], null);
            }

            for (int i = 0; i < n_irr; ++i)
            {
                model.AddConstr(supp_irr[i], GRB.LESS_EQUAL, x_irr[i], null);
                expr = new GRBLinExpr();
                expr.AddTerm(max_value, supp_irr[i]);
                model.AddConstr(expr, GRB.GREATER_EQUAL, x_irr[i], null);

                model.AddConstr(supp_irr[i], GRB.EQUAL, supp[irr_indices[i]], null);
                model.AddConstr(x_irr[i], GRB.EQUAL, x[irr_indices[i]], null);
            }


            model.Update();

        }

        private FluxPattern GetSolutionPattern()
        {

            if (model.Get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE)
                throw new CalculationException();

            bool[] pattern = new bool[n];

            bool is_rev = true;
            for (int j = 0; j < n; ++j)
            {
                if (Math.Abs(x[j].Get(GRB.DoubleAttr.X)) > TOLERANCE)
                {
                    pattern[j] = true;
                    is_rev &= rev[j];
                }
                //Console.Write("{0:0.##} ", x[j].Get(GRB.DoubleAttr.X));
            }
            //Console.WriteLine();
            //Console.WriteLine("{0:0.0} = {1:0.0} - {2:0.0} --> {3} = {4}+{5}", x[22].Get(GRB.DoubleAttr.X), x_pos[22].Get(GRB.DoubleAttr.X), x_neg[22].Get(GRB.DoubleAttr.X), supp[22].Get(GRB.DoubleAttr.X), supp_pos[22].Get(GRB.DoubleAttr.X), supp_neg[22].Get(GRB.DoubleAttr.X));
            //Console.ReadLine();

            return new FluxPattern(pattern, is_rev);
        }


        public FluxPattern CalculateMax(FluxPattern lower_bound, FluxPattern upper_bound, out LinkedList<FluxPattern> witnesses)
        {

            witnesses = new LinkedList<FluxPattern>();
            FluxPattern res;

            for (int i = 0; i < n; ++i)
            {
                if (lower_bound[i])
                {
                    supp[i].Set(GRB.DoubleAttr.LB, 1);
                    if (!rev[i])
                        x[i].Set(GRB.DoubleAttr.LB, 1);
                }
                if (!upper_bound[i])
                {
                    supp[i].Set(GRB.DoubleAttr.UB, 0);
                    x[i].Set(GRB.DoubleAttr.UB, 0);
                    if (!rev[i])
                        x[i].Set(GRB.DoubleAttr.LB, 0);
                }
            }
            model.Update();

            model.Optimize();
            res = GetSolutionPattern();

            witnesses.AddLast(res);
            //Console.WriteLine(res);
            //Console.ReadLine();


            for (int i = 0; i < n; ++i)
            {
                if (lower_bound[i])
                {
                    supp[i].Set(GRB.DoubleAttr.LB, 0);
                    if (!rev[i])
                        x[i].Set(GRB.DoubleAttr.LB, 0);
                }
                if (!upper_bound[i])
                {
                    supp[i].Set(GRB.DoubleAttr.UB, 1);
                    x[i].Set(GRB.DoubleAttr.UB, MAXVALUE);
                    if (!rev[i])
                        x[i].Set(GRB.DoubleAttr.LB, -MAXVALUE);
                }
            }
            model.Update();

            return res;
        }


        public FluxPattern CalculateMax(List<int> R, int n, ICollection<FluxPattern> witnesses, out LinkedList<FluxPattern> new_witnesses, out FluxPattern calcs, FluxPattern MAX = null, FluxPattern FREV = null, IIntCoupling C = null)
        {

            FluxPattern lb = new FluxPattern(new bool[n], true), ub;
            bool keep;
            foreach (FluxPattern w in witnesses)
            {
                keep = true;
                for (int i = 0; keep && i < R.Count; ++i)
                    keep = !w[R[i]];
                if (keep)
                    lb += w;
            }

            bool[] l, u;

            if (MAX != null && FREV != null)
            {
                HashSet<int> frev = new HashSet<int>(FREV.Values), nfrev = new HashSet<int>(MAX.Values);
                foreach (int r in frev)
                    nfrev.Remove(r);
                if (frev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in nfrev)
                        l[i] = true;
                    lb += new FluxPattern(l, false);
                }
                if (nfrev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in frev)
                        l[i] = true;
                    lb += new FluxPattern(l, true);
                }

            }

            if (MAX == null)
            {
                u = new bool[n];
                for (int i = 0; i < n; ++i)
                    u[i] = true;
                ub = new FluxPattern(u, false);
            }
            else
                u = MAX.Incidence;

            ICollection<int> c;
            if (C == null)
                foreach (int r in R)
                    u[r] = false;
            else
                foreach (int r in R)
                {
                    c = C[r];
                    foreach (int i in c)
                        u[i] = false;
                }
            ub = new FluxPattern(u, false);

            calcs = null;
            return this.CalculateMax(lb, ub, out new_witnesses);
        }


        public int SolverCalls
        {
            get { throw new NotImplementedException(); }
        }
    }


    public class LPIntFCACalculator : IIntFCACalculator<FluxPattern>
    {

        protected GRBModel model;

        protected GRBVar[] x;

        public readonly double TOLERANCE, MAXVALUE;
        public readonly int m, n;

        private readonly bool[] rev;

        int solver_calls;


        public LPIntFCACalculator(double[,] matrix, bool[] rev, double max_value, double tolerance)
        {

            this.solver_calls = 0;

            this.rev = rev;

            this.n = rev.Length;
            this.m = matrix.Length / n;

            this.TOLERANCE = tolerance;
            this.MAXVALUE = max_value;


            // Create variables

            GRBEnv env = new GRBEnv();
            model = new GRBModel(env);
            model.GetEnv().Set(GRB.IntParam.OutputFlag, 0);

            double[] lb = new double[n],
            ub = new double[n],
            obj = new double[n];
            char[] type = new char[n];

            for (int i = 0; i < n; ++i)
            {
                lb[i] = rev[i] ? -GRB.INFINITY : 0;
                ub[i] = GRB.INFINITY;
                obj[i] = 0;
                type[i] = GRB.CONTINUOUS;
            }
            x = model.AddVars(lb, ub, obj, type, null);


            // Integrate new variables

            model.Update();

            GRBLinExpr expr;


            // flux conservation
            for (int i = 0; i < m; ++i)
            {
                expr = new GRBLinExpr();


                for (int j = 0; j < n; ++j)
                    if (Math.Abs(matrix[i, j]) > tolerance)
                        expr.AddTerm(matrix[i, j], x[j]);

                model.AddConstr(expr, GRB.EQUAL, 0, null);
            }

            model.Update();
        }



        public FluxPattern CalculateMax(FluxPattern lower_bound, FluxPattern upper_bound, IIntCoupling coupling, out LinkedList<FluxPattern> witnesses, out FluxPattern calcs)
        {

            int n = lower_bound.Length;
            bool allrev = true;
            for (int i = 0; allrev && i < n; ++i)
                allrev = rev[i] || !upper_bound[i];


            witnesses = new LinkedList<FluxPattern>();
            FluxPattern a;

            bool[] ub = upper_bound.Incidence, calced = new bool[n];
            for (int i = 0; i < n; ++i)
            {
                if (!ub[i])
                {
                    x[i].Set(GRB.DoubleAttr.UB, 0);
                    if (rev[i])
                        x[i].Set(GRB.DoubleAttr.LB, 0);
                }
            }
            model.Update();

            ICollection<int> c;
            bool[] sol = new bool[n];
            for (int i = 0; i < n; ++i)
                if (calced[i] = ub[i] && !lower_bound[i])
                {

                    x[i].Set(GRB.DoubleAttr.LB, 1);
                    model.Update();
                    model.Optimize();
                    this.solver_calls++;

                    if (rev[i] && !allrev && model.Get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE)
                    {
                        x[i].Set(GRB.DoubleAttr.LB, -GRB.INFINITY);
                        x[i].Set(GRB.DoubleAttr.UB, -1);
                        model.Update();
                        model.Optimize();
                        this.solver_calls++;
                    }

                    if (model.Get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE)
                    {
                        if (coupling == null)
                        {
                            ub[i] = false;
                            x[i].Set(GRB.DoubleAttr.UB, 0);
                            if (rev[i])
                                x[i].Set(GRB.DoubleAttr.LB, 0);
                        }
                        else
                        {
                            c = coupling[i];
                            foreach (int r in c)
                            {
                                ub[r] = false;
                                x[r].Set(GRB.DoubleAttr.UB, 0);
                                if (rev[r])
                                    x[r].Set(GRB.DoubleAttr.LB, 0);
                            }
                        }
                    }
                    else
                    {
                        for (int j = 0; j < n; ++j)
                            sol[j] = Math.Abs(x[j].Get(GRB.DoubleAttr.X)) > TOLERANCE;

                        a = new FluxPattern(sol, false);
                        lower_bound += a;
                        witnesses.AddLast(a);
                    }

                    x[i].Set(GRB.DoubleAttr.UB, GRB.INFINITY);
                    x[i].Set(GRB.DoubleAttr.LB, rev[i] ? -GRB.INFINITY : 0);
                    model.Update();

                }

            for (int i = 0; i < n; ++i)
            {
                x[i].Set(GRB.DoubleAttr.UB, GRB.INFINITY);
                if (rev[i])
                    x[i].Set(GRB.DoubleAttr.LB, -GRB.INFINITY);

            }
            model.Update();

            //Console.WriteLine(lower_bound);
            //Console.ReadLine();
            calcs = new FluxPattern(calced, false);
            return lower_bound;
        }


        public FluxPattern CalculateMax(List<int> R, int n, ICollection<FluxPattern> witnesses, out LinkedList<FluxPattern> new_witnesses, out FluxPattern calcs, FluxPattern MAX = null, FluxPattern FREV = null, IIntCoupling C = null)
        {

            FluxPattern lb = new FluxPattern(new bool[n], true), ub;
            bool keep;
            foreach (FluxPattern w in witnesses)
            {
                keep = true;
                for (int i = 0; keep && i < R.Count; ++i)
                    keep = !w[R[i]];
                if (keep)
                    lb += w;
            }

            bool[] l, u;

            if (MAX != null && FREV != null)
            {
                HashSet<int> frev = new HashSet<int>(FREV.Values), nfrev = new HashSet<int>(MAX.Values);
                foreach (int r in frev)
                    nfrev.Remove(r);
                if (R.Count == 1 && frev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in nfrev)
                        l[i] = true;
                    lb += new FluxPattern(l, false);
                }
                if (nfrev.IsSupersetOf(R))
                {
                    l = new bool[n];
                    foreach (int i in frev)
                        l[i] = true;
                    lb += new FluxPattern(l, true);
                }

            }

            if (MAX == null)
            {
                u = new bool[n];
                for (int i = 0; i < n; ++i)
                    u[i] = true;
                ub = new FluxPattern(u, false);
            }
            else
                u = MAX.Incidence;

            ICollection<int> c;
            if (C == null)
                foreach (int r in R)
                    u[r] = false;
            else
                foreach (int r in R)
                {
                    c = C[r];
                    foreach (int i in c)
                        u[i] = false;
                }
            ub = new FluxPattern(u, false);


            return this.CalculateMax(lb, ub, C, out new_witnesses, out calcs);
        }

        public int SolverCalls
        {
            get { return this.solver_calls; }
        }
    }
}