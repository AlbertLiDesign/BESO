using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using System.IO;

namespace BESO
{
    public class BESO2D
    {
        #region Resolution
        public int nelx;
        public int nely;
        #endregion

        #region FE varialbes
        private double[] U;
        private int[] ik;
        private int[] jk;
        private double[] vk;
        #endregion
        #region BESO Parameters
        /// <summary>
        /// Filter radius
        /// </summary>
        public double rmin;

        /// <summary>
        /// Volume fraction
        /// </summary>
        public double vf;

        /// <summary>
        /// Penalty exponent
        /// </summary>
        public double p;

        /// <summary>
        /// Evolution rate
        /// </summary>
        public double ert;

        /// <summary>
        /// The maximum iteration
        /// </summary>
        public int maxIter;
        #endregion

        #region BESO variables
        /// <summary>
        /// Design variables
        /// </summary>
        public double[] Xe;

        public double[] dc;
        public double[] dc_old;

        /// <summary>
        /// Elemental stiffness matrix
        /// </summary>
        private double[] Ke;

        /// <summary>
        /// The minimum design variable
        /// </summary>
        private double Xmin = 0.001;

        /// <summary>
        /// The isovalue for extracting isosurface.
        /// </summary>
        public List<double> isovalues = new List<double>();

        /// <summary>
        /// The iterative history of the global compliance
        /// </summary>
        private List<double> HistoryC = new List<double>();

        public int[] free_dofs;

        public double Compliance = 0.0;

        public double vol = 1.0;
        public double delta = 1.0;
        public int iter = 0;
        public bool convergence = false;
        public string info = null;
        public StringBuilder initInfo;
        public StringBuilder optInfo;
        #endregion

        #region Flt variables
        public int[] ih;
        public int[] jh;
        public double[] vh;
        public double[] sh;
        #endregion

        private Stopwatch stopwatch;

        public bool OutputK = false;
        public bool changeSupports = true;
        public BESO2D() { }
        public BESO2D(double rmin, double vf, double ert = 0.02, double p = 3.0, int maxIter = 100)
        {
            if (rmin <= 0.0)
                throw new Exception("Rmin must be large than 0.");
            if (!(vf > 0.0 && vf < 1.0))
                throw new Exception("Vt must be large than 0 and be less than 1.");

            this.vf = vf;
            this.p = p;
            this.ert = ert;
            this.maxIter = maxIter;
            this.rmin = rmin;

            stopwatch = new Stopwatch();
            optInfo = new StringBuilder("====================== Optimization ======================" + '\n');
        }
        public BESO2D(BESO2D ibeso)
        {
            vf = ibeso.vf;
            p = ibeso.p;
            ert = ibeso.ert;
            maxIter = ibeso.maxIter;
            rmin = ibeso.rmin;

            nelx = ibeso.nelx;
            nely = ibeso.nely;

            Xe = (double[])ibeso.Xe.Clone();
            dc = (double[])ibeso.dc.Clone();
            dc_old = (double[])ibeso.dc_old.Clone();
            Ke = (double[])ibeso.Ke.Clone();

            ih = (int[])ibeso.ih.Clone();
            jh = (int[])ibeso.jh.Clone();
            vh = (double[])ibeso.vh.Clone();
            sh = (double[])ibeso.sh.Clone();

            ik = (int[])ibeso.ik.Clone();
            jk = (int[])ibeso.jk.Clone();

            stopwatch = new Stopwatch();

            optInfo = new StringBuilder("====================== Optimization ======================" + '\n');
        }


        public void Initialize(int nelx, int nely)
        {
            initInfo = new StringBuilder("====================== Launch BESO ======================" + '\n');

            this.nelx = nelx;
            this.nely = nely;

            dc = new double[nely * nelx];
            dc_old = new double[nely * nelx];
            Xe = new double[nely* nelx];
            Array.Fill(Xe, 1.0);

            stopwatch.Start();
            GetKe();
            ik = new int[nelx * nely * 8 * 8];
            jk = new int[nelx * nely * 8 * 8];
            PreFE(nelx, nely, ik, jk);
            stopwatch.Stop();
            initInfo.Append("PreFE: " + stopwatch.Elapsed.TotalMilliseconds + '\n');

            stopwatch.Restart();
            PreFlt();
            stopwatch.Stop();
            initInfo.Append("PreFlt: " + stopwatch.Elapsed.TotalMilliseconds + '\n');
        }
        public void Optimize()
        {
            if (delta > 0.001 && iter < maxIter)
            {
                optInfo.Append("====================== Iter: " + iter.ToString() + " ======================" + '\n');
                iter += 1;
                vol = Math.Max(vf, vol * (1.0 - ert));

                #region FEA
                stopwatch.Restart();
                FE();
                stopwatch.Stop();
                #endregion
                #region Prepare report
                optInfo.Append("FEA:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion

                #region Get DC
                stopwatch.Restart();
                GetDc();
                HistoryC.Add(Compliance);
                stopwatch.Stop();
                #endregion
                #region Prepare report
                optInfo.Append("Getting Sensitivity:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion

                #region Flt
                stopwatch.Restart();
                Flt(dc.Length, dc, sh);

                if (iter > 1)
                    for (int j = 0; j < nely; j++)
                    {
                        for (int i = 0; i < nelx; i++)
                        {
                            dc[i * nely + j] = (dc[i * nely + j] + dc_old[i * nely + j]) * 0.5;
                        }
                    }

                // Record the sensitiveies in each step
                dc_old = (double[])dc.Clone();
                stopwatch.Stop();
                #endregion
                #region Prepare report
                optInfo.Append("Flt:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion

                #region ADD & DEL
                stopwatch.Restart();
                ADD_DEL(vol);
                stopwatch.Stop();
                #endregion
                #region Prepare report
                optInfo.Append("ADD & DEL:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion


                #region Checking Convergence
                stopwatch.Restart();

                // Check convergence 
                if (iter > 10)
                {
                    var newV = 0.0;
                    var lastV = 0.0;
                    for (int i = 1; i < 6; i++)
                    {
                        newV += HistoryC[HistoryC.Count - i];
                        lastV += HistoryC[HistoryC.Count - 5 - i];
                    }
                    delta = Math.Abs((newV - lastV) / lastV);
                }
                #endregion
                #region Prepare report
                optInfo.Append("Checking Convergence:" + stopwatch.Elapsed.TotalMilliseconds + '\n');

                optInfo.Append("Volume: " + vol.ToString() + '\n');
                optInfo.Append("Compliance: " + Compliance.ToString() + '\n');
                optInfo.Append("Change: " + delta.ToString() + '\n');
                #endregion

                info = "Iter: " + iter.ToString() + ", Volume: " + vol.ToString()
                    + ", Compliance: " + Compliance.ToString() + ", Change: " + delta.ToString();
            }
            else
            {
                convergence = true;
            }
        }

        private void PreFlt()
        {
            int rminf = (int)Math.Floor(rmin);

            ih = new int[(int)(nelx * nely * Math.Pow((2 * rminf + 1), 2))];
            jh = new int[ih.Length];
            vh = new double[ih.Length];
            sh = new double[nelx * nely];

            int sum = 0;
            for (int i = 0; i < nelx; i++)
            {
                for (int j = 0; j < nely; j++)
                {
                    var e1 = i * nely + j + 1;
                    for (int k = Math.Max(i - rminf, 0); k < Math.Min(i + rminf + 1, nelx); k++)
                    {
                        for (int l = Math.Max(j - rminf, 0); l < Math.Min(j + rminf + 1, nely); l++)
                        {
                            var e2 = k * nely + l + 1;
                            ih[sum] = e1 - 1;
                            jh[sum] = e2 - 1;
                            vh[sum] = Math.Max(0.0, rminf - Math.Sqrt((i - k) * (i - k) + (j - l) * (j - l)));
                            sum++;
                        }
                    }
                }
            }

            GetRowSum(sum, nelx * nely, ih, jh, vh, sh);
        }
        private void ADD_DEL(double volfra)
        {
            double lowest = dc.Min();
            double highest = dc.Max();
            double th = 0.0;
            double vol = volfra * nelx * nely;
            while (((highest - lowest) / highest) > 1e-5)
            {
                th = (highest + lowest) * 0.5;
                double sum = 0.0;
                for (int j = 0; j < nely; j++)
                {
                    for (int i = 0; i < nelx; i++)
                    {
                        Xe[j*nelx+ i] = dc[i * nely + j] > th ? 1.0 : Xmin;
                        sum += Xe[j * nelx + i];
                    }
                }
                if (sum - vol > 0.0) lowest = th;
                else highest = th;
            }
            isovalues.Add(th);
        }
        private void GetDc()
        {
            Compliance = 0.0;
            for (int ely = 0; ely < nely; ely++)
            {
                for (int elx = 0; elx < nelx; elx++)
                {
                    var n1 = (nely + 1) * elx + ely + 1;
                    var n2 = (nely + 1) * (elx + 1) + ely + 1;

                    double[] Ue = { U[2 * n1 - 2], U[2 * n1 - 1], U[2 * n2 - 2], U[2 * n2 - 1],
                    U[2 * n2], U[2 * n2 + 1], U[2 * n1], U[2 * n1 + 1]};

                    double v = TransposeMultiply(8, 8, Ke, Ue);

                    Compliance += 0.5 * Math.Pow(Xe[ely*nelx + elx], p) * v;
                    dc[elx * nely + ely] = 0.5 * Math.Pow(Xe[ely* nelx + elx], p - 1) * v;
                }
            }
        }
        private void FE()
        {
            int num_allDofs = 2 * (nelx + 1) * (nely + 1);
            int num_fixedDofs = 2 * (nely + 1);
            int num_freeDofs = num_allDofs - num_fixedDofs;

            // Assemble stiffness matrix with all DOFs
            vk = new double[64 * nelx * nely];
            for (int i = 0; i < nelx; i++)
            {
                for (int j = 0; j < nely; j++)
                {
                    var ex = Math.Pow(Xe[j * nelx + i], p);
                    for (int a = 0; a < 8; a++)
                    {
                        for (int b = 0; b < 8; b++)
                        {
                            vk[i * nely * 64 + j * 64 + a * 8 + b] = ex * Ke[a * 8 + b];
                        }
                    }
                }
            }

            var F = new double[num_allDofs];
            U = new double[num_allDofs];

            // Define force vector
            F[2 * (nelx + 1) * (nely + 1) - nely - 1] = -1.0;

            if (changeSupports)
            {
                // Define fixed dofs
                var fixed_dofs = new int[num_fixedDofs];
                for (int i = 0; i < num_fixedDofs; i++)
                    fixed_dofs[i] = i;

                var all_dofs = new int[num_allDofs];
                for (int i = 0; i < num_allDofs; i++)
                    all_dofs[i] = i;

                // Obtain free dofs
                free_dofs = all_dofs.Except(fixed_dofs).ToArray();
                changeSupports = false;
            }

            var U_freedof = new double[num_freeDofs];
            Assembly_Solve(num_freeDofs, num_allDofs, ik.Length, free_dofs, ik, jk, vk, F, U_freedof);

            for (int i = 0; i < num_freeDofs; i++)
            {
                U[free_dofs[i]] = U_freedof[i];
            }
        }

        private void GetKe()
        {
            var E = 1.0;
            var nu = 0.3;
            var w = E / (1 - nu * nu);
            var k = new double[8]
            {
            w * (0.5-nu/6.0), w * (0.125+nu/8.0), w * (-0.25-nu/12.0), w * (-0.125+3*nu/8.0),
            w * (-0.25+nu/12.0), w * (-0.125-nu/8.0), w * (nu/6.0), w * (0.125-3*nu/8.0)
            };
            Ke = new double[64]{
                k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],
                k[1],k[0],k[7],k[6],k[5],k[4],k[3],k[2],
                k[2],k[7],k[0],k[5],k[6],k[3],k[4],k[1],
                k[3],k[6],k[5],k[0],k[7],k[2],k[1],k[4],
                k[4],k[5],k[6],k[7],k[0],k[1],k[2],k[3],
                k[5],k[4],k[3],k[2],k[1],k[0],k[7],k[6],
                k[6],k[3],k[4],k[1],k[2],k[7],k[0],k[5],
                k[7],k[2],k[1],k[4],k[3],k[6],k[5],k[0]};
        }

        [DllImport("Solver.dll")]
        private static extern void PreFE(int nelx, int nely, int[] ik, int[] jk);
        [DllImport("Solver.dll")]
        private static extern void Assembly_Solve(int num_freeDofs, int num_allDofs, int num_triplets, int[] free_dofs, int[] ik, int[] jk, double[] vk, double[] F, double[] U);

        [DllImport("Solver.dll")]
        private static extern double TransposeMultiply(int rows, int cols, double[] A, double[] U);

        [DllImport("Solver.dll")]
        private static extern void Flt(int dc_length, double[] dc, double[] sh);
        [DllImport("Solver.dll")]
        private static extern void GetRowSum(int coo_length, int rows, int[] ih, int[] jh, double[] vh, double[] sh);

        #region Debug Methods
        public StringBuilder ModelInfo()
        {
            StringBuilder report = new StringBuilder("=================== Model Info ===================" + '\n');
            report.Append("Disp Nodes: " + ((nelx + 1) * (nely + 1)).ToString() + '\n');
            report.Append("Disp Elements: " + (nelx * nely).ToString() + '\n');
            report.Append('\n');
            report.Append("Opt Nodes: " + ((nelx + 1) * (nely + 1)).ToString() + '\n');
            report.Append("Opt Elements: " + (nelx * nely).ToString() + '\n');
            report.Append('\n');
            report.Append("=================== Parameters Info ===================" + '\n');
            report.Append("xCount: " + nelx.ToString() + '\n');
            report.Append("yCount: " + nely.ToString() + '\n');
            return report;
        }

        public void WriteXe(string path)
        {
            string output = path + '\\' + "Xe2.txt";
            StreamWriter sw = new StreamWriter(output);

            for (int i = 0; i < nelx; i++)
            {
                for (int j = 0; j < nely; j++)
                {
                    sw.WriteLine(Xe[j * nelx + i].ToString());
                }
            }
                
            sw.Flush();
            sw.Close();
            sw.Dispose();
        }
        #endregion
    }
}
