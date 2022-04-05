using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BESO
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Diagnostics;
    using System.Text;
    using System.Threading.Tasks;
    using System.Runtime.InteropServices;
    using Random = System.Random;
    using System.IO;

    public class iBESO
    {
        #region 解析度参数
        public int nelx;
        public int nely;
        public int nels;
        #endregion
        public int[,] subElems;

        #region FE varialbes
        private double[] U;
        private double[] F;
        private int[] ik;
        private int[] jk;
        private double[] vk;
        private double area_w;
        #endregion
        #region BESO Parameters
        public int N;
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
        /// <summary>
        /// Design variables
        /// </summary>
        public double[] sub_Xe;

        public double[] sub_dc;
        public double[] dc;
        private double[] dc_old;
        public double[] dc_nd;
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

        private int sol_id = 0;
        private Random random;

        public bool OutputK = false;
        public bool changeSupports = true;

        public iBESO() { }
        public iBESO(int ID, double rmin, double ert = 0.02, double p = 3.0, double vf = 0.5, int maxIter = 100)
        {
            this.sol_id = ID;
            random = new Random(ID);
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
        public iBESO(int ID, iBESO ibeso)
        {
            this.sol_id = ID;
            random = new Random(ID);
            vf = ibeso.vf;
            p = ibeso.p;
            ert = ibeso.ert;
            maxIter = ibeso.maxIter;
            rmin = ibeso.rmin;

            nelx = ibeso.nelx;
            nely = ibeso.nely;
            nels = ibeso.nels;
            N = nely * nelx * nels * nels;

            sub_Xe = (double[])ibeso.sub_Xe.Clone();
            Xe = (double[])ibeso.Xe.Clone();
            dc = (double[])ibeso.dc.Clone();
            sub_dc = (double[])ibeso.sub_dc.Clone();
            dc_old = (double[])ibeso.dc_old.Clone();
            Ke = (double[])ibeso.Ke.Clone();
            subElems = (int[,])ibeso.subElems.Clone();
            area_w = 1.0 / nels / nels;

            ih = (int[])ibeso.ih.Clone();
            jh = (int[])ibeso.jh.Clone();
            vh = (double[])ibeso.vh.Clone();
            sh = (double[])ibeso.sh.Clone();

            ik = (int[])ibeso.ik.Clone();
            jk = (int[])ibeso.jk.Clone();

            stopwatch = new Stopwatch();

            optInfo = new StringBuilder("====================== Optimization ======================" + '\n');
        }

        public void Initialize(int nelx, int nely, int nels)
        {
            initInfo = new StringBuilder("====================== Launch BESO ======================" + '\n');

            this.nelx = nelx;
            this.nely = nely;
            this.nels = nels;
            area_w = 1.0 / nels / nels;

            N = nely * nelx * nels * nels;
            dc = new double[N];
            dc_old = new double[N];
            sub_dc = new double[N];
            dc = new double[nely * nelx];
            Xe = new double[nely * nelx];
            for (int j = 0; j < nely; j++)
            {
                for (int i = 0; i < nelx; i++)
                {
                    Xe[j * nelx + i] = 1.0;
                }
            }
            sub_Xe = new double[N];
            for (int j = 0; j < nely * nels; j++)
            {
                for (int i = 0; i < nelx * nels; i++)
                {
                    sub_Xe[j * nelx * nels + i] = 1.0;
                }
            }
            AssociateMeshes();

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
                stopwatch.Stop();
                #endregion
                #region Prepare report
                stopwatch.Restart();
                optInfo.Append("Flt:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion

                #region Cal nodal sensitivity
                stopwatch.Restart();
                //GetNodalSensitivity();
                stopwatch.Stop();
                #endregion
                #region Prepare report
                optInfo.Append("Cal dc_node:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion

                #region Bilinear interpolation
                stopwatch.Restart();
                if (nels != 1) BilinearInterpolation();
                else sub_dc = dc;

                stopwatch.Stop();
                #endregion
                #region Prepare report
                optInfo.Append("Bilinear interpolation:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion

                if (iter > 1)
                    for (int j = 0; j < nely * nels; j++)
                    {
                        for (int i = 0; i < nelx * nels; i++)
                        {
                            sub_dc[i * nely * nels + j] = (sub_dc[i * nely * nels + j] + dc_old[i * nely * nels + j]) * 0.5;
                        }
                    }

                // Record the sensitiveies in each step
                dc_old = (double[])sub_dc.Clone();

                #region ADD & DEL
                stopwatch.Restart(); // 计时
                ADD_DEL(vol);
                stopwatch.Stop();
                #endregion
                #region Prepare report
                optInfo.Append("ADD & DEL:" + stopwatch.Elapsed.TotalMilliseconds + '\n');
                #endregion


                #region Checking Convergence
                stopwatch.Restart(); // 计时

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
        /// <summary>
        /// Get nodal sensitivity numbers
        /// </summary>
        public void GetNodalSensitivity()
        {
            dc_nd = new double[(nelx + 1) * (nely + 1)];
            int[] sum = new int[dc_nd.Length];
            for (int i = 0; i < nelx; i++)
            {
                for (int j = 0; j < nely; j++)
                {
                    dc_nd[i * (nely + 1) + j] += dc[i * nely + j];
                    dc_nd[i * (nely + 1) + j + 1] += dc[i * nely + j];
                    dc_nd[(i + 1) * (nely + 1) + j + 1] += dc[i * nely + j];
                    dc_nd[(i + 1) * (nely + 1) + j] += dc[i * nely + j];

                    sum[i * (nely + 1) + j] += 1;
                    sum[i * (nely + 1) + j + 1] += 1;
                    sum[(i + 1) * (nely + 1) + j + 1] += 1;
                    sum[(i + 1) * (nely + 1) + j] += 1;

                }
            }

            for (int i = 0; i < dc_nd.Length; i++)
            {
                dc_nd[i] /= sum[i];
            }
        }

        public void BilinearInterpolation()
        {
            for (int elx = 0; elx < nelx; elx++)
            {
                for (int ely = 0; ely < nely; ely++)
                {
                    for (int a = 0; a < nels; a++)
                    {
                        for (int b = 0; b < nels; b++)
                        {
                            var x = (a + 0.5) / nels;
                            var y = (b + 0.5) / nels;
                            int id = subElems[elx * nely + ely, a * nels + b];
                            sub_dc[id] = dc_nd[elx * (nely + 1) + ely] * (1 - x) * (1 - y) +
                              dc_nd[elx * (nely + 1) + ely + 1] * (1 - x) * y +
                              dc_nd[(elx + 1) * (nely + 1) + ely + 1] * x * y +
                              dc_nd[(elx + 1) * (nely + 1) + ely] * x * (1 - y);
                        }
                    }

                }
            }
        }
        private void PreFlt()
        {
            int rminf = (int)Math.Floor(rmin) - 1;

            ih = new int[(int)(nely * nelx * Math.Pow((2 * rminf + 1), 2))];
            jh = new int[ih.Length];
            vh = new double[ih.Length];
            sh = new double[nely * nelx];

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
                            vh[sum] = Math.Max(0.0, rmin - Math.Sqrt((i - k) * (i - k) + (j - l) * (j - l)));
                            sum++;
                        }
                    }
                }
            }

            GetRowSum(sum, nely * nelx, ih, jh, vh, sh);
        }
        private void ADD_DEL(double volfra)
        {
            double lowest = sub_dc.Min();
            double highest = sub_dc.Max();
            double th = 0.0;
            double vol = volfra * nely * nelx;
            while (((highest - lowest) / highest) > 1e-5)
            {
                th = (highest + lowest) * 0.5;
                double sum = 0.0;
                for (int ely = 0; ely < nely; ely++)
                {
                    for (int elx = 0; elx < nelx; elx++)
                    {
                        for (int a = 0; a < nels; a++)
                        {
                            for (int b = 0; b < nels; b++)
                            {
                                int id = subElems[elx * nely + ely, a * nels + b];
                                sub_Xe[id] = sub_dc[id] > th ? 1.0 : Xmin;
                                sum += sub_Xe[id] * area_w;
                            }
                        }
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
                    int id = elx * nely + ely;
                    Compliance += 0.5 * Math.Pow(sub_Xe[id], p) * v * area_w;

                    dc[id] = 0.5 * Math.Pow(sub_Xe[id], p - 1) * v * area_w;
                    if (sol_id != 0)
                    {
                        dc[id] *= random.Next(75, 100) * 0.01;
                    }
                    //for (int a = 0; a < nels; a++)
                    //{
                    //    for (int b = 0; b < nels; b++)
                    //    {
                    //        int id = subElems[elx * nely + ely, a * nels + b];
                    //        Compliance += 0.5 * Math.Pow(sub_Xe[id], p) * v * area_w;

                    //        dc[id] = 0.5 * Math.Pow(sub_Xe[id], p - 1) * v * area_w * scored_w[id];
                    //        if (sol_id != 0)
                    //        {
                    //            dc[id] *= random.Next(75, 100) * 0.01;
                    //        }

                    //    }
                    //}
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
                    var ex = Math.Pow(Xe[j * nely + i], p);
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

        private void AssociateMeshes()
        {
            subElems = new int[nelx * nely, nels * nels];

            for (int i = 0; i < nelx; i++)
            {
                for (int j = 0; j < nely; j++)
                    for (int a = 0; a < nels; a++)
                        for (int b = 0; b < nels; b++)
                        {
                            int id = i * nely * nels * nels + a * nels * nely + j * nels + b;
                            subElems[nely * i + j, a * nels + b] = id;
                        }
            }
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
            report.Append("subCount: " + nels.ToString() + '\n');

            return report;
        }

        public void WriteXe(string path)
        {
            string output = path + '\\' + "Xe1.txt";
            StreamWriter sw = new StreamWriter(output);

            for (int i = 0; i < nelx; i++)
            {
                for (int j = 0; j < nely; j++)
                {
                    sw.WriteLine(sub_Xe[i * nely + j].ToString());
                }
            }

            sw.Flush();
            sw.Close();
            sw.Dispose();
        }
        #endregion
    }


}
