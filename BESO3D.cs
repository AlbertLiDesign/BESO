using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using System.IO;
using MathNet.Numerics.LinearAlgebra.Double;

namespace BESO
{
    public class BESO3D
    {
        #region Resolution
        public int nelx;
        public int nely;
        public int nelz;
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
        public int[] kh;
        public double[] vh;
        public double[] sh;
        #endregion

        private Stopwatch stopwatch;

        public bool OutputK = false;
        public bool changeSupports = true;
        public BESO3D() { }
        public BESO3D(double rmin, double vf, double ert = 0.02, double p = 3.0, int maxIter = 100)
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
        public BESO3D(BESO3D beso)
        {
            vf = beso.vf;
            p = beso.p;
            ert = beso.ert;
            maxIter = beso.maxIter;
            rmin = beso.rmin;

            nelx = beso.nelx;
            nely = beso.nely;

            Xe = (double[])beso.Xe.Clone();
            dc = (double[])beso.dc.Clone();
            dc_old = (double[])beso.dc_old.Clone();
            Ke = (double[])beso.Ke.Clone();

            ih = (int[])beso.ih.Clone();
            jh = (int[])beso.jh.Clone();
            vh = (double[])beso.vh.Clone();
            sh = (double[])beso.sh.Clone();

            ik = (int[])beso.ik.Clone();
            jk = (int[])beso.jk.Clone();

            stopwatch = new Stopwatch();

            optInfo = new StringBuilder("====================== Optimization ======================" + '\n');
        }


        public void Initialize(int nelx, int nely, int nelz)
        {
            initInfo = new StringBuilder("====================== Launch BESO ======================" + '\n');

            this.nelx = nelx;
            this.nely = nely;
            this.nelz = nelz;

            dc = new double[nelx * nely * nelz];
            dc_old = new double[nelx * nely * nelz];
            Xe = new double[nelx * nely * nelz];

            Array.Fill(Xe, 1);

            stopwatch.Start();
            GetKe();
            ik = new int[nelx * nely * nelz * 24 * 24];
            jk = new int[nelx * nely * nelz * 24 * 24];

            PreFE3D(nelx, nely, nelz, ik, jk);
            stopwatch.Stop();
            initInfo.Append("PreFE: " + stopwatch.Elapsed.TotalMilliseconds + '\n');

            stopwatch.Restart();
            //PreFlt();
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

                // Record the previous sensitivity numbers
                if(iter > 1) dc_old = (double[])dc.Clone();

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
                //Flt(dc.Length, dc, sh);

                if (iter > 1)
                    for (int j = 0; j < nely; j++)
                    {
                        for (int i = 0; i < nelx; i++)
                        {
                            dc[i * nely + j] = (dc[i * nely + j] + dc_old[i * nely + j]) * 0.5;
                        }
                    }

                // Record the sensitiveies in each step
                
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
            //int rminf = (int)Math.Floor(rmin);

            //ih = new int[(int)(nelx * nely * Math.Pow((2 * rminf + 1), 2))];
            //jh = new int[ih.Length];
            //vh = new double[ih.Length];
            //sh = new double[nelx * nely];

            //int sum = 0;
            //for (int i = 0; i < nelx; i++)
            //{
            //    for (int j = 0; j < nely; j++)
            //    {
            //        var e1 = i * nely + j + 1;
            //        for (int k = Math.Max(i - rminf, 0); k < Math.Min(i + rminf + 1, nelx); k++)
            //        {
            //            for (int l = Math.Max(j - rminf, 0); l < Math.Min(j + rminf + 1, nely); l++)
            //            {
            //                var e2 = k * nely + l + 1;
            //                ih[sum] = e1 - 1;
            //                jh[sum] = e2 - 1;
            //                vh[sum] = Math.Max(0.0, rminf - Math.Sqrt((i - k) * (i - k) + (j - l) * (j - l)));
            //                sum++;
            //            }
            //        }
            //    }
            //}

            //GetRowSum(sum, nelx * nely, ih, jh, vh, sh);
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
                        Xe[j * nelx + i] = dc[i * nely + j] > th ? 1.0 : Xmin;
                        sum += Xe[j * nelx + i];
                    }
                }
                if (sum - vol > 0.0) lowest = th;
                else highest = th;
            }
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

                    //double v = TransposeMultiply(8, 8, Ke, Ue);

                    //Compliance += 0.5 * Math.Pow(Xe[ely * nelx + elx], p) * v;
                    //dc[elx * nely + ely] = 0.5 * Math.Pow(Xe[ely * nelx + elx], p - 1) * v;
                }
            }
        }

        private void PreFE3D(int nelx, int nely, int nelz, int[] ik, int[] jk)
        {
            int[,,] nodeNrs = new int[nelz + 1, nely + 1, nelx + 1];
            int[] cVec = new int[nelx * nely * nelz];
            double[,] cMat = new double[nelx* nely *nelz, 24];
            int[] edofs = new int[24]{ 
                0, 1, 2, 3 * (nely + 1) * (nelz + 1), 3 * (nely + 1) * (nelz + 1) + 1,
                3 * (nely + 1) * (nelz + 1) + 2, 3 * (nely + 1) * (nelz + 1) - 3, 3 * (nely + 1) * (nelz + 1) - 2,
                3 * (nely + 1) * (nelz + 1) - 1, -3, -2, -1, 3 * (nely + 1), 3 * (nely + 1) + 1, 3 * (nely + 1) + 2,
                3 * (nely + 1) * (nelz + 2), 3 * (nely + 1) * (nelz + 2) + 1, 3 * (nely + 1) * (nelz + 2) + 2,
                3 * (nely + 1) * (nelz + 2) -3, 3 * (nely + 1) * (nelz + 2) - 2, 3 * (nely + 1) * (nelz + 2) - 1,
                3 * (nely + 1) -3, 3 * (nely + 1) - 2, 3 * (nely + 1) -1
                };

            for (int z = 0; z < nelz + 1; z++)
            {
                for (int y = 0; y < nely + 1; y++)
                {
                    for (int x = 0; x < nelx + 1; x++)
                    {
                        nodeNrs[z, y, x] = x * (nely + 1) * (nelz + 1) + y * (nelz + 1) + z;
                    }
                }
            }

            for (int z = 0; z < nelz; z++)
            {
                for (int y = 0; y < nely; y++)
                {
                    for (int x = 0; x < nelx; x++)
                    {
                        cVec[x * nely * nelz + y * nelz + z] = 3 * nodeNrs[z, y, x] + 1;
                    }
                }
            }

            for (int i = 0; i < nelx * nely * nelz; i++)
            {
                for (int j = 0; j < 24; j++)
                {
                    cMat[i, j] = cVec[i] + edofs[j];
                }
            }

            var sI = new double[300];
            var sII = new double[300];
            int num = 0;
            int num2 = 0;
            for (int i = 0; i < 24; i++)
            {
                for (int j = i; j < 24; j++)
                {
                    sI[num] = j;
                    sII[num] = num2;
                    num++;
                }
                num2++;
            }


            //var mat = DenseMatrix.OfArray(cMat);
            //Console.WriteLine(mat);
        }
        private void FE()
        {
            int num_allDofs = 3 * (nelx + 1) * (nely + 1) * (nelz + 1);
            int num_fixedDofs = 3 * (nely + 1) * (nelx + 1);
            int num_freeDofs = num_allDofs - num_fixedDofs;

            var F = new double[num_allDofs];
            U = new double[num_allDofs];

            // Define force vector
            int forceID = (int)Math.Floor((nely + 1) * (nelx + 1) * 0.5);
            F[forceID * 3] = -1.0;

            if (changeSupports)
            {
                // Define fixed dofs
                var fixed_dofs = new int[num_fixedDofs];
                for (int i = 0; i < num_fixedDofs; i++)
                    fixed_dofs[i] = num_allDofs - 1 - i;

                var all_dofs = new int[num_allDofs];
                for (int i = 0; i < num_allDofs; i++)
                    all_dofs[i] = i;

                // Obtain free dofs
                free_dofs = all_dofs.Except(fixed_dofs).ToArray();
                changeSupports = false;
            }

            var U_freedof = new double[num_freeDofs];
            //Assembly_Solve(num_freeDofs, num_allDofs, ik.Length, free_dofs, ik, jk, vk, F, U_freedof);

            for (int i = 0; i < num_freeDofs; i++)
            {
                U[free_dofs[i]] = U_freedof[i];
            }
        }

        private void GetKe()
        {
            var E = 1.0;
            var nu = 0.3;

            var w = 1 / (1 + nu) / (2 * nu - 1) / 144;
            var p1 = new double[300]{-32, -6, -6, 8, 6, 6, 10, 6, 3, -4, -6, -3, -4, -3, -6, 10,
            3, 6, 8, 3, 3, 4, -3, -3, -32, -6, -6, -4, -3, 6, 10, 3, 6, 8, 6, -3, -4, -6, -3, 4, -3, 3, 8, 3,
            3, 10, 6, -32, -6, -3, -4, -3, -3, 4, -3, -6, -4, 6, 6, 8, 6, 3, 10, 3, 3, 8, 3, 6, 10, -32, 6, 6,
            -4, 6, 3, 10, -6, -3, 10, -3, -6, -4, 3, 6, 4, 3, 3, 8, -3, -3, -32, -6, -6, 8, 6, -6, 10, 3, 3, 4,
            -3, 3, -4, -6, -3, 10, 6, -3, 8, 3, -32, 3, -6, -4, 3, -3, 4, -6, 3, 10, -6, 6, 8, -3, 6, 10, -3,
            3, 8, -32, -6, 6, 8, 6, -6, 8, 3, -3, 4, -3, 3, -4, -3, 6, 10, 3, -6, -32, 6, -6, -4, 3, 3, 8, -3,
            3, 10, -6, -3, -4, 6, -3, 4, 3, -32, 6, 3, -4, -3, -3, 8, -3, -6, 10, -6, -6, 8, -6, -3, 10, -32,
            6, -6, 4, 3, -3, 8, -3, 3, 10, -3, 6, -4, 3, -6, -32, 6, -3, 10, -6, -3, 8, -3, 3, 4, 3, 3, -4, 6,
            -32, 3, -6, 10, 3, -3, 8, 6, -3, 10, 6, -6, 8, -32, -6, 6, 8, 6, -6, 10, 6, -3, -4, -6, 3, -32, 6,
            -6, -4, 3, 6, 10, -3, 6, 8, -6, -32, 6, 3, -4, 3, 3, 4, 3, 6, -4, -32, 6, -6, -4, 6, -3, 10, -6, 3,
            -32, 6, -6, 8, -6, -6, 10, -3, -32, -3, 6, -4, -3, 3, 4, -32, -6, -6, 8, 6, 6, -32, -6, -6, -4,
            -3, -32, -6, -3, -4, -32, 6, 6, -32, -6, -32};

            var p2 = new double[300]{48, 0, 0, 0, -24, -24, -12, 0, -12, 0,
            24, 0, 0, 0, 24, -12, -12, 0, -12, 0, 0, -12, 12, 12, 48, 0, 24, 0, 0, 0, -12, -12, -24, 0, -24,
            0, 0, 24, 12, -12, 12, 0, -12, 0, -12, -12, 0, 48, 24, 0, 0, 12, 12, -12, 0, 24, 0, -24, -24, 0,
            0, -12, -12, 0, 0, -12, -12, 0, -12, 48, 0, 0, 0, -24, 0, -12, 0, 12, -12, 12, 0, 0, 0, -24,
            -12, -12, -12, -12, 0, 0, 48, 0, 24, 0, -24, 0, -12, -12, -12, -12, 12, 0, 0, 24, 12, -12, 0,
            0, -12, 0, 48, 0, 24, 0, -12, 12, -12, 0, -12, -12, 24, -24, 0, 12, 0, -12, 0, 0, -12, 48, 0, 0,
            0, -24, 24, -12, 0, 0, -12, 12, -12, 0, 0, -24, -12, -12, 0, 48, 0, 24, 0, 0, 0, -12, 0, -12,
            -12, 0, 0, 0, -24, 12, -12, -12, 48, -24, 0, 0, 0, 0, -12, 12, 0, -12, 24, 24, 0, 0, 12, -12,
            48, 0, 0, -12, -12, 12, -12, 0, 0, -12, 12, 0, 0, 0, 24, 48, 0, 12, -12, 0, 0, -12, 0, -12, -12,
            -12, 0, 0, -24, 48, -12, 0, -12, 0, 0, -12, 0, 12, -12, -24, 24, 0, 48, 0, 0, 0, -24, 24, -12,
            0, 12, 0, 24, 0, 48, 0, 24, 0, 0, 0, -12, 12, -24, 0, 24, 48, -24, 0, 0, -12, -12, -12, 0, -24,
            0, 48, 0, 0, 0, -24, 0, -12, 0, -12, 48, 0, 24, 0, 24, 0, -12, 12, 48, 0, -24, 0, 12, -12, -12,
            48, 0, 0, 0, -24, -24, 48, 0, 24, 0, 0, 48, 24, 0, 0, 48, 0, 0, 48, 0, 48};

            var p = new double[300];
            for (int i = 0; i < 300; i++)
                p[i] = w * (p1[i] + nu * p2[i]);


            Ke = new double[24 * 24];

            int num = 0;
            for (int i = 0; i < 24; i++)
            {
                for (int j = 0; j < 24; j++)
                {
                    if (i<= j)
                    {
                        Ke[i * 24 + j] = p[num];
                        num++;
                    }
                }
            }

            for (int i = 0; i < 24; i++)
            {
                for (int j = 0; j < 24; j++)
                {
                    if (i > j)
                    {
                        Ke[i * 24 + j] = Ke[j * 24 + i];
                    }
                }
            }
        }

        [DllImport("Solver.dll")]
        private static extern void WriteMatrix(int row, int col, double[] val);

        //[DllImport("Solver.dll")]
        //private static extern void PreFE3D(int nelx, int nely, int nelz, int[] ik, int[] jk);

        #region Debug Methods
        public StringBuilder ModelInfo()
        {
            StringBuilder report = new StringBuilder("=================== Model Info ===================" + '\n');
            report.Append("Opt Nodes: " + ((nelx + 1) * (nely + 1) * (nelz + 1)).ToString() + '\n');
            report.Append("Opt Elements: " + (nelx * nely * nelz).ToString() + '\n');
            report.Append('\n');
            report.Append("=================== Parameters Info ===================" + '\n');
            report.Append("xCount: " + nelx.ToString() + '\n');
            report.Append("yCount: " + nely.ToString() + '\n');
            report.Append("zCount: " + nelz.ToString() + '\n');
            return report;
        }

        //public void WriteXe(string path)
        //{
        //    string output = path + '\\' + "Xe2.txt";
        //    StreamWriter sw = new StreamWriter(output);

        //    for (int i = 0; i < nelx; i++)
        //    {
        //        for (int j = 0; j < nely; j++)
        //        {
        //            sw.WriteLine(Xe[j * nelx + i].ToString());
        //        }
        //    }

        //    sw.Flush();
        //    sw.Close();
        //    sw.Dispose();
        //}
        #endregion
    }
}
