using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace BESO
{
    public class BESO3D_time
    {
        public bool parallel = true;
        #region Resolution
        public int nelx;
        public int nely;
        public int nelz;
        #endregion

        #region FE varialbes
        private int nEl;
        private double[] U;
        private int[] ik;
        private int[] jk;
        private double[] sk;
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
        /// Evolution rate
        /// </summary>
        public double ert;

        /// <summary>
        /// The maximum iteration
        /// </summary>
        public int maxIter;

        /// <summary>
        /// Penalty
        /// </summary>
        public int p;
        public double Xmin;
        #endregion

        #region BESO variables
        /// <summary>
        /// Design variables
        /// </summary>
        public double[] Xe;

        public double[] dc;
        public double[] dc_old;

        private int[,] cMat;
        /// <summary>
        /// Elemental stiffness matrix
        /// </summary>
        private double[] Ke0;
        private double[] Ke;
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
        private double prefetime;
        private double preflttime;
        private double featime;
        private double flttime;
        private double othertime;

        public bool OutputK = false;

        #region Constructors
        public BESO3D_time() { }
        public BESO3D_time(double rmin, double vf, double ert = 0.02, int p = 3, int maxIter = 100)
        {
            if (rmin <= 0.0)
                throw new Exception("Rmin must be large than 0.");
            if (!(vf > 0.0 && vf < 1.0))
                throw new Exception("Vt must be large than 0 and be less than 1.");

            this.vf = vf;
            this.ert = ert;
            this.maxIter = maxIter;
            this.rmin = rmin;
            this.p = p;
            Xmin = 1e-3;

            stopwatch = new Stopwatch();
            optInfo = new StringBuilder("====================== Optimization ======================" + '\n');
        }
        #endregion

        public void Initialize(int nelx, int nely, int nelz, bool parallel = false)
        {
            this.parallel = parallel;

            nEl = nelx * nely * nelz;
            initInfo = new StringBuilder("====================== Launch BESO ======================" + '\n');

            this.nelx = nelx;
            this.nely = nely;
            this.nelz = nelz;

            dc = new double[nEl];
            dc_old = new double[nEl];
            Xe = new double[nEl];

            Array.Fill(Xe, 1);

            stopwatch.Start();
            GetKe();

            PreFE3D(nelx, nely, nelz);
            stopwatch.Stop();
            initInfo.Append("PreFE: " + stopwatch.Elapsed.TotalMilliseconds + '\n');
            prefetime = stopwatch.Elapsed.TotalMilliseconds;

            stopwatch.Restart();
            PreFlt();
            stopwatch.Stop();
            initInfo.Append("PreFlt: " + stopwatch.Elapsed.TotalMilliseconds + '\n');
            preflttime = stopwatch.Elapsed.TotalMilliseconds;
        }
        public void Optimize()
        {
            if (delta > 0.001 && iter < maxIter)
            {
                iter += 1;
                vol = Math.Max(vf, vol * (1.0 - ert));

                // Record the previous sensitivity numbers
                if (iter > 1) dc_old = (double[])dc.Clone();

                #region FEA
                stopwatch.Restart();
                FE();
                stopwatch.Stop();
                featime += stopwatch.Elapsed.TotalMilliseconds;
                #endregion

                #region Get DC
                stopwatch.Restart();
                GetDc();
                HistoryC.Add(Compliance);
                stopwatch.Stop();
                othertime += stopwatch.Elapsed.TotalMilliseconds;
                #endregion

                #region Flt
                stopwatch.Restart();
                Wrapper.Flt(nEl, dc, sh);

                // Record the sensitiveies in each step
                if (iter > 1)
                    for (int i = 0; i < nEl; i++)
                        dc[i] = (dc[i] + dc_old[i]) * 0.5;

                stopwatch.Stop();
                flttime += stopwatch.Elapsed.TotalMilliseconds;
                #endregion

                #region ADD & DEL
                stopwatch.Restart();
                ADD_DEL(vol);
                stopwatch.Stop();
                othertime += stopwatch.Elapsed.TotalMilliseconds;
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
                othertime += stopwatch.Elapsed.TotalMilliseconds;
                #endregion

                info = "Iter: " + iter.ToString() + ", Volume: " + Math.Round(vol, p).ToString()
                    + ", Compliance: " + Math.Round(Compliance, p).ToString() + ", Change: " + Math.Round(delta, p).ToString();
            }
            else
            {
                convergence = true;
            }
        }
        private void PreFlt()
        {
            int rminf = (int)Math.Floor(rmin);

            ih = new int[(int)(nEl * Math.Pow((2 * rminf - 1), 3))];
            jh = new int[ih.Length];
            vh = new double[ih.Length];
            sh = new double[nEl];
            int sum = 0;
            for (int k1 = 1; k1 <= nelz; k1++)
            {
                for (int j1 = 1; j1 <= nely; j1++)
                {
                    for (int i1 = 1; i1 <= nelx; i1++)
                    {
                        var e1 = (k1 - 1) * nelx * nely + (j1 - 1) * nelx + i1;
                        for (int k2 = Math.Max(k1 - rminf + 1, 1); k2 <= Math.Min(k1 + rminf - 1, nelz); k2++)
                        {
                            for (int j2 = Math.Max(j1 - rminf + 1, 1); j2 <= Math.Min(j1 + rminf - 1, nely); j2++)
                            {
                                for (int i2 = Math.Max(i1 - rminf + 1, 1); i2 <= Math.Min(i1 + rminf - 1, nelx); i2++)
                                {
                                    var e2 = (k2 - 1) * nelx * nely + (j2 - 1) * nelx + i2;
                                    ih[sum] = e1 - 1;
                                    jh[sum] = e2 - 1;
                                    vh[sum] = Math.Max(0.0, rmin - Math.Sqrt((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2) + (k1 - k2) * (k1 - k2)));
                                    sum++;
                                }
                            }
                        }
                    }
                }
            }

            Wrapper.GetRowSum(sum, nEl, ih, jh, vh, sh);
        }
        private void ADD_DEL(double volfra)
        {
            double lowest = dc.Min();
            double highest = dc.Max();
            double th = 0.0;
            double vol = volfra * nEl;
            while (((highest - lowest) / highest) > 1e-5)
            {
                th = (highest + lowest) * 0.5;
                double sum = 0.0;
                for (int i = 0; i < nEl; i++)
                {
                    Xe[i] = dc[i] > th ? 1.0 : Xmin;
                    sum += Xe[i];
                }
                if (sum - vol > 0.0) lowest = th;
                else highest = th;
            }
        }
        private void GetDc()
        {
            Compliance = 0.0;
            if (parallel)
            {
                double[] C = new double[nEl];
                Parallel.For(0, nEl, i =>
                {
                    double[] Ue = new double[24]
                    {
                        U[cMat[i, 0]], U[cMat[i, 1]], U[cMat[i, 2]],
                        U[cMat[i, 3]], U[cMat[i, 4]], U[cMat[i, 5]],
                        U[cMat[i, 6]], U[cMat[i, 7]], U[cMat[i, 8]],
                        U[cMat[i, 9]], U[cMat[i, 10]], U[cMat[i, 11]],
                        U[cMat[i, 12]], U[cMat[i, 13]], U[cMat[i, 14]],
                        U[cMat[i, 15]], U[cMat[i, 16]], U[cMat[i, 17]],
                        U[cMat[i, 18]], U[cMat[i, 19]], U[cMat[i, 20]],
                        U[cMat[i, 21]], U[cMat[i, 22]], U[cMat[i, 23]]
                    };
                    double v = Wrapper.TransposeMultiply(24, 24, Ke0, Ue);

                    var p1 = Xe[i] == 1 ? 1.0 : 1e-9;
                    var p2 = Xe[i] == 1 ? 1.0 : 1e-6;
                    C[i] = 0.5 * p1 * v;
                    dc[i] = 0.5 * p2 * v;
                });
                Compliance = C.Sum();
            }
            else
            {
                for (int i = 0; i < nEl; i++)
                {
                    double[] Ue = new double[24]
                    {
                        U[cMat[i, 0]], U[cMat[i, 1]], U[cMat[i, 2]],
                        U[cMat[i, 3]], U[cMat[i, 4]], U[cMat[i, 5]],
                        U[cMat[i, 6]], U[cMat[i, 7]], U[cMat[i, 8]],
                        U[cMat[i, 9]], U[cMat[i, 10]], U[cMat[i, 11]],
                        U[cMat[i, 12]], U[cMat[i, 13]], U[cMat[i, 14]],
                        U[cMat[i, 15]], U[cMat[i, 16]], U[cMat[i, 17]],
                        U[cMat[i, 18]], U[cMat[i, 19]], U[cMat[i, 20]],
                        U[cMat[i, 21]], U[cMat[i, 22]], U[cMat[i, 23]]
                    };
                    double v = Wrapper.TransposeMultiply(24, 24, Ke0, Ue);

                    var p1 = Xe[i] == 1 ? 1.0 : 1e-9;
                    var p2 = Xe[i] == 1 ? 1.0 : 1e-6;
                    Compliance += 0.5 * p1 * v;
                    dc[i] = 0.5 * p2 * v;
                }
            }
        }

        #region Finite element analysis
        private void PreFE3D(int nelx, int nely, int nelz)
        {
            int[,,] nodeNrs = new int[nelz + 1, nely + 1, nelx + 1];
            int nEl = nelx * nely * nelz;
            int[] cVec = new int[nEl];
            cMat = new int[nEl, 24];
            int[] edofs = new int[24]
            {
                -3, -2, -1,
                0, 1, 2,
                3 * (nelx + 1), 3 * (nelx + 1) + 1, 3 * (nelx + 1) + 2,
                3 * (nelx + 1) - 3, 3 * (nelx + 1) - 2, 3 * (nelx + 1) - 1,

                3 * (nelx + 1) * (nely + 1) - 3, 3 * (nelx + 1) * (nely + 1) - 2, 3 * (nelx + 1) * (nely + 1) - 1,
                3 * (nelx + 1) * (nely + 1), 3 * (nelx + 1) * (nely + 1) + 1, 3 * (nelx + 1) * (nely + 1) + 2,
                3 * (nelx + 1) * (nely + 2), 3 * (nelx + 1) * (nely + 2) + 1, 3 * (nelx + 1) * (nely + 2) + 2,
                3 * (nelx + 1) * (nely + 2) - 3, 3 * (nelx + 1) * (nely + 2) - 2, 3 * (nelx + 1) * (nely + 2) - 1,
            };

            if (parallel)
            {
                Parallel.For(0, nelx, x =>
                {
                    for (int z = 0; z < nelz + 1; z++)
                        for (int y = 0; y < nely + 1; y++)
                            nodeNrs[z, y, x] = z * (nelx + 1) * (nely + 1) + y * (nelx + 1) + x;
                });
                Parallel.For(0, nelx, x =>
                {
                    for (int z = 0; z < nelz; z++)
                        for (int y = 0; y < nely; y++)
                            cVec[z * nelx * nely + y * nelx + x] = 3 * (nodeNrs[z, y, x] + 1) + 1;
                });

                Parallel.For(0, nEl, i =>
                {
                    for (int j = 0; j < 24; j++)
                        cMat[i, j] = cVec[i] + edofs[j] - 1;
                });

                var sI = new int[300];
                var sII = new int[300];
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

                ik = new int[300 * nEl];
                jk = new int[300 * nEl];

                Parallel.For(0, nEl, i =>
                {
                    for (int j = 0; j < 300; j++)
                    {
                        var a = cMat[i, sI[j]];
                        var b = cMat[i, sII[j]];
                        if (a >= b)
                        {
                            ik[i * 300 + j] = a;
                            jk[i * 300 + j] = b;
                        }
                        else
                        {
                            ik[i * 300 + j] = b;
                            jk[i * 300 + j] = a;
                        }
                    }

                });
            }
            else
            {
                for (int z = 0; z < nelz + 1; z++)
                    for (int y = 0; y < nely + 1; y++)
                        for (int x = 0; x < nelx + 1; x++)
                            nodeNrs[z, y, x] = z * (nelx + 1) * (nely + 1) + y * (nelx + 1) + x;

                for (int z = 0; z < nelz; z++)
                    for (int y = 0; y < nely; y++)
                        for (int x = 0; x < nelx; x++)
                            cVec[z * nelx * nely + y * nelx + x] = 3 * (nodeNrs[z, y, x] + 1) + 1;

                for (int i = 0; i < nEl; i++)
                    for (int j = 0; j < 24; j++)
                        cMat[i, j] = cVec[i] + edofs[j] - 1;

                var sI = new int[300];
                var sII = new int[300];
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

                ik = new int[300 * nEl];
                jk = new int[300 * nEl];

                for (int j = 0; j < nEl; j++)
                {
                    for (int i = 0; i < 300; i++)
                    {
                        var a = cMat[j, sI[i]];
                        var b = cMat[j, sII[i]];
                        if (a >= b)
                        {
                            ik[j * 300 + i] = a;
                            jk[j * 300 + i] = b;
                        }
                        else
                        {
                            ik[j * 300 + i] = b;
                            jk[j * 300 + i] = a;
                        }
                    }
                }
            }
        }
        private void FE()
        {
            // Assemble stiffness matrix with all DOFs
            sk = new double[300 * nEl];
            for (int i = 0; i < nEl; i++)
            {
                var ex = Xe[i] == 1 ? 1.0 : 1e-9;
                for (int j = 0; j < 300; j++)
                {
                    sk[300 * i + j] = ex * Ke[j];
                }
            }

            int num_allDofs = 3 * (nelx + 1) * (nely + 1) * (nelz + 1);
            var F = new double[num_allDofs];
            U = new double[num_allDofs];

            // Define force vector
            int forceID = (int)Math.Floor((nelz + 1) * 0.5) * (nelx + 1) * (nely + 1) + (int)Math.Floor((nelx + 1) * 0.5) + 1;
            F[forceID * 3 - 1] = -1.0;

            // Define fixed dofs
            int num_fixedDofs = 3 * (nelx + 1) * (nelz + 1);
            int num_freeDofs = num_allDofs - num_fixedDofs;
            var fixed_dofs = new int[num_fixedDofs];

            // my order
            for (int z = 0; z < nelz + 1; z++)
            {
                for (int x = 0; x < nelx + 1; x++)
                {
                    int id = nely * (nelx + 1) + z * (nelx + 1) * (nely + 1) + x;
                    fixed_dofs[3 * (z * (nelx + 1) + x)] = 3 * id;
                    fixed_dofs[3 * (z * (nelx + 1) + x) + 1] = 3 * id + 1;
                    fixed_dofs[3 * (z * (nelx + 1) + x) + 2] = 3 * id + 2;
                }
            }

            var all_dofs = new int[num_allDofs];
            for (int i = 0; i < num_allDofs; i++)
                all_dofs[i] = i;

            // Obtain free dofs
            free_dofs = all_dofs.Except(fixed_dofs).ToArray();

            var U_freedof = new double[num_freeDofs];

            Wrapper.Assembly_Solve(parallel, num_freeDofs, num_allDofs, ik.Length, free_dofs, ik, jk, sk, F, U_freedof);

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

            Ke = new double[300];
            Ke0 = new double[24 * 24];

            if (parallel)
            {
                Parallel.For(0, 300, i =>
                {
                    Ke[i] = E * w * (p1[i] + nu * p2[i]);
                });

                int num = 0;
                for (int i = 0; i < 24; i++)
                {
                    for (int j = 0; j < 24; j++)
                    {
                        if (i <= j)
                        {
                            Ke0[i * 24 + j] = Ke[num];
                            num++;
                        }
                    }
                }
                Parallel.For(0, 24, i =>
                {
                    for (int j = 0; j < 24; j++)
                        if (i > j)
                            Ke0[i * 24 + j] = Ke0[j * 24 + i];
                });
            }
            else
            {
                for (int i = 0; i < 300; i++)
                    Ke[i] = E * w * (p1[i] + nu * p2[i]);

                int num = 0;
                for (int i = 0; i < 24; i++)
                {
                    for (int j = 0; j < 24; j++)
                    {
                        if (i <= j)
                        {
                            Ke0[i * 24 + j] = Ke[num];
                            num++;
                        }
                    }
                }

                for (int i = 0; i < 24; i++)
                    for (int j = 0; j < 24; j++)
                        if (i > j)
                            Ke0[i * 24 + j] = Ke0[j * 24 + i];
            }
        }
        #endregion

        #region Debug Methods
        public void PrintTime()
        {
            Console.WriteLine("PreFE time:" + '\t' + prefetime.ToString());
            Console.WriteLine("PreFlt time:" + '\t' + preflttime.ToString());
            Console.WriteLine("FEA time:" + '\t' + (featime / (iter + 1)).ToString());
            Console.WriteLine("Flt time:" + '\t' + (flttime / (iter + 1)).ToString());
            Console.WriteLine("Other time:" + '\t' + (othertime / (iter + 1)).ToString());
        }
        public StringBuilder ModelInfo()
        {
            StringBuilder report = new StringBuilder("=================== Model Info ===================" + '\n');
            report.Append("Nodes: " + ((nelx + 1) * (nely + 1) * (nelz + 1)).ToString() + '\n');
            report.Append("Elements: " + (nelx * nely * nelz).ToString() + '\n');
            report.Append("Parallel mode: " + parallel.ToString() + '\n');
            report.Append('\n');
            report.Append("=================== Parameters Info ===================" + '\n');
            report.Append("xCount: " + nelx.ToString() + '\n');
            report.Append("yCount: " + nely.ToString() + '\n');
            report.Append("zCount: " + nelz.ToString() + '\n');
            return report;
        }

        private void WriteValue(string path, double[] a)
        {
            StreamWriter sw = new StreamWriter(path);

            for (int i = 0; i < a.Length; i++)
            {
                sw.WriteLine(a[i].ToString());
            }

            sw.Flush();
            sw.Close();
            sw.Dispose();

        }
        private void WriteValue(string path, int[] a)
        {
            StreamWriter sw = new StreamWriter(path);

            for (int i = 0; i < a.Length; i++)
            {
                sw.WriteLine(a[i].ToString());
            }

            sw.Flush();
            sw.Close();
            sw.Dispose();
        }
        #endregion
    }
}
