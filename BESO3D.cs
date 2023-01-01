using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using System.IO;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Globalization;

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
        #endregion

        #region BESO variables
        /// <summary>
        /// Design variables
        /// </summary>
        public double[] Xe;

        public double[] dc;
        public double[] dc_old;

        private int[] edofMat;
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

        public bool OutputK = false;

        #region Constructors
        public BESO3D() { }
        public BESO3D(double rmin, double vf, double ert = 0.02, int maxIter = 100)
        {
            if (rmin <= 0.0)
                throw new Exception("Rmin must be large than 0.");
            if (!(vf > 0.0 && vf < 1.0))
                throw new Exception("Vt must be large than 0 and be less than 1.");

            this.vf = vf;
            this.ert = ert;
            this.maxIter = maxIter;
            this.rmin = rmin;

            stopwatch = new Stopwatch();
            optInfo = new StringBuilder("====================== Optimization ======================" + '\n');
        }
        public BESO3D(BESO3D beso)
        {
            vf = beso.vf;
            ert = beso.ert;
            maxIter = beso.maxIter;
            rmin = beso.rmin;

            nelx = beso.nelx;
            nely = beso.nely;

            Xe = (double[])beso.Xe.Clone();
            dc = (double[])beso.dc.Clone();
            dc_old = (double[])beso.dc_old.Clone();
            Ke0 = (double[])beso.Ke0.Clone();

            ih = (int[])beso.ih.Clone();
            jh = (int[])beso.jh.Clone();
            vh = (double[])beso.vh.Clone();
            sh = (double[])beso.sh.Clone();

            ik = (int[])beso.ik.Clone();
            jk = (int[])beso.jk.Clone();

            stopwatch = new Stopwatch();

            optInfo = new StringBuilder("====================== Optimization ======================" + '\n');
        }
        #endregion

        public void Initialize(int nelx, int nely, int nelz)
        {
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
                Wrapper.Flt(dc.Length, dc, sh);

                if (iter > 1)
                    for (int k = 0; k < nelz; k++)
                        for (int j = 0; j < nely; j++)
                            for (int i = 0; i < nelx; i++)
                            {
                                int id = i * nely * nelz + j * nelz + k;
                                dc[id] = (dc[id] + dc_old[id]) * 0.5;
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
            int rminf = (int)Math.Floor(rmin);

            ih = new int[(int)(nEl * Math.Pow((3 * rminf - 1), 3))];
            jh = new int[ih.Length];
            vh = new double[ih.Length];
            sh = new double[nEl];

            int sum = 0;
            for (int k1 = 1; k1 <= nelz; k1++)
            {
                for (int i1 = 1; i1 <= nelx; i1++)
                {
                    for (int j1 = 1; j1 <= nely; j1++)
                    {
                        var e1 = (k1 - 1) * nelx * nely + (i1 - 1) * nely + j1;
                        //var e1 = (i1 - 1) * nely * nelz + (j1 - 1) * nelz + k1;
                        for (int k2 = Math.Max(k1 - rminf + 1, 1); k2 <= Math.Min(k1 + rminf - 1, nelz); k2++)
                        {
                            for (int i2 = Math.Max(i1 - rminf + 1, 1); i2 <= Math.Min(i1 + rminf - 1, nelx); i2++)
                            {
                                for (int j2 = Math.Max(j1 - rminf + 1, 1); j2 <= Math.Min(j1 + rminf - 1, nely); j2++)
                                {
                                    var e2 = (k2 - 1) * nelx * nely + (i2 - 1) * nely + j2;
                                    //var e2 = (i2 - 1) * nely * nelz + (j2 - 1) * nelz + k2;
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
                for (int k = 0; k < nelz; k++)
                {
                    for (int j = 0; j < nely; j++)
                    {
                        for (int i = 0; i < nelx; i++)
                        {
                            int id = i * nely * nelz + j * nelz + k;
                            Xe[id] = dc[id] > th ? 1.0 : 1e-3;
                            sum += Xe[j * nelx + i];
                        }
                    }
                }
                if (sum - vol > 0.0) lowest = th;
                else highest = th;
            }
        }
        private void GetDc()
        {
            Compliance = 0.0;
            for (int i = 0; i < nEl; i++)
            {
                double[] Ue = new double[24]
                {
                   U[edofMat[i + nEl * 0]], U[edofMat[i + nEl * 1]], U[edofMat[i + nEl * 2]],
                   U[edofMat[i + nEl * 3]], U[edofMat[i + nEl * 4]], U[edofMat[i + nEl * 5]],
                   U[edofMat[i + nEl * 6]], U[edofMat[i + nEl * 7]], U[edofMat[i + nEl * 8]],
                   U[edofMat[i + nEl * 9]], U[edofMat[i + nEl * 10]], U[edofMat[i + nEl * 11]],
                   U[edofMat[i + nEl * 12]], U[edofMat[i + nEl * 13]], U[edofMat[i + nEl * 14]],
                   U[edofMat[i + nEl * 15]], U[edofMat[i + nEl * 16]], U[edofMat[i + nEl * 17]],
                   U[edofMat[i + nEl * 18]], U[edofMat[i + nEl * 19]], U[edofMat[i + nEl * 20]],
                   U[edofMat[i + nEl * 21]], U[edofMat[i + nEl * 22]], U[edofMat[i + nEl * 23]]
                };
                double v = Wrapper.TransposeMultiply(24, 24, Ke0, Ue);

                var p1 = Xe[i] == 1 ? 1.0 : 1e-9;
                var p2 = Xe[i] == 1 ? 1.0 : 1e-6;
                Compliance += 0.5 * p1 * v;
                dc[i] = 0.5 * p2 * v;
            }
        }

        #region Finite element analysis
        private void PreFE3D(int nelx, int nely, int nelz)
        {
            var nodegrd = new int[nelx + 1, nely + 1];
            for (int i = 0; i < nelx + 1; i++)
                for (int j = 0; j < nely + 1; j++)
                    nodegrd[i, j] = i * (nely + 1) + j;

            var nodeids = new int[nelx * nely];
            for (int i = 0; i < nelx; i++)
                for (int j = 0; j < nely; j++)
                    nodeids[i * nely + j] = nodegrd[i, j];

            var nodeidz = new int[nelz];
            for (int i = 0; i < nelz; i++)
                nodeidz[i] = i * (nelx + 1) * (nely + 1);

            var edofVec = new int[nEl];
            for (int i = 0; i < nelx * nely; i++)
            {
                for (int j = 0; j < nelz; j++)
                {
                    edofVec[j * nelx * nely + i] = (nodeids[i] + 1 + nodeidz[j]) * 3;
                }
            }

            int a = 3 * nely;
            int b = 3 * (nely + 1) * (nelx + 1);
            int[] edofs = new int[24]
            {
                0, 1, 2,
                a + 3, a + 4, a + 5,
                a, a + 1, a + 2,
                -3, -2, -1,

                b, b + 1,b + 2,
                b + a + 3, b + a + 4, b + a + 5,
                b + a, b + a + 1, b + a + 2,
                b - 3, b - 2, b - 1
            };

            edofMat = new int[nEl * 24];
            for (int i = 0; i < 24; i++)
            {
                for (int j = 0; j < nEl; j++)
                {
                    edofMat[i * nEl + j] = edofVec[j] + edofs[i];
                }
            }

            ik = new int[24 * 24 * nEl];
            jk = new int[24 * 24 * nEl];

            Wrapper.Cal_ik_jk(nEl, edofMat, ik, jk);
        }
        private void PreFE3D2(int nelx, int nely, int nelz)
        {
            int[,,] nodeNrs = new int[nelx + 1, nely + 1, nelz + 1];

            for (int z = 0; z < nelz + 1; z++)
            {
                for (int y = 0; y < nely + 1; y++)
                {
                    for (int x = 0; x < nelx + 1; x++)
                    {
                        nodeNrs[x, y, z] = y * (nelx + 1) * (nelz + 1) + x * (nelz + 1) + z;
                        //nodeNrs[x, y, z] = x * (nely + 1) * (nelz + 1) + y * (nelz + 1) + z;
                    }
                }
            }

            int[] cVec = new int[nEl];
            for (int z = 0; z < nelz; z++)
            {
                for (int y = 0; y < nely; y++)
                {
                    for (int x = 0; x < nelx; x++)
                    {
                        cVec[y * nelx * nelz + x * nelz + z] = 3 * (nodeNrs[x, y, z] + 1) + 1;
                        //cVec[x * nely * nelz + y * nelz + z] = 3 * (nodeNrs[x, y, z] + 1) + 1;
                    }
                }
            }

            int[,] cMat = new int[nEl, 24];

            int[] edofs = new int[24]{
                -3, -2, -1,
                3 * (nelz + 1) -3, 3 * (nelz + 1) - 2, 3 * (nelz + 1) -1,
                3 * (nelz + 1) * (nelx + 2) -3, 3 * (nelz + 1) * (nelx + 2) - 2, 3 * (nelz + 1) * (nelx + 2) - 1,
                3 * (nelx + 1) * (nelz + 1) - 3, 3 * (nelx + 1) * (nelz + 1) - 2, 3 * (nelx + 1) * (nelz + 1) - 1,

                0, 1, 2,
                3 * (nelz + 1), 3 * (nelz + 1) + 1, 3 * (nelz + 1) + 2,
                3 * (nelz + 1) * (nelx + 2), 3 * (nelz + 1) * (nelx + 2) + 1, 3 * (nelz + 1) * (nelx + 2) + 2,
                3 * (nelz + 1) * (nelx + 1), 3 * (nelz + 1) * (nelx + 1) + 1, 3 * (nelz + 1) * (nelx + 1) + 2
                };
          

            for (int i = 0; i < nEl; i++)
            {
                for (int j = 0; j < 24; j++)
                {
                    cMat[i, j] = cVec[i] + edofs[j];
                }
            }

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
                    var a = cMat[j, sI[i]] - 1;
                    var b = cMat[j, sII[i]] - 1;
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
        private void FE()
        {
            // Assemble stiffness matrix with all DOFs
            sk = new double[24 * 24 * nEl];
            for (int z = 0; z < nelz; z++)
            {
                for (int y = 0; y < nely; y++)
                {
                    for (int x = 0; x < nelx; x++)
                    {
                        var id = x * nely * nelz + y * nelz + z;
                        var ex = Xe[id] == 1 ? 1.0 : 1e-9;
                        for (int a = 0; a < 24; a++)
                            for (int b = 0; b < 24; b++)
                                sk[24 * (24 * id + a) + b] = ex * Ke0[b * 24 + a];
                    }
                }
            }

            int num_allDofs = 3 * (nelx + 1) * (nely + 1) * (nelz + 1);
            int num_fixedDofs = 3 * (nely + 1) * (nelz + 1);
            int num_freeDofs = num_allDofs - num_fixedDofs;

            var F = new double[num_allDofs];
            U = new double[num_allDofs];

            // Define force vector
            //int forceID = (int)Math.Floor((nely + 1) * (nelx + 1) * 0.5) + 1;
            for (int z = 0; z < nelz + 1; z++)
            {
                var loadnid = z * (nelx + 1) * (nely + 1) + nelx * (nely + 1) + (nely + 1 - 0);
                F[loadnid * 3 - 1 - 1] = -1.0;
            }

            // Define fixed dofs
            var fixed_dofs = new int[num_fixedDofs];
            //for (int i = 0; i < num_fixedDofs; i++)
            //    fixed_dofs[i] = num_allDofs - 1 - i;
            int sum = 0;
            for (int z = 0; z < nelz + 1; z++)
            {
                for (int y = 0; y < nely + 1; y++)
                {
                    fixed_dofs[sum * 3] = 3 * (z * (nelx + 1) * (nely + 1) + (nely + 1 - y)) - 1;
                    fixed_dofs[sum * 3 + 1] = 3 * (z * (nelx + 1) * (nely + 1) + (nely + 1 - y)) - 2;
                    fixed_dofs[sum * 3 + 2] = 3 * (z * (nelx + 1) * (nely + 1) + (nely + 1 - y)) - 3;
                    sum++;
                }
            }

            var all_dofs = new int[num_allDofs];
            for (int i = 0; i < num_allDofs; i++)
                all_dofs[i] = i;

            // Obtain free dofs
            free_dofs = all_dofs.Except(fixed_dofs).ToArray();

            var U_freedof = new double[num_freeDofs];

            Wrapper.Assembly_Solve(num_freeDofs, num_allDofs, ik.Length, free_dofs, ik, jk, sk, F, U_freedof);

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
            for (int i = 0; i < 300; i++)
                Ke[i] = E * w * (p1[i] + nu * p2[i]);


            Ke0 = new double[24 * 24];

            int num = 0;
            for (int i = 0; i < 24; i++)
            {
                for (int j = 0; j < 24; j++)
                {
                    if (i<= j)
                    {
                        Ke0[i * 24 + j] = Ke[num];
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
                        Ke0[i * 24 + j] = Ke0[j * 24 + i];
                    }
                }
            }
        }
        #endregion

        private static int[,] KroneckerProduct(int[,] A , int[,] B)
        {
            int bHeight = B.GetLength(0);
            int bWidth = B.GetLength(1);
            int height = A.GetLength(0) * bHeight;
            int width = A.GetLength(1) * bWidth;
            var result = new int[height, width];
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    result[y, x] = (A[y / bHeight, x / bWidth] - 1) * B[y % bHeight, x % bWidth];
                }
            }
            return result;
        }

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
