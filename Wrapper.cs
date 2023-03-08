using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace BESO
{
    public class Wrapper
    {
        [DllImport("Solver.dll")]
        public static extern void PrintMatrix(int row, int col, double[] val);
        [DllImport("Solver.dll")]
        public static extern void Assembly_Solve_AMG(bool parallel, int num_freeDofs, int num_allDofs, int num_triplets, int[] free_dofs, int[] ik, int[] jk, double[] vk, double[] F, double[] U);

        [DllImport("Solver.dll")]
        public static extern void Assembly_Solve(bool parallel, int num_freeDofs, int num_allDofs, int num_triplets, int[] free_dofs, int[] ik, int[] jk, double[] vk, double[] F, double[] U);
        [DllImport("Solver.dll")]
        public static extern double TransposeMultiply(int rows, int cols, double[] A, double[] U);

        [DllImport("Solver.dll")]
        public static extern void PreFE(int nelx, int nely, int[] ik, int[] jk);

        [DllImport("Solver.dll")]
        public static extern void Flt(int dc_length, double[] dc, double[] sh);
        [DllImport("Solver.dll")]
        public static extern void Flt3D(int nEl, double[] dc, double[] sh);
        [DllImport("Solver.dll")]
        public static extern void GetRowSum(int coo_length, int rows, int[] ih, int[] jh, double[] vh, double[] sh);

        [DllImport("Solver.dll")]
        public static extern void Cal_ik_jk(int nEl, int[] edofMat, int[] ik, int[] jk);
    }

}
