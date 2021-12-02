//using CSparse.Double;
//using CSparse.Storage;
//using System;
//using System.Collections.Generic;
//using System.IO;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;


//public class OutputMatrix
//{
//    public static void WriteMatrix(SparseMatrix A, string path)
//    {
//        StreamWriter sw = new StreamWriter(path);
//        sw.WriteLine("%%MatrixMarket matrix coordinate real symmetric");
//        sw.WriteLine(A.RowCount.ToString() + ' ' + A.ColumnCount.ToString() + ' ' + A.NonZerosCount.ToString());
//        int id = 0;
//        for (int i = 0; i < A.ColumnCount; i++)
//        {
//            int dif = A.ColumnPointers[i + 1] - A.ColumnPointers[i];
//            for (int j = 0; j < dif; j++)
//            {
//                if (A.RowIndices[id] >= i)
//                {
//                    sw.WriteLine((A.RowIndices[id] + 1).ToString() + ' ' + (i + 1).ToString() + ' ' + A.Values[id].ToString());
//                }
//                id++;
//            }
//        }
//        sw.Flush();
//        sw.Close();
//        sw.Dispose();
//    }
//    public static void WriteMatrix(CoordinateStorage<double> A, string path)
//    {
//        StreamWriter sw = new StreamWriter(path);
//        sw.WriteLine("%%MatrixMarket matrix coordinate real symmetric");
//        sw.WriteLine(A.RowCount.ToString() + ' ' + A.ColumnCount.ToString() + ' ' + A.NonZerosCount.ToString());
//        for (int i = 0; i < A.Values.Length; i++)
//        {
//            sw.WriteLine((A.RowIndices[i] + 1).ToString() + ' ' + (A.ColumnIndices[i] + 1).ToString() + ' ' + A.Values[i].ToString());
//        }
//        sw.Flush();
//        sw.Close();
//        sw.Dispose();
//    }

//}

