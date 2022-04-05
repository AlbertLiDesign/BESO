using System;

namespace BESO
{
    class Program
    {
        static void Main(string[] args)
        {
            BESO beso = new BESO(3.0, 0.5);
            beso.Initialize(1000, 1000);

            Console.WriteLine(beso.ModelInfo());
            beso.Optimize();
            Console.WriteLine(beso.info);

            //while (!beso.convergence)
            //{
            //    beso.Optimize();
            //    Console.WriteLine(beso.info);

            //    //beso.WriteXe(@"E:\TestData");
            //}
            Console.WriteLine(beso.optInfo);
            Console.ReadKey();
        }
    }
}
