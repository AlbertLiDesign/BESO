using System;
using System.Diagnostics;

namespace BESO
{
    class Program
    {
        static void Main(string[] args)
        {
            testBESO3D();
        }

        private static void testBESO3D()
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            BESO3D beso = new BESO3D(3.0, 0.5);
            beso.Initialize(40, 30, 20, false);
            Console.WriteLine(beso.ModelInfo());

            sw.Stop();
            Console.WriteLine(
                "======================= Init. time: "
                + sw.ElapsedMilliseconds.ToString()
                + " =======================");

            while (!beso.convergence)
            {
                sw.Restart();
                beso.Optimize();
                sw.Stop();
                Console.WriteLine(beso.info);
                Console.WriteLine(
                    "======================= It. time: " 
                    + sw.ElapsedMilliseconds.ToString() 
                    + " =======================");
                //beso.WriteXe(@"E:\TestData");
            }
            //Console.WriteLine(beso.optInfo);
            Console.ReadKey();
        }

        private static void testBESO2D()
        {
            BESO2D beso = new BESO2D(3.0, 0.5);
            beso.Initialize(8, 5);

            Console.WriteLine(beso.ModelInfo());

            while (!beso.convergence)
            {
                beso.Optimize();
                Console.WriteLine(beso.info);

                //beso.WriteXe(@"E:\TestData");
            }
            Console.WriteLine(beso.optInfo);
            Console.ReadKey();
        }
    }
}
