using System;
using System.Diagnostics;

namespace BESO
{
    class Program
    {
        static void Main(string[] args)
        {
            Stopwatch stopwatch= new Stopwatch();
            stopwatch.Start();
            testBESO2DwithTime();
            stopwatch.Stop();
            Console.WriteLine("Total time:" + '\t' + stopwatch.ElapsedMilliseconds.ToString());
            Console.ReadKey();
        }

        private static void testBESO3D()
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            BESO3D beso = new BESO3D(3.0, 0.5);
            beso.Initialize(80, 50, 20, true);
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

        private static void testBESO2DwithTime()
        {
            Stopwatch stopwatch = new Stopwatch();

            BESO2D_time beso = new BESO2D_time(3.0, 0.5, 0.02, 3, 200);
            beso.Initialize(640, 400);

            Console.WriteLine(beso.ModelInfo());

            while (!beso.convergence)
            {
                stopwatch.Restart();
                beso.Optimize();
                Console.WriteLine(beso.info);
                stopwatch.Stop();
                Console.WriteLine(stopwatch.ElapsedMilliseconds.ToString());
                beso.PrintTime();
            }
            beso.PrintTime();
            //Console.WriteLine(beso.optInfo);
        }
        private static void testBESO2D()
        {
            Stopwatch stopwatch= new Stopwatch();

            BESO2D beso = new BESO2D(3.0, 0.5, 0.02, 3, 200);
            beso.Initialize(1000, 1000);

            Console.WriteLine(beso.ModelInfo());

            while (!beso.convergence)
            {
                stopwatch.Restart();
                beso.Optimize();
                Console.WriteLine(beso.info);
                stopwatch.Stop();
                Console.WriteLine(stopwatch.ElapsedMilliseconds.ToString());
            }
            //Console.WriteLine(beso.optInfo);
        }
    }
}
