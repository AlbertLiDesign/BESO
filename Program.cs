using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BESO
{
    class Program
    {
        static void Main(string[] args)
        {
            Stopwatch stopwatch= new Stopwatch();
            stopwatch.Start();
            testBESO3D();
            stopwatch.Stop();
            Console.WriteLine("Total time:" + '\t' + stopwatch.ElapsedMilliseconds.ToString());
            Console.ReadKey();
        }

        private static void testBESO3D()
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            BESO3D beso = new BESO3D(3.0, 0.5);
            beso.Initialize(10, 8, 6, true);
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
            }

            List<int> xeNum= new List<int>();
            for (int i = 0; i < beso.Xe.Length; i++)
            {
                if (beso.Xe[i] == 1)
                {
                    xeNum.Add(i);
                }
                
            }
            beso.WriteValue("Xe.txt", xeNum.ToArray());

            //Console.WriteLine(beso.optInfo);
            Console.ReadKey();
        }
        private static void testiBESOVR()
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            BESO3D beso = new BESO3D(3.0, 0.5);
            beso.Initialize(40, 30, 20, true);
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
        private static void testBESO3DwithTime()
        {
            Stopwatch stopwatch = new Stopwatch();

            BESO3D_time beso = new BESO3D_time(3.0, 0.5, 0.02, 3, 200);
            beso.Initialize(40, 20, 30);
            beso.parallel = true;

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
        }
        private static void testBESO2DwithTime()
        {
            Stopwatch stopwatch = new Stopwatch();

            BESO2D_time beso = new BESO2D_time(3.0, 0.5, 0.02, 3, 200);
            beso.Initialize(1000, 1000);
            beso.parallel = true;

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
