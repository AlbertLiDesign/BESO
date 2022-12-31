using System;

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
            BESO3D beso = new BESO3D(3.0, 0.5);
            beso.Initialize(10, 10, 10);

            Console.WriteLine(beso.ModelInfo());
            beso.Optimize();
            Console.WriteLine(beso.info);

            //while (!beso.convergence)
            //{
            //    beso.Optimize();
            //    Console.WriteLine(beso.info);

            //    //beso.WriteXe(@"E:\TestData");
            //}
            //Console.WriteLine(beso.optInfo);
            Console.ReadKey();
        }

        private static void testBESO2D()
        {
            BESO2D beso = new BESO2D(3.0, 0.5);
            beso.Initialize(80, 50);

            Console.WriteLine(beso.ModelInfo());
            //beso.Optimize();
            //Console.WriteLine(beso.info);

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
