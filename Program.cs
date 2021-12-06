using System;

namespace BESO
{
    class Program
    {
        static void Main(string[] args)
        {
            BESO beso = new BESO(3.0, 0.5);
            beso.Initialize(80, 50);

            Console.WriteLine(beso.ModelInfo());
            while (!beso.convergence)
            {
                beso.Optimize();
                Console.WriteLine(beso.info);
            }

            Console.ReadKey();
        }
    }
}
