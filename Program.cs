using System;

namespace BESO
{
    class Program
    {
        static void Main(string[] args)
        {
            iBESO ibeso = new iBESO(3.0);
            ibeso.Initialize(80, 50, 1);

            Console.WriteLine(ibeso.initInfo);

            //ibeso.OutputK = true;

            //ibeso.Optimize();
            //Console.WriteLine(ibeso.info);

            while (!ibeso.convergence)
            {
                ibeso.Optimize();
                Console.WriteLine(ibeso.info);
            }

            Console.WriteLine(ibeso.optInfo);
            Console.ReadKey();
        }
    }
}
