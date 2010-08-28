using System;
using System.Text;

namespace NPB {
    [Serializable()]
    public class MPIResults{
        public static void print_results(string name, 
                                         char Class, 
                                         int n1, 
                                         int n2, 
                                         int n3, 
                                         int niter, 
                                         int nprocs_compiled, 
                                         int nprocs_total, 
                                         double t, 
                                         double mops, 
                                         string optype, 
                                         bool verified, 
                                         string npbversion) {
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine(" " + name + " Benchmark Completed.");
            Console.WriteLine(" Class           = " + Class);
            //c   If this is not a grid-based problem (EP, FT, CG), then
            //c   we only print n1, which contains some measure of the
            //c   problem size. In that case, n2 and n3 are both zero.
            //c   Otherwise, we print the grid size n1xn2xn3
            if ((n2 == 0) && (n3 == 0)) {
                Console.WriteLine(" Size            = " + n1);
            }
            else {
                Console.WriteLine(" Size            = " + n1 + " x " + n2 + " x " + n3);
            }
            Console.WriteLine(" Iterations      = " + niter);
            Console.WriteLine(" Time in seconds = " + t);
            Console.WriteLine(" Total processes = " + nprocs_total);
            Console.WriteLine(" Compiled procs  = " + nprocs_compiled);
            Console.WriteLine(" Mop/s total     = " + mops);
            Console.WriteLine(" Mop/s/process   = " + mops / nprocs_total);
            Console.WriteLine(" Operation type  = " + optype);
            if (verified) {
                Console.WriteLine(" Verification    = SUCCESSFUL");
            }
            else {
                Console.WriteLine(" Verification    = UNSUCCESSFUL");
            }
            Console.WriteLine(" Version         = " + npbversion);
            Console.WriteLine("--------------------------------------------------");
            Console.WriteLine(" Please send the results of this run to:");
            Console.WriteLine(" NPB Development Team ");
            Console.WriteLine(" Internet: npb@nas.nasa.gov");
            Console.WriteLine(" If email is not available, send this to:");
            Console.WriteLine(" MS T27A-1");
            Console.WriteLine(" NASA Ames Research Center");
            Console.WriteLine(" Moffett Field, CA  94035-1000");
            Console.WriteLine(" Fax: 650-604-3957");
        }
    }
}