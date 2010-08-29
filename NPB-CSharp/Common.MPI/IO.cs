using System;
using System.Text;

namespace NPB {
    [Serializable()]
    public class IO{
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
        public static string[] readFileData(string path, int[] conf) {
            String s, temp = ""; int size = 0;
            for (int i = 0; i < conf.Length; i++) { if (conf[i]!=0) size = size + conf[i]; }
            string[] vet = new string[size];
            //conf is the quantity of variable in each line.
            //LU int[] conf = {0,0,2,0,0,1,0,0,1,0,0,1,0,0,5,0,0,3}; // conf numbers: first line have 1 elements, second line 1 elements, third line 3 elements... 
            //BT int[] conf = { 1, 1, 3, 2, 2 };                     // conf numbers: first line have 1 elements, second line 1 elements, third line 3 elements... 
            //FT int[] conf = {1,1,2};
            int count = 0;
            System.IO.StreamReader file = new System.IO.StreamReader(path);
            for (int k = 0; k < conf.Length; k++) {
                s = file.ReadLine(); s = s.Trim(' ');
                for (int j = 0; j < conf[k]; j++) {
                    for (int i = 0; i <= s.Length; i++) {
                        if (s.Length == 0) {
                            break;
                        }
                        else { if (s[0] == ' ') break; }
                        temp = temp + s[0];
                        s = s.Substring(1);
                        i--;
                    }
                    vet[count++] = temp; temp = "";
                    s = s.Trim(' ');
                }
            }
            for (int i = 0; i < vet.Length; i++) {
                string tTemp = "";
                for (int j = 0; j < vet[i].Length; j++) {
                    if (vet[i][j] == '.') {
                        tTemp = tTemp + ',';
                    }
                    else tTemp = tTemp + vet[i][j];
                }
                vet[i] = tTemp;
            }
            return vet;
        }
    }
}