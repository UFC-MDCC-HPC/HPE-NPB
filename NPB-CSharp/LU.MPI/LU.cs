using System;
using System.IO;
using NPB.Lub;
using NPB3_0_JAV;
using NPB3_0_JAV.BMInOut;

namespace NPB {
    public class LU:LUBase {

        public LU(char c):base(c){}

        static void Main(String[] argv) {

            LU lu = null;
            LUBase.debug = false;

            try {
                string param = argv[0];
            }
            catch (Exception) {
                argv = new String[1];
                argv[0] = "CLASS=S"; // CLASS DEFAULT, IF USER NOT TYPE CLASS=S IN COMMAND-LINE ARGUMENT
            }
            char paramClass;
            if (!LUBase.debug) {
                BMArgs.ParseCmdLineArgs(argv, BMName);
                paramClass = BMArgs.CLASS;
            }
            else {
                paramClass = 'K';  //DEBUG: CHANGE TO [K=(S and 4 PROCESSORS)] OR [S=(S and 1 PROCESSOR)]
            }                      //DEBUG: OR [T=(A and 4 PROCESSORS)] OR [I=(B and 4 PROCESSORS)]

            try {
                lu = new LU(paramClass);
            }
            catch (OutOfMemoryException e) {
                Console.WriteLine(e.ToString());
                Environment.Exit(0);
            }
            lu.runBenchMark();
        }

        public void runBenchMark() {
        }

        public static double mod(double a, double b) { return (a % b); }

        public double min(int n1, int n2) { return n1<n2?n1:n2; }

        public double max(double n1, double n2) { return n1>n2?n1:n2; }

        public double pow2(double p) { return p * p; }

    }
}


