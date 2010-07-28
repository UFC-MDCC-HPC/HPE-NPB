using System;
using NPB;
using NPB3_0_JAV;

namespace NPB.Lub {
    public class LUBase {
      //******************************************** Attributes *******************************************************/
            public const int le = 1;
        //npbparams.h
            protected static int maxcells, problem_size, niter_default;
            protected static double dt_default;
            protected static int wr_default;
            protected static int iotype;
            protected static bool convertdouble = false;
            protected static string compiletime;
            protected static string npbversion = "3.3";
        //end npbparans.h

        //mpinpb.h
            protected MPI.Environment mpi = null;
            protected MPI.Intracommunicator worldcomm, comm_setup, comm_solve, comm_rhs = null;
            protected static int node, no_nodes, total_nodes, dp_type;
            protected static bool active;
        //end mpinpb.h

        //Suporte
            protected int npDebug = 0, root = 0;
            protected Timer timer = new Timer();
            protected static bool debug = false;
            protected static String BMName = "LU";
            protected char clss;
        //Suporte

      //***************************************************************************************************************/

        public LUBase(char c){
            DateTime nowTime = DateTime.Now;
            compiletime = nowTime.Day + "/" + nowTime.Month + "/" + nowTime.Year;
            this.clss = c;

            
            
            
            mpi_start();
            initVars();
        }

        private void mpi_start() {
        }

        private void initVars(){
        }

        public static string[] readInputLuData(string path) {
            String s, temp = "";
            string[] vet = new string[9];
            int[] conf = {1,1,3,2,2}; // conf numbers: first line have 1 elements, second line 1 elements, third line 3 elements... 
            int count = 0;
            System.IO.StreamReader file = new System.IO.StreamReader(path);
            for (int k = 0; k < conf.Length; k++) {
                s = file.ReadLine(); s = s.Trim(' ');
                for (int j = 0; j < conf[k]; j++) {
                    for (int i = 0; i <= s.Length; i++) {
                        if (s[0] == ' ') {
                            break;
                        }
                        temp = temp + s[0];
                        s = s.Substring(1);
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
