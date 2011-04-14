using System;
using NPB;

namespace NPB {
    public class FTBase:Base {

      //******************************************** Attributes *******************************************************/
        //npbparams.h
            protected int nx, ny, nz, maxdim, niter_default, ntdivnp;
            protected bool convertdouble = false;
        //end.h

        //global.h
            protected static int np1, np2, np;
            protected static int layout_type;
        //constant
            protected int layout_0D = 0, layout_1D = 1, layout_2D = 2;
            protected int fftblock_default=16, fftblockpad_default=18;
            protected int transblock=32, transblockpad=34;
            protected int T_total = 1, T_setup = 2, T_fft = 3, T_evolve = 4, T_checksum = 5, T_fftlow = 6, 
                          T_fftcopy = 7, T_transpose = 8, T_transxzloc = 9, T_transxzglo = 10, T_transxzfin = 11, 
                          T_transxyloc = 12, T_transxyglo = 13, T_transxyfin = 14,  T_synch = 15, T_max = 15;
            protected double seed = 314159265, a = 1220703125, pi = Math.PI, alpha=.000001;
            protected static int REAL=0,IMAG=1;
        //end constant
            protected static int fftblock, fftblockpad;
            protected static int node, me1, me2;
            protected static int[,] dims = new int[3+1, 3+1];
            protected static int[] xstart = new int[3+1];
            protected static int[] ystart = new int[3+1];
            protected static int[] zstart = new int[3+1];
            protected static int[] xend = new int[3+1];
            protected static int[] yend = new int[3+1];
            protected static int[] zend = new int[3+1];
            protected bool timers_enabled = false;

            protected static double[,] u;           //complex u(nx);
            protected static int niter;
        //end.h

        //MPI
            protected MPI.Environment mpi = null;
            protected static MPI.Intracommunicator worldcomm, commslice1, commslice2 = null;
        //endMpi

        //ft.f
            protected static double[,,,] u0, u1, u2;                       //complex u0(ntdivnp), u1(ntdivnp), u2(ntdivnp)
            protected static double[] twiddle;                            // twiddle(ntdivnp)
            protected char CLSS; 
        //end.f

        //Support
            protected int root=0;
            public static String BMName = "FT";
        //endSupport

        public Timer timer = new Timer();

      //***************************************************************************************************************/

        public FTBase(char c){
            this.CLSS = c;
            switch (this.CLSS) {
                case ('S'):
                    nx=ny=nz=64;
                    niter_default=6;
                    maxdim=64;
                    break;
                case ('W'):
                    nx=ny=128;
                    nz=32;
                    niter_default=6;
                    maxdim = 128;
                    break;     
                case ('A'):
                    nx=256;
                    ny=256;
                    nz=128;
                    niter_default=6;
                    maxdim = 256;
                    break;      
                case ('B'):
                    nx=512;
                    ny=nz=256;
                    niter_default=20;
                    maxdim = 512;
                    break;
                case ('C'):
                    nx=ny=nz=512;
                    niter_default=20;
                    maxdim = 512;
                    break;
            }
            mpi_start();
        }

        private void mpi_start() {
            string[] args = System.Environment.GetCommandLineArgs();
            mpi = new MPI.Environment(ref args);
            worldcomm = MPI.Communicator.world;
            np = worldcomm.Size;
            node = worldcomm.Rank;
            ntdivnp = ((nx*ny)/np)*nz;
        }
    }
}
