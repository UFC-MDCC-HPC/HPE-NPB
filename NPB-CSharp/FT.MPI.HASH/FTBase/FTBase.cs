using System;
using NPB;

namespace NPB {
    public class FTBase: Base {

        //******************************************** Attributes *******************************************************/
        //public const int sizeArrayAdd = 1;
        //npbparams.h
        protected int nx, ny, nz, niter_default, maxdim, ntdivnp;//, np_min;
        //protected double ntotal_f;
        //protected string cs1, cs2, cs3, cs4, cs5, cs6, cs7;//compiletime, npbversion, 
        //end.h

        //global.h
        protected static int np1, np2, np;
        protected static int layout_type;
        protected int layout_0D = 0, layout_1D = 1, layout_2D = 2;
        protected int fftblock_default=16, fftblockpad_default=18;
        protected int transblock=32, transblockpad=34;
        protected static int fftblock, fftblockpad;
        protected static int node, me1, me2;

        protected static int[,] dims = new int[3, 3];
        protected static int[] xstart = new int[3];
        protected static int[] ystart = new int[3];
        protected static int[] zstart = new int[3];
        protected static int[] xend = new int[3];
        protected static int[] yend = new int[3];
        protected static int[] zend = new int[3];
        protected int T_total = 1, T_setup = 2, T_fft = 3, T_evolve = 4, T_checksum = 5, T_fftlow = 6, 
                          T_fftcopy = 7, T_transpose = 8, T_transxzloc = 9, T_transxzglo = 10, T_transxzfin = 11, 
                          T_transxyloc = 12, T_transxyglo = 13, T_transxyfin = 14,  T_synch = 15, T_max = 15;
        protected bool timers_enabled = false;
        protected int ilog2;
        protected double randlc;
        protected static bool debug, debugsynch, debugX=true;
        protected double seed = 314159265, a = 1220703125, pi = Math.PI, alpha=.000001;
        protected static int REAL=0,IMAG=1;

        protected static double[,] u;           //complex u(nx);
        protected static double[,] sums;        //complex sums(0:niter_default)
        protected static int niter;
        //end.h

        //MPI
        protected MPI.Environment mpi = null;
        protected static MPI.Intracommunicator worldcomm, commslice1, commslice2 = null;
        //ft.f
        protected int i, ierr;
        protected static double[,,,] u0, u1, u2;                       //complex u0(ntdivnp), u1(ntdivnp), u2(ntdivnp)
        protected static double[] twiddle;                            // twiddle(ntdivnp)
        protected static double[,] pad1 = new double[2, 3];
        protected static double[,] pad2 = new double[2, 3];
        protected static double[,] pad3 = new double[2, 3]; //complex pad1(3), pad2(3), pad3(3)
        protected char CLSS;
        //end.f

        //timer.f
        protected static double[] start = new double[64];
        protected static double[] elapsed = new double[64];
        //enf.f

        //Support
        protected int root=0;
        public static String BMName = "FT";
        //endSupport

        public Timer timer = new Timer();

        //***************************************************************************************************************/

        public FTBase(char c) {
            this.CLSS = c;
            switch(this.CLSS) {
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
