using System;
using NPB;

namespace NPB {
    public class FTBase:Base {

      //******************************************** Attributes *******************************************************/
        public const int sizeArrayAdd = 1;
        //npbparams.h
            protected int nx, ny, nz, maxdim, niter_default, ntdivnp, np_min;
            protected double ntotal_f;
            protected bool convertdouble;
            protected string compiletime, npbversion, cs1, cs2, cs3, cs4, cs5, cs6, cs7;
        //end.h

        //global.h
            protected static int np1, np2, np;
            protected static int layout_type;
            protected int layout_0D = 0, layout_1D = 1, layout_2D = 2;
            protected int fftblock_default=16, fftblockpad_default=18;
            protected int transblock=32, transblockpad=34;
            protected static int fftblock, fftblockpad;
            protected static int me, me1, me2;
            // Declared in block MPI: protected int commslice1, commslice2;
            protected static int[,] dims = new int[3+sizeArrayAdd, 3+sizeArrayAdd];
            protected static int[] xstart = new int[3+sizeArrayAdd];
            protected static int[] ystart = new int[3+sizeArrayAdd];
            protected static int[] zstart = new int[3+sizeArrayAdd];
            protected static int[] xend = new int[3+sizeArrayAdd];
            protected static int[] yend = new int[3+sizeArrayAdd];
            protected static int[] zend = new int[3+sizeArrayAdd];
            protected int T_total = 1, T_setup = 2, T_fft = 3, T_evolve = 4, T_checksum = 5, T_fftlow = 6, 
                          T_fftcopy = 7, T_transpose = 8, T_transxzloc = 9, T_transxzglo = 10, T_transxzfin = 11, 
                          T_transxyloc = 12, T_transxyglo = 13, T_transxyfin = 14,  T_synch = 15, T_max = 15;
            protected bool timers_enabled = false;
            protected double timer_read;
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
            //mpif.h
                protected int mpi_comm_world = 0;
                protected int mpi_max = 1, mpi_sum = 2, mpi_min = 3;
                protected int mpi_double_precision = 1, mpi_integer = 2, mpi_byte = 3, mpi_real= 4, 
                              mpi_complex = 5, mpi_double_complex = 6;
                protected int mpi_any_source = -1;
                protected int mpi_err_other = -1;
                protected double mpi_wtime;
                protected int mpi_status_size = 3;
            //end.h

            //mpinpb.h
                protected static int dc_type;
            //
        //endMpi

        //ft.f
            protected int i, ierr;
            protected static double[,,,] u0, u1, u2;                       //complex u0(ntdivnp), u1(ntdivnp), u2(ntdivnp)
            protected static double[] twiddle;                            // twiddle(ntdivnp)
            protected static double[,] pad1 = new double[3+sizeArrayAdd,2+sizeArrayAdd];
            protected static double[,] pad2 = new double[3+sizeArrayAdd,2+sizeArrayAdd];
            protected static double[,] pad3 = new double[3+sizeArrayAdd,2+sizeArrayAdd]; //complex pad1(3), pad2(3), pad3(3)
            protected int iter;
            protected double total_time, mflops;
            //protected bool verified;
            protected char clss; 
        //end.f

        //timer.f
            protected static double[] start = new double[64+sizeArrayAdd]; 
            protected static double[] elapsed = new double[64+sizeArrayAdd];
        //enf.f

        //Support
            protected int npDebug=0, root=0;
            public static String BMName = "FT";
        //endSupport

        public Timer timer = new Timer();

      //***************************************************************************************************************/

        public FTBase(char c){
            this.clss = c;
            string dayOfCompile="28 April 2010";
            switch (this.clss) {
                case ('S'):
                    nx=ny=nz=64;
                    niter_default=6;
                    maxdim=64;
                    ntotal_f=(double)(nx*ny*nz);
                    convertdouble = false;
                    compiletime = dayOfCompile;
                    npbversion="3.3";
                    break;
                case ('W'):
                    nx=ny=128;
                    nz=32;
                    niter_default=6;
                    maxdim = 128;
                    ntotal_f = (double)(nx * ny * nz);
                    convertdouble = false;
                    compiletime = dayOfCompile;
                    npbversion = "3.3";
                    break;     
                case ('A'):
                    nx=256;
                    ny=256;
                    nz=128;
                    niter_default=6;
                    maxdim = 256;
                    ntotal_f = (double)(nx * ny * nz);
                    convertdouble = false;
                    compiletime = dayOfCompile;
                    npbversion = "3.3";
                    break;      
                case ('B'):
                    nx=512;
                    ny=nz=256;
                    niter_default=20;
                    maxdim = 512;
                    ntotal_f = (double)(nx * ny * nz);
                    convertdouble = false;
                    compiletime = dayOfCompile;
                    npbversion = "3.3";
                    break;
                case ('C'):
                    nx=ny=nz=512;
                    niter_default=20;
                    maxdim = 512;
                    ntotal_f = (double)(nx * ny * nz);
                    convertdouble = false;
                    compiletime = dayOfCompile;
                    npbversion = "3.3";
                    break;
                //case ('K'):
                //    nx = ny = nz = 64;
                //    niter_default = 6;
                //    maxdim = 64;
                //    ntotal_f = (double)(nx * ny * nz);
                //    convertdouble = false;
                //    compiletime = dayOfCompile;
                //    npbversion = "3.3";
                //    npDebug = 3;
                //    this.clss = 'S';
                //    break;
                //case ('I'):
                //    nx = 512;
                //    ny = nz = 256;
                //    niter_default = 20;
                //    maxdim = 512;
                //    ntotal_f = (double)(nx * ny * nz);
                //    convertdouble = false;
                //    compiletime = dayOfCompile;
                //    npbversion = "3.3";
                //    npDebug = 3;
                //    this.clss = 'B';
                //    break;
                //case ('T'):
                //    nx = 256;
                //    ny = 256;
                //    nz = 128;
                //    niter_default = 6;
                //    maxdim = 256;
                //    ntotal_f = (double)(nx * ny * nz);
                //    convertdouble = false;
                //    compiletime = dayOfCompile;
                //    npbversion = "3.3";
                //    npDebug = 3;
                //    this.clss = 'A';
                //    break; 
            }
            mpi_start();
        }

        private void mpi_start() {
            string[] args = System.Environment.GetCommandLineArgs();
            mpi = new MPI.Environment(ref args);
            worldcomm = MPI.Communicator.world;
            np = worldcomm.Size + npDebug;
            me = worldcomm.Rank;
            np_min = np;
            ntdivnp = ((nx*ny)/np_min)*nz;
        }
    }
}
