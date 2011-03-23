using System;
using NPB;

namespace NPB {
    public class LUBase: Base {
        //******************************************** Attributes *******************************************************/
        //npbparams.h
        protected static int nnodes_compiled, isiz01, isiz02, isiz03, isiz1, isiz2, isiz3, itmax_default, inorm_default;
        protected static double dt_default;
        protected static bool convertdouble = false;
        protected static string compiletime, npbversion="3.3";
        //end npbparans.h

        //applu.incl
        protected static int ipr_default=1;
        protected static double omega_default=1.2d;
        protected static double tolrsd1_def=1.0E-08, tolrsd2_def=1.0E-08, tolrsd3_def=1.0E-08, tolrsd4_def=1.0E-08, tolrsd5_def=1.0E-08;
        protected static double c1=1.40d, c2=0.40d, c3=1.00E-01, c4=1.00d, c5=1.40d;
        //-- grid -------------------------------------------------------------
        protected static int nx, ny, nz, nx0, ny0, nz0, ipt, ist, iend, jpt, jst, jend, ii1, ii2, ji1, ji2, ki1, ki2;
        protected static double dxi, deta, dzeta, tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3;
        //---------------------------------------------------------------------
        //   dissipation ------------------------------------------------------            
        protected static double dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, dz1, dz2, dz3, dz4, dz5, dssp;
        //---------------------------------------------------------------------
        protected static double[,,,] u, rsd, frct, flux;        //       u[5, -1:isiz1+2, -1:isiz2+2, isiz3],
        //     rsd[5, -1:isiz1+2, -1:isiz2+2, isiz3],
        //    frct[5, -1:isiz1+2, -1:isiz2+2, isiz3],
        //    flux[5,  0:isiz1+1,  0:isiz2+1, isiz3];
        protected static int ipr, inorm;
        protected static int itmax, invert;
        protected static double dt, omega, frc, ttotal;
        /*tolrsd*/
        protected static double[] tolrsd = new double[5];     //tolrsd[5]
        /*rsdnm*/
        protected static double[] rsdnm  = new double[5];     //rsdnm[5]
        /*errnm*/
        protected static double[] errnm  = new double[5];     //errnm[5];
        protected static double[, , ,] a, b, c, d;              //a[5,5,isiz1,isiz2], b[5,5,isiz1,isiz2], c[5,5,isiz1,isiz2], d[5,5,isiz1,isiz2];
        /*ce*/
        //protected static double[,] ce = new double[13, 5];   //ce[5,13]
        protected static double[,] ce = 
                                        {{   2.0d,    1.0d,    2.0d,    2.0d,    5.0d},
                                         {   0.0d,    0.0d,    2.0d,    2.0d,    4.0d},
                                         {   0.0d,    0.0d,    0.0d,    0.0d,    3.0d},
                                         {   4.0d,    0.0d,    0.0d,    0.0d,    2.0d},
                                         {   5.0d,    1.0d,    0.0d,    0.0d, 1.0E-01},
                                         {   3.0d,    2.0d,    2.0d,    2.0d, 4.0E-01},
                                         {5.0E-01,    3.0d,    3.0d,    3.0d, 3.0E-01},
                                         {2.0E-02, 1.0E-02, 4.0E-02, 3.0E-02, 5.0E-02},
                                         {1.0E-02, 3.0E-02, 3.0E-02, 5.0E-02, 4.0E-02},
                                         {3.0E-02, 2.0E-02, 5.0E-02, 4.0E-02, 3.0E-02},
                                         {5.0E-01, 4.0E-01, 3.0E-01, 2.0E-01, 1.0E-01},
                                         {4.0E-01, 3.0E-01, 5.0E-01, 1.0E-01, 3.0E-01},
                                         {3.0E-01, 5.0E-01, 4.0E-01, 3.0E-01, 2.0E-01}};

        protected static int node, ndim, num, xdim, ydim, row, col;
        protected static int north, south, east, west;
        protected static int from_s = 1, from_n = 2, from_e = 3, from_w = 4;
        protected static int npmax;

        protected static bool[] icommn, icomms, icomme, icommw; //Fortran: icommn[npmax+1],icomms[npmax+1], icomme[npmax+1],icommw[npmax+1]
        //protected static double[,] buf, buf1;                   //Fortran: buf[5,2*isiz2*isiz3], buf1[5,2*isiz2*isiz3]
        protected static double maxtime;
        //end applu.incl

        //lu.f
        //protected static bool verified;
        protected static double mflops;
        //end lu.f

        //mpinpb.h
        protected static int no_nodes;// dp_type, node
        protected MPI.Environment mpi = null;
        protected MPI.Intracommunicator worldcomm;// comm_setup, comm_solve, comm_rhs = null;
        //end mpinpb.h

        //Suporte
        protected int root = 0;
        protected Timer timer = new Timer();
        protected static bool debug = false;
        protected static String BMName = "LU";
        protected char clss;
        //Suporte

        //***************************************************************************************************************/

        public LUBase(char c) {
            DateTime nowTime = DateTime.Now;
            compiletime = nowTime.Day + "/" + nowTime.Month + "/" + nowTime.Year;
            this.clss = c;
            switch(c) {
                case 'S':
                    isiz01 = isiz02 = isiz03 =12;
                    isiz3  = isiz03;
                    itmax_default = inorm_default=50;
                    dt_default = 0.5d;
                    break;
                case 'W':
                    isiz01 = isiz02 = isiz03 = 33;
                    isiz3  = isiz03;
                    itmax_default = inorm_default = 300;
                    dt_default = 1.5E-3; // dt_default = .0015;
                    break;
                case 'A':
                    isiz01 = isiz02 = isiz03 = 64;
                    isiz3  = isiz03;
                    itmax_default = inorm_default = 250;
                    dt_default = 2.0d;
                    break;
                case 'B':
                    isiz01 = isiz02 = isiz03 = 102;
                    isiz3  = isiz03;
                    itmax_default = inorm_default = 250;
                    dt_default = 2.0d;
                    break;
                case 'C':
                    isiz01 = isiz02 = isiz03 = 162;
                    isiz3  = isiz03;
                    itmax_default = inorm_default = 250;
                    dt_default = 2.0d;
                    break;
                case 'D':
                    isiz01 = isiz02 = isiz03 = 408;
                    isiz3  = isiz03;
                    itmax_default = inorm_default = 300;
                    dt_default = 1.0d;
                    break;
                case 'E':
                    isiz01 = isiz02 = isiz03 = 1020;
                    isiz3  = isiz03;
                    itmax_default = inorm_default = 300;
                    dt_default = 0.5d;
                    break;
            }
            mpi_start();
            initVars();
        }

        private void mpi_start() {
            string[] args = System.Environment.GetCommandLineArgs();
            mpi = new MPI.Environment(ref args);
            worldcomm = MPI.Communicator.world;//call MPI_INIT(IERROR)
            num = worldcomm.Size;//call MPI_COMM_SIZE(MPI_COMM_WORLD, num, IERROR)
            node = worldcomm.Rank;//call MPI_COMM_RANK(MPI_COMM_WORLD, id, IERROR)            
            ndim   = nodedim(num);
        }

        private void initVars() {
            npmax = isiz01 + isiz02;
            nnodes_compiled = num;
            int ydiv = ilog2(num) / 2;
            int xdiv = ydiv;
            if(xdiv + ydiv != ilog2(num))
                xdiv += 1;
            xdiv = ipow2(xdiv);
            ydiv = ipow2(ydiv);
            isiz1 = isiz01 / xdiv;
            if(isiz1 * xdiv < isiz01)
                isiz1++;
            isiz2 = isiz01 / ydiv;
            if(isiz2 * ydiv < isiz01)
                isiz2++;

            /**/
            icommn = new bool[npmax+1];//icommn[npmax+1]
            /**/
            icomms = new bool[npmax+1];//icomms[npmax+1]
            /**/
            icomme = new bool[npmax+1];//icomme[npmax+1]
            /**/
            icommw = new bool[npmax+1];//icommw[npmax+1]

            a = new double[isiz2, isiz1, 5, 5];//a[5,5,isiz1,isiz2];
            b = new double[isiz2, isiz1, 5, 5];//b[5,5,isiz1,isiz2];
            c = new double[isiz2, isiz1, 5, 5];//c[5,5,isiz1,isiz2];
            d = new double[isiz2, isiz1, 5, 5];//d[5,5,isiz1,isiz2];

            u    = new double[isiz3, isiz2+4, isiz1+4, 5];//       u[5, -1:isiz1+2, -1:isiz2+2, isiz3];
            rsd  = new double[isiz3, isiz2+4, isiz1+4, 5];//     rsd[5, -1:isiz1+2, -1:isiz2+2, isiz3];
            frct = new double[isiz3, isiz2+4, isiz1+4, 5];//    frct[5, -1:isiz1+2, -1:isiz2+2, isiz3];
            flux = new double[isiz3, isiz2+2, isiz1+2, 5];//    flux[5,  0:isiz1+1,  0:isiz2+1, isiz3];
        }
        public static int nodedim(double n) {
            return (int)(Math.Log(n) / Math.Log(2.0d) + 0.00001);
        }
        public static int ilog2(int i) {
            int log2;
            int exp2 = 1;
            if(i <= 0)
                return (-1);
            for(log2 = 0; log2 < 20; log2++) {
                if(exp2 == i)
                    return (log2);
                exp2 *= 2;
            }
            return (-1);
        }
        public static int ipow2(int i) {
            int pow2 = 1;
            if(i < 0)
                return (-1);
            if(i == 0)
                return (1);
            while(i-->0)
                pow2 *= 2;
            return (pow2);
        }
    }
}
