using System;
using NPB;
using NPB3_0_JAV;

namespace NPB.Btb {
    public class BTBase {
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

        //header.h
            protected int aa = 1, bb = 2, cc = 3, BLOCK_SIZE = 5;
            protected static int ncells;
            protected static int[] grid_points = new int[3+le];
            protected static double elapsed_time;
            protected static double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4, dy5, 
                                 dz1, dz2, dz3, dz4, dz5, dssp, dt, dxmax, dymax, dzmax, xxcon1, xxcon2, xxcon3, xxcon4, xxcon5, 
                                 dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4, yycon5, dy1ty1, dy2ty1, 
                                 dy3ty1, dy4ty1, dy5ty1, zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1, dz2tz1, dz3tz1, dz4tz1, 
                                 dz5tz1, dnxm1, dnym1, dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2, c3, c4, c5, c4dssp, c5dssp, 
                                 dtdssp, dttx1, bt, dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1, c2dtty1, c2dttz1, comz1, comz4, 
                                 comz5, comz6, c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;
            protected static double[,] ce = new double[5+le, 13+le];
            protected int EAST = 2000, WEST = 3000, NORTH = 4000, SOUTH = 5000, BOTTOM = 6000, TOP = 7000;
            protected static int[,] cell_coord;
            protected static int[,] cell_low;
            protected static int[,] cell_high;
            protected static int[,] cell_size;
            protected static int[] predecessor;
            protected static int[,] slice;
            protected static int[] grid_size;
            protected static int[] successor;
            protected static int[,] start;
            protected static int[,] end;
            protected int IMAX, JMAX, KMAX, MAX_CELL_DIM, BUF_SIZE;
            protected static double[, , ,] us, vs, ws, qs, rho_i, square, backsub_info;
            protected static double[, , , ,] forcing, u, rhs;
            protected static double[, , , , ,] lhsc;
            //protected static double[] in_buffer, out_buffer;
            protected static double[] sum;

            protected static double[] cv, rhon, rhos, rhoq, cuf, q;
            protected static double[,] ue, buf; //

            protected static int west_size, east_size, bottom_size, top_size, north_size, south_size, start_send_west, start_send_east,
                              start_send_south, start_send_north, start_send_bottom, start_send_top, start_recv_west, start_recv_east,
                              start_recv_south, start_recv_north, start_recv_bottom, start_recv_top;
            protected static double[,] tmp_block = new double[5, 5];
            protected static double[,] b_inverse = new double[5, 5];
            protected static double[] tmp_vec = new double[5];
            protected static double[] xce_sub = new double[5];
            protected static int collbuf_nodes, collbuf_size, iosize, eltext, combined_btype, fp, idump, record_length, idump_sub, rd_interval;
            protected static int iseek, element, combined_ftype;
        //end header.h

        //mpinpb.h
            protected MPI.Environment mpi = null;
            protected MPI.Intracommunicator worldcomm, comm_setup, comm_solve, comm_rhs = null;
            protected static int node, no_nodes, total_nodes, dp_type;
            protected static bool active;
        //end mpinpb.h

        // work_lhs.h
            protected double tmp1, tmp2, tmp3;
            protected double[,,] fjac; //[5, 5, -2:MAX_CELL_DIM+1];
            protected double[,,] njac; //[5, 5, -2:MAX_CELL_DIM+1];
            protected double[,,] lhsa; //[5, 5, -1:MAX_CELL_DIM]; 
            protected double[,,] lhsb; //[5, 5, -1:MAX_CELL_DIM];
        // end work_lhs.h

        //bt.f
            protected int i, niter, step, c, error, fstatus;
            protected double navg, mflops, mbytes, n3;
            protected double t, tmax, tiominv, tpc;
            protected bool verified;
            protected string cbuff;
            protected int wr_interval;
        //end bt.f

        //Suporte
            protected int npDebug = 0, root = 0;
            protected Timer timer = new Timer();
            protected static bool debug = false;
            protected static String BMName = "BT";
            protected char clss;
        //Suporte

      //***************************************************************************************************************/

        public BTBase(char c){
            DateTime nowTime = DateTime.Now;
            compiletime = nowTime.Day + "/" + nowTime.Month + "/" + nowTime.Year;
            this.clss = c;
            switch (this.clss) {
                case ('S'):
                    problem_size=12;
                    niter_default=60;
                    dt_default = 0.010;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    break;
                case ('W'):
                    problem_size=24;
                    niter_default=200;
                    dt_default = 0.0008;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    break;     
                case ('A'):
                    problem_size=64;
                    niter_default=200;
                    dt_default = 0.0008;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    break;      
                case ('B'):
                    problem_size=102;
                    niter_default=200;
                    dt_default = 0.0003;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    break;
                case ('C'):
                    problem_size=162;
                    niter_default=200;
                    dt_default = 0.0001;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    break;
                case ('K'):
                    problem_size = 12;
                    niter_default = 60;
                    dt_default = 0.010;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    npDebug = 3;
                    this.clss = 'S';
                    break;
                case ('T'):
                    problem_size = 64;
                    niter_default = 200;
                    dt_default = 0.0008;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    npDebug = 3;
                    this.clss = 'A';
                    break; 
                case ('I'):
                    problem_size = 102;
                    niter_default = 200;
                    dt_default = 0.0003;
                    wr_default = 5;
                    iotype = 0;
                    convertdouble = false;
                    npDebug = 3;
                    this.clss = 'B';
                    break;
            }
            mpi_start();
            initVars();
        }

        private void mpi_start() {
            int color, nc;

            string[] args = System.Environment.GetCommandLineArgs();
            mpi = new MPI.Environment(ref args);
            worldcomm = MPI.Communicator.world;
            total_nodes = worldcomm.Size + npDebug;
            node = worldcomm.Rank;

            nc = Convert.ToInt32(Math.Sqrt(total_nodes) + 0.00001d);
            maxcells = Convert.ToInt32(Math.Sqrt(total_nodes));

            if (nc > maxcells) nc = maxcells;
            if (node >= nc*nc) {
               active = false;
               color = 1;
            } else {
               active = true;
               color = 0;
            }

            comm_setup = (MPI.Intracommunicator)worldcomm.Split(color, node);//call mpi_comm_split[MPI_COMM_WORLD,color,node,comm_setup,error];

            if (!active) return;

            no_nodes = comm_setup.Size + npDebug;//call mpi_comm_size[comm_setup, no_nodes, error];
            comm_solve = (MPI.Intracommunicator)comm_setup.Clone();//call mpi_comm_dup[comm_setup, comm_solve, error];
            comm_rhs = (MPI.Intracommunicator)comm_setup.Clone();//call mpi_comm_dup[comm_setup, comm_rhs, error];

            root = 0;
        }

        private void initVars(){
            cell_coord   = new int[3+le,maxcells+le];
            cell_low     = new int[3+le,maxcells+le];
            cell_high    = new int[3+le,maxcells+le];
            cell_size    = new int[3+le,maxcells+le];
            predecessor  = new int[3+le];
            slice        = new int[3+le,maxcells+le];
            grid_size    = new int[3+le];
            successor    = new int[3+le];
            start        = new int[3+le,maxcells+le];
            end          = new int[3+le,maxcells+le];
            MAX_CELL_DIM = (problem_size/maxcells)+1;
            IMAX         = MAX_CELL_DIM;
            JMAX         = MAX_CELL_DIM;
            KMAX         = MAX_CELL_DIM;
            BUF_SIZE     = MAX_CELL_DIM*MAX_CELL_DIM*(maxcells-1)*60+1;
            sum          = new double[niter_default];

            us           = new double[maxcells+le,KMAX+2,JMAX+2,IMAX+2];            //us     [    -1:IMAX,  -1:JMAX,  -1:KMAX,   maxcells]
            vs           = new double[maxcells+le,KMAX+2,JMAX+2,IMAX+2];            //vs     [    -1:IMAX,  -1:JMAX,  -1:KMAX,   maxcells]
            ws           = new double[maxcells+le,KMAX+2,JMAX+2,IMAX+2];            //ws     [    -1:IMAX,  -1:JMAX,  -1:KMAX,   maxcells]
            qs           = new double[maxcells+le,KMAX+2,JMAX+2,IMAX+2];            //qs     [    -1:IMAX,  -1:JMAX,  -1:KMAX,   maxcells]
            rho_i        = new double[maxcells+le,KMAX+2,JMAX+2,IMAX+2];            //rho_i  [    -1:IMAX,  -1:JMAX,  -1:KMAX,   maxcells]
            square       = new double[maxcells+le,KMAX+2,JMAX+2,IMAX+2];            //square [    -1:IMAX,  -1:JMAX,  -1:KMAX,   maxcells]
            
            forcing      = new double[maxcells+le,KMAX  ,JMAX  ,IMAX  ,5+le];       //forcing[5, 0:IMAX-1, 0:JMAX-1, 0:KMAX-1, maxcells]
            u            = new double[maxcells+le,KMAX+4,JMAX+4,IMAX+4,5+le];       //u[5,-2:IMAX+1,-2:JMAX+1,-2:KMAX+1, maxcells]
            rhs          = new double[maxcells+le,KMAX+1,JMAX+1,IMAX+1,5+le];       //rhs    [5,  -1:IMAX-1,-1:JMAX-1,-1:KMAX-1, maxcells]
            lhsc         = new double[maxcells+le,(KMAX-1)+2,(JMAX-1)+2,(IMAX-1)+2,5+le,5+le];//lhsc[5,5,-1:IMAX-1,-1:JMAX-1,-1:KMAX-1, maxcells]
            backsub_info = new double[maxcells+le,MAX_CELL_DIM+1,MAX_CELL_DIM+1,5+le];//backsub_info[5, 0:MAX_CELL_DIM, 0:MAX_CELL_DIM, maxcells]
            //in_buffer    = new double[BUF_SIZE+le];                                    //in_buffer[BUF_SIZE]
            //out_buffer   = new double[BUF_SIZE+le];                                    //out_buffer[BUF_SIZE]

            cv           = new double[MAX_CELL_DIM + 4];//[-2:MAX_CELL_DIM+1];
            rhon         = new double[MAX_CELL_DIM + 4];//[-2:MAX_CELL_DIM+1];
            rhos         = new double[MAX_CELL_DIM + 4];//[-2:MAX_CELL_DIM+1];
            rhoq         = new double[MAX_CELL_DIM + 4];//[-2:MAX_CELL_DIM+1];
            cuf          = new double[MAX_CELL_DIM + 4];//[-2:MAX_CELL_DIM+1];
            q            = new double[MAX_CELL_DIM + 4];//[-2:MAX_CELL_DIM+1];
            ue           = new double[MAX_CELL_DIM+4,5+le];//[-2:MAX_CELL_DIM+1,5];
            buf          = new double[MAX_CELL_DIM+4,5+le];//[-2:MAX_CELL_DIM+1,5];

            fjac = new double[MAX_CELL_DIM+4, 5+le, 5+le]; //[5, 5, -2:MAX_CELL_DIM+1];
            njac = new double[MAX_CELL_DIM+4, 5+le, 5+le]; //[5, 5, -2:MAX_CELL_DIM+1];
            lhsa = new double[MAX_CELL_DIM+2, 5+le, 5+le]; //[5, 5, -1:MAX_CELL_DIM]; 
            lhsb = new double[MAX_CELL_DIM+2, 5+le, 5+le]; //[5, 5, -1:MAX_CELL_DIM];            
        }

        public static string[] readInputBtData(string path) {
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