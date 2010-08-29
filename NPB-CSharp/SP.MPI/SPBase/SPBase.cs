/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                               S P B a s e                               !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    SPbase implements base class for SP benchmark.                       !
!									  !
!    Permission to use, copy, distribute and modify this software	  !
!    for any purpose with or without fee is hereby granted.  We 	  !
!    request, however, that all derived work reference the NAS  	  !
!    Parallel Benchmarks 3.0. This software is provided "as is" 	  !
!    without express or implied warranty.				  !
!									  !
!    Information on NPB 3.0, including the Technical Report NAS-02-008	  !
!    "Implementation of the NAS Parallel Benchmarks in Java",		  !
!    original specifications, source code, results and information	  !
!    on how to submit new results, is available at:			  !
!									  !
!	    http://www.nas.nasa.gov/Software/NPB/			  !
!									  !
!    Send comments or suggestions to  npb@nas.nasa.gov  		  !
!									  !
!	   NAS Parallel Benchmarks Group				  !
!	   NASA Ames Research Center					  !
!	   Mail Stop: T27A-1						  !
!	   Moffett Field, CA   94035-1000				  !
!									  !
!	   E-mail:  npb@nas.nasa.gov					  !
!	   Fax:     (650) 604-3957					  !
!									  !
!-------------------------------------------------------------------------!
!     Translation to Java and to MultiThreaded Code:			  !
!     Michael A. Frumkin					          !
!     Mathew Schultz	   					          !
!-------------------------------------------------------------------------!
*/
using System;
//using System.Threading;
using NPB;
using MPI;

namespace NPB3_0_JAV.SPThreads {

public class SPBase /* : Thread*/
{
  
	public static String BMName = "SP";
	public char CLASS = 'S';

    protected int node, no_nodes, total_nodes, root, dp_type;
    protected bool active;
    protected const int DEFAULT_TAG = 0;

    protected MPI.Environment mpi = null;
    protected Intracommunicator worldcomm, comm_setup, comm_solve, comm_rhs = null;


	protected int IMAX = 0, JMAX = 0, KMAX = 0, MAX_CELL_DIM = 0, BUF_SIZE = 0, IMAXP, JMAXP,
				  problem_size = 0, maxcells = 0, ncells = 0;
	protected int[] grid_points = { 0, 0, 0 };
	protected int niter_default = 0;
	protected double dt_default = 0.0d;

	protected double[,,,,] u, rhs, forcing, lhs;

	protected double[,,,] us, vs, ws, qs, ainv, rho_i, speed, square;

	protected double[,] ue, buf;

	protected double[] cv, rhon, rhos, rhoq, cuf, q, in_buffer, out_buffer;

	protected static double tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3,
			                dx1, dx2, dx3, dx4, dx5, dy1, dy2, dy3, dy4,
			                dy5, dz1, dz2, dz3, dz4, dz5, dssp, dt,
			                dxmax, dymax, dzmax, xxcon1, xxcon2,
			                xxcon3, xxcon4, xxcon5, dx1tx1, dx2tx1, dx3tx1,
			                dx4tx1, dx5tx1, yycon1, yycon2, yycon3, yycon4,
			                yycon5, dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
			                zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, dz1tz1,
			                dz2tz1, dz3tz1, dz4tz1, dz5tz1, dnxm1, dnym1,
			                dnzm1, c1c2, c1c5, c3c4, c1345, conz1, c1, c2,
			                c3, c4, c5, c4dssp, c5dssp, dtdssp, dttx1, bt,
			                dttx2, dtty1, dtty2, dttz1, dttz2, c2dttx1,
			                c2dtty1, c2dttz1, comz1, comz4, comz5, comz6,
			                c3c4tx3, c3c4ty3, c3c4tz3, c2iv, con43, con16;

    protected const int EAST = 2000, WEST = 3000, NORTH = 4000, SOUTH = 5000, BOTTOM = 6000, TOP = 7000;

    protected int[,] cell_coord, cell_low, cell_high, cell_size, slice, start, end;
    protected int[] predecessor, grid_size, successor;


    protected static int west_size, east_size, bottom_size, top_size,
                    north_size, south_size, start_send_west,
                    start_send_east, start_send_south, start_send_north,
                    start_send_bottom, start_send_top, start_recv_west,
                    start_recv_east, start_recv_south, start_recv_north,
                    start_recv_bottom, start_recv_top;


    protected double[,] ce ={{2.0d,1.0d,2.0d,2.0d,5.0d},
		                     {0.0d,0.0d,2.0d,2.0d,4.0d},
		                     {0.0d,0.0d,0.0d,0.0d,3.0d},
		                     {4.0d,0.0d,0.0d,0.0d,2.0d},
		                     {5.0d,1.0d,0.0d,0.0d,0.1d},
		                     {3.0d,2.0d,2.0d,2.0d,0.4d},
		                     {0.5d,3.0d,3.0d,3.0d,0.3d},
		                     {0.02d,0.01d,0.04d,0.03d,0.05d},
		                     {0.01d,0.03d,0.03d,0.05d,0.04d},
		                     {0.03d,0.02d,0.05d,0.04d,0.03d},
		                     {0.5d,0.4d,0.3,0.2d,0.1d},
		                     {0.4d,0.3d,0.5,0.1d,0.3d},
		                     {0.3d,0.5d,0.4,0.3d,0.2d}};

	public bool timeron = false;
	public Timer timer = new Timer();
/*	public static int t_total = 1, t_rhsx = 2,
					 t_rhsy = 3, t_rhsz = 4, t_rhs = 5,
					 t_xsolve = 6, t_ysolve = 7, t_zsolve = 8, t_rdis1 = 9,
					 t_rdis2 = 10, t_txinvr = 11, t_pinvr = 12, t_ninvr = 13,
					 t_tzetar = 14, t_add = 15, t_last = 15; */

	public SPBase() { }

	public SPBase(char clss)
	{
        setup_mpi();

		CLASS = clss;
        
		switch (clss)
		{
			case 'S':
                problem_size = 12;
				dt_default = .015d;
				niter_default = 100;
				break;
			case 'W':
				problem_size = 36;
				dt_default = .0015d;
				niter_default = 400;
				break;
			case 'A':
				problem_size = 64;
				dt_default = .0015d;
				niter_default = 400;
				break;
			case 'B':
				problem_size = 102;
				dt_default = .001d;
				niter_default = 400;
				break;
			case 'C':
				problem_size = 162;
				dt_default = .00067d;
				niter_default = 400;
				break;
            default :
                Console.WriteLine(node + ": NO PROBLEM INSTANCE !!!");
                break;
		}

        MAX_CELL_DIM = (problem_size/maxcells)+1;                
		IMAX = JMAX = KMAX = grid_points[0] = grid_points[1] = grid_points[2] = MAX_CELL_DIM;
        IMAXP = IMAX / 2 * 2 + 1;
        JMAXP = JMAX / 2 * 2 + 1;
        BUF_SIZE = MAX_CELL_DIM*MAX_CELL_DIM*(maxcells-1)*60*2+1;
        u = new double[maxcells, 5, KMAX + 1 + 3, JMAXP + 1 + 3, IMAXP + 1 + 3];   /* x:-2 y:-2 z:-2 */

        rhs = new double[maxcells, 5, KMAX + 2, JMAXP + 2, IMAXP + 2];                     /* x:0  y:0  z:0 */
        forcing = new double[maxcells, 5, KMAX + 2, JMAXP + 2, IMAXP + 2];                 /* x:0  y:0  z:0 */
        lhs = new double[maxcells, 15, KMAX + 2, JMAXP + 2, IMAXP + 2];                    /* x:0  y:0  z:0 */

        us = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];               /* x:-1  y:-1  z:-1 */
        vs = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];               /* x:-1  y:-1  z:-1 */
        ws = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];               /* x:-1  y:-1  z:-1 */
        qs = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];               /* x:-1  y:-1  z:-1 */
        ainv = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];               /* x:-1  y:-1  z:-1 */
        rho_i = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];            /* x:-1  y:-1  z:-1 */
        speed = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];            /* x:-1  y:-1  z:-1 */
        square = new double[maxcells, KMAX + 3, JMAX + 3, IMAX + 3];           /* x:-1  y:-1  z:-1 */

		ue = new double[MAX_CELL_DIM + 4, 5];   /* -2 */
		buf = new double[MAX_CELL_DIM + 4, 5];  /* -2 */

		cv = new double[MAX_CELL_DIM + 4];     /* -2 */
		rhon = new double[MAX_CELL_DIM + 4];   /* -2 */
		rhos = new double[MAX_CELL_DIM + 4];   /* -2 */
		rhoq = new double[MAX_CELL_DIM + 4];   /* -2 */
		cuf = new double[MAX_CELL_DIM + 4];    /* -2 */
		q = new double[MAX_CELL_DIM + 4];      /* -2 */

        cell_coord = new int[maxcells,3];
        cell_high = new int[maxcells,3];
        cell_low = new int[maxcells,3];
        cell_size = new int[maxcells,3];
        slice = new int[maxcells,3];
        start = new int[maxcells,3];
        end = new int[maxcells,3];

        in_buffer = new double[BUF_SIZE];


        out_buffer = new double[BUF_SIZE];

        predecessor = new int[3];
        grid_size = new int[3];
        successor = new int[3];

	}

	//protected Thread master = null;

    private void setup_mpi()
    {
        int error, nc, color;

        string[] args = System.Environment.GetCommandLineArgs();
        mpi = new MPI.Environment(ref args);
        
        worldcomm = Communicator.world;

        total_nodes = worldcomm.Size;
        node = worldcomm.Rank;


        //---------------------------------------------------------------------
        //     compute square root; add small number to allow for roundoff
        //---------------------------------------------------------------------
        nc = Convert.ToInt32(Math.Sqrt(total_nodes) + 0.00001d);
        maxcells = Convert.ToInt32(Math.Sqrt(total_nodes));

        //---------------------------------------------------------------------
        // We handle a non-square number of nodes by making the excess nodes
        // inactive. However, we can never handle more cells than were compiled
        // in. 
        //---------------------------------------------------------------------

        if (nc > maxcells) nc = maxcells;

        if (node >= nc * nc)
        {
            active = false;
            color = 1;
        }
        else
        {
            active = true;
            color = 0;
        }

        comm_setup = (Intracommunicator)worldcomm.Split(color, node);

        if (!active) return;

        no_nodes = comm_setup.Size;
        comm_solve = (Intracommunicator)comm_setup.Clone();
        comm_rhs = (Intracommunicator)comm_setup.Clone();


        //---------------------------------------------------------------------
        //     let node 0 be the root for the group (there is only one)
        //---------------------------------------------------------------------
        root = 0;
    }

	public double dmax1(double a, double b)
	{
		if (a < b) return b; else return a;
	}

	public double dmax1(double a, double b, double c, double d)
	{
		return dmax1(dmax1(a, b), dmax1(c, d));
	}


	public void exact_solution(double xi, double eta, double zeta, double[] dtemp, int offset)
	{
		for (int m = 0; m < 5; m++)
		{
			dtemp[m + offset] = ce[0,m] +
			xi * (ce[1,m] + xi * (ce[4,m] + xi * (ce[7,m] + xi * ce[10,m]))) +
			eta * (ce[2,m] + eta * (ce[5,m] + eta * (ce[8,m] + eta * ce[11,m]))) +
			zeta * (ce[3,m] + zeta * (ce[6,m] + zeta * (ce[9,m] +
			zeta * ce[12,m])));
		}
	}

	public void lhsinit()
	{

        //---------------------------------------------------------------------
        // loop over all cells                                       
        //---------------------------------------------------------------------
        for (int c = 0; c < ncells; c++)
        {
            //---------------------------------------------------------------------
            // first, initialize the start and end arrays
            //---------------------------------------------------------------------

            for (int d = 0; d < 3; d++)
            {
                if (cell_coord[c, d] == 0)
                    start[c, d] = 3;
                else
                    start[c, d] = 2;

                if (cell_coord[c, d] == ncells - 1)
                    end[c, d] = 1;
                else
                    end[c, d] = 0;
            }


            //---------------------------------------------------------------------
            // zap the whole left hand side for starters
            //---------------------------------------------------------------------

            for (int n = 0; n < 15; n++)
            {
                for (int k = 2; k < 2 + cell_size[c, 2]; k++)
                {
                    for (int j = 2; j < 2 + cell_size[c, 1]; j++)
                    {
                        for (int i = 2; i < 2 + cell_size[c, 0]; i++)
                        {
                            lhs[c, n, k, j, i] = 0.0d;
                        }
                    }
                }
            }

            //---------------------------------------------------------------------
            // next, set all diagonal values to 1. This is overkill, but convenient
            //---------------------------------------------------------------------

            for (int n = 0; n < 3; n++)
            {
                for (int k = 2; k < 2 + cell_size[c, 2]; k++)
                {
                    for (int j = 2; j < 2 + cell_size[c, 1]; j++)
                    {
                        for (int i = 2; i < 2 + cell_size[c, 0]; i++)
                        {
                            lhs[c, 5 * n + 2, k, j, i] = 1.0d;
                        }
                    }
                }
            }
        }
   	}

	public void initialize()
	{
		int c, i, j, k, m, ii, jj, kk, ix, iy, iz;
		double xi, eta, zeta, Pxi, Peta, Pzeta; 
        double[] Pface = new double[5*3*2];
		double[] temp = new double[5];

		//---------------------------------------------------------------------
		//  Later (in compute_rhs) we compute 1/u for every element. A few of 
		//  the corner elements are not used, but it convenient (and faster) 
		//  to compute the whole thing with a simple loop. Make sure those 
		//  values are nonzero by initializing the whole thing here. 
		//---------------------------------------------------------------------
        for (c = 0; c < ncells; c++)
        {

            for (k = 1; k <= IMAX + 2; k++)
            {
                for (j = 1; j <= IMAX + 2; j++)
                {
                    for (i = 1; i <= IMAX + 2; i++)
                    {
                        u[c, 0, k, j, i] = 1.0d;
                        u[c, 1, k, j, i] = 0.0d;
                        u[c, 2, k, j, i] = 0.0d;
                        u[c, 3, k, j, i] = 0.0d;
                        u[c, 4, k, j, i] = 1.0d;
                    }
                }
            }
        }


		//---------------------------------------------------------------------
		// first store the "interpolated" values everywhere on the grid    
		//---------------------------------------------------------------------
        for (c = 0; c < maxcells; c++)
        {
            kk = 2;
            for (k = cell_low[c,2]; k <= cell_high[c,2]; k++)
            {
                zeta = k * dnzm1;
                jj = 2;
                for (j = cell_low[c,1]; j <= cell_high[c,1]; j++)
                {
                    eta = j * dnym1;
                    ii = 2;
                    for (i = cell_low[c,0]; i <= cell_high[c,0]; i++)
                    {
                        xi = i * dnxm1;

                        for (ix = 0; ix <= 1; ix++)
                        {
                            exact_solution(ix, eta, zeta, Pface, 0 + 0 * 5 + ix * 15);
                        }
                        for (iy = 0; iy <= 1; iy++)
                        {
                            exact_solution(xi, iy, zeta, Pface, 0 + 1 * 5 + iy * 15);
                        }

                        for (iz = 0; iz <= 1; iz++)
                        {
                            exact_solution(xi, eta, iz, Pface, 0 + 2 * 5 + iz * 15);
                        }

                        for (m = 0; m < 5; m++)
                        {
                            Pxi = xi * Pface[m + 0 * 5 + 1 * 15] +
                                (1.0d - xi) * Pface[m + 0 * 5 + 0 * 15];
                            Peta = eta * Pface[m + 1 * 5 + 1 * 15] +
                                    (1.0d - eta) * Pface[m + 1 * 5 + 0 * 15];
                            Pzeta = zeta * Pface[m + 2 * 5 + 1 * 15] +
                                    (1.0d - zeta) * Pface[m + 2 * 5 + 0 * 15];

                            u[c, m, kk, jj, ii] =
                              Pxi + Peta + Pzeta -
                                      Pxi * Peta - Pxi * Pzeta - Peta * Pzeta +
                                      Pxi * Peta * Pzeta;
//                            if (node==0) Console.WriteLine(u[c,m,kk,jj,ii]); 
                        }
                        ii++;
                    }
                    jj++;
                }
                kk++;
            }
        }
        

		//---------------------------------------------------------------------
		// now store the exact values on the boundaries        
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		// west face                                                  
		//---------------------------------------------------------------------

        c = slice[0, 0];
        ii = 2;
		xi = 0.0d;
		kk = 2;
		for (k = cell_low[c,2]; k <= cell_high[c,2]; k++)
		{
			zeta = k * dnzm1;
            jj = 2;
			for (j = cell_low[c,1]; j <= cell_high[c,1]; j++)
			{
				eta = j * dnym1;
				exact_solution(xi, eta, zeta, temp, 0);
				for (m = 0; m < 5 ; m++)
				{
                    u[c, m, kk, jj, ii] = temp[m];
                }
                jj++;
			}
            kk++;
		}

		//---------------------------------------------------------------------
		// east face                                                      
		//---------------------------------------------------------------------

        c = slice[ncells - 1, 0];
        ii = 2 + cell_size[c, 0] - 1;
		xi = 1.0d;
        kk = 2;
		for (k = cell_low[c,2]; k <= cell_high[c,2]; k++)
		{
			zeta = k * dnzm1;
            jj = 2;
			for (j = cell_low[c,1]; j <= cell_high[c,1]; j++)
			{
                eta = j * dnym1;
				exact_solution(xi, eta, zeta, temp, 0);
				for (m = 0; m <= 4; m++)
				{
                    u[c, m, kk, jj, ii] = temp[m];
                }
                jj++;
			}
            kk++;
		}

		//---------------------------------------------------------------------
		// south face                                                 
		//---------------------------------------------------------------------

        c = slice[0,1];
        jj = 2;
		eta = 0.0d;
		kk = 2;
		for (k = cell_low[c,2]; k <= cell_high[c,2]; k++)
		{
			zeta = k * dnzm1;
            ii = 2;
			for (i = cell_low[c,0]; i <= cell_high[c,0]; i++)
			{
				xi = i * dnxm1;
				exact_solution(xi, eta, zeta, temp, 0);
				for (m = 0; m <= 4; m++)
				{
                    u[c, m, kk, jj, ii] = temp[m];
                }
                ii++;
			}
            kk++;
		}

		//---------------------------------------------------------------------
		// north face                                    
		//---------------------------------------------------------------------

        c = slice[ncells - 1, 1];
        jj = 2 + cell_size[c, 1] - 1;
		eta = 1.0d;
        kk = 2;
		for (k = cell_low[c,2]; k <= cell_high[c,2]; k++)
		{
			zeta = k * dnzm1;
            ii = 2;
            for (i = cell_low[c,0]; i <= cell_high[c,0]; i++)
            {
                xi = i * dnxm1;
                exact_solution(xi, eta, zeta, temp, 0);
                for (m = 0; m <= 4; m++)
                {
                    u[c, m, kk, jj, ii] = temp[m];
                }
                ii++;
            }
            kk++;
		}

		//---------------------------------------------------------------------
		// bottom face                                       
		//---------------------------------------------------------------------

        c = slice[0, 2];
        kk = 2;
		zeta = 0.0d;
		jj = 2;
        for (j = cell_low[c, 1]; j <= cell_high[c, 1]; j++)
        {
            eta = j * dnym1;
            ii = 2;
            for (i = cell_low[c, 0]; i <= cell_high[c, 0]; i++)
            {
                xi = i * dnxm1;
                exact_solution(xi, eta, zeta, temp, 0);
                for (m = 0; m <= 4; m++)
                {
                    u[c, m, kk, jj, ii] = temp[m];
                }
                ii++;
            }
            jj++;
        }

		//---------------------------------------------------------------------
		// top face     
		//---------------------------------------------------------------------

        c = slice[ncells - 1, 2];
        kk = 2 + cell_size[c, 2] - 1;
		zeta = 1.0d;
        jj = 2;
		for (j = cell_low[c,1]; j <= cell_high[c,1]; j++)
		{
            eta = j * dnym1;
            ii = 2;
			for (i = cell_low[c,0]; i <= cell_high[c,0]; i++)
			{
                xi = i * dnxm1;
                exact_solution(xi, eta, zeta, temp, 0);
				for (m = 0; m <= 4; m++)
				{
					u[c,m,kk,jj,ii] = temp[m];
                }
                ii++;
			}
            jj++;
		}
	}
	public void set_constants(int ndid)
	{
		ce[0,0] = 2.0d * (1.0d + ((double)ndid) * 0.01d);
		//    ce[0]=2.0;

		c1 = 1.4d;
		c2 = 0.4d;
		c3 = 0.1d;
		c4 = 1.0d;
		c5 = 1.4d;

		bt = Math.Sqrt(0.5d);

		dnxm1 = 1.0d / (grid_points[0] - 1);
		dnym1 = 1.0d / (grid_points[1] - 1);
		dnzm1 = 1.0d / (grid_points[2] - 1);

		c1c2 = c1 * c2;
		c1c5 = c1 * c5;
		c3c4 = c3 * c4;
		c1345 = c1c5 * c3c4;

		conz1 = (1.0d - c1c5);

		tx1 = 1.0d / (dnxm1 * dnxm1);
		tx2 = 1.0d / (2.0d * dnxm1);
		tx3 = 1.0d / dnxm1;

		ty1 = 1.0d / (dnym1 * dnym1);
		ty2 = 1.0d / (2.0d * dnym1);
		ty3 = 1.0d / dnym1;

		tz1 = 1.0d / (dnzm1 * dnzm1);
		tz2 = 1.0d / (2.0d * dnzm1);
		tz3 = 1.0d / dnzm1;

		dx1 = 0.75d;
		dx2 = 0.75d;
		dx3 = 0.75d;
		dx4 = 0.75d;
		dx5 = 0.75d;

		dy1 = 0.75d;
		dy2 = 0.75d;
		dy3 = 0.75d;
		dy4 = 0.75d;
		dy5 = 0.75d;

		dz1 = 1.0d;
		dz2 = 1.0d;
		dz3 = 1.0d;
		dz4 = 1.0d;
		dz5 = 1.0d;

		dxmax = dmax1(dx3, dx4);
		dymax = dmax1(dy2, dy4);
		dzmax = dmax1(dz2, dz3);

		dssp = 0.25d * dmax1(dx1, dmax1(dy1, dz1));

		c4dssp = 4.0d * dssp;
		c5dssp = 5.0d * dssp;

		dttx1 = dt * tx1;
		dttx2 = dt * tx2;
		dtty1 = dt * ty1;
		dtty2 = dt * ty2;
		dttz1 = dt * tz1;
		dttz2 = dt * tz2;

		c2dttx1 = 2.0d * dttx1;
		c2dtty1 = 2.0d * dtty1;
		c2dttz1 = 2.0d * dttz1;

		dtdssp = dt * dssp;

		comz1 = dtdssp;
		comz4 = 4.0d * dtdssp;
		comz5 = 5.0d * dtdssp;
		comz6 = 6.0d * dtdssp;

		c3c4tx3 = c3c4 * tx3;
		c3c4ty3 = c3c4 * ty3;
		c3c4tz3 = c3c4 * tz3;

		dx1tx1 = dx1 * tx1;
		dx2tx1 = dx2 * tx1;
		dx3tx1 = dx3 * tx1;
		dx4tx1 = dx4 * tx1;
		dx5tx1 = dx5 * tx1;

		dy1ty1 = dy1 * ty1;
		dy2ty1 = dy2 * ty1;
		dy3ty1 = dy3 * ty1;
		dy4ty1 = dy4 * ty1;
		dy5ty1 = dy5 * ty1;

		dz1tz1 = dz1 * tz1;
		dz2tz1 = dz2 * tz1;
		dz3tz1 = dz3 * tz1;
		dz4tz1 = dz4 * tz1;
		dz5tz1 = dz5 * tz1;

		c2iv = 2.5d;
		con43 = 4.0d / 3.0d;
		con16 = 1.0d / 6.0d;

		xxcon1 = c3c4tx3 * con43 * tx3;
		xxcon2 = c3c4tx3 * tx3;
		xxcon3 = c3c4tx3 * conz1 * tx3;
		xxcon4 = c3c4tx3 * con16 * tx3;
		xxcon5 = c3c4tx3 * c1c5 * tx3;

		yycon1 = c3c4ty3 * con43 * ty3;
		yycon2 = c3c4ty3 * ty3;
		yycon3 = c3c4ty3 * conz1 * ty3;
		yycon4 = c3c4ty3 * con16 * ty3;
		yycon5 = c3c4ty3 * c1c5 * ty3;

		zzcon1 = c3c4tz3 * con43 * tz3;
		zzcon2 = c3c4tz3 * tz3;
		zzcon3 = c3c4tz3 * conz1 * tz3;
		zzcon4 = c3c4tz3 * con16 * tz3;
		zzcon5 = c3c4tz3 * c1c5 * tz3;
	}
}
}





