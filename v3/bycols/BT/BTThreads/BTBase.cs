﻿/*
!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.0         !
!                                                                         !
!                       J A V A         V E R S I O N                     !
!                                                                         !
!                               B T B a s e                               !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    BTbase implements base class for BT benchmark.                       !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.0. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.0, including the Technical Report NAS-02-008    !
!    "Implementation of the NAS Parallel Benchmarks in Java",             !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!
!     Translation to Java and to MultiThreaded Code:                      !
!     Michael A. Frumkin                                                  !
!     Mathew Schultz                                                      !
!-------------------------------------------------------------------------!
*/

using System;
using System.Threading;
using NPB3_0_JAV;

namespace NPB3_0_JAV.BTThreads
{
    public class BTBase /*extends Thread */ {

        public static String BMName = "BT";
        public char CLASS = 'S';
        //npb class parameters  

        protected int IMAX = 0, JMAX = 0, KMAX = 0, problem_size = 0;

        protected int[] grid_points = { 0, 0, 0 };
        protected int niter_default = 0;
        protected double dt_default = 0.0;

        //array declarations from header.h
        protected double[,,] us, vs, ws, qs, rho_i, square;
        protected int jsize1, ksize1;

        protected double[,,,] forcing, u, rhs;
        protected double[] cv, cuf, q;
        protected int isize2, jsize2, ksize2;

        protected double[,] ue, buf;
        protected int jsize3;

        //here 5 are the dimensions of the CFD vector
        //(density,x_impuls,y_impuls,z_impuls,energy)
        protected const int isize4 = 5, jsize4 = 5 * 5, ksize4 = 5 * 5 * 3;
        protected const int aa = 0, bb = 1, cc = 2, BLOCK_SIZE = 5;

        // constants
        protected static double tx1, tx2, tx3, dt,
               ty1, ty2, ty3,
               tz1, tz2, tz3,
               dx1, dx2, dx3, dx4, dx5,
               dy1, dy2, dy3, dy4, dy5,
               dz1, dz2, dz3, dz4, dz5,
               dssp, dxmax, dymax, dzmax,
               xxcon1, xxcon2, xxcon3, xxcon4, xxcon5,
               dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1,
               yycon1, yycon2, yycon3, yycon4, yycon5,
               dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
               zzcon1, zzcon2, zzcon3, zzcon4, zzcon5,
               dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1,
               dnxm1, dnym1, dnzm1,
               c1c2, c1c5, c3c4, c1345, conz1,
               c1, c2, c3, c4, c5, c4dssp, c5dssp, dtdssp,
               dttx1, dttx2, dtty1, dtty2, dttz1, dttz2,
               c2dttx1, c2dtty1, c2dttz1,
               comz1, comz4, comz5, comz6, c3c4tx3, c3c4ty3,
               c3c4tz3, c2iv, con43, con16;

        protected static double[,] ce = {
               	{2.0, 1.0, 2.0, 2.0, 5.0},
               	{0.0, 0.0, 2.0, 2.0, 4.0},
               	{0.0, 0.0, 0.0, 0.0, 3.0},
               	{4.0, 0.0, 0.0, 0.0, 2.0},
               	{5.0, 1.0, 0.0, 0.0, 0.1},
               	{3.0, 2.0, 2.0, 2.0, 0.4},
               	{0.5, 3.0, 3.0, 3.0, 0.3},
               	{0.02, 0.01, 0.04, 0.03, 0.05},
               	{0.01, 0.03, 0.03, 0.05, 0.04},
               	{0.03, 0.02, 0.05, 0.04, 0.03},
               	{0.5, 0.4, 0.3, 0.2, 0.1},
               	{0.4, 0.3, 0.5, 0.1, 0.3},
               	{0.3, 0.5, 0.4, 0.3, 0.2}};

        //timer constants
        public bool timeron = false;
        public const int t_rhsx = 2, t_rhsy = 3, t_rhsz = 4,
                         t_xsolve = 6, t_ysolve = 7, t_zsolve = 8,
                         t_rdis1 = 9, t_rdis2 = 10, t_add = 11,
                         t_rhs = 5, t_last = 11, t_total = 1;

        public Timer timer = new Timer();

        public BTBase() { }

        public BTBase(char clss, int np)
        {
            CLASS = clss;
            num_threads = np;
            switch (clss)
            {
                case 'S':
                    problem_size = IMAX = JMAX = KMAX =
                                   grid_points[0] = grid_points[1] = grid_points[2] = 12;
                    dt_default = .01;
                    niter_default = 60;
                    CLASS = 'S';
                    break;
                case 'W':
                    problem_size = IMAX = JMAX = KMAX =
                                   grid_points[0] = grid_points[1] = grid_points[2] = 24;
                    dt_default = .0008;
                    niter_default = 200;
                    CLASS = 'W';
                    break;
                case 'A':
                    problem_size = IMAX = JMAX = KMAX =
                                   grid_points[0] = grid_points[1] = grid_points[2] = 64;
                    dt_default = .0008;
                    niter_default = 200;
                    CLASS = 'A';
                    break;
                case 'B':
                    problem_size = IMAX = JMAX = KMAX =
                                   grid_points[0] = grid_points[1] = grid_points[2] = 102;
                    dt_default = .0003;
                    niter_default = 200;
                    CLASS = 'B';
                    break;
                case 'C':
                    problem_size = IMAX = JMAX = KMAX =
                                   grid_points[0] = grid_points[1] = grid_points[2] = 162;
                    dt_default = .0001;
                    niter_default = 200;
                    CLASS = 'C';
                    break;
            }

            // set up arrays and array strides

            jsize1 = IMAX / 2 * 2 + 1;
            ksize1 = (JMAX / 2 * 2 + 1) * (IMAX / 2 * 2 + 1);
            us = new double[IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];
            vs = new double[IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];
            ws = new double[IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];
            qs = new double[IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];
            rho_i = new double[IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];
            square = new double[IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];

            isize2 = 5;
            jsize2 = 5 * (IMAX / 2 * 2 + 1);
            ksize2 = 5 * (IMAX / 2 * 2 + 1) * (JMAX / 2 * 2 + 1);

            forcing = new double[5, IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];
            u = new double[5, IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];
            rhs = new double[5, IMAX / 2 * 2 + 1, JMAX / 2 * 2 + 1, KMAX];

            cv = new double[problem_size + 2];
            cuf = new double[problem_size + 2];
            q = new double[problem_size + 2];

            jsize3 = (problem_size + 2);
            ue = new double[ 5,(problem_size + 2)];
            buf = new double[ 5,(problem_size + 2)];
        }


        // thread variables
        protected Thread master = null;
        protected int num_threads;

        /*
        protected RHSCompute[] rhscomputer;
        protected XSolver[] xsolver;
        protected YSolver[] ysolver;
        protected ZSolver[] zsolver;
        protected RHSAdder[] rhsadder;

        public void setupThreads(BT bt)
        {
            master = bt;
            if (num_threads > problem_size - 2)
                num_threads = problem_size - 2;


            int[] interval1 = new int[num_threads];
            int[] interval2 = new int[num_threads];
            set_interval(problem_size, interval1);
            set_interval(problem_size - 2, interval2);
            int[][] partition1 = new int[interval1.Length][2];
            int[][] partition2 = new int[interval2.Length][2];
            set_partition(0, interval1, partition1);
            set_partition(1, interval2, partition2);

            rhscomputer = new RHSCompute[num_threads];
            xsolver = new XSolver[num_threads];
            ysolver = new YSolver[num_threads];
            zsolver = new ZSolver[num_threads];
            rhsadder = new RHSAdder[num_threads];

            // create and start threads  
            for (int ii = 0; ii < num_threads; ii++)
            {
                rhscomputer[ii] = new RHSCompute(bt, partition1[ii][0], partition1[ii][1],
                                                  partition2[ii][0], partition2[ii][1]);
                rhscomputer[ii].id = ii;
                rhscomputer[ii].start();

                xsolver[ii] = new XSolver(bt, partition2[ii][0], partition2[ii][1]);
                xsolver[ii].id = ii;
                xsolver[ii].start();

                ysolver[ii] = new YSolver(bt, partition2[ii][0], partition2[ii][1]);
                ysolver[ii].id = ii;
                ysolver[ii].start();

                zsolver[ii] = new ZSolver(bt, partition2[ii][0], partition2[ii][1]);
                zsolver[ii].id = ii;
                zsolver[ii].start();

                rhsadder[ii] = new RHSAdder(bt, partition2[ii][0], partition2[ii][1]);
                rhsadder[ii].id = ii;
                rhsadder[ii].start();
            }
        }

        public void set_interval(int problem_size, int[] interval)
        {
            interval[0] = problem_size / num_threads;
            for (int i = 1; i < num_threads; i++) interval[i] = interval[0];
            int remainder = problem_size % num_threads;
            for (int i = 0; i < remainder; i++) interval[i]++;
        }

        public void set_partition(int start, int[] interval, int[][] array)
        {
            array[0][0] = start;
            if (start == 0) array[0][1] = interval[0] - 1;
            else array[0][1] = interval[0];

            for (int i = 1; i < interval.Length; i++)
            {
                array[i][0] = array[i - 1][1] + 1;
                array[i][1] = array[i - 1][1] + interval[i];
            }
        }
		*/
		
        public double dmax1(double a, double b)
        {
            if (a < b) return b; else return a;
        }

        public void set_constants()
        {
            c1 = 1.4;
            c2 = 0.4;
            c3 = 0.1;
            c4 = 1.0;
            c5 = 1.4;
            dnxm1 = 1.0 / (grid_points[0] - 1);
            dnym1 = 1.0 / (grid_points[1] - 1);
            dnzm1 = 1.0 / (grid_points[2] - 1);

            c1c2 = c1 * c2;
            c1c5 = c1 * c5;
            c3c4 = c3 * c4;
            c1345 = c1c5 * c3c4;

            conz1 = (1.0 - c1c5);

            tx1 = 1.0 / (dnxm1 * dnxm1);
            tx2 = 1.0 / (2.0 * dnxm1);
            tx3 = 1.0 / dnxm1;

            ty1 = 1.0 / (dnym1 * dnym1);
            ty2 = 1.0 / (2.0 * dnym1);
            ty3 = 1.0 / dnym1;

            tz1 = 1.0 / (dnzm1 * dnzm1);
            tz2 = 1.0 / (2.0 * dnzm1);
            tz3 = 1.0 / dnzm1;

            dx1 = 0.75;
            dx2 = 0.75;
            dx3 = 0.75;
            dx4 = 0.75;
            dx5 = 0.75;

            dy1 = 0.75;
            dy2 = 0.75;
            dy3 = 0.75;
            dy4 = 0.75;
            dy5 = 0.75;

            dz1 = 1.0;
            dz2 = 1.0;
            dz3 = 1.0;
            dz4 = 1.0;
            dz5 = 1.0;

            dxmax = dmax1(dx3, dx4);
            dymax = dmax1(dy2, dy4);
            dzmax = dmax1(dz2, dz3);

            dssp = 0.25 * dmax1(dx1, dmax1(dy1, dz1));

            c4dssp = 4.0 * dssp;
            c5dssp = 5.0 * dssp;

            dttx1 = dt * tx1;
            dttx2 = dt * tx2;
            dtty1 = dt * ty1;
            dtty2 = dt * ty2;
            dttz1 = dt * tz1;
            dttz2 = dt * tz2;

            c2dttx1 = 2.0 * dttx1;
            c2dtty1 = 2.0 * dtty1;
            c2dttz1 = 2.0 * dttz1;

            dtdssp = dt * dssp;

            comz1 = dtdssp;
            comz4 = 4.0 * dtdssp;
            comz5 = 5.0 * dtdssp;
            comz6 = 6.0 * dtdssp;

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

            c2iv = 2.5;
            con43 = 4.0 / 3.0;
            con16 = 1.0 / 6.0;

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
            dt = dt_default;
        }

        public void exact_solution(double xi, double eta, double zeta,
                                   double[] dtemp, int dtmpoffst)
        {
            for (int m = 0; m < 5; m++)
            {
                dtemp[m + dtmpoffst] = ce[0,m]
                                     + xi * (ce[1,m] + xi * (ce[4,m]
                                     + xi * (ce[7,m] + xi * ce[10,m])))
                                     + eta * (ce[2,m] + eta * (ce[5,m]
                                     + eta * (ce[8,m] + eta * ce[11,m])))
                                     + zeta * (ce[3,m] + zeta * (ce[6,m]
                                     + zeta * (ce[9,m] + zeta * ce[12,m])));
            }
        }
        public void initialize()
        {
            int i, j, k, m, ix, iy, iz;
            double xi, eta, zeta;
            double[] Pface = new double[5 * 3 * 2];
            double Pxi, Peta, Pzeta;
            double[] temp = new double[5];


            //---------------------------------------------------------------------
            //  Later (in compute_rhs) we compute 1/u for every element. A few of
            //  the corner elements are not used, but it convenient (and faster)
            //  to compute the whole thing with a simple loop. Make sure those
            //  values are nonzero by initializing the whole thing here.
            //---------------------------------------------------------------------
            for (k = 0; k <= grid_points[2] - 1; k++)
            {
                for (j = 0; j <= grid_points[1] - 1; j++)
                {
                    for (i = 0; i <= grid_points[0] - 1; i++)
                    {
                        for (m = 0; m <= 4; m++)
                        {
                            u[m, i, j, k] = 1.0;
                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //     first store the "interpolated" values everywhere on the grid    
            //---------------------------------------------------------------------

            for (k = 0; k <= grid_points[2] - 1; k++)
            {
                zeta = k * dnzm1;
                for (j = 0; j <= grid_points[1] - 1; j++)
                {
                    eta = j * dnym1;
                    for (i = 0; i <= grid_points[0] - 1; i++)
                    {
                        xi = i * dnxm1;

                        for (ix = 0; ix <= 1; ix++)
                        {
                            exact_solution(ix, eta, zeta,
                                           Pface, 0 + 0 * 5 + ix * 15);

                        }
                        for (iy = 0; iy <= 1; iy++)
                        {
                            exact_solution(xi, iy, zeta,
                                           Pface, 0 + 1 * 5 + iy * 15);
                        }

                        for (iz = 0; iz <= 1; iz++)
                        {
                            exact_solution(xi, eta, iz,
                                           Pface, 0 + 2 * 5 + iz * 15);
                        }

                        for (m = 0; m <= 4; m++)
                        {
                            Pxi = xi * Pface[m + 0 * 5 + 1 * 15] +
                              (1.0 - xi) * Pface[m + 0 * 5 + 0 * 15];
                            Peta = eta * Pface[m + 1 * 5 + 1 * 15] +
                              (1.0 - eta) * Pface[m + 1 * 5 + 0 * 15];
                            Pzeta = zeta * Pface[m + 2 * 5 + 1 * 15] +
                              (1.0 - zeta) * Pface[m + 2 * 5 + 0 * 15];

                            u[m, i, j, k] = Pxi + Peta + Pzeta -
                              Pxi * Peta - Pxi * Pzeta - Peta * Pzeta +
                              Pxi * Peta * Pzeta;

                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //     now store the exact values on the boundaries        
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //     west face                                                  
            //---------------------------------------------------------------------
            i = 0;
            xi = 0.0;
            for (k = 0; k <= grid_points[2] - 1; k++)
            {
                zeta = k * dnzm1;
                for (j = 0; j <= grid_points[1] - 1; j++)
                {
                    eta = j * dnym1;
                    exact_solution(xi, eta, zeta, temp, 0);
                    for (m = 0; m <= 4; m++)
                    {
                        u[m, i, j, k] = temp[m];
                    }
                }
            }

            //---------------------------------------------------------------------
            //     east face                                                      
            //---------------------------------------------------------------------

            i = grid_points[0] - 1;
            xi = 1.0;
            for (k = 0; k <= grid_points[2] - 1; k++)
            {
                zeta = k * dnzm1;
                for (j = 0; j <= grid_points[1] - 1; j++)
                {
                    eta = j * dnym1;
                    exact_solution(xi, eta, zeta, temp, 0);
                    for (m = 0; m <= 4; m++)
                    {
                        u[m, i, j, k] = temp[m];
                    }
                }
            }

            //---------------------------------------------------------------------
            //     south face                                                
            //---------------------------------------------------------------------
            j = 0;
            eta = 0.0;
            for (k = 0; k <= grid_points[2] - 1; k++)
            {
                zeta = k * dnzm1;
                for (i = 0; i <= grid_points[0] - 1; i++)
                {
                    xi = i * dnxm1;
                    exact_solution(xi, eta, zeta, temp, 0);
                    for (m = 0; m <= 4; m++)
                    {
                        u[m, i, j, k] = temp[m];
                    }
                }
            }

            //---------------------------------------------------------------------
            //     north face                                    
            //---------------------------------------------------------------------
            j = grid_points[1] - 1;
            eta = 1.0;
            for (k = 0; k <= grid_points[2] - 1; k++)
            {
                zeta = k * dnzm1;
                for (i = 0; i <= grid_points[0] - 1; i++)
                {
                    xi = i * dnxm1;
                    exact_solution(xi, eta, zeta, temp, 0);
                    for (m = 0; m <= 4; m++)
                    {
                        u[m, i, j, k] = temp[m];
                    }
                }
            }

            //---------------------------------------------------------------------
            //     bottom face                                      
            //---------------------------------------------------------------------
            k = 0;
            zeta = 0.0;
            for (i = 0; i <= grid_points[0] - 1; i++)
            {
                xi = i * dnxm1;
                for (j = 0; j <= grid_points[1] - 1; j++)
                {
                    eta = j * dnym1;
                    exact_solution(xi, eta, zeta, temp, 0);
                    for (m = 0; m <= 4; m++)
                    {
                        u[m, i, j, k] = temp[m];
                    }
                }
            }

            //---------------------------------------------------------------------
            //     top face    
            //---------------------------------------------------------------------
            k = grid_points[2] - 1;
            zeta = 1.0;
            for (i = 0; i <= grid_points[0] - 1; i++)
            {
                xi = i * dnxm1;
                for (j = 0; j <= grid_points[1] - 1; j++)
                {
                    eta = j * dnym1;
                    exact_solution(xi, eta, zeta, temp, 0);
                    for (m = 0; m <= 4; m++)
                    {
                        u[m, i, j, k] = temp[m];
                    }
                }
            }
        }

        public void lhsinit(double[,,,] lhs, int size)
        {
            int i, m, n;
            //---------------------------------------------------------------------
            //     zero the whole left hand side for starters
            //---------------------------------------------------------------------
            for (i = 0; i <= size; i += size)
            {
                for (m = 0; m <= 4; m++)
                {
                    for (n = 0; n <= 4; n++)
                    {
                        lhs[m, n, 0, i] = 0.0;
                        lhs[m, n, 1, i] = 0.0;
                        lhs[m, n, 2, i] = 0.0;
                    }
                }
            }
            //---------------------------------------------------------------------
            //     next, set all diagonal values to 1. This is overkill, but convenient
            //---------------------------------------------------------------------
            for (i = 0; i <= size; i += size)
            {
                for (m = 0; m <= 4; m++)
                {
                    lhs[m, m, 1, i] = 1.0;
                }
            }
        }

        public void matvec_sub(double[,,,] ablock, int blkoffst1, int blkoffst2,
                               double[,,,] avect, int avcoffst1, int avcoffst2, int avcoffst3,
                               double[,,,] bvect, int bvcoffst1, int bvcoffst2, int bvcoffst3)
        {

            for (int i = 0; i < 5; i++)
            {
                bvect[i, bvcoffst1, bvcoffst2, bvcoffst3] += -ablock[i, 0, blkoffst1, blkoffst2] * avect[0, avcoffst1, avcoffst2, avcoffst3]
                                     - ablock[i, 1, blkoffst1, blkoffst2] * avect[1, avcoffst1, avcoffst2, avcoffst3]
                                     - ablock[i, 2, blkoffst1, blkoffst2] * avect[2, avcoffst1, avcoffst2, avcoffst3]
                                     - ablock[i, 3, blkoffst1, blkoffst2] * avect[3, avcoffst1, avcoffst2, avcoffst3]
                                     - ablock[i, 4, blkoffst1, blkoffst2] * avect[4, avcoffst1, avcoffst2, avcoffst3];
            }
        }

        public void matmul_sub(double[,,,] ablock, int ablkoffst1, int ablkoffst2,
		                       double[,,,] bblock, int bblkoffst1, int bblkoffst2,
		                       double[,,,] cblock, int cblkoffst1, int cblkoffst2)
        {
            for (int j = 0; j < 5; j++)
            {


                cblock[0, j, cblkoffst1, cblkoffst2] += -ablock[0, 0, ablkoffst1, ablkoffst2] * bblock[0, j, bblkoffst1, bblkoffst2]
                                          - ablock[0, 1, ablkoffst1, ablkoffst2] * bblock[1, j, bblkoffst1, bblkoffst2]
                                          - ablock[0, 2, ablkoffst1, ablkoffst2] * bblock[2, j, bblkoffst1, bblkoffst2]
                                          - ablock[0, 3, ablkoffst1, ablkoffst2] * bblock[3, j, bblkoffst1, bblkoffst2]
                                          - ablock[0, 4, ablkoffst1, ablkoffst2] * bblock[4, j, bblkoffst1, bblkoffst2];
				
				cblock[1, j, cblkoffst1, cblkoffst2] += -ablock[1, 0, ablkoffst1, ablkoffst2] * bblock[0, j, bblkoffst1, bblkoffst2]
                                          - ablock[1, 1, ablkoffst1, ablkoffst2] * bblock[1, j, bblkoffst1, bblkoffst2]
                                          - ablock[1, 2, ablkoffst1, ablkoffst2] * bblock[2, j, bblkoffst1, bblkoffst2]
                                          - ablock[1, 3, ablkoffst1, ablkoffst2] * bblock[3, j, bblkoffst1, bblkoffst2]
                                          - ablock[1, 4, ablkoffst1, ablkoffst2] * bblock[4, j, bblkoffst1, bblkoffst2];

				cblock[2, j, cblkoffst1, cblkoffst2] += -ablock[2, 0, ablkoffst1, ablkoffst2] * bblock[0, j, bblkoffst1, bblkoffst2]
                                          - ablock[2, 1, ablkoffst1, ablkoffst2] * bblock[1, j, bblkoffst1, bblkoffst2]
                                          - ablock[2, 2, ablkoffst1, ablkoffst2] * bblock[2, j, bblkoffst1, bblkoffst2]
                                          - ablock[2, 3, ablkoffst1, ablkoffst2] * bblock[3, j, bblkoffst1, bblkoffst2]
                                          - ablock[2, 4, ablkoffst1, ablkoffst2] * bblock[4, j, bblkoffst1, bblkoffst2];

				cblock[3, j, cblkoffst1, cblkoffst2] += -ablock[3, 0, ablkoffst1, ablkoffst2] * bblock[0, j, bblkoffst1, bblkoffst2]
                                          - ablock[3, 1, ablkoffst1, ablkoffst2] * bblock[1, j, bblkoffst1, bblkoffst2]
                                          - ablock[3, 2, ablkoffst1, ablkoffst2] * bblock[2, j, bblkoffst1, bblkoffst2]
                                          - ablock[3, 3, ablkoffst1, ablkoffst2] * bblock[3, j, bblkoffst1, bblkoffst2]
                                          - ablock[3, 4, ablkoffst1, ablkoffst2] * bblock[4, j, bblkoffst1, bblkoffst2];

				cblock[4, j, cblkoffst1, cblkoffst2] += -ablock[4, 0, ablkoffst1, ablkoffst2] * bblock[0, j, bblkoffst1, bblkoffst2]
                                          - ablock[4, 1, ablkoffst1, ablkoffst2] * bblock[1, j, bblkoffst1, bblkoffst2]
                                          - ablock[4, 2, ablkoffst1, ablkoffst2] * bblock[2, j, bblkoffst1, bblkoffst2]
                                          - ablock[4, 3, ablkoffst1, ablkoffst2] * bblock[3, j, bblkoffst1, bblkoffst2]
                                          - ablock[4, 4, ablkoffst1, ablkoffst2] * bblock[4, j, bblkoffst1, bblkoffst2];
            }
        }

        public void binvcrhs(double[,,,] lhss, int lhsoffst1, int lhsoffst2,
                             double[,,,] c, int coffst1, int coffst2,
		                     double[,,,] r, int roffst1, int roffst2, int roffst3)
        {
            double pivot;
            double coeff;

            pivot = 1.0 / lhss[0, 0, lhsoffst1, lhsoffst2];
            lhss[0, 1, lhsoffst1, lhsoffst2] *= pivot;
            lhss[0, 2, lhsoffst1, lhsoffst2] *= pivot;
            lhss[0, 3, lhsoffst1, lhsoffst2] *= pivot;
            lhss[0, 4, lhsoffst1, lhsoffst2] *= pivot;
            c[0, 0, coffst1, coffst2] *= pivot;
            c[0, 1, coffst1, coffst2] *= pivot;
            c[0, 2, coffst1, coffst2] *= pivot;
            c[0, 3, coffst1, coffst2] *= pivot;
            c[0, 4, coffst1, coffst2] *= pivot;
            r[0, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[1, 0, lhsoffst1, lhsoffst2];
            lhss[1, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[1, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[1, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[1, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            c[1, 0, coffst1, coffst2] -= coeff * c[0, 0, coffst1, coffst2];
            c[1, 1, coffst1, coffst2] -= coeff * c[0, 1, coffst1, coffst2];
            c[1, 2, coffst1, coffst2] -= coeff * c[0, 2, coffst1, coffst2];
            c[1, 3, coffst1, coffst2] -= coeff * c[0, 3, coffst1, coffst2];
            c[1, 4, coffst1, coffst2] -= coeff * c[0, 4, coffst1, coffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];

			coeff = lhss[2, 0, lhsoffst1, lhsoffst2];
            lhss[2, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[2, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[2, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[2, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            c[2, 0, coffst1, coffst2] -= coeff * c[0, 0, coffst1, coffst2];
            c[2, 1, coffst1, coffst2] -= coeff * c[0, 1, coffst1, coffst2];
            c[2, 2, coffst1, coffst2] -= coeff * c[0, 2, coffst1, coffst2];
            c[2, 3, coffst1, coffst2] -= coeff * c[0, 3, coffst1, coffst2];
            c[2, 4, coffst1, coffst2] -= coeff * c[0, 4, coffst1, coffst2];
            r[2, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];
			coeff = lhss[3, 0, lhsoffst1, lhsoffst2];
            lhss[3, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[3, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[3, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            c[3, 0, coffst1, coffst2] -= coeff * c[0, 0, coffst1, coffst2];
            c[3, 1, coffst1, coffst2] -= coeff * c[0, 1, coffst1, coffst2];
            c[3, 2, coffst1, coffst2] -= coeff * c[0, 2, coffst1, coffst2];
            c[3, 3, coffst1, coffst2] -= coeff * c[0, 3, coffst1, coffst2];
            c[3, 4, coffst1, coffst2] -= coeff * c[0, 4, coffst1, coffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];

            coeff = lhss[4, 0, lhsoffst1, lhsoffst2];
            lhss[4, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[4, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[4, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            c[4, 0, coffst1, coffst2] -= coeff * c[0, 0, coffst1, coffst2];
            c[4, 1, coffst1, coffst2] -= coeff * c[0, 1, coffst1, coffst2];
            c[4, 2, coffst1, coffst2] -= coeff * c[0, 2, coffst1, coffst2];
            c[4, 3, coffst1, coffst2] -= coeff * c[0, 3, coffst1, coffst2];
            c[4, 4, coffst1, coffst2] -= coeff * c[0, 4, coffst1, coffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];
			
			
            pivot = 1.0 / lhss[1, 1, lhsoffst1, lhsoffst2];
            lhss[1, 2, lhsoffst1, lhsoffst2] *= pivot;
            lhss[1, 3, lhsoffst1, lhsoffst2] *= pivot;
            lhss[1, 4, lhsoffst1, lhsoffst2] *= pivot;
            c[1, 0, coffst1, coffst2] *= pivot;
            c[1, 1, coffst1, coffst2] *= pivot;
            c[1, 2, coffst1, coffst2] *= pivot;
            c[1, 3, coffst1, coffst2] *= pivot;
            c[1, 4, coffst1, coffst2] *= pivot;
            r[1, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[0, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[0, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[0, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            c[0, 0, coffst1, coffst2] -= coeff * c[1, 0, coffst1, coffst2];
            c[0, 1, coffst1, coffst2] -= coeff * c[1, 1, coffst1, coffst2];
            c[0, 2, coffst1, coffst2] -= coeff * c[1, 2, coffst1, coffst2];
            c[0, 3, coffst1, coffst2] -= coeff * c[1, 3, coffst1, coffst2];
            c[0, 4, coffst1, coffst2] -= coeff * c[1, 4, coffst1, coffst2];
            r[0, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];			

            coeff = lhss[2, 1, lhsoffst1, lhsoffst2];
            lhss[2, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[2, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[2, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            c[2, 0, coffst1, coffst2] -= coeff * c[1, 0, coffst1, coffst2];
            c[2, 1, coffst1, coffst2] -= coeff * c[1, 1, coffst1, coffst2];
            c[2, 2, coffst1, coffst2] -= coeff * c[1, 2, coffst1, coffst2];
            c[2, 3, coffst1, coffst2] -= coeff * c[1, 3, coffst1, coffst2];
            c[2, 4, coffst1, coffst2] -= coeff * c[1, 4, coffst1, coffst2];
			r[2, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];

            coeff = lhss[3, 1, lhsoffst1, lhsoffst2];
            lhss[3, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[3, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            c[3, 0, coffst1, coffst2] -= coeff * c[1, 0, coffst1, coffst2];
            c[3, 1, coffst1, coffst2] -= coeff * c[1, 1, coffst1, coffst2];
            c[3, 2, coffst1, coffst2] -= coeff * c[1, 2, coffst1, coffst2];
            c[3, 3, coffst1, coffst2] -= coeff * c[1, 3, coffst1, coffst2];
            c[3, 4, coffst1, coffst2] -= coeff * c[1, 4, coffst1, coffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];
			
            coeff = lhss[4, 1, lhsoffst1, lhsoffst2];
            lhss[4, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[4, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            c[4, 0, coffst1, coffst2] -= coeff * c[1, 0, coffst1, coffst2];
            c[4, 1, coffst1, coffst2] -= coeff * c[1, 1, coffst1, coffst2];
            c[4, 2, coffst1, coffst2] -= coeff * c[1, 2, coffst1, coffst2];
            c[4, 3, coffst1, coffst2] -= coeff * c[1, 3, coffst1, coffst2];
			c[4, 4, coffst1, coffst2] -= coeff * c[1, 4, coffst1, coffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];

            pivot = 1.0 / lhss[2, 2, lhsoffst1, lhsoffst2];
            lhss[2, 3, lhsoffst1, lhsoffst2] *= pivot;
            lhss[2, 4, lhsoffst1, lhsoffst2] *= pivot;
            c[2, 0, coffst1, coffst2] *= pivot;
            c[2, 1, coffst1, coffst2] *= pivot;
            c[2, 2, coffst1, coffst2] *= pivot;
            c[2, 3, coffst1, coffst2] *= pivot;
            c[2, 4, coffst1, coffst2] *= pivot;
            r[2, roffst1, roffst2, roffst3] *= pivot;
            coeff = lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[0, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[0, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            c[0, 0, coffst1, coffst2] -= coeff * c[2, 0, coffst1, coffst2];
            c[0, 1, coffst1, coffst2] -= coeff * c[2, 1, coffst1, coffst2];
            c[0, 2, coffst1, coffst2] -= coeff * c[2, 2, coffst1, coffst2];
            c[0, 3, coffst1, coffst2] -= coeff * c[2, 3, coffst1, coffst2];
            c[0, 4, coffst1, coffst2] -= coeff * c[2, 4, coffst1, coffst2];
            r[0, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];
            coeff = lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[1, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[1, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            c[1, 0, coffst1, coffst2] -= coeff * c[2, 0, coffst1, coffst2];
            c[1, 1, coffst1, coffst2] -= coeff * c[2, 1, coffst1, coffst2];
            c[1, 2, coffst1, coffst2] -= coeff * c[2, 2, coffst1, coffst2];
            c[1, 3, coffst1, coffst2] -= coeff * c[2, 3, coffst1, coffst2];
            c[1, 4, coffst1, coffst2] -= coeff * c[2, 4, coffst1, coffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];
            coeff = lhss[3, 2, lhsoffst1, lhsoffst2];
            lhss[3, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            c[3, 0, coffst1, coffst2] -= coeff * c[2, 0, coffst1, coffst2];
            c[3, 1, coffst1, coffst2] -= coeff * c[2, 1, coffst1, coffst2];
            c[3, 2, coffst1, coffst2] -= coeff * c[2, 2, coffst1, coffst2];
            c[3, 3, coffst1, coffst2] -= coeff * c[2, 3, coffst1, coffst2];
            c[3, 4, coffst1, coffst2] -= coeff * c[2, 4, coffst1, coffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];
            coeff = lhss[4, 2, lhsoffst1, lhsoffst2];
            lhss[4, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            c[4, 0, coffst1, coffst2] -= coeff * c[2, 0, coffst1, coffst2];
            c[4, 1, coffst1, coffst2] -= coeff * c[2, 1, coffst1, coffst2];
            c[4, 2, coffst1, coffst2] -= coeff * c[2, 2, coffst1, coffst2];
            c[4, 3, coffst1, coffst2] -= coeff * c[2, 3, coffst1, coffst2];
            c[4, 4, coffst1, coffst2] -= coeff * c[2, 4, coffst1, coffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];


            pivot = 1.0 / lhss[3, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] *= pivot;
            c[3, 0, coffst1, coffst2] *= pivot;
            c[3, 1, coffst1, coffst2] *= pivot;
            c[3, 2, coffst1, coffst2] *= pivot;
            c[3, 3, coffst1, coffst2] *= pivot;
            c[3, 4, coffst1, coffst2] *= pivot;
            r[3, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[0, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            c[0, 0, coffst1, coffst2] -= coeff * c[3, 0, coffst1, coffst2];
            c[0, 1, coffst1, coffst2] -= coeff * c[3, 1, coffst1, coffst2];
            c[0, 2, coffst1, coffst2] -= coeff * c[3, 2, coffst1, coffst2];
            c[0, 3, coffst1, coffst2] -= coeff * c[3, 3, coffst1, coffst2];
            c[0, 4, coffst1, coffst2] -= coeff * c[3, 4, coffst1, coffst2];
            r[0, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];

            coeff = lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[1, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            c[1, 0, coffst1, coffst2] -= coeff * c[3, 0, coffst1, coffst2];
            c[1, 1, coffst1, coffst2] -= coeff * c[3, 1, coffst1, coffst2];
            c[1, 2, coffst1, coffst2] -= coeff * c[3, 2, coffst1, coffst2];
            c[1, 3, coffst1, coffst2] -= coeff * c[3, 3, coffst1, coffst2];
            c[1, 4, coffst1, coffst2] -= coeff * c[3, 4, coffst1, coffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];

            coeff = lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[2, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            c[2, 0, coffst1, coffst2] -= coeff * c[3, 0, coffst1, coffst2];
            c[2, 1, coffst1, coffst2] -= coeff * c[3, 1, coffst1, coffst2];
            c[2, 2, coffst1, coffst2] -= coeff * c[3, 2, coffst1, coffst2];
            c[2, 3, coffst1, coffst2] -= coeff * c[3, 3, coffst1, coffst2];
            c[2, 4, coffst1, coffst2] -= coeff * c[3, 4, coffst1, coffst2];
            r[2, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];
            coeff = lhss[4, 3, lhsoffst1, lhsoffst2];
			lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            c[4, 0, coffst1, coffst2] -= coeff * c[3, 0, coffst1, coffst2];
            c[4, 1, coffst1, coffst2] -= coeff * c[3, 1, coffst1, coffst2];
            c[4, 2, coffst1, coffst2] -= coeff * c[3, 2, coffst1, coffst2];
            c[4, 3, coffst1, coffst2] -= coeff * c[3, 3, coffst1, coffst2];
            c[4, 4, coffst1, coffst2] -= coeff * c[3, 4, coffst1, coffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];

            pivot = 1.0 / lhss[4, 4, lhsoffst1, lhsoffst2];
            c[4, 0, coffst1, coffst2] *= pivot;
            c[4, 1, coffst1, coffst2] *= pivot;
            c[4, 2, coffst1, coffst2] *= pivot;
            c[4, 3, coffst1, coffst2] *= pivot;
            c[4, 4, coffst1, coffst2] *= pivot;
			r[4, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[0, 4, lhsoffst1, lhsoffst2];
            c[0, 0, coffst1, coffst2] -= coeff * c[4, 0, coffst1, coffst2];
            c[0, 1, coffst1, coffst2] -= coeff * c[4, 1, coffst1, coffst2];
            c[0, 2, coffst1, coffst2] -= coeff * c[4, 2, coffst1, coffst2];
            c[0, 3, coffst1, coffst2] -= coeff * c[4, 3, coffst1, coffst2];
            c[0, 4, coffst1, coffst2] -= coeff * c[4, 4, coffst1, coffst2];
			r[0, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];

            coeff = lhss[1, 4, lhsoffst1, lhsoffst2];
            c[1, 0, coffst1, coffst2] -= coeff * c[4, 0, coffst1, coffst2];
            c[1, 1, coffst1, coffst2] -= coeff * c[4, 1, coffst1, coffst2];
            c[1, 2, coffst1, coffst2] -= coeff * c[4, 2, coffst1, coffst2];
            c[1, 3, coffst1, coffst2] -= coeff * c[4, 3, coffst1, coffst2];
            c[1, 4, coffst1, coffst2] -= coeff * c[4, 4, coffst1, coffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];

            coeff = lhss[2, 4, lhsoffst1, lhsoffst2];
            c[2, 0, coffst1, coffst2] -= coeff * c[4, 0, coffst1, coffst2];
            c[2, 1, coffst1, coffst2] -= coeff * c[4, 1, coffst1, coffst2];
            c[2, 2, coffst1, coffst2] -= coeff * c[4, 2, coffst1, coffst2];
            c[2, 3, coffst1, coffst2] -= coeff * c[4, 3, coffst1, coffst2];
            c[2, 4, coffst1, coffst2] -= coeff * c[4, 4, coffst1, coffst2];
            r[2, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];

            coeff = lhss[3, 4, lhsoffst1, lhsoffst2];
            c[3, 0, coffst1, coffst2] -= coeff * c[4, 0, coffst1, coffst2];
            c[3, 1, coffst1, coffst2] -= coeff * c[4, 1, coffst1, coffst2];
            c[3, 2, coffst1, coffst2] -= coeff * c[4, 2, coffst1, coffst2];
            c[3, 3, coffst1, coffst2] -= coeff * c[4, 3, coffst1, coffst2];
            c[3, 4, coffst1, coffst2] -= coeff * c[4, 4, coffst1, coffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];
        }

        public void binvrhs(double[,,,] lhss, int lhsoffst1, int lhsoffst2, double[,,,] r, int roffst1, int roffst2, int roffst3)
        {
            double pivot;
            double coeff;

            pivot = 1 / lhss[0, 0, lhsoffst1, lhsoffst2];
            lhss[0, 1, lhsoffst1, lhsoffst2] *= pivot;
            lhss[0, 2, lhsoffst1, lhsoffst2] *= pivot;
            lhss[0, 3, lhsoffst1, lhsoffst2] *= pivot;
            lhss[0, 4, lhsoffst1, lhsoffst2] *= pivot;
            r[0, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[1, 0, lhsoffst1, lhsoffst2];
            lhss[1, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[1, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[1, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[1, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];

            coeff = lhss[2, 0, lhsoffst1, lhsoffst2];
            lhss[2, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[2, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[2, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[2, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            r[2, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];

            coeff = lhss[3, 0, lhsoffst1, lhsoffst2];
            lhss[3, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[3, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[3, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];

            coeff = lhss[4, 0, lhsoffst1, lhsoffst2];
            lhss[4, 1, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[4, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[4, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[0, 4, lhsoffst1, lhsoffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[0, roffst1, roffst2, roffst3];
			

            pivot = 1 / lhss[1, 1, lhsoffst1, lhsoffst2];
            lhss[1, 2, lhsoffst1, lhsoffst2] *= pivot;
            lhss[1, 3, lhsoffst1, lhsoffst2] *= pivot;
            lhss[1, 4, lhsoffst1, lhsoffst2] *= pivot;
            r[1, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[0, 1, lhsoffst1, lhsoffst2];
            lhss[0, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[0, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[0, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            r[0, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];

            coeff = lhss[2, 1, lhsoffst1, lhsoffst2];
            lhss[2, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[2, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[2, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            r[2, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];

            coeff = lhss[3, 1, lhsoffst1, lhsoffst2];
            lhss[3, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[3, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];

            coeff = lhss[4, 1, lhsoffst1, lhsoffst2];
            lhss[4, 2, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[4, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[1, 4, lhsoffst1, lhsoffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[1, roffst1, roffst2, roffst3];
			

            pivot = 1 / lhss[2, 2, lhsoffst1, lhsoffst2];
            lhss[2, 3, lhsoffst1, lhsoffst2] *= pivot;
            lhss[2, 4, lhsoffst1, lhsoffst2] *= pivot;
            r[2, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[0, 2, lhsoffst1, lhsoffst2];
            lhss[0, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[0, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            r[0, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];

            coeff = lhss[1, 2, lhsoffst1, lhsoffst2];
            lhss[1, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[1, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];

            coeff = lhss[3, 2, lhsoffst1, lhsoffst2];
            lhss[3, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];

            coeff = lhss[4, 2, lhsoffst1, lhsoffst2];
            lhss[4, 3, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[2, 4, lhsoffst1, lhsoffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[2, roffst1, roffst2, roffst3];
			

            pivot = 1 / lhss[3, 3, lhsoffst1, lhsoffst2];
            lhss[3, 4, lhsoffst1, lhsoffst2] *= pivot;
            r[3, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[0, 3, lhsoffst1, lhsoffst2];
            lhss[0, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            r[0, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];

            coeff = lhss[1, 3, lhsoffst1, lhsoffst2];
            lhss[1, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];

            coeff = lhss[2, 3, lhsoffst1, lhsoffst2];
            lhss[2, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            r[2, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];

            coeff = lhss[4, 3, lhsoffst1, lhsoffst2];
            lhss[4, 4, lhsoffst1, lhsoffst2] -= coeff * lhss[3, 4, lhsoffst1, lhsoffst2];
            r[4, roffst1, roffst2, roffst3] -= coeff * r[3, roffst1, roffst2, roffst3];


            pivot = 1 / lhss[4, 4, lhsoffst1, lhsoffst2];
            r[4, roffst1, roffst2, roffst3] *= pivot;

            coeff = lhss[0, 4, lhsoffst1, lhsoffst2];
            r[0, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];

            coeff = lhss[1, 4, lhsoffst1, lhsoffst2];
            r[1, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];

            coeff = lhss[2, 4, lhsoffst1, lhsoffst2];
            r[2, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];

            coeff = lhss[3, 4, lhsoffst1, lhsoffst2];
            r[3, roffst1, roffst2, roffst3] -= coeff * r[4, roffst1, roffst2, roffst3];
        }
    }
}
 