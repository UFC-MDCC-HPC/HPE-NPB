/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                                  S P                                    !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is a serial/multithreaded version of                  !
!    the NPB3_0_JAV BT code.                                              !
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
!									  !
! Authors: R. Van der Wijngaart 					  !
!	   T. Harris							  !
!	   M. Yarrow							  !
! Modified for PBN (Programming Baseline for NPB):			  !
!	   H. Jin							  !
! Translation to Java and MultiThreaded Code				  !
!	   M. Frumkin							  !
!	   M. Schultz							  !
!-------------------------------------------------------------------------!
*/

using System;
using System.IO;
using NPB3_0_JAV.BTThreads;
using NPB3_0_JAV.BMInOut;
using MPI;

namespace NPB3_0_JAV{
    public class BT : BTBase{
        public int bid = -1;
        public static int t_total = 1;
        public BMResults results;
        public BT(char clss)
            : base(clss)
        {
        }
        public static void Main(String[] argv)
        {
            BT bt = null;

            BMArgs.ParseCmdLineArgs(argv, BMName);
            char CLSS = BMArgs.CLASS;

            try
            {
                bt = new BT(CLSS);
            }
            catch (OutOfMemoryException e)
            {
                Console.WriteLine(e.Message);
                BMArgs.outOfMemoryMessage();
                System.Environment.Exit(0);
            }
            bt.runBenchMark();

        }

        public void run() { runBenchMark(); }

        public void runBenchMark()
        {
            if (!active)
            {
                Console.WriteLine("not active !");
                System.Environment.Exit(0);
            }

            int niter = -1;
            if (node == root)
            {
                BMArgs.Banner(BMName, CLASS, false, total_nodes);

            }

            //---------------------------------------------------------------------
            //      Read input file (if it exists), else take
            //      defaults from parameters
            //---------------------------------------------------------------------
            niter = getInputPars();

            comm_setup.Broadcast<int>(ref niter, root);
            comm_setup.Broadcast<int>(ref dp_type, root);
            comm_setup.Broadcast<int>(ref grid_points, root);
            
            make_set();

            for (int c = 0; c < ncells; c++)
            {
                if ((cell_size[c, 0] > IMAX) ||
                    (cell_size[c, 1] > JMAX) ||
                    (cell_size[c, 2] > KMAX))
                {
                    Console.WriteLine("Problem size too big for compiled array sizes");
                    System.Environment.Exit(0);
                }
            }

            set_constants(0);
            initialize();
            lhsinit();
            exact_rhs();
            compute_buffer_size(5);

            //---------------------------------------------------------------------
            //      do one time step to touch all code, and reinitialize
            //---------------------------------------------------------------------
            adi();
            initialize();

            //---------------------------------------------------------------------
            //      Synchronize before placing time stamp
            //---------------------------------------------------------------------
            // mpi_barrier(comm_setup, error);
            comm_setup.Barrier();

            timer.resetAllTimers();
            timer.start(t_total);

            Console.WriteLine("STARTING"); Console.Out.Flush();

            for (int step = 1; step <= niter; step++)
            {
                if (node == 0 && (step % 20 == 0 || step == 1 || step == niter))
                {
                    Console.WriteLine("Time step " + step);
                }
                adi();
            }
            timer.stop(1);
            int verified = verify(niter);

            double tmax = comm_setup.Reduce<double>(timer.readTimer(t_total), Operation<double>.Max, root);

            if (node == root)
            {
                double time = timer.readTimer(t_total);
                results = new BMResults(BMName,
                              CLASS,
                              grid_points[0],
                              grid_points[1],
                              grid_points[2],
                              niter,
                              time,
                              getMFLOPS(time, niter),
                              "floating point",
                              verified,
                              true,
                              total_nodes,
                              bid);
                results.print();
            }

            worldcomm.Barrier();
            mpi.Dispose();

        }

        private void compute_buffer_size(int dim)
        {
            int c, face_size;

            if (ncells == 1) return;

            //---------------------------------------------------------------------
            //      compute the actual sizes of the buffers; note that there is 
            //      always one cell face that doesn't need buffer space, because it 
            //      is at the boundary of the grid
            //---------------------------------------------------------------------

            west_size = 0;
            east_size = 0;

            for (c = 0; c < ncells; c++)
            {
                face_size = cell_size[c, 1] * cell_size[c, 2] * dim * 2;
                if (cell_coord[c, 0] != 0) west_size = west_size + face_size;
                if (cell_coord[c, 0] != ncells - 1) east_size = east_size + face_size;
            }

            north_size = 0;
            south_size = 0;
            for (c = 0; c < ncells; c++)
            {
                face_size = cell_size[c, 0] * cell_size[c, 2] * dim * 2;
                if (cell_coord[c, 1] != 0) south_size = south_size + face_size;
                if (cell_coord[c, 1] != ncells - 1) north_size = north_size + face_size;
            }

            top_size = 0;
            bottom_size = 0;
            for (c = 0; c < ncells; c++)
            {
                face_size = cell_size[c, 0] * cell_size[c, 1] * dim * 2;
                if (cell_coord[c, 2] != 0) bottom_size = bottom_size + face_size;
                if (cell_coord[c, 2] != ncells - 1) top_size = top_size + face_size;
            }

            start_send_west = 1;
            start_send_east = start_send_west + west_size;
            start_send_south = start_send_east + east_size;
            start_send_north = start_send_south + south_size;
            start_send_bottom = start_send_north + north_size;
            start_send_top = start_send_bottom + bottom_size;
            start_recv_west = 1;
            start_recv_east = start_recv_west + west_size;
            start_recv_south = start_recv_east + east_size;
            start_recv_north = start_recv_south + south_size;
            start_recv_bottom = start_recv_north + north_size;
            start_recv_top = start_recv_bottom + bottom_size;
        }

        public double getMFLOPS(double total_time, int niter)
        {
            double mflops = 0.0d;
            if (total_time > 0)
            {
                mflops = (881.174d * Math.Pow(problem_size, 3)
                        - 4683.91d * Math.Pow(problem_size, 2)
                        + 11484.5d * problem_size
                        - 19272.4d) * niter / (1.0d * 1000000.0d);
            }
            return mflops;
        }

        public void adi()
        {
            copy_faces();
            x_solve();
            y_solve();
            z_solve();
            add();
        }

        public int getInputPars()
        {
            int niter = 0;
            if (File.Exists("inputbt.data"))
            {
                FileStream f2 = new FileStream("inputbt.data", System.IO.FileMode.Open);
                try
                {
                    StreamReader sss = new StreamReader(f2);
                    Console.WriteLine("Reading from input file inputbt.data");
                    ;
                    niter = int.Parse(sss.ReadLine());
                    dt = double.Parse(sss.ReadLine());
                    grid_points[0] = int.Parse(sss.ReadLine());
                    grid_points[1] = int.Parse(sss.ReadLine());
                    grid_points[2] = int.Parse(sss.ReadLine());
                    //			fis.close();
                }
                catch (Exception e)
                {
                    Console.Error.WriteLine(e.Message);
                }
            }
            else
            {
                Console.WriteLine("No input file inputbt.data," +
                                   "Using compiled defaults");
                niter = niter_default;
                dt = dt_default;
                grid_points[0] = problem_size;
                grid_points[1] = problem_size;
                grid_points[2] = problem_size;
            }
            Console.WriteLine("Size: " + grid_points[0]
                                 + " X " + grid_points[1]
                     + " X " + grid_points[2]);
            for (int c = 0; c < ncells; c++)
            {
                if ((cell_size[c, 0] > IMAX) ||
                    (cell_size[c, 1] > JMAX) ||
                    (cell_size[c, 2] > KMAX))
                {
                    Console.WriteLine("Problem size too big for compiled array sizes");
                    System.Environment.Exit(0);
                }
            }
            Console.WriteLine("Iterations: " + niter + " dt: " + dt);

            return niter;
        }

        public void add()
        {
            int c, i, j, k, m;

            for (c = 0; c < ncells; c++)
            {
                for (k = start[c, 2]; k < 2 + cell_size[c, 2] - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < 2 + cell_size[c, 1] - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < 2 + cell_size[c, 0] - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                u[c, k, j, i, m] += rhs[c, k, j, i, m];
                            }
                        }
                    }
                }
            }
        }

        public void error_norm(double[] rms)
        {
            int c, i, j, k, ii, jj, kk, m, d;
            double xi, eta, zeta;
            double[] u_exact = new double[5];
            double add;
            double[] rms_work = new double[5];

            for (m = 0; m < 5; m++)
            {
                rms_work[m] = 0.0d;
            }

            for (c = 0; c < ncells; c++)
            {
                kk = 2;
                for (k = cell_low[c, 2]; k <= cell_high[c, 2]; k++)
                {
                    zeta = k * dnzm1;
                    jj = 2;
                    for (j = cell_low[c, 1]; j <= cell_high[c, 1]; j++)
                    {
                        eta = j * dnym1;
                        ii = 2;
                        for (i = cell_low[c, 0]; i <= cell_high[c, 0]; i++)
                        {
                            xi = i * dnxm1;
                            exact_solution(xi, eta, zeta, u_exact, 0);

                            for (m = 0; m < 5; m++)
                            {
                                add = u[c, kk, jj, ii, m] - u_exact[m];
                                rms_work[m] += add * add;
                            }
                            ii++;
                        }
                        jj++;
                    }
                    kk++;
                }
            }

            comm_setup.Allreduce<double>(rms_work, Operation<double>.Add, ref rms);

            for (m = 0; m < 5; m++)
            {
                for (d = 0; d < 3; d++)
                {
                    rms[m] /= (grid_points[d] - 2);
                }
                rms[m] = Math.Sqrt(rms[m]);
            }
        }

        public void rhs_norm(double[] rms)
        {
            int c, i, j, k, d, m, ksize, jsize, isize;
            double add;
            double[] rms_work = new double[5];

            for (m = 0; m < 5; m++)
            {
                rms_work[m] = 0.0d;
            }

            for (c = 0; c < ncells; c++)
            {
                ksize = cell_size[c, 2] + 2;
                jsize = cell_size[c, 1] + 2;
                isize = cell_size[c, 0] + 2;

                for (m = 0; m < 5; m++)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                add = rhs[c, k, j, i, m];
                                rms_work[m] += add * add;
                            }
                        }
                    }
                }
            }

            //       mpi_allreduce(rms_work, rms, 5, dp_type, MPI_SUM, comm_setup, error)
            comm_setup.Allreduce<double>(rms_work, Operation<double>.Add, ref rms);


            for (m = 0; m < 5; m++)
            {
                for (d = 0; d < 3; d++)
                {
                    rms[m] = rms[m] / (grid_points[d] - 2);
                }
                rms[m] = Math.Sqrt(rms[m]);
            }
        }

        public void exact_rhs()
        {
            double[] dtemp = new double[5];
            double xi, eta, zeta, dtpp;
            int c, m, i, j, k, ip1, im1, jp1, jm1, km1, kp1, ksize, jsize, isize;

            // int ii, jj, kk;  /* +2 offset required by C# arrays for Fortran with -2 lower limit */


            for (c = 0; c < ncells; c++)
            {
                ksize = cell_size[c, 2] + 2;
                jsize = cell_size[c, 1] + 2;
                isize = cell_size[c, 0] + 2;

                //---------------------------------------------------------------------
                //      initialize                                  
                //---------------------------------------------------------------------
                for (m = 0; m < 5; m++)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                forcing[c, k, j, i, m] = 0.0d;
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      xi-direction flux differences                      
                //---------------------------------------------------------------------
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    zeta = (k - 2 + cell_low[c, 2]) * dnzm1;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        eta = (j - 2 + cell_low[c, 1]) * dnym1;
                        for (i = 2 * start[c, 0] - 4; i <= isize + 1 - 2 * end[c, 0]; i++)
                        {
                            xi = (i - 2 + cell_low[c, 0]) * dnxm1;

                            exact_solution(xi, eta, zeta, dtemp, 0);
                            for (m = 0; m < 5; m++)
                            {
                                ue[i, m] = dtemp[m]; // OK ue[i,m]
                            }

                            dtpp = 1.0d / dtemp[0];

                            for (m = 1; m < 5; m++)
                            {
                                buf[i, m] = dtpp * dtemp[m]; // OK buf[i,m]
                            }

                            cuf[i] = buf[i, 1] * buf[i, 1];
                            buf[i, 0] = cuf[i] + buf[i, 2] * buf[i, 2] +
                                                 buf[i, 3] * buf[i, 3];
                            q[i] = 0.5d * (buf[i, 1] * ue[i, 1] + buf[i, 2] * ue[i, 2] +
                                                                 buf[i, 3] * ue[i, 3]);
                            // OK cuf[i], buf[i,0], q[i]
                        }

                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            im1 = i - 1;
                            ip1 = i + 1;

                            double b = forcing[c, k, j, i, 0];

                            forcing[c, k, j, i, 0] = forcing[c, k, j, i, 0] -
                                             tx2 * (ue[ip1, 1] - ue[im1, 1]) +
                                             dx1tx1 * (ue[ip1, 0] - 2.0d * ue[i, 0] + ue[im1, 0]);

                            forcing[c, k, j, i, 1] = forcing[c, k, j, i, 1] - tx2 * (
                                            (ue[ip1, 1] * buf[ip1, 1] + c2 * (ue[ip1, 4] - q[ip1])) -
                                            (ue[im1, 1] * buf[im1, 1] + c2 * (ue[im1, 4] - q[im1]))) +
                                             xxcon1 * (buf[ip1, 1] - 2.0d * buf[i, 1] + buf[im1, 1]) +
                                             dx2tx1 * (ue[ip1, 1] - 2.0d * ue[i, 1] + ue[im1, 1]);

                            forcing[c, k, j, i, 2] = forcing[c, k, j, i, 2] - tx2 * (
                                             ue[ip1, 2] * buf[ip1, 1] - ue[im1, 2] * buf[im1, 1]) +
                                             xxcon2 * (buf[ip1, 2] - 2.0d * buf[i, 2] + buf[im1, 2]) +
                                             dx3tx1 * (ue[ip1, 2] - 2.0d * ue[i, 2] + ue[im1, 2]);


                            forcing[c, k, j, i, 3] = forcing[c, k, j, i, 3] - tx2 * (
                                             ue[ip1, 3] * buf[ip1, 1] - ue[im1, 3] * buf[im1, 1]) +
                                             xxcon2 * (buf[ip1, 3] - 2.0d * buf[i, 3] + buf[im1, 3]) +
                                             dx4tx1 * (ue[ip1, 3] - 2.0d * ue[i, 3] + ue[im1, 3]);

                            forcing[c, k, j, i, 4] = forcing[c, k, j, i, 4] - tx2 * (
                                             buf[ip1, 1] * (c1 * ue[ip1, 4] - c2 * q[ip1]) -
                                             buf[im1, 1] * (c1 * ue[im1, 4] - c2 * q[im1])) +
                                             0.5d * xxcon3 * (buf[ip1, 0] - 2.0d * buf[i, 0] +
                                                           buf[im1, 0]) +
                                             xxcon4 * (cuf[ip1] - 2.0d * cuf[i] + cuf[im1]) +
                                             xxcon5 * (buf[ip1, 4] - 2.0d * buf[i, 4] + buf[im1, 4]) +
                                             dx5tx1 * (ue[ip1, 4] - 2.0d * ue[i, 4] + ue[im1, 4]);
                        }

                        //---------------------------------------------------------------------
                        //            Fourth-order dissipation                         
                        //---------------------------------------------------------------------
                        if (start[c, 0] > 2)
                        {
                            
                            for (m = 0; m < 5; m++)
                            {
                                i = 3;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                                    (5.0d * ue[i, m] - 4.0d * ue[i + 1, m] + ue[i + 2, m]);
                                i = 4;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                                   (-4.0d * ue[i - 1, m] + 6.0d * ue[i, m] -
                                                     4.0d * ue[i + 1, m] + ue[i + 2, m]);
                            }
                        }

                        for (m = 0; m < 5; m++)
                        {
                            for (i = 3 * start[c, 0] - 4; i < isize - 3 * end[c, 0]; i++)
                            {
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                                 (ue[i - 2, m] - 4.0d * ue[i - 1, m] +
                                                  6.0d * ue[i, m] - 4.0d * ue[i + 1, m] + ue[i + 2, m]);
                            }
                        }

                        if (end[c, 0] > 0)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                i = isize - 3;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                                   (ue[i - 2, m] - 4.0d * ue[i - 1, m] +
                                                    6.0d * ue[i, m] - 4.0d * ue[i + 1, m]);
                                i = isize - 2;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                                   (ue[i - 2, m] - 4.0d * ue[i - 1, m] + 5.0d * ue[i, m]);
                            }
                        }
                    }
                } // end k

                //---------------------------------------------------------------------
                //  eta-direction flux differences             
                //---------------------------------------------------------------------
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    zeta = (k - 2 + cell_low[c, 2]) * dnzm1;
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        xi = (i - 2 + cell_low[c, 0]) * dnxm1;

                        for (j = 2 * start[c, 1] - 4; j <= jsize + 1 - 2 * end[c, 1]; j++)
                        {
                            eta = (j - 2 + cell_low[c, 1]) * dnym1;

                            exact_solution(xi, eta, zeta, dtemp, 0);
                            for (m = 0; m < 5; m++)
                            {
                                ue[j, m] = dtemp[m];
                            }
                            dtpp = 1.0d / dtemp[0];

                            for (m = 1; m < 5; m++)
                            {
                                buf[j, m] = dtpp * dtemp[m];
                            }

                            cuf[j] = buf[j, 2] * buf[j, 2];
                            buf[j, 0] = cuf[j] + buf[j, 1] * buf[j, 1] +
                                                 buf[j, 3] * buf[j, 3];
                            q[j] = 0.5d * (buf[j, 1] * ue[j, 1] +
                                          buf[j, 2] * ue[j, 2] +
                                          buf[j, 3] * ue[j, 3]);
                        }

                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            jm1 = j - 1;
                            jp1 = j + 1;

                            forcing[c, k, j, i, 0] = forcing[c, k, j, i, 0] -
                                  ty2 * (ue[jp1, 2] - ue[jm1, 2]) +
                                  dy1ty1 * (ue[jp1, 0] - 2.0d * ue[j, 0] + ue[jm1, 0]);

                            forcing[c, k, j, i, 1] = forcing[c, k, j, i, 1] - ty2 * (
                                  ue[jp1, 1] * buf[jp1, 2] - ue[jm1, 1] * buf[jm1, 2]) +
                                  yycon2 * (buf[jp1, 1] - 2.0d * buf[j, 1] + buf[jm1, 1]) +
                                  dy2ty1 * (ue[jp1, 1] - 2.0d * ue[j, 1] + ue[jm1, 1]);

                            forcing[c, k, j, i, 2] = forcing[c, k, j, i, 2] - ty2 * (
                                  (ue[jp1, 2] * buf[jp1, 2] + c2 * (ue[jp1, 4] - q[jp1])) -
                                  (ue[jm1, 2] * buf[jm1, 2] + c2 * (ue[jm1, 4] - q[jm1]))) +
                                  yycon1 * (buf[jp1, 2] - 2.0d * buf[j, 2] + buf[jm1, 2]) +
                                  dy3ty1 * (ue[jp1, 2] - 2.0d * ue[j, 2] + ue[jm1, 2]);

                            forcing[c, k, j, i, 3] = forcing[c, k, j, i, 3] - ty2 * (
                                  ue[jp1, 3] * buf[jp1, 2] - ue[jm1, 3] * buf[jm1, 2]) +
                                  yycon2 * (buf[jp1, 3] - 2.0d * buf[j, 3] + buf[jm1, 3]) +
                                  dy4ty1 * (ue[jp1, 3] - 2.0d * ue[j, 3] + ue[jm1, 3]);

                            forcing[c, k, j, i, 4] = forcing[c, k, j, i, 4] - ty2 * (
                                  buf[jp1, 2] * (c1 * ue[jp1, 4] - c2 * q[jp1]) -
                                  buf[jm1, 2] * (c1 * ue[jm1, 4] - c2 * q[jm1])) +
                                  0.5d * yycon3 * (buf[jp1, 0] - 2.0d * buf[j, 0] +
                                                buf[jm1, 0]) +
                                  yycon4 * (cuf[jp1] - 2.0d * cuf[j] + cuf[jm1]) +
                                  yycon5 * (buf[jp1, 4] - 2.0d * buf[j, 4] + buf[jm1, 4]) +
                                  dy5ty1 * (ue[jp1, 4] - 2.0d * ue[j, 4] + ue[jm1, 4]);
                        }

                        //---------------------------------------------------------------------
                        //            Fourth-order dissipation                      
                        //---------------------------------------------------------------------
                        if (start[c, 1] > 2)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                j = 3;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                          (5.0d * ue[j, m] - 4.0d * ue[j + 1, m] + ue[j + 2, m]);
                                j = 4;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                         (-4.0d * ue[j - 1, m] + 6.0d * ue[j, m] -
                                           4.0d * ue[j + 1, m] + ue[j + 2, m]);
                            }
                        }

                        for (m = 0; m < 5; m++)
                        {
                            for (j = 3 * start[c, 1] - 4; j < jsize - 3 * end[c, 1]; j++)
                            {
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                      (ue[j - 2, m] - 4.0d * ue[j - 1, m] +
                                       6.0d * ue[j, m] - 4.0d * ue[j + 1, m] + ue[j + 2, m]);
                            }
                        }

                        if (end[c, 1] > 0)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                j = jsize - 3;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                         (ue[j - 2, m] - 4.0d * ue[j - 1, m] +
                                          6.0d * ue[j, m] - 4.0d * ue[j + 1, m]);
                                j = jsize - 2;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                         (ue[j - 2, m] - 4.0d * ue[j - 1, m] + 5.0d * ue[j, m]);

                            }
                        }

                    }
                }

                //---------------------------------------------------------------------
                //      zeta-direction flux differences                      
                //---------------------------------------------------------------------
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    eta = (j - 2 + cell_low[c, 1]) * dnym1;
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        xi = (i - 2 + cell_low[c, 0]) * dnxm1;

                        for (k = 2 * start[c, 2] - 4; k <= ksize + 1 - 2 * end[c, 2]; k++)
                        {
                            zeta = (k - 2 + cell_low[c, 2]) * dnzm1;

                            exact_solution(xi, eta, zeta, dtemp, 0);
                            for (m = 0; m <= 4; m++)
                            {
                                ue[k, m] = dtemp[m];
                            }

                            dtpp = 1.0d / dtemp[0];

                            for (m = 1; m < 5; m++)
                            {
                                buf[k, m] = dtpp * dtemp[m];
                            }

                            cuf[k] = buf[k, 3] * buf[k, 3];
                            buf[k, 0] = cuf[k] + buf[k, 1] * buf[k, 1] +
                                       buf[k, 2] * buf[k, 2];
                            q[k] = 0.5d * (buf[k, 1] * ue[k, 1] + buf[k, 2] * ue[k, 2] +
                                          buf[k, 3] * ue[k, 3]);
                        }

                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            km1 = k - 1;
                            kp1 = k + 1;

                            forcing[c, k, j, i, 0] = forcing[c, k, j, i, 0] -
                                   tz2 * (ue[kp1, 3] - ue[km1, 3]) +
                                   dz1tz1 * (ue[kp1, 0] - 2.0d * ue[k, 0] + ue[km1, 0]);

                            forcing[c, k, j, i, 1] = forcing[c, k, j, i, 1] - tz2 * (
                                   ue[kp1, 1] * buf[kp1, 3] - ue[km1, 1] * buf[km1, 3]) +
                                   zzcon2 * (buf[kp1, 1] - 2.0d * buf[k, 1] + buf[km1, 1]) +
                                   dz2tz1 * (ue[kp1, 1] - 2.0d * ue[k, 1] + ue[km1, 1]);

                            forcing[c, k, j, i, 2] = forcing[c, k, j, i, 2] - tz2 * (
                                   ue[kp1, 2] * buf[kp1, 3] - ue[km1, 2] * buf[km1, 3]) +
                                   zzcon2 * (buf[kp1, 2] - 2.0d * buf[k, 2] + buf[km1, 2]) +
                                   dz3tz1 * (ue[kp1, 2] - 2.0d * ue[k, 2] + ue[km1, 2]);

                            forcing[c, k, j, i, 3] = forcing[c, k, j, i, 3] - tz2 * (
                                  (ue[kp1, 3] * buf[kp1, 3] + c2 * (ue[kp1, 4] - q[kp1])) -
                                  (ue[km1, 3] * buf[km1, 3] + c2 * (ue[km1, 4] - q[km1]))) +
                                  zzcon1 * (buf[kp1, 3] - 2.0d * buf[k, 3] + buf[km1, 3]) +
                                  dz4tz1 * (ue[kp1, 3] - 2.0d * ue[k, 3] + ue[km1, 3]);

                            forcing[c, k, j, i, 4] = forcing[c, k, j, i, 4] - tz2 * (
                                   buf[kp1, 3] * (c1 * ue[kp1, 4] - c2 * q[kp1]) -
                                   buf[km1, 3] * (c1 * ue[km1, 4] - c2 * q[km1])) +
                                   0.5d * zzcon3 * (buf[kp1, 0] - 2.0d * buf[k, 0]
                                                + buf[km1, 0]) +
                                   zzcon4 * (cuf[kp1] - 2.0d * cuf[k] + cuf[km1]) +
                                   zzcon5 * (buf[kp1, 4] - 2.0d * buf[k, 4] + buf[km1, 4]) +
                                   dz5tz1 * (ue[kp1, 4] - 2.0d * ue[k, 4] + ue[km1, 4]);
                        }

                        //---------------------------------------------------------------------
                        //            Fourth-order dissipation
                        //---------------------------------------------------------------------
                        if (start[c, 2] > 2)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                k = 3;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                          (5.0d * ue[k, m] - 4.0d * ue[k + 1, m] + ue[k + 2, m]);
                                k = 4;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                         (-4.0d * ue[k - 1, m] + 6.0d * ue[k, m] -
                                           4.0d * ue[k + 1, m] + ue[k + 2, m]);
                            }
                        }

                        for (m = 0; m < 5; m++)
                        {
                            for (k = 3 * start[c, 2] - 4; k < ksize - 3 * end[c, 2]; k++)
                            {
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                      (ue[k - 2, m] - 4.0d * ue[k - 1, m] +
                                       6.0d * ue[k, m] - 4.0d * ue[k + 1, m] + ue[k + 2, m]);
                            }
                        }

                        if (end[c, 2] > 0)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                k = ksize - 3;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                         (ue[k - 2, m] - 4.0d * ue[k - 1, m] +
                                          6.0d * ue[k, m] - 4.0d * ue[k + 1, m]);
                                k = ksize - 2;
                                forcing[c, k, j, i, m] = forcing[c, k, j, i, m] - dssp *
                                      (ue[k - 2, m] - 4.0d * ue[k - 1, m] + 5.0d * ue[k, m]);
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                // now change the sign of the forcing function, 
                //---------------------------------------------------------------------
                for (m = 0; m < 5; m++)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                forcing[c, k, j, i, m] = -1.0d * forcing[c, k, j, i, m];
                            }
                        }
                    }
                }
            } // cell loop
        }

        public void copy_faces()
        {
            int i, j, k, c, m, p0, p1, p2, p3, p4, p5, ksize, jsize, isize;
            Request[] requests;
            int[] b_size;

            double[][] out_buffer = new double[6][];
            double[][] in_buffer = new double[6][];

            requests = new Request[12];
            b_size = new int[6];

            if (no_nodes == 1)
            {
                compute_rhs();
                return;
            }

            b_size[0] = east_size;
            b_size[1] = west_size;
            b_size[2] = north_size;
            b_size[3] = south_size;
            b_size[4] = top_size;
            b_size[5] = bottom_size;

            for (i = 0; i < 6; i++)
            {
                out_buffer[i] = new double[b_size[i]];
                in_buffer[i] = new double[b_size[i]];
            }


            //---------------------------------------------------------------------
            // because the difference stencil for the diagonalized scheme is 
            // orthogonal, we do not have to perform the staged copying of faces, 
            // but can send all face information simultaneously to the neighboring 
            // cells in all directions          
            //---------------------------------------------------------------------
            p0 = 0;
            p1 = 0;
            p2 = 0;
            p3 = 0;
            p4 = 0;
            p5 = 0;

            for (c = 0; c < ncells; c++)
            {
                ksize = cell_size[c, 2] + 2;
                jsize = cell_size[c, 1] + 2;
                isize = cell_size[c, 0] + 2;

                //---------------------------------------------------------------------
                //            fill the buffer to be sent to eastern neighbors (i-dir)
                //---------------------------------------------------------------------
                if (cell_coord[c, 0] != ncells - 1)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = isize - 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    out_buffer[0][p0++] = u[c, k, j, i, m];
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //            fill the buffer to be sent to western neighbors 
                //---------------------------------------------------------------------
                if (cell_coord[c, 0] != 0)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = 2; i <= 3; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    out_buffer[1][p1++] = u[c, k, j, i, m];
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //            fill the buffer to be sent to northern neighbors (j_dir)
                //---------------------------------------------------------------------
                if (cell_coord[c, 1] != ncells - 1)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = jsize - 2; j < jsize; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    out_buffer[2][p2++] = u[c, k, j, i, m];
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //            fill the buffer to be sent to southern neighbors 
                //---------------------------------------------------------------------
                if (cell_coord[c, 1] != 0)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = 2; j <= 3; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    out_buffer[3][p3++] = u[c, k, j, i, m];
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //            fill the buffer to be sent to top neighbors (k-dir)
                //---------------------------------------------------------------------
                if (cell_coord[c, 2] != ncells - 1)
                {
                    for (k = ksize - 2; k < ksize; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    out_buffer[4][p4++] = u[c, k, j, i, m];
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //            fill the buffer to be sent to bottom neighbors
                //---------------------------------------------------------------------
                if (cell_coord[c, 2] != 0)
                {
                    for (k = 2; k <= 3; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    out_buffer[5][p5++] = u[c, k, j, i, m];
                                }
                            }
                        }
                    }
                }
            }


            RequestList requestList = new RequestList();

            requests[0] = comm_rhs.ImmediateReceive<double>(successor[0], WEST, in_buffer[0]);
            requests[1] = comm_rhs.ImmediateReceive<double>(predecessor[0], EAST, in_buffer[1]);
            requests[2] = comm_rhs.ImmediateReceive<double>(successor[1], SOUTH, in_buffer[2]);
            requests[3] = comm_rhs.ImmediateReceive<double>(predecessor[1], NORTH, in_buffer[3]);
            requests[4] = comm_rhs.ImmediateReceive<double>(successor[2], BOTTOM, in_buffer[4]);
            requests[5] = comm_rhs.ImmediateReceive<double>(predecessor[2], TOP, in_buffer[5]);
            requests[6] = comm_rhs.ImmediateSend<double>(out_buffer[0], successor[0], EAST);
            requests[7] = comm_rhs.ImmediateSend<double>(out_buffer[1], predecessor[0], WEST);
            requests[8] = comm_rhs.ImmediateSend<double>(out_buffer[2], successor[1], NORTH);
            requests[9] = comm_rhs.ImmediateSend<double>(out_buffer[3], predecessor[1], SOUTH);
            requests[10] = comm_rhs.ImmediateSend<double>(out_buffer[4], successor[2], TOP);
            requests[11] = comm_rhs.ImmediateSend<double>(out_buffer[5], predecessor[2], BOTTOM);

            foreach (Request request in requests)
            {
                requestList.Add(request);
            }

            requestList.WaitAll();

            //---------------------------------------------------------------------
            // unpack the data that has just been received;             
            //---------------------------------------------------------------------
            p0 = 0;
            p1 = 0;
            p2 = 0;
            p3 = 0;
            p4 = 0;
            p5 = 0;

            for (c = 0; c < ncells; c++)
            {
                isize = cell_size[c, 0] + 2;
                jsize = cell_size[c, 1] + 2;
                ksize = cell_size[c, 2] + 2;

                if (cell_coord[c, 0] != 0)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = 0; i <= 1; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    u[c, k, j, i, m] = in_buffer[1][p0++];
                                }
                            }
                        }
                    }
                }

                if (cell_coord[c, 0] != ncells - 1)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = isize; i <= isize + 1; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    u[c, k, j, i, m] = in_buffer[0][p1++];
                                }
                            }
                        }
                    }
                }

                if (cell_coord[c, 1] != 0)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = 0; j <= 1; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    u[c, k, j, i, m] = in_buffer[3][p2++];
                                }
                            }
                        }
                    }
                }

                if (cell_coord[c, 1] != ncells - 1)
                {
                    for (k = 2; k < ksize; k++)
                    {
                        for (j = jsize; j <= jsize + 1; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    u[c, k, j, i, m] = in_buffer[2][p3++];
                                }
                            }
                        }
                    }
                }

                if (cell_coord[c, 2] != 0)
                {
                    for (k = 0; k <= 1; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    u[c, k, j, i, m] = in_buffer[5][p4++];
                                }
                            }
                        }
                    }
                }

                if (cell_coord[c, 2] != ncells - 1)
                {
                    for (k = ksize; k <= ksize + 1; k++)
                    {
                        for (j = 2; j < jsize; j++)
                        {
                            for (i = 2; i < isize; i++)
                            {
                                for (m = 0; m < 5; m++)
                                {
                                    u[c, k, j, i, m] = in_buffer[4][p5++];
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      cells loop
                //---------------------------------------------------------------------
            }

            //---------------------------------------------------------------------
            // now that we have all the data, compute the rhs
            //---------------------------------------------------------------------
            compute_rhs();
        }

        //---------------------------------------------------------------------
        // This function allocates space for a set of cells and fills the set     
        // such that communication between cells on different nodes is only
        // nearest neighbor                                                   
        //---------------------------------------------------------------------
        public void make_set()
        {
            int p, i, j, c, dir, size, excess, ierrcode = 0;

            //---------------------------------------------------------------------
            //     compute square root; add small number to allow for roundoff
            //     (note: this is computed in setup_mpi.f also, but prefer to do
            //     it twice because of some include file problems).
            //---------------------------------------------------------------------
            ncells = Convert.ToInt32(Math.Sqrt(no_nodes));

            //---------------------------------------------------------------------
            //      this makes coding easier
            //---------------------------------------------------------------------
            p = ncells;

            //---------------------------------------------------------------------
            //      determine the location of the cell at the bottom of the 3D 
            //      array of cells
            //---------------------------------------------------------------------
            cell_coord[0, 0] = mod(node, p);   //mod(node,p);
            cell_coord[0, 1] = node / p;
            cell_coord[0, 2] = 0;

            //---------------------------------------------------------------------
            //      set the cell_coords for cells in the rest of the z-layers; 
            //      this comes down to a simple linear numbering in the z-direct-
            //      ion, and to the doubly-cyclic numbering in the other dirs     
            //---------------------------------------------------------------------
            for (c = 1; c < p; c++)
            {
                cell_coord[c, 0] = mod(cell_coord[c - 1, 0] + 1, p) ;                           // mod(cell_coord(1,c-1)+1,p);
                cell_coord[c, 1] = mod(cell_coord[c - 1, 1] - 1 + p, p) ;                           // mod(cell_coord(2,c-1)-1+p,p); 
                cell_coord[c, 2] = c ;
            }

            //---------------------------------------------------------------------
            //      slice(n,dir) contains the sequence number of the cell that is in
            //      coordinate plane n in the dir direction
            //---------------------------------------------------------------------
            for (dir = 0; dir < 3; dir++)
            {
                for (c = 0; c < p; c++)
                {
                    slice[cell_coord[c, dir], dir] = c;
                }
            }
            
            //---------------------------------------------------------------------
            //      fill the predecessor and successor entries, using the indices 
            //      of the bottom cells (they are the same at each level of k 
            //      anyway) acting as if full periodicity pertains; note that p is
            //      added to those arguments to the mod functions that might
            //      otherwise return wrong values when using the modulo function
            //---------------------------------------------------------------------
            i = cell_coord[0, 0];
            j = cell_coord[0, 1];

            predecessor[0] = mod(i - 1 + p, p) + p * j;
            predecessor[1] = i + p * mod(j - 1 + p, p);
            predecessor[2] = mod(i + 1, p) + p * mod(j - 1 + p, p);

            successor[0] = mod(i + 1, p) + p * j;
            successor[1] = i + p * mod(j + 1, p);
            successor[2] = mod(i - 1 + p, p) + p * mod(j + 1, p);

            //---------------------------------------------------------------------
            // now compute the sizes of the cells                                    
            //---------------------------------------------------------------------
            for (dir = 0; dir < 3; dir++)
            {
                //---------------------------------------------------------------------
                //         set cell_coord range for each direction                            
                //---------------------------------------------------------------------
                size = grid_points[dir] / p;
                excess = mod(grid_points[dir], p);

                for (c = 0; c < ncells; c++)
                {
                    if (cell_coord[c, dir] < excess)
                    {
                        cell_size[c, dir] = size + 1;
                        cell_low[c, dir] = (cell_coord[c, dir]) * (size + 1);
                        cell_high[c, dir] = cell_low[c, dir] + size;
                    }
                    else
                    {
                        cell_size[c, dir] = size;
                        cell_low[c, dir] = excess * (size + 1) + (cell_coord[c, dir] - excess) * size;
                        cell_high[c, dir] = cell_low[c, dir] + size - 1;
                    }
                    if (cell_size[c, dir] <= 2)
                    {
                        Console.WriteLine("Error: Cell size too small. Min size is 3");
                        worldcomm.Abort(ierrcode); 
                        System.Environment.Exit(0);
                    }
                }
            }
        }

        public int mod(int a, int b)
        {
            int r;
            Math.DivRem(a, b, out r);
            return r;
        }

        public void compute_rhs()
        {
            int c, i, j, k, m, ksize, jsize, isize;
            double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1,
                   wijk, wp1, wm1;

            //---------------------------------------------------------------------
            // loop over all cells owned by this node                           
            //---------------------------------------------------------------------
            for (c = 0; c < ncells; c++)
            {
                ksize = cell_size[c, 2] + 2;
                jsize = cell_size[c, 1] + 2;
                isize = cell_size[c, 0] + 2;

                //---------------------------------------------------------------------
                //      compute the reciprocal of density, and the kinetic energy, 
                //      and the speed of sound. 
                //---------------------------------------------------------------------

                for (k = 1; k <= ksize; k++)
                {
                    for (j = 1; j <= jsize; j++)
                    {
                        for (i = 1; i <= isize; i++)
                        {
                            rho_inv = 1.0d / u[c, k, j, i, 0];
                            rho_i[c, k, j, i] = rho_inv;
                            us[c, k, j, i] = u[c, k, j, i, 1] * rho_inv;
                            vs[c, k, j, i] = u[c, k, j, i, 2] * rho_inv;
                            ws[c, k, j, i] = u[c, k, j, i, 3] * rho_inv;
                            square[c, k, j, i] = 0.5d * (
                                          u[c, k, j, i, 1] * u[c, k, j, i, 1] +
                                          u[c, k, j, i, 2] * u[c, k, j, i, 2] +
                                          u[c, k, j, i, 3] * u[c, k, j, i, 3]) * rho_inv;
                            qs[c, k, j, i] = square[c, k, j, i] * rho_inv;
                            //---------------------------------------------------------------------
                            //               (don't need speed and ainx until the lhs computation)
                            //---------------------------------------------------------------------
                            aux = c1c2 * rho_inv * (u[c, k, j, i, 4] - square[c, k, j, i]);
                            aux = Math.Sqrt(aux);
                            speed[c, k, j, i] = aux;
                            ainv[c, k, j, i] = 1.0d / aux;
                        }
                    }
                }

                //---------------------------------------------------------------------
                // copy the exact forcing term to the right hand side;  because 
                // this forcing term is known, we can store it on the whole grid
                // including the boundary                   
                //---------------------------------------------------------------------

                for (k = 2; k < ksize; k++)
                {
                    for (j = 2; j < jsize; j++)
                    {
                        for (i = 2; i < isize; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = forcing[c, k, j, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      compute xi-direction fluxes 
                //---------------------------------------------------------------------
                //if (timeron) timer.start(t_rhsx);
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            uijk = us[c, k, j, i];
                            up1 = us[c, k, j, i + 1];
                            um1 = us[c, k, j, i - 1];

                            rhs[c, k, j, i, 0] = rhs[c, k, j, i, 0] + dx1tx1 *
                                      (u[c, k, j, i + 1, 0] - 2.0d * u[c, k, j, i, 0] +
                                       u[c, k, j, i - 1, 0]) -
                                      tx2 * (u[c, k, j, i + 1, 1] - u[c, k, j, i - 1, 1]);

                            rhs[c, k, j, i, 1] = rhs[c, k, j, i, 1] + dx2tx1 *
                                      (u[c, k, j, i + 1, 1] - 2.0d * u[c, k, j, i, 1] +
                                       u[c, k, j, i - 1, 1]) +
                                      xxcon2 * con43 * (up1 - 2.0d * uijk + um1) -
                                      tx2 * (u[c, k, j, i + 1, 1] * up1 -
                                             u[c, k, j, i - 1, 1] * um1 +
                                             (u[c, k, j, i + 1, 4] - square[c, k, j, i + 1] -
                                              u[c, k, j, i - 1, 4] + square[c, k, j, i - 1]) *
                                              c2);

                            rhs[c, k, j, i, 2] = rhs[c, k, j, i, 2] + dx3tx1 *
                                      (u[c, k, j, i + 1, 2] - 2.0d * u[c, k, j, i, 2] +
                                       u[c, k, j, i - 1, 2]) +
                                      xxcon2 * (vs[c, k, j, i + 1] - 2.0d * vs[c, k, j, i] +
                                                vs[c, k, j, i - 1]) -
                                      tx2 * (u[c, k, j, i + 1, 2] * up1 -
                                             u[c, k, j, i - 1, 2] * um1);

                            rhs[c, k, j, i, 3] = rhs[c, k, j, i, 3] + dx4tx1 *
                                      (u[c, k, j, i + 1, 3] - 2.0d * u[c, k, j, i, 3] +
                                       u[c, k, j, i - 1, 3]) +
                                      xxcon2 * (ws[c, k, j, i + 1] - 2.0d * ws[c, k, j, i] +
                                                ws[c, k, j, i - 1]) -
                                      tx2 * (u[c, k, j, i + 1, 3] * up1 -
                                             u[c, k, j, i - 1, 3] * um1);

                            rhs[c, k, j, i, 4] = rhs[c, k, j, i, 4] + dx5tx1 *
                                      (u[c, k, j, i + 1, 4] - 2.0d * u[c, k, j, i, 4] +
                                       u[c, k, j, i - 1, 4]) +
                                      xxcon3 * (qs[c, k, j, i + 1] - 2.0d * qs[c, k, j, i] +
                                                qs[c, k, j, i - 1]) +
                                      xxcon4 * (up1 * up1 - 2.0d * uijk * uijk +
                                                um1 * um1) +
                                      xxcon5 * (u[c, k, j, i + 1, 4] * rho_i[c, k, j, i + 1] -
                                                2.0d * u[c, k, j, i, 4] * rho_i[c, k, j, i] +
                                                u[c, k, j, i - 1, 4] * rho_i[c, k, j, i - 1]) -
                                      tx2 * ((c1 * u[c, k, j, i + 1, 4] -
                                               c2 * square[c, k, j, i + 1]) * up1 -
                                              (c1 * u[c, k, j, i - 1, 4] -
                                               c2 * square[c, k, j, i - 1]) * um1);
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      add fourth order xi-direction dissipation               
                //---------------------------------------------------------------------

                if (start[c, 0] > 2)
                {
                    i = 3;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                              (5.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j, i + 1, m] +
                                      u[c, k, j, i + 2, m]);
                            }
                        }
                    }

                    i = 4;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                              (-4.0d * u[c, k, j, i - 1, m] + 6.0d * u[c, k, j, i, m] -
                                4.0d * u[c, k, j, i + 1, m] + u[c, k, j, i + 2, m]);
                            }
                        }
                    }

                }


                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = 3 * start[c, 0] - 4; i < isize - 3 * end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                   (u[c, k, j, i - 2, m] - 4.0d * u[c, k, j, i - 1, m] +
                                    6.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j, i + 1, m] +
                                        u[c, k, j, i + 2, m]);
                            }
                        }
                    }
                }

                if (end[c, 0] > 0)
                {
                    i = ksize - 3;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j, i - 2, m] - 4.0d * u[c, k, j, i - 1, m] +
                                        6.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j, i + 1, m]);
                            }
                        }
                    }

                    i = ksize - 2;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j, i - 2, m] - 4.0d * u[c, k, j, i - 1, m] +
                                        5.0d * u[c, k, j, i, m]);
                            }
                        }
                    }

                }

                // if (timeron) timer.stop(t_rhsx);

                //---------------------------------------------------------------------
                //      compute eta-direction fluxes 
                //---------------------------------------------------------------------
                // if (timeron) timer.start(t_rhsy);

                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            vijk = vs[c, k, j, i];
                            vp1 = vs[c, k, j + 1, i];
                            vm1 = vs[c, k, j - 1, i];
                            rhs[c, k, j, i, 0] = rhs[c, k, j, i, 0] + dy1ty1 *
                                     (u[c, k, j + 1, i, 0] - 2.0d * u[c, k, j, i, 0] +
                                      u[c, k, j - 1, i, 0]) -
                                     ty2 * (u[c, k, j + 1, i, 2] - u[c, k, j - 1, i, 2]);
                            rhs[c, k, j, i, 1] = rhs[c, k, j, i, 1] + dy2ty1 *
                                     (u[c, k, j + 1, i, 1] - 2.0d * u[c, k, j, i, 1] +
                                      u[c, k, j - 1, i, 1]) +
                                     yycon2 * (us[c, k, j + 1, i] - 2.0d * us[c, k, j, i] +
                                               us[c, k, j - 1, i]) -
                                     ty2 * (u[c, k, j + 1, i, 1] * vp1 -
                                            u[c, k, j - 1, i, 1] * vm1);
                            rhs[c, k, j, i, 2] = rhs[c, k, j, i, 2] + dy3ty1 *
                                     (u[c, k, j + 1, i, 2] - 2.0d * u[c, k, j, i, 2] +
                                      u[c, k, j - 1, i, 2]) +
                                     yycon2 * con43 * (vp1 - 2.0d * vijk + vm1) -
                                     ty2 * (u[c, k, j + 1, i, 2] * vp1 -
                                            u[c, k, j - 1, i, 2] * vm1 +
                                            (u[c, k, j + 1, i, 4] - square[c, k, j + 1, i] -
                                             u[c, k, j - 1, i, 4] + square[c, k, j - 1, i])
                                            * c2);
                            rhs[c, k, j, i, 3] = rhs[c, k, j, i, 3] + dy4ty1 *
                                     (u[c, k, j + 1, i, 3] - 2.0d * u[c, k, j, i, 3] +
                                      u[c, k, j - 1, i, 3]) +
                                     yycon2 * (ws[c, k, j + 1, i] - 2.0d * ws[c, k, j, i] +
                                               ws[c, k, j - 1, i]) -
                                     ty2 * (u[c, k, j + 1, i, 3] * vp1 -
                                            u[c, k, j - 1, i, 3] * vm1);
                            rhs[c, k, j, i, 4] = rhs[c, k, j, i, 4] + dy5ty1 *
                                     (u[c, k, j + 1, i, 4] - 2.0d * u[c, k, j, i, 4] +
                                      u[c, k, j - 1, i, 4]) +
                                     yycon3 * (qs[c, k, j + 1, i] - 2.0d * qs[c, k, j, i] +
                                               qs[c, k, j - 1, i]) +
                                     yycon4 * (vp1 * vp1 - 2.0d * vijk * vijk +
                                               vm1 * vm1) +
                                     yycon5 * (u[c, k, j + 1, i, 4] * rho_i[c, k, j + 1, i] -
                                               2.0d * u[c, k, j, i, 4] * rho_i[c, k, j, i] +
                                               u[c, k, j - 1, i, 4] * rho_i[c, k, j - 1, i]) -
                                     ty2 * ((c1 * u[c, k, j + 1, i, 4] -
                                             c2 * square[c, k, j + 1, i]) * vp1 -
                                            (c1 * u[c, k, j - 1, i, 4] -
                                             c2 * square[c, k, j - 1, i]) * vm1);
                        }
                    }
                }


                //---------------------------------------------------------------------
                //      add fourth order eta-direction dissipation         
                //---------------------------------------------------------------------

                if (start[c, 1] > 2)
                {
                    j = 3;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (5.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j + 1, i, m] +
                                              u[c, k, j + 2, i, m]);
                            }
                        }
                    }


                    j = 4;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (-4.0d * u[c, k, j - 1, i, m] + 6.0d * u[c, k, j, i, m] -
                                        4.0d * u[c, k, j + 1, i, m] + u[c, k, j + 2, i, m]);
                            }
                        }
                    }
                }

                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = 3 * start[c, 1] - 4; j < jsize - 3 * end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                   (u[c, k, j - 2, i, m] - 4.0d * u[c, k, j - 1, i, m] +
                                    6.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j + 1, i, m] +
                                        u[c, k, j + 2, i, m]);
                            }
                        }
                    }
                }


                if (end[c, 1] > 0)
                {
                    j = jsize - 3;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j - 2, i, m] - 4.0d * u[c, k, j - 1, i, m] +
                                        6.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j + 1, i, m]);
                            }
                        }
                    }


                    j = jsize - 2;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j - 2, i, m] - 4.0d * u[c, k, j - 1, i, m] +
                                        5.0d * u[c, k, j, i, m]);
                            }
                        }
                    }
                }


                // if (timeron) timer.stop(t_rhsy);

                //---------------------------------------------------------------------
                //      compute zeta-direction fluxes 
                //---------------------------------------------------------------------
                // if (timeron) timer.start(t_rhsz);
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            wijk = ws[c, k, j, i];
                            wp1 = ws[c, k + 1, j, i];
                            wm1 = ws[c, k - 1, j, i];

                            rhs[c, k, j, i, 0] = rhs[c, k, j, i, 0] + dz1tz1 *
                                     (u[c, k + 1, j, i, 0] - 2.0d * u[c, k, j, i, 0] +
                                      u[c, k - 1, j, i, 0]) -
                                     tz2 * (u[c, k + 1, j, i, 3] - u[c, k - 1, j, i, 3]);
                            rhs[c, k, j, i, 1] = rhs[c, k, j, i, 1] + dz2tz1 *
                                     (u[c, k + 1, j, i, 1] - 2.0d * u[c, k, j, i, 1] +
                                      u[c, k - 1, j, i, 1]) +
                                     zzcon2 * (us[c, k + 1, j, i] - 2.0d * us[c, k, j, i] +
                                               us[c, k - 1, j, i]) -
                                     tz2 * (u[c, k + 1, j, i, 1] * wp1 -
                                            u[c, k - 1, j, i, 1] * wm1);
                            rhs[c, k, j, i, 2] = rhs[c, k, j, i, 2] + dz3tz1 *
                                     (u[c, k + 1, j, i, 2] - 2.0d * u[c, k, j, i, 2] +
                                      u[c, k - 1, j, i, 2]) +
                                     zzcon2 * (vs[c, k + 1, j, i] - 2.0d * vs[c, k, j, i] +
                                               vs[c, k - 1, j, i]) -
                                     tz2 * (u[c, k + 1, j, i, 2] * wp1 -
                                            u[c, k - 1, j, i, 2] * wm1);
                            rhs[c, k, j, i, 3] = rhs[c, k, j, i, 3] + dz4tz1 *
                                     (u[c, k + 1, j, i, 3] - 2.0d * u[c, k, j, i, 3] +
                                      u[c, k - 1, j, i, 3]) +
                                     zzcon2 * con43 * (wp1 - 2.0d * wijk + wm1) -
                                     tz2 * (u[c, k + 1, j, i, 3] * wp1 -
                                            u[c, k - 1, j, i, 3] * wm1 +
                                            (u[c, k + 1, j, i, 4] - square[c, k + 1, j, i] -
                                             u[c, k - 1, j, i, 4] + square[c, k - 1, j, i])
                                            * c2);
                            rhs[c, k, j, i, 4] = rhs[c, k, j, i, 4] + dz5tz1 *
                                     (u[c, k + 1, j, i, 4] - 2.0d * u[c, k, j, i, 4] +
                                      u[c, k - 1, j, i, 4]) +
                                     zzcon3 * (qs[c, k + 1, j, i] - 2.0d * qs[c, k, j, i] +
                                               qs[c, k - 1, j, i]) +
                                     zzcon4 * (wp1 * wp1 - 2.0d * wijk * wijk +
                                               wm1 * wm1) +
                                     zzcon5 * (u[c, k + 1, j, i, 4] * rho_i[c, k + 1, j, i] -
                                               2.0 * u[c, k, j, i, 4] * rho_i[c, k, j, i] +
                                               u[c, k - 1, j, i, 4] * rho_i[c, k - 1, j, i]) -
                                     tz2 * ((c1 * u[c, k + 1, j, i, 4] -
                                              c2 * square[c, k + 1, j, i]) * wp1 -
                                             (c1 * u[c, k - 1, j, i, 4] -
                                              c2 * square[c, k - 1, j, i]) * wm1);
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      add fourth order zeta-direction dissipation                
                //---------------------------------------------------------------------

                if (start[c, 2] > 2)
                {
                    k = 3;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (5.0d * u[c, k, j, i, m] - 4.0d * u[c, k + 1, j, i, m] +
                                              u[c, k + 2, j, i, m]);
                            }
                        }
                    }

                    k = 4;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (-4.0d * u[c, k - 1, j, i, m] + 6.0d * u[c, k, j, i, m] -
                                        4.0d * u[c, k + 1, j, i, m] + u[c, k + 2, j, i, m]);
                            }
                        }
                    }
                }

                for (k = 3 * start[c, 2] - 4; k < ksize - 3 * end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                   (u[c, k - 2, j, i, m] - 4.0d * u[c, k - 1, j, i, m] +
                                    6.0d * u[c, k, j, i, m] - 4.0d * u[c, k + 1, j, i, m] +
                                        u[c, k + 2, j, i, m]);
                            }
                        }
                    }
                }

                if (end[c, 2] > 0)
                {

                    k = ksize - 3;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k - 2, j, i, m] - 4.0d * u[c, k - 1, j, i, m] +
                                        6.0d * u[c, k, j, i, m] - 4.0d * u[c, k + 1, j, i, m]);
                            }
                        }
                    }

                    k = ksize - 2;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k - 2, j, i, m] - 4.0d * u[c, k - 1, j, i, m] +
                                        5.0d * u[c, k, j, i, m]);
                            }
                        }
                    }
                }

                //if (timeron) timer.stop(t_rhsz);


                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] * dt;
                            }
                        }
                    }
                }
            }
        }

        public int verify(int no_time_steps)
        {
            double[] xcrref = new double[5], xceref = new double[5],
                     xcrdif = new double[5], xcedif = new double[5],
                     xce = new double[5], xcr = new double[5];
            double dtref = 0;
            int m;
            int verified = -1;
            char clss = 'U';
            //---------------------------------------------------------------------
            //   compute the error norm and the residual norm, and exit if not printing
            //---------------------------------------------------------------------
            error_norm(xce);
            copy_faces();
            rhs_norm(xcr);

            for (m = 0; m < 5; m++) 
                xcr[m] /= dt;

            if(node != 0) 
                return 0;

            for (m = 0; m < 5; m++)
            {
                xcrref[m] = 1.0d;
                xceref[m] = 1.0d;
            }
            //---------------------------------------------------------------------
            //    reference data for 12X12X12 grids after 100 time steps, with DT = 1.50d-02
            //---------------------------------------------------------------------
            if((grid_points[0] == 12) && (grid_points[1] == 12) && (grid_points[2] == 12) && (no_time_steps == 60)) {
                clss = 'S';
                dtref = 0.01; // 1.0d / 100;//1.0d-2;
                //c---------------------------------------------------------------------
                //c  Reference values of RMS-norms of residual.
                //c---------------------------------------------------------------------
                xcrref[0] = 1.7034283709541311E-01d;//xcrref[1] = 1.7034283709541311d-01;
                xcrref[1] = 1.2975252070034097E-02d;
                xcrref[2] = 3.2527926989486055E-02d;
                xcrref[3] = 2.6436421275166801E-02d;
                xcrref[4] = 1.9211784131744430E-01d;
                //c---------------------------------------------------------------------
                //c  Reference values of RMS-norms of solution error.
                //c---------------------------------------------------------------------
                xceref[0] = 4.9976913345811579E-04d;
                xceref[1] = 4.5195666782961927E-05d;
                xceref[2] = 7.3973765172921357E-05d;
                xceref[3] = 7.3821238632439731E-05d;
                xceref[4] = 8.9269630987491446E-04d;
                //   c---------------------------------------------------------------------
                //   c    reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
                //   c---------------------------------------------------------------------
            }
            else if((grid_points[0] == 24) && (grid_points[1] == 24) && (grid_points[2] == 24) && (no_time_steps == 200)) {
                clss = 'W';
                dtref = 0.0008;//0.8d/1000; //0.8d-3
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[0] = 0.1125590409344E+03;
                xcrref[1] = 0.1180007595731E+02;
                xcrref[2] = 0.2710329767846E+02;
                xcrref[3] = 0.2469174937669E+02;
                xcrref[4] = 0.2638427874317E+03;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                xceref[0] = 0.4419655736008E+01;
                xceref[1] = 0.4638531260002;
                xceref[2] = 0.1011551749967E+01;
                xceref[3] = 0.9235878729944;
                xceref[4] = 0.1018045837718E+02;
                //   c---------------------------------------------------------------------
                //   c    reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
                //   c---------------------------------------------------------------------
            }
            else if((grid_points[0] == 64) && (grid_points[1] == 64) && (grid_points[2] == 64) && (no_time_steps == 200)) {
                clss = 'A';
                dtref = 0.0008;//0.8d/1000; //0.8d-3
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[0] = 1.0806346714637264E+02d;
                xcrref[1] = 1.1319730901220813E+01d;
                xcrref[2] = 2.5974354511582465E+01d;
                xcrref[3] = 2.3665622544678910E+01d;
                xcrref[4] = 2.5278963211748344E+02d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                xceref[0] = 4.2348416040525025E+00d;
                xceref[1] = 4.4390282496995698E-01d;
                xceref[2] = 9.6692480136345650E-01d;
                xceref[3] = 8.8302063039765474E-01d;
                xceref[4] = 9.7379901770829278E+00d;
                //   c---------------------------------------------------------------------
                //   c    reference data for 102X102X102 grids after 200 time steps,
                //   c    with DT = 3.0d-04
                //   c---------------------------------------------------------------------
            }
            else if((grid_points[0] == 102) && (grid_points[1] == 102) && (grid_points[2] == 102) && (no_time_steps == 200)) {
                clss = 'B';
                dtref = .0003;//3.0d/10000;  //3.0d-4
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[0] = 1.4233597229287254E+03d;
                xcrref[1] = 9.9330522590150238E+01d;
                xcrref[2] = 3.5646025644535285E+02d;
                xcrref[3] = 3.2485447959084092E+02d;
                xcrref[4] = 3.2707541254659363E+03d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                xceref[0] = 5.2969847140936856E+01d;
                xceref[1] = 4.4632896115670668E+00d;
                xceref[2] = 1.3122573342210174E+01d;
                xceref[3] = 1.2006925323559144E+01d;
                xceref[4] = 1.2459576151035986E+02d;
                //   c---------------------------------------------------------------------
                //   c    reference data for 162X162X162 grids after 200 time steps,
                //   c    with DT = 1.0d-04
                //   c---------------------------------------------------------------------
            }
            else if((grid_points[0] == 162) && (grid_points[1] == 162) && (grid_points[2] == 162) && (no_time_steps == 200)) {
                clss = 'C';
                dtref = .0001;//1.0d/10000; //1.0d-4
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[0] = 0.62398116551764615E+04d;
                xcrref[1] = 0.50793239190423964E+03d;
                xcrref[2] = 0.15423530093013596E+04d;
                xcrref[3] = 0.13302387929291190E+04d;
                xcrref[4] = 0.11604087428436455E+05d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                xceref[0] = 0.16462008369091265E+03d;
                xceref[1] = 0.11497107903824313E+02d;
                xceref[2] = 0.41207446207461508E+02d;
                xceref[3] = 0.37087651059694167E+02d;
                xceref[4] = 0.36211053051841265E+03d;
                //   c---------------------------------------------------------------------
                //   c    reference data for 408x408x408 grids after 250 time steps,
                //   c    with DT = 0.2d-04
                //   c---------------------------------------------------------------------
            }
            else if((grid_points[0] == 408) && (grid_points[1] == 408) && (grid_points[2] == 408) && (no_time_steps == 250)) {
                clss = 'D';
                dtref = 0.2 / 10000; //0.2d-4
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[0] = 0.2533188551738E+05d;
                xcrref[1] = 0.2346393716980E+04d;
                xcrref[2] = 0.6294554366904E+04d;
                xcrref[3] = 0.5352565376030E+04d;
                xcrref[4] = 0.3905864038618E+05d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                xceref[0] = 0.3100009377557E+03d;
                xceref[1] = 0.2424086324913E+02d;
                xceref[2] = 0.7782212022645E+02d;
                xceref[3] = 0.6835623860116E+02d;
                xceref[4] = 0.6065737200368E+03d;
                //   c---------------------------------------------------------------------
                //   c    reference data for 1020x1020x1020 grids after 250 time steps,
                //   c    with DT = 0.4d-05
                //   c---------------------------------------------------------------------
            }
            else if((grid_points[0] == 1020) && (grid_points[1] == 1020) && (grid_points[2] == 1020) && (no_time_steps == 250)) {
                clss = 'E';
                dtref = 0.4 / 100000; //0.4d-5;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[0] = 0.9795372484517E+05d;
                xcrref[1] = 0.9739814511521E+04d;
                xcrref[2] = 0.2467606342965E+05d;
                xcrref[3] = 0.2092419572860E+05d;
                xcrref[4] = 0.1392138856939E+06d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                xceref[0] = 0.4327562208414E+03d;
                xceref[1] = 0.3699051964887E+02d;
                xceref[2] = 0.1089845040954E+03d;
                xceref[3] = 0.9462517622043E+02d;
                xceref[4] = 0.7765512765309E+03d;
            }
            //---------------------------------------------------------------------
            //    verification test for residuals if gridsize is either 12X12X12 or 
            //    64X64X64 or 102X102X102 or 162X162X162
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //    Compute the difference of solution values and the known reference values.
            //---------------------------------------------------------------------
            for (m = 0; m < 5; m++)
            {
                xcrdif[m] = Math.Abs((xcr[m] - xcrref[m]) / xcrref[m]);
                xcedif[m] = Math.Abs((xce[m] - xceref[m]) / xceref[m]);
            }
            //---------------------------------------------------------------------
            //   tolerance level
            //---------------------------------------------------------------------
            double epsilon = 1.0E-8d;
            //---------------------------------------------------------------------
            //    Output the comparison of computed results to known cases.
            //---------------------------------------------------------------------
            if (clss != 'U')
            {
                Console.WriteLine(" Verification being performed for class " + clss);
                Console.WriteLine(" Accuracy setting for epsilon = " + epsilon);
                if (Math.Abs(dt - dtref) <= epsilon)
                {
                    if (verified == -1) verified = 1;
                }
                else
                {
                    verified = 0;
                    clss = 'U';
                    Console.WriteLine("DT does not match the reference value of " + dtref);
                }
                Console.WriteLine(" Comparison of RMS-norms of residual");
            }
            else
            {
                Console.WriteLine(" Unknown CLASS");
                Console.WriteLine(" RMS-norms of residual");
            }
            verified = BMResults.printComparisonStatus(clss, verified, epsilon,
                                                     xcr, xcrref, xcrdif);
            if (clss != 'U')
            {
                Console.WriteLine(" Comparison of RMS-norms of solution error");
            }
            else
            {
                Console.WriteLine(" RMS-norms of solution error");
            }
            verified = BMResults.printComparisonStatus(clss, verified, epsilon,
                                                     xce, xceref, xcedif);

            BMResults.printVerificationStatus(clss, verified, BMName);
            return verified;
        }


        public void x_solve() {
            //---------------------------------------------------------------------
            //     
            //     Performs line solves in X direction by first factoring
            //     the block-tridiagonal matrix into an upper triangular matrix, 
            //     and { performing back substitution to solve for the unknow
            //     vectors of each line.  
            //     
            //     Make sure we treat elements zero to cell_size in the direction
            //     of the sweep.
            //     
            //---------------------------------------------------------------------
            int c, stage, first, last, buffer_size; //r_status[MPI_STATUS_SIZE];
            buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * (5 * 5 + 5);
            //MPI.Request[] recv_id = new MPI.Request[1];
            //MPI.Request[] send_id = new MPI.Request[1];
            MPI.Request[] requests = new MPI.Request[2] { null, null };
            double[] out_buffer_x = new double[buffer_size];
            //istart = 0;
            //---------------------------------------------------------------------
            //     in our terminology stage is the number of the cell in the x-direction
            //     i.e. stage = 1 means the start of the line stage=ncells means end
            //---------------------------------------------------------------------
            for(stage = 0; stage < ncells; stage++) {
                c = slice[stage, 0];
                //isize = cell_size[c,0] - 1;
                //jsize = cell_size[c,1] - 1;
                //ksize = cell_size[c,2] - 1;
                //---------------------------------------------------------------------
                //     set last-cell flag
                //---------------------------------------------------------------------
                if(stage == ncells-1) {
                    last = 1;
                }
                else {
                    last = 0;
                }
                if(stage == 0) {
                    //---------------------------------------------------------------------
                    //     This is the first cell, so solve without receiving data
                    //---------------------------------------------------------------------
                    first = 1;
                    //          //c            call lhsx[c];
                    x_solve_cell(first, last, c);
                }
                else {
                    //---------------------------------------------------------------------
                    //     Not the first cell of this line, so receive info from
                    //     processor working on preceeding cell
                    //---------------------------------------------------------------------
                    first = 0;
                    int jp = cell_coord[c,1];
                    int kp = cell_coord[c,2];
                    requests[0] = comm_solve.ImmediateReceive<double>(predecessor[0], WEST + jp + kp * ncells, out_buffer_x); //x_receive_solve_info(out_buffer_x, requests, c); //recv_id, c);  //x_receive_solve_info(recv_id,c);
                    //---------------------------------------------------------------------
                    //     overlap computations and communications
                    //---------------------------------------------------------------------
                    //            call lhsx[c]
                    //---------------------------------------------------------------------
                    //     wait for completion
                    //---------------------------------------------------------------------
                    requests[1].Wait(); //send_id[0].Wait();
                    requests[0].Wait(); //recv_id[0].Wait();
                    //Fortran: call mpi_wait[send_id,r_status,error]
                    //Fortran: call mpi_wait[recv_id,r_status,error]
                    //---------------------------------------------------------------------
                    //     install C'[istart] and rhs'[istart] to be used in this cell
                    //---------------------------------------------------------------------
                    x_unpack_solve_info(out_buffer_x, c);
                    x_solve_cell(first, last, c);
                }
                if(last == 0) {
                    double[] in_buffer_x = new double[buffer_size];//buffer_size=(MAX_CELL_DIM*MAX_CELL_DIM*(5*5+5))

                    x_pack_solve_info(in_buffer_x, c); //send_id, c); //x_send_solve_info(send_id,c);
                    
                    int jp = cell_coord[c,1];
                    int kp = cell_coord[c,2];
                    requests[1] = comm_solve.ImmediateSend<double>(in_buffer_x, successor[0], WEST+jp+kp*ncells);
                }
            }
            out_buffer_x = null;
            buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * 5;
            out_buffer_x = new double[buffer_size];
            //---------------------------------------------------------------------
            //     now perform backsubstitution in reverse direction
            //---------------------------------------------------------------------
            for(stage = ncells-1; stage >= 0; stage--) {  //for(stage = ncells, 1, -1;
                c = slice[stage, 0];
                first = 0;
                last = 0;
                if(stage == 0)
                    first = 1;
                if(stage == (ncells-1)) {
                    last = 1; //---------------------------------------------------------------------
                    //     last cell, so perform back substitute without waiting
                    //---------------------------------------------------------------------
                    x_backsubstitute(first, last, c); //call x_backsubstitute[first, last,c];
                }
                else {
                    int jp = cell_coord[c,1];
                    int kp = cell_coord[c,2];
                    
                    requests[0] = comm_solve.ImmediateReceive<double>(successor[0], EAST+jp+kp*ncells, out_buffer_x); //x_receive_backsub_info(out_buffer_x, requests, c); //recv_id, c);  //      call x_receive_backsub_info[recv_id,c];

                    requests[1].Wait(); //send_id[0].Wait();//      call mpi_wait[send_id,r_status,error];
                    requests[0].Wait(); //recv_id[0].Wait();//      call mpi_wait[recv_id,r_status,error];
                    
                    x_unpack_backsub_info(out_buffer_x, c);
                    x_backsubstitute(first, last, c);   //      call x_backsubstitute[first,last,c];
                }
                if(first == 0) {
                    int jp = cell_coord[c,1];
                    int kp = cell_coord[c,2];
                    double[] in_buffer_x = new double[buffer_size];//buffer_size=(MAX_CELL_DIM * MAX_CELL_DIM * 5)
                    x_pack_backsub_info(in_buffer_x, c); //send_id, c);  //call x_send_backsub_info[send_id,c];
                    requests[1] = comm_solve.ImmediateSend<double>(in_buffer_x, predecessor[0], EAST+jp+kp*ncells);                    
                }
            }
        }

        public void x_unpack_solve_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack C'[-1] and rhs'[-1] for
            //     all j and k
            //---------------------------------------------------------------------
            int j, k, m, n, ptr, istart;
            istart = 2;
            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(j = 0; j <= JMAX - 1; j++) {
                    for(m = 1; m <= 5; m++) {
                        for(n = 0; n < 5; n++) {
                            lhsc[c, k+2, j+2, istart-1, n, m-1] = out_buffer_x[ptr + n]; //lhsc[m,n,istart-1,j,k,c] = out_buffer[ptr+n];
                        }
                        ptr = ptr + 5;
                    }
                    for(n = 0; n < 5; n++) {
                        rhs[c, k+2, j+2, istart-1, n] = out_buffer_x[ptr + n]; //rhs[n,istart-1,j,k,c] = out_buffer[ptr+n];
                    }
                    ptr = ptr + 5;
                }
            }
        }

        public void x_pack_solve_info(double[] in_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     pack up and send C'[iend] and rhs'[iend] for
            //     all j and k
            //---------------------------------------------------------------------
            int j, k, m, n, isize, ptr;
            isize = cell_size[c, 0] + 1;
            //---------------------------------------------------------------------
            //     pack up buffer
            //---------------------------------------------------------------------
            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(j = 0; j <= JMAX - 1; j++) {
                    for(m = 1; m <= 5; m++) {
                        for(n = 0; n < 5; n++) {
                            in_buffer_x[ptr + n] = lhsc[c, k+2, j+2, isize, n, m-1];  //in_buffer[ptr+n] = lhsc[m,n,isize,j,k,c];
                        }
                        ptr = ptr + 5;
                    }
                    for(n = 0; n < 5; n++) {
                        in_buffer_x[ptr + n] = rhs[c, k+2, j+2, isize, n];  //in_buffer[ptr+n] = rhs[n,isize,j,k,c];
                    }
                    ptr = ptr + 5;
                }
            }
            //---------------------------------------------------------------------
            //     send buffer 
            //---------------------------------------------------------------------
            //call mpi_isend[in_buffer, buffer_size,dp_type, successor[1],WEST+jp+kp*NCELLS, comm_solve, send_id,error];
            //send_id[0] =
        }

        public void x_pack_backsub_info(double[] in_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     pack up and send U[istart] for all j and k
            //---------------------------------------------------------------------
            int j, k, n, ptr, istart;
            int buffer_size;
            //---------------------------------------------------------------------
            //     Send element 0 to previous processor
            //---------------------------------------------------------------------
            istart = 2;

            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(j = 0; j <= JMAX - 1; j++) {
                    for(n = 0; n < 5; n++) {
                        in_buffer_x[ptr + n] = rhs[c, k+2, j+2, istart, n]; //in_buffer[ptr+n] = rhs[n,istart,j,k,c];
                    }
                    ptr = ptr + 5;
                }
            }
            //call mpi_isend[in_buffer, buffer_size,dp_type, predecessor[1], EAST+jp+kp*NCELLS, comm_solve, send_id,error];
            //send_id[0] =
        }

        public void x_unpack_backsub_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack U[isize] for all j and k
            //---------------------------------------------------------------------
            int j, k, n, ptr;
            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(j = 0; j <= JMAX - 1; j++) {
                    for(n = 0; n < 5; n++) {
                        backsub_info[c, k+2, j+2, n] = out_buffer_x[ptr + n];  //backsub_info[n,j,k,c] = out_buffer[ptr+n];
                    }
                    ptr = ptr + 5;
                }
            }
        }

        public void x_backsubstitute(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     back solve: if last cell, { generate U[isize]=rhs[isize]
            //     } else { assume U[isize] is loaded in un pack backsub_info
            //     so just use it
            //     after call u[istart] will be sent to next cell
            //---------------------------------------------------------------------
            int i, j, k;
            int m, n, isize, jsize, ksize, istart;

            istart = 2;
            isize = cell_size[c, 0] + 1;
            jsize = cell_size[c, 1] - end[c, 1] + 1;
            ksize = cell_size[c, 2] - end[c, 2] + 1;
            if(last == 0) {
                for(k = start[c, 2]; k <= ksize; k++) {
                    for(j = start[c, 1]; j <= jsize; j++) {
                        //---------------------------------------------------------------------
                        //     U[isize] uses info from previous cell if not last cell
                        //---------------------------------------------------------------------
                        for(m = 0; m < 5; m++) {
                            for(n = 0; n < 5; n++) {//rhs[m,isize,j,k,c] = rhs[m,isize,j,k,c] - lhsc[m,n,isize,j,k,c]*backsub_info[n,j,k,c]
                                rhs[c, k, j, isize, m] = rhs[c, k, j, isize, m] - lhsc[c, k, j, isize, n, m] * backsub_info[c, k, j, n];
                            }
                        }
                    }
                }
            }
            for(k = start[c, 2]; k <= ksize; k++) {
                for(j = start[c, 1]; j <= jsize; j++) {
                    for(i = isize-1; i >= istart; i--) {  //for(i=isize-1,istart,-1;
                        for(m = 0; m < 5; m++) {
                            for(n = 0; n < 5; n++) { //rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - lhsc[m,n,i,j,k,c]*rhs[n,i+1,j,k,c];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - lhsc[c, k, j, i, n, m] * rhs[c, k, j, i+1, n];
                            }
                        }
                    }
                }
            }
        }//isize ksize jsize start[ istart

        public void x_solve_cell(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     performs guaussian elimination on this cell.
            //     
            //     assumes that unpacking routines for non-first cells 
            //     preload C' and rhs' from previous cell.
            //     
            //     assumed send happens outside this routine, but that
            //     c'[IMAX] and rhs'[IMAX] will be sent to next cell
            //---------------------------------------------------------------------
            int i, j, k, isize, ksize, jsize, istart;
            double tmp1, tmp2, tmp3;
            istart = 2;
            isize = cell_size[c, 0] + 1;
            jsize = cell_size[c, 1] - end[c, 1] + 1;
            ksize = cell_size[c, 2] - end[c, 2] + 1;

            lhsabinit(lhsa, lhsb, isize);
            for(k = start[c, 2]; k <= ksize; k++) {
                for(j = start[c, 1]; j <= jsize; j++) {
                    //---------------------------------------------------------------------
                    //     This function computes the left hand side in the xi-direction
                    //---------------------------------------------------------------------
                    //---------------------------------------------------------------------
                    //     determine a [labeled f] and n jacobians for cell c
                    //---------------------------------------------------------------------
                    for(i = start[c, 0]-1; i <= cell_size[c, 0] - end[c, 0]+2; i++) {
                        tmp1 = rho_i[c, k, j, i];  //rho_i[i,j,k,c];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        fjac[i+1, 0, 0] = 0.0d; // fjac[1, 1, i] = 0.0d;
                        fjac[i+1, 1, 0] = 1.0d; // fjac[1, 2, i] = 1.0d;
                        fjac[i+1, 2, 0] = 0.0d; // fjac[1, 3, i] = 0.0d;
                        fjac[i+1, 3, 0] = 0.0d; // fjac[1, 4, i] = 0.0d;
                        fjac[i+1, 4, 0] = 0.0d; // fjac[1, 5, i] = 0.0d;

                        fjac[i+1, 0, 1] = -(u[c, k, j, i, 1] * tmp2 * u[c, k, j, i, 1]) + c2 * qs[c, k, j, i];
                        fjac[i+1, 1, 1] = (2.0d - c2) * (u[c, k, j, i, 1] * tmp1);
                        fjac[i+1, 2, 1] = -c2 * (u[c, k, j, i, 2] * tmp1);
                        fjac[i+1, 3, 1] = -c2 * (u[c, k, j, i, 3] * tmp1);
                        fjac[i+1, 4, 1] = c2;

                        fjac[i+1, 0, 2] = -(u[c, k, j, i, 1] * u[c, k, j, i, 2]) * tmp2;
                        fjac[i+1, 1, 2] = u[c, k, j, i, 2] * tmp1;
                        fjac[i+1, 2, 2] = u[c, k, j, i, 1] * tmp1;
                        fjac[i+1, 3, 2] = 0.0d;
                        fjac[i+1, 4, 2] = 0.0d;

                        fjac[i+1, 0, 3] = -(u[c, k, j, i, 1] * u[c, k, j, i, 3]) * tmp2;
                        fjac[i+1, 1, 3] = u[c, k, j, i, 3] * tmp1;
                        fjac[i+1, 2, 3] = 0.0d;
                        fjac[i+1, 3, 3] = u[c, k, j, i, 1] * tmp1;
                        fjac[i+1, 4, 3] = 0.0d;

                        fjac[i+1, 0, 4] = (c2 * 2.0d * qs[c, k, j, i] - c1 * (u[c, k, j, i, 4] * tmp1)) * (u[c, k, j, i, 1] * tmp1);
                        fjac[i+1, 1, 4] = c1 * u[c, k, j, i, 4] * tmp1 - c2 * (u[c, k, j, i, 1] * u[c, k, j, i, 1] * tmp2 + qs[c, k, j, i]);
                        fjac[i+1, 2, 4] = -c2 * (u[c, k, j, i, 2] * u[c, k, j, i, 1]) * tmp2;
                        fjac[i+1, 3, 4] = -c2 * (u[c, k, j, i, 3] * u[c, k, j, i, 1]) * tmp2;
                        fjac[i+1, 4, 4] = c1 * (u[c, k, j, i, 1] * tmp1);

                        njac[i+1, 0, 0] = 0.0d;
                        njac[i+1, 1, 0] = 0.0d;
                        njac[i+1, 2, 0] = 0.0d;
                        njac[i+1, 3, 0] = 0.0d;
                        njac[i+1, 4, 0] = 0.0d;

                        njac[i+1, 0, 1] = -con43 * c3c4 * tmp2 * u[c, k, j, i, 1];  //njac[2 + i, 1, 2] = -con43 * c3c4 * tmp2 * u[2, i, j, k, c];
                        njac[i+1, 1, 1] = con43 * c3c4 * tmp1;
                        njac[i+1, 2, 1] = 0.0d;
                        njac[i+1, 3, 1] = 0.0d;
                        njac[i+1, 4, 1] = 0.0d;

                        njac[i+1, 0, 2] = -c3c4 * tmp2 * u[c, k, j, i, 2];  //c3c4 * tmp2 * u[3,i,j,k,c];
                        njac[i+1, 1, 2] = 0.0d;
                        njac[i+1, 2, 2] = c3c4 * tmp1;
                        njac[i+1, 3, 2] = 0.0d;
                        njac[i+1, 4, 2] = 0.0d;

                        njac[i+1, 0, 3] = -c3c4 * tmp2 * u[c, k, j, i, 3];  //c3c4 * tmp2 * u[4,i,j,k,c];
                        njac[i+1, 1, 3] = 0.0d;
                        njac[i+1, 2, 3] = 0.0d;
                        njac[i+1, 3, 3] = c3c4 * tmp1;
                        njac[i+1, 4, 3] = 0.0d;

                        njac[i+1, 0, 4] = -(con43 * c3c4
                            - c1345) * tmp3 * (pow2(u[c, k, j, i, 1]))           /*- c1345 ) * tmp3 * (u[2,i,j,k,c]**2) */
                            - (c3c4 - c1345) * tmp3 * (pow2(u[c, k, j, i, 2]))  /*- ( c3c4 - c1345 ) * tmp3 * (u[3,i,j,k,c]**2) */
                            - (c3c4 - c1345) * tmp3 * (pow2(u[c, k, j, i, 3]))  /*- ( c3c4 - c1345 ) * tmp3 * (u[4,i,j,k,c]**2) */
                            - c1345 * tmp2 * u[c, k, j, i, 4];  //- c1345 * tmp2 * u[5,i,j,k,c];

                        njac[i+1, 1, 4] = (con43 * c3c4 - c1345) * tmp2 * u[c, k, j, i, 1];  // u[2,i,j,k,c];
                        njac[i+1, 2, 4] = (c3c4 - c1345) * tmp2 * u[c, k, j, i, 2];          // u[3,i,j,k,c];
                        njac[i+1, 3, 4] = (c3c4 - c1345) * tmp2 * u[c, k, j, i, 3];          // u[4,i,j,k,c];
                        njac[i+1, 4, 4] = (c1345) * tmp1;
                    }
                    //---------------------------------------------------------------------
                    //     now jacobians set, so form left hand side in x direction
                    //---------------------------------------------------------------------
                    for(i = start[c, 0]; i <= isize - end[c, 0]; i++) {

                        tmp1 = dt * tx1;
                        tmp2 = dt * tx2;

                        lhsa[i, 0, 0] = -tmp2 * fjac[i, 0, 0] - tmp1 * njac[i, 0, 0] - tmp1 * dx1;//2+i-1
                        lhsa[i, 1, 0] = -tmp2 * fjac[i, 1, 0] - tmp1 * njac[i, 1, 0];
                        lhsa[i, 2, 0] = -tmp2 * fjac[i, 2, 0] - tmp1 * njac[i, 2, 0];
                        lhsa[i, 3, 0] = -tmp2 * fjac[i, 3, 0] - tmp1 * njac[i, 3, 0];
                        lhsa[i, 4, 0] = -tmp2 * fjac[i, 4, 0] - tmp1 * njac[i, 4, 0];

                        lhsa[i, 0, 1] = -tmp2 * fjac[i, 0, 1] - tmp1 * njac[i, 0, 1];
                        lhsa[i, 1, 1] = -tmp2 * fjac[i, 1, 1] - tmp1 * njac[i, 1, 1] - tmp1 * dx2;
                        lhsa[i, 2, 1] = -tmp2 * fjac[i, 2, 1] - tmp1 * njac[i, 2, 1];
                        lhsa[i, 3, 1] = -tmp2 * fjac[i, 3, 1] - tmp1 * njac[i, 3, 1];
                        lhsa[i, 4, 1] = -tmp2 * fjac[i, 4, 1] - tmp1 * njac[i, 4, 1];

                        lhsa[i, 0, 2] = -tmp2 * fjac[i, 0, 2] - tmp1 * njac[i, 0, 2];
                        lhsa[i, 1, 2] = -tmp2 * fjac[i, 1, 2] - tmp1 * njac[i, 1, 2];
                        lhsa[i, 2, 2] = -tmp2 * fjac[i, 2, 2] - tmp1 * njac[i, 2, 2] - tmp1 * dx3;
                        lhsa[i, 3, 2] = -tmp2 * fjac[i, 3, 2] - tmp1 * njac[i, 3, 2];
                        lhsa[i, 4, 2] = -tmp2 * fjac[i, 4, 2] - tmp1 * njac[i, 4, 2];

                        lhsa[i, 0, 3] = -tmp2 * fjac[i, 0, 3] - tmp1 * njac[i, 0, 3];
                        lhsa[i, 1, 3] = -tmp2 * fjac[i, 1, 3] - tmp1 * njac[i, 1, 3];
                        lhsa[i, 2, 3] = -tmp2 * fjac[i, 2, 3] - tmp1 * njac[i, 2, 3];
                        lhsa[i, 3, 3] = -tmp2 * fjac[i, 3, 3] - tmp1 * njac[i, 3, 3] - tmp1 * dx4;
                        lhsa[i, 4, 3] = -tmp2 * fjac[i, 4, 3] - tmp1 * njac[i, 4, 3];

                        lhsa[i, 0, 4] = -tmp2 * fjac[i, 0, 4] - tmp1 * njac[i, 0, 4];
                        lhsa[i, 1, 4] = -tmp2 * fjac[i, 1, 4] - tmp1 * njac[i, 1, 4];
                        lhsa[i, 2, 4] = -tmp2 * fjac[i, 2, 4] - tmp1 * njac[i, 2, 4];
                        lhsa[i, 3, 4] = -tmp2 * fjac[i, 3, 4] - tmp1 * njac[i, 3, 4];
                        lhsa[i, 4, 4] = -tmp2 * fjac[i, 4, 4] - tmp1 * njac[i, 4, 4] - tmp1 * dx5;

                        lhsb[i, 0, 0] = 1.0d + tmp1 * 2.0d * njac[i+1, 0, 0] + tmp1 * 2.0d * dx1;
                        lhsb[i, 1, 0] = tmp1 * 2.0d * njac[i+1, 1, 0];
                        lhsb[i, 2, 0] = tmp1 * 2.0d * njac[i+1, 2, 0];
                        lhsb[i, 3, 0] = tmp1 * 2.0d * njac[i+1, 3, 0];
                        lhsb[i, 4, 0] = tmp1 * 2.0d * njac[i+1, 4, 0];

                        lhsb[i, 0, 1] = tmp1 * 2.0d * njac[i+1, 0, 1];
                        lhsb[i, 1, 1] = 1.0d + tmp1 * 2.0d * njac[i+1, 1, 1] + tmp1 * 2.0d * dx2;
                        lhsb[i, 2, 1] = tmp1 * 2.0d * njac[i+1, 2, 1];
                        lhsb[i, 3, 1] = tmp1 * 2.0d * njac[i+1, 3, 1];
                        lhsb[i, 4, 1] = tmp1 * 2.0d * njac[i+1, 4, 1];

                        lhsb[i, 0, 2] = tmp1 * 2.0d * njac[i+1, 0, 2];
                        lhsb[i, 1, 2] = tmp1 * 2.0d * njac[i+1, 1, 2];
                        lhsb[i, 2, 2] = 1.0d + tmp1 * 2.0d * njac[i+1, 2, 2] + tmp1 * 2.0d * dx3;
                        lhsb[i, 3, 2] = tmp1 * 2.0d * njac[i+1, 3, 2];
                        lhsb[i, 4, 2] = tmp1 * 2.0d * njac[i+1, 4, 2];

                        lhsb[i, 0, 3] = tmp1 * 2.0d * njac[i+1, 0, 3];
                        lhsb[i, 1, 3] = tmp1 * 2.0d * njac[i+1, 1, 3];
                        lhsb[i, 2, 3] = tmp1 * 2.0d * njac[i+1, 2, 3];
                        lhsb[i, 3, 3] = 1.0d + tmp1 * 2.0d * njac[i+1, 3, 3] + tmp1 * 2.0d * dx4;
                        lhsb[i, 4, 3] = tmp1 * 2.0d * njac[i+1, 4, 3];

                        lhsb[i, 0, 4] = tmp1 * 2.0d * njac[i+1, 0, 4];
                        lhsb[i, 1, 4] = tmp1 * 2.0d * njac[i+1, 1, 4];
                        lhsb[i, 2, 4] = tmp1 * 2.0d * njac[i+1, 2, 4];
                        lhsb[i, 3, 4] = tmp1 * 2.0d * njac[i+1, 3, 4];
                        lhsb[i, 4, 4] = 1.0d + tmp1 * 2.0d * njac[i+1, 4, 4] + tmp1 * 2.0d * dx5;

                        lhsc[c, k, j, i, 0, 0] = tmp2 * fjac[i+2, 0, 0] - tmp1 * njac[i+2, 0, 0] - tmp1 * dx1; //2+i+1
                        lhsc[c, k, j, i, 1, 0] = tmp2 * fjac[i+2, 1, 0] - tmp1 * njac[i+2, 1, 0];
                        lhsc[c, k, j, i, 2, 0] = tmp2 * fjac[i+2, 2, 0] - tmp1 * njac[i+2, 2, 0];
                        lhsc[c, k, j, i, 3, 0] = tmp2 * fjac[i+2, 3, 0] - tmp1 * njac[i+2, 3, 0];
                        lhsc[c, k, j, i, 4, 0] = tmp2 * fjac[i+2, 4, 0] - tmp1 * njac[i+2, 4, 0];

                        lhsc[c, k, j, i, 0, 1] = tmp2 * fjac[i+2, 0, 1] - tmp1 * njac[i+2, 0, 1];
                        lhsc[c, k, j, i, 1, 1] = tmp2 * fjac[i+2, 1, 1] - tmp1 * njac[i+2, 1, 1] - tmp1 * dx2;
                        lhsc[c, k, j, i, 2, 1] = tmp2 * fjac[i+2, 2, 1] - tmp1 * njac[i+2, 2, 1];
                        lhsc[c, k, j, i, 3, 1] = tmp2 * fjac[i+2, 3, 1] - tmp1 * njac[i+2, 3, 1];
                        lhsc[c, k, j, i, 4, 1] = tmp2 * fjac[i+2, 4, 1] - tmp1 * njac[i+2, 4, 1];

                        lhsc[c, k, j, i, 0, 2] = tmp2 * fjac[i+2, 0, 2] - tmp1 * njac[i+2, 0, 2];
                        lhsc[c, k, j, i, 1, 2] = tmp2 * fjac[i+2, 1, 2] - tmp1 * njac[i+2, 1, 2];
                        lhsc[c, k, j, i, 2, 2] = tmp2 * fjac[i+2, 2, 2] - tmp1 * njac[i+2, 2, 2] - tmp1 * dx3;
                        lhsc[c, k, j, i, 3, 2] = tmp2 * fjac[i+2, 3, 2] - tmp1 * njac[i+2, 3, 2];
                        lhsc[c, k, j, i, 4, 2] = tmp2 * fjac[i+2, 4, 2] - tmp1 * njac[i+2, 4, 2];

                        lhsc[c, k, j, i, 0, 3] = tmp2 * fjac[i+2, 0, 3] - tmp1 * njac[i+2, 0, 3];
                        lhsc[c, k, j, i, 1, 3] = tmp2 * fjac[i+2, 1, 3] - tmp1 * njac[i+2, 1, 3];
                        lhsc[c, k, j, i, 2, 3] = tmp2 * fjac[i+2, 2, 3] - tmp1 * njac[i+2, 2, 3];
                        lhsc[c, k, j, i, 3, 3] = tmp2 * fjac[i+2, 3, 3] - tmp1 * njac[i+2, 3, 3] - tmp1 * dx4;
                        lhsc[c, k, j, i, 4, 3] = tmp2 * fjac[i+2, 4, 3] - tmp1 * njac[i+2, 4, 3];

                        lhsc[c, k, j, i, 0, 4] = tmp2 * fjac[i+2, 0, 4] - tmp1 * njac[i+2, 0, 4];
                        lhsc[c, k, j, i, 1, 4] = tmp2 * fjac[i+2, 1, 4] - tmp1 * njac[i+2, 1, 4];
                        lhsc[c, k, j, i, 2, 4] = tmp2 * fjac[i+2, 2, 4] - tmp1 * njac[i+2, 2, 4];
                        lhsc[c, k, j, i, 3, 4] = tmp2 * fjac[i+2, 3, 4] - tmp1 * njac[i+2, 3, 4];
                        lhsc[c, k, j, i, 4, 4] = tmp2 * fjac[i+2, 4, 4] - tmp1 * njac[i+2, 4, 4] - tmp1 * dx5;
                    }
                    //---------------------------------------------------------------------
                    //     outer most for(loops - sweeping in i direction
                    //---------------------------------------------------------------------
                    if(first == 1) {
                        //---------------------------------------------------------------------
                        //     multiply c[istart,j,k] by b_inverse and copy back to c
                        //     multiply rhs[istart] by b_inverse[istart] and copy to rhs
                        //---------------------------------------------------------------------
                        /* call binvcrhs(lhsb[1,1,istart], lhsc[1,1,istart,j,k,c], rhs[1,istart,j,k,c]); */
                        //binvcrhs(lhsb, lhsc, rhs, istart-1, c+1, k-1, j-1, istart-1, c+1, k-1, j-1, istart-1);
                        binvcrhs(lhsb, lhsc, rhs, istart, c, k, j, istart, c, k, j, istart);
                        //binvcrhs(lhsb[istart+1, 1, 1], lhsc[c, k+1, j+1, istart+1, 1, 1], rhs[c, k+1, j+1, istart+1, 1]);  
                    }
                    //---------------------------------------------------------------------
                    //     begin inner most for(loop
                    //     for(all the elements of the cell unless last 
                    //---------------------------------------------------------------------
                    for(i = istart + first; i <= isize - last; i++) {
                        //---------------------------------------------------------------------
                        //     rhs[i] = rhs[i] - A*rhs[i-1]
                        //---------------------------------------------------------------------
                        /*Fortran:  matvec_sub(lhsa[1,1,i], rhs[1,i-1,j,k,c],rhs[1,i,j,k,c]);*/
                        /*Padrao c# matvec_sub(lhsa[i+1,1,1], rhs[c,k+1,j+1,i-1+1,1],rhs[c,k+1,j+1,i+1,1]);*/
                        //matvec_sub(ref lhsa, ref rhs, ref rhs, (i-1), c+1, (k-1), (j-1), i-2, c+1, (k-1), (j-1), (i-1));
                        matvec_sub(lhsa, rhs, rhs, i, c, k, j, i-1, c, k, j, i);
                        //---------------------------------------------------------------------
                        //     B[i] = B[i] - C[i-1]*A[i]
                        //---------------------------------------------------------------------
                        /*Fortran:   matmul_sub(lhsa[1,1,i], lhsc[1,1,i-1,j,k,c], lhsb[1,1,i]);*/
                        /*Padrao c#: matmul_sub(lhsa[i+1,1,1], lhsc[c,k+1,j+1,i-1+1,1,1], lhsb[i+1,1,1]);*/
                        //matmul_sub(lhsa, lhsc, lhsb, (i-1), c+1, (k-1), (j-1), i-2, (i-1));
                        matmul_sub(lhsa, lhsc, lhsb, i, c, k, j, i-1, i);
                        //---------------------------------------------------------------------
                        //     multiply c[i,j,k] by b_inverse and copy back to c
                        //     multiply rhs[1,j,k] by b_inverse[1,j,k] and copy to rhs
                        //---------------------------------------------------------------------
                        /*Fortran binvcrhs( lhsb[1,1,i], lhsc[1,1,i,j,k,c], rhs[1,i,j,k,c] );*/
                        //binvcrhs(lhsb, lhsc, rhs, i-1, c+1, k-1, j-1, i-1, c+1, k-1, j-1, i-1);
                        binvcrhs(lhsb, lhsc, rhs, i, c, k, j, i, c, k, j, i);
                    }
                    //---------------------------------------------------------------------
                    //     Now finish up special cases for last cell
                    //---------------------------------------------------------------------
                    if(last == 1) {
                        //---------------------------------------------------------------------
                        //     rhs[isize] = rhs[isize] - A*rhs[isize-1]
                        //---------------------------------------------------------------------
                        /*Fortran: matvec_sub[lhsa[1,1,isize], rhs[1,isize-1,j,k,c],rhs[1,isize,j,k,c]];*/
                        /*C#:      matvec_sub[lhsa[isize+1,1,1], rhs[c,k+1,j+1,isize-1+1,1],rhs[c,k+1,j+1,isize+1,1]];*/
                        //matvec_sub(ref lhsa, ref rhs, ref rhs, (isize-1), c+1, (k-1), (j-1), (isize-2), c+1, (k-1), (j-1), (isize-1));
                        matvec_sub(lhsa, rhs, rhs, isize, c, k, j, isize-1, c, k, j, isize);
                        //---------------------------------------------------------------------
                        //     B[isize] = B[isize] - C[isize-1]*A[isize]
                        //---------------------------------------------------------------------
                        /*Fortran: matmul_sub[lhsa[1,1,isize], lhsc[1,1,isize-1,j,k,c], lhsb[1,1,isize]];*/
                        /*C#:      matmul_sub(lhsa[isize+1,1,1], lhsc[c,k+1,j+1,isize-1+1,1,1], lhsb[isize+1,1,1]);*/
                        //matmul_sub(ref lhsa, ref lhsc, ref lhsb, (isize-1), c+1, (k-1), (j-1), isize-2, (isize-1));
                        matmul_sub(lhsa, lhsc, lhsb, isize, c, k, j, isize-1, isize);
                        //---------------------------------------------------------------------
                        //     multiply rhs[] by b_inverse[] and copy to rhs
                        //---------------------------------------------------------------------
                        /*Fortran: binvrhs[ lhsb[1,1,isize], rhs[1,isize,j,k,c] ];*/
                        /*C#:      binvrhs( lhsb[isize+1,1,1], rhs[c,k+1,j+1,isize+1,1] ); */
                        //binvrhs(ref lhsb, ref rhs, (isize-1), c+1, (k-1), (j-1), (isize-1));
                        binvrhs(lhsb, rhs, isize, c, k, j, isize);
                    }
                }
            }
        }//start[ isize jsize ksize istart

        /* x_solve **********************************************************************************************/


        public void y_solve() {
            //---------------------------------------------------------------------
            //     Performs line solves in Y direction by first factoring
            //     the block-tridiagonal matrix into an upper triangular matrix, 
            //     and { performing back substitution to solve for the unknow
            //     vectors of each line.  
            //     
            //     Make sure we treat elements zero to cell_size in the direction
            //     of the sweep.
            //---------------------------------------------------------------------
            int c, stage, first, last, buffer_size;//int r_status[MPI_STATUS_SIZE];
            buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * (5 * 5 + 5);
            //MPI.Request[] recv_id = new MPI.Request[1];
            //MPI.Request[] send_id = new MPI.Request[1];
            MPI.Request[] requests=new MPI.Request[2] { null, null };
            double[] out_buffer_y = new double[buffer_size];
            //---------------------------------------------------------------------
            //     in our terminology stage is the number of the cell in the y-direction
            //     i.e. stage = 1 means the start of the line stage=ncells means end
            //---------------------------------------------------------------------
            for(stage = 0; stage < ncells; stage++) {
                c = slice[stage, 1];
                //isize = cell_size[c, 0] - 1;
                //jsize = cell_size[c, 1] - 1;
                //ksize = cell_size[c, 2] - 1;
                //---------------------------------------------------------------------
                //     set last-cell flag
                //---------------------------------------------------------------------
                if(stage == ncells-1) {
                    last = 1;
                }
                else {
                    last = 0;
                }
                if(stage == 0) {
                    //---------------------------------------------------------------------
                    //     This is the first cell, so solve without receiving data
                    //---------------------------------------------------------------------
                    first = 1;
                    //c            call lhsy[c]
                    y_solve_cell(first, last, c); //call y_solve_cell[first,last,c];
                }
                else {
                    //c---------------------------------------------------------------------
                    //c     Not the first cell of this line, so receive info from
                    //c     processor working on preceeding cell
                    //c---------------------------------------------------------------------
                    first = 0;
                    int ip = cell_coord[c,0];
                    int kp = cell_coord[c,2];
                    requests[0] = comm_solve.ImmediateReceive<double>(predecessor[1], SOUTH+ip+kp*ncells, out_buffer_y);//y_receive_solve_info(out_buffer_y, requests, c); //recv_id, c); //call y_receive_solve_info[recv_id,c];
                    //      c---------------------------------------------------------------------
                    //      c     overlap computations and communications
                    //      c---------------------------------------------------------------------
                    //      c            call lhsy[c]
                    //      c---------------------------------------------------------------------
                    //      c     wait for completion
                    //      c---------------------------------------------------------------------
                    requests[1].Wait(); //send_id[0].Wait();
                    requests[0].Wait(); //recv_id[0].Wait();
                    //call mpi_wait[send_id,r_status,error];
                    //call mpi_wait[recv_id,r_status,error];
                    //      c---------------------------------------------------------------------
                    //      c     install C'[jstart+1] and rhs'[jstart+1] to be used in this cell
                    //      c---------------------------------------------------------------------
                    y_unpack_solve_info(out_buffer_y, c);
                    y_solve_cell(first, last, c); //call y_solve_cell[first,last,c];
                }
                if(last == 0) {
                    int ip = cell_coord[c,0];
                    int kp = cell_coord[c,2];
                    double[] in_buffer_y = new double[buffer_size];//buffer_size=MAX_CELL_DIM * MAX_CELL_DIM * (5 * 5 + 5)
                    y_pack_solve_info(in_buffer_y, c); //send_id, c);  //call y_send_solve_info[send_id,c];
                    requests[1] = comm_solve.ImmediateSend<double>(in_buffer_y, successor[1], SOUTH+ip+kp*ncells);
                }
            }
            out_buffer_y = null;
            buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * 5;
            out_buffer_y = new double[buffer_size];
            //---------------------------------------------------------------------
            //     now perform backsubstitution in reverse direction
            //---------------------------------------------------------------------
            for(stage = ncells-1; stage >= 0; stage--) {  //for(stage = ncells, 1, -1
                c = slice[stage, 1];
                first = 0;
                last = 0;
                if(stage == 0)
                    first = 1;
                if(stage == ncells-1) {
                    last = 1;
                    //---------------------------------------------------------------------
                    //     last cell, so perform back substitute without waiting
                    //---------------------------------------------------------------------
                    y_backsubstitute(first, last, c);     //call y_backsubstitute[first, last,c];
                }
                else {
                    int ip = cell_coord[c, 0];
                    int kp = cell_coord[c, 2];

                    requests[0] = comm_solve.ImmediateReceive<double>(successor[1], NORTH+ip+kp*ncells, out_buffer_y);//y_receive_backsub_info(out_buffer_y, requests, c); //recv_id, c);   //call y_receive_backsub_info[recv_id,c];

                    requests[1].Wait();//send_id[0].Wait();
                    requests[0].Wait();//recv_id[0].Wait();
                    //call mpi_wait[send_id,r_status,error];
                    //call mpi_wait[recv_id,r_status,error];
                    y_unpack_backsub_info(out_buffer_y, c);
                    y_backsubstitute(first, last, c);      //call y_backsubstitute[first,last,c];
                }
                if(first == 0) {
                    int ip = cell_coord[c,0];
                    int kp = cell_coord[c,2];                    
                    double[] in_buffer_y = new double[buffer_size];//buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * 5;                    
                    y_pack_backsub_info(in_buffer_y, c); //send_id, c);  //call y_send_backsub_info[send_id,c];
                    requests[1] = comm_solve.ImmediateSend<double>(in_buffer_y, predecessor[1], NORTH+ip+kp*ncells);
                }
            }
        }

        public void y_unpack_solve_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack C'[-1] and rhs'[-1] for
            //     all i and k
            //---------------------------------------------------------------------
            int i, k, m, n, ptr, jstart;
            jstart = 2;
            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(m = 1; m <= 5; m++) {
                        for(n = 0; n < 5; n++) {
                            lhsc[c, k+2, jstart-1, i+2, n, m-1] = out_buffer_x[ptr + n]; //lhsc[m,n,i,jstart-1,k,c] = out_buffer_x[ptr+n];
                        }
                        ptr = ptr + 5;
                    }
                    for(n = 0; n < 5; n++) {
                        rhs[c, k+2, jstart-1, i+2, n] = out_buffer_x[ptr + n];  // rhs[n,i,jstart-1,k,c] = out_buffer_x[ptr+n];
                    }
                    ptr = ptr + 5;
                }
            }
        }

        public void y_pack_solve_info(double[] in_buffer_y, int c) {
            //---------------------------------------------------------------------
            //     pack up and send C'[jend] and rhs'[jend] for
            //     all i and k
            //---------------------------------------------------------------------
            int i, k, m, n, jsize, ptr;
            int buffer_size;
            jsize = cell_size[c, 1] + 1;
            //---------------------------------------------------------------------
            //     pack up buffer
            //---------------------------------------------------------------------
            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(m = 1; m <= 5; m++) {
                        for(n = 0; n < 5; n++) {
                            in_buffer_y[ptr + n] = lhsc[c, k+2, jsize, i+2, n, m-1];  //in_buffer[ptr+n] = lhsc[m,n,i,jsize,k,c]
                        }
                        ptr = ptr + 5;
                    }
                    for(n = 0; n < 5; n++) {
                        in_buffer_y[ptr + n] = rhs[c, k+2, jsize, i+2, n];       //in_buffer[ptr+n] = rhs[n,i,jsize,k,c] 
                    }
                    ptr = ptr + 5;
                }
            }
            //---------------------------------------------------------------------
            //     send buffer 
            //---------------------------------------------------------------------
            //call mpi_isend[in_buffer, buffer_size, dp_type, successor[2], SOUTH+ip+kp*NCELLS, comm_solve, send_id,error];
            //send_id[0]=
        }//jsize

        public void y_pack_backsub_info(double[] in_buffer_y, int c) {
            //---------------------------------------------------------------------
            //     pack up and send U[jstart] for all i and k
            //---------------------------------------------------------------------
            int i, k, n, ptr, jstart;
            int buffer_size;
            //---------------------------------------------------------------------
            //     Send element 0 to previous processor
            //---------------------------------------------------------------------
            jstart = 2;
            
            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(n = 0; n < 5; n++) {
                        in_buffer_y[ptr + n] = rhs[c, k+2, jstart, i+2, n];  //in_buffer[ptr+n] = rhs[n,i,jstart,k,c];
                    }
                    ptr = ptr + 5;
                }
            }
            //call mpi_isend[in_buffer, buffer_size, dp_type, predecessor[2], NORTH+ip+kp*NCELLS, comm_solve, send_id,error];
            //send_id[0] =
        }

        public void y_unpack_backsub_info(double[] out_buffer_x, int c) {
            //  c---------------------------------------------------------------------
            //  c     unpack U[jsize] for all i and k
            //  c---------------------------------------------------------------------
            int i, k, n, ptr;

            ptr = 0;
            for(k = 0; k <= KMAX - 1; k++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(n = 0; n < 5; n++) {
                        backsub_info[c, k+2, i+2, n] = out_buffer_x[ptr + n];  //backsub_info[n,i,k,c] = out_buffer[ptr+n];
                    }
                    ptr = ptr + 5;
                }
            }
        }

        public void y_backsubstitute(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     back solve: if last cell, { generate U[jsize]=rhs[jsize]
            //     } else { assume U[jsize] is loaded in un pack backsub_info
            //     so just use it
            //     after call u[jstart] will be sent to next cell
            //---------------------------------------------------------------------
            int i, k;
            int m, n, j, jsize, isize, ksize, jstart;
            jstart = 2;
            isize = cell_size[c, 0] - end[c, 0] + 1;
            jsize = cell_size[c, 1] + 1;
            ksize = cell_size[c, 2] - end[c, 2] + 1;
            if(last == 0) {
                for(k = start[c, 2]; k <= ksize; k++) {
                    for(i = start[c, 0]; i <= isize; i++) {
                        //---------------------------------------------------------------------
                        //     U[jsize] uses info from previous cell if not last cell
                        //---------------------------------------------------------------------
                        for(m = 0; m < 5; m++) {
                            for(n = 0; n < 5; n++) { //rhs[m,i,jsize,k,c]=rhs[m,i,jsize,k,c]-lhsc[m,n,i,jsize,k,c]*backsub_info[n,i,k,c];
                                rhs[c, k, jsize, i, m] = rhs[c, k, jsize, i, m] - lhsc[c, k, jsize, i, n, m] * backsub_info[c, k, i, n];
                            }
                        }
                    }
                }
            }
            for(k = start[c, 2]; k <= ksize; k++) {
                for(j = jsize-1; j >= jstart; j--) {//for(j=jsize-1,jstart,-1;
                    for(i = start[c, 0]; i <= isize; i++) {
                        for(m = 0; m < 5; m++) {
                            for(n = 0; n < 5; n++) {  //rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - lhsc[m,n,i,j,k,c]*rhs[n,i,j+1,k,c];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - lhsc[c, k, j, i, n, m] * rhs[c, k, j+1, i, n];
                            }
                        }
                    }
                }
            }
        }//start[ isize jsize ksize jstart

        public void y_solve_cell(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     performs guaussian elimination on this cell.
            //     
            //     assumes that unpacking routines for non-first cells 
            //     preload C' and rhs' from previous cell.
            //     
            //     assumed send happens outside this routine, but that
            //     c'[JMAX] and rhs'[JMAX] will be sent to next cell
            //---------------------------------------------------------------------
            int i, j, k, isize, ksize, jsize, jstart;
            double tmp1, tmp2, tmp3;
            double[,] utmp = new double[JMAX + 4, 7];   //double utmp[6,-2:JMAX+1];

            jstart = 2;
            isize = cell_size[c, 0] - end[c, 0] + 1;
            jsize = cell_size[c, 1] + 1;
            ksize = cell_size[c, 2] - end[c, 2] + 1;

            lhsabinit(lhsa, lhsb, jsize);

            for(k = start[c, 2]; k <= ksize; k++) {
                for(i = start[c, 0]; i <= isize; i++) {
                    //---------------------------------------------------------------------
                    //     This function computes the left hand side for the three y-factors   
                    //---------------------------------------------------------------------
                    //---------------------------------------------------------------------
                    //     Compute the indices for storing the tri-diagonal matrix;
                    //     determine a [labeled f] and n jacobians for cell c
                    //---------------------------------------------------------------------
                    for(j = start[c, 1]-1; j <= cell_size[c, 1] - end[c, 1] + 2; j++) {
                        utmp[j, 1] = 1.0d / u[c, k, j, i, 0];  //u[1,i,j,k,c]
                        utmp[j, 2] = u[c, k, j, i, 1];         //u[2,i,j,k,c]; 
                        utmp[j, 3] = u[c, k, j, i, 2];         //u[3,i,j,k,c]; 
                        utmp[j, 4] = u[c, k, j, i, 3];         //u[4,i,j,k,c]; 
                        utmp[j, 5] = u[c, k, j, i, 4];         //u[5,i,j,k,c]; 
                        utmp[j, 6] = qs[c, k, j, i];           //qs[i,j,k,c]; 
                    }

                    for(j = start[c, 1]-1; j <= cell_size[c, 1] - end[c, 1] + 2; j++) {
                        tmp1 = utmp[j, 1];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        fjac[j+1, 0, 0] = 0.0d;
                        fjac[j+1, 1, 0] = 0.0d;
                        fjac[j+1, 2, 0] = 1.0d;
                        fjac[j+1, 3, 0] = 0.0d;
                        fjac[j+1, 4, 0] = 0.0d;

                        fjac[j+1, 0, 1] = -(utmp[j, 2] * utmp[j, 3]) * tmp2;
                        fjac[j+1, 1, 1] = utmp[j, 3] * tmp1;
                        fjac[j+1, 2, 1] = utmp[j, 2] * tmp1;
                        fjac[j+1, 3, 1] = 0.0d;
                        fjac[j+1, 4, 1] = 0.0d;

                        fjac[j+1, 0, 2] = -(utmp[j, 3] * utmp[j, 3] * tmp2) + c2 * utmp[j, 6];
                        fjac[j+1, 1, 2] = -c2 * utmp[j, 2] * tmp1;
                        fjac[j+1, 2, 2] = (2.0d - c2) * utmp[j, 3] * tmp1;
                        fjac[j+1, 3, 2] = -c2 * utmp[j, 4] * tmp1;
                        fjac[j+1, 4, 2] = c2;

                        fjac[j+1, 0, 3] = -(utmp[j, 3] * utmp[j, 4]) * tmp2;
                        fjac[j+1, 1, 3] = 0.0d;
                        fjac[j+1, 2, 3] = utmp[j, 4] * tmp1;
                        fjac[j+1, 3, 3] = utmp[j, 3] * tmp1;
                        fjac[j+1, 4, 3] = 0.0d;

                        fjac[j+1, 0, 4] = (c2 * 2.0d * utmp[j, 6] - c1 * utmp[j, 5] * tmp1) * utmp[j, 3] * tmp1;
                        fjac[j+1, 1, 4] = -c2 * utmp[j, 2] * utmp[j, 3] * tmp2;
                        fjac[j+1, 2, 4] = c1 * utmp[j, 5] * tmp1 - c2 * (utmp[j, 6] + utmp[j, 3] * utmp[j, 3] * tmp2);
                        fjac[j+1, 3, 4] = -c2 * (utmp[j, 3] * utmp[j, 4]) * tmp2;
                        fjac[j+1, 4, 4] = c1 * utmp[j, 3] * tmp1;

                        njac[j+1, 0, 0] = 0.0d;
                        njac[j+1, 1, 0] = 0.0d;
                        njac[j+1, 2, 0] = 0.0d;
                        njac[j+1, 3, 0] = 0.0d;
                        njac[j+1, 4, 0] = 0.0d;

                        njac[j+1, 0, 1] = -c3c4 * tmp2 * utmp[j, 2];
                        njac[j+1, 1, 1] = c3c4 * tmp1;
                        njac[j+1, 2, 1] = 0.0d;
                        njac[j+1, 3, 1] = 0.0d;
                        njac[j+1, 4, 1] = 0.0d;

                        njac[j+1, 0, 2] = -con43 * c3c4 * tmp2 * utmp[j, 3];
                        njac[j+1, 1, 2] = 0.0d;
                        njac[j+1, 2, 2] = con43 * c3c4 * tmp1;
                        njac[j+1, 3, 2] = 0.0d;
                        njac[j+1, 4, 2] = 0.0d;

                        njac[j+1, 0, 3] = -c3c4 * tmp2 * utmp[j, 4];
                        njac[j+1, 1, 3] = 0.0d;
                        njac[j+1, 2, 3] = 0.0d;
                        njac[j+1, 3, 3] = c3c4 * tmp1;
                        njac[j+1, 4, 3] = 0.0d;

                        njac[j+1, 0, 4] = -(c3c4 - c1345) * tmp3 * (pow2(utmp[j, 2])) -
                                        (con43 * c3c4 - c1345) * tmp3 * (pow2(utmp[j, 3])) - (c3c4 - c1345) * tmp3 * (pow2(utmp[j, 4])) - c1345 * tmp2 * utmp[j, 5];

                        njac[j+1, 1, 4] = (c3c4 - c1345) * tmp2 * utmp[j, 2];
                        njac[j+1, 2, 4] = (con43 * c3c4 - c1345) * tmp2 * utmp[j, 3];
                        njac[j+1, 3, 4] = (c3c4 - c1345) * tmp2 * utmp[j, 4];
                        njac[j+1, 4, 4] = (c1345) * tmp1;
                    }
                    //---------------------------------------------------------------------
                    //     now joacobians set, so form left hand side in y direction
                    //---------------------------------------------------------------------
                    for(j = start[c, 1]; j <= jsize - end[c, 1]; j++) {
                        tmp1 = dt * ty1;
                        tmp2 = dt * ty2;//lhsa[1,1,j]=-tmp2*fjac[1,1,j-1]-tmp1*njac[1,1,j-1]-tmp1*dy1;
                        lhsa[j, 0, 0] = -tmp2 * fjac[j, 0, 0] - tmp1 * njac[j, 0, 0] - tmp1 * dy1;
                        lhsa[j, 1, 0] = -tmp2 * fjac[j, 1, 0] - tmp1 * njac[j, 1, 0];//lhsa[1,2,j]=-tmp2*fjac[1,2,j-1]-tmp1*njac[1,2,j-1];
                        lhsa[j, 2, 0] = -tmp2 * fjac[j, 2, 0] - tmp1 * njac[j, 2, 0];
                        lhsa[j, 3, 0] = -tmp2 * fjac[j, 3, 0] - tmp1 * njac[j, 3, 0];
                        lhsa[j, 4, 0] = -tmp2 * fjac[j, 4, 0] - tmp1 * njac[j, 4, 0]; //Obs: fjac[2+j-1,5,1] - tmp1 * njac[2+j-1,5,1];

                        lhsa[j, 0, 1] = -tmp2 * fjac[j, 0, 1] - tmp1 * njac[j, 0, 1];
                        lhsa[j, 1, 1] = -tmp2 * fjac[j, 1, 1] - tmp1 * njac[j, 1, 1] - tmp1 * dy2;
                        lhsa[j, 2, 1] = -tmp2 * fjac[j, 2, 1] - tmp1 * njac[j, 2, 1];
                        lhsa[j, 3, 1] = -tmp2 * fjac[j, 3, 1] - tmp1 * njac[j, 3, 1];
                        lhsa[j, 4, 1] = -tmp2 * fjac[j, 4, 1] - tmp1 * njac[j, 4, 1];

                        lhsa[j, 0, 2] = -tmp2 * fjac[j, 0, 2] - tmp1 * njac[j, 0, 2];
                        lhsa[j, 1, 2] = -tmp2 * fjac[j, 1, 2] - tmp1 * njac[j, 1, 2];
                        lhsa[j, 2, 2] = -tmp2 * fjac[j, 2, 2] - tmp1 * njac[j, 2, 2] - tmp1 * dy3;
                        lhsa[j, 3, 2] = -tmp2 * fjac[j, 3, 2] - tmp1 * njac[j, 3, 2];
                        lhsa[j, 4, 2] = -tmp2 * fjac[j, 4, 2] - tmp1 * njac[j, 4, 2];

                        lhsa[j, 0, 3] = -tmp2 * fjac[j, 0, 3] - tmp1 * njac[j, 0, 3];
                        lhsa[j, 1, 3] = -tmp2 * fjac[j, 1, 3] - tmp1 * njac[j, 1, 3];
                        lhsa[j, 2, 3] = -tmp2 * fjac[j, 2, 3] - tmp1 * njac[j, 2, 3];
                        lhsa[j, 3, 3] = -tmp2 * fjac[j, 3, 3] - tmp1 * njac[j, 3, 3] - tmp1 * dy4;
                        lhsa[j, 4, 3] = -tmp2 * fjac[j, 4, 3] - tmp1 * njac[j, 4, 3];

                        lhsa[j, 0, 4] = -tmp2 * fjac[j, 0, 4] - tmp1 * njac[j, 0, 4];
                        lhsa[j, 1, 4] = -tmp2 * fjac[j, 1, 4] - tmp1 * njac[j, 1, 4];
                        lhsa[j, 2, 4] = -tmp2 * fjac[j, 2, 4] - tmp1 * njac[j, 2, 4];
                        lhsa[j, 3, 4] = -tmp2 * fjac[j, 3, 4] - tmp1 * njac[j, 3, 4];
                        lhsa[j, 4, 4] = -tmp2 * fjac[j, 4, 4] - tmp1 * njac[j, 4, 4] - tmp1 * dy5;

                        lhsb[j, 0, 0] = 1.0d + tmp1 * 2.0d * njac[j+1, 0, 0] + tmp1 * 2.0d * dy1;
                        lhsb[j, 1, 0] = tmp1 * 2.0d * njac[j+1, 1, 0];  //lhsb[1,2,j] = tmp1 * 2.0d * njac[1,2,j];
                        lhsb[j, 2, 0] = tmp1 * 2.0d * njac[j+1, 2, 0];  //lhsb[1,3,j] = tmp1 * 2.0d * njac[1,3,j];
                        lhsb[j, 3, 0] = tmp1 * 2.0d * njac[j+1, 3, 0];  //lhsb[1,4,j] = tmp1 * 2.0d * njac[1,4,j];
                        lhsb[j, 4, 0] = tmp1 * 2.0d * njac[j+1, 4, 0];  //lhsb[1,5,j] = tmp1 * 2.0d * njac[1,5,j];

                        lhsb[j, 0, 1] = tmp1 * 2.0d * njac[j+1, 0, 1];
                        lhsb[j, 1, 1] = 1.0d + tmp1 * 2.0d * njac[j+1, 1, 1] + tmp1 * 2.0d * dy2;
                        lhsb[j, 2, 1] = tmp1 * 2.0d * njac[j+1, 2, 1];
                        lhsb[j, 3, 1] = tmp1 * 2.0d * njac[j+1, 3, 1];
                        lhsb[j, 4, 1] = tmp1 * 2.0d * njac[j+1, 4, 1];

                        lhsb[j, 0, 2] = tmp1 * 2.0d * njac[j+1, 0, 2];
                        lhsb[j, 1, 2] = tmp1 * 2.0d * njac[j+1, 1, 2];
                        lhsb[j, 2, 2] = 1.0d + tmp1 * 2.0d * njac[j+1, 2, 2] + tmp1 * 2.0d * dy3;
                        lhsb[j, 3, 2] = tmp1 * 2.0d * njac[j+1, 3, 2];
                        lhsb[j, 4, 2] = tmp1 * 2.0d * njac[j+1, 4, 2];

                        lhsb[j, 0, 3] = tmp1 * 2.0d * njac[j+1, 0, 3];
                        lhsb[j, 1, 3] = tmp1 * 2.0d * njac[j+1, 1, 3];
                        lhsb[j, 2, 3] = tmp1 * 2.0d * njac[j+1, 2, 3];
                        lhsb[j, 3, 3] = 1.0d + tmp1 * 2.0d * njac[j+1, 3, 3] + tmp1 * 2.0d * dy4;
                        lhsb[j, 4, 3] = tmp1 * 2.0d * njac[j+1, 4, 3];

                        lhsb[j, 0, 4] = tmp1 * 2.0d * njac[j+1, 0, 4];
                        lhsb[j, 1, 4] = tmp1 * 2.0d * njac[j+1, 1, 4];
                        lhsb[j, 2, 4] = tmp1 * 2.0d * njac[j+1, 2, 4];
                        lhsb[j, 3, 4] = tmp1 * 2.0d * njac[j+1, 3, 4];
                        lhsb[j, 4, 4] = 1.0d + tmp1 * 2.0d * njac[j+1, 4, 4] + tmp1 * 2.0d * dy5;

                        lhsc[c, k, j, i, 0, 0] = tmp2 * fjac[j+2, 0, 0] - tmp1 * njac[j+2, 0, 0] - tmp1 * dy1;
                        lhsc[c, k, j, i, 1, 0] = tmp2 * fjac[j+2, 1, 0] - tmp1 * njac[j+2, 1, 0];//lhsc[1,2,i,j,k,c]=tmp2*fjac[1,2,j+1]-tmp1*njac[1,2,j+1];
                        lhsc[c, k, j, i, 2, 0] = tmp2 * fjac[j+2, 2, 0] - tmp1 * njac[j+2, 2, 0];//lhsc[1,3,i,j,k,c]=tmp2*fjac[1,3,j+1]-tmp1*njac[1,3,j+1];
                        lhsc[c, k, j, i, 3, 0] = tmp2 * fjac[j+2, 3, 0] - tmp1 * njac[j+2, 3, 0];//lhsc[1,4,i,j,k,c]=tmp2*fjac[1,4,j+1]-tmp1*njac[1,4,j+1];
                        lhsc[c, k, j, i, 4, 0] = tmp2 * fjac[j+2, 4, 0] - tmp1 * njac[j+2, 4, 0];//lhsc[1,5,i,j,k,c]=tmp2*fjac[1,5,j+1]-tmp1*njac[1,5,j+1];

                        lhsc[c, k, j, i, 0, 1] = tmp2 * fjac[j+2, 0, 1] - tmp1 * njac[j+2, 0, 1];
                        lhsc[c, k, j, i, 1, 1] = tmp2 * fjac[j+2, 1, 1] - tmp1 * njac[j+2, 1, 1] - tmp1 * dy2;
                        lhsc[c, k, j, i, 2, 1] = tmp2 * fjac[j+2, 2, 1] - tmp1 * njac[j+2, 2, 1];
                        lhsc[c, k, j, i, 3, 1] = tmp2 * fjac[j+2, 3, 1] - tmp1 * njac[j+2, 3, 1];
                        lhsc[c, k, j, i, 4, 1] = tmp2 * fjac[j+2, 4, 1] - tmp1 * njac[j+2, 4, 1];

                        lhsc[c, k, j, i, 0, 2] = tmp2 * fjac[j+2, 0, 2] - tmp1 * njac[j+2, 0, 2];
                        lhsc[c, k, j, i, 1, 2] = tmp2 * fjac[j+2, 1, 2] - tmp1 * njac[j+2, 1, 2];
                        lhsc[c, k, j, i, 2, 2] = tmp2 * fjac[j+2, 2, 2] - tmp1 * njac[j+2, 2, 2] - tmp1 * dy3;
                        lhsc[c, k, j, i, 3, 2] = tmp2 * fjac[j+2, 3, 2] - tmp1 * njac[j+2, 3, 2];
                        lhsc[c, k, j, i, 4, 2] = tmp2 * fjac[j+2, 4, 2] - tmp1 * njac[j+2, 4, 2];

                        lhsc[c, k, j, i, 0, 3] = tmp2 * fjac[j+2, 0, 3] - tmp1 * njac[j+2, 0, 3];
                        lhsc[c, k, j, i, 1, 3] = tmp2 * fjac[j+2, 1, 3] - tmp1 * njac[j+2, 1, 3];
                        lhsc[c, k, j, i, 2, 3] = tmp2 * fjac[j+2, 2, 3] - tmp1 * njac[j+2, 2, 3];
                        lhsc[c, k, j, i, 3, 3] = tmp2 * fjac[j+2, 3, 3] - tmp1 * njac[j+2, 3, 3] - tmp1 * dy4;
                        lhsc[c, k, j, i, 4, 3] = tmp2 * fjac[j+2, 4, 3] - tmp1 * njac[j+2, 4, 3];

                        lhsc[c, k, j, i, 0, 4] = tmp2 * fjac[j+2, 0, 4] - tmp1 * njac[j+2, 0, 4];
                        lhsc[c, k, j, i, 1, 4] = tmp2 * fjac[j+2, 1, 4] - tmp1 * njac[j+2, 1, 4];
                        lhsc[c, k, j, i, 2, 4] = tmp2 * fjac[j+2, 2, 4] - tmp1 * njac[j+2, 2, 4];
                        lhsc[c, k, j, i, 3, 4] = tmp2 * fjac[j+2, 3, 4] - tmp1 * njac[j+2, 3, 4];
                        lhsc[c, k, j, i, 4, 4] = tmp2 * fjac[j+2, 4, 4] - tmp1 * njac[j+2, 4, 4] - tmp1 * dy5;
                    }
                    //---------------------------------------------------------------------
                    //     outer most for(loops - sweeping in i direction
                    //---------------------------------------------------------------------
                    if(first == 1) {
                        //c---------------------------------------------------------------------
                        //c     multiply c[i,jstart,k] by b_inverse and copy back to c
                        //c     multiply rhs[jstart] by b_inverse[jstart] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvcrhs[ lhsb[1,1,jstart],lhsc[1,1,i,jstart,k,c], rhs[1,i,jstart,k,c] ];
                        //C#:           binvcrhs(lhsb[jstart+1, 1, 1], lhsc[c, k+1, jstart+1, i+1, 1, 1], rhs[c, k+1, jstart+1, i+1, 1]);
                        //binvcrhs(lhsb, lhsc, rhs, jstart-1, c+1, (k-1), jstart-1, (i-1), c+1, (k-1), jstart-1, (i-1));
                        binvcrhs(lhsb, lhsc, rhs, jstart, c, k, jstart, i, c, k, jstart, i);
                    }
                    //---------------------------------------------------------------------
                    //     begin inner most for(loop
                    //     for(all the elements of the cell unless last 
                    //---------------------------------------------------------------------
                    for(j = jstart + first; j <= jsize - last; j++) {
                        //c---------------------------------------------------------------------
                        //c     subtract A*lhs_vector[j-1] from lhs_vector[j]
                        //c     
                        //c     rhs[j] = rhs[j] - A*rhs[j-1]
                        //c---------------------------------------------------------------------
                        //Fortran: call matvec_sub[lhsa[1,1,j], rhs[1,i,j-1,k,c],rhs[1,i,j,k,c]];
                        //C#:           matvec_sub(lhsa[j+1,1,1], rhs[c,k+1,j,i+1,1],rhs[c,k+1,j+1,i+1,1]);
                        //matvec_sub(lhsa, rhs, rhs, (j-1), c+1, (k-1), j-2, (i-1), c+1, (k-1), (j-1), (i-1));
                        matvec_sub(lhsa, rhs, rhs, j, c, k, j-1, i, c, k, j, i);
                        //c---------------------------------------------------------------------
                        //c     B[j] = B[j] - C[j-1]*A[j]
                        //c---------------------------------------------------------------------
                        //call matmul_sub[lhsa[1,1,j],lhsc[1,1,i,j-1,k,c],lhsb[1,1,j]];
                        //C#:  matmul_sub(lhsa[j+1,1,1],lhsc[c,k+1,j-1+1,i+1,1,1],lhsb[j+1,1,1]);
                        //matmul_sub(lhsa, lhsc, lhsb, (j-1), c+1, (k-1), j-2, (i-1), (j-1));
                        matmul_sub(lhsa, lhsc, lhsb, j, c, k, j-1, i, j);
                        //c---------------------------------------------------------------------
                        //c     multiply c[i,j,k] by b_inverse and copy back to c
                        //c     multiply rhs[i,1,k] by b_inverse[i,1,k] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvcrhs[ lhsb[1,1,j],lhsc[1,1,i,j,k,c], rhs[1,i,j,k,c] ];
                        //C#:           binvcrhs( lhsb[j+1,1,1],lhsc[c,k+1,j+1,i+1,1,1], rhs[c,k+1,j+1,i+1,1] );
                        //binvcrhs(lhsb, lhsc, rhs, (j-1), c+1, (k-1), (j-1), (i-1), c+1, (k-1), (j-1), (i-1));
                        binvcrhs(lhsb, lhsc, rhs, j, c, k, j, i, c, k, j, i);
                    }
                    //c---------------------------------------------------------------------
                    //c     Now finish up special cases for last cell
                    //c---------------------------------------------------------------------
                    if(last == 1) {
                        //c---------------------------------------------------------------------
                        //c     rhs[jsize] = rhs[jsize] - A*rhs[jsize-1]
                        //c---------------------------------------------------------------------
                        //Fortran: call matvec_sub[lhsa[1,1,jsize], rhs[1,i,jsize-1,k,c],rhs[1,i,jsize,k,c]];
                        //C#:           matvec_sub(lhsa[jsize+1,1,1], rhs[c,k+1,jsize-1+1,i+1,1],rhs[c,k+1,jsize+1,i+1,1]);
                        //matvec_sub(lhsa, rhs, rhs, (jsize-1), c+1, (k-1), jsize-2, (i-1), c+1, (k-1), (jsize-1), (i-1));
                        matvec_sub(lhsa, rhs, rhs, jsize, c, k, jsize-1, i, c, k, jsize, i);
                        //c---------------------------------------------------------------------
                        //c     B[jsize] = B[jsize] - C[jsize-1]*A[jsize]
                        //c     call matmul_sub[aa,i,jsize,k,c,
                        //c     $              cc,i,jsize-1,k,c,bb,i,jsize,k,c]
                        //c---------------------------------------------------------------------
                        //Fortran: call matmul_sub[lhsa[1,1,jsize], lhsc[1,1,i,jsize-1,k,c], lhsb[1,1,jsize]];
                        //C#:           matmul_sub(lhsa[jsize+1,1,1], lhsc[c,k+1,jsize-1+1,i+1,1,1], lhsb[jsize+1,1,1]);
                        //matmul_sub(lhsa, lhsc, lhsb, (jsize-1), c+1, (k-1), jsize-2, (i-1), (jsize-1));
                        matmul_sub(lhsa, lhsc, lhsb, jsize, c, k, jsize-1, i, jsize);
                        //c---------------------------------------------------------------------
                        //c     multiply rhs[jsize] by b_inverse[jsize] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvrhs[ lhsb[1,1,jsize],rhs[1,i,jsize,k,c] ];
                        //C#:           binvrhs(lhsb[jsize+1,1,1],rhs[c,k+1,jsize+1,i+1,1] );
                        binvrhs(lhsb, rhs, jsize, c, k, jsize, i);
                    }
                }
            }
        }// start[ ksize isize jsize jstart istart kstart


        public void z_solve() {
            //---------------------------------------------------------------------
            //     Performs line solves in Z direction by first factoring
            //     the block-tridiagonal matrix into an upper triangular matrix, 
            //     and { performing back substitution to solve for the unknow
            //     vectors of each line.  
            //     
            //     Make sure we treat elements zero to cell_size in the direction
            //     of the sweep.
            //---------------------------------------------------------------------
            int c, stage, first, last, buffer_size; //int[] r_status[MPI_STATUS_SIZE];
            buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * (5 * 5 + 5);
            //MPI.Request[] recv_id = new MPI.Request[1];
            //MPI.Request[] send_id = new MPI.Request[1];
            MPI.Request[] requests = new MPI.Request[2] { null, null };
            double[] out_buffer_z = new double[buffer_size];
            //---------------------------------------------------------------------
            //     in our terminology stage is the number of the cell in the y-direction
            //     i.e. stage = 1 means the start of the line stage=ncells means end
            //---------------------------------------------------------------------
            for(stage = 0; stage < ncells; stage++) {
                c = slice[stage, 2];
                //isize = cell_size[c, 0] - 1;
                //jsize = cell_size[c, 1] - 1;
                //ksize = cell_size[c, 2] - 1;
                //---------------------------------------------------------------------
                //     set last-cell flag
                //---------------------------------------------------------------------
                if(stage == ncells-1) {
                    last = 1;
                }
                else {
                    last = 0;
                }

                if(stage == 0) {
                    //---------------------------------------------------------------------
                    //     This is the first cell, so solve without receiving data
                    //---------------------------------------------------------------------
                    first = 1;
                    z_solve_cell(first, last, c);  //call z_solve_cell[first,last,c];
                }
                else {
                    //---------------------------------------------------------------------
                    //     Not the first cell of this line, so receive info from
                    //     processor working on preceeding cell
                    //---------------------------------------------------------------------
                    first = 0;
                    int ip = cell_coord[c, 0];
                    int jp = cell_coord[c, 1];
                    requests[0] = comm_solve.ImmediateReceive<double>(predecessor[2], BOTTOM+ip+jp*ncells, out_buffer_z);//z_receive_solve_info(out_buffer_z, requests, c); //recv_id, c);//call z_receive_solve_info[recv_id,c];
                    //  c---------------------------------------------------------------------
                    //  c     overlap computations and communications
                    //  c---------------------------------------------------------------------
                    //      c      call lhsz[c]
                    //  c---------------------------------------------------------------------
                    //  c     wait for completion
                    //  c---------------------------------------------------------------------
                    requests[1].Wait(); //send_id[0].Wait();
                    requests[0].Wait(); //recv_id[0].Wait();
                    //Fortran: call mpi_wait[send_id,r_status,error];
                    //Fortran: call mpi_wait[recv_id,r_status,error];
                    //  c---------------------------------------------------------------------
                    //  c     install C'[kstart+1] and rhs'[kstart+1] to be used in this cell
                    //  c---------------------------------------------------------------------
                    z_unpack_solve_info(out_buffer_z, c);
                    z_solve_cell(first, last, c);  //call z_solve_cell[first,last,c];
                }
                if(last == 0) {
                    int ip = cell_coord[c,0];
                    int jp = cell_coord[c,1];                    
                    double[] in_buffer_z = new double[buffer_size];//buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * (5 * 5 + 5);
                    z_pack_solve_info(in_buffer_z, c); //send_id, c);  //call z_send_solve_info[send_id,c];
                    requests[1] = comm_solve.ImmediateSend<double>(in_buffer_z, successor[2], BOTTOM+ip+jp*ncells);                    
                }
            }
            out_buffer_z = null;
            buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * 5;
            out_buffer_z = new double[buffer_size];
            //---------------------------------------------------------------------
            //     now perform backsubstitution in reverse direction
            //---------------------------------------------------------------------
            for(stage = ncells-1; stage >= 0; stage--) {  //for(stage = ncells, 1, -1;
                c = slice[stage, 2];
                first = 0;
                last = 0;
                if(stage == 0)
                    first = 1;
                if(stage == ncells-1) {
                    last = 1;
                    //---------------------------------------------------------------------
                    //     last cell, so perform back substitute without waiting
                    //---------------------------------------------------------------------
                    z_backsubstitute(first, last, c); //call z_backsubstitute[first, last,c];
                }
                else {
                    int ip = cell_coord[c, 0];
                    int jp = cell_coord[c, 1];
                    requests[0] = comm_solve.ImmediateReceive<double>(successor[2], TOP+ip+jp*ncells, out_buffer_z);// z_receive_backsub_info(out_buffer_z, requests, c); //recv_id, c);            //call z_receive_backsub_info[recv_id,c];
                    
                    requests[1].Wait(); //send_id[0].Wait();
                    requests[0].Wait(); //recv_id[0].Wait();
                    // Fortran: call mpi_wait[send_id,r_status,error]; 
                    // Fortran: call mpi_wait[recv_id,r_status,error];
                    z_unpack_backsub_info(out_buffer_z, c);
                    z_backsubstitute(first, last, c); //call z_backsubstitute[first,last,c];
                }
                if(first == 0) {
                    int ip = cell_coord[c,0];
                    int jp = cell_coord[c,1];                    
                    double[] in_buffer_z = new double[buffer_size];//buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * 5;
                    z_pack_backsub_info(in_buffer_z, c); //send_id, c);  //call z_send_backsub_info[send_id,c];
                    requests[1] = comm_solve.ImmediateSend<double>(in_buffer_z, predecessor[2], TOP+ip+jp*ncells);
                }
            }
        }

        public void z_unpack_solve_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack C'[-1] and rhs'[-1] for
            //     all i and j
            //---------------------------------------------------------------------
            int i, j, m, n, ptr, kstart;
            kstart = 2;
            ptr = 0;
            for(j = 0; j <= JMAX - 1; j++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(m = 1; m <= 5; m++) {
                        for(n = 0; n < 5; n++) {  //lhsc[m,n,i,j,kstart-1,c] = out_buffer[ptr+n];
                            lhsc[c, kstart-1, j+2, i+2, n, m-1] = out_buffer_x[ptr + n];
                        }
                        ptr = ptr + 5;
                    }
                    for(n = 0; n < 5; n++) { //rhs[n,i,j,kstart-1,c] = out_buffer[ptr+n];
                        rhs[c, kstart-1, j+2, i+2, n] = out_buffer_x[ptr + n];
                    }
                    ptr = ptr + 5;
                }
            }
        }

        public void z_pack_solve_info(double[] in_buffer_z, int c) {
            //---------------------------------------------------------------------
            //     pack up and send C'[kend] and rhs'[kend] for
            //     all i and j
            //---------------------------------------------------------------------
            int i, j, m, n, ksize, ptr;
            ksize = cell_size[c, 2] + 1;
            //c---------------------------------------------------------------------
            //c     pack up buffer
            //c---------------------------------------------------------------------
            ptr = 0;
            for(j = 0; j <= JMAX - 1; j++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(m = 1; m <= 5; m++) {
                        for(n = 0; n < 5; n++) {  //in_buffer[ptr+n] = lhsc[m,n,i,j,ksize,c]
                            in_buffer_z[ptr + n] = lhsc[c, ksize, j+2, i+2, n, m-1];
                        }
                        ptr = ptr + 5;
                    }
                    for(n = 0; n < 5; n++) {  //in_buffer[ptr+n] = rhs[n,i,j,ksize,c];
                        in_buffer_z[ptr + n] = rhs[c, ksize, j+2, i+2, n];
                    }
                    ptr = ptr + 5;
                }
            }
            //---------------------------------------------------------------------
            //     send buffer 
            //---------------------------------------------------------------------
            //call mpi_isend[in_buffer, buffer_size, dp_type, successor[3], BOTTOM+ip+jp*NCELLS, comm_solve, send_id,error];
            //send_id[0] =
        }

        public void z_pack_backsub_info(double[] in_buffer_z, int c) {
            //---------------------------------------------------------------------
            //     pack up and send U[jstart] for all i and j
            //---------------------------------------------------------------------
            int i, j, n, ptr, kstart;
            //---------------------------------------------------------------------
            //     Send element 0 to previous processor
            //---------------------------------------------------------------------
            kstart = 2;
            ptr = 0;
            for(j = 0; j <= JMAX - 1; j++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(n = 0; n < 5; n++) {  //in_buffer[ptr+n] = rhs[n,i,j,kstart,c];
                        in_buffer_z[ptr + n] = rhs[c, kstart, j+2, i+2, n];
                    }
                    ptr = ptr + 5;
                }
            }
            //call mpi_isend[in_buffer, buffer_size, dp_type, predecessor[3], TOP+ip+jp*NCELLS, comm_solve, send_id,error];
            //send_id[0] =
        }

        public void z_unpack_backsub_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack U[ksize] for all i and j
            //---------------------------------------------------------------------
            int i, j, n, ptr;
            ptr = 0;
            for(j = 0; j <= JMAX - 1; j++) {
                for(i = 0; i <= IMAX - 1; i++) {
                    for(n = 0; n < 5; n++) {  //backsub_info[n,i,j,c] = out_buffer[ptr+n];
                        backsub_info[c, j+2, i+2, n] = out_buffer_x[ptr + n];
                    }
                    ptr = ptr + 5;
                }
            }
        }

        public void z_backsubstitute(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     back solve: if last cell, { generate U[ksize]=rhs[ksize]
            //     } else { assume U[ksize] is loaded in un pack backsub_info
            //     so just use it
            //     after call u[kstart] will be sent to next cell
            //---------------------------------------------------------------------
            int i, k;
            int m, n, j, jsize, isize, ksize, kstart;

            kstart = 2;
            isize = cell_size[c, 0] - end[c, 0] + 1;
            jsize = cell_size[c, 1] - end[c, 1] + 1;
            ksize = cell_size[c, 2] + 1;
            if(last == 0) {
                for(j = start[c, 1]; j <= jsize; j++) {
                    for(i = start[c, 0]; i <= isize; i++) {
                        //---------------------------------------------------------------------
                        //     U[jsize] uses info from previous cell if not last cell
                        //---------------------------------------------------------------------
                        for(m = 0; m < 5; m++) {
                            for(n = 0; n < 5; n++) { //rhs[m,i,j,ksize,c] = rhs[m,i,j,ksize,c] - lhsc[m,n,i,j,ksize,c]*backsub_info[n,i,j,c]
                                rhs[c, ksize, j, i, m] = rhs[c, ksize, j, i, m] - lhsc[c, ksize, j, i, n, m] * backsub_info[c, j, i, n];
                            }
                        }
                    }
                }
            }
            for(k = ksize-1; k >= kstart; k--) {  //for(k=ksize-1,kstart,-1;
                for(j = start[c, 1]; j <= jsize; j++) {
                    for(i = start[c, 0]; i <= isize; i++) {
                        for(m = 0; m < 5; m++) {
                            for(n = 0; n < 5; n++) {  //rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - lhsc[m,n,i,j,k,c]*rhs[n,i,j,k+1,c];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - lhsc[c, k, j, i, n, m] * rhs[c, k+1, j, i, n];
                            }
                        }
                    }
                }
            }
        }//start[ ksize jsize isize istart jstart kstart

        public void z_solve_cell(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     performs guaussian elimination on this cell.
            //     
            //     assumes that unpacking routines for non-first cells 
            //     preload C' and rhs' from previous cell.
            //     
            //     assumed send happens outside this routine, but that
            //     c'[KMAX] and rhs'[KMAX] will be sent to next cell.
            //---------------------------------------------------------------------
            int i, j, k, isize, ksize, jsize, kstart;
            double tmp1, tmp2, tmp3;
            double[,] utmp = new double[KMAX + 4, 7];   //double utmp[6,-2:KMAX+1];
            kstart = 2;
            isize = cell_size[c, 0] - end[c, 0] + 1;
            jsize = cell_size[c, 1] - end[c, 1] + 1;
            ksize = cell_size[c, 2] + 1;

            lhsabinit(lhsa, lhsb, ksize); //call lhsabinit[lhsa, lhsb, ksize];

            for(j = start[c, 1]; j <= jsize; j++) {
                for(i = start[c, 0]; i <= isize; i++) {
                    //---------------------------------------------------------------------
                    //     This function computes the left hand side for the three z-factors   
                    //---------------------------------------------------------------------
                    //---------------------------------------------------------------------
                    //     Compute the indices for storing the block-diagonal matrix;
                    //     determine c [labeled f] and s jacobians for cell c
                    //---------------------------------------------------------------------
                    for(k = start[c, 2]-1; k <= cell_size[c, 2] - end[c, 2]+2; k++) {
                        utmp[k, 1] = 1.0d / u[c, k, j, i, 0];  //utmp[1,k] = 1.0d0 / u[1,i,j,k,c]
                        utmp[k, 2] = u[c, k, j, i, 1];  //utmp[2,k] =         u[2,i,j,k,c]
                        utmp[k, 3] = u[c, k, j, i, 2];  //utmp[3,k] =         u[3,i,j,k,c]
                        utmp[k, 4] = u[c, k, j, i, 3];  //utmp[4,k] =         u[4,i,j,k,c]
                        utmp[k, 5] = u[c, k, j, i, 4];  //utmp[5,k] =         u[5,i,j,k,c]
                        utmp[k, 6] = qs[c, k, j, i];  //utmp[6,k] =          qs[i,j,k,c]
                    }
                    for(k = start[c, 2]-1; k <= cell_size[c, 2] - end[c, 2]+2; k++) {

                        tmp1 = utmp[k, 1];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        fjac[k+1, 0, 0] = 0.0d;   //fjac[1,1,k] = 0.0d+00
                        fjac[k+1, 1, 0] = 0.0d;
                        fjac[k+1, 2, 0] = 0.0d;
                        fjac[k+1, 3, 0] = 1.0d;
                        fjac[k+1, 4, 0] = 0.0d;

                        fjac[k+1, 0, 1] = -(utmp[k, 2] * utmp[k, 4]) * tmp2;
                        fjac[k+1, 1, 1] = utmp[k, 4] * tmp1;
                        fjac[k+1, 2, 1] = 0.0d;
                        fjac[k+1, 3, 1] = utmp[k, 2] * tmp1;
                        fjac[k+1, 4, 1] = 0.0d;

                        fjac[k+1, 0, 2] = -(utmp[k, 3] * utmp[k, 4]) * tmp2;
                        fjac[k+1, 1, 2] = 0.0d;
                        fjac[k+1, 2, 2] = utmp[k, 4] * tmp1;
                        fjac[k+1, 3, 2] = utmp[k, 3] * tmp1;
                        fjac[k+1, 4, 2] = 0.0d;

                        fjac[k+1, 0, 3] = -(utmp[k, 4] * utmp[k, 4] * tmp2) + c2 * utmp[k, 6];
                        fjac[k+1, 1, 3] = -c2 * utmp[k, 2] * tmp1;
                        fjac[k+1, 2, 3] = -c2 * utmp[k, 3] * tmp1;
                        fjac[k+1, 3, 3] = (2.0d - c2) * utmp[k, 4] * tmp1;
                        fjac[k+1, 4, 3] = c2;

                        fjac[k+1, 0, 4] = (c2 * 2.0d * utmp[k, 6] - c1 * (utmp[k, 5] * tmp1)) * (utmp[k, 4] * tmp1);
                        fjac[k+1, 1, 4] = -c2 * (utmp[k, 2] * utmp[k, 4]) * tmp2;
                        fjac[k+1, 2, 4] = -c2 * (utmp[k, 3] * utmp[k, 4]) * tmp2;
                        fjac[k+1, 3, 4] = c1 * (utmp[k, 5] * tmp1) - c2 * (utmp[k, 6] + utmp[k, 4] * utmp[k, 4] * tmp2);
                        fjac[k+1, 4, 4] = c1 * utmp[k, 4] * tmp1;

                        njac[k+1, 0, 0] = 0.0d;
                        njac[k+1, 1, 0] = 0.0d;
                        njac[k+1, 2, 0] = 0.0d;
                        njac[k+1, 3, 0] = 0.0d;
                        njac[k+1, 4, 0] = 0.0d;

                        njac[k+1, 0, 1] = -c3c4 * tmp2 * utmp[k, 2];
                        njac[k+1, 1, 1] = c3c4 * tmp1;
                        njac[k+1, 2, 1] = 0.0d;
                        njac[k+1, 3, 1] = 0.0d;
                        njac[k+1, 4, 1] = 0.0d;

                        njac[k+1, 0, 2] = -c3c4 * tmp2 * utmp[k, 3];
                        njac[k+1, 1, 2] = 0.0d;
                        njac[k+1, 2, 2] = c3c4 * tmp1;
                        njac[k+1, 3, 2] = 0.0d;
                        njac[k+1, 4, 2] = 0.0d;

                        njac[k+1, 0, 3] = -con43 * c3c4 * tmp2 * utmp[k, 4];
                        njac[k+1, 1, 3] = 0.0d;
                        njac[k+1, 2, 3] = 0.0d;
                        njac[k+1, 3, 3] = con43 * c3 * c4 * tmp1;
                        njac[k+1, 4, 3] = 0.0d;

                        njac[k+1, 0, 4] = -(c3c4 - c1345) * tmp3 * (pow2(utmp[k, 2])) -
                                        (c3c4 - c1345) * tmp3 * (pow2(utmp[k, 3])) - (con43 * c3c4 - c1345) * tmp3 * (pow2(utmp[k, 4])) - c1345 * tmp2 * utmp[k, 5];

                        njac[k+1, 1, 4] = (c3c4 - c1345) * tmp2 * utmp[k, 2];
                        njac[k+1, 2, 4] = (c3c4 - c1345) * tmp2 * utmp[k, 3];
                        njac[k+1, 3, 4] = (con43 * c3c4 - c1345) * tmp2 * utmp[k, 4];
                        njac[k+1, 4, 4] = (c1345) * tmp1;
                    }
                    //---------------------------------------------------------------------
                    //     now joacobians set, so form left hand side in z direction
                    //---------------------------------------------------------------------
                    for(k = start[c, 2]; k <= ksize - end[c, 2]; k++) {
                        tmp1 = dt * tz1;
                        tmp2 = dt * tz2;
                        lhsa[k, 0, 0] = -tmp2 * fjac[k, 0, 0] - tmp1 * njac[k, 0, 0] - tmp1 * dz1;//lhsa[1,1,k]=-tmp2*fjac[1,1,k-1]-tmp1*njac[1,1,k-1]
                        lhsa[k, 1, 0] = -tmp2 * fjac[k, 1, 0] - tmp1 * njac[k, 1, 0];    //fjac[k-1+2,2,1] njac[k-1+2,2,1]
                        lhsa[k, 2, 0] = -tmp2 * fjac[k, 2, 0] - tmp1 * njac[k, 2, 0];
                        lhsa[k, 3, 0] = -tmp2 * fjac[k, 3, 0] - tmp1 * njac[k, 3, 0];
                        lhsa[k, 4, 0] = -tmp2 * fjac[k, 4, 0] - tmp1 * njac[k, 4, 0];

                        lhsa[k, 0, 1] = -tmp2 * fjac[k, 0, 1] - tmp1 * njac[k, 0, 1];
                        lhsa[k, 1, 1] = -tmp2 * fjac[k, 1, 1] - tmp1 * njac[k, 1, 1] - tmp1 * dz2;
                        lhsa[k, 2, 1] = -tmp2 * fjac[k, 2, 1] - tmp1 * njac[k, 2, 1];
                        lhsa[k, 3, 1] = -tmp2 * fjac[k, 3, 1] - tmp1 * njac[k, 3, 1];
                        lhsa[k, 4, 1] = -tmp2 * fjac[k, 4, 1] - tmp1 * njac[k, 4, 1];

                        lhsa[k, 0, 2] = -tmp2 * fjac[k, 0, 2] - tmp1 * njac[k, 0, 2];
                        lhsa[k, 1, 2] = -tmp2 * fjac[k, 1, 2] - tmp1 * njac[k, 1, 2];
                        lhsa[k, 2, 2] = -tmp2 * fjac[k, 2, 2] - tmp1 * njac[k, 2, 2] - tmp1 * dz3;
                        lhsa[k, 3, 2] = -tmp2 * fjac[k, 3, 2] - tmp1 * njac[k, 3, 2];
                        lhsa[k, 4, 2] = -tmp2 * fjac[k, 4, 2] - tmp1 * njac[k, 4, 2];

                        lhsa[k, 0, 3] = -tmp2 * fjac[k, 0, 3] - tmp1 * njac[k, 0, 3];
                        lhsa[k, 1, 3] = -tmp2 * fjac[k, 1, 3] - tmp1 * njac[k, 1, 3];
                        lhsa[k, 2, 3] = -tmp2 * fjac[k, 2, 3] - tmp1 * njac[k, 2, 3];
                        lhsa[k, 3, 3] = -tmp2 * fjac[k, 3, 3] - tmp1 * njac[k, 3, 3] - tmp1 * dz4;
                        lhsa[k, 4, 3] = -tmp2 * fjac[k, 4, 3] - tmp1 * njac[k, 4, 3];

                        lhsa[k, 0, 4] = -tmp2 * fjac[k, 0, 4] - tmp1 * njac[k, 0, 4];
                        lhsa[k, 1, 4] = -tmp2 * fjac[k, 1, 4] - tmp1 * njac[k, 1, 4];
                        lhsa[k, 2, 4] = -tmp2 * fjac[k, 2, 4] - tmp1 * njac[k, 2, 4];
                        lhsa[k, 3, 4] = -tmp2 * fjac[k, 3, 4] - tmp1 * njac[k, 3, 4];
                        lhsa[k, 4, 4] = -tmp2 * fjac[k, 4, 4] - tmp1 * njac[k, 4, 4] - tmp1 * dz5;

                        lhsb[k, 0, 0] = 1.0d + tmp1 * 2.0d * njac[k+1, 0, 0] + tmp1 * 2.0d * dz1; //lhsb[1,1,k]=1.0d+tmp1*2.0d+00*njac[1,1,k]
                        lhsb[k, 1, 0] = tmp1 * 2.0d * njac[k+1, 1, 0];
                        lhsb[k, 2, 0] = tmp1 * 2.0d * njac[k+1, 2, 0];
                        lhsb[k, 3, 0] = tmp1 * 2.0d * njac[k+1, 3, 0];
                        lhsb[k, 4, 0] = tmp1 * 2.0d * njac[k+1, 4, 0];

                        lhsb[k, 0, 1] = tmp1 * 2.0d * njac[k+1, 0, 1];
                        lhsb[k, 1, 1] = 1.0d + tmp1 * 2.0d * njac[k+1, 1, 1] + tmp1 * 2.0d * dz2;
                        lhsb[k, 2, 1] = tmp1 * 2.0d * njac[k+1, 2, 1];
                        lhsb[k, 3, 1] = tmp1 * 2.0d * njac[k+1, 3, 1];
                        lhsb[k, 4, 1] = tmp1 * 2.0d * njac[k+1, 4, 1];

                        lhsb[k, 0, 2] = tmp1 * 2.0d * njac[k+1, 0, 2];
                        lhsb[k, 1, 2] = tmp1 * 2.0d * njac[k+1, 1, 2];
                        lhsb[k, 2, 2] = 1.0d + tmp1 * 2.0d * njac[k+1, 2, 2] + tmp1 * 2.0d * dz3;
                        lhsb[k, 3, 2] = tmp1 * 2.0d * njac[k+1, 3, 2];
                        lhsb[k, 4, 2] = tmp1 * 2.0d * njac[k+1, 4, 2];

                        lhsb[k, 0, 3] = tmp1 * 2.0d * njac[k+1, 0, 3];
                        lhsb[k, 1, 3] = tmp1 * 2.0d * njac[k+1, 1, 3];
                        lhsb[k, 2, 3] = tmp1 * 2.0d * njac[k+1, 2, 3];
                        lhsb[k, 3, 3] = 1.0d + tmp1 * 2.0d * njac[k+1, 3, 3] + tmp1 * 2.0d * dz4;
                        lhsb[k, 4, 3] = tmp1 * 2.0d * njac[k+1, 4, 3];

                        lhsb[k, 0, 4] = tmp1 * 2.0d * njac[k+1, 0, 4];
                        lhsb[k, 1, 4] = tmp1 * 2.0d * njac[k+1, 1, 4];
                        lhsb[k, 2, 4] = tmp1 * 2.0d * njac[k+1, 2, 4];
                        lhsb[k, 3, 4] = tmp1 * 2.0d * njac[k+1, 3, 4];
                        lhsb[k, 4, 4] = 1.0d + tmp1 * 2.0d * njac[k+1, 4, 4] + tmp1 * 2.0d * dz5;

                        lhsc[c, k, j, i, 0, 0] = tmp2 * fjac[k+2, 0, 0] - tmp1 * njac[k+2, 0, 0] - tmp1 * dz1;//lhsc[1,1,i,j,k,c] fjac[1,1,k+1] njac[1,1,k+1]
                        lhsc[c, k, j, i, 1, 0] = tmp2 * fjac[k+2, 1, 0] - tmp1 * njac[k+2, 1, 0];//fjac[k+1+2,2,1] njac[k+1+2,2,1]
                        lhsc[c, k, j, i, 2, 0] = tmp2 * fjac[k+2, 2, 0] - tmp1 * njac[k+2, 2, 0];
                        lhsc[c, k, j, i, 3, 0] = tmp2 * fjac[k+2, 3, 0] - tmp1 * njac[k+2, 3, 0];
                        lhsc[c, k, j, i, 4, 0] = tmp2 * fjac[k+2, 4, 0] - tmp1 * njac[k+2, 4, 0];

                        lhsc[c, k, j, i, 0, 1] = tmp2 * fjac[k+2, 0, 1] - tmp1 * njac[k+2, 0, 1];
                        lhsc[c, k, j, i, 1, 1] = tmp2 * fjac[k+2, 1, 1] - tmp1 * njac[k+2, 1, 1] - tmp1 * dz2;
                        lhsc[c, k, j, i, 2, 1] = tmp2 * fjac[k+2, 2, 1] - tmp1 * njac[k+2, 2, 1];
                        lhsc[c, k, j, i, 3, 1] = tmp2 * fjac[k+2, 3, 1] - tmp1 * njac[k+2, 3, 1];
                        lhsc[c, k, j, i, 4, 1] = tmp2 * fjac[k+2, 4, 1] - tmp1 * njac[k+2, 4, 1];

                        lhsc[c, k, j, i, 0, 2] = tmp2 * fjac[k+2, 0, 2] - tmp1 * njac[k+2, 0, 2];
                        lhsc[c, k, j, i, 1, 2] = tmp2 * fjac[k+2, 1, 2] - tmp1 * njac[k+2, 1, 2];
                        lhsc[c, k, j, i, 2, 2] = tmp2 * fjac[k+2, 2, 2] - tmp1 * njac[k+2, 2, 2] - tmp1 * dz3;
                        lhsc[c, k, j, i, 3, 2] = tmp2 * fjac[k+2, 3, 2] - tmp1 * njac[k+2, 3, 2];
                        lhsc[c, k, j, i, 4, 2] = tmp2 * fjac[k+2, 4, 2] - tmp1 * njac[k+2, 4, 2];

                        lhsc[c, k, j, i, 0, 3] = tmp2 * fjac[k+2, 0, 3] - tmp1 * njac[k+2, 0, 3];
                        lhsc[c, k, j, i, 1, 3] = tmp2 * fjac[k+2, 1, 3] - tmp1 * njac[k+2, 1, 3];
                        lhsc[c, k, j, i, 2, 3] = tmp2 * fjac[k+2, 2, 3] - tmp1 * njac[k+2, 2, 3];
                        lhsc[c, k, j, i, 3, 3] = tmp2 * fjac[k+2, 3, 3] - tmp1 * njac[k+2, 3, 3] - tmp1 * dz4;
                        lhsc[c, k, j, i, 4, 3] = tmp2 * fjac[k+2, 4, 3] - tmp1 * njac[k+2, 4, 3];

                        lhsc[c, k, j, i, 0, 4] = tmp2 * fjac[k+2, 0, 4] - tmp1 * njac[k+2, 0, 4];
                        lhsc[c, k, j, i, 1, 4] = tmp2 * fjac[k+2, 1, 4] - tmp1 * njac[k+2, 1, 4];
                        lhsc[c, k, j, i, 2, 4] = tmp2 * fjac[k+2, 2, 4] - tmp1 * njac[k+2, 2, 4];
                        lhsc[c, k, j, i, 3, 4] = tmp2 * fjac[k+2, 3, 4] - tmp1 * njac[k+2, 3, 4];
                        lhsc[c, k, j, i, 4, 4] = tmp2 * fjac[k+2, 4, 4] - tmp1 * njac[k+2, 4, 4] - tmp1 * dz5;
                    }
                    //---------------------------------------------------------------------
                    //     outer most for(loops - sweeping in i direction
                    //---------------------------------------------------------------------
                    if(first == 1) {
                        //---------------------------------------------------------------------
                        //     multiply c[i,j,kstart] by b_inverse and copy back to c
                        //     multiply rhs[kstart] by b_inverse[kstart] and copy to rhs
                        //---------------------------------------------------------------------
                        //Fortran: call binvcrhs[ lhsb[1,1,kstart], lhsc[1,1,i,j,kstart,c], rhs[1,i,j,kstart,c] ];
                        //C#:           binvcrhs( lhsb[kstart+1,1,1], lhsc[c,kstart+1,j+1,i+1,1,1], rhs[c,kstart+1,j+1,i+1,1] );
                        //binvcrhs(lhsb, lhsc, rhs, kstart-1, c+1, kstart-1, (j-1), (i-1), c+1, kstart-1, (j-1), (i-1));
                        binvcrhs(lhsb, lhsc, rhs, kstart, c, kstart, j, i, c, kstart, j, i);
                    }
                    //c---------------------------------------------------------------------
                    //c     begin inner most for(loop
                    //c     for(all the elements of the cell unless last 
                    //c---------------------------------------------------------------------
                    for(k = kstart + first; k <= ksize - last; k++) {
                        //c---------------------------------------------------------------------
                        //c     subtract A*lhs_vector[k-1] from lhs_vector[k]
                        //c     
                        //c     rhs[k] = rhs[k] - A*rhs[k-1]
                        //c---------------------------------------------------------------------
                        //Fortran: call matvec_sub[lhsa[1,1,k], rhs[1,i,j,k-1,c],rhs[1,i,j,k,c]];
                        //C#: matvec_sub(lhsa[k+1,1,1], rhs[c,k-1+1,j+1,i+1,1],rhs[c,k+1,j+1,i+1,1]);
                        //matvec_sub(lhsa, rhs, rhs, (k-1), c+1, k-2, (j-1), (i-1), c+1, (k-1), (j-1), (i-1));
                        matvec_sub(lhsa, rhs, rhs, k, c, k-1, j, i, c, k, j, i);
                        //c---------------------------------------------------------------------
                        //c     B[k] = B[k] - C[k-1]*A[k]
                        //c     call matmul_sub[aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k,c]
                        //c---------------------------------------------------------------------
                        //Fortran: call matmul_sub[lhsa[1,1,k], lhsc[1,1,i,j,k-1,c], lhsb[1,1,k]];
                        //C#:           matmul_sub(lhsa[k+1,1,1], lhsc[c,k-1+1,j+1,i+1,1,1], lhsb[k+1,1,1]);
                        //matmul_sub(lhsa, lhsc, lhsb, (k-1), c+1, k-2, (j-1), (i-1), (k-1));
                        matmul_sub(lhsa, lhsc, lhsb, k, c, k-1, j, i, k);
                        //c---------------------------------------------------------------------
                        //c     multiply c[i,j,k] by b_inverse and copy back to c
                        //c     multiply rhs[i,j,1] by b_inverse[i,j,1] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvcrhs[ lhsb[1,1,k], lhsc[1,1,i,j,k,c], rhs[1,i,j,k,c] ];
                        //C#:           binvcrhs( lhsb[k+1,1,1],    lhsc[c,k+1,j+1,i+1,1,1],    rhs[c,k+1,j+1,i+1,1] );
                        //binvcrhs(lhsb, lhsc, rhs, (k-1), c+1, (k-1), (j-1), (i-1), c+1, (k-1), (j-1), (i-1));
                        binvcrhs(lhsb, lhsc, rhs, k, c, k, j, i, c, k, j, i);
                    }
                    //c---------------------------------------------------------------------
                    //c     Now finish up special cases for last cell
                    //c---------------------------------------------------------------------
                    if(last == 1) {
                        //c---------------------------------------------------------------------
                        //c     rhs[ksize] = rhs[ksize] - A*rhs[ksize-1]
                        //c---------------------------------------------------------------------
                        //Fortran: call matvec_sub[lhsa[1,1,ksize], rhs[1,i,j,ksize-1,c],rhs[1,i,j,ksize,c]];
                        //C#:           matvec_sub(lhsa[ksize+1,1,1], rhs[c,ksize-1+1,j+1,i+1,1],rhs[c,ksize+1,j+1,i+1,1]);
                        //matvec_sub(lhsa, rhs, rhs, (ksize-1), c+1, (ksize-2), (j-1), (i-1), c+1, (ksize-1), (j-1), (i-1));
                        matvec_sub(lhsa, rhs, rhs, ksize, c, ksize-1, j, i, c, ksize, j, i);
                        //c---------------------------------------------------------------------
                        //c     B[ksize] = B[ksize] - C[ksize-1]*A[ksize]
                        //c     call matmul_sub[aa,i,j,ksize,c,
                        //c     $              cc,i,j,ksize-1,c,bb,i,j,ksize,c]
                        //c---------------------------------------------------------------------
                        //Fortran: call matmul_sub[lhsa[1,1,ksize],lhsc[1,1,i,j,ksize-1,c],lhsb[1,1,ksize]];
                        //C#:           matmul_sub[lhsa[ksize+1,1,1],lhsc[c,ksize-1+1,j+1,i+1,1,1],lhsb[ksize+1,1,1]];
                        //matmul_sub(lhsa, lhsc, lhsb, (ksize-1), c+1, (ksize-2), (j-1), (i-1), (ksize-1));
                        matmul_sub(lhsa, lhsc, lhsb, ksize, c, ksize-1, j, i, ksize);
                        //c---------------------------------------------------------------------
                        //c     multiply rhs[ksize] by b_inverse[ksize] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvrhs[ lhsb[1,1,ksize], rhs[1,i,j,ksize,c] ];
                        //C#:           binvrhs(lhsb[ksize+1,1,1], rhs[c,ksize+1,j+1,i+1,1] );
                        //binvrhs(lhsb, rhs, (ksize-1), c+1, (ksize-1), (j-1), (i-1));
                        binvrhs(lhsb, rhs, ksize, c, ksize, j, i);
                    }
                }
            }
        }//start[ ksize isize jsize kstart istart jstart


        public void binvcrhs(double[, ,] lhs, double[, , , , ,] c, double[, , , ,] r, int l1, int c1, int c2, int c3, int c4, int r1, int r2, int r3, int r4) {
            double pivot, coeff; //dimension lhs[5,5]; //double c[5,5], r[5];
            pivot = 1.00d / lhs[l1, 0, 0];
            lhs[l1, 1, 0] = lhs[l1, 1, 0] * pivot;
            lhs[l1, 2, 0] = lhs[l1, 2, 0] * pivot;
            lhs[l1, 3, 0] = lhs[l1, 3, 0] * pivot;
            lhs[l1, 4, 0] = lhs[l1, 4, 0] * pivot;

            c[c1, c2, c3, c4, 0, 0] = c[c1, c2, c3, c4, 0, 0] * pivot;
            c[c1, c2, c3, c4, 1, 0] = c[c1, c2, c3, c4, 1, 0] * pivot;
            c[c1, c2, c3, c4, 2, 0] = c[c1, c2, c3, c4, 2, 0] * pivot;
            c[c1, c2, c3, c4, 3, 0] = c[c1, c2, c3, c4, 3, 0] * pivot;
            c[c1, c2, c3, c4, 4, 0] = c[c1, c2, c3, c4, 4, 0] * pivot;

            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] * pivot;
            coeff = lhs[l1, 0, 1];
            lhs[l1, 1, 1] = lhs[l1, 1, 1] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 1] = lhs[l1, 2, 1] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 1] = lhs[l1, 3, 1] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 0];

            c[c1, c2, c3, c4, 0, 1] = c[c1, c2, c3, c4, 0, 1] - coeff * c[c1, c2, c3, c4, 0, 0];
            c[c1, c2, c3, c4, 1, 1] = c[c1, c2, c3, c4, 1, 1] - coeff * c[c1, c2, c3, c4, 1, 0];
            c[c1, c2, c3, c4, 2, 1] = c[c1, c2, c3, c4, 2, 1] - coeff * c[c1, c2, c3, c4, 2, 0];
            c[c1, c2, c3, c4, 3, 1] = c[c1, c2, c3, c4, 3, 1] - coeff * c[c1, c2, c3, c4, 3, 0];
            c[c1, c2, c3, c4, 4, 1] = c[c1, c2, c3, c4, 4, 1] - coeff * c[c1, c2, c3, c4, 4, 0];

            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 0];
            coeff = lhs[l1, 0, 2];
            lhs[l1, 1, 2] = lhs[l1, 1, 2] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 2] = lhs[l1, 2, 2] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 0];

            c[c1, c2, c3, c4, 0, 2] = c[c1, c2, c3, c4, 0, 2] - coeff * c[c1, c2, c3, c4, 0, 0];
            c[c1, c2, c3, c4, 1, 2] = c[c1, c2, c3, c4, 1, 2] - coeff * c[c1, c2, c3, c4, 1, 0];
            c[c1, c2, c3, c4, 2, 2] = c[c1, c2, c3, c4, 2, 2] - coeff * c[c1, c2, c3, c4, 2, 0];
            c[c1, c2, c3, c4, 3, 2] = c[c1, c2, c3, c4, 3, 2] - coeff * c[c1, c2, c3, c4, 3, 0];
            c[c1, c2, c3, c4, 4, 2] = c[c1, c2, c3, c4, 4, 2] - coeff * c[c1, c2, c3, c4, 4, 0];

            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 0];
            coeff = lhs[l1, 0, 3];
            lhs[l1, 1, 3] = lhs[l1, 1, 3] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 3] = lhs[l1, 2, 3] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 0];

            c[c1, c2, c3, c4, 0, 3] = c[c1, c2, c3, c4, 0, 3] - coeff * c[c1, c2, c3, c4, 0, 0];
            c[c1, c2, c3, c4, 1, 3] = c[c1, c2, c3, c4, 1, 3] - coeff * c[c1, c2, c3, c4, 1, 0];
            c[c1, c2, c3, c4, 2, 3] = c[c1, c2, c3, c4, 2, 3] - coeff * c[c1, c2, c3, c4, 2, 0];
            c[c1, c2, c3, c4, 3, 3] = c[c1, c2, c3, c4, 3, 3] - coeff * c[c1, c2, c3, c4, 3, 0];
            c[c1, c2, c3, c4, 4, 3] = c[c1, c2, c3, c4, 4, 3] - coeff * c[c1, c2, c3, c4, 4, 0];

            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 0];
            coeff = lhs[l1, 0, 4];
            lhs[l1, 1, 4] = lhs[l1, 1, 4] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 4] = lhs[l1, 2, 4] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 0];

            c[c1, c2, c3, c4, 0, 4] = c[c1, c2, c3, c4, 0, 4] - coeff * c[c1, c2, c3, c4, 0, 0];
            c[c1, c2, c3, c4, 1, 4] = c[c1, c2, c3, c4, 1, 4] - coeff * c[c1, c2, c3, c4, 1, 0];
            c[c1, c2, c3, c4, 2, 4] = c[c1, c2, c3, c4, 2, 4] - coeff * c[c1, c2, c3, c4, 2, 0];
            c[c1, c2, c3, c4, 3, 4] = c[c1, c2, c3, c4, 3, 4] - coeff * c[c1, c2, c3, c4, 3, 0];
            c[c1, c2, c3, c4, 4, 4] = c[c1, c2, c3, c4, 4, 4] - coeff * c[c1, c2, c3, c4, 4, 0];

            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 0];
            pivot = 1.00d / lhs[l1, 1, 1];
            lhs[l1, 2, 1] = lhs[l1, 2, 1] * pivot;
            lhs[l1, 3, 1] = lhs[l1, 3, 1] * pivot;
            lhs[l1, 4, 1] = lhs[l1, 4, 1] * pivot;

            c[c1, c2, c3, c4, 0, 1] = c[c1, c2, c3, c4, 0, 1] * pivot;
            c[c1, c2, c3, c4, 1, 1] = c[c1, c2, c3, c4, 1, 1] * pivot;
            c[c1, c2, c3, c4, 2, 1] = c[c1, c2, c3, c4, 2, 1] * pivot;
            c[c1, c2, c3, c4, 3, 1] = c[c1, c2, c3, c4, 3, 1] * pivot;
            c[c1, c2, c3, c4, 4, 1] = c[c1, c2, c3, c4, 4, 1] * pivot;

            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] * pivot;
            coeff = lhs[l1, 1, 0];
            lhs[l1, 2, 0] = lhs[l1, 2, 0] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 0] = lhs[l1, 3, 0] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 0] = lhs[l1, 4, 0] - coeff * lhs[l1, 4, 1];

            c[c1, c2, c3, c4, 0, 0] = c[c1, c2, c3, c4, 0, 0] - coeff * c[c1, c2, c3, c4, 0, 1];
            c[c1, c2, c3, c4, 1, 0] = c[c1, c2, c3, c4, 1, 0] - coeff * c[c1, c2, c3, c4, 1, 1];
            c[c1, c2, c3, c4, 2, 0] = c[c1, c2, c3, c4, 2, 0] - coeff * c[c1, c2, c3, c4, 2, 1];
            c[c1, c2, c3, c4, 3, 0] = c[c1, c2, c3, c4, 3, 0] - coeff * c[c1, c2, c3, c4, 3, 1];
            c[c1, c2, c3, c4, 4, 0] = c[c1, c2, c3, c4, 4, 0] - coeff * c[c1, c2, c3, c4, 4, 1];

            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 1];
            coeff = lhs[l1, 1, 2];
            lhs[l1, 2, 2] = lhs[l1, 2, 2] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 1];

            c[c1, c2, c3, c4, 0, 2] = c[c1, c2, c3, c4, 0, 2] - coeff * c[c1, c2, c3, c4, 0, 1];
            c[c1, c2, c3, c4, 1, 2] = c[c1, c2, c3, c4, 1, 2] - coeff * c[c1, c2, c3, c4, 1, 1];
            c[c1, c2, c3, c4, 2, 2] = c[c1, c2, c3, c4, 2, 2] - coeff * c[c1, c2, c3, c4, 2, 1];
            c[c1, c2, c3, c4, 3, 2] = c[c1, c2, c3, c4, 3, 2] - coeff * c[c1, c2, c3, c4, 3, 1];
            c[c1, c2, c3, c4, 4, 2] = c[c1, c2, c3, c4, 4, 2] - coeff * c[c1, c2, c3, c4, 4, 1];

            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 1];
            coeff = lhs[l1, 1, 3];
            lhs[l1, 2, 3] = lhs[l1, 2, 3] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 1];

            c[c1, c2, c3, c4, 0, 3] = c[c1, c2, c3, c4, 0, 3] - coeff * c[c1, c2, c3, c4, 0, 1];
            c[c1, c2, c3, c4, 1, 3] = c[c1, c2, c3, c4, 1, 3] - coeff * c[c1, c2, c3, c4, 1, 1];
            c[c1, c2, c3, c4, 2, 3] = c[c1, c2, c3, c4, 2, 3] - coeff * c[c1, c2, c3, c4, 2, 1];
            c[c1, c2, c3, c4, 3, 3] = c[c1, c2, c3, c4, 3, 3] - coeff * c[c1, c2, c3, c4, 3, 1];
            c[c1, c2, c3, c4, 4, 3] = c[c1, c2, c3, c4, 4, 3] - coeff * c[c1, c2, c3, c4, 4, 1];

            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 1];
            coeff = lhs[l1, 1, 4];
            lhs[l1, 2, 4] = lhs[l1, 2, 4] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 1];

            c[c1, c2, c3, c4, 0, 4] = c[c1, c2, c3, c4, 0, 4] - coeff * c[c1, c2, c3, c4, 0, 1];
            c[c1, c2, c3, c4, 1, 4] = c[c1, c2, c3, c4, 1, 4] - coeff * c[c1, c2, c3, c4, 1, 1];
            c[c1, c2, c3, c4, 2, 4] = c[c1, c2, c3, c4, 2, 4] - coeff * c[c1, c2, c3, c4, 2, 1];
            c[c1, c2, c3, c4, 3, 4] = c[c1, c2, c3, c4, 3, 4] - coeff * c[c1, c2, c3, c4, 3, 1];
            c[c1, c2, c3, c4, 4, 4] = c[c1, c2, c3, c4, 4, 4] - coeff * c[c1, c2, c3, c4, 4, 1];

            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 1];
            pivot = 1.00d / lhs[l1, 2, 2];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] * pivot;
            lhs[l1, 4, 2] = lhs[l1, 4, 2] * pivot;

            c[c1, c2, c3, c4, 0, 2] = c[c1, c2, c3, c4, 0, 2] * pivot;
            c[c1, c2, c3, c4, 1, 2] = c[c1, c2, c3, c4, 1, 2] * pivot;
            c[c1, c2, c3, c4, 2, 2] = c[c1, c2, c3, c4, 2, 2] * pivot;
            c[c1, c2, c3, c4, 3, 2] = c[c1, c2, c3, c4, 3, 2] * pivot;
            c[c1, c2, c3, c4, 4, 2] = c[c1, c2, c3, c4, 4, 2] * pivot;

            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] * pivot;
            coeff = lhs[l1, 2, 0];
            lhs[l1, 3, 0] = lhs[l1, 3, 0] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 0] = lhs[l1, 4, 0] - coeff * lhs[l1, 4, 2];

            c[c1, c2, c3, c4, 0, 0] = c[c1, c2, c3, c4, 0, 0] - coeff * c[c1, c2, c3, c4, 0, 2];
            c[c1, c2, c3, c4, 1, 0] = c[c1, c2, c3, c4, 1, 0] - coeff * c[c1, c2, c3, c4, 1, 2];
            c[c1, c2, c3, c4, 2, 0] = c[c1, c2, c3, c4, 2, 0] - coeff * c[c1, c2, c3, c4, 2, 2];
            c[c1, c2, c3, c4, 3, 0] = c[c1, c2, c3, c4, 3, 0] - coeff * c[c1, c2, c3, c4, 3, 2];
            c[c1, c2, c3, c4, 4, 0] = c[c1, c2, c3, c4, 4, 0] - coeff * c[c1, c2, c3, c4, 4, 2];

            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 2];
            coeff = lhs[l1, 2, 1];
            lhs[l1, 3, 1] = lhs[l1, 3, 1] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 2];

            c[c1, c2, c3, c4, 0, 1] = c[c1, c2, c3, c4, 0, 1] - coeff * c[c1, c2, c3, c4, 0, 2];
            c[c1, c2, c3, c4, 1, 1] = c[c1, c2, c3, c4, 1, 1] - coeff * c[c1, c2, c3, c4, 1, 2];
            c[c1, c2, c3, c4, 2, 1] = c[c1, c2, c3, c4, 2, 1] - coeff * c[c1, c2, c3, c4, 2, 2];
            c[c1, c2, c3, c4, 3, 1] = c[c1, c2, c3, c4, 3, 1] - coeff * c[c1, c2, c3, c4, 3, 2];
            c[c1, c2, c3, c4, 4, 1] = c[c1, c2, c3, c4, 4, 1] - coeff * c[c1, c2, c3, c4, 4, 2];

            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 2];
            coeff = lhs[l1, 2, 3];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 2];

            c[c1, c2, c3, c4, 0, 3] = c[c1, c2, c3, c4, 0, 3] - coeff * c[c1, c2, c3, c4, 0, 2];
            c[c1, c2, c3, c4, 1, 3] = c[c1, c2, c3, c4, 1, 3] - coeff * c[c1, c2, c3, c4, 1, 2];
            c[c1, c2, c3, c4, 2, 3] = c[c1, c2, c3, c4, 2, 3] - coeff * c[c1, c2, c3, c4, 2, 2];
            c[c1, c2, c3, c4, 3, 3] = c[c1, c2, c3, c4, 3, 3] - coeff * c[c1, c2, c3, c4, 3, 2];
            c[c1, c2, c3, c4, 4, 3] = c[c1, c2, c3, c4, 4, 3] - coeff * c[c1, c2, c3, c4, 4, 2];

            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 2];
            coeff = lhs[l1, 2, 4];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 2];

            c[c1, c2, c3, c4, 0, 4] = c[c1, c2, c3, c4, 0, 4] - coeff * c[c1, c2, c3, c4, 0, 2];
            c[c1, c2, c3, c4, 1, 4] = c[c1, c2, c3, c4, 1, 4] - coeff * c[c1, c2, c3, c4, 1, 2];
            c[c1, c2, c3, c4, 2, 4] = c[c1, c2, c3, c4, 2, 4] - coeff * c[c1, c2, c3, c4, 2, 2];
            c[c1, c2, c3, c4, 3, 4] = c[c1, c2, c3, c4, 3, 4] - coeff * c[c1, c2, c3, c4, 3, 2];
            c[c1, c2, c3, c4, 4, 4] = c[c1, c2, c3, c4, 4, 4] - coeff * c[c1, c2, c3, c4, 4, 2];

            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 2];
            pivot = 1.00d / lhs[l1, 3, 3];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] * pivot;

            c[c1, c2, c3, c4, 0, 3] = c[c1, c2, c3, c4, 0, 3] * pivot;
            c[c1, c2, c3, c4, 1, 3] = c[c1, c2, c3, c4, 1, 3] * pivot;
            c[c1, c2, c3, c4, 2, 3] = c[c1, c2, c3, c4, 2, 3] * pivot;
            c[c1, c2, c3, c4, 3, 3] = c[c1, c2, c3, c4, 3, 3] * pivot;
            c[c1, c2, c3, c4, 4, 3] = c[c1, c2, c3, c4, 4, 3] * pivot;

            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] * pivot;
            coeff = lhs[l1, 3, 0];
            lhs[l1, 4, 0] = lhs[l1, 4, 0] - coeff * lhs[l1, 4, 3];

            c[c1, c2, c3, c4, 0, 0] = c[c1, c2, c3, c4, 0, 0] - coeff * c[c1, c2, c3, c4, 0, 3];
            c[c1, c2, c3, c4, 1, 0] = c[c1, c2, c3, c4, 1, 0] - coeff * c[c1, c2, c3, c4, 1, 3];
            c[c1, c2, c3, c4, 2, 0] = c[c1, c2, c3, c4, 2, 0] - coeff * c[c1, c2, c3, c4, 2, 3];
            c[c1, c2, c3, c4, 3, 0] = c[c1, c2, c3, c4, 3, 0] - coeff * c[c1, c2, c3, c4, 3, 3];
            c[c1, c2, c3, c4, 4, 0] = c[c1, c2, c3, c4, 4, 0] - coeff * c[c1, c2, c3, c4, 4, 3];

            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 3];
            coeff = lhs[l1, 3, 1];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 3];

            c[c1, c2, c3, c4, 0, 1] = c[c1, c2, c3, c4, 0, 1] - coeff * c[c1, c2, c3, c4, 0, 3];
            c[c1, c2, c3, c4, 1, 1] = c[c1, c2, c3, c4, 1, 1] - coeff * c[c1, c2, c3, c4, 1, 3];
            c[c1, c2, c3, c4, 2, 1] = c[c1, c2, c3, c4, 2, 1] - coeff * c[c1, c2, c3, c4, 2, 3];
            c[c1, c2, c3, c4, 3, 1] = c[c1, c2, c3, c4, 3, 1] - coeff * c[c1, c2, c3, c4, 3, 3];
            c[c1, c2, c3, c4, 4, 1] = c[c1, c2, c3, c4, 4, 1] - coeff * c[c1, c2, c3, c4, 4, 3];

            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 3];
            coeff = lhs[l1, 3, 2];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 3];

            c[c1, c2, c3, c4, 0, 2] = c[c1, c2, c3, c4, 0, 2] - coeff * c[c1, c2, c3, c4, 0, 3];
            c[c1, c2, c3, c4, 1, 2] = c[c1, c2, c3, c4, 1, 2] - coeff * c[c1, c2, c3, c4, 1, 3];
            c[c1, c2, c3, c4, 2, 2] = c[c1, c2, c3, c4, 2, 2] - coeff * c[c1, c2, c3, c4, 2, 3];
            c[c1, c2, c3, c4, 3, 2] = c[c1, c2, c3, c4, 3, 2] - coeff * c[c1, c2, c3, c4, 3, 3];
            c[c1, c2, c3, c4, 4, 2] = c[c1, c2, c3, c4, 4, 2] - coeff * c[c1, c2, c3, c4, 4, 3];

            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 3];
            coeff = lhs[l1, 3, 4];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 3];

            c[c1, c2, c3, c4, 0, 4] = c[c1, c2, c3, c4, 0, 4] - coeff * c[c1, c2, c3, c4, 0, 3];
            c[c1, c2, c3, c4, 1, 4] = c[c1, c2, c3, c4, 1, 4] - coeff * c[c1, c2, c3, c4, 1, 3];
            c[c1, c2, c3, c4, 2, 4] = c[c1, c2, c3, c4, 2, 4] - coeff * c[c1, c2, c3, c4, 2, 3];
            c[c1, c2, c3, c4, 3, 4] = c[c1, c2, c3, c4, 3, 4] - coeff * c[c1, c2, c3, c4, 3, 3];
            c[c1, c2, c3, c4, 4, 4] = c[c1, c2, c3, c4, 4, 4] - coeff * c[c1, c2, c3, c4, 4, 3];

            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 3];
            pivot = 1.00d / lhs[l1, 4, 4];
            c[c1, c2, c3, c4, 0, 4] = c[c1, c2, c3, c4, 0, 4] * pivot;
            c[c1, c2, c3, c4, 1, 4] = c[c1, c2, c3, c4, 1, 4] * pivot;
            c[c1, c2, c3, c4, 2, 4] = c[c1, c2, c3, c4, 2, 4] * pivot;
            c[c1, c2, c3, c4, 3, 4] = c[c1, c2, c3, c4, 3, 4] * pivot;
            c[c1, c2, c3, c4, 4, 4] = c[c1, c2, c3, c4, 4, 4] * pivot;

            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] * pivot;
            coeff = lhs[l1, 4, 0];
            c[c1, c2, c3, c4, 0, 0] = c[c1, c2, c3, c4, 0, 0] - coeff * c[c1, c2, c3, c4, 0, 4];
            c[c1, c2, c3, c4, 1, 0] = c[c1, c2, c3, c4, 1, 0] - coeff * c[c1, c2, c3, c4, 1, 4];
            c[c1, c2, c3, c4, 2, 0] = c[c1, c2, c3, c4, 2, 0] - coeff * c[c1, c2, c3, c4, 2, 4];
            c[c1, c2, c3, c4, 3, 0] = c[c1, c2, c3, c4, 3, 0] - coeff * c[c1, c2, c3, c4, 3, 4];
            c[c1, c2, c3, c4, 4, 0] = c[c1, c2, c3, c4, 4, 0] - coeff * c[c1, c2, c3, c4, 4, 4];

            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 4];
            coeff = lhs[l1, 4, 1];
            c[c1, c2, c3, c4, 0, 1] = c[c1, c2, c3, c4, 0, 1] - coeff * c[c1, c2, c3, c4, 0, 4];
            c[c1, c2, c3, c4, 1, 1] = c[c1, c2, c3, c4, 1, 1] - coeff * c[c1, c2, c3, c4, 1, 4];
            c[c1, c2, c3, c4, 2, 1] = c[c1, c2, c3, c4, 2, 1] - coeff * c[c1, c2, c3, c4, 2, 4];
            c[c1, c2, c3, c4, 3, 1] = c[c1, c2, c3, c4, 3, 1] - coeff * c[c1, c2, c3, c4, 3, 4];
            c[c1, c2, c3, c4, 4, 1] = c[c1, c2, c3, c4, 4, 1] - coeff * c[c1, c2, c3, c4, 4, 4];

            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 4];
            coeff = lhs[l1, 4, 2];
            c[c1, c2, c3, c4, 0, 2] = c[c1, c2, c3, c4, 0, 2] - coeff * c[c1, c2, c3, c4, 0, 4];
            c[c1, c2, c3, c4, 1, 2] = c[c1, c2, c3, c4, 1, 2] - coeff * c[c1, c2, c3, c4, 1, 4];
            c[c1, c2, c3, c4, 2, 2] = c[c1, c2, c3, c4, 2, 2] - coeff * c[c1, c2, c3, c4, 2, 4];
            c[c1, c2, c3, c4, 3, 2] = c[c1, c2, c3, c4, 3, 2] - coeff * c[c1, c2, c3, c4, 3, 4];
            c[c1, c2, c3, c4, 4, 2] = c[c1, c2, c3, c4, 4, 2] - coeff * c[c1, c2, c3, c4, 4, 4];

            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 4];
            coeff = lhs[l1, 4, 3];
            c[c1, c2, c3, c4, 0, 3] = c[c1, c2, c3, c4, 0, 3] - coeff * c[c1, c2, c3, c4, 0, 4];
            c[c1, c2, c3, c4, 1, 3] = c[c1, c2, c3, c4, 1, 3] - coeff * c[c1, c2, c3, c4, 1, 4];
            c[c1, c2, c3, c4, 2, 3] = c[c1, c2, c3, c4, 2, 3] - coeff * c[c1, c2, c3, c4, 2, 4];
            c[c1, c2, c3, c4, 3, 3] = c[c1, c2, c3, c4, 3, 3] - coeff * c[c1, c2, c3, c4, 3, 4];
            c[c1, c2, c3, c4, 4, 3] = c[c1, c2, c3, c4, 4, 3] - coeff * c[c1, c2, c3, c4, 4, 4];

            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 4];
        }//lhs[ c[ r[ c[ r[

        public void matvec_sub(double[, ,] ablock, double[, , , ,] avec, double[, , , ,] bvec, int ab1, int av1, int av2, int av3, int av4, int bv1, int bv2, int bv3, int bv4) {
            //---------------------------------------------------------------------
            //     subtracts bvec=bvec - ablock*avec
            //---------------------------------------------------------------------
            //double ablock,avec,bvec
            //dimension ablock[5,5],avec[5],bvec[5]
            //---------------------------------------------------------------------
            //            rhs[i,ic,jc,kc,ccell] = rhs[i,ic,jc,kc,ccell] 
            //     $           - lhs[i,1,ablock,ia,ja,ka,acell]*
            //---------------------------------------------------------------------
            bvec[bv1, bv2, bv3, bv4, 0] = bvec[bv1, bv2, bv3, bv4, 0] - ablock[ab1, 0, 0] * avec[av1, av2, av3, av4, 0] - ablock[ab1, 1, 0] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 0] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 0] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 0] * avec[av1, av2, av3, av4, 4];
            bvec[bv1, bv2, bv3, bv4, 1] = bvec[bv1, bv2, bv3, bv4, 1] - ablock[ab1, 0, 1] * avec[av1, av2, av3, av4, 0] - ablock[ab1, 1, 1] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 1] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 1] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 1] * avec[av1, av2, av3, av4, 4];
            bvec[bv1, bv2, bv3, bv4, 2] = bvec[bv1, bv2, bv3, bv4, 2] - ablock[ab1, 0, 2] * avec[av1, av2, av3, av4, 0] - ablock[ab1, 1, 2] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 2] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 2] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 2] * avec[av1, av2, av3, av4, 4];
            bvec[bv1, bv2, bv3, bv4, 3] = bvec[bv1, bv2, bv3, bv4, 3] - ablock[ab1, 0, 3] * avec[av1, av2, av3, av4, 0] - ablock[ab1, 1, 3] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 3] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 3] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 3] * avec[av1, av2, av3, av4, 4];
            bvec[bv1, bv2, bv3, bv4, 4] = bvec[bv1, bv2, bv3, bv4, 4] - ablock[ab1, 0, 4] * avec[av1, av2, av3, av4, 0] - ablock[ab1, 1, 4] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 4] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 4] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 4] * avec[av1, av2, av3, av4, 4];
        }

        public void matmul_sub(double[, ,] ablock, double[, , , , ,] bblock, double[, ,] cblock, int a1, int b1, int b2, int b3, int b4, int c1) {
            //---------------------------------------------------------------------
            //     subtracts a[i,j,k] X b[i,j,k] from c[i,j,k]
            //---------------------------------------------------------------------
            //      double ablock, bblock, cblock
            //      dimension ablock[5,5], bblock[5,5], cblock[5,5]
            cblock[c1, 0, 0] = cblock[c1, 0, 0] - ablock[a1, 0, 0] * bblock[b1, b2, b3, b4, 0, 0]
                                     - ablock[a1, 1, 0] * bblock[b1, b2, b3, b4, 0, 1]
                                     - ablock[a1, 2, 0] * bblock[b1, b2, b3, b4, 0, 2]
                                     - ablock[a1, 3, 0] * bblock[b1, b2, b3, b4, 0, 3]
                                     - ablock[a1, 4, 0] * bblock[b1, b2, b3, b4, 0, 4];

            cblock[c1, 0, 1] = cblock[c1, 0, 1] - ablock[a1, 0, 1] * bblock[b1, b2, b3, b4, 0, 0]
                                     - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 0, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 0, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 0, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 0, 4];

            cblock[c1, 0, 2] = cblock[c1, 0, 2] - ablock[a1, 0, 2] * bblock[b1, b2, b3, b4, 0, 0]
                                     - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 0, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 0, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 0, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 0, 4];

            cblock[c1, 0, 3] = cblock[c1, 0, 3] - ablock[a1, 0, 3] * bblock[b1, b2, b3, b4, 0, 0]
                                     - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 0, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 0, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 0, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 0, 4];

            cblock[c1, 0, 4] = cblock[c1, 0, 4] - ablock[a1, 0, 4] * bblock[b1, b2, b3, b4, 0, 0]
                                     - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 0, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 0, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 0, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 0, 4];

            cblock[c1, 1, 0] = cblock[c1, 1, 0] - ablock[a1, 0, 0] * bblock[b1, b2, b3, b4, 1, 0]
                                     - ablock[a1, 1, 0] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 0] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 0] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 0] * bblock[b1, b2, b3, b4, 1, 4];

            cblock[c1, 1, 1] = cblock[c1, 1, 1] - ablock[a1, 0, 1] * bblock[b1, b2, b3, b4, 1, 0]
                                     - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 1, 4];

            cblock[c1, 1, 2] = cblock[c1, 1, 2] - ablock[a1, 0, 2] * bblock[b1, b2, b3, b4, 1, 0]
                                     - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 1, 4];

            cblock[c1, 1, 3] = cblock[c1, 1, 3] - ablock[a1, 0, 3] * bblock[b1, b2, b3, b4, 1, 0]
                                     - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 1, 4];

            cblock[c1, 1, 4] = cblock[c1, 1, 4] - ablock[a1, 0, 4] * bblock[b1, b2, b3, b4, 1, 0]
                                     - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 1, 4];

            cblock[c1, 2, 0] = cblock[c1, 2, 0] - ablock[a1, 0, 0] * bblock[b1, b2, b3, b4, 2, 0]
                                     - ablock[a1, 1, 0] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 0] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 0] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 0] * bblock[b1, b2, b3, b4, 2, 4];

            cblock[c1, 2, 1] = cblock[c1, 2, 1] - ablock[a1, 0, 1] * bblock[b1, b2, b3, b4, 2, 0]
                                     - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 2, 4];

            cblock[c1, 2, 2] = cblock[c1, 2, 2] - ablock[a1, 0, 2] * bblock[b1, b2, b3, b4, 2, 0]
                                     - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 2, 4];

            cblock[c1, 2, 3] = cblock[c1, 2, 3] - ablock[a1, 0, 3] * bblock[b1, b2, b3, b4, 2, 0]
                                     - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 2, 4];

            cblock[c1, 2, 4] = cblock[c1, 2, 4] - ablock[a1, 0, 4] * bblock[b1, b2, b3, b4, 2, 0]
                                     - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 2, 4];

            cblock[c1, 3, 0] = cblock[c1, 3, 0] - ablock[a1, 0, 0] * bblock[b1, b2, b3, b4, 3, 0]
                                     - ablock[a1, 1, 0] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 0] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 0] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 0] * bblock[b1, b2, b3, b4, 3, 4];

            cblock[c1, 3, 1] = cblock[c1, 3, 1] - ablock[a1, 0, 1] * bblock[b1, b2, b3, b4, 3, 0]
                                     - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 3, 4];

            cblock[c1, 3, 2] = cblock[c1, 3, 2] - ablock[a1, 0, 2] * bblock[b1, b2, b3, b4, 3, 0]
                                     - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 3, 4];

            cblock[c1, 3, 3] = cblock[c1, 3, 3] - ablock[a1, 0, 3] * bblock[b1, b2, b3, b4, 3, 0]
                                     - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 3, 4];

            cblock[c1, 3, 4] = cblock[c1, 3, 4] - ablock[a1, 0, 4] * bblock[b1, b2, b3, b4, 3, 0]
                                     - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 3, 4];

            cblock[c1, 4, 0] = cblock[c1, 4, 0] - ablock[a1, 0, 0] * bblock[b1, b2, b3, b4, 4, 0]
                                     - ablock[a1, 1, 0] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 0] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 0] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 0] * bblock[b1, b2, b3, b4, 4, 4];

            cblock[c1, 4, 1] = cblock[c1, 4, 1] - ablock[a1, 0, 1] * bblock[b1, b2, b3, b4, 4, 0]
                                     - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 4, 4];

            cblock[c1, 4, 2] = cblock[c1, 4, 2] - ablock[a1, 0, 2] * bblock[b1, b2, b3, b4, 4, 0]
                                     - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 4, 4];

            cblock[c1, 4, 3] = cblock[c1, 4, 3] - ablock[a1, 0, 3] * bblock[b1, b2, b3, b4, 4, 0]
                                     - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 4, 4];

            cblock[c1, 4, 4] = cblock[c1, 4, 4] - ablock[a1, 0, 4] * bblock[b1, b2, b3, b4, 4, 0]
                                     - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 4, 4];
        }

        public void binvrhs(double[, ,] lhs, double[, , , ,] r, int l1, int r1, int r2, int r3, int r4) {
            double pivot, coeff; // dimension lhs[5,5]; r[5];

            pivot = 1.00d / lhs[l1, 0, 0];
            lhs[l1, 1, 0] = lhs[l1, 1, 0] * pivot;
            lhs[l1, 2, 0] = lhs[l1, 2, 0] * pivot;
            lhs[l1, 3, 0] = lhs[l1, 3, 0] * pivot;
            lhs[l1, 4, 0] = lhs[l1, 4, 0] * pivot;
            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] * pivot;

            coeff =         lhs[l1, 0, 1];
            lhs[l1, 1, 1] = lhs[l1, 1, 1] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 1] = lhs[l1, 2, 1] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 1] = lhs[l1, 3, 1] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 0];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 0];

            coeff =         lhs[l1, 0, 2];
            lhs[l1, 1, 2] = lhs[l1, 1, 2] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 2] = lhs[l1, 2, 2] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 0];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 0];

            coeff =         lhs[l1, 0, 3];
            lhs[l1, 1, 3] = lhs[l1, 1, 3] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 3] = lhs[l1, 2, 3] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 0];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 0];

            coeff =         lhs[l1, 0, 4];
            lhs[l1, 1, 4] = lhs[l1, 1, 4] - coeff * lhs[l1, 1, 0];
            lhs[l1, 2, 4] = lhs[l1, 2, 4] - coeff * lhs[l1, 2, 0];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 0];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 0];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 0];

            pivot = 1.00d / lhs[l1, 1, 1];
            lhs[l1, 2, 1] = lhs[l1, 2, 1] * pivot;
            lhs[l1, 3, 1] = lhs[l1, 3, 1] * pivot;
            lhs[l1, 4, 1] = lhs[l1, 4, 1] * pivot;
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] * pivot;

            coeff =         lhs[l1, 1, 0];
            lhs[l1, 2, 0] = lhs[l1, 2, 0] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 0] = lhs[l1, 3, 0] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 0] = lhs[l1, 4, 0] - coeff * lhs[l1, 4, 1];
            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 1];

            coeff =         lhs[l1, 1, 2];
            lhs[l1, 2, 2] = lhs[l1, 2, 2] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 1];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 1];

            coeff =         lhs[l1, 1, 3];
            lhs[l1, 2, 3] = lhs[l1, 2, 3] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 1];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 1];

            coeff =         lhs[l1, 1, 4];
            lhs[l1, 2, 4] = lhs[l1, 2, 4] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 1];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 1];

            pivot = 1.00d / lhs[l1, 2, 2];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] * pivot;
            lhs[l1, 4, 2] = lhs[l1, 4, 2] * pivot;
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] * pivot;

            coeff =         lhs[l1, 2, 0];
            lhs[l1, 3, 0] = lhs[l1, 3, 0] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 0] = lhs[l1, 4, 0] - coeff * lhs[l1, 4, 2];
            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 2];

            coeff =         lhs[l1, 2, 1];
            lhs[l1, 3, 1] = lhs[l1, 3, 1] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 2];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 2];

            coeff =         lhs[l1, 2, 3];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 2];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 2];

            coeff =         lhs[l1, 2, 4];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 2];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 2];

            pivot = 1.00d / lhs[l1, 3, 3];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] * pivot;
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] * pivot;

            coeff = lhs[l1, 3, 0];
            lhs[l1, 4, 0] = lhs[l1, 4, 0] - coeff * lhs[l1, 4, 3];
            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 3];

            coeff = lhs[l1, 3, 1];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 3];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 3];

            coeff = lhs[l1, 3, 2];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 3];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 3];

            coeff = lhs[l1, 3, 4];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 3];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 3];

            pivot = 1.00d / lhs[l1, 4, 4];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] * pivot;

            coeff = lhs[l1, 4, 0];
            r[r1, r2, r3, r4, 0] = r[r1, r2, r3, r4, 0] - coeff * r[r1, r2, r3, r4, 4];

            coeff = lhs[l1, 4, 1];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 4];

            coeff = lhs[l1, 4, 2];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 4];

            coeff = lhs[l1, 4, 3];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 4];
        }//lhs[ lhs[ r[ r[ c[ c[
    }
}




