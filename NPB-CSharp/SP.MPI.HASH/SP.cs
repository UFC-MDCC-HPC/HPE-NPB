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
!    the NPB3_0_JAV SP code.                                              !
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
using NPB3_0_JAV.SPThreads;
using NPB3_0_JAV.BMInOut;
using MPI;

namespace NPB3_0_JAV
{

    public class SP : SPBase
    {
        public int bid = -1;
        public static int t_total = 1;
        public BMResults results;
        public SP(char clss)
            : base(clss)
        {
        }
        public static void Main(String[] argv)
        {
            SP sp = null;

            BMArgs.ParseCmdLineArgs(argv, BMName);
            char CLSS = BMArgs.CLASS;

            try
            {
                sp = new SP(CLSS);
            }
            catch (OutOfMemoryException e)
            {
                Console.WriteLine(e.Message);
                BMArgs.outOfMemoryMessage();
                System.Environment.Exit(0);
            }
            sp.runBenchMark();

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
           
            Console.WriteLine("Rank " + Communicator.world.Rank
                                   + " (running on " + MPI.Environment.ProcessorName + ")");
              
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
			
			create_buffers();
			
            set_constants(0);
            initialize();
            lhsinit();
            exact_rhs();
            compute_buffer_size(5);
			
			init_copy_faces();
			
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

            if(node==root) Console.WriteLine("STARTING"); Console.Out.Flush();

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
            txinvr();
            x_solve();
            y_solve();
            z_solve();
            add();
        }

        public int getInputPars()
        {
            int niter = 0;
            if (File.Exists("inputsp.data"))
            {
                FileStream f2 = new FileStream("inputsp.data", System.IO.FileMode.Open);
                try
                {
                    StreamReader sss = new StreamReader(f2);
                    Console.WriteLine("Reading from input file inputsp.data");
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
                if(node==root) Console.WriteLine("No input file inputsp.data," + "Using compiled defaults");
                niter = niter_default;
                dt = dt_default;
                grid_points[0] = problem_size;
                grid_points[1] = problem_size;
                grid_points[2] = problem_size;
            }
            if(node==root) Console.WriteLine("Size: " + grid_points[0] + " X " + grid_points[1] + " X " + grid_points[2]);
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
            if(node==root) Console.WriteLine("Iterations: " + niter + " dt: " + dt);

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

        public void ninvr(int c)
        {
            int i, j, k;
            double r1, r2, r3, r4, r5, t1, t2;

            int ksize = cell_size[c, 2] + 2;
            int jsize = cell_size[c, 1] + 2;
            int isize = cell_size[c, 0] + 2;

            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {

                        r1 = rhs[c, k, j, i, 0];
                        r2 = rhs[c, k, j, i, 1];
                        r3 = rhs[c, k, j, i, 2];
                        r4 = rhs[c, k, j, i, 3];
                        r5 = rhs[c, k, j, i, 4];

                        t1 = bt * r3;
                        t2 = 0.5d * (r4 + r5);

                        rhs[c, k, j, i, 0] = -r2;
                        rhs[c, k, j, i, 1] = r1;
                        rhs[c, k, j, i, 2] = bt * (r4 - r5);
                        rhs[c, k, j, i, 3] = -t1 + t2;
                        rhs[c, k, j, i, 4] = t1 + t2;
                    }
                }
            }
        }

        public void pinvr(int c)
        {
            int i, j, k;
            double r1, r2, r3, r4, r5, t1, t2;

            int ksize = cell_size[c, 2] + 2;
            int jsize = cell_size[c, 1] + 2;
            int isize = cell_size[c, 0] + 2;

            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {

                        r1 = rhs[c, k, j, i, 0];
                        r2 = rhs[c, k, j, i, 1];
                        r3 = rhs[c, k, j, i, 2];
                        r4 = rhs[c, k, j, i, 3];
                        r5 = rhs[c, k, j, i, 4];

                        t1 = bt * r1;
                        t2 = 0.5d * (r4 + r5);

                        rhs[c, k, j, i, 0] = bt * (r4 - r5);
                        rhs[c, k, j, i, 1] = -r3;
                        rhs[c, k, j, i, 2] = r2;
                        rhs[c, k, j, i, 3] = -t1 + t2;
                        rhs[c, k, j, i, 4] = t1 + t2;
                    }
                }
            }
        }
		
		private double[][] out_buffer_copy_faces;
		private double[][] in_buffer_copy_faces;
		private Request[] requests_copy_faces;
        private int[] b_size_copy_faces;
		
		public void init_copy_faces()
		{
            out_buffer_copy_faces = new double[6][];
            in_buffer_copy_faces = new double[6][];
			
            requests_copy_faces = new Request[12];
            b_size_copy_faces = new int[6];			
			
            b_size_copy_faces[0] = east_size;
            b_size_copy_faces[1] = west_size;
            b_size_copy_faces[2] = north_size;
            b_size_copy_faces[3] = south_size;
            b_size_copy_faces[4] = top_size;
            b_size_copy_faces[5] = bottom_size;

            for (int i = 0; i < 6; i++)
            {
                out_buffer_copy_faces[i] = new double[b_size_copy_faces[i]];
                in_buffer_copy_faces[i] = new double[b_size_copy_faces[i]];
            }
		}

        public void copy_faces()
        {
            int i, j, k, c, m, p0, p1, p2, p3, p4, p5, ksize, jsize, isize;
            Request[] requests;
            
            double[][] out_buffer = out_buffer_copy_faces; // new double[6][];
            double[][] in_buffer = in_buffer_copy_faces; // new double[6][];

            requests = requests_copy_faces; // new Request[12];

            if (no_nodes == 1)
            {
                compute_rhs();
                return;
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
								//if (node == 0)
								//	Console.WriteLine("forcing = " + forcing[c, k, j, i, m]);
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
							
//							for (int mm=0; mm<5;mm++)
//							{
//							    if (node==0) Console.WriteLine("TRACE XI 1: " + rhs[c,k,j,i,mm]);	
//							}
							
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      add fourth order xi-direction dissipation               
                //---------------------------------------------------------------------

                if (start[c, 0] > 2)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                    		i = 3;
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                              (5.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j, i + 1, m] +
                                      u[c, k, j, i + 2, m]);
								
							 //   if (node==0) Console.WriteLine("TRACE XI 2: " + rhs[c,k,j,i,m]);	
                            }
							
                    		i = 4;
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                              (-4.0d * u[c, k, j, i - 1, m] + 6.0d * u[c, k, j, i, m] -
                                4.0d * u[c, k, j, i + 1, m] + u[c, k, j, i + 2, m]);
								
							   // if (node==0) Console.WriteLine("TRACE XI 3: " + rhs[c,k,j,i,m]);	
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
								
							 //   if (node==0) Console.WriteLine("TRACE XI 4: " + rhs[c,k,j,i,m]);	
                            }
                        }
                    }
                }

                if (end[c, 0] > 0)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                    		i = isize - 3;
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j, i - 2, m] - 4.0d * u[c, k, j, i - 1, m] +
                                        6.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j, i + 1, m]);
								
	//						    if (node==0) Console.WriteLine("TRACE XI 5: " + rhs[c,k,j,i,m]);	
                            }
							
                    		i = isize - 2;
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j, i - 2, m] - 4.0d * u[c, k, j, i - 1, m] +
                                        5.0d * u[c, k, j, i, m]);
								
	//						    if (node==0) Console.WriteLine("TRACE XI 6: " + rhs[c,k,j,i,m]);	
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
							
			//				for (int mm=0; mm<5;mm++)
			//				{
			//				    if (node==0) Console.WriteLine("TRACE ETA 1: " + rhs[c,k,j,i,mm]);	
			//				}
                        }
                    }
                }


                //---------------------------------------------------------------------
                //      add fourth order eta-direction dissipation         
                //---------------------------------------------------------------------

                if (start[c, 1] > 2)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                   		j = 3;
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (5.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j + 1, i, m] +
                                              u[c, k, j + 2, i, m]);
								
							   // if (node==0) Console.WriteLine("TRACE ETA 2: " + rhs[c,k,j,i,m]);	
                            }
                        }
						
                   		j = 4;
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (-4.0d * u[c, k, j - 1, i, m] + 6.0d * u[c, k, j, i, m] -
                                        4.0d * u[c, k, j + 1, i, m] + u[c, k, j + 2, i, m]);
								
							 //   if (node==0) Console.WriteLine("TRACE ETA 3: " + rhs[c,k,j,i,m]);	
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
								
//							    if (node==0) Console.WriteLine("TRACE ETA 4: " + rhs[c,k,j,i,m]);	
                            }
                        }
                    }
                }


                if (end[c, 1] > 0)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                   		j = jsize - 3;
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j - 2, i, m] - 4.0d * u[c, k, j - 1, i, m] +
                                        6.0d * u[c, k, j, i, m] - 4.0d * u[c, k, j + 1, i, m]);
								
			//				    if (node==0) Console.WriteLine("TRACE ETA 5: " + rhs[c,k,j,i,m]);	
                            }
                        }
						
                        j = jsize - 2;
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m <= 4; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k, j - 2, i, m] - 4.0d * u[c, k, j - 1, i, m] +
                                        5.0d * u[c, k, j, i, m]);
								
			//				    if (node==0) Console.WriteLine("TRACE ETA 6: " + rhs[c,k,j,i,m]);	
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
							
	//						for (int mm=0; mm<5;mm++)
	//						{
	//						    if (node==0) Console.WriteLine("TRACE ZETA 1: " + rhs[c,k,j,i,mm]);	
	//						}
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
								
				//			    if (node==0) Console.WriteLine("TRACE ZETA 2: " + rhs[c,k,j,i,m]);	
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
								
							   // if (node==0) Console.WriteLine("TRACE ZETA 3: " + rhs[c,k,j,i,m]);	
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
								
							 //   if (node==0) Console.WriteLine("TRACE ZETA 4: " + rhs[c,k,j,i,m]);	
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
								
	//						    if (node==0) Console.WriteLine("TRACE ZETA 5: " + rhs[c,k,j,i,m]);	
                            }
                        }
                    }
					
                    k = ksize - 2;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (m = 0; m < 5; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - dssp *
                                      (u[c, k - 2, j, i, m] - 4.0d * u[c, k - 1, j, i, m] +
                                        5.0d * u[c, k, j, i, m]);
								
			//				    if (node==0) Console.WriteLine("TRACE ZETA 6: " + rhs[c,k,j,i,m]);	
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
							    	
			//				    if (node==0) Console.WriteLine("TRACE DT 0: " + rhs[c,k,j,i,m]);	
                            }
                        }
                    }
                }
            }
        }

        public void txinvr()
        {
            int c, i, j, k, isize, jsize, ksize;
            double t1, t2, t3, ac, ru1, xvel, yvel, zvel,
                   r1, r2, r3, r4, r5, ac2inv;
            
            for (c = 0; c < ncells; c++)
            {
                ksize = cell_size[c, 2] + 2;
                jsize = cell_size[c, 1] + 2;
                isize = cell_size[c, 0] + 2;

                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            ru1 = rho_i[c, k, j, i];
                            xvel = us[c, k, j, i];
                            yvel = vs[c, k, j, i];
                            zvel = ws[c, k, j, i];
                            ac = speed[c, k, j, i];
                            ac2inv = ainv[c, k, j, i] * ainv[c, k, j, i];

                            r1 = rhs[c, k, j, i, 0];
                            r2 = rhs[c, k, j, i, 1];
                            r3 = rhs[c, k, j, i, 2];
                            r4 = rhs[c, k, j, i, 3];
                            r5 = rhs[c, k, j, i, 4];
                                
                            t1 = c2 * ac2inv * (qs[c, k, j, i] * r1 - xvel * r2 -
                                yvel * r3 - zvel * r4 + r5);
                            t2 = bt * ru1 * (xvel * r1 - r2);
                            t3 = (bt * ru1 * ac) * t1;

                            rhs[c, k, j, i, 0] = r1 - t1;
                            rhs[c, k, j, i, 1] = -ru1 * (zvel * r1 - r4);
                            rhs[c, k, j, i, 2] = ru1 * (yvel * r1 - r3);
                            rhs[c, k, j, i, 3] = -t2 + t3;
                            rhs[c, k, j, i, 4] = t2 + t3;
                        }
                    }
                }
            }
        }

        public void tzetar(int c)
        {   
            int i, j, k;
            int ksize, jsize, isize;
            double t1, t2, t3, ac, xvel, yvel, zvel,
                   r1, r2, r3, r4, r5, btuz, acinv, ac2u, uzik1;

            ksize = cell_size[c, 2] + 2;
            jsize = cell_size[c, 1] + 2;
            isize = cell_size[c, 0] + 2;

            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        xvel = us[c, k, j, i];
                        yvel = vs[c, k, j, i];
                        zvel = ws[c, k, j, i];
                        ac = speed[c, k, j, i];
                        acinv = ainv[c, k, j, i];

                        ac2u = ac * ac;

                        r1 = rhs[c, k, j, i, 0];
                        r2 = rhs[c, k, j, i, 1];
                        r3 = rhs[c, k, j, i, 2];
                        r4 = rhs[c, k, j, i, 3];
                        r5 = rhs[c, k, j, i, 4];

                        uzik1 = u[c, k, j, i, 0];
                        btuz = bt * uzik1;

                        t1 = btuz * acinv * (r4 + r5);
                        t2 = r3 + t1;
                        t3 = btuz * (r4 - r5);

                        rhs[c, k, j, i, 0] = t2;
                        rhs[c, k, j, i, 1] = -uzik1 * r2 + xvel * t2;
                        rhs[c, k, j, i, 2] = uzik1 * r1 + yvel * t2;
                        rhs[c, k, j, i, 3] = zvel * t2 + t3;
                        rhs[c, k, j, i, 4] = uzik1 * (-xvel * r2 + yvel * r1) +
                              qs[c, k, j, i] * t2 + c2iv * ac2u * t1 + zvel * t3;

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
            if (grid_points[0] == 12
                  && grid_points[1] == 12
                  && grid_points[2] == 12
                  && no_time_steps == 100
            )
            {

                clss = 'S';
                dtref = .015d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 2.7470315451339479E-2d;
                xcrref[1] = 1.0360746705285417E-2d;
                xcrref[2] = 1.6235745065095532E-2d;
                xcrref[3] = 1.5840557224455615E-2d;
                xcrref[4] = 3.4849040609362460E-2d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 2.7289258557377227E-5d;
                xceref[1] = 1.0364446640837285E-5d;
                xceref[2] = 1.6154798287166471E-5d;
                xceref[3] = 1.5750704994480102E-5d;
                xceref[4] = 3.4177666183390531E-5d;


                //---------------------------------------------------------------------
                //    reference data for 36X36X36 grids after 400 time steps, with DT = 
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 36) &&
                   (grid_points[1] == 36) &&
                   (grid_points[2] == 36) &&
                   (no_time_steps == 400))
            {

                clss = 'W';
                dtref = .0015;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 0.1893253733584E-2d;
                xcrref[1] = 0.1717075447775E-3d;
                xcrref[2] = 0.2778153350936E-3d;
                xcrref[3] = 0.2887475409984E-3d;
                xcrref[4] = 0.3143611161242E-2d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 0.7542088599534E-4d;
                xceref[1] = 0.6512852253086E-5d;
                xceref[2] = 0.1049092285688E-4d;
                xceref[3] = 0.1128838671535E-4d;
                xceref[4] = 0.1212845639773E-3d;

                //---------------------------------------------------------------------
                //    reference data for 64X64X64 grids after 400 time steps, with DT = 1.5d-03
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 64) &&
                  (grid_points[1] == 64) &&
                  (grid_points[2] == 64) &&
                  (no_time_steps == 400))
            {

                clss = 'A';
                dtref = .0015d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 2.4799822399300195d;
                xcrref[1] = 1.1276337964368832d;
                xcrref[2] = 1.5028977888770491d;
                xcrref[3] = 1.4217816211695179d;
                xcrref[4] = 2.1292113035138280d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 1.0900140297820550E-4d;
                xceref[1] = 3.7343951769282091E-5d;
                xceref[2] = 5.0092785406541633E-5d;
                xceref[3] = 4.7671093939528255E-5d;
                xceref[4] = 1.3621613399213001E-4d;

                //---------------------------------------------------------------------
                //    reference data for 102X102X102 grids after 400 time steps,
                //    with DT = 1.0d-03
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 102) &&
                  (grid_points[1] == 102) &&
                  (grid_points[2] == 102) &&
                  (no_time_steps == 400))
            {

                clss = 'B';
                dtref = .001d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 0.6903293579998E+02d;
                xcrref[1] = 0.3095134488084E+02d;
                xcrref[2] = 0.4103336647017E+02d;
                xcrref[3] = 0.3864769009604E+02d;
                xcrref[4] = 0.5643482272596E+02d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 0.9810006190188E-02d;
                xceref[1] = 0.1022827905670E-02d;
                xceref[2] = 0.1720597911692E-02d;
                xceref[3] = 0.1694479428231E-02d;
                xceref[4] = 0.1847456263981E-01d;

                //---------------------------------------------------------------------
                //    reference data for 162X162X162 grids after 400 time steps,
                //    with DT = 0.67d-03
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 162) &&
                  (grid_points[1] == 162) &&
                  (grid_points[2] == 162) &&
                  (no_time_steps == 400))
            {

                clss = 'C';
                dtref = .00067d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 0.5881691581829E+03d;
                xcrref[1] = 0.2454417603569E+03d;
                xcrref[2] = 0.3293829191851E+03d;
                xcrref[3] = 0.3081924971891E+03d;
                xcrref[4] = 0.4597223799176E+03d;

                //---------------------------------------------------------------------
                //    Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 0.2598120500183d;
                xceref[1] = 0.2590888922315E-01d;
                xceref[2] = 0.5132886416320E-01d;
                xceref[3] = 0.4806073419454E-01d;
                xceref[4] = 0.5483377491301d;
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
		
		protected double[,,][] in_buffer_solver, out_buffer_solver;
		
		void create_buffers() 		
		{
			int c, stage, isize, jsize, ksize, buffer_size;
			
			in_buffer_solver = new double[3,2,ncells][];			
			out_buffer_solver = new double[3,2,ncells][];
			
			for (stage = 0; stage < ncells; stage++)
            {
                c = slice[stage, 0];
				
                jsize = cell_size[c, 1] + 2;
                ksize = cell_size[c, 2] + 2;
				
		buffer_size = (jsize - start[c, 1] - end[c, 1]) * (ksize - start[c, 2] - end[c, 2]);
				
                in_buffer_solver[0,0,stage] = new double[22*buffer_size];
                out_buffer_solver[0,0,stage] = new double[22*buffer_size];				
				
                in_buffer_solver[0,1,stage] = new double[10 * buffer_size];
                out_buffer_solver[0,1,stage] = new double[10 * buffer_size];
				
                c = slice[stage, 1];

                isize = cell_size[c, 0] + 2;
                ksize = cell_size[c, 2] + 2;

                buffer_size = (isize - start[c, 0] - end[c, 0]) * (ksize - start[c, 2] - end[c, 2]);

                in_buffer_solver[1,0,stage] = new double[22*buffer_size];
                out_buffer_solver[1,0,stage] = new double[22 * buffer_size];
				
                in_buffer_solver[1,1,stage] = new double[10 * buffer_size];
                out_buffer_solver[1,1,stage] = new double[10 * buffer_size];
				
                c = slice[stage, 2];

                isize = cell_size[c, 0] + 2;
                jsize = cell_size[c, 1] + 2;

                buffer_size = (isize - start[c, 0] - end[c, 0]) * (jsize - start[c, 1] - end[c, 1]);

		in_buffer_solver[2,0,stage] = new double[22*buffer_size];
                out_buffer_solver[2,0,stage] = new double[22 * buffer_size];
				
                in_buffer_solver[2,1,stage] = new double[10 * buffer_size];
                out_buffer_solver[2,1,stage] = new double[10 * buffer_size];
				
			}
			
		}
		
		
        public void x_solve()
        {
            int c, i, j, k, n, iend, jsize, ksize, i1, i2, m, buffer_size, p, istart, stage;
            double r1, r2, d, e, sm1, sm2, fac1, fac2;
            double[] s = new double[5];
            Request[] requests = new Request[2] { null, null};
            double[] in_buffer_x;
            double[] out_buffer_x;

            //---------------------------------------------------------------------
            //---------------------------------------------------------------------

            //if (timeron) timer.start(t_xsolve);

            for (stage = 0; stage < ncells; stage++)
            {
                c = slice[stage, 0];

                istart = 2;
                iend = 2 + cell_size[c, 0] - 1;

                jsize = cell_size[c, 1] + 2;
                ksize = cell_size[c, 2] + 2;
								
                // buffer_size = (jsize - start[c, 1] - end[c, 1]) * (ksize - start[c, 2] - end[c, 2]);

                 in_buffer_x = in_buffer_solver[0,0,stage]; // new double[22*buffer_size];
                 out_buffer_x = out_buffer_solver[0,0,stage]; // new double[22*buffer_size];
                
                if (stage != 0)
                {
                    //---------------------------------------------------------------------
                    //            if this is not the first processor in this row of cells, 
                    //            receive data from predecessor containing the right hand
                    //            sides and the upper diagonal elements of the previous two rows
                    //---------------------------------------------------------------------

                    //if (requests[0] != null)
                    //    requestList.Remove(requests[0]);

                    requests[0] = comm_solve.ImmediateReceive<double>(predecessor[0], DEFAULT_TAG, in_buffer_x);
              
                    //requestList.Add(requests[0]);

                    //---------------------------------------------------------------------
                    //            communication has already been started. 
                    //            compute the left hand side while waiting for the msg
                    //---------------------------------------------------------------------
                    lhsx(c);

                    //---------------------------------------------------------------------
                    //            wait for pending communication to complete
                    //            This waits on the current receive and on the send
                    //            from the previous stage. They always come in pairs. 
                    //---------------------------------------------------------------------

                    //requestList.WaitAll();
                    requests[1].Wait();
                    requests[0].Wait();

                    //---------------------------------------------------------------------
                    //            unpack the buffer                                 
                    //---------------------------------------------------------------------
                    i = istart;
                    i1 = istart + 1;
                    n = -1;

                    //---------------------------------------------------------------------
                    //            create a running pointer
                    //---------------------------------------------------------------------

                    p = 0;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                                     in_buffer_x[p    ] * lhs[c, k, j, i, n + 1];
                            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                                     in_buffer_x[p + 1] * lhs[c, k, j, i, n + 1];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                      in_buffer_x[p + 2 + m] * lhs[c, k, j, i, n + 1];
                            }
                            d = in_buffer_x[p + 5];
                            e = in_buffer_x[p + 6];
                            for (m = 0; m <= 2; m++)
                            {
                                s[m] = in_buffer_x[p + 7 + m];
                            }
                            r1 = lhs[c, k, j, i, n + 2];
                            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                            lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                            }
                            r2 = lhs[c, k, j, i1, n + 1];
                            lhs[c, k, j, i1, n + 2] = lhs[c, k, j, i1, n + 2] - d * r2;
                            lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] - e * r2;
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] - s[m] * r2;
                            }
                            p = p + 10;
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                            {
                                lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                                         in_buffer_x[p] * lhs[c, k, j, i, n + 1];
                                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                                         in_buffer_x[p + 1] * lhs[c, k, j, i, n + 2];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         in_buffer_x[p + 2] * lhs[c, k, j, i, n + 1];
                                d = in_buffer_x[p + 3];
                                e = in_buffer_x[p + 4];
                                s[m] = in_buffer_x[p + 5];
                                r1 = lhs[c, k, j, i, n + 2];
                                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                                lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                                r2 = lhs[c, k, j, i1, n + 1];
                                lhs[c, k, j, i1, n + 2] = lhs[c, k, j, i1, n + 2] - d * r2;
                                lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] - e * r2;
                                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] - s[m] * r2;
                                p = p + 6;
                            }
                        }
                    }
                }
                else
                {
                    lhsx(c);
                }

                //---------------------------------------------------------------------
                //                          FORWARD ELIMINATION  
                //---------------------------------------------------------------------

                //---------------------------------------------------------------------
                //      perform the Thomas algorithm; first, FORWARD ELIMINATION     
                //---------------------------------------------------------------------
                n = -1;

                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = istart; i <= iend - 2; i++)
                        {
                            i1 = i + 1;
                            i2 = i + 2;
                            fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                            lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                            }
                            lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] -
                                           lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i1, n + 4] = lhs[c, k, j, i1, n + 4] -
                                           lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                                            lhs[c, k, j, i1, n + 2] * rhs[c, k, j, i, m];
                            }
                            lhs[c, k, j, i2, n + 2] = lhs[c, k, j, i2, n + 2] -
                                           lhs[c, k, j, i2, n + 1] * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i2, n + 3] = lhs[c, k, j, i2, n + 3] -
                                           lhs[c, k, j, i2, n + 1] * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i2, m] = rhs[c, k, j, i2, m] -
                                            lhs[c, k, j, i2, n + 1] * rhs[c, k, j, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      The last two rows in this grid block are a bit different, 
                //      since they do not have two more rows available for the
                //      elimination of off-diagonal entries
                //---------------------------------------------------------------------

                i = iend - 1;
                i1 = iend;
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                        lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                        lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                        }
                        lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] -
                                       lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 4];
                        lhs[c, k, j, i1, n + 4] = lhs[c, k, j, i1, n + 4] -
                                       lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                                        lhs[c, k, j, i1, n + 2] * rhs[c, k, j, i, m];
                        }
                        //---------------------------------------------------------------------
                        //            scale the last row immediately 
                        //---------------------------------------------------------------------
                        fac2 = 1.0d / lhs[c, k, j, i1, n + 3];
                        lhs[c, k, j, i1, n + 4] = fac2 * lhs[c, k, j, i1, n + 4];
                        lhs[c, k, j, i1, n + 5] = fac2 * lhs[c, k, j, i1, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k, j, i1, m] = fac2 * rhs[c, k, j, i1, m];
                        }
                    }
                }

                //---------------------------------------------------------------------
                //      do the u+c and the u-c factors                 
                //---------------------------------------------------------------------

                for (m = 3; m <= 4; m++)
                {
                    n = (m - 2) * 5 - 1;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = istart; i <= iend - 2; i++)
                            {
                                i1 = i + 1;
                                i2 = i + 2;
                                fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                                lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                                lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                                rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                                lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] -
                                              lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 4];
                                lhs[c, k, j, i1, n + 4] = lhs[c, k, j, i1, n + 4] -
                                              lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 5];
                                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                                              lhs[c, k, j, i1, n + 2] * rhs[c, k, j, i, m];
                                lhs[c, k, j, i2, n + 2] = lhs[c, k, j, i2, n + 2] -
                                              lhs[c, k, j, i2, n + 1] * lhs[c, k, j, i, n + 4];
                                lhs[c, k, j, i2, n + 3] = lhs[c, k, j, i2, n + 3] -
                                              lhs[c, k, j, i2, n + 1] * lhs[c, k, j, i, n + 5];
                                rhs[c, k, j, i2, m] = rhs[c, k, j, i2, m] -
                                              lhs[c, k, j, i2, n + 1] * rhs[c, k, j, i, m];
                            }
                        }
                    }


//                    Console.WriteLine(node + ": X-SOLVE - PASS 3" + stage);

                    //---------------------------------------------------------------------
                    //         And again the last two rows separately
                    //---------------------------------------------------------------------
                    i = iend - 1;
                    i1 = iend;

                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                            lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                            rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];                        //*
                            lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] -
                                           lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 4];       //*
                            lhs[c, k, j, i1, n + 4] = lhs[c, k, j, i1, n + 4] -
                                           lhs[c, k, j, i1, n + 2] * lhs[c, k, j, i, n + 5];
                            rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                                           lhs[c, k, j, i1, n + 2] * rhs[c, k, j, i, m];           //*
                            //---------------------------------------------------------------------
                            //               Scale the last row immediately
                            //---------------------------------------------------------------------
                            fac2 = 1.0d / lhs[c, k, j, i1, n + 3];
                            lhs[c, k, j, i1, n + 4] = fac2 * lhs[c, k, j, i1, n + 4];
                            lhs[c, k, j, i1, n + 5] = fac2 * lhs[c, k, j, i1, n + 5];
                            rhs[c, k, j, i1, m] = fac2 * rhs[c, k, j, i1, m];                     //*

                        }
                    }
                }

                //---------------------------------------------------------------------
                //         send information to the next processor, except when this
                //         is the last grid block
                //---------------------------------------------------------------------
                if (stage != ncells - 1)
                {

                    //---------------------------------------------------------------------
                    //            create a running pointer for the send buffer  
                    //---------------------------------------------------------------------
                    p = 0;
                    n = -1;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = iend - 1; i <= iend; i++)
                            {
                                out_buffer_x[p] = lhs[c, k, j, i, n + 4];
                                out_buffer_x[p + 1] = lhs[c, k, j, i, n + 5];
                                for (m = 0; m <= 2; m++)
                                {
                                    out_buffer_x[p + 2 + m] = rhs[c, k, j, i, m];
                                }
                                p = p + 5;
                            }
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                            {
                                for (i = iend - 1; i <= iend; i++)
                                {
                                    out_buffer_x[p] = lhs[c, k, j, i, n + 4];
                                    out_buffer_x[p + 1] = lhs[c, k, j, i, n + 5];
                                    out_buffer_x[p + 2] = rhs[c, k, j, i, m];
                                    p = p + 3;
                                }
                            }
                        }
                    }

                    //if (requests[1] != null)
                    //    requestList.Remove(requests[1]);

                    requests[1] = comm_solve.ImmediateSend<double>(out_buffer_x, successor[0], DEFAULT_TAG);                    

                    //requestList.Add(requests[1]);

                }


            } // cells loop

            //---------------------------------------------------------------------
            //                         BACKSUBSTITUTION 
            //---------------------------------------------------------------------

            for (stage = ncells - 1; stage >= 0; stage--)
            {
                c = slice[stage, 0];

                istart = 2;
                iend = 2 + cell_size[c, 0] - 1;

                jsize = cell_size[c, 1] + 2;
                ksize = cell_size[c, 2] + 2;

                // buffer_size = (jsize - start[c, 1] - end[c, 1]) * (ksize - start[c, 2] - end[c, 2]);

                in_buffer_x = in_buffer_solver[0,1,stage]; // new double[10 * buffer_size];
                out_buffer_x = out_buffer_solver[0,1,stage]; // new double[10 * buffer_size];

                if (stage != ncells - 1)
                {
                    //---------------------------------------------------------------------
                    //            if this is not the starting cell in this row of cells, 
                    //            wait for a message to be received, containing the 
                    //            solution of the previous two stations     
                    //---------------------------------------------------------------------

                    //if (requests[0] != null)
                    //    requestList.Remove(requests[0]);

                    requests[0] = comm_solve.ImmediateReceive<double>(successor[0], DEFAULT_TAG, in_buffer_x);

                    //requestList.Add(requests[0]);

                    //---------------------------------------------------------------------
                    //            communication has already been started
                    //            while waiting, do the block-diagonal inversion for the 
                    //            cell that was just finished                
                    //---------------------------------------------------------------------

                    ninvr(slice[stage + 1, 0]);

                    //---------------------------------------------------------------------
                    //            wait for pending communication to complete
                    //---------------------------------------------------------------------

                    //requestList.WaitAll();
                    requests[1].Wait();
                    requests[0].Wait();

                    //             mpi_waitall(2, requests, statuses, error)

                    //---------------------------------------------------------------------
                    //            unpack the buffer for the first three factors         
                    //---------------------------------------------------------------------
                    n = -1;
                    p = 0;
                    i = iend;
                    i1 = i - 1;
                    for (m = 0; m <= 2; m++)
                    {
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                            {
                                sm1 = in_buffer_x[p];
                                sm2 = in_buffer_x[p + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i, n + 4] * sm1 -
                                       lhs[c, k, j, i, n + 5] * sm2;
                                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                                       lhs[c, k, j, i1, n + 4] * rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i1, n + 5] * sm1;
                                p = p + 2;
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            now unpack the buffer for the remaining two factors
                    //---------------------------------------------------------------------
                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                            {
                                sm1 = in_buffer_x[p];
                                sm2 = in_buffer_x[p + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i, n + 4] * sm1 -
                                       lhs[c, k, j, i, n + 5] * sm2;
                                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                                       lhs[c, k, j, i1, n + 4] * rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i1, n + 5] * sm1;
                                p = p + 2;
                            }
                        }
                    }
                }
                else
                {
                    //---------------------------------------------------------------------
                    //            now we know this is the first grid block on the back sweep,
                    //            so we don't need a message to start the substitution. 
                    //---------------------------------------------------------------------
                    i = iend - 1;
                    i1 = iend;
                    n = -1;
                    for (m = 0; m <= 2; m++)
                    {
                        for (k = start[c, 2]; k < ksize - end[c, 2] - 1; k++)
                        {
                            for (j = start[c, 1]; j < jsize - end[c, 1] - 1; j++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                            lhs[c, k, j, i, n + 4] * rhs[c, k, j, i1, m];
                            }
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                             lhs[c, k, j, i, n + 4] * rhs[c, k, j, i1, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         Whether or not this is the last processor, we always have
                //         to complete the back-substitution 
                //---------------------------------------------------------------------

                //---------------------------------------------------------------------
                //         The first three factors
                //---------------------------------------------------------------------
                n = -1;
                for (m = 0; m <= 2; m++)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = iend - 2; i >= istart; i--)
                            {
                                i1 = i + 1;
                                i2 = i + 2;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         lhs[c, k, j, i, n + 4] * rhs[c, k, j, i1, m] -
                                         lhs[c, k, j, i, n + 5] * rhs[c, k, j, i2, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         And the remaining two
                //---------------------------------------------------------------------
                for (m = 3; m <= 4; m++)
                {
                    n = (m - 2) * 5 - 1;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = iend - 2; i >= istart; i--)
                            {
                                i1 = i + 1;
                                i2 = i + 2;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         lhs[c, k, j, i, n + 4] * rhs[c, k, j, i1, m] -
                                         lhs[c, k, j, i, n + 5] * rhs[c, k, j, i2, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         send on information to the previous processor, if needed
                //---------------------------------------------------------------------
                if (stage != 0)
                {
                    i = istart;
                    i1 = istart + 1;
                    p = 0;
                    for (m = 0; m <= 4; m++)
                    {
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                            {
                                out_buffer_x[p] = rhs[c, k, j, i, m];
                                out_buffer_x[p + 1] = rhs[c, k, j, i1, m];
                                p = p + 2;
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            pack and send the buffer
                    //---------------------------------------------------------------------

                    //if (requests[1] != null)
                    //    requestList.Remove(requests[1]);

                    requests[1] = comm_solve.ImmediateSend<double>(out_buffer_x, predecessor[0], DEFAULT_TAG);

                    //requestList.Add(requests[1]);

                }

                //if (timeron) timer.stop(t_xsolve);

                //---------------------------------------------------------------------
                //         If this was the last stage, do the block-diagonal inversion          
                //---------------------------------------------------------------------
                if (stage == 0)
                {
                   // if (timeron) timer.start(t_ninvr);
                    ninvr(c);
                   // if (timeron) timer.stop(t_ninvr);
                }
            }
        }

        public void y_solve()
        {
            int i, j, k, stage, n, isize, jend, ksize, j1, j2, buffer_size, c, m, p, jstart; /* requests(2), statuses(MPI_STATUS_SIZE, 2);*/
            double r1, r2, d, e, sm1, sm2, fac1, fac2;
            double[] s = new double[5];
            Request[] requests = new Request[2] { null, null };
            double[] in_buffer_y;
            double[] out_buffer_y;

            //---------------------------------------------------------------------
            //---------------------------------------------------------------------

            // if (timeron) timer.start(t_ysolve);

            //---------------------------------------------------------------------
            // now do a sweep on a layer-by-layer basis, i.e. sweeping through cells
            // on this node in the direction of increasing i for the forward sweep,
            // and after that reversing the direction for the backsubstitution  
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //                          FORWARD ELIMINATION  
            //---------------------------------------------------------------------
            for (stage = 0; stage < ncells; stage++)
            {
                c = slice[stage, 1];

                jstart = 2;
                jend = 2 + cell_size[c, 1] - 1;

                isize = cell_size[c, 0] + 2;
                ksize = cell_size[c, 2] + 2;

//                buffer_size = (isize - start[c, 0] - end[c, 0]) * (ksize - start[c, 2] - end[c, 2]);

                in_buffer_y = in_buffer_solver[1,0,stage]; // new double[22*buffer_size];
                out_buffer_y = out_buffer_solver[1,0,stage]; // new double[22 * buffer_size];

                if (stage != 0)
                {

                    //---------------------------------------------------------------------
                    //            if this is not the first processor in this row of cells, 
                    //            receive data from predecessor containing the right hand
                    //            sides and the upper diagonal elements of the previous two rows
                    //---------------------------------------------------------------------

                    //if (requests[0] != null)
                    //    requestList.Remove(requests[0]);

                    requests[0] = comm_solve.ImmediateReceive<double>(predecessor[1], DEFAULT_TAG, in_buffer_y);

                    //requestList.Add(requests[0]);

                    //             mpi_irecv(in_buffer, 22*buffer_size, 
                    //     >                      dp_type, predecessor(2), 
                    //     >                      DEFAULT_TAG, comm_solve, 
                    //     >                      requests(1), error)

                    //---------------------------------------------------------------------
                    //            communication has already been started. 
                    //            compute the left hand side while waiting for the msg
                    //---------------------------------------------------------------------
                    lhsy(c);

                    //---------------------------------------------------------------------
                    //            wait for pending communication to complete
                    //            This waits on the current receive and on the send
                    //            from the previous stage. They always come in pairs. 
                    //---------------------------------------------------------------------

                    //requestList.WaitAll();
                    requests[1].Wait();
                    requests[0].Wait();

                    //             mpi_waitall(2, requests, statuses, error)

                    //---------------------------------------------------------------------
                    //            unpack the buffer                                 
                    //---------------------------------------------------------------------
                    j = jstart;
                    j1 = jstart + 1;
                    n = -1;
                    //---------------------------------------------------------------------
                    //            create a running pointer
                    //---------------------------------------------------------------------
                    p = 0;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                                    in_buffer_y[p] * lhs[c, k, j, i, n + 1];
                            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                                    in_buffer_y[p + 1] * lhs[c, k, j, i, n + 1];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                      in_buffer_y[p + 2 + m] * lhs[c, k, j, i, n + 1];
                            }
                            d = in_buffer_y[p + 5]; ;
                            e = in_buffer_y[p + 6];
                            for (m = 0; m <= 2; m++)
                            {
                                s[m] = in_buffer_y[p + 7 + m];
                            }
                            r1 = lhs[c, k, j, i, n + 2];
                            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                            lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                            }
                            r2 = lhs[c, k, j1, i, n+1];
                            lhs[c, k, j1, i, n + 2] = lhs[c, k, j1, i, n + 2] - d * r2;
                            lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] - e * r2;
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] - s[m] * r2;
                            }
                            p = p + 10;
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                                         in_buffer_y[p] * lhs[c, k, j, i, n + 1];
                                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                                         in_buffer_y[p + 1] * lhs[c, k, j, i, n + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         in_buffer_y[p + 2] * lhs[c, k, j, i, n + 1];
                                d = in_buffer_y[p + 3];
                                e = in_buffer_y[p + 4];
                                s[m] = in_buffer_y[p + 5];
                                r1 = lhs[c, k, j, i, n + 2];
                                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                                lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                                r2 = lhs[c, k, j1, i, n + 1];
                                lhs[c, k, j1, i, n + 2] = lhs[c, k, j1, i, n + 2] - d * r2;
                                lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] - e * r2;
                                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] - s[m] * r2;
                                p = p + 6;
                            }
                        }
                    }
                }
                else
                {

                    //---------------------------------------------------------------------
                    //            if this IS the first cell, we still compute the lhs
                    //---------------------------------------------------------------------
                    lhsy(c);
                }


                //---------------------------------------------------------------------
                //         perform the Thomas algorithm; first, FORWARD ELIMINATION     
                //---------------------------------------------------------------------
                n = -1;

                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = jstart; j <= jend - 2; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            j1 = j + 1;
                            j2 = j + 2;
                            fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                            lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                            }
                            lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] -
                                       lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j1, i, n + 4] = lhs[c, k, j1, i, n + 4] -
                                       lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                                        lhs[c, k, j1, i, n + 2] * rhs[c, k, j, i, m];
                            }
                            lhs[c, k, j2, i, n + 2] = lhs[c, k, j2, i, n + 2] -
                                       lhs[c, k, j2, i, n + 1] * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j2, i, n + 3] = lhs[c, k, j2, i, n + 3] -
                                       lhs[c, k, j2, i, n + 1] * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j2, i, m] = rhs[c, k, j2, i, m] -
                                        lhs[c, k, j2, i, n + 1] * rhs[c, k, j, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         The last two rows in this grid block are a bit different, 
                //         since they do not have two more rows available for the
                //         elimination of off-diagonal entries
                //---------------------------------------------------------------------

                j = jend - 1;
                j1 = jend;
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                        lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                        lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                        }
                        lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] -
                                   lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 4];
                        lhs[c, k, j1, i, n + 4] = lhs[c, k, j1, i, n + 4] -
                                   lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                                    lhs[c, k, j1, i, n + 2] * rhs[c, k, j, i, m];
                        }
                        //---------------------------------------------------------------------
                        //               scale the last row immediately (some of this is
                        //               overkill in case this is the last cell)
                        //---------------------------------------------------------------------
                        fac2 = 1.0d / lhs[c, k, j1, i, n + 3];
                        lhs[c, k, j1, i, n + 4] = fac2 * lhs[c, k, j1, i, n + 4];
                        lhs[c, k, j1, i, n + 5] = fac2 * lhs[c, k, j1, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k, j1, i, m] = fac2 * rhs[c, k, j1, i, m];
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         do the u+c and the u-c factors                 
                //---------------------------------------------------------------------
                for (m = 3; m <= 4; m++)
                {
                    n = (m - 2) * 5 - 1;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = jstart; j <= jend - 2; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                j1 = j + 1;
                                j2 = j + 2;
                                fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                                lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                                lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                                rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                                lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] -
                                           lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 4];
                                lhs[c, k, j1, i, n + 4] = lhs[c, k, j1, i, n + 4] -
                                          lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 5];
                                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                                           lhs[c, k, j1, i, n + 2] * rhs[c, k, j, i, m];
                                lhs[c, k, j2, i, n + 2] = lhs[c, k, j2, i, n + 2] -
                                           lhs[c, k, j2, i, n + 1] * lhs[c, k, j, i, n + 4];
                                lhs[c, k, j2, i, n + 3] = lhs[c, k, j2, i, n + 3] -
                                           lhs[c, k, j2, i, n + 1] * lhs[c, k, j, i, n + 5];
                                rhs[c, k, j2, i, m] = rhs[c, k, j2, i, m] -
                                          lhs[c, k, j2, i, n + 1] * rhs[c, k, j, i, m];
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            And again the last two rows separately
                    //---------------------------------------------------------------------
                    j = jend - 1;
                    j1 = jend;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                            lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                            rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                            lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] -
                                       lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j1, i, n + 4] = lhs[c, k, j1, i, n + 4] -
                                       lhs[c, k, j1, i, n + 2] * lhs[c, k, j, i, n + 5];
                            rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                                       lhs[c, k, j1, i, n + 2] * rhs[c, k, j, i, m];
                            //---------------------------------------------------------------------
                            //              Scale the last row immediately (some of this is overkill
                            //               if this is the last cell)
                            //---------------------------------------------------------------------
                            fac2 = 1.0d / lhs[c, k, j1, i, n + 3];
                            lhs[c, k, j1, i, n + 4] = fac2 * lhs[c, k, j1, i, n + 4];
                            lhs[c, k, j1, i, n + 5] = fac2 * lhs[c, k, j1, i, n + 5];
                            rhs[c, k, j1, i, m] = fac2 * rhs[c, k, j1, i, m];

  
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         send information to the next processor, except when this
                //         is the last grid block;
                //---------------------------------------------------------------------

                if (stage != ncells - 1)
                {

                    //---------------------------------------------------------------------
                    //            create a running pointer for the send buffer  
                    //---------------------------------------------------------------------
                    p = 0;
                    n = -1;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (j = jend - 1; j <= jend; j++)
                            {
                                out_buffer_y[p    ] = lhs[c, k, j, i, n + 4];
                                out_buffer_y[p + 1] = lhs[c, k, j, i, n + 5];
                                for (m = 0; m <= 2; m++)
                                {
                                    out_buffer_y[p + 2 + m] = rhs[c, k, j, i, m];
                                }
                                p = p + 5;
                            }
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                for (j = jend - 1; j <= jend; j++)
                                {
                                    out_buffer_y[p] = lhs[c, k, j, i, n + 4];
                                    out_buffer_y[p + 1] = lhs[c, k, j, i, n + 5];
                                    out_buffer_y[p + 2] = rhs[c, k, j, i, m];
                                    p = p + 3;
                                }
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            pack and send the buffer
                    //---------------------------------------------------------------------

                    //if (requests[1] != null)
                    //    requestList.Remove(requests[1]);

                    requests[1] = comm_solve.ImmediateSend<double>(out_buffer_y, successor[1], DEFAULT_TAG);

                    //requestList.Add(requests[1]);

                }
            }

            //---------------------------------------------------------------------
            //      now go in the reverse direction                      
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //                         BACKSUBSTITUTION 
            //---------------------------------------------------------------------
            for (stage = ncells - 1; stage >= 0; stage--)
            {
                c = slice[stage, 1];

                jstart = 2;
                jend = 2 + cell_size[c, 1] - 1;

                isize = cell_size[c, 0] + 2;
                ksize = cell_size[c, 2] + 2;

//                buffer_size = (isize - start[c, 0] - end[c, 0]) * (ksize - start[c, 2] - end[c, 2]);

                in_buffer_y = in_buffer_solver[1,1,stage]; // new double[10 * buffer_size];
                out_buffer_y = out_buffer_solver[1,1,stage]; // new double[10 * buffer_size];

                if (stage != ncells - 1)
                {

                    //---------------------------------------------------------------------
                    //            if this is not the starting cell in this row of cells, 
                    //            wait for a message to be received, containing the 
                    //            solution of the previous two stations     
                    //---------------------------------------------------------------------

                    //if (requests[0] != null)
                    //    requestList.Remove(requests[0]);

                    requests[0] = comm_solve.ImmediateReceive<double>(successor[1], DEFAULT_TAG, in_buffer_y);

                    //requestList.Add(requests[0]);

                    //             mpi_irecv(in_buffer, 10*buffer_size, 
                    //     >                      dp_type, successor(2), 
                    //     >                      DEFAULT_TAG, comm_solve, 
                    //     >                      requests(1), error)


                    //---------------------------------------------------------------------
                    //            communication has already been started
                    //            while waiting, do the block-diagonal inversion for the 
                    //            cell that was just finished                
                    //---------------------------------------------------------------------

                    pinvr(slice[stage + 1, 1]);

                    //---------------------------------------------------------------------
                    //            wait for pending communication to complete
                    //---------------------------------------------------------------------

                    //requestList.WaitAll();
                    requests[1].Wait();
                    requests[0].Wait();

                    //             mpi_waitall(2, requests, statuses, error);

                    //---------------------------------------------------------------------
                    //            unpack the buffer for the first three factors         
                    //---------------------------------------------------------------------
                    n = -1;
                    p = 0;
                    j = jend;
                    j1 = j - 1;
                    for (m = 0; m <= 2; m++)
                    {
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                sm1 = in_buffer_y[p];
                                sm2 = in_buffer_y[p + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i, n + 4] * sm1 -
                                       lhs[c, k, j, i, n + 5] * sm2;
                                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                                       lhs[c, k, j1, i, n + 4] * rhs[c, k, j, i, m] -
                                       lhs[c, k, j1, i, n + 5] * sm1;
                                p = p + 2;
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            now unpack the buffer for the remaining two factors
                    //---------------------------------------------------------------------
                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                sm1 = in_buffer_y[p];
                                sm2 = in_buffer_y[p + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i, n + 4] * sm1 -
                                       lhs[c, k, j, i, n + 5] * sm2;
                                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                                       lhs[c, k, j1, i, n + 4] * rhs[c, k, j, i, m] -
                                       lhs[c, k, j1, i, n + 5] * sm1;
                                p = p + 2;
                            }
                        }
                    }
                }
                else
                {
                    //---------------------------------------------------------------------
                    //            now we know this is the first grid block on the back sweep,
                    //            so we don't need a message to start the substitution. 
                    //---------------------------------------------------------------------

                    j = jend - 1;
                    j1 = jend;
                    n = -1;
                    for (m = 0; m <= 2; m++)
                    {
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                            lhs[c, k, j, i, n + 4] * rhs[c, k, j1, i, m];
                            }
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                            lhs[c, k, j, i, n + 4] * rhs[c, k, j1, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         Whether or not this is the last processor, we always have
                //         to complete the back-substitution 
                //---------------------------------------------------------------------

                //---------------------------------------------------------------------
                //         The first three factors
                //---------------------------------------------------------------------
                n = -1;
                for (m = 0; m <= 2; m++)
                {
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = jend - 2; j >= jstart; j--)   
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                j1 = j + 1;
                                j2 = j + 2;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         lhs[c, k, j, i, n + 4] * rhs[c, k, j1, i, m] -
                                         lhs[c, k, j, i, n + 5] * rhs[c, k, j2, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         And the remaining two
                //---------------------------------------------------------------------
                for (m = 3; m <= 4; m++)
                {
                    n = (m - 2) * 5 - 1;
                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        for (j = jend - 2; j >= jstart; j--)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                j1 = j + 1;
                                j2 = j1 + 1;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         lhs[c, k, j, i, n + 4] * rhs[c, k, j1, i, m] -
                                         lhs[c, k, j, i, n + 5] * rhs[c, k, j2, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         send on information to the previous processor, if needed
                //---------------------------------------------------------------------
                if (stage != 0)
                {
                    j = jstart;
                    j1 = jstart + 1;
                    p = 0;
                    for (m = 0; m <= 4; m++)
                    {
                        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                out_buffer_y[p] = rhs[c, k, j, i, m];
                                out_buffer_y[p + 1] = rhs[c, k, j1, i, m];
                                p = p + 2;
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            pack and send the buffer
                    //---------------------------------------------------------------------

                    //if (requests[1] != null)
                    //    requestList.Remove(requests[1]);

                    requests[1] = comm_solve.ImmediateSend<double>(out_buffer_y, predecessor[1], DEFAULT_TAG);

                    //requestList.Add(requests[1]);
                }

                // if (timeron) timer.stop(t_ysolve);

                //---------------------------------------------------------------------
                //         If this was the last stage, do the block-diagonal inversion          
                //---------------------------------------------------------------------
                if (stage == 0)
                {
                 //   if (timeron) timer.start(t_pinvr);
                    pinvr(c);
                 //   if (timeron) timer.stop(t_pinvr);
                }

            }

        }

        public void z_solve()
        {
            int i, j, k, stage, n, isize, jsize, kend, k1, k2, buffer_size, c, m, p, kstart;
            double r1, r2, d, e, sm1, sm2, fac1, fac2;
            double[] s = new double[5];
            double[] rtmp = new double[5 * (KMAX + 1)];
            Request[] requests = new Request[2] { null, null };
            double[] in_buffer_z;
            double[] out_buffer_z;            

            //---------------------------------------------------------------------
            // now do a sweep on a layer-by-layer basis, i.e. sweeping through cells
            // on this node in the direction of increasing i for the forward sweep,
            // and after that reversing the direction for the backsubstitution  
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //                          FORWARD ELIMINATION  
            //---------------------------------------------------------------------
            for (stage = 0; stage < ncells; stage++)
            {
                c = slice[stage, 2];

                kstart = 2;
                kend = 2 + cell_size[c, 2] - 1;

                isize = cell_size[c, 0] + 2;
                jsize = cell_size[c, 1] + 2;

//                buffer_size = (isize - start[c, 0] - end[c, 0]) * (jsize - start[c, 1] - end[c, 1]);

                in_buffer_z = in_buffer_solver[2,0,stage];// new double[22*buffer_size];
                out_buffer_z = out_buffer_solver[2,0,stage];// new double[22 * buffer_size];

                if (stage != 0)
                {
                    //---------------------------------------------------------------------
                    //            if this is not the first processor in this row of cells, 
                    //            receive data from predecessor containing the right hand
                    //            sides and the upper diagonal elements of the previous two rows
                    //---------------------------------------------------------------------

                    //if (requests[0] != null)
                    //    requestList.Remove(requests[0]);

                    requests[0] = comm_solve.ImmediateReceive<double>(predecessor[2], DEFAULT_TAG, in_buffer_z);

                    //requestList.Add(requests[0]);

                    //---------------------------------------------------------------------
                    //            communication has already been started. 
                    //            compute the left hand side while waiting for the msg
                    //---------------------------------------------------------------------
                    lhsz(c);

                    //---------------------------------------------------------------------
                    //            wait for pending communication to complete
                    //---------------------------------------------------------------------

                    //requestList.WaitAll();
                    requests[1].Wait();
                    requests[0].Wait();

                    //---------------------------------------------------------------------
                    //            unpack the buffer                                 
                    //---------------------------------------------------------------------
                    k = kstart;
                    k1 = kstart + 1;
                    n = -1;

                    //---------------------------------------------------------------------
                    //            create a running pointer
                    //---------------------------------------------------------------------
                    p = 0;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                                     in_buffer_z[p] * lhs[c, k, j, i, n + 1];
                            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                                     in_buffer_z[p + 1] * lhs[c, k, j, i, n + 1];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                      in_buffer_z[p + 2 + m] * lhs[c, k, j, i, n + 1];
                            }
                            d = in_buffer_z[p + 5];
                            e = in_buffer_z[p + 6];
                            for (m = 0; m <= 2; m++)
                            {
                                s[m] = in_buffer_z[p + 7 + m];
                            }
                            r1 = lhs[c, k, j, i, n + 2];
                            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                            lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                            }
                            r2 = lhs[c, k1, j, i, n + 1];
                            lhs[c, k1, j, i, n + 2] = lhs[c, k1, j, i, n + 2] - d * r2;
                            lhs[c, k1, j, i, n + 3] = lhs[c, k1, j, i, n + 3] - e * r2;
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] - s[m] * r2;
                            }
                            p = p + 10;
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                                         in_buffer_z[p] * lhs[c, k, j, i, n + 1];
                                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                                         in_buffer_z[p + 1] * lhs[c, k, j, i, n + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         in_buffer_z[p + 2] * lhs[c, k, j, i, n + 1];
                                d = in_buffer_z[p + 3];
                                e = in_buffer_z[p + 4];
                                s[m] = in_buffer_z[p + 5];
                                r1 = lhs[c, k, j, i, n + 2];
                                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                                lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                                r2 = lhs[c, k1, j, i, n + 1];
                                lhs[c, k1, j, i, n + 2] = lhs[c, k1, j, i, n + 2] - d * r2;
                                lhs[c, k1, j, i, n + 3] = lhs[c, k1, j, i, n + 3] - e * r2;
                                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] - s[m] * r2;
                                p = p + 6;
                            }
                        }
                    }
                }
                else
                {

                    //---------------------------------------------------------------------
                    //            if this IS the first cell, we still compute the lhs
                    //---------------------------------------------------------------------
                    lhsz(c);
                }

                //---------------------------------------------------------------------
                //         perform the Thomas algorithm; first, FORWARD ELIMINATION     
                //---------------------------------------------------------------------
                n = -1;

                for (k = kstart; k <= kend - 2; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            k1 = k + 1;
                            k2 = k + 2;
                            fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                            lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                            }
                            lhs[c, k1, j, i, n + 3] = lhs[c, k1, j, i, n + 3] -
                                       lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 4];
                            lhs[c, k1, j, i, n + 4] = lhs[c, k1, j, i, n + 4] -
                                       lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                                        lhs[c, k1, j, i, n + 2] * rhs[c, k, j, i, m];
                            }
                            lhs[c, k2, j, i, n + 2] = lhs[c, k2, j, i, n + 2] -
                                       lhs[c, k2, j, i, n + 1] * lhs[c, k, j, i, n + 4];
                            lhs[c, k2, j, i, n + 3] = lhs[c, k2, j, i, n + 3] -
                                       lhs[c, k2, j, i, n + 1] * lhs[c, k, j, i, n + 5];
                            for (m = 0; m <= 2; m++)
                            {
                                rhs[c, k2, j, i, m] = rhs[c, k2, j, i, m] -
                                        lhs[c, k2, j, i, n + 1] * rhs[c, k, j, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         The last two rows in this grid block are a bit different, 
                //         since they do not have two more rows available for the
                //         elimination of off-diagonal entries
                //---------------------------------------------------------------------
                k = kend - 1;
                k1 = kend;
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                        lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                        lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                        }
                        lhs[c, k1, j, i, n + 3] = lhs[c, k1, j, i, n + 3] -
                                   lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 4];
                        lhs[c, k1, j, i, n + 4] = lhs[c, k1, j, i, n + 4] -
                                   lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                                    lhs[c, k1, j, i, n + 2] * rhs[c, k, j, i, m];
                        }
                        //---------------------------------------------------------------------
                        //               scale the last row immediately (some of this is
                        //               overkill in case this is the last cell)
                        //---------------------------------------------------------------------
                        fac2 = 1.0d / lhs[c, k1, j, i, n + 3];
                        lhs[c, k1, j, i, n + 4] = fac2 * lhs[c, k1, j, i, n + 4];
                        lhs[c, k1, j, i, n + 5] = fac2 * lhs[c, k1, j, i, n + 5];
                        for (m = 0; m <= 2; m++)
                        {
                            rhs[c, k1, j, i, m] = fac2 * rhs[c, k1, j, i, m];
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         do the u+c and the u-c factors               
                //---------------------------------------------------------------------
                for (m = 3; m <= 4; m++)
                {
                    n = (m - 2) * 5 - 1;
                    for (k = kstart; k <= kend - 2; k++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                k1 = k + 1;
                                k2 = k + 2;
                                fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                                lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                                lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                                rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                                lhs[c, k1, j, i, n + 3] = lhs[c, k1, j, i, n + 3] -
                                           lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 4];
                                lhs[c, k1, j, i, n + 4] = lhs[c, k1, j, i, n + 4] -
                                           lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 5];
                                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                                           lhs[c, k1, j, i, n + 2] * rhs[c, k, j, i, m];
                                lhs[c, k2, j, i, n + 2] = lhs[c, k2, j, i, n + 2] -
                                           lhs[c, k2, j, i, n + 1] * lhs[c, k, j, i, n + 4];
                                lhs[c, k2, j, i, n + 3] = lhs[c, k2, j, i, n + 3] -
                                           lhs[c, k2, j, i, n + 1] * lhs[c, k, j, i, n + 5];
                                rhs[c, k2, j, i, m] = rhs[c, k2, j, i, m] -
                                           lhs[c, k2, j, i, n + 1] * rhs[c, k, j, i, m];
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            And again the last two rows separately
                    //---------------------------------------------------------------------
                    k = kend - 1;
                    k1 = kend;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            fac1 = 1.0d / lhs[c, k, j, i, n + 3];
                            lhs[c, k, j, i, n + 4] = fac1 * lhs[c, k, j, i, n + 4];
                            lhs[c, k, j, i, n + 5] = fac1 * lhs[c, k, j, i, n + 5];
                            rhs[c, k, j, i, m] = fac1 * rhs[c, k, j, i, m];
                            lhs[c, k1, j, i, n + 3] = lhs[c, k1, j, i, n + 3] -
                                       lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 4];
                            lhs[c, k1, j, i, n + 4] = lhs[c, k1, j, i, n + 4] -
                                       lhs[c, k1, j, i, n + 2] * lhs[c, k, j, i, n + 5];
                            rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                                      lhs[c, k1, j, i, n + 2] * rhs[c, k, j, i, m];
                            //---------------------------------------------------------------------
                            //               Scale the last row immediately (some of this is overkill
                            //               if this is the last cell)
                            //---------------------------------------------------------------------
                            fac2 = 1.0d / lhs[c, k1, j, i, n + 3];
                            lhs[c, k1, j, i, n + 4] = fac2 * lhs[c, k1, j, i, n + 4];
                            lhs[c, k1, j , i, n + 5] = fac2 * lhs[c, k1, j , i, n + 5];
                            rhs[c, k1, j, i, m] = fac2 * rhs[c, k1, j, i, m];

                        }
                    }
                }

                //---------------------------------------------------------------------
                //         send information to the next processor, except when this
                //         is the last grid block,
                //---------------------------------------------------------------------

                if (stage != ncells - 1)
                {

                    //---------------------------------------------------------------------
                    //            create a running pointer for the send buffer  
                    //---------------------------------------------------------------------
                    p = 0;
                    n = -1;
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            for (k = kend - 1; k <= kend; k++)
                            {
                                out_buffer_z[p    ] = lhs[c, k, j, i, n + 4];
                                out_buffer_z[p + 1] = lhs[c, k, j, i, n + 5];
                                for (m = 0; m <= 2; m++)
                                {
                                    out_buffer_z[p + 2 + m] = rhs[c, k, j, i, m];
                                }
                                p = p + 5;
                            }
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                for (k = kend - 1; k <= kend; k++)
                                {
                                    out_buffer_z[p] = lhs[c, k, j, i, n + 4];
                                    out_buffer_z[p + 1] = lhs[c, k, j, i, n + 5];
                                    out_buffer_z[p + 2] = rhs[c, k, j, i, m];
                                    p = p + 3;
                                }
                            }
                        }
                    }

                    //if (requests[1] != null)
                    //    requestList.Remove(requests[1]);

                    requests[1] = comm_solve.ImmediateSend<double>(out_buffer_z, successor[2], DEFAULT_TAG);

                    //requestList.Add(requests[1]);

                }
            }

            //---------------------------------------------------------------------
            //      now go in the reverse direction                      
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //                         BACKSUBSTITUTION 
            //---------------------------------------------------------------------
            for (stage = ncells - 1; stage >= 0; stage--)
            {
                c = slice[stage, 2];

                kstart = 2;
                kend = 2 + cell_size[c, 2] - 1;

                isize = cell_size[c, 0] + 2;
                jsize = cell_size[c, 1] + 2;
				
//                buffer_size = (isize - start[c, 0] - end[c, 0]) * (jsize - start[c, 1] - end[c, 1]);

                in_buffer_z = in_buffer_solver[2,1,stage]; // new double[10 * buffer_size];
               out_buffer_z = out_buffer_solver[2,1,stage]; // new double[10 * buffer_size];

                if (stage != ncells - 1)
                {

                    //---------------------------------------------------------------------
                    //            if this is not the starting cell in this row of cells, 
                    //            wait for a message to be received, containing the 
                    //            solution of the previous two stations     
                    //---------------------------------------------------------------------

                    //if (requests[0] != null)
                    //    requestList.Remove(requests[0]);

                    requests[0] = comm_solve.ImmediateReceive<double>(successor[2], DEFAULT_TAG, in_buffer_z);

                    //requestList.Add(requests[0]);

                    //---------------------------------------------------------------------
                    //            communication has already been started
                    //            while waiting, do the  block-diagonal inversion for the 
                    //            cell that was just finished                
                    //---------------------------------------------------------------------

                    tzetar(slice[stage + 1, 2]);

                    //---------------------------------------------------------------------
                    //            wait for pending communication to complete
                    //---------------------------------------------------------------------

                    //requestList.WaitAll();
                    requests[1].Wait();
                    requests[0].Wait();

                    //---------------------------------------------------------------------
                    //            unpack the buffer for the first three factors         
                    //---------------------------------------------------------------------
                    n = -1;
                    p = 0;
                    k = kend;
                    k1 = k - 1;
                    for (m = 0; m <= 2; m++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                sm1 = in_buffer_z[p];
                                sm2 = in_buffer_z[p + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i, n + 4] * sm1 -
                                       lhs[c, k, j, i, n + 5] * sm2;
                                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                                       lhs[c, k1, j, i, n + 4] * rhs[c, k, j, i, m] -
                                       lhs[c, k1, j , i, n + 5] * sm1;
                                p = p + 2;
                            }
                        }
                    }

                    //---------------------------------------------------------------------
                    //            now unpack the buffer for the remaining two factors
                    //---------------------------------------------------------------------
                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                sm1 = in_buffer_z[p];
                                sm2 = in_buffer_z[p + 1];
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                       lhs[c, k, j, i, n + 4] * sm1 -
                                       lhs[c, k, j, i, n + 5] * sm2;
                                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                                       lhs[c, k1, j, i, n + 4] * rhs[c, k, j, i, m] -
                                       lhs[c, k1, j , i, n + 5] * sm1;
                                p = p + 2;
                            }
                        }
                    }
                }
                else
                {

                    //---------------------------------------------------------------------
                    //            now we know this is the first grid block on the back sweep,
                    //            so we don't need a message to start the substitution. 
                    //---------------------------------------------------------------------

                    k = kend - 1;
                    k1 = kend;
                    n = -1;
                    for (m = 0; m <= 2; m++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                            lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m];
                            }
                        }
                    }

                    for (m = 3; m <= 4; m++)
                    {
                        n = (m - 2) * 5 - 1;
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                            lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m];
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         Whether or not this is the last processor, we always have
                //         to complete the back-substitution 
                //---------------------------------------------------------------------

                //---------------------------------------------------------------------
                //         The first three factors
                //---------------------------------------------------------------------
                n = -1;
                for (m = 0; m <= 2; m++)
                {
                    for (k = kend - 2; k >= kstart; k--)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                k1 = k + 1;
                                k2 = k + 2;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m] -
                                         lhs[c, k, j, i, n + 5] * rhs[c, k2, j, i, m];
                                if (node == 0)
                                {
                                    // Console.WriteLine(i + " " + j + " " + k + " z BACK 1 : " + rhs[c, k, j, i, m]);
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         And the remaining two
                //---------------------------------------------------------------------
                for (m = 3; m <= 4; m++)
                {
                    n = (m - 2) * 5 - 1;
                    for (k = kend - 2; k >= kstart; k--)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                k1 = k + 1;
                                k2 = k + 2;
                                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                         lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m] -
                                         lhs[c, k, j, i, n + 5] * rhs[c, k2, j, i, m];
                                if (node == 0)
                                {
                                    // Console.WriteLine(i + " " + j + " " + k + " z BACK 2 : " + rhs[c, k, j, i, m]);
                                }
                            }
                        }
                    }
                }

                //---------------------------------------------------------------------
                //         send on information to the previous processor, if needed
                //---------------------------------------------------------------------
                if (stage != 0)
                {
                    k = kstart;
                    k1 = kstart + 1;
                    p = 0;
                    for (m = 0; m <= 4; m++)
                    {
                        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                        {
                            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                            {
                                out_buffer_z[p] = rhs[c, k, j, i, m];
                                out_buffer_z[p + 1] = rhs[c, k1, j, i, m];
                                p = p + 2;
                            }
                        }
                    }

                    //if (requests[1] != null)
                    //    requestList.Remove(requests[1]);

                    requests[1] = comm_solve.ImmediateSend<double>(out_buffer_z, predecessor[2], DEFAULT_TAG);

                    //requestList.Add(requests[1]);

                }

                // if (timeron) timer.stop(t_zsolve);

                //---------------------------------------------------------------------
                //         If this was the last stage, do the block-diagonal inversion
                //---------------------------------------------------------------------
                if (stage == 0)
                {
                 //   if (timeron) timer.start(t_tzetar);
                    tzetar(c);
                 //   if (timeron) timer.stop(t_tzetar);
                }

            }

        }

        public void lhsx(int c)
        {
            double ru1;
            int i, j, k;

            int ksize = cell_size[c, 2] + 2;
            int jsize = cell_size[c, 1] + 2;
            int isize = cell_size[c, 0] + 2;

            //---------------------------------------------------------------------
            //      treat only cell c             
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //      first fill the lhs for the u-eigenvalue                   
            //---------------------------------------------------------------------
            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0] - 1; i < isize - end[c, 0] + 1; i++)
                    {
                        ru1 = c3c4 * rho_i[c, k, j, i];
                        cv[i] = us[c, k, j, i];
                        rhon[i] = dmax1(dx2 + con43 * ru1, dx5 + c1c5 * ru1, dxmax + ru1, dx1);
                    }

                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0] = 0.0d;
                        lhs[c, k, j, i, 1] = -dttx2 * cv[i - 1] - dttx1 * rhon[i - 1];
                        lhs[c, k, j, i, 2] = 1.0d + c2dttx1 * rhon[i];
                        lhs[c, k, j, i, 3] = dttx2 * cv[i + 1] - dttx1 * rhon[i + 1];
                        lhs[c, k, j, i, 4] = 0.0d;
                    }
                }
            }

            //---------------------------------------------------------------------
            //      add fourth order dissipation                             
            //---------------------------------------------------------------------
            if (start[c, 0] > 2)
            {
                i = 3;
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz5;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;

                        lhs[c, k, j, i + 1, 1] = lhs[c, k, j, i + 1, 1] - comz4;
                        lhs[c, k, j, i + 1, 2] = lhs[c, k, j, i + 1, 2] + comz6;
                        lhs[c, k, j, i + 1, 3] = lhs[c, k, j, i + 1, 3] - comz4;
                        lhs[c, k, j, i + 1, 4] = lhs[c, k, j, i + 1, 4] + comz1;
                    }
                }
            }

            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = 3 * start[c, 0] - 4; i < isize - 3 * end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;
                    }
                }
            }

            if (end[c, 0] > 0)
            {
                i = isize - 3;
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;

                        lhs[c, k, j, i + 1, 0] = lhs[c, k, j, i + 1, 0] + comz1;
                        lhs[c, k, j, i + 1, 1] = lhs[c, k, j, i + 1, 1] - comz4;
                        lhs[c, k, j, i + 1, 2] = lhs[c, k, j, i + 1, 2] + comz5;
                    }
                }
            }

            //---------------------------------------------------------------------
            //      subsequently, fill the other factors (u+c), (u-c) by a4ing to 
            //      the first  
            //---------------------------------------------------------------------
            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0 + 5] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 5] = lhs[c, k, j, i, 1] -
                                         dttx2 * speed[c, k, j, i - 1];
                        lhs[c, k, j, i, 2 + 5] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 5] = lhs[c, k, j, i, 3] +
                                         dttx2 * speed[c, k, j, i + 1];
                        lhs[c, k, j, i, 4 + 5] = lhs[c, k, j, i, 4];
                        lhs[c, k, j, i, 0 + 10] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 10] = lhs[c, k, j, i, 1] +
                                         dttx2 * speed[c, k, j, i - 1];
                        lhs[c, k, j, i, 2 + 10] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 10] = lhs[c, k, j, i, 3] -
                                         dttx2 * speed[c, k, j, i + 1];
                        lhs[c, k, j, i, 4 + 10] = lhs[c, k, j, i, 4];
                    }
                }
            }
        }

        public void lhsy(int c)
        {
            double ru1;
            int i, j, k;

            int ksize = cell_size[c, 2] + 2;
            int jsize = cell_size[c, 1] + 2;
            int isize = cell_size[c, 0] + 2;

            //---------------------------------------------------------------------
            //      treat only cell c
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //      first fill the lhs for the u-eigenvalue         
            //---------------------------------------------------------------------
            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                {

                    for (j = start[c, 1] - 1; j <= jsize - end[c, 1]; j++)
                    {
                        ru1 = c3c4 * rho_i[c, k, j, i];
                        cv[j] = vs[c, k, j, i];
                        rhoq[j] = dmax1(dy3 + con43 * ru1,
                                        dy5 + c1c5 * ru1,
                                        dymax + ru1,
                                        dy1);
                    }

                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        lhs[c, k, j, i, 0] = 0.0d;
                        lhs[c, k, j, i, 1] = -dtty2 * cv[j - 1] - dtty1 * rhoq[j - 1];
                        lhs[c, k, j, i, 2] = 1.0d + c2dtty1 * rhoq[j];
                        lhs[c, k, j, i, 3] = dtty2 * cv[j + 1] - dtty1 * rhoq[j + 1];
                        lhs[c, k, j, i, 4] = 0.0d;
                    }
                }
            }

            //---------------------------------------------------------------------
            //      add fourth order dissipation                             
            //---------------------------------------------------------------------
            if (start[c, 1] > 2)
            {
                j = 3;
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {

                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz5;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;

                        lhs[c, k, j + 1, i, 1] = lhs[c, k, j + 1, i, 1] - comz4;
                        lhs[c, k, j + 1, i, 2] = lhs[c, k, j + 1, i, 2] + comz6;
                        lhs[c, k, j + 1, i, 3] = lhs[c, k, j + 1, i, 3] - comz4;
                        lhs[c, k, j + 1, i, 4] = lhs[c, k, j + 1, i, 4] + comz1;
                    }
                }
            }

            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = 3 * start[c, 1] - 4; j < jsize - 3 * end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;
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
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;

                        lhs[c, k, j + 1, i, 0] = lhs[c, k, j + 1, i, 0] + comz1;
                        lhs[c, k, j + 1, i, 1] = lhs[c, k, j + 1, i, 1] - comz4;
                        lhs[c, k, j + 1, i, 2] = lhs[c, k, j + 1, i, 2] + comz5;
                    }
                }
            }

            //---------------------------------------------------------------------
            //      subsequently, do the other two factors                    
            //---------------------------------------------------------------------
            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0 + 5] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 5] = lhs[c, k, j, i, 1] -
                                         dtty2 * speed[c, k, j - 1, i];
                        lhs[c, k, j, i, 2 + 5] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 5] = lhs[c, k, j, i, 3] +
                                         dtty2 * speed[c, k, j + 1, i];
                        lhs[c, k, j, i, 4 + 5] = lhs[c, k, j, i, 4];
                        lhs[c, k, j, i, 0 + 10] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 10] = lhs[c, k, j, i, 1] +
                                         dtty2 * speed[c, k, j - 1, i];
                        lhs[c, k, j, i, 2 + 10] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 10] = lhs[c, k, j, i, 3] -
                                         dtty2 * speed[c, k, j + 1, i];
                        lhs[c, k, j, i, 4 + 10] = lhs[c, k, j, i, 4];

                    }
                }
            }
        }

        public void lhsz(int c)
        {
            double ru1;
            int i, j, k;

            int ksize = cell_size[c, 2] + 2;
            int jsize = cell_size[c, 1] + 2;
            int isize = cell_size[c, 0] + 2;

            //---------------------------------------------------------------------
            //      treat only cell c                                         
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            // first fill the lhs for the u-eigenvalue                          
            //---------------------------------------------------------------------
            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
            {
                for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                {

                    for (k = start[c, 2] - 1; k <= ksize - end[c, 2]; k++)
                    {
                        ru1 = c3c4 * rho_i[c, k, j, i];
                        cv[k] = ws[c, k, j, i];
                        rhos[k] = dmax1(dz4 + con43 * ru1,
                                       dz5 + c1c5 * ru1,
                                       dzmax + ru1,
                                      dz1);
                    }

                    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                    {
                        lhs[c, k, j, i, 0] = 0.0d;
                        lhs[c, k, j, i, 1] = -dttz2 * cv[k - 1] - dttz1 * rhos[k - 1];
                        lhs[c, k, j, i, 2] = 1.0d + c2dttz1 * rhos[k];
                        lhs[c, k, j, i, 3] = dttz2 * cv[k + 1] - dttz1 * rhos[k + 1];
                        lhs[c, k, j, i, 4] = 0.0d;
                    }
                }
            }

            //---------------------------------------------------------------------
            //      add fourth order dissipation                                  
            //---------------------------------------------------------------------
            if (start[c, 2] > 2)
            {
                k = 3;
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz5;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;

                        lhs[c, k + 1, j, i, 1] = lhs[c, k + 1, j, i, 1] - comz4;
                        lhs[c, k + 1, j, i, 2] = lhs[c, k + 1, j, i, 2] + comz6;
                        lhs[c, k + 1, j, i, 3] = lhs[c, k + 1, j, i, 3] - comz4;
                        lhs[c, k + 1, j, i, 4] = lhs[c, k + 1, j, i, 4] + comz1;
                    }
                }
            }

            for (k = 3 * start[c, 2] - 4; k < ksize - 3 * end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;
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
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;

                        lhs[c, k + 1, j, i, 0] = lhs[c, k + 1, j, i, 0] + comz1;
                        lhs[c, k + 1, j, i, 1] = lhs[c, k + 1, j, i, 1] - comz4;
                        lhs[c, k + 1, j, i, 2] = lhs[c, k + 1, j, i, 2] + comz5;
                    }
                }
            }


            //---------------------------------------------------------------------
            //      subsequently, fill the other factors (u+c), (u-c) 
            //---------------------------------------------------------------------
            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0 + 5] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 5] = lhs[c, k, j, i, 1] -
                                         dttz2 * speed[c, k - 1, j, i];
                        lhs[c, k, j, i, 2 + 5] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 5] = lhs[c, k, j, i, 3] +
                                         dttz2 * speed[c, k + 1, j, i];
                        lhs[c, k, j, i, 4 + 5] = lhs[c, k, j, i, 4];
                        lhs[c, k, j, i, 0 + 10] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 10] = lhs[c, k, j, i, 1] +
                                         dttz2 * speed[c, k - 1, j, i];
                        lhs[c, k, j, i, 2 + 10] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 10] = lhs[c, k, j, i, 3] -
                                        dttz2 * speed[c, k + 1, j, i];
                        lhs[c, k, j, i, 4 + 10] = lhs[c, k, j, i, 4];
                    }
                }
            }
        }

        public double getTime() { return timer.readTimer(1); }

        public void finalize() 
        {
            Console.WriteLine("LU: is about to be garbage collected");
            //base.finalize();
        }
    }

}




