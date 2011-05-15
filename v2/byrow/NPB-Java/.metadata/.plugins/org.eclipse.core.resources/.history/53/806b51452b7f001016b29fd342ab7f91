

/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                                  B T                                    !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is a serial/multithreaded version of the BT code.     !
!								          !
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
! Translation to Java and to MultiThreaded Code				  !
!	   M. Frumkin							  !
!	   M. Schultz							  !
!-------------------------------------------------------------------------!
*/

using System;
using System.IO;
using NPB3_0_JAV.BTThreads;
using NPB3_0_JAV.BMInOut;

namespace NPB3_0_JAV
{

    public class BT : BTBase
    {

        public int bid = -1;
        public BMResults results;
        public bool serial = true;
        double[][][] fjac;
        double[][][] njac;
        double[][][][] lhs;

        double tmp1;
        double tmp2;
        double tmp3;

        public BT(char clss, int threads, bool ser) : base(clss, threads)
        {
            //super(clss, threads);
            serial = ser;
            fjac =  instantiate_jagged_array_3(problem_size + 1, 5,5);
            njac = instantiate_jagged_array_3(problem_size + 1, 5,5);
            lhs = instantiate_jagged_array_4(problem_size + 1, 3, 5,5);
        }

        public static void Main(String[] argv)
        {
            BT bt = null;

            BMArgs.ParseCmdLineArgs(argv, BMName);
            char CLSS = BMArgs.CLASS;
            int np = BMArgs.num_threads;
            bool serial = BMArgs.serial;
            try
            {
                bt = new BT(CLSS, np, serial);
            }
            catch (OutOfMemoryException e)
            {
				Console.Error.WriteLine(e.Message);
                BMArgs.outOfMemoryMessage();
                Environment.Exit(0);
            }
            bt.runBenchMark();
        }

        public void run() { runBenchMark(); }

        public void runBenchMark()
        {
            BMArgs.Banner(BMName, CLASS, serial, num_threads);

            int numTimers = t_last + 1;
            String[] t_names = new String[numTimers];
            double[] trecs = new double[numTimers];
            setTimers(t_names);
            int niter = getInputPars();

            set_constants();
            initialize();
            exact_rhs();

            //if (!serial) setupThreads(this);
            //---------------------------------------------------------------------
            //      do one time step to touch all code, and reinitialize
            //---------------------------------------------------------------------
            if (serial) adi_serial();
            //else adi();
            initialize();

            timer.resetAllTimers();
            timer.start(t_total);

            for (int step = 1; step <= niter; step++)
            {   //niter
                if (step % 20 == 0 || step == 1 || step == niter)
                {
                    Console.WriteLine("Time step " + step);
                }
                if (serial) adi_serial();
                //else adi();
            }

            timer.stop(t_total);
            int verified = verify(niter);

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
                          serial,
                          num_threads,
                          bid);
            results.print();
            if (timeron) printTimers(t_names, trecs, time);
            //Console.ReadKey(true);
        }

        public double getMFLOPS(double total_time, int niter)
        {
            double mflops = 0.0;
            if (total_time > 0)
            {
                double n3 = grid_points[0] * grid_points[1] * grid_points[2];
                double navg = (grid_points[0] + grid_points[1] + grid_points[2]) / 3.0;
                mflops = 3478.8 * n3 - 17655.7 * Math.Pow(navg, 2) + 28023.7 * navg;
                mflops *= niter / (total_time * 1000000.0);
            }
            return mflops;
        }

        public void adi_serial()
        {
            compute_rhs();
            x_solve();
            y_solve();
            z_solve();
            add();
        }

        /*
              public void adi(){
                if(timeron)timer.start(t_rhs);
                doRHS();

                if(timeron)timer.start(t_rhsx);
                doRHS();
                if(timeron)timer.stop(t_rhsx);   

                if(timeron)timer.start(t_rhsy);
                doRHS(); 
                if(timeron)timer.stop(t_rhsy);

                if(timeron)timer.start(t_rhsz);
                doRHS();
                if(timeron)timer.stop(t_rhsz);

                doRHS();
                if(timeron)timer.stop(t_rhs);

                if(timeron)timer.start(t_xsolve);
                synchronized(this){
                  for(int m=0;m<num_threads;m++)
                synchronized(xsolver[m]){
                      xsolver[m].done=false;
                      xsolver[m].notify();
                    }
                  for(int m=0;m<num_threads;m++)
                      while(!xsolver[m].done){
                    try{wait();}catch(InterruptedException e){} 
                        notifyAll();
                  }
                }
                if(timeron)timer.stop(t_xsolve);

                if(timeron)timer.start(t_ysolve);
                synchronized(this){
                  for(int m=0;m<num_threads;m++)
                synchronized(ysolver[m]){
                      ysolver[m].done=false;
                      ysolver[m].notify();
                    }
                  for(int m=0;m<num_threads;m++)
                      while(!ysolver[m].done){
                    try{wait();}catch(InterruptedException e){}
                        notifyAll();
                  }
                 }
                 if(timeron)timer.stop(t_ysolve);

                if(timeron)timer.start(t_zsolve);
                synchronized(this){
                  for(int m=0;m<num_threads;m++)
                synchronized(zsolver[m]){
                      zsolver[m].done=false;
                      zsolver[m].notify();
                    }
                  for(int m=0;m<num_threads;m++)
                      while(!zsolver[m].done){
                    try{wait();}catch(InterruptedException e){} 
                        notifyAll();
                  }
                }
                if(timeron)timer.stop(t_zsolve);

                if(timeron)timer.start(t_add);
                synchronized(this){
                  for(int m=0;m<num_threads;m++)
                synchronized(rhsadder[m]){
                      rhsadder[m].done=false;
                      rhsadder[m].notify();
                    }
                  for(int m=0;m<num_threads;m++)
                      while(!rhsadder[m].done){
                    try{wait();}catch(InterruptedException e){} 
                        notifyAll();
                  }
                }
                if(timeron)timer.stop(t_add);    
              } 

              synchronized void doRHS(){
                int m;
                for(m=0;m<num_threads;m++)
                synchronized(rhscomputer[m]){
                      rhscomputer[m].done=false;
                      rhscomputer[m].notify();
                    }
                for(m=0;m<num_threads;m++)
                while(!rhscomputer[m].done){
                  try{wait();}catch(InterruptedException e){}
                      notifyAll();
                }
              }
        */

        public void printTimers(String[] t_names, double[] trecs, double tmax)
        {
            //DecimalFormat fmt = new DecimalFormat("0.000");
            
            double t;
            Console.WriteLine("SECTION  Time           (secs)");
            for (int i = 1; i <= t_last; i++) trecs[i] = timer.readTimer(i);
            if (tmax == 0.0) tmax = 1.0;
            for (int i = 1; i <= t_last; i++)
            {
                Console.WriteLine(t_names[i] + ":" + String.Format("{0:0.000}", trecs[i]) + ":" +
                                   "  (" + String.Format("{0:0.000}", trecs[i] * 100 / tmax) + "%)");

                if (i == t_rhs)
                {
                    t = trecs[t_rhsx] + trecs[t_rhsy] + trecs[t_rhsz];
                    Console.WriteLine("    --> total ");
                    Console.WriteLine("sub-rhs ");
                    Console.WriteLine(String.Format("{0:0.000}", t));
                    Console.WriteLine("  (");
                    Console.WriteLine(String.Format("{0:0.000}",t * 100 / tmax));
                    Console.WriteLine("%)");
                    t = trecs[t_rhs] - trecs[t_rhsx] + trecs[t_rhsy] + trecs[t_rhsz];
                    Console.WriteLine("    --> total ");
                    Console.WriteLine("rest-rhs ");
                    Console.WriteLine(String.Format("{0:0.000}",t));
                    Console.WriteLine("  (");
                    Console.WriteLine(String.Format("{0:0.000}",t * 100 / tmax));
                    Console.WriteLine("%)");
                }
                else if (i == t_zsolve)
                {
                    t = trecs[t_zsolve] - trecs[t_rdis1] - trecs[t_rdis2];
                    Console.WriteLine("    --> total ");
                    Console.WriteLine("sub-zsol ");
                    Console.WriteLine(String.Format("{0:0.000}",t));
                    Console.WriteLine("  ");
                    Console.WriteLine(String.Format("{0:0.000}",t * 100 / tmax));
                    Console.WriteLine();
                }
                else if (i == t_rdis2)
                {
                    t = trecs[t_rdis1] + trecs[t_rdis2];
                    Console.WriteLine("    --> total ");
                    Console.WriteLine("redist ");
                    Console.WriteLine(String.Format("{0:0.000}",t));
                    Console.WriteLine("  ");
                    Console.WriteLine(String.Format("{0:0.000}",t * 100 / tmax));
                }
            }
        }


        public int getInputPars()
        {
            int niter = 0;
            //File f2 = new File("inputbt.data");
            if (File.Exists("inputbt.data"))
            {
            	FileStream f2 = new FileStream("inputbt.data", System.IO.FileMode.Open);
                try
                {
                    //FileInputStream fis = new FileInputStream(f2);
                    StreamReader datafile = new StreamReader(f2);
                    //StreamReader datafile = new StreamReader(fis);
                    Console.WriteLine("Reading from input file inputbt.data");
                    niter = int.Parse(datafile.ReadLine());
                    dt = double.Parse(datafile.ReadLine());
                    grid_points[0] = int.Parse(datafile.ReadLine());
                    grid_points[1] = int.Parse(datafile.ReadLine());
                    grid_points[2] = int.Parse(datafile.ReadLine());
                    datafile.Close();
                }
                catch (Exception e)
                {
                    Console.Error.WriteLine("exception caught! " + e.Message);
                }
            }
            else
            {
                Console.WriteLine("No input file inputbt.data, Using compiled defaults");
                niter = niter_default;
                dt = dt_default;
                grid_points[0] = problem_size;
                grid_points[1] = problem_size;
                grid_points[2] = problem_size;
            }
            Console.WriteLine("Size: " + grid_points[0]
                       + " X " + grid_points[1]
                       + " X " + grid_points[2]);
            if ((grid_points[0] > IMAX) ||
             (grid_points[1] > JMAX) ||
             (grid_points[2] > KMAX))
            {
                Console.WriteLine("Problem size too big for array");
                Environment.Exit(0);
            }
            Console.WriteLine("Iterations: " + niter + " dt: " + dt);
            return niter;
        }

        public void setTimers(String[] t_names)
        {
            //File f1 = new File("timer.flag");
            timeron = false;
            if (File.Exists("timer.flag"))
            {
                timeron = true;
                t_names[t_total] = "total    ";
                t_names[t_rhsx] = "rhsx     ";
                t_names[t_rhsy] = "rhsy     ";
                t_names[t_rhsz] = "rhsz     ";
                t_names[t_rhs] = "rhs      ";
                t_names[t_xsolve] = "xsolve   ";
                t_names[t_ysolve] = "ysolve   ";
                t_names[t_zsolve] = "zsolve   ";
                t_names[t_rdis1] = "redist1  ";
                t_names[t_rdis2] = "redist2  ";
                t_names[t_add] = "add      ";
            }
        }

        public void rhs_norm(double[] rms, int rmsoffst)
        {
            int i, j, k, d, m;
            double add;

            for (m = 0; m < rms.Length; m++) rms[m + rmsoffst] = 0.0;

            for (k = 1; k <= grid_points[2] - 2; k++)
            {
                for (j = 1; j <= grid_points[1] - 2; j++)
                {
                    for (i = 1; i <= grid_points[0] - 2; i++)
                    {
                        for (m = 0; m < rms.Length; m++)
                        {
                            add = rhs[k][j][i][m];
                            rms[m] += add * add;
                        }
                    }
                }
            }

            for (m = 0; m < rms.Length; m++)
            {
                for (d = 0; d <= 2; d++)
                {
                    rms[m] /= grid_points[d] - 2;
                }
                rms[m] = Math.Sqrt(rms[m + rmsoffst]);
            }
        }

        public void error_norm(double[] rms, int rmsoffst)
        {
            int i, j, k, m, d;
            double[] u_exact = new double[5];
            double xi, eta, zeta, add;

            for (m = 0; m < rms.Length; m++) rms[m + rmsoffst] = 0.0;

            for (k = 0; k <= grid_points[2] - 1; k++)
            {
                zeta = k * dnzm1;
                for (j = 0; j <= grid_points[1] - 1; j++)
                {
                    eta = j * dnym1;
                    for (i = 0; i <= grid_points[0] - 1; i++)
                    {
                        xi = i * dnxm1;
                        exact_solution(xi, eta, zeta, u_exact, 0);
                        for (m = 0; m < rms.Length; m++)
                        {
                            add = u[k][j][i][m] - u_exact[m];
                            rms[m] += add * add;
                        }
                    }
                }
            }

            for (m = 0; m < rms.Length; m++)
            {
                for (d = 0; d <= 2; d++)
                {
                    rms[m] /= grid_points[d] - 2;
                }
                rms[m] = Math.Sqrt(rms[m]);
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
            error_norm(xce, 0);
            compute_rhs();
            rhs_norm(xcr, 0);

            for (m = 0; m < xcr.Length; m++) xcr[m] = xcr[m] / dt;

            for (m = 1; m < xcrref.Length; m++)
            {
                xcrref[m] = 1.0;
                xceref[m] = 1.0;
            }

            //---------------------------------------------------------------------
            //    reference data for 12X12X12 grids after 100 time steps, with DT = 1.0d-02
            //---------------------------------------------------------------------
            if (grid_points[0] == 12
                  && grid_points[1] == 12
                  && grid_points[2] == 12
                  && no_time_steps == 60
            )
            {

                clss = 'S';
                dtref = 0.01;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 1.7034283709541311E-01;
                xcrref[1] = 1.2975252070034097E-02;
                xcrref[2] = 3.2527926989486055E-02;
                xcrref[3] = 2.6436421275166801E-02;
                xcrref[4] = 1.9211784131744430E-01;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 4.9976913345811579E-04;
                xceref[1] = 4.5195666782961927E-05;
                xceref[2] = 7.3973765172921357E-05;
                xceref[3] = 7.3821238632439731E-05;
                xceref[4] = 8.9269630987491446E-04;

                //---------------------------------------------------------------------
                //    reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 24) &&
                  (grid_points[1] == 24) &&
                  (grid_points[2] == 24) &&
                  (no_time_steps == 200))
            {

                clss = 'W';
                dtref = 0.0008;
                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 0.1125590409344E+03;
                xcrref[1] = 0.1180007595731E+02;
                xcrref[2] = 0.2710329767846E+02;
                xcrref[3] = 0.2469174937669E+02;
                xcrref[4] = 0.2638427874317E+03;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 0.4419655736008E+01;
                xceref[1] = 0.4638531260002;
                xceref[2] = 0.1011551749967E+01;
                xceref[3] = 0.9235878729944;
                xceref[4] = 0.1018045837718E+02;


                //---------------------------------------------------------------------
                //    reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 64) &&
                 (grid_points[1] == 64) &&
                 (grid_points[2] == 64) &&
                 (no_time_steps == 200))
            {

                clss = 'A';
                dtref = 0.0008;
                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 1.0806346714637264E+02;
                xcrref[1] = 1.1319730901220813E+01;
                xcrref[2] = 2.5974354511582465E+01;
                xcrref[3] = 2.3665622544678910E+01;
                xcrref[4] = 2.5278963211748344E+02;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 4.2348416040525025;
                xceref[1] = 4.4390282496995698E-01;
                xceref[2] = 9.6692480136345650E-01;
                xceref[3] = 8.8302063039765474E-01;
                xceref[4] = 9.7379901770829278;

                //---------------------------------------------------------------------
                //    reference data for 102X102X102 grids after 200 time steps,
                //    with DT = 3.0d-04
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 102) &&
              (grid_points[1] == 102) &&
              (grid_points[2] == 102) &&
              (no_time_steps == 200))
            {

                clss = 'B';
                dtref = .0003;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 1.4233597229287254E+03;
                xcrref[1] = 9.9330522590150238E+01;
                xcrref[2] = 3.5646025644535285E+02;
                xcrref[3] = 3.2485447959084092E+02;
                xcrref[4] = 3.2707541254659363E+03;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 5.2969847140936856E+01;
                xceref[1] = 4.4632896115670668;
                xceref[2] = 1.3122573342210174E+01;
                xceref[3] = 1.2006925323559144E+01;
                xceref[4] = 1.2459576151035986E+02;

                //---------------------------------------------------------------------
                //    reference data for 162X162X162 grids after 200 time steps,
                //    with DT = .0001
                //---------------------------------------------------------------------
            }
            else if ((grid_points[0] == 162) &&
                 (grid_points[1] == 162) &&
                 (grid_points[2] == 162) &&
                 (no_time_steps == 200))
            {

                clss = 'C';
                dtref = .0001;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of residual.
                //---------------------------------------------------------------------
                xcrref[0] = 0.62398116551764615E+04;
                xcrref[1] = 0.50793239190423964E+03;
                xcrref[2] = 0.15423530093013596E+04;
                xcrref[3] = 0.13302387929291190E+04;
                xcrref[4] = 0.11604087428436455E+05;

                //---------------------------------------------------------------------
                //  Reference values of RMS-norms of solution error.
                //---------------------------------------------------------------------
                xceref[0] = 0.16462008369091265E+03;
                xceref[1] = 0.11497107903824313E+02;
                xceref[2] = 0.41207446207461508E+02;
                xceref[3] = 0.37087651059694167E+02;
                xceref[4] = 0.36211053051841265E+03;
            }
            //---------------------------------------------------------------------
            //    Compute the difference of solution values and the known reference values.
            //---------------------------------------------------------------------
            for (m = 0; m < xcr.Length; m++)
            {
                xcrdif[m] = Math.Abs((xcr[m] - xcrref[m]) / xcrref[m]);
                xcedif[m] = Math.Abs((xce[m] - xceref[m]) / xceref[m]);
            }
            //---------------------------------------------------------------------
            //   tolerance level
            //---------------------------------------------------------------------
            double epsilon = 1.0 * Math.Pow(.1, 8);
            //---------------------------------------------------------------------
            //    Output the comparison of computed results to known cases.
            //---------------------------------------------------------------------
            if (clss != 'U')
            {
                Console.WriteLine("Verification being performed for class " + clss);
                Console.WriteLine("accuracy setting for epsilon = " + epsilon);
                if (Math.Abs(dt - dtref) <= epsilon)
                {
                    verified = 1;
                }
                else
                {
                    verified = 0;
                    clss = 'U';
                    Console.WriteLine("DT does not match the reference value of " + dtref);
                }
            }
            else
            {
                Console.WriteLine("Unknown class");
            }

            if (clss != 'U') Console.WriteLine("Comparison of RMS-norms of residual");
            else Console.WriteLine("RMS-norms of residual");
            verified = BMResults.printComparisonStatus(clss, verified, epsilon,
                                                     xcr, xcrref, xcrdif);

            if (clss != 'U')
            {
                Console.WriteLine("Comparison of RMS-norms of solution error");
            }
            else
            {
                Console.WriteLine("RMS-norms of solution error");
            }
            verified = BMResults.printComparisonStatus(clss, verified, epsilon,
                                                     xce, xceref, xcedif);

            BMResults.printVerificationStatus(clss, verified, BMName);
            return verified;
        }

        public void add()
        {
            int i, j, k, m;
            if (timeron) timer.start(t_add);
            for (k = 1; k <= grid_points[2] - 2; k++)
            {
                for (j = 1; j <= grid_points[1] - 2; j++)
                {
                    for (i = 1; i <= grid_points[0] - 2; i++)
                    {
                        for (m = 0; m <= 4; m++)
                        {
                            u[k][j][i][m] += rhs[k][j][i][m];
                        }
                    }
                }
            }
            if (timeron) timer.stop(t_add);
        }

	public void exact_rhs()
	{
		double[] dtemp = new double[5];
        double xi, eta, zeta, dtpp;
		int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;

		//---------------------------------------------------------------------
		//      initialize                                  
		//---------------------------------------------------------------------
		for (k = 0; k <= grid_points[2] - 1; k++)
		{
			for (j = 0; j <= grid_points[1] - 1; j++)
			{
				for (i = 0; i <= grid_points[0] - 1; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						forcing[k][j][i][m] = 0.0;
					}
				}
			}
		}
		//---------------------------------------------------------------------
		//      xi-direction flux differences                      
		//---------------------------------------------------------------------
		for (k = 1; k <= grid_points[2] - 2; k++)
		{
			zeta = k * dnzm1;
			for (j = 1; j <= grid_points[1] - 2; j++)
			{
				eta = j * dnym1;
				for (i = 0; i <= grid_points[0] - 1; i++)
				{
					xi = i * dnxm1;

					exact_solution(xi, eta, zeta, dtemp, 0);
					for (m = 0; m <= 4; m++)
					{
						ue[i][m] = dtemp[m];
					}

					dtpp = 1.0 / dtemp[0];

					for (m = 1; m <= 4; m++)
					{
						buf[i][m] = dtpp * dtemp[m];
					}

					cuf[i] = buf[i][1] * buf[i][1];
					buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] +
										 buf[i][3] * buf[i][3];
					q[i] = 0.5 * (buf[i][1] * ue[i][1] + buf[i][2] * ue[i][2] +
											buf[i][3] * ue[i][3]);

				}

				for (i = 1; i <= grid_points[0] - 2; i++)
				{
					im1 = i - 1;
					ip1 = i + 1;

					forcing[k][j][i][0] = forcing[k][j][i][0] -
									 tx2 * (ue[ip1][1] - ue[im1][1]) +
									 dx1tx1 * (ue[ip1][0] - 2.0 * ue[i][0] + ue[im1][0]);

					forcing[k][j][i][1] = forcing[k][j][i][1] - tx2 * (
									(ue[ip1][1] * buf[ip1][1] + c2 * (ue[ip1][4] - q[ip1])) -
									(ue[im1][1] * buf[im1][1] + c2 * (ue[im1][4] - q[im1]))) +
									 xxcon1 * (buf[ip1][1] - 2.0 * buf[i][1] + buf[im1][1]) +
									 dx2tx1 * (ue[ip1][1] - 2.0 * ue[i][1] + ue[im1][1]);

					forcing[k][j][i][2] = forcing[k][j][i][2] - tx2 * (
									 ue[ip1][2] * buf[ip1][1] - ue[im1][2] * buf[im1][1]) +
									 xxcon2 * (buf[ip1][2] - 2.0 * buf[i][2] + buf[im1][2]) +
									 dx3tx1 * (ue[ip1][2] - 2.0 * ue[i][2] + ue[im1][2]);


					forcing[k][j][i][3] = forcing[k][j][i][3] - tx2 * (
									 ue[ip1][3] * buf[ip1][1] - ue[im1][3] * buf[im1][1]) +
									 xxcon2 * (buf[ip1][3] - 2.0 * buf[i][3] + buf[im1][3]) +
									 dx4tx1 * (ue[ip1][3] - 2.0 * ue[i][3] + ue[im1][3]);

					forcing[k][j][i][4] = forcing[k][j][i][4] - tx2 * (
									 buf[ip1][1] * (c1 * ue[ip1][4] - c2 * q[ip1]) -
									 buf[im1][1] * (c1 * ue[im1][4] - c2 * q[im1])) +
									 0.5 * xxcon3 * (buf[ip1][0] - 2.0 * buf[i][0] +
												   buf[im1][0]) +
									 xxcon4 * (cuf[ip1] - 2.0 * cuf[i] + cuf[im1]) +
									 xxcon5 * (buf[ip1][4] - 2.0 * buf[i][4] + buf[im1][4]) +
									 dx5tx1 * (ue[ip1][4] - 2.0 * ue[i][4] + ue[im1][4]);
                }

				//---------------------------------------------------------------------
				//            Fourth-order dissipation                         
				//---------------------------------------------------------------------
				for (m = 0; m <= 4; m++)
				{
					i = 1;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
										(5.0 * ue[i][m] - 4.0 * ue[i+1][m] + ue[i+2][m]);
					i = 2;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
									   (-4.0 * ue[i-1][m] + 6.0 * ue[i][m] -
										 4.0 * ue[i+1][m] + ue[i+2][m]);
				}

				for (m = 0; m <= 4; m++)
				{
					for (i = 3; i <= grid_points[0] - 4; i++)
					{
						forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
										 (ue[i-2][m] - 4.0 * ue[i-1][m] +
										  6.0 * ue[i][m] - 4.0 * ue[i+1][m] + ue[i+2][m]);
					}
				}

				for (m = 0; m <= 4; m++)
				{
					i = grid_points[0] - 3;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
									   (ue[i-2][m] - 4.0 * ue[i-1][m] +
										6.0 * ue[i][m] - 4.0 * ue[i+1][m]);
					i = grid_points[0] - 2;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
									   (ue[i-2][m] - 4.0 * ue[i-1][m] + 5.0 * ue[i][m]);
				}
			}
		}

		//---------------------------------------------------------------------
		//  eta-direction flux differences             
		//---------------------------------------------------------------------
		for (k = 1; k <= grid_points[2] - 2; k++)
		{
			zeta = k * dnzm1;
			for (i = 1; i <= grid_points[0] - 2; i++)
			{
				xi = i * dnxm1;

				for (j = 0; j <= grid_points[1] - 1; j++)
				{
					eta = j * dnym1;

					exact_solution(xi, eta, zeta, dtemp, 0);
					for (m = 0; m <= 4; m++)
					{
						ue[j][m] = dtemp[m];
					}
					dtpp = 1.0 / dtemp[0];

					for (m = 1; m <= 4; m++)
					{
						buf[j][m] = dtpp * dtemp[m];
					}

					cuf[j] = buf[j][2] * buf[j][2];
					buf[j][0] = cuf[j] + buf[j][1] * buf[j][1] +
							   buf[j][3] * buf[j][3];
					q[j] = 0.5 * (buf[j][1] * ue[j][1] + buf[j][2] * ue[j][2] +
								  buf[j][3] * ue[j][3]);
				}

				for (j = 1; j <= grid_points[1] - 2; j++)
				{
					jm1 = j - 1;
					jp1 = j + 1;

					forcing[k][j][i][0] = forcing[k][j][i][0] -
						  ty2 * (ue[jp1][2] - ue[jm1][2]) +
						  dy1ty1 * (ue[jp1][0] - 2.0 * ue[j][0] + ue[jm1][0]);

					forcing[k][j][i][1] = forcing[k][j][i][1] - ty2 * (
						  ue[jp1][1] * buf[jp1][2] - ue[jm1][1] * buf[jm1][2]) +
						  yycon2 * (buf[jp1][1] - 2.0 * buf[j][1] + buf[jm1][1]) +
						  dy2ty1 * (ue[jp1][1] - 2.0 * ue[j][1] + ue[jm1][1]);

					forcing[k][j][i][2] = forcing[k][j][i][2] - ty2 * (
						  (ue[jp1][2] * buf[jp1][2] + c2 * (ue[jp1][4] - q[jp1])) -
						  (ue[jm1][2] * buf[jm1][2] + c2 * (ue[jm1][4] - q[jm1]))) +
						  yycon1 * (buf[jp1][2] - 2.0 * buf[j][2] + buf[jm1][2]) +
						  dy3ty1 * (ue[jp1][2] - 2.0 * ue[j][2] + ue[jm1][2]);

					forcing[k][j][i][3] = forcing[k][j][i][3] - ty2 * (
						  ue[jp1][3] * buf[jp1][2] - ue[jm1][3] * buf[jm1][2]) +
						  yycon2 * (buf[jp1][3] - 2.0 * buf[j][3] + buf[jm1][3]) +
						  dy4ty1 * (ue[jp1][3] - 2.0 * ue[j][3] + ue[jm1][3]);

					forcing[k][j][i][4] = forcing[k][j][i][4] - ty2 * (
						  buf[jp1][2] * (c1 * ue[jp1][4] - c2 * q[jp1]) -
						  buf[jm1][2] * (c1 * ue[jm1][4] - c2 * q[jm1])) +
						  0.5 * yycon3 * (buf[jp1][0] - 2.0 * buf[j][0] +
										buf[jm1][0]) +
						  yycon4 * (cuf[jp1] - 2.0 * cuf[j] + cuf[jm1]) +
						  yycon5 * (buf[jp1][4] - 2.0 * buf[j][4] + buf[jm1][4]) +
						  dy5ty1 * (ue[jp1][4] - 2.0 * ue[j][4] + ue[jm1][4]);
				}

				//---------------------------------------------------------------------
				//            Fourth-order dissipation                      
				//---------------------------------------------------------------------
				for (m = 0; m <= 4; m++)
				{
					j = 1;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							  (5.0 * ue[j][m] - 4.0 * ue[j+1][m] + ue[j+2][m]);
					j = 2;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							 (-4.0 * ue[j-1][m] + 6.0 * ue[j][m] -
							   4.0 * ue[j+1][m] + ue[j+2][m]);
				}

				for (m = 0; m <= 4; m++)
				{
					for (j = 3; j <= grid_points[1] - 4; j++)
					{
						forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							  (ue[j-2][m] - 4.0 * ue[j-1][m] +
							   6.0 * ue[j][m] - 4.0 * ue[j+1][m] + ue[j+2][m]);
					}
				}

				for (m = 0; m <= 4; m++)
				{
					j = grid_points[1] - 3;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							 (ue[j-2][m] - 4.0 * ue[j-1][m] +
							  6.0 * ue[j][m] - 4.0 * ue[j+1][m]);
					j = grid_points[1] - 2;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							 (ue[j-2][m] - 4.0 * ue[j-1][m] + 5.0 * ue[j][m]);

				}

			}
		}

		//---------------------------------------------------------------------
		//      zeta-direction flux differences                      
		//---------------------------------------------------------------------
		for (j = 1; j <= grid_points[1] - 2; j++)
		{
			eta = j * dnym1;
			for (i = 1; i <= grid_points[0] - 2; i++)
			{
				xi = i * dnxm1;

				for (k = 0; k <= grid_points[2] - 1; k++)
				{
					zeta = k * dnzm1;

					exact_solution(xi, eta, zeta, dtemp, 0);
					for (m = 0; m <= 4; m++)
					{
						ue[k][m] = dtemp[m];
					}

					dtpp = 1.0 / dtemp[0];

					for (m = 1; m <= 4; m++)
					{
						buf[k][m] = dtpp * dtemp[m];
					}

					cuf[k] = buf[k][3] * buf[k][3];
					buf[k][0] = cuf[k] + buf[k][1] * buf[k][1] +
							   buf[k][2] * buf[k][2];
					q[k] = 0.5 * (buf[k][1] * ue[k][1] + buf[k][2] * ue[k][2] +
								  buf[k][3] * ue[k][3]);
				}

				for (k = 1; k <= grid_points[2] - 2; k++)
				{
					km1 = k - 1;
					kp1 = k + 1;

					forcing[k][j][i][0] = forcing[k][j][i][0] -
						   tz2 * (ue[kp1][3] - ue[km1][3]) +
						   dz1tz1 * (ue[kp1][0] - 2.0 * ue[k][0] + ue[km1][0]);

					forcing[k][j][i][1] = forcing[k][j][i][1] - tz2 * (
						   ue[kp1][1] * buf[kp1][3] - ue[km1][1] * buf[km1][3]) +
						   zzcon2 * (buf[kp1][1] - 2.0 * buf[k][1] + buf[km1][1]) +
						   dz2tz1 * (ue[kp1][1] - 2.0 * ue[k][1] + ue[km1][1]);

					forcing[k][j][i][2] = forcing[k][j][i][2] - tz2 * (
						   ue[kp1][2] * buf[kp1][3] - ue[km1][2] * buf[km1][3]) +
						   zzcon2 * (buf[kp1][2] - 2.0 * buf[k][2] + buf[km1][2]) +
						   dz3tz1 * (ue[kp1][2] - 2.0 * ue[k][2] + ue[km1][2]);

					forcing[k][j][i][3] = forcing[k][j][i][3] - tz2 * (
						  (ue[kp1][3] * buf[kp1][3] + c2 * (ue[kp1][4] - q[kp1])) -
						  (ue[km1][3] * buf[km1][3] + c2 * (ue[km1][4] - q[km1]))) +
						  zzcon1 * (buf[kp1][3] - 2.0 * buf[k][3] + buf[km1][3]) +
						  dz4tz1 * (ue[kp1][3] - 2.0 * ue[k][3] + ue[km1][3]);

					forcing[k][j][i][4] = forcing[k][j][i][4] - tz2 * (
						   buf[kp1][3] * (c1 * ue[kp1][4] - c2 * q[kp1]) -
						   buf[km1][3] * (c1 * ue[km1][4] - c2 * q[km1])) +
						   0.5 * zzcon3 * (buf[kp1][0] - 2.0 * buf[k][0]
										+ buf[km1][0]) +
						   zzcon4 * (cuf[kp1] - 2.0 * cuf[k] + cuf[km1]) +
						   zzcon5 * (buf[kp1][4] - 2.0 * buf[k][4] + buf[km1][4]) +
						   dz5tz1 * (ue[kp1][4] - 2.0 * ue[k][4] + ue[km1][4]);
				}

				//---------------------------------------------------------------------
				//            Fourth-order dissipation
				//---------------------------------------------------------------------
				for (m = 0; m <= 4; m++)
				{
					k = 1;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							  (5.0 * ue[k][m] - 4.0 * ue[k+1][m] + ue[k+2][m]);
					k = 2;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							 (-4.0 * ue[k-1][m] + 6.0 * ue[k][m] -
							   4.0 * ue[k+1][m] + ue[k+2][m]);
				}

				for (m = 0; m <= 4; m++)
				{
					for (k = 3; k <= grid_points[2] - 4; k++)
					{
						forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							  (ue[k-2][m] - 4.0 * ue[k-1][m] +
							   6.0 * ue[k][m] - 4.0 * ue[k+1][m] + ue[k+2][m]);
					}
				}

				for (m = 0; m <= 4; m++)
				{
					k = grid_points[2] - 3;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
							 (ue[k-2][m] - 4.0 * ue[k-1][m] +
							  6.0 * ue[k][m] - 4.0 * ue[k+1][m]);
					k = grid_points[2] - 2;
					forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
						  (ue[k-2][m] - 4.0 * ue[k-1][m] + 5.0 * ue[k][m]);
				}
			}
		}

		//---------------------------------------------------------------------
		// now change the sign of the forcing function, 
		//---------------------------------------------------------------------
		for (k = 1; k <= grid_points[2] - 2; k++)
		{
			for (j = 1; j <= grid_points[1] - 2; j++)
			{
				for (i = 1; i <= grid_points[0] - 2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						forcing[k][j][i][m] = -1.0 * forcing[k][j][i][m];
					}
				}
			}
		}
	}

        public void x_solve()
        {
            int i, j, k, m, n, isize;

            if (timeron) timer.start(t_xsolve);
            //---------------------------------------------------------------------
            //     This function computes the left hand side in the xi-direction
            //---------------------------------------------------------------------

            isize = grid_points[0] - 1;

            //---------------------------------------------------------------------
            //     determine a (labeled f) and n jacobians
            //---------------------------------------------------------------------
            for (k = 1; k <= grid_points[2] - 2; k++)
            {
                for (j = 1; j <= grid_points[1] - 2; j++)
                {
                    for (i = 0; i <= isize; i++)
                    {

                        tmp1 = rho_i[k][j][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;
                        //---------------------------------------------------------------------
                        //---------------------------------------------------------------------
                        fjac[i][0][0] = 0.0;
                        fjac[i][1][0] = 1.0;
                        fjac[i][2][0] = 0.0;
                        fjac[i][3][0] = 0.0;
                        fjac[i][4][0] = 0.0;

                        fjac[i][0][1] = -(u[k][j][i][1]) * tmp2 *
                             u[k][j][i][1]
                             + c2 * qs[k][j][i];
                        fjac[i][1][1] = (2.0 - c2)
                             * (u[k][j][i][1] / u[k][j][i][0]);
                        fjac[i][2][1] = -c2 * (u[k][j][i][2] * tmp1);
                        fjac[i][3][1] = -c2 * (u[k][j][i][3] * tmp1);
                        fjac[i][4][1] = c2;

                        fjac[i][0][2] = -(u[k][j][i][1] * u[k][j][i][2]) * tmp2;
                        fjac[i][1][2] = u[k][j][i][2] * tmp1;
                        fjac[i][2][2] = u[k][j][i][1] * tmp1;
                        fjac[i][3][2] = 0.0;
                        fjac[i][4][2] = 0.0;

                        fjac[i][0][3] = -(u[k][j][i][1] * u[k][j][i][3]) * tmp2;
                        fjac[i][1][3] = u[k][j][i][3] * tmp1;
                        fjac[i][2][3] = 0.0;
                        fjac[i][3][3] = u[k][j][i][1] * tmp1;
                        fjac[i][4][3] = 0.0;

                        fjac[i][0][4] = (c2 * 2.0 * square[k][j][i]
                             - c1 * u[k][j][i][4])
                             * (u[k][j][i][1] * tmp2);
                        fjac[i][1][4] = c1 * u[k][j][i][4] * tmp1
                             - c2
                             * (u[k][j][i][1] * u[k][j][i][1] * tmp2
                             + qs[k][j][i]);
                        fjac[i][2][4] = -c2 * (u[k][j][i][2] * u[k][j][i][1])
                             * tmp2;
                        fjac[i][3][4] = -c2 * (u[k][j][i][3] * u[k][j][i][1])
                             * tmp2;
                        fjac[i][4][4] = c1 * (u[k][j][i][1] * tmp1);

                        njac[i][0][0] = 0.0;
                        njac[i][1][0] = 0.0;
                        njac[i][2][0] = 0.0;
                        njac[i][3][0] = 0.0;
                        njac[i][4][0] = 0.0;

                        njac[i][0][1] = -con43 * c3c4 * tmp2 * u[k][j][i][1];
                        njac[i][1][1] = con43 * c3c4 * tmp1;
                        njac[i][2][1] = 0.0;
                        njac[i][3][1] = 0.0;
                        njac[i][4][1] = 0.0;

                        njac[i][0][2] = -c3c4 * tmp2 * u[k][j][i][2];
                        njac[i][1][2] = 0.0;
                        njac[ i][2][2] = c3c4 * tmp1;
                        njac[i][3][2] = 0.0;
                        njac[i][4][2] = 0.0;

                        njac[i][0][3] = -c3c4 * tmp2 * u[k][j][i][3];
                        njac[i][1][3] = 0.0;
                        njac[i][2][3] = 0.0;
                        njac[i][3][3] = c3c4 * tmp1;
                        njac[i][4][3] = 0.0;

                        njac[i][0][4] = -(con43 * c3c4
                             - c1345) * tmp3 * (Math.Pow(u[k][j][i][1], 2))
                             - (c3c4 - c1345) * tmp3 * (Math.Pow(u[k][j][i][2], 2))
                             - (c3c4 - c1345) * tmp3 * (Math.Pow(u[k][j][i][3], 2))
                             - c1345 * tmp2 * u[k][j][i][4];

                        njac[i][1][4] = (con43 * c3c4
                             - c1345) * tmp2 * u[k][j][i][1];
                        njac[i][2][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][2];
                        njac[i][3][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][3];
                        njac[i][4][4] = (c1345) * tmp1;

                    }
                    //---------------------------------------------------------------------
                    //     now jacobians set, so form left hand side in x direction
                    //---------------------------------------------------------------------
                    lhsinit(lhs, isize);

                    for (i = 1; i <= isize - 1; i++)
                    {

                        tmp1 = dt * tx1;
                        tmp2 = dt * tx2;

                        lhs[i][aa][0][0] = -tmp2 * fjac[i - 1][0][0]
                             - tmp1 * njac[i - 1][0][0]
                             - tmp1 * dx1;
                        lhs[i][aa][1][0] = -tmp2 * fjac[i - 1][1][0]
                             - tmp1 * njac[i - 1][1][0];
                        lhs[i][aa][2][0] = -tmp2 * fjac[i - 1][2][0]
                             - tmp1 * njac[i - 1][2][0];
                        lhs[i][aa][3][0] = -tmp2 * fjac[i - 1][3][0]
                             - tmp1 * njac[i - 1][3][0];
                        lhs[i][aa][4][0] = -tmp2 * fjac[i - 1][4][0]
                             - tmp1 * njac[i - 1][4][0];

                        lhs[i][aa][0][1] = -tmp2 * fjac[i - 1][0][1]
                             - tmp1 * njac[i - 1][0][1];
                        lhs[i][aa][1][1] = -tmp2 * fjac[i - 1][1][1]
                             - tmp1 * njac[i - 1][1][1]
                             - tmp1 * dx2;
                        lhs[i][aa][2][1] = -tmp2 * fjac[i - 1][2][1]
                             - tmp1 * njac[i - 1][2][1];
                        lhs[i][aa][3][1] = -tmp2 * fjac[i - 1][3][1]
                             - tmp1 * njac[i - 1][3][1];
                        lhs[i][aa][4][1] = -tmp2 * fjac[i - 1][4][1]
                             - tmp1 * njac[i - 1][4][1];

                        lhs[i][aa][0][2] = -tmp2 * fjac[i - 1][0][2]
                             - tmp1 * njac[i - 1][0][2];
                        lhs[i][aa][1][2] = -tmp2 * fjac[i - 1][1][2]
                             - tmp1 * njac[i - 1][1][2];
                        lhs[i][aa][2][2] = -tmp2 * fjac[i - 1][2][2]
                             - tmp1 * njac[i - 1][2][2]
                             - tmp1 * dx3;
                        lhs[i][aa][3][2] = -tmp2 * fjac[i - 1][3][2]
                             - tmp1 * njac[i - 1][3][2];
                        lhs[i][aa][4][2] = -tmp2 * fjac[i - 1][4][2]
                             - tmp1 * njac[i - 1][4][2];

                        lhs[i][aa][0][3] = -tmp2 * fjac[i - 1][0][3]
                             - tmp1 * njac[i - 1][0][3];
                        lhs[i][aa][1][3] = -tmp2 * fjac[i - 1][1][3]
                             - tmp1 * njac[i - 1][1][3];
                        lhs[i][aa][2][3] = -tmp2 * fjac[i - 1][2][3]
                             - tmp1 * njac[i - 1][2][3];
                        lhs[i][aa][3][3] = -tmp2 * fjac[i - 1][3][3]
                             - tmp1 * njac[i - 1][3][3]
                             - tmp1 * dx4;
                        lhs[i][aa][4][3] = -tmp2 * fjac[i - 1][4][3]
                             - tmp1 * njac[i - 1][4][3];

                        lhs[i][aa][0][4] = -tmp2 * fjac[i - 1][0][4]
                             - tmp1 * njac[i - 1][0][4];
                        lhs[i][aa][1][4] = -tmp2 * fjac[i - 1][1][4]
                             - tmp1 * njac[i - 1][1][4];
                        lhs[i][aa][2][4] = -tmp2 * fjac[i - 1][2][4]
                             - tmp1 * njac[i - 1][2][4];
                        lhs[i][aa][3][4] = -tmp2 * fjac[i - 1][3][4]
                             - tmp1 * njac[i - 1][3][4];
                        lhs[i][aa][4][4] = -tmp2 * fjac[i - 1][4][4]
                             - tmp1 * njac[i - 1][4][4]
                             - tmp1 * dx5;

                        lhs[ i][bb][0][0] = 1.0
                             + tmp1 * 2.0 * njac[i][0][0]
                             + tmp1 * 2.0 * dx1;
                        lhs[i][bb][1][0] = tmp1 * 2.0 * njac[i][1][0];
                        lhs[i][bb][2][0] = tmp1 * 2.0 * njac[i][2][0];
                        lhs[i][bb][3][0] = tmp1 * 2.0 * njac[i][3][0];
                        lhs[i][bb][4][0] = tmp1 * 2.0 * njac[i][4][0];



                        lhs[ i][bb][0][1] = tmp1 * 2.0 * njac[i][0][1];
                        lhs[i][bb][1][1] = 1.0
                             + tmp1 * 2.0 * njac[i][1][1]
                             + tmp1 * 2.0 * dx2;
                        lhs[i][bb][2][1] = tmp1 * 2.0 * njac[i][2][1];
                        lhs[i][bb][3][1] = tmp1 * 2.0 * njac[i][3][1];
                        lhs[i][bb][4][1] = tmp1 * 2.0 * njac[i][4][1];

                        lhs[ i][bb][0][2] = tmp1 * 2.0 * njac[i][0][2];
                        lhs[i][bb][1][2] = tmp1 * 2.0 * njac[i][1][2];
                        lhs[i][bb][2][2] = 1.0
                             + tmp1 * 2.0 * njac[ i][2][2]
                             + tmp1 * 2.0 * dx3;
                        lhs[i][bb][3][2] = tmp1 * 2.0 * njac[i][3][2];
                        lhs[i][bb][4][2] = tmp1 * 2.0 * njac[i][4][2];

                        lhs[ i][bb][0][3] = tmp1 * 2.0 * njac[i][0][3];
                        lhs[i][bb][1][3] = tmp1 * 2.0 * njac[i][1][3];
                        lhs[i][bb][2][3] = tmp1 * 2.0 * njac[i][2][3];
                        lhs[i][bb][3][3] = 1.0
                             + tmp1 * 2.0 * njac[i][3][3]
                             + tmp1 * 2.0 * dx4;
                        lhs[i][bb][4][3] = tmp1 * 2.0 * njac[i][4][3];

                        lhs[ i][bb][0][4] = tmp1 * 2.0 * njac[i][0][4];
                        lhs[ i][bb][1][4] = tmp1 * 2.0 * njac[i][1][4];
                        lhs[ i][bb][2][4] = tmp1 * 2.0 * njac[i][2][4];
                        lhs[ i][bb][3][4] = tmp1 * 2.0 * njac[i][3][4];
                        lhs[ i][bb][4][4] = 1.0
                             + tmp1 * 2.0 * njac[i][4][4]
                             + tmp1 * 2.0 * dx5;


                        lhs[i][cc][0][0] = tmp2 * fjac[i+1][0][0]
                                 - tmp1 * njac[i + 1][0][0]
                                 - tmp1 * dx1;
                        lhs[i][cc][1][0] = tmp2 * fjac[i+1][1][0]
                             - tmp1 * njac[i + 1][1][0];
                        lhs[i][cc][2][0] = tmp2 * fjac[i+1][2][0]
                             - tmp1 * njac[i + 1][2][0];
                        lhs[i][cc][3][0] = tmp2 * fjac[i+1][3][0]
                             - tmp1 * njac[i + 1][3][0];
                        lhs[i][cc][4][0] = tmp2 * fjac[i+1][4][0]
                             - tmp1 * njac[i + 1][4][0];

                        lhs[i][cc][0][1] = tmp2 * fjac[i+1][0][1]
                             - tmp1 * njac[i + 1][0][1];
                        lhs[i][cc][1][1] = tmp2 * fjac[i+1][1][1]
                             - tmp1 * njac[i + 1][1][1]
                             - tmp1 * dx2;
                        lhs[i][cc][2][1] = tmp2 * fjac[i+1][2][1]
                             - tmp1 * njac[i + 1][2][1];
                        lhs[i][cc][3][1] = tmp2 * fjac[i+1][3][1]
                             - tmp1 * njac[i + 1][3][1];
                        lhs[i][cc][4][1] = tmp2 * fjac[i+1][4][1]
                             - tmp1 * njac[i + 1][4][1];

                        lhs[i][cc][0][2] = tmp2 * fjac[i+1][0][2]
                             - tmp1 * njac[i + 1][0][2];
                        lhs[i][cc][1][2] = tmp2 * fjac[i+1][1][2]
                             - tmp1 * njac[i + 1][1][2];
                        lhs[i][cc][2][2] = tmp2 * fjac[i+1][2][2]
                             - tmp1 * njac[i + 1][2][2]
                             - tmp1 * dx3;
                        lhs[i][cc][3][2] = tmp2 * fjac[i+1][3][2]
                             - tmp1 * njac[i + 1][3][2];
                        lhs[i][cc][4][2] = tmp2 * fjac[i+1][4][2]
                             - tmp1 * njac[i + 1][4][2];

                        lhs[i][cc][0][3] = tmp2 * fjac[i + 1][0][3]
                             - tmp1 * njac[i + 1][0][3];
                        lhs[i][cc][1][3] = tmp2 * fjac[i + 1][1][3]
                             - tmp1 * njac[i + 1][1][3];
                        lhs[i][cc][2][3] = tmp2 * fjac[i + 1][2][3]
                             - tmp1 * njac[i + 1][2][3];
                        lhs[i][cc][3][3] = tmp2 * fjac[i + 1][3][3]
                             - tmp1 * njac[i + 1][3][3]
                             - tmp1 * dx4;
                        lhs[i][cc][4][3] = tmp2 * fjac[i+1][4][3]
                             - tmp1 * njac[i + 1][4][3];

                        lhs[i][cc][0][4] = tmp2 * fjac[i+1][0][4]
                             - tmp1 * njac[i + 1][0][4];
                        lhs[i][cc][1][4] = tmp2 * fjac[i+1][1][4]
                             - tmp1 * njac[i + 1][1][4];
                        lhs[i][cc][2][4] = tmp2 * fjac[i+1][2][4]
                             - tmp1 * njac[i + 1][2][4];
                        lhs[i][cc][3][4] = tmp2 * fjac[i+1][3][4]
                             - tmp1 * njac[i + 1][3][4];
                        lhs[i][cc][4][4] = tmp2 * fjac[i+1][4][4]
                             - tmp1 * njac[i + 1][4][4]
                             - tmp1 * dx5;

                    }
                    //---------------------------------------------------------------------
                    //     performs guaussian elimination on this cell.
                    //     
                    //     assumes that unpacking routines for non-first cells 
                    //     preload C' and rhs' from previous cell.
                    //     
                    //     assumed send happens outside this routine, but that
                    //     c'(IMAX) and rhs'(IMAX) will be sent to next cell
                    //---------------------------------------------------------------------

                    //---------------------------------------------------------------------
                    //     outer most do loops - sweeping in i direction
                    //---------------------------------------------------------------------

                    //---------------------------------------------------------------------
                    //     multiply c(0,j,k) by b_inverse and copy back to c
                    //     multiply rhs(0) by b_inverse(0) and copy to rhs
                    //---------------------------------------------------------------------
					binvcrhs(lhs, bb, 0, lhs, cc, 0, rhs, 0, j, k);

                    //---------------------------------------------------------------------
                    //     begin inner most do loop
                    //     do all the elements of the cell unless last 
                    //---------------------------------------------------------------------
                    for (i = 1; i <= isize - 1; i++)
                    {

                        //---------------------------------------------------------------------
                        //     rhs(i) = rhs(i) - A*rhs(i-1)
                        //---------------------------------------------------------------------
                        matvec_sub(lhs, aa, i,
                                   rhs, (i - 1), j, k,
                                   rhs, i, j, k);

                        //---------------------------------------------------------------------
                        //     B(i) = B(i) - C(i-1)*A(i)
                        //---------------------------------------------------------------------
                        matmul_sub(lhs, aa, i, lhs, cc, (i - 1), lhs, bb, i);


                        //---------------------------------------------------------------------
                        //     multiply c(i,j,k) by b_inverse and copy back to c
                        //     multiply rhs(1,j,k) by b_inverse(1,j,k) and copy to rhs
                        //---------------------------------------------------------------------
                        binvcrhs(lhs, bb, i,
						         lhs, cc, i,
						         rhs, i, j, k);

                    }
                    //---------------------------------------------------------------------
                    //     rhs(isize) = rhs(isize) - A*rhs(isize-1)
                    //---------------------------------------------------------------------
                    matvec_sub(lhs, aa, isize,
					           rhs, (isize - 1), j, k,
					           rhs, isize, j, k);

                    //---------------------------------------------------------------------
                    //     B(isize) = B(isize) - C(isize-1)*A(isize)
                    //---------------------------------------------------------------------
                    matmul_sub(lhs, aa, isize,
					           lhs, cc, (isize - 1),
					           lhs, bb, isize);

                    //---------------------------------------------------------------------
                    //     multiply rhs() by b_inverse() and copy to rhs
                    //---------------------------------------------------------------------
                    binvrhs(lhs, bb, isize, rhs, isize, j, k);

                    //---------------------------------------------------------------------
                    //     back solve: if last cell, then generate U(isize)=rhs(isize)
                    //     else assume U(isize) is loaded in un pack backsub_info
                    //     so just use it
                    //     after call u(istart) will be sent to next cell
                    //---------------------------------------------------------------------

                    for (i = isize - 1; i >= 0; i--)
                    {
                        for (m = 0; m <= BLOCK_SIZE - 1; m++)
                        {
                            for (n = 0; n <= BLOCK_SIZE - 1; n++)
                            {
                                rhs[k][j][i][m] = rhs[k][j][i][m]
                                     - lhs[i][cc][n][m] * rhs[ k][ j][ i + 1][n];
                            }
                        }
                    }
                }
            }
            if (timeron) timer.stop(t_xsolve);
        }

	public void compute_rhs()
	{
		int nx2 = grid_points[0] - 2;
		int ny2 = grid_points[1] - 2;
		int nz2 = grid_points[2] - 2;

			
		int i, j, k, m;
		double rho_inv, uijk, up1, um1, vijk, vp1, vm1,
			   wijk, wp1, wm1;
		//---------------------------------------------------------------------
		//      compute the reciprocal of density, and the kinetic energy, 
		//      and the speed of sound. 
		//---------------------------------------------------------------------

		for (k = 0; k <= grid_points[2] - 1; k++)
		{
			for (j = 0; j <= grid_points[1] - 1; j++)
			{
				for (i = 0; i <= grid_points[0] - 1; i++)
				{
					rho_inv = 1.0 / u[k][j][i][0];
					rho_i[k][j][i] = rho_inv;
					us[k][j][i] = u[k][j][i][1] * rho_inv;
					vs[k][j][i] = u[k][j][i][2] * rho_inv;
					ws[k][j][i] = u[k][j][i][3] * rho_inv;
					square[k][j][i] = 0.5 * (
								  u[k][j][i][1] * u[k][j][i][1] +
								  u[k][j][i][2] * u[k][j][i][2] +
								  u[k][j][i][3] * u[k][j][i][3]) * rho_inv;
					qs[k][j][i] = square[k][j][i] * rho_inv;
				}
			}
		}

		//---------------------------------------------------------------------
		// copy the exact forcing term to the right hand side;  because 
		// this forcing term is known, we can store it on the whole grid
		// including the boundary                   
		//---------------------------------------------------------------------

		for (k = 0; k <= grid_points[2] - 1; k++)
		{
			for (j = 0; j <= grid_points[1] - 1; j++)
			{
				for (i = 0; i <= grid_points[0] - 1; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						rhs[k][j][i][m] = forcing[k][j][i][m];
					}
				}
			}
		}

		//---------------------------------------------------------------------
		//      compute xi-direction fluxes 
		//---------------------------------------------------------------------
		if (timeron) timer.start(t_rhsx);
		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					uijk = us[k][j][i];
					up1 = us[k][j][i+1];
					um1 = us[k][j][i-1];

					rhs[k][j][i][0] = rhs[k][j][i][0] + dx1tx1 *
							  (u[k][j][i+1][0] - 2.0 * u[k][j][i][0] +
							   u[k][j][i-1][0]) -
							  tx2 * (u[k][j][i+1][1] - u[k][j][i-1][1]);

					rhs[k][j][i][1] = rhs[k][j][i][1] + dx2tx1 *
							  (u[k][j][i+1][1] - 2.0 * u[k][j][i][1] +
							   u[k][j][i-1][1]) +
							  xxcon2 * con43 * (up1 - 2.0 * uijk + um1) -
							  tx2 * (u[k][j][i+1][1] * up1 -
									 u[k][j][i-1][1] * um1 +
									 (u[k][j][i+1][4] - square[k][j][i+1] -
									  u[k][j][i-1][4] + square[k][j][i-1]) *
									  c2);

					rhs[k][j][i][2] = rhs[k][j][i][2] + dx3tx1 *
							  (u[k][j][i+1][2] - 2.0 * u[k][j][i][2] +
							   u[k][j][i-1][2]) +
							  xxcon2 * (vs[k][j][i+1] - 2.0 * vs[k][j][i] +
										vs[k][j][i-1]) -
							  tx2 * (u[k][j][i+1][2] * up1 -
									 u[k][j][i-1][2] * um1);

					rhs[k][j][i][3] = rhs[k][j][i][3] + dx4tx1 *
							  (u[k][j][i+1][3] - 2.0 * u[k][j][i][3] +
							   u[k][j][i-1][3]) +
							  xxcon2 * (ws[k][j][i+1] - 2.0 * ws[k][j][i] +
										ws[k][j][i-1]) -
							  tx2 * (u[k][j][i+1][3] * up1 -
									 u[k][j][i-1][3] * um1);

					rhs[k][j][i][4] = rhs[k][j][i][4] + dx5tx1 *
							  (u[k][j][i+1][4] - 2.0 * u[k][j][i][4] +
							   u[k][j][i-1][4]) +
							  xxcon3 * (qs[k][j][i+1] - 2.0 * qs[k][j][i] +
										qs[k][j][i-1]) +
							  xxcon4 * (up1 * up1 - 2.0 * uijk * uijk +
										um1 * um1) +
							  xxcon5 * (u[k][j][i+1][4] * rho_i[k][j][i+1] -
										2.0 * u[k][j][i][4] * rho_i[k][j][i] +
										u[k][j][i-1][4] * rho_i[k][j][i-1]) -
							  tx2 * ((c1 * u[k][j][i+1][4] -
									   c2 * square[k][j][i+1]) * up1 -
									  (c1 * u[k][j][i-1][4] -
									   c2 * square[k][j][i-1]) * um1);
				}

				//---------------------------------------------------------------------
				//      add fourth order xi-direction dissipation               
				//---------------------------------------------------------------------

				i = 1;
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (5.0 * u[k][j][i][m] - 4.0 * u[k][j][i+1][m] +
									  u[k][j][i+2][m]);
				}

				i = 2;
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (-4.0 * u[k][j][i-1][m] + 6.0 * u[k][j][i][m] -
								4.0 * u[k][j][i+1][m] + u[k][j][i+2][m]);
				}

				for (i = 3; i <= nx2 - 2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							   (u[k][j][i-2][m] - 4.0 * u[k][j][i-1][m] +
								6.0 * u[k][j][i][m] - 4.0 * u[k][j][i+1][m] +
									u[k][j][i+2][m]);
					}
				}

				i = nx2 - 1;
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (u[k][j][i-2][m] - 4.0 * u[k][j][i-1][m] +
								6.0 * u[k][j][i][m] - 4.0 * u[k][j][i+1][m]);
				}

				i = nx2;
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (u[k][j][i-2][m] - 4.0 * u[k][j][i-1][m] +
								5.0 * u[k][j][i][m]);
				}
			}
		}
		if (timeron) timer.stop(t_rhsx);

		//---------------------------------------------------------------------
		//      compute eta-direction fluxes 
		//---------------------------------------------------------------------
		if (timeron) timer.start(t_rhsy);
		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					vijk = vs[k][j][i];
					vp1 = vs[k][j+1][i];
					vm1 = vs[k][j-1][i];
					rhs[k][j][i][0] = rhs[k][j][i][0] + dy1ty1 *
							 (u[k][j+1][i][0] - 2.0 * u[k][j][i][0] +
							  u[k][j-1][i][0]) -
							 ty2 * (u[k][j+1][i][2] - u[k][j-1][i][2]);
					rhs[k][j][i][1] = rhs[k][j][i][1] + dy2ty1 *
							 (u[k][j+1][i][1] - 2.0 * u[k][j][i][1] +
							  u[k][j-1][i][1]) +
							 yycon2 * (us[k][j+1][i] - 2.0 * us[k][j][i] +
									   us[k][j-1][i]) -
							 ty2 * (u[k][j+1][i][1] * vp1 -
									u[k][j-1][i][1] * vm1);
					rhs[k][j][i][2] = rhs[k][j][i][2] + dy3ty1 *
							 (u[k][j+1][i][2] - 2.0 * u[k][j][i][2] +
							  u[k][j-1][i][2]) +
							 yycon2 * con43 * (vp1 - 2.0 * vijk + vm1) -
							 ty2 * (u[k][j+1][i][2] * vp1 -
									u[k][j-1][i][2] * vm1 +
									(u[k][j+1][i][4] - square[k][j+1][i] -
									 u[k][j-1][i][4] + square[k][j-1][i])
									* c2);
					rhs[k][j][i][3] = rhs[k][j][i][3] + dy4ty1 *
							 (u[k][j+1][i][3] - 2.0 * u[k][j][i][3] +
							  u[k][j-1][i][3]) +
							 yycon2 * (ws[k][j+1][i] - 2.0 * ws[k][j][i] +
									   ws[k][j-1][i]) -
							 ty2 * (u[k][j+1][i][3] * vp1 -
									u[k][j-1][i][3] * vm1);
					rhs[k][j][i][4] = rhs[k][j][i][4] + dy5ty1 *
							 (u[k][j+1][i][4] - 2.0 * u[k][j][i][4] +
							  u[k][j-1][i][4]) +
							 yycon3 * (qs[k][j+1][i] - 2.0 * qs[k][j][i] +
									   qs[k][j-1][i]) +
							 yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk +
									   vm1 * vm1) +
							 yycon5 * (u[k][j+1][i][4] * rho_i[k][j+1][i] -
									   2.0 * u[k][j][i][4] * rho_i[k][j][i] +
									   u[k][j-1][i][4] * rho_i[k][j-1][i]) -
							 ty2 * ((c1 * u[k][j+1][i][4] -
									 c2 * square[k][j+1][i]) * vp1 -
									(c1 * u[k][j-1][i][4] -
									 c2 * square[k][j-1][i]) * vm1);
				}
			}

			//---------------------------------------------------------------------
			//      add fourth order eta-direction dissipation         
			//---------------------------------------------------------------------

			j = 1;
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (5.0 * u[k][j][i][m] - 4.0 * u[k][j+1][i][m] +
									  u[k][j+2][i][m]);
				}
			}

			j = 2;
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (-4.0 * u[k][j-1][i][m] + 6.0 * u[k][j][i][m] -
								4.0 * u[k][j+1][i][m] + u[k][j+2][i][m]);
				}
			}

			for (j = 3; j <= ny2 - 2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							   (u[k][j-2][i][m] - 4.0 * u[k][j-1][i][m] +
								6.0 * u[k][j][i][m] - 4.0 * u[k][j+1][i][m] +
									u[k][j+2][i][m]);
					}
				}
			}

			j = ny2 - 1;
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (u[k][j-2][i][m] - 4.0 * u[k][j-1][i][m] +
								6.0 * u[k][j][i][m] - 4.0 * u[k][j+1][i][m]);
				}
			}

			j = ny2;
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (u[k][j-2][i][m] - 4.0 * u[k][j-1][i][m] +
								5.0 * u[k][j][i][m]);
				}
			}
		}
		if (timeron) timer.stop(t_rhsy);

		//---------------------------------------------------------------------
		//      compute zeta-direction fluxes 
		//---------------------------------------------------------------------
		if (timeron) timer.start(t_rhsz);
		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					wijk = ws[k][j][i];
					wp1 = ws[k+1][j][i];
					wm1 = ws[k-1][j][i];

					rhs[k][j][i][0] = rhs[k][j][i][0] + dz1tz1 *
							 (u[k+1][j][i][0] - 2.0 * u[k][j][i][0] +
							  u[k-1][j][i][0]) -
							 tz2 * (u[k+1][j][i][3] - u[k-1][j][i][3]);
					rhs[k][j][i][1] = rhs[k][j][i][1] + dz2tz1 *
							 (u[k+1][j][i][1] - 2.0 * u[k][j][i][1] +
							  u[k-1][j][i][1]) +
							 zzcon2 * (us[k+1][j][i] - 2.0 * us[k][j][i] +
									   us[k-1][j][i]) -
							 tz2 * (u[k+1][j][i][1] * wp1 -
									u[k-1][j][i][1] * wm1);
					rhs[k][j][i][2] = rhs[k][j][i][2] + dz3tz1 *
							 (u[k+1][j][i][2] - 2.0 * u[k][j][i][2] +
							  u[k-1][j][i][2]) +
							 zzcon2 * (vs[k+1][j][i] - 2.0 * vs[k][j][i] +
									   vs[k-1][j][i]) -
							 tz2 * (u[k+1][j][i][2] * wp1 -
									u[k-1][j][i][2] * wm1);
					rhs[k][j][i][3] = rhs[k][j][i][3] + dz4tz1 *
							 (u[k+1][j][i][3] - 2.0 * u[k][j][i][3] +
							  u[k-1][j][i][3]) +
							 zzcon2 * con43 * (wp1 - 2.0 * wijk + wm1) -
							 tz2 * (u[k+1][j][i][3] * wp1 -
									u[k-1][j][i][3] * wm1 +
									(u[k+1][j][i][4] - square[k+1][j][i] -
									 u[k-1][j][i][4] + square[k-1][j][i])
									* c2);
					rhs[k][j][i][4] = rhs[k][j][i][4] + dz5tz1 *
							 (u[k+1][j][i][4] - 2.0 * u[k][j][i][4] +
							  u[k-1][j][i][4]) +
							 zzcon3 * (qs[k+1][j][i] - 2.0 * qs[k][j][i] +
									   qs[k-1][j][i]) +
							 zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk +
									   wm1 * wm1) +
							 zzcon5 * (u[k+1][j][i][4] * rho_i[k+1][j][i] -
									   2.0 * u[k][j][i][4] * rho_i[k][j][i] +
									   u[k-1][j][i][4] * rho_i[k-1][j][i]) -
							 tz2 * ((c1 * u[k+1][j][i][4] -
									  c2 * square[k+1][j][i]) * wp1 -
									 (c1 * u[k-1][j][i][4] -
									  c2 * square[k-1][j][i]) * wm1);
				}
			}
		}

		//---------------------------------------------------------------------
		//      add fourth order zeta-direction dissipation                
		//---------------------------------------------------------------------

		k = 1;
		for (j = 1; j <= ny2; j++)
		{
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (5.0 * u[k][j][i][m] - 4.0 * u[k+1][j][i][m] +
									  u[k+2][j][i][m]);
				}
			}
		}

		k = 2;
		for (j = 1; j <= ny2; j++)
		{
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (-4.0 * u[k-1][j][i][m] + 6.0 * u[k][j][i][m] -
								4.0 * u[k+1][j][i][m] + u[k+2][j][i][m]);
				}
			}
		}

		for (k = 3; k <= nz2 - 2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							   (u[k-2][j][i][m] - 4.0 * u[k-1][j][i][m] +
								6.0 * u[k][j][i][m] - 4.0 * u[k+1][j][i][m] +
									u[k+2][j][i][m]);
					}
				}
			}
		}

		k = nz2 - 1;
		for (j = 1; j <= ny2; j++)
		{
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (u[k-2][j][i][m] - 4.0 * u[k-1][j][i][m] +
								6.0 * u[k][j][i][m] - 4.0 * u[k+1][j][i][m]);
				}
			}
		}

		k = nz2;
		for (j = 1; j <= ny2; j++)
		{
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] - dssp *
							  (u[k-2][j][i][m] - 4.0 * u[k-1][j][i][m] +
								5.0 * u[k][j][i][m]);
				}
			}
		}
		if (timeron) timer.stop(t_rhsz);


        for (m = 0; m <= 4; m++)
        {
            for (k = 1; k <= nz2; k++)
            {

                for (j = 1; j <= ny2; j++)
                {
                    for (i = 1; i <= nx2; i++)
                    {
                        rhs[k][j][i][m] = rhs[k][j][i][m] * dt;

                    }
                }
            }
        }
		
	}

        public void print_lhs()
        {
            double count1 = 0, count2 = 0, count3 = 0;
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    for (int m = 0; m < problem_size + 1; m++)
                    {
                        count2 += njac[m][j][i];
                        count3 += fjac[m][j][i];
                        for (int k = 0; k < 3; k++)
                        {
                            count1 += lhs[m][k][j][i];
                        }
                    }
                }
            }
            Console.WriteLine("lhs checksum is: ");
            Console.WriteLine(count1);
            Console.WriteLine("fjac checksum is: ");
            Console.WriteLine(count3);
            Console.WriteLine("njac checksum is: ");
            Console.WriteLine(count2);
        }

        public void y_solve()
        {
            int i, j, k, m, n, jsize;

            if (timeron) timer.start(t_ysolve);

            //---------------------------------------------------------------------
            //     This function computes the left hand side for the three y-factors   
            //---------------------------------------------------------------------

            jsize = grid_points[1] - 1;

            //---------------------------------------------------------------------
            //     Compute the indices for storing the tri-diagonal matrix;
            //     determine a (labeled f) and n jacobians for cell c
            //---------------------------------------------------------------------
            for (k = 1; k <= grid_points[2] - 2; k++)
            {
                for (i = 1; i <= grid_points[0] - 2; i++)
                {
                    for (j = 0; j <= jsize; j++)
                    {

                        tmp1 = rho_i[k][j][i];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        fjac[ j][0][0] = 0.0;
                        fjac[ j][1][0] = 0.0;
                        fjac[ j][2][0] = 1.0;
                        fjac[ j][3][0] = 0.0;
                        fjac[ j][4][0] = 0.0;

                        fjac[j][0][1] = -(u[k][j][i][1] * u[k][j][i][2])
                             * tmp2;
                        fjac[j][1][1] = u[k][j][i][2] * tmp1;
                        fjac[j][2][1] = u[k][j][i][1] * tmp1;
                        fjac[j][3][1] = 0.0;
                        fjac[ j][4][1] = 0.0;

                        fjac[j][0][2] = -(u[k][j][i][2] * u[k][j][i][2] * tmp2)
                             + c2 * qs[k][j][i];
                        fjac[j][1][2] = -c2 * u[k][j][i][1] * tmp1;
                        fjac[j][2][2] = (2.0 - c2)
                             * u[k][j][i][2] * tmp1;
                        fjac[j][3][2] = -c2 * u[k][j][i][3] * tmp1;
                        fjac[j][4][2] = c2;

                        fjac[j][0][3] = -(u[k][j][i][2] * u[k][j][i][3])
                             * tmp2;
                        fjac[j][1][3] = 0.0;
                        fjac[j][2][3] = u[k][j][i][3] * tmp1;
                        fjac[j][3][3] = u[k][j][i][2] * tmp1;
                        fjac[j][4][3] = 0.0;

                        fjac[j][0][4] = (c2 * 2.0 * square[k][j][i]
                             - c1 * u[k][j][i][4])
                             * u[k][j][i][2] * tmp2;
                        fjac[j][1][4] = -c2 * u[k][j][i][1] * u[k][j][i][2]
                             * tmp2;
                        fjac[j][2][4] = c1 * u[k][j][i][4] * tmp1
                             - c2
                             * (qs[k][j][i]
                             + u[k][j][i][2] * u[k][j][i][2] * tmp2);
                        fjac[j][3][4] = -c2 * (u[k][j][i][2] * u[k][j][i][3])
                             * tmp2;
                        fjac[j][4][4] = c1 * u[k][j][i][2] * tmp1;

                        njac[j][0][0] = 0.0;
                        njac[j][1][0] = 0.0;
                        njac[j][2][0] = 0.0;
                        njac[j][3][0] = 0.0;
                        njac[j][4][0] = 0.0;

                        njac[j][0][1] = -c3c4 * tmp2 * u[k][j][i][1];
                        njac[j][1][1] = c3c4 * tmp1;
                        njac[j][2][1] = 0.0;
                        njac[j][3][1] = 0.0;
                        njac[j][4][1] = 0.0;

                        njac[j][0][2] = -con43 * c3c4 * tmp2 * u[k][j][i][2];
                        njac[j][1][2] = 0.0;
                        njac[j][2][2] = con43 * c3c4 * tmp1;
                        njac[j][3][2] = 0.0;
                        njac[j][4][2] = 0.0;

                        njac[j][0][3] = -c3c4 * tmp2 * u[k][j][i][3];
                        njac[j][1][3] = 0.0;
                        njac[j][2][3] = 0.0;
                        njac[j][3][3] = c3c4 * tmp1;
                        njac[j][4][3] = 0.0;

                        njac[j][0][4] = -(c3c4
                             - c1345) * tmp3 * (Math.Pow(u[k][j][i][1], 2))
                             - (con43 * c3c4
                             - c1345) * tmp3 * (Math.Pow(u[k][j][i][2], 2))
                             - (c3c4 - c1345) * tmp3 * (Math.Pow(u[k][j][i][3], 2))
                             - c1345 * tmp2 * u[k][j][i][4];

                        njac[j][1][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][1];
                        njac[j][2][4] = (con43 * c3c4
                             - c1345) * tmp2 * u[k][j][i][2];
                        njac[j][3][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][3];
                        njac[j][4][4] = (c1345) * tmp1;
                    }

                    //---------------------------------------------------------------------
                    //     now joacobians set, so form left hand side in y direction
                    //---------------------------------------------------------------------
                    lhsinit(lhs, jsize);
                    for (j = 1; j <= jsize - 1; j++)
                    {

                        tmp1 = dt * ty1;
                        tmp2 = dt * ty2;


                        lhs[ j][aa][0][0] = -tmp2 * fjac[j-1][0][0]
                             - tmp1 * njac[j-1][0][0]
                             - tmp1 * dy1;
                        lhs[j][aa][1][0] = -tmp2 * fjac[j-1][1][0]
                             - tmp1 * njac[j-1][1][0];
                        lhs[j][aa][2][0] = -tmp2 * fjac[j-1][2][0]
                             - tmp1 * njac[j-1][2][0];
                        lhs[j][aa][3][0] = -tmp2 * fjac[j-1][3][0]
                             - tmp1 * njac[j-1][3][0];
                        lhs[j][aa][4][0] = -tmp2 * fjac[j-1][4][0]
                             - tmp1 * njac[j-1][4][0];

                        lhs[j][aa][0][1] = -tmp2 * fjac[j-1][0][1]
                             - tmp1 * njac[j-1][0][1];
                        lhs[j][aa][1][1] = -tmp2 * fjac[j-1][1][1]
                             - tmp1 * njac[j-1][1][1]
                             - tmp1 * dy2;
                        lhs[j][aa][2][1] = -tmp2 * fjac[j-1][2][1]
                             - tmp1 * njac[j-1][2][1];
                        lhs[j][aa][3][1] = -tmp2 * fjac[j-1][3][1]
                             - tmp1 * njac[j-1][3][1];
                        lhs[j][aa][4][1] = -tmp2 * fjac[j-1][4][1]
                             - tmp1 * njac[j-1][4][1];


                        lhs[j][aa][0][2] = -tmp2 * fjac[j-1][0][2]
                             - tmp1 * njac[j-1][0][2];
                        lhs[j][aa][1][2] = -tmp2 * fjac[j-1][1][2]
                             - tmp1 * njac[j-1][1][2];
                        lhs[j][aa][2][2] = -tmp2 * fjac[j-1][2][2]
                             - tmp1 * njac[j-1][2][2]
                             - tmp1 * dy3;
                        lhs[j][aa][3][2] = -tmp2 * fjac[j-1][3][2]
                             - tmp1 * njac[j-1][3][2];
                        lhs[j][aa][4][2] = -tmp2 * fjac[j-1][4][2]
                             - tmp1 * njac[j-1][4][2];


                        lhs[j][aa][0][3] = -tmp2 * fjac[j-1][0][3]
                             - tmp1 * njac[j-1][0][3];
                        lhs[j][aa][1][3] = -tmp2 * fjac[j-1][1][3]
                             - tmp1 * njac[j-1][1][3];
                        lhs[j][aa][2][3] = -tmp2 * fjac[j-1][2][3]
                             - tmp1 * njac[j-1][2][3];
                        lhs[j][aa][3][3] = -tmp2 * fjac[j-1][3][3]
                             - tmp1 * njac[j-1][3][3]
                             - tmp1 * dy4;
                        lhs[j][aa][4][3] = -tmp2 * fjac[j-1][4][3]
                             - tmp1 * njac[j-1][4][3];



                        lhs[j][aa][0][4] = -tmp2 * fjac[j-1][0][4]
                             - tmp1 * njac[j-1][0][4];
                        lhs[j][aa][1][4] = -tmp2 * fjac[j-1][1][4]
                             - tmp1 * njac[j-1][1][4];
                        lhs[j][aa][2][4] = -tmp2 * fjac[j-1][2][4]
                             - tmp1 * njac[j-1][2][4];
                        lhs[j][aa][3][4] = -tmp2 * fjac[j-1][3][4]
                             - tmp1 * njac[j-1][3][4];
                        lhs[j][aa][4][4] = -tmp2 * fjac[j-1][4][4]
                             - tmp1 * njac[j-1][4][4]
                             - tmp1 * dy5;

                        lhs[j][bb][0][0] = 1.0
                             + tmp1 * 2.0 * njac[j][0][0]
                             + tmp1 * 2.0 * dy1;
                        lhs[j][bb][1][0] = tmp1 * 2.0 * njac[j][1][0];
                        lhs[j][bb][2][0] = tmp1 * 2.0 * njac[j][2][0];
                        lhs[j][bb][3][0] = tmp1 * 2.0 * njac[j][3][0];
                        lhs[j][bb][4][0] = tmp1 * 2.0 * njac[j][4][0];

                        lhs[j][bb][0][1] = tmp1 * 2.0 * njac[j][0][1];
                        lhs[j][bb][1][1] = 1.0
                             + tmp1 * 2.0 * njac[j][1][1]
                             + tmp1 * 2.0 * dy2;
                        lhs[j][bb][2][1] = tmp1 * 2.0 * njac[j][2][1];
                        lhs[j][bb][3][1] = tmp1 * 2.0 * njac[j][3][1];
                        lhs[j][bb][4][1] = tmp1 * 2.0 * njac[j][4][1];

                        lhs[j][bb][0][2] = tmp1 * 2.0 * njac[j][0][2];
                        lhs[j][bb][1][2] = tmp1 * 2.0 * njac[j][1][2];
                        lhs[j][bb][2][2] = 1.0
                             + tmp1 * 2.0 * njac[j][2][2]
                             + tmp1 * 2.0 * dy3;
                        lhs[j][bb][3][2] = tmp1 * 2.0 * njac[j][3][2];
                        lhs[j][bb][4][2] = tmp1 * 2.0 * njac[j][4][2];

                        lhs[j][bb][0][3] = tmp1 * 2.0 * njac[j][0][3];
                        lhs[j][bb][1][3] = tmp1 * 2.0 * njac[j][1][3];
                        lhs[j][bb][2][3] = tmp1 * 2.0 * njac[j][2][3];
                        lhs[j][bb][3][3] = 1.0
                             + tmp1 * 2.0 * njac[j][3][3]
                             + tmp1 * 2.0 * dy4;
                        lhs[j][bb][4][3] = tmp1 * 2.0 * njac[j][4][3];

                        lhs[j][bb][0][4] = tmp1 * 2.0 * njac[j][0][4];
                        lhs[j][bb][1][4] = tmp1 * 2.0 * njac[j][1][4];
                        lhs[j][bb][2][4] = tmp1 * 2.0 * njac[j][2][4];
                        lhs[j][bb][3][4] = tmp1 * 2.0 * njac[j][3][4];
                        lhs[j][bb][4][4] = 1.0
                             + tmp1 * 2.0 * njac[j][4][4]
                             + tmp1 * 2.0 * dy5;

                        lhs[j][cc][0][0] = tmp2 * fjac[j + 1][0][0]
                             - tmp1 * njac[j + 1][0][0]
                             - tmp1 * dy1;
                        lhs[j][cc][1][0] = tmp2 * fjac[j + 1][1][0]
                             - tmp1 * njac[j + 1][1][0];
                        lhs[j][cc][2][0] = tmp2 * fjac[j + 1][2][0]
                             - tmp1 * njac[j + 1][2][0];
                        lhs[j][cc][3][0] = tmp2 * fjac[j + 1][3][0]
                             - tmp1 * njac[j + 1][3][0];
                        lhs[j][cc][4][0] = tmp2 * fjac[j + 1][4][0]
                             - tmp1 * njac[j + 1][4][0];

                        lhs[j][cc][0][1] = tmp2 * fjac[j + 1][0][1]
                             - tmp1 * njac[j + 1][0][1];
                        lhs[j][cc][1][1] = tmp2 * fjac[j + 1][1][1]
                             - tmp1 * njac[j + 1][1][1]
                             - tmp1 * dy2;
                        lhs[j][cc][2][1] = tmp2 * fjac[j + 1][2][1]
                             - tmp1 * njac[j + 1][2][1];
                        lhs[j][cc][3][1] = tmp2 * fjac[j + 1][3][1]
                             - tmp1 * njac[j + 1][3][1];
                        lhs[j][cc][4][1] = tmp2 * fjac[j + 1][4][1]
                             - tmp1 * njac[j + 1][4][1];

                        lhs[j][cc][0][2] = tmp2 * fjac[j + 1][0][2]
                             - tmp1 * njac[j + 1][0][2];
                        lhs[j][cc][1][2] = tmp2 * fjac[j + 1][1][2]
                             - tmp1 * njac[j + 1][1][2];
                        lhs[j][cc][2][2] = tmp2 * fjac[j + 1][2][2]
                             - tmp1 * njac[j + 1][2][2]
                             - tmp1 * dy3;
                        lhs[j][cc][3][2] = tmp2 * fjac[j + 1][3][2]
                             - tmp1 * njac[j + 1][3][2];
                        lhs[j][cc][4][2] = tmp2 * fjac[j + 1][4][2]
                             - tmp1 * njac[j + 1][4][2];

                        lhs[j][cc][0][3] = tmp2 * fjac[j + 1][0][3]
                             - tmp1 * njac[j + 1][0][3];
                        lhs[j][cc][1][3] = tmp2 * fjac[j + 1][1][3]
                             - tmp1 * njac[j + 1][1][3];
                        lhs[j][cc][2][3] = tmp2 * fjac[j + 1][2][3]
                             - tmp1 * njac[j + 1][2][3];
                        lhs[j][cc][3][3] = tmp2 * fjac[j + 1][3][3]
                             - tmp1 * njac[j + 1][3][3]
                             - tmp1 * dy4;
                        lhs[j][cc][4][3] = tmp2 * fjac[j + 1][4][3]
                             - tmp1 * njac[j + 1][4][3];

                        lhs[j][cc][0][4] = tmp2 * fjac[j + 1][0][4]
                             - tmp1 * njac[j + 1][0][4];
                        lhs[j][cc][1][4] = tmp2 * fjac[j + 1][1][4]
                             - tmp1 * njac[j + 1][1][4];
                        lhs[j][cc][2][4] = tmp2 * fjac[j + 1][2][4]
                             - tmp1 * njac[j + 1][2][4];
                        lhs[j][cc][3][4] = tmp2 * fjac[j + 1][3][4]
                             - tmp1 * njac[j + 1][3][4];
                        lhs[j][cc][4][4] = tmp2 * fjac[j + 1][4][4]
                             - tmp1 * njac[j + 1][4][4]
                             - tmp1 * dy5;
                    }
                    //---------------------------------------------------------------------
                    //     performs guaussian elimination on this cell.
                    //     
                    //     assumes that unpacking routines for non-first cells 
                    //     preload C' and rhs' from previous cell.
                    //     
                    //     assumed send happens outside this routine, but that
                    //     c'(JMAX) and rhs'(JMAX) will be sent to next cell
                    //---------------------------------------------------------------------

                    //---------------------------------------------------------------------
                    //     multiply c(i,0,k) by b_inverse and copy back to c
                    //     multiply rhs(0) by b_inverse(0) and copy to rhs
                    //---------------------------------------------------------------------
                    binvcrhs(lhs, bb, 0,
					         lhs, cc, 0,
					         rhs, i, 0, k);

                    //---------------------------------------------------------------------
                    //     begin inner most do loop
                    //     do all the elements of the cell unless last 
                    //---------------------------------------------------------------------
                    for (j = 1; j <= jsize - 1; j++)
                    {

                        //---------------------------------------------------------------------
                        //     subtract A*lhs_vector(j-1) from lhs_vector(j)
                        //     
                        //     rhs(j) = rhs(j) - A*rhs(j-1)
                        //---------------------------------------------------------------------
                        matvec_sub(lhs, aa, j,
                                   rhs, i, (j - 1), k,
						           rhs, i, j, k);

                        //---------------------------------------------------------------------
                        //     B(j) = B(j) - C(j-1)*A(j)
                        //---------------------------------------------------------------------
                        matmul_sub(lhs, aa, j,  lhs, cc, (j - 1), lhs, bb, j);

                        //---------------------------------------------------------------------
                        //     multiply c(i,j,k) by b_inverse and copy back to c
                        //     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
                        //---------------------------------------------------------------------
                        binvcrhs(lhs, bb, j,
						         lhs, cc, j,
						         rhs, i, j, k);
                    }

                    //---------------------------------------------------------------------
                    //     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
                    //---------------------------------------------------------------------
                    matvec_sub(lhs, aa, jsize,
					           rhs, i, (jsize - 1), k,
					           rhs, i, jsize, k);

                    //---------------------------------------------------------------------
                    //     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
                    //       matmul_sub(aa,i,jsize,k,c,
                    //     $              cc,i,jsize-1,k,c,bb,i,jsize,k)
                    //---------------------------------------------------------------------
                    matmul_sub(lhs, aa, jsize,
					           lhs, cc, (jsize - 1),
					           lhs, bb, jsize);

                    //---------------------------------------------------------------------
                    //     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
                    //---------------------------------------------------------------------
                    binvrhs(lhs, bb, jsize,
					        rhs, i, jsize, k);

                    //---------------------------------------------------------------------
                    //     back solve: if last cell, then generate U(jsize)=rhs(jsize)
                    //     else assume U(jsize) is loaded in un pack backsub_info
                    //     so just use it
                    //     after   u(jstart) will be sent to next cell
                    //---------------------------------------------------------------------

                    for (j = jsize - 1; j >= 0; j--)
                    {
                        for (m = 0; m <= BLOCK_SIZE - 1; m++)
                        {
                            for (n = 0; n <= BLOCK_SIZE - 1; n++)
                            {
                                rhs[k][j][i][m] = rhs[k][j][i][m]
                                     - lhs[ j][cc][n][m] * rhs[k][j+1][i][n];
                            }
                        }
                    }
                }
            }
            if (timeron) timer.stop(t_ysolve);
        }

        public void z_solve()
        {
            int i, j, k, m, n, ksize;

            if (timeron) timer.start(t_zsolve);

            //---------------------------------------------------------------------
            //     This function computes the left hand side for the three z-factors   
            //---------------------------------------------------------------------

            ksize = grid_points[2] - 1;

            //---------------------------------------------------------------------
            //     Compute the indices for storing the block-diagonal matrix;
            //     determine c (labeled f) and s jacobians
            //---------------------------------------------------------------------
            for (j = 1; j <= grid_points[1] - 2; j++)
            {
                for (i = 1; i <= grid_points[0] - 2; i++)
                {
                    for (k = 0; k <= ksize; k++)
                    {

                        tmp1 = 1.0 / u[k][j][i][0];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        fjac[ k][0][0] = 0.0;
                        fjac[k][1][0] = 0.0;
                        fjac[k][2][0] = 0.0;
                        fjac[k][3][0] = 1.0;
                        fjac[k][4][0] = 0.0;

                        fjac[k][0][1] = -(u[k][j][i][1] * u[k][j][i][3])
                             * tmp2;
                        fjac[k][1][1] = u[k][j][i][3] * tmp1;
                        fjac[k][2][1] = 0.0;
                        fjac[k][3][1] = u[k][j][i][1] * tmp1;
                        fjac[k][4][1] = 0.0;

                        fjac[k][0][2] = -(u[k][j][i][2] * u[k][j][i][3])
                             * tmp2;
                        fjac[k][1][2] = 0.0;
                        fjac[k][2][2] = u[k][j][i][3] * tmp1;
                        fjac[k][3][2] = u[k][j][i][2] * tmp1;
                        fjac[k][4][2] = 0.0;

                        fjac[k][0][3] = -(u[k][j][i][3] * u[k][j][i][3] * tmp2)
                             + c2 * qs[k][j][i];
                        fjac[k][1][3] = -c2 * u[k][j][i][1] * tmp1;
                        fjac[k][2][3] = -c2 * u[k][j][i][2] * tmp1;
                        fjac[k][3][3] = (2.0 - c2)
                             * u[k][j][i][3] * tmp1;
                        fjac[k][4][3] = c2;

                        fjac[k][0][4] = (c2 * 2.0 * square[k][j][i]
                                 - c1 * u[k][j][i][4])
                                 * u[k][j][i][3] * tmp2;
                        fjac[k][1][4] = -c2 * (u[k][j][i][1] * u[k][j][i][3])
                             * tmp2;
                        fjac[k][2][4] = -c2 * (u[k][j][i][2] * u[k][j][i][3])
                             * tmp2;
                        fjac[k][3][4] = c1 * (u[k][j][i][4] * tmp1)
                             - c2
                             * (qs[k][j][i]
                             + u[k][j][i][3] * u[k][j][i][3] * tmp2);
                        fjac[k][4][4] = c1 * u[k][j][i][3] * tmp1;

                        njac[k][0][0] = 0.0;
                        njac[k][1][0] = 0.0;
                        njac[k][2][0] = 0.0;
                        njac[k][3][0] = 0.0;
                        njac[k][4][0] = 0.0;

                        njac[k][0][1] = -c3c4 * tmp2 * u[k][j][i][1];
                        njac[k][1][1] = c3c4 * tmp1;
                        njac[k][2][1] = 0.0;
                        njac[k][3][1] = 0.0;
                        njac[k][4][1] = 0.0;

                        njac[k][0][2] = -c3c4 * tmp2 * u[k][j][i][2];
                        njac[k][1][2] = 0.0;
                        njac[k][2][2] = c3c4 * tmp1;
                        njac[k][3][2] = 0.0;
                        njac[k][4][2] = 0.0;

                        njac[k][0][3] = -con43 * c3c4 * tmp2 * u[k][j][i][3];
                        njac[k][1][3] = 0.0;
                        njac[k][2][3] = 0.0;
                        njac[k][3][3] = con43 * c3 * c4 * tmp1;
                        njac[k][4][3] = 0.0;

                        njac[k][0][4] = -(c3c4
                             - c1345) * tmp3 * (Math.Pow(u[k][j][i][1], 2))
                             - (c3c4 - c1345) * tmp3 * (Math.Pow(u[k][j][i][2], 2))
                             - (con43 * c3c4
                             - c1345) * tmp3 * (Math.Pow(u[k][j][i][3], 2))
                             - c1345 * tmp2 * u[k][j][i][4];

                        njac[k][1][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][1];
                        njac[k][2][4] = (c3c4 - c1345) * tmp2 * u[k][j][i][2];
                        njac[k][3][4] = (con43 * c3c4
                             - c1345) * tmp2 * u[k][j][i][3];
                        njac[k][4][4] = (c1345) * tmp1;
                    }

                    //---------------------------------------------------------------------
                    //     now jacobians set, so form left hand side in z direction
                    //---------------------------------------------------------------------
                    lhsinit(lhs, ksize);
                    for (k = 1; k <= ksize - 1; k++)
                    {

                        tmp1 = dt * tz1;
                        tmp2 = dt * tz2;

                        lhs[k][aa][0][0] = -tmp2 * fjac[k-1][0][0]
                                 - tmp1 * njac[k-1][0][0]
                                 - tmp1 * dz1;
                        lhs[k][aa][1][0] = -tmp2 * fjac[k-1][1][0]
                             - tmp1 * njac[k-1][1][0];
                        lhs[k][aa][2][0] = -tmp2 * fjac[k-1][2][0]
                             - tmp1 * njac[k-1][2][0];
                        lhs[k][aa][3][0] = -tmp2 * fjac[k-1][3][0]
                             - tmp1 * njac[k-1][3][0];
                        lhs[k][aa][4][0] = -tmp2 * fjac[k-1][4][0]
                             - tmp1 * njac[k-1][4][0];

                        lhs[k][aa][0][1] = -tmp2 * fjac[k-1][0][1]
                             - tmp1 * njac[k-1][0][1];
                        lhs[k][aa][1][1] = -tmp2 * fjac[k-1][1][1]
                             - tmp1 * njac[k-1][1][1]
                             - tmp1 * dz2;
                        lhs[k][aa][2][1] = -tmp2 * fjac[k-1][2][1]
                             - tmp1 * njac[k-1][2][1];
                        lhs[k][aa][3][1] = -tmp2 * fjac[k-1][3][1]
                             - tmp1 * njac[k-1][3][1];
                        lhs[k][aa][4][1] = -tmp2 * fjac[k-1][4][1]
                             - tmp1 * njac[k-1][4][1];

                        lhs[k][aa][0][2] = -tmp2 * fjac[k-1][0][2]
                             - tmp1 * njac[k-1][0][2];
                        lhs[k][aa][1][2] = -tmp2 * fjac[k-1][1][2]
                             - tmp1 * njac[k-1][1][2];
                        lhs[k][aa][2][2] = -tmp2 * fjac[k-1][2][2]
                             - tmp1 * njac[k-1][2][2]
                             - tmp1 * dz3;
                        lhs[k][aa][3][2] = -tmp2 * fjac[k-1][3][2]
                             - tmp1 * njac[k-1][3][2];
                        lhs[k][aa][4][2] = -tmp2 * fjac[k-1][4][2]
                             - tmp1 * njac[k-1][4][2];

                        lhs[k][aa][0][3] = -tmp2 * fjac[k-1][0][3]
                             - tmp1 * njac[k-1][0][3];
                        lhs[k][aa][1][3] = -tmp2 * fjac[k-1][1][3]
                             - tmp1 * njac[k-1][1][3];
                        lhs[k][aa][2][3] = -tmp2 * fjac[k-1][2][3]
                             - tmp1 * njac[k-1][2][3];
                        lhs[k][aa][3][3] = -tmp2 * fjac[k-1][3][3]
                             - tmp1 * njac[k-1][3][3]
                             - tmp1 * dz4;
                        lhs[k][aa][4][3] = -tmp2 * fjac[k-1][4][3]
                             - tmp1 * njac[k-1][4][3];

                        lhs[k][aa][0][4] = -tmp2 * fjac[k-1][0][4]
                             - tmp1 * njac[k-1][0][4];
                        lhs[k][aa][1][4] = -tmp2 * fjac[k-1][1][4]
                             - tmp1 * njac[k-1][1][4];
                        lhs[k][aa][2][4] = -tmp2 * fjac[k-1][2][4]
                             - tmp1 * njac[k-1][2][4];
                        lhs[k][aa][3][4] = -tmp2 * fjac[k-1][3][4]
                             - tmp1 * njac[k-1][3][4];
                        lhs[k][aa][4][4] = -tmp2 * fjac[k-1][4][4]
                             - tmp1 * njac[k-1][4][4]
                             - tmp1 * dz5;

                        lhs[k][bb][0][0] = 1.0
                             + tmp1 * 2.0 * njac[k][0][0]
                             + tmp1 * 2.0 * dz1;
                        lhs[k][bb][1][0] = tmp1 * 2.0 * njac[k][1][0];
                        lhs[k][bb][2][0] = tmp1 * 2.0 * njac[k][2][0];
                        lhs[k][bb][3][0] = tmp1 * 2.0 * njac[k][3][0];
                        lhs[k][bb][4][0] = tmp1 * 2.0 * njac[k][4][0];

                        lhs[k][bb][0][1] = tmp1 * 2.0 * njac[k][0][1];
                        lhs[k][bb][1][1] = 1.0
                             + tmp1 * 2.0 * njac[k][1][1]
                             + tmp1 * 2.0 * dz2;
                        lhs[k][bb][2][1] = tmp1 * 2.0 * njac[k][2][1];
                        lhs[k][bb][3][1] = tmp1 * 2.0 * njac[k][3][1];
                        lhs[k][bb][4][1] = tmp1 * 2.0 * njac[k][4][1];

                        lhs[k][bb][0][2] = tmp1 * 2.0 * njac[k][0][2];
                        lhs[k][bb][1][2] = tmp1 * 2.0 * njac[k][1][2];
                        lhs[k][bb][2][2] = 1.0
                             + tmp1 * 2.0 * njac[k][2][2]
                             + tmp1 * 2.0 * dz3;
                        lhs[k][bb][3][2] = tmp1 * 2.0 * njac[k][3][2];
                        lhs[k][bb][4][2] = tmp1 * 2.0 * njac[k][4][2];

                        lhs[k][bb][0][3] = tmp1 * 2.0 * njac[k][0][3];
                        lhs[k][bb][1][3] = tmp1 * 2.0 * njac[k][1][3];
                        lhs[k][bb][2][3] = tmp1 * 2.0 * njac[k][2][3];
                        lhs[k][bb][3][3] = 1.0
                             + tmp1 * 2.0 * njac[k][3][3]
                             + tmp1 * 2.0 * dz4;
                        lhs[k][bb][4][3] = tmp1 * 2.0 * njac[k][4][3];

                        lhs[k][bb][0][4] = tmp1 * 2.0 * njac[k][0][4];
                        lhs[k][bb][1][4] = tmp1 * 2.0 * njac[k][1][4];
                        lhs[k][bb][2][4] = tmp1 * 2.0 * njac[k][2][4];
                        lhs[k][bb][3][4] = tmp1 * 2.0 * njac[k][3][4];
                        lhs[k][bb][4][4] = 1.0
                             + tmp1 * 2.0 * njac[k][4][4]
                             + tmp1 * 2.0 * dz5;

                        lhs[k][cc][0][0] = tmp2 * fjac[k + 1][0][0]
                             - tmp1 * njac[k+1][0][0]
                             - tmp1 * dz1;
                        lhs[k][cc][1][0] = tmp2 * fjac[k + 1][1][0]
                             - tmp1 * njac[k+1][1][0];
                        lhs[k][cc][2][0] = tmp2 * fjac[k + 1][2][0]
                             - tmp1 * njac[k+1][2][0];
                        lhs[k][cc][3][0] = tmp2 * fjac[k + 1][3][0]
                             - tmp1 * njac[k+1][3][0];
                        lhs[k][cc][4][0] = tmp2 * fjac[k + 1][4][0]
                             - tmp1 * njac[k+1][4][0];

                        lhs[k][cc][0][1] = tmp2 * fjac[k + 1][0][1]
                             - tmp1 * njac[k+1][0][1];
                        lhs[k][cc][1][1] = tmp2 * fjac[k + 1][1][1]
                             - tmp1 * njac[k+1][1][1]
                             - tmp1 * dz2;
                        lhs[k][cc][2][1] = tmp2 * fjac[k + 1][2][1]
                             - tmp1 * njac[k+1][2][1];
                        lhs[k][cc][3][1] = tmp2 * fjac[k + 1][3][1]
                             - tmp1 * njac[k+1][3][1];
                        lhs[k][cc][4][1] = tmp2 * fjac[k + 1][4][1]
                             - tmp1 * njac[k+1][4][1];

                        lhs[k][cc][0][2] = tmp2 * fjac[k + 1][0][2]
                             - tmp1 * njac[k+1][0][2];
                        lhs[k][cc][1][2] = tmp2 * fjac[k + 1][1][2]
                             - tmp1 * njac[k+1][1][2];
                        lhs[k][cc][2][2] = tmp2 * fjac[k + 1][2][2]
                             - tmp1 * njac[k+1][2][2]
                             - tmp1 * dz3;
                        lhs[k][cc][3][2] = tmp2 * fjac[k + 1][3][2]
                             - tmp1 * njac[k+1][3][2];
                        lhs[k][cc][4][2] = tmp2 * fjac[k + 1][4][2]
                             - tmp1 * njac[k+1][4][2];

                        lhs[k][cc][0][3] = tmp2 * fjac[k + 1][0][3]
                             - tmp1 * njac[k+1][0][3];
                        lhs[k][cc][1][3] = tmp2 * fjac[k + 1][1][3]
                             - tmp1 * njac[k+1][1][3];
                        lhs[k][cc][2][3] = tmp2 * fjac[k + 1][2][3]
                             - tmp1 * njac[k+1][2][3];
                        lhs[k][cc][3][3] = tmp2 * fjac[k + 1][3][3]
                             - tmp1 * njac[k+1][3][3]
                             - tmp1 * dz4;
                        lhs[k][cc][4][3] = tmp2 * fjac[k + 1][4][3]
                             - tmp1 * njac[k+1][4][3];

                        lhs[k][cc][0][4] = tmp2 * fjac[k + 1][0][4]
                             - tmp1 * njac[k+1][0][4];
                        lhs[k][cc][1][4] = tmp2 * fjac[k + 1][1][4]
                             - tmp1 * njac[k+1][1][4];
                        lhs[k][cc][2][4] = tmp2 * fjac[k + 1][2][4]
                             - tmp1 * njac[k+1][2][4];
                        lhs[k][cc][3][4] = tmp2 * fjac[k + 1][3][4]
                             - tmp1 * njac[k+1][3][4];
                        lhs[k][cc][4][4] = tmp2 * fjac[k + 1][4][4]
                             - tmp1 * njac[k+1][4][4]
                             - tmp1 * dz5;
                    }

                    //---------------------------------------------------------------------
                    //     performs guaussian elimination on this cell.
                    //     
                    //     assumes that unpacking routines for non-first cells 
                    //     preload C' and rhs' from previous cell.
                    //     
                    //     assumed send happens outside this routine, but that
                    //     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
                    //---------------------------------------------------------------------

                    //---------------------------------------------------------------------
                    //     outer most do loops - sweeping in i direction
                    //---------------------------------------------------------------------

                    //---------------------------------------------------------------------
                    //     multiply c(i,j,0) by b_inverse and copy back to c
                    //     multiply rhs(0) by b_inverse(0) and copy to rhs
                    //---------------------------------------------------------------------
                    binvcrhs(lhs, bb, 0,
					         lhs, cc, 0,
					         rhs, i, j, 0);

                    //---------------------------------------------------------------------
                    //     begin inner most do loop
                    //     do all the elements of the cell unless last 
                    //---------------------------------------------------------------------
                    for (k = 1; k <= ksize - 1; k++)
                    {

                        //---------------------------------------------------------------------
                        //     subtract A*lhs_vector(k-1) from lhs_vector(k)
                        //     
                        //     rhs(k) = rhs(k) - A*rhs(k-1)
                        //---------------------------------------------------------------------
                        matvec_sub(lhs, aa, k,
                                   rhs, i, j, (k - 1),
						           rhs, i, j, k);

                        //---------------------------------------------------------------------
                        //     B(k) = B(k) - C(k-1)*A(k)
                        //     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k)
                        //---------------------------------------------------------------------
                        matmul_sub(lhs, aa, k,
						           lhs, cc, (k - 1),
						           lhs, bb, k);

                        //---------------------------------------------------------------------
                        //     multiply c(i,j,k) by b_inverse and copy back to c
                        //     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
                        //---------------------------------------------------------------------
                        binvcrhs(lhs, bb, k,
						         lhs, cc, k,
						         rhs, i, j, k);
                    }

                    //---------------------------------------------------------------------
                    //     Now finish up special cases for last cell
                    //---------------------------------------------------------------------

                    //---------------------------------------------------------------------
                    //     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
                    //---------------------------------------------------------------------
                    matvec_sub(lhs, aa, ksize,
					           rhs, i, j, (ksize - 1),
					           rhs, i, j, ksize);

                    //---------------------------------------------------------------------
                    //     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
                    //     call matmul_sub(aa,i,j,ksize,c,
                    //     $              cc,i,j,ksize-1,c,bb,i,j,ksize)
                    //---------------------------------------------------------------------
                    matmul_sub(lhs, aa, ksize,
					           lhs, cc, (ksize - 1),
					           lhs, bb, ksize);

                    //---------------------------------------------------------------------
                    //     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
                    //---------------------------------------------------------------------
                    binvrhs(lhs, bb, ksize,
					        rhs, i, j, ksize);

                    //---------------------------------------------------------------------
                    //     back solve: if last cell, then generate U(ksize)=rhs(ksize)
                    //     else assume U(ksize) is loaded in un pack backsub_info
                    //     so just use it
                    //     after call u(kstart) will be sent to next cell
                    //---------------------------------------------------------------------

                    for (k = ksize - 1; k >= 0; k--)
                    {
                        for (m = 0; m <= BLOCK_SIZE - 1; m++)
                        {
                            for (n = 0; n <= BLOCK_SIZE - 1; n++)
                            {
                                rhs[k][j][i][m] += -lhs[k][cc][n][m]* rhs[k+1][j][i][n];
                            }
                        }
                    }
                }
            }
            if (timeron) timer.stop(t_zsolve);
        }
        public double checkSum(double[] arr)
        {
            double csum = 0.0;
            for (int k = 0; k <= grid_points[2] - 1; k++)
            {
                for (int j = 0; j <= grid_points[1] - 1; j++)
                {
                    for (int i = 0; i <= grid_points[0] - 1; i++)
                    {
                        for (int m = 0; m <= 4; m++)
                        {
                            int offset = m + i * isize2 + j * jsize2 + k * ksize2;
                            csum += (arr[offset] * arr[offset]) /
                                 (double)(grid_points[2] * grid_points[1] * grid_points[0] * 5);
                        }
                    }
                }
            }
            return csum;
        }

        public double getTime() { return timer.readTimer(1); }
        public void finalize() /*throws Throwable */ {
            Console.WriteLine("BT: is about to be garbage collected");
            //super.finalize();
        }
    }
}
