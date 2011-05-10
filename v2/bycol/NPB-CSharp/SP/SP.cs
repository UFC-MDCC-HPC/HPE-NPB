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

namespace NPB3_0_JAV {

public class SP : SPBase
{
	public int bid = -1;
	public BMResults results;
	public bool serial = false;
	public SP(char clss, int np, bool ser) : base(clss,np)
	{
		serial = ser;
	}
	public static void Main(String[] argv)
	{
		SP sp = null;
        Console.WriteLine("testw");
		BMArgs.ParseCmdLineArgs(argv, BMName);
		char CLSS = BMArgs.CLASS;
		int np = BMArgs.num_threads;
		bool serial = BMArgs.serial;
		try
		{
			sp = new SP(CLSS, np, serial);
		}
		catch (OutOfMemoryException e)
		{
			BMArgs.outOfMemoryMessage();
			Environment.Exit(0);
		}
		sp.runBenchMark();

	}

	public void run() { runBenchMark(); }

	public void runBenchMark()
	{
		BMArgs.Banner(BMName, CLASS, serial, num_threads);

		int numTimers = t_last + 1;
		String[] t_names = new String[numTimers];
		double[] trecs = new double[numTimers];
		setTimers(t_names);
		//---------------------------------------------------------------------
		//      Read input file (if it exists), else take
		//      defaults from parameters
		//---------------------------------------------------------------------
		int niter = getInputPars();
		set_constants(0);
		initialize();
		exact_rhs();

//		if (!serial) setupThreads(this);
		//---------------------------------------------------------------------
		//      do one time step to touch all code, and reinitialize
		//---------------------------------------------------------------------
		if (serial) adi_serial();
//		else adi();
		initialize();

		timer.resetAllTimers();
		timer.start(t_total);
		for (int step = 1; step <= niter; step++)
		{
			if (step % 20 == 0 || step == 1 || step == niter)
			{
				Console.WriteLine("Time step " + step);
			}
			if (serial) adi_serial();
//			else adi();
		}
		timer.stop(1);
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
	}

	public double getMFLOPS(double total_time, int niter)
	{
		double mflops = 0.0;
		if (total_time > 0)
		{
			int n3 = grid_points[0] * grid_points[1] * grid_points[2];
			double t = (grid_points[0] + grid_points[1] + grid_points[2]) / 3.0;
			mflops = 881.174 * n3
					 - 4683.91 * t * t
					 + 11484.5 * t - 19272.4;
			mflops *= niter / (total_time * 1000000.0);
		}
		return mflops;
	}

	public void adi_serial()
	{
		if (timeron) timer.start(t_rhs);
		compute_rhs();
		if (timeron) timer.stop(t_rhs);
		if (timeron) timer.start(t_txinvr);
		txinvr();
		if (timeron) timer.stop(t_txinvr);
		x_solve();
		y_solve();
		z_solve();
		if (timeron) timer.start(t_add);
		add();
		if (timeron) timer.stop(t_add);
	}

/*	public void adi()
	{
		if (timeron) timer.start(t_rhs);
		doRHS();
		doRHS();

		if (timeron) timer.start(t_rhsx);
		doRHS();
		if (timeron) timer.stop(t_rhsx);

		if (timeron) timer.start(t_rhsy);
		doRHS();
		if (timeron) timer.stop(t_rhsy);

		if (timeron) timer.start(t_rhsz);
		doRHS();
		if (timeron) timer.stop(t_rhsz);

		doRHS();
		if (timeron) timer.stop(t_rhs);

		if (timeron) timer.start(t_txinvr);
		synchronized (this)
		{
			for (int m = 0; m < num_threads; m++)
				synchronized (txinverse[m])
				{
					txinverse[m].done = false;
					txinverse[m].notify();
				}
			for (int m = 0; m < num_threads; m++)
				while (!txinverse[m].done)
				{
					try { wait(); }
					catch (InterruptedException e) { }
					notifyAll();
				}
		}
		if (timeron) timer.stop(t_txinvr);

		if (timeron) timer.start(t_xsolve);
		doXsolve();
		if (timeron) timer.stop(t_xsolve);

		if (timeron) timer.start(t_ninvr);
		doXsolve();
		if (timeron) timer.stop(t_ninvr);

		if (timeron) timer.start(t_ysolve);
		doYsolve();
		if (timeron) timer.stop(t_ysolve);

		if (timeron) timer.start(t_pinvr);
		doYsolve();
		if (timeron) timer.stop(t_pinvr);

		if (timeron) timer.start(t_zsolve);
		doZsolve();
		if (timeron) timer.stop(t_zsolve);

		if (timeron) timer.start(t_tzetar);
		doZsolve();
		if (timeron) timer.stop(t_tzetar);

		if (timeron) timer.start(t_add);
		synchronized (this)
		{
			for (int m = 0; m < num_threads; m++)
				synchronized (rhsadder[m])
				{
					rhsadder[m].done = false;
					rhsadder[m].notify();
				}
			for (int m = 0; m < num_threads; m++)
				while (!rhsadder[m].done)
				{
					try { wait(); }
					catch (InterruptedException e) { }
					notifyAll();
				}
		}
		if (timeron) timer.stop(t_add);
	}
  
    synchronized void doRHS()
	{
		int m;
		for (m = 0; m < num_threads; m++)
			synchronized (rhscomputer[m])
			{
				rhscomputer[m].done = false;
				rhscomputer[m].notify();
			}
		for (m = 0; m < num_threads; m++)
			while (!rhscomputer[m].done)
			{
				try { wait(); }
				catch (InterruptedException e) { }
				notifyAll();
			}
	}

	synchronized void doXsolve()
	{
		int m;
		for (m = 0; m < num_threads; m++)
			synchronized (xsolver[m])
			{
				xsolver[m].done = false;
				xsolver[m].notify();
			}
		for (m = 0; m < num_threads; m++)
			while (!xsolver[m].done)
			{
				try { wait(); }
				catch (InterruptedException e) { }
				notifyAll();
			}
	}

	synchronized void doYsolve()
	{
		int m;
		for (m = 0; m < num_threads; m++)
			synchronized (ysolver[m])
			{
				ysolver[m].done = false;
				ysolver[m].notify();
			}
		for (m = 0; m < num_threads; m++)
			while (!ysolver[m].done)
			{
				try { wait(); }
				catch (InterruptedException e) { }
				notifyAll();
			}
	}

	synchronized void doZsolve()
	{
		int m;
		for (m = 0; m < num_threads; m++)
			synchronized (zsolver[m])
			{
				zsolver[m].done = false;
				zsolver[m].notify();
			}
		for (m = 0; m < num_threads; m++)
			while (!zsolver[m].done)
			{
				try { wait(); }
				catch (InterruptedException e) { }
				notifyAll();
			}
	}
 */

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
				Console.Error.WriteLine("exception caught!");
			}
		}
		else
		{
			Console.WriteLine("No input file inputsp.data," +
							   "Using compiled defaults");
			niter = niter_default;
			dt = dt_default;
			grid_points[0] = problem_size;
			grid_points[1] = problem_size;
			grid_points[2] = problem_size;
		}
		if ((grid_points[0] > IMAX) ||
		 (grid_points[1] > JMAX) ||
		 (grid_points[2] > KMAX))
		{
			Console.WriteLine("Problem size too big for array");
			Environment.Exit(0);
		}
		Console.WriteLine("Iterations: " + niter + " dt: " + dt);

		nx2 = grid_points[0] - 2;
		ny2 = grid_points[1] - 2;
		nz2 = grid_points[2] - 2;
		return niter;
	}
	public void setTimers(String[] t_names)
	{
		timeron = false;
		if (File.Exists("timer.flag"))
		{
			timeron = true;
			t_names[t_total] = "total";
			t_names[t_rhsx] = "rhsx";
			t_names[t_rhsy] = "rhsy";
			t_names[t_rhsz] = "rhsz";
			t_names[t_rhs] = "rhs";
			t_names[t_xsolve] = "xsolve";
			t_names[t_ysolve] = "ysolve";
			t_names[t_zsolve] = "zsolve";
			t_names[t_rdis1] = "redist1";
			t_names[t_rdis2] = "redist2";
			t_names[t_tzetar] = "tzetar";
			t_names[t_ninvr] = "ninvr";
			t_names[t_pinvr] = "pinvr";
			t_names[t_txinvr] = "txinvr";
			t_names[t_add] = "add";
		}
	}
	public void printTimers(String[] t_names, double[] trecs, double tmax)
	{
		double t;
		Console.WriteLine("  SECTION   Time (secs)");
		for (int i = 1; i <= t_last; i++)
		{
			trecs[i] = timer.readTimer(i);
		}
		if (tmax == 0.0) tmax = 1.0;
		for (int i = 1; i < t_last; i++)
		{
            double dbl = trecs[i] * 100 / tmax;
			Console.WriteLine(t_names[i] + ":" + trecs[i].ToString("N3") +
						   "  (" + dbl.ToString("N3") + "%)");
			if (i == t_rhs)
			{
				t = trecs[t_rhsx] + trecs[t_rhsy] + trecs[t_rhsz];
                dbl = (t * 100.0 / tmax);
				Console.WriteLine("    --> total " + "sub-rhs" + ":" + (t.ToString("N3")) +
							   "  (" + dbl.ToString("N3") + "%)");
				t = trecs[t_rhs] - t;
                dbl = (t * 100.0 / tmax);
				Console.WriteLine("    --> total " + "rest-rhs" + ":" + (t.ToString("N3")) +
							   "  (" + dbl.ToString("N3") + "%)");
			}
			else if (i == t_zsolve)
			{
				t = trecs[t_zsolve] - trecs[t_rdis1] - trecs[t_rdis2];
                dbl = (t * 100.0 / tmax);
				Console.WriteLine("    --> total " + "sub-zsol" + ":" + (t.ToString("N3")) +
							   "  (" + dbl.ToString("N3") + "%)");
			}
			else if (i == t_rdis2)
			{
				t = trecs[t_rdis1] + trecs[t_rdis2];
                dbl = (t * 100.0 / tmax);
				Console.WriteLine("    --> total " + "redist" + ":" + (t.ToString("N3")) +
							   "  (" + dbl.ToString("N3") + "%)");
			}
		}
	}


	public void add()
	{
		int i, j, k, m;
		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						u[m][i][j][k] +=
								   rhs[m][i][j][k];
					}
				}
			}
		}
	}

	public void error_norm(double[] rms)
	{
		int i, j, k, m, d;
		double xi, eta, zeta;
        double[] u_exact = new double[5];
        double add;

		for (m = 0; m <= 4; m++)
		{
			rms[m] = 0.0;
		}

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

					for (m = 0; m <= 4; m++)
					{
						add = u[m][i][j][k] - u_exact[m];
						rms[m] = rms[m] + add * add;
					}
				}
			}
		}

        foreach (double x in rms) {

        }

		for (m = 0; m <= 4; m++)
		{
			for (d = 0; d <= 2; d++)
			{
				rms[m] = rms[m] / (grid_points[d] - 2);
			}
			rms[m] = Math.Sqrt(rms[m]);
		}
	}

	public void rhs_norm(double[] rms)
	{

		int i, j, k, d, m;
		double add;

		for (m = 0; m <= 4; m++)
		{
			rms[m] = 0.0;
		}

		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						add = rhs[m][i][j][k];
						rms[m] = rms[m] + add * add;
					}
				}
			}
		}
		for (m = 0; m <= 4; m++)
		{
			for (d = 0; d <= 2; d++)
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
						forcing[m][i][j][k] = 0.0;
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
						ue[m][i] = dtemp[m];
					}

					dtpp = 1.0 / dtemp[0];

					for (m = 1; m <= 4; m++)
					{
						buf[m][i] = dtpp * dtemp[m];
					}

					cuf[i] = buf[1][i] * buf[1][i];
					buf[0][i] = cuf[i] + buf[2][i] * buf[2][i] +
										 buf[3][i] * buf[3][i];
					q[i] = 0.5 * (buf[1][i] * ue[1][i] + buf[2][i] * ue[2][i] +
											buf[3][i] * ue[3][i]);

				}

				for (i = 1; i <= grid_points[0] - 2; i++)
				{
					im1 = i - 1;
					ip1 = i + 1;

					forcing[0][i][j][k] = forcing[0][i][j][k] -
									 tx2 * (ue[1][ip1] - ue[1][im1]) +
									 dx1tx1 * (ue[0][ip1] - 2.0 * ue[0][i] + ue[0][im1]);

					forcing[1][i][j][k] = forcing[1][i][j][k] - tx2 * (
									(ue[1][ip1] * buf[1][ip1] + c2 * (ue[4][ip1] - q[ip1])) -
									(ue[1][im1] * buf[1][im1] + c2 * (ue[4][im1] - q[im1]))) +
									 xxcon1 * (buf[1][ip1] - 2.0 * buf[1][i] + buf[1][im1]) +
									 dx2tx1 * (ue[1][ip1] - 2.0 * ue[1][i] + ue[1][im1]);

					forcing[2][i][j][k] = forcing[2][i][j][k] - tx2 * (
									 ue[2][ip1] * buf[1][ip1] - ue[2][im1] * buf[1][im1]) +
									 xxcon2 * (buf[2][ip1] - 2.0 * buf[2][i] + buf[2][im1]) +
									 dx3tx1 * (ue[2][ip1] - 2.0 * ue[2][i] + ue[2][im1]);


					forcing[3][i][j][k] = forcing[3][i][j][k] - tx2 * (
									 ue[3][ip1] * buf[1][ip1] - ue[3][im1] * buf[1][im1]) +
									 xxcon2 * (buf[3][ip1] - 2.0 * buf[3][i] + buf[3][im1]) +
									 dx4tx1 * (ue[3][ip1] - 2.0 * ue[3][i] + ue[3][im1]);

					forcing[4][i][j][k] = forcing[4][i][j][k] - tx2 * (
									 buf[1][ip1] * (c1 * ue[4][ip1] - c2 * q[ip1]) -
									 buf[1][im1] * (c1 * ue[4][im1] - c2 * q[im1])) +
									 0.5 * xxcon3 * (buf[0][ip1] - 2.0 * buf[0][i] +
												   buf[0][im1]) +
									 xxcon4 * (cuf[ip1] - 2.0 * cuf[i] + cuf[im1]) +
									 xxcon5 * (buf[4][ip1] - 2.0 * buf[4][i] + buf[4][im1]) +
									 dx5tx1 * (ue[4][ip1] - 2.0 * ue[4][i] + ue[4][im1]);
                }

				//---------------------------------------------------------------------
				//            Fourth-order dissipation                         
				//---------------------------------------------------------------------
				for (m = 0; m <= 4; m++)
				{
					i = 1;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
										(5.0 * ue[m][i] - 4.0 * ue[m][i+1] + ue[m][i+2]);
					i = 2;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
									   (-4.0 * ue[m][i-1] + 6.0 * ue[m][i] -
										 4.0 * ue[m][i+1] + ue[m][i+2]);
				}

				for (m = 0; m <= 4; m++)
				{
					for (i = 3; i <= grid_points[0] - 4; i++)
					{
						forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
										 (ue[m][i-2] - 4.0 * ue[m][i-1] +
										  6.0 * ue[m][i] - 4.0 * ue[m][i+1] + ue[m][i+2]);
					}
				}

				for (m = 0; m <= 4; m++)
				{
					i = grid_points[0] - 3;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
									   (ue[m][i-2] - 4.0 * ue[m][i-1] +
										6.0 * ue[m][i] - 4.0 * ue[m][i+1]);
					i = grid_points[0] - 2;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
									   (ue[m][i-2] - 4.0 * ue[m][i-1] + 5.0 * ue[m][i]);
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
						ue[m][j] = dtemp[m];
					}
					dtpp = 1.0 / dtemp[0];

					for (m = 1; m <= 4; m++)
					{
						buf[m][j] = dtpp * dtemp[m];
					}

					cuf[j] = buf[2][j] * buf[2][j];
					buf[0][j] = cuf[j] + buf[1][j] * buf[1][j] +
							   buf[3][j] * buf[3][j];
					q[j] = 0.5 * (buf[1][j] * ue[1][j] + buf[2][j] * ue[2][j] +
								  buf[3][j] * ue[3][j]);
				}

				for (j = 1; j <= grid_points[1] - 2; j++)
				{
					jm1 = j - 1;
					jp1 = j + 1;

					forcing[0][i][j][k] = forcing[0][i][j][k] -
						  ty2 * (ue[2][jp1] - ue[2][jm1]) +
						  dy1ty1 * (ue[0][jp1] - 2.0 * ue[0][j] + ue[0][jm1]);

					forcing[1][i][j][k] = forcing[1][i][j][k] - ty2 * (
						  ue[1][jp1] * buf[2][jp1] - ue[1][jm1] * buf[2][jm1]) +
						  yycon2 * (buf[1][jp1] - 2.0 * buf[1][j] + buf[1][jm1]) +
						  dy2ty1 * (ue[1][jp1] - 2.0 * ue[1][j] + ue[1][jm1]);

					forcing[2][i][j][k] = forcing[2][i][j][k] - ty2 * (
						  (ue[2][jp1] * buf[2][jp1] + c2 * (ue[4][jp1] - q[jp1])) -
						  (ue[2][jm1] * buf[2][jm1] + c2 * (ue[4][jm1] - q[jm1]))) +
						  yycon1 * (buf[2][jp1] - 2.0 * buf[2][j] + buf[2][jm1]) +
						  dy3ty1 * (ue[2][jp1] - 2.0 * ue[2][j] + ue[2][jm1]);

					forcing[3][i][j][k] = forcing[3][i][j][k] - ty2 * (
						  ue[3][jp1] * buf[2][jp1] - ue[3][jm1] * buf[2][jm1]) +
						  yycon2 * (buf[3][jp1] - 2.0 * buf[3][j] + buf[3][jm1]) +
						  dy4ty1 * (ue[3][jp1] - 2.0 * ue[3][j] + ue[3][jm1]);

					forcing[4][i][j][k] = forcing[4][i][j][k] - ty2 * (
						  buf[2][jp1] * (c1 * ue[4][jp1] - c2 * q[jp1]) -
						  buf[2][jm1] * (c1 * ue[4][jm1] - c2 * q[jm1])) +
						  0.5 * yycon3 * (buf[0][jp1] - 2.0 * buf[0][j] +
										buf[0][jm1]) +
						  yycon4 * (cuf[jp1] - 2.0 * cuf[j] + cuf[jm1]) +
						  yycon5 * (buf[4][jp1] - 2.0 * buf[4][j] + buf[4][jm1]) +
						  dy5ty1 * (ue[4][jp1] - 2.0 * ue[4][j] + ue[4][jm1]);
				}

				//---------------------------------------------------------------------
				//            Fourth-order dissipation                      
				//---------------------------------------------------------------------
				for (m = 0; m <= 4; m++)
				{
					j = 1;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							  (5.0 * ue[m][j] - 4.0 * ue[m][j+1] + ue[m][j+2]);
					j = 2;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							 (-4.0 * ue[m][j-1] + 6.0 * ue[m][j] -
							   4.0 * ue[m][j+1] + ue[m][j+2]);
				}

				for (m = 0; m <= 4; m++)
				{
					for (j = 3; j <= grid_points[1] - 4; j++)
					{
						forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							  (ue[m][j-2] - 4.0 * ue[m][j-1] +
							   6.0 * ue[m][j] - 4.0 * ue[m][j+1] + ue[m][j+2]);
					}
				}

				for (m = 0; m <= 4; m++)
				{
					j = grid_points[1] - 3;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							 (ue[m][j-2] - 4.0 * ue[m][j-1] +
							  6.0 * ue[m][j] - 4.0 * ue[m][j+1]);
					j = grid_points[1] - 2;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							 (ue[m][j-2] - 4.0 * ue[m][j-1] + 5.0 * ue[m][j]);

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
						ue[m][k] = dtemp[m];
					}

					dtpp = 1.0 / dtemp[0];

					for (m = 1; m <= 4; m++)
					{
						buf[m][k] = dtpp * dtemp[m];
					}

					cuf[k] = buf[3][k] * buf[3][k];
					buf[0][k] = cuf[k] + buf[1][k] * buf[1][k] +
							   buf[2][k] * buf[2][k];
					q[k] = 0.5 * (buf[1][k] * ue[1][k] + buf[2][k] * ue[2][k] +
								  buf[3][k] * ue[3][k]);
				}

				for (k = 1; k <= grid_points[2] - 2; k++)
				{
					km1 = k - 1;
					kp1 = k + 1;

					forcing[0][i][j][k] = forcing[0][i][j][k] -
						   tz2 * (ue[3][kp1] - ue[3][km1]) +
						   dz1tz1 * (ue[0][kp1] - 2.0 * ue[0][k] + ue[0][km1]);

					forcing[1][i][j][k] = forcing[1][i][j][k] - tz2 * (
						   ue[1][kp1] * buf[3][kp1] - ue[1][km1] * buf[3][km1]) +
						   zzcon2 * (buf[1][kp1] - 2.0 * buf[1][k] + buf[1][km1]) +
						   dz2tz1 * (ue[1][kp1] - 2.0 * ue[1][k] + ue[1][km1]);

					forcing[2][i][j][k] = forcing[2][i][j][k] - tz2 * (
						   ue[2][kp1] * buf[3][kp1] - ue[2][km1] * buf[3][km1]) +
						   zzcon2 * (buf[2][kp1] - 2.0 * buf[2][k] + buf[2][km1]) +
						   dz3tz1 * (ue[2][kp1] - 2.0 * ue[2][k] + ue[2][km1]);

					forcing[3][i][j][k] = forcing[3][i][j][k] - tz2 * (
						  (ue[3][kp1] * buf[3][kp1] + c2 * (ue[4][kp1] - q[kp1])) -
						  (ue[3][km1] * buf[3][km1] + c2 * (ue[4][km1] - q[km1]))) +
						  zzcon1 * (buf[3][kp1] - 2.0 * buf[3][k] + buf[3][km1]) +
						  dz4tz1 * (ue[3][kp1] - 2.0 * ue[3][k] + ue[3][km1]);

					forcing[4][i][j][k] = forcing[4][i][j][k] - tz2 * (
						   buf[3][kp1] * (c1 * ue[4][kp1] - c2 * q[kp1]) -
						   buf[3][km1] * (c1 * ue[4][km1] - c2 * q[km1])) +
						   0.5 * zzcon3 * (buf[0][kp1] - 2.0 * buf[0][k]
										+ buf[0][km1]) +
						   zzcon4 * (cuf[kp1] - 2.0 * cuf[k] + cuf[km1]) +
						   zzcon5 * (buf[4][kp1] - 2.0 * buf[4][k] + buf[4][km1]) +
						   dz5tz1 * (ue[4][kp1] - 2.0 * ue[4][k] + ue[4][km1]);
				}

				//---------------------------------------------------------------------
				//            Fourth-order dissipation
				//---------------------------------------------------------------------
				for (m = 0; m <= 4; m++)
				{
					k = 1;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							  (5.0 * ue[m][k] - 4.0 * ue[m][k+1] + ue[m][k+2]);
					k = 2;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							 (-4.0 * ue[m][k-1] + 6.0 * ue[m][k] -
							   4.0 * ue[m][k+1] + ue[m][k+2]);
				}

				for (m = 0; m <= 4; m++)
				{
					for (k = 3; k <= grid_points[2] - 4; k++)
					{
						forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							  (ue[m][k-2] - 4.0 * ue[m][k-1] +
							   6.0 * ue[m][k] - 4.0 * ue[m][k+1] + ue[m][k+2]);
					}
				}

				for (m = 0; m <= 4; m++)
				{
					k = grid_points[2] - 3;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
							 (ue[m][k-2] - 4.0 * ue[m][k-1] +
							  6.0 * ue[m][k] - 4.0 * ue[m][k+1]);
					k = grid_points[2] - 2;
					forcing[m][i][j][k] = forcing[m][i][j][k] - dssp *
						  (ue[m][k-2] - 4.0 * ue[m][k-1] + 5.0 * ue[m][k]);
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
						forcing[m][i][j][k] = -1.0 * forcing[m][i][j][k];
					}
				}
			}
		}
	}

	public void ninvr()
	{
		int i, j, k;
		double r1, r2, r3, r4, r5, t1, t2;

		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{

					r1 = rhs[0][i][j][k];
					r2 = rhs[1][i][j][k];
					r3 = rhs[2][i][j][k];
					r4 = rhs[3][i][j][k];
					r5 = rhs[4][i][j][k];

					t1 = bt * r3;
					t2 = 0.5 * (r4 + r5);

					rhs[0][i][j][k] = -r2;
					rhs[1][i][j][k] = r1;
					rhs[2][i][j][k] = bt * (r4 - r5);
					rhs[3][i][j][k] = -t1 + t2;
					rhs[4][i][j][k] = t1 + t2;
				}
			}
		}
	}

	public void pinvr()
	{
		int i, j, k;
		double r1, r2, r3, r4, r5, t1, t2;

		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{

					r1 = rhs[0][i][j][k];
					r2 = rhs[1][i][j][k];
					r3 = rhs[2][i][j][k];
					r4 = rhs[3][i][j][k];
					r5 = rhs[4][i][j][k];

					t1 = bt * r1;
					t2 = 0.5 * (r4 + r5);

					rhs[0][i][j][k] = bt * (r4 - r5);
					rhs[1][i][j][k] = -r3;
					rhs[2][i][j][k] = r2;
					rhs[3][i][j][k] = -t1 + t2;
					rhs[4][i][j][k] = t1 + t2;
				}
			}
		}
	}

	public void compute_rhs()
	{
		int i, j, k, m;
		double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1,
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
					rho_inv = 1.0 / u[0][i][j][k];
					rho_i[i][j][k] = rho_inv;
					us[i][j][k] = u[1][i][j][k] * rho_inv;
					vs[i][j][k] = u[2][i][j][k] * rho_inv;
					ws[i][j][k] = u[3][i][j][k] * rho_inv;
					square[i][j][k] = 0.5 * (
								  u[1][i][j][k] * u[1][i][j][k] +
								  u[2][i][j][k] * u[2][i][j][k] +
								  u[3][i][j][k] * u[3][i][j][k]) * rho_inv;
					qs[i][j][k] = square[i][j][k] * rho_inv;
					//---------------------------------------------------------------------
					//               (don't need speed and ainx until the lhs computation)
					//---------------------------------------------------------------------
					aux = c1c2 * rho_inv * (u[4][i][j][k] - square[i][j][k]);
					speed[i][j][k] = Math.Sqrt(aux);
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
						rhs[m][i][j][k] = forcing[m][i][j][k];
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
					uijk = us[i][j][k];
					up1 = us[i+1][j][k];
					um1 = us[i-1][j][k];

					rhs[0][i][j][k] = rhs[0][i][j][k] + dx1tx1 *
							  (u[0][i+1][j][k] - 2.0 * u[0][i][j][k] +
							   u[0][i-1][j][k]) -
							  tx2 * (u[1][i+1][j][k] - u[1][i-1][j][k]);

					rhs[1][i][j][k] = rhs[1][i][j][k] + dx2tx1 *
							  (u[1][i+1][j][k] - 2.0 * u[1][i][j][k] +
							   u[1][i-1][j][k]) +
							  xxcon2 * con43 * (up1 - 2.0 * uijk + um1) -
							  tx2 * (u[1][i+1][j][k] * up1 -
									 u[1][i-1][j][k] * um1 +
									 (u[4][i+1][j][k] - square[i+1][j][k] -
									  u[4][i-1][j][k] + square[i-1][j][k]) *
									  c2);

					rhs[2][i][j][k] = rhs[2][i][j][k] + dx3tx1 *
							  (u[2][i+1][j][k] - 2.0 * u[2][i][j][k] +
							   u[2][i-1][j][k]) +
							  xxcon2 * (vs[i+1][j][k] - 2.0 * vs[i][j][k] +
										vs[i-1][j][k]) -
							  tx2 * (u[2][i+1][j][k] * up1 -
									 u[2][i-1][j][k] * um1);

					rhs[3][i][j][k] = rhs[3][i][j][k] + dx4tx1 *
							  (u[3][i+1][j][k] - 2.0 * u[3][i][j][k] +
							   u[3][i-1][j][k]) +
							  xxcon2 * (ws[i+1][j][k] - 2.0 * ws[i][j][k] +
										ws[i-1][j][k]) -
							  tx2 * (u[3][i+1][j][k] * up1 -
									 u[3][i-1][j][k] * um1);

					rhs[4][i][j][k] = rhs[4][i][j][k] + dx5tx1 *
							  (u[4][i+1][j][k] - 2.0 * u[4][i][j][k] +
							   u[4][i-1][j][k]) +
							  xxcon3 * (qs[i+1][j][k] - 2.0 * qs[i][j][k] +
										qs[i-1][j][k]) +
							  xxcon4 * (up1 * up1 - 2.0 * uijk * uijk +
										um1 * um1) +
							  xxcon5 * (u[4][i+1][j][k] * rho_i[i+1][j][k] -
										2.0 * u[4][i][j][k] * rho_i[i][j][k] +
										u[4][i-1][j][k] * rho_i[i-1][j][k]) -
							  tx2 * ((c1 * u[4][i+1][j][k] -
									   c2 * square[i+1][j][k]) * up1 -
									  (c1 * u[4][i-1][j][k] -
									   c2 * square[i-1][j][k]) * um1);
				}

				//---------------------------------------------------------------------
				//      add fourth order xi-direction dissipation               
				//---------------------------------------------------------------------

				i = 1;
				for (m = 0; m <= 4; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (5.0 * u[m][i][j][k] - 4.0 * u[m][i+1][j][k] +
									  u[m][i+2][j][k]);
				}

				i = 2;
				for (m = 0; m <= 4; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (-4.0 * u[m][i-1][j][k] + 6.0 * u[m][i][j][k] -
								4.0 * u[m][i+1][j][k] + u[m][i+2][j][k]);
				}

				for (i = 3; i <= nx2 - 2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							   (u[m][i-2][j][k] - 4.0 * u[m][i-1][j][k] +
								6.0 * u[m][i][j][k] - 4.0 * u[m][i+1][j][k] +
									u[m][i+2][j][k]);
					}
				}

				i = nx2 - 1;
				for (m = 0; m <= 4; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (u[m][i-2][j][k] - 4.0 * u[m][i-1][j][k] +
								6.0 * u[m][i][j][k] - 4.0 * u[m][i+1][j][k]);
				}

				i = nx2;
				for (m = 0; m <= 4; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (u[m][i-2][j][k] - 4.0 * u[m][i-1][j][k] +
								5.0 * u[m][i][j][k]);
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
					vijk = vs[i][j][k];
					vp1 = vs[i][j+1][k];
					vm1 = vs[i][j-1][k];
					rhs[0][i][j][k] = rhs[0][i][j][k] + dy1ty1 *
							 (u[0][i][j+1][k] - 2.0 * u[0][i][j][k] +
							  u[0][i][j-1][k]) -
							 ty2 * (u[2][i][j+1][k] - u[2][i][j-1][k]);
					rhs[1][i][j][k] = rhs[1][i][j][k] + dy2ty1 *
							 (u[1][i][j+1][k] - 2.0 * u[1][i][j][k] +
							  u[1][i][j-1][k]) +
							 yycon2 * (us[i][j+1][k] - 2.0 * us[i][j][k] +
									   us[i][j-1][k]) -
							 ty2 * (u[1][i][j+1][k] * vp1 -
									u[1][i][j-1][k] * vm1);
					rhs[2][i][j][k] = rhs[2][i][j][k] + dy3ty1 *
							 (u[2][i][j+1][k] - 2.0 * u[2][i][j][k] +
							  u[2][i][j-1][k]) +
							 yycon2 * con43 * (vp1 - 2.0 * vijk + vm1) -
							 ty2 * (u[2][i][j+1][k] * vp1 -
									u[2][i][j-1][k] * vm1 +
									(u[4][i][j+1][k] - square[i][j+1][k] -
									 u[4][i][j-1][k] + square[i][j-1][k])
									* c2);
					rhs[3][i][j][k] = rhs[3][i][j][k] + dy4ty1 *
							 (u[3][i][j+1][k] - 2.0 * u[3][i][j][k] +
							  u[3][i][j-1][k]) +
							 yycon2 * (ws[i][j+1][k] - 2.0 * ws[i][j][k] +
									   ws[i][j-1][k]) -
							 ty2 * (u[3][i][j+1][k] * vp1 -
									u[3][i][j-1][k] * vm1);
					rhs[4][i][j][k] = rhs[4][i][j][k] + dy5ty1 *
							 (u[4][i][j+1][k] - 2.0 * u[4][i][j][k] +
							  u[4][i][j-1][k]) +
							 yycon3 * (qs[i][j+1][k] - 2.0 * qs[i][j][k] +
									   qs[i][j-1][k]) +
							 yycon4 * (vp1 * vp1 - 2.0 * vijk * vijk +
									   vm1 * vm1) +
							 yycon5 * (u[4][i][j+1][k] * rho_i[i][j+1][k] -
									   2.0 * u[4][i][j][k] * rho_i[i][j][k] +
									   u[4][i][j-1][k] * rho_i[i][j-1][k]) -
							 ty2 * ((c1 * u[4][i][j+1][k] -
									 c2 * square[i][j+1][k]) * vp1 -
									(c1 * u[4][i][j-1][k] -
									 c2 * square[i][j-1][k]) * vm1);
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
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (5.0 * u[m][i][j][k] - 4.0 * u[m][i][j+1][k] +
									  u[m][i][j+2][k]);
				}
			}

			j = 2;
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (-4.0 * u[m][i][j-1][k] + 6.0 * u[m][i][j][k] -
								4.0 * u[m][i][j+1][k] + u[m][i][j+2][k]);
				}
			}

			for (j = 3; j <= ny2 - 2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{
					for (m = 0; m <= 4; m++)
					{
						rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							   (u[m][i][j-2][k] - 4.0 * u[m][i][j-1][k] +
								6.0 * u[m][i][j][k] - 4.0 * u[m][i][j+1][k] +
									u[m][i][j+2][k]);
					}
				}
			}

			j = ny2 - 1;
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (u[m][i][j-2][k] - 4.0 * u[m][i][j-1][k] +
								6.0 * u[m][i][j][k] - 4.0 * u[m][i][j+1][k]);
				}
			}

			j = ny2;
			for (i = 1; i <= nx2; i++)
			{
				for (m = 0; m <= 4; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (u[m][i][j-2][k] - 4.0 * u[m][i][j-1][k] +
								5.0 * u[m][i][j][k]);
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
					wijk = ws[i][j][k];
					wp1 = ws[i][j][k+1];
					wm1 = ws[i][j][k-1];

					rhs[0][i][j][k] = rhs[0][i][j][k] + dz1tz1 *
							 (u[0][i][j][k+1] - 2.0 * u[0][i][j][k] +
							  u[0][i][j][k-1]) -
							 tz2 * (u[3][i][j][k+1] - u[3][i][j][k-1]);
					rhs[1][i][j][k] = rhs[1][i][j][k] + dz2tz1 *
							 (u[1][i][j][k+1] - 2.0 * u[1][i][j][k] +
							  u[1][i][j][k-1]) +
							 zzcon2 * (us[i][j][k+1] - 2.0 * us[i][j][k] +
									   us[i][j][k-1]) -
							 tz2 * (u[1][i][j][k+1] * wp1 -
									u[1][i][j][k-1] * wm1);
					rhs[2][i][j][k] = rhs[2][i][j][k] + dz3tz1 *
							 (u[2][i][j][k+1] - 2.0 * u[2][i][j][k] +
							  u[2][i][j][k-1]) +
							 zzcon2 * (vs[i][j][k+1] - 2.0 * vs[i][j][k] +
									   vs[i][j][k-1]) -
							 tz2 * (u[2][i][j][k+1] * wp1 -
									u[2][i][j][k-1] * wm1);
					rhs[3][i][j][k] = rhs[3][i][j][k] + dz4tz1 *
							 (u[3][i][j][k+1] - 2.0 * u[3][i][j][k] +
							  u[3][i][j][k-1]) +
							 zzcon2 * con43 * (wp1 - 2.0 * wijk + wm1) -
							 tz2 * (u[3][i][j][k+1] * wp1 -
									u[3][i][j][k-1] * wm1 +
									(u[4][i][j][k+1] - square[i][j][k+1] -
									 u[4][i][j][k-1] + square[i][j][k-1])
									* c2);
					rhs[4][i][j][k] = rhs[4][i][j][k] + dz5tz1 *
							 (u[4][i][j][k+1] - 2.0 * u[4][i][j][k] +
							  u[4][i][j][k-1]) +
							 zzcon3 * (qs[i][j][k+1] - 2.0 * qs[i][j][k] +
									   qs[i][j][k-1]) +
							 zzcon4 * (wp1 * wp1 - 2.0 * wijk * wijk +
									   wm1 * wm1) +
							 zzcon5 * (u[4][i][j][k+1] * rho_i[i][j][k+1] -
									   2.0 * u[4][i][j][k] * rho_i[i][j][k] +
									   u[4][i][j][k-1] * rho_i[i][j][k-1]) -
							 tz2 * ((c1 * u[4][i][j][k+1] -
									  c2 * square[i][j][k+1]) * wp1 -
									 (c1 * u[4][i][j][k-1] -
									  c2 * square[i][j][k-1]) * wm1);
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
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (5.0 * u[m][i][j][k] - 4.0 * u[m][i][j][k+1] +
									  u[m][i][j][k+2]);
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
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (-4.0 * u[m][i][j][k-1] + 6.0 * u[m][i][j][k] -
								4.0 * u[m][i][j][k+1] + u[m][i][j][k+2]);
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
						rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							   (u[m][i][j][k-2] - 4.0 * u[m][i][j][k-1] +
								6.0 * u[m][i][j][k] - 4.0 * u[m][i][j][k+1] +
									u[m][i][j][k+2]);
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
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (u[m][i][j][k-2] - 4.0 * u[m][i][j][k-1] +
								6.0 * u[m][i][j][k] - 4.0 * u[m][i][j][k+1]);
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
					rhs[m][i][j][k] = rhs[m][i][j][k] - dssp *
							  (u[m][i][j][k-2] - 4.0 * u[m][i][j][k-1] +
								5.0 * u[m][i][j][k]);
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
                        rhs[m][i][j][k] = rhs[m][i][j][k] * dt;

                    }
                }
            }
        }
		
	}
    
	public void txinvr()
	{
		int i, j, k;
		double t1, t2, t3, ac, ru1, uu, vv, ww,
			   r1, r2, r3, r4, r5, ac2inv;

		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{

					ru1 = rho_i[i][j][k];
					uu = us[i][j][k];
					vv = vs[i][j][k];
					ww = ws[i][j][k];
					ac = speed[i][j][k];
					ac2inv = 1.0 / (ac * ac);

					r1 = rhs[0][i][j][k];
					r2 = rhs[1][i][j][k];
					r3 = rhs[2][i][j][k];
					r4 = rhs[3][i][j][k];
					r5 = rhs[4][i][j][k];

					t1 = c2 * ac2inv * (qs[i][j][k] * r1 - uu * r2 -
						vv * r3 - ww * r4 + r5);
					t2 = bt * ru1 * (uu * r1 - r2);
					t3 = (bt * ru1 * ac) * t1;

					rhs[0][i][j][k] = r1 - t1;
					rhs[1][i][j][k] = -ru1 * (ww * r1 - r4);
					rhs[2][i][j][k] = ru1 * (vv * r1 - r3);
					rhs[3][i][j][k] = -t2 + t3;
					rhs[4][i][j][k] = t2 + t3;

				}
			}
		}
	}

	public void tzetar()
	{
		int i, j, k;
		double t1, t2, t3, ac, xvel, yvel, zvel,
				r1, r2, r3, r4, r5, btuz, acinv, ac2u, uzik1;

		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{
				for (i = 1; i <= nx2; i++)
				{

					xvel = us[i][j][k];
					yvel = vs[i][j][k];
					zvel = ws[i][j][k];
					ac = speed[i][j][k];

					ac2u = ac * ac;

					r1 = rhs[0][i][j][k];
					r2 = rhs[1][i][j][k];
					r3 = rhs[2][i][j][k];
					r4 = rhs[3][i][j][k];
					r5 = rhs[4][i][j][k];

					uzik1 = u[0][i][j][k];
					btuz = bt * uzik1;

					t1 = btuz / ac * (r4 + r5);
					t2 = r3 + t1;
					t3 = btuz * (r4 - r5);

					rhs[0][i][j][k] = t2;
					rhs[1][i][j][k] = -uzik1 * r2 + xvel * t2;
					rhs[2][i][j][k] = uzik1 * r1 + yvel * t2;
					rhs[3][i][j][k] = zvel * t2 + t3;
					rhs[4][i][j][k] = uzik1 * (-xvel * r2 + yvel * r1) +
						  qs[i][j][k] * t2 + c2iv * ac2u * t1 + zvel * t3;

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
		compute_rhs();
		rhs_norm(xcr);

		for (m = 0; m <= 4; m++) xcr[m] = xcr[m] / dt;

		for (m = 0; m <= 4; m++)
		{
			xcrref[m] = 1.0;
			xceref[m] = 1.0;
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
			dtref = .015;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of residual.
			//---------------------------------------------------------------------
			xcrref[0] = 2.7470315451339479E-2;
			xcrref[1] = 1.0360746705285417E-2;
			xcrref[2] = 1.6235745065095532E-2;
			xcrref[3] = 1.5840557224455615E-2;
			xcrref[4] = 3.4849040609362460E-2;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of solution error.
			//---------------------------------------------------------------------
			xceref[0] = 2.7289258557377227E-5;
			xceref[1] = 1.0364446640837285E-5;
			xceref[2] = 1.6154798287166471E-5;
			xceref[3] = 1.5750704994480102E-5;
			xceref[4] = 3.4177666183390531E-5;


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
			xcrref[0] = 0.1893253733584E-2;
			xcrref[1] = 0.1717075447775E-3;
			xcrref[2] = 0.2778153350936E-3;
			xcrref[3] = 0.2887475409984E-3;
			xcrref[4] = 0.3143611161242E-2;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of solution error.
			//---------------------------------------------------------------------
			xceref[0] = 0.7542088599534E-4;
			xceref[1] = 0.6512852253086E-5;
			xceref[2] = 0.1049092285688E-4;
			xceref[3] = 0.1128838671535E-4;
			xceref[4] = 0.1212845639773E-3;

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
			dtref = .0015;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of residual.
			//---------------------------------------------------------------------
			xcrref[0] = 2.4799822399300195;
			xcrref[1] = 1.1276337964368832;
			xcrref[2] = 1.5028977888770491;
			xcrref[3] = 1.4217816211695179;
			xcrref[4] = 2.1292113035138280;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of solution error.
			//---------------------------------------------------------------------
			xceref[0] = 1.0900140297820550E-4;
			xceref[1] = 3.7343951769282091E-5;
			xceref[2] = 5.0092785406541633E-5;
			xceref[3] = 4.7671093939528255E-5;
			xceref[4] = 1.3621613399213001E-4;

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
			dtref = .001;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of residual.
			//---------------------------------------------------------------------
			xcrref[0] = 0.6903293579998E+02;
			xcrref[1] = 0.3095134488084E+02;
			xcrref[2] = 0.4103336647017E+02;
			xcrref[3] = 0.3864769009604E+02;
			xcrref[4] = 0.5643482272596E+02;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of solution error.
			//---------------------------------------------------------------------
			xceref[0] = 0.9810006190188E-02;
			xceref[1] = 0.1022827905670E-02;
			xceref[2] = 0.1720597911692E-02;
			xceref[3] = 0.1694479428231E-02;
			xceref[4] = 0.1847456263981E-01;

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
			dtref = .00067;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of residual.
			//---------------------------------------------------------------------
			xcrref[0] = 0.5881691581829E+03;
			xcrref[1] = 0.2454417603569E+03;
			xcrref[2] = 0.3293829191851E+03;
			xcrref[3] = 0.3081924971891E+03;
			xcrref[4] = 0.4597223799176E+03;

			//---------------------------------------------------------------------
			//    Reference values of RMS-norms of solution error.
			//---------------------------------------------------------------------
			xceref[0] = 0.2598120500183;
			xceref[1] = 0.2590888922315E-01;
			xceref[2] = 0.5132886416320E-01;
			xceref[3] = 0.4806073419454E-01;
			xceref[4] = 0.5483377491301;
		}
		//---------------------------------------------------------------------
		//    verification test for residuals if gridsize is either 12X12X12 or 
		//    64X64X64 or 102X102X102 or 162X162X162
		//---------------------------------------------------------------------

		//---------------------------------------------------------------------
		//    Compute the difference of solution values and the known reference values.
		//---------------------------------------------------------------------
		for (m = 0; m <= 4; m++)
		{
			xcrdif[m] = Math.Abs((xcr[m] - xcrref[m]) / xcrref[m]);
			xcedif[m] = Math.Abs((xce[m] - xceref[m]) / xceref[m]);
		}
		//---------------------------------------------------------------------
		//   tolerance level
		//---------------------------------------------------------------------
		double epsilon = 1.0E-8;
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

	public void x_solve()
	{
		int i, j, k, n, i1, i2, m;
		double ru1, fac1, fac2;

		//---------------------------------------------------------------------
		//---------------------------------------------------------------------

		if (timeron) timer.start(t_xsolve);
		for (k = 1; k <= nz2; k++)
		{
			for (j = 1; j <= ny2; j++)
			{

				//---------------------------------------------------------------------
				// Computes the left hand side for the three x-factors  
				//---------------------------------------------------------------------

				//---------------------------------------------------------------------
				//      first fill the lhs for the u-eigenvalue                   
				//---------------------------------------------------------------------
				for (i = 0; i <= grid_points[0] - 1; i++)
				{
					ru1 = c3c4 * rho_i[i][j][k];
					cv[i] = us[i][j][k];
					rhon[i] = dmax1(dx2 + con43 * ru1,
									dx5 + c1c5 * ru1,
									dxmax + ru1,
									dx1);
				}

				lhsinit(grid_points[0] - 1);
				for (i = 1; i <= nx2; i++)
				{
					lhs[0][i] = 0.0;
					lhs[1][i] = -dttx2 * cv[(i - 1)] - dttx1 * rhon[i - 1];
					lhs[2][i] = 1.0 + c2dttx1 * rhon[i];
					lhs[3][i] = dttx2 * cv[i + 1] - dttx1 * rhon[i + 1];
					lhs[4][i] = 0.0;
				}

				//---------------------------------------------------------------------
				//      add fourth order dissipation                             
				//---------------------------------------------------------------------

				i = 1;
				lhs[2][i] = lhs[2][i] + comz5;
				lhs[3][i] = lhs[3][i] - comz4;
				lhs[4][i] = lhs[4][i] + comz1;

				lhs[1][i+1] = lhs[1][i+1] - comz4;
				lhs[2][i+1] = lhs[2][i+1] + comz6;
				lhs[3][i+1] = lhs[3][i+1] - comz4;
				lhs[4][i+1] = lhs[4][i+1] + comz1;

				for (i = 3; i <= grid_points[0] - 4; i++)
				{
					lhs[0][i] = lhs[0][i] + comz1;
					lhs[1][i] = lhs[1][i] - comz4;
					lhs[2][i] = lhs[2][i] + comz6;
					lhs[3][i] = lhs[3][i] - comz4;
					lhs[4][i] = lhs[4][i] + comz1;
				}


				i = grid_points[0] - 3;
				lhs[0][i] = lhs[0][i] + comz1;
				lhs[1][i] = lhs[1][i] - comz4;
				lhs[2][i] = lhs[2][i] + comz6;
				lhs[3][i] = lhs[3][i] - comz4;

				lhs[0][i+1] = lhs[0][i+1] + comz1;
				lhs[1][i+1] = lhs[1][i+1] - comz4;
				lhs[2][i+1] = lhs[2][i+1] + comz5;

				//---------------------------------------------------------------------
				//      subsequently, fill the other factors (u+c), (u-c) by adding to 
				//      the first  
				//---------------------------------------------------------------------
				for (i = 1; i <= nx2; i++)
				{
					lhsp[0][i] = lhs[0][i];
					lhsp[1][i] = lhs[1][i] -
									  dttx2 * speed[i-1][j][k];
					lhsp[2][i] = lhs[2][i];
					lhsp[3][i] = lhs[3][i] +
									  dttx2 * speed[i+1][j][k];
					lhsp[4][i] = lhs[4][i];
					lhsm[0][i] = lhs[0][i];
					lhsm[1][i] = lhs[1][i] +
									  dttx2 * speed[i-1][j][k];
					lhsm[2][i] = lhs[2][i];
					lhsm[3][i] = lhs[3][i] -
									  dttx2 * speed[i+1][j][k];
					lhsm[4][i] = lhs[4][i];
				}

				//---------------------------------------------------------------------
				//                          FORWARD ELIMINATION  
				//---------------------------------------------------------------------

				//---------------------------------------------------------------------
				//      perform the Thomas algorithm; first, FORWARD ELIMINATION     
				//---------------------------------------------------------------------

				for (i = 0; i <= grid_points[0] - 3; i++)
				{
					i1 = i + 1;
					i2 = i + 2;
					fac1 = 1.0 / lhs[2][i];
					lhs[3][i] = fac1 * lhs[3][i];
					lhs[4][i] = fac1 * lhs[4][i];
					for (m = 0; m <= 2; m++)
					{
						rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
					}
					lhs[2][i1] = lhs[2][i1] -
								   lhs[1][i1] * lhs[3][i];
					lhs[3][i1] = lhs[3][i1] -
								   lhs[1][i1] * lhs[4][i];
					for (m = 0; m <= 2; m++)
					{
						rhs[m][i1][j][k] = rhs[m][i1][j][k] -
									lhs[1][i1] * rhs[m][i][j][k];
					}
					lhs[1][i2] = lhs[1][i2] -
								   lhs[0][i2] * lhs[3][i];
					lhs[2][i2] = lhs[2][i2] -
								   lhs[0][i2] * lhs[4][i];
					for (m = 0; m <= 2; m++)
					{
						rhs[m][i2][j][k] = rhs[m][i2][j][k] -
									lhs[0][i2] * rhs[m][i][j][k];
					}
				}

				//---------------------------------------------------------------------
				//      The last two rows in this grid block are a bit different, 
				//      since they do not have two more rows available for the
				//      elimination of off-diagonal entries
				//---------------------------------------------------------------------

				i = grid_points[0] - 2;
				i1 = grid_points[0] - 1;
				fac1 = 1.0 / lhs[2][i];
				lhs[3][i] = fac1 * lhs[3][i];
				lhs[4][i] = fac1 * lhs[4][i];
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
				}
				lhs[2][i1] = lhs[2][i1] -
							   lhs[1][i1] * lhs[3][i];
				lhs[3][i1] = lhs[3][i1] -
							   lhs[1][i1] * lhs[4][i];
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i1][j][k] = rhs[m][i1][j][k] -
								lhs[1][i1] * rhs[m][i][j][k];
				}
				//---------------------------------------------------------------------
				//            scale the last row immediately 
				//---------------------------------------------------------------------
				fac2 = 1.0 / lhs[2][i1];
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i1][j][k] = fac2 * rhs[m][i1][j][k];
				}

				//---------------------------------------------------------------------
				//      do the u+c and the u-c factors                 
				//---------------------------------------------------------------------

				for (i = 0; i <= grid_points[0] - 3; i++)
				{
					i1 = i + 1;
					i2 = i + 2;
					m = 3;
					fac1 = 1.0 / lhsp[2][i];
					lhsp[3][i] = fac1 * lhsp[3][i];
					lhsp[4][i] = fac1 * lhsp[4][i];
					rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
					lhsp[2][i1] = lhsp[2][i1] -
								  lhsp[1][i1] * lhsp[3][i];
					lhsp[3][i1] = lhsp[3][i1] -
								  lhsp[1][i1] * lhsp[4][i];
					rhs[m][i1][j][k] = rhs[m][i1][j][k] -
								  lhsp[1][i1] * rhs[m][i][j][k];
					lhsp[1][i2] = lhsp[1][i2] -
								  lhsp[0][i2] * lhsp[3][i];
					lhsp[2][i2] = lhsp[2][i2] -
								  lhsp[0][i2] * lhsp[4][i];
					rhs[m][i2][j][k] = rhs[m][i2][j][k] -
								  lhsp[0][i2] * rhs[m][i][j][k];
					m = 4;
					fac1 = 1.0 / lhsm[2][i];
					lhsm[3][i] = fac1 * lhsm[3][i];
					lhsm[4][i] = fac1 * lhsm[4][i];
					rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
					lhsm[2][i1] = lhsm[2][i1] -
								  lhsm[1][i1] * lhsm[3][i];
					lhsm[3][i1] = lhsm[3][i1] -
								  lhsm[1][i1] * lhsm[4][i];
					rhs[m][i1][j][k] = rhs[m][i1][j][k] -
								  lhsm[1][i1] * rhs[m][i][j][k];
					lhsm[1][i2] = lhsm[1][i2] -
								  lhsm[0][i2] * lhsm[3][i];
					lhsm[2][i2] = lhsm[2][i2] -
								  lhsm[0][i2] * lhsm[4][i];
					rhs[m][i2][j][k] = rhs[m][i2][j][k] -
								  lhsm[0][i2] * rhs[m][i][j][k];
				}

				//---------------------------------------------------------------------
				//         And again the last two rows separately
				//---------------------------------------------------------------------
				i = grid_points[0] - 2;
				i1 = grid_points[0] - 1;
				m = 3;
				fac1 = 1.0 / lhsp[2][i];
				lhsp[3][i] = fac1 * lhsp[3][i];
				lhsp[4][i] = fac1 * lhsp[4][i];
				rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
				lhsp[2][i1] = lhsp[2][i1] -
							   lhsp[1][i1] * lhsp[3][i];
				lhsp[3][i1] = lhsp[3][i1] -
							   lhsp[1][i1] * lhsp[4][i];
				rhs[m][i1][j][k] = rhs[m][i1][j][k] -
							   lhsp[1][i1] * rhs[m][i][j][k];
				m = 4;
				fac1 = 1.0 / lhsm[2][i];
				lhsm[3][i] = fac1 * lhsm[3][i];
				lhsm[4][i] = fac1 * lhsm[4][i];
				rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
				lhsm[2][i1] = lhsm[2][i1] -
							   lhsm[1][i1] * lhsm[3][i];
				lhsm[3][i1] = lhsm[3][i1] -
							   lhsm[1][i1] * lhsm[4][i];
				rhs[m][i1][j][k] = rhs[m][i1][j][k] -
							   lhsm[1][i1] * rhs[m][i][j][k];
				//---------------------------------------------------------------------
				//               Scale the last row immediately
				//---------------------------------------------------------------------
				rhs[3][i1][j][k] = rhs[3][i1][j][k] / lhsp[2][i1];
				rhs[4][i1][j][k] = rhs[4][i1][j][k] / lhsm[2][i1];

				//---------------------------------------------------------------------
				//                         BACKSUBSTITUTION 
				//---------------------------------------------------------------------

				i = grid_points[0] - 2;
				i1 = grid_points[0] - 1;
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] -
									   lhs[3][i] * rhs[m][i1][j][k];
				}

				rhs[3][i][j][k] = rhs[3][i][j][k] -
								   lhsp[3][i] * rhs[3][i1][j][k];
				rhs[4][i][j][k] = rhs[4][i][j][k] -
								   lhsm[3][i] * rhs[4][i1][j][k];


				//---------------------------------------------------------------------
				//      The first three factors
				//---------------------------------------------------------------------
				for (i = grid_points[0] - 3; i >= 0; i--)
				{
					i1 = i + 1;
					i2 = i + 2;

					for (m = 0; m <= 2; m++)
					{
						rhs[m][i][j][k] = rhs[m][i][j][k] -
									 lhs[3][i] * rhs[m][i1][j][k] -
									 lhs[4][i] * rhs[m][i2][j][k];
					}
					//---------------------------------------------------------------------
					//      And the remaining two
					//---------------------------------------------------------------------
					rhs[3][i][j][k] = rhs[3][i][j][k] -
									lhsp[3][i] * rhs[3][i1][j][k] -
									lhsp[4][i] * rhs[3][i2][j][k];
					rhs[4][i][j][k] = rhs[4][i][j][k] -
									lhsm[3][i] * rhs[4][i1][j][k] -
									lhsm[4][i] * rhs[4][i2][j][k];
				}
			}
		}
		if (timeron) timer.stop(t_xsolve);

		//---------------------------------------------------------------------
		//      Do the block-diagonal inversion          
		//---------------------------------------------------------------------
		if (timeron) timer.start(t_ninvr);
		ninvr();
		if (timeron) timer.stop(t_ninvr);
	}

	public void y_solve()
	{
		int i, j, k, n, j1, j2, m;
		double ru1, fac1, fac2;

		//---------------------------------------------------------------------
		//---------------------------------------------------------------------

		if (timeron) timer.start(t_ysolve);
		for (k = 1; k <= grid_points[2] - 2; k++)
		{
			for (i = 1; i <= grid_points[0] - 2; i++)
			{

				//---------------------------------------------------------------------
				// Computes the left hand side for the three y-factors   
				//---------------------------------------------------------------------

				//---------------------------------------------------------------------
				//      first fill the lhs for the u-eigenvalue         
				//---------------------------------------------------------------------

				for (j = 0; j <= grid_points[1] - 1; j++)
				{
					ru1 = c3c4 * rho_i[i][j][k];
					cv[j] = vs[i][j][k];
					rhoq[j] = dmax1(dy3 + con43 * ru1,
									 dy5 + c1c5 * ru1,
									 dymax + ru1,
									 dy1);
				}

				lhsinit(grid_points[1] - 1);

				for (j = 1; j <= grid_points[1] - 2; j++)
				{
					lhs[0][j] = 0.0;
					lhs[1][j] = -dtty2 * cv[j - 1] - dtty1 * rhoq[j - 1];
					lhs[2][j] = 1.0 + c2dtty1 * rhoq[j];
					lhs[3][j] = dtty2 * cv[j + 1] - dtty1 * rhoq[j + 1];
					lhs[4][j] = 0.0;
				}

				//---------------------------------------------------------------------
				//      add fourth order dissipation                             
				//---------------------------------------------------------------------

				j = 1;

				lhs[2][j] = lhs[2][j] + comz5;
				lhs[3][j] = lhs[3][j] - comz4;
				lhs[4][j] = lhs[4][j] + comz1;

				lhs[1][j+1] = lhs[1][j+1] - comz4;
				lhs[2][j+1] = lhs[2][j+1] + comz6;
				lhs[3][j+1] = lhs[3][j+1] - comz4;
				lhs[4][j+1] = lhs[4][j+1] + comz1;

				for (j = 3; j <= grid_points[1] - 4; j++)
				{

					lhs[0][j] = lhs[0][j] + comz1;
					lhs[1][j] = lhs[1][j] - comz4;
					lhs[2][j] = lhs[2][j] + comz6;
					lhs[3][j] = lhs[3][j] - comz4;
					lhs[4][j] = lhs[4][j] + comz1;
				}

				j = grid_points[1] - 3;
				lhs[0][j] = lhs[0][j] + comz1;
				lhs[1][j] = lhs[1][j] - comz4;
				lhs[2][j] = lhs[2][j] + comz6;
				lhs[3][j] = lhs[3][j] - comz4;

				lhs[0][j+1] = lhs[0][j+1] + comz1;
				lhs[1][j+1] = lhs[1][j+1] - comz4;
				lhs[2][j+1] = lhs[2][j+1] + comz5;

				//---------------------------------------------------------------------
				//      subsequently, do the other two factors                    
				//---------------------------------------------------------------------
				for (j = 1; j <= grid_points[1] - 2; j++)
				{

					lhsp[0][j] = lhs[0][j];
					lhsp[1][j] = lhs[1][j] -
										   dtty2 * speed[i][j-1][k];
					lhsp[2][j] = lhs[2][j];
					lhsp[3][j] = lhs[3][j] +
								   dtty2 * speed[i][j+1][k];
					lhsp[4][j] = lhs[4][j];

					lhsm[0][j] = lhs[0][j];
					lhsm[1][j] = lhs[1][j] +
										   dtty2 * speed[i][j-1][k];
					lhsm[2][j] = lhs[2][j];
					lhsm[3][j] = lhs[3][j] -
										   dtty2 * speed[i][j+1][k];
					lhsm[4][j] = lhs[4][j];

				}

				//---------------------------------------------------------------------
				//                          FORWARD ELIMINATION  
				//---------------------------------------------------------------------

				for (j = 0; j <= grid_points[1] - 3; j++)
				{
					j1 = j + 1;
					j2 = j + 2;
					fac1 = 1.0 / lhs[2][j];
					lhs[3][j] = fac1 * lhs[3][j];
					lhs[4][j] = fac1 * lhs[4][j];
					for (m = 0; m <= 2; m++)
					{
						rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
					}
					lhs[2][j1] = lhs[2][j1] -
								   lhs[1][j1] * lhs[3][j];
					lhs[3][j1] = lhs[3][j1] -
								   lhs[1][j1] * lhs[4][j];
					for (m = 0; m <= 2; m++)
					{
						rhs[m][i][j1][k] = rhs[m][i][j1][k] -
									lhs[1][j1] * rhs[m][i][j][k];
					}
					lhs[1][j2] = lhs[1][j2] -
								   lhs[0][j2] * lhs[3][j];
					lhs[2][j2] = lhs[2][j2] -
								   lhs[0][j2] * lhs[4][j];
					for (m = 0; m <= 2; m++)
					{
						rhs[m][i][j2][k] = rhs[m][i][j2][k] -
									lhs[0][j2] * rhs[m][i][j][k];
					}
				}

				//---------------------------------------------------------------------
				//      The last two rows in this grid block are a bit different, 
				//      since they do not have two more rows available for the
				//      elimination of off-diagonal entries
				//---------------------------------------------------------------------

				j = grid_points[1] - 2;
				j1 = grid_points[1] - 1;
				fac1 = 1.0 / lhs[2][j];
				lhs[3][j] = fac1 * lhs[3][j];
				lhs[4][j] = fac1 * lhs[4][j];
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
				}
				lhs[2][j1] = lhs[2][j1] -
							   lhs[1][j1] * lhs[3][j];
				lhs[3][j1] = lhs[3][j1] -
							   lhs[1][j1] * lhs[4][j];
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i][j1][k] = rhs[m][i][j1][k] -
								lhs[1][j1] * rhs[m][i][j][k];
				}
				//---------------------------------------------------------------------
				//            scale the last row immediately 
				//---------------------------------------------------------------------
				fac2 = 1.0 / lhs[2][j1];
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i][j1][k] = fac2 * rhs[m][i][j1][k];
				}

				//---------------------------------------------------------------------
				//      do the u+c and the u-c factors                 
				//---------------------------------------------------------------------
				for (j = 0; j <= grid_points[1] - 3; j++)
				{
					j1 = j + 1;
					j2 = j + 2;
					m = 3;
					fac1 = 1.0 / lhsp[2][j];
					lhsp[3][j] = fac1 * lhsp[3][j];
					lhsp[4][j] = fac1 * lhsp[4][j];
					rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
					lhsp[2][j1] = lhsp[2][j1] -
							lhsp[1][j1] * lhsp[3][j];
					lhsp[3][j1] = lhsp[3][j1] -
							lhsp[1][j1] * lhsp[4][j];
					rhs[m][i][j1][k] = rhs[m][i][j1][k] -
							lhsp[1][j1] * rhs[m][i][j][k];
					lhsp[1][j2] = lhsp[1][j2] -
							lhsp[0][j2] * lhsp[3][j];
					lhsp[2][j2] = lhsp[2][j2] -
							lhsp[0][j2] * lhsp[4][j];
					rhs[m][i][j2][k] = rhs[m][i][j2][k] -
							lhsp[0][j2] * rhs[m][i][j][k];
					m = 4;
					fac1 = 1.0 / lhsm[2][j];
					lhsm[3][j] = fac1 * lhsm[3][j];
					lhsm[4][j] = fac1 * lhsm[4][j];
					rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
					lhsm[2][j1] = lhsm[2][j1] -
							lhsm[1][j1] * lhsm[3][j];
					lhsm[3][j1] = lhsm[3][j1] -
							lhsm[1][j1] * lhsm[4][j];
					rhs[m][i][j1][k] = rhs[m][i][j1][k] -
							lhsm[1][j1] * rhs[m][i][j][k];
					lhsm[1][j2] = lhsm[1][j2] -
							lhsm[0][j2] * lhsm[3][j];
					lhsm[2][j2] = lhsm[2][j2] -
							lhsm[0][j2] * lhsm[4][j];
					rhs[m][i][j2][k] = rhs[m][i][j2][k] -
							lhsm[0][j2] * rhs[m][i][j][k];
				}

				//---------------------------------------------------------------------
				//         And again the last two rows separately
				//---------------------------------------------------------------------
				j = grid_points[1] - 2;
				j1 = grid_points[1] - 1;
				m = 3;
				fac1 = 1.0 / lhsp[2][j];
				lhsp[3][j] = fac1 * lhsp[3][j];
				lhsp[4][j] = fac1 * lhsp[4][j];
				rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
				lhsp[2][j1] = lhsp[2][j1] -
						lhsp[1][j1] * lhsp[3][j];
				lhsp[3][j1] = lhsp[3][j1] -
						lhsp[1][j1] * lhsp[4][j];
				rhs[m][i][j1][k] = rhs[m][i][j1][k] -
						lhsp[1][j1] * rhs[m][i][j][k];
				m = 4;
				fac1 = 1.0 / lhsm[2][j];
				lhsm[3][j] = fac1 * lhsm[3][j];
				lhsm[4][j] = fac1 * lhsm[4][j];
				rhs[m][i][j][k] = fac1 * rhs[m][i][j][k];
				lhsm[2][j1] = lhsm[2][j1] -
						lhsm[1][j1] * lhsm[3][j];
				lhsm[3][j1] = lhsm[3][j1] -
						lhsm[1][j1] * lhsm[4][j];
				rhs[m][i][j1][k] = rhs[m][i][j1][k] -
						lhsm[1][j1] * rhs[m][i][j][k];
				//---------------------------------------------------------------------
				//               Scale the last row immediately 
				//---------------------------------------------------------------------
				rhs[3][i][j1][k] = rhs[3][i][j1][k] / lhsp[2][j1];
				rhs[4][i][j1][k] = rhs[4][i][j1][k] / lhsm[2][j1];

				//---------------------------------------------------------------------
				//                         BACKSUBSTITUTION 
				//---------------------------------------------------------------------

				j = grid_points[1] - 2;
				j1 = grid_points[1] - 1;
				for (m = 0; m <= 2; m++)
				{
					rhs[m][i][j][k] = rhs[m][i][j][k] -
									 lhs[3][j] * rhs[m][i][j1][k];
				}

				rhs[3][i][j][k] = rhs[3][i][j][k] -
									lhsp[3][j] * rhs[3][i][j1][k];
				rhs[4][i][j][k] = rhs[4][i][j][k] -
									lhsm[3][j] * rhs[4][i][j1][k];

				//---------------------------------------------------------------------
				//      The first three factors
				//---------------------------------------------------------------------
				for (j = grid_points[1] - 3; j >= 0; j--)
				{
					j1 = j + 1;
					j2 = j + 2;
					for (m = 0; m <= 2; m++)
					{
						rhs[m][i][j][k] = rhs[m][i][j][k] -
									 lhs[3][j] * rhs[m][i][j1][k] -
									 lhs[4][j] * rhs[m][i][j2][k];
					}

					//---------------------------------------------------------------------
					//      And the remaining two
					//---------------------------------------------------------------------
					rhs[3][i][j][k] = rhs[3][i][j][k] -
									lhsp[3][j] * rhs[3][i][j1][k] -
									lhsp[4][j] * rhs[3][i][j2][k];
					rhs[4][i][j][k] = rhs[4][i][j][k] -
									lhsm[3][j] * rhs[4][i][j1][k] -
									lhsm[4][j] * rhs[4][i][j2][k];
				}
			}
		}
		if (timeron) timer.stop(t_ysolve);

		if (timeron) timer.start(t_pinvr);
		pinvr();
		if (timeron) timer.stop(t_pinvr);
	}

	public void z_solve()
	{
		int i, j, k, n, k1, k2, m;
		double ru1, fac1, fac2;
        double[] rtmp = new double[5*(KMAX+1)];


		//---------------------------------------------------------------------
		//---------------------------------------------------------------------

		if (timeron) timer.start(t_zsolve);
		for (j = 1; j <= ny2; j++)
		{
			for (i = 1; i <= nx2; i++)
			{

				//---------------------------------------------------------------------
				// Computes the left hand side for the three z-factors   
				//---------------------------------------------------------------------

				//---------------------------------------------------------------------
				// first fill the lhs for the u-eigenvalue                          
				//---------------------------------------------------------------------

				for (k = 0; k <= nz2 + 1; k++)
				{
					ru1 = c3c4 * rho_i[i][j][k];
					cv[k] = ws[i][j][k];
					rhos[k] = dmax1(dz4 + con43 * ru1,
									dz5 + c1c5 * ru1,
									dzmax + ru1,
									dz1);
				}

				lhsinit(grid_points[2] - 1);
				for (k = 1; k <= nz2; k++)
				{
					lhs[0][k] = 0.0;
					lhs[1][k] = -dttz2 * cv[k - 1] - dttz1 * rhos[k - 1];
					lhs[2][k] = 1.0 + c2dttz1 * rhos[k];
					lhs[3][k] = dttz2 * cv[k + 1] - dttz1 * rhos[k + 1];
					lhs[4][k] = 0.0;
				}

				//---------------------------------------------------------------------
				//      add fourth order dissipation                                  
				//---------------------------------------------------------------------

				k = 1;
				lhs[2][k] = lhs[2][k] + comz5;
				lhs[3][k] = lhs[3][k] - comz4;
				lhs[4][k] = lhs[4][k] + comz1;

				k = 2;
				lhs[1][k] = lhs[1][k] - comz4;
				lhs[2][k] = lhs[2][k] + comz6;
				lhs[3][k] = lhs[3][k] - comz4;
				lhs[4][k] = lhs[4][k] + comz1;

				for (k = 3; k <= nz2 - 2; k++)
				{
					lhs[0][k] = lhs[0][k] + comz1;
					lhs[1][k] = lhs[1][k] - comz4;
					lhs[2][k] = lhs[2][k] + comz6;
					lhs[3][k] = lhs[3][k] - comz4;
					lhs[4][k] = lhs[4][k] + comz1;
				}

				k = nz2 - 1;
				lhs[0][k] = lhs[0][k] + comz1;
				lhs[1][k] = lhs[1][k] - comz4;
				lhs[2][k] = lhs[2][k] + comz6;
				lhs[3][k] = lhs[3][k] - comz4;

				k = nz2;
				lhs[0][k] = lhs[0][k] + comz1;
				lhs[1][k] = lhs[1][k] - comz4;
				lhs[2][k] = lhs[2][k] + comz5;

				//---------------------------------------------------------------------
				//      subsequently, fill the other factors (u+c), (u-c) 
				//---------------------------------------------------------------------
				for (k = 1; k <= nz2; k++)
				{
					lhsp[0][k] = lhs[0][k];
					lhsp[1][k] = lhs[1][k] - dttz2 * speed[i][j][k-1];
					lhsp[2][k] = lhs[2][k];
					lhsp[3][k] = lhs[3][k] + dttz2 * speed[i][j][k+1];
					lhsp[4][k] = lhs[4][k];

					lhsm[0][k] = lhs[0][k];
					lhsm[1][k] = lhs[1][k] + dttz2 * speed[i][j][k-1];
					lhsm[2][k] = lhs[2][k];
					lhsm[3][k] = lhs[3][k] - dttz2 * speed[i][j][k+1];
					lhsm[4][k] = lhs[4][k];
				}

				//---------------------------------------------------------------------
				// Load a row of K data
				//---------------------------------------------------------------------

				for (k = 0; k <= nz2 + 1; k++)
				{
					rtmp[0 + k * 5] = rhs[0][i][j][k];
					rtmp[1 + k * 5] = rhs[1][i][j][k];
					rtmp[2 + k * 5] = rhs[2][i][j][k];
					rtmp[3 + k * 5] = rhs[3][i][j][k];
					rtmp[4 + k * 5] = rhs[4][i][j][k];
				}

				//---------------------------------------------------------------------
				//                          FORWARD ELIMINATION  
				//---------------------------------------------------------------------

				for (k = 0; k <= grid_points[2] - 3; k++)
				{
					k1 = k + 1;
					k2 = k + 2;
					fac1 = 1.0 / lhs[2][k];
					lhs[3][k] = fac1 * lhs[3][k];
					lhs[4][k] = fac1 * lhs[4][k];
					for (m = 0; m <= 2; m++)
					{
						rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
					}
					lhs[2][k1] = lhs[2][k1] -
								   lhs[1][k1] * lhs[3][k];
					lhs[3][k1] = lhs[3][k1] -
								   lhs[1][k1] * lhs[4][k];
					for (m = 0; m <= 2; m++)
					{
						rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
									lhs[1][k1] * rtmp[m + k * 5];
					}
					lhs[1][k2] = lhs[1][k2] -
								   lhs[0][k2] * lhs[3][k];
					lhs[2][k2] = lhs[2][k2] -
								   lhs[0][k2] * lhs[4][k];
					for (m = 0; m <= 2; m++)
					{
						rtmp[m + k2 * 5] = rtmp[m + k2 * 5] -
									lhs[0][k2] * rtmp[m + k * 5];
					}
				}

				//---------------------------------------------------------------------
				//      The last two rows in this grid block are a bit different, 
				//      since they do not have two more rows available for the
				//      elimination of off-diagonal entries
				//---------------------------------------------------------------------
				k = grid_points[2] - 2;
				k1 = grid_points[2] - 1;
				fac1 = 1.0 / lhs[2][k];
				lhs[3][k] = fac1 * lhs[3][k];
				lhs[4][k] = fac1 * lhs[4][k];
				for (m = 0; m <= 2; m++)
				{
					rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
				}
				lhs[2][k1] = lhs[2][k1] -
							   lhs[1][k1] * lhs[3][k];
				lhs[3][k1] = lhs[3][k1] -
							   lhs[1][k1] * lhs[4][k];
				for (m = 0; m <= 2; m++)
				{
					rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
								lhs[1][k1] * rtmp[m + k * 5];
				}
				//---------------------------------------------------------------------
				//               scale the last row immediately
				//---------------------------------------------------------------------
				fac2 = 1.0 / lhs[2][k1];
				for (m = 0; m <= 2; m++)
				{
					rtmp[m + k1 * 5] = fac2 * rtmp[m + k1 * 5];
				}

				//---------------------------------------------------------------------
				//      do the u+c and the u-c factors               
				//---------------------------------------------------------------------
				for (k = 0; k <= grid_points[2] - 3; k++)
				{
					k1 = k + 1;
					k2 = k + 2;
					m = 3;
					fac1 = 1.0 / lhsp[2][k];
					lhsp[3][k] = fac1 * lhsp[3][k];
					lhsp[4][k] = fac1 * lhsp[4][k];
					rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
					lhsp[2][k1] = lhsp[2][k1] -
							lhsp[1][k1] * lhsp[3][k];
					lhsp[3][k1] = lhsp[3][k1] -
							lhsp[1][k1] * lhsp[4][k];
					rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
							lhsp[1][k1] * rtmp[m + k * 5];
					lhsp[1][k2] = lhsp[1][k2] -
							lhsp[0][k2] * lhsp[3][k];
					lhsp[2][k2] = lhsp[2][k2] -
							lhsp[0][k2] * lhsp[4][k];
					rtmp[m + k2 * 5] = rtmp[m + k2 * 5] -
							lhsp[0][k2] * rtmp[m + k * 5];
					m = 4;
					fac1 = 1.0 / lhsm[2][k];
					lhsm[3][k] = fac1 * lhsm[3][k];
					lhsm[4][k] = fac1 * lhsm[4][k];
					rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
					lhsm[2][k1] = lhsm[2][k1] -
							lhsm[1][k1] * lhsm[3][k];
					lhsm[3][k1] = lhsm[3][k1] -
							lhsm[1][k1] * lhsm[4][k];
					rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
							lhsm[1][k1] * rtmp[m + k * 5];
					lhsm[1][k2] = lhsm[1][k2] -
							lhsm[0][k2] * lhsm[3][k];
					lhsm[2][k2] = lhsm[2][k2] -
							lhsm[0][k2] * lhsm[4][k];
					rtmp[m + k2 * 5] = rtmp[m + k2 * 5] -
							lhsm[0][k2] * rtmp[m + k * 5];
				}

				//---------------------------------------------------------------------
				//         And again the last two rows separately
				//---------------------------------------------------------------------
				k = grid_points[2] - 2;
				k1 = grid_points[2] - 1;
				m = 3;
				fac1 = 1.0 / lhsp[2][k];
				lhsp[3][k] = fac1 * lhsp[3][k];
				lhsp[4][k] = fac1 * lhsp[4][k];
				rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
				lhsp[2][k1] = lhsp[2][k1] -
						lhsp[1][k1] * lhsp[3][k];
				lhsp[3][k1] = lhsp[3][k1] -
						lhsp[1][k1] * lhsp[4][k];
				rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
						lhsp[1][k1] * rtmp[m + k * 5];
				m = 4;
				fac1 = 1.0 / lhsm[2][k];
				lhsm[3][k] = fac1 * lhsm[3][k];
				lhsm[4][k] = fac1 * lhsm[4][k];
				rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
				lhsm[2][k1] = lhsm[2][k1] -
						lhsm[1][k1] * lhsm[3][k];
				lhsm[3][k1] = lhsm[3][k1] -
						lhsm[1][k1] * lhsm[4][k];
				rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
						lhsm[1][k1] * rtmp[m + k * 5];
				//---------------------------------------------------------------------
				//               Scale the last row immediately (some of this is overkill
				//               if this is the last cell)
				//---------------------------------------------------------------------
				rtmp[3 + k1 * 5] = rtmp[3 + k1 * 5] / lhsp[2][k1];
				rtmp[4 + k1 * 5] = rtmp[4 + k1 * 5] / lhsm[2][k1];

				//---------------------------------------------------------------------
				//                         BACKSUBSTITUTION 
				//---------------------------------------------------------------------

				k = grid_points[2] - 2;
				k1 = grid_points[2] - 1;
				for (m = 0; m <= 2; m++)
				{
					rtmp[m + k * 5] = rtmp[m + k * 5] -
									   lhs[3][k] * rtmp[m + k1 * 5];
				}

				rtmp[3 + k * 5] = rtmp[3 + k * 5] -
									  lhsp[3][k] * rtmp[3 + k1 * 5];
				rtmp[4 + k * 5] = rtmp[4 + k * 5] -
									  lhsm[3][k] * rtmp[4 + k1 * 5];

				//---------------------------------------------------------------------
				//      Whether or not this is the last processor, we always have
				//      to complete the back-substitution 
				//---------------------------------------------------------------------

				//---------------------------------------------------------------------
				//      The first three factors
				//---------------------------------------------------------------------
				for (k = grid_points[2] - 3; k >= 0; k--)
				{
					k1 = k + 1;
					k2 = k + 2;
					for (m = 0; m <= 2; m++)
					{
						rtmp[m + k * 5] = rtmp[m + k * 5] -
									 lhs[3][k] * rtmp[m + k1 * 5] -
									 lhs[4][k] * rtmp[m + k2 * 5];
					}

					//---------------------------------------------------------------------
					//      And the remaining two
					//---------------------------------------------------------------------
					rtmp[3 + k * 5] = rtmp[3 + k * 5] -
									lhsp[3][k] * rtmp[3 + k1 * 5] -
									lhsp[4][k] * rtmp[3 + k2 * 5];
					rtmp[4 + k * 5] = rtmp[4 + k * 5] -
									lhsm[3][k] * rtmp[4 + k1 * 5] -
									lhsm[4][k] * rtmp[4 + k2 * 5];
				}

				//---------------------------------------------------------------------
				//      Store result
				//---------------------------------------------------------------------
				for (k = 0; k <= nz2 + 1; k++)
				{
					rhs[0][i][j][k] = rtmp[0 + k * 5];
					rhs[1][i][j][k] = rtmp[1 + k * 5];
					rhs[2][i][j][k] = rtmp[2 + k * 5];
					rhs[3][i][j][k] = rtmp[3 + k * 5];
					rhs[4][i][j][k] = rtmp[4 + k * 5];
				}

			}
		}
		if (timeron) timer.stop(t_zsolve);
		if (timeron) timer.start(t_tzetar);
		tzetar();
		if (timeron) timer.stop(t_tzetar);
	}

	public double getTime() { return timer.readTimer(1); }

	public void finalize() 
	{
	    Console.WriteLine("SP: is about to be garbage collected");
		//base.finalize();
	}
}

}




