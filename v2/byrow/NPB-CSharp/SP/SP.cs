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
						u[k][j][i][m] +=
								   rhs[k][j][i][m];
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
						add = u[k][j][i][m] - u_exact[m];
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
						add = rhs[k][j][i][m];
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

					r1 = rhs[k][j][i][0];
					r2 = rhs[k][j][i][1];
					r3 = rhs[k][j][i][2];
					r4 = rhs[k][j][i][3];
					r5 = rhs[k][j][i][4];

					t1 = bt * r3;
					t2 = 0.5 * (r4 + r5);

					rhs[k][j][i][0] = -r2;
					rhs[k][j][i][1] = r1;
					rhs[k][j][i][2] = bt * (r4 - r5);
					rhs[k][j][i][3] = -t1 + t2;
					rhs[k][j][i][4] = t1 + t2;
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

					r1 = rhs[k][j][i][0];
					r2 = rhs[k][j][i][1];
					r3 = rhs[k][j][i][2];
					r4 = rhs[k][j][i][3];
					r5 = rhs[k][j][i][4];

					t1 = bt * r1;
					t2 = 0.5 * (r4 + r5);

					rhs[k][j][i][0] = bt * (r4 - r5);
					rhs[k][j][i][1] = -r3;
					rhs[k][j][i][2] = r2;
					rhs[k][j][i][3] = -t1 + t2;
					rhs[k][j][i][4] = t1 + t2;
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
					//---------------------------------------------------------------------
					//               (don't need speed and ainx until the lhs computation)
					//---------------------------------------------------------------------
					aux = c1c2 * rho_inv * (u[k][j][i][4] - square[k][j][i]);
					speed[k][j][i] = Math.Sqrt(aux);
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

					ru1 = rho_i[k][j][i];
					uu = us[k][j][i];
					vv = vs[k][j][i];
					ww = ws[k][j][i];
					ac = speed[k][j][i];
					ac2inv = 1.0 / (ac * ac);

					r1 = rhs[k][j][i][0];
					r2 = rhs[k][j][i][1];
					r3 = rhs[k][j][i][2];
					r4 = rhs[k][j][i][3];
					r5 = rhs[k][j][i][4];

					t1 = c2 * ac2inv * (qs[k][j][i] * r1 - uu * r2 -
						vv * r3 - ww * r4 + r5);
					t2 = bt * ru1 * (uu * r1 - r2);
					t3 = (bt * ru1 * ac) * t1;

					rhs[k][j][i][0] = r1 - t1;
					rhs[k][j][i][1] = -ru1 * (ww * r1 - r4);
					rhs[k][j][i][2] = ru1 * (vv * r1 - r3);
					rhs[k][j][i][3] = -t2 + t3;
					rhs[k][j][i][4] = t2 + t3;

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

					xvel = us[k][j][i];
					yvel = vs[k][j][i];
					zvel = ws[k][j][i];
					ac = speed[k][j][i];

					ac2u = ac * ac;

					r1 = rhs[k][j][i][0];
					r2 = rhs[k][j][i][1];
					r3 = rhs[k][j][i][2];
					r4 = rhs[k][j][i][3];
					r5 = rhs[k][j][i][4];

					uzik1 = u[k][j][i][0];
					btuz = bt * uzik1;

					t1 = btuz / ac * (r4 + r5);
					t2 = r3 + t1;
					t3 = btuz * (r4 - r5);

					rhs[k][j][i][0] = t2;
					rhs[k][j][i][1] = -uzik1 * r2 + xvel * t2;
					rhs[k][j][i][2] = uzik1 * r1 + yvel * t2;
					rhs[k][j][i][3] = zvel * t2 + t3;
					rhs[k][j][i][4] = uzik1 * (-xvel * r2 + yvel * r1) +
						  qs[k][j][i] * t2 + c2iv * ac2u * t1 + zvel * t3;

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
					ru1 = c3c4 * rho_i[k][j][i];
					cv[i] = us[k][j][i];
					rhon[i] = dmax1(dx2 + con43 * ru1,
									dx5 + c1c5 * ru1,
									dxmax + ru1,
									dx1);
				}

				lhsinit(grid_points[0] - 1);
				for (i = 1; i <= nx2; i++)
				{
					lhs[i][0] = 0.0;
					lhs[i][1] = -dttx2 * cv[(i - 1)] - dttx1 * rhon[i - 1];
					lhs[i][2] = 1.0 + c2dttx1 * rhon[i];
					lhs[i][3] = dttx2 * cv[i + 1] - dttx1 * rhon[i + 1];
					lhs[i][4] = 0.0;
				}

				//---------------------------------------------------------------------
				//      add fourth order dissipation                             
				//---------------------------------------------------------------------

				i = 1;
				lhs[i][2] = lhs[i][2] + comz5;
				lhs[i][3] = lhs[i][3] - comz4;
				lhs[i][4] = lhs[i][4] + comz1;

				lhs[i+1][1] = lhs[i+1][1] - comz4;
				lhs[i+1][2] = lhs[i+1][2] + comz6;
				lhs[i+1][3] = lhs[i+1][3] - comz4;
				lhs[i+1][4] = lhs[i+1][4] + comz1;

				for (i = 3; i <= grid_points[0] - 4; i++)
				{
					lhs[i][0] = lhs[i][0] + comz1;
					lhs[i][1] = lhs[i][1] - comz4;
					lhs[i][2] = lhs[i][2] + comz6;
					lhs[i][3] = lhs[i][3] - comz4;
					lhs[i][4] = lhs[i][4] + comz1;
				}


				i = grid_points[0] - 3;
				lhs[i][0] = lhs[i][0] + comz1;
				lhs[i][1] = lhs[i][1] - comz4;
				lhs[i][2] = lhs[i][2] + comz6;
				lhs[i][3] = lhs[i][3] - comz4;

				lhs[i+1][0] = lhs[i+1][0] + comz1;
				lhs[i+1][1] = lhs[i+1][1] - comz4;
				lhs[i+1][2] = lhs[i+1][2] + comz5;

				//---------------------------------------------------------------------
				//      subsequently, fill the other factors (u+c), (u-c) by adding to 
				//      the first  
				//---------------------------------------------------------------------
				for (i = 1; i <= nx2; i++)
				{
					lhsp[i][0] = lhs[i][0];
					lhsp[i][1] = lhs[i][1] -
									  dttx2 * speed[k][j][i-1];
					lhsp[i][2] = lhs[i][2];
					lhsp[i][3] = lhs[i][3] +
									  dttx2 * speed[k][j][i+1];
					lhsp[i][4] = lhs[i][4];
					lhsm[i][0] = lhs[i][0];
					lhsm[i][1] = lhs[i][1] +
									  dttx2 * speed[k][j][i-1];
					lhsm[i][2] = lhs[i][2];
					lhsm[i][3] = lhs[i][3] -
									  dttx2 * speed[k][j][i+1];
					lhsm[i][4] = lhs[i][4];
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
					fac1 = 1.0 / lhs[i][2];
					lhs[i][3] = fac1 * lhs[i][3];
					lhs[i][4] = fac1 * lhs[i][4];
					for (m = 0; m <= 2; m++)
					{
						rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
					}
					lhs[i1][2] = lhs[i1][2] -
								   lhs[i1][1] * lhs[i][3];
					lhs[i1][3] = lhs[i1][3] -
								   lhs[i1][1] * lhs[i][4];
					for (m = 0; m <= 2; m++)
					{
						rhs[k][j][i1][m] = rhs[k][j][i1][m] -
									lhs[i1][1] * rhs[k][j][i][m];
					}
					lhs[i2][1] = lhs[i2][1] -
								   lhs[i2][0] * lhs[i][3];
					lhs[i2][2] = lhs[i2][2] -
								   lhs[i2][0] * lhs[i][4];
					for (m = 0; m <= 2; m++)
					{
						rhs[k][j][i2][m] = rhs[k][j][i2][m] -
									lhs[i2][0] * rhs[k][j][i][m];
					}
				}

				//---------------------------------------------------------------------
				//      The last two rows in this grid block are a bit different, 
				//      since they do not have two more rows available for the
				//      elimination of off-diagonal entries
				//---------------------------------------------------------------------

				i = grid_points[0] - 2;
				i1 = grid_points[0] - 1;
				fac1 = 1.0 / lhs[i][2];
				lhs[i][3] = fac1 * lhs[i][3];
				lhs[i][4] = fac1 * lhs[i][4];
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
				}
				lhs[i1][2] = lhs[i1][2] -
							   lhs[i1][1] * lhs[i][3];
				lhs[i1][3] = lhs[i1][3] -
							   lhs[i1][1] * lhs[i][4];
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j][i1][m] = rhs[k][j][i1][m] -
								lhs[i1][1] * rhs[k][j][i][m];
				}
				//---------------------------------------------------------------------
				//            scale the last row immediately 
				//---------------------------------------------------------------------
				fac2 = 1.0 / lhs[i1][2];
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j][i1][m] = fac2 * rhs[k][j][i1][m];
				}

				//---------------------------------------------------------------------
				//      do the u+c and the u-c factors                 
				//---------------------------------------------------------------------

				for (i = 0; i <= grid_points[0] - 3; i++)
				{
					i1 = i + 1;
					i2 = i + 2;
					m = 3;
					fac1 = 1.0 / lhsp[i][2];
					lhsp[i][3] = fac1 * lhsp[i][3];
					lhsp[i][4] = fac1 * lhsp[i][4];
					rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
					lhsp[i1][2] = lhsp[i1][2] -
								  lhsp[i1][1] * lhsp[i][3];
					lhsp[i1][3] = lhsp[i1][3] -
								  lhsp[i1][1] * lhsp[i][4];
					rhs[k][j][i1][m] = rhs[k][j][i1][m] -
								  lhsp[i1][1] * rhs[k][j][i][m];
					lhsp[i2][1] = lhsp[i2][1] -
								  lhsp[i2][0] * lhsp[i][3];
					lhsp[i2][2] = lhsp[i2][2] -
								  lhsp[i2][0] * lhsp[i][4];
					rhs[k][j][i2][m] = rhs[k][j][i2][m] -
								  lhsp[i2][0] * rhs[k][j][i][m];
					m = 4;
					fac1 = 1.0 / lhsm[i][2];
					lhsm[i][3] = fac1 * lhsm[i][3];
					lhsm[i][4] = fac1 * lhsm[i][4];
					rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
					lhsm[i1][2] = lhsm[i1][2] -
								  lhsm[i1][1] * lhsm[i][3];
					lhsm[i1][3] = lhsm[i1][3] -
								  lhsm[i1][1] * lhsm[i][4];
					rhs[k][j][i1][m] = rhs[k][j][i1][m] -
								  lhsm[i1][1] * rhs[k][j][i][m];
					lhsm[i2][1] = lhsm[i2][1] -
								  lhsm[i2][0] * lhsm[i][3];
					lhsm[i2][2] = lhsm[i2][2] -
								  lhsm[i2][0] * lhsm[i][4];
					rhs[k][j][i2][m] = rhs[k][j][i2][m] -
								  lhsm[i2][0] * rhs[k][j][i][m];
				}

				//---------------------------------------------------------------------
				//         And again the last two rows separately
				//---------------------------------------------------------------------
				i = grid_points[0] - 2;
				i1 = grid_points[0] - 1;
				m = 3;
				fac1 = 1.0 / lhsp[i][2];
				lhsp[i][3] = fac1 * lhsp[i][3];
				lhsp[i][4] = fac1 * lhsp[i][4];
				rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
				lhsp[i1][2] = lhsp[i1][2] -
							   lhsp[i1][1] * lhsp[i][3];
				lhsp[i1][3] = lhsp[i1][3] -
							   lhsp[i1][1] * lhsp[i][4];
				rhs[k][j][i1][m] = rhs[k][j][i1][m] -
							   lhsp[i1][1] * rhs[k][j][i][m];
				m = 4;
				fac1 = 1.0 / lhsm[i][2];
				lhsm[i][3] = fac1 * lhsm[i][3];
				lhsm[i][4] = fac1 * lhsm[i][4];
				rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
				lhsm[i1][2] = lhsm[i1][2] -
							   lhsm[i1][1] * lhsm[i][3];
				lhsm[i1][3] = lhsm[i1][3] -
							   lhsm[i1][1] * lhsm[i][4];
				rhs[k][j][i1][m] = rhs[k][j][i1][m] -
							   lhsm[i1][1] * rhs[k][j][i][m];
				//---------------------------------------------------------------------
				//               Scale the last row immediately
				//---------------------------------------------------------------------
				rhs[k][j][i1][3] = rhs[k][j][i1][3] / lhsp[i1][2];
				rhs[k][j][i1][4] = rhs[k][j][i1][4] / lhsm[i1][2];

				//---------------------------------------------------------------------
				//                         BACKSUBSTITUTION 
				//---------------------------------------------------------------------

				i = grid_points[0] - 2;
				i1 = grid_points[0] - 1;
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] -
									   lhs[i][3] * rhs[k][j][i1][m];
				}

				rhs[k][j][i][3] = rhs[k][j][i][3] -
								   lhsp[i][3] * rhs[k][j][i1][3];
				rhs[k][j][i][4] = rhs[k][j][i][4] -
								   lhsm[i][3] * rhs[k][j][i1][4];


				//---------------------------------------------------------------------
				//      The first three factors
				//---------------------------------------------------------------------
				for (i = grid_points[0] - 3; i >= 0; i--)
				{
					i1 = i + 1;
					i2 = i + 2;

					for (m = 0; m <= 2; m++)
					{
						rhs[k][j][i][m] = rhs[k][j][i][m] -
									 lhs[i][3] * rhs[k][j][i1][m] -
									 lhs[i][4] * rhs[k][j][i2][m];
					}
					//---------------------------------------------------------------------
					//      And the remaining two
					//---------------------------------------------------------------------
					rhs[k][j][i][3] = rhs[k][j][i][3] -
									lhsp[i][3] * rhs[k][j][i1][3] -
									lhsp[i][4] * rhs[k][j][i2][3];
					rhs[k][j][i][4] = rhs[k][j][i][4] -
									lhsm[i][3] * rhs[k][j][i1][4] -
									lhsm[i][4] * rhs[k][j][i2][4];
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
					ru1 = c3c4 * rho_i[k][j][i];
					cv[j] = vs[k][j][i];
					rhoq[j] = dmax1(dy3 + con43 * ru1,
									 dy5 + c1c5 * ru1,
									 dymax + ru1,
									 dy1);
				}

				lhsinit(grid_points[1] - 1);

				for (j = 1; j <= grid_points[1] - 2; j++)
				{
					lhs[j][0] = 0.0;
					lhs[j][1] = -dtty2 * cv[j - 1] - dtty1 * rhoq[j - 1];
					lhs[j][2] = 1.0 + c2dtty1 * rhoq[j];
					lhs[j][3] = dtty2 * cv[j + 1] - dtty1 * rhoq[j + 1];
					lhs[j][4] = 0.0;
				}

				//---------------------------------------------------------------------
				//      add fourth order dissipation                             
				//---------------------------------------------------------------------

				j = 1;

				lhs[j][2] = lhs[j][2] + comz5;
				lhs[j][3] = lhs[j][3] - comz4;
				lhs[j][4] = lhs[j][4] + comz1;

				lhs[j+1][1] = lhs[j+1][1] - comz4;
				lhs[j+1][2] = lhs[j+1][2] + comz6;
				lhs[j+1][3] = lhs[j+1][3] - comz4;
				lhs[j+1][4] = lhs[j+1][4] + comz1;

				for (j = 3; j <= grid_points[1] - 4; j++)
				{

					lhs[j][0] = lhs[j][0] + comz1;
					lhs[j][1] = lhs[j][1] - comz4;
					lhs[j][2] = lhs[j][2] + comz6;
					lhs[j][3] = lhs[j][3] - comz4;
					lhs[j][4] = lhs[j][4] + comz1;
				}

				j = grid_points[1] - 3;
				lhs[j][0] = lhs[j][0] + comz1;
				lhs[j][1] = lhs[j][1] - comz4;
				lhs[j][2] = lhs[j][2] + comz6;
				lhs[j][3] = lhs[j][3] - comz4;

				lhs[j+1][0] = lhs[j+1][0] + comz1;
				lhs[j+1][1] = lhs[j+1][1] - comz4;
				lhs[j+1][2] = lhs[j+1][2] + comz5;

				//---------------------------------------------------------------------
				//      subsequently, do the other two factors                    
				//---------------------------------------------------------------------
				for (j = 1; j <= grid_points[1] - 2; j++)
				{

					lhsp[j][0] = lhs[j][0];
					lhsp[j][1] = lhs[j][1] -
										   dtty2 * speed[k][j-1][i];
					lhsp[j][2] = lhs[j][2];
					lhsp[j][3] = lhs[j][3] +
								   dtty2 * speed[k][j+1][i];
					lhsp[j][4] = lhs[j][4];

					lhsm[j][0] = lhs[j][0];
					lhsm[j][1] = lhs[j][1] +
										   dtty2 * speed[k][j-1][i];
					lhsm[j][2] = lhs[j][2];
					lhsm[j][3] = lhs[j][3] -
										   dtty2 * speed[k][j+1][i];
					lhsm[j][4] = lhs[j][4];

				}

				//---------------------------------------------------------------------
				//                          FORWARD ELIMINATION  
				//---------------------------------------------------------------------

				for (j = 0; j <= grid_points[1] - 3; j++)
				{
					j1 = j + 1;
					j2 = j + 2;
					fac1 = 1.0 / lhs[j][2];
					lhs[j][3] = fac1 * lhs[j][3];
					lhs[j][4] = fac1 * lhs[j][4];
					for (m = 0; m <= 2; m++)
					{
						rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
					}
					lhs[j1][2] = lhs[j1][2] -
								   lhs[j1][1] * lhs[j][3];
					lhs[j1][3] = lhs[j1][3] -
								   lhs[j1][1] * lhs[j][4];
					for (m = 0; m <= 2; m++)
					{
						rhs[k][j1][i][m] = rhs[k][j1][i][m] -
									lhs[j1][1] * rhs[k][j][i][m];
					}
					lhs[j2][1] = lhs[j2][1] -
								   lhs[j2][0] * lhs[j][3];
					lhs[j2][2] = lhs[j2][2] -
								   lhs[j2][0] * lhs[j][4];
					for (m = 0; m <= 2; m++)
					{
						rhs[k][j2][i][m] = rhs[k][j2][i][m] -
									lhs[j2][0] * rhs[k][j][i][m];
					}
				}

				//---------------------------------------------------------------------
				//      The last two rows in this grid block are a bit different, 
				//      since they do not have two more rows available for the
				//      elimination of off-diagonal entries
				//---------------------------------------------------------------------

				j = grid_points[1] - 2;
				j1 = grid_points[1] - 1;
				fac1 = 1.0 / lhs[j][2];
				lhs[j][3] = fac1 * lhs[j][3];
				lhs[j][4] = fac1 * lhs[j][4];
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
				}
				lhs[j1][2] = lhs[j1][2] -
							   lhs[j1][1] * lhs[j][3];
				lhs[j1][3] = lhs[j1][3] -
							   lhs[j1][1] * lhs[j][4];
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j1][i][m] = rhs[k][j1][i][m] -
								lhs[j1][1] * rhs[k][j][i][m];
				}
				//---------------------------------------------------------------------
				//            scale the last row immediately 
				//---------------------------------------------------------------------
				fac2 = 1.0 / lhs[j1][2];
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j1][i][m] = fac2 * rhs[k][j1][i][m];
				}

				//---------------------------------------------------------------------
				//      do the u+c and the u-c factors                 
				//---------------------------------------------------------------------
				for (j = 0; j <= grid_points[1] - 3; j++)
				{
					j1 = j + 1;
					j2 = j + 2;
					m = 3;
					fac1 = 1.0 / lhsp[j][2];
					lhsp[j][3] = fac1 * lhsp[j][3];
					lhsp[j][4] = fac1 * lhsp[j][4];
					rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
					lhsp[j1][2] = lhsp[j1][2] -
							lhsp[j1][1] * lhsp[j][3];
					lhsp[j1][3] = lhsp[j1][3] -
							lhsp[j1][1] * lhsp[j][4];
					rhs[k][j1][i][m] = rhs[k][j1][i][m] -
							lhsp[j1][1] * rhs[k][j][i][m];
					lhsp[j2][1] = lhsp[j2][1] -
							lhsp[j2][0] * lhsp[j][3];
					lhsp[j2][2] = lhsp[j2][2] -
							lhsp[j2][0] * lhsp[j][4];
					rhs[k][j2][i][m] = rhs[k][j2][i][m] -
							lhsp[j2][0] * rhs[k][j][i][m];
					m = 4;
					fac1 = 1.0 / lhsm[j][2];
					lhsm[j][3] = fac1 * lhsm[j][3];
					lhsm[j][4] = fac1 * lhsm[j][4];
					rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
					lhsm[j1][2] = lhsm[j1][2] -
							lhsm[j1][1] * lhsm[j][3];
					lhsm[j1][3] = lhsm[j1][3] -
							lhsm[j1][1] * lhsm[j][4];
					rhs[k][j1][i][m] = rhs[k][j1][i][m] -
							lhsm[j1][1] * rhs[k][j][i][m];
					lhsm[j2][1] = lhsm[j2][1] -
							lhsm[j2][0] * lhsm[j][3];
					lhsm[j2][2] = lhsm[j2][2] -
							lhsm[j2][0] * lhsm[j][4];
					rhs[k][j2][i][m] = rhs[k][j2][i][m] -
							lhsm[j2][0] * rhs[k][j][i][m];
				}

				//---------------------------------------------------------------------
				//         And again the last two rows separately
				//---------------------------------------------------------------------
				j = grid_points[1] - 2;
				j1 = grid_points[1] - 1;
				m = 3;
				fac1 = 1.0 / lhsp[j][2];
				lhsp[j][3] = fac1 * lhsp[j][3];
				lhsp[j][4] = fac1 * lhsp[j][4];
				rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
				lhsp[j1][2] = lhsp[j1][2] -
						lhsp[j1][1] * lhsp[j][3];
				lhsp[j1][3] = lhsp[j1][3] -
						lhsp[j1][1] * lhsp[j][4];
				rhs[k][j1][i][m] = rhs[k][j1][i][m] -
						lhsp[j1][1] * rhs[k][j][i][m];
				m = 4;
				fac1 = 1.0 / lhsm[j][2];
				lhsm[j][3] = fac1 * lhsm[j][3];
				lhsm[j][4] = fac1 * lhsm[j][4];
				rhs[k][j][i][m] = fac1 * rhs[k][j][i][m];
				lhsm[j1][2] = lhsm[j1][2] -
						lhsm[j1][1] * lhsm[j][3];
				lhsm[j1][3] = lhsm[j1][3] -
						lhsm[j1][1] * lhsm[j][4];
				rhs[k][j1][i][m] = rhs[k][j1][i][m] -
						lhsm[j1][1] * rhs[k][j][i][m];
				//---------------------------------------------------------------------
				//               Scale the last row immediately 
				//---------------------------------------------------------------------
				rhs[k][j1][i][3] = rhs[k][j1][i][3] / lhsp[j1][2];
				rhs[k][j1][i][4] = rhs[k][j1][i][4] / lhsm[j1][2];

				//---------------------------------------------------------------------
				//                         BACKSUBSTITUTION 
				//---------------------------------------------------------------------

				j = grid_points[1] - 2;
				j1 = grid_points[1] - 1;
				for (m = 0; m <= 2; m++)
				{
					rhs[k][j][i][m] = rhs[k][j][i][m] -
									 lhs[j][3] * rhs[k][j1][i][m];
				}

				rhs[k][j][i][3] = rhs[k][j][i][3] -
									lhsp[j][3] * rhs[k][j1][i][3];
				rhs[k][j][i][4] = rhs[k][j][i][4] -
									lhsm[j][3] * rhs[k][j1][i][4];

				//---------------------------------------------------------------------
				//      The first three factors
				//---------------------------------------------------------------------
				for (j = grid_points[1] - 3; j >= 0; j--)
				{
					j1 = j + 1;
					j2 = j + 2;
					for (m = 0; m <= 2; m++)
					{
						rhs[k][j][i][m] = rhs[k][j][i][m] -
									 lhs[j][3] * rhs[k][j1][i][m] -
									 lhs[j][4] * rhs[k][j2][i][m];
					}

					//---------------------------------------------------------------------
					//      And the remaining two
					//---------------------------------------------------------------------
					rhs[k][j][i][3] = rhs[k][j][i][3] -
									lhsp[j][3] * rhs[k][j1][i][3] -
									lhsp[j][4] * rhs[k][j2][i][3];
					rhs[k][j][i][4] = rhs[k][j][i][4] -
									lhsm[j][3] * rhs[k][j1][i][4] -
									lhsm[j][4] * rhs[k][j2][i][4];
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
					ru1 = c3c4 * rho_i[k][j][i];
					cv[k] = ws[k][j][i];
					rhos[k] = dmax1(dz4 + con43 * ru1,
									dz5 + c1c5 * ru1,
									dzmax + ru1,
									dz1);
				}

				lhsinit(grid_points[2] - 1);
				for (k = 1; k <= nz2; k++)
				{
					lhs[k][0] = 0.0;
					lhs[k][1] = -dttz2 * cv[k - 1] - dttz1 * rhos[k - 1];
					lhs[k][2] = 1.0 + c2dttz1 * rhos[k];
					lhs[k][3] = dttz2 * cv[k + 1] - dttz1 * rhos[k + 1];
					lhs[k][4] = 0.0;
				}

				//---------------------------------------------------------------------
				//      add fourth order dissipation                                  
				//---------------------------------------------------------------------

				k = 1;
				lhs[k][2] = lhs[k][2] + comz5;
				lhs[k][3] = lhs[k][3] - comz4;
				lhs[k][4] = lhs[k][4] + comz1;

				k = 2;
				lhs[k][1] = lhs[k][1] - comz4;
				lhs[k][2] = lhs[k][2] + comz6;
				lhs[k][3] = lhs[k][3] - comz4;
				lhs[k][4] = lhs[k][4] + comz1;

				for (k = 3; k <= nz2 - 2; k++)
				{
					lhs[k][0] = lhs[k][0] + comz1;
					lhs[k][1] = lhs[k][1] - comz4;
					lhs[k][2] = lhs[k][2] + comz6;
					lhs[k][3] = lhs[k][3] - comz4;
					lhs[k][4] = lhs[k][4] + comz1;
				}

				k = nz2 - 1;
				lhs[k][0] = lhs[k][0] + comz1;
				lhs[k][1] = lhs[k][1] - comz4;
				lhs[k][2] = lhs[k][2] + comz6;
				lhs[k][3] = lhs[k][3] - comz4;

				k = nz2;
				lhs[k][0] = lhs[k][0] + comz1;
				lhs[k][1] = lhs[k][1] - comz4;
				lhs[k][2] = lhs[k][2] + comz5;

				//---------------------------------------------------------------------
				//      subsequently, fill the other factors (u+c), (u-c) 
				//---------------------------------------------------------------------
				for (k = 1; k <= nz2; k++)
				{
					lhsp[k][0] = lhs[k][0];
					lhsp[k][1] = lhs[k][1] - dttz2 * speed[k-1][j][i];
					lhsp[k][2] = lhs[k][2];
					lhsp[k][3] = lhs[k][3] + dttz2 * speed[k+1][j][i];
					lhsp[k][4] = lhs[k][4];

					lhsm[k][0] = lhs[k][0];
					lhsm[k][1] = lhs[k][1] + dttz2 * speed[k-1][j][i];
					lhsm[k][2] = lhs[k][2];
					lhsm[k][3] = lhs[k][3] - dttz2 * speed[k+1][j][i];
					lhsm[k][4] = lhs[k][4];
				}

				//---------------------------------------------------------------------
				// Load a row of K data
				//---------------------------------------------------------------------

				for (k = 0; k <= nz2 + 1; k++)
				{
					rtmp[0 + k * 5] = rhs[k][j][i][0];
					rtmp[1 + k * 5] = rhs[k][j][i][1];
					rtmp[2 + k * 5] = rhs[k][j][i][2];
					rtmp[3 + k * 5] = rhs[k][j][i][3];
					rtmp[4 + k * 5] = rhs[k][j][i][4];
				}

				//---------------------------------------------------------------------
				//                          FORWARD ELIMINATION  
				//---------------------------------------------------------------------

				for (k = 0; k <= grid_points[2] - 3; k++)
				{
					k1 = k + 1;
					k2 = k + 2;
					fac1 = 1.0 / lhs[k][2];
					lhs[k][3] = fac1 * lhs[k][3];
					lhs[k][4] = fac1 * lhs[k][4];
					for (m = 0; m <= 2; m++)
					{
						rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
					}
					lhs[k1][2] = lhs[k1][2] -
								   lhs[k1][1] * lhs[k][3];
					lhs[k1][3] = lhs[k1][3] -
								   lhs[k1][1] * lhs[k][4];
					for (m = 0; m <= 2; m++)
					{
						rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
									lhs[k1][1] * rtmp[m + k * 5];
					}
					lhs[k2][1] = lhs[k2][1] -
								   lhs[k2][0] * lhs[k][3];
					lhs[k2][2] = lhs[k2][2] -
								   lhs[k2][0] * lhs[k][4];
					for (m = 0; m <= 2; m++)
					{
						rtmp[m + k2 * 5] = rtmp[m + k2 * 5] -
									lhs[k2][0] * rtmp[m + k * 5];
					}
				}

				//---------------------------------------------------------------------
				//      The last two rows in this grid block are a bit different, 
				//      since they do not have two more rows available for the
				//      elimination of off-diagonal entries
				//---------------------------------------------------------------------
				k = grid_points[2] - 2;
				k1 = grid_points[2] - 1;
				fac1 = 1.0 / lhs[k][2];
				lhs[k][3] = fac1 * lhs[k][3];
				lhs[k][4] = fac1 * lhs[k][4];
				for (m = 0; m <= 2; m++)
				{
					rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
				}
				lhs[k1][2] = lhs[k1][2] -
							   lhs[k1][1] * lhs[k][3];
				lhs[k1][3] = lhs[k1][3] -
							   lhs[k1][1] * lhs[k][4];
				for (m = 0; m <= 2; m++)
				{
					rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
								lhs[k1][1] * rtmp[m + k * 5];
				}
				//---------------------------------------------------------------------
				//               scale the last row immediately
				//---------------------------------------------------------------------
				fac2 = 1.0 / lhs[k1][2];
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
					fac1 = 1.0 / lhsp[k][2];
					lhsp[k][3] = fac1 * lhsp[k][3];
					lhsp[k][4] = fac1 * lhsp[k][4];
					rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
					lhsp[k1][2] = lhsp[k1][2] -
							lhsp[k1][1] * lhsp[k][3];
					lhsp[k1][3] = lhsp[k1][3] -
							lhsp[k1][1] * lhsp[k][4];
					rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
							lhsp[k1][1] * rtmp[m + k * 5];
					lhsp[k2][1] = lhsp[k2][1] -
							lhsp[k2][0] * lhsp[k][3];
					lhsp[k2][2] = lhsp[k2][2] -
							lhsp[k2][0] * lhsp[k][4];
					rtmp[m + k2 * 5] = rtmp[m + k2 * 5] -
							lhsp[k2][0] * rtmp[m + k * 5];
					m = 4;
					fac1 = 1.0 / lhsm[k][2];
					lhsm[k][3] = fac1 * lhsm[k][3];
					lhsm[k][4] = fac1 * lhsm[k][4];
					rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
					lhsm[k1][2] = lhsm[k1][2] -
							lhsm[k1][1] * lhsm[k][3];
					lhsm[k1][3] = lhsm[k1][3] -
							lhsm[k1][1] * lhsm[k][4];
					rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
							lhsm[k1][1] * rtmp[m + k * 5];
					lhsm[k2][1] = lhsm[k2][1] -
							lhsm[k2][0] * lhsm[k][3];
					lhsm[k1][2] = lhsm[k1][2] -
							lhsm[k2][0] * lhsm[k][4];
					rtmp[m + k2 * 5] = rtmp[m + k2 * 5] -
							lhsm[k2][0] * rtmp[m + k * 5];
				}

				//---------------------------------------------------------------------
				//         And again the last two rows separately
				//---------------------------------------------------------------------
				k = grid_points[2] - 2;
				k1 = grid_points[2] - 1;
				m = 3;
				fac1 = 1.0 / lhsp[k][2];
				lhsp[k][3] = fac1 * lhsp[k][3];
				lhsp[k][4] = fac1 * lhsp[k][4];
				rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
				lhsp[k1][2] = lhsp[k1][2] -
						lhsp[k1][1] * lhsp[k][3];
				lhsp[k1][3] = lhsp[k1][3] -
						lhsp[k1][1] * lhsp[k][4];
				rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
						lhsp[k1][1] * rtmp[m + k * 5];
				m = 4;
				fac1 = 1.0 / lhsm[k][2];
				lhsm[k][3] = fac1 * lhsm[k][3];
				lhsm[k][4] = fac1 * lhsm[k][4];
				rtmp[m + k * 5] = fac1 * rtmp[m + k * 5];
				lhsm[k1][2] = lhsm[k1][2] -
						lhsm[k1][1] * lhsm[k][3];
				lhsm[k1][3] = lhsm[k1][3] -
						lhsm[k1][1] * lhsm[k][4];
				rtmp[m + k1 * 5] = rtmp[m + k1 * 5] -
						lhsm[k1][1] * rtmp[m + k * 5];
				//---------------------------------------------------------------------
				//               Scale the last row immediately (some of this is overkill
				//               if this is the last cell)
				//---------------------------------------------------------------------
				rtmp[3 + k1 * 5] = rtmp[3 + k1 * 5] / lhsp[k1][2];
				rtmp[4 + k1 * 5] = rtmp[4 + k1 * 5] / lhsm[k1][2];

				//---------------------------------------------------------------------
				//                         BACKSUBSTITUTION 
				//---------------------------------------------------------------------

				k = grid_points[2] - 2;
				k1 = grid_points[2] - 1;
				for (m = 0; m <= 2; m++)
				{
					rtmp[m + k * 5] = rtmp[m + k * 5] -
									   lhs[k][3] * rtmp[m + k1 * 5];
				}

				rtmp[3 + k * 5] = rtmp[3 + k * 5] -
									  lhsp[k][3] * rtmp[3 + k1 * 5];
				rtmp[4 + k * 5] = rtmp[4 + k * 5] -
									  lhsm[k][3] * rtmp[4 + k1 * 5];

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
									 lhs[k][3] * rtmp[m + k1 * 5] -
									 lhs[k][4] * rtmp[m + k2 * 5];
					}

					//---------------------------------------------------------------------
					//      And the remaining two
					//---------------------------------------------------------------------
					rtmp[3 + k * 5] = rtmp[3 + k * 5] -
									lhsp[k][3] * rtmp[3 + k1 * 5] -
									lhsp[k][4] * rtmp[3 + k2 * 5];
					rtmp[4 + k * 5] = rtmp[4 + k * 5] -
									lhsm[k][3] * rtmp[4 + k1 * 5] -
									lhsm[k][4] * rtmp[4 + k2 * 5];
				}

				//---------------------------------------------------------------------
				//      Store result
				//---------------------------------------------------------------------
				for (k = 0; k <= nz2 + 1; k++)
				{
					rhs[k][j][i][0] = rtmp[0 + k * 5];
					rhs[k][j][i][1] = rtmp[1 + k * 5];
					rhs[k][j][i][2] = rtmp[2 + k * 5];
					rhs[k][j][i][3] = rtmp[3 + k * 5];
					rhs[k][j][i][4] = rtmp[4 + k * 5];
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




