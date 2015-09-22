// LU.cs created with MonoDevelop
// User: diane at 12:56 10/4/2009
//
// To change standard headers go to Edit->Preferences->Coding->Standard Headers
//
/*
!-------------------------------------------------------------------------!
!								          !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0         !
!								          !
!			C# 	V E R S I O N			  !
!									  !
!                                  LU                                     !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is a serial/multithreaded version of                  !
!    the NPB3_0_JAV LU code.                                              !
!                                                                         !
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
! Authors: R. Van der Wijngaart 					  !
!	   T. Harris							  !
!	   M. Yarrow							  !
! Translation to Java and MultiThreaded Code				  !
!	   M. Frumkin							  !
!	   M. Schultz							  !
! Translation to C#				  !
!      Elidiane Martins
!-------------------------------------------------------------------------!
*/


using NPB3_0_JAV.LUThreads;
using NPB3_0_JAV.BMInOut;
using System;
using System.IO;

namespace NPB3_0_JAV{
	

	public class LU : LUBase{
		public int bid = -1;
		public BMResults results;
        
		public LU(char clss, int np) : base(clss, np){
		}
		
		public static void Main(String[] argv){
			LU lu = null;

			BMArgs.ParseCmdLineArgs(argv, BMName);
			char CLSS = BMArgs.CLASS;
			int np = BMArgs.num_threads;
			
			try{ 
				lu = new LU(CLSS, np);
			} catch(OutOfMemoryException e){
				Console.WriteLine(e.Message);
				BMArgs.outOfMemoryMessage();
				Environment.Exit(0);
			}      
			lu.runBenchMark();
		}
  
		public void run() {runBenchMark();}
  
		public void runBenchMark(){
			BMArgs.Banner(BMName, CLASS, true, num_threads);
    
			int numTimers = t_last + 1;
			String[] t_names = new String[numTimers];
			double[] trecs = new double[numTimers];
			setTimers(t_names);

			getInputPars();

//---------------------------------------------------------------------
//   set up domain sizes
//---------------------------------------------------------------------
    
			domain();

//---------------------------------------------------------------------
//   set up coefficients
//---------------------------------------------------------------------
    
			setcoeff();
 //   if(!serial) setupThreads(this);
//---------------------------------------------------------------------
//   set the boundary values for dependent variables
//---------------------------------------------------------------------
    
			setbv();
  
//---------------------------------------------------------------------
//   set the initial values for dependent variables
//---------------------------------------------------------------------
   
			setiv();

//---------------------------------------------------------------------
//   compute the forcing term based on prescribed exact solution
//---------------------------------------------------------------------
    
			erhs();

//---------------------------------------------------------------------
//   perform the SSOR iterations
//---------------------------------------------------------------------
    
			double tm;
			tm = sssor();

//---------------------------------------------------------------------
//   compute the solution error
//---------------------------------------------------------------------
    
			error();

//---------------------------------------------------------------------
//   compute the surface integral
//---------------------------------------------------------------------
    
			pintgr();

//---------------------------------------------------------------------
//   verification test
//---------------------------------------------------------------------
    
			int verified = verify(rsdnm, errnm, frc);
			results = new BMResults(BMName,
			                        CLASS,
			                        nx0,
			                        ny0,
			                        nz0,
			                        itmax,
			                        tm,
			                        getMFLOPS(itmax, tm),
			                        "floating point",
			                        verified,
			                        true,
			                        num_threads,
			                        bid);
			results.print();				
//---------------------------------------------------------------------
//      More timers
//---------------------------------------------------------------------
    
			if(timeron) printTimers(t_names, trecs, tm);
		}
		public double getMFLOPS(int itmax, double tm){
			double mflops = 0.0;
			if( tm > 0 ){
				mflops = 1984.77 * nx0 * ny0 * nz0
					-10923.3 * pow2((nx0 + ny0 + nz0)/ 3.0) 
						+27770.9 * (nx0 + ny0 + nz0) / 3.0-144010.0;
				mflops *= itmax / (tm * 1000000.0);
			}
			return mflops;
		}
		public void printTimers(String[] t_names, double[] trecs, double tm){
			//DecimalFormat fmt = new DecimalFormat("0.000");
			Console.WriteLine("  SECTION     Time (secs)");
			for(int i = 0; i < t_last; i++) trecs[i] = timer.readTimer(i);
			if ( tm == 0.0 ) tm = 1.0;
			for(int i = 1; i < t_last; i++){
				double dbl =(trecs[i] * 100.0 / tm);
				Console.WriteLine("  " + t_names[i] + ":" + trecs[i].ToString("N3") + 
				                  "(" + dbl.ToString("N3") + "%)" );
				if (i == t_rhs) {
					double t = trecs[t_rhsx] + trecs[t_rhsy] + trecs[t_rhsz];
					dbl = (t * 100.0 / tm);
					Console.WriteLine("     " + "--> total " + "sub-rhs" + 
					                  ":" + t.ToString("N3") +
					                  "  (" + dbl.ToString("N3") + "%)");
					t = trecs[i] - t;
					dbl = (t * 100.0 / tm);
					Console.WriteLine("     " + "--> total " + "rest-rhs" +
					                  ":" + t.ToString("N3") + 
					                  "  (" + dbl.ToString("N3") + "%)");
				}
			}
		}
  
		public void setTimers(String[] t_names){
			//File f1 = new File("timer.flag");
			timeron = false;
			if(File.Exists("timer.flag")){
				timeron = true;
				t_names[t_total] = "total";
				t_names[t_rhsx] = "rhsx";
				t_names[t_rhsy] = "rhsy";
				t_names[t_rhsz] = "rhsz";
				t_names[t_rhs] = "rhs";
				t_names[t_jacld] = "jacld";
				t_names[t_blts] = "blts";
				t_names[t_jacu] = "jacu";
				t_names[t_buts] = "buts";
				t_names[t_add] = "add";
				t_names[t_l2norm] = "l2norm";
			}
		}

		public void blts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, 
		                 int k, double omega, 
		                 double[][][][] v, double[][][] tv, double[][][][] ldz, double[][][][] ldy, 
		                 double[][][][] ldx, double[][][][] d, 
		                 int ist, int iend, int jst, int jend, int nx0, int ny0){
			int i, j, m;
			double  tmp, tmp1;
			double[][]  tmat = instantiate_jagged_array_2(5, 5);

			for(j = jst-1; j <= jend-1; j++){
				for(i = ist-1; i <= iend-1; i++){
					for(m = 0; m <= 4; m++){

						tv[m][i][j] =  v[m][i][j][k]
						- omega * (  ldz[ m][ 0][ i][j] * v[ 0][ i][ j][k-1]
						           + ldz[ m][ 1][ i][j] * v[ 1][ i][ j][k-1]
						           + ldz[ m][ 2][ i][j] * v[ 2][ i][ j][k-1]
						           + ldz[ m][ 3][ i][j] * v[ 3][ i][ j][k-1]
						           + ldz[ m][ 4][ i][j] * v[ 4][ i][ j][k-1]  );
					}
				}
			}

			for(j=jst-1;j<=jend-1;j++){
				for(i=ist-1;i<=iend-1;i++){
					for(m=0;m<=4;m++){

						tv[ m][ i][j] =  tv[ m][ i ][j]
						- omega * (  ldy[ m][ 0][ i][j] * v[ 0][ i][ j-1][k]
						           + ldx[ m][ 0][ i][j] * v[ 0][ i-1][ j][k]
						           + ldy[ m][ 1][ i][j] * v[ 1][ i][ j-1][k]
						           + ldx[ m][ 1][ i][j] * v[ 1][ i-1][ j][k]
						           + ldy[ m][ 2][ i][j] * v[ 2][ i][ j-1][k]
						           + ldx[ m][ 2][ i][j] * v[ 2][ i-1][ j][k]
						           + ldy[ m][ 3][ i][j] * v[ 3][ i][ j-1][k]
						           + ldx[ m][ 3][ i][j] * v[ 3][ i-1][ j][k]
						           + ldy[ m][ 4][ i][j] * v[ 4][ i][ j-1][k]
						           + ldx[ m][ 4][ i][j] * v[ 4][ i-1][ j][k] );
					}
       
//---------------------------------------------------------------------
//   diagonal block inversion
//   forward elimination
//---------------------------------------------------------------------
            
					for(m=0;m<=4;m++){
						tmat[ m][0] = d[ m][ 0][ i][j];
						tmat[ m][1] = d[ m][ 1][ i][j];
						tmat[ m][2] = d[ m][ 2][ i][j];
						tmat[ m][3] = d[ m][ 3][ i][j];
						tmat[ m][4] = d[ m][ 4][ i][j];
					}

					tmp1 = 1.0 / tmat[ 0][0];
					tmp = tmp1 * tmat[ 1][0];
					tmat[ 1][1] =  tmat[ 1][1]
					- tmp * tmat[ 0][1];
					tmat[ 1][2] =  tmat[ 1][2]
					- tmp * tmat[ 0][2];
					tmat[ 1][3] =  tmat[ 1][3]
					- tmp * tmat[ 0][3];
					tmat[ 1][4] =  tmat[ 1][4]
					- tmp * tmat[ 0][4];
					tv[ 1][ i][j] = tv[ 1][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp = tmp1 * tmat[ 2][0];
					tmat[ 2][1] =  tmat[ 2][1]
					- tmp * tmat[ 0][1];
					tmat[ 2][2] =  tmat[ 2][2]
					- tmp * tmat[ 0][2];
					tmat[ 2][3] =  tmat[ 2][3]
					- tmp * tmat[ 0][3];
					tmat[ 2][4] =  tmat[ 2][4]
					- tmp * tmat[ 0][4];
					tv[ 2][ i][j] = tv[ 2][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp = tmp1 * tmat[ 3][0];
					tmat[ 3][1] =  tmat[ 3][1]
					- tmp * tmat[ 0][1];
					tmat[ 3][2] =  tmat[ 3][2]
					- tmp * tmat[ 0][2];
					tmat[ 3][3] =  tmat[ 3][3]
					- tmp * tmat[ 0][3];
					tmat[ 3][4] =  tmat[ 3][4]
					- tmp * tmat[ 0][4];
					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp = tmp1 * tmat[ 4][0];
					tmat[ 4][1] =  tmat[ 4][1]
					- tmp * tmat[ 0][1];
					tmat[ 4][2] =  tmat[ 4][2]
					- tmp * tmat[ 0][2];
					tmat[ 4][3] =  tmat[ 4][3]
					- tmp * tmat[ 0][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 0][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp1 = 1.0 / tmat[ 1][1];
					tmp = tmp1 * tmat[ 2][1];
					tmat[ 2][2] =  tmat[ 2][2]
					- tmp * tmat[ 1][2];
					tmat[ 2][3] =  tmat[ 2][3]
					- tmp * tmat[ 1][3];
					tmat[ 2][4] =  tmat[ 2][4]
					- tmp * tmat[ 1][4];
					tv[ 2][ i][j] = tv[ 2][ i][j]
					- tv[ 1][ i][j] * tmp;

					tmp = tmp1 * tmat[ 3][1];
					tmat[ 3][2] =  tmat[ 3][2]
					- tmp * tmat[ 1][2];
					tmat[ 3][3] =  tmat[ 3][3]
					- tmp * tmat[ 1][3];
					tmat[ 3][4] =  tmat[ 3][4]
					- tmp * tmat[ 1][4];
					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tv[ 1][ i][j] * tmp;

					tmp = tmp1 * tmat[ 4][1];
					tmat[ 4][2] =  tmat[ 4][2]
					- tmp * tmat[ 1][2];
					tmat[ 4][3] =  tmat[ 4][3]
					- tmp * tmat[ 1][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 1][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 1][ i][j] * tmp;

					tmp1 = 1.0 / tmat[ 2][2];
					tmp = tmp1 * tmat[ 3][2];
					tmat[ 3][3] =  tmat[ 3][3]
					- tmp * tmat[ 2][3];
					tmat[ 3][4] =  tmat[ 3][4]
					- tmp * tmat[ 2][4];
					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tv[ 2][ i][j] * tmp;

					tmp = tmp1 * tmat[ 4][2];
					tmat[ 4][3] =  tmat[ 4][3]
					- tmp * tmat[ 2][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 2][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 2][ i][j] * tmp;

					tmp1 = 1.0 / tmat[ 3][3];
					tmp = tmp1 * tmat[ 4][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 3][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 3][ i][j] * tmp;

//---------------------------------------------------------------------
//   back substitution
//---------------------------------------------------------------------
            
					v[ 4][ i][ j][k] = tv[ 4][ i][j]
					/ tmat[ 4][4];

					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tmat[ 3][4] * v[ 4][ i][ j][k];
					v[ 3][ i][ j][k] = tv[ 3][ i][j]
					/ tmat[ 3][3];

					tv[ 2][ i][j] = tv[ 2][ i][j]
					- tmat[ 2][3] * v[ 3][ i][ j][k]
					- tmat[ 2][4] * v[ 4][ i][ j][k];
					v[ 2][ i][ j][k] = tv[ 2][ i][j]
					/ tmat[ 2][2];

					tv[ 1][ i][j] = tv[ 1][ i][j]
					- tmat[ 1][2] * v[ 2][ i][ j][k]
					- tmat[ 1][3] * v[ 3][ i][ j][k]
					- tmat[ 1][4] * v[ 4][ i][ j][k];
					v[ 1][ i][ j][k] = tv[ 1][ i][j]
					/ tmat[ 1][1];

					tv[ 0][ i][j] = tv[ 0][ i][j]
					- tmat[ 0][1] * v[ 1][ i][ j][k]
					- tmat[ 0][2] * v[ 2][ i][ j][k]
					- tmat[ 0][3] * v[ 3][ i][ j][k]
					- tmat[ 0][4] * v[ 4][ i][ j][k];
					v[ 0][ i][ j][k] = tv[ 0][ i][j]
					/ tmat[ 0][0];
				}
			}
		}

		public void buts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, 
		                 int k, double omega, 
		                 double[][][][] v, double[][][] tv, double[][][][] d, double[][][][] udx, 
		                 double[][][][] udy, double[][][][] udz, 
		                 int ist, int iend, int jst, int jend, int nx0, int ny0){
			int i, j, m;
			double  tmp, tmp1;
			double[][]  tmat =  instantiate_jagged_array_2(5, 5);


			for(j=jend-1;j>=jst-1;j--){
				for(i=iend-1;i>=ist-1;i--){
					for(m=0;m<=4;m++){
						tv[ m][ i][j] = 
							omega * (  udz[ m][ 0][ i][j] * v[ 0][ i][ j][k+1]
							         + udz[ m][ 1][ i][j] * v[ 1][ i][ j][k+1]
							         + udz[ m][ 2][ i][j] * v[ 2][ i][ j][k+1]
							         + udz[ m][ 3][ i][j] * v[ 3][ i][ j][k+1]
							         + udz[ m][ 4][ i][j] * v[ 4][ i][ j][k+1] );
					}
				}
			}

      
			for(j=jend-1;j>=jst-1;j--){
				for(i=iend-1;i>=ist-1;i--){
					for(m=0;m<=4;m++){
						tv[ m][ i][j] = tv[ m][ i][j]
						+ omega * ( udy[ m][ 0][ i][j] * v[ 0][ i][ j+1][k]
						           + udx[ m][ 0][ i][j] * v[ 0][ i+1][ j][k]
						           + udy[ m][ 1][ i][j] * v[ 1][ i][ j+1][k]
						           + udx[ m][ 1][ i][j] * v[ 1][ i+1][ j][k]
						           + udy[ m][ 2][ i][j] * v[ 2][ i][ j+1][k]
						           + udx[ m][ 2][ i][j] * v[ 2][ i+1][ j][k]
						           + udy[ m][ 3][ i][j] * v[ 3][ i][ j+1][k]
						           + udx[ m][ 3][ i][j] * v[ 3][ i+1][ j][k]
						           + udy[ m][ 4][ i][j] * v[ 4][ i][ j+1][k]
						           + udx[ m][ 4][ i][j] * v[ 4][ i+1][ j][k] );
					}

//---------------------------------------------------------------------
//   diagonal block inversion
//---------------------------------------------------------------------
            
					for(m=0;m<=4;m++){
						tmat[ m][0] = d[ m][ 0][ i][j];
						tmat[ m][1] = d[ m][ 1][ i][j];
						tmat[ m][2] = d[ m][ 2][ i][j];
						tmat[ m][3] = d[ m][ 3][ i][j];
						tmat[ m][4] = d[ m][ 4][ i][j];
					}

					tmp1 = 1.0 / tmat[ 0][0];
					tmp = tmp1 * tmat[ 1][0];
					tmat[ 1][1] =  tmat[ 1][1]
					- tmp * tmat[ 0][1];
					tmat[ 1][2] =  tmat[ 1][2]
					- tmp * tmat[ 0][2];
					tmat[ 1][3] =  tmat[ 1][3]
					- tmp * tmat[ 0][3];
					tmat[ 1][4] =  tmat[ 1][4]
					- tmp * tmat[ 0][4];
					tv[ 1][ i][j] = tv[ 1][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp = tmp1 * tmat[ 2][0];
					tmat[ 2][1] =  tmat[ 2][1]
					- tmp * tmat[ 0][1];
					tmat[ 2][2] =  tmat[ 2][2]
					- tmp * tmat[ 0][2];
					tmat[ 2][3] =  tmat[ 2][3]
					- tmp * tmat[ 0][3];
					tmat[ 2][4] =  tmat[ 2][4]
					- tmp * tmat[ 0][4];
					tv[ 2][ i][j] = tv[ 2][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp = tmp1 * tmat[ 3][0];
					tmat[ 3][1] =  tmat[ 3][1]
					- tmp * tmat[ 0][1];
					tmat[ 3][2] =  tmat[ 3][2]
					- tmp * tmat[ 0][2];
					tmat[ 3][3] =  tmat[ 3][3]
					- tmp * tmat[ 0][3];
					tmat[ 3][4] =  tmat[ 3][4]
					- tmp * tmat[ 0][4];
					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp = tmp1 * tmat[ 4][0];
					tmat[ 4][1] =  tmat[ 4][1]
					- tmp * tmat[ 0][1];
					tmat[ 4][2] =  tmat[ 4][2]
					- tmp * tmat[ 0][2];
					tmat[ 4][3] =  tmat[ 4][3]
					- tmp * tmat[ 0][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 0][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 0][ i][j] * tmp;

					tmp1 = 1.0 / tmat[ 1][1];
					tmp = tmp1 * tmat[ 2][1];
					tmat[ 2][2] =  tmat[ 2][2]
					- tmp * tmat[ 1][2];
					tmat[ 2][3] =  tmat[ 2][3]
					- tmp * tmat[ 1][3];
					tmat[ 2][4] =  tmat[ 2][4]
					- tmp * tmat[ 1][4];
					tv[ 2][ i][j] = tv[ 2][ i][j]
					- tv[ 1][ i][j] * tmp;

					tmp = tmp1 * tmat[ 3][1];
					tmat[ 3][2] =  tmat[ 3][2]
					- tmp * tmat[ 1][2];
					tmat[ 3][3] =  tmat[ 3][3]
					- tmp * tmat[ 3][3];
					tmat[ 3][4] =  tmat[ 3][4]
					- tmp * tmat[ 1][4];
					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tv[ 1][ i][j] * tmp;

					tmp = tmp1 * tmat[ 4][1];
					tmat[ 4][2] =  tmat[ 4][2]
					- tmp * tmat[ 1][2];
					tmat[ 4][3] =  tmat[ 4][3]
					- tmp * tmat[ 1][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 1][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 1][ i][j] * tmp;

					tmp1 = 1.0 / tmat[ 2][2];
					tmp = tmp1 * tmat[ 3][2];
					tmat[ 3][3] =  tmat[ 3][3]
					- tmp * tmat[ 2][3];
					tmat[ 3][4] =  tmat[ 3][4]
					- tmp * tmat[ 2][4];
					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tv[ 2][ i][j] * tmp;

					tmp = tmp1 * tmat[ 4][2];
					tmat[ 4][3] =  tmat[ 4][3]
					- tmp * tmat[ 2][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 2][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 2][ i][j] * tmp;

					tmp1 = 1.0 / tmat[ 3][3];
					tmp = tmp1 * tmat[ 4][3];
					tmat[ 4][4] =  tmat[ 4][4]
					- tmp * tmat[ 3][4];
					tv[ 4][ i][j] = tv[ 4][ i][j]
					- tv[ 3][ i][j] * tmp;

//---------------------------------------------------------------------
//   back substitution
//---------------------------------------------------------------------
					tv[ 4][ i][j] = tv[ 4][ i][j]
					/ tmat[ 4][4];

					tv[ 3][ i][j] = tv[ 3][ i][j]
					- tmat[ 3][4] * tv[ 4][ i][j];
					tv[ 3][ i][j] = tv[ 3][ i][j]
					/ tmat[ 3][3];

					tv[ 2][ i][j] = tv[ 2][ i][j]
					- tmat[ 2][3] * tv[ 3][ i][j]
					- tmat[ 2][4] * tv[ 4][ i][j];
					tv[ 2][ i][j] = tv[ 2][ i][j]
					/ tmat[ 2][2];

					tv[ 1][ i][j] = tv[ 1][ i][j]
					- tmat[ 1][2] * tv[ 2][ i][j]
					- tmat[ 1][3] * tv[ 3][ i][j]
					- tmat[ 1][4] * tv[ 4][ i][j];
					tv[ 1][ i][j] = tv[ 1][ i][j]
					/ tmat[ 1][1];

					tv[ 0][ i][j] = tv[ 0][ i][j]
					- tmat[ 0][1] * tv[ 1][ i][j]
					- tmat[ 0][2] * tv[ 2][ i][j]
					- tmat[ 0][3] * tv[ 3][ i][j]
					- tmat[ 0][4] * tv[ 4][ i][j];
					tv[ 0][ i][j] = tv[ 0][ i][j]
					/ tmat[ 0][0];

					v[ 0][ i][j ][k] = v[ 0][ i][j ][k] 
					- tv[ 0][ i][j];
					v[ 1][ i][ j][k] = v[ 1][ i][ j][k] 
					- tv[ 1][ i][j];
					v[ 2][ i][ j][k] = v[ 2][ i][ j][k] 
					- tv[ 2][ i][j];
					v[ 3][ i][ j][k] = v[ 3][ i][ j][k] 
					- tv[ 3][ i][j];
					v[ 4][ i][ j][k] = v[ 4][ i][ j][k] 
					- tv[ 4][ i][j];	    
				}
			}
		}

  
		public void domain(){
			nx = nx0;
			ny = ny0;
			nz = nz0;

//---------------------------------------------------------------------
//   check the sub-domain size
//---------------------------------------------------------------------
			if ( ( nx < 4 ) || ( ny < 4 ) || ( nz < 4 ) ) {
				Console.WriteLine("     " + "SUBDOMAIN SIZE IS TOO SMALL - ");
				Console.WriteLine("     " + "ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS");
				Console.WriteLine("     " + "SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL");
				Console.WriteLine("     " + "TO 4 THEY ARE CURRENTLY "+ nx + " "+ny + " " + nz);
				Environment.Exit(0);
			}

			if ( ( nx > isiz1 ) || ( ny > isiz2 ) || ( nz > isiz3 ) ) {
				Console.WriteLine("     " + "SUBDOMAIN SIZE IS TOO LARGE - ");
				Console.WriteLine("     " + "ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS");
				Console.WriteLine("     " + "SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL");
				Console.WriteLine( "     " + "TO ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY."
				                  +" THEY ARE");
				Console.WriteLine("     " + " CURRENTLY "+ nx + " "+ny + " " + nz);
				Environment.Exit(0);
			}

//---------------------------------------------------------------------
//   set up the start and end in i and j extents for all processors
//---------------------------------------------------------------------
      
			ist = 2;
			iend = nx - 1;

			jst = 2;
			jend = ny - 1;

			ii1 = 2;
			ii2 = nx0 - 1;
			ji1 = 2;
			ji2 = ny0 - 2;
			ki1 = 3;
			ki2 = nz0 - 1;
		}

  public void erhs(){
        int i, j, k, m;
      double  xi, eta, zeta;
      double  q;
      double  u21, u31, u41;
      double  tmp;
      double  u21i, u31i, u41i, u51i;
      double  u21j, u31j, u41j, u51j;
      double  u21k, u31k, u41k, u51k;
      double  u21im1, u31im1, u41im1, u51im1;
      double  u21jm1, u31jm1, u41jm1, u51jm1;
      double  u21km1, u31km1, u41km1, u51km1;

      for(k=0;k<=nz-1;k++){
         for(j=0;j<=ny-1;j++){
            for(i=0;i<=nx-1;i++){
               for(m=0;m<=4;m++){
                  frct[ m ][ i][ j][k] = 0.0;
               }
            }
         }
      }

      for(k=0;k<=nz-1;k++){
         zeta = (double) k/(nz-1);
         for(j=0;j<=ny-1;j++){
            eta = (double) j/(ny0-1);
            for(i=0;i<=nx-1;i++){
               xi = (double) i/(nx0-1);
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] =  ce[0][m]
                       + (ce[1][m]
                       + (ce[4][m]
                       + (ce[7][m]
                       +  ce[10][m] * xi) * xi) * xi) * xi
                       + (ce[2][m]
                       + (ce[5][m]
                       + (ce[8][m]
                       +  ce[11][m] * eta) * eta) * eta) * eta
                       + (ce[3][m]
                       + (ce[6][m]
                       + (ce[9][m]
                       +  ce[12][m] * zeta) * zeta) * zeta) * zeta;
              }
            }
         }
      }

//---------------------------------------------------------------------
//   xi-direction flux differences
//---------------------------------------------------------------------

      for(k=1;k<=nz - 2;k++){
         for(j=jst-1;j<=jend-1;j++){
            for(i=0;i<=nx-1;i++){
               flux[0][i] = rsd[ 1][ i][ j][k];
               u21 = rsd[ 1][ i][ j][k] / rsd[ 0][ i][ j][k];
               q = 0.50 * (  rsd[ 1][ i][ j][k] * rsd[ 1][ i][ j][k]
                               + rsd[ 2][ i][ j][k] * rsd[ 2][ i][ j][k]
                               + rsd[ 3][ i][ j][k] * rsd[ 3][ i][ j][k] )
                            / rsd[ 0][ i][ j][k];
               flux[ 1][i] = rsd[ 1][ i][ j][k] * u21 + c2 * 
                               ( rsd[ 4][ i][ j][k] - q );
               flux[ 2][i] = rsd[ 2][ i][ j][k] * u21;
               flux[ 3][i] = rsd[ 3][ i][ j][k] * u21;
               flux[ 4][i] = ( c1 * rsd[ 4][ i][ j][k] - c2 * q ) * u21;
            }

            for(i=ist-1;i<=iend-1;i++){
               for(m=0;m<=4;m++){
                  frct[ m][ i][ j][k] =  frct[ m][ i][ j][k]
                         - tx2 * ( flux[ m][i+1] - flux[ m][i-1] );
               }
            }
            for(i=ist-1;i<=nx-1;i++){
               tmp = 1.0 / rsd[ 0][ i][ j][k];

               u21i = tmp * rsd[ 1][ i][ j][k];
               u31i = tmp * rsd[ 2][ i][ j][k];
               u41i = tmp * rsd[ 3][ i][ j][k];
               u51i = tmp * rsd[ 4][ i][ j][k];

               tmp = 1.0 / rsd[0][ i-1][ j][k];

               u21im1 = tmp * rsd[ 1][ i-1][ j][k];
               u31im1 = tmp * rsd[ 2][ i-1][ j][k];
               u41im1 = tmp * rsd[ 3][ i-1][ j][k];
               u51im1 = tmp * rsd[ 4][ i-1][ j][k];

               flux[1][i] = (4.0/3.0) * tx3 * 
                              ( u21i - u21im1 );
               flux[2][i] = tx3 * ( u31i - u31im1 );
               flux[3][i] = tx3 * ( u41i - u41im1 );
               flux[4][i] = 0.50 * ( 1.0 - c1*c5 )
                    * tx3 * ( ( pow2(u21i) + pow2(u31i) + pow2(u41i) )
                            - ( pow2(u21im1) + pow2(u31im1) + pow2(u41im1) ) )
                    + (1.0/6.0)
                    * tx3 * ( pow2(u21i) - pow2(u21im1) )
                    + c1 * c5 * tx3 * ( u51i - u51im1 );
            }

            for(i=ist-1;i<=iend-1;i++){
               frct[ 0][ i][ j][k] = frct[ 0][ i][ j][k]
                    + dx1 * tx1 * (            rsd[ 0][ i-1][ j][k]
                                   - 2.0 * rsd[ 0][ i][ j][k]
                                   +           rsd[ 0][ i+1][ j][k] );
               frct[ 1][ i][ j][k] = frct[ 1][ i][ j][k]
                 + tx3 * c3 * c4 * ( flux[ 1][i+1] - flux[ 1][i] )
                    + dx2 * tx1 * (            rsd[ 1][ i-1][ j][k]
                                   - 2.0 * rsd[ 1][ i][ j][k]
                                   +           rsd[ 1][ i+1][ j][k] );
               frct[ 2][ i][ j][k] = frct[ 2][ i][ j][k]
                 + tx3 * c3 * c4 * ( flux[ 2][i+1] - flux[ 2][i] )
                    + dx3 * tx1 * (            rsd[ 2][ i-1][ j][k]
                                   - 2.0 * rsd[ 2][ i][ j][k]
                                   +           rsd[ 2][ i+1][ j][k] );
               frct[ 3][ i][ j][k] = frct[ 3][ i][ j][k]
                  + tx3 * c3 * c4 * ( flux[ 3][i+1] - flux[ 3][i] )
                    + dx4 * tx1 * (            rsd[ 3][ i-1][ j][k]
                                   - 2.0 * rsd[ 3][ i][ j][k]
                                   +           rsd[ 3][ i+1][ j][k] );
               frct[ 4][ i][ j][k] = frct[ 4][ i][ j][k]
                 + tx3 * c3 * c4 * ( flux[ 4][i+1] - flux[ 4][i] )
                    + dx5 * tx1 * (            rsd[ 4][ i-1][ j][k]
                                   - 2.0 * rsd[ 4][ i][ j][k]
                                   +           rsd[ 4][ i+1][ j][k] );
            }
					
					
					

//---------------------------------------------------------------------
//   Fourth-order dissipation
//---------------------------------------------------------------------
            for(m=0;m<=4;m++){
               frct[ m][ 1][ j][k] = frct[ m][ 1][ j][k]
                 - dssp * ( + 5.0 * rsd[ m][ 1][ j][k]
                             - 4.0 * rsd[ m][ 2][ j][k]
                             +           rsd[ m][ 3][ j][k] );
               frct[ m][ 2][ j][k] = frct[ m][ 2][ j][k]
                 - dssp * ( - 4.0 * rsd[ m][ 1][ j][k]
                             + 6.0 * rsd[ m][ 2][ j][k]
                             - 4.0 * rsd[ m][ 3][ j][k]
                             +           rsd[ m][ 4][ j][k] );
            }

            for(i=3;i<=nx - 4;i++){
               for(m=0;m<=4;m++){
                  frct[ m][ i][ j][k] = frct[ m][ i][ j][k]
                    - dssp * (            rsd[ m][ i-2][ j][k]
                               - 4.0 * rsd[ m][ i-1][ j][k]
                               + 6.0 * rsd[ m][ i][ j][k]
                               - 4.0 * rsd[ m][ i+1][ j][k]
                               +           rsd[ m][ i+2][ j][k] );
               }
            }

            for(m=0;m<=4;m++){
               frct[ m][ nx-3][ j][k] = frct[ m][ nx-3][ j][k]
                 - dssp * (             rsd[ m][ nx-5][ j][k]
                             - 4.0 * rsd[ m][ nx-4][ j][k]
                             + 6.0 * rsd[ m][ nx-3][ j][k]
                             - 4.0 * rsd[ m][ nx-2][ j][k]  );
               frct[ m][ nx-2][ j][k] = frct[ m][ nx-2][ j][k]
                 - dssp * (             rsd[ m][ nx-4][ j][k]
                             - 4.0 * rsd[ m][ nx-3][ j][k]
                             + 5.0 * rsd[ m][ nx-2][ j][k] );
            }
         }
      }

//---------------------------------------------------------------------
//   eta-direction flux differences
//---------------------------------------------------------------------

      for(k=1;k<=nz - 2;k++){
         for(i=ist-1;i<=iend-1;i++){
            for(j=0;j<=ny-1;j++){
               flux[ 0][j] = rsd[ 2][ i][ j][k];
               u31 = rsd[ 2][ i][ j][k] / rsd[ 0][ i][ j][k];
               q = 0.50 * (  rsd[ 1][ i][ j][k] * rsd[ 1][ i][ j][k]
                               + rsd[ 2][ i][ j][k] * rsd[ 2][ i][ j][k]
                               + rsd[ 3][ i][ j][k] * rsd[ 3][ i][ j][k] )
                            / rsd[ 0][ i][ j][k];
               flux[ 1][j] = rsd[ 1][ i][ j][k] * u31 ;
               flux[ 2][j] = rsd[ 2][ i][ j][k] * u31 + c2 * 
                             ( rsd[ 4][ i][ j][k] - q );
               flux[ 3][j] = rsd[ 3][ i][ j][k] * u31;
               flux[ 4][j] = ( c1 * rsd[ 4][ i][ j][k] - c2 * q ) * u31;
            }

            for(j=jst-1;j<=jend-1;j++){
               for(m=0;m<=4;m++){
                  frct[ m][ i][ j][k] =  frct[ m][ i][ j][k]
                       - ty2 * ( flux[m][j+1] - flux[m][j-1] );
               }
            }

            for(j=jst-1;j<=ny-1;j++){
               tmp = 1.0 / rsd[ 0][ i][ j][k];

               u21j = tmp * rsd[ 1][ i][ j][k];
               u31j = tmp * rsd[ 2][ i][ j][k];
               u41j = tmp * rsd[ 3][ i][ j][k];
               u51j = tmp * rsd[ 4][ i][ j][k];

               tmp = 1.0 / rsd[ 0][ i][ j-1][k];

               u21jm1 = tmp * rsd[ 1][ i][ j-1][k];
               u31jm1 = tmp * rsd[ 2][ i][ j-1][k];
               u41jm1 = tmp * rsd[ 3][ i][ j-1][k];
               u51jm1 = tmp * rsd[ 4][ i][ j-1][k];

               flux[ 1][j] = ty3 * ( u21j - u21jm1 );
               flux[ 2][j] = (4.0/3.0) * ty3 * 
                             ( u31j - u31jm1 );
               flux[ 3][j] = ty3 * ( u41j - u41jm1 );
               flux[ 4][j] = 0.50 * ( 1.0 - c1*c5 )
                    * ty3 * ( ( pow2(u21j) + pow2(u31j) + pow2(u41j) )
                            - ( pow2(u21jm1) + pow2(u31jm1) + pow2(u41jm1) ) )
                    + (1.0/6.0)
                    * ty3 * ( pow2(u31j) - pow2(u31jm1) )
                    + c1 * c5 * ty3 * ( u51j - u51jm1 );
            }

            for(j=jst-1;j<=jend-1;j++){
               frct[ 0][ i][ j][k] = frct[ 0][ i][ j][k]
                    + dy1 * ty1 * (            rsd[ 0][ i][ j-1][k]
                                   - 2.0 * rsd[ 0][ i][ j][k]
                                   +           rsd[ 0][ i][ j+1][k] );
               frct[ 1][ i][ j][k] = frct[ 1][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 1][j+1] - flux[ 1][j] )
                    + dy2 * ty1 * (            rsd[ 1][ i][ j-1][k]
                                   - 2.0 * rsd[ 1][ i][ j][k]
                                   +           rsd[ 1][ i][ j+1][k] );
               frct[ 2][ i][ j][k] = frct[ 2][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 2][j+1] - flux[ 2][j] )
                    + dy3 * ty1 * (            rsd[ 2][ i][ j-1][k]
                                   - 2.0 * rsd[ 2][ i][ j][k]
                                   +           rsd[ 2][ i][ j+1][k] );
               frct[ 3][ i][ j][k] = frct[ 3][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 3][j+1] - flux[ 3][j] )
                    + dy4 * ty1 * (            rsd[ 3][ i][ j-1][k]
                                   - 2.0 * rsd[ 3][ i][ j][k]
                                   +           rsd[ 3][ i][ j+1][k] );
               frct[ 4][ i][ j][k] = frct[ 4][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 4][j+1] - flux[ 4][j] )
                    + dy5 * ty1 * (            rsd[ 4][ i][ j-1][k]
                                   - 2.0 * rsd[ 4][ i][ j][k]
                                   +           rsd[ 4][ i][ j+1][k] );
            }

//---------------------------------------------------------------------
//   fourth-order dissipation
//---------------------------------------------------------------------
            for(m=0;m<=4;m++){
               frct[ m][ i][ 1][k] = frct[ m][ i][ 1][k]
                 - dssp * ( + 5.0 * rsd[ m][ i][ 1][k]
                             - 4.0 * rsd[ m][ i][ 2][k]
                             +           rsd[ m][ i][ 3][k] );
               frct[ m][ i][ 2][k] = frct[ m][ i][ 2][k]
                 - dssp * ( - 4.0 * rsd[ m][ i][ 1][k]
                             + 6.0 * rsd[ m][ i][ 2][k]
                             - 4.0 * rsd[ m][ i][ 3][k]
                             +           rsd[ m][ i][ 4][k] );
            }

            for(j=3;j<=ny - 4;j++){
               for(m=0;m<=4;m++){
                  frct[ m][ i][ j][k] = frct[ m][ i][ j][k]
                    - dssp * (            rsd[ m][ i][ j-2][k]
                              - 4.0 * rsd[ m][ i][ j-1][k]
                              + 6.0 * rsd[ m][ i][ j][k]
                              - 4.0 * rsd[ m][ i][ j+1][k]
                              +           rsd[ m][ i][ j+2][k] );
               }
            }

            for(m=0;m<=4;m++){
               frct[ m][ i][ ny-3][k] = frct[ m][ i][ ny-3][k]
                 - dssp * (             rsd[ m][ i][ ny-5][k]
                             - 4.0 * rsd[ m][ i][ ny-4][k]
                             + 6.0 * rsd[ m][ i][ ny-3][k]
                             - 4.0 * rsd[ m][ i][ ny-2][k]  );
               frct[ m][ i][ ny-2][k] = frct[ m][ i][ ny-2][k]
                 - dssp * (             rsd[ m][ i][ ny-4][k]
                             - 4.0 * rsd[ m][ i][ ny-3][k]
                             + 5.0 * rsd[ m][ i][ ny-2][k]  );
            }
         }
      }

//---------------------------------------------------------------------
//   zeta-direction flux differences
//---------------------------------------------------------------------
      for(j=jst-1;j<=jend-1;j++){
         for(i=ist-1;i<=iend-1;i++){
            for(k=0;k<=nz-1;k++){
               flux[ 0][k] = rsd[ 3][ i][ j][k];
               u41 = rsd[ 3][ i][ j][k] / rsd[ 0][ i][ j][k];
               q = 0.50 * (  rsd[ 1][ i][ j][k] * rsd[ 1][ i][ j][k]
                               + rsd[ 2][ i][ j][k] * rsd[ 2][ i][ j][k]
                               + rsd[ 3][ i][ j][k] * rsd[ 3][ i][ j][k] )
                            / rsd[ 0][ i][ j][k];
               flux[ 1][k] = rsd[ 1][ i][ j][k] * u41 ;
               flux[ 2][k] = rsd[ 2][ i][ j][k] * u41 ;
               flux[ 3][k] = rsd[ 3][ i][ j][k] * u41 + c2 * 
                               ( rsd[ 4][ i][ j][k] - q );
               flux[ 4][k] = ( c1 * rsd[ 4][ i][ j][k] - c2 * q ) * u41;
            }

            for(k=1;k<=nz - 2;k++){
               for(m=0;m<=4;m++){
                  frct[ m][ i][ j][k] =  frct[ m][ i][ j][k]
                        - tz2 * ( flux[ m][k+1] - flux[ m][k-1] );
               }
            }

            for(k=1;k<=nz-1;k++){
               tmp = 1.0 / rsd[ 0][ i][ j][k];

               u21k = tmp * rsd[ 1][ i][ j][k];
               u31k = tmp * rsd[ 2][ i][ j][k];
               u41k = tmp * rsd[ 3][ i][ j][k];
               u51k = tmp * rsd[ 4][ i][ j][k];

               tmp = 1.0 / rsd[ 0][ i][ j][k-1];

               u21km1 = tmp * rsd[ 1][ i][ j][k-1];
               u31km1 = tmp * rsd[ 2][ i][ j][k-1];
               u41km1 = tmp * rsd[ 3][ i][ j][k-1];
               u51km1 = tmp * rsd[ 4][ i][ j][k-1];

               flux[ 1][k] = tz3 * ( u21k - u21km1 );
               flux[ 2][k] = tz3 * ( u31k - u31km1 );
               flux[ 3][k] = (4.0/3.0) * tz3 * ( u41k 
                             - u41km1 );
               flux[ 4][k] = 0.50 * ( 1.0 - c1*c5 )
                    * tz3 * ( ( pow2(u21k) + pow2(u31k) + pow2(u41k) )
                            - ( pow2(u21km1) + pow2(u31km1) + pow2(u41km1) ) )
                    + (1.0/6.0)
                    * tz3 * ( pow2(u41k) - pow2(u41km1) )
                    + c1 * c5 * tz3 * ( u51k - u51km1 );
            }

            for(k=1;k<=nz - 2;k++){
               frct[ 0][ i][ j][k] = frct[ 0][ i][ j][k]
                    + dz1 * tz1 * (            rsd[ 0][ i][ j][k+1]
                                   - 2.0 * rsd[ 0][ i][ j][k]
                                   +           rsd[ 0][ i][ j][k-1] );
               frct[ 1][ i][ j][k] = frct[ 1][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 1][k+1] - flux[ 1][k] )
                    + dz2 * tz1 * (            rsd[ 1][ i][ j][k+1]
                                   - 2.0 * rsd[ 1][ i][ j][k]
                                   +           rsd[ 1][ i][ j][k-1] );
               frct[ 2][ i][ j][k] = frct[ 2][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 2][k+1] - flux[ 2][k] )
                    + dz3 * tz1 * (            rsd[ 2][ i][ j][k+1]
                                   - 2.0 * rsd[ 2][ i][ j][k]
                                   +           rsd[ 2][ i][ j][k-1] );
               frct[ 3][ i][ j][k] = frct[ 3][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 3][k+1] - flux[ 3][k] )
                    + dz4 * tz1 * (            rsd[ 3][ i][ j][k+1]
                                   - 2.0 * rsd[ 3][ i][ j][k]
                                   +           rsd[ 3][ i][ j][k-1] );
               frct[ 4][ i][ j][k] = frct[ 4][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 4][k+1] - flux[ 4][k] )
                    + dz5 * tz1 * (            rsd[ 4][ i][ j][k+1]
                                   - 2.0 * rsd[ 4][ i][ j][k]
                                   +           rsd[ 4][ i][ j][k-1] );
            }

//---------------------------------------------------------------------
//   fourth-order dissipation
//---------------------------------------------------------------------
            for(m=0;m<=4;m++){
               frct[ m][ i][ j][1] = frct[ m][ i][ j][1]
                 - dssp * ( + 5.0 * rsd[ m][ i][ j][1]
                             - 4.0 * rsd[ m][ i][ j][2]
                             +           rsd[ m][ i][ j][3] );
               frct[ m][ i][ j][2] = frct[ m][ i][ j][2]
                 - dssp * (- 4.0 * rsd[ m][ i][ j][1]
                            + 6.0 * rsd[ m][ i][ j][2]
                            - 4.0 * rsd[ m][ i][ j][3]
                            +           rsd[ m][ i][ j][4] );
            }

            for(k=3;k<=nz - 4;k++){
               for(m=0;m<=4;m++){
                  frct[ m][ i][ j][k] = frct[ m][ i][ j][k]
                    - dssp * (           rsd[ m][ i][ j][k-2]
                              - 4.0 * rsd[ m][ i][ j][k-1]
                              + 6.0 * rsd[ m][ i][ j][k]
                              - 4.0 * rsd[ m][ i][ j][k+1]
                              +           rsd[ m][ i][ j][k+2] );
               }
            }

            for(m=0;m<=4;m++){
               frct[ m][ i][ j][nz-3] = frct[ m][ i][ j][nz-3]
                 - dssp * (            rsd[ m][ i][ j][nz-5]
                            - 4.0 * rsd[ m][ i][ j][nz-4]
                            + 6.0 * rsd[ m][ i][ j][nz-3]
                            - 4.0 * rsd[ m][ i][ j][nz-2]  );
               frct[ m][ i][ j][nz-2] = frct[ m][ i][ j][nz-2]
                 - dssp * (             rsd[ m][ i][ j][nz-4]
                             - 4.0 * rsd[ m][ i][ j][nz-3]
                             + 5.0 * rsd[ m][ i][ j][nz-2]  );
            }
         }
      }
  }

  public void error(){
    int i, j, k, m;
    double  tmp;
    double[]  u000ijk = new double[5];

    for(m=0;m<=4;m++){
      errnm[m] = 0.0;
    }

    for(k=1;k<=nz-2;k++){
      for(j=jst-1;j<=jend-1;j++){
	for(i=ist-1;i<=iend-1;i++){
	  exact( i+1, j+1, k+1, u000ijk );
	  for(m=0;m<=4;m++){
	    tmp = ( u000ijk[m] - u[ m][ i][ j][k] );
	    errnm[m] = errnm[m] + pow2(tmp);
	  }
	}
      }
    }

    for(m=0;m<=4;m++){
      errnm[m] = Math.Sqrt( errnm[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
    }
  }
    
  public void jacld(int k){
      int i, j;
      double  r43;
      double  c1345;
      double  c34;
      double  tmp1, tmp2, tmp3;

      r43 = ( 4.0 / 3.0 );
      c1345 = c1 * c3 * c4 * c5;
      c34 = c3 * c4;

         for(j=jst-1;j<=jend-1;j++){
            for(i=ist-1;i<=iend-1;i++){
//---------------------------------------------------------------------
//   form the block daigonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i][ j][k];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               d[ 0][ 0][ i][j] =  1.0
                             + dt * 2.0 * (   tx1 * dx1
                                            + ty1 * dy1
                                            + tz1 * dz1 );
               d[ 0][ 1][ i][j] =  0.0;
               d[ 0][ 2][ i][j] =  0.0;
               d[ 0][ 3][ i][j] =  0.0;
               d[ 0][ 4][ i][j] =  0.0;

               d[ 1][ 0][ i][j] = -dt * 2.0
                * (  tx1 * r43 + ty1 + tz1  )
                * c34 * tmp2 * u[ 1][ i][ j][k];
               d[ 1][ 1][ i][j] =  1.0
                + dt * 2.0 * c34 * tmp1 
                * (  tx1 * r43 + ty1 + tz1 )
                + dt * 2.0 * (   tx1 * dx2
                                   + ty1 * dy2
                                   + tz1 * dz2  );
               d[ 1][ 2][ i][j] = 0.0;
               d[ 1][ 3][ i][j] = 0.0;
               d[ 1][ 4][ i][j] = 0.0;

               d[ 2][ 0][ i][j] = -dt * 2.0
                 * (  tx1 + ty1 * r43 + tz1  )
                 * c34 * tmp2 * u[ 2][ i][ j][k];
               d[ 2][ 1][ i][j] = 0.0;
               d[ 2][ 2][ i][j] = 1.0
               + dt * 2.0 * c34 * tmp1
                    * (  tx1 + ty1 * r43 + tz1 )
               + dt * 2.0 * (  tx1 * dx3
                                 + ty1 * dy3
                                 + tz1 * dz3 );
               d[ 2][ 3][ i][j] = 0.0;
               d[ 2][ 4][ i][j] = 0.0;

               d[ 3][ 0][ i][j] = -dt * 2.0
                 * (  tx1 + ty1 + tz1 * r43  )
                 * c34 * tmp2 * u[ 3][ i][ j][k];
               d[ 3][ 1][ i][j] = 0.0;
               d[ 3][ 2][ i][j] = 0.0;
               d[ 3][ 3][ i][j] = 1.0
               + dt * 2.0 * c34 * tmp1
                    * (  tx1 + ty1 + tz1 * r43 )
               + dt * 2.0 * (  tx1 * dx4
                                 + ty1 * dy4
                                 + tz1 * dz4 );
               d[ 3][ 4][ i][j] = 0.0;

               d[ 4][ 0][ i][j] = -dt * 2.0
        * ( ( ( tx1 * ( r43*c34 - c1345 )
           + ty1 * ( c34 - c1345 )
           + tz1 * ( c34 - c1345 ) ) * ( pow2(u[ 1][ i][ j][k]) )
         + ( tx1 * ( c34 - c1345 )
           + ty1 * ( r43*c34 - c1345 )
           + tz1 * ( c34 - c1345 ) ) * ( pow2(u[ 2][ i][ j][k]) )
         + ( tx1 * ( c34 - c1345 )
           + ty1 * ( c34 - c1345 )
           + tz1 * ( r43*c34 - c1345 ) ) * ( pow2(u[ 3][ i][ j][k]) )
            ) * tmp3
         + ( tx1 + ty1 + tz1 ) * c1345 * tmp2 * u[ 4][ i][ j][k] );

               d[ 4][ 1][ i][j] = dt * 2.0 * tmp2 * u[ 1][ i][ j][k]
       * ( tx1 * ( r43*c34 - c1345 )
         + ty1 * (     c34 - c1345 )
         + tz1 * (     c34 - c1345 ) );
               d[ 4][ 2][ i][j] = dt * 2.0 * tmp2 * u[ 2][ i][ j][k]
       * ( tx1 * ( c34 - c1345 )
         + ty1 * ( r43*c34 -c1345 )
         + tz1 * ( c34 - c1345 ) );
               d[ 4][ 3][ i][j] = dt * 2.0 * tmp2 * u[ 3][ i][ j][k]
       * ( tx1 * ( c34 - c1345 )
         + ty1 * ( c34 - c1345 )
         + tz1 * ( r43*c34 - c1345 ) );
               d[ 4][ 4][ i][j] = 1.0
         + dt * 2.0 * ( tx1  + ty1 + tz1 ) * c1345 * tmp1
         + dt * 2.0 * (  tx1 * dx5
                          +  ty1 * dy5
                          +  tz1 * dz5 );

//---------------------------------------------------------------------
//   form the first block sub-diagonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i][ j][k-1];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               a[ 0][ 0][ i][j] = - dt * tz1 * dz1;
               a[ 0][ 1][ i][j] =   0.0;
               a[ 0][ 2][ i][j] =   0.0;
               a[ 0][ 3][ i][j] = - dt * tz2;
               a[ 0][ 4][ i][j] =   0.0;

               a[ 1][ 0][ i][j] = - dt * tz2
                 * ( - ( u[ 1][ i][ j][k-1]*u[ 3][ i][ j][k-1] ) * tmp2 )
                 - dt * tz1 * ( - c34 * tmp2 * u[ 1][ i][ j][k-1] );
               a[ 1][ 1][ i][j] = - dt * tz2 * ( u[ 3][ i][ j][k-1] * tmp1 )
                 - dt * tz1 * c34 * tmp1
                 - dt * tz1 * dz2 ;
               a[ 1][ 2][ i][j] = 0.0;
               a[ 1][ 3][ i][j] = - dt * tz2 * ( u[ 1][ i][ j][k-1] * tmp1 );
               a[ 1][ 4][ i][j] = 0.0;

               a[ 2][ 0][ i][j] = - dt * tz2
                 * ( - ( u[ 2][ i][ j][k-1]*u[ 3][ i][ j][k-1] ) * tmp2 )
                 - dt * tz1 * ( - c34 * tmp2 * u[ 2][ i][ j][k-1] );
               a[ 2][ 1][ i][j] = 0.0;
               a[ 2][ 2][ i][j] = - dt * tz2 * ( u[ 3][ i][ j][k-1] * tmp1 )
                 - dt * tz1 * ( c34 * tmp1 )
                 - dt * tz1 * dz3;
               a[ 2][ 3][ i][j] = - dt * tz2 * ( u[ 2][ i][ j][k-1] * tmp1 );
               a[ 2][ 4][ i][j] = 0.0;

               a[ 3][ 0][ i][j] = - dt * tz2
              * ( - pow2(( u[ 3][ i][ j][k-1] * tmp1 ))
                  + c2 * qs[ i][ j][k-1] * tmp1 )
              - dt * tz1 * ( - r43 * c34 * tmp2 * u[ 3][ i][ j][k-1] );
               a[ 3][ 1][ i][j] = - dt * tz2
                   * ( - c2 * ( u[ 1][ i][ j][k-1] * tmp1 ) );
               a[ 3][ 2][ i][j] = - dt * tz2
                   * ( - c2 * ( u[ 2][ i][ j][k-1] * tmp1 ) );
               a[ 3][ 3][ i][j] = - dt * tz2 * ( 2.0 - c2 )
                   * ( u[ 3][ i][ j][k-1] * tmp1 )
                   - dt * tz1 * ( r43 * c34 * tmp1 )
                   - dt * tz1 * dz4;
               a[ 3][ 4][ i][j] = - dt * tz2 * c2;

               a[ 4][ 0][ i][j] = - dt * tz2
             * ( ( c2 * 2.0 * qs[ i][ j][k-1]
             - c1 * u[ 4][ i][ j][k-1] )
                  * u[ 3][ i][ j][k-1] * tmp2 )
             - dt * tz1
             * ( - ( c34 - c1345 ) * tmp3 * pow2(u[ 1][ i][ j][k-1])
                 - ( c34 - c1345 ) * tmp3 * pow2(u[ 2][ i][ j][k-1])
                 - ( r43*c34 - c1345 )* tmp3 * pow2(u[ 3][ i][ j][k-1])
                - c1345 * tmp2 * u[ 4][ i][ j][k-1] );
               a[ 4][ 1][ i][j] = - dt * tz2
             * ( - c2 * ( u[ 1][ i][ j][k-1]*u[ 3][ i][ j][k-1] ) * tmp2 )
             - dt * tz1 * ( c34 - c1345 ) * tmp2 * u[ 1][ i][ j][k-1];
               a[ 4][ 2][ i][j] = - dt * tz2
             * ( - c2 * ( u[ 2][ i][ j][k-1]*u[ 3][ i][ j][k-1] ) * tmp2 )
             - dt * tz1 * ( c34 - c1345 ) * tmp2 * u[ 2][ i][ j][k-1];
               a[ 4][ 3][ i][j] = - dt * tz2
             * ( c1 * ( u[ 4][ i][ j][k-1] * tmp1 )
             - c2
             * ( qs[ i][ j][k-1] * tmp1
                  + u[ 3][ i][ j][k-1]*u[ 3][ i][ j][k-1] * tmp2 ) )
             - dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u[ 3][ i][ j][k-1];
               a[ 4][ 4][ i][j] = - dt * tz2
             * ( c1 * ( u[ 3][ i][ j][k-1] * tmp1 ) )
             - dt * tz1 * c1345 * tmp1
             - dt * tz1 * dz5;

//---------------------------------------------------------------------
//   form the second block sub-diagonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i][ j-1][k];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               b[ 0][ 0][ i][j] = - dt * ty1 * dy1;
               b[ 0][ 1][ i][j] =   0.0;
               b[ 0][ 2][ i][j] = - dt * ty2;
               b[ 0][ 3][ i][j] =   0.0;
               b[ 0][ 4][ i][j] =   0.0;

               b[ 1][ 0][ i][j] = - dt * ty2
                 * ( - ( u[ 1][ i][ j-1][k]*u[ 2][ i][ j-1][k] ) * tmp2 )
                 - dt * ty1 * ( - c34 * tmp2 * u[ 1][ i][ j-1][k] );
               b[ 1][ 1][ i][j] = - dt * ty2 * ( u[ 2][ i][ j-1][k] * tmp1 )
                - dt * ty1 * ( c34 * tmp1 )
                - dt * ty1 * dy2;
               b[ 1][ 2][ i][j] = - dt * ty2 * ( u[ 1][ i][ j-1][k] * tmp1 );
               b[ 1][ 3][ i][j] = 0.0;
               b[ 1][ 4][ i][j] = 0.0;

               b[ 2][ 0][ i][j] = - dt * ty2
                 * ( - pow2(( u[ 2][ i][ j-1][k] * tmp1 ))
             + c2 * ( qs[ i][ j-1][k] * tmp1 ) )
             - dt * ty1 * ( - r43 * c34 * tmp2 * u[ 2][ i][ j-1][k] );
               b[ 2][ 1][ i][j] = - dt * ty2
                         * ( - c2 * ( u[ 1][ i][ j-1][k] * tmp1 ) );
               b[ 2][ 2][ i][j] = - dt * ty2 * ( ( 2.0 - c2 )
                         * ( u[ 2][ i][ j-1][k] * tmp1 ) )
             - dt * ty1 * ( r43 * c34 * tmp1 )
             - dt * ty1 * dy3;
               b[ 2][ 3][ i][j] = - dt * ty2
                         * ( - c2 * ( u[ 3][ i][ j-1][k] * tmp1 ) );
               b[ 2][ 4][ i][j] = - dt * ty2 * c2;

               b[ 3][ 0][ i][j] = - dt * ty2
                    * ( - ( u[ 2][ i][ j-1][k]*u[ 3][ i][ j-1][k] ) * tmp2 )
             - dt * ty1 * ( - c34 * tmp2 * u[ 3][ i][ j-1][k] );
               b[ 3][ 1][ i][j] = 0.0;
               b[ 3][ 2][ i][j] = - dt * ty2 * ( u[ 3][ i][ j-1][k] * tmp1 );
               b[ 3][ 3][ i][j] = - dt * ty2 * ( u[ 2][ i][ j-1][k] * tmp1 )
                              - dt * ty1 * ( c34 * tmp1 )
                              - dt * ty1 * dy4;
               b[ 3][ 4][ i][j] = 0.0;

               b[ 4][ 0][ i][j] = - dt * ty2
                * ( ( c2 * 2.0 * qs[ i][ j-1][k]
                     - c1 * u[ 4][ i][ j-1][k] )
                * ( u[ 2][ i][ j-1][k] * tmp2 ) )
                - dt * ty1
                * ( - (     c34 - c1345 )*tmp3*pow2(u[ 1][ i][ j-1][k])
                    - ( r43*c34 - c1345 )*tmp3*pow2(u[ 2][ i][ j-1][k])
                    - (     c34 - c1345 )*tmp3*pow2(u[ 3][ i][ j-1][k])
                    - c1345*tmp2*u[ 4][ i][ j-1][k] );
               b[ 4][ 1][ i][j] = - dt * ty2
                * ( - c2 * ( u[ 1][ i][ j-1][k]*u[ 2][ i][ j-1][k] ) * tmp2 )
                - dt * ty1
                * ( c34 - c1345 ) * tmp2 * u[ 1][ i][ j-1][k];
               b[ 4][ 2][ i][j] = - dt * ty2
                * ( c1 * ( u[ 4][ i][ j-1][k] * tmp1 )
                - c2 
                * ( qs[ i][ j-1][k] * tmp1
                     + u[ 2][ i][ j-1][k]*u[ 2][ i][ j-1][k] * tmp2 ) )
                - dt * ty1
                * ( r43*c34 - c1345 ) * tmp2 * u[ 2][ i][ j-1][k];
               b[ 4][ 3][ i][j] = - dt * ty2
                * ( - c2 * ( u[ 2][ i][ j-1][k]*u[ 3][ i][ j-1][k] ) * tmp2 )
                - dt * ty1 * ( c34 - c1345 ) * tmp2 * u[ 3][ i][ j-1][k];
               b[ 4][ 4][ i][j] = - dt * ty2
                * ( c1 * ( u[ 2][ i][ j-1][k] * tmp1 ) )
                - dt * ty1 * c1345 * tmp1
                - dt * ty1 * dy5;
	    
//---------------------------------------------------------------------
//   form the third block sub-diagonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i-1][ j][k];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               c[ 0][ 0][ i][j] = - dt * tx1 * dx1;
               c[ 0][ 1][ i][j] = - dt * tx2;
               c[ 0][ 2][ i][j] =   0.0;
               c[ 0][ 3][ i][j] =   0.0;
               c[ 0][ 4][ i][j] =   0.0;

               c[ 1][ 0][ i][j] = - dt * tx2
                * ( - pow2(( u[ 1][ i-1][ j][k] * tmp1 ))
             + c2 * qs[ i-1][ j][k] * tmp1 )
                - dt * tx1 * ( - r43 * c34 * tmp2 * u[ 1][ i-1][ j][k] );
               c[ 1][ 1][ i][j] = - dt * tx2
                * ( ( 2.0 - c2 ) * ( u[ 1][ i-1][ j][k] * tmp1 ) )
                - dt * tx1 * ( r43 * c34 * tmp1 )
                - dt * tx1 * dx2;
               c[ 1][ 2][ i][j] = - dt * tx2
                    * ( - c2 * ( u[ 2][ i-1][ j][k] * tmp1 ) );
               c[ 1][ 3][ i][j] = - dt * tx2
                    * ( - c2 * ( u[ 3][ i-1][ j][k] * tmp1 ) );
               c[ 1][ 4][ i][j] = - dt * tx2 * c2 ;

               c[ 2][ 0][ i][j] = - dt * tx2
                    * ( - ( u[ 1][ i-1][ j][k] * u[ 2][ i-1][ j][k] ) * tmp2 )
               - dt * tx1 * ( - c34 * tmp2 * u[ 2][ i-1][ j][k] );
               c[ 2][ 1][ i][j] = - dt * tx2 * ( u[ 2][ i-1][ j][k] * tmp1 );
               c[ 2][ 2][ i][j] = - dt * tx2 * ( u[ 1][ i-1][ j][k] * tmp1 )
                - dt * tx1 * ( c34 * tmp1 )
                - dt * tx1 * dx3;
               c[ 2][ 3][ i][j] = 0.0;
               c[ 2][ 4][ i][j] = 0.0;

               c[ 3][ 0][ i][j] = - dt * tx2
                * ( - ( u[ 1][ i-1][ j][k]*u[ 3][ i-1][ j][k] ) * tmp2 )
                - dt * tx1 * ( - c34 * tmp2 * u[ 3][ i-1][ j][k] );
               c[ 3][ 1][ i][j] = - dt * tx2 * ( u[ 3][ i-1][ j][k] * tmp1 );
               c[ 3][ 2][ i][j] = 0.0;
               c[ 3][ 3][ i][j] = - dt * tx2 * ( u[ 1][ i-1][ j][k] * tmp1 )
                - dt * tx1 * ( c34 * tmp1 )
                - dt * tx1 * dx4;
               c[ 3][ 4][ i][j] = 0.0;

               c[ 4][ 0][ i][j] = - dt * tx2
                * ( ( c2 * 2.0 * qs[ i-1][ j][k]
                    - c1 * u[ 4][ i-1][ j][k] )
                * u[ 1][ i-1][ j][k] * tmp2 )
                - dt * tx1
                * ( - ( r43*c34 - c1345 ) * tmp3 * pow2( u[ 1][ i-1][ j][k])
                    - (     c34 - c1345 ) * tmp3 * pow2( u[ 2][ i-1][ j][k])
                    - (     c34 - c1345 ) * tmp3 * pow2( u[ 3][ i-1][ j][k])
                    - c1345 * tmp2 * u[ 4][ i-1][ j][k] );
               c[ 4][ 1][ i][j] = - dt * tx2
                * ( c1 * ( u[ 4][ i-1][ j][k] * tmp1 )
                   - c2
                   * ( u[ 1][ i-1][ j][k]*u[ 1][ i-1][ j][k] * tmp2
                        + qs[ i-1][ j][k] * tmp1 ) )
                 - dt * tx1
                 * ( r43*c34 - c1345 ) * tmp2 * u[ 1][ i-1][ j][k];
               c[ 4][ 2][ i][j] = - dt * tx2
                 * ( - c2 * ( u[ 2][ i-1][ j][k]*u[ 1][ i-1][ j][k] ) * tmp2 )
                 - dt * tx1
                 * (  c34 - c1345 ) * tmp2 * u[ 2][ i-1][ j][k];
               c[ 4][ 3][ i][j] = - dt * tx2
                 * ( - c2 * ( u[ 3][ i-1][ j][k]*u[ 1][ i-1][ j][k] ) * tmp2 )
                 - dt * tx1
                 * (  c34 - c1345 ) * tmp2 * u[ 3][ i-1][ j][k];
               c[ 4][ 4][ i][j] = - dt * tx2
                 * ( c1 * ( u[ 1][ i-1][ j][k] * tmp1 ) )
                 - dt * tx1 * c1345 * tmp1
                 - dt * tx1 * dx5;
            }
         }
  }

  public void jacu(int k){
      int i, j;    
      double  r43;
      double  c1345;
      double  c34;
      double  tmp1, tmp2, tmp3;

      r43 = ( 4.0 / 3.0 );
      c1345 = c1 * c3 * c4 * c5;
      c34 = c3 * c4;

         for(j=jend-1;j>=jst-1;j--){
            for(i=iend-1;i>=ist-1;i--){

//---------------------------------------------------------------------
//   form the block daigonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i][ j][k];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               d[ 0][ 0][ i][j] =  1.0
                             + dt * 2.0 * (   tx1 * dx1
                                                + ty1 * dy1
                                                + tz1 * dz1 );
               d[ 0][ 1][ i][j] =  0.0;
               d[ 0][ 2][ i][j] =  0.0;
               d[ 0][ 3][ i][j] =  0.0;
               d[ 0][ 4][ i][j] =  0.0;

               d[ 1][ 0][ i][j] =  dt * 2.0
                 * ( - tx1 * r43 - ty1 - tz1 )
                 * ( c34 * tmp2 * u[ 1][ i][ j][k] );
               d[ 1][ 1][ i][j] =  1.0
                + dt * 2.0 * c34 * tmp1 
                * (  tx1 * r43 + ty1 + tz1 )
                + dt * 2.0 * (   tx1 * dx2
                                   + ty1 * dy2
                                   + tz1 * dz2  );
               d[ 1][ 2][ i][j] = 0.0;
               d[ 1][ 3][ i][j] = 0.0;
               d[ 1][ 4][ i][j] = 0.0;

               d[ 2][ 0][ i][j] = dt * 2.0
                 * ( - tx1 - ty1 * r43 - tz1 )
                 * ( c34 * tmp2 * u[ 2][ i][ j][k] );
               d[ 2][ 1][ i][j] = 0.0;
               d[ 2][ 2][ i][j] = 1.0
               + dt * 2.0 * c34 * tmp1
                    * (  tx1 + ty1 * r43 + tz1 )
               + dt * 2.0 * (  tx1 * dx3
                                 + ty1 * dy3
                                 + tz1 * dz3 );
               d[ 2][ 3][ i][j] = 0.0;
               d[ 2][ 4][ i][j] = 0.0;

               d[ 3][ 0][ i][j] = dt * 2.0
                 * ( - tx1 - ty1 - tz1 * r43 )
                 * ( c34 * tmp2 * u[ 3][ i][ j][k] );
               d[ 3][ 1][ i][j] = 0.0;
               d[ 3][ 2][ i][j] = 0.0;
               d[ 3][ 3][ i][j] = 1.0
               + dt * 2.0 * c34 * tmp1
                    * (  tx1 + ty1 + tz1 * r43 )
               + dt * 2.0 * (  tx1 * dx4
                                 + ty1 * dy4
                                 + tz1 * dz4 );
               d[ 3][ 4][ i][j] = 0.0;

               d[ 4][ 0][ i][j] = -dt * 2.0
        * ( ( ( tx1 * ( r43*c34 - c1345 )
           + ty1 * ( c34 - c1345 )
           + tz1 * ( c34 - c1345 ) ) * pow2( u[ 1][ i][ j][k])
         + ( tx1 * ( c34 - c1345 )
           + ty1 * ( r43*c34 - c1345 )
           + tz1 * ( c34 - c1345 ) ) * pow2( u[ 2][ i][ j][k])
         + ( tx1 * ( c34 - c1345 )
           + ty1 * ( c34 - c1345 )
           + tz1 * ( r43*c34 - c1345 ) ) * pow2( u[ 3][ i][ j][k])
            ) * tmp3
         + ( tx1 + ty1 + tz1 ) * c1345 * tmp2 * u[ 4][ i][ j][k] );

               d[ 4][ 1][ i][j] = dt * 2.0
       * ( tx1 * ( r43*c34 - c1345 )
         + ty1 * (     c34 - c1345 )
         + tz1 * (     c34 - c1345 ) ) * tmp2 * u[ 1][ i][ j][k];
               d[ 4][ 2][ i][j] = dt * 2.0
       * ( tx1 * ( c34 - c1345 )
         + ty1 * ( r43*c34 -c1345 )
         + tz1 * ( c34 - c1345 ) ) * tmp2 * u[ 2][ i][ j][k];
               d[ 4][ 3][ i][j] = dt * 2.0
       * ( tx1 * ( c34 - c1345 )
         + ty1 * ( c34 - c1345 )
         + tz1 * ( r43*c34 - c1345 ) ) * tmp2 * u[ 3][ i][ j][k];
               d[ 4][ 4][ i][j] = 1.0
         + dt * 2.0 * ( tx1 + ty1 + tz1 ) * c1345 * tmp1
         + dt * 2.0 * (  tx1 * dx5
                          +  ty1 * dy5
                          +  tz1 * dz5 );
//---------------------------------------------------------------------
//   form the first block sub-diagonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i+1][ j][k];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               a[ 0][ 0][ i][j] = - dt * tx1 * dx1;
               a[ 0][ 1][ i][j] =   dt * tx2;
               a[ 0][ 2][ i][j] =   0.0;
               a[ 0][ 3][ i][j] =   0.0;
               a[ 0][ 4][ i][j] =   0.0;

               a[ 1][ 0][ i][j] =  dt * tx2
                * ( - pow2(( u[ 1][ i+1][ j][k] * tmp1 ))
           + c2 * qs[ i+1][ j][k] * tmp1 )
                - dt * tx1 * ( - r43 * c34 * tmp2 * u[ 1][ i+1][ j][k] );
               a[ 1][ 1][ i][j] =  dt * tx2
                * ( ( 2.0 - c2 ) * ( u[ 1][ i+1][ j][k] * tmp1 ) )
                - dt * tx1 * ( r43 * c34 * tmp1 )
                - dt * tx1 * dx2;
               a[ 1][ 2][ i][j] =  dt * tx2
                    * ( - c2 * ( u[ 2][ i+1][ j][k] * tmp1 ) );
               a[ 1][ 3][ i][j] =  dt * tx2
                    * ( - c2 * ( u[ 3][ i+1][ j][k] * tmp1 ) );
               a[ 1][ 4][ i][j] =  dt * tx2 * c2 ;

               a[ 2][ 0][ i][j] =  dt * tx2
                    * ( - ( u[ 1][ i+1][ j][k] * u[ 2][ i+1][ j][k] ) * tmp2 )
               - dt * tx1 * ( - c34 * tmp2 * u[ 2][ i+1][ j][k] );
               a[ 2][ 1][ i][j] =  dt * tx2 * ( u[ 2][ i+1][ j][k] * tmp1 );
               a[ 2][ 2][ i][j] =  dt * tx2 * ( u[ 1][ i+1][ j][k] * tmp1 )
                - dt * tx1 * ( c34 * tmp1 )
                - dt * tx1 * dx3;
               a[ 2][ 3][ i][j] = 0.0;
               a[ 2][ 4][ i][j] = 0.0;

               a[ 3][ 0][ i][j] = dt * tx2
                * ( - ( u[ 1][ i+1][ j][k]*u[ 3][ i+1][ j][k] ) * tmp2 )
                - dt * tx1 * ( - c34 * tmp2 * u[ 3][ i+1][ j][k] );
               a[ 3][ 1][ i][j] = dt * tx2 * ( u[ 3][ i+1][ j][k] * tmp1 );
               a[ 3][ 2][ i][j] = 0.0;
               a[ 3][ 3][ i][j] = dt * tx2 * ( u[ 1][ i+1][ j][k] * tmp1 )
                - dt * tx1 * ( c34 * tmp1 )
                - dt * tx1 * dx4;
               a[ 3][ 4][ i][j] = 0.0;

               a[ 4][ 0][ i][j] = dt * tx2
                * ( ( c2 * 2.0 * qs[ i+1][ j][k]
                    - c1 * u[ 4][ i+1][ j][k] )
                * ( u[ 1][ i+1][ j][k] * tmp2 ) )
                - dt * tx1
                * ( - ( r43*c34 - c1345 ) * tmp3 * pow2( u[ 1][ i+1][ j][k])
                    - (     c34 - c1345 ) * tmp3 * pow2( u[ 2][ i+1][ j][k])
                    - (     c34 - c1345 ) * tmp3 * pow2( u[ 3][ i+1][ j][k])
                    - c1345 * tmp2 * u[ 4][ i+1][ j][k] );
               a[ 4][ 1][ i][j] = dt * tx2
                * ( c1 * ( u[ 4][ i+1][ j][k] * tmp1 )
                   - c2
                   * (  u[ 1][ i+1][ j][k]*u[ 1][ i+1][ j][k] * tmp2
                        + qs[ i+1][ j][k] * tmp1 ) )
                 - dt * tx1
                 * ( r43*c34 - c1345 ) * tmp2 * u[ 1][ i+1][ j][k];
               a[ 4][ 2][ i][j] = dt * tx2
                 * ( - c2 * ( u[ 2][ i+1][ j][k]*u[ 1][ i+1][ j][k] ) * tmp2 )
                 - dt * tx1
                 * (  c34 - c1345 ) * tmp2 * u[ 2][ i+1][ j][k];
               a[ 4][ 3][ i][j] = dt * tx2
                 * ( - c2 * ( u[ 3][ i+1][ j][k]*u[ 1][ i+1][ j][k] ) * tmp2 )
                 - dt * tx1
                 * (  c34 - c1345 ) * tmp2 * u[ 3][ i+1][ j][k];
               a[ 4][ 4][ i][j] = dt * tx2
                 * ( c1 * ( u[ 1][ i+1][ j][k] * tmp1 ) )
                 - dt * tx1 * c1345 * tmp1
                 - dt * tx1 * dx5;

//---------------------------------------------------------------------
//   form the second block sub-diagonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i][ j+1][k];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               b[ 0][ 0][ i][j] = - dt * ty1 * dy1;
               b[ 0][ 1][ i][j] =   0.0;
               b[ 0][ 2][ i][j] =  dt * ty2;
               b[ 0][ 3][ i][j] =   0.0;
               b[ 0][ 4][ i][j] =   0.0;

               b[ 1][ 0][ i][j] =  dt * ty2
                 * ( - ( u[ 1][ i][ j+1][k]*u[ 2][ i][ j+1][k] ) * tmp2 )
                 - dt * ty1 * ( - c34 * tmp2 * u[ 1][ i][ j+1][k] );
               b[ 1][ 1][ i][j] =  dt * ty2 * ( u[ 2][ i][ j+1][k] * tmp1 )
                - dt * ty1 * ( c34 * tmp1 )
                - dt * ty1 * dy2;
               b[ 1][ 2][ i][j] =  dt * ty2 * ( u[ 1][ i][ j+1][k] * tmp1 );
               b[ 1][ 3][ i][j] = 0.0;
               b[ 1][ 4][ i][j] = 0.0;

               b[ 2][ 0][ i][j] =  dt * ty2
                 * ( - pow2(( u[ 2][ i][ j+1][k] * tmp1 ))
            + c2 * ( qs[ i][ j+1][k] * tmp1 ) )
             - dt * ty1 * ( - r43 * c34 * tmp2 * u[ 2][ i][ j+1][k] );
               b[ 2][ 1][ i][j] =  dt * ty2
                         * ( - c2 * ( u[ 1][ i][ j+1][k] * tmp1 ) );
               b[ 2][ 2][ i][j] =  dt * ty2 * ( ( 2.0 - c2 )
                         * ( u[ 2][ i][ j+1][k] * tmp1 ) )
             - dt * ty1 * ( r43 * c34 * tmp1 )
             - dt * ty1 * dy3;
               b[ 2][ 3][ i][j] =  dt * ty2
                         * ( - c2 * ( u[ 3][ i][ j+1][k] * tmp1 ) );
               b[ 2][ 4][ i][j] =  dt * ty2 * c2;

               b[ 3][ 0][ i][j] =  dt * ty2
                    * ( - ( u[ 2][ i][ j+1][k]*u[ 3][ i][ j+1][k] ) * tmp2 )
             - dt * ty1 * ( - c34 * tmp2 * u[ 3][ i][ j+1][k] );
               b[ 3][ 1][ i][j] = 0.0;
               b[ 3][ 2][ i][j] =  dt * ty2 * ( u[ 3][ i][ j+1][k] * tmp1 );
               b[ 3][ 3][ i][j] =  dt * ty2 * ( u[ 2][ i][ j+1][k] * tmp1 )
                              - dt * ty1 * ( c34 * tmp1 )
                              - dt * ty1 * dy4;
               b[ 3][ 4][ i][j] = 0.0;

               b[ 4][ 0][ i][j] =  dt * ty2
                * ( ( c2 * 2.0 * qs[ i][ j+1][k]
                     - c1 * u[ 4][ i][ j+1][k] )
                * ( u[ 2][ i][ j+1][k] * tmp2 ) )
                - dt * ty1
                * ( - (     c34 - c1345 )*tmp3*pow2(u[ 1][ i][ j+1][k])
                    - ( r43*c34 - c1345 )*tmp3*pow2(u[ 2][ i][ j+1][k])
                    - (     c34 - c1345 )*tmp3*pow2(u[ 3][ i][ j+1][k])
                    - c1345*tmp2*u[ 4][ i][ j+1][k] );
               b[ 4][ 1][ i][j] =  dt * ty2
                * ( - c2 * ( u[ 1][ i][ j+1][k]*u[ 2][ i][ j+1][k] ) * tmp2 )
                - dt * ty1
                * ( c34 - c1345 ) * tmp2 * u[ 1][ i][ j+1][k];
               b[ 4][ 2][ i][j] =  dt * ty2
                * ( c1 * ( u[ 4][ i][ j+1][k] * tmp1 )
                - c2 
                * ( qs[ i][ j+1][k] * tmp1
                     + u[ 2][ i][ j+1][k]*u[ 2][ i][ j+1][k] * tmp2 ) )
                - dt * ty1
                * ( r43*c34 - c1345 ) * tmp2 * u[ 2][ i][ j+1][k];
               b[ 4][ 3][ i][j] =  dt * ty2
                * ( - c2 * ( u[ 2][ i][ j+1][k]*u[ 3][ i][ j+1][k] ) * tmp2 )
                - dt * ty1 * ( c34 - c1345 ) * tmp2 * u[ 3][ i][ j+1][k];
               b[ 4][ 4][ i][j] =  dt * ty2
                * ( c1 * ( u[ 2][ i][ j+1][k] * tmp1 ) )
                - dt * ty1 * c1345 * tmp1
                - dt * ty1 * dy5;

//---------------------------------------------------------------------
//   form the third block sub-diagonal
//---------------------------------------------------------------------
               tmp1 = rho_i[ i][ j][k+1];
               tmp2 = tmp1 * tmp1;
               tmp3 = tmp1 * tmp2;

               c[ 0][ 0][ i][j] = - dt * tz1 * dz1;
               c[ 0][ 1][ i][j] =   0.0;
               c[ 0][ 2][ i][j] =   0.0;
               c[ 0][ 3][ i][j] = dt * tz2;
               c[ 0][ 4][ i][j] =   0.0;

               c[ 1][ 0][ i][j] = dt * tz2
                 * ( - ( u[ 1][ i][ j][k+1]*u[ 3][ i][ j][k+1] ) * tmp2 )
                 - dt * tz1 * ( - c34 * tmp2 * u[ 1][ i][ j][k+1] );
               c[ 1][ 1][ i][j] = dt * tz2 * ( u[ 3][ i][ j][k+1] * tmp1 )
                 - dt * tz1 * c34 * tmp1
                 - dt * tz1 * dz2 ;
               c[ 1][ 2][ i][j] = 0.0;
               c[ 1][ 3][ i][j] = dt * tz2 * ( u[ 1][ i][ j][k+1] * tmp1 );
               c[ 1][ 4][ i][j] = 0.0;

               c[ 2][ 0][ i][j] = dt * tz2
                 * ( - ( u[ 2][ i][ j][k+1]*u[ 3][ i][ j][k+1] ) * tmp2 )
                 - dt * tz1 * ( - c34 * tmp2 * u[ 2][ i][ j][k+1] );
               c[ 2][ 1][ i][j] = 0.0;
               c[ 2][ 2][ i][j] = dt * tz2 * ( u[ 3][ i][ j][k+1] * tmp1 )
                 - dt * tz1 * ( c34 * tmp1 )
                 - dt * tz1 * dz3;
               c[ 2][ 3][ i][j] = dt * tz2 * ( u[ 2][ i][ j][k+1] * tmp1 );
               c[ 2][ 4][ i][j] = 0.0;

               c[ 3][ 0][ i][j] = dt * tz2
              * ( - pow2(( u[ 3][ i][ j][k+1] * tmp1 ))
                  + c2 * ( qs[ i][ j][k+1] * tmp1 ) )
              - dt * tz1 * ( - r43 * c34 * tmp2 * u[ 3][ i][ j][k+1] );
               c[ 3][ 1][ i][j] = dt * tz2
                   * ( - c2 * ( u[ 1][ i][ j][k+1] * tmp1 ) );
               c[ 3][ 2][ i][j] = dt * tz2
                   * ( - c2 * ( u[ 2][ i][ j][k+1] * tmp1 ) );
               c[ 3][ 3][ i][j] = dt * tz2 * ( 2.0 - c2 )
                   * ( u[ 3][ i][ j][k+1] * tmp1 )
                   - dt * tz1 * ( r43 * c34 * tmp1 )
                   - dt * tz1 * dz4;
               c[ 3][ 4][ i][j] = dt * tz2 * c2;

               c[ 4][ 0][ i][j] = dt * tz2
           * ( ( c2 * 2.0 * qs[ i][ j][k+1]
             - c1 * u[ 4][ i][ j][k+1] )
                  * ( u[ 3][ i][ j][k+1] * tmp2 ) )
             - dt * tz1
             * ( - ( c34 - c1345 ) * tmp3 * pow2(u[ 1][ i][ j][k+1])
                 - ( c34 - c1345 ) * tmp3 * pow2(u[ 2][ i][ j][k+1])
                 - ( r43*c34 - c1345 )* tmp3 * pow2(u[ 3][ i][ j][k+1])
                - c1345 * tmp2 * u[ 4][ i][ j][k+1] );
               c[ 4][ 1][ i][j] = dt * tz2
             * ( - c2 * ( u[ 1][ i][ j][k+1]*u[ 3][ i][ j][k+1] ) * tmp2 )
             - dt * tz1 * ( c34 - c1345 ) * tmp2 * u[ 1][ i][ j][k+1];
               c[ 4][ 2][ i][j] = dt * tz2
             * ( - c2 * ( u[ 2][ i][ j][k+1]*u[ 3][ i][ j][k+1] ) * tmp2 )
             - dt * tz1 * ( c34 - c1345 ) * tmp2 * u[ 2][ i][ j][k+1];
               c[ 4][ 3][ i][j] = dt * tz2
             * ( c1 * ( u[ 4][ i][ j][k+1] * tmp1 )
             - c2
             * ( qs[ i][ j][k+1] * tmp1
                  + u[ 3][ i][ j][k+1]*u[ 3][ i][ j][k+1] * tmp2 ) )
             - dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u[ 3][ i][ j][k+1];
               c[ 4][ 4][ i][j] = dt * tz2
             * ( c1 * ( u[ 3][ i][ j][k+1] * tmp1 ) )
             - dt * tz1 * c1345 * tmp1
             - dt * tz1 * dz5;

            }
         } 
  }
  
  public void l2norm(int ldx, int ldy, int ldz, int nx0, int ny0, int nz0, 
    		     int ist, int iend, 
        	     int jst, int jend, double[][][][] v, double[] sum){
    int i, j, k, m;

    for(m=0;m<=4;m++){
       sum[m] = 0.0;
    }

    for(k=1;k<=nz0-2;k++){
       for(j=jst-1;j<=jend-1;j++){
    	  for(i=ist-1;i<=iend-1;i++){
    	     for(m=0;m<=4;m++){
    		sum[m] = sum[m] + v[ m][ i][ j][k] 
        			  * v[ m][ i][ j][k];
    	     }
    	  }
       }
    }

    for(m=0;m<=4;m++){
       sum[m] = Math.Sqrt( sum[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
    }
  }

  public void pintgr(){
      int i, j, k;
      int ibeg, ifin, ifin1;
      int jbeg, jfin, jfin1;
      double[][]  phi1 = instantiate_jagged_array_2(isiz2+2, isiz3+2); 
      double[][]  phi2 = instantiate_jagged_array_2(isiz2+2, isiz3+2);
      double  frc1, frc2, frc3;
      //int isize5 = (isiz2+2);

//---------------------------------------------------------------------
//   set up the sub-domains for intation in each processor
//---------------------------------------------------------------------
      ibeg = ii1;
      ifin = ii2;
      jbeg = ji1;
      jfin = ji2;
      ifin1 = ifin - 1;
      jfin1 = jfin - 1;

//---------------------------------------------------------------------
//   initialize
//---------------------------------------------------------------------
      for(i=0;i<=isiz2+1;i++){
        for(k=0;k<=isiz3+1;k++){
          phi1[k][i] = 0;
          phi2[k][i] = 0;
        }
      }

      for(j=jbeg-1;j<=jfin-1;j++){
         for(i=ibeg-1;i<=ifin-1;i++){

            k = ki1-1;

            phi1[j][i] = c2*(  u[ 4][ i][ j][k]
                 - 0.50 * (  pow2(u[ 1][ i][ j][k])
                               + pow2(u[ 2][ i][ j][k])
                               + pow2(u[ 3][ i][ j][k]) )
                              / u[ 0][ i][ j][k] );

            k = ki2-1;

            phi2[j][i] = c2*(  u[ 4][ i][ j][k]
                 - 0.50 * (  pow2(u[ 1][ i][ j][k])
                               + pow2(u[ 2][ i][ j][k])
                               + pow2(u[ 3][ i][ j][k]) )
                              / u[ 0][ i][ j][k] );
         }
      }

      frc1 = 0.0;

      for(j=jbeg-1;j<=jfin1-1;j++){
         for(i=ibeg-1;i<=ifin1-1;i++){
            frc1 = frc1 + (  phi1[j][i]
                           + phi1[j][i+1]
                           + phi1[j+1][i]
                           + phi1[j+1][i+1]
                           + phi2[j][i]
                           + phi2[j][i+1]
                           + phi2[j+1][i]
                           + phi2[j+1][i+1] );
         }
      }

      frc1 = dxi * deta * frc1;

//---------------------------------------------------------------------
//   initialize
//---------------------------------------------------------------------
      for(i=0;i<=isiz2+1;i++){
        for(k=0;k<=isiz3+1;k++){
          phi1[k][i] = 0;
          phi2[k][i] = 0;
        }
      }
      if (jbeg==ji1) {
        for(k=ki1-1;k<=ki2-1;k++){
           for(i=ibeg-1;i<=ifin-1;i++){
              phi1[k][i] = c2*(  u[ 4][ i][ jbeg-1][k]
                   - 0.50 * (  pow2(u[ 1][ i][ jbeg-1][k])
                                 + pow2(u[ 2][ i][ jbeg-1][k])
                                 + pow2(u[ 3][ i][ jbeg-1][k]) )
                                / u[ 0][ i][ jbeg-1][k] );
           }
        }
      }

      if (jfin==ji2) {
        for(k=ki1-1;k<=ki2-1;k++){
           for(i=ibeg-1;i<=ifin-1;i++){
              phi2[k][i] = c2*(  u[ 4][ i][ jfin-1][k]
                   - 0.50 * (  pow2(u[ 1][ i][ jfin-1][k])
                                 + pow2(u[ 2][ i][ jfin-1][k])
                                 + pow2(u[ 3][ i][ jfin-1][k]) )
                                / u[ 0][ i][ jfin-1][k] );
           }
        }
      }

      frc2 = 0.0;
      for(k=ki1-1;k<=ki2-2;k++){
         for(i=ibeg-1;i<=ifin1-1;i++){
            frc2 = frc2 + (  phi1[k][i]
                           + phi1[k][i+1]
                           + phi1[k+1][i]
                           + phi1[k+1][i+1]
                           + phi2[k][i]
                           + phi2[k][i+1]
                           + phi2[k+1][i]
                           + phi2[k+1][i+1] );
         }
      }

      frc2 = dxi * dzeta * frc2;

//---------------------------------------------------------------------
//   initialize
//---------------------------------------------------------------------
      for(i=0;i<=isiz2+1;i++){
        for(k=0;k<=isiz3+1;k++){
          phi1[k][i] = 0;
          phi2[k][i] = 0;
        }
      }
      if (ibeg==ii1) {
        for(k=ki1-1;k<=ki2-1;k++){
           for(j=jbeg-1;j<=jfin-1;j++){
              phi1[k][j] = c2*(  u[ 4][ ibeg-1][ j][k]
                   - 0.50 * ( pow2(u[ 1][ ibeg-1][ j][k])
                                 + pow2(u[ 2][ ibeg-1][ j][k])
                                 + pow2(u[ 3][ ibeg-1][ j][k]) )
                                / u[ 0][ ibeg-1][ j][k] );
           }
        }
      }

      if (ifin==ii2) {
        for(k=ki1-1;k<=ki2-1;k++){
           for(j=jbeg-1;j<=jfin-1;j++){
              phi2[k][j] = c2*(  u[ 4][ ifin-1][ j][k]
                   - 0.50 * (  pow2(u[ 1][ ifin-1][ j][k] )
                                 + pow2(u[ 2][ ifin-1][ j][k])
                                 + pow2(u[ 3][ ifin-1][ j][k]) )
                                / u[ 0][ ifin-1][ j][k] );
           }
        }
      }

      frc3 = 0.0;

      for(k=ki1-1;k<=ki2-2;k++){
         for(j=jbeg-1;j<=jfin1-1;j++){
            frc3 = frc3 + (  phi1[k][j]
                           + phi1[k][j+1]
                           + phi1[k+1][j]
                           + phi1[k+1][j+1]
                           + phi2[k][j]
                           + phi2[k][j+1]
                           + phi2[k+1][j]
                           + phi2[k+1][j+1] );
         }
      }
      frc3 = deta * dzeta * frc3;
      frc = 0.25 * ( frc1 + frc2 + frc3 );
  }
  
  public void getInputPars(){
      
//---------------------------------------------------------------------
//    if input file does not exist, it uses defaults
//       ipr = 1 for detailed progress output
//       inorm = how often the norm is printed (once every inorm iterations)
//       itmax = number of pseudo time steps
//       dt = time step
//       omega 1 over-relaxation factor for SSOR
//       tolrsd = steady state residual tolerance levels
//       nx, ny, nz = number of grid points in x, y, z directions
//---------------------------------------------------------------------
     // File f2 = new File("inputlu.data");
      if (File.Exists("inputlu.data")){
				
		FileStream f2 = new FileStream("inputlu.data", System.IO.FileMode.Open);
		try{  
  	  StreamReader datafile = new StreamReader(f2);
	  Console.WriteLine("Reading from input file inputlu.data");
	  
	  ipr = int.Parse(datafile.ReadLine());
	  inorm = int.Parse(datafile.ReadLine());
	  itmax = int.Parse(datafile.ReadLine());
	  dt = double.Parse(datafile.ReadLine());
	  omega = double.Parse(datafile.ReadLine());
	  tolrsd[0] = double.Parse(datafile.ReadLine());
	  tolrsd[1] = double.Parse(datafile.ReadLine());
	  tolrsd[2] = double.Parse(datafile.ReadLine());
	  tolrsd[3] = double.Parse(datafile.ReadLine());
	  tolrsd[4] = double.Parse(datafile.ReadLine());
	  nx0 = int.Parse(datafile.ReadLine());
	  ny0 = int.Parse(datafile.ReadLine()); 
	  nz0 = int.Parse(datafile.ReadLine());
	  
	 // fis.close();
        }catch(Exception e){  
	       Console.WriteLine("exception caught! " + e.Message);
        } 
      }else{
        ipr = ipr_default;
        inorm = inorm_default;
        itmax = itmax_default;
        dt = dt_default;
        omega = omega_default;
        tolrsd[0] = tolrsd1_def;
        tolrsd[1] = tolrsd2_def;
        tolrsd[2] = tolrsd3_def;
        tolrsd[3] = tolrsd4_def;
        tolrsd[4] = tolrsd5_def;
        nx0 = isiz1;
        ny0 = isiz2;
        nz0 = isiz3;
      }

//---------------------------------------------------------------------
//   check problem size
//---------------------------------------------------------------------
      if ( ( nx0 < 4 ) || ( ny0 < 4 ) || ( nz0 < 4 ) ) {
	Console.WriteLine("     PROBLEM SIZE IS TOO SMALL - ");
	Console.WriteLine("     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5");
	Environment.Exit(0);
      }
      if ( ( nx0 > isiz1 ) || ( ny0 > isiz2 ) || ( nz0 > isiz3 ) ) {
	Console.WriteLine("     PROBLEM SIZE IS TOO LARGE - ");
	Console.WriteLine("     NX, NY AND NZ SHOULD BE EQUAL TO");
	Console.WriteLine("     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY");
	Environment.Exit(0);
      }
      Console.WriteLine("LU: Iterations="+itmax+" dt="+dt);
  }
  
		
  public void rhs()
		{
      int i, j, k, m;
      double  q;
      double  tmp;
      double  u21, u31, u41;
      double  u21i, u31i, u41i, u51i;
      double  u21j, u31j, u41j, u51j;
      double  u21k, u31k, u41k, u51k;
      double  u21im1, u31im1, u41im1, u51im1;
      double  u21jm1, u31jm1, u41jm1, u51jm1;
      double  u21km1, u31km1, u41km1, u51km1;

      for(k=0;k<=nz-1;k++){
         for(j=0;j<=ny-1;j++){
            for(i=0;i<=nx-1;i++){
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] = - frct[ m][ i][ j][k];
               }
               tmp = 1.0 / u[ 0][ i][ j][k];
               rho_i[ i][ j][k] = tmp;
               qs[ i][ j][k] = 0.50 * (  
	                         u[ 1][ i][ j][k] * u[ 1][ i][ j][k]
                               + u[ 2][ i][ j][k] * u[ 2][ i][ j][k]
                               + u[ 3][ i][ j][k] * u[ 3][ i][ j][k] )
                            * tmp;
            }
         }
      }
      if (timeron)  timer.start(t_rhsx);
//---------------------------------------------------------------------
//   xi-direction flux differences
//---------------------------------------------------------------------

      for(k=1;k<=nz - 2;k++){
         for(j=jst-1;j<=jend-1;j++){
            for(i=0;i<=nx-1;i++){
               flux[ 0][i] = u[ 1][ i][ j][k];
               u21 = u[ 1][ i][ j][k] * rho_i[ i][ j][k];

               q = qs[ i][ j][k];

               flux[ 1][i] = u[ 1][ i][ j][k] * u21 + c2 * 
                              ( u[ 4][ i][ j][k] - q );
               flux[ 2][i] = u[ 2][ i][ j][k] * u21;
               flux[ 3][i] = u[ 3][ i][ j][k] * u21;
               flux[ 4][i] = ( c1 * u[ 4][ i][ j][k] - c2 * q ) * u21;
            }

            for(i=ist-1;i<=iend-1;i++){
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] =  rsd[ m][ i][ j][k]
                       - tx2 * ( flux[ m][i+1] - flux[ m][i-1] );
               }
            }

            for(i=ist-1;i<=nx-1;i++){
               tmp = rho_i[ i][ j][k];

               u21i = tmp * u[ 1][ i][ j][k];
               u31i = tmp * u[ 2][ i][ j][k];
               u41i = tmp * u[ 3][ i][ j][k];
               u51i = tmp * u[ 4][ i][ j][k];

               tmp = rho_i[ i-1][ j][k];

               u21im1 = tmp * u[ 1][ i-1][ j][k];
               u31im1 = tmp * u[ 2][ i-1][ j][k];
               u41im1 = tmp * u[ 3][ i-1][ j][k];
               u51im1 = tmp * u[ 4][ i-1][ j][k];

               flux[ 1][i] = (4.0/3.0) * tx3 * (u21i-u21im1);
               flux[ 2][i] = tx3 * ( u31i - u31im1 );
               flux[ 3][i] = tx3 * ( u41i - u41im1 );
               flux[ 4][i] = 0.50 * ( 1.0 - c1*c5 )
                    * tx3 * ( ( pow2(u21i) + pow2(u31i) + pow2(u41i) )
                            - ( pow2(u21im1) + pow2(u31im1) + pow2(u41im1) ) )
                    + (1.0/6.0)
                    * tx3 * ( pow2(u21i) - pow2(u21im1) )
                    + c1 * c5 * tx3 * ( u51i - u51im1 );
            }

            for(i=ist-1;i<=iend-1;i++){
               rsd[ 0][ i][ j][k] = rsd[ 0][ i][ j][k]
                    + dx1 * tx1 * (            u[ 0][ i-1][ j][k]
                                   - 2.0 * u[ 0][ i][ j][k]
                                   +           u[ 0][ i+1][ j][k] );
               rsd[ 1][ i][ j][k] = rsd[ 1][ i][ j][k]
                + tx3 * c3 * c4 * ( flux[ 1][i+1] - flux[ 1][i] )
                    + dx2 * tx1 * (            u[ 1][ i-1][ j][k]
                                   - 2.0 * u[ 1][ i][ j][k]
                                   +           u[ 1][ i+1][ j][k] );
               rsd[ 2][ i][ j][k] = rsd[ 2][ i][ j][k]
                + tx3 * c3 * c4 * ( flux[ 2][i+1] - flux[ 2][i] )
                    + dx3 * tx1 * (            u[ 2][ i-1][ j][k]
                                   - 2.0 * u[ 2][ i][ j][k]
                                   +           u[ 2][ i+1][ j][k] );
               rsd[ 3][ i][ j][k] = rsd[ 3][ i][ j][k]
                + tx3 * c3 * c4 * ( flux[ 3][i+1] - flux[ 3][i] )
                    + dx4 * tx1 * (            u[ 3][ i-1][ j][k]
                                   - 2.0 * u[ 3][ i][ j][k]
                                   +           u[ 3][ i+1][ j][k] );
               rsd[ 4][ i][ j][k] = rsd[ 4][ i][ j][k]
                + tx3 * c3 * c4 * ( flux[ 4][i+1] - flux[ 4][i] )
                    + dx5 * tx1 * (            u[ 4][ i-1][ j][k]
                                   - 2.0 * u[ 4][ i][ j][k]
                                   +           u[ 4][ i+1][ j][k] );
            }

//---------------------------------------------------------------------
//   Fourth-order dissipation
//---------------------------------------------------------------------
            for(m=0;m<=4;m++){
               rsd[ m][ 1][ j][k] = rsd[ m][ 1][ j][k]
                 - dssp * ( + 5.0 * u[ m][ 1][ j][k]
                            - 4.0 * u[ m][ 2][ j][k]
                            +           u[ m][ 3][ j][k] );
               rsd[ m][ 2][ j][k] = rsd[ m][ 2][ j][k]
                 - dssp * ( - 4.0 * u[ m][ 1][ j][k]
                            + 6.0 * u[ m][ 2][ j][k]
                            - 4.0 * u[ m][ 3][ j][k]
                            +           u[ m][ 4][ j][k] );
            }

            for(i=3;i<=nx - 4;i++){
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] = rsd[ m][ i][ j][k]
                    - dssp * (            u[ m][ i-2][ j][k]
                              - 4.0 * u[ m][ i-1][ j][k]
                              + 6.0 * u[ m][ i][ j][k]
                              - 4.0 * u[ m][ i+1][ j][k]
                              +           u[ m][ i+2][ j][k] );
               }
            }

	    
            for(m=0;m<=4;m++){
               rsd[ m][ nx-3][ j][k] = rsd[ m][ nx-3][ j][k]
                 - dssp * (             u[ m][ nx-5][ j][k]
                            - 4.0 * u[ m][ nx-4][ j][k]
                            + 6.0 * u[ m][ nx-3][ j][k]
                            - 4.0 * u[ m][ nx-2][ j][k]  );
               rsd[ m][ nx-2][ j][k] = rsd[ m][ nx-2][ j][k]
                 - dssp * (             u[ m][ nx-4][ j][k]
                            - 4.0 * u[ m][ nx-3][ j][k]
                            + 5.0 * u[ m][ nx-2][ j][k] );
            }
	    
         }
      }
      if (timeron)  timer.stop(t_rhsx);

      if (timeron)  timer.start(t_rhsy);

//---------------------------------------------------------------------
//   eta-direction flux differences
//---------------------------------------------------------------------
      for(k=1;k<=nz - 2;k++){
         for(i=ist-1;i<=iend-1;i++){
            for(j=0;j<=ny-1;j++){
               flux[ 0][j] = u[ 2][ i][ j][k];
               u31 = u[ 2][ i][ j][k] * rho_i[ i][ j][k];

               q = qs[ i][ j][k];

               flux[ 1][j] = u[ 1][ i][ j][k] * u31 ;
               flux[ 2][j] = u[ 2][ i][ j][k] * u31 + c2 * (u[ 4][ i][ j][k]-q);
               flux[ 3][j] = u[ 3][ i][ j][k] * u31;
               flux[ 4][j] = ( c1 * u[ 4][ i][ j][k] - c2 * q ) * u31;
            }

            for(j=jst-1;j<=jend-1;j++){
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] =  rsd[ m][ i][ j][k]
                         - ty2 * ( flux[ m][j+1] - flux[ m][j-1] );
               }
            }

            for(j=jst-1;j<=ny-1;j++){
               tmp = rho_i[ i][ j][k];

               u21j = tmp * u[ 1][ i][ j][k];
               u31j = tmp * u[ 2][ i][ j][k];
               u41j = tmp * u[ 3][ i][ j][k];
               u51j = tmp * u[ 4][ i][ j][k];

               tmp = rho_i[ i][ j-1][k];
               u21jm1 = tmp * u[ 1][ i][ j-1][k];
               u31jm1 = tmp * u[ 2][ i][ j-1][k];
               u41jm1 = tmp * u[ 3][ i][ j-1][k];
               u51jm1 = tmp * u[ 4][ i][ j-1][k];

               flux[ 1][j] = ty3 * ( u21j - u21jm1 );
               flux[ 2][j] = (4.0/3.0) * ty3 * (u31j-u31jm1);
               flux[ 3][j] = ty3 * ( u41j - u41jm1 );
               flux[ 4][j] = 0.50 * ( 1.0 - c1*c5 )
                    * ty3 * ( ( pow2(u21j) + pow2(u31j) + pow2(u41j) )
                            - ( pow2(u21jm1) + pow2(u31jm1) + pow2(u41jm1) ) )
                    + (1.0/6.0)
                    * ty3 * ( pow2(u31j) - pow2(u31jm1) )
                    + c1 * c5 * ty3 * ( u51j - u51jm1 );
            }

            for(j=jst-1;j<=jend-1;j++){

               rsd[ 0][ i][ j][k] = rsd[ 0][ i][ j][k]
                    + dy1 * ty1 * (            u[ 0][ i][ j-1][k]
                                   - 2.0 * u[ 0][ i][ j][k]
                                   +           u[ 0][ i][ j+1][k] );

               rsd[ 1][ i][ j][k] = rsd[ 1][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 1][j+1] - flux[ 1][j] )
                    + dy2 * ty1 * (            u[ 1][ i][ j-1][k]
                                   - 2.0 * u[ 1][ i][ j][k]
                                   +           u[ 1][ i][ j+1][k] );

               rsd[ 2][ i][ j][k] = rsd[ 2][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 2][j+1] - flux[ 2][j] )
                    + dy3 * ty1 * (            u[ 2][ i][ j-1][k]
                                   - 2.0 * u[ 2][ i][ j][k]
                                   +           u[ 2][ i][ j+1][k] );

               rsd[ 3][ i][ j][k] = rsd[ 3][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 3][j+1] - flux[ 3][j] )
                    + dy4 * ty1 * (            u[ 3][ i][ j-1][k]
                                   - 2.0 * u[ 3][ i][ j][k]
                                   +           u[ 3][ i][ j+1][k] );

               rsd[ 4][ i][ j][k] = rsd[ 4][ i][ j][k]
                + ty3 * c3 * c4 * ( flux[ 4][j+1] - flux[ 4][j] )
                    + dy5 * ty1 * (            u[ 4][ i][ j-1][k]
                                   - 2.0 * u[ 4][ i][ j][k]
                                   +           u[ 4][ i][ j+1][k] );

            }

//---------------------------------------------------------------------
//   fourth-order dissipation
//---------------------------------------------------------------------
            for(m=0;m<=4;m++){
               rsd[ m][ i][ 1][k] = rsd[ m][ i][ 1][k]
                 - dssp * ( + 5.0 * u[ m][ i][ 1][k]
                            - 4.0 * u[ m][ i][ 2][k]
                            +           u[ m][ i][ 3][k] );
               rsd[ m][ i][ 2][k] = rsd[ m][ i][ 2][k]
                 - dssp * ( - 4.0 * u[ m][ i][ 1][k]
                            + 6.0 * u[ m][ i][ 2][k]
                            - 4.0 * u[ m][ i][ 3][k]
                            +           u[ m][ i][ 4][k] );
            }

            for(j=3;j<=ny - 4;j++){
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] = rsd[ m][ i][ j][k]
                    - dssp * (            u[ m][ i][ j-2][k]
                              - 4.0 * u[ m][ i][ j-1][k]
                              + 6.0 * u[ m][ i][ j][k]
                              - 4.0 * u[ m][ i][ j+1][k]
                              +           u[ m][ i][ j+2][k] );
               }
            }

            for(m=0;m<=4;m++){
               rsd[ m][ i][ ny-3][k] = rsd[ m][ i][ ny-3][k]
                 - dssp * (             u[ m][ i][ ny-5][k]
                            - 4.0 * u[ m][ i][ ny-4][k]
                            + 6.0 * u[ m][ i][ ny-3][k]
                            - 4.0 * u[ m][ i][ ny-2][k]  );
               rsd[ m][ i][ ny-2][k] = rsd[ m][ i][ ny-2][k]
                 - dssp * (             u[ m][ i][ ny-4][k]
                            - 4.0 * u[ m][ i][ ny-3][k]
                            + 5.0 * u[ m][ i][ ny-2][k] );
            }

         }
      }
      if (timeron)  timer.stop(t_rhsy);

      if (timeron)  timer.start(t_rhsz);

//---------------------------------------------------------------------
//   zeta-direction flux differences
//---------------------------------------------------------------------
      for(j=jst-1;j<=jend-1;j++){
         for(i=ist-1;i<=iend-1;i++){
      	    for(k=0;k<=nz-1;k++){
               flux[ 0][k] = u[ 3][ i][ j][k];
               u41 = u[ 3][ i][ j][k] * rho_i[ i][ j][k];

               q = qs[ i][ j][k];

               flux[ 1][k] = u[ 1][ i][ j][k] * u41 ;
               flux[ 2][k] = u[ 2][ i][ j][k] * u41 ;
               flux[ 3][k] = u[ 3][ i][ j][k] * u41 + c2 * (u[ 4][ i][ j][k]-q);
               flux[ 4][k] = ( c1 * u[ 4][ i][ j][k] - c2 * q ) * u41;
            }

            for(k=1;k<=nz - 2;k++){
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] =  rsd[ m][ i][ j][k]
                      - tz2 * ( flux[ m][k+1] - flux[ m][k-1] );
               }
            }

            for(k=1;k<=nz-1;k++){
               tmp = rho_i[ i][ j][k];

               u21k = tmp * u[ 1][ i][ j][k];
               u31k = tmp * u[ 2][ i][ j][k];
               u41k = tmp * u[ 3][ i][ j][k];
               u51k = tmp * u[ 4][ i][ j][k];

               tmp = rho_i[ i][ j][k-1];

               u21km1 = tmp * u[ 1][ i][ j][k-1];
               u31km1 = tmp * u[ 2][ i][ j][k-1];
               u41km1 = tmp * u[ 3][ i][ j][k-1];
               u51km1 = tmp * u[ 4][ i][ j][k-1];

               flux[ 1][k] = tz3 * ( u21k - u21km1 );
               flux[ 2][k] = tz3 * ( u31k - u31km1 );
               flux[ 3][k] = (4.0/3.0) * tz3 * (u41k-u41km1);
               flux[ 4][k] = 0.50 * ( 1.0 - c1*c5 )
                    * tz3 * ( (pow2(u21k) + pow2(u31k) +pow2(u41k) )
                            - ( pow2(u21km1) + pow2(u31km1) +pow2(u41km1) ) )
                    + (1.0/6.0)
                    * tz3 * ( pow2(u41k) - pow2(u41km1) )
                    + c1 * c5 * tz3 * ( u51k - u51km1 );
            }

            for(k=1;k<=nz - 2;k++){
               rsd[ 0][ i][ j][k] = rsd[ 0][ i][ j][k]
                    + dz1 * tz1 * (            u[ 0][ i][ j][k-1]
                                   - 2.0 * u[ 0][ i][ j][k]
                                   +           u[ 0][ i][ j][k+1] );
               rsd[ 1][ i][ j][k] = rsd[ 1][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 1][k+1] - flux[ 1][k] )
                    + dz2 * tz1 * (            u[ 1][ i][ j][k-1]
                                   - 2.0 * u[ 1][ i][ j][k]
                                   +           u[ 1][ i][ j][k+1] );
               rsd[ 2][ i][ j][k] = rsd[ 2][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 2][k+1] - flux[ 2][k] )
                    + dz3 * tz1 * (            u[ 2][ i][ j][k-1]
                                   - 2.0 * u[ 2][ i][ j][k]
                                   +           u[ 2][ i][ j][k+1] );
               rsd[ 3][ i][ j][k] = rsd[ 3][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 3][k+1] - flux[ 3][k] )
                    + dz4 * tz1 * (            u[ 3][ i][ j][k-1]
                                   - 2.0 * u[ 3][ i][ j][k]
                                   +           u[ 3][ i][ j][k+1] );
               rsd[ 4][ i][ j][k] = rsd[ 4][ i][ j][k]
                + tz3 * c3 * c4 * ( flux[ 4][k+1] - flux[ 4][k] )
                    + dz5 * tz1 * (            u[ 4][ i][ j][k-1]
                                   - 2.0 * u[ 4][ i][ j][k]
                                   +           u[ 4][ i][ j][k+1] );
            }

//---------------------------------------------------------------------
//   fourth-order dissipation
//---------------------------------------------------------------------
            for(m=0;m<=4;m++){
               rsd[ m][ i][ j][1] = rsd[ m][ i][ j][1]
                 - dssp * ( + 5.0 * u[ m][ i][ j][1]
                            - 4.0 * u[ m][ i][ j][2]
                            +           u[ m][ i][ j][3] );
               rsd[ m][ i][ j][2] = rsd[ m][ i][ j][2]
                 - dssp * ( - 4.0 * u[ m][ i][ j][1]
                            + 6.0 * u[ m][ i][ j][2]
                            - 4.0 * u[ m][ i][ j][3]
                            +       u[ m][ i][ j][4] );
            }

            for(k=3;k<=nz - 4;k++){
               for(m=0;m<=4;m++){
                  rsd[ m][ i][ j][k] = rsd[ m][ i][ j][k]
                    - dssp * (            u[ m][ i][ j][k-2]
                              - 4.0 * u[ m][ i][ j][k-1]
                              + 6.0 * u[ m][ i][ j][k]
                              - 4.0 * u[ m][ i][ j][k+1]
                              +           u[ m][ i][ j][k+2] );
               }
            }

            for(m=0;m<=4;m++){
               rsd[ m][ i][ j][nz-3] = rsd[ m][ i][ j][nz-3]
                 - dssp * (             u[ m][ i][ j][nz-5]
                            - 4.0 * u[ m][ i][ j][nz-4]
                            + 6.0 * u[ m][ i][ j][nz-3]
                            - 4.0 * u[ m][ i][ j][nz-2]  );
               rsd[ m][ i][ j][nz-2] = rsd[ m][ i][ j][nz-2]
                 - dssp * (             u[ m][ i][ j][nz-4]
                            - 4.0 * u[ m][ i][ j][nz-3]
                            + 5.0 * u[ m][ i][ j][nz-2] );
            }
         }
      }
      if (timeron)  timer.stop(t_rhsz);  
  }
  
  public void setcoeff(){
      dxi = 1.0 / ( nx0 - 1 );
      deta = 1.0 / ( ny0 - 1 );
      dzeta = 1.0 / ( nz0 - 1 );

      tx1 = 1.0 / ( dxi * dxi );
      tx2 = 1.0 / ( 2.0 * dxi );
      tx3 = 1.0 / dxi;

      ty1 = 1.0 / ( deta * deta );
      ty2 = 1.0 / ( 2.0 * deta );
      ty3 = 1.0 / deta;

      tz1 = 1.0 / ( dzeta * dzeta );
      tz2 = 1.0 / ( 2.0 * dzeta );
      tz3 = 1.0 / dzeta;
      
//---------------------------------------------------------------------
//   diffusion coefficients
//---------------------------------------------------------------------
      dx1 = 0.75;
      dx2 = dx1;
      dx3 = dx1;
      dx4 = dx1;
      dx5 = dx1;

      dy1 = 0.75;
      dy2 = dy1;
      dy3 = dy1;
      dy4 = dy1;
      dy5 = dy1;

      dz1 = 1.00;
      dz2 = dz1;
      dz3 = dz1;
      dz4 = dz1;
      dz5 = dz1;

//---------------------------------------------------------------------
//   fourth difference dissipation
//---------------------------------------------------------------------
      dssp = ( max (dx1, dy1, dz1 ) ) / 4.0;
  }
   
  public void setbv(){
      int i, j, k, m;
      double[] temp1 = new double[5], temp2 = new double[5];

//---------------------------------------------------------------------
//   set the dependent variable values along the top and bottom faces
//---------------------------------------------------------------------
      for(j=0;j<=ny-1;j++){
         for(i=0;i<=nx-1;i++){
	    exact( i+1, j+1, 1, temp1 );
            exact( i+1, j+1, nz, temp2 );
	    for(m=0;m<=4;m++){
               u[ m][ i][ j][0] = temp1[m];
               u[ m][ i][ j][nz-1] = temp2[m];
           }
         }
      }

//---------------------------------------------------------------------
//   set the dependent variable values along north and south faces
//---------------------------------------------------------------------
      for(k=0;k<=nz-1;k++){
         for(i=0;i<=nx-1;i++){
            exact( i+1, 1, k+1, temp1 );
            exact( i+1, ny, k+1, temp2 );
	    for(m=0;m<=4;m++){
               u[ m][ i][ 0][k] = temp1[m];
               u[ m][ i][ ny-1][k] = temp2[m];
	    }
         }
      }
//---------------------------------------------------------------------
//   set the dependent variable values along east and west faces
//---------------------------------------------------------------------
      for(k=0;k<=nz-1;k++){
         for(j=0;j<=ny-1;j++){
             exact( 1, j+1, k+1, temp1 );
             exact( nx, j+1, k+1, temp2 );
	    for(m=0;m<=4;m++){
               u[ m][ 0][ j][k] = temp1[m];
               u[ m][ nx-1][ j][k] = temp2[m];
           }
         }
      }
  }
  
  public void setiv(){
    int i, j, k, m;
    double  xi, eta, zeta;
    double  pxi, peta, pzeta;
    double[]  ue_1jk = new double[5],
            ue_nx0jk = new double[5],
	    ue_i1k = new double[5],
    	    ue_iny0k = new double[5],
	    ue_ij1 = new double[5],
	    ue_ijnz = new double[5];

    for(k=1;k<=nz-2;k++){
       zeta = (double) k / (nz-1);
       for(j=1;j<=ny-2;j++){
         eta = (double) j / (ny0-1);
    	  for(i=1;i<=nx-2;i++){
            xi = (double) i / (nx0-1);
    	      exact (1,j+1,k+1,ue_1jk);
    	      exact (nx0,j+1,k+1,ue_nx0jk);
    	      exact (i+1,1,k+1,ue_i1k);
    	      exact (i+1,ny0,k+1,ue_iny0k);
    	      exact (i+1,j+1,1,ue_ij1);
    	      exact (i+1,j+1,nz,ue_ijnz);
    	     for(m=0;m<=4;m++){
    		pxi =	( 1.0 - xi ) * ue_1jk[m]
    				  + xi   * ue_nx0jk[m];
    		peta =  ( 1.0 - eta ) * ue_i1k[m]
    				  + eta   * ue_iny0k[m];
    		pzeta = ( 1.0 - zeta ) * ue_ij1[m]
    				  + zeta   * ue_ijnz[m];

    		u[ m][ i][ j][k] = pxi + peta + pzeta
    		     - pxi * peta - peta * pzeta - pzeta * pxi
    		     + pxi * peta * pzeta;

    	     }
    	  }
       }
    }
  }
  
  public double sssor(){
    int i, j, k, m, n;
    int istep;
    double  tmp;
    double[]  delunm = new double[5]; double[][][] tv = instantiate_jagged_array_3(5, isiz1, isiz2);
 
//---------------------------------------------------------------------
//   begin pseudo-time stepping iterations
//---------------------------------------------------------------------
    tmp = 1.0 / ( omega * ( 2.0 - omega ) ); 
//---------------------------------------------------------------------
//   initialize a,b,c,d to zero (guarantees that page tables have been
//   formed, if applicable on given architecture, before timestepping).
//---------------------------------------------------------------------
    for(j=0;j<=isiz2-1;j++){
      for(i=0;i<=isiz1-1;i++){
	for(n=0;n<=4;n++){
	  for(m=0;m<=4;m++){
	    a[m][ n ][ i][j] = 0;
	    b[m][ n ][ i][j] = 0;
	    c[m][ n ][ i][j] = 0;
	    d[m][ n ][ i][j] = 0;
	  }
	}
      }
    }

    timer.resetAllTimers();
//---------------------------------------------------------------------
//   compute the steady-state residuals
//---------------------------------------------------------------------
    rhs();

//---------------------------------------------------------------------
//   compute the L2 norms of newton iteration residuals
//---------------------------------------------------------------------
    l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
	    ist, iend, jst, jend,
	    rsd, rsdnm ); 
 
    timer.resetAllTimers();
    
    timer.start(1);
 
//---------------------------------------------------------------------
//   the timestep loop
//---------------------------------------------------------------------
    for(istep=1;istep<=itmax;istep++){         
      if (istep % 20 == 0 || istep == itmax || istep == 1) {
	Console.WriteLine(" Time step " + istep);
      }

//---------------------------------------------------------------------
//   perform SSOR iteration
//---------------------------------------------------------------------
      if (timeron)  timer.start(t_rhs);
      for(k=1;k<=nz - 2;k++){
	for(j=jst-1;j<=jend-1;j++){
	  for(i=ist-1;i<=iend-1;i++){
	    for(m=0;m<=4;m++){
	      rsd[ m][ i][ j][k] = 
	                      dt * rsd[ m][ i][ j][k];
	    }
	  }
	}
      }
      if (timeron)  timer.stop(t_rhs);
 
      for(k=1;k<=nz -2 ;k++){
//---------------------------------------------------------------------
//   form the lower triangular part of the jacobian matrix
//---------------------------------------------------------------------
	if (timeron)  timer.start(t_jacld);
	jacld(k);
	if (timeron)  timer.stop(t_jacld);

//---------------------------------------------------------------------
//   perform the lower triangular solution
//---------------------------------------------------------------------;
	if (timeron)  timer.start(t_blts);
	blts( isiz1, isiz2, isiz3,
	      nx, ny, nz, k,
	      omega,
	      rsd, tv,
	      a, b, c, d,
	      ist, iend, jst, jend, 
	      nx0, ny0 );
	if (timeron)  timer.stop(t_blts);
      }
      for(k = nz-2;k>=1; k--){
//---------------------------------------------------------------------
//   form the strictly upper triangular part of the jacobian matrix
//---------------------------------------------------------------------
	if (timeron)  timer.start(t_jacu);
	jacu(k);
	if (timeron)  timer.stop(t_jacu);

//---------------------------------------------------------------------
//   perform the upper triangular solution
//---------------------------------------------------------------------
	if (timeron)  timer.start(t_buts);
	buts( isiz1, isiz2, isiz3,
	      nx, ny, nz, k,
	      omega,
	      rsd, tv,
	      d, a, b, c,
	      ist, iend, jst, jend,
	      nx0, ny0 );
	if (timeron)  timer.stop(t_buts);
      }
 
//---------------------------------------------------------------------
//   update the variables
//---------------------------------------------------------------------

      if (timeron)  timer.start(t_add);
      for(k=1;k<=nz-2;k++){
	for(j=jst-1;j<=jend-1;j++){
	  for(i=ist-1;i<=iend-1;i++){
	    for(m=0;m<=4;m++){
	      u[ m][ i][ j][k] += 
		   + tmp * rsd[ m][ i][ j][k];
	    }
	  }
	}
      }
      if (timeron)  timer.stop(t_add);
 
//---------------------------------------------------------------------
//   compute the max-norms of newton iteration corrections
//---------------------------------------------------------------------
      if ( istep % inorm  == 0 ){
	if (timeron)timer.start(t_l2norm);
	l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
		ist, iend, jst, jend,
		rsd, delunm );
	if (timeron)  timer.stop(t_l2norm);
      }
 
//---------------------------------------------------------------------
//   compute the steady-state residuals
//---------------------------------------------------------------------
      if (timeron)  timer.start(t_rhs);
      rhs();    
      if (timeron)  timer.stop(t_rhs);
 
//---------------------------------------------------------------------
//   compute the max-norms of newton iteration residuals
//---------------------------------------------------------------------
      if ( istep % inorm == 0 || istep == itmax ){
	if (timeron)timer.start(t_l2norm);
	l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
		ist, iend, jst, jend,
		rsd, rsdnm );
	if (timeron)  timer.stop(t_l2norm);
      }

//---------------------------------------------------------------------
//   check the newton-iteration residuals against the tolerance levels
//---------------------------------------------------------------------
      if ( ( rsdnm[0] < tolrsd[0] ) &&
	   ( rsdnm[1] < tolrsd[1] ) &&
	   ( rsdnm[2] < tolrsd[2] ) &&
	   ( rsdnm[3] < tolrsd[3] ) &&
	   ( rsdnm[4] < tolrsd[4] ) ) {
           timer.stop(1);
           return timer.readTimer(1);
        }
    } 
    timer.stop(1);
    return timer.readTimer(1);
  }
  
  public int verify(double[] xcr, double[] xce, double xci){
    
			double[] xcrref = new double[5], xceref = new double[5]; double xciref=0; 
           double[] xcrdif = new double[5], xcedif = new double[5]; double xcidif=0,
           dtref=0;
    int m;
    int verified=-1;
    char clss = 'U';

    for(m=0;m<=4;m++){
      xcrref[m] = 1.0;
      xceref[m] = 1.0;
    }
    xciref = 1.0;

    if ( (nx0 == 12) && (ny0 == 12) &&(nz0 == 12) && (itmax == 50)){
      clss = 'S';
      dtref = .5;
//---------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (12X12X12) grid,
//   after 50 time steps, with  DT = 5.0d-01
//---------------------------------------------------------------------
      xcrref[0] = 1.6196343210976702E-2;
      xcrref[1] = 2.1976745164821318E-3;
      xcrref[2] = 1.5179927653399185E-3;
      xcrref[3] = 1.5029584435994323E-3;
      xcrref[4] = 3.4264073155896461E-2;

//---------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (12X12X12) grid,
//   after 50 time steps, with  DT = 5.0d-01
//---------------------------------------------------------------------
      xceref[0] = 6.4223319957960924E-4;
      xceref[1] = 8.4144342047347926E-5;
      xceref[2] = 5.8588269616485186E-5;
      xceref[3] = 5.8474222595157350E-5;
      xceref[4] = 1.3103347914111294E-3;

//---------------------------------------------------------------------
//   Reference value of surface integral, for the (12X12X12) grid,
//   after 50 time steps, with DT = 5.0d-01
//---------------------------------------------------------------------
      xciref = 7.8418928865937083;
    }else if ( (nx0 == 33) && (ny0 == 33) &&(nz0 == 33) &&(itmax == 300) ) {

      clss = 'W';
      dtref = 1.5E-3;
//---------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (33x33x33) grid,
//   after 300 time steps, with  DT = 1.5d-3
//---------------------------------------------------------------------
      xcrref[0] =   0.1236511638192E+02;
      xcrref[1] =   0.1317228477799E+01;
      xcrref[2] =   0.2550120713095E+01;
      xcrref[3] =   0.2326187750252E+01;
      xcrref[4] =   0.2826799444189E+02;

//---------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (33X33X33) grid,
//---------------------------------------------------------------------
      xceref[0] =   0.4867877144216;
      xceref[1] =   0.5064652880982E-1;
      xceref[2] =   0.9281818101960E-1;
      xceref[3] =   0.8570126542733E-1;
      xceref[4] =   0.1084277417792E+01;

//---------------------------------------------------------------------
//   Reference value of surface integral, for the (33X33X33) grid,
//   after 300 time steps, with  DT = 1.5d-3
//---------------------------------------------------------------------
      xciref    =   0.1161399311023E+02;
    }else if ( (nx0 == 64) && (ny0 == 64) &&(nz0 == 64) &&(itmax == 250) ) {
      clss = 'A';
      dtref = 2.0;
//---------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (64X64X64) grid,
//   after 250 time steps, with  DT = 2.0
//---------------------------------------------------------------------
      xcrref[0] = 7.7902107606689367E+02;
      xcrref[1] = 6.3402765259692870E+01;
      xcrref[2] = 1.9499249727292479E+02;
      xcrref[3] = 1.7845301160418537E+02;
      xcrref[4] = 1.8384760349464247E+03;

//---------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (64X64X64) grid,
//   after 250 time steps, with  DT = 2.0
//---------------------------------------------------------------------;
      xceref[0] = 2.9964085685471943E+01;
      xceref[1] = 2.8194576365003349;
      xceref[2] = 7.3473412698774742;
      xceref[3] = 6.7139225687777051;
      xceref[4] = 7.0715315688392578E+01;

//---------------------------------------------------------------------
//   Reference value of surface integral, for the (64X64X64) grid,
//   after 250 time steps, with DT = 2.0
//---------------------------------------------------------------------
      xciref = 2.6030925604886277E+01;
    }else if ( (nx0 == 102) && (ny0 == 102) && (nz0 == 102) && (itmax == 250) ) {
      clss = 'B';
      dtref = 2.0;

//---------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (102X102X102) grid,
//   after 250 time steps, with  DT = 2.0
//---------------------------------------------------------------------
      xcrref[0] = 3.5532672969982736E+03;
      xcrref[1] = 2.6214750795310692E+02;
      xcrref[2] = 8.8333721850952190E+02;
      xcrref[3] = 7.7812774739425265E+02;
      xcrref[4] = 7.3087969592545314E+03;

//---------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (102X102X102) 
//   grid, after 250 time steps, with  DT = 2.0
//---------------------------------------------------------------------
      xceref[0] = 1.1401176380212709E+02;
      xceref[1] = 8.1098963655421574;
      xceref[2] = 2.8480597317698308E+01;
      xceref[3] = 2.5905394567832939E+01;
      xceref[4] = 2.6054907504857413E+02;

//---------------------------------------------------------------------
//   Reference value of surface integral, for the (102X102X102) grid,
//   after 250 time steps, with DT = 2.0
//---------------------------------------------------------------------
      xciref = 4.7887162703308227E+01;
    }else if ( (nx0 == 162) && (ny0 == 162) && (nz0 == 162) && (itmax == 250) ) {

      clss = 'C';
      dtref = 2.0;

//---------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (162X162X162) grid,
//   after 250 time steps, with  DT = 2.0
//---------------------------------------------------------------------
      xcrref[0] = 1.03766980323537846E+04;
      xcrref[1] = 8.92212458801008552E+02;
      xcrref[2] = 2.56238814582660871E+03;
      xcrref[3] = 2.19194343857831427E+03;
      xcrref[4] = 1.78078057261061185E+04;

//---------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (162X162X162) 
//   grid, after 250 time steps, with  DT = 2.0
//---------------------------------------------------------------------
      xceref[0] = 2.15986399716949279E+02;
      xceref[1] = 1.55789559239863600E+01;
      xceref[2] = 5.41318863077207766E+01;
      xceref[3] = 4.82262643154045421E+01;
      xceref[4] = 4.55902910043250358E+02;

//---------------------------------------------------------------------
//   Reference value of surface integral, for the (162X162X162) grid,
//   after 250 time steps, with DT = 2.0
//---------------------------------------------------------------------
      xciref = 6.66404553572181300E+01;      
    }

//---------------------------------------------------------------------
//    verification test for residuals if gridsize is either 12X12X12 or 
//    64X64X64 or 102X102X102 or 162X162X162
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//    Compute the difference of solution values and the known reference values.
//---------------------------------------------------------------------
    for(m=0;m<=4;m++){
      xcrdif[m] = Math.Abs((xcr[m]-xcrref[m])/xcrref[m]) ;
      xcedif[m] = Math.Abs((xce[m]-xceref[m])/xceref[m]);       
    }
    xcidif = Math.Abs((xci - xciref)/xciref);

//---------------------------------------------------------------------
//   tolerance level
//---------------------------------------------------------------------
    double epsilon = 1.0E-8;
//---------------------------------------------------------------------
//    Output the comparison of computed results to known cases.
//---------------------------------------------------------------------

    if(clss != 'U') {
      Console.WriteLine(" Verification being performed for class " + clss);
      Console.WriteLine(" Accuracy setting for epsilon = " + epsilon);
      if (Math.Abs(dt-dtref) <= epsilon) {  
	if(verified==-1) verified = 1;
      }else{
	verified = 0;
	clss = 'U';
	Console.WriteLine(" DT= "+dt+
	                   " does not match the reference value of "+ dtref);
      }
    }else{ 
      Console.WriteLine(" Unknown class");
      verified = -1;
    }
    if (clss != 'U') {
      Console.WriteLine(" Comparison of RMS-norms of residual");
    }else{
      Console.WriteLine(" RMS-norms of residual");
      verified = -1;
    }
    verified=BMResults.printComparisonStatus(clss,verified,epsilon,
                                             xcr,xcrref,xcrdif);

    if (clss != 'U') {
      Console.WriteLine(" Comparison of RMS-norms of solution error");
    }else{
      Console.WriteLine(" RMS-norms of solution error");
    }
    verified=BMResults.printComparisonStatus(clss,verified,epsilon,
                                             xce,xceref,xcedif);

    if (clss != 'U') {
      Console.WriteLine(" Comparison of surface integral");
    }else{
      Console.WriteLine(" Surface integral");
    }
    verified=BMResults.printComparisonStatus(clss,verified,epsilon,
                                             xci,xciref,xcidif);
    BMResults.printVerificationStatus(clss,verified,BMName); 
    return verified;
  }
  public void checksum(double[] array, int size, 
                       String[] arrayname, bool stop){
    double sum = 0;
    for(int i=0; i<size; i++){
      sum += array[i];
    }
    Console.WriteLine("array:"+arrayname + " checksum is: " + sum);
    if(stop)Environment.Exit(0);
  }
  public double checkSum(double[] arr){
    double csum=0.0;
    for(int k=0;k<=nz-1;k++){
      for(int j=0;j<=ny-1;j++){
	for(int i=0;i<=nx-1;i++){
	  for(int m=0;m<=4;m++){
	  int offset=m+i*isize1+j*jsize1+k*ksize1;
	    csum+=(arr[offset]*arr[offset])/
	         (double)(nx*ny*nz*5);
	  }
	}
      }
    }
    return csum;
  }  
  public double getTime(){ return timer.readTimer(1);} 
  public void finalize() /*throws Throwable*/{
    Console.WriteLine("LU: is about to be garbage collected"); 
    //super.finalize();
  }
}
}