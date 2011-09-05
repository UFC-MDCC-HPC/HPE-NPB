// LUBase.cs created with MonoDevelop
// User: diane at 21:47 10/4/2009
//
// To change standard headers go to Edit->Preferences->Coding->Standard Headers
//
	/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                               L U B a s e                               !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    LUbase implements base class for LU benchmark.                       !
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
using System.Threading;
using NPB3_0_JAV;

namespace NPB3_0_JAV.LUThreads
{

public class LUBase /* : Thread*/{

  public static String BMName="LU";
  public char CLASS = 'S';
  
  protected int isiz1, isiz2, isiz3;
  protected int itmax_default, inorm_default;
  protected double dt_default;
  protected int ipr_default = 1;
  protected static double 
                   omega_default = 1.2,tolrsd1_def=.00000001, 
                   tolrsd2_def=.00000001, tolrsd3_def=.00000001, 
                   tolrsd4_def=.00000001, tolrsd5_def=.00000001,
		   c1 = 1.4, c2 = 0.4, c3 = .1, c4 = 1, c5 = 1.4;
//---------------------------------------------------------------------
//   grid
//---------------------------------------------------------------------
  protected int nx, ny, nz;
  protected int nx0, ny0, nz0;
  protected int ist, iend;
  protected int jst, jend;
  protected int ii1, ii2;
  protected int ji1, ji2;
  protected int ki1, ki2;
  protected double    dxi, deta, dzeta;
  protected double    tx1, tx2, tx3;
  protected double    ty1, ty2, ty3;
  protected double    tz1, tz2, tz3;

//---------------------------------------------------------------------
//   dissipation
//---------------------------------------------------------------------
  protected double   dx1, dx2, dx3, dx4, dx5;
  protected double   dy1, dy2, dy3, dy4, dy5;
  protected double   dz1, dz2, dz3, dz4, dz5;
  protected double   dssp;

//---------------------------------------------------------------------
//   field variables and residuals
//   to improve cache performance, second two dimensions padded by 1 
//   for even number sizes only.
//   Note: corresponding array (called "v") in routines blts, buts, 
//   and l2norm are similarly padded
//---------------------------------------------------------------------

  protected double[][][][] u, rsd, frct;
  protected int isize1, jsize1, ksize1; 

  protected double[][] flux; 
  protected int isize2;

  protected double[][][] qs, rho_i;
  protected int jsize3, ksize3;
//---------------------------------------------------------------------
//   output control parameters
//---------------------------------------------------------------------
  protected static int ipr, inorm;

//---------------------------------------------------------------------
//   newton-raphson iteration control parameters
//---------------------------------------------------------------------
  protected int itmax;
  protected double dt, omega, frc, ttotal;
  protected double[] tolrsd = new double[5],rsdnm = new double[5], errnm = new double[5];

  protected double[][][][]   a, b, c, d;
  protected int isize4, jsize4, ksize4;
//---------------------------------------------------------------------
//   coefficients of the exact solution
//---------------------------------------------------------------------
   private double[,] ce_ = 			
		{   {2.0, 1.0, 2.0, 2.0, 5.0},
			{0.0, 0.0, 2.0, 2.0, 4.0},
			{0.0, 0.0, 0.0, 0.0, 3.0},
			{4.0, 0.0, 0.0, 0.0, 2.0},
			{5.0, 1.0, 0.0, 0.0, .1 },
			{3.0, 2.0, 2.0, 2.0, .4 },
			{.5 , 3.0, 3.0, 3.0, .3 },
			{.02, .01, .04, .03, .05},
			{.01, .03, .03, .05, .04},
			{.03, .02, .05, .04, .03},
			{.5 , .4 , .3 , .2 , .1} ,
			{.4 , .3 , .5 , .1 , .3} ,
			{.3 , .5 , .4 , .3 , .2}
			};

   protected static double[][] ce = instantiate_jagged_array_2(13,5);

   protected void init_ce() 
   {	
			for (int n=0; n < 13; n++) 
			{
				for (int m=0; m < 5; m++) 
				{
					ce[n][m] = ce_[n,m];
				}	
			}
	}
//---------------------------------------------------------------------
//   timers
//---------------------------------------------------------------------
  public static int t_total = 1, t_rhsx = 2,  t_rhsy = 3, t_rhsz = 4,
               t_rhs = 5, t_jacld = 6, t_blts = 7,t_jacu = 8,
               t_buts = 9, t_add = 10, t_l2norm = 11,t_last = 11;
  public bool timeron;
  public Timer timer = new Timer();

  public LUBase(){ init_ce(); }

  public LUBase(char cls, int np){
    init_ce();
    CLASS=cls;
    num_threads = np;
   switch(cls){
    case 'S':
      isiz1=isiz2=isiz3=12;
      itmax_default=inorm_default=50;
      dt_default=.5;
      break;
    case 'W':
      isiz1=isiz2=isiz3=33;
      itmax_default=inorm_default=300;
      dt_default=.0015;
      break;
    case 'A':
      isiz1=isiz2=isiz3=64;
      itmax_default=inorm_default=250;
      dt_default=2;
      break;
    case 'B':
      isiz1=isiz2=isiz3=102;
      itmax_default=inorm_default=250;
      dt_default=2;
      break;
    case 'C':
      isiz1=isiz2=isiz3=162;
      itmax_default=inorm_default=250;
      dt_default=2;
      break;
    }

    u = instantiate_jagged_array_4(5, isiz1/2*2+1, isiz2/2*2+1, isiz3);
    rsd = instantiate_jagged_array_4(5, isiz1/2*2+1, isiz2/2*2+1, isiz3);
    frct = instantiate_jagged_array_4(5, isiz1/2*2+1, isiz2/2*2+1, isiz3);
    isize1=5;
    jsize1=5*(isiz1/2*2+1);
    ksize1=5*(isiz1/2*2+1)*(isiz2/2*2+1);
    
    flux= instantiate_jagged_array_2(5, isiz1);
    isize2=5;

    qs = instantiate_jagged_array_3(isiz1/2*2+1,isiz2/2*2+1,isiz3);
    rho_i = instantiate_jagged_array_3(isiz1/2*2+1,isiz2/2*2+1,isiz3);
    jsize3 = (isiz1/2*2+1);
    ksize3 = (isiz1/2*2+1)*(isiz2/2*2+1);

    a = instantiate_jagged_array_4(5, 5, isiz1/2*2+1, isiz2); 
    b = instantiate_jagged_array_4(5, 5, isiz1/2*2+1, isiz2);
    c = instantiate_jagged_array_4(5, 5, isiz1/2*2+1, isiz2); 
    d = instantiate_jagged_array_4(5, 5, isiz1/2*2+1, isiz2);

    isize4=5;
    jsize4=5*5;
    ksize4=5*5*(isiz1/2*2+1);
  }
		
	public static double[][][][] instantiate_jagged_array_4(int N1, int N2, int N3, int N4)	
	{
	    double[][][][] r = new double[N1][][][];
		for (int i=0; i < N1; i++) 
		{
			r[i] = new double[N2][][];
			for (int j=0; j<N2; j++) 
			{
				r[i][j] = new double[N3][];
				for (int k=0; k<N3; k++) 
				{
					r[i][j][k] = new double[N4];
				}
			}
		}
			
		return r;
	}
		
	public static double[][][] instantiate_jagged_array_3(int N1, int N2, int N3)	
	{
	    double[][][] r = new double[N1][][];
		for (int i=0; i < N1; i++) 
		{
			r[i] = new double[N2][];
			for (int j=0; j<N2; j++) 
			{
				r[i][j] = new double[N3];
			}
		}
			
		return r;
	}
		
	public static double[][] instantiate_jagged_array_2(int N1, int N2)	
	{
	    double[][] r = new double[N1][];
		for (int i=0; i < N1; i++) 
		{
			r[i] = new double[N2];
		}
			
		return r;
	}
		
  
  protected Thread master=null;
  protected int num_threads;
		
  /*protected RHSCompute[] rhscomputer = null;
  protected Scale[] scaler = null;
  protected Adder[] adder =null;
  protected LowerJac[] lowerjac =null;
  protected UpperJac[] upperjac =null;
  
  public void setupThreads(LU lu){
    master = lu;
   if(num_threads>isiz1-2)
   num_threads=isiz1-2;
		
		
    int[] interval1=new int[num_threads];
    int[] interval2=new int[num_threads];    
	set_interval(num_threads, isiz1, interval1);
    set_interval(num_threads, isiz1-2, interval2);
    int[,] partition1 = new int[interval1.Length][2];
    int[,] partition2 = new int[interval2.Length][2];
    set_partition(0,interval1,partition1);
    set_partition(1,interval2,partition2);
   
    rhscomputer = new RHSCompute[num_threads];
    scaler = new Scale[num_threads];
    adder = new Adder[num_threads];
    lowerjac = new LowerJac[num_threads];
    upperjac = new UpperJac[num_threads];
    for(int ii=0;ii<num_threads;ii++){
      rhscomputer[ii] =  new RHSCompute(lu,partition1[ii][0],partition1[ii][1],
					partition2[ii][0],partition2[ii][1]);
      rhscomputer[ii].id=ii;
      rhscomputer[ii].start();

      scaler[ii] =  new Scale(lu,partition2[ii][0],partition2[ii][1]);
      scaler[ii].id=ii;
      scaler[ii].start();

      adder[ii] =  new Adder(lu,partition2[ii][0],partition2[ii][1]);
      adder[ii].id=ii;
      adder[ii].start();

      lowerjac[ii] =  new LowerJac(lu,partition2[ii][0],partition2[ii][1]);
      lowerjac[ii].id=ii;
      lowerjac[ii].neighbor=lowerjac;
      lowerjac[ii].start();

      upperjac[ii] =  new UpperJac(lu,partition2[num_threads-ii-1][0],
                                      partition2[num_threads-ii-1][1]);
      upperjac[ii].id=ii;
      upperjac[ii].neighbor=upperjac;
      upperjac[ii].start();
    }
 // }*/

  public void checksum(double[] array, int size, 
                       String arrayname, bool stop){
    double sum = 0;
    for(int i=0; i<size; i++) sum += array[i];
    Console.WriteLine("array:"+arrayname + " checksum is: " + sum);
    if(stop) Environment.Exit(0);
  }
  
  public void set_interval(int threads, int problem_size, int[] interval ){
    interval[0]= problem_size/threads;
    for(int i=1;i<threads;i++)
      interval[i]=interval[0];
    int remainder = problem_size%threads;
    for(int i=0;i<remainder;i++)
      interval[i]++;
  } 
  
  public void set_partition(int start, int[] interval, int[][] array){
    array[0][0]=start;
    if(start==0) array[0][1]=interval[0]-1;
    else array[0][1]=interval[0];
    
    for(int i=1;i<interval.Length;i++){
      array[i][0]=array[i-1][1]+1;
      array[i][1]=array[i-1][1]+interval[i];
    }
  }

  public void exact(double i, double j, double k, double[] u0 ){
    double xi  =  (i-1) / ( nx0 - 1 );
    double eta  = (j-1) / ( ny0 - 1 );
    double zeta = (k-1) / ( nz  - 1 );

    for(int m=0;m<=4;m++){
      u0[m] =  ce[0][m]
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
  protected double max(double a, double b){
    if(a<b) return b;
    else return a;
  }
  protected double max(double a, double b, double c){
    return max( a, max( b , c) );
  }
		public double pow2(double p) { return p * p; }
}

}
