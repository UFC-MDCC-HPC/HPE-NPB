/*
!-------------------------------------------------------------------------!
!            P R O G R A M M I N G     B A S E L I N E S                  !
!                                                                         !
!                                 F O R                                   !
!                                                                         !
!      N  A  S     P A R A L L E L     B E N C H M A R K S   P B N 1.0    !
!                                                                         !
!                               F T B a s e                               !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    FTBase implements base class for FT benchmark.                       !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.               !
!    We request, however, that all derived work reference the             !
!    NAS Grid Benchmarks 1.0. This software is provided                   !
!    "as is" without express or implied warranty.                         ! 
!                                                                         !
!    Information on NPB 3.0, including the Technical Report NAS-02-008	  !
!    "Implementation of the NAS Parallel Benchmarks in Java",		  !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!         http://www.nas.nasa.gov/Software/NPB/                           !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!
! Authors: M. Frumkin           					  !
!          M. Schultz           					  !
!-------------------------------------------------------------------------!
*/
using System;
using System.Threading;
using NPB3_0_JAV;

namespace NPB3_0_JAV.FTThreads {

public class FTBase /* : Thread*/
{
  public static String BMName = "FT";
  public char CLASS = 'S';
  
  protected int nx, ny, nz, maxdim, niter_default;

  //complex arrays
  protected double[,,] scr;
  protected double[,,,] plane;
  //protected int isize2;
  //protected int isize3, jsize3, ksize3;
  //protected int isize4, jsize4, ksize4;

  protected double[,] checksum; 
  protected double[,,,] xtr;  //isize3=2;jsize3=2*(ny+1);ksize3=2*(ny+1)*nx;
  protected double[,,,] xnt;  //isize4=2;jsize4=2*(ny+1);ksize4=2*(ny+1)*nz;
  protected double[,] exp1,exp2,exp3; 

  public bool timeron=false;
  public Timer timer = new Timer();

  //constants
  protected static int REAL=0,IMAG=1;
  protected static double pi=Math.PI,alpha=.000001;

  public FTBase(){}
  
  public FTBase(char clss, int np, bool serial){
    CLASS = clss;
    num_threads=np;
    switch (CLASS){
    case 'S':
      nx=ny=nz=64;
      niter_default=6;
      break;
    case 'W':
      nx=ny=128;
      nz=32;
      niter_default=6;
      break;     
    case 'A':
      nx=256;
      ny=256;
      nz=128;
      niter_default=6;
      break;      
    case 'B':
      nx=512;
      ny=nz=256;
      niter_default=20;
      break;
    case 'C':
      nx=ny=nz=512;
      niter_default=20;      
      break;
    }    
    maxdim = max( nx , max( ny , nx ) );
    scr = new double[maxdim,(maxdim+1),2];       // scr = new double[2 * (maxdim + 1) * maxdim];
    plane = new double[1,maxdim,(maxdim+1),2];  // plane = new double[2 * (maxdim + 1) * maxdim];
    //isize2=2;
    //isize3=2;
    //jsize3=2*(ny+1);
    //ksize3=2*(ny+1)*nx;
    //isize4=2;
    //jsize4=2*(ny+1);
    //ksize4=2*(ny+1)*nz;
    //complex values
    checksum = new double[niter_default,2]; //isize2=2;

    xtr = new double[nz,nx,(ny+1),2]; 
    xnt = new double[nz,nx,(ny+1),2]; 
    exp1 = new double[nx,2];
    exp2 = new double[ny,2];
    exp3 = new double[nz,2];
  }

  // thread variables
  //protected Thread master = null;
  protected int num_threads;

  public int max(int a, int b){if(a>b)return a; else return b;}

  public void CompExp (int n, double[,] exponent){     
    int nu = n;
    int m = ilog2(n);
    exponent[0,0] = m;
    
    double eps=1.0E-16;
    int ku = 1;
    int ln = 1;
    for(int j=1;j<=m;j++){
        double t = pi / ln;
        for(int i=0;i<=ln-1;i++){
	        double ti = i * t;
            int idx = (i + ku); // obs  int idx=(i+ku)*2;
	        exponent[idx,REAL] = Math.Cos(ti);
            exponent[idx,IMAG] = Math.Sin(ti);
	        if(Math.Abs(exponent[idx,REAL]) < eps) exponent[idx,REAL]=0;
	        if(Math.Abs(exponent[idx,IMAG]) < eps) exponent[idx,IMAG]=0;
        }
        ku = ku + ln;
        ln = 2 * ln;
    }
  }

  public void initial_conditions(double[,,,] u0, int d1, int d2, int d3){
      double[,] tmp = new double[maxdim,2];
      double[] RanStarts = new double[maxdim];
      //seed has to be init here since
      //is called 2 times 
      double seed = 314159265;
      double a = Math.Pow(5.0, 13);
      double start = seed;
      //---------------------------------------------------------------------
      // Jump to the starting element for our first plane.
      //---------------------------------------------------------------------
      Random rng = new Random(seed);
      double an = rng.ipow46(a, 0);
      rng.randlc(seed, an);
      an = rng.ipow46(a, 2 * d1 * d2);
      //---------------------------------------------------------------------
      // Go through by z planes filling in one square at a time.
      //---------------------------------------------------------------------
      RanStarts[0] = start;
      for (int k = 1; k < d3; k++)
      {
          seed = rng.randlc(start, an);
          RanStarts[k] = start = seed;
      }
      for (int k = 0; k < d3; k++){
          double x0 = RanStarts[k];
          for (int j = 0; j < d1; j++){
              x0 = rng.vranlc2(d2, x0, a, tmp, 0);
              for (int i = 0; i < d2; i++){
                  u0[k,i,j,REAL] = tmp[i,REAL];
                  u0[k,i,j,IMAG] = tmp[i,IMAG];
              }
          }
      }
  }

  public int ilog2(int n){
    int nn, lg;
    if (n == 1) return 0;
    lg = 1;
    nn = 2;
    while (nn < n){
      nn = nn*2;
      lg = lg+1;
    }
   return lg;
  }

  protected static int fftblock_default=4*4096, //Size of L1 cache on SGI O2K
                       fftblock;
  public void Swarztrauber(int is0, int m, int len, int n, double[,,,] x, int xoffst, int xd1, double[,] exponent, double[,,] scr)
  {
      int i, j = 0, l, mx;
      int k, n1, li, lj, lk, ku, i11, i12, i21, i22;
      int BlockStart, BlockEnd;

      //complex values
      double[] u1 = new double[2];
      double[] x11 = new double[2];
      double[] x21 = new double[2];
      if (timeron) timer.start(4);

      //---------------------------------------------------------------------
      //   Perform one variant of the Stockham FFT.
      //---------------------------------------------------------------------

      fftblock = fftblock_default / n;
      if (fftblock < 8) fftblock = 8;
      for (BlockStart = 0; BlockStart < len; BlockStart += fftblock)
      {
          BlockEnd = BlockStart + fftblock - 1;
          if (BlockEnd >= len) BlockEnd = len - 1;
          for (l = 1; l <= m; l += 2)
          {
              n1 = n / 2;
              lk = (int)Math.Pow(2, l - 1);
              li = (int)Math.Pow(2, m - l);
              lj = 2 * lk;
              ku = li;

              for (i = 0; i <= li - 1; i++)
              {
                  i11 = i * lk;
                  i12 = i11 + n1;
                  i21 = i * lj;
                  i22 = i21 + lk;

                  u1[REAL] = exponent[(ku+i), REAL];
                  if (is0 >= 1)
                  {
                      u1[IMAG] = exponent[(ku+i),IMAG];
                  }
                  else
                  {
                      u1[IMAG] = -exponent[(ku+i),IMAG];
                  }
                  for (k = 0; k <= lk - 1; k++)
                  {
                      for (j = BlockStart; j <= BlockEnd; j++)
                      {
                          x11[REAL] = x[xoffst,(i11+k),j,REAL];
                          x11[IMAG] = x[xoffst,(i11+k),j,IMAG];
                          x21[REAL] = x[xoffst,(i12+k),j,REAL];
                          x21[IMAG] = x[xoffst,(i12+k),j,IMAG];
                          scr[(i21+k),j,REAL] = x11[REAL] + x21[REAL];
                          scr[(i21+k),j,IMAG] = x11[IMAG] + x21[IMAG];
                          scr[(i22+k),j,REAL] = u1[REAL] * (x11[REAL] - x21[REAL]) - u1[IMAG] * (x11[IMAG] - x21[IMAG]);
                          scr[(i22+k),j,IMAG] = u1[IMAG] * (x11[REAL] - x21[REAL]) + u1[REAL] * (x11[IMAG] - x21[IMAG]);
                      }
                  }
              }
              if (l == m)
              {
                  for (k = 0; k < n; k++)
                  {
                      for (j = BlockStart; j <= BlockEnd; j++)
                      {
                          x[xoffst,k,j,REAL] = scr[k,j,REAL];
                          x[xoffst,k,j,IMAG] = scr[k,j,IMAG];
                      }
                  }
              }
              else
              {
                  n1 = n / 2;
                  lk = (int)Math.Pow(2, l);
                  li = (int)Math.Pow(2, m - l - 1);
                  lj = 2 * lk;
                  ku = li;

                  for (i = 0; i <= li - 1; i++)
                  {
                      i11 = i * lk;
                      i12 = i11 + n1;
                      i21 = i * lj;
                      i22 = i21 + lk;

                      u1[REAL] = exponent[(ku+i),REAL];
                      if (is0 >= 1)
                      {
                          u1[IMAG] = exponent[(ku+i),IMAG];
                      }
                      else
                      {
                          u1[IMAG] = -exponent[(ku+i),IMAG];
                      }
                      for (k = 0; k <= lk - 1; k++)
                      {
                          for (j = BlockStart; j <= BlockEnd; j++)
                          {
                              x11[REAL] = scr[(i11+k),j,REAL];
                              x11[IMAG] = scr[(i11+k),j,IMAG];

                              x21[REAL] = scr[(i12+k),j,REAL];
                              x21[IMAG] = scr[(i12+k),j,IMAG];

                              x[xoffst,(i21+k),j,REAL] = x11[REAL] + x21[REAL];
                              x[xoffst,(i21+k),j,IMAG] = x11[IMAG] + x21[IMAG];
                              x[xoffst,(i22+k),j,REAL] = u1[REAL] * (x11[REAL] - x21[REAL]) - u1[IMAG] * (x11[IMAG] - x21[IMAG]);
                              x[xoffst,(i22+k),j,IMAG] = u1[IMAG] * (x11[REAL] - x21[REAL]) + u1[REAL] * (x11[IMAG] - x21[IMAG]);
                          }
                      }
                  }
              }
          }
      }
      if (timeron) timer.stop(4);
  }
}
}
