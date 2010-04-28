using System;

namespace NPB {

public class Random{

  protected static int REAL = 0, IMAG = 1;

  //default seed
  public double tran = 314159265.0;   //First 9 digits of PI
  //Random Number Multiplier
  public double amult = 1220703125.0; //=Math.pow(5.0,13);
  public int KS=0;
  public double R23, R46, T23, T46;
  //constants
  public static double d2m46=Math.Pow(0.5,46);
  protected static long i246m1 = (long)Math.Pow(2,46)-1;
  
  public Random(){}
  public Random(double sd){seed=sd;}
  //Random number generator with an external seed
  public double randlc(double x, double a) { // versao original serial
      //double[] y;
      double r23, r46, t23, t46, t1, t2, t3, t4, a1, a2, x1, x2, z;
      r23 = Math.Pow(0.5, 23);
      r46 = Math.Pow(r23, 2);
      t23 = Math.Pow(2.0, 23);
      t46 = Math.Pow(t23, 2);
      //---------------------------------------------------------------------
      //   Break A into two parts such that A = 2^23 * A1 + A2.
      //---------------------------------------------------------------------
      t1 = r23 * a;
      a1 = (int)t1;
      a2 = a - t23 * a1;
      //---------------------------------------------------------------------
      //   Break X into two parts such that X = 2^23 * X1 + X2, compute
      //   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
      //   X = 2^23 * Z + A2 * X2  (mod 2^46).
      //---------------------------------------------------------------------
      t1 = r23 * x;
      x1 = (int)t1;
      x2 = x - t23 * x1;
      t1 = a1 * x2 + a2 * x1;
      t2 = (int)(r23 * t1);
      z = t1 - t23 * t2;
      t3 = t23 * z + a2 * x2;
      t4 = (int)(r46 * t3);
      x = t3 - t46 * t4;
      return x;
  }
  //Random number generator with an internal seed
  public double randlc(double a){
    //double[] y; 
    double r23,r46,t23,t46,t1,t2,t3,t4,a1,a2,x1,x2,z;
    r23 = Math.Pow(0.5,23); 
    r46 = Math.Pow(r23, 2); 
    t23 = Math.Pow(2.0,23);
    t46 = Math.Pow(t23, 2);
//---------------------------------------------------------------------
//   Break A into two parts such that A = 2^23 * A1 + A2.
//---------------------------------------------------------------------
    t1 = r23 * a;
    a1 = (int) t1;
    a2 = a - t23 * a1;
//---------------------------------------------------------------------
//   Break X into two parts such that X = 2^23 * X1 + X2, compute
//   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
//   X = 2^23 * Z + A2 * X2  (mod 2^46).
//---------------------------------------------------------------------
    t1 = r23 * tran;
    x1 = (int) t1;
    x2 = tran - t23 * x1;
    t1 = a1 * x2 + a2 * x1;
    t2 = (int) (r23 * t1);
    z = t1 - t23 * t2;
    t3 = t23 * z + a2 * x2;
    t4 = (int) (r46 * t3);
    tran = t3 - t46 * t4;
    return(r46 * tran);
  }
  public double vranlc(double n, double x, double a, double[] y,int offset){ 
    long Lx = (long)x;
    long La = (long)a;

    for(int i=0;i<n;i++){
      Lx   = (Lx*La) & (i246m1);
      y[offset+i] = (double)(d2m46* Lx);
    }
    return (double) Lx;
  }
  public double vranlc2(double n, double x, double a, double[,] y,int offset){ 
    long Lx = (long)x;
    long La = (long)a;

    for(int i=0;i<n;i++){
      Lx   = (Lx*La) & (i246m1);
      y[i+offset,REAL] = (double)(d2m46* Lx);
      Lx   = (Lx*La) & (i246m1);
      y[i+offset,IMAG] = (double)(d2m46* Lx);				
    }
    return (double) Lx;
  }

  public double seed;
  public double ipow46(double a, int exponent ){
      int n, n2;
      double q, r;
//---------------------------------------------------------------------
// Use
//   a^n = a^(n/2)*a^(n/2) if n even else
//   a^n = a*a^(n-1)       if n odd
//---------------------------------------------------------------------
      if (exponent == 0) return seed;
      q = a;
      r = 1;
      n = exponent;

      while(n>1){
         n2 = n/2;
         if (n2*2==n){
            seed = randlc(q,q);
	    q=seed;
            n = n2;
         }else{
            seed = randlc(r,q);
	    r=seed;
            n = n-1;
	 }
      }
      seed = randlc(r,q);
      return seed;
  }
  public double power( double a, int n ){
//c---------------------------------------------------------------------
//c     power  raises an integer, disguised as a double
//c     precision real, to an integer power
//c---------------------------------------------------------------------
      double aj,ajj,pow;
      int nj;

      pow = 1.0;
      nj = n;
      aj = a;
      while( nj != 0 ){
        if( nj%2==1 ) {
	  seed=randlc(pow, aj );
	  pow=seed;
	}
	ajj=aj;
        seed=randlc(aj, ajj );
	aj=seed;
        nj = nj/2;
      }
      return pow;
   }
  public double ipow46P(double a, int exp_1, int exp_2) {

      //c---------------------------------------------------------------------
      //c---------------------------------------------------------------------

      //c---------------------------------------------------------------------
      //c compute a^exponent mod 2^46
      //c---------------------------------------------------------------------
      //Random rd = new Random();
      double dummy, q, r;
      int n, n2; // ierr;
      //      external randlc
      //      double precision randlc
      bool two_pow;
      //---------------------------------------------------------------------
      //c Use
      //c   a^n = a^(n/2)*a^(n/2) if n even else
      //c   a^n = a*a^(n-1)       if n odd
      //c---------------------------------------------------------------------
      r = 1;
      if ((exp_2 == 0) || (exp_1 == 0)) return r;
      q = a;
      n = exp_1;
      two_pow = true;

      while (two_pow) {
          n2 = n / 2;
          if (n2 * 2 == n) {
              dummy = this.randlcP(ref q, ref q);
              n = n2;
          }
          else {
              n = n * exp_2;
              two_pow = false;
          }
      }

      while (n > 1) {
          n2 = n / 2;
          if (n2 * 2 == n) {
              dummy = this.randlcP(ref q, ref q);
              n = n2;
          }
          else {
              dummy = this.randlcP(ref r, ref q);
              n = n - 1;
          }
      }
      dummy = this.randlcP(ref r, ref q);
      return r;
  }
  public double randlcP(ref double x, ref double a) {
      //c---------------------------------------------------------------------
      //c
      //c   This routine returns a uniform pseudorandom double precision number in the
      //c   range (0, 1) by using the linear congruential generator
      //c
      //c   x_{k+1} = a x_k  (mod 2^46)
      //c
      //c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
      //c   before repeating.  The argument A is the same as 'a' in the above formula,
      //c   and X is the same as x_0.  A and X must be odd double precision integers
      //c   in the range (1, 2^46).  The returned value RANDLC is normalized to be
      //c   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
      //c   the new seed x_1, so that subsequent calls to RANDLC using the same
      //c   arguments will generate a continuous sequence.
      //c
      //c   This routine should produce the same results on any computer with at least
      //c   48 mantissa bits in double precision floating point data.  On 64 bit
      //c   systems, double precision should be disabled.
      //c
      //c   David H. Bailey     October 26, 1990
      //c
      //c---------------------------------------------------------------------
      double r23, r46, t23, t46, t1, t2, t3, t4, a1, a2, x1, x2, z;
      r23 = Math.Pow(0.5, 23);
      r46 = Math.Pow(r23, 2);
      t23 = Math.Pow(2.0, 23);
      t46 = Math.Pow(t23, 2);

      //c---------------------------------------------------------------------
      //c   Break A into two parts such that A = 2^23 * A1 + A2.
      //c---------------------------------------------------------------------
      t1 = r23 * a;
      a1 = (int)(t1);
      a2 = a - t23 * a1;

      //c---------------------------------------------------------------------
      //c   Break X into two parts such that X = 2^23 * X1 + X2, compute
      //c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
      //c   X = 2^23 * Z + A2 * X2  (mod 2^46).
      //c---------------------------------------------------------------------
      t1 = r23 * x;
      x1 = (int)(t1);
      x2 = x - t23 * x1;

      t1 = a1 * x2 + a2 * x1;
      t2 = (int)(r23 * t1);
      z = t1 - t23 * t2;
      t3 = t23 * z + a2 * x2;
      t4 = (int)(r46 * t3);
      x = t3 - t46 * t4;
      return (r46 * x);
  }
  public void vranlcP(int k, int nn, int nnn, int nnnn, double x, double a, double[,,,] y) {
      //---------------------------------------------------------------------
      //c   This routine generates N uniform pseudorandom double precision numbers in
      //c   the range (0, 1) by using the linear congruential generator
      //c   
      //c   x_{k+1} = a x_k  (mod 2^46)
      //c   
      //c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
      //c   before repeating.  The argument A is the same as 'a' in the above formula,
      //c   and X is the same as x_0.  A and X must be odd double precision integers
      //c   in the range (1, 2^46).  The N results are placed in Y and are normalized
      //c   to be between 0 and 1.  X is updated to contain the new seed, so that
      //c   subsequent calls to RANDLC using the same arguments will generate a
      //c   continuous sequence.
      //c   
      //c   This routine generates the output sequence in batches of length NV, for
      //c   convenience on vector computers.  This routine should produce the same
      //c   results on any computer with at least 48 mantissa bits in double precision
      //c   floating point data.  On Cray systems, double precision should be disabled.
      //c   
      //c   David H. Bailey    August 30, 1990
      //c---------------------------------------------------------------------
      // u1 = y = new double[dims[1, 3] + 1, dims[1, 2] + 1, dims[1, 1] + 1];
      double r23, r46, t23, t46;
      int nv, n;
      int di1 = y.GetLength(0), di2 = y.GetLength(1), di3 = y.GetLength(2), di4 = y.GetLength(3);
      n = nn * nnn * nnnn;
      r23 = Math.Pow(2, (-23));
      r46 = r23 * r23;
      t23 = Math.Pow(2, 23);
      t46 = t23 * t23;
      nv = 64;
      double[] xv = new double[nv + 1];
      double t1, t2, t3, t4, an, a1 = 0, a2 = 0, x1, x2, yy;
      int n1, i, j;
      //      external randlc
      //      double precision randlc
      //c---------------------------------------------------------------------
      //c     Compute the first NV elements of the sequence using RANDLC.
      //c---------------------------------------------------------------------
      t1 = x;
      n1 = min(n, nv);
      for (i = 1; i <= n1; i++) {
          xv[i] = t46 * this.randlcP(ref t1, ref a);
      }
      //c---------------------------------------------------------------------
      //c     It is not necessary to compute AN, A1 or A2 unless N is greater than NV.
      //c---------------------------------------------------------------------
      if (n > nv) {
          //c---------------------------------------------------------------------
          //c     Compute AN = AA ^ NV (mod 2^46) using successive calls to RANDLC.
          //c---------------------------------------------------------------------
          t1 = a;
          t2 = r46 * a;

          for (i = 1; i <= nv - 1; i++) {
              t2 = randlcP(ref t1, ref a);
          }

          an = t46 * t2;

          //c---------------------------------------------------------------------
          //c     Break AN into two parts such that AN = 2^23 * A1 + A2.
          //c---------------------------------------------------------------------
          t1 = r23 * an;
          a1 = (int)(t1);
          a2 = an - t23 * a1;
      }

      //c---------------------------------------------------------------------
      //c     Compute N pseudorandom results in batches of size NV.
      //c---------------------------------------------------------------------
      int idx = 0; int xi = 0; int count = -1;
      int div = 0;
      for (j = 0; j < nn; j++) {
          n1 = min(nv, n - j*nnn);
          refazer:
          for (i=0+div; i < ((int) (nnn/2)+div); i++) {
              for (xi = 0; xi < nnnn; xi++) {
                  count++;
                  int kuu = xi % 2 == 0 ? 0 : 1;
                  idx = ((i)*nnnn)+xi+1;
                  idx = (i >=(int)(nnn/2))?idx-(nnn):idx;
                  y[k-1, j, i, xi] = r46 * xv[idx]; // y[i+j]
                  if (kuu == 1) count--;
              }
          }
          count = -1;
          //c---------------------------------------------------------------------
          //c     If this is the last pass through the 140 loop, it is not necessary to
          //c     update the XV vector.
          //c---------------------------------------------------------------------
          if ((j + n1) == n) goto desvio;

          //c---------------------------------------------------------------------
          //c     Update the XV vector by multiplying each element by AN (mod 2^46).
          //c---------------------------------------------------------------------
          for (i = 1; i <= nv; i++) {
              t1 = r23 * xv[i];      
              x1 = (int)(t1);         
              x2 = xv[i] - t23 * x1;  
              t1 = a1 * x2 + a2 * x1; 
              t2 = (int)(r23 * t1);   
              yy = t1 - t23 * t2;     
              t3 = t23 * yy + a2 * x2;
              t4 = (int)(r46 * t3);   
              xv[i] = t3 - t46 * t4;                                                
          }
          if (div == 0) {
              div = 32;
              goto refazer;
          } else div = 0;
      }

  //c---------------------------------------------------------------------
  //c     Save the last seed in X so that subsequent calls to VRANLC will generate
  //c     a continuous sequence.
  //c---------------------------------------------------------------------
  desvio:
      x = xv[n1];
  }
  //  //**
  public int min(int a, int b) {
      if (a < b) return a; else return b;
  }

  public void imprime(double[,,,] u1) {

      System.IO.TextWriter arquivo = System.IO.File.AppendText("D:/u1.txt");
      int p1, p2, p3, p4; bool flag = true;
      for (p1 = 0; p1 < u1.GetLength(0); p1++)
          for (p2 = 0; p2 < u1.GetLength(1); p2++)
              for (p3 = 0; p3 < u1.GetLength(2); p3++) {
                  for (p4 = 0; p4 < u1.GetLength(3); p4++) {
                      if (u1[p1, p2, p3, p4] != 0) {
                          if (flag) {
                              arquivo.Write("[" + p1 + "," + p2 + "," + p3 + "," + p4 + "]=" + u1[p1, p2, p3, p4] + " :: ");
                              flag = !flag;
                          }
                          else {
                              arquivo.WriteLine("[" + p1 + "," + p2 + "," + p3 + "," + p4 + "]=" + u1[p1, p2, p3, p4]);
                              flag = !flag;
                          }
                          //arquivo.WriteLine(" <---> [" + p1 + "," + p2 + "," + p3 + "]=" + u1[p1, p2+1, 1]);
                      }
                  }
              }

      //arquivo.WriteLine("Posicao: " + p1 + " " + p2 + " " + p3 + " " + p4 + " " + " Valor: " + x[p1, p2, p3, p4]);
      arquivo.Close();

  }

}

}
