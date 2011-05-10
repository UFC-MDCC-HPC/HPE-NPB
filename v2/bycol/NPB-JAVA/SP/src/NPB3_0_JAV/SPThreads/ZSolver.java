/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                             Z S o l v e r                               !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    ZSolver implements thread for z_solve subroutine                     !
!    of the SP benchmark.                                                 !
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
package NPB3_0_JAV.SPThreads;

import NPB3_0_JAV.SP;

public class ZSolver extends SPBase {
	public int id;
	public boolean done = true;

	// private arrays and data
	int lower_bound;
	int upper_bound;
	int state = 1;
	double lhs[][], lhsm[][], lhsp[][], cv[], rhos[];

	public ZSolver(SP sp, int low, int high) {
		Init(sp);
		lower_bound = low;
		upper_bound = high;
		setPriority(Thread.MAX_PRIORITY);
		setDaemon(true);
		master = sp;
		lhs = new double[5][(problem_size + 1)];
		lhsp = new double[5][(problem_size + 1)];
		lhsm = new double[5][(problem_size + 1)];
		cv = new double[problem_size];
		rhos = new double[problem_size];
	}

	void Init(SP sp) {
		// initialize shared data
		IMAX = sp.IMAX;
		JMAX = sp.JMAX;
		KMAX = sp.KMAX;
		problem_size = sp.problem_size;
		nx2 = sp.nx2;
		ny2 = sp.ny2;
		nz2 = sp.nz2;
		grid_points = sp.grid_points;
		niter_default = sp.niter_default;
		dt_default = sp.dt_default;
		u = sp.u;
		rhs = sp.rhs;
		forcing = sp.forcing;
		isize1 = sp.isize1;
		jsize1 = sp.jsize1;
		ksize1 = sp.ksize1;
		us = sp.us;
		vs = sp.vs;
		ws = sp.ws;
		qs = sp.qs;
		rho_i = sp.rho_i;
		speed = sp.speed;
		square = sp.square;
		jsize2 = sp.jsize2;
		ksize2 = sp.ksize2;
		ue = sp.ue;
		buf = sp.buf;
		jsize3 = sp.jsize3;
		lhs = sp.lhs;
		lhsp = sp.lhsp;
		lhsm = sp.lhsm;
		jsize4 = sp.jsize4;
		cv = sp.cv;
		rhon = sp.rhon;
		rhos = sp.rhos;
		rhoq = sp.rhoq;
		cuf = sp.cuf;
		q = sp.q;
		ce = sp.ce;
	}

	public void run() {
		for (;;) {
			synchronized (this) {
				while (done == true) {
					try {
						wait();
						synchronized (master) {
							master.notify();
						}
					} catch (InterruptedException ie) {
					}
				}
				step();
				synchronized (master) {
					done = true;
					master.notify();
				}
			}
		}
	}

	public void step(){
    int i, j, k, n, k1, k2, m;
    double  t1, t2, t3, ac, xvel, yvel, zvel, r1, r2, r3,r4, r5, 
            btuz, acinv, ac2u, uzik1;
    double ru1, fac1, fac2, rtmp[] = new double[5*(KMAX+1)];

    switch(state){
    case 1:
       for(j=lower_bound;j<=upper_bound;j++){
          for(i=1;i<=nx2;i++){

//---------------------------------------------------------------------
// Computes the left hand side for the three z-factors   
//---------------------------------------------------------------------

//---------------------------------------------------------------------
// first fill the lhs for the u-eigenvalue                          
//---------------------------------------------------------------------

             for(k=0;k<=nz2+1;k++){
                ru1 = c3c4*rho_i[i][j][k];
                cv[k] = ws[i][j][k];
                rhos[k] = dmax1(dz4 + con43 * ru1,
                                dz5 + c1c5 * ru1,
                                dzmax + ru1,
                                dz1);
             }

             lhsinit(grid_points[2]-1);
             for(k=1;k<=nz2;k++){
                lhs[0][k] =  0.0;
                lhs[1][k] = -dttz2 * cv[k-1] - dttz1 * rhos[k-1];
                lhs[2][k] =  1.0 + c2dttz1 * rhos[k];
                lhs[3][k] =  dttz2 * cv[k+1] - dttz1 * rhos[k+1];
                lhs[4][k] =  0.0;
             }

//---------------------------------------------------------------------
//      add fourth order dissipation                                  
//---------------------------------------------------------------------

             k = 1;
             lhs[2][k] =  lhs[2][k] + comz5;
             lhs[3][k] =  lhs[3][k] - comz4;
             lhs[4][k] =  lhs[4][k] + comz1;
	       
             k = 2;
             lhs[1][k] =  lhs[1][k] - comz4;
             lhs[2][k] = 
            	 lhs[2][k] + comz6;
             lhs[3][k] =  lhs[3][k] - comz4;
             lhs[4][k] =  lhs[4][k] + comz1;

             for(k=3;k<=nz2-2;k++){
                lhs[0][k] = 
                	lhs[0][k] + comz1;
                lhs[1][k] =   lhs[1][k] - comz4;
                lhs[2][k] =  lhs[2][k] + comz6;
                lhs[3][k] =	lhs[3][k] - comz4;
                lhs[4][k] = lhs[4][k] + comz1;
             }

             k = nz2-1;
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
             for(k=1;k<=nz2;k++){
                lhsp[0][k] = lhs[0][k];
                lhsp[1][k] = lhs[1][k] - 
                                  dttz2 * speed[i][j][(k-1)];
                lhsp[2][k] = lhs[2][k];
                lhsp[3][k] = lhs[3][k] + dttz2 * speed[i][j][(k+1)];
                lhsp[4][k] = lhs[4][k];
                lhsm[0][k] = lhs[0][k];
                lhsm[1][k] = lhs[1][k] +  dttz2 * speed[i][j][(k-1)];
                lhsm[2][k] = lhs[2][k];
                lhsm[3][k] = lhs[3][k] -   dttz2 * speed[i][j][(k+1)];
                lhsm[4][k] = lhs[4][k];
             }

//---------------------------------------------------------------------
// Load a row of K data
//---------------------------------------------------------------------

             for(k=0;k<=nz2+1;k++){
	        rtmp[0][k] = rhs[0][i][j][k];
	        rtmp[1][k] = rhs[1][i][j][k];
	        rtmp[2][k] = rhs[2][i][j][k];
	        rtmp[3][k] = rhs[3][i][j][k];
	        rtmp[4][k] = rhs[4][i][j][k];
	     }

//---------------------------------------------------------------------
//                          FORWARD ELIMINATION  
//---------------------------------------------------------------------

             for(k=0;k<=grid_points[2]-3;k++){
                k1 = k  + 1;
                k2 = k  + 2;
                fac1      = 1./lhs[2][k];
                lhs[3][k]  = fac1*lhs[3][k];
                lhs[4][k]  = fac1*lhs[4][k];
                for(m=0;m<=2;m++){
                   rtmp[m][k] = fac1*rtmp[m][k];
                }
                lhs[2][k1] = lhs[2][k1] - lhs[1][k1]*lhs[3][k];
                lhs[3][k1] = lhs[3][k1] - lhs[1][k1]*lhs[4][k];
                for(m=0;m<=2;m++){
                   rtmp[m][k1] = rtmp[m][k1] - lhs[1][k1]*rtmp[m][k];
                }
                lhs[1][k2] = lhs[1][k2] -lhs[0][k2]*lhs[3][k];
                lhs[2][k2] = lhs[2][k2] -lhs[0][k2]*lhs[4][k];
                for(m=0;m<=2;m++){
                   rtmp[m][k2] = rtmp[m][k2] - lhs[0][k2]*rtmp[m][k];
                }
             }

//---------------------------------------------------------------------
//      The last two rows in this grid block are a bit different, 
//      since they do not have two more rows available for the
//      elimination of off-diagonal entries
//---------------------------------------------------------------------
             k  = grid_points[2]-2;
             k1 = grid_points[2]-1;
             fac1      = 1./lhs[2][k];
             lhs[3][k]  = fac1*lhs[3][k];
             lhs[4][k]  = fac1*lhs[4][k];
             for(m=0;m<=2;m++){
                rtmp[m][k] = fac1*rtmp[m][k];
             }
             lhs[2][k1] = lhs[2][k1] - lhs[1][k1]*lhs[3][k];
             lhs[3][k1] = lhs[3][k1] - lhs[1][k1]*lhs[4][k];
             for(m=0;m<=2;m++){
                rtmp[m][k1] = rtmp[m][k1] - lhs[1][k1]*rtmp[m][k];
             }
//---------------------------------------------------------------------
//               scale the last row immediately
//---------------------------------------------------------------------
             fac2      = 1./lhs[2][k1];
             for(m=0;m<=2;m++){
                rtmp[m][k1] = fac2*rtmp[m][k1];
             }

//---------------------------------------------------------------------
//      do the u+c and the u-c factors               
//---------------------------------------------------------------------
             for(k=0;k<=grid_points[2]-3;k++){
             	k1 = k  + 1;
             	k2 = k  + 2;
	     	m = 3;
             	fac1	   = 1./lhsp[2][k];
             	lhsp[3][k]  = fac1*lhsp[3][k];
             	lhsp[4][k]  = fac1*lhsp[4][k];
             	rtmp[m][k]  = fac1*rtmp[m][k];
             	lhsp[2][k1] = lhsp[2][k1] - lhsp[1][k1]*lhsp[3][k];
             	lhsp[3][k1] = lhsp[3][k1] -lhsp[1][k1]*lhsp[4][k];
             	rtmp[m][k1] = rtmp[m][k1] - lhsp[1][k1]*rtmp[m][k];
             	lhsp[1][k2] = lhsp[1][k2] -lhsp[0][k2]*lhsp[3][k];
             	lhsp[2][k2] = lhsp[2][k2] - lhsp[0][k2]*lhsp[4][k];
             	rtmp[m][k2] = rtmp[m][k2] - lhsp[0][k2]*rtmp[m][k];
	     	m = 4;
             	fac1	   = 1./lhsm[2][k];
             	lhsm[3][k]  = fac1*lhsm[3][k];
             	lhsm[4][k]  = fac1*lhsm[4][k];
             	rtmp[m][k]  = fac1*rtmp[m][k];
             	lhsm[2][k1] = lhsm[2][k1] -lhsm[1][k1]*lhsm[3][k];
             	lhsm[3][k1] = lhsm[3][k1] -lhsm[1][k1]*lhsm[4][k];
             	rtmp[m][k1] = rtmp[m][k1] - lhsm[1][k1]*rtmp[m][k];
             	lhsm[1][k2] = lhsm[1][k2] - lhsm[0][k2]*lhsm[3][k];
             	lhsm[2][k2] = lhsm[2][k2] - lhsm[0][k2]*lhsm[4][k];
             	rtmp[m][k2] = rtmp[m][k2] -lhsm[0][k2]*rtmp[m][k];
             }

//---------------------------------------------------------------------
//         And again the last two rows separately
//---------------------------------------------------------------------
             k  = grid_points[2]-2;
             k1 = grid_points[2]-1;
	     m = 3;
             fac1	= 1./lhsp[2][k];
             lhsp[3][k]  = fac1*lhsp[3][k];
             lhsp[4][k]  = fac1*lhsp[4][k];
             rtmp[m][k]  = fac1*rtmp[m][k];
             lhsp[2][k1] = lhsp[2][k1] - lhsp[1][k1]*lhsp[3][k];
             lhsp[3][k1] = lhsp[3][k1] -lhsp[1][k1]*lhsp[4][k];
             rtmp[m][k1] = rtmp[m][k1] - lhsp[1][k1]*rtmp[m][k];
	     m = 4;
             fac1	= 1.0/lhsm[2][k];
             lhsm[3][k]  = fac1*lhsm[3][k];
             lhsm[4][k]  = fac1*lhsm[4][k];
             rtmp[m][k]  = fac1*rtmp[m][k];
             lhsm[2][k1] = lhsm[2][k1] - lhsm[1][k1]*lhsm[3][k];
             lhsm[3][k1] = lhsm[3][k1] -lhsm[1][k1]*lhsm[4][k];
             rtmp[m][k1] = rtmp[m][k1] - lhsm[1][k1]*rtmp[m][k];
//---------------------------------------------------------------------
//               Scale the last row immediately (some of this is overkill
//               if this is the last cell)
//---------------------------------------------------------------------
             rtmp[3][k1] = rtmp[3][k1]/lhsp[2][k1];
             rtmp[4][k1] = rtmp[4][k1]/lhsm[2][k1];

//---------------------------------------------------------------------
//                         BACKSUBSTITUTION 
//---------------------------------------------------------------------

             k  = grid_points[2]-2;
             k1 = grid_points[2]-1;
             for(m=0;m<=2;m++){
                rtmp[m][k] = rtmp[m][k] -  lhs[3][k]*rtmp[m][k1];
             }

             rtmp[3][k] = rtmp[3][k] -  lhsp[3][k]*rtmp[3][k1];
             rtmp[4][k] = rtmp[4][k] -  lhsm[3][k]*rtmp[4][k1];

//---------------------------------------------------------------------
//      Whether or not this is the last processor, we always have
//      to complete the back-substitution 
//---------------------------------------------------------------------

//---------------------------------------------------------------------
//      The first three factors
//---------------------------------------------------------------------
             for(k=grid_points[2]-3;k>=0;k--){
                k1 = k  + 1;
                k2 = k  + 2;
                for(m=0;m<=2;m++){
                   rtmp[m][k] = rtmp[m][k] -  lhs[3][k]*rtmp[m][k1] - lhs[4][k]*rtmp[m][k2];
                }

//---------------------------------------------------------------------
//      And the remaining two
//---------------------------------------------------------------------
                rtmp[3][k] = rtmp[3][k] - 
                                lhsp[3][k]*rtmp[3][k1] - lhsp[4][k]*rtmp[3][k2];
                rtmp[4][k] = rtmp[4][k] - lhsm[3][k]*rtmp[4][k1] - lhsm[4][k]*rtmp[4][k2];
             }

//---------------------------------------------------------------------
//      Store result
//---------------------------------------------------------------------
             for(k=0;k<=nz2+1;k++){
	        rhs[0][i][j][k] = rtmp[0][k];
	        rhs[1][i][j][k] = rtmp[1][k];
	        rhs[2][i][j][k] = rtmp[2][k];
	        rhs[3][i][j][k] = rtmp[3][k];
	        rhs[4][i][j][k] = rtmp[4][k];
	     }
          }
       }

      break;
    case 2:

       for(k=lower_bound;k<=upper_bound;k++){
          for(j=1;j<=ny2;j++){
             for(i=1;i<=nx2;i++){

                xvel = us[i][j][k];
                yvel = vs[i][j][k];
                zvel = ws[i][j][k];
                ac   = speed[i][j][k];

                ac2u = ac*ac;

                r1 = rhs[0][i][j][k];
                r2 = rhs[1][i][j][k];
                r3 = rhs[2][i][j][k];
                r4 = rhs[3][i][j][k];
                r5 = rhs[4][i][j][k];

                uzik1 = u[0][i][j][k];
                btuz  = bt * uzik1;

                t1 = btuz/ac * (r4 + r5);
                t2 = r3 + t1;
                t3 = btuz * (r4 - r5);

                rhs[0][i][j][k] = t2;
                rhs[1][i][j][k] = -uzik1*r2 + xvel*t2;
                rhs[2][i][j][k] =  uzik1*r1 + yvel*t2;
                rhs[3][i][j][k] =  zvel*t2  + t3;
                rhs[4][i][j][k] =  uzik1*(-xvel*r2 + yvel*r1) +
                          qs[i+j*jsize2+k*ksize2]*t2 + c2iv*ac2u*t1 + zvel*t3;

             }
          }
       }

      break;      
    };
    state++;
    if(state==3)state=1;
  }

	public void lhsinit(int size) {
		int i, n;

		// ---------------------------------------------------------------------
		// zap the whole left hand side for starters
		// ---------------------------------------------------------------------
		for (i = 0; i <= size; i += size) {
			for (n = 0; n <= 4; n++) {
				lhs[n][i] = 0.0;
				lhsp[n][i] = 0.0;
				lhsm[n][i] = 0.0;
			}
		}

		// ---------------------------------------------------------------------
		// next, set all diagonal values to 1. This is overkill, but
		// convenient
		// ---------------------------------------------------------------------
		for (i = 0; i <= size; i += size) {
			lhs[2][i] = 1.0;
			lhsp[2][i] = 1.0;
			lhsm[2][i] = 1.0;
		}
	}
}
