/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                            R H S C o m p u t e                          !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    RHSCompute implements thread for compute_rhs subroutine of SP        !
!    benchmark.                                                           !
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
!     Translation to Java and to MultiThreaded Code:			  !
!     Michael A. Frumkin					          !
!     Mathew Schultz	   					          !
!-------------------------------------------------------------------------!
*/
package NPB3_0_JAV.SPThreads;
import NPB3_0_JAV.SP;

public class RHSCompute extends SPBase{
  public int id;
  public boolean done = true;
 
  //private arrays and data
  int lower_bound1;
  int upper_bound1;
  int lower_bound2;
  int upper_bound2;
  int state;

  public RHSCompute(SP sp,int low1, int high1, int low2, int high2){
    Init(sp);
    lower_bound1=low1;
    upper_bound1=high1;
    lower_bound2=low2;
    upper_bound2=high2;
    state=1;
    setPriority(Thread.MAX_PRIORITY);
    setDaemon(true);
    master=sp;
  }
  void Init(SP sp){
    //initialize shared data
    IMAX=sp.IMAX;
    JMAX=sp.JMAX; 
    KMAX=sp.KMAX; 
    problem_size=sp.problem_size; 
    nx2=sp.nx2;
    ny2=sp.ny2;
    nz2=sp.nz2;
    grid_points=sp.grid_points;
    niter_default=sp.niter_default;
    dt_default=sp.dt_default;
    u=sp.u;
    rhs=sp.rhs;
    forcing=sp.forcing;
    isize1=sp.isize1;
    jsize1=sp.jsize1;
    ksize1=sp.ksize1;
    us=sp.us;
    vs=sp.vs;
    ws=sp.ws;
    qs=sp.qs;
    rho_i=sp.rho_i;
    speed=sp.speed;
    square=sp.square;
    jsize2=sp.jsize2;
    ksize2=sp.ksize2;
    ue=sp.ue;
    buf=sp.buf;
    jsize3=sp.jsize3;
    lhs=sp.lhs;
    lhsp=sp.lhsp;
    lhsm=sp.lhsm;
    jsize4=sp.jsize4;
    cv=sp.cv;
    rhon=sp.rhon;
    rhos=sp.rhos;
    rhoq=sp.rhoq;
    cuf=sp.cuf;
    q=sp.q;
    ce=sp.ce;
  }
  public void run(){    
    for(;;){
      synchronized(this){ 
      while(done==true){
	try{
	  wait();
	synchronized(master){ master.notify();}
	}catch(InterruptedException ie){}
      }
      step();
      synchronized(master){done=true;master.notify();}
      }
    }
  }
  
  public void step(){
    int i, j, k, m;
    double aux, rho_inv, uijk, up1, um1, vijk, vp1, vm1,wijk, wp1, wm1;
    switch(state){
    case 1:
//---------------------------------------------------------------------
//      compute the reciprocal of density, and the kinetic energy, 
//      and the speed of sound. 
//---------------------------------------------------------------------

      for(k=lower_bound1;k<=upper_bound1;k++){
	for(j=0;j<=grid_points[1]-1;j++){
	    for(i=0;i<=grid_points[0]-1;i++){
	      rho_inv = 1.0/u[0][i][j][k];
	      rho_i[i][j][k] = rho_inv;
	      us[i][j][k] = 
	    	  u[1][i][j][k] * rho_inv;
	      vs[i][j][k] = 
	    	  u[2][i][j][k] * rho_inv;
	      ws[i][j][k] = 
	    	  u[3][i][j][k] * rho_inv;
	      square[i][j][k]     = 
	    	  0.5* (         u[1][i][j][k]*u[1][i][j][k] +           
	    			  		u[2][i][j][k]*u[2][i][j][k] +
                              u[3][i][j][k]* u[3][i][j][k] ) * rho_inv;
	      qs[i][j][k] = 
	    	  square[i][j][k] * rho_inv;
//---------------------------------------------------------------------
//               (don't need speed and ainx until the lhs computation)
//---------------------------------------------------------------------
	      aux = c1c2*rho_inv* (u[4][i][j][k] -  square[i][j][k]);
	      speed[i][j][k] = Math.sqrt(aux);
	    }
	}
      }      
      break;
    case 2:
//---------------------------------------------------------------------
// copy the exact forcing term to the right hand side;  because 
// this forcing term is known, we can store it on the whole grid
// including the boundary                   
//---------------------------------------------------------------------

       for(k=lower_bound1;k<=upper_bound1;k++){
          for(j=0;j<=grid_points[1]-1;j++){
             for(i=0;i<=grid_points[0]-1;i++){
                for(m=0;m<=4;m++){
                   rhs[m][i][j][k] = 
                	   forcing[m][i][j][k];
                }
             }
          }
       }

      break;
    case 3:
//---------------------------------------------------------------------
//      compute xi-direction fluxes 
//---------------------------------------------------------------------
       for(k=lower_bound2;k<=upper_bound2;k++){
          for(j=1;j<=ny2;j++){
             for(i=1;i<=nx2;i++){
                uijk = us[i][j][k];
                up1  = us[i][1+j][k];
                um1  = us[i-1][j][k];

                rhs[0][i][j][k] =
                	rhs[0][i][j][k] + dx1tx1 * 
                          (u[0][(i+1)][j][k] - 
                        		  2.0*u[0][i][j][k] + 
                           u[0][(i-1)][j][k]) -
                          tx2 * (u[1][(i+1)][j][k] - 
                        		  u[1][(i-1)][j][k]);

                rhs[1][i][j][k] = 
                	rhs[1][i][j][k] + dx2tx1 * 
                          (u[1][(i+1)][j][k] - 
                        		  2.0*u[1][i][j][k] + 
                           u[1][(i-1)][j][k]) +
                          xxcon2*con43 * (up1 - 2.0*uijk + um1) -
                          tx2 * (u[1][(i+1)][j][k]*up1 - 
                                 u[1][(i-1)][j][k]*um1 +
                                 (u[4][(i+1)][j][k]- 
                                		 square[(i+1)][j][k]-
                                  u[4][(i-1)][j][k]+ 
                                  square[(i-1)][j][k])*
                                  c2);

                rhs[2][i][j][k] = 
                	rhs[2][i][j][k] + dx3tx1 * 
                          (u[2][(i+1)][j][k] - 
                        		  2.0*u[2][i][j][k] +
                           u[2][(i-1)][j][k]) +
                          xxcon2 * (vs[i][1+j][k] -
                        		  2.0*vs[i][j][k] +
                                    vs[i-1][j][k]) -
                          tx2 * (u[2][(i+1)][j][k]*up1 - 
                                 u[2][(i-1)][j][k]*um1);

                rhs[3][i][j][k] =
                	rhs[3][i][j][k] + dx4tx1 * 
                          (u[3][(i+1)][j][k] - 
                        		  2.0*u[3][i][j][k] +
                           u[3][(i-1)][j][k]) +
                          xxcon2 * (ws[i][1+j][k] -
                        		  2.0*ws[i][j][k] +
                                    ws[i-1][j][k]) -
                          tx2 * (u[3][(i+1)][j][k]*up1 - 
                                 u[3][(i-1)][j][k]*um1);

                rhs[4][i][j][k] = 
                	rhs[4][i][j][k] + dx5tx1 * 
                          (u[4][(i+1)][j][k] - 
                        		  2.0*u[4][i][j][k] +
                           u[4][(i-1)][j][k]) +
                          xxcon3 * (qs[i][1+j][k] - 
                        		  2.0*qs[i][j][k] +
                                    qs[i-1][j][k]) +
                          xxcon4 * (up1*up1 -       2.0*uijk*uijk + 
                                    um1*um1) +
                          xxcon5 * (u[4][(i+1)][j][k]*rho_i[i][1+j][k] - 
                                    2.0*u[4][i][j][k]*rho_i[i][j][k] +
                                    u[4][(i-1)][j][k]*
                                    rho_i[i-1][j][k]) -
                          tx2 * ( (c1*u[4][(i+1)][j][k] - 
                                   c2*square[i][1+j][k])*up1 -
                                  (c1*u[4][(i-1)][j][k] - 
                                   c2*square[i-1][j][k])*um1 );
             }

//---------------------------------------------------------------------
//      add fourth order xi-direction dissipation               
//---------------------------------------------------------------------

             i = 1;
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] =
                	rhs[m][i][j][k]- dssp * 
                          ( 5.0*u[m][i][j][k] 
                                  - 4.0*u[m][(i+1)][j][k] +
                                  u[m][(i+2)][j][k]);
             }

             i = 2;
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] = 
                	rhs[m][i][j][k] - dssp * 
                          (-4.0*u[m][(i-1)][j][k] + 
                        		  6.0*u[m][i][j][k] -
                            4.0*u[m][(i+1)][j][k] + 
                            u[m][(i+2)][j][k]);
             }

             for(i=3;i<=nx2-2;i++){
                for(m=0;m<=4;m++){
                   rhs[m][i][j][k] = 
                	   rhs[m][i][j][k] - dssp * 
                          (  u[m][(i-2)][j][k] 
                               - 4.0*u[m][(i-1)][j][k] + 
                           6.0*u[m][i][j][k] - 
                           4.0*u[m][(i+1)][j][k] + 
                               u[m][(i+2)][j][k] );
                }
             }

             i = nx2-1;
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] = 
                	rhs[m][i][j][k] - dssp *
                          ( u[m][(i-2)][j][k] - 
                        		  4.0*u[m][(i-1)][j][k] + 
                            6.0*u[m][i][j][k] -
                            4.0*u[m][(i+1)][j][k] );
             }

             i = nx2;
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] = 
                	rhs[m][i][j][k] - dssp *
                          ( u[m][(i-2)][j][k] -
                        		  4.*u[m][(i-1)][j][k] +
                            5.*u[m][i][j][k] );
             }
          }
       }

      break;
    case 4:
//---------------------------------------------------------------------
//      compute eta-direction fluxes 
//---------------------------------------------------------------------
       for(k=lower_bound2;k<=upper_bound2;k++){
          for(j=1;j<=ny2;j++){
             for(i=1;i<=nx2;i++){
                vijk = vs[i][j][k];
                vp1  = vs[i][(j+1)][k];
                vm1  = vs[i][(j-1)][k];
                rhs[0][i][j][k] = 
                	rhs[0][i][j][k] + dy1ty1 * 
                         (u[0][i][(j+1)][k] - 
                        		 2.0*u[0][i][j][k] + 
                          u[0][i][(j-1)][k]) -
                         ty2 * (u[2][i][(j+1)][k] - 
                        		 u[2][i][(j-1)][k]);
                rhs[1][i][j][k] = 
                	rhs[1][i][j][k] + dy2ty1 * 
                         (u[1][i][(j+1)][k] - 
                        		 2.0*u[1][i][j][k] + 
                          u[1][i][(j-1)][k]) +
                         yycon2 * (us[i][(j+1)][k] - 
                        		 2.0*us[i][j][k] + 
                                   us[i][(j-1)][k]) -
                         ty2 * (u[1][i][(j+1)][k]*vp1 - 
                                u[1][i][(j-1)][k]*vm1);
                rhs[2][i][j][k] = 
                	rhs[2][i][j][k] + dy3ty1 * 
                         (u[2][i][(j+1)][k] - 
                        		 2.0*u[2][i][j][k] + 
                          u[2][i][(j-1)][k]) +
                         yycon2*con43 * (vp1 - 2.0*vijk + vm1) -
                         ty2 * (u[2][i][(j+1)][k]*vp1 - 
                                u[2][i][(j-1)][k]*vm1 +
                                (u[4][i][(j+1)][k] - 
                                		square[i][(j+1)][k] - 
                                 u[4][i][(j-1)][k] + 
                                 square[i][(j-1)][k])
                                *c2);
                rhs[3][i][j][k] = 
                	rhs[3][i][j][k] + dy4ty1 * 
                         (u[3][i][(j+1)][k] - 
                        		 2.0*u[3][i][j][k] + 
                          u[3][i][(j-1)][k]) +
                         yycon2 * (ws[i][(j+1)][k] - 
                        		 2.0*ws[i][j][k] + 
                                   ws[i][(j-1)][k]) -
                         ty2 * (u[3][i][(j+1)][k]*vp1 - 
                                u[3][i][(j-1)][k]*vm1);
                rhs[4][i][j][k] = 
                	rhs[4][i][j][k] + dy5ty1 * 
                         (u[4][i][(j+1)][k] - 
                        		 2.0*u[4][i][j][k] + 
                          u[4][i][(j-1)][k]) +
                         yycon3 * (qs[i][(j+1)][k] - 
                        		 2.0*qs[i][j][k] + 
                                   qs[i][(j-1)][k]) +
                         yycon4 * (vp1*vp1       - 2.0*vijk*vijk + 
                                   vm1*vm1) +
                         yycon5 * (u[4][i][(j+1)][k]*
                        		 rho_i[i][(j+1)][k] - 
                                   2.0*u[4][i][j][k]*
                                   rho_i[i][j][k] +
                                   u[4][i][(j-1)][k]*
                                   rho_i[i][(j-1)][k]) -
                         ty2 * ((c1*u[4][i][(j+1)][k] - 
                                 c2*square[i][(j+1)][k]) * vp1 -
                                (c1*u[4][i][(j-1)][k] - 
                                 c2*square[i][(j-1)][k]) * vm1);
             }
          }

//---------------------------------------------------------------------
//      add fourth order eta-direction dissipation         
//---------------------------------------------------------------------

          j = 1;
          for(i=1;i<=nx2;i++){
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] = 
                	rhs[m][i][j][k]- dssp * 
                          ( 5.0*u[m][i][j][k] - 
                        		  4.0*u[m][i][(j+1)][k] +
                                  u[m][i][(j+2)][k]);
             }
          }

          j = 2;
          for(i=1;i<=nx2;i++){
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] =
                	rhs[m][i][j][k] - dssp * 
                          (-4.0*u[m][i][(j-1)][k] +
                        		  6.0*u[m][i][j][k] -
                            4.0*u[m][i][(j+1)][k] + 
                            u[m][i][(j+2)][k]);
             }
          }

          for(j=3;j<=ny2-2;j++){
             for(i=1;i<=nx2;i++){
                for(m=0;m<=4;m++){
                   rhs[m][i][j][k] =
                	   rhs[m][i][j][k] - dssp * 
                          (  u[m][i][(j-2)][k] - 
                        		  4.0*u[m][i][(j-1)][k] + 
                           6.0*u[m][i][j][k] -
                           4.0*u[m][i][(j+1)][k] + 
                               u[m][i][(j+2)][k] );
                }
             }
          }
 
          j = ny2-1;
          for(i=1;i<=nx2;i++){
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] = 
                	rhs[m][i][j][k] - dssp *
                          ( u[m][i][(j-2)][k] - 
                        		  4.0*u[m][i][(j-1)][k] + 
                            6.0*u[m][i][j][k] - 
                            4.0*u[m][i][(j+1)][k] );
             }
          }

          j = ny2;
          for(i=1;i<=nx2;i++){
             for(m=0;m<=4;m++){
                rhs[m][i][j][k] = 
                	rhs[m][i][j][k] - dssp *
                          ( u[m][i][(j-2)][k] - 
                        		  4.*u[m][i][(j-1)][k] +
                            5.*u[m][i][j][k] );
             }
          }
       }
      break;
    case 5:
//---------------------------------------------------------------------
//      compute zeta-direction fluxes 
//---------------------------------------------------------------------
       for(j=lower_bound2;j<=upper_bound2;j++){
          for(k=1;k<=ny2;k++){
             for(i=1;i<=nx2;i++){
                wijk = ws[i][j][k];
                wp1  = ws[i][j][(k+1)];
                wm1  = ws[i][j][(k-1)];

                rhs[0][i][j][k] = 
                	rhs[0][i][j][k] + dz1tz1 * 
                         (u[0][i][j][(k+1)] - 
                        		 2.0*u[0][i][j][k] + 
                          u[0][i][j][(k-1)]) -
                         tz2 * (u[3][i][j][(k+1)] - 
                        		 u[3][i][j][(k-1)]);
                rhs[1][i][j][k] = 
                	rhs[1][i][j][k] + dz2tz1 * 
                         (u[1][i][j][(k+1)] - 
                        		 2.0*u[1][i][j][k] + 
                          u[1][i][j][(k-1)]) +
                         zzcon2 * (us[i][j][(k+1)] - 
                        		 2.0*us[i][j][k] + 
                                   us[i][j][(k-1)]) -
                         tz2 * (u[1][i][j][(k+1)]*wp1 - 
                                u[1][i][j][(k-1)]*wm1);
                rhs[2][i][j][k] = 
                	rhs[2][i][j][k] + dz3tz1 * 
                         (u[2][i][j][(k+1)] - 
                        		 2.0*u[2][i][j][k] + 
                          u[2][i][j][(k-1)]) +
                         zzcon2 * (vs[i][j][(k+1)] - 
                        		 2.0*vs[i][j][k] + 
                                   vs[i][j][(k-1)]) -
                         tz2 * (u[2][i][j][(k+1)]*wp1 - 
                                u[2][i][j][(k-1)]*wm1);
                rhs[3][i][j][k] = 
                	rhs[3][i][j][k] + dz4tz1 * 
                         (u[3][i][j][(k+1)] - 
                        		 2.0*u[3][i][j][k] + 
                          u[3][i][j][(k-1)]) +
                         zzcon2*con43 * (wp1 - 2.0*wijk + wm1) -
                         tz2 * (u[3][i][j][(k+1)]*wp1 - 
                                u[3][i][j][(k-1)]*wm1 +
                                (u[4][i][j][(k+1)] - 
                                		square[i][j][(k+1)] - 
                                 u[4][i][j][(k-1)] + 
                                 square[i][j][(k-1)])
                                *c2);
                rhs[4][i][j][k] =
                	rhs[4][i][j][k] + dz5tz1 * 
                         (u[4][i][j][(k+1)] - 
                        		 2.0*u[4][i][j][k] + 
                          u[4][i][j][(k-1)]) +
                         zzcon3 * (qs[i][j][(k+1)] - 
                        		 2.0*qs[i][j][k] + 
                                   qs[i][j][(k-1)]) +
                         zzcon4 * (wp1*wp1 - 2.0*wijk*wijk + 
                                   wm1*wm1) +
                         zzcon5 * (u[4][i][j][(k+1)]*
                        		 rho_i[i][j][(k+1)] - 
                                   2.0*u[4][i][j][k]*
                                   rho_i[i][j][k] +
                                   u[4][i][j][(k-1)]*
                                   rho_i[i][j][(k-1)]) -
                         tz2 * ( (c1*u[4][i][j][(k+1)] - 
                                  c2*square[i][j][(k+1)])*wp1 -
                                 (c1*u[4][i][j][(k-1)] - 
                                  c2*square[i][j][(k-1)])*wm1);
             }
          }

//---------------------------------------------------------------------
//      add fourth order zeta-direction dissipation                
//---------------------------------------------------------------------
	   k = 1;
	       for(i=1;i<=nx2;i++){
		   for(m=0;m<=4;m++){
		       rhs[m][i][j][k] = 
		    	   rhs[m][i][j][k]- dssp * 
			   ( 5.0*u[m][i][j][k] -
					   4.0*u[m][i][j][(k+1)] +
                                  u[m][i][j][(k+2)]);
		   }
	       }

	   k = 2;
	       for(i=1;i<=nx2;i++){
		   for(m=0;m<=4;m++){
		       rhs[m][i][j][k] = 
		    	   rhs[m][i][j][k] - dssp * 
			   (-4.0*u[m][i][j][(k-1)] + 
					   6.0*u[m][i][j][k] -
                            4.0*u[m][i][j][(k+1)] +
                            u[m][i][j][(k+2)]);
		   }
	       }

	   for(k=3;k<=nz2-2;k++){
		   for(i=1;i<=nx2;i++){
		       for(m=0;m<=4;m++){
			   rhs[m][i][j][k] = 
				   rhs[m][i][j][k] - dssp * 
			       (  u[m][i][j][(k-2)] -
			    		   4.0*u[m][i][j][(k-1)] + 
				  6.0*u[m][i][j][k] - 
				  4.0*u[m][i][j][(k+1)] + 
				  u[m][i][j][(k+2)] );
		       }
		   }
	   }
 
	   k = nz2-1;
	       for(i=1;i<=nx2;i++){
		   for(m=0;m<=4;m++){
		       rhs[m][i][j][k] = 
		    	   rhs[m][i][j][k] - dssp *
                          ( u[m][i][j][(k-2)] - 
                        		  4.0*u[m][i][j][(k-1)] + 
                            6.0*u[m][i][j][k] - 
                            4.0*u[m][i][j][(k+1)] );
		   }
	       }

	   k = nz2;
	       for(i=1;i<=nx2;i++){
		   for(m=0;m<=4;m++){
		       rhs[m][i][j][k] = 
		    	   rhs[m][i][j][k] - dssp *
                          ( u[m][i][j][(k-2)] -
                        		  4.*u[m][i][j][(k-1)] +
                            5.*u[m][i][j][k] );
		   }
	       }
       }
     break;    
    case 6:
       for(k=lower_bound2;k<=upper_bound2;k++){
          for(j=1;j<=ny2;j++){
             for(i=1;i<=nx2;i++){
                for(m=0;m<=4;m++){
                   rhs[m][i][j][k] = 
                	   rhs[m][i][j][k] * dt;
                }
             }
          }
       }

      break;
   }
    state++;
    if(state==7)state=1;
  }
}
