using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.direction.Backward;
using lu.problem_size.Instance;
using common.problem_size.Class;
using lu.ssor.BlockTriangularSolutionWFWork;

namespace impl.lu.ssor.BlockUpperTriangularSolutionWFWork { 

public class IBTSUpperWork<DIR, I, C> : BaseIBTSUpperWork<DIR, I, C>, IBTSWork<DIR, I, C>
	where DIR:IBackwardDirection
	where I:IInstance<C>
	where C:IClass
   {

	    private int k;
	    private double omega;
	    
		public override int go() 
		{
            int i, j, m;
            double tmp, tmp1;
            double[,] tmat = new double[5, 5];
            double[,,] tv = new double[isiz2, isiz1, 5];
			
			//Console.WriteLine("BTS UPPER - n=" + north + " s=" + south + " e=" + east + " w=" + west);
			
            readBuffer(rsd, k);
			
            for(j = jend-1; j>= jst-1; j--) 
			{
                for(i = iend-1; i>= ist-1; i--) 
				{
                    for(m = 0; m< 5; m++) 
					{
                        tv[j, i, m] = omega*( udz[j, i, 0, m] * rsd[k, j+2, i+2, 0]
                                            + udz[j, i, 1, m] * rsd[k, j+2, i+2, 1]
                                            + udz[j, i, 2, m] * rsd[k, j+2, i+2, 2]
                                            + udz[j, i, 3, m] * rsd[k, j+2, i+2, 3]
                                            + udz[j, i, 4, m] * rsd[k, j+2, i+2, 4]);
                    }
                }
            }
			
            for(j = jend-1; j>=jst-1; j--) 
			{
                for(i = iend-1; i>=ist-1; i--) 
				{
                    for(m = 0; m< 5; m++) 
					{
                        tv[j, i, m] = tv[j, i, m]
                                        + omega * ( udy[j, i, 0, m] * rsd[k-1, j+3, i+2, 0]
                                                  + udx[j, i, 0, m] * rsd[k-1, j+2, i+3, 0]
                                                  + udy[j, i, 1, m] * rsd[k-1, j+3, i+2, 1]
                                                  + udx[j, i, 1, m] * rsd[k-1, j+2, i+3, 1]
                                                  + udy[j, i, 2, m] * rsd[k-1, j+3, i+2, 2]
                                                  + udx[j, i, 2, m] * rsd[k-1, j+2, i+3, 2]
                                                  + udy[j, i, 3, m] * rsd[k-1, j+3, i+2, 3]
                                                  + udx[j, i, 3, m] * rsd[k-1, j+2, i+3, 3]
                                                  + udy[j, i, 4, m] * rsd[k-1, j+3, i+2, 4]
                                                  + udx[j, i, 4, m] * rsd[k-1, j+2, i+3, 4]);
                    }
					
                    for(m = 0; m< 5; m++) 
					{
                        tmat[0, m] = d[j, i, 0, m];
                        tmat[1, m] = d[j, i, 1, m];
                        tmat[2, m] = d[j, i, 2, m];
                        tmat[3, m] = d[j, i, 3, m];
                        tmat[4, m] = d[j, i, 4, m];
                    }

                    tmp1 = 1.0d /tmat[0, 0];
                    tmp = tmp1 * tmat[0, 1];
                    tmat[1, 1] =  tmat[1, 1] - tmp * tmat[1, 0];
                    tmat[2, 1] =  tmat[2, 1] - tmp * tmat[2, 0];
                    tmat[3, 1] =  tmat[3, 1] - tmp * tmat[3, 0];
                    tmat[4, 1] =  tmat[4, 1] - tmp * tmat[4, 0];
                    tv[j, i, 1] = tv[j, i, 1]- tv[j, i, 0] * tmp;

                    tmp = tmp1 * tmat[0, 2];
                    tmat[1, 2] =  tmat[1, 2] - tmp * tmat[1, 0];
                    tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 0];
                    tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 0];
                    tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 0];
                    tv[j, i, 2] = tv[j, i, 2] - tv[j, i, 0] * tmp;

                    tmp = tmp1 * tmat[0, 3];
                    tmat[1, 3] =  tmat[1, 3] - tmp * tmat[1, 0];
                    tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 0];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 0];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 0];
                    tv[j, i, 3] = tv[j, i, 3] - tv[j, i, 0] * tmp;

                    tmp = tmp1 * tmat[0, 4];
                    tmat[1, 4] =  tmat[1, 4] - tmp * tmat[1, 0];
                    tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 0];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 0];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 0];
                    tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 0] * tmp;

                    tmp1 = 1.0d /tmat[1, 1];
                    tmp = tmp1 * tmat[1, 2];
                    tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 1];
                    tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 1];
                    tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 1];
                    tv[j, i, 2] = tv[j, i, 2] - tv[j, i, 1] * tmp;

                    tmp = tmp1 * tmat[1, 3];
                    tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 1];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 1];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 1];
                    tv[j, i, 3] = tv[j, i, 3] - tv[j, i, 1] * tmp;

                    tmp = tmp1 * tmat[1, 4];
                    tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 1];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 1];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 1];
                    tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 1] * tmp;

                    tmp1 = 1.0d /tmat[2, 2];
                    tmp = tmp1 * tmat[2, 3];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 2];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 2];
                    tv[j, i, 3] = tv[j, i, 3] - tv[j, i, 2] * tmp;

                    tmp = tmp1 * tmat[2, 4];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 2];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 2];
                    tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 2] * tmp;

                    tmp1 = 1.0d /tmat[3, 3];
                    tmp = tmp1 * tmat[3, 4];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 3];
                    tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 3] * tmp;

                    tv[j, i, 4] = tv[j, i, 4]/tmat[4,4];
                    tv[j, i, 3] = tv[j, i, 3]-tmat[4,3]*tv[j, i, 4];
                    tv[j, i, 3] = tv[j, i, 3]/tmat[3,3];
                    tv[j, i, 2] = tv[j, i, 2]-tmat[3,2]*tv[j, i, 3]-tmat[4, 2]*tv[j, i, 4];
                    tv[j, i, 2] = tv[j, i, 2]/tmat[2,2];
                    tv[j, i, 1] = tv[j, i, 1]-tmat[2,1]*tv[j, i, 2]-tmat[3, 1]*tv[j, i, 3]-tmat[4, 1]*tv[j, i, 4];
                    tv[j, i, 1] = tv[j, i, 1]/tmat[1,1];
                    tv[j, i, 0] = tv[j, i, 0]-tmat[1,0]*tv[j, i, 1]-tmat[2, 0]*tv[j, i, 2]-tmat[3, 0]*tv[j, i, 3] - tmat[4, 0]*tv[j, i, 4];
                    tv[j, i, 0] = tv[j, i, 0]/tmat[0,0];

                    rsd[k-1, j+2, i+2, 0] = rsd[k-1, j+2, i+2, 0] - tv[j, i, 0];
                    rsd[k-1, j+2, i+2, 1] = rsd[k-1, j+2, i+2, 1] - tv[j, i, 1];
                    rsd[k-1, j+2, i+2, 2] = rsd[k-1, j+2, i+2, 2] - tv[j, i, 2];
                    rsd[k-1, j+2, i+2, 3] = rsd[k-1, j+2, i+2, 3] - tv[j, i, 3];
                    rsd[k-1, j+2, i+2, 4] = rsd[k-1, j+2, i+2, 4] - tv[j, i, 4];
                }
            }

            writeBuffer(rsd, k);

			return 0;
		}
		
		
		
		public void setParameters(int k, double omega)
		{
		   this.k     = k;
		   this.omega = omega;
		}
		
    // From IExchance11.cs	
	private void readBuffer(double[,,,] g, int k) 
	{
	   int i, j, idx;
	
	    double[] dum1;
	    
		if (south != -1)	    
	    {
	        dum1 = Y_buffer.Array;       
	        
	        idx = 0;
            for(j=jst; j<=jend; j++) 
			{
                g[k-1, j+1, nx+2, 0] = dum1[0+idx];
                g[k-1, j+1, nx+2, 1] = dum1[1+idx];
                g[k-1, j+1, nx+2, 2] = dum1[2+idx];
                g[k-1, j+1, nx+2, 3] = dum1[3+idx];
                g[k-1, j+1, nx+2, 4] = dum1[4+idx];
                idx = idx + 5;
            }
	    }
       
		if (east != -1)	    
        {
	        dum1 = X_buffer.Array;
	       
	        idx = 0;
            for(i=ist; i<=iend; i++) {
                g[k-1, ny+2, i+1, 0] = dum1[0+idx];
                g[k-1, ny+2, i+1, 1] = dum1[1+idx];
                g[k-1, ny+2, i+1, 2] = dum1[2+idx];
                g[k-1, ny+2, i+1, 3] = dum1[3+idx];
                g[k-1, ny+2, i+1, 4] = dum1[4+idx];
                idx = idx + 5;
            }
        }
       
	}
	 
	// From IExchange13.cs 
	private void writeBuffer(double[,,,] g, int k) 
	{
        int i, j, idx;
        double[] dum;

	    if (north != -1)
        {        
	        dum = Y_buffer.Array;
	        
	        idx = 0;
	        for(j=jst; j<=jend; j++) {
	            dum[0+idx] = g[k-1, j+1, 2, 0];
	            dum[1+idx] = g[k-1, j+1, 2, 1];
	            dum[2+idx] = g[k-1, j+1, 2, 2];
	            dum[3+idx] = g[k-1, j+1, 2, 3];
	            dum[4+idx] = g[k-1, j+1, 2, 4];
	            idx = idx + 5;
	        }
        }

        if (west != -1)
        {	   
	        dum = X_buffer.Array;
	        
	        idx = 0;
	        for(i=ist; i<=iend; i++) {
	            dum[0+idx] = g[k-1, 2, i+1, 0];
	            dum[1+idx] = g[k-1, 2, i+1, 1];
	            dum[2+idx] = g[k-1, 2, i+1, 2];
	            dum[3+idx] = g[k-1, 2, i+1, 3];
	            dum[4+idx] = g[k-1, 2, i+1, 4];
	            idx = idx + 5;
	        }
        }
	}

   }

}
