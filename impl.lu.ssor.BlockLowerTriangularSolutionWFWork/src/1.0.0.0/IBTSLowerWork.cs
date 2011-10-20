using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using lu.problem_size.Instance;
using common.problem_size.Class;
using common.direction.Forward;
using lu.ssor.BlockTriangularSolutionWFWork;

namespace impl.lu.ssor.BlockLowerTriangularSolutionWFWork { 

public class IBTSLowerWork<DIR, I, C> : BaseIBTSLowerWork<DIR, I, C>, IBTSWork<DIR, I, C>
	where DIR:IForwardDirection
	where I:IInstance<C>
	where C:IClass
{
	public override int go() 
	{ 
        int i, j, m;
        double  tmp, tmp1;
        double[,] tmat = new double[5, 5];
			
		//Console.WriteLine("BTS LOWER - n=" + north + " s=" + south + " e=" + east + " w=" + west);
			
        readBuffer(rsd ,k);
		
        for(j = jst; j<= jend; j++) 
		{
            for(i = ist; i<= iend; i++) 
			{
                for(m = 1; m<= 5; m++) 
				{
                    rsd[k-1, j+1, i+1, m-1] =  rsd[k-1, j+1, i+1, m -1]
         - omega * (  ldz[j-1, i-1, 0, m-1] * rsd[k-2, j+1, i+1, 0]
                    + ldz[j-1, i-1, 1, m-1] * rsd[k-2, j+1, i+1, 1]
                    + ldz[j-1, i-1, 2, m-1] * rsd[k-2, j+1, i+1, 2]
                    + ldz[j-1, i-1, 3, m-1] * rsd[k-2, j+1, i+1, 3]
                    + ldz[j-1, i-1, 4, m-1] * rsd[k-2, j+1, i+1, 4]);
                }
            }
        }
		
        for(j=jst; j<=jend; j++) 
		{
            for(i = ist; i<= iend; i++) 
			{
                for(m = 1; m<= 5; m++) 
				{
                    rsd[k-1, j+1, i+1, m-1] =  rsd[k-1, j+1, i+1, m-1]
                    - omega * ( ldy[j-1, i-1, 0, m-1] * rsd[k-1, j, i+1, 0]
                              + ldx[j-1, i-1, 0, m-1] * rsd[k-1, j+1, i, 0]
                              + ldy[j-1, i-1, 1, m-1] * rsd[k-1, j, i+1, 1]
                              + ldx[j-1, i-1, 1, m-1] * rsd[k-1, j+1, i, 1]
                              + ldy[j-1, i-1, 2, m-1] * rsd[k-1, j, i+1, 2]
                              + ldx[j-1, i-1, 2, m-1] * rsd[k-1, j+1, i, 2]
                              + ldy[j-1, i-1, 3, m-1] * rsd[k-1, j, i+1, 3]
                              + ldx[j-1, i-1, 3, m-1] * rsd[k-1, j+1, i, 3]
                              + ldy[j-1, i-1, 4, m-1] * rsd[k-1, j, i+1, 4]
                              + ldx[j-1, i-1, 4, m-1] * rsd[k-1, j+1, i, 4]);
                }
				
                //---------------------------------------------------------------------
                //   diagonal block inversion
                //
                //   forward elimination
                //---------------------------------------------------------------------
				for(m = 0; m< 5; m++) 
				{
                    tmat[0, m] = d[j-1, i-1, 0, m];
                    tmat[1, m] = d[j-1, i-1, 1, m];
                    tmat[2, m] = d[j-1, i-1, 2, m];
                    tmat[3, m] = d[j-1, i-1, 3, m];
                    tmat[4, m] = d[j-1, i-1, 4, m];
                }
				
                tmp1 = 1.0d /tmat[0, 0];
                tmp = tmp1 * tmat[0, 1];
                tmat[1, 1] =  tmat[1, 1] - tmp * tmat[1, 0];
                tmat[2, 1] =  tmat[2, 1] - tmp * tmat[2, 0];
                tmat[3, 1] =  tmat[3, 1] - tmp * tmat[3, 0];
                tmat[4, 1] =  tmat[4, 1] - tmp * tmat[4, 0];
                rsd[k-1, j+1, i+1, 1] = rsd[k-1, j+1, i+1, 1] - rsd[k-1, j+1, i+1, 0] * tmp;

                tmp = tmp1 * tmat[0, 2];
                tmat[1, 2] =  tmat[1, 2] - tmp * tmat[1, 0];
                tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 0];
                tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 0];
                tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 0];
                rsd[k-1, j+1, i+1, 2] = rsd[k-1, j+1, i+1, 2] - rsd[k-1, j+1, i+1, 0] * tmp;

                tmp = tmp1 * tmat[0, 3];
                tmat[1, 3] =  tmat[1, 3] - tmp * tmat[1, 0];
                tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 0];
                tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 0];
                tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 0];
                rsd[k-1, j+1, i+1, 3] = rsd[k-1, j+1, i+1, 3] - rsd[k-1, j+1, i+1, 0] * tmp;

                tmp = tmp1 * tmat[0, 4];
                tmat[1, 4] =  tmat[1, 4] - tmp * tmat[1, 0];
                tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 0];
                tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 0];
                tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 0];
                rsd[k-1, j+1, i+1, 4] = rsd[k-1, j+1, i+1, 4] - rsd[k-1, j+1, i+1, 0] * tmp;

                tmp1 = 1.0d /tmat[1, 1];
                tmp = tmp1 * tmat[1, 2];
                tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 1];
                tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 1];
                tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 1];
                rsd[k-1, j+1, i+1, 2] = rsd[k-1, j+1, i+1, 2] - rsd[k-1, j+1, i+1, 1] * tmp;

                tmp = tmp1 * tmat[1, 3];
                tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 1];
                tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 1];
                tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 1];
                rsd[k-1, j+1, i+1, 3] = rsd[k-1, j+1, i+1, 3] - rsd[k-1, j+1, i+1, 1] * tmp;

                tmp = tmp1 * tmat[1, 4];
                tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 1];
                tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 1];
                tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 1];
                rsd[k-1, j+1, i+1, 4] = rsd[k-1, j+1, i+1, 4] - rsd[k-1, j+1, i+1, 1] * tmp;

                tmp1 = 1.0d /tmat[2, 2];
                tmp = tmp1 * tmat[2, 3];
                tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 2];
                tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 2];
                rsd[k-1, j+1, i+1, 3] = rsd[k-1, j+1, i+1, 3] - rsd[k-1, j+1, i+1, 2] * tmp;

                tmp = tmp1 * tmat[2, 4];
                tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 2];
                tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 2];
                rsd[k-1, j+1, i+1, 4] = rsd[k-1, j+1, i+1, 4] - rsd[k-1, j+1, i+1, 2] * tmp;

                tmp1 = 1.0d /tmat[3, 3];
                tmp = tmp1 * tmat[3, 4];
                tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 3];
                rsd[k-1, j+1, i+1, 4] = rsd[k-1, j+1, i+1, 4] - rsd[k-1, j+1, i+1, 3] * tmp;

                //---------------------------------------------------------------------
                //   back substitution
                //---------------------------------------------------------------------

				rsd[k-1, j+1, i+1, 4] = rsd[k-1, j+1, i+1, 4]/ tmat[4, 4];
                rsd[k-1, j+1, i+1, 3] = rsd[k-1, j+1, i+1, 3]- tmat[4, 3] * rsd[k-1, j+1, i+1, 4];
                rsd[k-1, j+1, i+1, 3] = rsd[k-1, j+1, i+1, 3]/ tmat[3, 3];
                rsd[k-1, j+1, i+1, 2] = rsd[k-1, j+1, i+1, 2] -tmat[3, 2] * rsd[k-1, j+1, i+1, 3] - tmat[4, 2] * 
                                                                   rsd[k-1, j+1, i+1, 4];
                rsd[k-1, j+1, i+1, 2] = rsd[k-1, j+1, i+1, 2] /tmat[2, 2];
                rsd[k-1, j+1, i+1, 1] = rsd[k-1, j+1, i+1, 1]- tmat[2, 1] * rsd[k-1, j+1, i+1, 2]-tmat[3, 1]*
                                                                   rsd[k-1, j+1, i+1, 3]-tmat[4, 1]*
                                                                   rsd[k-1, j+1, i+1, 4];
                rsd[k-1, j+1, i+1, 1] = rsd[k-1, j+1, i+1, 1] /tmat[1, 1];
                rsd[k-1, j+1, i+1, 0] = rsd[k-1, j+1, i+1, 0] -tmat[1, 0] * rsd[k-1, j+1, i+1, 1]-
                                                       tmat[2, 0] * rsd[k-1, j+1, i+1, 2]-
                                                       tmat[3, 0] * rsd[k-1, j+1, i+1, 3]-
                                                       tmat[4, 0] * rsd[k-1, j+1, i+1, 4];
                rsd[k-1, j+1, i+1, 0] = rsd[k-1, j+1, i+1, 0] /tmat[0, 0];
				
			
			}
        }
		
		writeBuffer(rsd, k);
					
		return 0;
	
	} // end activate method
	
    private int k;
    private double omega;
	
	public void setParameters(int k, double omega)
	{
	   this.k     = k;
	   this.omega = omega;
	}	

    // From IExchance10.cs	
	private void readBuffer(double[,,,] g, int k) 
	{
	   int i, j, idx;
	
	    double[] dum1;
	    
	    if (north != -1)
	    {
	        dum1 = Y_buffer.Array;       
	        
	        idx = 0;
	        for(j=jst; j<=jend; j++) 
			{
	            g[k-1, j+1, 1, 0] = dum1[0+idx];
	            g[k-1, j+1, 1, 1] = dum1[1+idx];
	            g[k-1, j+1, 1, 2] = dum1[2+idx];
	            g[k-1, j+1, 1, 3] = dum1[3+idx];
	            g[k-1, j+1, 1, 4] = dum1[4+idx];
					
/*					if (this.GlobalRank == 1)
								Console.WriteLine(
									              dum1[0+idx] + " " + 
									              dum1[1+idx] + " " + 
									              dum1[2+idx] + " " + 
									              dum1[3+idx] + " " +
									              dum1[4+idx] 
									); */
	            idx = idx + 5;
					
	        }
	    }
       
        if (west != -1)
        {
	        dum1 = X_buffer.Array;
	       
	        idx = 0;
	        for(i=ist; i<=iend; i++) 
			{
	            g[k-1, 1, i+1, 0] = dum1[0+idx];
	            g[k-1, 1, i+1, 1] = dum1[1+idx];
	            g[k-1, 1, i+1, 2] = dum1[2+idx];
	            g[k-1, 1, i+1, 3] = dum1[3+idx];
	            g[k-1, 1, i+1, 4] = dum1[4+idx];
	            idx = idx + 5;
	        }
        }
       
	}
	 
	// From IExchange12.cs 
	private void writeBuffer(double[,,,] g, int k) 
	{
        int i, j, idx;
        double[] dum;
        
	    if (south != -1)
        {
	        dum = Y_buffer.Array;
	        
//						Console.WriteLine("$$$$$$");
	        idx = 0;
	        for(j=jst; j<=jend; j++) 
	        {
	            dum[0+idx] = g[k-1, j+1, nx+1, 0];
	            dum[1+idx] = g[k-1, j+1, nx+1, 1];
	            dum[2+idx] = g[k-1, j+1, nx+1, 2];
	            dum[3+idx] = g[k-1, j+1, nx+1, 3];
	            dum[4+idx] = g[k-1, j+1, nx+1, 4];
					
/*					if (this.GlobalRank == 0)
								Console.WriteLine(
									              dum[0+idx] + " " + 
									              dum[1+idx] + " " + 
									              dum[2+idx] + " " + 
									              dum[3+idx] + " " +
									              dum[4+idx] 
									); */
	            idx = idx + 5;
	        }
        }
	   
	    if (east != -1)
        {
	        dum = X_buffer.Array;
	        
	        idx = 0;
	        for(i=ist; i<=iend; i++) 
	        {
	            dum[0+idx] = g[k-1, ny+1, i+1, 0];
	            dum[1+idx] = g[k-1, ny+1, i+1, 1];
	            dum[2+idx] = g[k-1, ny+1, i+1, 2];
	            dum[3+idx] = g[k-1, ny+1, i+1, 3];
	            dum[4+idx] = g[k-1, ny+1, i+1, 4];
	            idx = idx + 5;
	        }
        }
	}
		
}

}


