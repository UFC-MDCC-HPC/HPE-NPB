using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using lu.problem_size.Instance_LU;
using common.problem_size.Class;
using lu.exchange.ExchangePattern5;
using lu.exchange.Exchange2D;
using MPI;


namespace impl.lu.exchange.Exchange5 
{ 
	public class IExchange5<I, C, E> : BaseIExchange5<I, C, E>, IExchange2D<I, C, E>
		where I:IInstance_LU<C>
		where C:IClass 
		where E:IExchangePattern5
	{
	   
	    protected static int from_s = 1, from_n = 2, from_e = 3, from_w = 4;

		private double[,] g;
	    private int beg;
	    private int fin1; 
	    
 		public override int go() 
		{ 
            int k;            
			
            if(fin1==nx) 
			{
                double[] dum = new double[nz];
                MPI.Request msgid1 = worldcomm.ImmediateReceive<double>(south, from_s, dum); 
				msgid1.Wait();
				
                for(k = 1; k<=nz; k++) 
				{
                    g[k, nx+1] = dum[k-1];
                }
            }
			
            if(beg==1) 
			{
                double[] dum = new double[nz];
                for(k = 1; k<=nz; k++) 
				{
                    dum[k-1] = g[k, 1];
                }
                worldcomm.Send<double>(dum, north, from_s);				
            }
            
			return 0;
		}
		
		
		public void setParameters(double[,] g, int beg, int fin1)
		{
		   this.g = g;
		   this.beg = beg;
		   this.fin1 = fin1;
		}
		
		public void setParameters(double[,,,] g) {}
		public void setParameters(double[,,,] g, int k) {}
		public void setParameters(double[,] g, double[,] h, int ibeg, int ifin1, int jbeg, int jfin1) {}

	}
}
