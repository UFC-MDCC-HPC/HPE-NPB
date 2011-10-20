using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.axis.XAxis;
using lu.exchange.Exchange3;

namespace impl.lu.exchange.Exchange3XImpl { 

	public class IExchange3XImpl<O> : BaseIExchange3XImpl<O>, IExchange3<O>
   	   where O:IX
	{	
		protected static int from_s = 1, from_n = 2, from_e = 3, from_w = 4;
		
		public override int go() 
		{ 
			Shift_forward.setTag(from_w);
			Shift_forward.go();
			
			Shift_backward.setTag(from_e);
			Shift_backward.go();
		
			return 0;
		} 
	
		private double[,,,] g;
		
 	    public void setParameters(double[,,,] g, int nx, int ny, int nz) 
		{
			this.g = g;
		    Buffer_reader_backward.setParameters(g, nx, ny, nz);
			Buffer_reader_forward.setParameters(g, nx, ny, nz);
			Buffer_writer_backward.setParameters(g, nx, ny, nz);
			Buffer_writer_forward.setParameters(g, nx, ny, nz);
		}
	}

}
