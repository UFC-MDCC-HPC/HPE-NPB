using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.interactionpattern.wavefront.WavefrontWork;
using common.direction.Forward;
using common.interactionpattern.Wavefront;

namespace impl.common.interactionpattern.WavefrontForwardImpl 
{ 
	public class IWavefrontForwardImpl<W, DIR> : BaseIWavefrontForwardImpl<W, DIR>, IWavefront<W, DIR>
		where W:IWavefrontWork<DIR>
		where DIR:IForwardDirection
	{	
	    protected static int from_s = 1, from_n = 2, from_e = 3, from_w = 4;
		
		public override int go() 
		{ 
		    receive();
			Wavefront_work.go();
			send();
			
			return 0;
		}  
		
		private void receive() 
		{
		   if(north != -1) 
		   {
			  double[] buf = Y_buffer.Array;
			  comm.Receive<double>(north, from_n, ref buf);
		   }
		   
		   if(west != -1) 
		   {  
			  double[] buf = X_buffer.Array;
		      comm.Receive<double>(west, from_w, ref buf);
		   }
		}
		
		private void send() 
		{
		   if(south != -1) 
		   { 
		     double[] buf = Y_buffer.Array;
		     comm.Send<double>(buf, south, from_n);
		   }
		     
		   if(east != -1) 
		   { 
		     double[] buf = X_buffer.Array; 
		     comm.Send<double>(buf, east, from_w);
		   }		     
		}
	
	}

}
