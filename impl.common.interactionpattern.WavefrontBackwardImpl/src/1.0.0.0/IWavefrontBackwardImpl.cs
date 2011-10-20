using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.interactionpattern.wavefront.WavefrontWork;
using common.direction.Backward;
using common.interactionpattern.Wavefront;

namespace impl.common.interactionpattern.WavefrontBackwardImpl 
{ 
	public class IWavefrontBackwardImpl<W, DIR> : BaseIWavefrontBackwardImpl<W, DIR>, IWavefront<W, DIR>
		where W:IWavefrontWork<DIR>
		where DIR:IBackwardDirection
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
		   if(south != -1) 
		   {
			  double[] buf = Y_buffer.Array;
		      comm.Receive<double>(south, from_s, ref buf);
		   }
		   
		   if(east != -1) 
		   {  
			  double[] buf = X_buffer.Array;
		      comm.Receive<double>(east, from_e, ref buf);
		   }
		}
		
		private void send() 
		{
		   if(north != -1) 
		   { 
		     double[] buf = Y_buffer.Array;
		     comm.Send<double>(buf, north, from_s);
		   }
		     
		   if(west != -1) 
		   { 
		     double[] buf = X_buffer.Array;  
		     comm.Send<double>(buf, west, from_e);
		   }		     
		}
	
	}

}
