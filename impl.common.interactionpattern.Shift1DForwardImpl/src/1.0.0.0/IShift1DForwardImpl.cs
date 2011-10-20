using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.interactionpattern.shift1D.BufferReader;
using common.direction.Forward;
using common.interactionpattern.shift1D.BufferWriter;
using common.interactionpattern.Shift1D;
using MPI;

namespace impl.common.interactionpattern.Shift1DForwardImpl { 

	public class IShift1DForwardImpl<BR, DIR, BW> : BaseIShift1DForwardImpl<BR, DIR, BW>, IShift1D<DIR, BR, BW>
		where BR:IBufferReader<DIR>
		where DIR:IForwardDirection
		where BW:IBufferWriter<DIR>
	{	
	    private static int DEFAULT_TAG = 0;		
	
		private int previous, next;	
				
		private Intracommunicator comm;
		
		override public void initialize()
		{
		   comm = (Intracommunicator) this.WorldComm; //Mpi.localComm(this);	   
		}

		override public void post_initialize()
		{
		   next = Axis.successor;
		   previous = Axis.predecessor;
		}
		
		public override int go() 
		{ 
		   Request handle = null;
			
		   if (previous != -1) 
		   {
			   double[] buf = Input_buffer.Array;
		       handle = comm.ImmediateReceive<double>(previous, tag, buf);
		   }
			
		   if (next != -1) 
		   {
		      Buffer_writer.go();
//			  Console.Error.WriteLine("FORWARD BEGIN RANK #" + this.GlobalRank + " SEND TO " + next + " - tag=" + tag);
			  double[] buf = Output_buffer.Array;
		      comm.Send<double>(buf, next, tag);
//			  Console.Error.WriteLine("FORWARD END RANK #" + this.GlobalRank + " SEND TO " + next + " - tag=" + tag);
		   }
		      
		   if (previous != -1) 
			{
//			  Console.Error.WriteLine("FORWARD BEGIN RANK #" + this.GlobalRank + " RECEIVE FROM " + previous + " - tag=" + tag);
			  handle.Wait();
//			  Console.Error.WriteLine("FORWARD END RANK #" + this.GlobalRank + " RECEIVE FROM " + previous + " - tag=" + tag);
		      Buffer_reader.go();			
			}
		
		   return 0;
		} // end activate method 
	
	    private int tag = DEFAULT_TAG;	

		public void setTag(int tag)
		{
		   this.tag = tag;
		}
		
	}

}
