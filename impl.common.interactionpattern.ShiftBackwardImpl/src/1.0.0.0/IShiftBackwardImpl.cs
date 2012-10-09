using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.interactionpattern.Shift;
using MPI;
using common.direction.Backward;

namespace impl.common.interactionpattern.ShiftBackwardImpl { 

public class IShiftBackwardImpl<DIR> : BaseIShiftBackwardImpl<DIR>, IShift<DIR>
		where DIR:IBackwardDirection
{

	private RequestList requestList = new RequestList();
			
	public IShiftBackwardImpl() 
	{  
	} 
			
	private static int DEFAULT_TAG = 0;		
			
	public void initiate_send()
	{
		if (handle_left != null)
				requestList.Remove(handle_left);
		// int rank = comm.Rank;
		// Console.WriteLine(rank + ": shift-to-left : initiate_send : " + Cell.predecessor);
		handle_left = comm.ImmediateSend<double>(Output_buffer.Array, Cell.predecessor, DEFAULT_TAG);      
		requestList.Add(handle_left);
	}
	
	public void initiate_recv() 
	{
		if (handle_right != null)
				requestList.Remove(handle_right);
		//int rank = comm.Rank;
		//Console.WriteLine(rank + ": shift-to-left : initiate_recv : " + Cell.successor);
		handle_right = comm.ImmediateReceive<double>(Cell.successor, DEFAULT_TAG, Input_buffer.Array);			
		requestList.Add(handle_right);
	}
			
	private Request handle_left;
	private Request handle_right;
			
	public Request HandleLeft { get { return handle_left; } }
	public Request HandleRight { get { return handle_right; } }
			
	public override int go() 
	{ 	   
	   int rank = comm.Rank;
	   //Console.WriteLine(rank + ": shift to left begin - from " + Cell.successor + ", to " + Cell.predecessor + ", count=" + requestList.Count);
	   requestList.WaitAll();
	   //clearRequests();
	   //Console.WriteLine(rank + ": shift to left end");
			
	   return 0;		
	} // end activate method 

	private void clearRequests() 
	{
	   if (handle_left != null) requestList.Remove(handle_left);
	   if (handle_right != null) requestList.Remove(handle_right);
	   handle_left = null;
	   handle_right = null;
	}
}

}
