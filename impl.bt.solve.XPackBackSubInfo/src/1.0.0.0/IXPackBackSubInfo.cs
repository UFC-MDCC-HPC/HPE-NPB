using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.XAxis;
using bt.solve.BTMethod;
using bt.solve.PackBackSubInfo;

namespace impl.bt.solve.XPackBackSubInfo 
{ 
	public class IXPackBackSubInfo<I, C, DIR, MTH> : BaseIXPackBackSubInfo<I, C, DIR, MTH>, IPackBackSubInfo<I, C, DIR, MTH>
		where I:IInstance_BT<C>
		where C:IClass
		where DIR:IX
		where MTH:IBTMethod 
	{	    
		private int c;
		
		
		private bool buffers_created = false;
		
		private void create_buffers()
		{
		}
		
		public override int go()   
		{
			if (!buffers_created)
			{
				this.create_buffers();
				buffers_created = true;
			}
			
			int stage = this.Iteration_control.getCurrentStage();
				
			c = slice[stage,0];		
			
			int ptr = 0, istart = 2;
			for(int k = 0; k <= KMAX - 1; k++) 
			{
			    for(int j = 0; j <= JMAX - 1; j++) 
			    {
			        for(int n = 0; n < 5; n++) 
			        {
			            in_buffer_x[ptr + n] = rhs[c, k+2, j+2, istart, n];
			        }
			        ptr = ptr + 5;
			    }
			}
		
		    return 0;	
		}		
	}
}
