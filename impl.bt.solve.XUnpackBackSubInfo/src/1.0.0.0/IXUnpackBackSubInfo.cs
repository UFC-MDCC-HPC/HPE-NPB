using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.XAxis;
using bt.solve.BTMethod;
using bt.solve.UnpackBackSubInfo;

namespace impl.bt.solve.XUnpackBackSubInfo 
{ 
	public class IXUnpackBackSubInfo<I, C, DIR, MTH> : BaseIXUnpackBackSubInfo<I, C, DIR, MTH>, IUnpackBackSubInfo<I, C, DIR, MTH>
		where I:IInstance_BT<C>
		where C:IClass
		where DIR:IX
		where MTH:IBTMethod 
	{	   
		private int c;		

		public override int go()   
		{
			backsub_info = Backsub_info.Field;
						
			int stage = this.Iteration_control.getCurrentStage();
				
			c = slice[stage,0];
			
			int ptr = 0;
			for(int k = 0; k <= KMAX - 1; k++) 
			{
			    for(int j = 0; j <= JMAX - 1; j++) 
			    {
			        for(int n = 0; n < 5; n++) 
			        {
			            backsub_info[c, k+2, j+2, n, 0] = out_buffer_x[ptr + n];
			        }
			        ptr = ptr + 5;
			    }
			}
			
			return 0;
		}
	}
}
