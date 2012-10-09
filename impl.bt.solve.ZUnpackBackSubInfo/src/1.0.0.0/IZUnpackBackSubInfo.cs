using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.ZAxis;
using bt.solve.BTMethod;
using bt.solve.UnpackBackSubInfo;

namespace impl.bt.solve.ZUnpackBackSubInfo { 
	public class IZUnpackBackSubInfo<I, C, DIR, MTH> : BaseIZUnpackBackSubInfo<I, C, DIR, MTH>, IUnpackBackSubInfo<I, C, DIR, MTH>
	where I:IInstance_BT<C>
	where C:IClass
	where DIR:IZ
	where MTH:IBTMethod {
	   
		private int c;
		
		public override int go()   
		{
			backsub_info = Backsub_info.Field;
			
			int stage = this.Iteration_control.getCurrentStage();
				
			c = slice[stage,2];
			
            int ptr = 0;
            for(int j = 0; j <= JMAX - 1; j++) 
            {
                for(int i = 0; i <= IMAX - 1; i++) 
                {
                    for(int n = 0; n < 5; n++) 
                    {
                        backsub_info[c, j+2, i+2, n, 0] = out_buffer_z[ptr + n];
                    }
                    ptr = ptr + 5;
                }
            }
            return 0;
		}
	}
}
