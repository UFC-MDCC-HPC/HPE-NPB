using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.YAxis;
using bt.solve.BTMethod;
using bt.solve.PackSolveInfo;

namespace impl.bt.solve.YPackSolveInfo 
{
	public class IYPackSolveInfo<I, C, DIR, MTH> : BaseIYPackSolveInfo<I, C, DIR, MTH>, IPackSolveInfo<I, C, DIR, MTH>
		where I:IInstance_BT<C>
		where C:IClass
		where DIR:IY
		where MTH:IBTMethod 
	{
		private int c;
		
		public override int go()   
		{
			lhsc = Lhsc.Field6;
			
			int stage = this.Iteration_control.getCurrentStage();
				
			c = slice[stage,1];
			
            int jsize, ptr=0;
            jsize = cell_size[c, 1] + 1;
            
            for(int k = 0; k <= KMAX - 1; k++) 
            {
                for(int i = 0; i <= IMAX - 1; i++) 
                {
                    for(int m = 1; m <= 5; m++) 
                    {
                        for(int n = 0; n < 5; n++) 
                        {
                            in_buffer_y[ptr + n] = lhsc[c, k+2, jsize, i+2, n, m-1];
                        }
                        ptr = ptr + 5;
                    }
                    
                    for(int n = 0; n < 5; n++) 
                    {
                        in_buffer_y[ptr + n] = rhs[c, k+2, jsize, i+2, n];
                    }
                    ptr = ptr + 5;
                }
            }
            
            return 0;
		}
	}
}
