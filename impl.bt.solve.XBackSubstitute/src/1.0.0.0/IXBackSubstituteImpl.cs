using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.XAxis;
using bt.solve.BTMethod;
using bt.solve.BackSubstitute;

namespace impl.bt.solve.XBackSubstitute 
{ 

public class IXBackSubstituteImpl<I, C, DIR, MTH> : BaseIXBackSubstituteImpl<I, C, DIR, MTH>, IBackSubstitute<I, C, DIR, MTH>
where I:IInstance_BT<C>
where C:IClass
where DIR:IX
where MTH:IBTMethod
{
	private int last, c;
	
	public override int go()   
	{ 
		lhsc = Lhsc.Field6;
		backsub_info = Backsub_info.Field;
			
		int i, j, k;
		int m, n, isize, jsize, ksize, istart, stage;
			
		stage = this.Iteration_control.getCurrentStage();
		last = this.Iteration_control.is_last_stage() ? 1 : 0;
			
		c = slice[stage,0];		
		
		istart = 2;
		isize = cell_size[c, 0] + 1;
		jsize = cell_size[c, 1] - end[c, 1] + 1;
		ksize = cell_size[c, 2] - end[c, 2] + 1;
		
		if(last == 0) 
		{
		    for(k = start[c, 2]; k <= ksize; k++) 
		    {
		        for(j = start[c, 1]; j <= jsize; j++) 
		        {
		            //---------------------------------------------------------------------
		            //     U[isize] uses info from previous cell if not last cell
		            //---------------------------------------------------------------------
		            for(m = 0; m < 5; m++) 
		            {
		                for(n = 0; n < 5; n++) 
		                {//rhs[m,isize,j,k,c] = rhs[m,isize,j,k,c] - lhsc[m,n,isize,j,k,c]*backsub_info[n,j,k,c]
		                    rhs[c, k, j, isize, m] = rhs[c, k, j, isize, m] - lhsc[c, k, j, isize, n, m] * backsub_info[c, k, j, n, 0];
		                }
		            }
		        }
		    }
		}
		
		for(k = start[c, 2]; k <= ksize; k++) 
		{
		    for(j = start[c, 1]; j <= jsize; j++) 
		    {
		        for(i = isize-1; i >= istart; i--) 
		        {  //for(i=isize-1,istart,-1;
		            for(m = 0; m < 5; m++) 
		            {
		                for(n = 0; n < 5; n++) 
		                { //rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - lhsc[m,n,i,j,k,c]*rhs[n,i+1,j,k,c];
		                    rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - lhsc[c, k, j, i, n, m] * rhs[c, k, j, i+1, n];
		                }
		            }
		        }
		    }
		}
	
	   return 0;	   
	} // end activate method 

}

}
