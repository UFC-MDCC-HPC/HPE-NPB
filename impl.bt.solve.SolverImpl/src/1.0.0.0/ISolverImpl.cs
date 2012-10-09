using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.Axis;
using bt.solve.BTMethod;
using bt.solve.Solver;

namespace impl.bt.solve.SolverImpl 
{ 
	public class ISolverImpl<I, C, MTH, DIR> : BaseISolverImpl<I, C, MTH, DIR>, IBTSolver<I, C, MTH, DIR>
		where I:IInstance_BT<C>
		where C:IClass
		where DIR:IAxis
		where MTH:IBTMethod 
	{
		
		public override int go() 
        {		    
			Iteration_control_forward.start();		
			
	        Solve_cell.go();
						
			while (!Iteration_control_forward.is_last_stage())
			{
		        Pack_solve_info.go();				
		        Shift_forward.initiate_send();				
				
				Iteration_control_forward.advance();				
				
		        Shift_forward.initiate_recv();
		        Shift_forward.go();				
		        Unpack_solve_info.go();				
		        Solve_cell.go();
			}
			
			Iteration_control_forward.end();
			
			Iteration_control_backward.start();			
			Back_substitute.go();			
			
			while (!Iteration_control_backward.is_first_stage())
			{
		        Pack_back_sub_info.go();				
		        Shift_backward.initiate_send();				
				Iteration_control_backward.advance();			        
			    Shift_backward.initiate_recv();
		        Shift_backward.go();		        				
		        Unpack_back_sub_info.go();				
		        Back_substitute.go();			    
			}
			
			Iteration_control_backward.end ();
			
			return 0;
		}
	}
}

