using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.ZAxis;
using bt.solve.BTMethod;
using bt.solve.Solver;

namespace impl.bt.solve.ZSolver 
{ 
	public class IZSolverImpl<I, C, MTH, DIR> : BaseIZSolverImpl<I, C, MTH, DIR>, IBTSolver<I, C, MTH, DIR>
	where I:IInstance_BT<C>
	where C:IClass
	where DIR:IZ
	where MTH:IBTMethod 
	{
		private double[][] in_buffer_solver;
		private double[][] out_buffer_solver;
		
		private void create_buffers_solvers()
		{
		    in_buffer_solver = new double[2][];
			out_buffer_solver = new double[2][];
			int buffer_size;
			
			buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * (5 * 5 + 5);
			in_buffer_solver[0] = new double[buffer_size];
			out_buffer_solver[0] = new double[buffer_size];
			
			buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * 5;
			in_buffer_solver[1] = new double[buffer_size];
			out_buffer_solver[1] = new double[buffer_size];
		}

		private bool buffers_ok = false;
		
		public override int go() 
        {		    
			if (!buffers_ok)
			{
				this.create_buffers_solvers();
				buffers_ok = true;
			}

			//int c, stage, first, last;
						
			Input_buffer.Array = out_buffer_solver[0] ;
			
			Iteration_control_forward.start();
			
	        Solve_cell.go();
						
			while (!Iteration_control_forward.is_last_stage())
			{
		        Output_buffer.Array = in_buffer_solver[0]; // new double[buffer_size];
				
		        Pack_solve_info.go();
				
		        Shift_lr.initiate_send();
				
				Iteration_control_forward.advance();
				
		        Shift_lr.initiate_recv();
		        Shift_lr.go();
				
		        Unpack_solve_info.go();
				
		        Solve_cell.go();
			}
			
			Iteration_control_forward.end();
			
			Input_buffer.Array = out_buffer_solver[1]; // new double[buffer_size];
			
			Iteration_control_backward.start();
			
			Back_substitute.go();			
			
			while (!Iteration_control_backward.is_first_stage())
			{
		        Output_buffer.Array = in_buffer_solver[1]; // new double[buffer_size];
				
		        Pack_back_sub_info.go();
				
		        Shift_rl.initiate_send();
				
				Iteration_control_backward.advance();
			        
			    Shift_rl.initiate_recv();
		        Shift_rl.go();		        
				
		        Unpack_back_sub_info.go();
				
		        Back_substitute.go();			    
			}
			
			Iteration_control_backward.end ();
			
			return 0;
		}
	}
}
