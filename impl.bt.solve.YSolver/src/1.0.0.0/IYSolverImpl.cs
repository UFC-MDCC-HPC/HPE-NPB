using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.YAxis;
using bt.solve.BTMethod;
using bt.solve.Solver;

namespace impl.bt.solve.YSolver 
{ 
	public class IYSolverImpl<I, C, MTH, DIR> : BaseIYSolverImpl<I, C, MTH, DIR>, IBTSolver<I, C, MTH, DIR>
		where I:IInstance_BT<C>
		where C:IClass
		where DIR:IY
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

			int c, stage, first, last;
            //buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * (5 * 5 + 5);
            double[] out_buffer_y;// = new double[buffer_size];
            
            Input_buffer.Array = out_buffer_y = out_buffer_solver[0]; // new double[buffer_size];
            
            for(stage = 0; stage < ncells; stage++) 
            {
                c = slice[stage, 1];
                if(stage == ncells-1) 
                {
                    last = 1;
                }
                else 
                {
                    last = 0;
                }
                
                if(stage == 0) 
                {
                    first = 1;
                    Solve_cell.setParameters(lhsc, first, last, c);
                    Solve_cell.go();
                }                
                else 
                {
                    first = 0;
                   // int ip = cell_coord[c,0];
                   // int kp = cell_coord[c,2];
			        Shift_lr.initiate_recv();
			        Shift_lr.go();                    
                    //requests[0] = comm_solve.ImmediateReceive<double>(predecessor[1], SOUTH+ip+kp*ncells, out_buffer_y);
                    //requests[1].Wait();
                    //requests[0].Wait();
                    Unpack_solve_info.setParameters(lhsc, out_buffer_y, c);
                    Unpack_solve_info.go();
                    Solve_cell.setParameters(lhsc, first, last, c); 
                    Solve_cell.go();
                }
                
                if(last == 0) 
                {
                   // int ip = cell_coord[c,0];
                   // int kp = cell_coord[c,2];
                    double[] in_buffer_y = Output_buffer.Array = in_buffer_solver[0]; // new double[buffer_size];
                    Pack_solve_info.setParameters(lhsc, in_buffer_y, c);
                    Pack_solve_info.go();
			        Shift_lr.initiate_send();                    
                    //requests[1] = comm_solve.ImmediateSend<double>(in_buffer_y, successor[1], SOUTH+ip+kp*ncells);
                }
            }
            
            //out_buffer_y = null;
            Input_buffer.Array = out_buffer_y = null;
            //buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * 5;
            //out_buffer_y = new double[buffer_size];
            Input_buffer.Array = out_buffer_y = out_buffer_solver[1]; // new double[buffer_size];
            
            for(stage = ncells-1; stage >= 0; stage--) 
            {  
                c = slice[stage, 1];
                first = 0;
                last = 0;
                if(stage == 0)
                    first = 1;
                if(stage == ncells-1) 
                {
                    last = 1;
                    Back_substitute.setParameters(lhsc, backsub_info, first, last, c);
                    Back_substitute.go();
                }
                else 
                {
                    //int ip = cell_coord[c, 0];
                    //int kp = cell_coord[c, 2];
			        Shift_rl.initiate_recv();
			        Shift_rl.go();                    
                    //requests[0] = comm_solve.ImmediateReceive<double>(successor[1], NORTH+ip+kp*ncells, out_buffer_y);
                    //requests[1].Wait();
                    //requests[0].Wait();
                    Unpack_back_sub_info.setParameters(backsub_info, out_buffer_y, c);
                    Unpack_back_sub_info.go();
                    Back_substitute.setParameters(lhsc, backsub_info, first, last, c);
                    Back_substitute.go();
                }
                
                if(first == 0) 
                {
                   // int ip = cell_coord[c,0];
                    //int kp = cell_coord[c,2];                    
                    double[] in_buffer_y = Output_buffer.Array = in_buffer_solver[1]; // new double[buffer_size];
                    Pack_back_sub_info.setParameters(in_buffer_y, c);
                    Pack_back_sub_info.go(); 
			        Shift_rl.initiate_send();
                    //requests[1] = comm_solve.ImmediateSend<double>(in_buffer_y, predecessor[1], NORTH+ip+kp*ncells);
                }
            }
            
            return 0;
		}
	}
}