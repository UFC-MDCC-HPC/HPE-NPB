using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using sp.solve.SPMethod;
using common.axis.YAxis;
using sp.solve.Solver;

namespace impl.sp.solve.YSolver { 

	public class IYSolverImpl<I, C, MTH, DIR> : BaseIYSolverImpl<I, C, MTH, DIR>, ISPSolver<I, C, MTH, DIR>
		where I:IInstance_SP<C>
		where C:IClass
		where MTH:ISPMethod
		where DIR:IY
	{		
		public override int go() 
		{		
            //---------------------------------------------------------------------
            //                          FORWARD ELIMINATION  
            //---------------------------------------------------------------------
			
			Forward.begin();
			Lhs.begin ();
			Read_buffer_forward.begin();
			Write_buffer_forward.begin();
			
			Lhs.go(); 
			Lhs.advance ();
			
			Forward.go(); 
			Forward.advance();
				
			int k=0;
			while (!Forward.finished())
		    {
				Console.WriteLine("{0}: BEGIN LOOP FORWARD-Y {1}", this.Rank, k);
				Write_buffer_forward.go();
				
				Read_buffer_forward.advance ();
				
				Shift_forward.initiate_send();
				Shift_forward.initiate_recv();					
				
				Lhs.go();					
				Lhs.advance ();		
				
				Shift_forward.go();		            
				
				Write_buffer_forward.advance();
				
				Read_buffer_forward.go ();
				
				Forward.go();
				Forward.advance();
				Console.WriteLine("{0}: END LOOP FORWARD-Y {1}", this.Rank, k++);
				
		    } // cells loop
		    
		
		    //---------------------------------------------------------------------
		    //                         BACKWARD SUBSTITUTION 
		    //---------------------------------------------------------------------
		
			Backward.begin();
			Matvecproduct.begin();
			Read_buffer_backward.begin();
			Write_buffer_backward.begin();
			
			Backward.init();
			Backward.go();
			Backward.advance();
			
			k=0;
			while (!Backward.finished())
		    {				
				Console.WriteLine("{0}: BEGIN LOOP BACKWARD-Y {1}", this.Rank, k);
				Write_buffer_backward.go();					
				
				Read_buffer_backward.advance();
				
				Shift_backward.initiate_send();
				Shift_backward.initiate_recv();
	
				Matvecproduct.go();
			    Matvecproduct.advance();
	
				Shift_backward.go();
				
				Write_buffer_backward.advance();
	
				Read_buffer_backward.go();
				
				Backward.go();
		        Backward.advance ();								
				Console.WriteLine("{0}: END LOOP BACKWARD-Y {1}", this.Rank, k++);
		    }
	        //---------------------------------------------------------------------
	        //         If this was the last stage, do the block-diagonal inversion          
	        //---------------------------------------------------------------------
            //Matvecproduct.advance();
		    Matvecproduct.go();
					
			return 0;
		}  
		
	}

}
