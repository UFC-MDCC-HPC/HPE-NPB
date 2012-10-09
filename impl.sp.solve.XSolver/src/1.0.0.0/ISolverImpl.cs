using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using sp.solve.SPMethod;
using common.axis.Axis;
using sp.solve.Solver;

namespace impl.sp.solve.SolverImpl { 
	
	public class ISolverImpl<I, C, MTH, DIR> : BaseISolverImpl<I, C, MTH, DIR>, ISPSolver<I, C, MTH, DIR>
		where I:IInstance_SP<C>
		where C:IClass
		where MTH:ISPMethod
		where DIR:IAxis
	{		
		public override int go() 
		{		 
//		    Console.WriteLine("{0}: sucessor = {1}, predecessor = {2}",this.Rank, Topology.successor, Topology.predecessor);
			
            //---------------------------------------------------------------------
            //                          FORWARD ELIMINATION  
            //---------------------------------------------------------------------
			
			Forward.begin();
			Lhs.begin ();
			Read_buffer_forward.begin();
			Write_buffer_forward.begin();
			
			//Console.WriteLine("{0}: FIRST LHS !!!!", this.Rank);
			Lhs.go();       // lhs.stage = 0
			Lhs.advance (); // lhs.stage = 1
			
			Forward.go();       // lhs.stage=0
			Forward.advance();  // lhs.stage = 1
				
		int k=0;
			while (!Forward.finished())
		    {				
				//Console.WriteLine("{0}: BEGIN LOOP FORWARD {1}", this.Rank, k);
				
				Write_buffer_forward.go();         // write_buffer_foward.stage = 0
				Write_buffer_forward.advance();    // write_buffer_forward.stage = 1
								
				Read_buffer_forward.advance ();   // read_buffer_forward.stage = 1
				
				Shift_forward.initiate_send();    
				Shift_forward.initiate_recv();					
				
				Lhs.go();					      // lhs.stage = 1
				Lhs.advance ();		              // lhs.stage = 2
				
				Shift_forward.go();		            
				
				Read_buffer_forward.go ();		  // read_buffer_forward.stage = 1		
				
				Forward.go();                     // forward.stage = 1
				Forward.advance();                // forward.stage = 2
				//Console.WriteLine("{0}: END LOOP FORWARD {1}", this.Rank, k++);
				
		    } // cells loop
		    
		
		    //---------------------------------------------------------------------
		    //                         BACKWARD SUBSTITUTION 
		    //---------------------------------------------------------------------
		
			Backward.begin();
			Matvecproduct.begin();
			Read_buffer_backward.begin();
			Write_buffer_backward.begin();
			
			Backward.init();        // backward.stage = 1
			Backward.go();          // backward.stage = 1
			Backward.advance();     // backward.stage = 0
			
			k=0;
			while (!Backward.finished())
		    {				
				//Console.WriteLine("{0}: BEGIN LOOP BACKWARD {1}", this.Rank, k);
				Write_buffer_backward.go();		  // write_buffer_backward.stage = 1

				Write_buffer_backward.advance();  // write_buffer_backward.stage = 0
				Read_buffer_backward.advance();   // read_buffer_backward = 0
				
				Shift_backward.initiate_send();   
				Shift_backward.initiate_recv();
	
				Matvecproduct.go();               // matvecproduct.stage = 1
			    Matvecproduct.advance();          // matvecproduct.stage = 0
	
				Shift_backward.go();
	
				Read_buffer_backward.go();        // read_buffer_backward = 0
				
				Backward.go();                    // backward.stage = 0
		        Backward.advance ();              // backward.stage = -1
				//Console.WriteLine("{0}: END LOOP BACKWARD {1}", this.Rank, k++);
		    }	
	        //---------------------------------------------------------------------
	        //         If this was the last stage, do the block-diagonal inversion          
	        //---------------------------------------------------------------------
        //    Matvecproduct.advance();
		    Matvecproduct.go(); // matvecproduct.stage = 0
					
			return 0;
		}  
	}
}
