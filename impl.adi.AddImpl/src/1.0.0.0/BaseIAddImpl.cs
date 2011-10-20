/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using adi.Add;
using common.problem_size.Class;
using common.problem_size.Instance;


namespace impl.adi.AddImpl { 

	public abstract class BaseIAddImpl<I,C>: Computation, BaseIAdd<I,C>
			where I:IInstance<C>
			where C:IClass
	{
			
	#region data
			
	protected int ncells;
	protected double[,,,,] u;
	protected double[,,,,] rhs;
	protected int[,] start;
	protected int[,] end;
	protected int[,] cell_size;
			
	override public void initialize() 
	{
		start = Cells.cell_start;
		end = Cells.cell_end;
		cell_size = Cells.cell_size;   
				
	    ncells = Problem.NCells;			
		u = Problem.Field_u;
		rhs = Problem.Field_rhs;
	}
			
	#endregion		
			
	private ICells cells = null;
	
	public ICells Cells {
		get {
			if (this.cells == null) 
			{
			    this.cells = (ICells) Services.getPort("cells_info");
			}
			return this.cells;
		}
	}
	
	private IProblemDefinition<I,C> problem = null;
	
	public IProblemDefinition<I,C> Problem {
		get {
			if (this.problem == null)
			{
			   this.problem = (IProblemDefinition<I,C>) Services.getPort("problem_data");					
			}
			return this.problem;
		}
	}
	
	
	abstract public int go(); 
	
	
	}

}
