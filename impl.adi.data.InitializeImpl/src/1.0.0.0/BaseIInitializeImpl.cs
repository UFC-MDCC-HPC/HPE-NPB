/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using adi.data.ExactSolution;
using adi.data.Initialize;

namespace impl.adi.data.InitializeImpl { 

public abstract class BaseIInitializeImpl<I, C>: Computation, BaseIInitialize<I, C>
where I:IInstance<C>
where C:IClass
{
		
#region data		
		
protected int[,] cell_size, cell_low, cell_high, start, end, slice;		
protected int ncells, IMAX, JMAX, KMAX, maxcells;		
protected double[,,,,] u;		
protected double dnzm1, dnym1, dnxm1; 

		
override public void initialize()
{
	cell_size = Cells.cell_size;
	cell_low = Cells.cell_low;
	cell_high = Cells.cell_high;
	start = Cells.cell_start;
	end = Cells.cell_end;
	slice = Cells.cell_slice;
	ncells = Problem.NCells;
	
	u = Problem.Field_u;
				
	dnzm1 = Constants.dnzm1;
	dnym1 = Constants.dnym1;
	dnxm1 = Constants.dnxm1;
	
	IMAX = Problem.IMAX;
	JMAX = Problem.JMAX;
	KMAX = Problem.KMAX;
	maxcells = Problem.maxcells;
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

private IProblemDefinition<I, C> problem = null;

public IProblemDefinition<I, C> Problem {
	get {
		if (this.problem == null)
		{
			this.problem = (IProblemDefinition<I, C>) Services.getPort("problem_data");
		}	
		return this.problem;
	}
}

private IExactSolution exact_solution = null;

protected IExactSolution Exact_solution {
	get {
		if (this.exact_solution == null)
			this.exact_solution = (IExactSolution) Services.getPort("exact_solution");
		return this.exact_solution;
	}
}


abstract public int go(); 


}

}
