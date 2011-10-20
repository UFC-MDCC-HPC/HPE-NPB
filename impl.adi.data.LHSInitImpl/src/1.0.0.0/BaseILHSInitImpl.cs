/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using adi.data.LHSInit;

namespace impl.adi.data.LHSInitImpl { 

public abstract class BaseILHSInitImpl<I, C>: Computation, BaseILHSInit<I, C>
where I:IInstance<C>
where C:IClass
{
		
#region data		
		
protected int[,] cell_size, cell_low, cell_high, start, end, slice, cell_coord;		
protected int ncells, IMAX, JMAX, KMAX, maxcells;		
protected double[,,,,] lhs;

override public void initialize()
{	
	cell_size = Cells.cell_size;
	cell_low = Cells.cell_low;
	cell_high = Cells.cell_high;
	cell_coord = Cells.cell_coord;
	start = Cells.cell_start;
	end = Cells.cell_end;
	slice = Cells.cell_slice;
	
	ncells = Problem.NCells;
	lhs = Problem.Field_lhs;			
				
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


abstract public int go(); 


}

}
