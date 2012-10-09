/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.axis.Axis;
using common.axis.ZAxis;
using sp.solve.SPMethod;
using sp.solve.Backward;

namespace impl.sp.solve.ZBackwardImpl { 

public abstract class BaseIZBackwardImpl<I, C, DIR, MTH>: Computation, BaseIBackward<I, C, DIR, MTH>
where I:IInstance_SP<C>
where C:IClass
where DIR:IZ
where MTH:ISPMethod
{
		
#region data
		
protected int[,] start, end, slice, cell_size;
protected double[,,,,] lhs, rhs;
protected int ncells;
		
override public void initialize()
{	
	start = Cells.cell_start;
	end = Cells.cell_end;
	slice = Cells.cell_slice;
	cell_size = Cells.cell_size;
	
	lhs = Problem.Field_lhs;
	rhs = Problem.Field_rhs;
	ncells = Problem.NCells;
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

private DIR axis = default(DIR);

protected DIR Axis {
	get {
		if (this.axis == null)
			this.axis = (DIR) Services.getPort("axis");
		return this.axis;
	}
}

private MTH method = default(MTH);

protected MTH Method {
	get {
		if (this.method == null)
			this.method = (MTH) Services.getPort("method");
		return this.method;
	}
}


abstract public int go(); 


}

}
