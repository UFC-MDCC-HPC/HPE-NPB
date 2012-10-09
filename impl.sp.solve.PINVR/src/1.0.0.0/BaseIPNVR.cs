/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using sp.solve.BlockDiagonalMatVecProduct;
using common.axis.Axis;
using sp.solve.SPMethod;
using common.axis.YAxis;

namespace impl.sp.solve.PINVR { 

public abstract class BaseIPNVR<I,C,DIR,MTH>: Computation, BaseIBlockDiagonalMatVecProduct<I,C,DIR,MTH>
	where I:IInstance<C>
	where C:IClass
	where DIR:IY
	where MTH:ISPMethod
{
		
#region data		
		
protected int[,] start, end, cell_size, slice;
protected double[,,,,] rhs;
protected double bt;
protected int ncells;		
		
override public void initialize()
{
	start = Cells.cell_start;
	end = Cells.cell_end;
	cell_size = Cells.cell_size;
	slice = Cells.cell_slice;					
	
	bt = Constants.bt;					
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
