/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using adi.data.ProblemDefinition;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.Buffer;
using common.axis.Axis;
using common.direction.Backward;
using common.axis.YAxis;
using sp.solve.shift.WriteBuffer;
using common.datapartition.MultiPartitionCells;

namespace impl.sp.solve.shift.YWriteBufferBackward { 

public abstract class BaseIYWriteBufferBackwardImpl<I, C, D, DIR>: Computation, BaseIWriteBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IBackwardDirection
where DIR:IY
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
	
	ncells = Problem.NCells;
	lhs = Problem.Field_lhs;
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
		
private D direction = default(D);

protected D Direction {
	get {
		if (this.direction == null)
			this.direction = (D) Services.getPort("direction");
		return this.direction;
	}
}

private IProblemDefinition<I, C> problem_data_2 = null;

public IProblemDefinition<I, C> Problem {
	get {
		if (this.problem_data_2 == null)
			this.problem_data_2 = (IProblemDefinition<I, C>) Services.getPort("problem_data_2");
		return this.problem_data_2;
	}
}

private IBuffer buffer = null;

public IBuffer Buffer {
	get {
		if (this.buffer == null)
			this.buffer = (IBuffer) Services.getPort("buffer");
		return this.buffer;
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


abstract public int go(); 


}

}
