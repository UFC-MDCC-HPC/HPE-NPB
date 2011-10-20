/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using common.topology.Ring;
using sp.solve.LHS;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using adi.data.ProblemDefinition;
using sp.solve.BlockDiagonalMatVecProduct;
using common.Buffer;
using sp.solve.Forward;
using sp.solve.SPMethod;
using common.axis.YAxis;
using common.interactionpattern.Shift;
using common.direction.Forward;
using common.direction.Backward;
using environment.MPIDirect;
using common.problem_size.Instance;
using sp.solve.Backward;
using sp.solve.Solver;

namespace impl.sp.solve.YSolver { 

public abstract class BaseIYSolverImpl<I, C, MTH, DIR>: Computation, BaseISPSolver<I, C, MTH, DIR>
where I:IInstance_SP<C>
where C:IClass
where MTH:ISPMethod
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

private ICell cell = null;

public ICell Cell {
	get {
		if (this.cell == null)
		{
			this.cell = (ICell) Services.getPort("topology");					
		}
		return this.cell;
	}
}

private ILHS<I, C, DIR, MTH> lhs_ = null;

protected ILHS<I, C, DIR, MTH> Lhs {
	get {
		if (this.lhs_ == null)
			this.lhs_ = (ILHS<I, C, DIR, MTH>) Services.getPort("lhs");
		return this.lhs_;
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

private IBlockDiagonalMatVecProduct<I, C, DIR, MTH> matvecproduct = null;

protected IBlockDiagonalMatVecProduct<I, C, DIR, MTH> Matvecproduct {
	get {
		if (this.matvecproduct == null)
			this.matvecproduct = (IBlockDiagonalMatVecProduct<I, C, DIR, MTH>) Services.getPort("matvecproduct");
		return this.matvecproduct;
	}
}

private IBuffer output_buffer_forward = null;

protected IBuffer Output_buffer_forward {
	get {
		if (this.output_buffer_forward == null)
			this.output_buffer_forward = (IBuffer) Services.getPort("output_buffer_forward");
		return this.output_buffer_forward;
	}
}

private IBuffer input_buffer_forward = null;

protected IBuffer Input_buffer_forward {
	get {
		if (this.input_buffer_forward == null)
			this.input_buffer_forward = (IBuffer) Services.getPort("input_buffer_forward");
		return this.input_buffer_forward;
	}
}

private IBuffer output_buffer_backward = null;

protected IBuffer Output_buffer_backward {
	get {
		if (this.output_buffer_backward == null)
			this.output_buffer_backward = (IBuffer) Services.getPort("output_buffer_backward");
		return this.output_buffer_backward;
	}
}

private IBuffer input_buffer_backward = null;

protected IBuffer Input_buffer_backward {
	get {
		if (this.input_buffer_backward == null)
			this.input_buffer_backward = (IBuffer) Services.getPort("input_buffer_backward");
		return this.input_buffer_backward;
	}
}

private IForward<I, C, MTH, DIR> forward = null;

protected IForward<I, C, MTH, DIR> Forward {
	get {
		if (this.forward == null)
			this.forward = (IForward<I, C, MTH, DIR>) Services.getPort("forward");
		return this.forward;
	}
}

private IShift<IForwardDirection> shift_forward = null;

protected IShift<IForwardDirection> Shift_forward {
	get {
		if (this.shift_forward == null)
			this.shift_forward = (IShift<IForwardDirection>) Services.getPort("shift_forward");
		return this.shift_forward;
	}
}

private IShift<IBackwardDirection> shift_backward = null;

protected IShift<IBackwardDirection> Shift_backward {
	get {
		if (this.shift_backward == null)
			this.shift_backward = (IShift<IBackwardDirection>) Services.getPort("shift_backward");
		return this.shift_backward;
	}
}

private IMPIDirect mpi = null;

public IMPIDirect Mpi {
	get {
		if (this.mpi == null)
			this.mpi = (IMPIDirect) Services.getPort("mpi");
		return this.mpi;
	}
}


private IBackward<I, C, DIR, MTH> backward = null;

protected IBackward<I, C, DIR, MTH> Backward {
	get {
		if (this.backward == null)
			this.backward = (IBackward<I, C, DIR, MTH>) Services.getPort("backward");
		return this.backward;
	}
}


abstract public int go(); 


}

}
