/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.axis.Axis;
using common.axis.YAxis;
using bt.solve.BTMethod;
using bt.solve.BackSubstitute;
using adi.solve.IterationControl;
using common.direction.Backward;
using common.data.Field;

namespace impl.bt.solve.YBackSubstitute { 

public abstract class BaseIYBackSubstituteImpl<I, C, DIR, MTH>: Computation, BaseIBackSubstitute<I, C, DIR, MTH>
where I:IInstance_BT<C>
where C:IClass
where DIR:IY
where MTH:IBTMethod
{
#region data
		
protected int[,] start, end, cell_size, slice;
protected double[,,,,] rhs;
protected double[,,,,,] lhsc;
protected double[,,,,] backsub_info;

override public void initialize()		
{
	start = Cells.cell_start;
	end = Cells.cell_end;
	cell_size = Cells.cell_size;
	slice = Cells.cell_slice;
			
	rhs = Problem.Field_rhs;
			
	Iteration_control.setNumberOfStages(Cells.ncells);		
			
			
	int MAX_CELL_DIM = Problem.MAX_CELL_DIM;
	int maxcells = Problem.maxcells;			
	Backsub_info.initialize_field("backsub_info", maxcells, MAX_CELL_DIM+3, MAX_CELL_DIM+3, 5);						
	
}
		
#endregion  

private IField backsub_info_ = null;

protected IField Backsub_info {
	get {
		if (backsub_info_ == null) 
			backsub_info_ = (IField) Services.getPort("backsub_info");
		return backsub_info_;
	}
}
		
private IField lhsc_ = null;

protected IField Lhsc {
	get {
		if (lhsc_ == null) 
			lhsc_ = (IField) Services.getPort("lhsc");
		return lhsc_;
	}
}
		
		
private IIterationControl<IBackwardDirection> iteration_control = null;

public IIterationControl<IBackwardDirection> Iteration_control {
	get {
		if (this.iteration_control == null)
			this.iteration_control = (IIterationControl<IBackwardDirection>) Services.getPort("iteration_control");
		return this.iteration_control;
	}
}

private ICells cells = null;

public ICells Cells {
	get {
		if (this.cells == null)
			this.cells = (ICells) Services.getPort("cells_info");
		return this.cells;
	}
}

private IProblemDefinition<I, C> problem = null;

public IProblemDefinition<I, C> Problem {
	get {
		if (this.problem == null)
			this.problem = (IProblemDefinition<I, C>) Services.getPort("problem_data");
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
