/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.datapartition.MultiPartition;
using environment.MPIDirect;
using adi.Verify;
using common.datapartition.MultiPartitionCells;
using common.topology.Ring;
using adi.data.Initialize;
using adi.data.ExactRHS;
using adi.data.LHSInit;
using common.benchmarking.Timer;
using adi.ADI_Solver3D;
using common.solve.SolvingMethod;
using adi.ADI;

namespace impl.adi.ADI_Solver3DImpl { 

public abstract class BaseIADI_Solver3DImpl<MTH, CLASS, I>: Application, BaseIADI_Solver3D<MTH, CLASS, I>
	where MTH:ISolvingMethod
	where CLASS:IClass
	where I:IInstance<CLASS>
{
					
#region MyRegion
		
public PROBLEM_CLASS problem_class;
		
protected int ncells;
protected int[,] cell_size;
protected int[] grid_points;		
protected int problem_size;
protected int niter;
		
override public void initialize()
{	
	cell_size = Cells.cell_size;		
	ncells = Cells.ncells;			
	
	problem_size = Instance.problem_size;			
	problem_class = Instance.CLASS;	
	niter = Instance.niter_default;
	
	grid_points = Problem.grid_points;			
}
		
#endregion
		
		
private ITimer timer = null;

protected ITimer Timer {
	get {
		if (this.timer == null)
			this.timer = (ITimer) Services.getPort("timer");
		return this.timer;
	}
}

private IInitialize<I, CLASS> initialize_ = null;

protected IInitialize<I, CLASS> Initialize {
	get {
		if (this.initialize_ == null)
			this.initialize_ = (IInitialize<I, CLASS>) Services.getPort("initialize");
		return this.initialize_;
	}
}

private IExactRHS<I, CLASS> exact_rhs = null;

protected IExactRHS<I, CLASS> Exact_rhs {
	get {
		if (this.exact_rhs == null)
			this.exact_rhs = (IExactRHS<I, CLASS>) Services.getPort("exact_rhs");
		return this.exact_rhs;
	}
}

private ILHSInit<I, CLASS> lhsinit = null;

protected ILHSInit<I, CLASS> Lhsinit {
	get {
		if (this.lhsinit == null)
			this.lhsinit = (ILHSInit<I, CLASS>) Services.getPort("lhsinit");
		return this.lhsinit;
	}
}
		
		
private ICells cells = null;

protected ICells Cells {
	get {
		if (this.cells == null) 
		{
			this.cells = (ICells) Services.getPort("cells_info");
		}
		return this.cells;
	}
}

private ICell x = null;

protected ICell X {
	get {
		if (this.x == null)
			this.x = (ICell) Services.getPort("x");
		return this.x;
	}
}

private ICell y = null;

protected ICell Y {
	get {
		if (this.y == null)
			this.y = (ICell) Services.getPort("y");
		return this.y;
	}
}

private ICell z = null;

protected ICell Z {
	get {
		if (this.z == null)
			this.z = (ICell) Services.getPort("z");
		return this.z;
	}
}
		
		
private I instance = default(I);

protected I Instance {
	get {
		if (instance==null) 
		{
			this.instance = (I) Services.getPort("instance");
		}
		return instance;
	}
}		
		

private IProblemDefinition<I, CLASS> problem = null;

protected IProblemDefinition<I, CLASS> Problem {
	get {
		if (this.problem == null) 
		{
			this.problem = (IProblemDefinition<I, CLASS>) Services.getPort("problem_data");
		}
		return this.problem;
	}
}

private IMultiPartition<I, CLASS> process = null;

protected IMultiPartition<I, CLASS> Process {
	get {
		if (this.process == null)
			this.process = (IMultiPartition<I, CLASS>) Services.getPort("data_partition");
		return this.process;
	}
}

private IMPIDirect mpi = null;

protected IMPIDirect Mpi {
	get {
		if (this.mpi == null)
			this.mpi = (IMPIDirect) Services.getPort("mpi");
		return this.mpi;
	}
}

private IVerify<I, CLASS> verify = null;

protected IVerify<I, CLASS> Verify {
	get {
		if (this.verify == null)
			this.verify = (IVerify<I, CLASS>) Services.getPort("verify");
		return this.verify;
	}
}

private IADI<MTH, CLASS, I> adi = null;

protected IADI<MTH, CLASS, I> Adi {
	get {
		if (this.adi == null)
			this.adi = (IADI<MTH, CLASS, I>) Services.getPort("adi");
		return this.adi;
	}
}


abstract public int go(); 


}

}
