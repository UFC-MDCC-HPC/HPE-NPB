/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using common.topology.Ring;
using adi.CopyFaces;
using common.problem_size.Class;
using adi.ComputeRHS;
using adi.Solver;
using common.axis.ZAxis;
using common.direction.Direction;
using environment.MPIDirect;
using common.axis.XAxis;
using common.axis.YAxis;
using adi.Add;
using adi.data.ProblemDefinition;
using adi.ADI;
using common.solve.SolvingMethod;
using common.problem_size.Instance;

namespace impl.adi.ADIImpl { 

public abstract class BaseIADIImpl<I, CLASS, MTH>: Computation, BaseIADI<I, CLASS, MTH>
	where MTH:ISolvingMethod
	where CLASS:IClass
	where I:IInstance<CLASS>
{

private ICells cells = null;

public ICells Cells {
	get {
		if (this.cells == null)
			this.cells = (ICells) Services.getPort("cells_info");
		return this.cells;
	}
}

private ICell cell = null;

public ICell Cell {
	get {
		if (this.cell == null)
			this.cell = (ICell) Services.getPort("topology");
		return this.cell;
	}
}

private ICell x = null;

public ICell X {
	get {
		if (this.x == null)
			this.x = (ICell) Services.getPort("x");
		return this.x;
	}
}

private ICell z = null;

public ICell Z {
	get {
		if (this.z == null)
			this.z = (ICell) Services.getPort("z");
		return this.z;
	}
}

private ICell y = null;

public ICell Y {
	get {
		if (this.y == null)
			this.y = (ICell) Services.getPort("y");
		return this.y;
	}
}

private ICopyFaces<I, CLASS> copy_faces = null;

protected ICopyFaces<I, CLASS> Copy_faces {
	get {
		if (this.copy_faces == null)
			this.copy_faces = (ICopyFaces<I, CLASS>) Services.getPort("copy_faces");
		return this.copy_faces;
	}
}

private IComputeRHS<I, CLASS> compute_rhs = null;

protected IComputeRHS<I, CLASS> Compute_rhs {
	get {
		if (this.compute_rhs == null)
			this.compute_rhs = (IComputeRHS<I, CLASS>) Services.getPort("compute_rhs");
		return this.compute_rhs;
	}
}

private ISolver<I, CLASS, MTH, IZ> z_solve = null;

protected ISolver<I, CLASS, MTH, IZ> Z_solve {
	get {
		if (this.z_solve == null)
			this.z_solve = (ISolver<I, CLASS, MTH, IZ>) Services.getPort("z_solve");
		return this.z_solve;
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

private ISolver<I, CLASS, MTH, IX> x_solve = null;

protected ISolver<I, CLASS, MTH, IX> X_solve {
	get {
		if (this.x_solve == null)
			this.x_solve = (ISolver<I, CLASS, MTH, IX>) Services.getPort("x_solve");
		return this.x_solve;
	}
}

private ISolver<I, CLASS, MTH, IY> y_solve = null;

protected ISolver<I, CLASS, MTH, IY> Y_solve {
	get {
		if (this.y_solve == null)
			this.y_solve = (ISolver<I, CLASS, MTH, IY>) Services.getPort("y_solve");
		return this.y_solve;
	}
}

private IAdd<I, CLASS> add = null;

protected IAdd<I, CLASS> Add {
	get {
		if (this.add == null)
			this.add = (IAdd<I, CLASS>) Services.getPort("add");
		return this.add;
	}
}

private IProblemDefinition<I, CLASS> problem = null;

public IProblemDefinition<I, CLASS> Problem {
	get {
		if (this.problem == null)
			this.problem = (IProblemDefinition<I, CLASS>) Services.getPort("problem_data");
		return this.problem;
	}
}


abstract public int go(); 


}

}
