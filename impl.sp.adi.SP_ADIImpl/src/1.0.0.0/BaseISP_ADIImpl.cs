/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using common.topology.Ring;
using adi.CopyFaces;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using adi.Solver;
using sp.solve.SPMethod;
using common.axis.YAxis;
using common.axis.ZAxis;
using sp.solve.BlockDiagonalMatVecProduct;
using common.axis.XYZAxes;
using environment.MPIDirect;
using adi.Add;
using adi.data.ProblemDefinition;
using adi.ComputeRHS;
using common.axis.XAxis;
using sp.adi.SP_ADI;


namespace impl.sp.adi.SP_ADIImpl { 

public abstract class BaseISP_ADIImpl<MTH, C, I>: Computation, BaseISP_ADI<MTH, C, I>
where MTH:ISPMethod
where C:IClass
where I:IInstance_SP<C>
{

private ICells cells = null;

public ICells Cells {
	get {
		if (this.cells == null)
			this.cells = (ICells) Services.getPort("cells_info");
		return this.cells;
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

private ICell y = null;

public ICell Y {
	get {
		if (this.y == null)
			this.y = (ICell) Services.getPort("y");
		return this.y;
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

private ICopyFaces<I, C> copy_faces = null;

protected ICopyFaces<I, C> Copy_faces {
	get {
		if (this.copy_faces == null)
			this.copy_faces = (ICopyFaces<I, C>) Services.getPort("copy_faces");
		return this.copy_faces;
	}
}

private ISolver<I, C, MTH, IX> x_solve = null;

protected ISolver<I, C, MTH, IX> X_solve {
	get {
		if (this.x_solve == null)
			this.x_solve = (ISolver<I, C, MTH, IX>) Services.getPort("x_solve");
		return this.x_solve;
	}
}

private ISolver<I, C, MTH, IY> y_solve = null;

protected ISolver<I, C, MTH, IY> Y_solve {
	get {
		if (this.y_solve == null)
			this.y_solve = (ISolver<I, C, MTH, IY>) Services.getPort("y_solve");
		return this.y_solve;
	}
}

private ISolver<I, C, MTH, IZ> z_solve = null;

protected ISolver<I, C, MTH, IZ> Z_solve {
	get {
		if (this.z_solve == null)
			this.z_solve = (ISolver<I, C, MTH, IZ>) Services.getPort("z_solve");
		return this.z_solve;
	}
}

private IBlockDiagonalMatVecProduct<I, C, IXYZ, MTH> txinvr = null;

protected IBlockDiagonalMatVecProduct<I, C, IXYZ, MTH> Txinvr {
	get {
		if (this.txinvr == null)
			this.txinvr = (IBlockDiagonalMatVecProduct<I, C, IXYZ, MTH>) Services.getPort("txinvr");
		return this.txinvr;
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

private IAdd<I, C> add = null;

protected IAdd<I, C> Add {
	get {
		if (this.add == null)
			this.add = (IAdd<I, C>) Services.getPort("add");
		return this.add;
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

private IComputeRHS<I, C> compute_rhs = null;

protected IComputeRHS<I, C> Compute_rhs {
	get {
		if (this.compute_rhs == null)
			this.compute_rhs = (IComputeRHS<I, C>) Services.getPort("compute_rhs");
		return this.compute_rhs;
	}
}


abstract public int go(); 


}

}
