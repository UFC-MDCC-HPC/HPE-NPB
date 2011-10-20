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
using common.axis.ZAxis;
using bt.solve.BTMethod;
using bt.solve.UnpackBackSubInfo;

namespace impl.bt.solve.ZUnpackBackSubInfo { 

public abstract class BaseIZUnpackBackSubInfo<I, C, DIR, MTH>: Computation, BaseIUnpackBackSubInfo<I, C, DIR, MTH>
where I:IInstance_BT<C>
where C:IClass
where DIR:IZ
where MTH:IBTMethod
{
#region data
	protected int JMAX;
	protected int IMAX;

	override public void initialize(){
		JMAX = Problem.JMAX;
		IMAX = Problem.IMAX;
	}
#endregion

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
