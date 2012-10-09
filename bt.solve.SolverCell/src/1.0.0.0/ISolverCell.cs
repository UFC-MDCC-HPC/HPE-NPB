using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.axis.Axis;
using common.solve.SolvingMethod;

namespace bt.solve.SolverCell { 
	public interface ISolverCell<I, C, DIR, MTH> : BaseISolverCell<I, C, DIR, MTH>
	where I:IInstance<C>
	where C:IClass
	where DIR:IAxis
	where MTH:ISolvingMethod 
	{
	}
}
