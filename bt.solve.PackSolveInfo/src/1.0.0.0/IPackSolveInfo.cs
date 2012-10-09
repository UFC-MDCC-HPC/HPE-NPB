using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.axis.Axis;
using common.solve.SolvingMethod;

namespace bt.solve.PackSolveInfo { 
	public interface IPackSolveInfo<I, C, DIR, MTH> : BaseIPackSolveInfo<I, C, DIR, MTH>
	where I:IInstance<C>
	where C:IClass
	where DIR:IAxis
	where MTH:ISolvingMethod 
	{
	} // end main interface 
} // end namespace 
