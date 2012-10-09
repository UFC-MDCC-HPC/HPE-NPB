using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.axis.Axis;
using common.solve.SolvingMethod;

namespace bt.solve.UnpackBackSubInfo { 
	public interface IUnpackBackSubInfo<I, C, DIR, MTH> : BaseIUnpackBackSubInfo<I, C, DIR, MTH>
	where I:IInstance<C>
	where C:IClass
	where DIR:IAxis
	where MTH:ISolvingMethod 
	{
	} // end main interface 
} // end namespace 
