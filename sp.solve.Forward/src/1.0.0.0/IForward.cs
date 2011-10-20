using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.solve.SolvingMethod;
using common.axis.Axis;

namespace sp.solve.Forward { 

public interface IForward<I, C, MTH, DIR> : BaseIForward<I, C, MTH, DIR>
where I:IInstance<C>
where C:IClass
where MTH:ISolvingMethod
where DIR:IAxis
{
	void enterStage(int stage);

} // end main interface 

} // end namespace 
