using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.axis.Axis;
using common.solve.SolvingMethod;

namespace sp.solve.Backward { 

public interface IBackward<I, C, DIR, MTH> : BaseIBackward<I, C, DIR, MTH>
where I:IInstance<C>
where C:IClass
where DIR:IAxis
where MTH:ISolvingMethod
{
	void init();
	void enterStage(int stage);
} // end main interface 

} // end namespace 
