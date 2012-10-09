using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.axis.Axis;
using common.solve.SolvingMethod;

namespace sp.solve.LHS { 

public interface ILHS<I, C, DIR, MTH> : BaseILHS<I, C, DIR, MTH>
where I:IInstance<C>
where C:IClass
where DIR:IAxis
where MTH:ISolvingMethod
{
	void begin();
	void advance();
	bool finished();
} // end main interface 

} // end namespace 
