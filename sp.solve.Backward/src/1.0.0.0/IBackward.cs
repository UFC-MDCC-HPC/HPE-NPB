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
	void init();         // action
	void begin();        // action
	void advance();      // action
	bool finished();     // condition
	bool first_stage();  // condition
	bool last_stage();   // condition
		
} // end main interface 

} // end namespace 
