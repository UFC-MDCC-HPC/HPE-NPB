using br.ufc.pargo.hpe.kinds;
using common.solve.SolvingMethod;
using common.axis.Axis;
using common.problem_size.Class;
using common.problem_size.Instance;

namespace adi.Solver { 

	public interface ISolver<I, C, MTH, DIR> : BaseISolver<I, C, MTH, DIR>
		where MTH:ISolvingMethod
		where DIR:IAxis
		where C:IClass
		where I:IInstance<C>
	{
	
	
	}  

}  
