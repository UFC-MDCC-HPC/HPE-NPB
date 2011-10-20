using br.ufc.pargo.hpe.kinds;
using common.problem_size.Class;
using common.problem_size.Instance;

namespace adi.Add { 

public interface IAdd<I,C> : BaseIAdd<I,C>
		where I:IInstance<C>
		where C:IClass
{

} // end main interface 

} // end namespace 
