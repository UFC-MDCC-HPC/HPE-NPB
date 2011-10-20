using br.ufc.pargo.hpe.kinds;
using common.problem_size.Class;
using common.solve.SolvingMethod;
using common.problem_size.Instance; 

namespace adi.ADI { 

public interface IADI<MTH, CLASS, I> : BaseIADI<MTH, CLASS, I>
where MTH:ISolvingMethod
where CLASS:IClass
where I:IInstance<CLASS>
{


} // end main interface 

} // end namespace 
