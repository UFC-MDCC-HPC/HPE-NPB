using br.ufc.pargo.hpe.kinds;
using common.problem_size.Class;
using common.solve.SolvingMethod;
using common.problem_size.Instance; 

namespace adi.ADI { 

public interface IADI<I, CLASS, MTH> : BaseIADI<I, CLASS, MTH>
where MTH:ISolvingMethod
where CLASS:IClass
where I:IInstance<CLASS>
{


} // end main interface 

} // end namespace 
