using br.ufc.pargo.hpe.kinds;
using common.problem_size.Class;
using common.problem_size.Instance;
using common.solve.SolvingMethod;

namespace adi.ADI_Solver3D { 

public interface IADI_Solver3D<MTH, CLASS, I> : BaseIADI_Solver3D<MTH, CLASS, I>
where MTH:ISolvingMethod
where CLASS:IClass
where I:IInstance<CLASS>
{


} // end main interface 

} // end namespace 
