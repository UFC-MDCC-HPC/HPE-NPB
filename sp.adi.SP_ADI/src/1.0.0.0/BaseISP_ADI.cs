/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using sp.solve.SPMethod;
using common.problem_size.Class;
using sp.problem_size.Instance_SP;
using adi.ADI;

namespace sp.adi.SP_ADI { 

public interface BaseISP_ADI<I, C, MTH> : BaseIADI<I, C, MTH>, IComputationKind 
where MTH:ISPMethod
where C:IClass
where I:IInstance_SP<C>
{



} // end main interface 

} // end namespace 
