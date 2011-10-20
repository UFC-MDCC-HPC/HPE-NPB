using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.solve.SolvingMethod;
using common.axis.Axis;
using adi.Solver;

namespace sp.solve.Solver { 

public interface ISPSolver<I, C, MTH, DIR> : BaseISPSolver<I, C, MTH, DIR>, ISolver<I, C, MTH, DIR>
where I:IInstance_SP<C>
where C:IClass
where MTH:ISolvingMethod
where DIR:IAxis
{


} // end main interface 

} // end namespace 
