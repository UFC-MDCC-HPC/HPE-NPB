using br.ufc.pargo.hpe.kinds;
using common.solve.SolvingMethod;
using common.axis.Axis;
using common.problem_size.Class;
using common.problem_size.Instance;
using adi.Solver;

namespace bt.solve.Solver { 

public interface IBTSolver<I, C, MTH, DIR> : BaseIBTSolver<I, C, MTH, DIR>, ISolver<I, C, MTH, DIR>
where MTH:ISolvingMethod
where DIR:IAxis
where C:IClass
where I:IInstance<C>
{


} // end main interface 

} // end namespace 
