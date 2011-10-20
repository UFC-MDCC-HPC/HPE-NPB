/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using common.topology.Ring;
using adi.data.ProblemDefinition;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.Buffer;
using environment.MPIDirect;
using common.solve.SolvingMethod;
using common.axis.Axis;
using adi.Solver;

namespace sp.solve.Solver { 

public interface BaseISPSolver<I, C, MTH, DIR> : BaseISolver<I, C, MTH, DIR>, IComputationKind 
where I:IInstance_SP<C>
where C:IClass
where MTH:ISolvingMethod
where DIR:IAxis
{

	ICells Cells {get;}
	ICell Cell {get;}
	IProblemDefinition<I, C> Problem {get;}
	IMPIDirect Mpi {get;}


} // end main interface 

} // end namespace 
