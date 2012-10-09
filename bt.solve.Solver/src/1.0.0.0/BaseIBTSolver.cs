/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.Buffer;
using environment.MPIDirect;
using common.axis.Axis;
using common.solve.SolvingMethod;
using adi.Solver;

namespace bt.solve.Solver { 

public interface BaseIBTSolver<I, C, MTH, DIR> : BaseISolver<I, C, MTH, DIR>,IComputationKind 
where I:IInstance<C>
where C:IClass
where DIR:IAxis
where MTH:ISolvingMethod
{

	ICell Topology {get;}
	ICells Cells_info {get;}
	IProblemDefinition<I, C> Problem_data {get;}
	IMPIDirect Mpi {get;}


} // end main interface 

} // end namespace 
