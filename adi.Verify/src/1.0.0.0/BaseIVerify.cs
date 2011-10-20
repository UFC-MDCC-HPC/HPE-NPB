/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.datapartition.MultiPartitionCells;
using common.problem_size.Instance;
using common.problem_size.Class;
using environment.MPIDirect;
using adi.data.ProblemDefinition;

namespace adi.Verify { 

public interface BaseIVerify<I, C> : IComputationKind
where C:IClass 
where I:IInstance<C>
{
	IMPIDirect Mpi {get;}
	IProblemDefinition<I,C> Problem {get;}
	ICell X {get;}
	ICell Y {get;}
	ICell Z {get;}
	ICells Cells_info {get;}
} // end main interface 

} // end namespace 
