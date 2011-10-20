/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using sp.problem_size.Instance_SP; 
using common.problem_size.Class;
using environment.MPIDirect;

namespace sp.ADI { 

public interface BaseIADI<C> : IComputationKind 
where C:IClass
{

	IProblemDefinition<IInstance_SP<C>, C> Problem {get;}
	IMPIDirect Mpi {get;}
	ICell X {get;}
	ICell Y {get;}
	ICell Z {get;}
	ICells Cells {get;}


} // end main interface 

} // end namespace 
