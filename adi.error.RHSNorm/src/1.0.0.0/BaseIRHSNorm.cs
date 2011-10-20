/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.datapartition.MultiPartitionCells;
using environment.MPIDirect;

namespace adi.error.RHSNorm { 

public interface BaseIRHSNorm<I,C> : IComputationKind 
where I:IInstance<C>
where C:IClass
{
		
	ICells Cells {get;}
	IProblemDefinition<I,C> Problem {get;}
	IMPIDirect Mpi {get;}



} // end main interface 

} // end namespace 
