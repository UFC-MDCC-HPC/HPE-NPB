/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Class;
using common.problem_size.Instance;


namespace adi.ComputeRHS { 

public interface BaseIComputeRHS<I,C> : IComputationKind 
		where I:IInstance<C>
		where C:IClass
{

	ICells Cells {get;}
	IProblemDefinition<I,C> Problem {get;}


} // end main interface 

} // end namespace 
