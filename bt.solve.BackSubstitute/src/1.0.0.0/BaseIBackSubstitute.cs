/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.solve.SolvingMethod;
using common.axis.Axis;

namespace bt.solve.BackSubstitute { 

public interface BaseIBackSubstitute<I, C, DIR, MTH> : IComputationKind 
where I:IInstance<C>
where C:IClass
where DIR:IAxis
where MTH:ISolvingMethod
{

	ICells Cells {get;}
	IProblemDefinition<I, C> Problem {get;}


} // end main interface 

} // end namespace 
