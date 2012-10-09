/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using adi.data.ProblemDefinition;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.Buffer;
using common.direction.Direction;
using common.axis.Axis;
using common.datapartition.MultiPartitionCells;

namespace sp.solve.shift.ReadBuffer { 

public interface BaseIReadBuffer<I, C, D, DIR> : IComputationKind 
where I:IInstance<C>
where C:IClass
where D:IDirection
where DIR:IAxis
{

	IProblemDefinition<I, C> Problem {get;}
	IBuffer Buffer {get;}
	ICells Cells {get;}


} // end main interface 

} // end namespace 
