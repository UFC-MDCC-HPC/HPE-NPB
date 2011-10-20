/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.datapartition.MultiPartitionCells;
using common.problem_size.Instance;
using common.problem_size.Class;

namespace common.datapartition.MultiPartition { 

public interface BaseIMultiPartition<I,C> : IEnvironmentKind 
where I:IInstance<C>
where C:IClass
{

	ICell Z {get;}
	ICell Y {get;}
	ICell X {get;}
	ICells Cells {get;}


} // end main interface 

} // end namespace 
