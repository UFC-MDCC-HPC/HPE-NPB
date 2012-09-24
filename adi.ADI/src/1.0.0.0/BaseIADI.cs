/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Instance; 
using common.problem_size.Class;
using common.solve.SolvingMethod;
using environment.MPIDirect;

namespace adi.ADI { 

public interface BaseIADI<I, CLASS, MTH> : IComputationKind 
where MTH:ISolvingMethod
where CLASS:IClass
where I:IInstance<CLASS>
{

	IProblemDefinition<I, CLASS> Problem {get;}
	IMPIDirect Mpi {get;}
	ICell X {get;}
	ICell Y {get;}
	ICell Z {get;}
	ICells Cells {get;}


} // end main interface 

} // end namespace 
