/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using lu.data.ProblemDefinition;
using lu.problem_size.Instance;
using common.problem_size.Class;
using lu.datapartition.BlocksInfo;
using common.direction.Direction;

namespace lu.ssor.TriangularPartJacobianMatrix { 
	public interface BaseITriangularPartJacobianMatrix<I, C, DIS> : IComputationKind 
	where I:IInstance<C>
	where C:IClass
	where DIS:IDirection{
		IProblemDefinition<I, C> Problem {get;}
		IBlocksInfo<I, C> Blocks {get;}
	}
}
