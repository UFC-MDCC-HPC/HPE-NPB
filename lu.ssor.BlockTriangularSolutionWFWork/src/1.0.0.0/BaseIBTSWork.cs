/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using lu.data.ProblemDefinition;
using lu.problem_size.Instance;
using common.problem_size.Class;
using lu.datapartition.BlocksInfo;
using common.direction.Direction;
using common.interactionpattern.wavefront.WavefrontWork;

namespace lu.ssor.BlockTriangularSolutionWFWork { 

public interface BaseIBTSWork<DIR, I, C> : BaseIWavefrontWork<DIR>, IComputationKind 
where DIR:IDirection
where I:IInstance<C>
where C:IClass
{

	IProblemDefinition<I, C> Problem {get;}
	IBlocksInfo<I, C> Blocks {get;}


} // end main interface 

} // end namespace 
