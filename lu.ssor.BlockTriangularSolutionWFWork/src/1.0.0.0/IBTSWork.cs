using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using lu.problem_size.Instance;
using common.problem_size.Class;
using common.interactionpattern.wavefront.WavefrontWork;

namespace lu.ssor.BlockTriangularSolutionWFWork { 

public interface IBTSWork<DIR, I, C> : BaseIBTSWork<DIR, I, C>, IWavefrontWork<DIR>
where DIR:IDirection
where I:IInstance<C>
where C:IClass
{
   void setParameters(int k, double omega);

} // end main interface 

} // end namespace 
