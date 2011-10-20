/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.Buffer;
using common.interactionpattern.wavefront.WavefrontWork;
using common.direction.Direction;
using environment.MPIDirect;

namespace common.interactionpattern.Wavefront { 

public interface BaseIWavefront<W, DIR> : IComputationKind 
where W:IWavefrontWork<DIR>
where DIR:IDirection
{
	ICell X {get;}
	ICell Y {get;}
	W Wavefront_work {get;}
	IMPIDirect Mpi {get;}


} // end main interface 

} // end namespace 
