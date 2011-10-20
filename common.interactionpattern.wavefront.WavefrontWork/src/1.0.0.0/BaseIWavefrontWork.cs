/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.Buffer;
using common.direction.Direction;

namespace common.interactionpattern.wavefront.WavefrontWork { 

public interface BaseIWavefrontWork<DIR> : IComputationKind 
where DIR:IDirection
{

	IBuffer X_buffer {get;}
	IBuffer Y_buffer {get;}


} // end main interface 

} // end namespace 
