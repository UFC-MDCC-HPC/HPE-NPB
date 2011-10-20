using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;

namespace common.interactionpattern.wavefront.WavefrontWork { 

public interface IWavefrontWork<DIR> : BaseIWavefrontWork<DIR>
where DIR:IDirection
{


} // end main interface 

} // end namespace 
