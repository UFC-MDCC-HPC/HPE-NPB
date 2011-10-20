using br.ufc.pargo.hpe.kinds;
using common.interactionpattern.wavefront.WavefrontWork;
using common.direction.Direction;

namespace common.interactionpattern.Wavefront { 

public interface IWavefront<W, DIR> : BaseIWavefront<W, DIR>
where W:IWavefrontWork<DIR>
where DIR:IDirection
{


} // end main interface 

} // end namespace 
