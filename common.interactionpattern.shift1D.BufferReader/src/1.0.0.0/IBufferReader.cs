using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;

namespace common.interactionpattern.shift1D.BufferReader { 

public interface IBufferReader<DIR> : BaseIBufferReader<DIR>
where DIR:IDirection
{


} // end main interface 

} // end namespace 
