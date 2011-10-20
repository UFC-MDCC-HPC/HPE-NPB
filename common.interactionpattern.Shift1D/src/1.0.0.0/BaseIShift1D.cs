/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using common.interactionpattern.shift1D.BufferReader;
using common.interactionpattern.shift1D.BufferWriter;

namespace common.interactionpattern.Shift1D { 

public interface BaseIShift1D<DIR, BR, BW> : ISynchronizerKind 
where DIR:IDirection
where BR:IBufferReader<DIR>
where BW:IBufferWriter<DIR>
{



} // end main interface 

} // end namespace 
