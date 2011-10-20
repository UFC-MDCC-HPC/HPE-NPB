using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using common.interactionpattern.shift1D.BufferReader;
using common.interactionpattern.shift1D.BufferWriter;

namespace common.interactionpattern.Shift1D { 

public interface IShift1D<DIR, BR, BW> : BaseIShift1D<DIR, BR, BW>
where DIR:IDirection
where BR:IBufferReader<DIR>
where BW:IBufferWriter<DIR>
{
	void setTag(int tag);

} // end main interface 

} // end namespace 
