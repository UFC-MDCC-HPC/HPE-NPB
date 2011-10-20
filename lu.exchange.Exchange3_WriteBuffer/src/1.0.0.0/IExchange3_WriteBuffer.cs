using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using common.axis.Axis;
using common.interactionpattern.shift1D.BufferWriter;

namespace lu.exchange.Exchange3_WriteBuffer { 

public interface IExchange3_WriteBuffer<DIR, O> : BaseIExchange3_WriteBuffer<DIR, O>, IBufferWriter<DIR>
where DIR:IDirection
where O:IAxis
{
   void setParameters(double[,,,] g, int nx, int ny, int nz);

} // end main interface 

} // end namespace 
