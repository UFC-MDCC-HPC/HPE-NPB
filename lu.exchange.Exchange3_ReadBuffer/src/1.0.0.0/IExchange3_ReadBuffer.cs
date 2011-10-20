using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using common.axis.Axis;
using common.interactionpattern.shift1D.BufferReader;

namespace lu.exchange.Exchange3_ReadBuffer { 

public interface IExchange3_ReadBuffer<DIR, O> : BaseIExchange3_ReadBuffer<DIR, O>, IBufferReader<DIR>
where DIR:IDirection
where O:IAxis
{
   void setParameters(double[,,,] g, int nx, int ny, int nz);

} // end main interface 

} // end namespace 
