using br.ufc.pargo.hpe.kinds;
using common.axis.Axis;

namespace lu.exchange.Exchange3 { 

public interface IExchange3<O> : BaseIExchange3<O>
where O:IAxis
{
   void setParameters(double[,,,] g, int nx, int ny, int nz);

} // end main interface 

} // end namespace 
