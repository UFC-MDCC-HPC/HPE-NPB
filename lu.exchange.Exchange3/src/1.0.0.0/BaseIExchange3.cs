/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.interactionpattern.shift1D.BufferReader;
using common.direction.Backward;
using common.interactionpattern.shift1D.BufferWriter;
using common.direction.Forward;
using environment.MPIDirect;
using common.axis.Axis;

namespace lu.exchange.Exchange3 { 

public interface BaseIExchange3<O> : IComputationKind 
where O:IAxis
{

	ICell Cell {get;}
	IMPIDirect Mpi {get;}


} // end main interface 

} // end namespace 
