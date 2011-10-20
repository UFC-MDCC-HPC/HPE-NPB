/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.Buffer;
using common.direction.Direction;
using common.axis.Axis;
using common.direction.Forward;
using common.axis.YAxis;
using lu.exchange.Exchange3_ReadBuffer;

namespace impl.lu.exchange.Exchange3_YForwardReadBufferImpl { 

public abstract class BaseIExchange3_YForwardReadBufferImpl<DIR, O>: Computation, BaseIExchange3_ReadBuffer<DIR, O>
where DIR:IForwardDirection
where O:IY
{

private IBuffer buffer = null;

public IBuffer Buffer {
	get {
		if (this.buffer == null)
			this.buffer = (IBuffer) Services.getPort("buffer");
		return this.buffer;
	}
}

private DIR direction = default(DIR);

protected DIR Direction {
	get {
		if (this.direction == null)
			this.direction = (DIR) Services.getPort("direction");
		return this.direction;
	}
}

private O axis = default(O);

protected O Axis {
	get {
		if (this.axis == null)
			this.axis = (O) Services.getPort("axis");
		return this.axis;
	}
}


abstract public int go(); 


}

}
