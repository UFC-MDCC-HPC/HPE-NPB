/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.Buffer;
using common.direction.Direction;
using common.axis.Axis;
using common.direction.Backward;
using common.axis.YAxis;
using lu.exchange.Exchange3_WriteBuffer;

namespace impl.lu.exchange.Exchange3_YBackwardWriteBufferImpl { 

public abstract class BaseIExchange3_YBackwardWriteBufferImpl<DIR, O>: Computation, BaseIExchange3_WriteBuffer<DIR, O>
where DIR:IBackwardDirection
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
