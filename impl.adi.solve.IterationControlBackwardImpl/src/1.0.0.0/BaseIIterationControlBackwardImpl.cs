/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using common.direction.Backward;
using adi.solve.IterationControl;

namespace impl.adi.solve.IterationControlBackwardImpl { 

public abstract class BaseIIterationControlBackwardImpl<D>: Computation, BaseIIterationControl<D>
where D:IBackwardDirection
{

private D direction = default(D);

protected D Direction {
	get {
		if (this.direction == null)
			this.direction = (D) Services.getPort("direction");
		return this.direction;
	}
}


abstract public int go(); 


}

}
