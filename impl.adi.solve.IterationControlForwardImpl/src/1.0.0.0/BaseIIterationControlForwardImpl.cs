/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;
using common.direction.Forward;
using adi.solve.IterationControl;

namespace impl.adi.solve.IterationControlForwardImpl { 

public abstract class BaseIIterationControlForwardImpl<D>: Computation, BaseIIterationControl<D>
where D:IForwardDirection
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
