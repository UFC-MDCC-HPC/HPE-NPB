/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.Buffer;
using common.interactionpattern.wavefront.WavefrontWork;
using common.direction.Forward;
using environment.MPIDirect;
using common.direction.Direction;
using common.interactionpattern.Wavefront;
using MPI;

namespace impl.common.interactionpattern.WavefrontForwardImpl { 

public abstract class BaseIWavefrontForwardImpl<W, DIR>: Computation, BaseIWavefront<W, DIR>
where W:IWavefrontWork<DIR>
where DIR:IForwardDirection
{

protected Intracommunicator comm;			
	
protected int north, west, south, east;	
		
override public void initialize()
{
   comm = this.WorldComm; //Mpi.localComm(this);
}		
		
override public void post_initialize()
{
   east = X.successor;
   west = X.predecessor;
   south = Y.successor;
   north = Y.predecessor;   
}		
		
private ICell x = null;

public ICell X {
	get {
		if (this.x == null)
			this.x = (ICell) Services.getPort("x");
		return this.x;
	}
}

private ICell y = null;

public ICell Y {
	get {
		if (this.y == null)
			this.y = (ICell) Services.getPort("y");
		return this.y;
	}
}

private IBuffer x_buffer = null;

protected IBuffer X_buffer {
	get {
		if (this.x_buffer == null)
			this.x_buffer = (IBuffer) Services.getPort("x_buffer");
		return this.x_buffer;
	}
}

private IBuffer y_buffer = null;

protected IBuffer Y_buffer {
	get {
		if (this.y_buffer == null)
			this.y_buffer = (IBuffer) Services.getPort("y_buffer");
		return this.y_buffer;
	}
}

private W wavefront_work = default(W);

public W Wavefront_work {
	get {
		if (this.wavefront_work == null)
			this.wavefront_work = (W) Services.getPort("work");
		return this.wavefront_work;
	}
}

private IMPIDirect mpi = null;

public IMPIDirect Mpi {
	get {
		if (this.mpi == null)
			this.mpi = (IMPIDirect) Services.getPort("mpi");
		return this.mpi;
	}
}

private IDirection direction = null;

protected IDirection Direction {
	get {
		if (this.direction == null)
			this.direction = (IDirection) Services.getPort("direction");
		return this.direction;
	}
}


abstract public int go(); 


}

}
