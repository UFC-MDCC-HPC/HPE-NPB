/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using common.interactionpattern.shift1D.BufferReader;
using common.direction.Backward;
using common.interactionpattern.shift1D.BufferWriter;
using common.direction.Forward;
using environment.MPIDirect;
using common.interactionpattern.Shift1D;
using lu.exchange.Exchange3_ReadBuffer;
using lu.exchange.Exchange3_WriteBuffer;
using common.axis.Axis;
using common.axis.YAxis;
using lu.exchange.Exchange3;

namespace impl.lu.exchange.Exchange3YImpl { 

public abstract class BaseIExchange3YImpl<O>: Computation, BaseIExchange3<O>
where O:IY
{

private ICell cell = null;

public ICell Cell {
	get {
		if (this.cell == null)
			this.cell = (ICell) Services.getPort("topology");
		return this.cell;
	}
}

private IExchange3_ReadBuffer<IBackwardDirection, O> buffer_reader_backward = null;

protected IExchange3_ReadBuffer<IBackwardDirection, O> Buffer_reader_backward {
	get {
		if (this.buffer_reader_backward == null)
			this.buffer_reader_backward = (IExchange3_ReadBuffer<IBackwardDirection, O>) Services.getPort("buffer_reader_backward");
		return this.buffer_reader_backward;
	}
}

private IExchange3_WriteBuffer<IForwardDirection, O> buffer_writer_forward = null;

protected IExchange3_WriteBuffer<IForwardDirection, O> Buffer_writer_forward {
	get {
		if (this.buffer_writer_forward == null)
			this.buffer_writer_forward = (IExchange3_WriteBuffer<IForwardDirection, O>) Services.getPort("buffer_writer_forward");
		return this.buffer_writer_forward;
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

private IShift1D<IBackwardDirection, IExchange3_ReadBuffer<IBackwardDirection, O>, IExchange3_WriteBuffer<IBackwardDirection, O>> shift_backward = null;

protected IShift1D<IBackwardDirection, IExchange3_ReadBuffer<IBackwardDirection, O>, IExchange3_WriteBuffer<IBackwardDirection, O>> Shift_backward {
	get {
		if (this.shift_backward == null)
			this.shift_backward = (IShift1D<IBackwardDirection, IExchange3_ReadBuffer<IBackwardDirection, O>, IExchange3_WriteBuffer<IBackwardDirection, O>>) Services.getPort("shift_backward");
		return this.shift_backward;
	}
}

private IExchange3_WriteBuffer<IBackwardDirection, O> buffer_writer_backward = null;

protected IExchange3_WriteBuffer<IBackwardDirection, O> Buffer_writer_backward {
	get {
		if (this.buffer_writer_backward == null)
			this.buffer_writer_backward = (IExchange3_WriteBuffer<IBackwardDirection, O>) Services.getPort("buffer_writer_backward");
		return this.buffer_writer_backward;
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

private IExchange3_ReadBuffer<IForwardDirection, O> buffer_reader_forward = null;

protected IExchange3_ReadBuffer<IForwardDirection, O> Buffer_reader_forward {
	get {
		if (this.buffer_reader_forward == null)
			this.buffer_reader_forward = (IExchange3_ReadBuffer<IForwardDirection, O>) Services.getPort("buffer_reader_forward");
		return this.buffer_reader_forward;
	}
}

private IShift1D<IForwardDirection, IExchange3_ReadBuffer<IForwardDirection, O>, IExchange3_WriteBuffer<IForwardDirection, O>> shift_forward = null;

protected IShift1D<IForwardDirection, IExchange3_ReadBuffer<IForwardDirection, O>, IExchange3_WriteBuffer<IForwardDirection, O>> Shift_forward {
	get {
		if (this.shift_forward == null)
			this.shift_forward = (IShift1D<IForwardDirection, IExchange3_ReadBuffer<IForwardDirection, O>, IExchange3_WriteBuffer<IForwardDirection, O>>) Services.getPort("shift_forward");
		return this.shift_forward;
	}
}


abstract public int go(); 


}

}
