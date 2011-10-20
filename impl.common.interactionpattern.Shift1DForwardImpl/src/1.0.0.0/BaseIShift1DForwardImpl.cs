/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.interactionpattern.shift1D.BufferWriter;
using common.direction.Forward;
using common.topology.Ring;
using common.Buffer;
using common.interactionpattern.shift1D.BufferReader;
using common.direction.Direction;
using common.interactionpattern.Shift1D;

namespace impl.common.interactionpattern.Shift1DForwardImpl { 

public abstract class BaseIShift1DForwardImpl<BR, DIR, BW>: Synchronizer, BaseIShift1D<DIR, BR, BW>
where BR:IBufferReader<DIR>
where DIR:IForwardDirection
where BW:IBufferWriter<DIR>
{

private BW buffer_writer = default(BW);

public BW Buffer_writer {
	get {
		if (this.buffer_writer == null)
			this.buffer_writer = (BW) Services.getPort("buffer_writer");
		return this.buffer_writer;
	}
}

private ICell axis = null;

public ICell Axis {
	get {
		if (this.axis == null)
			this.axis = (ICell) Services.getPort("axis");
		return this.axis;
	}
}

private IBuffer input_buffer = null;

protected IBuffer Input_buffer {
	get {
		if (this.input_buffer == null)
			this.input_buffer = (IBuffer) Services.getPort("input_buffer");
		return this.input_buffer;
	}
}

private IBuffer output_buffer = null;

protected IBuffer Output_buffer {
	get {
		if (this.output_buffer == null)
			this.output_buffer = (IBuffer) Services.getPort("output_buffer");
		return this.output_buffer;
	}
}

private BR buffer_reader = default(BR);

public BR Buffer_reader {
	get {
		if (this.buffer_reader == null)
			this.buffer_reader = (BR) Services.getPort("buffer_reader");
		return this.buffer_reader;
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


abstract public int go(); 


}

}
