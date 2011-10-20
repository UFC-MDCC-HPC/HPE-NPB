/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.topology.Mesh3D;
using common.topology.Ring;
using common.datapartition.MultiPartitionCells;
using common.datapartition.MultiPartition;
using common.problem_size.Class;
using common.problem_size.Instance;


namespace impl.common.datapartition.MultiPartitionImpl { 

public abstract class BaseIMultiPartitionImpl<I, C>: br.ufc.pargo.hpe.kinds.Environment, BaseIMultiPartition<I, C>
where I:IInstance<C>
where C:IClass
{
		
		
private I instance = default(I);

protected I Instance {
	get {
		if (instance==null) 
			instance = (I) Services.getPort("instance");
		return instance;
	}
}		
		
private ICell3D cell = null;

protected ICell3D Cell {
	get {
		if (this.cell == null)
			this.cell = (ICell3D) Services.getPort("topology");
		return this.cell;
	}
}

private ICell z = null;

public ICell Z {
	get {
		if (this.z == null)
			this.z = (ICell) Services.getPort("z");
		return this.z;
	}
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

private ICells cells = null;

public ICells Cells {
	get {
		if (this.cells == null)
			this.cells = (ICells) Services.getPort("cells_info");
		return this.cells;
	}
}


}

}
