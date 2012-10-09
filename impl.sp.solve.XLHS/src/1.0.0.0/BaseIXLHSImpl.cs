/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.datapartition.MultiPartitionCells;
using adi.data.ProblemDefinition;
using common.problem_size.Class;
using sp.solve.LHS;
using sp.problem_size.Instance_SP;
using common.axis.XAxis;
using sp.solve.SPMethod;

namespace impl.sp.solve.XLHS { 

public abstract class BaseIXLHSImpl<I,C,DIR,MTH>: Computation, BaseILHS<I,C,DIR,MTH>
	where I:IInstance_SP<C>
	where C:IClass
	where DIR:IX
	where MTH:ISPMethod
{
		
#region data
		
protected int[,] start, end, cell_size, slice;
protected double[,,,,] lhs, rho_i, speed, us;
protected double c3c4, dttx2, c2dttx1, dttx1, con43, dx5, dx1,
               comz5, comz4, comz1, comz6, dx2, c1c5, dxmax;
protected int MAX_CELL_DIM;
protected double[] cv, rhon;
protected int ncells;

override public void initialize()		
{
	start = Cells.cell_start;
	end = Cells.cell_end;
	cell_size = Cells.cell_size;
	slice = Cells.cell_slice;
			
			
	MAX_CELL_DIM = Problem.MAX_CELL_DIM;
	
	lhs = Problem.Field_lhs;
	rho_i = Problem.Field_rho;
	speed = Problem.Field_speed;
	us = Problem.Field_us;
	ncells = Problem.NCells;
			
	c3c4 = Constants.c3c4;
	dttx2 = Constants.dttx2;
	c2dttx1 = Constants.c2dttx1;
	dttx1 = Constants.dttx1;
	con43 = Constants.con43;
	dx5 = Constants.dx5;
	dx1 = Constants.dx1;
	comz4 = Constants.comz4;
	comz1 = Constants.comz1;
	comz6 = Constants.comz6;
	comz5 = Constants.comz5;
	dx2 = Constants.dx2;
	c1c5 = Constants.c1c5;
	dxmax = Constants.dxmax;
			
    cv = new double[MAX_CELL_DIM + 4];     /* -2 */   // lhsx, lhsy, lhsz (def/use)
    rhon = new double[MAX_CELL_DIM + 4];   /* -2 */   // lhsx (local)
}
		
#endregion
		
private ICells cells = null;

public ICells Cells {
	get {
		if (this.cells == null)
		{
			this.cells = (ICells) Services.getPort("cells_info");
		}
		return this.cells;
	}
}

private IProblemDefinition<I,C> problem = null;

public IProblemDefinition<I,C> Problem {
	get {
		if (this.problem == null)
			this.problem = (IProblemDefinition<I,C>) Services.getPort("problem_data");
		return this.problem;
	}
}
		
private DIR axis = default(DIR);

protected DIR Axis {
	get {
		if (this.axis == null)
			this.axis = (DIR) Services.getPort("axis");
		return this.axis;
	}
}

private MTH method = default(MTH);

protected MTH Method {
	get {
		if (this.method == null)
			this.method = (MTH) Services.getPort("method");
		return this.method;
	}
}
		

public BaseIXLHSImpl() { 

} 

abstract public int go(); 


}

}
