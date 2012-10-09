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
using common.axis.YAxis;
using sp.solve.SPMethod;

namespace impl.sp.solve.YLHS { 

public abstract class BaseIYLHSImpl<I,C,DIR,MTH>: Computation, BaseILHS<I,C,DIR,MTH>
	where I:IInstance_SP<C>
	where C:IClass
	where DIR:IY
	where MTH:ISPMethod
{
		
		
#region data
		
protected int[,] start, end, cell_size, slice;
protected double[,,,,] lhs, rho_i, speed, vs;
protected double c3c4, dtty2, c2dtty1, dtty1, con43, dy5, dy1,
               comz5, comz4, comz1, comz6, dy3, c1c5, dymax;
protected int ncells;
		
protected double[] cv, rhoq;

protected int MAX_CELL_DIM;
		
override public void initialize()
{	
	start = Cells.cell_start;
	end = Cells.cell_end;
	cell_size = Cells.cell_size;
	slice = Cells.cell_slice;
	
	MAX_CELL_DIM = Problem.MAX_CELL_DIM;
	
    cv = new double[MAX_CELL_DIM + 4];     /* -2 */   // lhsx, lhsy, lhsz (def/use)
    rhoq = new double[MAX_CELL_DIM + 4];   /* -2 */   // lhsx (local)
	
			
	lhs = Problem.Field_lhs;
	rho_i = Problem.Field_rho;
	speed = Problem.Field_speed;
	vs = Problem.Field_vs;
	ncells = Problem.NCells;
			
	c3c4 = Constants.c3c4;
	dtty2 = Constants.dtty2;
	c2dtty1 = Constants.c2dtty1;
	dtty1 = Constants.dtty1;
	con43 = Constants.con43;
	dy5 = Constants.dy5;
	dy1 = Constants.dy1;
	comz4 = Constants.comz4;
	comz1 = Constants.comz1;
	comz6 = Constants.comz6;
	comz5 = Constants.comz5;
	dy3 = Constants.dy3;
	c1c5 = Constants.c1c5;
	dymax = Constants.dymax;
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
		{
			this.problem = (IProblemDefinition<I,C>) Services.getPort("problem_data");
		}
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
		
abstract public int go(); 


}

}
