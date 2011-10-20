using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.problem_size.Class;
using adi.ADI; 
using common.solve.SolvingMethod;
using common.problem_size.Instance;

namespace impl.adi.ADIImpl { 

public class IADIImpl<MTH, CLASS, I> : BaseIADIImpl<MTH, CLASS, I>, IADI<MTH, CLASS, I>
	where MTH:ISolvingMethod
	where CLASS:IClass
	where I:IInstance<CLASS>
{

public override int go() 
{
	int no_nodes = Ranks.Length;
	if (no_nodes > 1) 
	{
	   Copy_faces.go();
	}
	Compute_rhs.go();
	
	X_solve.go();
	Y_solve.go();
	Z_solve.go();
	Add.go();
	
	return 0;
} // end activate method 

}

}
