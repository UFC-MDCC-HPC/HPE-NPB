using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.problem_size.Class;
using sp.problem_size.Instance_SP;
using sp.solve.SPMethod;
using sp.adi.SP_ADI;

namespace impl.sp.adi.SP_ADIImpl { 

public class ISP_ADIImpl<I,C,MTH> : BaseISP_ADIImpl<I,C,MTH>, ISP_ADI<I,C,MTH>
where MTH:ISPMethod
where C:IClass
where I:IInstance_SP<C>
{

public override int go() 
{ 
	int no_nodes = this.Size;
	if (no_nodes > 1) 
	{
	   Copy_faces.go();
	}
	Compute_rhs.go();
	Txinvr.go();
	X_solve.go();
	Y_solve.go();
	Z_solve.go();
	Add.go();
			
	return 0;
} // end activate method 

} 

}

