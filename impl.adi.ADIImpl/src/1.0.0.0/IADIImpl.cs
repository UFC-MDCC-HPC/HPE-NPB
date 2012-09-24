using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using adi.ADI; 
using common.problem_size.Class;
using common.solve.SolvingMethod;
using common.problem_size.Instance;

namespace impl.adi.ADIImpl 
{ 
	public class IADIImpl<I, CLASS, MTH> : BaseIADIImpl<I, CLASS, MTH>, IADI<I, CLASS, MTH>
		where MTH:ISolvingMethod
		where CLASS:IClass
		where I:IInstance<CLASS>
	{	
		public override int go() 
		{
			int no_nodes = Size;
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
		} 
	
	}

}
