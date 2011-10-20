using br.ufc.pargo.hpe.kinds;
using common.solve.SolvingMethod;

namespace bt.solve.MatMulSub { 
	public interface IMatMulSub<MTH> : BaseIMatMulSub<MTH>
	where MTH:ISolvingMethod {
	   void setParameters(double[,,] ablock, 
	                      double[,,,,,] bblock, 
	                      double[,,] cblock, 
	                      int a1, 
	                      int b1, int b2, int b3, int b4, 
	                      int c1);
	}
}
