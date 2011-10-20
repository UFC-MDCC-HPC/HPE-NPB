using br.ufc.pargo.hpe.kinds;
using common.solve.SolvingMethod;

namespace bt.solve.BinvRhs { 
	public interface IBinvRhs<MTH> : BaseIBinvRhs<MTH>
	where MTH:ISolvingMethod {
	   void setParameters(double[,,] lhs, 
	                      double[,,,,] r, 
	                      int l1, 
	                      int r1, int r2, int r3, int r4);
	}
}
