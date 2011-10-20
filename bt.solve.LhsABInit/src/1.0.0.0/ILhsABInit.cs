using br.ufc.pargo.hpe.kinds;

namespace bt.solve.LhsABInit { 
	public interface ILhsABInit : BaseILhsABInit {
	   void setParameters(double[,,] lhsa, double[,,] lhsb, int size);
	} // end main interface 
} // end namespace 
