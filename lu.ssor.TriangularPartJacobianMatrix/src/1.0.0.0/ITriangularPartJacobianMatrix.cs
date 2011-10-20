using br.ufc.pargo.hpe.kinds;
using lu.problem_size.Instance;
using common.problem_size.Class;
using common.direction.Direction;

namespace lu.ssor.TriangularPartJacobianMatrix { 
	public interface ITriangularPartJacobianMatrix<I, C, DIS> : BaseITriangularPartJacobianMatrix<I, C, DIS>
	where I:IInstance<C>
	where C:IClass
	where DIS:IDirection{
	   void setParameters(int k);
	}
}
