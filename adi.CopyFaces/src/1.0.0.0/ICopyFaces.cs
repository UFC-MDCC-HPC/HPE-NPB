using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;

namespace adi.CopyFaces { 

public interface ICopyFaces<I,C> : BaseICopyFaces<I,C>
		where I:IInstance<C>
		where C:IClass
{
		bool is_multiple();
} // end main interface 

} // end namespace 
