/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using ft.problem_size.Instance;
using common.problem_size.Class;
using ft.datapartition.BlocksInfo;

namespace ft.data.ProblemDefinition 
{ 
	public interface BaseIProblemDefinition<I, C> : IDataStructureKind 
		where I:IInstance<C>
		where C:IClass 
	{
		IBlocks<I,C> Blocks {get;}
	}
}
