/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using common.topology.Ring;
using lu.data.ProblemDefinition;
using lu.problem_size.Instance;
using common.problem_size.Class;
using lu.datapartition.BlocksInfo;
using environment.MPIDirect;
using lu.exchange.ExchangePattern;

namespace lu.exchange.Exchange2D { 

public interface BaseIExchange2D<I, C, E> : ISynchronizerKind 
where I:IInstance<C>
where C:IClass
where E:IExchangePattern
{

	ICell X {get;}
	ICell Y {get;}
	IProblemDefinition<I, C> Problem {get;}
	IBlocksInfo<I, C> Blocks {get;}
	IMPIDirect Mpi {get;}


} // end main interface 

} // end namespace 
