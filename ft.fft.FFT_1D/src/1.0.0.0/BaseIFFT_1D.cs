/* AUTOMATICALLY GENERATE CODE */

using br.ufc.pargo.hpe.kinds;
using ft.datapartition.BlocksInfo;
using ft.data.ProblemDefinition;
using ft.problem_size.Instance;
using common.problem_size.Class;
using common.solve.SolvingMethod;

namespace ft.fft.FFT_1D { 

public interface BaseIFFT_1D<I, C, M> : IComputationKind 
where I:IInstance<C>
where C:IClass
where M:ISolvingMethod
{

	IBlocks<I,C> Blocks {get;}
	IProblemDefinition<I, C> Problem {get;}


} // end main interface 

} // end namespace 
