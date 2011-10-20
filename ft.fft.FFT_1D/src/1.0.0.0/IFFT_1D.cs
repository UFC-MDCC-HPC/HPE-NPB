using br.ufc.pargo.hpe.kinds;
using ft.problem_size.Instance;
using common.problem_size.Class;
using common.solve.SolvingMethod;

namespace ft.fft.FFT_1D { 
	public interface IFFT_1D<I, C, M> : BaseIFFT_1D<I, C, M>
	where I:IInstance<C>
	where C:IClass
	where M:ISolvingMethod
	{
	   void setParameters(int dir, int m, int n, double[, , ,] y);
	}
}
