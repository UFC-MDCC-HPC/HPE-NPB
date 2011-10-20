/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using ft.datapartition.BlocksInfo;
using ft.data.ProblemDefinition;
using ft.problem_size.Instance_FT;
using common.problem_size.Class;
using ft.fft.FFT_1D;
using ft.fft.method.SwarztrauberMethod;

namespace impl.ft.fft.SwarztrauberFFT1DImpl { 
	public abstract class BaseISwarztrauberFFT1DImpl<I, C, M>: Computation, BaseIFFT_1D<I, C, M>
	where I:IInstance_FT<C>
	where C:IClass
	where M:ISwarztrauberMethod
	{
	   
		#region data
			protected int REAL, IMAG, fftblock, fftblockpad;
			protected double[,] u;
			override public void initialize(){
				REAL = Constants.REAL;
				IMAG = Constants.IMAG;
			    fftblock = Blocks.fftblock;
                fftblockpad = Blocks.fftblockpad;
				u           = Problem.u;
			}
		#endregion
	
		private IBlocks<I,C> blocks = null;
		
		public IBlocks<I,C> Blocks {
			get {
				if (this.blocks == null)
					this.blocks = (IBlocks<I,C>) Services.getPort("blocks_info");
				return this.blocks;
			}
		}
		
		private IProblemDefinition<I, C> problem = null;
		
		public IProblemDefinition<I, C> Problem {
			get {
				if (this.problem == null)
					this.problem = (IProblemDefinition<I, C>) Services.getPort("problem_data");
				return this.problem;
			}
		}
		abstract public int go(); 
	}
}
