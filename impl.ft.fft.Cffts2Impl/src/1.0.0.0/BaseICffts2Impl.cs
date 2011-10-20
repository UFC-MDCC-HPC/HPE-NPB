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
using common.axis.Axis;
using common.axis.YAxis;
using ft.fft.Cffts;
using common.solve.SolvingMethod;


namespace impl.ft.fft.Cffts2Impl { 
	public abstract class BaseICffts2Impl<I, C, DIR>: Computation, BaseICffts<I, C, DIR>
	where I:IInstance_FT<C>
	where C:IClass
	where DIR:IY{
	   
		#region data
			protected int size1, size2, REAL, IMAG, fftblock, fftblockpad;
			override public void initialize(){
				REAL = Constants.REAL;
				IMAG = Constants.IMAG;
			    fftblock = Blocks.fftblock;
                fftblockpad = Blocks.fftblockpad;
                size1 = Blocks.size1;
                size2 = Blocks.size2;
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
		
		private IFFT_1D<I, C, ISolvingMethod> cfftz = null;
		
		protected IFFT_1D<I, C, ISolvingMethod> Cfftz {
			get {
				if (this.cfftz == null)
					this.cfftz = (IFFT_1D<I, C, ISolvingMethod>) Services.getPort("cfftz");
				return this.cfftz;
			}
		}
		
		private DIR axis = default(DIR);
		
		protected DIR Axis {
			get {
				if (this.axis == null)
					this.axis = (DIR) Services.getPort("axis");
				return this.axis;
			}
		}
		abstract public int go(); 
	}
}
