/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using ft.datapartition.BlocksInfo;
using ft.data.ProblemDefinition;
using ft.problem_size.Instance_FT;
using common.problem_size.Class;
using common.axis.Axis;
using common.axis.XY_Z_Axes;
using ft.fft.TransposeFinish;

namespace impl.ft.fft.TransposeXY_ZFinish { 
	public abstract class BaseITransposeXY_ZFinish<I, C, DIR>: Computation, BaseITransposeFinish<I, C, DIR>
	where I:IInstance_FT<C>
	where C:IClass
	where DIR:IXY_Z_Axes
	{
	
		#region data
			protected int np1, np2, size1, size2;
			protected int REAL;
			protected int IMAG;
			
			override public void initialize(){
				np1 = Blocks.np1;
				np2 = Blocks.np2;
				size1 = Blocks.size1;
				size2 = Blocks.size2;
				REAL = Constants.REAL;
				IMAG = Constants.IMAG;
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
