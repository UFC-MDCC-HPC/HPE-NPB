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
using common.axis.X_YZ_Axis;
using ft.fft.TransposeGlobal;
using MPI;

namespace impl.ft.fft.TransposeX_YZGlobal { 
	public abstract class BaseITransposeX_YZGlobal<I, C, DIR>: Computation, BaseITransposeGlobal<I, C, DIR>
	where I:IInstance_FT<C>
	where C:IClass
	where DIR:IX_YZ_Axis
	{
	   
		#region data
		protected int me1,me2;
		
		protected int np;
		
		override public void initialize()
		{
				me1 = Blocks.me1;
				me2 = Blocks.me2;
			    np  = Blocks.np;
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
