/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using ft.datapartition.BlocksInfo;
using ft.fft.Transpose;
using ft.problem_size.Instance_FT;
using common.problem_size.Class;
using common.axis.XY_Z_Axes;
using ft.data.ProblemDefinition;
using common.axis.X_YZ_Axis;
using common.axis.X_Y_Axis;
using common.axis.X_Z_Axis;
using ft.fft.Cffts;
using common.axis.YAxis;
using common.axis.ZAxis;
using common.axis.XAxis;
using ft.fft.Fft;

namespace impl.ft.fft.FftImpl { 
	public abstract class BaseIFftImpl<I, C>: Computation, BaseIFft<I, C>
	where I:IInstance_FT<C>
	where C:IClass{
	   
		#region data
			protected int layout_0D, layout_1D, layout_2D, layout_type;
			protected int[,] dims;
			override public void initialize(){
				layout_0D = Constants.layout_0D;
				layout_1D = Constants.layout_1D;
				layout_2D = Constants.layout_2D;
				layout_type = Blocks.layout_type;
				
				dims = Blocks.dims;
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
		
		private ITranspose<I, C, IXY_Z_Axes> transpose_xy_z = null;
		
		protected ITranspose<I, C, IXY_Z_Axes> Transpose_xy_z {
			get {
				if (this.transpose_xy_z == null)
					this.transpose_xy_z = (ITranspose<I, C, IXY_Z_Axes>) Services.getPort("transpose_xy_z");
				return this.transpose_xy_z;
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
		
		private ITranspose<I, C, IX_YZ_Axis> transpose_x_yz = null;
		
		protected ITranspose<I, C, IX_YZ_Axis> Transpose_x_yz {
			get {
				if (this.transpose_x_yz == null)
					this.transpose_x_yz = (ITranspose<I, C, IX_YZ_Axis>) Services.getPort("transpose_x_yz");
				return this.transpose_x_yz;
			}
		}
		
		private ITranspose<I, C, IX_Y_Axis> transpose_x_y = null;
		
		protected ITranspose<I, C, IX_Y_Axis> Transpose_x_y {
			get {
				if (this.transpose_x_y == null)
					this.transpose_x_y = (ITranspose<I, C, IX_Y_Axis>) Services.getPort("transpose_x_y");
				return this.transpose_x_y;
			}
		}
		
		private ITranspose<I, C, IX_Z_Axis> transpose_x_z = null;
		
		protected ITranspose<I, C, IX_Z_Axis> Transpose_x_z {
			get {
				if (this.transpose_x_z == null)
					this.transpose_x_z = (ITranspose<I, C, IX_Z_Axis>) Services.getPort("transpose_x_z");
				return this.transpose_x_z;
			}
		}
		
		private ICffts<I, C, IY> cffts2 = null;
		
		protected ICffts<I, C, IY> Cffts2 {
			get {
				if (this.cffts2 == null)
					this.cffts2 = (ICffts<I, C, IY>) Services.getPort("cffts2");
				return this.cffts2;
			}
		}
		
		private ICffts<I, C, IZ> cffts3 = null;
		
		protected ICffts<I, C, IZ> Cffts3 {
			get {
				if (this.cffts3 == null)
					this.cffts3 = (ICffts<I, C, IZ>) Services.getPort("cffts3");
				return this.cffts3;
			}
		}
		
		private ICffts<I, C, IX> cffts1 = null;
		
		protected ICffts<I, C, IX> Cffts1 {
			get {
				if (this.cffts1 == null)
					this.cffts1 = (ICffts<I, C, IX>) Services.getPort("cffts1");
				return this.cffts1;
			}
		}
		abstract public int go(); 
	}
}
