/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using adi.data.ProblemDefinition;
using bt.problem_size.Instance_BT;
using common.problem_size.Class;
using common.Buffer;
using bt.solve.BackSubstitute;
using common.axis.Axis;
using bt.solve.BTMethod;
using common.interactionpattern.Shift;
using common.direction.Forward;
using bt.solve.UnpackSolveInfo;
using bt.solve.SolverCell;
using common.direction.Backward;
using environment.MPIDirect;
using bt.solve.PackBackSubInfo;
using common.datapartition.MultiPartitionCells;
using common.topology.Ring;
using bt.solve.PackSolveInfo;
using bt.solve.UnpackBackSubInfo;
using bt.solve.Solver;
using adi.solve.IterationControl;

namespace impl.bt.solve.SolverImpl { 
	public abstract class BaseISolverImpl<I, C, MTH, DIR>: Computation, BaseIBTSolver<I, C, MTH, DIR>
	where I:IInstance_BT<C>
	where C:IClass
	where DIR:IAxis
	where MTH:IBTMethod {
		#region data
			protected int KMAX, JMAX, IMAX, ncells, MAX_CELL_DIM, maxcells;
			protected int[,] slice, cell_coord;
			protected double[, , , , ,] lhsc;
		    protected double[, , ,] backsub_info;
			
			override public void initialize(){
				KMAX = Problem_data.KMAX;
				JMAX = Problem_data.JMAX;
				IMAX = Problem_data.IMAX;
				ncells = Problem_data.NCells;
				maxcells = Problem_data.maxcells;
				MAX_CELL_DIM = Problem_data.MAX_CELL_DIM;				
	            slice = Cells_info.cell_slice;
	            cell_coord = Cells_info.cell_coord;
                lhsc = new double[maxcells, KMAX+2, JMAX+2, IMAX+2, 5, 5];
                backsub_info = new double[maxcells, MAX_CELL_DIM+3, MAX_CELL_DIM+3, 5];
			}
		#endregion
		
		private IIterationControl<IForwardDirection> iteration_control_forward = null;
		
		public IIterationControl<IForwardDirection> Iteration_control_forward {
			get {
				if (this.iteration_control_forward == null)
					this.iteration_control_forward = (IIterationControl<IForwardDirection>) Services.getPort("iteration_control_forward");
				return this.iteration_control_forward;
			}
		}
		
		private IIterationControl<IBackwardDirection> iteration_control_backward = null;
		
		public IIterationControl<IBackwardDirection> Iteration_control_backward {
			get {
				if (this.iteration_control_backward == null)
					this.iteration_control_backward = (IIterationControl<IBackwardDirection>) Services.getPort("iteration_control_backward");
				return this.iteration_control_backward;
			}
		}
		
		private IProblemDefinition<I, C> problem_data = null;
		
		public IProblemDefinition<I, C> Problem_data {
			get {
				if (this.problem_data == null)
					this.problem_data = (IProblemDefinition<I, C>) Services.getPort("problem_data");
				return this.problem_data;
			}
		}
		
		private IBuffer input_buffer_forward = null;
		
		protected IBuffer Input_buffer_forward {
			get {
				if (this.input_buffer_forward == null)
					this.input_buffer_forward = (IBuffer) Services.getPort("input_buffer_forward");
				return this.input_buffer_forward;
			}
		}
		
		private IBuffer input_buffer_backward = null;
		
		protected IBuffer Input_buffer_backward {
			get {
				if (this.input_buffer_backward == null)
					this.input_buffer_backward = (IBuffer) Services.getPort("input_buffer_backward");
				return this.input_buffer_backward;
			}
		}
		
		private IBuffer output_buffer_forward = null;
		
		protected IBuffer Output_buffer_forward {
			get {
				if (this.output_buffer_forward == null)
					this.output_buffer_forward = (IBuffer) Services.getPort("output_buffer_forward");
				return this.output_buffer_forward;
			}
		}
		
		private IBuffer output_buffer_backward = null;
		
		protected IBuffer Output_buffer_backward {
			get {
				if (this.output_buffer_backward == null)
					this.output_buffer_backward = (IBuffer) Services.getPort("output_buffer_backward");
				return this.output_buffer_backward;
			}
		}
		
		private IBackSubstitute<I, C, DIR, MTH> back_substitute = null;
		
		protected IBackSubstitute<I, C, DIR, MTH> Back_substitute {
			get {
				if (this.back_substitute == null)
					this.back_substitute = (IBackSubstitute<I, C, DIR, MTH>) Services.getPort("back_substitute");
				return this.back_substitute;
			}
		}
		
		private IShift<IForwardDirection> shift_forward = null;
		
		protected IShift<IForwardDirection> Shift_forward {
			get {
				if (this.shift_forward == null)
					this.shift_forward = (IShift<IForwardDirection>) Services.getPort("shift_forward");
				return this.shift_forward;
			}
		}
		
		private IUnpackSolveInfo<I, C, DIR, MTH> unpack_solve_info = null;
		
		protected IUnpackSolveInfo<I, C, DIR, MTH> Unpack_solve_info {
			get {
				if (this.unpack_solve_info == null)
					this.unpack_solve_info = (IUnpackSolveInfo<I, C, DIR, MTH>) Services.getPort("unpack_solve_info");
				return this.unpack_solve_info;
			}
		}
		
		private ISolverCell<I, C, DIR, MTH> solve_cell = null;
		
		protected ISolverCell<I, C, DIR, MTH> Solve_cell {
			get {
				if (this.solve_cell == null)
					this.solve_cell = (ISolverCell<I, C, DIR, MTH>) Services.getPort("solve_cell");
				return this.solve_cell;
			}
		}
		
		private IShift<IBackwardDirection> shift_backward = null;
		
		protected IShift<IBackwardDirection> Shift_backward {
			get {
				if (this.shift_backward == null)
					this.shift_backward = (IShift<IBackwardDirection>) Services.getPort("shift_backward");
				return this.shift_backward;
			}
		}
		
		private IMPIDirect mpi = null;
		
		public IMPIDirect Mpi {
			get {
				if (this.mpi == null)
					this.mpi = (IMPIDirect) Services.getPort("mpi");
				return this.mpi;
			}
		}
		
		private IPackBackSubInfo<I, C, DIR, MTH> pack_back_sub_info = null;
		
		protected IPackBackSubInfo<I, C, DIR, MTH> Pack_back_sub_info {
			get {
				if (this.pack_back_sub_info == null)
					this.pack_back_sub_info = (IPackBackSubInfo<I, C, DIR, MTH>) Services.getPort("pack_back_sub_info");
				return this.pack_back_sub_info;
			}
		}
		
		private ICells cells_info = null;
		
		public ICells Cells_info {
			get {
				if (this.cells_info == null)
					this.cells_info = (ICells) Services.getPort("cells_info");
				return this.cells_info;
			}
		}
		
		private ICell topology = null;
		
		public ICell Topology {
			get {
				if (this.topology == null)
					this.topology = (ICell) Services.getPort("topology");
				return this.topology;
			}
		}
		
		private IPackSolveInfo<I, C, DIR, MTH> pack_solve_info = null;
		
		protected IPackSolveInfo<I, C, DIR, MTH> Pack_solve_info {
			get {
				if (this.pack_solve_info == null)
					this.pack_solve_info = (IPackSolveInfo<I, C, DIR, MTH>) Services.getPort("pack_solve_info");
				return this.pack_solve_info;
			}
		}
		
		private IUnpackBackSubInfo<I, C, DIR, MTH> unpack_back_sub_info = null;
		
		protected IUnpackBackSubInfo<I, C, DIR, MTH> Unpack_back_sub_info {
			get {
				if (this.unpack_back_sub_info == null)
					this.unpack_back_sub_info = (IUnpackBackSubInfo<I, C, DIR, MTH>) Services.getPort("unpack_back_sub_info");
				return this.unpack_back_sub_info;
			}
		}
		abstract public int go(); 
	}
}
