using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using ft.problem_size.Instance_FT;
using common.problem_size.Class;
using ft.data.ProblemDefinition;
using MPI;

namespace impl.ft.data.ProblemDefinitionImpl 
{ 
	public class IProblemDefinitionImpl<I, C> : BaseIProblemDefinitionImpl<I, C>, IProblemDefinition<I, C>
	where I:IInstance_FT<C>
	where C:IClass
	{				
      	protected double[]  _twiddle;
		protected double[,] _u;

		public double [,,,] Field_u0     { get { return U0.Field; } }
		public double [,,,] Field_u1     { get { return U1.Field; } }
		public double [,,,] Field_u2     { get { return U2.Field; } }
		public double[] twiddle          { get { return _twiddle; } }
		public double[,] u               { get { return _u;       } }

		override public void initialize()
		{           
			int[,] _dims = Blocks.dims;
			int _ntdivnp = Blocks.ntdivnp;
			int nx = Instance.nx;
			
			
            U0.initialize_field("u0", _dims[1, 0], _dims[2, 0], _dims[0, 0], 2);
            U1.initialize_field("u1", _dims[1, 0], _dims[2, 0], _dims[0, 0], 2);
            U2.initialize_field("u2", _dims[1, 0], _dims[2, 0], _dims[0, 0], 2);

            _u = new double[nx, 2];
            _twiddle = new double[_ntdivnp];
			
        }
		
	}
}
