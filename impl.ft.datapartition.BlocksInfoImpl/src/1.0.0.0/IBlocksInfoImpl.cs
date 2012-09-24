using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using ft.datapartition.BlocksInfo;
using common.problem_size.Class;
using ft.problem_size.Instance_FT;
using MPI;

namespace impl.ft.datapartition.BlocksInfoImpl 
{ 
	public class IBlocksInfoImpl<I, C> : BaseIBlocksInfoImpl<I, C>, IBlocks<I, C>
		where I:IInstance_FT<C>
		where C:IClass	
	{
		public static int layout_0D = 0, layout_1D = 1, layout_2D = 2;
		
		override public void initialize()
		{			
            niter_default = Instance.niter_default;
		    nx = Instance.nx;
			ny = Instance.ny;
			nz = Instance.nz;			
			
			initialConfig();
			problemConfig();			
			blocksConfig();			
		}
		
		protected int niter_default;
		private int root=0;

		public static double mod(double a, double b) { return (a % b); }
		
		private int _fftblock, _fftblockpad, _node, _me1, _me2, _size1, _size2;
		private int[] _xstart = {0,0,0};
		private int[] _ystart = {0,0,0};
		private int[] _zstart = {0,0,0};
		private int[] _xend = {0,0,0};
		private int[] _yend = {0,0,0};
		private int[] _zend = {0,0,0};
		
		public int fftblock    {get { return _fftblock; }    set { _fftblock = value; } }
		public int fftblockpad {get { return _fftblockpad; } set { _fftblockpad = value; } }
		public int node        {get { return _node; }        set { _node = value; } }
		public int me1         {get { return _me1; }         set { _me1 = value; } }
		public int me2         {get { return _me2; }         set { _me2 = value; } }
		
		public int size1       { get { return _size1;       } set { _size1 = value;       } }
		public int size2       { get { return _size2;       } set { _size2 = value;       } }
		
		public int[] xstart {get { return _xstart; } }
		public int[] ystart {get { return _ystart; } }
		public int[] zstart {get { return _zstart; } }
		public int[] xend   {get { return _xend; } }
		public int[] yend   {get { return _yend; } }
		public int[] zend   {get { return _zend; } }
		
		public int[,] dims            { get { return _dims;    } }		
		public int np          { get { return _np;          } set { _np = value;          } }
		public int np1         { get { return _np1;         } set { _np1 = value;         } }
		public int np2         { get { return _np2;         } set { _np2 = value;         } }
		public int layout_type { get { return _layout_type; } set { _layout_type = value; } }
		public int ntdivnp     { get { return _ntdivnp;     } set { _ntdivnp = value;     } }
		
		protected int[,] _dims = new int[3, 3];
		protected int _np,_np1,_np2,_layout_type,_ntdivnp, nx, ny, nz;
		

		public void blocksConfig()
		{
		    int layout_0D = 0, layout_1D = 1, layout_2D = 2;
		    int fftblock_default=16, fftblockpad_default=18;
		    int node = this.GlobalRank;
            _me2 = (int)mod(node, np2);
            _me1 = node/np2;
            if(layout_type == layout_0D) {
                for(int i = 0; i < 3; i++) {
                    _xstart[i] = 1;
                    _xend[i]   = nx;
                    _ystart[i] = 1;
                    _yend[i]   = ny;
                    _zstart[i] = 1;
                    _zend[i]   = nz;
                }

            }
            else if(layout_type == layout_1D) {
                _xstart[0] = 1;
                _xend[0]   = nx;
                _ystart[0] = 1;
                _yend[0]   = ny;
                _zstart[0] = 1 + _me2 * nz/np2;
                _zend[0]   = (_me2+1) * nz/np2;

                _xstart[1] = 1;
                _xend[1]   = nx;
                _ystart[1] = 1;
                _yend[1]   = ny;
                _zstart[1] = 1 + _me2 * nz/np2;
                _zend[1]   = (_me2+1) * nz/np2;

                _xstart[2] = 1;
                _xend[2]   = nx;
                _ystart[2] = 1 + _me2 * ny/np2;
                _yend[2]   = (_me2+1) * ny/np2;
                _zstart[2] = 1;
                _zend[2] = nz;

            }
            else if(layout_type == layout_2D) {

                _xstart[0] = 1;
                _xend[0]   = nx;
                _ystart[0] = 1 + _me1 * ny/np1;
                _yend[0]   = (_me1+1) * ny/np1;
                _zstart[0] = 1 + _me2 * nz/np2;
                _zend[0]   = (_me2+1) * nz/np2;

                _xstart[1] = 1 + _me1 * nx/np1;
                _xend[1]   = (_me1+1)*nx/np1;
                _ystart[1] = 1;
                _yend[1]   = ny;
                _zstart[1] = _zstart[0];
                _zend[1]   = _zend[0];

                _xstart[2] = _xstart[1];
                _xend[2]   = _xend[1];
                _ystart[2] = 1 + _me2 *ny/np2;
                _yend[2]   = (_me2+1)*ny/np2;
                _zstart[2] = 1;
                _zend[2] = nz;
            }
            _fftblock = fftblock_default;
            _fftblockpad = fftblockpad_default;

            int dim1 = ny/np1;
            int dim2 = nx/np1;
            int dim3 = nx/np1;
            if(layout_type == layout_2D) {
                if(dim1 < _fftblock)
                    _fftblock = dim1;
                if(dim2 < _fftblock)
                    _fftblock = dim2;
                if(dim3 < _fftblock)
                    _fftblock = dim3;
            }

            if(_fftblock != fftblock_default)
                _fftblockpad = _fftblock + 3;
                
            _size1 = ((int)(nz/np2))*nx*2;
            _size2 = nx*2;
			
     		Console.Error.WriteLine("size1=" + _size1 + ", size2="+ _size2);
		}
		
		
		
        public void problemConfig()
		{
            _np = this.Size;
			_ntdivnp = ((nx*ny)/_np)*nz;

            if(_layout_type == layout_0D) 
			{
                for(int i = 0; i < 3; i++) 
				{
                    _dims[0, i] = nx;
                    _dims[1, i] = ny;
                    _dims[2, i] = nz;
                }
            }
            else if(_layout_type == layout_1D) 
			{
                _dims[0, 0] = nx;
                _dims[1, 0] = ny;
                _dims[2, 0] = nz;

                _dims[0, 1] = nx;
                _dims[1, 1] = ny;
                _dims[2, 1] = nz;

                _dims[0, 2] = nz;
                _dims[1, 2] = nx;
                _dims[2, 2] = ny;
            }
            else if(_layout_type == layout_2D) 
			{
                _dims[0, 0] = nx;
                _dims[1, 0] = ny;
                _dims[2, 0] = nz;

                _dims[0, 1] = ny;
                _dims[1, 1] = nx;
                _dims[2, 1] = nz;

                _dims[0, 2] = nz;
                _dims[1, 2] = nx;
                _dims[2, 2] = ny;

            }
			
            _dims[1, 0] = _dims[1, 0] / _np1;
            _dims[2, 0] = _dims[2, 0] / _np2;
            _dims[1, 1] = _dims[1, 1] / _np1;
            _dims[2, 1] = _dims[2, 1] / _np2;
            _dims[1, 2] = _dims[1, 2] / _np1;
            _dims[2, 2] = _dims[2, 2] / _np2;
			
        }
		
		

        public int initialConfig() 
		{            
			Intracommunicator worldcomm = this.WorldComm;
			
            np = worldcomm.Size;
            int node = this.GlobalRank;
            int niter = niter_default;
			
            if(node == 0) 
			{
                Console.WriteLine(" NAS Parallel Benchmarks 3.3 -- FT Benchmark ");

                if(np == 1) {
                    np1 = 1;
                    np2 = 1;
                    layout_type = layout_0D;
                }
                else if(np <= nz) {
                    np1 = 1;
                    np2 = np;
                    layout_type = layout_1D;
                }
                else {
                    np1 = nz;
                    np2 = np/nz;
                    layout_type = layout_2D;
                }
                
                Console.WriteLine(" Size: " + nx + "x" + ny + "x" + nz);
                Console.WriteLine(" Iterations: "+ niter);
                Console.WriteLine(" Number of processes : "+ np);
                Console.WriteLine(" Processor array     : "+np1+"x"+np2);
				
                if(layout_type == layout_0D) 
				{
                    Console.WriteLine(" Layout type: OD");
                }
                else if(layout_type == layout_1D) 
				{
                    Console.WriteLine(" Layout type: 1D");
                }
                else 
				{
                    Console.WriteLine(" Layout type: 2D");
                }
            }
			
            worldcomm.Broadcast<int>(ref _np1, root);
            worldcomm.Broadcast<int>(ref niter, root);
            worldcomm.Broadcast<int>(ref _np2, root);
			
            if(np1 == 1 && np2 == 1) 
			{
                layout_type = layout_0D;
            }
            else if(np1 == 1) 
			{
                layout_type = layout_1D;
            }
            else 
			{
                layout_type = layout_2D;
            }
			
            return niter;
        }
		
	}
}
