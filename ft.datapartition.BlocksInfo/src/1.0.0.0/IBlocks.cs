using br.ufc.pargo.hpe.kinds;
using ft.problem_size.Instance;
using common.problem_size.Class;

namespace ft.datapartition.BlocksInfo 
{ 
	public interface IBlocks<I, C> : BaseIBlocks<I, C>
		where I:IInstance<C>
		where C:IClass 
	{
		int fftblock {get;set;}
		int fftblockpad {get;set;}
		int node {get;set;}
		int me1 {get;set;}
		int me2 {get;set;}
		int size1 {get;set;}
		int size2 {get;set;}
		int[] xstart {get;}
		int[] ystart {get;}
		int[] zstart {get;}
		int[] xend {get;}
		int[] yend {get;}
		int[] zend {get;}
		
		int[,] dims     { get; }		
		int np          { get; }
		int np1         { get; }
		int np2         { get; }
		int layout_type { get; }
		int ntdivnp     { get; }
		
		
		//void blocksConfig();
	}
}
