using br.ufc.pargo.hpe.kinds;

namespace common.datapartition.MultiPartitionCells { 

public interface ICells : BaseICells
{
  int[,] cell_coord {get;}
  int[,] cell_size {get;}
  int[,] cell_low {get;}
  int[,] cell_high {get;}
  int[,] cell_slice {get;}
  
  int[,] cell_start {get;}
  int[,] cell_end {get;}

  int ncells {get;set;} // sqrt(nodes)
  int ndirs {get;set;}  // Assumed to be 3 ?
		
  int[] grid_points { get; }
		
 // void initialize();
		
} // end main interface 

} // end namespace 
