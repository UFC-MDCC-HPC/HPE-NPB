using br.ufc.pargo.hpe.kinds;

namespace common.data.Field { 

public interface IField : BaseIField
{

	double[,,,,] Field { get; }
	double[,,,,,] Field6 { get; }
		
	void initialize_field(string fieldName, int max_cells, int IMAX, int JMAX, int KMAX);	
	void initialize_field(string fieldName, int max_cells, int IMAX, int JMAX, int KMAX, int neq);	
	void initialize_field(string fieldName, int max_cells, int IMAX, int JMAX, int KMAX, int neq1, int neq2);	
		
	string FieldName { get; }
    int IMAX { get; }
	int JMAX { get; }
	int KMAX { get; }
	int MaxCells { get; }

} // end main interface 

} // end namespace 
