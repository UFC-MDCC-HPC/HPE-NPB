using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.data.Field;

namespace impl.common.data.FieldImpl { 

public class IFieldImpl : BaseIFieldImpl, IField
{
		
    private double[,,,,,] field6 = null;
		
    public double[,,,,,] Field6 { 
			get 
			{ 
				if (field6 == null) 
				{
					throw new NotInitializedFieldException(this.FieldName);
				}
				return field6; 
			} 
		}
		
    private double[,,,,] field = null;
		
    public double[,,,,] Field { 
			get 
			{ 
				if (field == null) 
				{
					throw new NotInitializedFieldException(this.FieldName);
				}
				return field; 
			} 
		}
		
	public void initialize_field(string fieldName, int maxcells, int KMAX, int JMAX, int IMAX, int neq1, int neq2) 
	{
		this.fieldName = fieldName;
		this._IMAX = IMAX;
		this._JMAX = JMAX;
		this._KMAX = KMAX;
		this.maxcells = maxcells;
		this.neq1 = neq1;
		this.neq2 = neq2;
			
		field6 = new double[maxcells, KMAX, JMAX, IMAX, neq1, neq2];
	}
		
	public void initialize_field(string fieldName, int maxcells, int KMAX, int JMAX, int IMAX, int neq) 
	{
		this.fieldName = fieldName;
		this._IMAX = IMAX;
		this._JMAX = JMAX;
		this._KMAX = KMAX;
		this.maxcells = maxcells;
		this.neq1 = neq;
			
		field = new double[maxcells, KMAX, JMAX, IMAX, neq1];
	}
		
	public void initialize_field(string fieldName, int maxcells, int KMAX, int JMAX, int IMAX) 
	{
		this.initialize_field(fieldName, maxcells, KMAX, JMAX, IMAX, 1);
	}

	private string fieldName;
	private int maxcells;
	private int _IMAX, _JMAX, _KMAX;
	private int neq1 = 1;
	private int neq2 = 1;
		
	public string FieldName { get {return fieldName; } }
    public int IMAX { get { return _IMAX; } }
	public int JMAX { get { return _JMAX; } }
	public int KMAX { get { return _KMAX; } }
	public int MaxCells { get { return maxcells; } }
	public int NEQ { get {return neq1; } }
		

}

public class NotInitializedFieldException : Exception	
{
    public NotInitializedFieldException(string fieldName) : base("Field " + fieldName + " not initialized.")
	{
			
	}
}
	
}