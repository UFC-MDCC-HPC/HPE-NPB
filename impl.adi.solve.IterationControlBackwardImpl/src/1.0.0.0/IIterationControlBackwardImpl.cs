using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.direction.Backward;
using adi.solve.IterationControl;

namespace impl.adi.solve.IterationControlBackwardImpl { 

public class IIterationControlBackwardImpl<D> : BaseIIterationControlBackwardImpl<D>, IIterationControl<D>
where D:IBackwardDirection
{		
	private int stage = -2;
	private int number_of_stages = -1;	
		
	public void setNumberOfStages(int number_of_stages)
	{
		this.number_of_stages = number_of_stages;
	}
			
	public bool has_started()
	{
		return stage > -2;
	}
				
	public bool is_first_stage()
	{
		return stage == 0;
	}
					
	public int getCurrentStage()
	{
	    return this.stage;
	}
		
	public bool is_last_stage()
	{
		return stage == this.number_of_stages - 1;
	}
					
	public bool has_finished()
	{
		return stage == -1;
	}
	
	public void start()
	{
		stage = number_of_stages - 1;
	}
			
	public void advance()
	{
		if (stage >= 0) stage--;	
			
	}
			
	public void end()
	{
		stage = -2;
	}
		
	public override int go() { 
	
		 return 0;
	}
	
}

}
