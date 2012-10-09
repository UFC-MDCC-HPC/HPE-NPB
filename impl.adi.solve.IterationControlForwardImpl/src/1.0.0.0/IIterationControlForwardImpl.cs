using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.direction.Forward;
using adi.solve.IterationControl;

namespace impl.adi.solve.IterationControlForwardImpl { 

public class IIterationControlForwardImpl<D> : BaseIIterationControlForwardImpl<D>, IIterationControl<D>
where D:IForwardDirection
{
	private int stage = -2;
	private int number_of_stages = -1;	
		
	public void setNumberOfStages(int number_of_stages)
	{
		this.number_of_stages = number_of_stages;
	}
		
	public int getCurrentStage()
	{
	    return this.stage;
	}

	public bool has_started()
	{
		return stage > -2;
	}
				
	public bool is_first_stage()
	{
		return stage == 0;
	}
					
	public bool is_last_stage()
	{
		return stage == this.number_of_stages - 1;
	}
					
	public bool has_finished()
	{
		return stage == this.number_of_stages;
	}
	
	public void start()
	{
		stage = 0;
	}
			
	public void advance()
	{
		if (this.stage < this.number_of_stages)
				stage++;
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
