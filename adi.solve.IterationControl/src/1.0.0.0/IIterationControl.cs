using br.ufc.pargo.hpe.kinds;
using common.direction.Direction;

namespace adi.solve.IterationControl { 

public interface IIterationControl<DIR> : BaseIIterationControl<DIR>
where DIR:IDirection
{
   void setNumberOfStages(int max);
   int getCurrentStage();
		
   bool has_started();
   bool is_first_stage();
   bool is_last_stage();
   bool has_finished();

   void start();
   void advance();
   void end();
   
} // end main interface 

} // end namespace 
