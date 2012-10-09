using br.ufc.pargo.hpe.kinds;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.direction.Direction;
using common.axis.Axis;

namespace sp.solve.shift.WriteBuffer { 

public interface IWriteBuffer<I, C, D, DIR> : BaseIWriteBuffer<I, C, D, DIR>
where I:IInstance<C>
where C:IClass
where D:IDirection
where DIR:IAxis
{
   void begin();
   void advance();
   bool finished();

} // end main interface 

} // end namespace 
