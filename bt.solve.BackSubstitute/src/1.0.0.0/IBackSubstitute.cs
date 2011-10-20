using br.ufc.pargo.hpe.kinds;
using common.solve.SolvingMethod;
using common.problem_size.Instance;
using common.problem_size.Class;
using common.axis.Axis;

namespace bt.solve.BackSubstitute { 

public interface IBackSubstitute<I, C, DIR, MTH> : BaseIBackSubstitute<I, C, DIR, MTH>
where I:IInstance<C>
where C:IClass
where DIR:IAxis
where MTH:ISolvingMethod
{
   void setParameters(double[,,,,,] lhsc, double[,,,] backsub_info, int first, int last, int c);

} // end main interface 

} // end namespace 
