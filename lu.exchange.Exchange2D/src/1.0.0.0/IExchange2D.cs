using br.ufc.pargo.hpe.kinds;
using lu.problem_size.Instance;
using common.problem_size.Class;
using lu.exchange.ExchangePattern;

namespace lu.exchange.Exchange2D { 

public interface IExchange2D<I, C, E> : BaseIExchange2D<I, C, E>
where I:IInstance<C>
where C:IClass
where E:IExchangePattern
{

  void setParameters(double[,] g, double[,] h, int ibeg, int ifin1, int jbeg, int jfin1);
  void setParameters(double[,] g, int beg, int fin1);

} // end main interface 

} // end namespace 
