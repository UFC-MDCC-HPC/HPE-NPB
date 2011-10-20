/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using lu.datapartition.BlocksInfo;
using lu.problem_size.Instance_LU;
using common.problem_size.Class;

namespace impl.lu.datapartition.BlocksInfoImpl 
{ 
	public abstract class BaseIBlocksInfoImpl<I, C>: br.ufc.pargo.hpe.kinds.Environment, BaseIBlocksInfo<I, C> 
		where I:IInstance_LU<C>
		where C:IClass
	{
		
		private I instance = default(I);
		
		protected I Instance {
			get {
				if (this.instance == null)
					this.instance = (I) Services.getPort("instance_type");
				return this.instance;
			}
		}
	
	}
}
