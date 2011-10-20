/* Automatically Generated Code */

using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using ft.datapartition.BlocksInfo;
using common.problem_size.Class;
using ft.problem_size.Instance_FT;

namespace impl.ft.datapartition.BlocksInfoImpl 
{ 
	public abstract class BaseIBlocksInfoImpl<I, C>: br.ufc.pargo.hpe.kinds.Environment, BaseIBlocks<I, C>
		where I:IInstance_FT<C>
		where C:IClass	
	{
		public BaseIBlocksInfoImpl()
		{
	
		}
		
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
