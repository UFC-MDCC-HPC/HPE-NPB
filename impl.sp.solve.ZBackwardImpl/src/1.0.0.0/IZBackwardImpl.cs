using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.axis.ZAxis;
using sp.solve.SPMethod;
using sp.solve.Backward;

namespace impl.sp.solve.ZBackwardImpl { 

	public class IZBackwardImpl<I, C, DIR, MTH> : BaseIZBackwardImpl<I, C, DIR, MTH>, IBackward<I, C, DIR, MTH>
		where I:IInstance_SP<C>
		where C:IClass
		where DIR:IZ
		where MTH:ISPMethod
	{	
		public IZBackwardImpl()
		{
		}
				
		public void init()
		{
			int i, j, k, n, k1, m;
			
			int c, kend, isize, jsize;
			
            c = slice[stage, 2];

            kend = 2 + cell_size[c, 2] - 1;

            isize = cell_size[c, 0] + 2;
            jsize = cell_size[c, 1] + 2;
			
            #region backward init
            //---------------------------------------------------------------------
            //            now we know this is the first grid block on the back sweep,
            //            so we do not need a message to start the substitution. 
            //---------------------------------------------------------------------

            k = kend - 1;
            k1 = kend;
            n = -1;
            for (m = 0; m <= 2; m++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                    lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m];
                    }
                }
            }

            for (m = 3; m <= 4; m++)
            {
                n = (m - 2) * 5 - 1;
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                    lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m];
                    }
                }
            }
            #endregion backward_init
		}		
		
		
		private int stage = -1;
		
		public void begin()
		{
			stage = ncells-1;
		}
		
		
		public bool first_stage()
		{
			return (stage == ncells - 1);
		}
		
		public bool last_stage()
		{
			return (stage == 0);
		}

		public bool finished()
		{
			return stage < 0;
		}
		
		public void advance()			
		{
			stage--;
		}

		public override int go() 
		{ 
			int c, kstart, kend, isize, jsize;
			
            c = slice[stage, 2];

            kstart = 2;
            kend = 2 + cell_size[c, 2] - 1;

            isize = cell_size[c, 0] + 2;
            jsize = cell_size[c, 1] + 2;
			
			int i, j, k, n, k1, k2, m;
			//double fac1, fac2;

			#region backward
            //---------------------------------------------------------------------
            //         Whether or not this is the last processor, we always have
            //         to complete the back-substitution 
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //         The first three factors
            //---------------------------------------------------------------------
            n = -1;
            for (m = 0; m <= 2; m++)
            {
                for (k = kend - 2; k >= kstart; k--)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            k1 = k + 1;
                            k2 = k + 2;
                            rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                     lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m] -
                                     lhs[c, k, j, i, n + 5] * rhs[c, k2, j, i, m];
                        }
                    }
                }
            }

            //---------------------------------------------------------------------
            //         And the remaining two
            //---------------------------------------------------------------------
            for (m = 3; m <= 4; m++)
            {
                n = (m - 2) * 5 - 1;
                for (k = kend - 2; k >= kstart; k--)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                        {
                            k1 = k + 1;
                            k2 = k + 2;
                            rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                                     lhs[c, k, j, i, n + 4] * rhs[c, k1, j, i, m] -
                                     lhs[c, k, j, i, n + 5] * rhs[c, k2, j, i, m];
                        }
                    }
                }
            }
            #endregion
			
			return 0;
		} // end activate method 
		
	}

}
