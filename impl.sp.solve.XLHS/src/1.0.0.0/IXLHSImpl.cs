using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using adi.data.ProblemDefinition;
using sp.solve.LHS;
using common.problem_size.Class;
using sp.problem_size.Instance_SP;
using common.axis.XAxis;
using sp.solve.SPMethod;

namespace impl.sp.solve.XLHS { 

	public class IXLHSImpl<I,C,DIR,MTH> : BaseIXLHSImpl<I,C,DIR,MTH>, ILHS<I,C,DIR,MTH>
	where I:IInstance_SP<C>
	where C:IClass
	where DIR:IX
	where MTH:ISPMethod
	{
	
		public IXLHSImpl() 
		{ 
		} 
			
		private int c;

		private int stage = -1;
		
		public void begin()
		{
			stage = 0;
		}
		
		public bool finished()
		{
			return stage >= ncells;
		}
		
		public void advance()			
		{
			stage++;            	
		}
			
						
		public override int go() 
		{ 
            double ru1;
            int i, j, k;
			
			c = slice[stage, 0];
			
            int ksize = cell_size[c, 2] + 2;
            int jsize = cell_size[c, 1] + 2;
            int isize = cell_size[c, 0] + 2;
			
			if (this.Rank==1)
			Console.WriteLine("c={0}, start[c,0]={1}, start[c,1]={2}, start[c,2]={3}", c, start[c,0], start[c,1], start[c,2]);
			
            //---------------------------------------------------------------------
            //      treat only cell c             
            //---------------------------------------------------------------------

            //---------------------------------------------------------------------
            //      first fill the lhs for the u-eigenvalue                   
            //---------------------------------------------------------------------
            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0] - 1; i < isize - end[c, 0] + 1; i++)
                    {
                        ru1 = c3c4 * rho_i[c, k, j, i, 0];
                        cv[i] = us[c, k, j, i, 0];
                        rhon[i] = dmax1(dx2 + con43 * ru1, dx5 + c1c5 * ru1, dxmax + ru1, dx1);
                    }

                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0] = 0.0d;
                        lhs[c, k, j, i, 1] = -dttx2 * cv[i - 1] - dttx1 * rhon[i - 1];
                        lhs[c, k, j, i, 2] = 1.0d + c2dttx1 * rhon[i];
                        lhs[c, k, j, i, 3] = dttx2 * cv[i + 1] - dttx1 * rhon[i + 1];
                        lhs[c, k, j, i, 4] = 0.0d;
						
					//	if (this.Rank == 1) 
					//		Console.WriteLine("LHS - {0}, {1}, {2}, {3}, {4}, {5}",  lhs[c, k, j, i, 0], lhs[c, k, j, i, 1], lhs[c, k, j, i, 2], lhs[c, k, j, i, 3], lhs[c, k, j, i, 0], lhs[c, k, j, i, 4]);
				/*		if (this.Rank == 1) 
							Console.WriteLine("dttx2={0}, " +
								"cv[{8}-1]={1}, " +
								"dttx1={2}, " +
								"rhon[{8}-1]={3}, " +
								"c2dttx1={4}, " +
								"rhon[{8}]={5}, " +
								"cv[{8}+1]={6}, " +
								"rhon[{8}+1]={7}",  
							                  dttx2, 
							                  cv[i-1], 
							                  dttx1, 
							                  rhon[i-1], 
							                  c2dttx1, 
							                  rhon[i], 
							                  cv[i+1], 
							                  rhon[i+1],
							                  i);*/
                    }
                }
            }

            //---------------------------------------------------------------------
            //      add fourth order dissipation                             
            //---------------------------------------------------------------------
            if (start[c, 0] > 2)
            {
                i = 3;
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz5;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;

                        lhs[c, k, j, i + 1, 1] = lhs[c, k, j, i + 1, 1] - comz4;
                        lhs[c, k, j, i + 1, 2] = lhs[c, k, j, i + 1, 2] + comz6;
                        lhs[c, k, j, i + 1, 3] = lhs[c, k, j, i + 1, 3] - comz4;
                        lhs[c, k, j, i + 1, 4] = lhs[c, k, j, i + 1, 4] + comz1;
                    }
                }
            }

            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = 3 * start[c, 0] - 4; i < isize - 3 * end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;
                        lhs[c, k, j, i, 4] = lhs[c, k, j, i, 4] + comz1;
                    }
                }
            }

            if (end[c, 0] > 0)
            {
                i = isize - 3;
                for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
                {
                    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                    {
                        lhs[c, k, j, i, 0] = lhs[c, k, j, i, 0] + comz1;
                        lhs[c, k, j, i, 1] = lhs[c, k, j, i, 1] - comz4;
                        lhs[c, k, j, i, 2] = lhs[c, k, j, i, 2] + comz6;
                        lhs[c, k, j, i, 3] = lhs[c, k, j, i, 3] - comz4;

                        lhs[c, k, j, i + 1, 0] = lhs[c, k, j, i + 1, 0] + comz1;
                        lhs[c, k, j, i + 1, 1] = lhs[c, k, j, i + 1, 1] - comz4;
                        lhs[c, k, j, i + 1, 2] = lhs[c, k, j, i + 1, 2] + comz5;
                    }
                }
            }

            //---------------------------------------------------------------------
            //      subsequently, fill the other factors (u+c), (u-c) by a4ing to 
            //      the first  
            //---------------------------------------------------------------------
            for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
            {
                for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
                {
                    for (i = start[c, 0]; i < isize - end[c, 0]; i++)
                    {
                        lhs[c, k, j, i, 0 + 5] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 5] = lhs[c, k, j, i, 1] -
                                         dttx2 * speed[c, k, j, i - 1, 0];
                        lhs[c, k, j, i, 2 + 5] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 5] = lhs[c, k, j, i, 3] +
                                         dttx2 * speed[c, k, j, i + 1, 0];
                        lhs[c, k, j, i, 4 + 5] = lhs[c, k, j, i, 4];
                        lhs[c, k, j, i, 0 + 10] = lhs[c, k, j, i, 0];
                        lhs[c, k, j, i, 1 + 10] = lhs[c, k, j, i, 1] +
                                         dttx2 * speed[c, k, j, i - 1, 0];
                        lhs[c, k, j, i, 2 + 10] = lhs[c, k, j, i, 2];
                        lhs[c, k, j, i, 3 + 10] = lhs[c, k, j, i, 3] -
                                         dttx2 * speed[c, k, j, i + 1, 0];
                        lhs[c, k, j, i, 4 + 10] = lhs[c, k, j, i, 4];
                    }
                }
            }
			
			return 0;
		} // end activate method 
		
		public double dmax1(double a, double b)
		{
			if (a < b) return b; else return a;
		}
	
		public double dmax1(double a, double b, double c, double d)
		{
			return dmax1(dmax1(a, b), dmax1(c, d));
		}
	
	}

}
