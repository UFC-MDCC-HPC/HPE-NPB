using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Forward;
using common.axis.YAxis;
using sp.solve.shift.ReadBuffer;

namespace impl.sp.solve.shift.YReadBufferForward { 

public class IYReadBufferForwardImpl<I, C, D, DIR> : BaseIYReadBufferForwardImpl<I, C, D, DIR>, IReadBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IForwardDirection
where DIR:IY
{

private double[][] in_buffer_solver;
double[] in_buffer_y;

private int stage = -1;

#pragma action
public void begin()
{
	if (!buffer_created)		
	{
		create_buffers();
		buffer_created = true;
	}
			
	stage = 0;
    Buffer.Array = in_buffer_y = in_buffer_solver[stage];
}

#pragma action
public void advance()
{
    stage++;
    Buffer.Array = in_buffer_y = in_buffer_solver[stage];
}
		
#pragma condition		
public bool finished()
{
	return stage >= ncells;
}
		
void create_buffers() 		
{
	int c, stage, isize, ksize, buffer_size;
	
	in_buffer_solver = new double[ncells][];			
	
	for (stage = 0; stage < ncells; stage++)
	{
		c = slice[stage, 1];
		
		isize = cell_size[c, 0] + 2;
		ksize = cell_size[c, 2] + 2;
						
		buffer_size = (isize - start[c, 0] - end[c, 0]) * (ksize - start[c, 2] - end[c, 2]);

		in_buffer_solver[stage] = new double[22*buffer_size];
		
	}
	
}

bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 			
	int i, j, k, n, isize, ksize, j1, c, m, p, jstart; /* requests(2), statuses(MPI_STATUS_SIZE, 2);*/
	double r1, r2, d, e;
	double[] s = new double[5];
	
	c = slice[stage, 1];
	
	jstart = 2;
	
	isize = cell_size[c, 0] + 2;
	ksize = cell_size[c, 2] + 2;
			
    //---------------------------------------------------------------------
    //            unpack the buffer                                 
    //---------------------------------------------------------------------
    j = jstart;
    j1 = jstart + 1;
    n = -1;
    //---------------------------------------------------------------------
    //            create a running pointer
    //---------------------------------------------------------------------
    p = 0;
    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
    {
        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
        {
      //      Console.WriteLine("in_buffer_y_f["+ p + "] = " + in_buffer_y[p]);
      //      Console.WriteLine("in_buffer_y_f["+ (p+1) + "] = " + in_buffer_y[p+1]);
            lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                    in_buffer_y[p] * lhs[c, k, j, i, n + 1];
            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                    in_buffer_y[p + 1] * lhs[c, k, j, i, n + 1];
            for (m = 0; m <= 2; m++)
            {
      //          Console.WriteLine("in_buffer_y_f["+ (p+2+m) + "] = " + in_buffer_y[p+2+m]);
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                      in_buffer_y[p + 2 + m] * lhs[c, k, j, i, n + 1];
            }
      //      Console.WriteLine("in_buffer_y_f["+ (p+5) + "] = " + in_buffer_y[p+5]);
      //      Console.WriteLine("in_buffer_y_f["+ (p+6) + "] = " + in_buffer_y[p+6]);
            d = in_buffer_y[p + 5]; ;
            e = in_buffer_y[p + 6];
            for (m = 0; m <= 2; m++)
            {
       //        Console.WriteLine("in_buffer_y_f["+ (p+7+m) + "] = " + in_buffer_y[p+7+m]);
                s[m] = in_buffer_y[p + 7 + m];
            }
            r1 = lhs[c, k, j, i, n + 2];
            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
            lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
            for (m = 0; m <= 2; m++)
            {
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
            }
            r2 = lhs[c, k, j1, i, n+1];
            lhs[c, k, j1, i, n + 2] = lhs[c, k, j1, i, n + 2] - d * r2;
            lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] - e * r2;
            for (m = 0; m <= 2; m++)
            {
                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] - s[m] * r2;
            }
            p = p + 10;
        }
    }

    for (m = 3; m <= 4; m++)
    {
        n = (m - 2) * 5 - 1;
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
                lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                         in_buffer_y[p] * lhs[c, k, j, i, n + 1];
                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                         in_buffer_y[p + 1] * lhs[c, k, j, i, n + 1];
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                         in_buffer_y[p + 2] * lhs[c, k, j, i, n + 1];
                d = in_buffer_y[p + 3];
                e = in_buffer_y[p + 4];
                s[m] = in_buffer_y[p + 5];
                r1 = lhs[c, k, j, i, n + 2];
                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                r2 = lhs[c, k, j1, i, n + 1];
                lhs[c, k, j1, i, n + 2] = lhs[c, k, j1, i, n + 2] - d * r2;
                lhs[c, k, j1, i, n + 3] = lhs[c, k, j1, i, n + 3] - e * r2;
                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] - s[m] * r2;
                p = p + 6;
            }
        }
    }

			return 0;
}

}

}
