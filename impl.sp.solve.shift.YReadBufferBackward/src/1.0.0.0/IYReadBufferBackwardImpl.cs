using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Backward;
using common.axis.YAxis;
using sp.solve.shift.ReadBuffer;

namespace impl.sp.solve.shift.YReadBufferBackward { 

public class IYReadBufferBackwardImpl<I, C, D, DIR> : BaseIYReadBufferBackwardImpl<I, C, D, DIR>, IReadBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IBackwardDirection
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
			
	stage = ncells - 1;
    Buffer.Array = in_buffer_y = in_buffer_solver[stage];
}

#pragma action
public void advance()
{
    Buffer.Array = in_buffer_y = in_buffer_solver[--stage];
}
		
#pragma condition		
public bool finished()
{
	return stage < 0;
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

		in_buffer_solver[stage] = new double[10 * buffer_size];
	}
}

bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
	int i, j, k, n, isize, jend, ksize, j1, c, m, p; /* requests(2), statuses(MPI_STATUS_SIZE, 2);*/
	double sm1, sm2;
	
	c = slice[stage, 1];

	jend = 2 + cell_size[c, 1] - 1;
	
	isize = cell_size[c, 0] + 2;
	ksize = cell_size[c, 2] + 2;

    //---------------------------------------------------------------------
    //            unpack the buffer for the first three factors         
    //---------------------------------------------------------------------
    n = -1;
    p = 0;
    j = jend;
    j1 = j - 1;
    for (m = 0; m <= 2; m++)
    {
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
          //      Console.WriteLine("in_buffer_y["+ p + "] = " + in_buffer_y[p]);
          //      Console.WriteLine("in_buffer_y["+ (p+1) + "] = " + in_buffer_y[p+1]);
                sm1 = in_buffer_y[p];
                sm2 = in_buffer_y[p + 1];
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                       lhs[c, k, j, i, n + 4] * sm1 -
                       lhs[c, k, j, i, n + 5] * sm2;
                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                       lhs[c, k, j1, i, n + 4] * rhs[c, k, j, i, m] -
                       lhs[c, k, j1, i, n + 5] * sm1;
                p = p + 2;
            }
        }
    }

    //---------------------------------------------------------------------
    //            now unpack the buffer for the remaining two factors
    //---------------------------------------------------------------------
    for (m = 3; m <= 4; m++)
    {
        n = (m - 2) * 5 - 1;
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
           //     Console.WriteLine("in_buffer_y["+ p + "] = " + in_buffer_y[p]);
           //     Console.WriteLine("in_buffer_y["+ (p+1) + "] = " + in_buffer_y[p+1]);
                sm1 = in_buffer_y[p];
                sm2 = in_buffer_y[p + 1];
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                       lhs[c, k, j, i, n + 4] * sm1 -
                       lhs[c, k, j, i, n + 5] * sm2;
                rhs[c, k, j1, i, m] = rhs[c, k, j1, i, m] -
                       lhs[c, k, j1, i, n + 4] * rhs[c, k, j, i, m] -
                       lhs[c, k, j1, i, n + 5] * sm1;
                p = p + 2;
            }
        }
    }
			
	 return 0;
}

}

}
