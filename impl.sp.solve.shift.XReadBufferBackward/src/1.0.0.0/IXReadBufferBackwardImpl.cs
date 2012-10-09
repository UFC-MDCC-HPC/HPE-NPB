using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Backward;
using common.axis.XAxis;
using sp.solve.shift.ReadBuffer;

namespace impl.sp.solve.shift.XReadBufferBackward { 

public class IXReadBufferBackwardImpl<I, C, D, DIR> : BaseIXReadBufferBackwardImpl<I, C, D, DIR>, IReadBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IBackwardDirection
where DIR:IX
{
		
private double[][] in_buffer_solver;
double[] in_buffer_x;

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
    Buffer.Array = in_buffer_x = in_buffer_solver[stage];
}

#pragma action
public void advance()
{
    stage--;
    Buffer.Array = in_buffer_x = in_buffer_solver[stage];
}
		
#pragma condition		
public bool finished()
{
	return stage < 0;
}
		
void create_buffers() 		
{
	int c, stage, jsize, ksize, buffer_size;
	
	in_buffer_solver = new double[ncells][];			
	
	for (stage = 0; stage < ncells; stage++)
    {
        c = slice[stage, 0];
		
        jsize = cell_size[c, 1] + 2;
        ksize = cell_size[c, 2] + 2;
		
		buffer_size = (jsize - start[c, 1] - end[c, 1]) * (ksize - start[c, 2] - end[c, 2]);
		
        in_buffer_solver[stage] = new double[10*buffer_size];
	}			
}
		
bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 						
    int c, i, j, k, n, iend, jsize, ksize, i1, m, p;
    double sm1, sm2;
	
    c = slice[stage, 0];

    iend = 2 + cell_size[c, 0] - 1;

    jsize = cell_size[c, 1] + 2;
    ksize = cell_size[c, 2] + 2;
			
	//---------------------------------------------------------------------
    //            unpack the buffer for the first three factors         
    //---------------------------------------------------------------------
    n = -1;
    p = 0;
    i = iend;
    i1 = i - 1;
    for (m = 0; m <= 2; m++)
    {
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
            {
                sm1 = in_buffer_x[p];
                sm2 = in_buffer_x[p + 1];
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                       lhs[c, k, j, i, n + 4] * sm1 -
                       lhs[c, k, j, i, n + 5] * sm2;
                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                       lhs[c, k, j, i1, n + 4] * rhs[c, k, j, i, m] -
                       lhs[c, k, j, i1, n + 5] * sm1;
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
            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
            {
                sm1 = in_buffer_x[p];
                sm2 = in_buffer_x[p + 1];
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                       lhs[c, k, j, i, n + 4] * sm1 -
                       lhs[c, k, j, i, n + 5] * sm2;
                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] -
                       lhs[c, k, j, i1, n + 4] * rhs[c, k, j, i, m] -
                       lhs[c, k, j, i1, n + 5] * sm1;
                p = p + 2;
            }
        }
    }

	return 0;
}

}

}
