using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Backward;
using common.axis.YAxis;
using sp.solve.shift.WriteBuffer;

namespace impl.sp.solve.shift.YWriteBufferBackward { 

public class IYWriteBufferBackwardImpl<I, C, D, DIR> : BaseIYWriteBufferBackwardImpl<I, C, D, DIR>, IWriteBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IBackwardDirection
where DIR:IY
{

private double[][] out_buffer_solver;
double[] out_buffer_y;

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
}

#pragma action
public void advance()
{
    stage--;
}
		
#pragma condition		
public bool finished()
{
	return stage < 0;
}
		
void create_buffers() 		
{
	int c, stage, isize, ksize, buffer_size;
	
	out_buffer_solver = new double[ncells][];
	
	for (stage = 0; stage < ncells; stage++)
	{
		c = slice[stage, 1];
		
		isize = cell_size[c, 0] + 2;
		ksize = cell_size[c, 2] + 2;
						
		buffer_size = (isize - start[c, 0] - end[c, 0]) * (ksize - start[c, 2] - end[c, 2]);

		out_buffer_solver[stage] = new double[10 * buffer_size];
						
	}
	
}

bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
    Buffer.Array = out_buffer_y = out_buffer_solver[stage];
			
	int i, j, k, isize, ksize, j1, c, m, p, jstart; /* requests(2), statuses(MPI_STATUS_SIZE, 2);*/
	
	c = slice[stage, 1];
	
	jstart = 2;
	
	isize = cell_size[c, 0] + 2;
	ksize = cell_size[c, 2] + 2;
			
    j = jstart;
    j1 = jstart + 1;
    p = 0;
    for (m = 0; m <= 4; m++)
    {
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
                out_buffer_y[p] = rhs[c, k, j, i, m];
                out_buffer_y[p + 1] = rhs[c, k, j1, i, m];
                p = p + 2;
            }
        }
    }
			
	return 0;
}

}

}
