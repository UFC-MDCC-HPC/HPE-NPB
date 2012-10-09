using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Backward;
using common.axis.XAxis;
using sp.solve.shift.WriteBuffer;

namespace impl.sp.solve.shift.XWriteBufferBackward { 

public class IXWriteBufferBackwardImpl<I, C, D, DIR> : BaseIXWriteBufferBackwardImpl<I, C, D, DIR>, IWriteBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IBackwardDirection
where DIR:IX
{

private double[][] out_buffer_solver;
double[] out_buffer_x;

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
	int c, stage, jsize, ksize, buffer_size;
	
	out_buffer_solver = new double[ncells][];			
	
	for (stage = 0; stage < ncells; stage++)
    {
        c = slice[stage, 0];
		
        jsize = cell_size[c, 1] + 2;
        ksize = cell_size[c, 2] + 2;
		
		buffer_size = (jsize - start[c, 1] - end[c, 1]) * (ksize - start[c, 2] - end[c, 2]);
		
        out_buffer_solver[stage] = new double[10*buffer_size];								
	}			
}
		
bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
    Buffer.Array = out_buffer_x = out_buffer_solver[stage];
			
    int c, i, j, k, jsize, ksize, i1, m, p, istart;

    c = slice[stage, 0];

    istart = 2;

    jsize = cell_size[c, 1] + 2;
    ksize = cell_size[c, 2] + 2;

	i = istart;
    i1 = istart + 1;
    p = 0;
    for (m = 0; m <= 4; m++)
    {
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
            {
                out_buffer_x[p] = rhs[c, k, j, i, m];
                out_buffer_x[p + 1] = rhs[c, k, j, i1, m];
                p = p + 2;
            }
        }
    }
    return 0;
}

}

}
