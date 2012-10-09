using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Backward;
using common.axis.ZAxis;
using sp.solve.shift.WriteBuffer;

namespace impl.sp.solve.shift.ZWriteBufferBackward { 

public class IZWriteBufferBackwardImpl<I, C, D, DIR> : BaseIZWriteBufferBackwardImpl<I, C, D, DIR>, IWriteBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IBackwardDirection
where DIR:IZ
{
private double[][] out_buffer_solver;
double[] out_buffer_z;

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
	int c, stage, isize, jsize, buffer_size;
	
	out_buffer_solver = new double[ncells][];			
	
	for (stage = 0; stage < ncells; stage++)
	{
		c = slice[stage, 2];
 		
 		isize = cell_size[c, 0] + 2;
 		jsize = cell_size[c, 1] + 2;
 		
 		buffer_size = (isize - start[c, 0] - end[c, 0]) * (jsize - start[c, 1] - end[c, 1]);						

		out_buffer_solver[stage] = new double[10 * buffer_size];
	}
	
}

bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
    Buffer.Array = out_buffer_z = out_buffer_solver[stage];
			
	int i, j, k, isize, jsize, k1, c, m, p, kstart;
	
	c = slice[stage, 2];
	
	kstart = 2;
	
	isize = cell_size[c, 0] + 2;
	jsize = cell_size[c, 1] + 2;
			
    k = kstart;
    k1 = kstart + 1;
    p = 0;
    for (m = 0; m <= 4; m++)
    {
        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
                out_buffer_z[p] = rhs[c, k, j, i, m];
                out_buffer_z[p + 1] = rhs[c, k1, j, i, m];
                p = p + 2;
            }
        }
    }

	return 0;
}

}

}
