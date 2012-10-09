using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Forward;
using common.axis.ZAxis;
using sp.solve.shift.WriteBuffer;

namespace impl.sp.solve.shift.ZWriteBufferForward { 

public class IZWriteBufferForwardImpl<I, C, D, DIR> : BaseIZWriteBufferForwardImpl<I, C, D, DIR>, IWriteBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IForwardDirection
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
			
	stage = 0;
}

#pragma action
public void advance()
{
    stage++;
}

#pragma condition		
public bool finished()
{
	return stage >= ncells;
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

		out_buffer_solver[stage] = new double[22 * buffer_size];
	}
	
}

bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
	Buffer.Array = out_buffer_z = out_buffer_solver[stage];
			
	int i, j, k, n, isize, jsize, kend, c, m, p;
	
	c = slice[stage, 2];
	
	kend = 2 + cell_size[c, 2] - 1;
	
	isize = cell_size[c, 0] + 2;
	jsize = cell_size[c, 1] + 2;
			
    //---------------------------------------------------------------------
    //            create a running pointer for the send buffer  
    //---------------------------------------------------------------------
    p = 0;
    n = -1;
    for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
    {
        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
        {
            for (k = kend - 1; k <= kend; k++)
            {
                out_buffer_z[p    ] = lhs[c, k, j, i, n + 4];
                out_buffer_z[p + 1] = lhs[c, k, j, i, n + 5];
                for (m = 0; m <= 2; m++)
                {
                    out_buffer_z[p + 2 + m] = rhs[c, k, j, i, m];
                }
                p = p + 5;
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
                for (k = kend - 1; k <= kend; k++)
                {
                    out_buffer_z[p] = lhs[c, k, j, i, n + 4];
                    out_buffer_z[p + 1] = lhs[c, k, j, i, n + 5];
                    out_buffer_z[p + 2] = rhs[c, k, j, i, m];
                    p = p + 3;
                }
            }
        }
    }

	return 0;
}

}

}
