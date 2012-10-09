using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Forward;
using common.axis.YAxis;
using sp.solve.shift.WriteBuffer;

namespace impl.sp.solve.shift.YWriteBufferForward { 

public class IYWriteBufferForwardImpl<I, C, D, DIR> : BaseIYWriteBufferForwardImpl<I, C, D, DIR>, IWriteBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IForwardDirection
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
	int c, isize, ksize, buffer_size;
	
	out_buffer_solver = new double[ncells][];
	
	for (int stage = 0; stage < ncells; stage++)
	{
		c = slice[stage, 1];
		
		isize = cell_size[c, 0] + 2;
		ksize = cell_size[c, 2] + 2;
						
		buffer_size = (isize - start[c, 0] - end[c, 0]) * (ksize - start[c, 2] - end[c, 2]);

		out_buffer_solver[stage] = new double[22 * buffer_size];
						
	}
	
}

bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
    Buffer.Array = out_buffer_y = out_buffer_solver[stage];
			
	int i, j, k, n, isize, jend, ksize, j1, c, m, p, jstart; /* requests(2), statuses(MPI_STATUS_SIZE, 2);*/
	double r1, r2, d, e, sm1, sm2;
	double[] s = new double[5];
	
	c = slice[stage, 1];
	
	jstart = 2;
	jend = 2 + cell_size[c, 1] - 1;
	
	isize = cell_size[c, 0] + 2;
	ksize = cell_size[c, 2] + 2;
	
	//---------------------------------------------------------------------
    //            create a running pointer for the send buffer  
    //---------------------------------------------------------------------
    p = 0;
    n = -1;
    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
    {
        for (i = start[c, 0]; i < isize - end[c, 0]; i++)
        {
            for (j = jend - 1; j <= jend; j++)
            {
                out_buffer_y[p    ] = lhs[c, k, j, i, n + 4];
                out_buffer_y[p + 1] = lhs[c, k, j, i, n + 5];
           //     Console.WriteLine("out_buffer_y["+ p + "] = " + out_buffer_y[p]);
            //    Console.WriteLine("out_buffer_y["+ (p+1) + "] = " + out_buffer_y[p+1]);
                for (m = 0; m <= 2; m++)
                {
                    out_buffer_y[p + 2 + m] = rhs[c, k, j, i, m];
             //       Console.WriteLine("out_buffer_y["+ (p+2+m) + "] = " + out_buffer_y[p+2+m]);
                }
                p = p + 5;
            }
        }
    }

    for (m = 3; m <= 4; m++)
    {
        n = (m - 2) * 5 - 1;
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
                for (j = jend - 1; j <= jend; j++)
                {
                    out_buffer_y[p] = lhs[c, k, j, i, n + 4];
                    out_buffer_y[p + 1] = lhs[c, k, j, i, n + 5];
                    out_buffer_y[p + 2] = rhs[c, k, j, i, m];
               //     Console.WriteLine("out_buffer_y["+ (p) + "] = " + out_buffer_y[p]);
                //    Console.WriteLine("out_buffer_y["+ (p+1) + "] = " + out_buffer_y[p+1]);
               //     Console.WriteLine("out_buffer_y["+ (p+2) + "] = " + out_buffer_y[p+2]);
                    p = p + 3;
                }
            }
        }
    }

	return 0;
}

}

}
