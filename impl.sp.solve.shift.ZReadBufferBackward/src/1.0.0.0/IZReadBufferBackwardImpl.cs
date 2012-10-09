using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Backward;
using common.axis.ZAxis;
using sp.solve.shift.ReadBuffer;

namespace impl.sp.solve.shift.ZReadBufferBackward { 

public class IZReadBufferBackwardImpl<I, C, D, DIR> : BaseIZReadBufferBackwardImpl<I, C, D, DIR>, IReadBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IBackwardDirection
where DIR:IZ
{

private double[][] in_buffer_solver;
double[] in_buffer_z;

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
    Buffer.Array = in_buffer_z = in_buffer_solver[stage];
}

#pragma action
public void advance()
{
    stage--;
    Buffer.Array = in_buffer_z = in_buffer_solver[stage];
}
		
#pragma condition		
public bool finished()
{
	return stage < 0;
}
		
void create_buffers() 		
{
	int c, stage, isize, jsize, buffer_size;
	
	in_buffer_solver = new double[ncells][];			
	
	for (stage = 0; stage < ncells; stage++)
	{
		c = slice[stage, 2];
 		
 		isize = cell_size[c, 0] + 2;
 		jsize = cell_size[c, 1] + 2;
 		
 		buffer_size = (isize - start[c, 0] - end[c, 0]) * (jsize - start[c, 1] - end[c, 1]);						

		in_buffer_solver[stage] = new double[10 * buffer_size];						
	}
	
}

bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 

	int i, j, k, n, isize, jsize, kend, k1, c, m, p;
	double sm1, sm2;
	
	c = slice[stage, 2];
	
	kend = 2 + cell_size[c, 2] - 1;
	
	isize = cell_size[c, 0] + 2;
	jsize = cell_size[c, 1] + 2;

	//---------------------------------------------------------------------
    //            unpack the buffer for the first three factors         
    //---------------------------------------------------------------------
    n = -1;
    p = 0;
    k = kend;
    k1 = k - 1;
    for (m = 0; m <= 2; m++)
    {
        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
         //       Console.WriteLine("in_buffer_z["+ p + "] = " + in_buffer_z[p]);
         //       Console.WriteLine("in_buffer_z["+ (p+1) + "] = " + in_buffer_z[p+1]);
                sm1 = in_buffer_z[p];
                sm2 = in_buffer_z[p + 1];
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                       lhs[c, k, j, i, n + 4] * sm1 -
                       lhs[c, k, j, i, n + 5] * sm2;
                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                       lhs[c, k1, j, i, n + 4] * rhs[c, k, j, i, m] -
                       lhs[c, k1, j , i, n + 5] * sm1;
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
        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
        {
            for (i = start[c, 0]; i < isize - end[c, 0]; i++)
            {
//                Console.WriteLine("in_buffer_z["+ p + "] = " + in_buffer_z[p]);
 //               Console.WriteLine("in_buffer_z["+ (p+1) + "] = " + in_buffer_z[p+1]);
                sm1 = in_buffer_z[p];
                sm2 = in_buffer_z[p + 1];
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                       lhs[c, k, j, i, n + 4] * sm1 -
                       lhs[c, k, j, i, n + 5] * sm2;
                rhs[c, k1, j, i, m] = rhs[c, k1, j, i, m] -
                       lhs[c, k1, j, i, n + 4] * rhs[c, k, j, i, m] -
                       lhs[c, k1, j , i, n + 5] * sm1;
                p = p + 2;
            }
        }
    }
    return 0;
}

}

}
