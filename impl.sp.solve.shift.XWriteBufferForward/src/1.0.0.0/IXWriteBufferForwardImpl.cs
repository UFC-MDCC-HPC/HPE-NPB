using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Forward;
using common.axis.XAxis;
using sp.solve.shift.WriteBuffer;

namespace impl.sp.solve.shift.XWriteBufferForward { 

public class IXWriteBufferForwardImpl<I, C, D, DIR> : BaseIXWriteBufferForwardImpl<I, C, D, DIR>, IWriteBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IForwardDirection
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
	int c, stage, jsize, ksize, buffer_size;
	
	out_buffer_solver = new double[ncells][];			
	
	for (stage = 0; stage < ncells; stage++)
    {
        c = slice[stage, 0];
		
        jsize = cell_size[c, 1] + 2;
        ksize = cell_size[c, 2] + 2;
		
		buffer_size = (jsize - start[c, 1] - end[c, 1]) * (ksize - start[c, 2] - end[c, 2]);
		
        out_buffer_solver[stage] = new double[22*buffer_size];				
	}			
}
		
bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
    Buffer.Array = out_buffer_x = out_buffer_solver[stage];
			
	int c, i, j, k, n, iend, jsize, ksize, m, p;
	
    c = slice[stage, 0];

    iend = 2 + cell_size[c, 0] - 1;

    jsize = cell_size[c, 1] + 2;
    ksize = cell_size[c, 2] + 2;
						
    //--------------------------------------------------------------------
    //            create a running pointer for the send buffer  
    //---------------------------------------------------------------------
    p = 0;
    n = -1;
    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
    {
        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
        {
            for (i = iend - 1; i <= iend; i++)
            {
                out_buffer_x[p] = lhs[c, k, j, i, n + 4];  //if (this.Rank == 1) Console.WriteLine("out_buffer[{0}]={1}", p, out_buffer_x[p]);
                out_buffer_x[p + 1] = lhs[c, k, j, i, n + 5]; //if (this.Rank == 1) Console.WriteLine("out_buffer[{0}]={1}", p+1, out_buffer_x[p+1]);
                for (m = 0; m <= 2; m++)
                {
                    out_buffer_x[p + 2 + m] = rhs[c, k, j, i, m]; //if (this.Rank == 1) Console.WriteLine("out_buffer[{0}]={1}", p+2+m, out_buffer_x[p+2+m]);
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
            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
            {
                for (i = iend - 1; i <= iend; i++)
                {
                    out_buffer_x[p] = lhs[c, k, j, i, n + 4];      //if (this.Rank == 1) Console.WriteLine("out_buffer[{0}]={1}", p, out_buffer_x[p]);
                    out_buffer_x[p + 1] = lhs[c, k, j, i, n + 5];  //if (this.Rank == 1) Console.WriteLine("out_buffer[{0}]={1}", p+1, out_buffer_x[p+1]);
                    out_buffer_x[p + 2] = rhs[c, k, j, i, m];      //if (this.Rank == 1) Console.WriteLine("out_buffer[{0}]={1}", p+2, out_buffer_x[p+2]);
                    p = p + 3;
                }
            }
        }
    }
			
	return 0;
}

}

}
