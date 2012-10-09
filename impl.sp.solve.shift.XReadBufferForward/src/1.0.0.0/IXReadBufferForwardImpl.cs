using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using sp.problem_size.Instance_SP;
using common.problem_size.Class;
using common.direction.Forward;
using common.axis.XAxis;
using sp.solve.shift.ReadBuffer;

namespace impl.sp.solve.shift.XReadBufferForward { 

public class IXReadBufferForwardImpl<I, C, D, DIR> : BaseIXReadBufferForwardImpl<I, C, D, DIR>, IReadBuffer<I, C, D, DIR>
where I:IInstance_SP<C>
where C:IClass
where D:IForwardDirection
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
	stage = 0;
    Buffer.Array = in_buffer_x = in_buffer_solver[stage];
}

#pragma action
public void advance()
{
    stage++;
    Buffer.Array = in_buffer_x = in_buffer_solver[stage];
}
		
#pragma condition		
public bool finished()
{
	return stage >= ncells;
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
		
        in_buffer_solver[stage] = new double[22*buffer_size];
	}			
}
		
bool buffer_created = false;		
		
#pragma action
public override int go() 
{ 
			
    int c, i, j, k, n, jsize, ksize, i1, m, p, istart;
    double r1, r2, d, e;

    c = slice[stage, 0];

    istart = 2;

    jsize = cell_size[c, 1] + 2;
    ksize = cell_size[c, 2] + 2;
	
    double[] s = new double[5];
			
    //---------------------------------------------------------------------
    //            unpack the buffer                                 
    //---------------------------------------------------------------------
    i = istart;
    i1 = istart + 1;
    n = -1;

    //---------------------------------------------------------------------
    //            create a running pointer
    //---------------------------------------------------------------------

    p = 0;
    for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
    {
        for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
        {
            lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                     in_buffer_x[p    ] * lhs[c, k, j, i, n + 1];   //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p, in_buffer_x[p]);
            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                     in_buffer_x[p + 1] * lhs[c, k, j, i, n + 1];   //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+1, in_buffer_x[p+1]);
            for (m = 0; m <= 2; m++)
            {		                       
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                      in_buffer_x[p + 2 + m] * lhs[c, k, j, i, n + 1]; //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+2+m, in_buffer_x[p+2+m]);
            }
            
            d = in_buffer_x[p + 5]; //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+5, in_buffer_x[p+5]);
            e = in_buffer_x[p + 6]; //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+6, in_buffer_x[p+6]);
            for (m = 0; m <= 2; m++)
            {
                s[m] = in_buffer_x[p + 7 + m]; //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+7+m, in_buffer_x[p+7+m]);
            }
            r1 = lhs[c, k, j, i, n + 2];
            lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
            lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
            for (m = 0; m <= 2; m++)
            {
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
            }
            r2 = lhs[c, k, j, i1, n + 1];
            lhs[c, k, j, i1, n + 2] = lhs[c, k, j, i1, n + 2] - d * r2;
            lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] - e * r2;
            for (m = 0; m <= 2; m++)
            {
                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] - s[m] * r2;
            }
            p = p + 10;
        }
    }

    for (m = 3; m <= 4; m++)
    {
        n = (m - 2) * 5 - 1;
        for (k = start[c, 2]; k < ksize - end[c, 2]; k++)
        {
            for (j = start[c, 1]; j < jsize - end[c, 1]; j++)
            {
                lhs[c, k, j, i, n + 2] = lhs[c, k, j, i, n + 2] -
                         in_buffer_x[p] * lhs[c, k, j, i, n + 1]; //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p, in_buffer_x[p]);
                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] -
                         in_buffer_x[p + 1] * lhs[c, k, j, i, n + 2]; //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+1, in_buffer_x[p+1]);
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] -
                         in_buffer_x[p + 2] * lhs[c, k, j, i, n + 1]; //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+2, in_buffer_x[p+2]);
                d = in_buffer_x[p + 3];  //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+3, in_buffer_x[p+3]);
                e = in_buffer_x[p + 4];  //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+4, in_buffer_x[p+4]);
                s[m] = in_buffer_x[p + 5];  //if (this.Rank == 0) Console.WriteLine("in_buffer[{0}]={1}", p+5, in_buffer_x[p+5]);
                r1 = lhs[c, k, j, i, n + 2];
                lhs[c, k, j, i, n + 3] = lhs[c, k, j, i, n + 3] - d * r1;
                lhs[c, k, j, i, n + 4] = lhs[c, k, j, i, n + 4] - e * r1;
                rhs[c, k, j, i, m] = rhs[c, k, j, i, m] - s[m] * r1;
                r2 = lhs[c, k, j, i1, n + 1];
                lhs[c, k, j, i1, n + 2] = lhs[c, k, j, i1, n + 2] - d * r2;
                lhs[c, k, j, i1, n + 3] = lhs[c, k, j, i1, n + 3] - e * r2;
                rhs[c, k, j, i1, m] = rhs[c, k, j, i1, m] - s[m] * r2;
                p = p + 6;
            }
        }
    }

	return 0;
}

}

}
