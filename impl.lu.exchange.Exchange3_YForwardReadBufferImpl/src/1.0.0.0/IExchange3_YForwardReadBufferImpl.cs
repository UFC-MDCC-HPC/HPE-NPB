using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.direction.Forward;
using common.axis.YAxis;
using lu.exchange.Exchange3_ReadBuffer;

namespace impl.lu.exchange.Exchange3_YForwardReadBufferImpl { 

	public class IExchange3_YForwardReadBufferImpl<DIR, O> : BaseIExchange3_YForwardReadBufferImpl<DIR, O>, IExchange3_ReadBuffer<DIR, O>
		where DIR:IForwardDirection
		where O:IY
	{
	
		public override int go() 
		{ 			
            int j, k0;
            int ipos1, ipos2;
						
            for(k0 = 1; k0<=nz; k0++) 
			{
                for(j = 1; j<=ny; j++) 
				{
                    ipos1 = (k0-1)*ny + j - 1;
                    ipos2 = ipos1 + ny*nz;
                    g[k0-1, j+1, 0, 0] = buf1[0*size2+ipos1];
                    g[k0-1, j+1, 0, 1] = buf1[1*size2+ipos1];
                    g[k0-1, j+1, 0, 2] = buf1[2*size2+ipos1];
                    g[k0-1, j+1, 0, 3] = buf1[3*size2+ipos1];
                    g[k0-1, j+1, 0, 4] = buf1[4*size2+ipos1];

                    g[k0-1, j+1, 1, 0] = buf1[0*size2+ipos2];
                    g[k0-1, j+1, 1, 1] = buf1[1*size2+ipos2];
                    g[k0-1, j+1, 1, 2] = buf1[2*size2+ipos2];
                    g[k0-1, j+1, 1, 3] = buf1[3*size2+ipos2];
                    g[k0-1, j+1, 1, 4] = buf1[4*size2+ipos2];
					
                }
            }
		
		    return 0;
		} // end activate method 

		private double[,,,] g;
		private int  nx, ny, nz;
		
        private int size2;
        private double[] buf1;
		
		public void setParameters(double[,,,] g, int nx, int ny, int nz)
		{
			this.g = g;
			this.nx = nx;
			this.ny = ny;
			this.nz = nz;
			
            int bsize = 10*ny*nz;
			
            size2 = bsize / 5;
            //Console.Error.WriteLine(this.GlobalRank + ": *** BUFFER SIZE IS " + bsize + " ny=" + ny + ", nz=" + nz);
            buf1 = Buffer.Array = new double[bsize];
		}
	}

}
