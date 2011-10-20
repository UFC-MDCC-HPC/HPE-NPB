using System;
using br.ufc.pargo.hpe.backend.DGAC;
using br.ufc.pargo.hpe.basic;
using br.ufc.pargo.hpe.kinds;
using common.problem_size.Class;
using ft.FT;
using MPI;
using NPB3_0_JAV.BMInOut;

namespace impl.ft.FTImpl { 
	public class IFTImpl<C> : BaseIFTImpl<C>, IFT<C>
	where C:IClass{
	   
	    private Intracommunicator worldcomm;
	    private int node, np, /* np1, np2, layout_type,*/ root=0;
	    public static String BMName = "FT";
		protected int[,] dims;
		protected double[] twiddle;
			
		public override int go() 
		{ 
			runBenchMark();	
			return 0;
		}
				
        public void runBenchMark() 
		{
            worldcomm = this.WorldComm;
            np = worldcomm.Size;
            node = worldcomm.Rank;
            for(int i = 1; i <= T_max; i++) Timer.resetTimer(i);
                        
			dims = Blocks.dims;
			twiddle = Problem.twiddle;
            
            Fft.setParameters(1, u1, u0);

            Compute_index_map.go();            
            Compute_initial_conditions.go();
            Fftinit.go();            
            Fft.go();

            for(int i = 1; i <= T_max; i++) Timer.resetTimer(i);
            worldcomm.Barrier();

            Timer.start(T_total);

            Compute_index_map.go();
            Compute_initial_conditions.go();            
            Fftinit.go();
            Fft.go();

            double[] sums = new double[niter_default*2];
			
            Fft.setParameters(-1, u1, u2);
                
            for(int iter = 0; iter < niter; iter++) {
                Evolve.go();
                Fft.go();
                Checksum.setParameters(iter, sums);
                Checksum.go();
            }
            
            Verify.setParameters(niter, sums);
            int verified = Verify.go();
            Timer.stop(T_total);
            double total_time = Timer.readTimer(T_total);

            double ntotal_f = (double)(nx*ny*nz);
            double mflops=0.0;
            if(total_time != 0) {
                mflops = 0.000001*ntotal_f * (14.8157+7.19641*Math.Log(ntotal_f) +  (5.23518+7.21113*Math.Log(ntotal_f))*niter)/total_time;
            }
            else {
                mflops = 0.0;
            }
            if(node == 0) {
                BMResults results = new BMResults(BMName,
                                        problem_class.ToString()[0],
                                        nx,
                                        ny,
                                        nz,
                                        niter,
                                        total_time,
                                        mflops,
                                        "floating point",
                                        verified,
                                        true,
                                        np,
                                        -1);
                results.print();
            }
        }		
	}
}
