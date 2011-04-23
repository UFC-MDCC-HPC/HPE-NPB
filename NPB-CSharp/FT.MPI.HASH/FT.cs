/*
-------------------------------------------------------------------------
        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         
                                   F T                                   
-------------------------------------------------------------------------
    This benchmark is part of the NAS Parallel Benchmark 3.3 suite.      
    It is described in NAS Technical Reports 95-020 and 02-007           
                                                                         
    Permission to use, copy, distribute and modify this software         
    for any purpose with or without fee is hereby granted.  We           
    request, however, that all derived work reference the NAS            
    Parallel Benchmarks 3.3. This software is provided "as is"           
    without express or implied warranty.                                 
                                                                         
    Information on NPB 3.3, including the technical report, the          
    original specifications, source code, results and information        
    on how to submit new results, is available at:                       
                                                                         
           http://www.nas.nasa.gov/Software/NPB/                         
                                                                         
    Send comments or suggestions to  npb@nas.nasa.gov                    
                                                                         
          NAS Parallel Benchmarks Group                                  
          NASA Ames Research Center                                      
          Mail Stop: T27A-1                                              
          Moffett Field, CA   94035-1000                                 
                                                                         
          E-mail:  npb@nas.nasa.gov                                      
          Fax:     (650) 604-3957                                        
-------------------------------------------------------------------------
TO REDUCE THE AMOUNT OF MEMORY REQUIRED BY THE BENCHMARK WE NO LONGER
STORE THE ENTIRE TIME EVOLUTION ARRAY "EX" FOR ALL TIME STEPS, BUT
JUST FOR THE FIRST. ALSO, IT IS STORED ONLY FOR THE PART OF THE GRID
FOR WHICH THE CALLING PROCESSOR IS RESPONSIBLE, SO THAT THE MEMORY 
USAGE BECOMES SCALABLE. THIS NEW ARRAY IS CALLED "TWIDDLE" (SEE
NPB3.0-SER)

TO AVOID PROBLEMS WITH VERY LARGE ARRAY SIZES THAT ARE COMPUTED BY
MULTIPLYING GRID DIMENSIONS (CAUSING INTEGER OVERFLOW IN THE VARIABLE
NTOTAL) AND SUBSEQUENTLY DIVIDING BY THE NUMBER OF PROCESSORS, WE
COMPUTE THE SIZE OF ARRAY PARTITIONS MORE CONSERVATIVELY AS
((NX*NY)/NP)*NZ, WHERE NX, NY, AND NZ ARE GRID DIMENSIONS AND NP IS
THE NUMBER OF PROCESSORS, THE RESULT IS STORED IN "NTDIVNP". FOR THE 
PERFORMANCE CALCULATION WE STORE THE TOTAL NUMBER OF GRID POINTS IN A 
FLOATING POINT NUMBER "NTOTAL_F" INSTEAD OF AN INTEGER.
THIS FIX WILL FAIL IF THE NUMBER OF PROCESSORS IS SMALL.

UGLY HACK OF SUBROUTINE IPOW46: FOR VERY LARGE GRIDS THE SINGLE EXPONENT
FROM NPB2.3 MAY NOT FIT IN A 32-BIT INTEGER. HOWEVER, WE KNOW THAT THE
"EXPONENT" ARGUMENT OF THIS ROUTINE CAN ALWAYS BE FACTORED INTO A TERM 
DIVISIBLE BY NX (EXP_1) AND ANOTHER TERM (EXP_2). NX IS USUALLY A POWER
OF TWO, SO WE CAN KEEP HALVING IT UNTIL THE PRODUCT OF EXP_1
AND EXP_2 IS SMALL ENOUGH (NAMELY EXP_2 ITSELF). THIS UPDATED VERSION
OF IPWO46, WHICH NOW TAKES THE TWO FACTORS OF "EXPONENT" AS SEPARATE
ARGUMENTS, MAY BREAK DOWN IF EXP_1 DOES NOT CONTAIN A LARGE POWER OF TWO.
-------------------------------------------------------------------------
 Authors: D. Bailey
          W. Saphir
          R. F. Van der Wijngaart
 Translation to C# and MPI.NET Code
          Cenez Araújo de Rezende, MDCC/UFC
          Francisco Heron de Carvalho Junior (MDCC/UFC)
-------------------------------------------------------------------------
*/
using System;
using System.IO;

namespace NPB {
    public class FT: FTBase {

        public FT(char c) : base(c) {
        }

        static void Main(String[] argv) {

            FT ft = null;
            bool debug = false;

            try {
                string param = argv[0];
            }
            catch(Exception) {
                argv = new String[1];
                argv[0] = "CLASS=S"; // CLASS DEFAULT, IF USER NOT TYPE CLASS=S IN COMMAND-LINE ARGUMENT
            }
            char paramClass;
            if(!debug) {
                IO.parseCmdLineArgs(argv);
                paramClass = IO.CLASS;
            }
            else {
                paramClass = 'S';  //DEBUG: CHANGE TO [K=(S and 4 PROCESSORS)] OR [S=(S and 1 PROCESSOR)]
            }                      //DEBUG: OR [T=(A and 4 PROCESSORS)] OR [I=(B and 4 PROCESSORS)]

            try {
                ft = new FT(paramClass);
            }
            catch(OutOfMemoryException e) {
                Console.WriteLine(e.ToString());
                Environment.Exit(0);
            }

            ft.runBenchMark();

        }

        public void runBenchMark() {
            for(int i = 1; i <= T_max; i++) timer.resetTimer(i);
            initialConfig();
            problemDefination();
            blocksInfo();
            if(timers_enabled)
                synchup();
            compute_indexmap(twiddle);
            compute_initial_conditions(u1);
            fft_init(dims[0, 0]);  //control u
            fft(1, u1, u0); // fft(1, u1, u0);
            //c---------------------------------------------------------------------
            //c Start over from the beginning. Note that all operations must
            //c be timed, in contrast to other benchmarks. 
            //c---------------------------------------------------------------------
            for(int i = 1; i <= T_max; i++) timer.resetTimer(i);
            worldcomm.Barrier();

            timer.start(T_total);
            if(timers_enabled) timer.start(T_setup);

            compute_indexmap(twiddle);
            compute_initial_conditions(u1);
            fft_init(dims[0, 0]);

            if(timers_enabled) synchup();
            if(timers_enabled) timer.stop(T_setup);

            if(timers_enabled) timer.start(T_fft);
            fft(1, u1, u0);
            if(timers_enabled) timer.stop(T_fft);

            double[] sums = new double[niter_default*2];
            for(int iter = 0; iter < niter; iter++) {
                if(timers_enabled) timer.start(T_evolve);
                evolve(u0, u1, twiddle, dims[0, 0], dims[1, 0], dims[2, 0]);
                if(timers_enabled) timer.stop(T_evolve);
                if(timers_enabled) timer.start(T_fft);
                fft(-1, u1, u2);
                if(timers_enabled) timer.stop(T_fft);
                if(timers_enabled) synchup();
                if(timers_enabled) timer.start(T_checksum);
                checksum(iter, sums, u2, dims[0, 0], dims[1, 0], dims[2, 0]);
                if(timers_enabled) timer.stop(T_checksum);
            }

            int verified = verify(nx, ny, nz, niter, sums);
            timer.stop(T_total);
            double total_time = timer.readTimer(T_total); //total_time = timer_read(t_total);

            double ntotal_f = (double)(nx*ny*nz);
            double mflops=0.0;
            if(total_time != 0) {
                mflops = 0.000001*ntotal_f * (14.8157+7.19641*Math.Log(ntotal_f) +  (5.23518+7.21113*Math.Log(ntotal_f))*niter)/total_time;
            }
            else {
                mflops = 0.0;
            }
            if(node == 0) {
                IO.print_results(BMName, CLSS, nx, ny, nz, niter, np, np, total_time, mflops, "floating point", verified, "3.3");
            }
            if(timers_enabled)
                print_timers();
            mpi.Dispose();
        }

        public void initialConfig() {
            int fstatus=0;
            if(node == 0) {
                Console.WriteLine(" NAS Parallel Benchmarks "+ "3.3" +" -- FT Benchmark ");
                try {
                    Console.Write("Trying Read from input file inputft.data: ");
                    int[] conf = { 1, 1, 2 };
                    string[] vetTemp = IO.readFileData("inputft.data", conf);
                    niter = int.Parse(vetTemp[0]);
                    layout_type = int.Parse(vetTemp[1]);
                    np1 = int.Parse(vetTemp[2]);
                    np2 = int.Parse(vetTemp[3]);
                }
                catch /*(System.IO.FileNotFoundException e)*/ {
                    Console.WriteLine("inputft.data not found");
                    fstatus = 1;
                    //Console.WriteLine(e.ToString());
                }

                if(fstatus == 0) {
                    Console.WriteLine("inputft.data found");

                    //c---------------------------------------------------------------------
                    //c check to make sure input data is consistent
                    //c---------------------------------------------------------------------
                    //c---------------------------------------------------------------------
                    //c 1. product of processor grid dims must equal number of processors
                    //c---------------------------------------------------------------------

                    if(np1 * np2 != np) {
                        Console.WriteLine(" np1 and np2 given in input file are not valid.");
                        Console.WriteLine("Product is "+ np1*np2+" and should be "+np);
                        mpi.Dispose();
                        Environment.Exit(0);
                    }

                    //c---------------------------------------------------------------------
                    //c 2. layout type must be valid
                    //c---------------------------------------------------------------------

                    if(layout_type != layout_0D && layout_type != layout_1D && layout_type != layout_2D) {
                        Console.WriteLine(" Layout type specified in inputft.data is invalid ");
                        mpi.Dispose();
                        Environment.Exit(0);
                    }

                    //c---------------------------------------------------------------------
                    //c 3. 0D layout must be 1x1 grid
                    //c---------------------------------------------------------------------

                    if(layout_type == layout_0D && (np1 != 1 || np2 != 1)) {
                        Console.WriteLine(" For 0D layout, both np1 and np2 must be 1 ");
                        mpi.Dispose();
                        Environment.Exit(0);
                    }
                    //c---------------------------------------------------------------------
                    //c 4. 1D layout must be 1xN grid
                    //c---------------------------------------------------------------------

                    if(layout_type == layout_1D && np1 != 1) {
                        Console.WriteLine(" For 1D layout, np1 must be 1 ");
                        mpi.Dispose();
                        Environment.Exit(0);
                    }
                }
                else {
                    Console.WriteLine(" No input file inputft.data. Using compiled defaults");
                    niter = niter_default;
                    if(np == 1) {
                        np1 = 1;
                        np2 = 1;
                        layout_type = layout_0D;
                    }
                    else if(np <= nz) {
                        np1 = 1;
                        np2 = np;
                        layout_type = layout_1D;
                    }
                    else {
                        np1 = nz;
                        np2 = np/nz;
                        layout_type = layout_2D;
                    }
                }

                Console.WriteLine(" Size: " + nx + "x" + ny + "x" + nz);
                Console.WriteLine(" Iterations: "+ niter);
                Console.WriteLine(" Number of processes : "+ np);
                Console.WriteLine(" Processor array     : "+np1+"x"+np2);
                if(layout_type == layout_0D) {
                    Console.WriteLine(" Layout type: OD");
                }
                else if(layout_type == layout_1D) {
                    Console.WriteLine(" Layout type: 1D");
                }
                else {
                    Console.WriteLine(" Layout type: 2D");
                }
            }
            //c---------------------------------------------------------------------
            //c Since np1, np2 and layout_type are in a common block, 
            //c this sends all three. 
            //c---------------------------------------------------------------------

            worldcomm.Broadcast<int>(ref np1, root);   //call MPI_BCAST(np1, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            worldcomm.Broadcast<int>(ref niter, root); //call MPI_BCAST(niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            worldcomm.Broadcast<int>(ref np2, root);
            //worldcomm.Broadcast<int>(ref layout_type, root);
            if(np1 == 1 && np2 == 1) {
                layout_type = layout_0D;
            }
            else if(np1 == 1) {
                layout_type = layout_1D;
            }
            else {
                layout_type = layout_2D;
            }
        }

        public void problemDefination() {
            if(layout_type == layout_0D) {
                for(int i = 0; i < 3; i++) {
                    dims[0, i] = nx;
                    dims[1, i] = ny;
                    dims[2, i] = nz;
                }
            }
            else if(layout_type == layout_1D) {
                dims[0, 0] = nx;
                dims[1, 0] = ny;
                dims[2, 0] = nz;

                dims[0, 1] = nx;
                dims[1, 1] = ny;
                dims[2, 1] = nz;

                dims[0, 2] = nz;
                dims[1, 2] = nx;
                dims[2, 2] = ny;
            }
            else if(layout_type == layout_2D) {
                dims[0, 0] = nx;
                dims[1, 0] = ny;
                dims[2, 0] = nz;

                dims[0, 1] = ny;
                dims[1, 1] = nx;
                dims[2, 1] = nz;

                dims[0, 2] = nz;
                dims[1, 2] = nx;
                dims[2, 2] = ny;

            }
            //for(i = 0; i < 3; i++) {
            //    dims[1, i] = dims[1, i] / np1;
            //    dims[2, i] = dims[2, i] / np2;
            //}
            dims[1, 0] = dims[1, 0] / np1;
            dims[2, 0] = dims[2, 0] / np2;
            dims[1, 1] = dims[1, 1] / np1;
            dims[2, 1] = dims[2, 1] / np2;
            dims[1, 2] = dims[1, 2] / np1;
            dims[2, 2] = dims[2, 2] / np2;

            //complex u0(ntdivnp), u1(ntdivnp), u2(ntdivnp)
            u1 = new double[dims[1, 0], dims[2, 0], dims[0, 0], 2];//u1 = new double[dims[2, 0], dims[1, 0], dims[0, 0], 2];
            u0 = new double[dims[1, 0], dims[2, 0], dims[0, 0], 2];
            u2 = new double[dims[1, 0], dims[2, 0], dims[0, 0], 2];
            u = new double[nx, 2];
            // twiddle(ntdivnp)
            twiddle = new double[ntdivnp];
        }

        public void blocksInfo() {
            //c---------------------------------------------------------------------
            //c Determine processor coordinates of this processor
            //c Processor grid is np1xnp2. 
            //c Arrays are always (n1, n2/np1, n3/np2)
            //c Processor coords are zero-based. 
            //c---------------------------------------------------------------------

            me2 = (int)mod(node, np2);  // goes from 0...np2-1
            me1 = node/np2;        // goes from 0...np1-1

            //c---------------------------------------------------------------------
            //c Communicators for rows/columns of processor grid. 
            //c commslice1 is communicator of all procs with same me1, ranked as me2
            //c commslice2 is communicator of all procs with same me2, ranked as me1
            //c mpi_comm_split(comm, color, key, ...)
            //c---------------------------------------------------------------------

            commslice1=(MPI.Intracommunicator)worldcomm.Split(me1, me2);//MPI_Comm_split(MPI_COMM_WORLD,me1,me2,commslice1,ierr)
            commslice2=(MPI.Intracommunicator)worldcomm.Split(me2, me1);//MPI_Comm_split(MPI_COMM_WORLD,me2,me1,commslice2,ierr)

            //if(timers_enabled)
            //    synchup();

            //c---------------------------------------------------------------------
            //c Determine which section of the grid is owned by this
            //c processor. 
            //c---------------------------------------------------------------------
            if(layout_type == layout_0D) {

                for(int i = 0; i < 3; i++) {
                    xstart[i] = 1;
                    xend[i]   = nx;
                    ystart[i] = 1;
                    yend[i]   = ny;
                    zstart[i] = 1;
                    zend[i]   = nz;
                }

            }
            else if(layout_type == layout_1D) {
                xstart[0] = 1;
                xend[0]   = nx;
                ystart[0] = 1;
                yend[0]   = ny;
                zstart[0] = 1 + me2 * nz/np2;
                zend[0]   = (me2+1) * nz/np2;

                xstart[1] = 1;
                xend[1]   = nx;
                ystart[1] = 1;
                yend[1]   = ny;
                zstart[1] = 1 + me2 * nz/np2;
                zend[1]   = (me2+1) * nz/np2;

                xstart[2] = 1;
                xend[2]   = nx;
                ystart[2] = 1 + me2 * ny/np2;
                yend[2]   = (me2+1) * ny/np2;
                zstart[2] = 1;
                zend[2] = nz;

            }
            else if(layout_type == layout_2D) {

                xstart[0] = 1;
                xend[0]   = nx;
                ystart[0] = 1 + me1 * ny/np1;
                yend[0]   = (me1+1) * ny/np1;
                zstart[0] = 1 + me2 * nz/np2;
                zend[0]   = (me2+1) * nz/np2;

                xstart[1] = 1 + me1 * nx/np1;
                xend[1]   = (me1+1)*nx/np1;
                ystart[1] = 1;
                yend[1]   = ny;
                zstart[1] = zstart[0];
                zend[1]   = zend[0];

                xstart[2] = xstart[1];
                xend[2]   = xend[1];
                ystart[2] = 1 + me2 *ny/np2;
                yend[2]   = (me2+1)*ny/np2;
                zstart[2] = 1;
                zend[2] = nz;
            }

            //c---------------------------------------------------------------------
            //c Set up info for blocking of ffts and transposes.  This improves
            //c performance on cache-based systems. Blocking involves
            //c working on a chunk of the problem at a time, taking chunks
            //c along the first, second, or third dimension. 
            //c
            //c - In cffts1 blocking is on 2nd dimension (with fft on 1st dim)
            //c - In cffts2/3 blocking is on 1st dimension (with fft on 2nd and 3rd dims)

            //c Since 1st dim is always in processor, we'll assume it's long enough 
            //c (default blocking factor is 16 so min size for 1st dim is 16)
            //c The only case we have to worry about is cffts1 in a 2d decomposition. 
            //c so the blocking factor should not be larger than the 2nd dimension. 
            //c---------------------------------------------------------------------

            fftblock = fftblock_default;
            fftblockpad = fftblockpad_default;

            int dim1 = ny/np1;
            int dim2 = nx/np1;
            int dim3 = nx/np1;
            //if(layout_type == layout_2D) {
            //    if(dims[1, 0] < fftblock)
            //        fftblock = dims[1, 0];
            //    if(dims[1, 1] < fftblock)
            //        fftblock = dims[1, 1];
            //    if(dims[1, 2] < fftblock)
            //        fftblock = dims[1, 2];
            //}
            if(layout_type == layout_2D) {
                if(dim1 < fftblock)
                    fftblock = dim1;
                if(dim2 < fftblock)
                    fftblock = dim2;
                if(dim3 < fftblock)
                    fftblock = dim3;
            }

            if(fftblock != fftblock_default)
                fftblockpad = fftblock + 3;
            size1 = ((int)(nz/np2))*nx*2;
            size2 = nx*2;
        }

        public void synchup() {
            timer.start(T_synch);
            worldcomm.Barrier();
            timer.stop(T_synch);
        }

        public void compute_indexmap(double[] twiddle) {
            int i, j, k, ii, ii2, jj, ij2, kk;
            double ap;
            double alpha=.000001, pi = Math.PI;
            //Fortran: twiddle(d1, d2, d3)
            //C#     : twiddle[d3, d2, d1]

            //c---------------------------------------------------------------------
            //c this function is very different depending on whether 
            //c we are in the 0d, 1d or 2d layout. Compute separately. 
            //c basically we want to convert the fortran indices 
            //c   1 2 3 4 5 6 7 8 
            //c to 
            //c   0 1 2 3 -4 -3 -2 -1
            //c The following magic formula does the trick:
            //c mod(i-1+n/2, n) - n/2
            //c---------------------------------------------------------------------
            int d1 = dims[0, 2];
            int d2 = dims[1, 2];
            int d3 = dims[2, 2];

            ap = -4.0 * alpha * pi * pi;

            int idx;
            if(layout_type == layout_0D) { //xyz layout
                for(i = 1; i <= d1; i++) {
                    ii =  (int)mod(i+xstart[2]-2+nx/2, nx) - nx/2;
                    ii2 = ii*ii;
                    for(j = 1; j <= d2; j++) {
                        jj = (int)mod(j+ystart[2]-2+ny/2, ny) - ny/2;
                        ij2 = jj * jj + ii2;
                        for(k = 1; k <= d3; k++) {
                            kk = (int)mod(k+zstart[2]-2+nz/2, nz) - nz/2;
                            idx = (((k-1)*d2+(j-1))*d1+(i-1));//twiddle[k, j, i]
                            twiddle[idx] = Math.Exp(ap * (double)(kk * kk + ij2));
                        }
                    }
                }
            }
            else if(layout_type == layout_1D) { // zxy layout 
                for(i = 1; i <= d2; i++) {
                    ii =  (int)mod(i+xstart[2]-2+nx/2, nx) - nx/2;
                    ii2 = ii*ii;
                    for(j = 1; j <= d3; j++) {
                        jj = (int)mod(j+ystart[2]-2+ny/2, ny) - ny/2;
                        ij2 = jj*jj+ii2;
                        for(k = 1; k <= d1; k++) {
                            kk = (int)mod(k+zstart[2]-2+nz/2, nz) - nz/2;
                            idx = (((j-1)*d2+(i-1))*d1+(k-1)); //twiddle[j, i, k] 
                            twiddle[idx] = Math.Exp(ap * (double)(kk * kk + ij2));
                        }
                    }
                }
            }
            else if(layout_type == layout_2D) { // zxy layout
                for(i = 1; i <= d2; i++) {
                    ii =  (int)mod(i+xstart[2]-2+nx/2, nx) - nx/2;
                    ii2 = ii*ii;
                    for(j = 1; j <= d3; j++) {
                        jj = (int)mod(j+ystart[2]-2+ny/2, ny) - ny/2;
                        ij2 = jj*jj+ii2;
                        for(k = 1; k <= d1; k++) {
                            kk = (int)mod(k+zstart[2]-2+nz/2, nz) - nz/2;
                            idx = (((j-1)*d2+(i-1))*d1+(k-1)); // twiddle[j,i,k]
                            twiddle[idx] = Math.Exp(ap * (double)(kk*kk+ij2));
                        }
                    }
                }
            }
            else {
                Console.WriteLine(" Unknown layout type " + layout_type);
            }
        }

        public void compute_initial_conditions(double[, , ,] u1) {
            int k;
            double x0, start, an, dummy;
            double seed = 314159265, a = 1220703125;
            //c---------------------------------------------------------------------
            //c 0-D and 1-D layouts are easy because each processor gets a contiguous
            //c chunk of the array, in the Fortran ordering sense. 
            //c For a 2-D layout, it's a bit more complicated. We always
            //c have entire x-lines (contiguous) in processor. 
            //c We can do ny/np1 of them at a time since we have
            //c ny/np1 contiguous in y-direction. But then we jump
            //c by z-planes (nz/np2 of them, total). 
            //c For the 0-D and 1-D layouts we could do larger chunks, but
            //c this turns out to have no measurable impact on performance. 
            //c---------------------------------------------------------------------

            start = seed;

            //c---------------------------------------------------------------------
            //c Jump to the starting element for our first plane.
            //c---------------------------------------------------------------------

            an = ipow46(a, 2*nx, (zstart[0]-1)*ny + (ystart[0]-1));
            dummy = randlcGet(ref start, an);
            an = ipow46(a, 2*nx, ny);

            //c---------------------------------------------------------------------
            //c Go through by z planes filling in one square at a time.
            //c---------------------------------------------------------------------
            for(k = 0; k < dims[2, 0]; k++) { // nz/np2;
                x0 = start;
                //call vranlc(2*nx*dims(2, 1), x0, a, u0(1, 1, k)) : call native
                vranlc(2 * nx * dims[1, 0], x0, a, u1, k); //(1, 1, k));
                if((k+1) != dims[2, 0])
                    dummy = randlcGet(ref start, an);
            }
        }

        public double ipow46(double a, int exp_1, int exp_2) {
            //c---------------------------------------------------------------------
            //c compute a^exponent mod 2^46
            //c---------------------------------------------------------------------
            //double precision a, 
            double result;
            double dummy, q, r;
            int n, n2; //ierr;
            //external randlc;
            //double precision randlc;
            bool  two_pow;
            //c---------------------------------------------------------------------
            //c Use
            //c   a^n = a^(n/2)*a^(n/2) if n even else
            //c   a^n = a*a^(n-1)       if n odd
            //c---------------------------------------------------------------------
            result = 1;
            if(exp_2 == 0 || exp_1 == 0)
                return result;
            q = a;
            r = 1;
            n = exp_1;
            two_pow = true;

            while(two_pow) {
                n2 = n/2;
                if(n2 * 2 == n) {
                    dummy = randlcGet(ref q, q);
                    n = n2;
                }
                else {
                    n = n * exp_2;
                    two_pow = false;
                }
            }

            while(n > 1) {
                n2 = n/2;
                if(n2 * 2 == n) {
                    dummy = randlcGet(ref q, q);
                    n = n2;
                }
                else {
                    dummy = randlcGet(ref r, q);
                    n = n-1;
                }
            }
            dummy = randlcGet(ref r, q);
            result = r;
            return result;
        }

        public double randlcGet(ref double x, double a) {
            //c---------------------------------------------------------------------
            //c
            //c   This routine returns a uniform pseudorandom double precision number in the
            //c   range (0, 1) by using the linear congruential generator
            //c
            //c   x_{k+1} = a x_k  (mod 2^46)
            //c
            //c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
            //c   before repeating.  The argument A is the same as 'a' in the above formula,
            //c   and X is the same as x_0.  A and X must be odd double precision integers
            //c   in the range (1, 2^46).  The returned value RANDLC is normalized to be
            //c   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
            //c   the new seed x_1, so that subsequent calls to RANDLC using the same
            //c   arguments will generate a continuous sequence.
            //c
            //c   This routine should produce the same results on any computer with at least
            //c   48 mantissa bits in double precision floating point data.  On 64 bit
            //c   systems, double precision should be disabled.
            //c
            //c   David H. Bailey     October 26, 1990
            //c
            //c---------------------------------------------------------------------

            double r23, r46, t23, t46, t1, t2, t3, t4, a1, a2, x1, x2, z;

            r23 = Math.Pow(0.5, 23); //r23 = 0.5d0 ** 23
            r46 = Math.Pow(r23, 2);  //r46 = r23 ** 2
            t23 = Math.Pow(2.0, 23); //t23 = 2.d0 ** 23
            t46 = Math.Pow(t23, 2);  //t46 = t23 ** 2

            //c---------------------------------------------------------------------
            //c   Break A into two parts such that A = 2^23 * A1 + A2.
            //c---------------------------------------------------------------------
            t1 = r23 * a;
            a1 = (int)(t1);
            a2 = a - t23 * a1;

            //c---------------------------------------------------------------------
            //c   Break X into two parts such that X = 2^23 * X1 + X2, compute
            //c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
            //c   X = 2^23 * Z + A2 * X2  (mod 2^46).
            //c---------------------------------------------------------------------
            t1 = r23 * x;
            x1 = (int)(t1);
            x2 = x - t23 * x1;


            t1 = a1 * x2 + a2 * x1;
            t2 = (int)(r23 * t1);
            z = t1 - t23 * t2;
            t3 = t23 * z + a2 * x2;
            t4 = (int)(r46 * t3);
            x = t3 - t46 * t4;
            return (r46*x);
        }

        public void vranlc(int n, double x, double a, double[, , ,] u1, int k) {
            //---------------------------------------------------------------------
            //   This routine generates N uniform pseudorandom double precision numbers in
            //   the range (0, 1) by using the linear congruential generator
            //   
            //   x_{k+1} = a x_k  (mod 2^46)
            //   
            //   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
            //   before repeating.  The argument A is the same as 'a' in the above formula,
            //   and X is the same as x_0.  A and X must be odd double precision integers
            //   in the range (1, 2^46).  The N results are placed in Y and are normalized
            //   to be between 0 and 1.  X is updated to contain the new seed, so that
            //   subsequent calls to RANDLC using the same arguments will generate a
            //   continuous sequence.
            //   
            //   This routine generates the output sequence in batches of length NV, for
            //   convenience on vector computers.  This routine should produce the same
            //   results on any computer with at least 48 mantissa bits in double precision
            //   floating point data.  On Cray systems, double precision should be disabled.
            //   
            //   David H. Bailey    August 30, 1990
            //---------------------------------------------------------------------
            double r23, r46, t23, t46;
            int nv;
            r23 = Math.Pow(2.0, (-23));
            r46 = r23 * r23;
            t23 = Math.Pow(2.0, 23);
            t46 = t23 * t23;
            nv = 64;
            double[] xv = new double[nv];
            double t1, t2, t3, t4, an, a1=0, a2=0, x1, x2, yy;
            int n1, i, j;
            //---------------------------------------------------------------------
            //     Compute the first NV elements of the sequence using RANDLC.
            //---------------------------------------------------------------------
            t1 = x;
            n1 = (int)min(n, nv);
            for(i = 1; i <= n1; i++) {
                xv[i-1] = t46 * randlcGet(ref t1, a);
            }
            // ---------------------------------------------------------------------
            // It is not necessary to compute AN, A1 or A2 unless N is greater than NV.
            // ---------------------------------------------------------------------
            if(n > nv) {
                //---------------------------------------------------------------------
                //     Compute AN = AA ^ NV (mod 2^46) using successive calls to RANDLC.
                //---------------------------------------------------------------------
                t1 = a;
                t2 = r46 * a;

                for(i = 1; i <= nv - 1; i++) {
                    t2 = randlcGet(ref t1, a);
                }
                an = t46 * t2;
                //---------------------------------------------------------------------
                //     Break AN into two parts such that AN = 2^23 * A1 + A2.
                //---------------------------------------------------------------------
                t1 = r23 * an;
                a1 = (int)(t1);
                a2 = an - t23 * a1;
            }
            //---------------------------------------------------------------------
            //     Compute N pseudorandom results in batches of size NV.
            //---------------------------------------------------------------------
            for(j = 0; j <= n - 1; j = j + nv) {
                n1 = (int)min(nv, n - j);
                //---------------------------------------------------------------------
                //     Compute up to NV results based on the current seed vector XV.
                //---------------------------------------------------------------------
                int io;
                for(i = 0; i < n1; i++) { //y(i+j) = r46 * xv(i)
                    io = (i + j) + (2 * nx * dims[1, 0] * k);//idx=(i+j)+(2*nx*dims[1,0]*k); 
                    //Point.setValue(u1, idx, r46*xv[i]);  //u0[d3, d2, d1] //y[idx] = r46 * xv[i];
                    int m1 = (io % size1);
                    int m2 = (m1 % size2);
                    int _i = io/size1;
                    int _j = m1/size2;
                    int _k = m2/2;
                    int _t = (m2 % 2);

                    u1[_i, _j, _k, _t] = r46*xv[i];
                }
                //---------------------------------------------------------------------
                //     If this is the last pass through the 140 loop, it is not necessary to
                //     update the XV vector.
                //---------------------------------------------------------------------
                if(j + n1 == n) {
                    x = xv[n1-1];
                    return;
                }
                //---------------------------------------------------------------------
                //     Update the XV vector by multiplying each element by AN (mod 2^46).
                //---------------------------------------------------------------------
                for(i = 0; i < nv; i++) {
                    t1 = r23 * xv[i];
                    x1 = (int)(t1);
                    x2 = xv[i] - t23 * x1;
                    t1 = a1 * x2 + a2 * x1;
                    t2 = (int)(r23 * t1);
                    yy = t1 - t23 * t2;
                    t3 = t23 * yy + a2 * x2;
                    t4 = (int)(r46 * t3);
                    xv[i] = t3 - t46 * t4;
                }
            }
            //---------------------------------------------------------------------
            //     Save the last seed in X so that subsequent calls to VRANLC will generate
            //     a continuous sequence.
            //---------------------------------------------------------------------
        }

        public void fft_init(int n) {
            //c---------------------------------------------------------------------
            //c compute the roots-of-unity array that will be used for subsequent FFTs. 
            //c---------------------------------------------------------------------
            int m,nu,ku,i,j,ln;
            double t, ti;
            double pi = Math.PI;

            //c---------------------------------------------------------------------
            //c   Initialize the U array with sines and cosines in a manner that permits
            //c   stride one access at each FFT iteration.
            //c---------------------------------------------------------------------

            nu = n;
            m = ilog2(n);
            u[0, 0] = m; // u(1)
            u[0, 1] = 0.0;
            ku = 2;
            ln = 1;
            int idx;
            for(j = 1; j <= m; j++) {
                t = pi / ln;
                for(i = 0; i <= ln - 1; i++) {
                    ti = i * t;
                    //u[i+ku] = dcmplx (cos (ti), sin(ti));
                    idx = (i+ku)-1;
                    u[idx, REAL] = Math.Cos(ti);  //u[i+ku]
                    u[idx, IMAG] = Math.Sin(ti);  //u[i+ku]
                }
                ku = ku + ln;
                ln = 2 * ln;
            }
        }

        public int ilog2(int n) {
            int nn, lg;
            if(n == 1) {
                return 0;
            }
            lg = 1;
            nn = 2;
            while(nn < n) {
                nn = nn * 2;
                lg = lg + 1;
            }
            return lg;
        }

        public void fft(int dir, double[, , ,] u1, double[, , ,] u02) {
            //double complex u1(ntdivnp), u02(ntdivnp)
            //double complex scratch(fftblockpad_default*maxdim*2) !scratch equal y in cfftsX
            //c---------------------------------------------------------------------
            //c note: args u1, u02 must be different arrays
            //c note: args for cfftsx are (direction, layout, xin, xout, scratch)
            //c       xin/xout may be the same and it can be somewhat faster
            //c       if they are
            //c note: args for transpose are (layout1, layout2, xin, xout)
            //c       xin/xout must be different
            //c---------------------------------------------------------------------

            if(dir == 1) {
                if(layout_type == layout_0D) {
                    cffts1(dir, dims[0, 0], dims[1, 0], dims[2, 0], u1, u1);
                    cffts2(dir, dims[0, 1], dims[1, 1], dims[2, 1], u1, u1);
                    cffts3(dir, dims[0, 2], dims[1, 2], dims[2, 2], u1, u02);
                }
                else if(layout_type == layout_1D) {
                    cffts1(dir, dims[0, 0], dims[1, 0], dims[2, 0], u1, u1);
                    cffts2(dir, dims[0, 1], dims[1, 1], dims[2, 1], u1, u1);
                    //if(timers_enabled)
                    //    timer.start(T_transpose);
                    transpose_xy_z(1, 2, u1, u02);
                    //if(timers_enabled)
                    //    timer.stop(T_transpose);
                    cffts1(dir, dims[0, 2], dims[1, 2], dims[2, 2], u02, u02);
                }
                else if(layout_type == layout_2D) {
                    cffts1(dir, dims[0, 0], dims[1, 0], dims[2, 0], u1, u1);
                    //if(timers_enabled)
                    //    timer.start(T_transpose);
                    transpose_x_y(0, 1, u1, u02);
                    //if(timers_enabled)
                    //    timer.stop(T_transpose);
                    cffts1(dir, dims[0, 1], dims[1, 1], dims[2, 1], u02, u02);
                    //if(timers_enabled)
                    //    timer.start(T_transpose);
                    transpose_x_z(1, 2, u02, u1);
                    //if(timers_enabled)
                    //    timer.stop(T_transpose);
                    cffts1(dir, dims[0, 2], dims[1, 2], dims[2, 2], u1, u02);
                }
            }
            else {
                if(layout_type == layout_0D) {
                    cffts3(dir, dims[0, 2], dims[1, 2], dims[2, 2], u1, u1);
                    cffts2(dir, dims[0, 1], dims[1, 1], dims[2, 1], u1, u1);
                    cffts1(dir, dims[0, 0], dims[1, 0], dims[2, 0], u1, u02);
                }
                else if(layout_type == layout_1D) {
                    cffts1(dir, dims[0, 2], dims[1, 2], dims[2, 2], u1, u1);
                    //if(timers_enabled)
                    //    timer.start(T_transpose);
                    transpose_x_yz(2, 1, u1, u02);
                    //if(timers_enabled)
                    //    timer.stop(T_transpose);
                    cffts2(dir, dims[0, 1], dims[1, 1], dims[2, 1], u02, u02);
                    cffts1(dir, dims[0, 0], dims[1, 0], dims[2, 0], u02, u02);
                }
                else if(layout_type == layout_2D) {
                    cffts1(dir, dims[0, 2], dims[1, 2], dims[2, 2], u1, u1);
                    //if(timers_enabled)
                    //    timer.start(T_transpose);
                    transpose_x_z(2, 1, u1, u02);
                    //if(timers_enabled)
                    //    timer.stop(T_transpose);
                    cffts1(dir, dims[0, 1], dims[1, 1], dims[2, 1], u02, u02);
                    //if(timers_enabled)
                    //    timer.start(T_transpose);
                    transpose_x_y(1, 0, u02, u1);
                    //if(timers_enabled)
                    //    timer.stop(T_transpose);
                    cffts1(dir, dims[0, 0], dims[1, 0], dims[2, 0], u1, u02);
                }
            }
        }

        public void cffts1(int dir, int d1, int d2, int d3, double[, , ,] x, double[, , ,] xout) {
            int logd1;
            //Fortran
            //double complex x(d1,d2,d3);
            //double complex xout(d1,d2,d3);
            //double complex y(fftblockpad, d1, 2) ;
            //C#
            //y   [2, d1, fftblockpad]
            //x   [d3,d2,d1];
            //xout[d3,d2,d1];
            double[,,,] y = new double[2, d1, fftblockpad, 2];

            int i, j, k, jj, io;
            logd1 = ilog2(d1);
            for(k = 0; k < d3; k++) {
                for(jj = 0; jj <= (d2-fftblock); jj = jj + fftblock) {
                    //if(timers_enabled)
                    //    timer.start(T_fftcopy);
                    for(j = 0; j < fftblock; j++) {
                        for(i = 0; i < d1; i++) {//y(j,i,1) = x(i,j+jj,k)
                            io = ((k*d2+(j+jj))*d1+i)*2;
                            //y[0, i, j, REAL] = Point.getValue(x, io+REAL); //y[1,i,j,real] = x[k,j+jj,i,real]
                            //y[0, i, j, IMAG] = Point.getValue(x, io+IMAG);
                            int m1 = (io % size1);
                            int m2 = (m1 % size2);
                            int _i = io/size1;
                            int _j = m1/size2;
                            int _k = m2/2;
                            //int _t = (int)(m2 % dm4);
                            y[0, i, j, REAL] = x[_i, _j, _k, REAL];
                            y[0, i, j, IMAG] = x[_i, _j, _k, IMAG];
                        }
                    }
                    //if(timers_enabled)
                    //    timer.stop(T_fftcopy);
                    //if(timers_enabled)
                    //    timer.start(T_fftlow);
                    cfftz(dir, logd1, d1, y); //cfftz (iis, logd1, d1, y, y(1,1,2)); 
                    //if(timers_enabled)
                    //    timer.stop(T_fftlow);

                    //if(timers_enabled)
                    //    timer.start(T_fftcopy);

                    for(j = 0; j < fftblock; j++) {
                        for(i = 0; i < d1; i++) {
                            //iin   = ((0*d1+i)*fftblockpad+j)*2;
                            io  = (((k*d2+(j+jj))*d1+i)*2);
                            //Point.setAddress(y, iin+REAL, xout, io+REAL);
                            //Point.setAddress(y, iin+IMAG, xout, io+IMAG);
                            int m1 = (io % size1);
                            int m2 = (m1 % size2);
                            int _i = io/size1;
                            int _j = m1/size2;
                            int _k = m2/2;
                            //int _t = (int)(m2 % dm4);
                            xout[_i, _j, _k, REAL] = y[0, i, j, REAL];//xout(i,j+jj,k) = y(j,i,1)
                            xout[_i, _j, _k, IMAG] = y[0, i, j, IMAG];
                        }
                    }
                    //if(timers_enabled)
                    //    timer.stop(T_fftcopy);
                }
            }
        }

        public void cffts2(int dir, int d1, int d2, int d3, double[, , ,] x, double[, , ,] xout) {
            int logd2;
            //Fortran: double complex x(d1,d2,d3);
            //         double complex xout(d1,d2,d3);
            //         double complex y(fftblockpad, d2, 2);
            //C#:  x   [d3,d2,d1];
            //     xout[d3,d2,d1];
            //     y   [2, d2, fftblockpad];
            double[,,,] y = new double[2, d2, fftblockpad, 2];

            int i, j, k, ii, io;
            logd2 = ilog2(d2);
            for(k = 0; k < d3; k++) {
                for(ii = 0; ii <= d1 - fftblock; ii = ii + fftblock) {
                    //if(timers_enabled)
                    //    timer.start(T_fftcopy);
                    for(j = 0; j < d2; j++) {
                        for(i = 0; i < fftblock; i++) {
                            io = ((k*d2+j)*d1+(i+ii))*2;
                            //y[0, j, i, REAL] = Point.getValue(x, io+REAL);
                            //y[0, j, i, IMAG] = Point.getValue(x, io+IMAG);
                            int m1 = (io % size1);
                            int m2 = (m1 % size2);
                            int _i = io/size1;
                            int _j = m1/size2;
                            int _k = m2/2;
                            //int _t = (int)(m2 % dm4);
                            y[0, j, i, REAL] = x[_i, _j, _k, REAL];
                            y[0, j, i, IMAG] = x[_i, _j, _k, IMAG];
                        }
                    }
                    //if(timers_enabled)
                    //    timer.stop(T_fftcopy);

                    //if(timers_enabled)
                    //    timer.start(T_fftlow);
                    cfftz(dir, logd2, d2, y); //y(1, 1, 2));
                    //if(timers_enabled)
                    //    timer.stop(T_fftlow);

                    //if(timers_enabled)
                    //    timer.start(T_fftcopy);
                    for(j = 0; j < d2; j++) {
                        for(i = 0; i < fftblock; i++) {
                            //iin = ((0*d2+j)*fftblockpad+i)*2;
                            io = ((k * d2 + j) * d1 + (i + ii)) * 2;
                            //Point.setAddress(y, iin+REAL, xout, io+REAL);
                            //Point.setAddress(y, iin+IMAG, xout, io+IMAG);
                            int m1 = (io % size1);
                            int m2 = (m1 % size2);
                            int _i = io/size1;
                            int _j = m1/size2;
                            int _k = m2/2;
                            //int _t = (int)(m2 % dm4);
                            xout[_i, _j, _k, REAL] = y[0, j, i, REAL];
                            xout[_i, _j, _k, IMAG] = y[0, j, i, IMAG];
                        }
                    }
                    //if(timers_enabled)
                    //    timer.stop(T_fftcopy);
                }
            }
        }

        public void cffts3(int dir, int d1, int d2, int d3, double[, , ,] x, double[, , ,] xout) {
            int logd3;
            //Fortran: double complex x(d1,d2,d3);
            //         double complex xout(d1,d2,d3);
            //         double complex y(fftblockpad, d3, 2); 
            //C#:  x    [d3,d2,d1];
            //     xout [d3,d2,d1];
            //     y    [2, d3, fftblockpad]; 
            double[,,,] y = new double[2, d3, fftblockpad, 2];

            int i, j, k, ii, io;

            logd3 = ilog2(d3);

            for(j = 0; j < d2; j++) {
                for(ii = 0; ii <= d1 - fftblock; ii = ii + fftblock) {
                    //if(timers_enabled)
                    //    timer.start(T_fftcopy);
                    for(k = 0; k < d3; k++) {
                        for(i = 0; i < fftblock; i++) {
                            io = ((k*d2+j)*d1+(i+ii))*2;
                            //y[0, k, i, REAL] = Point.getValue(x, io+REAL);
                            //y[0, k, i, IMAG] = Point.getValue(x, io+IMAG);
                            int m1 = (io % size1);
                            int m2 = (m1 % size2);
                            int _i = io/size1;
                            int _j = m1/size2;
                            int _k = m2/2;
                            //int _t = (int)(m2 % dm4);
                            y[0, k, i, REAL] = x[_i, _j, _k, REAL];
                            y[0, k, i, IMAG] = x[_i, _j, _k, IMAG];
                        }
                    }
                    //if(timers_enabled)
                    //    timer.stop(T_fftcopy);

                    //if(timers_enabled)
                    //    timer.start(T_fftlow);
                    cfftz(dir, logd3, d3, y); //y(1, 1, 2));
                    //if(timers_enabled)
                    //    timer.stop(T_fftlow);

                    //if(timers_enabled)
                    //    timer.start(T_fftcopy);
                    for(k = 0; k < d3; k++) {
                        for(i = 0; i < fftblock; i++) {
                            //iin = ((0*d3+k)*fftblockpad+i)*2;
                            io  = (((k*d2+j)*d1+(i+ii))*2);
                            //Point.setAddress(y, iin+REAL, xout, io+REAL);
                            //Point.setAddress(y, iin+IMAG, xout, io+IMAG);
                            int m1 = (io % size1);
                            int m2 = (m1 % size2);
                            int _i = io/size1;
                            int _j = m1/size2;
                            int _k = m2/2;
                            //int _t = (int)(m2 % dm4);
                            xout[_i, _j, _k, REAL] = y[0, k, i, REAL];
                            xout[_i, _j, _k, IMAG] = y[0, k, i, IMAG];
                        }
                    }
                    //if(timers_enabled)
                    //    timer.stop(T_fftcopy);
                }
            }
        }

        public void cfftz(int dir, int m, int n, double[, , ,] y) {
            //c---------------------------------------------------------------------
            //c   Computes NY N-point complex-to-complex FFTs of X using an algorithm due
            //c   to Swarztrauber.  X is both the input and the output array, while Y is a 
            //c   scratch array.  It is assumed that N = 2^M.  Before calling CFFTZ to 
            //c   perform FFTs, the array U must be initialized by calling CFFTZ with IS 
            //c   set to 0 and M set to MX, where MX is the maximum value of M for any 
            //c   subsequent call.
            //c---------------------------------------------------------------------
            //Fortran: dimension x(fftblockpad,n), y(fftblockpad,n);
            int i,j,l,mx;
            //C#: dimension x[n,fftblockpad], y[n,fftblockpad];
            //c---------------------------------------------------------------------
            //c   Check if input parameters are invalid.
            //c---------------------------------------------------------------------
            mx = (int)u[0, 0]; //mx = u(1);
            if((dir != 1 && dir != -1) || m < 1 || m > mx) {
                Console.WriteLine("CFFTZ: Either U has not been initialized, or else one of the input parameters iis invalid " + dir + " " + m + " " + mx);
            }
            //c---------------------------------------------------------------------
            //c   Perform one variant of the Stockham FFT.
            //c---------------------------------------------------------------------
            for(l = 1; l <= m; l = l + 2) {
                fftz2(dir, l, m, n, fftblock, fftblockpad, u, y, 0, 1);
                if(l == m) {
                    for(j = 0; j < n; j++) {
                        for(i = 0; i < fftblock; i++) { //x(i,j) = y(i,j);   //C#: dimension x[n,fftblockpad], y[n,fftblockpad];
                            y[0, j, i, REAL] = y[1, j, i, REAL];
                            y[0, j, i, IMAG] = y[1, j, i, IMAG];
                        }
                    }
                    return;
                }
                fftz2(dir, l + 1, m, n, fftblock, fftblockpad, u, y, 1, 0);
            }
        }

        public void fftz2(int dir, int l, int m, int n, int ny, int ny1, double[,] u, double[, , ,] y, int iread, int iwrite) { //u=u x=ytemp y = ytemp
            //c---------------------------------------------------------------------
            //c   Performs the L-th iteration of the second variant of the Stockham FFT.
            //c---------------------------------------------------------------------
            int k, n1, li, lj, lk, ku, i, j, i11, i12, i21, i22;
            //double complex u,x,y,u1,x11,x21;
            double[] u1 = new double[2];
            double[] x11= new double[2];
            double[] x21= new double[2];
            //dimension u(n), x(ny1,n), y(ny1,n);
            //c---------------------------------------------------------------------
            //c   Set initial parameters.
            //c---------------------------------------------------------------------

            n1 = n / 2;
            lk = (int)Math.Pow(2, (l - 1));
            li = (int)Math.Pow(2, (m - l));
            lj = 2 * lk;
            ku = li;// +1;

            for(i = 0; i <= li - 1; i++) {
                i11 = i * lk;
                i12 = i11 + n1;
                i21 = i * lj;
                i22 = i21 + lk;

                u1[REAL] = u[(ku+i), REAL];
                if(dir >= 1) {
                    //    u1 = u(ku+i);
                    u1[1] = u[ku+i, IMAG];
                }
                else {
                    //    u1 = dconjg (u(ku+i));
                    u1[1] = -1*u[ku+i, IMAG];
                }

                //  c---------------------------------------------------------------------
                //  c   This loop is vectorizable.
                //  c---------------------------------------------------------------------
                for(k = 0; k <= lk - 1; k++) {
                    for(j = 0; j < ny; j++) {
                        x11[REAL] = y[iread, i11 + k, j, REAL]; //x11[0] = x[j,i11+k];
                        x11[IMAG] = y[iread, i11 + k, j, IMAG];
                        x21[REAL] = y[iread, i12 + k, j, REAL]; //x21 = x(j,i12+k);
                        x21[IMAG] = y[iread, i12 + k, j, IMAG];
                        y[iwrite, i21 + k, j, REAL] = x11[REAL] + x21[REAL]; //y(j,i21+k) = x11 + x21;
                        y[iwrite, i21 + k, j, IMAG] = x11[IMAG] + x21[IMAG];
                        y[iwrite, i22 + k, j, REAL] = u1[REAL] * (x11[REAL] - x21[REAL]) - u1[IMAG] * (x11[IMAG] - x21[IMAG]); //y(j,i22+k) = u1 * (x11 - x21);
                        y[iwrite, i22 + k, j, IMAG] = u1[IMAG] * (x11[REAL] - x21[REAL]) + u1[REAL] * (x11[IMAG] - x21[IMAG]);
                    }
                }
            }
        }

        public void transpose_x_y(int l1, int l2, double[, , ,] xin, double[, , ,] xout) {  // x1, x2); x1=u1 x2 = u0
            /*  double complex xin(ntdivnp), xout(ntdivnp)
                ---------------------------------------------------------------------
                 xy transpose is a little tricky, since we don't want
                 to touch 3rd axis. But alltoall must involve 3rd axis (most 
                 slowly varying) to be efficient. So we do
                 (nx, ny/np1, nz/np2) -> (ny/np1, nz/np2, nx) (local)
                 (ny/np1, nz/np2, nx) -> ((ny/np1*nz/np2)*np1, nx/np1) (global)
                 then local finish. 
                --------------------------------------------------------------------- */
            //imprimir4(xin, 0, "xin");

            transpose_x_y_local(dims[0, l1], dims[1, l1], dims[2, l1], xin, xout);

            transpose_x_y_global(dims[0, l1], dims[1, l1], dims[2, l1], xout, xin);

            transpose_x_y_finish(dims[0, l2], dims[1, l2], dims[2, l2], xin, xout);

        }

        public void transpose_x_y_local(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            //implicit none
            //include 'global.h'
            //integer d1, d2, d3
            //double complex uxin(d1, d2, d3)  ===> [d3, d2, d1]       
            //double complex uxout(d2, d3, d1) ===> [d1, d3, d2]
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;

            //if(timers_enabled)
            //    timer.start(T_transxyloc);
            //if(timers_enabled)
            //    timer.start(T_transxzloc);
            int i, j, k, ii, io;
            for(k = 0; k < d3; k++) { //k=1 <= d3        
                for(i = 0; i < d1; i++) { //i=1 <= d1 
                    for(j = 0; j < d2; j++) { // j=1 <= d2                     //xout(j,k,i)    =    xin(i,j,k); 
                        ii = ((k*d2+j)*d1+i)*2;
                        io = ((i*d3+k)*d2+j)*2;
                        //Point.setAddress(xin, ii+REAL, xout, io+REAL);
                        //Point.setAddress(xin, ii+IMAG, xout, io+IMAG);
                        //xout[i, k, j, REAL] = xin[k, j, i, REAL];
                        //xout[i, k, j, IMAG] = xin[k, j, i, IMAG];

                        m1 = (ii % size1);
                        m2 = (m1 % size2);
                        _i = ii/size1;
                        _j = m1/size2;
                        _k = m2/2;

                        om1 = (io % size1);
                        om2 = (om1 % size2);
                        o_i = io/size1;
                        o_j = om1/size2;
                        o_k = om2/2;

                        xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                        xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                    }
                }
            }
            //if(timers_enabled)
            //    timer.stop(T_transxyloc);
        }

/**/    public void transpose_x_y_global(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            /* ---------------------------------------------------------------------
               array is in form (ny/np1, nz/np2, nx)
               ---------------------------------------------------------------------
               double complex uxin(d2,d3,d1)
               double complex uxout(d2,d3,d1) ! not real layout but right size      */
            //if(timers_enabled)
            //    synchup();
            /* ---------------------------------------------------------------------
               do transpose among all processes with same 1-coord (me1)
               --------------------------------------------------------------------- */
            //if(timers_enabled)
            //    timer.start(T_transxyglo);
            double[] src       = new double[d1*d2*d3*2];
            double[] dst       = new double[d1*d2*d3*2];
            setVetor(xin, src);
            commslice2.AlltoallFlattened<double>(src, d1*d2*d3*2/np1, ref dst);
            setVetor(dst, xout);
            //if(timers_enabled)
            //    timer.stop(T_transxyglo);
        }

/**/    public void transpose_x_y_finish(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            //Fortran
            //double complex xin(d1/np1, d3, d2, 0:np1-1); 
            //double complex xout(d1,d2,d3);
            //C#
            //uxin [np1, d2, d3, d1/np1, 2]; 
            //uxout[d3    , d2, d1,   2, 0];

            int i, j, k, p, ioff;
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;

            //if(timers_enabled)
            //    timer.start(T_transxyfin);
            int io=0, ii=0;
            for(p = 0; p <= (np1-1); p++) {
                ioff = p*d1/np1;
                for(k = 0; k < d3; k++) {
                    for(j = 0; j < d2; j++) {
                        for(i = 0; i < (d1/np1); i++) { //io = ((((k-1)*size3+(j-1))*size4*np1+(i+ioff-1))) *2;  //uxout[k,j,i+ioff] ii  = ((((p)*size2+(j-1))*size4+(k-1))*size3+(i-1))*2; //uxin[p,j,k,i]
                            ii  = (((p*d2+j)*d3+k)*(d1/np1)+i)*2;
                            io  = ((k*d2+j)*d1+i+ioff)*2;
                            //Point.setAddress(xin, ii+REAL, xout, io+REAL);    //xout(i+ioff,j,k) = xin(i,k,j,p);
                            //Point.setAddress(xin, ii+IMAG, xout, io+IMAG);
                            m1 = (ii % size1);
                            m2 = (m1 % size2);
                            _i = ii/size1;
                            _j = m1/size2;
                            _k = m2/2;

                            om1 = (io % size1);
                            om2 = (om1 % size2);
                            o_i = io/size1;
                            o_j = om1/size2;
                            o_k = om2/2;

                            xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                            xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                        }
                    }
                }
            }
            //if(timers_enabled)
            //    timer.stop(T_transxyfin);
        }

        public void transpose_xy_z(int l1, int l2, double[, , ,] xin, double[, , ,] xout) {
            //double complex xin(ntdivnp), xout(ntdivnp)
            transpose_xy_z_local(dims[0, l1], dims[1, l1], dims[2, l1], xin, xout);
            transpose_xy_z_global(dims[1, 0], dims[2, 0], dims[0, 0], xout, xin);
            transpose_xy_z_finish(dims[0, l1],dims[1, l1], dims[2, l1], xin, xout);
        }

        public void transpose_xy_z_local(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            //Fortran
            //double complex xin(n1, n2), xout(n2, n1)
            //double complex z(transblockpad, transblock)
            //C#
            //xin [n2,n1]
            //xout[n1,n2]
            //z[transblock, transblockpad]

            double[,,] z = new double[transblockpad, transblock, 2];
            int i, j, ii, jj, iin, io;
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;
            int n1 = d1*d2;
            int n2 = d3;

            //if(timers_enabled)
            //    timer.start(T_transxzloc);

            //  c---------------------------------------------------------------------
            //  c If possible, block the transpose for cache memory systems. 
            //  c How much does this help? Example: R8000 Power Challenge (90 MHz)
            //  c Blocked version decreases time spend in this routine 
            //  c from 14 seconds to 5.2 seconds on 8 nodes class A.
            //  c---------------------------------------------------------------------

            if(n1 < transblock || n2 < transblock) {
                if(n1 >= n2) {
                    for(j = 0; j < n2; j++) {
                        for(i = 0; i < n1; i++) {
                            //uxout[i, j] = uxin[j, i];                 //xout(j, i) = xin(i, j);
                            iin = (j * n1 + i)*2;
                            io = (i * n2 + j)*2;
                            //Point.setAddress(xin, iin+REAL, xout, io+REAL);
                            //Point.setAddress(xin, iin+IMAG, xout, io+IMAG);
                            m1 = (iin % size1);
                            m2 = (m1 % size2);
                            _i = iin/size1;
                            _j = m1/size2;
                            _k = m2/2;

                            om1 = (io % size1);
                            om2 = (om1 % size2);
                            o_i = io/size1;
                            o_j = om1/size2;
                            o_k = om2/2;

                            xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                            xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                        }
                    }
                }
                else {
                    for(i = 0; i < n1; i++) {
                        for(j = 0; j < n2; j++) {                   //xout(j, i) = xin(i, j);
                            iin = (j * n1 + i)*2;
                            io = (i * n2 + j)*2;
                            //Point.setAddress(xin, iin+REAL, xout, io+REAL);
                            //Point.setAddress(xin, iin+IMAG, xout, io+IMAG);
                            m1 = (iin % size1);
                            m2 = (m1 % size2);
                            _i = iin/size1;
                            _j = m1/size2;
                            _k = m2/2;

                            om1 = (io % size1);
                            om2 = (om1 % size2);
                            o_i = io/size1;
                            o_j = om1/size2;
                            o_k = om2/2;

                            xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                            xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                        }
                    }
                }
            }
            else {
                for(j = 0; j <= n2 - 1; j = j + transblock) {
                    for(i = 0; i <= n1 - 1; i = i + transblock) {
                        //c---------------------------------------------------------------------
                        //c Note: compiler should be able to take j+jj out of inner loop
                        //c---------------------------------------------------------------------
                        for(jj = 0; jj < transblock; jj++) {
                            for(ii = 0; ii < transblock; ii++) { //z(jj,ii) = xin(i+ii, j+jj);
                                iin = ((j+jj)*n1+(i+ii))*2;      //xin[j+jj, i+ii];
                                //io = ((ii)*transblockpad+jj)*2; //z[ii,jj]
                                //Point.setAddress(xin, iin + REAL, z, io + REAL);
                                //Point.setAddress(xin, iin + IMAG, z, io + IMAG);
                                m1 = (iin % size1);
                                m2 = (m1 % size2);
                                _i = iin/size1;
                                _j = m1/size2;
                                _k = m2/2;

                                z[ii,jj,REAL] = xin[_i, _j, _k, REAL];
                                z[ii,jj,IMAG] = xin[_i, _j, _k, IMAG];
                            }
                        }
                        for(ii = 0; ii < transblock; ii++) {
                            for(jj = 0; jj < transblock; jj++) {//xout(j+jj, i+ii) = z(jj,ii);
                                //iin = (ii*transblockpad+jj)*2;  //z[ii,jj];
                                io = ((i+ii)*n2+(j+jj))*2;//xout[i+ii, j+jj]
                                //Point.setAddress(z, iin + REAL, xout, io + REAL);
                                //Point.setAddress(z, iin + IMAG, xout, io + IMAG);
                                m1 = (io % size1);
                                m2 = (m1 % size2);
                                _i = io/size1;
                                _j = m1/size2;
                                _k = m2/2;
                                xout[_i, _j, _k, REAL] = z[ii, jj, REAL];
                                xout[_i, _j, _k, IMAG] = z[ii, jj, IMAG];
                            }
                        }

                    }
                }
            }
            //if(timers_enabled)
            //    timer.stop(T_transxzloc);

        }

        public void transpose_x_yz_local(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            //Fortran
            //double complex xin(n1, n2), xout(n2, n1)
            //double complex z(transblockpad, transblock)
            //C#
            //xin [n2,n1]
            //xout[n1,n2]
            //z[transblock, transblockpad]

            double[,,] z = new double[transblockpad, transblock, 2];
            int i, j, ii, jj, iin, io;
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;
            int n1 = d1;
            int n2 = d2*d3;

            //if(timers_enabled)
            //    timer.start(T_transxzloc);

            //  c---------------------------------------------------------------------
            //  c If possible, block the transpose for cache memory systems. 
            //  c How much does this help? Example: R8000 Power Challenge (90 MHz)
            //  c Blocked version decreases time spend in this routine 
            //  c from 14 seconds to 5.2 seconds on 8 nodes class A.
            //  c---------------------------------------------------------------------

            if(n1 < transblock || n2 < transblock) {
                if(n1 >= n2) {
                    for(j = 0; j < n2; j++) {
                        for(i = 0; i < n1; i++) {
                            //uxout[i, j] = uxin[j, i];                 //xout(j, i) = xin(i, j);
                            iin = (j * n1 + i)*2;
                            io = (i * n2 + j)*2;
                            //Point.setAddress(xin, iin+REAL, xout, io+REAL);
                            //Point.setAddress(xin, iin+IMAG, xout, io+IMAG);
                            m1 = (iin % size1);
                            m2 = (m1 % size2);
                            _i = iin/size1;
                            _j = m1/size2;
                            _k = m2/2;

                            om1 = (io % size1);
                            om2 = (om1 % size2);
                            o_i = io/size1;
                            o_j = om1/size2;
                            o_k = om2/2;

                            xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                            xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                        }
                    }
                }
                else {
                    for(i = 0; i < n1; i++) {
                        for(j = 0; j < n2; j++) {                   //xout(j, i) = xin(i, j);
                            iin = (j * n1 + i)*2;
                            io = (i * n2 + j)*2;
                            //Point.setAddress(xin, iin+REAL, xout, io+REAL);
                            //Point.setAddress(xin, iin+IMAG, xout, io+IMAG);
                            m1 = (iin % size1);
                            m2 = (m1 % size2);
                            _i = iin/size1;
                            _j = m1/size2;
                            _k = m2/2;

                            om1 = (io % size1);
                            om2 = (om1 % size2);
                            o_i = io/size1;
                            o_j = om1/size2;
                            o_k = om2/2;

                            xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                            xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                        }
                    }
                }
            }
            else {
                for(j = 0; j <= n2 - 1; j = j + transblock) {
                    for(i = 0; i <= n1 - 1; i = i + transblock) {
                        //c---------------------------------------------------------------------
                        //c Note: compiler should be able to take j+jj out of inner loop
                        //c---------------------------------------------------------------------
                        for(jj = 0; jj < transblock; jj++) {
                            for(ii = 0; ii < transblock; ii++) { //z(jj,ii) = xin(i+ii, j+jj);
                                iin = ((j+jj)*n1+(i+ii))*2;      //xin[j+jj, i+ii];
                                //io = ((ii)*transblockpad+jj)*2; //z[ii,jj]
                                //Point.setAddress(xin, iin + REAL, z, io + REAL);
                                //Point.setAddress(xin, iin + IMAG, z, io + IMAG);
                                m1 = (iin % size1);
                                m2 = (m1 % size2);
                                _i = iin/size1;
                                _j = m1/size2;
                                _k = m2/2;

                                z[ii, jj, REAL] = xin[_i, _j, _k, REAL];
                                z[ii, jj, IMAG] = xin[_i, _j, _k, IMAG];
                            }
                        }
                        for(ii = 0; ii < transblock; ii++) {
                            for(jj = 0; jj < transblock; jj++) {//xout(j+jj, i+ii) = z(jj,ii);
                                //iin = (ii*transblockpad+jj)*2;  //z[ii,jj];
                                io = ((i+ii)*n2+(j+jj))*2;//xout[i+ii, j+jj]
                                //Point.setAddress(z, iin + REAL, xout, io + REAL);
                                //Point.setAddress(z, iin + IMAG, xout, io + IMAG);
                                m1 = (io % size1);
                                m2 = (m1 % size2);
                                _i = io/size1;
                                _j = m1/size2;
                                _k = m2/2;
                                xout[_i, _j, _k, REAL] = z[ii, jj, REAL];
                                xout[_i, _j, _k, IMAG] = z[ii, jj, IMAG];
                            }
                        }

                    }
                }
            }
            //if(timers_enabled)
            //    timer.stop(T_transxzloc);

        }

/******/public void transpose_xy_z_global(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            // double complex xin(ntdivnp)
            // double complex xout(ntdivnp) 
            double[] src = new double[d1*d2*d3*2];//ntdivnp*2
            double[] dst = new double[d1*d2*d3*2];//ntdivnp*2
            //if(timers_enabled)
            //    synchup();
            //if(timers_enabled)
            //    timer.start(T_transxzglo);
            setVetor(xin, src);
            commslice1.AlltoallFlattened<double>(src, (d1*d2*d3*2)/np, ref dst);//ntdivnp*2
            setVetor(dst, xout);
            // call mpi_alltoall(xin, ntdivnp/np, dc_type, xout, ntdivnp/np, dc_type, commslice1, ierr);
            //if(timers_enabled)
            //    timer.stop(T_transxzglo);

        }

        public void transpose_x_yz_global(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            // double complex xin(ntdivnp)
            // double complex xout(ntdivnp) 
            double[] src = new double[d1*d2*d3*2];//ntdivnp*2
            double[] dst = new double[d1*d2*d3*2];//ntdivnp*2
            //if(timers_enabled)
            //    synchup();
            //if(timers_enabled)
            //    timer.start(T_transxzglo);
            setVetor(xin, src);
            commslice1.AlltoallFlattened<double>(src, d1*d2*d3*2 / np, ref dst);//ntdivnp*2
            setVetor(dst, xout);
            // call mpi_alltoall(xin, ntdivnp/np, dc_type, xout, ntdivnp/np, dc_type, commslice1, ierr);
            //if(timers_enabled)
            //    timer.stop(T_transxzglo);

        }

/**/    public void transpose_xy_z_finish(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {//Ok
            //Fortran
            //double complex xin(n2, n1/np2, 0:np2-1), 
            //               xout(n2*np2, n1/np2)
            //C#
            //uxin  [   np2, n1/np2, n2] 
            //uxout [n1/np2, n2*np2    ]
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;
            int n1 = d1*d2;
            int n2 = d3;
            int i, j, p, ioff, ii, io;
            //if(timers_enabled)
            //    timer.start(T_transxzfin);
            for(p = 0; p <= np2 - 1; p++) {
                ioff = p*n2;
                for(j = 0; j < n1 / np2; j++) {
                    for(i = 0; i < n2; i++) { //xout(i+ioff, j) = xin(i, j, p);
                        ii = ((p*(n1/np2)+j)*n2+i)*2; //uxin[p, j, i]
                        io = (j*n2*np2+(i+ioff))*2; //uxout[j, i+ioff]
                        //Point.setAddress(xin, ii+REAL, xout, io+REAL);
                        //Point.setAddress(xin, ii+IMAG, xout, io+IMAG);
                        m1 = (ii % size1);
                        m2 = (m1 % size2);
                        _i = ii/size1;
                        _j = m1/size2;
                        _k = m2/2;

                        om1 = (io % size1);
                        om2 = (om1 % size2);
                        o_i = io/size1;
                        o_j = om1/size2;
                        o_k = om2/2;

                        xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                        xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                    }
                }
            }
            //if(timers_enabled)
            //    timer.stop(T_transxzfin);

        }

/**/    public void transpose_x_yz_finish(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) { //ok
            //Fortran
            //double complex xin(n2, n1/np2, 0:np2-1), 
            //               xout(n2*np2, n1/np2)
            //C#
            //uxin  [   np2, n1/np2, n2] 
            //uxout [n1/np2, n2*np2    ]
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;
            int n1 = d1;
            int n2 = d2*d3;
            int i, j, p, ioff, ii, io;
            //if(timers_enabled)
            //    timer.start(T_transxzfin);
            for(p = 0; p <= np2 - 1; p++) {
                ioff = p*n2;
                for(j = 0; j < n1 / np2; j++) {
                    for(i = 0; i < n2; i++) { //xout(i+ioff, j) = xin(i, j, p);
                        ii = ((p*(n1/np2)+j)*n2+i)*2; //uxin[p, j, i]
                        io = (j*n2*np2+(i+ioff))*2; //uxout[j, i+ioff]
                        //Point.setAddress(xin, ii+REAL, xout, io+REAL);
                        //Point.setAddress(xin, ii+IMAG, xout, io+IMAG);
                        m1 = (ii % size1);
                        m2 = (m1 % size2);
                        _i = ii/size1;
                        _j = m1/size2;
                        _k = m2/2;

                        om1 = (io % size1);
                        om2 = (om1 % size2);
                        o_i = io/size1;
                        o_j = om1/size2;
                        o_k = om2/2;

                        xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                        xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                    }
                }
            }
            //if(timers_enabled)
            //    timer.stop(T_transxzfin);

        }

        public void transpose_x_z(int l1, int l2, double[, , ,] xin, double[, , ,] xout) {
            //double complex xin(ntdivnp), xout(ntdivnp);

            transpose_x_z_local(dims[0, l1], dims[1, l1], dims[2, l1], xin, xout);
            transpose_x_z_global(dims[0, l1], dims[1, l1], dims[2, l1], xout, xin);
            transpose_x_z_finish(dims[0, l2], dims[1, l2], dims[2, l2], xin, xout);
        }

        public void transpose_x_z_local(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            //Fortran
            //double complex xin(d1,d2,d3)
            //double complex xout(d3,d2,d1)
            //double complex buf(transblockpad, maxdim)
            //C#
            //xin [d3,d2,d1]
            //xout[d1,d2,d3]
            double[,,] buf = new double[maxdim, transblockpad, 2];
            int block1, block3;
            int i, j, k, kk, ii, i1, k1;
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;

            //if(timers_enabled)
            //    timer.start(T_transxzloc);
            //if(timers_enabled)
            //    timer.start(T_transxzloc);
            if(d1 < 32)
                goto G100;
            block3 = d3;
            if(block3 == 1)
                goto G100;
            if(block3 > transblock)
                block3 = transblock;
            block1 = d1;
            if(block1*block3 > transblock*transblock)
                block1 = transblock*transblock/block3;
            /*---------------------------------------------------------------------
              blocked transpose
              ---------------------------------------------------------------------*/
            int iin = 0, io = 0;
            for(j = 0; j < d2; j++) {
                for(kk = 0; kk <= d3 - block3; kk = kk + block3) {
                    for(ii = 0; ii <= d1 - block1; ii = ii + block1) {
                        for(k = 0; k < block3; k++) {
                            k1 = k + kk;
                            for(i = 0; i < block1; i++) {                                 //buf(k, i) = xin(i+ii, j, k1);
                                iin = ((k1*d2+j)*d1+(i+ii))*2; //xin[k1, j, i+ii];
                                //io  = (i*transblockpad+k)*2;   //buf[i, k]
                                //Point.setAddress(xin, iin + REAL, buf, io + REAL);
                                //Point.setAddress(xin, iin + IMAG, buf, io + IMAG);
                                m1 = (iin % size1);
                                m2 = (m1 % size2);
                                _i = iin/size1;
                                _j = m1/size2;
                                _k = m2/2;

                                buf[i, k, REAL] = xin[_i, _j, _k, REAL];
                                buf[i, k, IMAG] = xin[_i, _j, _k, IMAG];
                            }
                        }
                        for(i = 0; i < block1; i++) {
                            i1 = i + ii;
                            for(k = 0; k < block3; k++) {                                 //xout(k+kk, j, i1) = buf(k, i);
                                //iin = (i*transblockpad+k)*2; //buf[i, k];
                                io  = ((i1*d2+j)*d3+(k+kk))*2;//xout[i1, j, k+kk]
                                //Point.setAddress(buf, iin + REAL, xout, io + REAL);
                                //Point.setAddress(buf, iin + IMAG, xout, io + IMAG);
                                m1 = (io % size1);
                                m2 = (m1 % size2);
                                _i = io/size1;
                                _j = m1/size2;
                                _k = m2/2;

                                xout[_i, _j, _k, REAL] = buf[i, k, REAL];
                                xout[_i, _j, _k, IMAG] = buf[i, k, IMAG];
                            }
                        }
                    }
                }
            }
            goto G200;
        //---------------------------------------------------------------------
        // basic transpose
        //---------------------------------------------------------------------
        G100:  //continue;
            for(j = 0; j < d2; j++) {
                for(k = 0; k < d3; k++) {
                    for(i = 0; i < d1; i++) {                                             //xout(k, j, i) = xin(i, j, k);
                        iin = ((k*d2+j)*d1+i)*2; //xin[k, j, i];
                        io  = ((i*d2+j)*d3+k)*2; //xout[i, j, k]
                        //Point.setAddress(xin, iin + REAL, xout, io + REAL);
                        //Point.setAddress(xin, iin + IMAG, xout, io + IMAG);
                        m1 = (iin % size1);
                        m2 = (m1 % size2);
                        _i = iin/size1;
                        _j = m1/size2;
                        _k = m2/2;

                        om1 = (io % size1);
                        om2 = (om1 % size2);
                        o_i = io/size1;
                        o_j = om1/size2;
                        o_k = om2/2;

                        xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                        xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                    }
                }
            }
        //---------------------------------------------------------------------
        // all done
        //---------------------------------------------------------------------
        G200:{} //continue;
            //if(timers_enabled)
            //    timer.stop(T_transxzloc);
        }

        public void transpose_x_z_global(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            //Fortran
            //double complex xin(d3,d2,d1);
            //double complex xout(d3,d2,d1) ! not real layout, but right size;
            //C#
            //xin[d1,d2,d3];
            //xout[d1,d2,d3]
            //if(timers_enabled)
            //    synchup();
            //---------------------------------------------------------------------
            // do transpose among all  processes with same 1-coord (me1)
            //---------------------------------------------------------------------
            //if(timers_enabled)
            //    timer.start(T_transxzglo);
            double[] src = new double[ntdivnp * 2];
            double[] dst = new double[ntdivnp * 2];
            setVetor(xin, src);
            commslice1.AlltoallFlattened<double>(src, d1*d2*d3*2/np2, ref dst);
            setVetor(dst, xout);
            //call mpi_alltoall(xin, d1*d2*d3/np2, dc_type,xout, d1*d2*d3/np2, dc_type,commslice1, ierr);
            //if(timers_enabled)
            //    timer.stop(T_transxzglo);
        }

/**/    public void transpose_x_z_finish(int d1, int d2, int d3, double[, , ,] xin, double[, , ,] xout) {
            //Fortran
            //double complex xin(d1/np2, d2, d3, 0:np2-1)
            //double complex xout(d1,d2,d3)
            //C#
            //xin  [np2, d3, d2, d1/np2]
            //xout [d3,d2,d1]
            int i, j, k, p, ioff, iin, io;
            int m1, m2, _i, _j, _k, om1, om2, o_i, o_j, o_k;

            //if(timers_enabled)
            //    timer.start(T_transxzfin);
            for(p = 0; p <= np2 - 1; p++) {
                ioff = p*d1/np2;
                for(k = 0; k < d3; k++) {
                    for(j = 0; j < d2; j++) {
                        for(i = 0; i < d1 / np2; i++) {                                   //xout(i+ioff, j, k) = xin(i, j, k, p);
                            iin = (((p*d3+k)*d2+j)*(d1/np2)+i)*2; //xin[p, k, j, i];
                            io  = ((k*d2+j)*d1+(i+ioff))*2; //xout[k, j, i+ioff]
                            //Point.setAddress(xin, iin + REAL, xout, io + REAL);
                            //Point.setAddress(xin, iin + IMAG, xout, io + IMAG);
                            m1 = (iin % size1);
                            m2 = (m1 % size2);
                            _i = iin/size1;
                            _j = m1/size2;
                            _k = m2/2;

                            om1 = (io % size1);
                            om2 = (om1 % size2);
                            o_i = io/size1;
                            o_j = om1/size2;
                            o_k = om2/2;

                            xout[o_i, o_j, o_k, REAL] = xin[_i, _j, _k, REAL];
                            xout[o_i, o_j, o_k, IMAG] = xin[_i, _j, _k, IMAG];
                        }
                    }
                }
            }
            //if(timers_enabled)
            //    timer.stop(T_transxzfin);
        }

        public void transpose_x_yz(int l1, int l2, double[, , ,] xin, double[, , ,] xout) {
            //double complex xin(ntdivnp), xout(ntdivnp)
            transpose_x_yz_local(dims[0, l1], dims[1, l1], dims[2, l1], xin, xout);
            transpose_x_yz_global(dims[1, 0], dims[2, 0], dims[0, 0], xout, xin);
            transpose_x_yz_finish(dims[0, l1], dims[1, l1],dims[2, l1], xin, xout);
        }

        public void evolve(double[, , ,] u0, double[, , ,] u1, double[] twiddle, int d1, int d2, int d3) {
            //  c---------------------------------------------------------------------
            //  c evolve u0 -> u1 (t time steps) in fourier space
            //  c---------------------------------------------------------------------
            //Fortran:
            //double complex u00(d1,d2,d3)
            //double complex u11(d1,d2,d3)
            //double precision twiddle1(d1,d2,d3)
            //C#:
            //double complex u00[d3,d2,d1]
            //double complex u11[d3,d2,d1]
            //double precision twiddle1[d3,d2,d1]

            int i, j, k, idx, m1, m2, _i, _j, _k;
            for(k = 0; k < d3; k++) {
                for(j = 0; j < d2; j++) {
                    for(i = 0; i < d1; i++) {//u00(i,j,k) = u00(i,j,k)*(twiddle1(i,j,k)); u11(i,j,k) = u00(i,j,k);//u00[k,j,i] = u00[k,j,i]*(twiddle1[k,j,i]);//u11[k,j,i] = u00[k,j,i];
                        idx = ((k*d2+j)*d1+i);

                        m1 = (idx*2 % size1);
                        m2 = (m1 % size2);
                        _i = idx*2/size1;
                        _j = m1/size2;
                        _k = m2/2;

                        u0[_i,_j,_k,REAL] = u0[_i,_j,_k,REAL]*twiddle[idx];
                        u0[_i,_j,_k,IMAG] = u0[_i,_j,_k,IMAG]*twiddle[idx];
                        u1[_i,_j,_k,REAL] = u0[_i,_j,_k,REAL];
                        u1[_i,_j,_k,IMAG] = u0[_i,_j,_k,IMAG];
                    }
                }
            }

        }

        public void checksum(int iter, double[] sums, double[, , ,] u2, int d1, int d2, int d3) {
            //Fortran:     double complex u1(d1, d2, d3);
            //C#     :     double complex u1[d3, d2, d1];
            int j, q,r,s, m1,m2,_i,_j,_k;
            double chk_Real, chk_Imag;
            double allchk_Real=0, allchk_Imag=0;   //double complex chk,allchk;

            chk_Real = 0.0;
            chk_Imag = 0.0;

            int idx=0;
            for(j = 1; j <= 1024; j++) {
                q = (int)mod(j, nx)+1;
                if(q >= xstart[0] && q <= xend[0]) {
                    r = (int)mod(3*j, ny)+1;
                    if(r >= ystart[0] && r <= yend[0]) {
                        s = (int)mod(5*j, nz)+1;
                        if(s >= zstart[0] && s <= zend[0]) {//chk=chk+u11(q-xstart(1)+1,r-ystart(1)+1,s-zstart(1)+1);
                            //C#     :     double complex u11[d3, d2, d1];
                            idx = (((s-zstart[0])*d2+(r-ystart[0]))*d1+(q-xstart[0]))*2;
                            m1 = (idx % size1);
                            m2 = (m1 % size2);
                            _i = idx/size1;
                            _j = m1/size2;
                            _k = m2/2;
                            chk_Real=chk_Real+u2[_i,_j,_k,REAL]; //u11[i1, i2, i3];//chk_Real=chk_Real+Point.getValue(u2, idx+REAL); //u11[i1, i2, i3];
                            chk_Imag=chk_Imag+u2[_i,_j,_k,IMAG]; //u11[i1, i2, i3];//chk_Imag=chk_Imag+Point.getValue(u2, idx+IMAG); //u11[i1, i2, i3];
                        }
                    }
                }
            }
            chk_Real = chk_Real/((double)(nx*ny*nz));
            chk_Imag = chk_Imag/((double)(nx*ny*nz));

            allchk_Real = worldcomm.Reduce<double>(chk_Real, MPI.Operation<double>.Add, root);
            allchk_Imag = worldcomm.Reduce<double>(chk_Imag, MPI.Operation<double>.Add, root);

            if(node == 0) {
                Console.WriteLine(" T = " + iter + "  Checksum = (" + allchk_Real + ") (" + allchk_Imag + ")");//+1P2D22.12);
            }
            if(iter >= 0) {
                sums[iter*2+REAL] = allchk_Real;
                sums[iter*2+IMAG] = allchk_Imag;
            }

        }

        public int verify(int d1, int d2, int d3, int nt, double[] sums) {
            int i;
            double err, epsilon;
            double[] csum_ref = new double[25*2]; //     double complex csum_ref(25); even=Real odd=Imag
            char Class = 'U';
            if(node != 0)
                return 0;
            epsilon = 0.000000000001;//epsilon = 1.0D-12;
            int verified = 0;
            if(d1 == 64 && d2 == 64 && d3 == 64 && nt == 6) {
                Class = 'S';
                csum_ref[0] = 554.6087004964;
                csum_ref[1] = 484.5363331978;
                csum_ref[2] = 554.6385409189;
                csum_ref[3] = 486.5304269511;
                csum_ref[4] = 554.6148406171;
                csum_ref[5] = 488.3910722336;
                csum_ref[6] = 554.5423607415;
                csum_ref[7] = 490.1273169046;
                csum_ref[8] = 554.4255039624;
                csum_ref[9] = 491.7475857993;
                csum_ref[10] = 554.2683411902;
                csum_ref[11] = 493.2597244941;
            }
            else if(d1 == 128 && d2 == 128 && d3 == 32 && nt == 6) {
                Class = 'W';
                csum_ref[0] = 567.3612178944;
                csum_ref[1] = 529.3246849175;
                csum_ref[2] = 563.1436885271;
                csum_ref[3] = 528.2149986629;
                csum_ref[4] = 559.4024089970;
                csum_ref[5] = 527.0996558037;
                csum_ref[6] = 556.0698047020;
                csum_ref[7] = 526.0027904925;
                csum_ref[8] = 553.0898991250;
                csum_ref[9] = 524.9400845633;
                csum_ref[10] = 550.4159734538;
                csum_ref[11] = 523.9212247086;
            }
            else if(d1 == 256 && d2 == 256 && d3 == 128 && nt == 6) {
                Class = 'A';
                csum_ref[0] = 504.6735008193;
                csum_ref[1] = 511.4047905510;
                csum_ref[2] = 505.9412319734;
                csum_ref[3] = 509.8809666433;
                csum_ref[4] = 506.9376896287;
                csum_ref[5] = 509.8144042213;
                csum_ref[6] = 507.7892868474;
                csum_ref[7] = 510.1336130759;
                csum_ref[8] = 508.5233095391;
                csum_ref[9] = 510.4914655194;
                csum_ref[10] = 509.1487099959;
                csum_ref[11] = 510.7917842803;
            }
            else if(d1 == 512 && d2 == 256 && d3 == 256 && nt == 20) {
                Class = 'B';
                csum_ref[0] = 517.7643571579;
                csum_ref[1] = 507.7803458597;
                csum_ref[2] = 515.4521291263;
                csum_ref[3] = 508.8249431599;
                csum_ref[4] = 514.6409228649;
                csum_ref[5] = 509.6208912659;
                csum_ref[6] = 514.2378756213;
                csum_ref[7] = 510.1023387619;
                csum_ref[8] = 513.9626667737;
                csum_ref[9] = 510.3976610617;
                csum_ref[10] = 513.7423460082;
                csum_ref[11] = 510.5948019802;
                csum_ref[12] = 513.5547056878;
                csum_ref[13] = 510.7404165783;
                csum_ref[14] = 513.3910925466;
                csum_ref[15] = 510.8576573661;
                csum_ref[16] = 513.2470705390;
                csum_ref[17] = 510.9577278523;
                csum_ref[18] = 513.1197729984;
                csum_ref[19] = 511.0460304483;
                csum_ref[20] = 513.0070319283;
                csum_ref[21] = 511.1252433800;
                csum_ref[22] = 512.9070537032;
                csum_ref[23] = 511.1968077718;
                csum_ref[24] = 512.8182883502;
                csum_ref[25] = 511.2616233064;
                csum_ref[26] = 512.7393733383;
                csum_ref[27] = 511.3203605551;
                csum_ref[28] = 512.6691062020;
                csum_ref[29] = 511.3735928093;
                csum_ref[30] = 512.6064276004;
                csum_ref[31] = 511.4218460548;
                csum_ref[32] = 512.5504076570;
                csum_ref[33] = 511.4656139760;
                csum_ref[34] = 512.5002331720;
                csum_ref[35] = 511.5053595966;
                csum_ref[36] = 512.4551951846;
                csum_ref[37] = 511.5415130407;
                csum_ref[38] = 512.4146770029;
                csum_ref[39] = 511.5744692211;
            }
            else if(d1 == 512 && d2 == 512 && d3 == 512 && nt == 20) {
                Class = 'C';
                csum_ref[0] = 519.5078707457;
                csum_ref[1] = 514.9019699238;
                csum_ref[2] = 515.5422171134;
                csum_ref[3] = 512.7578201997;
                csum_ref[4] = 514.4678022222;
                csum_ref[5] = 512.2251847514;
                csum_ref[6] = 514.0150594328;
                csum_ref[7] = 512.1090289018;
                csum_ref[8] = 513.7550426810;
                csum_ref[9] = 512.1143685824;
                csum_ref[10] = 513.5811056728;
                csum_ref[11] = 512.1496764568;
                csum_ref[12] = 513.4569343165;
                csum_ref[13] = 512.1870921893;
                csum_ref[14] = 513.3651975661;
                csum_ref[15] = 512.2193250322;
                csum_ref[16] = 513.2955192805;
                csum_ref[17] = 512.2454735794;
                csum_ref[18] = 513.2410471738;
                csum_ref[19] = 512.2663649603;
                csum_ref[20] = 513.1971141679;
                csum_ref[21] = 512.2830879827;
                csum_ref[22] = 513.1605205716;
                csum_ref[23] = 512.2965869718;
                csum_ref[24] = 513.1290734194;
                csum_ref[25] = 512.3075927445;
                csum_ref[26] = 513.1012720314;
                csum_ref[27] = 512.3166486553;
                csum_ref[28] = 513.0760908195;
                csum_ref[29] = 512.3241541685;
                csum_ref[30] = 513.0528295923;
                csum_ref[31] = 512.3304037599;
                csum_ref[32] = 513.0310107773;
                csum_ref[33] = 512.3356167976;
                csum_ref[34] = 513.0103090133;
                csum_ref[35] = 512.3399592211;
                csum_ref[36] = 512.9905029333;
                csum_ref[37] = 512.3435588985;
                csum_ref[38] = 512.9714421109;
                csum_ref[39] = 512.3465164008;
            }
            //else if(d1 == 2048 && d2 == 1024 && d3 == 1024 && nt == 25) {
            //    Class = 'D';
            //    csum_ref[0*2+REAL] = (512.2230065252);
            //    csum_ref[1*2+REAL] = (512.0463975765);
            //    csum_ref[2*2+REAL] = (511.9865766760);
            //    csum_ref[3*2+REAL] = (511.9518799488);
            //    csum_ref[4*2+REAL] = (511.9269088223);
            //    csum_ref[5*2+REAL] = (511.9082416858);
            //    csum_ref[6*2+REAL] = (511.8943814638);
            //    csum_ref[7*2+REAL] = (511.8842385057);
            //    csum_ref[8*2+REAL] = (511.8769435632);
            //    csum_ref[9*2+REAL] = (511.8718203448);
            //    csum_ref[10*2+REAL] = (511.8683569061);
            //    csum_ref[11*2+REAL] = (511.8661708593);
            //    csum_ref[12*2+REAL] = (511.8649768950);
            //    csum_ref[13*2+REAL] = (511.8645605626);
            //    csum_ref[14*2+REAL] = (511.8647586618);
            //    csum_ref[15*2+REAL] = (511.8654451572);
            //    csum_ref[16*2+REAL] = (511.8665212451);
            //    csum_ref[17*2+REAL] = (511.8679083821);
            //    csum_ref[18*2+REAL] = (511.8695433664);
            //    csum_ref[19*2+REAL] = (511.8713748264);
            //    csum_ref[20*2+REAL] = (511.8733606701);
            //    csum_ref[21*2+REAL] = (511.8754661974);
            //    csum_ref[22*2+REAL] = (511.8776626738);
            //    csum_ref[23*2+REAL] = (511.8799262314);
            //    csum_ref[24*2+REAL] = (511.8822370068);

            //    csum_ref[0*2+IMAG] = (511.8534037109);
            //    csum_ref[1*2+IMAG] = (511.7061181082);
            //    csum_ref[2*2+IMAG] = (511.7096364601);
            //    csum_ref[3*2+IMAG] = (511.7373863950);
            //    csum_ref[4*2+IMAG] = (511.7680347632);
            //    csum_ref[5*2+IMAG] = (511.7967875532);
            //    csum_ref[6*2+IMAG] = (511.8225281841);
            //    csum_ref[7*2+IMAG] = (511.8451629348);
            //    csum_ref[8*2+IMAG] = (511.8649119387);
            //    csum_ref[9*2+IMAG] = (511.8820803844);
            //    csum_ref[10*2+IMAG] = (511.8969781011);
            //    csum_ref[11*2+IMAG] = (511.9098918835);
            //    csum_ref[12*2+IMAG] = (511.9210777066);
            //    csum_ref[13*2+IMAG] = (511.9307604484);
            //    csum_ref[14*2+IMAG] = (511.9391362671);
            //    csum_ref[15*2+IMAG] = (511.9463757241);
            //    csum_ref[16*2+IMAG] = (511.9526269238);
            //    csum_ref[17*2+IMAG] = (511.9580184108);
            //    csum_ref[18*2+IMAG] = (511.9626617538);
            //    csum_ref[19*2+IMAG] = (511.9666538138);
            //    csum_ref[20*2+IMAG] = (511.9700787219);
            //    csum_ref[21*2+IMAG] = (511.9730095953);
            //    csum_ref[22*2+IMAG] = (511.9755100241);
            //    csum_ref[23*2+IMAG] = (511.9776353561);
            //    csum_ref[24*2+IMAG] = (511.9794338060);

            //}
            //else if(d1 == 4096 && d2 == 2048 && d3 == 2048 && nt == 25) {
            //    Class = 'E';
            //    csum_ref[0*2+REAL] = 512.1601045346;
            //    csum_ref[1*2+REAL] = 512.0905403678;
            //    csum_ref[2*2+REAL] = 512.0623229306;
            //    csum_ref[3*2+REAL] = 512.0438418997;
            //    csum_ref[4*2+REAL] = 512.0311521872;
            //    csum_ref[5*2+REAL] = 512.0226088809;
            //    csum_ref[6*2+REAL] = 512.0169296534;
            //    csum_ref[7*2+REAL] = 512.0131225172;
            //    csum_ref[8*2+REAL] = 512.0104767108;
            //    csum_ref[9*2+REAL] = 512.0085127969;
            //    csum_ref[10*2+REAL] = 512.0069224127;
            //    csum_ref[11*2+REAL] = 512.0055158164;
            //    csum_ref[12*2+REAL] = 512.0041820159;
            //    csum_ref[13*2+REAL] = 512.0028605402;
            //    csum_ref[14*2+REAL] = 512.0015223011;
            //    csum_ref[15*2+REAL] = 512.0001570022;
            //    csum_ref[16*2+REAL] = 511.9987650555;
            //    csum_ref[17*2+REAL] = 511.9973525091;
            //    csum_ref[18*2+REAL] = 511.9959279472;
            //    csum_ref[19*2+REAL] = 511.9945006558;
            //    csum_ref[20*2+REAL] = 511.9930795911;
            //    csum_ref[21*2+REAL] = 511.9916728462;
            //    csum_ref[22*2+REAL] = 511.9902874185;
            //    csum_ref[23*2+REAL] = 511.9889291565;
            //    csum_ref[24*2+REAL] = 511.9876028049;

            //    csum_ref[0*2+IMAG] = 511.7395998266;
            //    csum_ref[1*2+IMAG] = 511.8614716182;
            //    csum_ref[2*2+IMAG] = 511.9074203747;
            //    csum_ref[3*2+IMAG] = 511.9345900733;
            //    csum_ref[4*2+IMAG] = 511.9551325550;
            //    csum_ref[5*2+IMAG] = 511.9720179919;
            //    csum_ref[6*2+IMAG] = 511.9861371665;
            //    csum_ref[7*2+IMAG] = 511.9979364402;
            //    csum_ref[8*2+IMAG] = 512.0077674092;
            //    csum_ref[9*2+IMAG] = 512.0159443121;
            //    csum_ref[10*2+IMAG] = 512.0227453670;
            //    csum_ref[11*2+IMAG] = 512.0284096041;
            //    csum_ref[12*2+IMAG] = 512.0331373793;
            //    csum_ref[13*2+IMAG] = 512.0370938679;
            //    csum_ref[14*2+IMAG] = 512.0404138831;
            //    csum_ref[15*2+IMAG] = 512.0432068837;
            //    csum_ref[16*2+IMAG] = 512.0455615860;
            //    csum_ref[17*2+IMAG] = 512.0475499442;
            //    csum_ref[18*2+IMAG] = 512.0492304629;
            //    csum_ref[19*2+IMAG] = 512.0506508902;
            //    csum_ref[20*2+IMAG] = 512.0518503782;
            //    csum_ref[21*2+IMAG] = 512.0528612016;
            //    csum_ref[22*2+IMAG] = 512.0537101195;
            //    csum_ref[23*2+IMAG] = 512.0544194514;
            //    csum_ref[24*2+IMAG] = 512.0550079284;

            //}
            double a, b, c, d, r1, r2;
            if(Class != 'U') {
                for(i = 0; i < nt; i++) {
                    a = sums[i*2+REAL] - csum_ref[i*2+REAL];
                    b = sums[i*2+IMAG] - csum_ref[i*2+IMAG];
                    c = csum_ref[i*2+REAL];
                    d = csum_ref[i*2+IMAG];
                    r1 = ((a * c + b * d) / ((c * c) + (d * d))) * ((a * c + b * d) / ((c * c) + (d * d)));
                    r2 = ((c * b - a * d) / (c * c + d * d)) * ((c * b - a * d) / (c * c + d * d));
                    err = Math.Sqrt(r1 + r2);
                    if(!(err <= epsilon))
                        goto Go100;
                }
                verified = 1;
            }
        Go100:
            if(worldcomm.Size != np) {
                Console.WriteLine(" Warning: benchmark was compiled for "+np+" processors");
                Console.WriteLine(" Must be run on this many processors for official verification");
                Console.WriteLine(" so memory access is repeatable");
                verified = 0;
            }
            if(Class != 'U') {
                if(verified==1) {
                    Console.WriteLine(" Result verification successful");
                }
                else {
                    Console.WriteLine(" Result verification failed");
                }
            }
            Console.WriteLine("Class = "+ Class);
            return verified;
        }

        public void print_timers() {
            int i;
            string[] tstrings = new string[T_max];
            tstrings[0]="          total ";
            tstrings[1]="          setup ";
            tstrings[2]="            fft ";
            tstrings[3]="         evolve ";
            tstrings[4]="       checksum ";
            tstrings[5]="         fftlow ";
            tstrings[6]="        fftcopy ";
            tstrings[7]="      transpose ";
            tstrings[8]=" transpose1_loc ";
            tstrings[9]=" transpose1_glo ";
            tstrings[10]=" transpose1_fin ";
            tstrings[11]=" transpose2_loc ";
            tstrings[12]=" transpose2_glo ";
            tstrings[13]=" transpose2_fin ";
            tstrings[14]="           sync ";

            if(node != 0)
                return;
            for(i = 1; i <= T_max; i++) {
                if(timer.readTimer(i) != 0.0) {
                    Console.WriteLine(" timer "+ i + tstrings[i-1] + timer.readTimer(i));
                }
            }
        }

        public static unsafe void setVetor(double[, , ,] s, double[] d) {
            int size = s.Length;
            if(size == d.Length) {
                fixed(double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    for(int n = 0; n < size/2; n++) {
                        *((decimal*)p2) = *((decimal*)p1);
                        p2 += 2;
                        p1 += 2;
                    }
                }
            }
            else {
                throw new IndexOutOfRangeException();
            }
        }

        public static unsafe void setVetor(double[] s, double[, , ,] d) {
            int size = s.Length;
            if(size == d.Length) {
                fixed(double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    for(int n = 0; n < size / 2; n++) {
                        *((decimal*)p2) = *((decimal*)p1);
                        p2 += 2;
                        p1 += 2;
                    }
                }
            }
            else {
                throw new IndexOutOfRangeException();
            }
        }
    }
}
