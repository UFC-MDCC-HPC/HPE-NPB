using System;
using System.IO;
using NPB.Ftb;
using NPB.Ftc;

namespace NPB {
    public class FT:FTBase {

        public FT(char c):base(c){}

        static void Main(String[] argv){

            FT ft = null;

            try { 
                string param = argv[0]; 
            } catch (Exception) {
                argv = new String[1];
                argv[0] = "DEBUG";
            }
            char paramClass;
            if (argv[0] != "DEBUG") {
                BMArgs.ParseCmdLineArgs(argv, BMName);
                paramClass = BMArgs.CLASS;
            } else {
                paramClass = 'S';  //k=(S e 4 para np) S=(S e 1 para np)
            }

            try {
                ft = new FT(paramClass);
            } catch (OutOfMemoryException e) {
                Console.WriteLine(e.ToString());
                Environment.Exit(0);
            }

            ft.runBenchMark();

        }

        //**** codigos para debugar************************************************************
        public void debugImprimeFile1(double[] vetor, int nodo) {
            if (me == nodo) {
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/CSharpTwiddle" + clss + np + ".txt";
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                int p1; //p2, p3, p4; bool flag = true;
                for (p1 = 0; p1 < vetor.GetLength(0); p1++) {
                    arquivo.WriteLine("[" + p1 + "]=" + vetor[p1]);
                }
                arquivo.Close();
            }
        }

        public void debugImprimeFile2(double[] vetor, int nodo) {
            if (me == nodo) {
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/CSharpU1" + clss + np + ".txt";
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                int p1; //p2, p3, p4; bool flag = true;
                for (p1 = 0; p1 < vetor.GetLength(0); p1++) {
                    if (p1 % 2 == 0) {
                        arquivo.Write("[" + p1 + "]=" + vetor[p1]);
                    }
                    else {
                        arquivo.WriteLine(" :: [" + p1 + "]=" + vetor[p1]);
                    }

                }
                arquivo.Close();
            }
        }

        public void debugImprimeFile3(double[] vetor, int nodo) {
            if (me == nodo) {
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/CSharp-U-" + clss + np + ".txt";
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                int p1; //p2, p3, p4; bool flag = true;
                for (p1 = 0; p1 < vetor.GetLength(0); p1++) {
                    arquivo.WriteLine("[" + p1 + "]=" + vetor[p1]);

                }
                arquivo.Close();
            }
        }

        //*******************************************************************************************

        public void startBigArrays() {
            //complex u0(ntdivnp), u1(ntdivnp), u2(ntdivnp)
            u0 = new double[ntdivnp*2];
            u1 = new double[ntdivnp*2];
            u2 = new double[ntdivnp*2];
            u  = new double[nx,2];
            uUnidimensional = new double[nx*2];

            // twiddle(ntdivnp)
            twiddle = new double[ntdivnp];    
        }

        public void runBenchMark() {
            for (i = 1; i <= T_max; i++) {
                timer_clear(i);
            }
            setup();
            startBigArrays();

            compute_indexmap(twiddle, dims[1,3], dims[2,3], dims[3,3]); 
            //debugImprimeFile1(twiddle,root);

            compute_initial_conditions(u1, dims[1, 1], dims[2, 1], dims[3, 1]);
            //debugImprimeFile2(u1, root);

            fft_init(dims[1, 1]);  //control u
            //call fft(1, u1, u0)

            //c---------------------------------------------------------------------
            //c Start over from the beginning. Note that all operations must
            //c be timed, in contrast to other benchmarks. 
            //c---------------------------------------------------------------------
            //do i = 1, t_max
            //   call timer_clear(i)
            //end do
            //call MPI_Barrier(MPI_COMM_WORLD, ierr)

            //call timer_start(T_total)
            //if (timers_enabled) call timer_start(T_setup)

            //call compute_indexmap(twiddle, dims(1,3), dims(2,3), dims(3,3))
            //call compute_initial_conditions(u1, dims(1,1), dims(2,1), dims(3,1))
            //call fft_init (dims(1,1))

            //if (timers_enabled) call synchup()
            //if (timers_enabled) call timer_stop(T_setup)

            //if (timers_enabled) call timer_start(T_fft)
            //call fft(1, u1, u0)
            //if (timers_enabled) call timer_stop(T_fft)

            //do iter = 1, niter
            //   if (timers_enabled) call timer_start(T_evolve)
            //   call evolve(u0, u1, twiddle, dims(1,1), dims(2,1), dims(3,1))
            //   if (timers_enabled) call timer_stop(T_evolve)
            //   if (timers_enabled) call timer_start(T_fft)
            //   call fft(-1, u1, u2)
            //   if (timers_enabled) call timer_stop(T_fft)
            //   if (timers_enabled) call synchup()
            //   if (timers_enabled) call timer_start(T_checksum)
            //   call checksum(iter, u2, dims(1,1), dims(2,1), dims(3,1))
            //   if (timers_enabled) call timer_stop(T_checksum)
            //end do

            //call verify(nx, ny, nz, niter, verified, class)
            //call timer_stop(t_total)
            //if (np .ne. np_min) verified = .false.
            //total_time = timer_read(t_total)

            //if( total_time .ne. 0. ) then
            //   mflops = 1.0d-6*ntotal_f * (14.8157+7.19641*log(ntotal_f) +  (5.23518+7.21113*log(ntotal_f))*niter)/total_time
            //else
            //   mflops = 0.0
            //endif
            //if (me .eq. 0) then
            //   call print_results('FT', class, nx, ny, nz, niter, np_min, np,total_time, mflops, '          floating point', verified, npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)
            //endif
            //if (timers_enabled) call print_timers()

            mpi.Dispose();
        }

        public void setup() {

            //common /procgrid/ np1, np2, layout_type, np
            //common /blockinfo/ fftblock, fftblockpad
            //common /coords/ me, me1, me2
            //common /comms/ commslice1, commslice2
            //common /layout/ dims,xstart, ystart, zstart, xend, yend, zend
            //common /ucomm/ u
            //common /dbg/ debug, debugsynch
            //common /sumcomm/ sums
            //common /iter/ niter

            //common /mpistuff/ dc_type

            int i, fstatus=0; //ierr, j;
            debug = false;

            if (!convertdouble) {
                dc_type = 1275072546; //MPI_DOUBLE_COMPLEX;
            } else {
                dc_type = 1275070494; //MPI_COMPLEX;
            }

            if (me == 0) {
                Console.WriteLine(" NAS Parallel Benchmarks "+ npbversion +" -- FT Benchmark ");
                //        open (unit=2,file='inputft.data',status='old', iostat=fstatus)
                try {
                    Console.Write("Trying Read from input file inputft.data: ");
                    int[] vetTemp = BMArgs.readInputFtData("inputft.data");
                    niter = vetTemp[0]; layout_type = vetTemp[1]; np1 = vetTemp[2]; np2 = vetTemp[3];
                }
                catch (System.IO.FileNotFoundException e) {
                    Console.WriteLine("inputft.data not found");
                    fstatus = 1;
                    Console.WriteLine(e.ToString());
                }

                if (fstatus == 0) {
                    Console.WriteLine("inputft.data found");

                    //c---------------------------------------------------------------------
                    //c check to make sure input data is consistent
                    //c---------------------------------------------------------------------
                    //c---------------------------------------------------------------------
                    //c 1. product of processor grid dims must equal number of processors
                    //c---------------------------------------------------------------------

                    if (np1 * np2 != np) {
                        Console.WriteLine(" np1 and np2 given in input file are not valid.");
                        Console.WriteLine("Product is "+ np1*np2+" and should be "+np);
                        mpi.Dispose();
                        Environment.Exit(0);
                    }

                    //c---------------------------------------------------------------------
                    //c 2. layout type must be valid
                    //c---------------------------------------------------------------------

                    if (layout_type != layout_0D && layout_type != layout_1D && layout_type != layout_2D) {
                        Console.WriteLine(" Layout type specified in inputft.data is invalid ");
                        mpi.Dispose();
                        Environment.Exit(0);
                    }

                    //c---------------------------------------------------------------------
                    //c 3. 0D layout must be 1x1 grid
                    //c---------------------------------------------------------------------

                    if (layout_type == layout_0D && (np1 != 1 || np2 != 1)) {
                        Console.WriteLine(" For 0D layout, both np1 and np2 must be 1 ");
                        mpi.Dispose();
                        Environment.Exit(0);
                    }
                    //c---------------------------------------------------------------------
                    //c 4. 1D layout must be 1xN grid
                    //c---------------------------------------------------------------------

                    if (layout_type == layout_1D && np1 != 1) {
                        Console.WriteLine(" For 1D layout, np1 must be 1 ");
                        mpi.Dispose();
                        Environment.Exit(0);
                    }
                }
                else {
                    Console.WriteLine(" No input file inputft.data. Using compiled defaults"); 
                    niter = niter_default;
                    if (np == 1) {
                        np1 = 1;
                        np2 = 1;
                        layout_type = layout_0D;
                    }
                    else if (np <= nz) {
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

                if (np < np_min) {
                    Console.WriteLine(" Error: Compiled for "+ np_min + " processors. ");
                    Console.WriteLine(" Only "+ np + " processors found ");
                    mpi.Dispose();
                    Environment.Exit(0);
                }
                Console.WriteLine(" Size: " + nx + "x" + ny + "x" + nz);
                Console.WriteLine(" Iterations: "+ niter);
                Console.WriteLine(" Number of processes : "+ np);
                Console.WriteLine(" Processor array     : "+np1+"x"+np2);
                if (np != np_min) Console.WriteLine(" WARNING: compiled for "+np_min+" processes. Will not verify. ");
                if (layout_type == layout_0D) {
                    Console.WriteLine(" Layout type: OD");
                }
                else if (layout_type == layout_1D) {
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
            worldcomm.Broadcast<int>(ref layout_type, root);

            if (np1 == 1 && np2 == 1) {
                layout_type = layout_0D;
            }
            else if (np1 == 1) {
                layout_type = layout_1D;
            }
            else {
                layout_type = layout_2D;
            }

            if (layout_type == layout_0D) {
                for (i = 1; i <= 3; i++) {
                    dims[1, i] = nx;
                    dims[2, i] = ny;
                    dims[3, i] = nz;
                }
            } else if (layout_type == layout_1D) {
                dims[1, 1] = nx;
                dims[2, 1] = ny;
                dims[3, 1] = nz;

                dims[1, 2] = nx;
                dims[2, 2] = ny;
                dims[3, 2] = nz;

                dims[1, 3] = nz;
                dims[2, 3] = nx;
                dims[3, 3] = ny;
            }
            else if (layout_type == layout_2D) {
                    dims[1, 1] = nx;
                    dims[2, 1] = ny;
                    dims[3, 1] = nz;

                    dims[1, 2] = ny;
                    dims[2, 2] = nx;
                    dims[3, 2] = nz;

                    dims[1, 3] = nz;
                    dims[2, 3] = nx;
                    dims[3, 3] = ny;

            }

            for (i = 1; i <= 3; i++) {
                dims[2, i] = dims[2, i] / np1;
                dims[3, i] = dims[3, i] / np2;
            }


            //c---------------------------------------------------------------------
            //c Determine processor coordinates of this processor
            //c Processor grid is np1xnp2. 
            //c Arrays are always (n1, n2/np1, n3/np2)
            //c Processor coords are zero-based. 
            //c---------------------------------------------------------------------

            me2 = (int) mod(me, np2);  // goes from 0...np2-1
            me1 = me/np2;        // goes from 0...np1-1

            //c---------------------------------------------------------------------
            //c Communicators for rows/columns of processor grid. 
            //c commslice1 is communicator of all procs with same me1, ranked as me2
            //c commslice2 is communicator of all procs with same me2, ranked as me1
            //c mpi_comm_split(comm, color, key, ...)
            //c---------------------------------------------------------------------

            commslice1=(MPI.Intracommunicator)worldcomm.Split(me1,me2);//MPI_Comm_split(MPI_COMM_WORLD,me1,me2,commslice1,ierr)
            commslice2=(MPI.Intracommunicator)worldcomm.Split(me2,me1);//MPI_Comm_split(MPI_COMM_WORLD,me2,me1,commslice2,ierr)

            if (timers_enabled) synchup();

            if (debug) Console.WriteLine("proc coords: " + me +" "+ me1 +" "+ me2);

            //c---------------------------------------------------------------------
            //c Determine which section of the grid is owned by this
            //c processor. 
            //c---------------------------------------------------------------------
            if (layout_type == layout_0D) {

                for (i = 1; i <= 3; i++) {
                    xstart[i] = 1;
                    xend[i]   = nx;
                    ystart[i] = 1;
                    yend[i]   = ny;
                    zstart[i] = 1;
                    zend[i]   = nz;
                }

            } else if (layout_type == layout_1D) {
                xstart[1] = 1;
                xend[1]   = nx;
                ystart[1] = 1;
                yend[1]   = ny;
                zstart[1] = 1 + me2 * nz/np2;
                zend[1]   = (me2+1) * nz/np2;

                xstart[2] = 1;
                xend[2]   = nx;
                ystart[2] = 1;
                yend[2]   = ny;
                zstart[2] = 1 + me2 * nz/np2;
                zend[2]   = (me2+1) * nz/np2;

                xstart[3] = 1;
                xend[3]   = nx;
                ystart[3] = 1 + me2 * ny/np2;
                yend[3]   = (me2+1) * ny/np2;
                zstart[3] = 1;
                zend[3] = nz;

            }
            else if (layout_type == layout_2D) {

                xstart[1] = 1;
                xend[1]   = nx;
                ystart[1] = 1 + me1 * ny/np1;
                yend[1]   = (me1+1) * ny/np1;
                zstart[1] = 1 + me2 * nz/np2;
                zend[1]   = (me2+1) * nz/np2;

                xstart[2] = 1 + me1 * nx/np1;
                xend[2]   = (me1+1)*nx/np1;
                ystart[2] = 1;
                yend[2]   = ny;
                zstart[2] = zstart[1];
                zend[2]   = zend[1];

                xstart[3] = xstart[2];
                xend[3]   = xend[2];
                ystart[3] = 1 + me2 *ny/np2;
                yend[3]   = (me2+1)*ny/np2;
                zstart[3] = 1;
                zend[3] = nz;
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

            if (layout_type == layout_2D) {
                if (dims[2, 1] < fftblock) fftblock = dims[2, 1];
                if (dims[2, 2] < fftblock) fftblock = dims[2, 2];
                if (dims[2, 3] < fftblock) fftblock = dims[2, 3];
            }

            if (fftblock != fftblock_default) fftblockpad = fftblock + 3;
        }

        public static double mod(double a, double b) {
            return (a % b);
        }

        public void synchup(){
            //common /procgrid/ np1, np2, layout_type, np
            //common /blockinfo/ fftblock, fftblockpad
            //common /coords/ me, me1, me2
            //common /comms/ commslice1, commslice2
            //common /layout/ dims,xstart, ystart, zstart, xend, yend, zend
            //common /ucomm/ u
            //common /dbg/ debug, debugsynch
            //common /sumcomm/ sums
            //common /iter/ niter

            //common /mpistuff/ dc_type
            timer_start(T_synch);
            worldcomm.Barrier();
            timer_stop(T_synch);
        }

        public static void timer_clear(int n) {
            elapsed[n] = 0.0;
        }

        public static void timer_start(int n) {
            start[n] = MPI.Unsafe.MPI_Wtime();
        }

        public static void timer_stop(int n) {
            double t, now;
            now = MPI.Unsafe.MPI_Wtime();
            t = now - start[n];
            elapsed[n] = elapsed[n] + t;
        }

        public static double timer_readGet(int n) { //Note: timer_read: There is a variable static with this name, so this
                                                    // function name renamed to timer_readGet.
            return elapsed[n];
        }

        public void compute_indexmap(double[] twiddle, int d1, int d2, int d3) {

            //c---------------------------------------------------------------------
            //c compute function from local (i,j,k) to ibar^2+jbar^2+kbar^2 
            //c for time evolution exponent. 
            //c---------------------------------------------------------------------

            //common /procgrid/ np1, np2, layout_type, np
            //common /blockinfo/ fftblock, fftblockpad
            //common /coords/ me, me1, me2
            //common /comms/ commslice1, commslice2
            //common /layout/ dims,xstart, ystart, zstart, xend, yend, zend
            //common /ucomm/ u
            //common /dbg/ debug, debugsynch
            //common /sumcomm/ sums
            //common /iter/ niter

            //common /mpistuff/ dc_type

            int i, j, k, ii, ii2, jj, ij2, kk;
            double ap; 
            //twiddle(d1, d2, d3)

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

            ap = -4.0 * alpha * pi * pi;

            if (layout_type == layout_0D) { //xyz layout
                int size1 = dims[1, 3]; int size2 = dims[2, 3]; int size3 = dims[3, 3]; int idx;
                for (i = 1; i <= size1; i++) {
                    ii =  (int) mod(i+xstart[3]-2+nx/2, nx) - nx/2;
                    ii2 = ii*ii;
                    for (j = 1; j <= size2; j++) {
                        jj = (int) mod(j+ystart[3]-2+ny/2, ny) - ny/2;
                        ij2 = jj * jj + ii2;
                        for (k = 1; k <= size3; k++) {
                            kk = (int) mod(k+zstart[3]-2+nz/2, nz) - nz/2;
                            //twiddle[i, j, k] = dexp(ap * dfloat(kk * kk + ij2));
                            idx = (((i-1)*size2+(j-1))*size3+(k-1));
                            twiddle[idx] = Math.Exp(ap * (double)(kk * kk + ij2));
                        }
                    }
                }
            } else if (layout_type == layout_1D) { // zxy layout 
                int size1 = dims[2, 3]; int size2 = dims[3, 3]; int size3 = dims[1, 3]; int idx;
                for (i = 1; i <= size1; i++) {
                    ii =  (int) mod(i+xstart[3]-2+nx/2, nx) - nx/2;
                    ii2 = ii*ii;
                    for (j = 1; j <= size2; j++) {
                        jj = (int) mod(j+ystart[3]-2+ny/2, ny) - ny/2;
                        ij2 = jj*jj+ii2;
                        for (k = 1; k <= size3; k++) {
                            kk = (int) mod(k+zstart[3]-2+nz/2, nz) - nz/2;
                            //twiddle(k, i, j) = dexp(ap * dfloat(kk * kk + ij2));
                            idx = (((i - 1) * size2 + (j - 1)) * size3 + (k - 1));
                            twiddle[idx] = Math.Exp(ap * (double)(kk * kk + ij2));
                        }
                    }
                }
            }
            else if (layout_type == layout_2D) { // zxy layout
                int size1 = dims[2, 3]; int size2 = dims[3, 3]; int size3 = dims[1, 3]; int idx;
                for (i = 1; i <= size1; i++) {
                    ii =  (int) mod(i+xstart[3]-2+nx/2, nx) - nx/2;
                    ii2 = ii*ii;
                    for (j = 1; j <= size2; j++) {
                        jj = (int) mod(j+ystart[3]-2+ny/2, ny) - ny/2;
                        ij2 = jj*jj+ii2;
                        for (k = 1; k <= size3; k++) {
                            kk = (int) mod(k+zstart[3]-2+nz/2, nz) - nz/2;
                            //twiddle(k,i,j) = dexp(ap*dfloat(kk*kk+ij2)); 
                            idx = (((i - 1) * size2 + (j - 1)) * size3 + (k - 1));
                            twiddle[idx] = Math.Exp(ap * (double)(kk*kk+ij2));
                        }
                    }
                }
            }
            else {
                Console.WriteLine(" Unknown layout type " + layout_type);
            }
        }

        public void compute_initial_conditions(double[] u0, int d1, int d2, int d3) {
            //c---------------------------------------------------------------------
            //c Fill in array u0 with initial conditions from 
            //c random number generator 
            //c---------------------------------------------------------------------

            //common /procgrid/ np1, np2, layout_type, np
            //common /blockinfo/ fftblock, fftblockpad
            //common /coords/ me, me1, me2
            //common /comms/ commslice1, commslice2
            //common /layout/ dims,xstart, ystart, zstart, xend, yend, zend
            //common /ucomm/ u
            //common /dbg/ debug, debugsynch
            //common /sumcomm/ sums
            //common /iter/ niter

            //double complex u0(d1, d2, d3)

            int k;
            double x0, start, an, dummy;

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

            an = ipow46(a, 2*nx, (zstart[1]-1)*ny + (ystart[1]-1));
            dummy = randlcGet(ref start, ref an);
            an = ipow46(a, 2*nx, ny);

            //c---------------------------------------------------------------------
            //c Go through by z planes filling in one square at a time.
            //c---------------------------------------------------------------------
            for (k = 1; k <= dims[3, 1]; k++) { // nz/np2;
                x0 = start;
                //call vranlc(2*nx*dims(2, 1), x0, a, u0(1, 1, k)) : call native
                vranlc(2 * nx * dims[2, 1], x0, a, u0, k); //(1, 1, k));
                if (k != dims[3, 1]) dummy = randlcGet(ref start, ref an);
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
            if (exp_2 == 0 || exp_1 == 0) return result;
            q = a;
            r = 1;
            n = exp_1;
            two_pow = true;

            while (two_pow) {
                n2 = n/2;
                if (n2 * 2 == n) {
                    dummy = randlcGet(ref q, ref q);
                    n = n2;
                } else {
                    n = n * exp_2;
                    two_pow = false;
                }
            }

            while (n > 1) {
                n2 = n/2;
                if (n2 * 2 == n) {
                    dummy = randlcGet(ref q, ref q);
                    n = n2;
                } else {
                    dummy = randlcGet(ref r, ref q);
                    n = n-1;
                }
            }
            dummy = randlcGet(ref r, ref q);
            result = r;
            return result;
        }


        public double randlcGet(ref double x, ref double a) {
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
            a1 = (int) (t1);
            a2 = a - t23 * a1;

            //c---------------------------------------------------------------------
            //c   Break X into two parts such that X = 2^23 * X1 + X2, compute
            //c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
            //c   X = 2^23 * Z + A2 * X2  (mod 2^46).
            //c---------------------------------------------------------------------
            t1 = r23 * x;
            x1 = (int) (t1);
            x2 = x - t23 * x1;


            t1 = a1 * x2 + a2 * x1;
            t2 = (int) (r23 * t1);
            z = t1 - t23 * t2;
            t3 = t23 * z + a2 * x2;
            t4 = (int) (r46 * t3);
            x = t3 - t46 * t4;
            return (r46*x);
        }

        public void vranlc(int n, double x, double a, double[] y, int k) {
            //     c---------------------------------------------------------------------
            //     c   This routine generates N uniform pseudorandom double precision numbers in
            //     c   the range (0, 1) by using the linear congruential generator
            //     c   
            //     c   x_{k+1} = a x_k  (mod 2^46)
            //     c   
            //     c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
            //     c   before repeating.  The argument A is the same as 'a' in the above formula,
            //     c   and X is the same as x_0.  A and X must be odd double precision integers
            //     c   in the range (1, 2^46).  The N results are placed in Y and are normalized
            //     c   to be between 0 and 1.  X is updated to contain the new seed, so that
            //     c   subsequent calls to RANDLC using the same arguments will generate a
            //     c   continuous sequence.
            //     c   
            //     c   This routine generates the output sequence in batches of length NV, for
            //     c   convenience on vector computers.  This routine should produce the same
            //     c   results on any computer with at least 48 mantissa bits in double precision
            //     c   floating point data.  On Cray systems, double precision should be disabled.
            //     c   
            //     c   David H. Bailey    August 30, 1990
            //     c---------------------------------------------------------------------
            double r23, r46, t23, t46;
            int nv;
            r23 = Math.Pow(2.0,(-23)); 
            r46 = r23 * r23;
            t23 = Math.Pow(2.0,23); 
            t46 = t23 * t23; 
            nv = 64;
            double[] xv = new double[nv];
            double t1, t2, t3, t4, an, a1=0, a2=0, x1, x2, yy;
            int n1, i, j;
            //     c---------------------------------------------------------------------
            //     c     Compute the first NV elements of the sequence using RANDLC.
            //     c---------------------------------------------------------------------
            t1 = x;
            n1 = (int) min(n, nv);

            for (i = 1; i <= n1; i++) {
                xv[i-1] = t46 * randlcGet(ref t1, ref a);
            }

            //     c---------------------------------------------------------------------
            //     c     It is not necessary to compute AN, A1 or A2 unless N is greater than NV.
            //     c---------------------------------------------------------------------
            if (n > nv) {

                //c---------------------------------------------------------------------
                //c     Compute AN = AA ^ NV (mod 2^46) using successive calls to RANDLC.
                //c---------------------------------------------------------------------
                t1 = a;
                t2 = r46 * a;

                for (i = 1; i <= nv - 1; i++) {
                    t2 = randlcGet(ref t1, ref a);
                }

                an = t46 * t2;

                //c---------------------------------------------------------------------
                //c     Break AN into two parts such that AN = 2^23 * A1 + A2.
                //c---------------------------------------------------------------------
                t1 = r23 * an;
                a1 = (int) (t1);
                a2 = an - t23 * a1;
            }

            //     c---------------------------------------------------------------------
            //     c     Compute N pseudorandom results in batches of size NV.
            //     c---------------------------------------------------------------------
            for (j = 0; j <= n - 1; j = j + nv) {
                n1 = (int) min(nv, n - j);

                //c---------------------------------------------------------------------
                //c     Compute up to NV results based on the current seed vector XV.
                //c---------------------------------------------------------------------
                int idx;
                for (i = 1; i <= n1; i++) {
                    idx = (i+j-1) + (2*nx*dims[2,1]*(k-1));
                    y[idx] = r46 * xv[i-1];
                }

                //c---------------------------------------------------------------------
                //c     If this is the last pass through the 140 loop, it is not necessary to
                //c     update the XV vector.
                //c---------------------------------------------------------------------
                if (j + n1 == n) goto goto150;

                //c---------------------------------------------------------------------
                //c     Update the XV vector by multiplying each element by AN (mod 2^46).
                //c---------------------------------------------------------------------
                for (i = 1; i <= nv; i++) {
                    t1 = r23 * xv[i-1];
                    x1 = (int) (t1);
                    x2 = xv[i-1] - t23 * x1;
                    t1 = a1 * x2 + a2 * x1;
                    t2 = (int) (r23 * t1);
                    yy = t1 - t23 * t2;
                    t3 = t23 * yy + a2 * x2;
                    t4 = (int) (r46 * t3);
                    xv[i-1] = t3 - t46 * t4;
                }

            }

            //     c---------------------------------------------------------------------
            //     c     Save the last seed in X so that subsequent calls to VRANLC will generate
            //     c     a continuous sequence.
            //     c---------------------------------------------------------------------
            goto150:  x = xv[n1-1];

        }

        public double min(int n1, int n2) {
            return n1<n2?n1:n2;
        }

        public void fft_init(int n) {

            //c---------------------------------------------------------------------
            //c compute the roots-of-unity array that will be used for subsequent FFTs. 
            //c---------------------------------------------------------------------

            //implicit none
            //include 'global.h'

            int m,nu,ku,i,j,ln;
            double t, ti;

            //c---------------------------------------------------------------------
            //c   Initialize the U array with sines and cosines in a manner that permits
            //c   stride one access at each FFT iteration.
            //c---------------------------------------------------------------------

            nu = n;
            m = ilog2Get(n);
            u[0,0] = m; // u(1)
            u[0,1] = 0.0;
            ku = 2;
            ln = 1;
            int idx;
            // Teste uUnidimensional 
                int count=1, idx2;
                uUnidimensional[0] = m;
                uUnidimensional[1] = 0.0;
            // fim
            for (j = 1; j <= m; j++) {
                t = pi / ln;

                for (i = 0; i <= ln - 1; i++) {
                    ti = i * t;
                    //u[i+ku] = dcmplx (cos (ti), sin(ti));
                    idx = (i+ku)-1;
                    u[idx,REAL] = Math.Cos(ti);  //u[i+ku]
                    u[idx,IMAG] = Math.Sin(ti);  //u[i+ku]
                    // Teste uUnidimensional:
                        idx2 = idx+count;
                        uUnidimensional[idx2]   = Math.Cos(ti);
                        uUnidimensional[idx2+1] = Math.Sin(ti);
                        count = count+1;
                    //fim Teste uUnidimensional
                }

                ku = ku + ln;
                ln = 2 * ln;
            }

        }

        public int ilog2Get(int n) {
            int nn, lg;
            if (n == 1) {
                return 0;
            }
            lg = 1;
            nn = 2;
            while (nn < n) {
                nn = nn * 2;
                lg = lg + 1;
            }
            return lg;
        }

    }
}
