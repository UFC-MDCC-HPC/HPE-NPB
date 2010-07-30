using System;
using System.IO;
using NPB.Lub;
using NPB3_0_JAV;
using NPB3_0_JAV.BMInOut;

namespace NPB {
    public class LU:LUBase {

        public LU(char c):base(c){}

        static void Main(String[] argv) {

            LU lu = null;
            LUBase.debug = false;

            try {
                string param = argv[0];
            }
            catch (Exception) {
                argv = new String[1];
                argv[0] = "CLASS=S"; // CLASS DEFAULT, IF USER NOT TYPE CLASS=S IN COMMAND-LINE ARGUMENT
            }
            char paramClass;
            if (!LUBase.debug) {
                BMArgs.ParseCmdLineArgs(argv, BMName);
                paramClass = BMArgs.CLASS;
            }
            else {
                paramClass = 'K';  //DEBUG: CHANGE TO [K=(S and 4 PROCESSORS)] OR [S=(S and 1 PROCESSOR)]
            }                      //DEBUG: OR [T=(A and 4 PROCESSORS)] OR [I=(B and 4 PROCESSORS)]

            try {
                lu = new LU(paramClass);
            }
            catch (OutOfMemoryException e) {
                Console.WriteLine(e.ToString());
                Environment.Exit(0);
            }
            lu.runBenchMark();
        }

        public void runBenchMark() {
            read_input();
            //c---------------------------------------------------------------------
            //c   set up processor grid
            //c---------------------------------------------------------------------
            //call proc_grid[];
            //c---------------------------------------------------------------------
            //c   determine the neighbors
            //c---------------------------------------------------------------------
            //call neighbors[];
            //c---------------------------------------------------------------------
            //c   set up sub-domain sizes
            //c---------------------------------------------------------------------
            //call subdomain[];
            //c---------------------------------------------------------------------
            //c   set up coefficients
            //c---------------------------------------------------------------------
            //call setcoeff[];
            //c---------------------------------------------------------------------
            //c   set the masks required for comm
            //c---------------------------------------------------------------------
            //call sethyper[];
            //c---------------------------------------------------------------------
            //c   set the boundary values for dependent variables
            //c---------------------------------------------------------------------
            //call setbv[];
            //c---------------------------------------------------------------------
            //c   set the initial values for dependent variables
            //c---------------------------------------------------------------------
            //call setiv[];
            //c---------------------------------------------------------------------
            //c   compute the forcing term based on prescribed exact solution
            //c---------------------------------------------------------------------
            //call erhs[];
            //c---------------------------------------------------------------------
            //c   perform one SSOR iteration to touch all data and program pages 
            //c---------------------------------------------------------------------
            //call ssor[1];
            //c---------------------------------------------------------------------
            //c   reset the boundary and initial values
            //c---------------------------------------------------------------------
            //call setbv[];
            //call setiv[];
            //c---------------------------------------------------------------------
            //c   perform the SSOR iterations
            //c---------------------------------------------------------------------
            //call ssor[itmax];
            //c---------------------------------------------------------------------
            //c   compute the solution error
            //c---------------------------------------------------------------------
            //call error[];
            //c---------------------------------------------------------------------
            //c   compute the surface integral
            //c---------------------------------------------------------------------
            //call pintgr[];
            //c---------------------------------------------------------------------
            //c   verification test
            //c---------------------------------------------------------------------
            //IF [id==0] {
                //call verify [ rsdnm, errnm, frc, class, verified ];
                //mflops = float[itmax]*[1984.77*float[ nx0 ] *float[ ny0 ] *float[ nz0 ] -10923.3*[float[ nx0+ny0+nz0 ]/3.]**2+27770.9* float[ nx0+ny0+nz0 ]/3.-144010.] / [maxtime*1000000.];
                //call print_results['LU', class, nx0, ny0, nz0, itmax, nnodes_compiled, num, maxtime, mflops, '          floating point', verified, npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6, '[none]'];
            //}
            //call mpi_finalize[ierr];
            //end;

            worldcomm.Barrier();
            mpi.Dispose();
        }

        public static double mod(double a, double b) { return (a % b); }

        public double min(int n1, int n2) { return n1<n2?n1:n2; }

        public double max(double n1, double n2) { return n1>n2?n1:n2; }

        public double pow2(double p) { return p * p; }

        public void read_input() {
            int fstatus=0, nnodes;
            //  c---------------------------------------------------------------------
            //  c    only root reads the input file
            //  c    if input file does not exist, it uses defaults
            //  c       ipr = 1 for detailed progress output
            //  c       inorm = how often the norm is printed [once every inorm iterations]
            //  c       itmax = number of pseufor(time steps
            //  c       dt = time step
            //  c       omega 1 over-relaxation factor for SSOR
            //  c       tolrsd = steady state residual tolerance levels
            //  c       nx, ny, nz = number of grid points in x, y, z directions
            //  c---------------------------------------------------------------------
            if (id == root) {
                string[] vetTemp = new string[13];
                try {
                    double kkkkk = 4.5;
                    Console.Write("Trying Reading from input file inputlu.data: ");
                    vetTemp = LUBase.readInputLuData("inputlu.data");//open [unit=3,file='inputlu.data',status='old', access='sequential',form='formatted', iostat=fstatus];
                }
                catch (System.IO.FileNotFoundException) {
                    Console.WriteLine("inputlu.data not found");
                    fstatus = 1;
                }
                Console.WriteLine(" NAS Parallel Benchmarks "+npbversion+" -- LU Benchmark ");
               if (fstatus == 0) {
                   Console.WriteLine("Reading from input file inputlu.data");
                   ipr       = int.Parse(vetTemp[0]);//read [3,*] ipr, inorm
                   inorm     = int.Parse(vetTemp[1]);
                   itmax     = int.Parse(vetTemp[2]);//read [3,*] itmax
                   dt        = double.Parse(vetTemp[3]);//read [3,*] dt
                   omega     = double.Parse(vetTemp[4]);//read [3,*] omega
                   tolrsd[1] = double.Parse(vetTemp[5]);
                   tolrsd[2] = double.Parse(vetTemp[6]);
                   tolrsd[3] = double.Parse(vetTemp[7]);
                   tolrsd[4] = double.Parse(vetTemp[8]);
                   tolrsd[5] = double.Parse(vetTemp[9]);//read [3,*] tolrsd[1],tolrsd[2],tolrsd[3],tolrsd[4],tolrsd[5]
                   nx0       = int.Parse(vetTemp[10]);
                   ny0       = int.Parse(vetTemp[11]);
                   nz0       = int.Parse(vetTemp[12]);//read [3,*] nx0, ny0, nz0
               } else {
                   ipr = ipr_default;
                   inorm = inorm_default;
                   itmax = itmax_default;
                   dt = dt_default;
                   omega = omega_default;
                   tolrsd[1] = tolrsd1_def;
                   tolrsd[2] = tolrsd2_def;
                   tolrsd[3] = tolrsd3_def;
                   tolrsd[4] = tolrsd4_def;
                   tolrsd[5] = tolrsd5_def;
                   nx0 = isiz01;
                   ny0 = isiz02;
                   nz0 = isiz03;
               }
               nnodes = num;//   call MPI_COMM_SIZE[MPI_COMM_WORLD, nnodes, ierror];
               //c---------------------------------------------------------------------
               //c   check problem size
               //c---------------------------------------------------------------------
               if (nnodes != nnodes_compiled) {
                   Console.WriteLine("Warning: program is running on"+nnodes+" processors, but was compiled for "+nnodes_compiled);
               }
               if ((nx0 < 4) || (ny0 < 4) || (nz0 < 4)) {
                   Console.WriteLine("PROBLEM SIZE IS TOO SMALL - "+" SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5");
                   worldcomm.Abort(0);
                   mpi.Dispose();//CALL MPI_ABORT[ MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR ];
                   Environment.Exit(0);
               }

               if ((nx0 > isiz01) || (ny0 > isiz02) || (nz0 > isiz03)){
                   Console.WriteLine("PROBLEM SIZE IS TOO LARGE - NX, NY AND NZ SHOULD BE LESS THAN OR EQUAL TO ISIZ01, ISIZ02 AND ISIZ03 RESPECTIVELY");
                   worldcomm.Abort(0);
                   mpi.Dispose();//      CALL MPI_ABORT[ MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR ];
                   Environment.Exit(0);
               }
               Console.WriteLine(" Size: " + nx0 + " x " + ny0 + " x " + nz0);
               Console.WriteLine(" Iterations: " + itmax);
               Console.WriteLine(" Number of processes: " + nnodes);
            }
            //call bcast_inputs; call bcast below:
            worldcomm.Broadcast<int>(ref ipr, root);          //call MPI_BCAST[ipr, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref inorm, root);        //call MPI_BCAST[inorm, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref itmax, root);        //call MPI_BCAST[itmax, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<double>(ref dt, root);        //call MPI_BCAST[dt, 1, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<double>(ref omega, root);     //call MPI_BCAST[omega, 1, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<double>(ref tolrsd[1], root); //call MPI_BCAST[tolrsd, 5, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<double>(ref tolrsd[2], root); //call MPI_BCAST[tolrsd, 5, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<double>(ref tolrsd[3], root); //call MPI_BCAST[tolrsd, 5, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<double>(ref tolrsd[4], root); //call MPI_BCAST[tolrsd, 5, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<double>(ref tolrsd[5], root); //call MPI_BCAST[tolrsd, 5, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref nx0, root);          //call MPI_BCAST[nx0, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref ny0, root);          //call MPI_BCAST[ny0, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref nz0, root);          //call MPI_BCAST[nz0, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
        }
    }
}


