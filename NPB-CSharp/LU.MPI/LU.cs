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

            //double pp = Math.Pow(2,((double)3)/2);
            //double ppp = Math.Pow(2,1.5);
            //int pare=0;

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

        public void imprimir1(int nd, double[] vetor, int i1, string nome) {
            if (node == nd) {
                int linha = 0;
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/" + nome + "-" + clss + "-Size-" + no_nodes + "-Prt-" + node + ".txt";
                System.IO.File.Delete(strPath);
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                int p1;
                for (p1 = i1; p1 < vetor.GetLength(0); p1++)
                    arquivo.WriteLine((linha++) + " [" + p1 + "]=" + vetor[p1]);
                arquivo.Close();
            }
        }

        public void imprimir2(int nd, double[,] vetor, int i1, int i2, string nome) {
            if (node == nd) {
                int linha = 0;
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/" + nome + "-" + clss + "-Size-" + no_nodes + "-Prt-" + node + ".txt";
                System.IO.File.Delete(strPath);
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                int p1, p2;
                for (p1 = i1; p1 < vetor.GetLength(0); p1++)
                    for (p2 = i2; p2 < vetor.GetLength(1); p2++)
                        arquivo.WriteLine((linha++) + " [" + p1 + "," + p2 + "]=" + vetor[p1, p2]);

                arquivo.Close();
            }
        }

        public void imprimir2Inverso(int nd, double[,] vetor, int i1, int i2, string nome) {
            if (node == nd) {
                int p1, p2, size0 = vetor.GetLength(0), size1 = vetor.GetLength(1);
                double[,] tmp2 = new double[size1, size0];
                for (int i = 0; i < size0; i++)
                    for (int j = 0; j < size1; j++)
                        tmp2[j, i] = vetor[i, j]; //tmp[i*size1+j];

                int linha = 0;
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/" + nome + "-" + clss + "-Size-" + no_nodes + "-Prt-" + node + ".txt";
                System.IO.File.Delete(strPath);
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                for (p1 = i2; p1 < tmp2.GetLength(0); p1++)
                    for (p2 = i1; p2 < tmp2.GetLength(1); p2++)
                        arquivo.WriteLine((linha++) + " [" + p1 + "," + p2 + "]=" + tmp2[p1, p2]);

                arquivo.Close();
            }
        }

        public void imprimir3(int nd, double[, ,] vetor, int i1, int i2, int i3, string nome) {
            if (node == nd) {
                int linha = 0;
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/" + nome + "-" + clss + "-Size-" + no_nodes + "-Prt-" + node + ".txt";
                System.IO.File.Delete(strPath);
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                int p1, p2, p3;
                for (p1 = i1; p1 < vetor.GetLength(0); p1++)
                    for (p2 = i2; p2 < vetor.GetLength(1); p2++)
                        for (p3 = i3; p3 < vetor.GetLength(2); p3++)
                            arquivo.WriteLine((linha++) + " [" + p1 + "," + p2 + "," + p3 + "]=" + vetor[p1, p2, p3]);

                arquivo.Close();
            }
        }

        public void imprimir4(int nd, double[, , ,] vetor, int i1, int i2, int i3, int i4, string nome) {
            if (id == nd) {
                int linha = 0;
                string strPath = "C:/Documents and Settings/Administrador/Desktop/Logs/" + nome + "-" + clss + "-Size-" + no_nodes + "-Prt-" + node + ".txt";
                System.IO.File.Delete(strPath);
                System.IO.TextWriter arquivo = System.IO.File.AppendText(strPath);
                int p1, p2, p3, p4;
                for (p1 = i1; p1 < vetor.GetLength(0); p1++)
                    for (p2 = i2; p2 < vetor.GetLength(1); p2++)
                        for (p3 = i3; p3 < vetor.GetLength(2); p3++)
                            for (p4 = i4; p4 < vetor.GetLength(3); p4++)
                                arquivo.WriteLine((linha++) + " [" + p1 + "," + p2 + "," + p3 + "," + p4 + "]=" + vetor[p1, p2, p3, p4]);

                arquivo.Close();
            }
        }

        public void imprimir(string s) {
            if(s=="u") imprimir4(root, u, 1, 0, 0, 1, "U");
            if(s=="rsd") imprimir4(root, rsd, 1, 0, 0, 1, "RSD");//rsd  = new double[isiz3+1, isiz2+4, isiz1+4, 5+1];//     rsd[5, -1:isiz1+2, -1:isiz2+2, isiz3];
            if(s=="frct") imprimir4(root, frct, 1, 0, 0, 1, "FRCT");
            if(s=="flux") imprimir4(root, flux, 1, 0, 0, 1, "FLUX");
            if(s=="a") imprimir4(root, a, 1, 1, 1, 1, "a");
            if(s=="b") imprimir4(root, b, 1, 1, 1, 1, "b");
            if(s=="c") imprimir4(root, c, 1, 1, 1, 1, "c");
            if(s=="d") imprimir4(root, d, 1, 1, 1, 1, "d");
            mpi.Dispose();
            Environment.Exit(0);
        }

        public void runBenchMark() {
            read_input();
            //c---------------------------------------------------------------------
            //c   set up processor grid
            //c---------------------------------------------------------------------
            proc_grid();
            //c---------------------------------------------------------------------
            //c   determine the neighbors
            //c---------------------------------------------------------------------
            neighbors();
            //c---------------------------------------------------------------------
            //c   set up sub-domain sizes
            //c---------------------------------------------------------------------
            subdomain();
            //c---------------------------------------------------------------------
            //c   set up coefficients
            //c---------------------------------------------------------------------
            setcoeff();
            //c---------------------------------------------------------------------
            //c   set the masks required for comm
            //c---------------------------------------------------------------------
            sethyper();
            //c---------------------------------------------------------------------
            //c   set the boundary values for dependent variables
            //c---------------------------------------------------------------------
            setbv();
            //c---------------------------------------------------------------------
            //c   set the initial values for dependent variables
            //c---------------------------------------------------------------------
            setiv();
            //c---------------------------------------------------------------------
            //c   compute the forcing term based on prescribed exact solution
            //c---------------------------------------------------------------------
            erhs();
            //c---------------------------------------------------------------------
            //c   perform one SSOR iteration to touch all data and program pages 
            //c---------------------------------------------------------------------
            ssor(1);
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

        public void proc_grid() {
            //c---------------------------------------------------------------------
            //c
            //c   set up a two-d grid for processors: column-major ordering of unknowns
            //c   NOTE: assumes a power-of-two number of processors
            //c
            //c---------------------------------------------------------------------
            xdim   = (int) Math.Pow(2,(ndim/2));//xdim   = 2**(ndim/2);
            if (mod(ndim,2)==1) xdim = xdim + xdim;
            ydim   = num/xdim;
            row    = (int) mod(id,xdim) + 1;
            col    = id/xdim + 1;
        }

        public void neighbors() {
            //  c---------------------------------------------------------------------
            //  c     figure out the neighbors and their wrap numbers for each processor
            //  c---------------------------------------------------------------------
            south = -1;
            east  = -1;
            north = -1;
            west  = -1;
            if (row>1) {
                    north = id -1;
            } else {
                    north = -1;
            }
            if (row < xdim) {
                south = id + 1;
            } else {
                south = -1;
            }

            if (col > 1) {
                west = id - xdim;
            } else {
                west = -1;
            }
            if (col < ydim) {
                east = id + xdim;
            } else {
                east = -1;
            }
        }

        public void subdomain() {
            int mm;
            //  c---------------------------------------------------------------------
            //  c
            //  c   set up the sub-domain sizes
            //  c
            //  c---------------------------------------------------------------------
            //  c---------------------------------------------------------------------
            //  c   x dimension
            //  c---------------------------------------------------------------------
            mm   = (int) mod(nx0,xdim);
            if (row<=mm) {
              nx = nx0/xdim + 1;
              ipt = (row-1)*nx;
            } else {
              nx = nx0/xdim;
              ipt = (row-1)*nx + mm;
            }
            //  c---------------------------------------------------------------------
            //  c   y dimension
            //  c---------------------------------------------------------------------
            mm   = (int) mod(ny0,ydim);
            if (col<=mm) {
              ny = ny0/ydim + 1;
              jpt = (col-1)*ny;
            } else {
              ny = ny0/ydim;
              jpt = (col-1)*ny + mm;
            }
            //  c---------------------------------------------------------------------
            //  c   z dimension
            //  c---------------------------------------------------------------------
            nz = nz0;
            //  c---------------------------------------------------------------------
            //  c   check the sub-domain size
            //  c---------------------------------------------------------------------
            if ( ( nx < 4 ) || ( ny < 4 ) || ( nz < 4 ) ) {
                Console.WriteLine("SUBDOMAIN SIZE IS TOO SMALL - ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS, "+
                    "SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL TO 4 THEY ARE CURRENTLY: "+nx+"x"+ny+"x"+nz);
                worldcomm.Abort(0);//CALL MPI_ABORT[ MPI_COMM_WORLD,ERRORCODE,IERROR ]
                mpi.Dispose();
                Environment.Exit(0);            
            }
            if ( ( nx > isiz1 ) || ( ny > isiz2 ) || ( nz > isiz3 ) ) {
                Console.WriteLine("SUBDOMAIN SIZE IS TOO LARGE - ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS" +
                    "SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY. THEY ARE CURRENTLY"+
                    " "+nx+"x"+ny+"x"+nz);
                worldcomm.Abort(0);//CALL MPI_ABORT[ MPI_COMM_WORLD,ERRORCODE, IERROR ]
                mpi.Dispose();
                Environment.Exit(0);
            }
            //  c---------------------------------------------------------------------
            //  c   set up the start and end in i and j extents for all processors
            //  c---------------------------------------------------------------------
            ist = 1;
            iend = nx;
            if (north==-1) ist = 2;
            if (south==-1) iend = nx - 1;
            jst = 1;
            jend = ny;
            if (west==-1) jst = 2;
            if (east==-1) jend = ny - 1;
        }

        public void setcoeff() {
            //---------------------------------------------------------------------
            //   set up coefficients
            //---------------------------------------------------------------------
            dxi   = 1.0d / ( nx0 - 1 );
            deta  = 1.0d / ( ny0 - 1 );
            dzeta = 1.0d / ( nz0 - 1 );

            tx1 = 1.0d/( dxi*dxi);
            tx2 = 1.0d/(2.0d*dxi);
            tx3 = 1.0d/dxi;

            ty1 = 1.0d/ ( deta * deta );
            ty2 = 1.0d/ ( 2.0d* deta );
            ty3 = 1.0d/ deta;

            tz1 = 1.0d/ ( dzeta * dzeta );
            tz2 = 1.0d/ ( 2.0d* dzeta );
            tz3 = 1.0d/ dzeta;

            ii1 = 2;
            ii2 = nx0 - 1;
            ji1 = 2;
            ji2 = ny0 - 2;
            ki1 = 3;
            ki2 = nz0 - 1;

            //  c---------------------------------------------------------------------
            //  c   diffusion coefficients
            //  c---------------------------------------------------------------------
            dx1 = 0.75d;
            dx2 = dx1;
            dx3 = dx1;
            dx4 = dx1;
            dx5 = dx1;

            dy1 = 0.75d;
            dy2 = dy1;
            dy3 = dy1;
            dy4 = dy1;
            dy5 = dy1;

            dz1 = 1.00d;
            dz2 = dz1;
            dz3 = dz1;
            dz4 = dz1;
            dz5 = dz1;
            //  c---------------------------------------------------------------------
            //  c   fourth difference dissipation
            //  c---------------------------------------------------------------------      
            dssp = (max(max(dx1,dy1),dz1))/4.0d; //dssp=(max(dx1, dy1, dz1))/4.0d
            //  c---------------------------------------------------------------------
            //  c   coefficients of the exact solution to the first pde
            //  c---------------------------------------------------------------------
            ce[1,1] = 2.0d;
            ce[1,2] = 0.0d;
            ce[1,3] = 0.0d;
            ce[1,4] = 4.0d;
            ce[1,5] = 5.0d;
            ce[1,6] = 3.0d;
            ce[1,7] = 5.0E-01;
            ce[1,8] = 2.0E-02;
            ce[1,9] = 1.0E-02;
            ce[1,10] = 3.0E-02;
            ce[1,11] = 5.0E-01;
            ce[1,12] = 4.0E-01;
            ce[1,13] = 3.0E-01;
            //  c---------------------------------------------------------------------
            //  c   coefficients of the exact solution to the second pde
            //  c---------------------------------------------------------------------
            ce[2,1] = 1.0d;
            ce[2,2] = 0.0d;
            ce[2,3] = 0.0d;
            ce[2,4] = 0.0d;
            ce[2,5] = 1.0d;
            ce[2,6] = 2.0d;
            ce[2,7] = 3.0d;
            ce[2,8] = 1.0E-02;
            ce[2,9] = 3.0E-02;
            ce[2,10] = 2.0E-02;
            ce[2,11] = 4.0E-01;
            ce[2,12] = 3.0E-01;
            ce[2,13] = 5.0E-01;
            //  c---------------------------------------------------------------------
            //  c   coefficients of the exact solution to the third pde
            //  c---------------------------------------------------------------------
            ce[3,1] = 2.0d;
            ce[3,2] = 2.0d;
            ce[3,3] = 0.0d;
            ce[3,4] = 0.0d;
            ce[3,5] = 0.0d;
            ce[3,6] = 2.0d;
            ce[3,7] = 3.0d;
            ce[3,8] = 4.0E-02;
            ce[3,9] = 3.0E-02;
            ce[3,10] = 5.0E-02;
            ce[3,11] = 3.0E-01;
            ce[3,12] = 5.0E-01;
            ce[3,13] = 4.0E-01;
            //  c---------------------------------------------------------------------
            //  c   coefficients of the exact solution to the fourth pde
            //  c---------------------------------------------------------------------
            ce[4,1] = 2.0d;
            ce[4,2] = 2.0d;
            ce[4,3] = 0.0d;
            ce[4,4] = 0.0d;
            ce[4,5] = 0.0d;
            ce[4,6] = 2.0d;
            ce[4,7] = 3.0d;
            ce[4,8] = 3.0E-02;
            ce[4,9] = 5.0E-02;
            ce[4,10] = 4.0E-02;
            ce[4,11] = 2.0E-01;
            ce[4,12] = 1.0E-01;
            ce[4,13] = 3.0E-01;
            //  c---------------------------------------------------------------------
            //  c   coefficients of the exact solution to the fifth pde
            //  c---------------------------------------------------------------------
            ce[5,1] = 5.0d;
            ce[5,2] = 4.0d;
            ce[5,3] = 3.0d;
            ce[5,4] = 2.0d;
            ce[5,5] = 1.0E-01;
            ce[5,6] = 4.0E-01;
            ce[5,7] = 3.0E-01;
            ce[5,8] = 5.0E-02;
            ce[5,9] = 4.0E-02;
            ce[5,10] = 3.0E-02;
            ce[5,11] = 1.0E-01;
            ce[5,12] = 3.0E-01;
            ce[5,13] = 2.0E-01;
        }

        public void sethyper() {
            //---------------------------------------------------------------------
            //    for each column in a hyperplane, istart = first row,
            //---------------------------------------------------------------------
            int i, j, iglob, jglob, kp;
            //---------------------------------------------------------------------
            // compute the pointers for hyperplanes
            //---------------------------------------------------------------------
              for(kp = 2; kp<=(nx0+ny0); kp++){
                icomms[kp] = false;
                icommn[kp] = false;
                icomme[kp] = false;
                icommw[kp] = false;
                //---------------------------------------------------------------------
                //  check to see if comm. to south is required
                //---------------------------------------------------------------------
                if (south!=-1) {
                  i     = iend;
                  iglob = ipt + i;
                  jglob = kp - iglob;
                  j     = jglob - jpt;
                  if (jglob>=2 && jglob<=ny0-1 && j>=jst && j<=jend) icomms[kp] = true;
                }
                //---------------------------------------------------------------------
                //  check to see if comm. to north is required
                //---------------------------------------------------------------------
                if (north!=-1) {
                  i     = ist;
                  iglob = ipt + i;
                  jglob = kp - iglob;
                  j     = jglob - jpt;
                  if (jglob>=2 && jglob<=ny0-1 && j>=jst && j<=jend) icommn[kp] = true;
                }
                //---------------------------------------------------------------------
                //  check to see if comm. to east is required
                //---------------------------------------------------------------------
                if (east!=-1) {
                  j     = jend;
                  jglob = jpt + j;
                  iglob = kp - jglob;
                  i     = iglob - ipt;
                  if (iglob>=2 && iglob<=nx0-1 && i>=ist && i<=iend) icomme[kp] = true;
                }
                //---------------------------------------------------------------------
                //  check to see if comm. to west is required
                //---------------------------------------------------------------------
                if (west!=-1) {
                  j = jst;
                  jglob = jpt + j;
                  iglob = kp - jglob;
                  i     = iglob - ipt;
                  if (iglob>=2 && iglob<=nx0-1 && i>=ist && i<=iend) icommw[kp] = true;
                }
              }
              icomms[1] = false;
              icommn[1] = false;
              icomme[1] = false;
              icommw[1] = false;
              icomms[nx0+ny0+1] = false;
              icommn[nx0+ny0+1] = false;
              icomme[nx0+ny0+1] = false;
              icommw[nx0+ny0+1] = false;
        }

        public void setbv() {
            //  c---------------------------------------------------------------------
            //  c   set the boundary values of dependent variables
            //  c---------------------------------------------------------------------
            int i, j, k, iglob, jglob;
            //  c---------------------------------------------------------------------
            //  c   set the dependent variable values along the top and bottom faces
            //  c---------------------------------------------------------------------
            for(j = 1; j<= ny; j++){
               jglob = jpt + j;
               for(i = 1; i<= nx; i++){
                 iglob = ipt + i;
                  exact4(iglob,jglob,1, u,      1, j+1, i+1);   //exact( iglob, jglob, 1, u[ 1, i, j, 1 ] );
                  exact4(iglob,jglob,nz,u,     nz, j+1, i+1);   //exact( iglob, jglob, nz, u[ 1, i, j, nz ] );
               }
            }
            //---------------------------------------------------------------------
            //   set the dependent variable values along north and south faces
            //---------------------------------------------------------------------
            if (west==-1) {
               for(k = 1; k<= nz; k++){
                  for(i = 1; i<= nx; i++){
                     iglob = ipt + i;
                     exact4(iglob,1,k,u,   k,1+1,i+1);   //call exact[ iglob, 1, k, u[ 1, i, 1, k ] ];
                  }
               }
            }
            if (east==-1) {
                for(k = 1; k<= nz; k++){
                   for(i = 1; i<= nx; i++){
                      iglob = ipt + i;
                      exact4(iglob,ny0,k,u,     k,ny+1,i+1); //call exact[ iglob, ny0, k, u[ 1, i, ny, k ] ];
                   }
                }
            }
            //---------------------------------------------------------------------
            //   set the dependent variable values along east and west faces
            //---------------------------------------------------------------------
            if (north==-1) {
               for(k = 1; k<= nz; k++){
                  for(j = 1; j<= ny; j++){
                     jglob = jpt + j;
                     exact4(1,jglob,k,u,   k,j+1,1+1);  //call exact[ 1, jglob, k, u[ 1, 1, j, k ] ];
                  }
               }
            }
            if (south==-1) {
               for(k = 1; k<= nz; k++){
                  for(j = 1; j<= ny; j++){
                    jglob = jpt + j;
                    exact4(nx0,jglob,k,u,    k,j+1,nx+1); // call exact[ nx0, jglob, k, u[ 1, nx, j, k ] ];
                  }
               }
            }
        }

        public void exact4(int i, int j, int k, double[,,,] u000ijk, int i1, int i2, int i3) {
            //---------------------------------------------------------------------
            //
            //   compute the exact solution at [i,j,k]
            //
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //int i, j, k;
            //double u000ijk[*];
            int m;
            double xi, eta, zeta;

            xi   = ((double) (i-1))/(nx0-1); //( dble ( i - 1 ) ) / ( nx0 - 1 );
            eta  = ((double) (j-1))/(ny0-1); //( dble ( j - 1 ) ) / ( ny0 - 1 );
            zeta = ((double) (k-1))/(nz -1);  //( dble ( k - 1 ) ) / ( nz - 1 );
            for(m = 1; m<= 5; m++){
               u000ijk[i1,i2,i3,m] = ce[m,1] + ce[m,2]*xi + ce[m,3]*eta + ce[m,4]*zeta + ce[m,5]*xi*xi + ce[m,6]*eta*eta
                   + ce[m,7]*zeta*zeta + ce[m,8]*xi*xi*xi + ce[m,9]*eta*eta*eta + ce[m,10]*zeta*zeta*zeta + ce[m,11]*xi*xi*xi*xi
                   + ce[m,12]*eta*eta*eta*eta + ce[m,13]*zeta*zeta*zeta*zeta;
            }
        }

        public void exact1(int i, int j, int k, double[] u000ijk) {
            //---------------------------------------------------------------------
            //
            //   compute the exact solution at [i,j,k]
            //
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //int i, j, k;
            //double u000ijk[*];
            int m;
            double xi, eta, zeta;
            xi   = ((double) (i-1))/(nx0-1); //( dble ( i - 1 ) ) / ( nx0 - 1 );
            eta  = ((double) (j-1))/(ny0-1); //( dble ( j - 1 ) ) / ( ny0 - 1 );
            zeta = ((double) (k-1))/(nz -1);  //( dble ( k - 1 ) ) / ( nz - 1 );
            for(m = 1; m<= 5; m++){
               u000ijk[m] = ce[m,1] + ce[m,2]*xi + ce[m,3]*eta + ce[m,4]*zeta + ce[m,5]*xi*xi + ce[m,6]*eta*eta
                   + ce[m,7]*zeta*zeta + ce[m,8]*xi*xi*xi + ce[m,9]*eta*eta*eta + ce[m,10]*zeta*zeta*zeta + ce[m,11]*xi*xi*xi*xi
                   + ce[m,12]*eta*eta*eta*eta + ce[m,13]*zeta*zeta*zeta*zeta;
            }
        }

        public void setiv() {
            //---------------------------------------------------------------------
            //   set the initial values of independent variables based on tri-linear
            //   interpolation of boundary values in the computational space.
            int i, j, k, m;
            int iglob, jglob;
            double  xi, eta, zeta;
            double  pxi, peta, pzeta;
            double[] ue_1jk   = new double[5+1];   //ue_1jk[5]
            double[] ue_nx0jk = new double[5+1];   //ue_nx0jk[5]
            double[] ue_i1k   = new double[5+1];   //ue_i1k[5]
            double[] ue_iny0k = new double[5+1];   //ue_iny0k[5]
            double[] ue_ij1   = new double[5+1];   //ue_ij1[5]
            double[] ue_ijnz  = new double[5+1];   //ue_ijnz[5]
            for(k = 2; k<=nz-1; k++){
                zeta = ((double)(k-1))/(nz-1);
                for(j = 1; j<= ny; j++){
                    jglob = jpt + j;
                    if (jglob!=1 && jglob!=ny0) {
                        eta = ((double)(jglob-1))/(ny0-1);
                        for(i = 1; i<= nx; i++){
                            iglob = ipt + i;
                            if (iglob!=1 && iglob!=nx0) {
                                xi = ((double)(iglob-1))/(nx0-1);
                                exact1(1,jglob,k,ue_1jk);
                                exact1(nx0,jglob,k,ue_nx0jk);
                                exact1(iglob,1,k,ue_i1k);
                                exact1(iglob,ny0,k,ue_iny0k);
                                exact1(iglob,jglob,1,ue_ij1);
                                exact1(iglob,jglob,nz,ue_ijnz);
                                for(m = 1; m<= 5; m++){
                                    pxi =   (1.0d-xi) * ue_1jk[m] + xi   * ue_nx0jk[m];
                                    peta =  (1.0d-eta) * ue_i1k[m] + eta   * ue_iny0k[m];
                                    pzeta = (1.0d-zeta) * ue_ij1[m] + zeta   * ue_ijnz[m];
                                    u[k,j+1,i+1,m] = pxi + peta + pzeta - pxi * peta - peta * pzeta - pzeta * pxi + pxi * peta * pzeta;
                                }
                            }
                        }
                    }
                }
            }
        }

        //erhs.f
        public void erhs() {
            //---------------------------------------------------------------------
            //   compute the right hand side based on exact solution
            //---------------------------------------------------------------------
            int i, j, k, m;
            int iglob, jglob;
            int iex;
            int L1, L2;
            int ist1, iend1;
            int jst1, jend1;
            double dsspm;
            double xi, eta, zeta;
            double q;
            double u21, u31, u41;
            double tmp;
            double u21i, u31i, u41i, u51i;
            double u21j, u31j, u41j, u51j;
            double u21k, u31k, u41k, u51k;
            double u21im1, u31im1, u41im1, u51im1;
            double u21jm1, u31jm1, u41jm1, u51jm1;
            double u21km1, u31km1, u41km1, u51km1;
            dsspm = dssp;
            for(k = 1; k<= nz; k++){
               for(j = 1; j<= ny; j++){
                  for(i = 1; i<= nx; i++){
                     for(m = 1; m<= 5; m++){
                        frct[k,j+1,i+1,m] = 0.0d; //frct[ m, i, j, k ] = 0.0d;
                     }
                  }
               }
            }
            for(k = 1; k<= nz; k++){
               zeta = ((double)(k-1))/(nz-1);
               for(j = 1; j<= ny; j++){
                  jglob = jpt + j;
                  eta = ((double)(jglob-1))/(ny0-1);
                  for(i = 1; i<= nx; i++){
                     iglob = ipt + i;
                     xi = ((double)(iglob-1))/(nx0-1);
                     for(m = 1; m<= 5; m++){  //rsd[m,i,j,k] =  ce[m,1]
                        rsd[k,j+1,i+1,m] =  ce[m,1]
                            + ce[m,2] * xi
                            + ce[m,3] * eta
                            + ce[m,4] * zeta
                            + ce[m,5] * xi * xi
                            + ce[m,6] * eta * eta
                            + ce[m,7] * zeta * zeta
                            + ce[m,8] * xi * xi * xi
                            + ce[m,9] * eta * eta * eta
                            + ce[m,10] * zeta * zeta * zeta
                            + ce[m,11] * xi * xi * xi * xi
                            + ce[m,12] * eta * eta * eta * eta
                            + ce[m,13] * zeta * zeta * zeta * zeta;
                     }
                  }
               }
            }
            //---------------------------------------------------------------------
            //   xi-direction flux differences
            //---------------------------------------------------------------------
            //c   iex = flag : iex = 0  north/south communication
            //c              : iex = 1  east/west communication
            //c---------------------------------------------------------------------
            iex   = 0;
            //---------------------------------------------------------------------
            //   communicate and receive/send two rows of data
            //---------------------------------------------------------------------
            exchange_3(rsd,iex);
            L1 = 0;
            if (north==-1) L1 = 1;
            L2 = nx + 1;
            if (south==-1) L2 = nx;

            ist1 = 1;
            iend1 = nx;
            if (north==-1) ist1 = 4;
            if (south==-1) iend1 = nx - 3;
            for(k = 2; k<= nz - 1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = L1; i<= L2; i++){
                     flux[k,j,i,1] = rsd[k,j+1,i+1,2]; //flux[1,i,j,k] = rsd[2,i,j,k];
                     u21           = rsd[k,j+1,i+1,2]/rsd[k,j+1,i+1,1]; //u21 = rsd[2,i,j,k] / rsd[1,i,j,k];
                     //c -- q = 0.50d*(rsd[2,i,j,k]*rsd[2,i,j,k] + rsd[3,i,j,k]*rsd[3,i,j,k] + rsd[4,i,j,k]*rsd[4,i,j,k])/rsd[1,i,j,k];
                     q=0.50d*(rsd[k,j+1,i+1,2]*rsd[k,j+1,i+1,2]+rsd[k,j+1,i+1,3]*rsd[k,j+1,i+1,3]+rsd[k,j+1,i+1,4]*rsd[k,j+1,i+1,4])/rsd[k,j+1,i+1,1];
                     flux[k,j,i,2] =     rsd[k,j+1,i+1,2]*u21 + c2*(rsd[k,j+1,i+1,5] - q);//flux[2,i,j,k]=rsd[2,i,j,k]*u21+c2*(rsd[5,i,j,k]-q);
                     flux[k,j,i,3] =     rsd[k,j+1,i+1,3] * u21;                          //flux[3,i,j,k]=rsd[3,i,j,k] * u21;
                     flux[k,j,i,4] =     rsd[k,j+1,i+1,4] * u21;                          //flux[4,i,j,k]=rsd[4,i,j,k] * u21;
                     flux[k,j,i,5] = (c1*rsd[k,j+1,i+1,5] - c2*q)*u21;                    //flux[5,i,j,k]=(c1*rsd[5,i,j,k] - c2*q)*u21;
                  }
               }
            }
            for(k = 2; k<= nz - 1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){ //frct[m,i,j,k] =  frct[m,i,j,k] - tx2 * (flux[m,i+1,j,k] - flux[m,i-1,j,k]);
                        frct[k,j+1,i+1,m] =  frct[k,j+1,i+1,m] - tx2 * (flux[k,j,i+1,m] - flux[k,j,i-1,m]);
                     }
                  }
                  for(i = ist; i<= L2; i++){
                     tmp   = 1.0d/rsd[k,j+1,i+1,1];
                     u21i = tmp * rsd[k,j+1,i+1,2];
                     u31i = tmp * rsd[k,j+1,i+1,3];
                     u41i = tmp * rsd[k,j+1,i+1,4];
                     u51i = tmp * rsd[k,j+1,i+1,5];
                     tmp   = 1.0d/rsd[k,j+1,i,1];
                     u21im1 = tmp*rsd[k,j+1,i,2];
                     u31im1 = tmp*rsd[k,j+1,i,3];
                     u41im1 = tmp*rsd[k,j+1,i,4];
                     u51im1 = tmp*rsd[k,j+1,i,5];

                     flux[k,j,i,2] = (4.0d/3.0d)*tx3*(u21i - u21im1);
                     flux[k,j,i,3] = tx3 * (u31i - u31im1);
                     flux[k,j,i,4] = tx3 * (u41i - u41im1);
                     flux[k,j,i,5] = 0.50d*(1.0d-c1*c5)*tx3*
                         ((pow2(u21i)+pow2(u31i)+pow2(u41i))
                         -(pow2(u21im1)+pow2(u31im1)+pow2(u41im1)))
                         + (1.0d/6.0d)*tx3*(pow2(u21i) - pow2(u21im1))+c1*c5*tx3*(u51i-u51im1);
                  }
                  for(i = ist; i<= iend; i++){
                     frct[k,j+1,i+1,1] = frct[k,j+1,i+1,1]+dx1*tx1*(rsd[k,j+1,i,1]-2.0d*rsd[k,j+1,i+1,1]+rsd[k,j+1,i+2,1]);
                     frct[k,j+1,i+1,2] = frct[k,j+1,i+1,2]+tx3*c3*c4*(flux[k,j,i+1,2]-flux[k,j,i,2])+dx2*tx1*(rsd[k,j+1,i,2]-2.0d*rsd[k,j+1,i+1,2]+rsd[k,j+1,i+2,2]);
                     frct[k,j+1,i+1,3] = frct[k,j+1,i+1,3]+tx3*c3*c4*(flux[k,j,i+1,3]-flux[k,j,i,3])+dx3*tx1*(rsd[k,j+1,i,3]-2.0d*rsd[k,j+1,i+1,3]+rsd[k,j+1,i+2,3]);
                     frct[k,j+1,i+1,4] = frct[k,j+1,i+1,4]+tx3*c3*c4*(flux[k,j,i+1,4]-flux[k,j,i,4])+dx4*tx1*(rsd[k,j+1,i,4]-2.0d*rsd[k,j+1,i+1,4]+rsd[k,j+1,i+2,4]);
                     frct[k,j+1,i+1,5] = frct[k,j+1,i+1,5]+tx3*c3*c4*(flux[k,j,i+1,5]-flux[k,j,i,5])+dx5*tx1*(rsd[k,j+1,i,5]-2.0d*rsd[k,j+1,i+1,5]+rsd[k,j+1,i+2,5]);
                  }
                  //c---------------------------------------------------------------------
                  //c   Fourth-order dissipation
                  //c---------------------------------------------------------------------
                  if (north==-1) {
                   for(m = 1; m<= 5; m++){
                     frct[k,j+1,3,m] = frct[k,j+1,3,m]-dsspm*(+5.0d*rsd[k,j+1,3,m]-4.0d*rsd[k,j+1,4,m]+rsd[k,j+1,5,m]);
                     frct[k,j+1,4,m] = frct[k,j+1,4,m]-dsspm*(-4.0d*rsd[k,j+1,3,m]+6.0d*rsd[k,j+1,4,m]-4.0d*rsd[k,j+1,5,m]+rsd[k,j+1,6,m]);
                   }
                  }
                  for(i = ist1; i<=iend1; i++){
                     for(m = 1; m<= 5; m++){
                        frct[k,j+1,i+1,m] = frct[k,j+1,i+1,m]-dsspm*(rsd[k,j+1,i-1,m]-
                            4.0d*rsd[k,j+1,i,m]+6.0d*rsd[k,j+1,i+1,m]-4.0d*rsd[k,j+1,i+2,m]+rsd[k,j+1,i+3,m]);
                     }
                  }
                  if (south==-1) {
                   for(m = 1; m<= 5; m++){
                     frct[k,j+1,nx-1,m] = frct[k,j+1,nx-1,m]-dsspm*(rsd[k,j+1,nx-3,m]-4.0d*rsd[k,j+1,nx-2,m]+6.0d*rsd[k,j+1,nx-1,m]-4.0d*rsd[k,j+1,nx,m]);
                     frct[k,j+1,nx,m]   = frct[k,j+1,nx,m]  -dsspm*(rsd[k,j+1,nx-2,m]-4.0d*rsd[k,j+1,nx-1,m]+5.0d*rsd[k,j+1,nx,m]);
                   }
                  }
               }
            }
            //---------------------------------------------------------------------
            //   eta-direction flux differences
            //---------------------------------------------------------------------
            //c   iex = flag : iex = 0  north/south communication
            //c              : iex = 1  east/west communication
            //c
            //c---------------------------------------------------------------------
            iex   = 1;
            //---------------------------------------------------------------------
            //   communicate and receive/send two rows of data
            //---------------------------------------------------------------------
            exchange_3(rsd,iex);
            L1 = 0;
            if (west==-1) L1 = 1;
            L2 = ny + 1;
            if (east==-1) L2 = ny;
            jst1 = 1;
            jend1 = ny;
            if (west==-1) jst1 = 4;
            if (east==-1) jend1 = ny - 3;
            for(k = 2; k<= nz - 1; k++){
               for(j = L1; j<= L2; j++){
                  for(i = ist; i<= iend; i++){
                     flux[k,j,i,1] = rsd[k,j+1,i+1,3];
                     u31 = rsd[k,j+1,i+1,3] / rsd[k,j+1,i+1,1];
                     q          = 0.50d*(rsd[k,j+1,i+1,2]*rsd[k,j+1,i+1,2]+rsd[k,j+1,i+1,3]*rsd[k,j+1,i+1,3]
                                        +rsd[k,j+1,i+1,4]*rsd[k,j+1,i+1,4])/rsd[k,j+1,i+1,1];
                     flux[k,j,i,2] =     rsd[k,j+1,i+1,2]*u31;
                     flux[k,j,i,3] =     rsd[k,j+1,i+1,3]*u31+c2*(rsd[k,j+1,i+1,5]-q);
                     flux[k,j,i,4] =     rsd[k,j+1,i+1,4]*u31;
                     flux[k,j,i,5] = (c1*rsd[k,j+1,i+1,5]-c2*q)*u31;
                  }
               }
            }
            for(k = 2; k<= nz - 1; k++){
               for(i = ist; i<= iend; i++){
                  for(j = jst; j<= jend; j++){
                     for(m = 1; m<= 5; m++){
                        frct[k,j+1,i+1,m] =  frct[k,j+1,i+1,m] - ty2 * ( flux[k,j+1,i,m] - flux[k,j-1,i,m] );
                     }
                  }
               }
                for(j = jst; j<= L2; j++){
                    for(i = ist; i<= iend; i++){
                        tmp = 1.0d / rsd[k,j+1,i+1,1];
                        u21j = tmp * rsd[k,j+1,i+1,2];
                        u31j = tmp * rsd[k,j+1,i+1,3];
                        u41j = tmp * rsd[k,j+1,i+1,4];
                        u51j = tmp * rsd[k,j+1,i+1,5];
                        tmp = 1.0d / rsd[k,j,i+1,1];
                        u21jm1 = tmp*rsd[k,j,i+1,2];
                        u31jm1 = tmp*rsd[k,j,i+1,3];
                        u41jm1 = tmp*rsd[k,j,i+1,4];
                        u51jm1 = tmp*rsd[k,j,i+1,5];
                        flux[k,j,i,2] = ty3*(u21j-u21jm1);
                        flux[k,j,i,3] = (4.0d/3.0d)*ty3*(u31j-u31jm1);
                        flux[k,j,i,4] = ty3*(u41j-u41jm1);
                        flux[k,j,i,5] = 0.50d*(1.0d-c1*c5)*ty3*((pow2(u21j)+pow2(u31j)+pow2(u41j))-(pow2(u21jm1)+pow2(u31jm1)+pow2(u41jm1)))+(1.0d/6.0d)*ty3*(pow2(u31j)-pow2(u31jm1))+c1*c5*ty3*(u51j-u51jm1);
                    }
                }
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     frct[k,j+1,i+1,1] = frct[k,j+1,i+1,1]+dy1*ty1*(rsd[k,j,i+1,1]-2.0d*rsd[k,j+1,i+1,1]+rsd[k,j+2,i+1,1]);
                     frct[k,j+1,i+1,2] = frct[k,j+1,i+1,2]+ty3*c3*c4*(flux[k,j+1,i,2]-flux[k,j,i,2])+dy2*ty1*(rsd[k,j,i+1,2]-2.0d*rsd[k,j+1,i+1,2]+rsd[k,j+2,i+1,2]);
                     frct[k,j+1,i+1,3] = frct[k,j+1,i+1,3]+ty3*c3*c4*(flux[k,j+1,i,3]-flux[k,j,i,3])+dy3*ty1*(rsd[k,j,i+1,3]-2.0d*rsd[k,j+1,i+1,3]+rsd[k,j+2,i+1,3]);
                     frct[k,j+1,i+1,4] = frct[k,j+1,i+1,4]+ty3*c3*c4*(flux[k,j+1,i,4]-flux[k,j,i,4])+dy4*ty1*(rsd[k,j,i+1,4]-2.0d*rsd[k,j+1,i+1,4]+rsd[k,j+2,i+1,4]);
                     frct[k,j+1,i+1,5] = frct[k,j+1,i+1,5]+ty3*c3*c4*(flux[k,j+1,i,5]-flux[k,j,i,5])+dy5*ty1*(rsd[k,j,i+1,5]-2.0d*rsd[k,j+1,i+1,5]+rsd[k,j+2,i+1,5]);
                  }
               }
            //      c---------------------------------------------------------------------
            //      c   fourth-order dissipation
            //      c---------------------------------------------------------------------
               if (west==-1) {
                  for(i = ist; i<= iend; i++){
                   for(m = 1; m<= 5; m++){
                     frct[k,3,i+1,m] = frct[k,3,i+1,m]-dsspm*(+5.0d*rsd[k,3,i+1,m]-4.0d*rsd[k,4,i+1,m]+rsd[k,5,i+1,m]);
                     frct[k,4,i+1,m] = frct[k,4,i+1,m]-dsspm*(-4.0d*rsd[k,3,i+1,m]+6.0d*rsd[k,4,i+1,m]-4.0d*rsd[k,5,i+1,m]+rsd[k,6,i+1,m]);
                   }
                  }
               }
               for(j = jst1; j<= jend1; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        frct[k,j+1,i+1,m]=frct[k,j+1,i+1,m]-dsspm*(rsd[k,j-1,i+1,m]-4.0d*rsd[k,j,i+1,m]+6.0d*rsd[k,j+1,i+1,m]-4.0d*rsd[k,j+2,i+1,m]+rsd[k,j+3,i+1,m]);
                     }
                  }
               }
               if (east==-1) {
                  for(i = ist; i<= iend; i++){
                   for(m = 1; m<= 5; m++){
                     frct[k,ny-1,i+1,m] = frct[k,ny-1,i+1,m]-dsspm*(rsd[k,ny-3,i+1,m]-4.0d*rsd[k,ny-2,i+1,m]+6.0d*rsd[k,ny-1,i+1,m]-4.0d*rsd[k,ny,i+1,m]);
                     frct[k,ny  ,i+1,m] = frct[k,ny  ,i+1,m]-dsspm*(rsd[k,ny-2,i+1,m]-4.0d*rsd[k,ny-1,i+1,m]+5.0d*rsd[k,ny  ,i+1,m]);
                   }
                  }
               }
            }
            //  c---------------------------------------------------------------------
            //  c   zeta-direction flux differences
            //  c---------------------------------------------------------------------
            for(k = 1; k<= nz; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     flux[k,j,i,1] = rsd[k,j+1,i+1,4];      //flux[1,i,j,k] = rsd[4,i,j,k];
                     u41 = rsd[k,j+1,i+1,4] / rsd[k,j+1,i+1,1]; //u41 = rsd[4,i,j,k] / rsd[1,i,j,k];
                     q = 0.50d*(rsd[k,j+1,i+1,2]*rsd[k,j+1,i+1,2]+rsd[k,j+1,i+1,3]*rsd[k,j+1,i+1,3]+rsd[k,j+1,i+1,4]*rsd[k,j+1,i+1,4])/rsd[k,j+1,i+1,1];
                     flux[k,j,i,2] =rsd[k,j+1,i+1,2] * u41;
                     flux[k,j,i,3] =rsd[k,j+1,i+1,3] * u41;
                     flux[k,j,i,4] =rsd[k,j+1,i+1,4] * u41 + c2*(rsd[k,j+1,i+1,5] - q);
                     flux[k,j,i,5] =(c1*rsd[k,j+1,i+1,5]-c2*q)*u41;
                  }
               }
            }
            for(k = 2; k<= nz - 1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        frct[k,j+1,i+1,m] =  frct[k,j+1,i+1,m] - tz2 * (flux[k+1,j,i,m] - flux[k-1,j,i,m]);
                     }
                  }
               }
            }
            for(k = 2; k<= nz; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     tmp = 1.0d / rsd[k,j+1,i+1,1];
                     u21k = tmp * rsd[k,j+1,i+1,2];
                     u31k = tmp * rsd[k,j+1,i+1,3];
                     u41k = tmp * rsd[k,j+1,i+1,4];
                     u51k = tmp * rsd[k,j+1,i+1,5];

                     tmp = 1.0d / rsd[k-1,j+1,i+1,1];
                     u21km1 = tmp*rsd[k-1,j+1,i+1,2];
                     u31km1 = tmp*rsd[k-1,j+1,i+1,3];
                     u41km1 = tmp*rsd[k-1,j+1,i+1,4];
                     u51km1 = tmp*rsd[k-1,j+1,i+1,5];

                     flux[k,j,i,2] = tz3 * (u21k - u21km1);
                     flux[k,j,i,3] = tz3 * (u31k - u31km1);
                     flux[k,j,i,4] = (4.0d/3.0d) * tz3 * (u41k - u41km1);
                     flux[k,j,i,5] = 0.50d*(1.0d-c1*c5)*tz3*((pow2(u21k)+pow2(u31k)+pow2(u41k))-(pow2(u21km1)+pow2(u31km1)+pow2(u41km1)))+(1.0d/6.0d)*tz3*(pow2(u41k)-pow2(u41km1))+c1*c5*tz3*(u51k-u51km1);
                  }
               }
            }
            for(k = 2; k<= nz - 1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     frct[k,j+1,i+1,1] = frct[k,j+1,i+1,1]+dz1*tz1*(rsd[k+1,j+1,i+1,1]-2.0d*rsd[k,j+1,i+1,1]+rsd[k-1,j+1,i+1,1]);
                     frct[k,j+1,i+1,2] = frct[k,j+1,i+1,2]+tz3*c3*c4*(flux[k+1,j,i,2]-flux[k,j,i,2])+dz2*tz1*(rsd[k+1,j+1,i+1,2]-2.0d*rsd[k,j+1,i+1,2]+rsd[k-1,j+1,i+1,2]);
                     frct[k,j+1,i+1,3] = frct[k,j+1,i+1,3]+tz3*c3*c4*(flux[k+1,j,i,3]-flux[k,j,i,3])+dz3*tz1*(rsd[k+1,j+1,i+1,3]-2.0d*rsd[k,j+1,i+1,3]+rsd[k-1,j+1,i+1,3]);
                     frct[k,j+1,i+1,4] = frct[k,j+1,i+1,4]+tz3*c3*c4*(flux[k+1,j,i,4]-flux[k,j,i,4])+dz4*tz1*(rsd[k+1,j+1,i+1,4]-2.0d*rsd[k,j+1,i+1,4]+rsd[k-1,j+1,i+1,4]);
                     frct[k,j+1,i+1,5] = frct[k,j+1,i+1,5]+tz3*c3*c4*(flux[k+1,j,i,5]-flux[k,j,i,5])+dz5*tz1*(rsd[k+1,j+1,i+1,5]-2.0d*rsd[k,j+1,i+1,5]+rsd[k-1,j+1,i+1,5]);
                  }
               }
            }
            //  c---------------------------------------------------------------------
            //  c   fourth-order dissipation
            //  c---------------------------------------------------------------------
            for(j = jst; j<= jend; j++){
               for(i = ist; i<= iend; i++){
                  for(m = 1; m<= 5; m++){
                     frct[2,j+1,i+1,m] = frct[2,j+1,i+1,m]-dsspm*(+5.0d*rsd[2,j+1,i+1,m]-4.0d*rsd[3,j+1,i+1,m]+rsd[4,j+1,i+1,m]);
                     frct[3,j+1,i+1,m] = frct[3,j+1,i+1,m]-dsspm*(-4.0d*rsd[2,j+1,i+1,m]+6.0d*rsd[3,j+1,i+1,m]-4.0d*rsd[4,j+1,i+1,m]+rsd[5,j+1,i+1,m]);
                  }
               }
            }
            for(k = 4; k<= nz - 3; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        frct[k,j+1,i+1,m]=frct[k,j+1,i+1,m]-dsspm*(rsd[k-2,j+1,i+1,m]-4.0d*rsd[k-1,j+1,i+1,m]+6.0d*rsd[k,j+1,i+1,m]-4.0d*rsd[k+1,j+1,i+1,m]+rsd[k+2,j+1,i+1,m]);
                     }
                  }
               }
            }
            for(j = jst; j<= jend; j++){
               for(i = ist; i<= iend; i++){
                  for(m = 1; m<= 5; m++){
                     frct[nz-2,j+1,i+1,m]=frct[nz-2,j+1,i+1,m]-dsspm*(rsd[nz-4,j+1,i+1,m]- 4.0d*rsd[nz-3,j+1,i+1,m]+6.0d*rsd[nz-2,j+1,i+1,m]-4.0d*rsd[nz-1,j+1,i+1,m]);
                     frct[nz-1,j+1,i+1,m]=frct[nz-1,j+1,i+1,m]-dsspm*(rsd[nz-3,j+1,i+1,m]- 4.0d*rsd[nz-2,j+1,i+1,m]+5.0d*rsd[nz-1,j+1,i+1,m]);
                  }
               }
            }
        }
           //Exchange_3.f
        public void exchange_3(double[, , ,] g, int iex) {
            //c---------------------------------------------------------------------
            //c   compute the right hand side based on exact solution
            //c---------------------------------------------------------------------
            //c---------------------------------------------------------------------
            //c  input parameters
            //c---------------------------------------------------------------------
            //double  g[5,-1:isiz1+2,-1:isiz2+2,isiz3];
            //int iex;
            int i, j, k;
            int ipos1, ipos2;
            //int mid;
            //int STATUS[MPI_STATUS_SIZE];
            //int IERROR;

            int bsize = 10*ny*nz;
            int size2 = bsize / 5;
            MPI.Request[] mid = new MPI.Request[1];
            double[] buf1 = new double[bsize];
            double[] buf  = new double[bsize];

            if (iex==0) {
                //---------------------------------------------------------------------
                //   communicate in the south and north directions
                //---------------------------------------------------------------------
                if (north!=-1) {
                    //call MPI_IRECV[ buf1, 10*ny*nz, dp_type, MPI_ANY_SOURCE, from_n, MPI_COMM_WORLD, mid, IERROR ];
                    mid[0] = worldcomm.ImmediateReceive<double>(MPI.Unsafe.MPI_ANY_SOURCE,from_n,buf1);
                }
                //---------------------------------------------------------------------
                //   send south
                //---------------------------------------------------------------------
                if (south!=-1) {
                    for(k = 1; k<=nz; k++){
                      for(j = 1; j<=ny; j++){
                        ipos1 = (k-1)*ny+j              -1;  //ipos1 = (k-1)*ny+j;
                        ipos2 = ipos1 + ny*nz;               //ipos2 = ipos1 + ny*nz;
                        buf[0*size2+ipos1] = g[k,j+1,nx,1];  //buf[1,ipos1] = g[1,nx-1,j,k];
                        buf[1*size2+ipos1] = g[k,j+1,nx,2];  //buf[2,ipos1] = g[2,nx-1,j,k];
                        buf[2*size2+ipos1] = g[k,j+1,nx,3];  //buf[3,ipos1] = g[3,nx-1,j,k];
                        buf[3*size2+ipos1] = g[k,j+1,nx,4];  //buf[4,ipos1] = g[4,nx-1,j,k]; 
                        buf[4*size2+ipos1] = g[k,j+1,nx,5];  //buf[5,ipos1] = g[5,nx-1,j,k];

                        buf[0*size2+ipos2] = g[k,j+1,nx+1,1];    //buf[1,ipos2] = g[1,nx,j,k];
                        buf[1*size2+ipos2] = g[k,j+1,nx+1,2];    //buf[2,ipos2] = g[2,nx,j,k];
                        buf[2*size2+ipos2] = g[k,j+1,nx+1,3];    //buf[3,ipos2] = g[3,nx,j,k];
                        buf[3*size2+ipos2] = g[k,j+1,nx+1,4];    //buf[4,ipos2] = g[4,nx,j,k];
                        buf[4*size2+ipos2] = g[k,j+1,nx+1,5];    //buf[5,ipos2] = g[5,nx,j,k];
                      }
                    }
                    //call MPI_SEND[ buf, 10*ny*nz, dp_type, south, from_n, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, south, from_n);
                }
                //---------------------------------------------------------------------
                //   receive from north
                //---------------------------------------------------------------------
                if (north!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++){
                        for(j = 1; j<=ny; j++){
                            ipos1 = (k-1)*ny + j           -1;     //ipos1 = (k-1)*ny + j;
                            ipos2 = ipos1 + ny*nz;                 //ipos2 = ipos1 + ny*nz; 
                            g[k,j+1,0,1] = buf1[0*size2+ipos1];     //g[1,-1,j,k] = buf1[1,ipos1];       
                            g[k,j+1,0,2] = buf1[1*size2+ipos1];     //g[2,-1,j,k] = buf1[2,ipos1];
                            g[k,j+1,0,3] = buf1[2*size2+ipos1];     //g[3,-1,j,k] = buf1[3,ipos1];
                            g[k,j+1,0,4] = buf1[3*size2+ipos1];     //g[4,-1,j,k] = buf1[4,ipos1];
                            g[k,j+1,0,5] = buf1[4*size2+ipos1];     //g[5,-1,j,k] = buf1[5,ipos1];

                            g[k,j+1,1,1] = buf1[0*size2+ipos2];      //g[1,0,j,k] = buf1[1,ipos2];
                            g[k,j+1,1,2] = buf1[1*size2+ipos2];      //g[2,0,j,k] = buf1[2,ipos2];
                            g[k,j+1,1,3] = buf1[2*size2+ipos2];      //g[3,0,j,k] = buf1[3,ipos2]; 
                            g[k,j+1,1,4] = buf1[3*size2+ipos2];      //g[4,0,j,k] = buf1[4,ipos2];
                            g[k,j+1,1,5] = buf1[4*size2+ipos2];      //g[5,0,j,k] = buf1[5,ipos2];
                        }
                    }
                }
                if (south!=-1) {
                    //call MPI_IRECV[buf1, 10*ny*nz, dp_type, MPI_ANY_SOURCE, from_s, MPI_COMM_WORLD, mid, IERROR];
                    mid[0] = worldcomm.ImmediateReceive<double>(MPI.Unsafe.MPI_ANY_SOURCE, from_s, buf1);
                }
                //  c---------------------------------------------------------------------
                //  c   send north
                //  c---------------------------------------------------------------------
                if (north!=-1) {
                    for(k = 1; k<=nz; k++){
                        for(j = 1; j<=ny; j++){
                            ipos1 = (k-1)*ny + j   -1;          //ipos1 = (k-1)*ny + j;
                            ipos2 = ipos1 + ny*nz;              //ipos2 = ipos1 + ny*nz;
                            buf[0*size2+ipos1] = g[k,j+1,3,1];  //buf[1,ipos1] = g[1,2,j,k];
                            buf[1*size2+ipos1] = g[k,j+1,3,2];  //buf[2,ipos1] = g[2,2,j,k];
                            buf[2*size2+ipos1] = g[k,j+1,3,3];  //buf[3,ipos1] = g[3,2,j,k];
                            buf[3*size2+ipos1] = g[k,j+1,3,4];  //buf[4,ipos1] = g[4,2,j,k];
                            buf[4*size2+ipos1] = g[k,j+1,3,5];  //buf[5,ipos1] = g[5,2,j,k];

                            buf[0*size2+ipos2] = g[k,j+1,2,1];  //buf[1,ipos2] = g[1,1,j,k];
                            buf[1*size2+ipos2] = g[k,j+1,2,2];  //buf[2,ipos2] = g[2,1,j,k];
                            buf[2*size2+ipos2] = g[k,j+1,2,3];  //buf[3,ipos2] = g[3,1,j,k];
                            buf[3*size2+ipos2] = g[k,j+1,2,4];  //buf[4,ipos2] = g[4,1,j,k];
                            buf[4*size2+ipos2] = g[k,j+1,2,5];  //buf[5,ipos2] = g[5,1,j,k];
                        }
                    }
                    //call MPI_SEND[ buf, 10*ny*nz, dp_type, north, from_s, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, north, from_s);
                }
                //  c---------------------------------------------------------------------
                //  c   receive from south
                //  c---------------------------------------------------------------------
                if (south!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++){
                        for(j = 1; j<=ny; j++){
                            ipos1 = (k-1)*ny + j                -1; //ipos1 = (k-1)*ny + j;
                            ipos2 = ipos1 + ny*nz;                  //ipos2 = ipos1 + ny*nz;
                            g[k,j+1,nx+3,1]  = buf1[0*size2+ipos1]; //g[1,nx+2,j,k]  = buf1[1,ipos1];
                            g[k,j+1,nx+3,2]  = buf1[1*size2+ipos1]; //g[2,nx+2,j,k]  = buf1[2,ipos1];
                            g[k,j+1,nx+3,3]  = buf1[2*size2+ipos1]; //g[3,nx+2,j,k]  = buf1[3,ipos1];
                            g[k,j+1,nx+3,4]  = buf1[3*size2+ipos1]; //g[4,nx+2,j,k]  = buf1[4,ipos1];
                            g[k,j+1,nx+3,5]  = buf1[4*size2+ipos1]; //g[5,nx+2,j,k]  = buf1[5,ipos1];

                            g[k,j+1,nx+2,1] = buf1[0*size2+ipos2];  //g[1,nx+1,j,k] = buf1[1,ipos2];
                            g[k,j+1,nx+2,2] = buf1[1*size2+ipos2];  //g[2,nx+1,j,k] = buf1[2,ipos2];
                            g[k,j+1,nx+2,3] = buf1[2*size2+ipos2];  //g[3,nx+1,j,k] = buf1[3,ipos2];
                            g[k,j+1,nx+2,4] = buf1[3*size2+ipos2];  //g[4,nx+1,j,k] = buf1[4,ipos2];
                            g[k,j+1,nx+2,5] = buf1[4*size2+ipos2];  //g[5,nx+1,j,k] = buf1[5,ipos2];
                        }
                    }
                }
            } else {
                bsize = 10*nx*nz;
                size2 = bsize/5;
                buf1 = new double[bsize];
                buf  = new double[bsize];
                //---------------------------------------------------------------------
                //   communicate in the east and west directions
                //---------------------------------------------------------------------
                if (west!=-1) {
                    //call MPI_IRECV[ buf1, 10*nx*nz, dp_type, MPI_ANY_SOURCE, from_w, MPI_COMM_WORLD, mid, IERROR ];
                    mid[0] = worldcomm.ImmediateReceive<double>(MPI.Unsafe.MPI_ANY_SOURCE, from_w, buf1);
                }
                //---------------------------------------------------------------------
                //   send east
                //---------------------------------------------------------------------
                if (east!=-1) {
                    for(k = 1; k<=nz; k++){
                        for(i = 1; i<=nx; i++){
                            ipos1 = (k-1)*nx+i        -1;         //ipos1 = (k-1)*nx+i;
                            ipos2 = ipos1+nx*nz;                  //ipos2 = ipos1+nx*nz;
                            buf[0*size2+ipos1] = g[k,ny,i+1,1];   //buf[1,ipos1] = g[1,i,ny-1,k];
                            buf[1*size2+ipos1] = g[k,ny,i+1,2];   //buf[2,ipos1] = g[2,i,ny-1,k];
                            buf[2*size2+ipos1] = g[k,ny,i+1,3];   //buf[3,ipos1] = g[3,i,ny-1,k];
                            buf[3*size2+ipos1] = g[k,ny,i+1,4];   //buf[4,ipos1] = g[4,i,ny-1,k];
                            buf[4*size2+ipos1] = g[k,ny,i+1,5];   //buf[5,ipos1] = g[5,i,ny-1,k];

                            buf[0*size2+ipos2] = g[k,ny+1,i+1,1];     //buf[1,ipos2] = g[1,i,ny,k];
                            buf[1*size2+ipos2] = g[k,ny+1,i+1,2];     //buf[2,ipos2] = g[2,i,ny,k];
                            buf[2*size2+ipos2] = g[k,ny+1,i+1,3];     //buf[3,ipos2] = g[3,i,ny,k];
                            buf[3*size2+ipos2] = g[k,ny+1,i+1,4];     //buf[4,ipos2] = g[4,i,ny,k];
                            buf[4*size2+ipos2] = g[k,ny+1,i+1,5];     //buf[5,ipos2] = g[5,i,ny,k];
                        }
                    }
                    //call MPI_SEND[ buf, 10*nx*nz, dp_type, east, from_w, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, east, from_w);
                }
                //---------------------------------------------------------------------
                //   receive from west
                //---------------------------------------------------------------------
                if (west!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++){
                        for(i = 1; i<=nx; i++){
                            ipos1 = (k-1)*nx + i     -1;           //ipos1 = (k-1)*nx + i;
                            ipos2 = ipos1 + nx*nz;                 //ipos2 = ipos1 + nx*nz;
                            g[k,0,i+1,1] = buf1[0*size2+ipos1];    //g[1,i,-1,k] = buf1[1,ipos1];
                            g[k,0,i+1,2] = buf1[1*size2+ipos1];    //g[2,i,-1,k] = buf1[2,ipos1];
                            g[k,0,i+1,3] = buf1[2*size2+ipos1];    //g[3,i,-1,k] = buf1[3,ipos1];
                            g[k,0,i+1,4] = buf1[3*size2+ipos1];    //g[4,i,-1,k] = buf1[4,ipos1];
                            g[k,0,i+1,5] = buf1[4*size2+ipos1];    //g[5,i,-1,k] = buf1[5,ipos1];

                            g[k,1,i+1,1] = buf1[0*size2+ipos2];     //g[1,i,0,k] = buf1[1,ipos2];
                            g[k,1,i+1,2] = buf1[1*size2+ipos2];     //g[2,i,0,k] = buf1[2,ipos2];
                            g[k,1,i+1,3] = buf1[2*size2+ipos2];     //g[3,i,0,k] = buf1[3,ipos2];
                            g[k,1,i+1,4] = buf1[3*size2+ipos2];     //g[4,i,0,k] = buf1[4,ipos2];
                            g[k,1,i+1,5] = buf1[4*size2+ipos2];     //g[5,i,0,k] = buf1[5,ipos2];
                        }
                    }
                }
                if (east!=-1) {
                    //call MPI_IRECV[ buf1, 10*nx*nz, dp_type, MPI_ANY_SOURCE, from_e, MPI_COMM_WORLD, mid, IERROR ];
                    mid[0] = worldcomm.ImmediateReceive<double>(MPI.Unsafe.MPI_ANY_SOURCE, from_e, buf1);
                }
                //---------------------------------------------------------------------
                //   send west
                //---------------------------------------------------------------------
                if (west!=-1) {
                    for(k = 1; k<=nz; k++){
                        for(i = 1; i<=nx; i++){
                            ipos1 = (k-1)*nx + i   -1;          //ipos1 = (k-1)*nx + i;
                            ipos2 = ipos1 + nx*nz;              //ipos2 = ipos1 + nx*nz;
                            buf[0*size2+ipos1] = g[k,3,i+1,1];  //buf[1,ipos1] = g[1,i,2,k];
                            buf[1*size2+ipos1] = g[k,3,i+1,2];  //buf[2,ipos1] = g[2,i,2,k];
                            buf[2*size2+ipos1] = g[k,3,i+1,3];  //buf[3,ipos1] = g[3,i,2,k];
                            buf[3*size2+ipos1] = g[k,3,i+1,4];  //buf[4,ipos1] = g[4,i,2,k];
                            buf[4*size2+ipos1] = g[k,3,i+1,5];  //buf[5,ipos1] = g[5,i,2,k];

                            buf[0*size2+ipos2] = g[k,2,i+1,1];  //buf[1,ipos2] = g[1,i,1,k];
                            buf[1*size2+ipos2] = g[k,2,i+1,2];  //buf[2,ipos2] = g[2,i,1,k];
                            buf[2*size2+ipos2] = g[k,2,i+1,3];  //buf[3,ipos2] = g[3,i,1,k];
                            buf[3*size2+ipos2] = g[k,2,i+1,4];  //buf[4,ipos2] = g[4,i,1,k];
                            buf[4*size2+ipos2] = g[k,2,i+1,5];  //buf[5,ipos2] = g[5,i,1,k];
                        }
                    }
                    //call MPI_SEND[buf, 10*nx*nz, dp_type, west, from_e, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, west, from_e);
                }
                //---------------------------------------------------------------------
                //   receive from east
                //---------------------------------------------------------------------
                if (east!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++){
                        for(i = 1; i<=nx; i++){
                            ipos1 = (k-1)*nx + i        -1;         //ipos1 = (k-1)*nx + i;
                            ipos2 = ipos1 + nx*nz;                  //ipos2 = ipos1 + nx*nz;
                            g[k,ny+3,i+1,1]  = buf1[0*size2+ipos1]; //g[1,i,ny+2,k]  = buf1[1,ipos1];
                            g[k,ny+3,i+1,2]  = buf1[1*size2+ipos1]; //g[2,i,ny+2,k]  = buf1[2,ipos1];
                            g[k,ny+3,i+1,3]  = buf1[2*size2+ipos1]; //g[3,i,ny+2,k]  = buf1[3,ipos1];
                            g[k,ny+3,i+1,4]  = buf1[3*size2+ipos1]; //g[4,i,ny+2,k]  = buf1[4,ipos1];
                            g[k,ny+3,i+1,5]  = buf1[4*size2+ipos1]; //g[5,i,ny+2,k]  = buf1[5,ipos1];

                            g[k,ny+2,i+1,1] = buf1[0*size2+ipos2];  //g[1,i,ny+1,k] = buf1[1,ipos2];
                            g[k,ny+2,i+1,2] = buf1[1*size2+ipos2];  //g[2,i,ny+1,k] = buf1[2,ipos2];
                            g[k,ny+2,i+1,3] = buf1[2*size2+ipos2];  //g[3,i,ny+1,k] = buf1[3,ipos2];
                            g[k,ny+2,i+1,4] = buf1[3*size2+ipos2];  //g[4,i,ny+1,k] = buf1[4,ipos2];
                            g[k,ny+2,i+1,5] = buf1[4*size2+ipos2];  //g[5,i,ny+1,k] = buf1[5,ipos2];
                        }
                    }
                }
            }
        }
           //end Exchange_3.f
        //End erhs.f

        //ssor.f
        public void ssor(int niter) {
            //c---------------------------------------------------------------------
            //c   to perform pseudo-time stepping SSOR iterations
            //c   for five nonlinear pde's.
            //c---------------------------------------------------------------------
            //int  niter;
            int i, j, k, m;
            int istep;
            double  tmp;
            double[] delunm = new double[5+1];//delunm[5];
            double[,,] tv = new double[isiz2+1,isiz1+1,5+1];//tv[5,isiz1,isiz2];
            //external timer_read;
            //double wtime, timer_read;
            double wtime;
            //int IERROR;
            //ROOT = 0;
            //---------------------------------------------------------------------
            //   begin pseudo-time stepping iterations
            //---------------------------------------------------------------------
            tmp = 1.0d/(omega*(2.0d-omega));
            //---------------------------------------------------------------------
            //   initialize a,b,c,d to zero [guarantees that page tables have been
            //   formed, if applicable on given architecture, before timestepping].
            //---------------------------------------------------------------------
            for(m=1; m<=isiz2; m++){
               for(k=1; k<=isiz1; k++){
                  for(j=1; j<=5; j++){
                     for(i=1; i<=5; i++){
                        a[m,k,j,i] = 0.0;//a[i,j,k,m] = 0.0;
                        b[m,k,j,i] = 0.0;//b[i,j,k,m] = 0.0;
                        c[m,k,j,i] = 0.0;//c[i,j,k,m] = 0.0;
                        d[m,k,j,i] = 0.0;//d[i,j,k,m] = 0.0;
                     }
                  }
               }
            }
            //---------------------------------------------------------------------
            //   compute the steady-state residuals
            //---------------------------------------------------------------------
            rhs();
            //---------------------------------------------------------------------
            //   compute the L2 norms of newton iteration residuals
            //---------------------------------------------------------------------
            l2norm(isiz1, isiz2, isiz3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
            worldcomm.Barrier(); //call MPI_BARRIER[ MPI_COMM_WORLD, IERROR ]
            timer.resetTimer(1); //call timer_clear[1]
            timer.start(1);      //call timer_start[1]
            //---------------------------------------------------------------------
            //   the timestep loop
            //---------------------------------------------------------------------
            for(istep = 1; istep<= niter; istep++){
               if (id == 0) {
                  if (mod(istep, 20) == 0 || istep == itmax || istep == 1) {
                     if (niter > 1) Console.WriteLine(" Time step "+istep); //write[ *, 200] istep       200           format[' Time step ', i4]
                  }
               }
               //---------------------------------------------------------------------
               //   perform SSOR iteration
               //---------------------------------------------------------------------
               for(k = 2; k<= nz - 1; k++){
                  for(j = jst; j<= jend; j++){
                     for(i = ist; i<= iend; i++){
                        for(m = 1; m<= 5; m++){
                            rsd[k, j+1, i+1, m] = dt * rsd[k, j+1, i+1, m];   //rsd[m, i, j, k] = dt * rsd[m, i, j, k];
                        }
                     }
                  }
               }
               for(k = 2; k<= nz -1; k++){
                    //---------------------------------------------------------------------
                    //   form the lower triangular part of the jacobian matrix
                    //---------------------------------------------------------------------
                    jacld(k);
                    //---------------------------------------------------------------------
                    //   perform the lower triangular solution
                    //---------------------------------------------------------------------
                    blts(isiz1, isiz2, isiz3,nx, ny, nz, k,omega,rsd,a, b, c, d,ist, iend, jst, jend, nx0, ny0, ipt, jpt);
                }
                for(k=nz-1; k>= 2; k--){ //for(k = nz - 1, 2, -1;
                  //c---------------------------------------------------------------------
                  //c   form the strictly upper triangular part of the jacobian matrix
                  //c---------------------------------------------------------------------
                  jacu(k);
                  //c---------------------------------------------------------------------
                  //c   perform the upper triangular solution
                  //c---------------------------------------------------------------------
                  buts(isiz1, isiz2, isiz3, nx, ny, nz, k, omega, rsd, tv, d, a, b, c, ist, iend, jst, jend, nx0, ny0, ipt, jpt);
                }
                //---------------------------------------------------------------------
                //   update the variables
                //---------------------------------------------------------------------
                for(k = 2; k<= nz-1; k++){
                    for(j = jst; j<= jend; j++){
                        for(i = ist; i<= iend; i++){
                            for(m = 1; m<= 5; m++){     //u[ m, i, j, k ] = u[ m, i, j, k ] + tmp * rsd[ m, i, j, k ];
                                u[k,j+1,i+1,m] = u[k,j+1,i+1,m] + tmp * rsd[k,j+1,i+1,m];
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //   compute the max-norms of newton iteration corrections
                //---------------------------------------------------------------------
                if (mod(istep, inorm) == 0) {
                   l2norm(isiz1, isiz2, isiz3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, delunm);
                    //c            if [ ipr == 1 && id == 0 ] {
                    //c                write [*,1006] [ delunm[m], m = 1, 5 ]
                    //c            } else if [ ipr == 2 && id == 0 ] {
                    //c                write [*,'[i5,f15.6]'] istep,delunm[5]
                    //c            }
                }
                //---------------------------------------------------------------------
                //   compute the steady-state residuals
                //---------------------------------------------------------------------
                rhs();
                //---------------------------------------------------------------------
                //   compute the max-norms of newton iteration residuals
                //---------------------------------------------------------------------
                if ((mod(istep,inorm)== 0) || (istep==itmax)) {
                    l2norm(isiz1, isiz2, isiz3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
                    //c            if [ ipr == 1&&id==0 ] {
                    //c                write [*,1007] [ rsdnm[m], m = 1, 5 ]
                    //c            }
                }
                //---------------------------------------------------------------------
                //   check the newton-iteration residuals against the tolerance levels
                //---------------------------------------------------------------------
                if ((rsdnm[1]<tolrsd[1]) && (rsdnm[2]<tolrsd[2]) && (rsdnm[3]<tolrsd[3]) && (rsdnm[4]<tolrsd[4]) && (rsdnm[5]<tolrsd[5])) {
                    //c            if [ipr == 1 && id==0] {
                    //c               write [*,1004] istep
                    //c            }
                    return;//   return;
                }
            }
            timer.stop(1); //call timer_stop[1];
            wtime = timer.readTimer(1); //wtime = timer_read[1];
            //call MPI_ALLREDUCE[wtime, maxtime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR ];
            maxtime = worldcomm.Allreduce<double>(wtime, MPI.Operation<double>.Max);
            //return;
        }
            //rhs.f
        public void rhs() {
            //---------------------------------------------------------------------
            //   compute the right hand sides
            //---------------------------------------------------------------------
            int i, j, k, m;
            int iex;
            int L1, L2;
            int ist1, iend1;
            int jst1, jend1;
            double q;
            double u21, u31, u41;
            double tmp;
            double u21i, u31i, u41i, u51i;
            double u21j, u31j, u41j, u51j;
            double u21k, u31k, u41k, u51k;
            double u21im1, u31im1, u41im1, u51im1;
            double u21jm1, u31jm1, u41jm1, u51jm1;
            double u21km1, u31km1, u41km1, u51km1;

            for(k = 1; k<= nz; k++){
               for(j = 1; j<= ny; j++){
                  for(i = 1; i<= nx; i++){
                     for(m = 1; m<= 5; m++){
                        rsd[k,j+1,i+1,m] = -frct[k,j+1,i+1,m];//rsd[m,i,j,k] = - frct[m,i,j,k];
                     }
                  }
               }
            }
            //---------------------------------------------------------------------
            //   xi-direction flux differences
            //---------------------------------------------------------------------
            //c---------------------------------------------------------------------
            //c   iex = flag : iex = 0  north/south communication
            //c              : iex = 1  east/west communication
            //c---------------------------------------------------------------------
            iex   = 0;
            //---------------------------------------------------------------------
            //   communicate and receive/send two rows of data
            //---------------------------------------------------------------------
            exchange_3(u,iex);
            L1 = 0;
            if (north==-1) L1 = 1;
            L2 = nx + 1;
            if (south==-1) L2 = nx;

            ist1 = 1;
            iend1 = nx;
            if (north==-1) ist1 = 4;
            if (south==-1) iend1 = nx - 3;

            for(k = 2; k<= nz - 1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = L1; i<= L2; i++){
                     flux[k,j,i,1] = u[k,j+1,i+1,2];  //flux[1,i,j,k] = u[2,i,j,k];
                     u21=u[k,j+1,i+1,2]/u[k,j+1,i+1,1];  //u21 = u[2,i,j,k] / u[1,i,j,k];

                     q = 0.50d*(u[k,j+1,i+1,2]*u[k,j+1,i+1,2]+u[k,j+1,i+1,3]*u[k,j+1,i+1,3]+u[k,j+1,i+1,4]*u[k,j+1,i+1,4])/u[k,j+1,i+1,1];//q = 0.50d*(u[2,i,j,k]*u[2,i,j,k]+u[3,i,j,k]*u[3,i,j,k]+u[4,i,j,k]*u[4,i,j,k])/u[1,i,j,k];

                     flux[k,j,i,2] =     u[k,j+1,i+1,2] * u21 + c2 *(u[k,j+1,i+1,5] - q);   //flux[2,i,j,k]=u[2,i,j,k]*u21+c2*(u[5,i,j,k]-q);
                     flux[k,j,i,3] =     u[k,j+1,i+1,3] * u21;                              //flux[3,i,j,k]=u[3,i,j,k]*u21;
                     flux[k,j,i,4] =     u[k,j+1,i+1,4] * u21;                              //flux[4,i,j,k]=u[4,i,j,k]*u21;
                     flux[k,j,i,5] = (c1*u[k,j+1,i+1,5]-c2*q)*u21;                          //flux[5,i,j,k]=(c1*u[5,i,j,k]-c2*q)*u21;
                  }
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        rsd[k,j+1,i+1,m]=rsd[k,j+1,i+1,m]-tx2*(flux[k,j,i+1,m] - flux[k,j,i-1,m]);  //rsd[m,i,j,k] =  rsd[m,i,j,k]- tx2 *(flux[m,i+1,j,k] - flux[m,i-1,j,k]);
                     }
                  }
                  for(i = ist; i<= L2; i++){
                     tmp = 1.0d/u[k,j+1,i+1,1];
                     u21i = tmp*u[k,j+1,i+1,2];
                     u31i = tmp*u[k,j+1,i+1,3];
                     u41i = tmp*u[k,j+1,i+1,4];
                     u51i = tmp*u[k,j+1,i+1,5];

                     tmp = 1.0d/ u[k,j+1,i,1];

                     u21im1 = tmp * u[k,j+1,i,2];
                     u31im1 = tmp * u[k,j+1,i,3];
                     u41im1 = tmp * u[k,j+1,i,4];
                     u51im1 = tmp * u[k,j+1,i,5];

                     flux[k,j,i,2] = (4.0d/3.0d)*tx3*(u21i-u21im1);
                     flux[k,j,i,3] = tx3 * ( u31i - u31im1 );
                     flux[k,j,i,4] = tx3 * ( u41i - u41im1 );
                     flux[k,j,i,5] = 0.50d*(1.0d-c1*c5)*tx3*((pow2(u21i)+pow2(u31i)+pow2(u41i))-(pow2(u21im1)+pow2(u31im1)+pow2(u41im1)))+(1.0d/6.0d)*tx3*(pow2(u21i)-pow2(u21im1))+c1*c5*tx3*(u51i-u51im1);
                  }
                  for(i = ist; i<= iend; i++){
                     rsd[k,j+1,i+1,1]=rsd[k,j+1,i+1,1]+dx1*tx1*(u[k,j+1,i,1]-2.0d*u[k,j+1,i+1,1]+u[k,j+1,i+2,1]);
                     rsd[k,j+1,i+1,2]=rsd[k,j+1,i+1,2]+tx3*c3*c4*(flux[k,j,i+1,2]-flux[k,j,i,2])+dx2*tx1*(u[k,j+1,i,2]-2.0d*u[k,j+1,i+1,2]+u[k,j+1,i+2,2]);
                     rsd[k,j+1,i+1,3]=rsd[k,j+1,i+1,3]+tx3*c3*c4*(flux[k,j,i+1,3]-flux[k,j,i,3])+dx3*tx1*(u[k,j+1,i,3]-2.0d*u[k,j+1,i+1,3]+u[k,j+1,i+2,3]);
                     rsd[k,j+1,i+1,4]=rsd[k,j+1,i+1,4]+tx3*c3*c4*(flux[k,j,i+1,4]-flux[k,j,i,4])+dx4*tx1*(u[k,j+1,i,4]-2.0d*u[k,j+1,i+1,4]+u[k,j+1,i+2,4]);
                     rsd[k,j+1,i+1,5]=rsd[k,j+1,i+1,5]+tx3*c3*c4*(flux[k,j,i+1,5]-flux[k,j,i,5])+dx5*tx1*(u[k,j+1,i,5]-2.0d*u[k,j+1,i+1,5]+u[k,j+1,i+2,5]);
                  }
                  //---------------------------------------------------------------------
                  //   Fourth-order dissipation
                  //---------------------------------------------------------------------
                  if (north==-1) {
                   for(m = 1; m<= 5; m++){
                     rsd[k,j+1,3,m] = rsd[k,j+1,3,m]-dssp*(+5.0d*u[k,j+1,3,m]-4.0d*u[k,j+1,4,m]+     u[k,j+1,5,m]);
                     rsd[k,j+1,4,m] = rsd[k,j+1,4,m]-dssp*(-4.0d*u[k,j+1,3,m]+6.0d*u[k,j+1,4,m]-4.0d*u[k,j+1,5,m]+u[k,j+1,6,m]);
                   }
                  }
                  for(i = ist1; i<=iend1; i++){
                     for(m = 1; m<= 5; m++){
                        rsd[k,j+1,i+1,m] = rsd[k,j+1,i+1,m]-dssp*(u[k,j+1,i-1,m]-4.0d*u[k,j+1,i,m]+6.0d*u[k,j+1,i+1,m]-4.0d*u[k,j+1,i+2,m]+u[k,j+1,i+3,m]);
                     }
                  }
                  if (south==-1) {
                   for(m = 1; m<= 5; m++){
                     rsd[k,j+1,nx-1,m] = rsd[k,j+1,nx-1,m]-dssp*(u[k,j+1,nx-3,m]-4.0d*u[k,j+1,nx-2,m]+6.0d*u[k,j+1,nx-1,m]-4.0d*u[k,j+1,nx,m]);
                     rsd[k,j+1,nx,m]   = rsd[k,j+1,nx,m]  -dssp*(u[k,j+1,nx-2,m]-4.0d*u[k,j+1,nx-1,m]+5.0d*u[k,j+1,nx,m]);
                   }
                  }
               }
            } 
            //---------------------------------------------------------------------
            //   eta-direction flux differences
            //---------------------------------------------------------------------
            //c---------------------------------------------------------------------
            //c   iex = flag : iex = 0  north/south communication
            //c---------------------------------------------------------------------
            iex   = 1;
            //---------------------------------------------------------------------
            //   communicate and receive/send two rows of data
            //---------------------------------------------------------------------
            exchange_3(u,iex);

            L1 = 0;
            if (west==-1) L1 = 1;
            L2 = ny + 1;
            if (east==-1) L2 = ny;

            jst1 = 1;
            jend1 = ny;
            if (west==-1) jst1 = 4;
            if (east==-1) jend1 = ny - 3;

            for(k = 2; k<= nz - 1; k++){
               for(j = L1; j<= L2; j++){
                  for(i = ist; i<= iend; i++){
                     flux[k,j,i,1]=u[k,j+1,i+1,3];          //flux[1,i,j,k] = u[3,i,j,k];
                     u31=u[k,j+1,i+1,3]/u[k,j+1,i+1,1];     //u31 = u[3,i,j,k] / u[1,i,j,k];

                     q = 0.50d*(u[k,j+1,i+1,2]*u[k,j+1,i+1,2]+u[k,j+1,i+1,3]*u[k,j+1,i+1,3]+u[k,j+1,i+1,4]*u[k,j+1,i+1,4])/u[k,j+1,i+1,1];      //q = 0.50d*(u[2,i,j,k]*u[2,i,j,k]+u[3,i,j,k]*u[3,i,j,k]+u[4,i,j,k]*u[4,i,j,k])/u[1,i,j,k];

                     flux[k,j,i,2] =     u[k,j+1,i+1,2] * u31;                              //flux[2,i,j,k]=u[2,i,j,k] * u31;
                     flux[k,j,i,3] =     u[k,j+1,i+1,3] * u31 + c2 * (u[k,j+1,i+1,5]-q);    //flux[3,i,j,k]=u[3,i,j,k]*u31+c2*(u[5,i,j,k]-q);
                     flux[k,j,i,4] =     u[k,j+1,i+1,4] * u31;                              //flux[4,i,j,k] =     u[4,i,j,k] * u31;
                     flux[k,j,i,5] = (c1*u[k,j+1,i+1,5]-c2*q)*u31;                          //flux[5,i,j,k] = (c1*u[5,i,j,k]-c2*q)*u31;
                  }
               }
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){  //rsd[m,i,j,k] =  rsd[m,i,j,k]- ty2 * [ flux[m,i,j+1,k] - flux[m,i,j-1,k] ];
                        rsd[k,j+1,i+1,m] =  rsd[k,j+1,i+1,m]- ty2 * (flux[k,j+1,i,m] - flux[k,j-1,i,m]);
                     }
                  }
               }
               for(j = jst; j<= L2; j++){
                  for(i = ist; i<= iend; i++){
                     tmp = 1.0d / u[k,j+1,i+1,1];
                     u21j = tmp * u[k,j+1,i+1,2];
                     u31j = tmp * u[k,j+1,i+1,3];
                     u41j = tmp * u[k,j+1,i+1,4];
                     u51j = tmp * u[k,j+1,i+1,5];

                     tmp = 1.0d / u[k,j,i+1,1];
                     u21jm1=tmp * u[k,j,i+1,2];
                     u31jm1=tmp * u[k,j,i+1,3];
                     u41jm1=tmp * u[k,j,i+1,4];
                     u51jm1=tmp * u[k,j,i+1,5];

                     flux[k,j,i,2] = ty3 * (u21j - u21jm1);
                     flux[k,j,i,3] = (4.0d/3.0d) * ty3 * (u31j-u31jm1);
                     flux[k,j,i,4] = ty3 * (u41j - u41jm1);
                     flux[k,j,i,5] = 0.50d*(1.0d-c1*c5)*ty3*((pow2(u21j)+pow2(u31j)+pow2(u41j))-(pow2(u21jm1)+pow2(u31jm1)+pow2(u41jm1)))+(1.0d/6.0d)*ty3*(pow2(u31j)-pow2(u31jm1))+c1*c5*ty3*(u51j-u51jm1);
                  }
               }
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     rsd[k,j+1,i+1,1] = rsd[k,j+1,i+1,1]+dy1*ty1*(u[k,j,i+1,1]-2.0d*u[k,j+1,i+1,1]+u[k,j+2,i+1,1]);
                     rsd[k,j+1,i+1,2] = rsd[k,j+1,i+1,2]+ty3*c3*c4*(flux[k,j+1,i,2]-flux[k,j,i,2])+dy2*ty1*(u[k,j,i+1,2]-2.0d*u[k,j+1,i+1,2]+u[k,j+2,i+1,2]);
                     rsd[k,j+1,i+1,3] = rsd[k,j+1,i+1,3]+ty3*c3*c4*(flux[k,j+1,i,3]-flux[k,j,i,3])+dy3*ty1*(u[k,j,i+1,3]-2.0d*u[k,j+1,i+1,3]+u[k,j+2,i+1,3]);
                     rsd[k,j+1,i+1,4] = rsd[k,j+1,i+1,4]+ty3*c3*c4*(flux[k,j+1,i,4]-flux[k,j,i,4])+dy4*ty1*(u[k,j,i+1,4]-2.0d*u[k,j+1,i+1,4]+u[k,j+2,i+1,4]);
                     rsd[k,j+1,i+1,5] = rsd[k,j+1,i+1,5]+ty3*c3*c4*(flux[k,j+1,i,5]-flux[k,j,i,5])+dy5*ty1*(u[k,j,i+1,5]-2.0d*u[k,j+1,i+1,5]+u[k,j+2,i+1,5]);
                  }
               }
               //---------------------------------------------------------------------
               //   fourth-order dissipation
               //---------------------------------------------------------------------
               if (west==-1) {
                  for(i = ist; i<= iend; i++){
                   for(m = 1; m<= 5; m++){
                     rsd[k,3,i+1,m] = rsd[k,3,i+1,m]-dssp*(+5.0d*u[k,3,i+1,m]-4.0d*u[k,4,i+1,m]+       u[k,5,i+1,m]);
                     rsd[k,4,i+1,m] = rsd[k,4,i+1,m]-dssp*(-4.0d*u[k,3,i+1,m]+6.0d*u[k,4,i+1,m]-4.0d * u[k,5,i+1,m]+u[k,6,i+1,m]);
                   }
                  }
               }
               for(j = jst1; j<= jend1; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        rsd[k,j+1,i+1,m]=rsd[k,j+1,i+1,m]-dssp*(u[k,j-1,i+1,m]-4.0d*u[k,j,i+1,m]+6.0d*u[k,j+1,i+1,m]-4.0d*u[k,j+2,i+1,m]+u[k,j+3,i+1,m]);
                     }
                  }
               }
               if (east==-1) {
                  for(i = ist; i<= iend; i++){
                   for(m = 1; m<= 5; m++){
                     rsd[k,ny-1,i+1,m]=rsd[k,ny-1,i+1,m]-dssp*(u[k,ny-3,i+1,m]-4.0d*u[k,ny-2,i+1,m]+6.0d*u[k,ny-1,i+1,m]-4.0d*u[k,ny,i+1,m]);
                     rsd[k,ny,  i+1,m]=rsd[k,ny,  i+1,m]-dssp*(u[k,ny-2,i+1,m]-4.0d*u[k,ny-1,i+1,m]+5.0d*u[k,ny  ,i+1,m]);
                   }
                  }
               }
            }
            //---------------------------------------------------------------------
            //   zeta-direction flux differences
            //---------------------------------------------------------------------
            for(k = 1; k<= nz; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     flux[k,j,i,1]=u[k,j+1,i+1,4];
                     u41=u[k,j+1,i+1,4]/u[k,j+1,i+1,1];

                     q = 0.50d * (u[k,j+1,i+1,2] * u[k,j+1,i+1,2]+ u[k,j+1,i+1,3] * u[k,j+1,i+1,3]+ u[k,j+1,i+1,4] * u[k,j+1,i+1,4])/u[k,j+1,i+1,1];

                     flux[k,j,i,2] =   u[k,j+1,i+1,2] * u41;
                     flux[k,j,i,3] =   u[k,j+1,i+1,3] * u41;
                     flux[k,j,i,4] =   u[k,j+1,i+1,4] * u41 + c2 * (u[k,j+1,i+1,5]-q);
                     flux[k,j,i,5]=(c1*u[k,j+1,i+1,5]-c2*q)*u41;
                  }
               }
            }
            for(k = 2; k<= nz - 1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        rsd[k,j+1,i+1,m] =  rsd[k,j+1,i+1,m]- tz2 * (flux[k+1,j,i,m] - flux[k-1,j,i,m]);
                     }
                  }
               }
            }
            for(k = 2; k<= nz; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     tmp = 1.0d / u[k,j+1,i+1,1];
                     u21k = tmp * u[k,j+1,i+1,2];
                     u31k = tmp * u[k,j+1,i+1,3];
                     u41k = tmp * u[k,j+1,i+1,4];
                     u51k = tmp * u[k,j+1,i+1,5];

                     tmp   = 1.0d / u[k-1,j+1,i+1,1];
                     u21km1 = tmp * u[k-1,j+1,i+1,2];
                     u31km1 = tmp * u[k-1,j+1,i+1,3];
                     u41km1 = tmp * u[k-1,j+1,i+1,4];
                     u51km1 = tmp * u[k-1,j+1,i+1,5];

                     flux[k,j,i,2] = tz3 * (u21k - u21km1);
                     flux[k,j,i,3] = tz3 * (u31k - u31km1);
                     flux[k,j,i,4] = (4.0d/3.0d) * tz3 * (u41k-u41km1);
                     flux[k,j,i,5] = 0.50d*(1.0d-c1*c5)*tz3*((pow2(u21k)+pow2(u31k)+pow2(u41k))-(pow2(u21km1)+pow2(u31km1)+pow2(u41km1)))+(1.0d/6.0d)*tz3*(pow2(u41k)-pow2(u41km1))+c1*c5*tz3*(u51k-u51km1);
                  }
               }
            }
            for(k = 2; k<= nz - 1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     rsd[k,j+1,i+1,1] = rsd[k,j+1,i+1,1]+dz1*tz1*(u[k-1,j+1,i+1,1]-2.0d*u[k,j+1,i+1,1]+u[k+1,j+1,i+1,1]);
                     rsd[k,j+1,i+1,2] = rsd[k,j+1,i+1,2]+tz3*c3*c4*(flux[k+1,j,i,2]-flux[k,j,i,2])+dz2*tz1*(u[k-1,j+1,i+1,2]-2.0d*u[k,j+1,i+1,2]+u[k+1,j+1,i+1,2]);
                     rsd[k,j+1,i+1,3] = rsd[k,j+1,i+1,3]+tz3*c3*c4*(flux[k+1,j,i,3]-flux[k,j,i,3])+dz3*tz1*(u[k-1,j+1,i+1,3]-2.0d*u[k,j+1,i+1,3]+u[k+1,j+1,i+1,3]);
                     rsd[k,j+1,i+1,4] = rsd[k,j+1,i+1,4]+tz3*c3*c4*(flux[k+1,j,i,4]-flux[k,j,i,4])+dz4*tz1*(u[k-1,j+1,i+1,4]-2.0d*u[k,j+1,i+1,4]+u[k+1,j+1,i+1,4]);
                     rsd[k,j+1,i+1,5] = rsd[k,j+1,i+1,5]+tz3*c3*c4*(flux[k+1,j,i,5]-flux[k,j,i,5])+dz5*tz1*(u[k-1,j+1,i+1,5]-2.0d*u[k,j+1,i+1,5]+u[k+1,j+1,i+1,5]);
                  }
               }
            }
            //---------------------------------------------------------------------
            //   fourth-order dissipation
            //---------------------------------------------------------------------
            for(j = jst; j<= jend; j++){
               for(i = ist; i<= iend; i++){
                  for(m = 1; m<= 5; m++){
                     rsd[2,j+1,i+1,m]=rsd[2,j+1,i+1,m]-dssp*(+ 5.0d*u[2,j+1,i+1,m]-4.0d*u[3,j+1,i+1,m]+     u[4,j+1,i+1,m]);
                     rsd[3,j+1,i+1,m]=rsd[3,j+1,i+1,m]-dssp*(- 4.0d*u[2,j+1,i+1,m]+6.0d*u[3,j+1,i+1,m]-4.0d*u[4,j+1,i+1,m]+u[5,j+1,i+1,m]);
                  }
               }
            }
            for(k = 4; k<= nz - 3; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        rsd[k,j+1,i+1,m] = rsd[k,j+1,i+1,m]-dssp*(u[k-2,j+1,i+1,m]-4.0d*u[k-1,j+1,i+1,m]+6.0d*u[k,j+1,i+1,m]-4.0d*u[k+1,j+1,i+1,m]+u[k+2,j+1,i+1,m]);
                     }
                  }
               }
            }
            for(j = jst; j<= jend; j++){
               for(i = ist; i<= iend; i++){
                  for(m = 1; m<= 5; m++){
                     rsd[nz-2,j+1,i+1,m]=rsd[nz-2,j+1,i+1,m]-dssp*(u[nz-4,j+1,i+1,m]-4.0d*u[nz-3,j+1,i+1,m]+6.0d*u[nz-2,j+1,i+1,m]-4.0d*u[nz-1,j+1,i+1,m]);
                     rsd[nz-1,j+1,i+1,m]=rsd[nz-1,j+1,i+1,m]-dssp*(u[nz-3,j+1,i+1,m]-4.0d*u[nz-2,j+1,i+1,m]+5.0d*u[nz-1,j+1,i+1,m]);
                  }
               }
            }
        }
            //end rhs.f
            //l2norm.f
        public void l2norm(int ldx, int ldy, int ldz, int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend, double[, , ,] v, double[] sum) {
            //---------------------------------------------------------------------
            //   to compute the l2-norm of vector v.
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //int ldx, ldy, ldz;
            //int nx0, ny0, nz0;
            //int ist, iend;
            //int jst, jend;
            //double  v[5,-1:ldx+2,-1:ldy+2,*], sum[5];

            //teste debug
            if ((isiz1 + 2) != (ldx + 2) || (isiz2 + 2) != (ldy + 2)) {
                throw new ArgumentException("Look this code: vetor v");
            }

            int i, j, k, m;
            double[] dummy = new double[5+1];//dummy[5];

            //int IERROR;

            for(m = 1; m<= 5; m++){
               dummy[m] = 0.0d;
            }
            for(k = 2; k<= nz0-1; k++){
               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                     for(m = 1; m<= 5; m++){
                        dummy[m] = dummy[m] + v[k,j+1,i+1,m] * v[k,j+1,i+1,m];  //dummy[m] = dummy[m] + v[m,i,j,k] * v[m,i,j,k];
                     }
                  }
               }
            }
            //---------------------------------------------------------------------
            //   compute the global sum of individual contributions to dot product.
            //---------------------------------------------------------------------
            worldcomm.Allreduce<double>(dummy, MPI.Operation<double>.Add, ref sum);//call MPI_ALLREDUCE[ dummy,sum,5,dp_type,MPI_SUM,MPI_COMM_WORLD,IERROR ]

            for(m = 1; m<= 5; m++){
               sum[m] = Math.Sqrt(sum[m]/((nx0-2)*(ny0-2)*(nz0-2))); //sum[m] = sqrt(sum[m]/((nx0-2)*(ny0-2)*(nz0-2)));
            }
        }
            //end l2norm.f
            //jacld.f
        public void jacld(int k) {
            //---------------------------------------------------------------------
            //   compute the lower triangular part of the jacobian matrix
            //---------------------------------------------------------------------
            int i, j;
            double  r43;
            double  c1345;
            double  c34;
            double  tmp1, tmp2, tmp3;

            r43 = (4.0d/3.0d);
            c1345 = c1*c3*c4*c5;
            c34 = c3*c4;

               for(j = jst; j<= jend; j++){
                  for(i = ist; i<= iend; i++){
                       //---------------------------------------------------------------------
                       //   form the block daigonal
                       //---------------------------------------------------------------------
                     tmp1 = 1.0d / u[k,j+1,i+1,1];  //tmp1 = 1.0d / u[1, i, j, k];
                     tmp2 = tmp1 * tmp1;
                     tmp3 = tmp1 * tmp2;

                     d[j,i,1,1] =  1.0d+ dt * 2.0d * (tx1 * dx1+ ty1 * dy1+ tz1 * dz1);
                     d[j,i,2,1] =  0.0d;
                     d[j,i,3,1] =  0.0d;
                     d[j,i,4,1] =  0.0d;
                     d[j,i,5,1] =  0.0d;

                     d[j,i,1,2] =  dt*2.0d*(tx1*(-r43*c34*tmp2*u[k,j+1,i+1,2])+ty1*(-c34*tmp2*u[k,j+1,i+1,2])+tz1*(-c34*tmp2*u[k,j+1,i+1,2]));
                     d[j,i,2,2] =  1.0d+dt*2.0d*(tx1*r43*c34*tmp1+ty1*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx2+ty1*dy2+tz1*dz2);
                     d[j,i,3,2] = 0.0d;
                     d[j,i,4,2] = 0.0d;
                     d[j,i,5,2] = 0.0d;

                     d[j,i,1,3] = dt*2.0d*(tx1*(-c34*tmp2*u[k,j+1,i+1,3])+ty1*(-r43*c34*tmp2*u[k,j+1,i+1,3])+tz1*(-c34*tmp2*u[k,j+1,i+1,3]));
                     d[j,i,2,3] = 0.0d;
                     d[j,i,3,3] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*r43*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx3+ty1*dy3+tz1*dz3);
                     d[j,i,4,3] = 0.0d;
                     d[j,i,5,3] = 0.0d;

                     d[j,i,1,4] = dt*2.0d*(tx1*(-c34*tmp2*u[k,j+1,i+1,4])+ty1*(-c34*tmp2*u[k,j+1,i+1,4])+tz1*(-r43*c34*tmp2*u[k,j+1,i+1,4]));
                     d[j,i,2,4] = 0.0d;
                     d[j,i,3,4] = 0.0d;
                     d[j,i,4,4] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*c34*tmp1+tz1*r43*c34*tmp1)+dt*2.0d*(tx1*dx4+ty1*dy4+tz1*dz4);
                     d[j,i,5,4] = 0.0d;

                     d[j,i,1,5] = dt*2.0d*(tx1*(-(r43*c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,2]))-(c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,3]))-(c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,4]))-(c1345)*tmp2*u[k,j+1,i+1,5])+ty1*(-(c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,2]))-(r43*c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,3]))-(c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,4]))-(c1345)*tmp2*u[k,j+1,i+1,5])+tz1*(-(c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,2]))-(c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,3]))-(r43*c34-c1345)*tmp3*(pow2(u[k,j+1,i+1,4]))-(c1345)*tmp2*u[k,j+1,i+1,5]));
                     d[j,i,2,5] = dt*2.0d*(tx1*(r43*c34-c1345)*tmp2*u[k,j+1,i+1,2]+ty1*(c34-c1345)*tmp2*u[k,j+1,i+1,2]+tz1*(c34-c1345)*tmp2*u[k,j+1,i+1,2]);
                     d[j,i,3,5] = dt*2.0d*(tx1*(c34-c1345)*tmp2*u[k,j+1,i+1,3]+ty1*(r43*c34-c1345)*tmp2*u[k,j+1,i+1,3]+tz1*(c34-c1345)*tmp2*u[k,j+1,i+1,3]);
                     d[j,i,4,5] = dt*2.0d*(tx1*(c34-c1345)*tmp2*u[k,j+1,i+1,4]+ty1*(c34-c1345)*tmp2*u[k,j+1,i+1,4]+tz1*(r43*c34-c1345)*tmp2*u[k,j+1,i+1,4]);
                     d[j,i,5,5] = 1.0d+dt*2.0d*(tx1*c1345*tmp1+ty1*c1345*tmp1+tz1*c1345*tmp1)+dt*2.0d*(tx1*dx5+ty1*dy5+tz1*dz5);
                     //---------------------------------------------------------------------
                     //   form the first block sub-diagonal
                     //---------------------------------------------------------------------
                     tmp1 = 1.0d/u[k-1,j+1,i+1,1];
                     tmp2 = tmp1 * tmp1;
                     tmp3 = tmp1 * tmp2;

                     a[j,i,1,1] = - dt * tz1 * dz1;
                     a[j,i,2,1] =   0.0d;
                     a[j,i,3,1] =   0.0d;
                     a[j,i,4,1] = - dt * tz2;
                     a[j,i,5,1] =   0.0d;

                     a[j,i,1,2] = - dt*tz2*(-(u[k-1,j+1,i+1,2]*u[k-1,j+1,i+1,4])*tmp2)-dt*tz1*(-c34*tmp2*u[k-1,j+1,i+1,2]);
                     a[j,i,2,2] = - dt*tz2*(u[k-1,j+1,i+1,4]*tmp1)-dt*tz1*c34*tmp1-dt*tz1*dz2;
                     a[j,i,3,2] = 0.0d;
                     a[j,i,4,2] = -dt*tz2*(u[k-1,j+1,i+1,2]*tmp1);
                     a[j,i,5,2] = 0.0d;

                     a[j,i,1,3] = - dt * tz2* ( - ( u[k-1,j+1,i+1,3]*u[k-1,j+1,i+1,4] ) * tmp2 )- dt * tz1 * ( - c34 * tmp2 * u[k-1,j+1,i+1,3] );
                     a[j,i,2,3] = 0.0d;
                     a[j,i,3,3] = - dt * tz2 * ( u[k-1,j+1,i+1,4] * tmp1 )- dt * tz1 * ( c34 * tmp1 )- dt * tz1 * dz3;
                     a[j,i,4,3] = - dt * tz2 * ( u[k-1,j+1,i+1,3] * tmp1 );
                     a[j,i,5,3] = 0.0d;

                     a[j,i,1,4] = -dt*tz2*(-pow2((u[k-1,j+1,i+1,4]*tmp1))+0.50d*c2*((u[k-1,j+1,i+1,2]*u[k-1,j+1,i+1,2]+u[k-1,j+1,i+1,3]*u[k-1,j+1,i+1,3]+u[k-1,j+1,i+1,4]*u[k-1,j+1,i+1,4])*tmp2))-dt*tz1*(-r43*c34*tmp2*u[k-1,j+1,i+1,4]);
                     a[j,i,2,4] = - dt * tz2* ( - c2 * ( u[k-1,j+1,i+1,2] * tmp1 ) );
                     a[j,i,3,4] = - dt * tz2* ( - c2 * ( u[k-1,j+1,i+1,3] * tmp1 ) );
                     a[j,i,4,4] = - dt * tz2 *(2.0d-c2)*(u[k-1,j+1,i+1,4] * tmp1 )- dt * tz1 * ( r43 * c34 * tmp1 )- dt * tz1 * dz4;
                     a[j,i,5,4] = - dt * tz2 * c2;

                     a[j,i,1,5] = -dt*tz2*((c2*(u[k-1,j+1,i+1,2]*u[k-1,j+1,i+1,2]+ u[k-1,j+1,i+1,3] * u[k-1,j+1,i+1,3]+ u[k-1,j+1,i+1,4] * u[k-1,j+1,i+1,4] ) * tmp2- c1 * ( u[k-1,j+1,i+1,5] * tmp1 ) )* ( u[k-1,j+1,i+1,4] * tmp1 ) )- dt * tz1* ( - ( c34 - c1345 ) * tmp3 * (pow2(u[k-1,j+1,i+1,2]))- ( c34 - c1345 ) * tmp3 * (pow2(u[k-1,j+1,i+1,3]))- ( r43*c34 - c1345)* tmp3 *(pow2(u[k-1,j+1,i+1,4]))- c1345 * tmp2 * u[k-1,j+1,i+1,5]);
                     a[j,i,2,5] = -dt*tz2*(-c2*(u[k-1,j+1,i+1,2]*u[k-1,j+1,i+1,4])*tmp2)-dt*tz1*(c34-c1345)*tmp2*u[k-1,j+1,i+1,2];
                     a[j,i,3,5] = -dt*tz2*(-c2*(u[k-1,j+1,i+1,3]*u[k-1,j+1,i+1,4])*tmp2)-dt*tz1*(c34-c1345)*tmp2*u[k-1,j+1,i+1,3];
                     a[j,i,4,5] = -dt*tz2*(c1*(u[k-1,j+1,i+1,5]*tmp1)-0.50d*c2*((u[k-1,j+1,i+1,2]*u[k-1,j+1,i+1,2]+u[k-1,j+1,i+1,3]*u[k-1,j+1,i+1,3]+3.0d*u[k-1,j+1,i+1,4]*u[k-1,j+1,i+1,4])*tmp2))-dt*tz1*(r43*c34-c1345)*tmp2*u[k-1,j+1,i+1,4];
                     a[j,i,5,5] = -dt*tz2*(c1*(u[k-1,j+1,i+1,4]*tmp1))-dt*tz1*c1345*tmp1-dt*tz1*dz5;
                     //---------------------------------------------------------------------
                     //   form the second block sub-diagonal
                     //---------------------------------------------------------------------
                     tmp1 = 1.0d / u[k,j,i+1,1];
                     tmp2 = tmp1 * tmp1;
                     tmp3 = tmp1 * tmp2;

                     b[j,i,1,1] = - dt * ty1 * dy1;
                     b[j,i,2,1] =   0.0d;
                     b[j,i,3,1] = - dt * ty2;
                     b[j,i,4,1] =   0.0d;
                     b[j,i,5,1] =   0.0d;

                     b[j,i,1,2] = -dt*ty2*(-(u[k,j,i+1,2]*u[k,j,i+1,3])*tmp2)-dt*ty1*(-c34*tmp2*u[k,j,i+1,2]);
                     b[j,i,2,2] = -dt*ty2*(u[k,j,i+1,3]*tmp1)-dt*ty1*(c34*tmp1)-dt*ty1*dy2;
                     b[j,i,3,2] = -dt*ty2*(u[k,j,i+1,2]*tmp1);
                     b[j,i,4,2] = 0.0d;
                     b[j,i,5,2] = 0.0d;

                     b[j,i,1,3] = -dt*ty2*(-pow2((u[k,j,i+1,3]*tmp1))+0.50d*c2*((u[k,j,i+1,2]* u[k,j,i+1,2]+ u[k,j,i+1,3]* u[k,j,i+1,3]+ u[k,j,i+1,4]* u[k,j,i+1,4])* tmp2))- dt * ty1 *(- r43 * c34 * tmp2 * u[k,j,i+1,3]);
                     b[j,i,2,3] = -dt*ty2*(-c2*(u[k,j,i+1,2]*tmp1));
                     b[j,i,3,3] = -dt*ty2*((2.0d-c2)*(u[k,j,i+1,3]*tmp1))-dt*ty1*(r43*c34*tmp1)-dt*ty1*dy3;
                     b[j,i,4,3] = -dt*ty2*(-c2*(u[k,j,i+1,4]*tmp1));
                     b[j,i,5,3] = -dt*ty2*c2;

                     b[j,i,1,4] = -dt*ty2*(-(u[k,j,i+1,3]*u[k,j,i+1,4]) * tmp2 )- dt * ty1 * (- c34 * tmp2 * u[k,j,i+1,4]);                     
                     b[j,i,2,4] = 0.0d;
                     b[j,i,3,4] = -dt*ty2* ( u[k,j,i+1,4] * tmp1 );
                     b[j,i,4,4] = -dt*ty2* ( u[k,j,i+1,3] * tmp1 )- dt * ty1 * ( c34 * tmp1 )- dt * ty1 * dy4;
                     b[j,i,5,4] = 0.0d;

                     b[j,i,1,5] = - dt * ty2* ( ( c2 * (  u[k,j,i+1,2] *u[k,j,i+1,2]+ u[k,j,i+1,3] *u[k,j,i+1,3]+ u[k,j,i+1,4] * u[k,j,i+1,4] ) * tmp2- c1 * ( u[k,j,i+1,5] * tmp1 ) )* ( u[k,j,i+1,3] * tmp1 ) )- dt * ty1* (-(c34 - c1345)*tmp3*(pow2(u[k,j,i+1,2]))- (r43*c34 - c1345)*tmp3*(pow2(u[k,j,i+1,3]))- (c34 - c1345)*tmp3*(pow2(u[k,j,i+1,4]))- c1345*tmp2*u[k,j,i+1,5]);
                     b[j,i,2,5] = - dt * ty2* ( - c2 * (u[k,j,i+1,2]*u[k,j,i+1,3]) * tmp2 )- dt * ty1* (c34 - c1345) * tmp2 * u[k,j,i+1,2];
                     b[j,i,3,5] = - dt * ty2* (c1 *(u[k,j,i+1,5] * tmp1)- 0.50d * c2 *((u[k,j,i+1,2]*u[k,j,i+1,2]+ 3.0d * u[k,j,i+1,3]*u[k,j,i+1,3]+u[k,j,i+1,4]*u[k,j,i+1,4]) * tmp2))- dt * ty1* (r43*c34 - c1345) * tmp2 * u[k,j,i+1,3];
                     b[j,i,4,5] = - dt * ty2* (- c2 *(u[k,j,i+1,3]*u[k,j,i+1,4])* tmp2)- dt * ty1 * (c34 - c1345)* tmp2 *u[k,j,i+1,4];
                     b[j,i,5,5] = - dt * ty2* (c1 *(u[k,j,i+1,3] * tmp1))- dt * ty1 * c1345 * tmp1- dt * ty1 * dy5;
                     //---------------------------------------------------------------------
                     //   form the third block sub-diagonal
                     //---------------------------------------------------------------------                      
                     tmp1 = 1.0d / u[k,j+1,i,1];//tmp1 = 1.0d / u[1,i-1,j,k];
                     tmp2 = tmp1 * tmp1;
                     tmp3 = tmp1 * tmp2;

                     c[j,i,1,1] = - dt * tx1 * dx1;
                     c[j,i,2,1] = - dt * tx2;
                     c[j,i,3,1] =   0.0d;
                     c[j,i,4,1] =   0.0d;
                     c[j,i,5,1] =   0.0d;

                     c[j,i,1,2] = -dt*tx2*(-pow2((u[k,j+1,i,2]*tmp1))+ c2 * 0.50d * ( u[k,j+1,i,2] * u[k,j+1,i,2]+ u[k,j+1,i,3] * u[k,j+1,i,3]+ u[k,j+1,i,4] * u[k,j+1,i,4] ) * tmp2 )- dt * tx1 * ( - r43 * c34 * tmp2 * u[k,j+1,i,2] );
                     c[j,i,2,2] = -dt*tx2*((2.0d-c2)*(u[k,j+1,i,2] * tmp1 ) )- dt * tx1 * ( r43 * c34 * tmp1 )- dt * tx1 * dx2;
                     c[j,i,3,2] = -dt*tx2*(-c2*(u[k,j+1,i,3] * tmp1 ));
                     c[j,i,4,2] = -dt*tx2*(-c2*(u[k,j+1,i,4] * tmp1 ));
                     c[j,i,5,2] = -dt*tx2*c2;

                     c[j,i,1,3] = -dt*tx2*(-( u[k,j+1,i,2] * u[k,j+1,i,3] ) * tmp2 )- dt * tx1 * ( - c34 * tmp2 * u[k,j+1,i,3] );
                     c[j,i,2,3] = -dt*tx2*(   u[k,j+1,i,3] * tmp1 );
                     c[j,i,3,3] = -dt*tx2*(   u[k,j+1,i,2] * tmp1 )- dt * tx1 * ( c34 * tmp1 )- dt * tx1 * dx3;
                     c[j,i,4,3] = 0.0d;
                     c[j,i,5,3] = 0.0d;

                     c[j,i,1,4] = - dt * tx2*(-( u[k,j+1,i,2]*u[k,j+1,i,4] ) * tmp2 )- dt * tx1 * ( - c34 * tmp2 * u[k,j+1,i,4] );
                     c[j,i,2,4] = - dt * tx2 * ( u[k,j+1,i,4] * tmp1 );                    
                     c[j,i,3,4] = 0.0d;
                     c[j,i,4,4] = - dt * tx2 * ( u[k,j+1,i,2] * tmp1 )- dt * tx1 * ( c34 * tmp1 )- dt * tx1 * dx4;
                     c[j,i,5,4] = 0.0d;

                     c[j,i,1,5] = - dt * tx2* ((c2 * (  u[k,j+1,i,2] * u[k,j+1,i,2]+ u[k,j+1,i,3] * u[k,j+1,i,3]+ u[k,j+1,i,4] * u[k,j+1,i,4] ) * tmp2- c1 * ( u[k,j+1,i,5] * tmp1 ))* ( u[k,j+1,i,2] * tmp1 ))- dt * tx1* ( - ( r43*c34 - c1345 ) * tmp3 * (pow2(u[k,j+1,i,2]))- (c34 - c1345) * tmp3 * (pow2(u[k,j+1,i,3]))- (c34 - c1345) * tmp3 * (pow2(u[k,j+1,i,4]))- c1345 * tmp2 * u[k,j+1,i,5]);
                     c[j,i,2,5] = - dt * tx2* (c1 * (u[k,j+1,i,5] * tmp1)- 0.50d * c2* ((  3.0d*u[k,j+1,i,2]*u[k,j+1,i,2]+u[k,j+1,i,3]*u[k,j+1,i,3]+u[k,j+1,i,4]*u[k,j+1,i,4] ) * tmp2 ))- dt * tx1* ( r43*c34 - c1345 ) * tmp2 * u[k,j+1,i,2];
                     c[j,i,3,5] = - dt * tx2* (- c2 * ( u[k,j+1,i,3]*u[k,j+1,i,2] ) * tmp2 )- dt * tx1* (  c34 - c1345 ) * tmp2 * u[k,j+1,i,3];
                     c[j,i,4,5] = - dt * tx2* (- c2 * ( u[k,j+1,i,4]*u[k,j+1,i,2] ) * tmp2 )- dt * tx1* (  c34 - c1345 ) * tmp2 * u[k,j+1,i,4];
                     c[j,i,5,5] = - dt * tx2* (c1 * ( u[k,j+1,i,2] * tmp1 ))- dt * tx1 * c1345 * tmp1- dt * tx1 * dx5;
                  }
               }
        }
            //end jacld.f
            // blts.f
        public void blts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k, double omega, double[, , ,] v, double[, , ,] ldz, double[, , ,] ldy, double[, , ,] ldx, double[, , ,] d, int ist, int iend, int jst, int jend, int nx0, int ny0, int ipt, int jpt) {
            //c---------------------------------------------------------------------
            //c
            //c   compute the regular-sparse, block lower triangular solution:
            //c
            //c                     v <-- [ L-inv ] * v
            //c
            //c---------------------------------------------------------------------
            //c---------------------------------------------------------------------
            //c  input parameters
            //c---------------------------------------------------------------------
            //int ldmx, ldmy, ldmz,nx, ny, nz,k;
            //double  omega;
            //double  v[ 5, -1:ldmx+2, -1:ldmy+2, *],ldz[ 5, 5, ldmx, ldmy],ldy[ 5, 5, ldmx, ldmy],ldx[ 5, 5, ldmx, ldmy],d[ 5, 5, ldmx, ldmy];
            //int ist, iend,jst, jend,nx0, ny0,ipt, jpt;
            //c---------------------------------------------------------------------
            //c  local variables
            //c---------------------------------------------------------------------
            int i, j, m, iex;
            double  tmp, tmp1;
            double[,] tmat = new double[5+1,5+1]; //tmat[5,5]  Obs: nao invertida a ordem
            //---------------------------------------------------------------------
            //   receive data from north and west
            //---------------------------------------------------------------------
            //Check Debug
            if ((isiz1 + 2) != (ldmx + 2) || (isiz2 + 2) != (ldmy + 2)) {
                throw new ArgumentException("Look this code: vetor v");
            }//end Check Debug

            iex = 0;
            exchange_1(v,k,iex);
            for(j = jst; j<= jend; j++){
               for(i = ist; i<= iend; i++){
                  for(m = 1; m<= 5; m++){
                      v[ k, j+1, i+1, m ] =  v[ k, j+1, i+1, m ]
             - omega * (  ldz[j,i,1,m] * v[k-1,j+1,i+1,1 ]
                        + ldz[j,i,2,m] * v[k-1,j+1,i+1,2 ]
                        + ldz[j,i,3,m] * v[k-1,j+1,i+1,3 ]
                        + ldz[j,i,4,m] * v[k-1,j+1,i+1,4 ]
                        + ldz[j,i,5,m] * v[k-1,j+1,i+1,5 ]  );
                  }
               }
            }
            for(j=jst; j<=jend; j++){
              for(i = ist; i<= iend; i++){
                  for(m = 1; m<= 5; m++){
                        v[ k, j+1, i+1, m ] =  v[ k, j+1, i+1, m ]
                        - omega * ( ldy[ j, i,1, m ] * v[k,j  , i+1, 1 ]
                                  + ldx[ j, i,1, m ] * v[k,j+1, i  , 1 ]
                                  + ldy[ j, i,2, m ] * v[k,j  , i+1, 2 ]
                                  + ldx[ j, i,2, m ] * v[k,j+1, i  , 2 ]
                                  + ldy[ j, i,3, m ] * v[k,j  , i+1, 3 ]
                                  + ldx[ j, i,3, m ] * v[k,j+1, i  , 3 ]
                                  + ldy[ j, i,4, m ] * v[k,j  , i+1, 4 ]
                                  + ldx[ j, i,4, m ] * v[k,j+1, i  , 4 ]
                                  + ldy[ j, i,5, m ] * v[k,j  , i+1, 5 ]
                                  + ldx[ j, i,5, m ] * v[k,j+1, i  , 5 ] );
                  }
                  //---------------------------------------------------------------------
                  //   diagonal block inversion
                  //
                  //   forward elimination
                  //---------------------------------------------------------------------
                  for(m = 1; m<= 5; m++){
                     tmat[ m, 1 ] = d[ j, i, 1, m ];
                     tmat[ m, 2 ] = d[ j, i, 2, m ];
                     tmat[ m, 3 ] = d[ j, i, 3, m ];
                     tmat[ m, 4 ] = d[ j, i, 4, m ];
                     tmat[ m, 5 ] = d[ j, i, 5, m ];
                  }
                  tmp1 = 1.0d / tmat[1, 1];
                  tmp = tmp1 * tmat[2, 1];
                  tmat[2, 2] = tmat[2, 2] - tmp * tmat[1, 2];
                  tmat[2, 3] = tmat[2, 3] - tmp * tmat[1, 3];
                  tmat[2, 4] = tmat[2, 4] - tmp * tmat[1, 4];
                  tmat[2, 5] = tmat[2, 5] - tmp * tmat[1, 5];
                  v[k, j+1, i+1, 2] = v[k, j+1, i+1, 2] - v[k, j+1, i+1, 1] * tmp;

                  tmp = tmp1 * tmat[3, 1];
                  tmat[3, 2] = tmat[3, 2] - tmp * tmat[1, 2];
                  tmat[3, 3] = tmat[3, 3] - tmp * tmat[1, 3];
                  tmat[3, 4] = tmat[3, 4] - tmp * tmat[1, 4];
                  tmat[3, 5] = tmat[3, 5] - tmp * tmat[1, 5];
                  v[k, j+1,i+1, 3] = v[k, j+1,i+1, 3] - v[k, j+1,i+1, 1] * tmp;

                  tmp = tmp1 * tmat[4, 1];
                  tmat[4, 2] = tmat[4, 2] - tmp * tmat[1, 2];
                  tmat[4, 3] = tmat[4, 3] - tmp * tmat[1, 3];
                  tmat[4, 4] = tmat[4, 4] - tmp * tmat[1, 4];
                  tmat[4, 5] = tmat[4, 5] - tmp * tmat[1, 5];
                  v[k, j+1,i+1, 4] = v[k, j+1,i+1, 4] - v[k, j+1,i+1, 1] * tmp;

                  tmp = tmp1 * tmat[5, 1];
                  tmat[5, 2] = tmat[5, 2] - tmp * tmat[1, 2];
                  tmat[5, 3] = tmat[5, 3] - tmp * tmat[1, 3];
                  tmat[5, 4] = tmat[5, 4] - tmp * tmat[1, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[1, 5];
                  v[k, j+1,i+1, 5] = v[k, j+1,i+1, 5] - v[k, j+1,i+1, 1] * tmp;


                  tmp1 = 1.0d / tmat[2, 2];
                  tmp = tmp1 * tmat[3, 2];
                  tmat[3, 3] = tmat[3, 3] - tmp * tmat[2, 3];
                  tmat[3, 4] = tmat[3, 4] - tmp * tmat[2, 4];
                  tmat[3, 5] = tmat[3, 5] - tmp * tmat[2, 5];
                  v[k, j+1,i+1, 3] = v[k, j+1,i+1, 3] - v[k, j+1,i+1, 2] * tmp;

                  tmp = tmp1 * tmat[4, 2];
                  tmat[4, 3] = tmat[4, 3] - tmp * tmat[2, 3];
                  tmat[4, 4] = tmat[4, 4] - tmp * tmat[2, 4];
                  tmat[4, 5] = tmat[4, 5] - tmp * tmat[2, 5];
                  v[k, j+1,i+1, 4] = v[k, j+1,i+1, 4] - v[k, j+1,i+1, 2] * tmp;

                  tmp = tmp1 * tmat[5, 2];
                  tmat[5, 3] = tmat[5, 3] - tmp * tmat[2, 3];
                  tmat[5, 4] = tmat[5, 4] - tmp * tmat[2, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[2, 5];
                  v[k, j+1,i+1, 5] = v[k, j+1,i+1, 5] - v[k, j+1,i+1, 2] * tmp;

                  tmp1 = 1.0d / tmat[3, 3];
                  tmp = tmp1 * tmat[4, 3];
                  tmat[4, 4] = tmat[4, 4] - tmp * tmat[3, 4];
                  tmat[4, 5] = tmat[4, 5] - tmp * tmat[3, 5];
                  v[k, j+1,i+1, 4] = v[k, j+1,i+1, 4] - v[k, j+1,i+1, 3] * tmp;

                  tmp = tmp1 * tmat[5, 3];
                  tmat[5, 4] = tmat[5, 4] - tmp * tmat[3, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[3, 5];
                  v[k, j+1,i+1, 5] = v[k, j+1,i+1, 5] - v[k, j+1,i+1, 3] * tmp;

                  tmp1 = 1.0d / tmat[4, 4];
                  tmp = tmp1 * tmat[5, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[4, 5];
                  v[k, j+1,i+1, 5] = v[k, j+1,i+1, 5] - v[k, j+1,i+1, 4] * tmp;
                  //---------------------------------------------------------------------
                  //   back substitution
                  //---------------------------------------------------------------------
                  v[k, j+1,i+1, 5] = v[k, j+1,i+1, 5]/ tmat[ 5, 5 ];
                  v[k, j+1,i+1, 4] = v[k, j+1,i+1, 4]- tmat[ 4, 5 ] * v[ k, j+1, i+1, 5 ];
                  v[k, j+1,i+1, 4] = v[k, j+1,i+1, 4]/ tmat[ 4, 4 ];
                  v[k, j+1,i+1, 3] = v[k, j+1,i+1, 3] - tmat[3, 4] * v[k, j+1, i+1, 4] - tmat[3, 5] * v[k, j+1, i+1, 5];
                  v[k, j+1,i+1, 3] = v[k, j+1,i+1, 3] / tmat[3, 3];
                  v[k, j+1,i+1, 2] = v[k,j+1,i+1,2] - tmat[2, 3] * v[k,j+1,i+1,3]-tmat[2,4]*v[k,j+1,i+1,4]-tmat[2,5]*v[k,j+1,i+1,5];                
                  v[k, j+1,i+1, 2] = v[k, j+1,i+1, 2] / tmat[2, 2];
                  v[k,j+1,i+1,1] = v[k,j+1,i+1,1]-tmat[1,2]*v[k,j+1,i+1,2]-tmat[1,3]*v[k,j+1,i+1,3]-tmat[1,4]*v[k,j+1,i+1,4]-tmat[1,5]*v[k,j+1,i+1,5];
                  v[k, j+1,i+1, 1] = v[k, j+1,i+1, 1] / tmat[1, 1];
              }
            }
            //---------------------------------------------------------------------
            //   send data to east and south
            //---------------------------------------------------------------------
            iex = 2;
            exchange_1(v,k,iex);
        }
            // end blts.f
            // Exchange_1.f
        public void exchange_1(double[, , ,] g, int k, int iex) {
            //double  g[5,-1:isiz1+2,-1:isiz2+2,isiz3];
            //int k,iex;
            int i, j;
            //double dum[5,isiz1+isiz2], dum1[5,isiz1+isiz2];
            //if (id==0) Console.WriteLine("isiz1+isiz2: "+(isiz1+isiz2)+" ||| 5*[jend-jst+1]: "+(5*(jend-jst+1))+" jst: "+jst+" jend: "+jend);

            //int STATUS[MPI_STATUS_SIZE];
            //int IERROR;

            if(iex == 0) {
                if(north != -1) {
                    double[] dum1 = new double[5*(jend-jst+1)];
                    int idx = 0;
                    //call MPI_RECV[ dum1[1,jst],5*[jend-jst+1],dp_type,north,from_n,MPI_COMM_WORLD,status,IERROR ]
                    worldcomm.Receive<double>(north,from_n,ref dum1);
                    for(j=jst; j<=jend; j++){
                        g[k,j+1,1,1] = dum1[0+idx];//g[1,0,j,k] = dum1[1,j];
                        g[k,j+1,1,2] = dum1[1+idx];//g[2,0,j,k] = dum1[2,j];
                        g[k,j+1,1,3] = dum1[2+idx];//g[3,0,j,k] = dum1[3,j];
                        g[k,j+1,1,4] = dum1[3+idx];//g[4,0,j,k] = dum1[4,j];
                        g[k,j+1,1,5] = dum1[4+idx];//g[5,0,j,k] = dum1[5,j];
                        idx = idx + 5;
                    }
                }
                if(west != -1) {
                    double[] dum1 = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    //call MPI_RECV[ dum1[1,ist],5*[iend-ist+1],dp_type,west,from_w,MPI_COMM_WORLD,status,IERROR ]
                    worldcomm.Receive<double>(west, from_w, ref dum1);
                    for(i=ist; i<=iend; i++){
                        g[k,1,i+1,1] = dum1[0+idx];//g[1,i,0,k] = dum1[1,i];
                        g[k,1,i+1,2] = dum1[1+idx];//g[2,i,0,k] = dum1[2,i];
                        g[k,1,i+1,3] = dum1[2+idx];//g[3,i,0,k] = dum1[3,i];
                        g[k,1,i+1,4] = dum1[3+idx];//g[4,i,0,k] = dum1[4,i];
                        g[k,1,i+1,5] = dum1[4+idx];//g[5,i,0,k] = dum1[5,i];
                        idx = idx + 5;
                    }
                }
            } else if(iex == 1){
                if(south != -1){
                    double[] dum1 = new double[(5*(jend-jst+1))];
                    int idx = 0;
                    //call MPI_RECV[ dum1[1,jst],5*[jend-jst+1],dp_type,south,from_s,MPI_COMM_WORLD,status,IERROR ]
                    worldcomm.Receive<double>(south, from_s,ref dum1);
                    for(j=jst; j<=jend; j++){
                        g[k,j+1,nx+2,1] = dum1[0+idx];//g[1,nx+1,j,k] = dum1[1,j];
                        g[k,j+1,nx+2,2] = dum1[1+idx];//g[2,nx+1,j,k] = dum1[2,j];
                        g[k,j+1,nx+2,3] = dum1[2+idx];//g[3,nx+1,j,k] = dum1[3,j];
                        g[k,j+1,nx+2,4] = dum1[3+idx];//g[4,nx+1,j,k] = dum1[4,j];
                        g[k,j+1,nx+2,5] = dum1[4+idx];//g[5,nx+1,j,k] = dum1[5,j];
                        idx = idx + 5;
                    }
                }
                if(east != -1) {
                    double[] dum1 = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    //call MPI_RECV[ dum1[1,ist],5*[iend-ist+1],dp_type,east,from_e,MPI_COMM_WORLD,status,IERROR ]
                    worldcomm.Receive<double>(east, from_e, ref dum1);
                    for(i=ist; i<=iend; i++){
                        g[k,ny+2,i+1,1] = dum1[0+idx];//g[1,i,ny+1,k] = dum1[1,i];
                        g[k,ny+2,i+1,2] = dum1[1+idx];//g[2,i,ny+1,k] = dum1[2,i];
                        g[k,ny+2,i+1,3] = dum1[2+idx];//g[3,i,ny+1,k] = dum1[3,i];
                        g[k,ny+2,i+1,4] = dum1[3+idx];//g[4,i,ny+1,k] = dum1[4,i];
                        g[k,ny+2,i+1,5] = dum1[4+idx];//g[5,i,ny+1,k] = dum1[5,i];
                        idx = idx + 5;
                    }
                }
            } else if(iex == 2) {
                if(south != -1){
                    double[] dum = new double[5*(jend-jst+1)];
                    int idx = 0;
                    for(j=jst; j<=jend; j++){
                        dum[0+idx] = g[k,j+1,nx+1,1];//dum[1,j] = g[1,nx,j,k];
                        dum[1+idx] = g[k,j+1,nx+1,2];//dum[2,j] = g[2,nx,j,k];
                        dum[2+idx] = g[k,j+1,nx+1,3];//dum[3,j] = g[3,nx,j,k];
                        dum[3+idx] = g[k,j+1,nx+1,4];//dum[4,j] = g[4,nx,j,k];
                        dum[4+idx] = g[k,j+1,nx+1,5];//dum[5,j] = g[5,nx,j,k];
                        idx = idx + 5;
                    }
                    //call MPI_SEND[ dum[1,jst], 5*[jend-jst+1], dp_type, south, from_n, MPI_COMM_WORLD, IERROR ]
                    worldcomm.Send<double>(dum, south, from_n);
                }
                if(east != -1){
                    double[] dum = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    for(i=ist; i<=iend; i++){
                        dum[0+idx] = g[k,ny+1,i+1,1];//dum[1,i] = g[1,i,ny,k];
                        dum[1+idx] = g[k,ny+1,i+1,2];//dum[2,i] = g[2,i,ny,k];
                        dum[2+idx] = g[k,ny+1,i+1,3];//dum[3,i] = g[3,i,ny,k];
                        dum[3+idx] = g[k,ny+1,i+1,4];//dum[4,i] = g[4,i,ny,k];
                        dum[4+idx] = g[k,ny+1,i+1,5];//dum[5,i] = g[5,i,ny,k];
                        idx = idx + 5;
                    }
                    //call MPI_SEND[ dum[1,ist], 5*[iend-ist+1], dp_type, east, from_w, MPI_COMM_WORLD, IERROR ]
                    worldcomm.Send<double>(dum, east, from_w);
                }
            } else {
                if(north != -1) {
                    double[] dum = new double[(5*(jend-jst+1))];
                    int idx = 0;
                    for(j=jst; j<=jend; j++){
                        dum[0+idx] = g[k,j+1,2,1];//dum[1,j] = g[1,1,j,k];
                        dum[1+idx] = g[k,j+1,2,2];//dum[2,j] = g[2,1,j,k];
                        dum[2+idx] = g[k,j+1,2,3];//dum[3,j] = g[3,1,j,k];
                        dum[3+idx] = g[k,j+1,2,4];//dum[4,j] = g[4,1,j,k];
                        dum[4+idx] = g[k,j+1,2,5];//dum[5,j] = g[5,1,j,k];
                        idx = idx + 5;
                    }
                    //call MPI_SEND[ dum[1,jst], 5*[jend-jst+1], dp_type, north, from_s, MPI_COMM_WORLD, IERROR ]
                    worldcomm.Send<double>(dum, north, from_s);
                }
                if(west != -1) {
                    double[] dum = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    for(i=ist; i<=iend; i++){
                        dum[0+idx] = g[k,2,i+1,1];//dum[1,i] = g[1,i,1,k];
                        dum[1+idx] = g[k,2,i+1,2];//dum[2,i] = g[2,i,1,k];
                        dum[2+idx] = g[k,2,i+1,3];//dum[3,i] = g[3,i,1,k];
                        dum[3+idx] = g[k,2,i+1,4];//dum[4,i] = g[4,i,1,k];
                        dum[4+idx] = g[k,2,i+1,5];//dum[5,i] = g[5,i,1,k];
                        idx = idx + 5;
                    }
                    //call MPI_SEND[ dum[1,ist], 5*[iend-ist+1], dp_type, west, from_e, MPI_COMM_WORLD, IERROR ]
                    worldcomm.Send<double>(dum, west, from_e);
                }
            }
        }
            //end Exchange_1.f
            // jacu.f
        public void jacu(int k) {
            //---------------------------------------------------------------------
            //   compute the upper triangular part of the jacobian matrix
            //---------------------------------------------------------------------
            int i, j;
            double  r43,c1345,c34,tmp1, tmp2, tmp3;

            r43 = ( 4.0d / 3.0d );
            c1345 = c1 * c3 * c4 * c5;
            c34 = c3 * c4;

            for(j = jst; j<= jend; j++){
                for(i = ist; i<= iend; i++){
                    //---------------------------------------------------------------------
                    //   form the block daigonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k, j+1, i+1, 1];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    d[j,i,1,1] =  1.0d+ dt * 2.0d * (tx1 * dx1+ ty1 * dy1+ tz1 * dz1);
                    d[j,i,2,1] =  0.0d;
                    d[j,i,3,1] =  0.0d;
                    d[j,i,4,1] =  0.0d;
                    d[j,i,5,1] =  0.0d;

                    d[j,i,1,2] =  dt*2.0d*(tx1*(-r43*c34*tmp2*u[k,j+1,i+1,2])+ty1*(-c34*tmp2*u[k,j+1,i+1,2])+tz1*(-c34*tmp2*u[k,j+1,i+1,2]));
                    d[j,i,2,2] =  1.0d+dt*2.0d*(tx1*r43*c34*tmp1+ty1*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx2+ty1*dy2+tz1*dz2);
                    d[j,i,3,2] = 0.0d;
                    d[j,i,4,2] = 0.0d;
                    d[j,i,5,2] = 0.0d;

                    d[j,i,1,3] = dt * 2.0d*(tx1*(-c34*tmp2*u[k,j+1,i+1,3])+ty1*(-r43*c34*tmp2*u[k,j+1,i+1,3])+tz1*(-c34*tmp2 * u[k,j+1,i+1,3]));
                    d[j,i,2,3] = 0.0d;
                    d[j,i,3,3] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*r43*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx3+ty1*dy3+tz1*dz3);
                    d[j,i,4,3] = 0.0d;
                    d[j,i,5,3] = 0.0d;

                    d[j,i,1,4] = dt * 2.0d*(tx1*(-c34*tmp2*u[k,j+1,i+1,4])+ty1*(-c34*tmp2*u[k,j+1,i+1,4])+tz1*(-r43*c34*tmp2*u[k,j+1,i+1,4]));
                    d[j,i,2,4] = 0.0d;
                    d[j,i,3,4] = 0.0d;
                    d[j,i,4,4] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*c34*tmp1+tz1*r43*c34*tmp1)+dt*2.0d*(tx1*dx4+ty1*dy4+tz1*dz4);
                    d[j,i,5,4] = 0.0d;

                    d[j,i,1,5] = dt * 2.0d*(tx1*  (- ( r43*c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,2]))- ( c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,3]))- ( c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,4]))- ( c1345 ) * tmp2 * u[k,j+1,i+1,5])+ ty1 * ( - ( c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,2]))- ( r43*c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,3]))- ( c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,4]))- ( c1345 ) * tmp2 * u[k,j+1,i+1,5])+ tz1 * ( - ( c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,2]))- ( c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,3]))- ( r43*c34 - c1345)*tmp3*(pow2(u[k,j+1,i+1,4]))- ( c1345 ) * tmp2 * u[k,j+1,i+1,5]));
                    d[j,i,2,5] = dt * 2.0d  * ( tx1 * ( r43*c34 - c1345 ) * tmp2 * u[k,j+1,i+1,2]+ ty1 * (     c34 - c1345 ) * tmp2 * u[k,j+1,i+1,2]+ tz1 * (     c34 - c1345 ) * tmp2 * u[k,j+1,i+1,2] );
                    d[j,i,3,5] = dt * 2.0d * ( tx1 * ( c34 - c1345 ) * tmp2 * u[k,j+1,i+1,3]+ ty1 * ( r43*c34 -c1345 ) * tmp2 * u[k,j+1,i+1,3]+ tz1 * ( c34 - c1345 ) * tmp2 * u[k,j+1,i+1,3] ); 
                    d[j,i,4,5] = dt * 2.0d* ( tx1 * ( c34 - c1345 ) * tmp2 * u[k,j+1,i+1,4]+ ty1 * ( c34 - c1345 ) * tmp2 * u[k,j+1,i+1,4]+ tz1 * ( r43*c34 - c1345 ) * tmp2 * u[k,j+1,i+1,4] );
                    d[j,i,5,5] = 1.0d+dt*2.0d*(tx1*c1345*tmp1+ty1*c1345*tmp1+tz1*c1345*tmp1)+dt*2.0d*(tx1*dx5+ty1*dy5+tz1*dz5);
                    //---------------------------------------------------------------------
                    //   form the first block sub-diagonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k,j+1,i+2,1]; //tmp1 = 1.0d / u[1,i+1,j,k];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    a[j,i,1,1] = - dt * tx1 * dx1;
                    a[j,i,2,1] =   dt * tx2;
                    a[j,i,3,1] =   0.0d;
                    a[j,i,4,1] =   0.0d;
                    a[j,i,5,1] =   0.0d;

                    a[j,i,1,2] =  dt * tx2*(-pow2((u[k,j+1,i+2,2]*tmp1))+ c2 * 0.50d * (  u[k,j+1,i+2,2] * u[k,j+1,i+2,2]+ u[k,j+1,i+2,3] * u[k,j+1,i+2,3]+ u[k,j+1,i+2,4] * u[k,j+1,i+2,4] ) * tmp2 )-dt*tx1*(-r43*c34*tmp2*u[k,j+1,i+2,2] );     
                    a[j,i,2,2] =  dt*tx2*((2.0d-c2)*(u[k,j+1,i+2,2]*tmp1))-dt*tx1*(r43*c34*tmp1)-dt*tx1*dx2;
                    a[j,i,3,2] =  dt* tx2*(-c2*(u[k,j+1,i+2,3] * tmp1));
                    a[j,i,4,2] =  dt* tx2*(-c2*(u[k,j+1,i+2,4] * tmp1));
                    a[j,i,5,2] =  dt * tx2 * c2; 

                    a[j,i,1,3] =  dt * tx2*(- (u[k,j+1,i+2,2] * u[k,j+1,i+2,3] ) * tmp2 )- dt * tx1 * (- c34 * tmp2 * u[k,j+1,i+2,3]);
                    a[j,i,2,3] =  dt * tx2 * ( u[k,j+1,i+2,3] * tmp1 );
                    a[j,i,3,3] =  dt * tx2 * ( u[k,j+1,i+2,2] * tmp1 )- dt * tx1 * ( c34 * tmp1 )- dt * tx1 * dx3;
                    a[j,i,4,3] = 0.0d;
                    a[j,i,5,3] = 0.0d;

                    a[j,i,1,4] = dt * tx2*(-( u[k,j+1,i+2,2]*u[k,j+1,i+2,4] ) * tmp2 )- dt * tx1 * ( - c34 * tmp2 * u[k,j+1,i+2,4] );
                    a[j,i,2,4] = dt * tx2 * ( u[k,j+1,i+2,4] * tmp1 );     
                    a[j,i,3,4] = 0.0d;
                    a[j,i,4,4] = dt * tx2 * ( u[k,j+1,i+2,2] * tmp1 )- dt * tx1 * ( c34 * tmp1 )- dt * tx1 * dx4;
                    a[j,i,5,4] = 0.0d;

                    a[j,i,1,5] = dt*tx2*(( c2 * (u[k,j+1,i+2,2] * u[k,j+1,i+2,2]+ u[k,j+1,i+2,3] * u[k,j+1,i+2,3]+ u[k,j+1,i+2,4] * u[k,j+1,i+2,4] ) * tmp2- c1 * ( u[k,j+1,i+2,5] * tmp1 ))* ( u[k,j+1,i+2,2] * tmp1 ))- dt * tx1* ( - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[k,j+1,i+2,2]) )- (     c34 - c1345 ) * tmp3 * ( pow2(u[k,j+1,i+2,3]) )- (c34 - c1345 ) * tmp3 * ( pow2(u[k,j+1,i+2,4]) )- c1345 * tmp2 * u[k,j+1,i+2,5] );
                    a[j,i,2,5] = dt * tx2*(c1 * (u[k,j+1,i+2,5]*tmp1)-0.50d*c2*((3.0d*u[k,j+1,i+2,2]*u[k,j+1,i+2,2]+ u[k,j+1,i+2,3]*u[k,j+1,i+2,3]+ u[k,j+1,i+2,4]*u[k,j+1,i+2,4] ) * tmp2 ))- dt * tx1*(r43*c34-c1345)*tmp2*u[k,j+1,i+2,2];
                    a[j,i,3,5] = dt*tx2*(-c2*(u[k,j+1,i+2,3]*u[k,j+1,i+2,2])*tmp2)-dt*tx1*(c34-c1345)*tmp2* u[k,j+1,i+2,3];
                    a[j,i,4,5] = dt*tx2*(-c2*(u[k,j+1,i+2,4]*u[k,j+1,i+2,2] ) * tmp2 )- dt * tx1* (  c34 - c1345 ) * tmp2 * u[k,j+1,i+2,4];
                    a[j,i,5,5] = dt*tx2*(c1*( u[k,j+1,i+2,2] * tmp1 ))- dt * tx1 * c1345 * tmp1- dt * tx1 * dx5;
                    // c---------------------------------------------------------------------
                    // c   form the second block sub-diagonal
                    // c---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k,j+2,i+1,1]; //tmp1 = 1.0d / u[1,i,j+1,k];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    b[j,i,1,1] = - dt * ty1 * dy1;
                    b[j,i,2,1] =   0.0d;
                    b[j,i,3,1] =  dt * ty2;
                    b[j,i,4,1] =   0.0d;
                    b[j,i,5,1] =   0.0d;

                    b[j,i,1,2] = dt*ty2*(-(u[k,j+2,i+1,2]*u[k,j+2,i+1,3] ) * tmp2 )- dt * ty1 * ( - c34 * tmp2 * u[k,j+2,i+1,2] );
                    b[j,i,2,2] =  dt*ty2*( u[k,j+2,i+1,3] * tmp1 )- dt * ty1 * ( c34 * tmp1 )- dt * ty1 * dy2;
                    b[j,i,3,2] =  dt*ty2*( u[k,j+2,i+1,2] * tmp1 );
                    b[j,i,4,2] = 0.0d;
                    b[j,i,5,2] = 0.0d;

                    b[j,i,1,3] =  dt * ty2* ( -pow2((u[k,j+2,i+1,3]*tmp1))+0.50d*c2*((u[k,j+2,i+1,2] * u[k,j+2,i+1,2]+ u[k,j+2,i+1,3] * u[k,j+2,i+1,3]+ u[k,j+2,i+1,4] * u[k,j+2,i+1,4] )* tmp2 ))-dt*ty1*(-r43*c34*tmp2* u[k,j+2,i+1,3] );
                    b[j,i,2,3] =  dt * ty2*(-c2*( u[k,j+2,i+1,2] * tmp1 )); 
                    b[j,i,3,3] =  dt*ty2*((2.0d-c2)*(u[k,j+2,i+1,3] * tmp1 ))- dt * ty1 * ( r43 * c34 * tmp1 )- dt * ty1 * dy3;
                    b[j,i,4,3] =  dt*ty2*(-c2*(   u[k,j+2,i+1,4] * tmp1 ));
                    b[j,i,5,3] =  dt * ty2 * c2;

                    b[j,i,1,4] =  dt * ty2* ( - ( u[k,j+2,i+1,3]*u[k,j+2,i+1,4] ) * tmp2 )- dt * ty1 * ( - c34 * tmp2 * u[k,j+2,i+1,4] );         
                    b[j,i,2,4] = 0.0d;
                    b[j,i,3,4] =  dt * ty2 *    ( u[k,j+2,i+1,4] * tmp1 );
                    b[j,i,4,4] =  dt * ty2 *    ( u[k,j+2,i+1,3] * tmp1 )- dt * ty1 * ( c34 * tmp1 )- dt * ty1 * dy4;
                    b[j,i,5,4] = 0.0d;

                    b[j,i,1,5] =  dt * ty2* ( ( c2*(u[k,j+2,i+1,2] * u[k,j+2,i+1,2]+ u[k,j+2,i+1,3] * u[k,j+2,i+1,3]+ u[k,j+2,i+1,4] *u[k,j+2,i+1,4] ) * tmp2- c1 * ( u[k,j+2,i+1,5]*tmp1))*(u[k,j+2,i+1,3] * tmp1 ))- dt * ty1* ( - (     c34 - c1345 )*tmp3*   (pow2(u[k,j+2,i+1,2]))- ( r43*c34 - c1345 )*tmp3*   (pow2(u[k,j+2,i+1,3]))- (     c34 - c1345 )*tmp3*   (pow2(u[k,j+2,i+1,4]))- c1345*tmp2*u[k,j+2,i+1,5] );
                    b[j,i,2,5] =  dt * ty2* ( - c2*(u[k,j+2,i+1,2]*u[k,j+2,i+1,3] ) * tmp2 )- dt * ty1* ( c34 - c1345 ) * tmp2 * u[k,j+2,i+1,2];
                    b[j,i,3,5] =  dt * ty2* ( c1 * ( u[k,j+2,i+1,5] * tmp1 )- 0.50d * c2 * ((u[k,j+2,i+1,2]*u[k,j+2,i+1,2]+ 3.0d * u[k,j+2,i+1,3]*u[k,j+2,i+1,3]+ u[k,j+2,i+1,4]*u[k,j+2,i+1,4] ) * tmp2 ))- dt * ty1* ( r43*c34 - c1345 ) * tmp2 * u[k,j+2,i+1,3];
                    b[j,i,4,5] =  dt * ty2* ( - c2*( u[k,j+2,i+1,3]*u[k,j+2,i+1,4] ) * tmp2 )- dt * ty1 * ( c34 - c1345 ) * tmp2 * u[k,j+2,i+1,4];
                    b[j,i,5,5] =  dt * ty2* ( c1 * ( u[k,j+2,i+1,3] * tmp1 ) )- dt * ty1 * c1345 * tmp1- dt * ty1 * dy5;
                    // c---------------------------------------------------------------------
                    // c   form the third block sub-diagonal
                    // c---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k+1,j+1,i+1,1];   //tmp1 = 1.0d / u[1,i,j,k+1];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    c[j,i,1,1] = - dt * tz1 * dz1;
                    c[j,i,2,1] =   0.0d;
                    c[j,i,3,1] =   0.0d;
                    c[j,i,4,1] = dt * tz2;
                    c[j,i,5,1] =   0.0d;

                    c[j,i,1,2] = dt * tz2* ( - ( u[k+1,j+1,i+1,2]*u[k+1,j+1,i+1,4] ) * tmp2 )- dt * tz1 * ( - c34 * tmp2 * u[k+1,j+1,i+1,2] );
                    c[j,i,2,2] = dt * tz2 *    ( u[k+1,j+1,i+1,4] * tmp1 )- dt * tz1 * c34 * tmp1- dt * tz1 * dz2;       
                    c[j,i,3,2] = 0.0d;
                    c[j,i,4,2] = dt * tz2 *    ( u[k+1,j+1,i+1,2] * tmp1 );
                    c[j,i,5,2] = 0.0d;

                    c[j,i,1,3] = dt * tz2* ( - ( u[k+1,j+1,i+1,3]*u[k+1,j+1,i+1,4] ) * tmp2 )- dt * tz1 * ( - c34 * tmp2 * u[k+1,j+1,i+1,3] );                    
                    c[j,i,2,3] = 0.0d;
                    c[j,i,3,3] = dt * tz2 *    ( u[k+1,j+1,i+1,4] * tmp1 )- dt * tz1 * ( c34 * tmp1 )- dt * tz1 * dz3;
                    c[j,i,4,3] = dt * tz2 *    ( u[k+1,j+1,i+1,3] * tmp1 );
                    c[j,i,5,3] = 0.0d;

                    c[j,i,1,4] = dt * tz2     * ( - pow2((u[k+1,j+1,i+1,4]*tmp1)) + 0.50d * c2 * ((u[k+1,j+1,i+1,2] * u[k+1,j+1,i+1,2]+ u[k+1,j+1,i+1,3] * u[k+1,j+1,i+1,3]+ u[k+1,j+1,i+1,4] * u[k+1,j+1,i+1,4] ) * tmp2 ))- dt * tz1 * ( - r43 * c34 * tmp2 * u[k+1,j+1,i+1,4] );
                    c[j,i,2,4] = dt * tz2* ( - c2 * ( u[k+1,j+1,i+1,2] * tmp1 ) );
                    c[j,i,3,4] = dt * tz2* ( - c2 * ( u[k+1,j+1,i+1,3] * tmp1 ) );
                    c[j,i,4,4] = dt * tz2*(2.0d-c2)*( u[k+1,j+1,i+1,4] * tmp1 )- dt * tz1 * ( r43 * c34 * tmp1 )- dt * tz1 * dz4;
                    c[j,i,5,4] = dt * tz2 * c2;

                    c[j,i,1,5] = dt * tz2*((c2 * (u[k+1,j+1,i+1,2] * u[k+1,j+1,i+1,2]+ u[k+1,j+1,i+1,3] * u[k+1,j+1,i+1,3]+ u[k+1,j+1,i+1,4] * u[k+1,j+1,i+1,4] ) * tmp2- c1 * ( u[k+1,j+1,i+1,5] * tmp1 ))* ( u[k+1,j+1,i+1,4] * tmp1 )) - dt * tz1*(-(c34-c1345)*tmp3*(pow2(u[k+1,j+1,i+1,2]))-(c34-c1345)*tmp3 * (pow2(u[k+1,j+1,i+1,3]))-(r43*c34-c1345)*tmp3*(pow2(u[k+1,j+1,i+1,4]))- c1345 * tmp2 * u[k+1,j+1,i+1,5] );
                    c[j,i,2,5] = dt * tz2* (-c2*( u[k+1,j+1,i+1,2]*u[k+1,j+1,i+1,4] ) * tmp2 )- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[k+1,j+1,i+1,2];
                    c[j,i,3,5] = dt * tz2*(-c2* ( u[k+1,j+1,i+1,3]*u[k+1,j+1,i+1,4] ) * tmp2 )- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[k+1,j+1,i+1,3];
                    c[j,i,4,5] = dt * tz2* (c1* ( u[k+1,j+1,i+1,5] * tmp1 ) - 0.50d * c2* ((  u[k+1,j+1,i+1,2]*u[k+1,j+1,i+1,2]+ u[k+1,j+1,i+1,3]*u[k+1,j+1,i+1,3]+ 3.0d*u[k+1,j+1,i+1,4]*u[k+1,j+1,i+1,4] ) * tmp2 ))- dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u[k+1,j+1,i+1,4];
                    c[j,i,5,5] = dt * tz2* (c1* ( u[k+1,j+1,i+1,4] * tmp1 ))- dt * tz1 * c1345 * tmp1- dt * tz1 * dz5;
               }
            }
        }
            // end jacu.f
            // buts.f
        public void buts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k, double omega,double[, , ,] v,double[, ,] tv,double[, , ,] d,double[, , ,] udx,double[, , ,] udy,double[, , ,] udz, int ist, int iend, int jst, int jend, int nx0, int ny0, int ipt, int jpt) {
            //c---------------------------------------------------------------------
            //c
            //c   compute the regular-sparse, block upper triangular solution:
            //c
            //c                     v <-- [ U-inv ] * v
            //c
            //c---------------------------------------------------------------------
            //c---------------------------------------------------------------------
            //c  input parameters
            //c---------------------------------------------------------------------
            //int ldmx, ldmy, ldmz, nx, ny, nz, k; double  omega;
            //double  v[ 5, -1:ldmx+2, -1:ldmy+2, *], tv[5, ldmx, ldmy], 
            //        d[ 5, 5, ldmx, ldmy], udx[ 5, 5, ldmx, ldmy], udy[ 5, 5, ldmx, ldmy], udz[ 5, 5, ldmx, ldmy ];
            //int ist, iend,jst, jend,nx0, ny0,ipt, jpt;
            //c---------------------------------------------------------------------
            //c  local variables
            //c---------------------------------------------------------------------
            int i, j, m, iex;
            double tmp, tmp1;
            double[,] tmat = new double[5 + 1, 5 + 1];//tmat[5,5]
            //---------------------------------------------------------------------
            //   receive data from south and east
            //---------------------------------------------------------------------
            iex = 1;
            exchange_1(v,k,iex);
            //Debug teste
            if ((isiz1 + 2) != (ldmx + 2) || (isiz2 + 2) != (ldmy + 2)) {
                throw new ArgumentException("Look this code: vetor v");
            } // end debug teste
            
            for(j = jend; j>= jst; j--){ //for(j = jend, jst, -1;
               for(i = iend; i>= ist; i--){ //for(i = iend, ist, -1;
                  for(m = 1; m<= 5; m++){//tv[ m, i, j ] = 
                        tv[ j, i, m ] = 
                                         omega * (  udz[j ,i, 1, m ] * v[k+1,j+1,i+1,1]
                                                  + udz[j ,i, 2, m ] * v[k+1,j+1,i+1,2]
                                                  + udz[j ,i, 3, m ] * v[k+1,j+1,i+1,3]
                                                  + udz[j ,i, 4, m ] * v[k+1,j+1,i+1,4]
                                                  + udz[j ,i, 5, m ] * v[k+1,j+1,i+1,5] );
                  }
               }
            }
            for(j = jend; j>=jst; j--){   //for(j = jend,jst,-1;
              for(i = iend; i>=ist; i--){ //for(i = iend,ist,-1;
                  for(m = 1; m<= 5; m++){
                        tv[j,i,m] = tv[j,i,m]
                                            + omega * ( udy[j,i,1,m] * v[k, j+2, i+1,  1]
                                                      + udx[j,i,1,m] * v[k, j+1, i+2,  1]
                                                      + udy[j,i,2,m] * v[k, j+2, i+1,  2]
                                                      + udx[j,i,2,m] * v[k, j+1, i+2,  2]
                                                      + udy[j,i,3,m] * v[k, j+2, i+1,  3]
                                                      + udx[j,i,3,m] * v[k, j+1, i+2,  3]
                                                      + udy[j,i,4,m] * v[k, j+2, i+1,  4]
                                                      + udx[j,i,4,m] * v[k, j+1, i+2,  4]
                                                      + udy[j,i,5,m] * v[k, j+2, i+1,  5]
                                                      + udx[j,i,5,m] * v[k, j+1, i+2,  5] );
                  }
                  //---------------------------------------------------------------------
                  //   diagonal block inversion
                  //---------------------------------------------------------------------
                  for(m = 1; m<= 5; m++){
                     tmat[ m, 1 ] = d[j, i,1, m];
                     tmat[ m, 2 ] = d[j, i,2, m];
                     tmat[ m, 3 ] = d[j, i,3, m];
                     tmat[ m, 4 ] = d[j, i,4, m];
                     tmat[ m, 5 ] = d[j, i,5, m];
                  }
                  tmp1 = 1.0d / tmat[1, 1];
                  tmp = tmp1 * tmat[ 2, 1 ];
                  tmat[ 2, 2 ] =  tmat[ 2, 2 ]- tmp * tmat[ 1, 2 ];
                  tmat[ 2, 3 ] =  tmat[ 2, 3 ]- tmp * tmat[ 1, 3 ];
                  tmat[ 2, 4 ] =  tmat[ 2, 4 ]- tmp * tmat[ 1, 4 ];
                  tmat[ 2, 5 ] =  tmat[ 2, 5 ]- tmp * tmat[ 1, 5 ];
                  tv[j,i,2] = tv[j,i,2]- tv[j,i,1] * tmp;
                  
                  tmp = tmp1 * tmat[3, 1];
                  tmat[3, 2] = tmat[3, 2] - tmp * tmat[1, 2];
                  tmat[3, 3] = tmat[3, 3] - tmp * tmat[1, 3];
                  tmat[3, 4] = tmat[3, 4] - tmp * tmat[1, 4];
                  tmat[3, 5] = tmat[3, 5] - tmp * tmat[1, 5];
                  tv[j, i, 3] = tv[j, i, 3] - tv[j, i, 1] * tmp;
                  
                  tmp = tmp1 * tmat[4, 1];
                  tmat[4, 2] = tmat[4, 2] - tmp * tmat[1, 2];
                  tmat[4, 3] = tmat[4, 3] - tmp * tmat[1, 3];
                  tmat[4, 4] = tmat[4, 4] - tmp * tmat[1, 4];
                  tmat[4, 5] = tmat[4, 5] - tmp * tmat[1, 5];
                  tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 1] * tmp;
                  
                  tmp = tmp1 * tmat[5, 1];
                  tmat[5, 2] = tmat[5, 2] - tmp * tmat[1, 2];
                  tmat[5, 3] = tmat[5, 3] - tmp * tmat[1, 3];
                  tmat[5, 4] = tmat[5, 4] - tmp * tmat[1, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[1, 5];
                  tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 1] * tmp;
                  
                  tmp1 = 1.0d / tmat[2, 2];
                  tmp = tmp1 * tmat[3, 2];
                  tmat[3, 3] = tmat[3, 3] - tmp * tmat[2, 3];
                  tmat[3, 4] = tmat[3, 4] - tmp * tmat[2, 4];
                  tmat[3, 5] = tmat[3, 5] - tmp * tmat[2, 5];
                  tv[j, i, 3] = tv[j, i, 3] - tv[j, i, 2] * tmp;
                  
                  tmp = tmp1 * tmat[4, 2];
                  tmat[4, 3] = tmat[4, 3] - tmp * tmat[2, 3];
                  tmat[4, 4] = tmat[4, 4] - tmp * tmat[2, 4];
                  tmat[4, 5] = tmat[4, 5] - tmp * tmat[2, 5];
                  tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 2] * tmp;
                  
                  tmp = tmp1 * tmat[5, 2];
                  tmat[5, 3] = tmat[5, 3] - tmp * tmat[2, 3];
                  tmat[5, 4] = tmat[5, 4] - tmp * tmat[2, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[2, 5];
                  tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 2] * tmp;

                  tmp1 = 1.0d / tmat[3, 3];
                  tmp = tmp1 * tmat[4, 3];
                  tmat[4, 4] = tmat[4, 4] - tmp * tmat[3, 4];
                  tmat[4, 5] = tmat[4, 5] - tmp * tmat[3, 5];
                  tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 3] * tmp;

                  tmp = tmp1 * tmat[5, 3];
                  tmat[5, 4] = tmat[5, 4] - tmp * tmat[3, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[3, 5];
                  tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 3] * tmp;

                  tmp1 = 1.0d / tmat[4, 4];
                  tmp = tmp1 * tmat[5, 4];
                  tmat[5, 5] = tmat[5, 5] - tmp * tmat[4, 5];
                  tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 4] * tmp;
                  //c---------------------------------------------------------------------
                  //c   back substitution
                  //c---------------------------------------------------------------------
                  tv[j, i, 5] = tv[j, i, 5]/ tmat[ 5, 5 ];
                  tv[j, i, 4] = tv[j, i, 4]- tmat[ 4, 5 ] * tv[j, i, 5];
                  tv[j, i, 4] = tv[j, i, 4]/ tmat[ 4, 4 ];
                  tv[j, i, 3] = tv[j, i, 3]- tmat[ 3, 4 ] * tv[j, i, 4]- tmat[ 3, 5 ] * tv[j, i, 5];
                  tv[j, i, 3] = tv[j, i, 3]/ tmat[ 3, 3 ];
                  tv[j, i, 2] = tv[j, i, 2]- tmat[ 2, 3 ] * tv[j, i, 3]- tmat[ 2, 4 ] * tv[j, i, 4]- tmat[ 2, 5 ] * tv[j, i, 5];
                  tv[j, i, 2] = tv[j, i, 2]/ tmat[ 2, 2 ];
                  tv[j, i, 1] = tv[j, i, 1]-tmat[1, 2]*tv[j,i,2]-tmat[1,3]*tv[j,i,3]-tmat[1,4]*tv[j,i,4]-tmat[1,5]*tv[j,i,5];
                  tv[j, i, 1] = tv[j, i, 1] / tmat[1, 1];
                  v[k,j+1,i+1,1] = v[k,j+1,i+1,1] - tv[j, i, 1];
                  v[k,j+1,i+1,2] = v[k,j+1,i+1,2] - tv[j, i, 2];
                  v[k,j+1,i+1,3] = v[k,j+1,i+1,3] - tv[j, i, 3];
                  v[k,j+1,i+1,4] = v[k,j+1,i+1,4] - tv[j, i, 4];
                  v[k,j+1,i+1,5] = v[k,j+1,i+1,5] - tv[j, i, 5];
              }
            }
            //---------------------------------------------------------------------
            //   send data to north and west
            //---------------------------------------------------------------------
            iex = 3;
            exchange_1(v,k,iex);
        }
            // end buts.f

        //end ssor.f
    }
}
