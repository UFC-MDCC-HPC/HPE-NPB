/*
-------------------------------------------------------------------------
        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         
                                   L U                                   
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
 Authors: S. Weeratunga
          V. Venkatakrishnan
          E. Barszcz
          M. Yarrow
 Translation to C# and MPI.NET Code
          Cenez Araújo de Rezende, MDCC/UFC
          Francisco Heron de Carvalho Junior (MDCC/UFC)
-------------------------------------------------------------------------
*/
using System;
using System.IO;
using NPB3_0_JAV.BMInOut;
using MPI;

namespace NPB {
    public class LU: LUBase {

        public LU(char c) : base(c) {
        }

        static void Main(String[] argv) {

            LU lu = null;
            LUBase.debug = false;

            try {
                string param = argv[0];
            }
            catch(Exception) {
                argv = new String[1];
                argv[0] = "CLASS=S"; 
            }
            char paramClass;
            if(!LUBase.debug) {
                IO.parseCmdLineArgs(argv);
                paramClass = IO.CLASS;
            }
            else {
                paramClass = 'K';  
            }                      

            try {
                lu = new LU(paramClass);
            }
            catch(OutOfMemoryException e) {
                Console.WriteLine(e.ToString());
                System.Environment.Exit(0);
            }
            lu.runBenchMark();
        }

        public void runBenchMark() {
            read_input();
            //---------------------------------------------------------------------
            //   set up processor grid
            //---------------------------------------------------------------------
            proc_grid();
            //---------------------------------------------------------------------
            //   determine the neighbors
            //---------------------------------------------------------------------
            neighbors();
            //---------------------------------------------------------------------
            //   set up sub-domain sizes
            //---------------------------------------------------------------------
            subdomain();
            //---------------------------------------------------------------------
            //   set up coefficients
            //---------------------------------------------------------------------
            setConstants();
            //---------------------------------------------------------------------
            //   set the masks required for comm
            //---------------------------------------------------------------------
            sethyper();
            //---------------------------------------------------------------------
            //   set the boundary values for dependent variables
            //---------------------------------------------------------------------
            setbv();
            //---------------------------------------------------------------------
            //   set the initial values for dependent variables
            //---------------------------------------------------------------------
            setiv();
            //---------------------------------------------------------------------
            //   compute the forcing term based on prescribed exact solution
            //---------------------------------------------------------------------
            erhs();
            //---------------------------------------------------------------------
            //   perform one SSOR iteration to touch all data and program pages 
            //---------------------------------------------------------------------
            ssor(1);
            //---------------------------------------------------------------------
            //   reset the boundary and initial values
            //---------------------------------------------------------------------
            setbv();
            setiv();
            //---------------------------------------------------------------------
            //   perform the SSOR iterations
            //---------------------------------------------------------------------
            ssor(itmax);
            //---------------------------------------------------------------------
            //   compute the solution error
            //---------------------------------------------------------------------
            error();
            //---------------------------------------------------------------------
            //   compute the surface integral
            //---------------------------------------------------------------------
            pintgr();
            //---------------------------------------------------------------------
            //   verification test
            //---------------------------------------------------------------------
            if(node==0) {
                mflops = ((double)(itmax))*(1984.77*((double)(nx0))*((double)(ny0))*((double)(nz0))-10923.3*pow2((((double)(nx0+ny0+nz0))/3.0))+27770.9*((double)(nx0+ny0+nz0))/3.0-144010.0) / (maxtime*1000000.0);
                //IO.print_results(BMName, clss, nx0, ny0, nz0, itmax, nnodes_compiled, num, maxtime, mflops, "floating point", verified, npbversion);//compiletime, cs1, cs2, cs3, cs4, cs5, cs6, '[none]');
                int verified = verify(rsdnm, errnm, frc);
                BMResults results = new BMResults(BMName,
                                        clss,
                                        nx0,
                                        ny0,
                                        nz0,
                                        itmax,
                                        maxtime,
                                        mflops,
                                        "floating point",
                                        verified,
                                        true,
                                        num,
                                        -1);
                results.print();            
            }
            mpi.Dispose();//call mpi_finalize[ierr];
        }

        public void read_input() {
            int fstatus=0, nnodes;
            //---------------------------------------------------------------------
            //    only root reads the input file
            //    if input file does not exist, it uses defaults
            //       ipr = 1 for detailed progress output
            //       inorm = how often the norm is printed [once every inorm iterations]
            //       itmax = number of pseufor(time steps
            //       dt = time step
            //       omega 1 over-relaxation factor for SSOR
            //       tolrsd = steady state residual tolerance levels
            //       nx, ny, nz = number of grid points in x, y, z directions
            //---------------------------------------------------------------------
            if(node == root) {
                string[] vetTemp = new string[13];
                try {
                    Console.Write("Trying Reading from input file inputlu.data: ");
                    int[] conf = { 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 5, 0, 0, 3 };
                    vetTemp = IO.readFileData("inputlu.data", conf);//open [unit=3,file='inputlu.data',status='old', access='sequential',form='formatted', iostat=fstatus];
                }
                catch(System.IO.FileNotFoundException) {
                    Console.WriteLine("inputlu.data not found");
                    fstatus = 1;
                }
                Console.WriteLine(" NAS Parallel Benchmarks "+npbversion+" -- LU Benchmark ");
                if(fstatus == 0) {
                    Console.WriteLine("Reading from input file inputlu.data");
                    ipr       = int.Parse(vetTemp[0]);//read [3,*] ipr, inorm
                    inorm     = int.Parse(vetTemp[1]);
                    itmax     = int.Parse(vetTemp[2]);//read [3,*] itmax
                    dt        = double.Parse(vetTemp[3]);//read [3,*] dt
                    omega     = double.Parse(vetTemp[4]);//read [3,*] omega
                    tolrsd[0] = double.Parse(vetTemp[5]);
                    tolrsd[1] = double.Parse(vetTemp[6]);
                    tolrsd[2] = double.Parse(vetTemp[7]);
                    tolrsd[3] = double.Parse(vetTemp[8]);
                    tolrsd[4] = double.Parse(vetTemp[9]);//read [3,*] tolrsd[1],tolrsd[2],tolrsd[3],tolrsd[4],tolrsd[5]
                    nx0       = int.Parse(vetTemp[10]);
                    ny0       = int.Parse(vetTemp[11]);
                    nz0       = int.Parse(vetTemp[12]);//read [3,*] nx0, ny0, nz0
                }
                else {
                    ipr = ipr_default;
                    inorm = inorm_default;
                    itmax = itmax_default;
                    dt = dt_default;
                    omega = omega_default;
                    tolrsd[0] = tolrsd1_def;
                    tolrsd[1] = tolrsd2_def;
                    tolrsd[2] = tolrsd3_def;
                    tolrsd[3] = tolrsd4_def;
                    tolrsd[4] = tolrsd5_def;
                    nx0 = isiz01;
                    ny0 = isiz02;
                    nz0 = isiz03;
                }
                nnodes = num;//   call MPI_COMM_SIZE[MPI_COMM_WORLD, nnodes, ierror];
                //---------------------------------------------------------------------
                //   check problem size
                //---------------------------------------------------------------------
                if(nnodes != nnodes_compiled) {
                    Console.WriteLine("Warning: program is running on"+nnodes+" processors, but was compiled for "+nnodes_compiled);
                }
                if((nx0 < 4) || (ny0 < 4) || (nz0 < 4)) {
                    Console.WriteLine("PROBLEM SIZE IS TOO SMALL - "+" SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5");
                    worldcomm.Abort(0);
                    mpi.Dispose();//CALL MPI_ABORT[ MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR ];
                    System.Environment.Exit(0);
                }

                if((nx0 > isiz01) || (ny0 > isiz02) || (nz0 > isiz03)) {
                    Console.WriteLine("PROBLEM SIZE IS TOO LARGE - NX, NY AND NZ SHOULD BE LESS THAN OR EQUAL TO ISIZ01, ISIZ02 AND ISIZ03 RESPECTIVELY");
                    worldcomm.Abort(0);
                    mpi.Dispose();//      CALL MPI_ABORT[ MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR ];
                    System.Environment.Exit(0);
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
            worldcomm.Broadcast<double>(ref tolrsd, root);    //call MPI_BCAST[tolrsd, 5, dp_type, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref nx0, root);          //call MPI_BCAST[nx0, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref ny0, root);          //call MPI_BCAST[ny0, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
            worldcomm.Broadcast<int>(ref nz0, root);          //call MPI_BCAST[nz0, 1, MPI_int, root, MPI_COMM_WORLD, ierr]
        }

        public void proc_grid() {
            //---------------------------------------------------------------------
            //
            //   set up a two-d grid for processors: column-major ordering of unknowns
            //   NOTE: assumes a power-of-two number of processors
            //
            //---------------------------------------------------------------------
            xdim   = (int)Math.Pow(2, (ndim/2));//xdim   = 2**(ndim/2);
            if(mod(ndim, 2)==1)
                xdim = xdim + xdim;
            ydim   = num/xdim;
            row    = (int)mod(node, xdim) + 1;
            col    = node/xdim + 1;
        }

        public void neighbors() {
            //---------------------------------------------------------------------
            //     figure out the neighbors and their wrap numbers for each processor
            //---------------------------------------------------------------------
            south = -1;
            east  = -1;
            north = -1;
            west  = -1;
            if(row>1) {
                north = node -1;
            }
            else {
                north = -1;
            }
            if(row < xdim) {
                south = node + 1;
            }
            else {
                south = -1;
            }

            if(col > 1) {
                west = node - xdim;
            }
            else {
                west = -1;
            }
            if(col < ydim) {
                east = node + xdim;
            }
            else {
                east = -1;
            }
        }

        public void subdomain() {
            int mm;
            //---------------------------------------------------------------------
            //   set up the sub-domain sizes
            //---------------------------------------------------------------------
            //   x dimension
            //---------------------------------------------------------------------
            mm   = (int)mod(nx0, xdim);
            if(row<=mm) {
                nx = nx0/xdim + 1;
                ipt = (row-1)*nx;
            }
            else {
                nx = nx0/xdim;
                ipt = (row-1)*nx + mm;
            }
            //---------------------------------------------------------------------
            //   y dimension
            //---------------------------------------------------------------------
            mm   = (int)mod(ny0, ydim);
            if(col<=mm) {
                ny = ny0/ydim + 1;
                jpt = (col-1)*ny;
            }
            else {
                ny = ny0/ydim;
                jpt = (col-1)*ny + mm;
            }
            //---------------------------------------------------------------------
            //   z dimension
            //---------------------------------------------------------------------
            nz = nz0;
            //---------------------------------------------------------------------
            //   check the sub-domain size
            //---------------------------------------------------------------------
            if((nx < 4) || (ny < 4) || (nz < 4)) {
                Console.WriteLine("SUBDOMAIN SIZE IS TOO SMALL - ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS, "+
                    "SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL TO 4 THEY ARE CURRENTLY: "+nx+"x"+ny+"x"+nz);
                worldcomm.Abort(0);//CALL MPI_ABORT[ MPI_COMM_WORLD,ERRORCODE,IERROR ]
                mpi.Dispose();
                System.Environment.Exit(0);
            }
            if((nx > isiz1) || (ny > isiz2) || (nz > isiz3)) {
                Console.WriteLine("SUBDOMAIN SIZE IS TOO LARGE - ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS" +
                    "SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY. THEY ARE CURRENTLY"+
                    " "+nx+"x"+ny+"x"+nz);
                worldcomm.Abort(0);//CALL MPI_ABORT[ MPI_COMM_WORLD,ERRORCODE, IERROR ]
                mpi.Dispose();
                System.Environment.Exit(0);
            }
            //---------------------------------------------------------------------
            //   set up the start and end in i and j extents for all processors
            //---------------------------------------------------------------------
            ist = 1;
            iend = nx;
            if(north==-1)
                ist = 2;
            if(south==-1)
                iend = nx - 1;
            jst = 1;
            jend = ny;
            if(west==-1)
                jst = 2;
            if(east==-1)
                jend = ny - 1;
        }

        public void setConstants() { // setcoeff()
            //---------------------------------------------------------------------
            //   set up coefficients
            //---------------------------------------------------------------------
            dxi   = 1.0d / (nx0 - 1);
            deta  = 1.0d / (ny0 - 1);
            dzeta = 1.0d / (nz0 - 1);

            tx1 = 1.0d/(dxi*dxi);
            tx2 = 1.0d/(2.0d*dxi);
            tx3 = 1.0d/dxi;

            ty1 = 1.0d/ (deta * deta);
            ty2 = 1.0d/ (2.0d* deta);
            ty3 = 1.0d/ deta;

            tz1 = 1.0d/ (dzeta * dzeta);
            tz2 = 1.0d/ (2.0d* dzeta);
            tz3 = 1.0d/ dzeta;

            ii1 = 2;
            ii2 = nx0 - 1;
            ji1 = 2;
            ji2 = ny0 - 2;
            ki1 = 3;
            ki2 = nz0 - 1;

            //---------------------------------------------------------------------
            //   diffusion coefficients
            //---------------------------------------------------------------------
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
            //---------------------------------------------------------------------
            //   fourth difference dissipation
            //---------------------------------------------------------------------      
            dssp = (max(max(dx1, dy1), dz1))/4.0d; //dssp=(max(dx1, dy1, dz1))/4.0d
         }

        public void sethyper() {
            //---------------------------------------------------------------------
            //    for each column in a hyperplane, istart = first row,
            //---------------------------------------------------------------------
            int i, j, iglob, jglob, kp;
            //---------------------------------------------------------------------
            // compute the pointers for hyperplanes
            //---------------------------------------------------------------------
            for(kp = 1; kp<(nx0+ny0); kp++) {
                icomms[kp] = false;
                icommn[kp] = false;
                icomme[kp] = false;
                icommw[kp] = false;
                //---------------------------------------------------------------------
                //  check to see if comm. to south is required
                //---------------------------------------------------------------------
                if(south!=-1) {
                    i     = iend;
                    iglob = ipt + i;
                    jglob = kp - iglob + 1;
                    j     = jglob - jpt;
                    if(jglob>=2 && jglob<=ny0-1 && j>=jst && j<=jend)
                        icomms[kp] = true;
                }
                //---------------------------------------------------------------------
                //  check to see if comm. to north is required
                //---------------------------------------------------------------------
                if(north!=-1) {
                    i     = ist;
                    iglob = ipt + i;
                    jglob = kp - iglob + 1;
                    j     = jglob - jpt;
                    if(jglob>=2 && jglob<=ny0-1 && j>=jst && j<=jend)
                        icommn[kp] = true;
                }
                //---------------------------------------------------------------------
                //  check to see if comm. to east is required
                //---------------------------------------------------------------------
                if(east!=-1) {
                    j     = jend;
                    jglob = jpt + j;
                    iglob = kp - jglob + 1;
                    i     = iglob - ipt;
                    if(iglob>=2 && iglob<=nx0-1 && i>=ist && i<=iend)
                        icomme[kp] = true;
                }
                //---------------------------------------------------------------------
                //  check to see if comm. to west is required
                //---------------------------------------------------------------------
                if(west!=-1) {
                    j = jst;
                    jglob = jpt + j;
                    iglob = kp - jglob + 1;
                    i     = iglob - ipt;
                    if(iglob>=2 && iglob<=nx0-1 && i>=ist && i<=iend)
                        icommw[kp] = true;
                }
            }
            icomms[0] = false;
            icommn[0] = false;
            icomme[0] = false;
            icommw[0] = false;
            icomms[nx0+ny0] = false;
            icommn[nx0+ny0] = false;
            icomme[nx0+ny0] = false;
            icommw[nx0+ny0] = false;
        }

        public void setbv() {
            //---------------------------------------------------------------------
            //   set the boundary values of dependent variables
            //---------------------------------------------------------------------
            int i, j, k, iglob, jglob;
            //---------------------------------------------------------------------
            //   set the dependent variable values along the top and bottom faces
            //---------------------------------------------------------------------
            for(j = 1; j<= ny; j++) {
                jglob = jpt + j;
                for(i = 1; i<= nx; i++) {
                    iglob = ipt + i;
                    exact(iglob, jglob, 1, u, 0, j+1, i+1);   //exact( iglob, jglob, 1, u[ 1, i, j, 1 ] );
                    exact(iglob, jglob, nz, u, nz-1, j+1, i+1);   //exact( iglob, jglob, nz, u[ 1, i, j, nz ] );
                }
            }
            //---------------------------------------------------------------------
            //   set the dependent variable values along north and south faces
            //---------------------------------------------------------------------
            if(west==-1) {
                for(k = 1; k<= nz; k++) {
                    for(i = 1; i<= nx; i++) {
                        iglob = ipt + i;
                        exact(iglob, 1, k, u, k-1, 1+1, i+1);   //call exact[ iglob, 1, k, u[ 1, i, 1, k ] ];
                    }
                }
            }
            if(east==-1) {
                for(k = 1; k<= nz; k++) {
                    for(i = 1; i<= nx; i++) {
                        iglob = ipt + i;
                        exact(iglob, ny0, k, u, k-1, ny+1, i+1); //call exact[ iglob, ny0, k, u[ 1, i, ny, k ] ];
                    }
                }
            }
            //---------------------------------------------------------------------
            //   set the dependent variable values along east and west faces
            //---------------------------------------------------------------------
            if(north==-1) {
                for(k = 1; k<= nz; k++) {
                    for(j = 1; j<= ny; j++) {
                        jglob = jpt + j;
                        exact(1, jglob, k, u, k-1, j+1, 1+1);  //call exact[ 1, jglob, k, u[ 1, 1, j, k ] ];
                    }
                }
            }
            if(south==-1) {
                for(k = 1; k<= nz; k++) {
                    for(j = 1; j<= ny; j++) {
                        jglob = jpt + j;
                        exact(nx0, jglob, k, u, k-1, j+1, nx+1); // call exact[ nx0, jglob, k, u[ 1, nx, j, k ] ];
                    }
                }
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
            double[,,,] ue_1jk   = new double[1, 1, 1, 5];   //ue_1jk[5]
            double[,,,] ue_nx0jk = new double[1, 1, 1, 5];   //ue_nx0jk[5]
            double[,,,] ue_i1k   = new double[1, 1, 1, 5];   //ue_i1k[5]
            double[,,,] ue_iny0k = new double[1, 1, 1, 5];   //ue_iny0k[5]
            double[,,,] ue_ij1   = new double[1, 1, 1, 5];   //ue_ij1[5]
            double[,,,] ue_ijnz  = new double[1, 1, 1, 5];   //ue_ijnz[5]
            for(k = 2; k<=nz-1; k++) {
                zeta = ((double)(k-1))/(nz-1);
                for(j = 1; j<= ny; j++) {
                    jglob = jpt + j;
                    if(jglob!=1 && jglob!=ny0) {
                        eta = ((double)(jglob-1))/(ny0-1);
                        for(i = 1; i<= nx; i++) {
                            iglob = ipt + i;
                            if(iglob!=1 && iglob!=nx0) {
                                xi = ((double)(iglob-1))/(nx0-1);
                                exact(1, jglob, k, ue_1jk, 0, 0, 0);
                                exact(nx0, jglob, k, ue_nx0jk, 0, 0, 0);
                                exact(iglob, 1, k, ue_i1k, 0, 0, 0);
                                exact(iglob, ny0, k, ue_iny0k, 0, 0, 0);
                                exact(iglob, jglob, 1, ue_ij1, 0, 0, 0);
                                exact(iglob, jglob, nz, ue_ijnz, 0, 0, 0);
                                for(m = 0; m< 5; m++) {
                                    pxi =   (1.0d-xi) * ue_1jk[0, 0, 0, m] + xi   * ue_nx0jk[0, 0, 0, m];
                                    peta =  (1.0d-eta) * ue_i1k[0, 0, 0, m] + eta   * ue_iny0k[0, 0, 0, m];
                                    pzeta = (1.0d-zeta) * ue_ij1[0, 0, 0, m] + zeta   * ue_ijnz[0, 0, 0, m];
                                    u[k-1, j+1, i+1, m] = pxi + peta + pzeta - pxi * peta - peta * pzeta - pzeta * pxi + pxi * peta * pzeta;
                                }
                            }
                        }
                    }
                }
            }
        }//ue_1jk    ue_nx0jk      ue_i1k      ue_iny0k     ue_ij1     ue_ijnz m

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
            for(k = 0; k< nz; k++) {
                for(j = 2; j< ny+2; j++) {
                    for(i = 2; i< nx+2; i++) {
                        for(m = 0; m< 5; m++) {
                            frct[k, j, i, m] = 0.0d; //frct[ m, i, j, k ] = 0.0d;
                        }
                    }
                }
            }
            for(k = 0; k< nz; k++) {
                zeta = ((double)(k))/(nz-1);
                for(j = 2; j< ny+2; j++) {
                    jglob = jpt + j - 2;
                    eta = ((double)(jglob))/(ny0-1);
                    for(i = 2; i< nx+2; i++) {
                        iglob = ipt + i - 2;
                        xi = ((double)(iglob))/(nx0-1);
                        for(m = 0; m< 5; m++) {  //rsd[m,i,j,k] =  ce[m,1]
                            rsd[k, j, i, m] =  ce[0, m]
                            + ce[1, m] * xi
                            + ce[2, m] * eta
                            + ce[3, m] * zeta
                            + ce[4, m] * xi * xi
                            + ce[5, m] * eta * eta
                            + ce[6, m] * zeta * zeta
                            + ce[7, m] * xi * xi * xi
                            + ce[8, m] * eta * eta * eta
                            + ce[9, m] * zeta * zeta * zeta
                            + ce[10, m] * xi * xi * xi * xi
                            + ce[11, m] * eta * eta * eta * eta
                            + ce[12, m] * zeta * zeta * zeta * zeta;
                        }
                    }
                }
            }//k i j jglob iglob
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
            exchange_3(rsd, iex);
            L1 = 0;
            if(north==-1)
                L1 = 1;
            L2 = nx + 1;
            if(south==-1)
                L2 = nx;

            ist1 = 1;
            iend1 = nx;
            if(north==-1)
                ist1 = 4;
            if(south==-1)
                iend1 = nx - 3;
            for(k = 2; k<= nz - 1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = L1; i<= L2; i++) {
                        flux[k-1, j, i, 0] = rsd[k-1, j+1, i+1, 1]; //flux[1,i,j,k] = rsd[2,i,j,k];
                        u21           = rsd[k-1, j+1, i+1, 1]/rsd[k-1, j+1, i+1, 0]; //u21 = rsd[2,i,j,k] / rsd[1,i,j,k];
                        //c -- q = 0.50d*(rsd[2,i,j,k]*rsd[2,i,j,k] + rsd[3,i,j,k]*rsd[3,i,j,k] + rsd[4,i,j,k]*rsd[4,i,j,k])/rsd[1,i,j,k];
                        q=0.50d*(rsd[k-1, j+1, i+1, 1]*rsd[k-1, j+1, i+1, 1]+rsd[k-1, j+1, i+1, 2]*rsd[k-1, j+1, i+1, 2]+rsd[k-1, j+1, i+1, 3]*rsd[k-1, j+1, i+1, 3])/rsd[k-1, j+1, i+1, 0];
                        flux[k-1, j, i, 1] =     rsd[k-1, j+1, i+1, 1]*u21 + c2*(rsd[k-1, j+1, i+1, 4] - q);//flux[2,i,j,k]=rsd[2,i,j,k]*u21+c2*(rsd[5,i,j,k]-q);
                        flux[k-1, j, i, 2] =     rsd[k-1, j+1, i+1, 2] * u21;                          //flux[3,i,j,k]=rsd[3,i,j,k] * u21;
                        flux[k-1, j, i, 3] =     rsd[k-1, j+1, i+1, 3] * u21;                          //flux[4,i,j,k]=rsd[4,i,j,k] * u21;
                        flux[k-1, j, i, 4] = (c1*rsd[k-1, j+1, i+1, 4] - c2*q)*u21;                    //flux[5,i,j,k]=(c1*rsd[5,i,j,k] - c2*q)*u21;
                    }
                }
            }
            for(k = 2; k<= nz - 1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) { //frct[m,i,j,k] =  frct[m,i,j,k] - tx2 * (flux[m,i+1,j,k] - flux[m,i-1,j,k]);
                            frct[k-1, j+1, i+1, m] =  frct[k-1, j+1, i+1, m] - tx2 * (flux[k-1, j, i+1, m] - flux[k-1, j, i-1, m]);
                        }
                    }
                    for(i = ist; i<= L2; i++) {
                        tmp   = 1.0d/rsd[k-1, j+1, i+1, 0];
                        u21i = tmp * rsd[k-1, j+1, i+1, 1];
                        u31i = tmp * rsd[k-1, j+1, i+1, 2];
                        u41i = tmp * rsd[k-1, j+1, i+1, 3];
                        u51i = tmp * rsd[k-1, j+1, i+1, 4];
                        tmp   = 1.0d/rsd[k-1, j+1, i, 0];
                        u21im1 = tmp*rsd[k-1, j+1, i, 1];
                        u31im1 = tmp*rsd[k-1, j+1, i, 2];
                        u41im1 = tmp*rsd[k-1, j+1, i, 3];
                        u51im1 = tmp*rsd[k-1, j+1, i, 4];

                        flux[k-1, j, i, 1] = (4.0d/3.0d)*tx3*(u21i - u21im1);
                        flux[k-1, j, i, 2] = tx3 * (u31i - u31im1);
                        flux[k-1, j, i, 3] = tx3 * (u41i - u41im1);
                        flux[k-1, j, i, 4] = 0.50d*(1.0d-c1*c5)*tx3*
                         ((pow2(u21i)+pow2(u31i)+pow2(u41i))
                         -(pow2(u21im1)+pow2(u31im1)+pow2(u41im1)))
                         + (1.0d/6.0d)*tx3*(pow2(u21i) - pow2(u21im1))+c1*c5*tx3*(u51i-u51im1);
                    }
                    for(i = ist; i<= iend; i++) {
                        frct[k-1, j+1, i+1, 0] = frct[k-1, j+1, i+1, 0]+dx1*tx1*(rsd[k-1, j+1, i, 0]-2.0d*rsd[k-1, j+1, i+1, 0]+rsd[k-1, j+1, i+2, 0]);
                        frct[k-1, j+1, i+1, 1] = frct[k-1, j+1, i+1, 1]+tx3*c3*c4*(flux[k-1, j, i+1, 1]-flux[k-1, j, i, 1])+dx2*tx1*(rsd[k-1, j+1, i, 1]-2.0d*rsd[k-1, j+1, i+1, 1]+rsd[k-1, j+1, i+2, 1]);
                        frct[k-1, j+1, i+1, 2] = frct[k-1, j+1, i+1, 2]+tx3*c3*c4*(flux[k-1, j, i+1, 2]-flux[k-1, j, i, 2])+dx3*tx1*(rsd[k-1, j+1, i, 2]-2.0d*rsd[k-1, j+1, i+1, 2]+rsd[k-1, j+1, i+2, 2]);
                        frct[k-1, j+1, i+1, 3] = frct[k-1, j+1, i+1, 3]+tx3*c3*c4*(flux[k-1, j, i+1, 3]-flux[k-1, j, i, 3])+dx4*tx1*(rsd[k-1, j+1, i, 3]-2.0d*rsd[k-1, j+1, i+1, 3]+rsd[k-1, j+1, i+2, 3]);
                        frct[k-1, j+1, i+1, 4] = frct[k-1, j+1, i+1, 4]+tx3*c3*c4*(flux[k-1, j, i+1, 4]-flux[k-1, j, i, 4])+dx5*tx1*(rsd[k-1, j+1, i, 4]-2.0d*rsd[k-1, j+1, i+1, 4]+rsd[k-1, j+1, i+2, 4]);
                    }
                    //---------------------------------------------------------------------
                    //   Fourth-order dissipation
                    //---------------------------------------------------------------------
                    if(north==-1) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, j+1, 3, m] = frct[k-1, j+1, 3, m]-dsspm*(+5.0d*rsd[k-1, j+1, 3, m]-4.0d*rsd[k-1, j+1, 4, m]+rsd[k-1, j+1, 5, m]);
                            frct[k-1, j+1, 4, m] = frct[k-1, j+1, 4, m]-dsspm*(-4.0d*rsd[k-1, j+1, 3, m]+6.0d*rsd[k-1, j+1, 4, m]-4.0d*rsd[k-1, j+1, 5, m]+rsd[k-1, j+1, 6, m]);
                        }
                    }
                    for(i = ist1; i<=iend1; i++) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, j+1, i+1, m] = frct[k-1, j+1, i+1, m]-dsspm*(rsd[k-1, j+1, i-1, m]-
                            4.0d*rsd[k-1, j+1, i, m]+6.0d*rsd[k-1, j+1, i+1, m]-4.0d*rsd[k-1, j+1, i+2, m]+rsd[k-1, j+1, i+3, m]);
                        }
                    }
                    if(south==-1) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, j+1, nx-1, m] = frct[k-1, j+1, nx-1, m]-dsspm*(rsd[k-1, j+1, nx-3, m]-4.0d*rsd[k-1, j+1, nx-2, m]+6.0d*rsd[k-1, j+1, nx-1, m]-4.0d*rsd[k-1, j+1, nx, m]);
                            frct[k-1, j+1, nx, m]   = frct[k-1, j+1, nx, m]  -dsspm*(rsd[k-1, j+1, nx-2, m]-4.0d*rsd[k-1, j+1, nx-1, m]+5.0d*rsd[k-1, j+1, nx, m]);
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
            exchange_3(rsd, iex);
            L1 = 0;
            if(west==-1)
                L1 = 1;
            L2 = ny + 1;
            if(east==-1)
                L2 = ny;
            jst1 = 1;
            jend1 = ny;
            if(west==-1)
                jst1 = 4;
            if(east==-1)
                jend1 = ny - 3;
            for(k = 2; k<= nz - 1; k++) {
                for(j = L1; j<= L2; j++) {
                    for(i = ist; i<= iend; i++) {
                        flux[k-1, j, i, 0] = rsd[k-1, j+1, i+1, 2];
                        u31 = rsd[k-1, j+1, i+1, 2] / rsd[k-1, j+1, i+1, 0];
                        q          = 0.50d*(rsd[k-1, j+1, i+1, 1]*rsd[k-1, j+1, i+1, 1]+rsd[k-1, j+1, i+1, 2]*rsd[k-1, j+1, i+1, 2]
                                        +rsd[k-1, j+1, i+1, 3]*rsd[k-1, j+1, i+1, 3])/rsd[k-1, j+1, i+1, 0];
                        flux[k-1, j, i, 1] =     rsd[k-1, j+1, i+1, 1]*u31;
                        flux[k-1, j, i, 2] =     rsd[k-1, j+1, i+1, 2]*u31+c2*(rsd[k-1, j+1, i+1, 4]-q);
                        flux[k-1, j, i, 3] =     rsd[k-1, j+1, i+1, 3]*u31;
                        flux[k-1, j, i, 4] = (c1*rsd[k-1, j+1, i+1, 4]-c2*q)*u31;
                    }
                }
            }
            for(k = 2; k<= nz - 1; k++) {
                for(i = ist; i<= iend; i++) {
                    for(j = jst; j<= jend; j++) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, j+1, i+1, m] =  frct[k-1, j+1, i+1, m] - ty2 * (flux[k-1, j+1, i, m] - flux[k-1, j-1, i, m]);
                        }
                    }
                }
                for(j = jst; j<= L2; j++) {
                    for(i = ist; i<= iend; i++) {
                        tmp = 1.0d / rsd[k-1, j+1, i+1, 0];
                        u21j = tmp * rsd[k-1, j+1, i+1, 1];
                        u31j = tmp * rsd[k-1, j+1, i+1, 2];
                        u41j = tmp * rsd[k-1, j+1, i+1, 3];
                        u51j = tmp * rsd[k-1, j+1, i+1, 4];
                        tmp = 1.0d / rsd[k-1, j, i+1, 0];
                        u21jm1 = tmp*rsd[k-1, j, i+1, 1];
                        u31jm1 = tmp*rsd[k-1, j, i+1, 2];
                        u41jm1 = tmp*rsd[k-1, j, i+1, 3];
                        u51jm1 = tmp*rsd[k-1, j, i+1, 4];
                        flux[k-1, j, i, 1] = ty3*(u21j-u21jm1);
                        flux[k-1, j, i, 2] = (4.0d/3.0d)*ty3*(u31j-u31jm1);
                        flux[k-1, j, i, 3] = ty3*(u41j-u41jm1);
                        flux[k-1, j, i, 4] = 0.50d*(1.0d-c1*c5)*ty3*((pow2(u21j)+pow2(u31j)+pow2(u41j))-(pow2(u21jm1)+pow2(u31jm1)+pow2(u41jm1)))+(1.0d/6.0d)*ty3*(pow2(u31j)-pow2(u31jm1))+c1*c5*ty3*(u51j-u51jm1);
                    }
                }
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        frct[k-1, j+1, i+1, 0] = frct[k-1, j+1, i+1, 0]+dy1*ty1*(rsd[k-1, j, i+1, 0]-2.0d*rsd[k-1, j+1, i+1, 0]+rsd[k-1, j+2, i+1, 0]);
                        frct[k-1, j+1, i+1, 1] = frct[k-1, j+1, i+1, 1]+ty3*c3*c4*(flux[k-1, j+1, i, 1]-flux[k-1, j, i, 1])+dy2*ty1*(rsd[k-1, j, i+1, 1]-2.0d*rsd[k-1, j+1, i+1, 1]+rsd[k-1, j+2, i+1, 1]);
                        frct[k-1, j+1, i+1, 2] = frct[k-1, j+1, i+1, 2]+ty3*c3*c4*(flux[k-1, j+1, i, 2]-flux[k-1, j, i, 2])+dy3*ty1*(rsd[k-1, j, i+1, 2]-2.0d*rsd[k-1, j+1, i+1, 2]+rsd[k-1, j+2, i+1, 2]);
                        frct[k-1, j+1, i+1, 3] = frct[k-1, j+1, i+1, 3]+ty3*c3*c4*(flux[k-1, j+1, i, 3]-flux[k-1, j, i, 3])+dy4*ty1*(rsd[k-1, j, i+1, 3]-2.0d*rsd[k-1, j+1, i+1, 3]+rsd[k-1, j+2, i+1, 3]);
                        frct[k-1, j+1, i+1, 4] = frct[k-1, j+1, i+1, 4]+ty3*c3*c4*(flux[k-1, j+1, i, 4]-flux[k-1, j, i, 4])+dy5*ty1*(rsd[k-1, j, i+1, 4]-2.0d*rsd[k-1, j+1, i+1, 4]+rsd[k-1, j+2, i+1, 4]);
                    }
                }
                //---------------------------------------------------------------------
                //   fourth-order dissipation
                //---------------------------------------------------------------------
                if(west==-1) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, 3, i+1, m] = frct[k-1, 3, i+1, m]-dsspm*(+5.0d*rsd[k-1, 3, i+1, m]-4.0d*rsd[k-1, 4, i+1, m]+rsd[k-1, 5, i+1, m]);
                            frct[k-1, 4, i+1, m] = frct[k-1, 4, i+1, m]-dsspm*(-4.0d*rsd[k-1, 3, i+1, m]+6.0d*rsd[k-1, 4, i+1, m]-4.0d*rsd[k-1, 5, i+1, m]+rsd[k-1, 6, i+1, m]);
                        }
                    }
                }
                for(j = jst1; j<= jend1; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, j+1, i+1, m]=frct[k-1, j+1, i+1, m]-dsspm*(rsd[k-1, j-1, i+1, m]-4.0d*rsd[k-1, j, i+1, m]+6.0d*rsd[k-1, j+1, i+1, m]-4.0d*rsd[k-1, j+2, i+1, m]+rsd[k-1, j+3, i+1, m]);
                        }
                    }
                }
                if(east==-1) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, ny-1, i+1, m] = frct[k-1, ny-1, i+1, m]-dsspm*(rsd[k-1, ny-3, i+1, m]-4.0d*rsd[k-1, ny-2, i+1, m]+6.0d*rsd[k-1, ny-1, i+1, m]-4.0d*rsd[k-1, ny, i+1, m]);
                            frct[k-1, ny, i+1, m] = frct[k-1, ny, i+1, m]-dsspm*(rsd[k-1, ny-2, i+1, m]-4.0d*rsd[k-1, ny-1, i+1, m]+5.0d*rsd[k-1, ny, i+1, m]);
                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //   zeta-direction flux differences
            //---------------------------------------------------------------------
            for(k = 1; k<= nz; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        flux[k-1, j, i, 0] = rsd[k-1, j+1, i+1, 3];      //flux[1,i,j,k] = rsd[4,i,j,k];
                        u41 = rsd[k-1, j+1, i+1, 3] / rsd[k-1, j+1, i+1, 0]; //u41 = rsd[4,i,j,k] / rsd[1,i,j,k];
                        q = 0.50d*(rsd[k-1, j+1, i+1, 1]*rsd[k-1, j+1, i+1, 1]+rsd[k-1, j+1, i+1, 2]*rsd[k-1, j+1, i+1, 2]+rsd[k-1, j+1, i+1, 3]*rsd[k-1, j+1, i+1, 3])/rsd[k-1, j+1, i+1, 0];
                        flux[k-1, j, i, 1] =rsd[k-1, j+1, i+1, 1] * u41;
                        flux[k-1, j, i, 2] =rsd[k-1, j+1, i+1, 2] * u41;
                        flux[k-1, j, i, 3] =rsd[k-1, j+1, i+1, 3] * u41 + c2*(rsd[k-1, j+1, i+1, 4] - q);
                        flux[k-1, j, i, 4] =(c1*rsd[k-1, j+1, i+1, 4]-c2*q)*u41;
                    }
                }
            }
            for(k = 2; k<= nz - 1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, j+1, i+1, m] =  frct[k-1, j+1, i+1, m] - tz2 * (flux[k, j, i, m] - flux[k-2, j, i, m]);
                        }
                    }
                }
            }
            for(k = 2; k<= nz; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        tmp = 1.0d / rsd[k-1, j+1, i+1, 0];
                        u21k = tmp * rsd[k-1, j+1, i+1, 1];
                        u31k = tmp * rsd[k-1, j+1, i+1, 2];
                        u41k = tmp * rsd[k-1, j+1, i+1, 3];
                        u51k = tmp * rsd[k-1, j+1, i+1, 4];

                        tmp = 1.0d / rsd[k-2, j+1, i+1, 0];
                        u21km1 = tmp*rsd[k-2, j+1, i+1, 1];
                        u31km1 = tmp*rsd[k-2, j+1, i+1, 2];
                        u41km1 = tmp*rsd[k-2, j+1, i+1, 3];
                        u51km1 = tmp*rsd[k-2, j+1, i+1, 4];

                        flux[k-1, j, i, 1] = tz3 * (u21k - u21km1);
                        flux[k-1, j, i, 2] = tz3 * (u31k - u31km1);
                        flux[k-1, j, i, 3] = (4.0d/3.0d) * tz3 * (u41k - u41km1);
                        flux[k-1, j, i, 4] = 0.50d*(1.0d-c1*c5)*tz3*((pow2(u21k)+pow2(u31k)+pow2(u41k))-(pow2(u21km1)+pow2(u31km1)+pow2(u41km1)))+(1.0d/6.0d)*tz3*(pow2(u41k)-pow2(u41km1))+c1*c5*tz3*(u51k-u51km1);
                    }
                }
            }
            for(k = 2; k<= nz - 1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        frct[k-1, j+1, i+1, 0] = frct[k-1, j+1, i+1, 0]+dz1*tz1*(rsd[k, j+1, i+1, 0]-2.0d*rsd[k-1, j+1, i+1, 0]+rsd[k-2, j+1, i+1, 0]);
                        frct[k-1, j+1, i+1, 1] = frct[k-1, j+1, i+1, 1]+tz3*c3*c4*(flux[k, j, i, 1]-flux[k-1, j, i, 1])+dz2*tz1*(rsd[k, j+1, i+1, 1]-2.0d*rsd[k-1, j+1, i+1, 1]+rsd[k-2, j+1, i+1, 1]);
                        frct[k-1, j+1, i+1, 2] = frct[k-1, j+1, i+1, 2]+tz3*c3*c4*(flux[k, j, i, 2]-flux[k-1, j, i, 2])+dz3*tz1*(rsd[k, j+1, i+1, 2]-2.0d*rsd[k-1, j+1, i+1, 2]+rsd[k-2, j+1, i+1, 2]);
                        frct[k-1, j+1, i+1, 3] = frct[k-1, j+1, i+1, 3]+tz3*c3*c4*(flux[k, j, i, 3]-flux[k-1, j, i, 3])+dz4*tz1*(rsd[k, j+1, i+1, 3]-2.0d*rsd[k-1, j+1, i+1, 3]+rsd[k-2, j+1, i+1, 3]);
                        frct[k-1, j+1, i+1, 4] = frct[k-1, j+1, i+1, 4]+tz3*c3*c4*(flux[k, j, i, 4]-flux[k-1, j, i, 4])+dz5*tz1*(rsd[k, j+1, i+1, 4]-2.0d*rsd[k-1, j+1, i+1, 4]+rsd[k-2, j+1, i+1, 4]);
                    }
                }
            }
            //---------------------------------------------------------------------
            //   fourth-order dissipation
            //---------------------------------------------------------------------
            for(j = jst; j<= jend; j++) {
                for(i = ist; i<= iend; i++) {
                    for(m = 0; m< 5; m++) {
                        frct[-1+2, j+1, i+1, m] = frct[-1+2, j+1, i+1, m]-dsspm*(+5.0d*rsd[-1+2, j+1, i+1, m]-4.0d*rsd[-1+3, j+1, i+1, m]+rsd[-1+4, j+1, i+1, m]);
                        frct[-1+3, j+1, i+1, m] = frct[-1+3, j+1, i+1, m]-dsspm*(-4.0d*rsd[-1+2, j+1, i+1, m]+6.0d*rsd[-1+3, j+1, i+1, m]-4.0d*rsd[-1+4, j+1, i+1, m]+rsd[-1+5, j+1, i+1, m]);
                    }
                }
            }
            for(k = 4; k<= nz - 3; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            frct[k-1, j+1, i+1, m]=frct[k-1, j+1, i+1, m]-dsspm*(rsd[k-3, j+1, i+1, m]-4.0d*rsd[k-2, j+1, i+1, m]+6.0d*rsd[k-1, j+1, i+1, m]-4.0d*rsd[k, j+1, i+1, m]+rsd[k+1, j+1, i+1, m]);
                        }
                    }
                }
            }
            for(j = jst; j<= jend; j++) {
                for(i = ist; i<= iend; i++) {
                    for(m = 0; m< 5; m++) {
                        frct[-1+nz-2, j+1, i+1, m]=frct[-1+nz-2, j+1, i+1, m]-dsspm*(rsd[-1+nz-4, j+1, i+1, m]- 4.0d*rsd[-1+nz-3, j+1, i+1, m]+6.0d*rsd[-1+nz-2, j+1, i+1, m]-4.0d*rsd[-1+nz-1, j+1, i+1, m]);
                        frct[-1+nz-1, j+1, i+1, m]=frct[-1+nz-1, j+1, i+1, m]-dsspm*(rsd[-1+nz-3, j+1, i+1, m]- 4.0d*rsd[-1+nz-2, j+1, i+1, m]+5.0d*rsd[-1+nz-1, j+1, i+1, m]);
                    }
                }
            }
        }
        //Exchange_3.f
        public void exchange_3(double[, , ,] g, int iex) {
            //---------------------------------------------------------------------
            //   compute the right hand side based on exact solution
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
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

            if(iex==0) {
                //---------------------------------------------------------------------
                //   communicate in the south and north directions
                //---------------------------------------------------------------------
                if(north!=-1) {
                    //call MPI_IRECV[ buf1, 10*ny*nz, dp_type, MPI_ANY_SOURCE, from_n, MPI_COMM_WORLD, mid, IERROR ];
                    mid[0] = worldcomm.ImmediateReceive<double>(north, from_n, buf1);
                }
                //---------------------------------------------------------------------
                //   send south
                //---------------------------------------------------------------------
                if(south!=-1) {
                    for(k = 1; k<=nz; k++) {
                        for(j = 1; j<=ny; j++) {
                            ipos1 = (k-1)*ny+j              -1;  //ipos1 = (k-1)*ny+j;
                            ipos2 = ipos1 + ny*nz;               //ipos2 = ipos1 + ny*nz;
                            buf[0*size2+ipos1] = g[k-1, j+1, nx, 0];  //buf[1,ipos1] = g[1,nx-1,j,k];
                            buf[1*size2+ipos1] = g[k-1, j+1, nx, 1];  //buf[2,ipos1] = g[2,nx-1,j,k];
                            buf[2*size2+ipos1] = g[k-1, j+1, nx, 2];  //buf[3,ipos1] = g[3,nx-1,j,k];
                            buf[3*size2+ipos1] = g[k-1, j+1, nx, 3];  //buf[4,ipos1] = g[4,nx-1,j,k]; 
                            buf[4*size2+ipos1] = g[k-1, j+1, nx, 4];  //buf[5,ipos1] = g[5,nx-1,j,k];

                            buf[0*size2+ipos2] = g[k-1, j+1, nx+1, 0];    //buf[1,ipos2] = g[1,nx,j,k];
                            buf[1*size2+ipos2] = g[k-1, j+1, nx+1, 1];    //buf[2,ipos2] = g[2,nx,j,k];
                            buf[2*size2+ipos2] = g[k-1, j+1, nx+1, 2];    //buf[3,ipos2] = g[3,nx,j,k];
                            buf[3*size2+ipos2] = g[k-1, j+1, nx+1, 3];    //buf[4,ipos2] = g[4,nx,j,k];
                            buf[4*size2+ipos2] = g[k-1, j+1, nx+1, 4];    //buf[5,ipos2] = g[5,nx,j,k];
                        }
                    }
                    //call MPI_SEND[ buf, 10*ny*nz, dp_type, south, from_n, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, south, from_n);
                }
                //---------------------------------------------------------------------
                //   receive from north
                //---------------------------------------------------------------------
                if(north!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++) {
                        for(j = 1; j<=ny; j++) {
                            ipos1 = (k-1)*ny + j           -1;     //ipos1 = (k-1)*ny + j;
                            ipos2 = ipos1 + ny*nz;                 //ipos2 = ipos1 + ny*nz; 
                            g[k-1, j+1, 0, 0] = buf1[0*size2+ipos1];     //g[1,-1,j,k] = buf1[1,ipos1];       
                            g[k-1, j+1, 0, 1] = buf1[1*size2+ipos1];     //g[2,-1,j,k] = buf1[2,ipos1];
                            g[k-1, j+1, 0, 2] = buf1[2*size2+ipos1];     //g[3,-1,j,k] = buf1[3,ipos1];
                            g[k-1, j+1, 0, 3] = buf1[3*size2+ipos1];     //g[4,-1,j,k] = buf1[4,ipos1];
                            g[k-1, j+1, 0, 4] = buf1[4*size2+ipos1];     //g[5,-1,j,k] = buf1[5,ipos1];

                            g[k-1, j+1, 1, 0] = buf1[0*size2+ipos2];      //g[1,0,j,k] = buf1[1,ipos2];
                            g[k-1, j+1, 1, 1] = buf1[1*size2+ipos2];      //g[2,0,j,k] = buf1[2,ipos2];
                            g[k-1, j+1, 1, 2] = buf1[2*size2+ipos2];      //g[3,0,j,k] = buf1[3,ipos2]; 
                            g[k-1, j+1, 1, 3] = buf1[3*size2+ipos2];      //g[4,0,j,k] = buf1[4,ipos2];
                            g[k-1, j+1, 1, 4] = buf1[4*size2+ipos2];      //g[5,0,j,k] = buf1[5,ipos2];
                        }
                    }
                }
                if(south!=-1) {
                    //call MPI_IRECV[buf1, 10*ny*nz, dp_type, MPI_ANY_SOURCE, from_s, MPI_COMM_WORLD, mid, IERROR];
                    mid[0] = worldcomm.ImmediateReceive<double>(south, from_s, buf1);
                }
                //---------------------------------------------------------------------
                //   send north
                //---------------------------------------------------------------------
                if(north!=-1) {
                    for(k = 1; k<=nz; k++) {
                        for(j = 1; j<=ny; j++) {
                            ipos1 = (k-1)*ny + j   -1;          //ipos1 = (k-1)*ny + j;
                            ipos2 = ipos1 + ny*nz;              //ipos2 = ipos1 + ny*nz;
                            buf[0*size2+ipos1] = g[k-1, j+1, 3, 0];  //buf[1,ipos1] = g[1,2,j,k];
                            buf[1*size2+ipos1] = g[k-1, j+1, 3, 1];  //buf[2,ipos1] = g[2,2,j,k];
                            buf[2*size2+ipos1] = g[k-1, j+1, 3, 2];  //buf[3,ipos1] = g[3,2,j,k];
                            buf[3*size2+ipos1] = g[k-1, j+1, 3, 3];  //buf[4,ipos1] = g[4,2,j,k];
                            buf[4*size2+ipos1] = g[k-1, j+1, 3, 4];  //buf[5,ipos1] = g[5,2,j,k];

                            buf[0*size2+ipos2] = g[k-1, j+1, 2, 0];  //buf[1,ipos2] = g[1,1,j,k];
                            buf[1*size2+ipos2] = g[k-1, j+1, 2, 1];  //buf[2,ipos2] = g[2,1,j,k];
                            buf[2*size2+ipos2] = g[k-1, j+1, 2, 2];  //buf[3,ipos2] = g[3,1,j,k];
                            buf[3*size2+ipos2] = g[k-1, j+1, 2, 3];  //buf[4,ipos2] = g[4,1,j,k];
                            buf[4*size2+ipos2] = g[k-1, j+1, 2, 4];  //buf[5,ipos2] = g[5,1,j,k];
                        }
                    }
                    //call MPI_SEND[ buf, 10*ny*nz, dp_type, north, from_s, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, north, from_s);
                }
                //---------------------------------------------------------------------
                //   receive from south
                //---------------------------------------------------------------------
                if(south!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++) {
                        for(j = 1; j<=ny; j++) {
                            ipos1 = (k-1)*ny + j                -1; //ipos1 = (k-1)*ny + j;
                            ipos2 = ipos1 + ny*nz;                  //ipos2 = ipos1 + ny*nz;
                            g[k-1, j+1, nx+3, 0]  = buf1[0*size2+ipos1]; //g[1,nx+2,j,k]  = buf1[1,ipos1];
                            g[k-1, j+1, nx+3, 1]  = buf1[1*size2+ipos1]; //g[2,nx+2,j,k]  = buf1[2,ipos1];
                            g[k-1, j+1, nx+3, 2]  = buf1[2*size2+ipos1]; //g[3,nx+2,j,k]  = buf1[3,ipos1];
                            g[k-1, j+1, nx+3, 3]  = buf1[3*size2+ipos1]; //g[4,nx+2,j,k]  = buf1[4,ipos1];
                            g[k-1, j+1, nx+3, 4]  = buf1[4*size2+ipos1]; //g[5,nx+2,j,k]  = buf1[5,ipos1];

                            g[k-1, j+1, nx+2, 0] = buf1[0*size2+ipos2];  //g[1,nx+1,j,k] = buf1[1,ipos2];
                            g[k-1, j+1, nx+2, 1] = buf1[1*size2+ipos2];  //g[2,nx+1,j,k] = buf1[2,ipos2];
                            g[k-1, j+1, nx+2, 2] = buf1[2*size2+ipos2];  //g[3,nx+1,j,k] = buf1[3,ipos2];
                            g[k-1, j+1, nx+2, 3] = buf1[3*size2+ipos2];  //g[4,nx+1,j,k] = buf1[4,ipos2];
                            g[k-1, j+1, nx+2, 4] = buf1[4*size2+ipos2];  //g[5,nx+1,j,k] = buf1[5,ipos2];
                        }
                    }
                }
            }
            else {
                bsize = 10*nx*nz;
                size2 = bsize/5;
                buf1 = new double[bsize];
                buf  = new double[bsize];
                //---------------------------------------------------------------------
                //   communicate in the east and west directions
                //---------------------------------------------------------------------
                if(west!=-1) {
                    //call MPI_IRECV[ buf1, 10*nx*nz, dp_type, MPI_ANY_SOURCE, from_w, MPI_COMM_WORLD, mid, IERROR ];
                    mid[0] = worldcomm.ImmediateReceive<double>(west, from_w, buf1);
                }
                //---------------------------------------------------------------------
                //   send east
                //---------------------------------------------------------------------
                if(east!=-1) {
                    for(k = 1; k<=nz; k++) {
                        for(i = 1; i<=nx; i++) {
                            ipos1 = (k-1)*nx+i        -1;         //ipos1 = (k-1)*nx+i;
                            ipos2 = ipos1+nx*nz;                  //ipos2 = ipos1+nx*nz;
                            buf[0*size2+ipos1] = g[k-1, ny, i+1, 0];   //buf[1,ipos1] = g[1,i,ny-1,k];
                            buf[1*size2+ipos1] = g[k-1, ny, i+1, 1];   //buf[2,ipos1] = g[2,i,ny-1,k];
                            buf[2*size2+ipos1] = g[k-1, ny, i+1, 2];   //buf[3,ipos1] = g[3,i,ny-1,k];
                            buf[3*size2+ipos1] = g[k-1, ny, i+1, 3];   //buf[4,ipos1] = g[4,i,ny-1,k];
                            buf[4*size2+ipos1] = g[k-1, ny, i+1, 4];   //buf[5,ipos1] = g[5,i,ny-1,k];

                            buf[0*size2+ipos2] = g[k-1, ny+1, i+1, 0];     //buf[1,ipos2] = g[1,i,ny,k];
                            buf[1*size2+ipos2] = g[k-1, ny+1, i+1, 1];     //buf[2,ipos2] = g[2,i,ny,k];
                            buf[2*size2+ipos2] = g[k-1, ny+1, i+1, 2];     //buf[3,ipos2] = g[3,i,ny,k];
                            buf[3*size2+ipos2] = g[k-1, ny+1, i+1, 3];     //buf[4,ipos2] = g[4,i,ny,k];
                            buf[4*size2+ipos2] = g[k-1, ny+1, i+1, 4];     //buf[5,ipos2] = g[5,i,ny,k];
                        }
                    }
                    //call MPI_SEND[ buf, 10*nx*nz, dp_type, east, from_w, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, east, from_w);
                }
                //---------------------------------------------------------------------
                //   receive from west
                //---------------------------------------------------------------------
                if(west!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++) {
                        for(i = 1; i<=nx; i++) {
                            ipos1 = (k-1)*nx + i     -1;           //ipos1 = (k-1)*nx + i;
                            ipos2 = ipos1 + nx*nz;                 //ipos2 = ipos1 + nx*nz;
                            g[k-1, 0, i+1, 0] = buf1[0*size2+ipos1];    //g[1,i,-1,k] = buf1[1,ipos1];
                            g[k-1, 0, i+1, 1] = buf1[1*size2+ipos1];    //g[2,i,-1,k] = buf1[2,ipos1];
                            g[k-1, 0, i+1, 2] = buf1[2*size2+ipos1];    //g[3,i,-1,k] = buf1[3,ipos1];
                            g[k-1, 0, i+1, 3] = buf1[3*size2+ipos1];    //g[4,i,-1,k] = buf1[4,ipos1];
                            g[k-1, 0, i+1, 4] = buf1[4*size2+ipos1];    //g[5,i,-1,k] = buf1[5,ipos1];

                            g[k-1, 1, i+1, 0] = buf1[0*size2+ipos2];     //g[1,i,0,k] = buf1[1,ipos2];
                            g[k-1, 1, i+1, 1] = buf1[1*size2+ipos2];     //g[2,i,0,k] = buf1[2,ipos2];
                            g[k-1, 1, i+1, 2] = buf1[2*size2+ipos2];     //g[3,i,0,k] = buf1[3,ipos2];
                            g[k-1, 1, i+1, 3] = buf1[3*size2+ipos2];     //g[4,i,0,k] = buf1[4,ipos2];
                            g[k-1, 1, i+1, 4] = buf1[4*size2+ipos2];     //g[5,i,0,k] = buf1[5,ipos2];
                        }
                    }
                }
                if(east!=-1) {
                    //call MPI_IRECV[ buf1, 10*nx*nz, dp_type, MPI_ANY_SOURCE, from_e, MPI_COMM_WORLD, mid, IERROR ];
                    mid[0] = worldcomm.ImmediateReceive<double>(east, from_e, buf1);
                }
                //---------------------------------------------------------------------
                //   send west
                //---------------------------------------------------------------------
                if(west!=-1) {
                    for(k = 1; k<=nz; k++) {
                        for(i = 1; i<=nx; i++) {
                            ipos1 = (k-1)*nx + i   -1;          //ipos1 = (k-1)*nx + i;
                            ipos2 = ipos1 + nx*nz;              //ipos2 = ipos1 + nx*nz;
                            buf[0*size2+ipos1] = g[k-1, 3, i+1, 0];  //buf[1,ipos1] = g[1,i,2,k];
                            buf[1*size2+ipos1] = g[k-1, 3, i+1, 1];  //buf[2,ipos1] = g[2,i,2,k];
                            buf[2*size2+ipos1] = g[k-1, 3, i+1, 2];  //buf[3,ipos1] = g[3,i,2,k];
                            buf[3*size2+ipos1] = g[k-1, 3, i+1, 3];  //buf[4,ipos1] = g[4,i,2,k];
                            buf[4*size2+ipos1] = g[k-1, 3, i+1, 4];  //buf[5,ipos1] = g[5,i,2,k];

                            buf[0*size2+ipos2] = g[k-1, 2, i+1, 0];  //buf[1,ipos2] = g[1,i,1,k];
                            buf[1*size2+ipos2] = g[k-1, 2, i+1, 1];  //buf[2,ipos2] = g[2,i,1,k];
                            buf[2*size2+ipos2] = g[k-1, 2, i+1, 2];  //buf[3,ipos2] = g[3,i,1,k];
                            buf[3*size2+ipos2] = g[k-1, 2, i+1, 3];  //buf[4,ipos2] = g[4,i,1,k];
                            buf[4*size2+ipos2] = g[k-1, 2, i+1, 4];  //buf[5,ipos2] = g[5,i,1,k];
                        }
                    }
                    //call MPI_SEND[buf, 10*nx*nz, dp_type, west, from_e, MPI_COMM_WORLD, IERROR ];
                    worldcomm.Send<double>(buf, west, from_e);
                }
                //---------------------------------------------------------------------
                //   receive from east
                //---------------------------------------------------------------------
                if(east!=-1) {
                    //call MPI_WAIT[ mid, STATUS, IERROR ];
                    mid[0].Wait();
                    for(k = 1; k<=nz; k++) {
                        for(i = 1; i<=nx; i++) {
                            ipos1 = (k-1)*nx + i        -1;         //ipos1 = (k-1)*nx + i;
                            ipos2 = ipos1 + nx*nz;                  //ipos2 = ipos1 + nx*nz;
                            g[k-1, ny+3, i+1, 0]  = buf1[0*size2+ipos1]; //g[1,i,ny+2,k]  = buf1[1,ipos1];
                            g[k-1, ny+3, i+1, 1]  = buf1[1*size2+ipos1]; //g[2,i,ny+2,k]  = buf1[2,ipos1];
                            g[k-1, ny+3, i+1, 2]  = buf1[2*size2+ipos1]; //g[3,i,ny+2,k]  = buf1[3,ipos1];
                            g[k-1, ny+3, i+1, 3]  = buf1[3*size2+ipos1]; //g[4,i,ny+2,k]  = buf1[4,ipos1];
                            g[k-1, ny+3, i+1, 4]  = buf1[4*size2+ipos1]; //g[5,i,ny+2,k]  = buf1[5,ipos1];

                            g[k-1, ny+2, i+1, 0] = buf1[0*size2+ipos2];  //g[1,i,ny+1,k] = buf1[1,ipos2];
                            g[k-1, ny+2, i+1, 1] = buf1[1*size2+ipos2];  //g[2,i,ny+1,k] = buf1[2,ipos2];
                            g[k-1, ny+2, i+1, 2] = buf1[2*size2+ipos2];  //g[3,i,ny+1,k] = buf1[3,ipos2];
                            g[k-1, ny+2, i+1, 3] = buf1[3*size2+ipos2];  //g[4,i,ny+1,k] = buf1[4,ipos2];
                            g[k-1, ny+2, i+1, 4] = buf1[4*size2+ipos2];  //g[5,i,ny+1,k] = buf1[5,ipos2];
                        }
                    }
                }
            }
        }
        //end Exchange_3.f
        //End erhs.f
        //ssor.f
        public void ssor(int niter) {
            //---------------------------------------------------------------------
            //   to perform pseudo-time stepping SSOR iterations
            //   for five nonlinear pde's.
            //---------------------------------------------------------------------
            //int  niter;
            int i, j, k, m;
            int istep;
            double  tmp;
            double[] delunm = new double[5];//delunm[5];
            double[,,] tv = new double[isiz2+1, isiz1+1, 5+1];//tv[5,isiz1,isiz2];
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
            for(m=0; m<isiz2; m++) {
                for(k=0; k<isiz1; k++) {
                    for(j=0; j<5; j++) {
                        for(i=0; i<5; i++) {
                            a[m, k, j, i] = 0.0;//a[i,j,k,m] = 0.0;
                            b[m, k, j, i] = 0.0;//b[i,j,k,m] = 0.0;
                            c[m, k, j, i] = 0.0;//c[i,j,k,m] = 0.0;
                            d[m, k, j, i] = 0.0;//d[i,j,k,m] = 0.0;
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
            for(istep = 1; istep<= niter; istep++) {
                if(node == 0) {
                    if(mod(istep, 20) == 0 || istep == itmax || istep == 1) {
                        if(niter > 1)
                            Console.WriteLine(" Time step "+istep); //write[ *, 200] istep       200           format[' Time step ', i4]
                    }
                }
                //---------------------------------------------------------------------
                //   perform SSOR iteration
                //---------------------------------------------------------------------
                for(k = 1; k< nz - 1; k++) {
                    for(j = jst; j<= jend; j++) {
                        for(i = ist; i<= iend; i++) {
                            for(m = 0; m< 5; m++) {
                                rsd[k, j+1, i+1, m] = dt * rsd[k, j+1, i+1, m];   //rsd[m, i, j, k] = dt * rsd[m, i, j, k];
                            }
                        }
                    }
                }
                for(k = 2; k<= nz -1; k++) {
                    //---------------------------------------------------------------------
                    //   form the lower triangular part of the jacobian matrix
                    //---------------------------------------------------------------------
                    jacld(k);
                    //---------------------------------------------------------------------
                    //   perform the lower triangular solution
                    //---------------------------------------------------------------------
                    blts(isiz1, isiz2, isiz3, nx, ny, nz, k, omega, rsd, a, b, c, d, ist, iend, jst, jend, nx0, ny0, ipt, jpt);
                }
                for(k=nz-1; k>= 2; k--) { //for(k = nz - 1, 2, -1;
                    //---------------------------------------------------------------------
                    //   form the strictly upper triangular part of the jacobian matrix
                    //---------------------------------------------------------------------
                    jacu(k);
                    //---------------------------------------------------------------------
                    //   perform the upper triangular solution
                    //---------------------------------------------------------------------
                    buts(isiz1, isiz2, isiz3, nx, ny, nz, k, omega, rsd, tv, d, a, b, c, ist, iend, jst, jend, nx0, ny0, ipt, jpt);
                }
                //---------------------------------------------------------------------
                //   update the variables
                //---------------------------------------------------------------------
                for(k = 1; k< nz-1; k++) {
                    for(j = jst; j<= jend; j++) {
                        for(i = ist; i<= iend; i++) {
                            for(m = 0; m< 5; m++) {     //u[ m, i, j, k ] = u[ m, i, j, k ] + tmp * rsd[ m, i, j, k ];
                                u[k, j+1, i+1, m] = u[k, j+1, i+1, m] + tmp * rsd[k, j+1, i+1, m];
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //   compute the max-norms of newton iteration corrections
                //---------------------------------------------------------------------
                if(mod(istep, inorm) == 0) {
                    l2norm(isiz1, isiz2, isiz3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, delunm);
                }
                //---------------------------------------------------------------------
                //   compute the steady-state residuals
                //---------------------------------------------------------------------
                rhs();
                //---------------------------------------------------------------------
                //   compute the max-norms of newton iteration residuals
                //---------------------------------------------------------------------
                if((mod(istep, inorm)== 0) || (istep==itmax)) {
                    l2norm(isiz1, isiz2, isiz3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
                }
                //---------------------------------------------------------------------
                //   check the newton-iteration residuals against the tolerance levels
                //---------------------------------------------------------------------
                if((rsdnm[0]<tolrsd[0]) && (rsdnm[1]<tolrsd[1]) && (rsdnm[2]<tolrsd[2]) && (rsdnm[3]<tolrsd[3]) && (rsdnm[4]<tolrsd[4])) {
                    return;//   return;
                }
            }
            timer.stop(1); //call timer_stop[1];
            wtime = timer.readTimer(1); //wtime = timer_read[1];
            //call MPI_ALLREDUCE[wtime, maxtime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR ];
            maxtime = worldcomm.Allreduce<double>(wtime, MPI.Operation<double>.Max);
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

            for(k = 0; k< nz; k++) {
                for(j = 1; j<= ny; j++) {
                    for(i = 1; i<= nx; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k, j+1, i+1, m] = -frct[k, j+1, i+1, m];//rsd[m,i,j,k] = - frct[m,i,j,k];
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
            exchange_3(u, iex);
            L1 = 0;
            if(north==-1)
                L1 = 1;
            L2 = nx + 1;
            if(south==-1)
                L2 = nx;

            ist1 = 1;
            iend1 = nx;
            if(north==-1)
                ist1 = 4;
            if(south==-1)
                iend1 = nx - 3;

            for(k = 2; k<= nz - 1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = L1; i<= L2; i++) {
                        flux[k-1, j, i, 0] = u[k-1, j+1, i+1, 1];  //flux[1,i,j,k] = u[2,i,j,k];
                        u21=u[k-1, j+1, i+1, 1]/u[k-1, j+1, i+1, 0];  //u21 = u[2,i,j,k] / u[1,i,j,k];

                        q = 0.50d*(u[k-1, j+1, i+1, 1]*u[k-1, j+1, i+1, 1]+u[k-1, j+1, i+1, 2]*u[k-1, j+1, i+1, 2]+u[k-1, j+1, i+1, 3]*u[k-1, j+1, i+1, 3])/u[k-1, j+1, i+1, 0];//q = 0.50d*(u[2,i,j,k]*u[2,i,j,k]+u[3,i,j,k]*u[3,i,j,k]+u[4,i,j,k]*u[4,i,j,k])/u[1,i,j,k];

                        flux[k-1, j, i, 1] =     u[k-1, j+1, i+1, 1] * u21 + c2 *(u[k-1, j+1, i+1, 4] - q);   //flux[2,i,j,k]=u[2,i,j,k]*u21+c2*(u[5,i,j,k]-q);
                        flux[k-1, j, i, 2] =     u[k-1, j+1, i+1, 2] * u21;                              //flux[3,i,j,k]=u[3,i,j,k]*u21;
                        flux[k-1, j, i, 3] =     u[k-1, j+1, i+1, 3] * u21;                              //flux[4,i,j,k]=u[4,i,j,k]*u21;
                        flux[k-1, j, i, 4] = (c1*u[k-1, j+1, i+1, 4]-c2*q)*u21;                          //flux[5,i,j,k]=(c1*u[5,i,j,k]-c2*q)*u21;
                    }
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, j+1, i+1, m]=rsd[k-1, j+1, i+1, m]-tx2*(flux[k-1, j, i+1, m] - flux[k-1, j, i-1, m]);  //rsd[m,i,j,k] =  rsd[m,i,j,k]- tx2 *(flux[m,i+1,j,k] - flux[m,i-1,j,k]);
                        }
                    }
                    for(i = ist; i<= L2; i++) {
                        tmp = 1.0d/u[k-1, j+1, i+1, 0];
                        u21i = tmp*u[k-1, j+1, i+1, 1];
                        u31i = tmp*u[k-1, j+1, i+1, 2];
                        u41i = tmp*u[k-1, j+1, i+1, 3];
                        u51i = tmp*u[k-1, j+1, i+1, 4];

                        tmp = 1.0d/ u[k-1, j+1, i, 0];

                        u21im1 = tmp * u[k-1, j+1, i, 1];
                        u31im1 = tmp * u[k-1, j+1, i, 2];
                        u41im1 = tmp * u[k-1, j+1, i, 3];
                        u51im1 = tmp * u[k-1, j+1, i, 4];

                        flux[k-1, j, i, 1] = (4.0d/3.0d)*tx3*(u21i-u21im1);
                        flux[k-1, j, i, 2] = tx3 * (u31i - u31im1);
                        flux[k-1, j, i, 3] = tx3 * (u41i - u41im1);
                        flux[k-1, j, i, 4] = 0.50d*(1.0d-c1*c5)*tx3*((pow2(u21i)+pow2(u31i)+pow2(u41i))-(pow2(u21im1)+pow2(u31im1)+pow2(u41im1)))+(1.0d/6.0d)*tx3*(pow2(u21i)-pow2(u21im1))+c1*c5*tx3*(u51i-u51im1);
                    }
                    for(i = ist; i<= iend; i++) {
                        rsd[k-1, j+1, i+1, 0]=rsd[k-1, j+1, i+1, 0]+dx1*tx1*(u[k-1, j+1, i, 0]-2.0d*u[k-1, j+1, i+1, 0]+u[k-1, j+1, i+2, 0]);
                        rsd[k-1, j+1, i+1, 1]=rsd[k-1, j+1, i+1, 1]+tx3*c3*c4*(flux[k-1, j, i+1, 1]-flux[k-1, j, i, 1])+dx2*tx1*(u[k-1, j+1, i, 1]-2.0d*u[k-1, j+1, i+1, 1]+u[k-1, j+1, i+2, 1]);
                        rsd[k-1, j+1, i+1, 2]=rsd[k-1, j+1, i+1, 2]+tx3*c3*c4*(flux[k-1, j, i+1, 2]-flux[k-1, j, i, 2])+dx3*tx1*(u[k-1, j+1, i, 2]-2.0d*u[k-1, j+1, i+1, 2]+u[k-1, j+1, i+2, 2]);
                        rsd[k-1, j+1, i+1, 3]=rsd[k-1, j+1, i+1, 3]+tx3*c3*c4*(flux[k-1, j, i+1, 3]-flux[k-1, j, i, 3])+dx4*tx1*(u[k-1, j+1, i, 3]-2.0d*u[k-1, j+1, i+1, 3]+u[k-1, j+1, i+2, 3]);
                        rsd[k-1, j+1, i+1, 4]=rsd[k-1, j+1, i+1, 4]+tx3*c3*c4*(flux[k-1, j, i+1, 4]-flux[k-1, j, i, 4])+dx5*tx1*(u[k-1, j+1, i, 4]-2.0d*u[k-1, j+1, i+1, 4]+u[k-1, j+1, i+2, 4]);
                    }
                    //---------------------------------------------------------------------
                    //   Fourth-order dissipation
                    //---------------------------------------------------------------------
                    if(north==-1) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, j+1, 3, m] = rsd[k-1, j+1, 3, m]-dssp*(+5.0d*u[k-1, j+1, 3, m]-4.0d*u[k-1, j+1, 4, m]+     u[k-1, j+1, 5, m]);
                            rsd[k-1, j+1, 4, m] = rsd[k-1, j+1, 4, m]-dssp*(-4.0d*u[k-1, j+1, 3, m]+6.0d*u[k-1, j+1, 4, m]-4.0d*u[k-1, j+1, 5, m]+u[k-1, j+1, 6, m]);
                        }
                    }
                    for(i = ist1; i<=iend1; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, j+1, i+1, m] = rsd[k-1, j+1, i+1, m]-dssp*(u[k-1, j+1, i-1, m]-4.0d*u[k-1, j+1, i, m]+6.0d*u[k-1, j+1, i+1, m]-4.0d*u[k-1, j+1, i+2, m]+u[k-1, j+1, i+3, m]);
                        }
                    }
                    if(south==-1) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, j+1, nx-1, m] = rsd[k-1, j+1, nx-1, m]-dssp*(u[k-1, j+1, nx-3, m]-4.0d*u[k-1, j+1, nx-2, m]+6.0d*u[k-1, j+1, nx-1, m]-4.0d*u[k-1, j+1, nx, m]);
                            rsd[k-1, j+1, nx, m]   = rsd[k-1, j+1, nx, m]  -dssp*(u[k-1, j+1, nx-2, m]-4.0d*u[k-1, j+1, nx-1, m]+5.0d*u[k-1, j+1, nx, m]);
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
            exchange_3(u, iex);

            L1 = 0;
            if(west==-1)
                L1 = 1;
            L2 = ny + 1;
            if(east==-1)
                L2 = ny;

            jst1 = 1;
            jend1 = ny;
            if(west==-1)
                jst1 = 4;
            if(east==-1)
                jend1 = ny - 3;

            for(k = 2; k<= nz - 1; k++) {
                for(j = L1; j<= L2; j++) {
                    for(i = ist; i<= iend; i++) {
                        flux[k-1, j, i, 0]=u[k-1, j+1, i+1, 2];          //flux[1,i,j,k] = u[3,i,j,k];
                        u31=u[k-1, j+1, i+1, 2]/u[k-1, j+1, i+1, 0];     //u31 = u[3,i,j,k] / u[1,i,j,k];

                        q = 0.50d*(u[k-1, j+1, i+1, 1]*u[k-1, j+1, i+1, 1]+u[k-1, j+1, i+1, 2]*u[k-1, j+1, i+1, 2]+u[k-1, j+1, i+1, 3]*u[k-1, j+1, i+1, 3])/u[k-1, j+1, i+1, 0];      //q = 0.50d*(u[2,i,j,k]*u[2,i,j,k]+u[3,i,j,k]*u[3,i,j,k]+u[4,i,j,k]*u[4,i,j,k])/u[1,i,j,k];

                        flux[k-1, j, i, 1] =     u[k-1, j+1, i+1, 1] * u31;                              //flux[2,i,j,k]=u[2,i,j,k] * u31;
                        flux[k-1, j, i, 2] =     u[k-1, j+1, i+1, 2] * u31 + c2 * (u[k-1, j+1, i+1, 4]-q);    //flux[3,i,j,k]=u[3,i,j,k]*u31+c2*(u[5,i,j,k]-q);
                        flux[k-1, j, i, 3] =     u[k-1, j+1, i+1, 3] * u31;                              //flux[4,i,j,k] =     u[4,i,j,k] * u31;
                        flux[k-1, j, i, 4] = (c1*u[k-1, j+1, i+1, 4]-c2*q)*u31;                          //flux[5,i,j,k] = (c1*u[5,i,j,k]-c2*q)*u31;
                    }
                }
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {  //rsd[m,i,j,k] =  rsd[m,i,j,k]- ty2 * [ flux[m,i,j+1,k] - flux[m,i,j-1,k] ];
                            rsd[k-1, j+1, i+1, m] =  rsd[k-1, j+1, i+1, m]- ty2 * (flux[k-1, j+1, i, m] - flux[k-1, j-1, i, m]);
                        }
                    }
                }
                for(j = jst; j<= L2; j++) {
                    for(i = ist; i<= iend; i++) {
                        tmp = 1.0d / u[k-1, j+1, i+1, 0];
                        u21j = tmp * u[k-1, j+1, i+1, 1];
                        u31j = tmp * u[k-1, j+1, i+1, 2];
                        u41j = tmp * u[k-1, j+1, i+1, 3];
                        u51j = tmp * u[k-1, j+1, i+1, 4];

                        tmp = 1.0d / u[k-1, j, i+1, 0];
                        u21jm1=tmp * u[k-1, j, i+1, 1];
                        u31jm1=tmp * u[k-1, j, i+1, 2];
                        u41jm1=tmp * u[k-1, j, i+1, 3];
                        u51jm1=tmp * u[k-1, j, i+1, 4];

                        flux[k-1, j, i, 1] = ty3 * (u21j - u21jm1);
                        flux[k-1, j, i, 2] = (4.0d/3.0d) * ty3 * (u31j-u31jm1);
                        flux[k-1, j, i, 3] = ty3 * (u41j - u41jm1);
                        flux[k-1, j, i, 4] = 0.50d*(1.0d-c1*c5)*ty3*((pow2(u21j)+pow2(u31j)+pow2(u41j))-(pow2(u21jm1)+pow2(u31jm1)+pow2(u41jm1)))+(1.0d/6.0d)*ty3*(pow2(u31j)-pow2(u31jm1))+c1*c5*ty3*(u51j-u51jm1);
                    }
                }
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        rsd[k-1, j+1, i+1, 0] = rsd[k-1, j+1, i+1, 0]+dy1*ty1*(u[k-1, j, i+1, 0]-2.0d*u[k-1, j+1, i+1, 0]+u[k-1, j+2, i+1, 0]);
                        rsd[k-1, j+1, i+1, 1] = rsd[k-1, j+1, i+1, 1]+ty3*c3*c4*(flux[k-1, j+1, i, 1]-flux[k-1, j, i, 1])+dy2*ty1*(u[k-1, j, i+1, 1]-2.0d*u[k-1, j+1, i+1, 1]+u[k-1, j+2, i+1, 1]);
                        rsd[k-1, j+1, i+1, 2] = rsd[k-1, j+1, i+1, 2]+ty3*c3*c4*(flux[k-1, j+1, i, 2]-flux[k-1, j, i, 2])+dy3*ty1*(u[k-1, j, i+1, 2]-2.0d*u[k-1, j+1, i+1, 2]+u[k-1, j+2, i+1, 2]);
                        rsd[k-1, j+1, i+1, 3] = rsd[k-1, j+1, i+1, 3]+ty3*c3*c4*(flux[k-1, j+1, i, 3]-flux[k-1, j, i, 3])+dy4*ty1*(u[k-1, j, i+1, 3]-2.0d*u[k-1, j+1, i+1, 3]+u[k-1, j+2, i+1, 3]);
                        rsd[k-1, j+1, i+1, 4] = rsd[k-1, j+1, i+1, 4]+ty3*c3*c4*(flux[k-1, j+1, i, 4]-flux[k-1, j, i, 4])+dy5*ty1*(u[k-1, j, i+1, 4]-2.0d*u[k-1, j+1, i+1, 4]+u[k-1, j+2, i+1, 4]);
                    }
                }
                //---------------------------------------------------------------------
                //   fourth-order dissipation
                //---------------------------------------------------------------------
                if(west==-1) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, 3, i+1, m] = rsd[k-1, 3, i+1, m]-dssp*(+5.0d*u[k-1, 3, i+1, m]-4.0d*u[k-1, 4, i+1, m]+       u[k-1, 5, i+1, m]);
                            rsd[k-1, 4, i+1, m] = rsd[k-1, 4, i+1, m]-dssp*(-4.0d*u[k-1, 3, i+1, m]+6.0d*u[k-1, 4, i+1, m]-4.0d * u[k-1, 5, i+1, m]+u[k-1, 6, i+1, m]);
                        }
                    }
                }
                for(j = jst1; j<= jend1; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, j+1, i+1, m]=rsd[k-1, j+1, i+1, m]-dssp*(u[k-1, j-1, i+1, m]-4.0d*u[k-1, j, i+1, m]+6.0d*u[k-1, j+1, i+1, m]-4.0d*u[k-1, j+2, i+1, m]+u[k-1, j+3, i+1, m]);
                        }
                    }
                }
                if(east==-1) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, ny-1, i+1, m]=rsd[k-1, ny-1, i+1, m]-dssp*(u[k-1, ny-3, i+1, m]-4.0d*u[k-1, ny-2, i+1, m]+6.0d*u[k-1, ny-1, i+1, m]-4.0d*u[k-1, ny, i+1, m]);
                            rsd[k-1, ny, i+1, m]=rsd[k-1, ny, i+1, m]-dssp*(u[k-1, ny-2, i+1, m]-4.0d*u[k-1, ny-1, i+1, m]+5.0d*u[k-1, ny, i+1, m]);
                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //   zeta-direction flux differences
            //---------------------------------------------------------------------
            for(k = 1; k<= nz; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        flux[k-1, j, i, 0]=u[k-1, j+1, i+1, 3];
                        u41=u[k-1, j+1, i+1, 3]/u[k-1, j+1, i+1, 0];

                        q = 0.50d * (u[k-1, j+1, i+1, 1] * u[k-1, j+1, i+1, 1]+ u[k-1, j+1, i+1, 2] * u[k-1, j+1, i+1, 2]+ u[k-1, j+1, i+1, 3] * u[k-1, j+1, i+1, 3])/u[k-1, j+1, i+1, 0];

                        flux[k-1, j, i, 1] =   u[k-1, j+1, i+1, 1] * u41;
                        flux[k-1, j, i, 2] =   u[k-1, j+1, i+1, 2] * u41;
                        flux[k-1, j, i, 3] =   u[k-1, j+1, i+1, 3] * u41 + c2 * (u[k-1, j+1, i+1, 4]-q);
                        flux[k-1, j, i, 4]=(c1*u[k-1, j+1, i+1, 4]-c2*q)*u41;
                    }
                }
            }
            for(k = 2; k<= nz - 1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, j+1, i+1, m] =  rsd[k-1, j+1, i+1, m]- tz2 * (flux[k, j, i, m] - flux[k-2, j, i, m]);
                        }
                    }
                }
            }
            for(k = 2; k<= nz; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        tmp = 1.0d / u[k-1, j+1, i+1, 0];
                        u21k = tmp * u[k-1, j+1, i+1, 1];
                        u31k = tmp * u[k-1, j+1, i+1, 2];
                        u41k = tmp * u[k-1, j+1, i+1, 3];
                        u51k = tmp * u[k-1, j+1, i+1, 4];

                        tmp   = 1.0d / u[k-2, j+1, i+1, 0];
                        u21km1 = tmp * u[k-2, j+1, i+1, 1];
                        u31km1 = tmp * u[k-2, j+1, i+1, 2];
                        u41km1 = tmp * u[k-2, j+1, i+1, 3];
                        u51km1 = tmp * u[k-2, j+1, i+1, 4];

                        flux[k-1, j, i, 1] = tz3 * (u21k - u21km1);
                        flux[k-1, j, i, 2] = tz3 * (u31k - u31km1);
                        flux[k-1, j, i, 3] = (4.0d/3.0d) * tz3 * (u41k-u41km1);
                        flux[k-1, j, i, 4] = 0.50d*(1.0d-c1*c5)*tz3*((pow2(u21k)+pow2(u31k)+pow2(u41k))-(pow2(u21km1)+pow2(u31km1)+pow2(u41km1)))+(1.0d/6.0d)*tz3*(pow2(u41k)-pow2(u41km1))+c1*c5*tz3*(u51k-u51km1);
                    }
                }
            }
            for(k = 2; k<= nz - 1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        rsd[k-1, j+1, i+1, 0] = rsd[k-1, j+1, i+1, 0]+dz1*tz1*(u[k-2, j+1, i+1, 0]-2.0d*u[k-1, j+1, i+1, 0]+u[k, j+1, i+1, 0]);
                        rsd[k-1, j+1, i+1, 1] = rsd[k-1, j+1, i+1, 1]+tz3*c3*c4*(flux[k, j, i, 1]-flux[k-1, j, i, 1])+dz2*tz1*(u[k-2, j+1, i+1, 1]-2.0d*u[k-1, j+1, i+1, 1]+u[k, j+1, i+1, 1]);
                        rsd[k-1, j+1, i+1, 2] = rsd[k-1, j+1, i+1, 2]+tz3*c3*c4*(flux[k, j, i, 2]-flux[k-1, j, i, 2])+dz3*tz1*(u[k-2, j+1, i+1, 2]-2.0d*u[k-1, j+1, i+1, 2]+u[k, j+1, i+1, 2]);
                        rsd[k-1, j+1, i+1, 3] = rsd[k-1, j+1, i+1, 3]+tz3*c3*c4*(flux[k, j, i, 3]-flux[k-1, j, i, 3])+dz4*tz1*(u[k-2, j+1, i+1, 3]-2.0d*u[k-1, j+1, i+1, 3]+u[k, j+1, i+1, 3]);
                        rsd[k-1, j+1, i+1, 4] = rsd[k-1, j+1, i+1, 4]+tz3*c3*c4*(flux[k, j, i, 4]-flux[k-1, j, i, 4])+dz5*tz1*(u[k-2, j+1, i+1, 4]-2.0d*u[k-1, j+1, i+1, 4]+u[k, j+1, i+1, 4]);
                    }
                }
            }
            //---------------------------------------------------------------------
            //   fourth-order dissipation
            //---------------------------------------------------------------------
            for(j = jst; j<= jend; j++) {
                for(i = ist; i<= iend; i++) {
                    for(m = 0; m< 5; m++) {
                        rsd[-1+2, j+1, i+1, m]=rsd[-1+2, j+1, i+1, m]-dssp*(+5.0d*u[-1+2, j+1, i+1, m]-4.0d*u[-1+3, j+1, i+1, m]+     u[-1+4, j+1, i+1, m]);
                        rsd[-1+3, j+1, i+1, m]=rsd[-1+3, j+1, i+1, m]-dssp*(-4.0d*u[-1+2, j+1, i+1, m]+6.0d*u[-1+3, j+1, i+1, m]-4.0d*u[-1+4, j+1, i+1, m]+u[-1+5, j+1, i+1, m]);
                    }
                }
            }
            for(k = 4; k<= nz - 3; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            rsd[k-1, j+1, i+1, m] = rsd[k-1, j+1, i+1, m]-dssp*(u[k-3, j+1, i+1, m]-4.0d*u[k-2, j+1, i+1, m]+6.0d*u[k-1, j+1, i+1, m]-4.0d*u[k, j+1, i+1, m]+u[k+1, j+1, i+1, m]);
                        }
                    }
                }
            }
            for(j = jst; j<= jend; j++) {
                for(i = ist; i<= iend; i++) {
                    for(m = 0; m< 5; m++) {
                        rsd[-1+nz-2, j+1, i+1, m]=rsd[-1+nz-2, j+1, i+1, m]-dssp*(u[-1+nz-4, j+1, i+1, m]-4.0d*u[-1+nz-3, j+1, i+1, m]+6.0d*u[-1+nz-2, j+1, i+1, m]-4.0d*u[-1+nz-1, j+1, i+1, m]);
                        rsd[-1+nz-1, j+1, i+1, m]=rsd[-1+nz-1, j+1, i+1, m]-dssp*(u[-1+nz-3, j+1, i+1, m]-4.0d*u[-1+nz-2, j+1, i+1, m]+5.0d*u[-1+nz-1, j+1, i+1, m]);
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

            //debug
            //if ((isiz1 + 2) != (ldx + 2) || (isiz2 + 2) != (ldy + 2)) {
            //throw new ArgumentException("Look this code: vetor v");
            //}

            int i, j, k, m;
            double[] dummy = new double[5];//dummy[5];

            //int IERROR;

            for(m = 0; m< 5; m++) {
                dummy[m] = 0.0d;
            }
            for(k = 2; k<= nz0-1; k++) {
                for(j = jst; j<= jend; j++) {
                    for(i = ist; i<= iend; i++) {
                        for(m = 0; m< 5; m++) {
                            dummy[m] = dummy[m] + v[k-1, j+1, i+1, m] * v[k-1, j+1, i+1, m];  //dummy[m] = dummy[m] + v[m,i,j,k] * v[m,i,j,k];
                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //   compute the global sum of individual contributions to dot product.
            //---------------------------------------------------------------------
            worldcomm.Allreduce<double>(dummy, MPI.Operation<double>.Add, ref sum);//call MPI_ALLREDUCE[ dummy,sum,5,dp_type,MPI_SUM,MPI_COMM_WORLD,IERROR ]

            for(m = 0; m< 5; m++) {
                sum[m] = Math.Sqrt(sum[m]/((nx0-2)*(ny0-2)*(nz0-2))); //sum[m] = sqrt(sum[m]/((nx0-2)*(ny0-2)*(nz0-2)));
            }
        }//sum  dummy
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

            for(j = jst; j<= jend; j++) {
                for(i = ist; i<= iend; i++) {
                    //---------------------------------------------------------------------
                    //   form the block daigonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k-1, j+1, i+1, 0];  //tmp1 = 1.0d / u[1, i, j, k];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    d[j-1, i-1, 0, 0] =  1.0d+ dt * 2.0d * (tx1 * dx1+ ty1 * dy1+ tz1 * dz1);
                    d[j-1, i-1, 1, 0] =  0.0d;
                    d[j-1, i-1, 2, 0] =  0.0d;
                    d[j-1, i-1, 3, 0] =  0.0d;
                    d[j-1, i-1, 4, 0] =  0.0d;

                    d[j-1, i-1, 0, 1] =  dt*2.0d*(tx1*(-r43*c34*tmp2*u[k-1, j+1, i+1, 1])+ty1*(-c34*tmp2*u[k-1, j+1, i+1, 1])+tz1*(-c34*tmp2*u[k-1, j+1, i+1, 1]));
                    d[j-1, i-1, 1, 1] =  1.0d+dt*2.0d*(tx1*r43*c34*tmp1+ty1*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx2+ty1*dy2+tz1*dz2);
                    d[j-1, i-1, 2, 1] = 0.0d;
                    d[j-1, i-1, 3, 1] = 0.0d;
                    d[j-1, i-1, 4, 1] = 0.0d;

                    d[j-1, i-1, 0, 2] = dt*2.0d*(tx1*(-c34*tmp2*u[k-1, j+1, i+1, 2])+ty1*(-r43*c34*tmp2*u[k-1, j+1, i+1, 2])+tz1*(-c34*tmp2*u[k-1, j+1, i+1, 2]));
                    d[j-1, i-1, 1, 2] = 0.0d;
                    d[j-1, i-1, 2, 2] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*r43*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx3+ty1*dy3+tz1*dz3);
                    d[j-1, i-1, 3, 2] = 0.0d;
                    d[j-1, i-1, 4, 2] = 0.0d;

                    d[j-1, i-1, 0, 3] = dt*2.0d*(tx1*(-c34*tmp2*u[k-1, j+1, i+1, 3])+ty1*(-c34*tmp2*u[k-1, j+1, i+1, 3])+tz1*(-r43*c34*tmp2*u[k-1, j+1, i+1, 3]));
                    d[j-1, i-1, 1, 3] = 0.0d;
                    d[j-1, i-1, 2, 3] = 0.0d;
                    d[j-1, i-1, 3, 3] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*c34*tmp1+tz1*r43*c34*tmp1)+dt*2.0d*(tx1*dx4+ty1*dy4+tz1*dz4);
                    d[j-1, i-1, 4, 3] = 0.0d;

                    d[j-1, i-1, 0, 4] = dt*2.0d*(tx1*(-(r43*c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 1]))-(c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 2]))-(c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 3]))-(c1345)*tmp2*u[k-1, j+1, i+1, 4])+ty1*(-(c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 1]))-(r43*c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 2]))-(c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 3]))-(c1345)*tmp2*u[k-1, j+1, i+1, 4])+tz1*(-(c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 1]))-(c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 2]))-(r43*c34-c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 3]))-(c1345)*tmp2*u[k-1, j+1, i+1, 4]));
                    d[j-1, i-1, 1, 4] = dt*2.0d*(tx1*(r43*c34-c1345)*tmp2*u[k-1, j+1, i+1, 1]+ty1*(c34-c1345)*tmp2*u[k-1, j+1, i+1, 1]+tz1*(c34-c1345)*tmp2*u[k-1, j+1, i+1, 1]);
                    d[j-1, i-1, 2, 4] = dt*2.0d*(tx1*(c34-c1345)*tmp2*u[k-1, j+1, i+1, 2]+ty1*(r43*c34-c1345)*tmp2*u[k-1, j+1, i+1, 2]+tz1*(c34-c1345)*tmp2*u[k-1, j+1, i+1, 2]);
                    d[j-1, i-1, 3, 4] = dt*2.0d*(tx1*(c34-c1345)*tmp2*u[k-1, j+1, i+1, 3]+ty1*(c34-c1345)*tmp2*u[k-1, j+1, i+1, 3]+tz1*(r43*c34-c1345)*tmp2*u[k-1, j+1, i+1, 3]);
                    d[j-1, i-1, 4, 4] = 1.0d+dt*2.0d*(tx1*c1345*tmp1+ty1*c1345*tmp1+tz1*c1345*tmp1)+dt*2.0d*(tx1*dx5+ty1*dy5+tz1*dz5);
                    //---------------------------------------------------------------------
                    //   form the first block sub-diagonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d/u[k-2, j+1, i+1, 0];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    a[j-1, i-1, 0, 0] = -dt * tz1 * dz1;
                    a[j-1, i-1, 1, 0] =   0.0d;
                    a[j-1, i-1, 2, 0] =   0.0d;
                    a[j-1, i-1, 3, 0] = -dt * tz2;
                    a[j-1, i-1, 4, 0] =   0.0d;

                    a[j-1, i-1, 0, 1] = -dt*tz2*(-(u[k-2, j+1, i+1, 1]*u[k-2, j+1, i+1, 3])*tmp2)-dt*tz1*(-c34*tmp2*u[k-2, j+1, i+1, 1]);
                    a[j-1, i-1, 1, 1] = -dt*tz2*(u[k-2, j+1, i+1, 3]*tmp1)-dt*tz1*c34*tmp1-dt*tz1*dz2;
                    a[j-1, i-1, 2, 1] = 0.0d;
                    a[j-1, i-1, 3, 1] = -dt*tz2*(u[k-2, j+1, i+1, 1]*tmp1);
                    a[j-1, i-1, 4, 1] = 0.0d;

                    a[j-1, i-1, 0, 2] = -dt * tz2* (-(u[k-2, j+1, i+1, 2]*u[k-2, j+1, i+1, 3]) * tmp2)- dt * tz1 * (-c34 * tmp2 * u[k-2, j+1, i+1, 2]);
                    a[j-1, i-1, 1, 2] = 0.0d;
                    a[j-1, i-1, 2, 2] = -dt * tz2 * (u[k-2, j+1, i+1, 3] * tmp1)- dt * tz1 * (c34 * tmp1)- dt * tz1 * dz3;
                    a[j-1, i-1, 3, 2] = -dt * tz2 * (u[k-2, j+1, i+1, 2] * tmp1);
                    a[j-1, i-1, 4, 2] = 0.0d;

                    a[j-1, i-1, 0, 3] = -dt*tz2*(-pow2((u[k-2, j+1, i+1, 3]*tmp1))+0.50d*c2*((u[k-2, j+1, i+1, 1]*u[k-2, j+1, i+1, 1]+u[k-2, j+1, i+1, 2]*u[k-2, j+1, i+1, 2]+u[k-2, j+1, i+1, 3]*u[k-2, j+1, i+1, 3])*tmp2))-dt*tz1*(-r43*c34*tmp2*u[k-2, j+1, i+1, 3]);
                    a[j-1, i-1, 1, 3] = -dt * tz2* (-c2 * (u[k-2, j+1, i+1, 1] * tmp1));
                    a[j-1, i-1, 2, 3] = -dt * tz2* (-c2 * (u[k-2, j+1, i+1, 2] * tmp1));
                    a[j-1, i-1, 3, 3] = -dt * tz2 *(2.0d-c2)*(u[k-2, j+1, i+1, 3] * tmp1)- dt * tz1 * (r43 * c34 * tmp1)- dt * tz1 * dz4;
                    a[j-1, i-1, 4, 3] = -dt * tz2 * c2;

                    a[j-1, i-1, 0, 4] = -dt*tz2*((c2*(u[k-2, j+1, i+1, 1]*u[k-2, j+1, i+1, 1]+ u[k-2, j+1, i+1, 2] * u[k-2, j+1, i+1, 2]+ u[k-2, j+1, i+1, 3] * u[k-2, j+1, i+1, 3]) * tmp2- c1 * (u[k-2, j+1, i+1, 4] * tmp1))* (u[k-2, j+1, i+1, 3] * tmp1))- dt * tz1* (-(c34 - c1345) * tmp3 * (pow2(u[k-2, j+1, i+1, 1]))- (c34 - c1345) * tmp3 * (pow2(u[k-2, j+1, i+1, 2]))- (r43*c34 - c1345)* tmp3 *(pow2(u[k-2, j+1, i+1, 3]))- c1345 * tmp2 * u[k-2, j+1, i+1, 4]);
                    a[j-1, i-1, 1, 4] = -dt*tz2*(-c2*(u[k-2, j+1, i+1, 1]*u[k-2, j+1, i+1, 3])*tmp2)-dt*tz1*(c34-c1345)*tmp2*u[k-2, j+1, i+1, 1];
                    a[j-1, i-1, 2, 4] = -dt*tz2*(-c2*(u[k-2, j+1, i+1, 2]*u[k-2, j+1, i+1, 3])*tmp2)-dt*tz1*(c34-c1345)*tmp2*u[k-2, j+1, i+1, 2];
                    a[j-1, i-1, 3, 4] = -dt*tz2*(c1*(u[k-2, j+1, i+1, 4]*tmp1)-0.50d*c2*((u[k-2, j+1, i+1, 1]*u[k-2, j+1, i+1, 1]+u[k-2, j+1, i+1, 2]*u[k-2, j+1, i+1, 2]+3.0d*u[k-2, j+1, i+1, 3]*u[k-2, j+1, i+1, 3])*tmp2))-dt*tz1*(r43*c34-c1345)*tmp2*u[k-2, j+1, i+1, 3];
                    a[j-1, i-1, 4, 4] = -dt*tz2*(c1*(u[k-2, j+1, i+1, 3]*tmp1))-dt*tz1*c1345*tmp1-dt*tz1*dz5;
                    //---------------------------------------------------------------------
                    //   form the second block sub-diagonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k-1, j, i+1, 0];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    b[j-1, i-1, 0, 0] = -dt * ty1 * dy1;
                    b[j-1, i-1, 1, 0] =   0.0d;
                    b[j-1, i-1, 2, 0] = -dt * ty2;
                    b[j-1, i-1, 3, 0] =   0.0d;
                    b[j-1, i-1, 4, 0] =   0.0d;

                    b[j-1, i-1, 0, 1] = -dt*ty2*(-(u[k-1, j, i+1, 1]*u[k-1, j, i+1, 2])*tmp2)-dt*ty1*(-c34*tmp2*u[k-1, j, i+1, 1]);
                    b[j-1, i-1, 1, 1] = -dt*ty2*(u[k-1, j, i+1, 2]*tmp1)-dt*ty1*(c34*tmp1)-dt*ty1*dy2;
                    b[j-1, i-1, 2, 1] = -dt*ty2*(u[k-1, j, i+1, 1]*tmp1);
                    b[j-1, i-1, 3, 1] = 0.0d;
                    b[j-1, i-1, 4, 1] = 0.0d;

                    b[j-1, i-1, 0, 2] = -dt*ty2*(-pow2((u[k-1, j, i+1, 2]*tmp1))+0.50d*c2*((u[k-1, j, i+1, 1]* u[k-1, j, i+1, 1]+ u[k-1, j, i+1, 2]* u[k-1, j, i+1, 2]+ u[k-1, j, i+1, 3]* u[k-1, j, i+1, 3])* tmp2))- dt * ty1 *(-r43 * c34 * tmp2 * u[k-1, j, i+1, 2]);
                    b[j-1, i-1, 1, 2] = -dt*ty2*(-c2*(u[k-1, j, i+1, 1]*tmp1));
                    b[j-1, i-1, 2, 2] = -dt*ty2*((2.0d-c2)*(u[k-1, j, i+1, 2]*tmp1))-dt*ty1*(r43*c34*tmp1)-dt*ty1*dy3;
                    b[j-1, i-1, 3, 2] = -dt*ty2*(-c2*(u[k-1, j, i+1, 3]*tmp1));
                    b[j-1, i-1, 4, 2] = -dt*ty2*c2;

                    b[j-1, i-1, 0, 3] = -dt*ty2*(-(u[k-1, j, i+1, 2]*u[k-1, j, i+1, 3]) * tmp2)- dt * ty1 * (-c34 * tmp2 * u[k-1, j, i+1, 3]);
                    b[j-1, i-1, 1, 3] = 0.0d;
                    b[j-1, i-1, 2, 3] = -dt*ty2* (u[k-1, j, i+1, 3] * tmp1);
                    b[j-1, i-1, 3, 3] = -dt*ty2* (u[k-1, j, i+1, 2] * tmp1)- dt * ty1 * (c34 * tmp1)- dt * ty1 * dy4;
                    b[j-1, i-1, 4, 3] = 0.0d;

                    b[j-1, i-1, 0, 4] = -dt * ty2* ((c2 * (u[k-1, j, i+1, 1] *u[k-1, j, i+1, 1]+ u[k-1, j, i+1, 2] *u[k-1, j, i+1, 2]+ u[k-1, j, i+1, 3] * u[k-1, j, i+1, 3]) * tmp2- c1 * (u[k-1, j, i+1, 4] * tmp1))* (u[k-1, j, i+1, 2] * tmp1))- dt * ty1* (-(c34 - c1345)*tmp3*(pow2(u[k-1, j, i+1, 1]))- (r43*c34 - c1345)*tmp3*(pow2(u[k-1, j, i+1, 2]))- (c34 - c1345)*tmp3*(pow2(u[k-1, j, i+1, 3]))- c1345*tmp2*u[k-1, j, i+1, 4]);
                    b[j-1, i-1, 1, 4] = -dt * ty2* (-c2 * (u[k-1, j, i+1, 1]*u[k-1, j, i+1, 2]) * tmp2)- dt * ty1* (c34 - c1345) * tmp2 * u[k-1, j, i+1, 1];
                    b[j-1, i-1, 2, 4] = -dt * ty2* (c1 *(u[k-1, j, i+1, 4] * tmp1)- 0.50d * c2 *((u[k-1, j, i+1, 1]*u[k-1, j, i+1, 1]+ 3.0d * u[k-1, j, i+1, 2]*u[k-1, j, i+1, 2]+u[k-1, j, i+1, 3]*u[k-1, j, i+1, 3]) * tmp2))- dt * ty1* (r43*c34 - c1345) * tmp2 * u[k-1, j, i+1, 2];
                    b[j-1, i-1, 3, 4] = -dt * ty2* (-c2 *(u[k-1, j, i+1, 2]*u[k-1, j, i+1, 3])* tmp2)- dt * ty1 * (c34 - c1345)* tmp2 *u[k-1, j, i+1, 3];
                    b[j-1, i-1, 4, 4] = -dt * ty2* (c1 *(u[k-1, j, i+1, 2] * tmp1))- dt * ty1 * c1345 * tmp1- dt * ty1 * dy5;
                    //---------------------------------------------------------------------
                    //   form the third block sub-diagonal
                    //---------------------------------------------------------------------                      
                    tmp1 = 1.0d / u[k-1, j+1, i, 0];//tmp1 = 1.0d / u[1,i-1,j,k];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    c[j-1, i-1, 0, 0] = -dt * tx1 * dx1;
                    c[j-1, i-1, 1, 0] = -dt * tx2;
                    c[j-1, i-1, 2, 0] =   0.0d;
                    c[j-1, i-1, 3, 0] =   0.0d;
                    c[j-1, i-1, 4, 0] =   0.0d;

                    c[j-1, i-1, 0, 1] = -dt*tx2*(-pow2((u[k-1, j+1, i, 1]*tmp1))+ c2 * 0.50d * (u[k-1, j+1, i, 1] * u[k-1, j+1, i, 1]+ u[k-1, j+1, i, 2] * u[k-1, j+1, i, 2]+ u[k-1, j+1, i, 3] * u[k-1, j+1, i, 3]) * tmp2)- dt * tx1 * (-r43 * c34 * tmp2 * u[k-1, j+1, i, 1]);
                    c[j-1, i-1, 1, 1] = -dt*tx2*((2.0d-c2)*(u[k-1, j+1, i, 1] * tmp1))- dt * tx1 * (r43 * c34 * tmp1)- dt * tx1 * dx2;
                    c[j-1, i-1, 2, 1] = -dt*tx2*(-c2*(u[k-1, j+1, i, 2] * tmp1));
                    c[j-1, i-1, 3, 1] = -dt*tx2*(-c2*(u[k-1, j+1, i, 3] * tmp1));
                    c[j-1, i-1, 4, 1] = -dt*tx2*c2;

                    c[j-1, i-1, 0, 2] = -dt*tx2*(-(u[k-1, j+1, i, 1] * u[k-1, j+1, i, 2]) * tmp2)- dt * tx1 * (-c34 * tmp2 * u[k-1, j+1, i, 2]);
                    c[j-1, i-1, 1, 2] = -dt*tx2*(u[k-1, j+1, i, 2] * tmp1);
                    c[j-1, i-1, 2, 2] = -dt*tx2*(u[k-1, j+1, i, 1] * tmp1)- dt * tx1 * (c34 * tmp1)- dt * tx1 * dx3;
                    c[j-1, i-1, 3, 2] = 0.0d;
                    c[j-1, i-1, 4, 2] = 0.0d;

                    c[j-1, i-1, 0, 3] = -dt * tx2*(-(u[k-1, j+1, i, 1]*u[k-1, j+1, i, 3]) * tmp2)- dt * tx1 * (-c34 * tmp2 * u[k-1, j+1, i, 3]);
                    c[j-1, i-1, 1, 3] = -dt * tx2 * (u[k-1, j+1, i, 3] * tmp1);
                    c[j-1, i-1, 2, 3] = 0.0d;
                    c[j-1, i-1, 3, 3] = -dt * tx2 * (u[k-1, j+1, i, 1] * tmp1)- dt * tx1 * (c34 * tmp1)- dt * tx1 * dx4;
                    c[j-1, i-1, 4, 3] = 0.0d;

                    c[j-1, i-1, 0, 4] = -dt * tx2* ((c2 * (u[k-1, j+1, i, 1] * u[k-1, j+1, i, 1]+ u[k-1, j+1, i, 2] * u[k-1, j+1, i, 2]+ u[k-1, j+1, i, 3] * u[k-1, j+1, i, 3]) * tmp2- c1 * (u[k-1, j+1, i, 4] * tmp1))* (u[k-1, j+1, i, 1] * tmp1))- dt * tx1* (-(r43*c34 - c1345) * tmp3 * (pow2(u[k-1, j+1, i, 1]))- (c34 - c1345) * tmp3 * (pow2(u[k-1, j+1, i, 2]))- (c34 - c1345) * tmp3 * (pow2(u[k-1, j+1, i, 3]))- c1345 * tmp2 * u[k-1, j+1, i, 4]);
                    c[j-1, i-1, 1, 4] = -dt * tx2* (c1 * (u[k-1, j+1, i, 4] * tmp1)- 0.50d * c2* ((3.0d*u[k-1, j+1, i, 1]*u[k-1, j+1, i, 1]+u[k-1, j+1, i, 2]*u[k-1, j+1, i, 2]+u[k-1, j+1, i, 3]*u[k-1, j+1, i, 3]) * tmp2))- dt * tx1* (r43*c34 - c1345) * tmp2 * u[k-1, j+1, i, 1];
                    c[j-1, i-1, 2, 4] = -dt * tx2* (-c2 * (u[k-1, j+1, i, 2]*u[k-1, j+1, i, 1]) * tmp2)- dt * tx1* (c34 - c1345) * tmp2 * u[k-1, j+1, i, 2];
                    c[j-1, i-1, 3, 4] = -dt * tx2* (-c2 * (u[k-1, j+1, i, 3]*u[k-1, j+1, i, 1]) * tmp2)- dt * tx1* (c34 - c1345) * tmp2 * u[k-1, j+1, i, 3];
                    c[j-1, i-1, 4, 4] = -dt * tx2* (c1 * (u[k-1, j+1, i, 1] * tmp1))- dt * tx1 * c1345 * tmp1- dt * tx1 * dx5;
                }
            }
        }
        //end jacld.f
        // blts.f
        public void blts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k, double omega, double[, , ,] v, 
                         double[, , ,] zld, double[, , ,] yld, double[, , ,] xld, double[, , ,] d, 
                         int ist, int iend, int jst, int jend, int nx0, int ny0, int ipt, int jpt) {
            //---------------------------------------------------------------------
            //   compute the regular-sparse, block lower triangular solution:
            //                     v <-- [ L-inv ] * v
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //int ldmx, ldmy, ldmz,nx, ny, nz,k;
            //double  omega;
            //double  v[ 5, -1:ldmx+2, -1:ldmy+2, *],ldz[ 5, 5, ldmx, ldmy],ldy[ 5, 5, ldmx, ldmy],ldx[ 5, 5, ldmx, ldmy],d[ 5, 5, ldmx, ldmy];
            //int ist, iend,jst, jend,nx0, ny0,ipt, jpt;
            //---------------------------------------------------------------------
            //  local variables
            //---------------------------------------------------------------------
            int i, j, m, iex;
            double  tmp, tmp1;
            double[,] tmat = new double[5, 5]; //tmat[5,5]
            //---------------------------------------------------------------------
            //   receive data from north and west
            //---------------------------------------------------------------------
            //Debug
            //if ((isiz1 + 2) != (ldmx + 2) || (isiz2 + 2) != (ldmy + 2)) {
            //    throw new ArgumentException("Look this code: vetor v");
            //}//end Debug

            iex = 0;
            exchange_1(v, k, iex);
            for(j = jst; j<= jend; j++) {
                for(i = ist; i<= iend; i++) {
                    for(m = 1; m<= 5; m++) {
                        v[k-1, j+1, i+1, m -1] =  v[k-1, j+1, i+1, m -1]
             - omega * (  zld[j-1, i-1, 0, m-1] * v[k-2, j+1, i+1, 0]
                        + zld[j-1, i-1, 1, m-1] * v[k-2, j+1, i+1, 1]
                        + zld[j-1, i-1, 2, m-1] * v[k-2, j+1, i+1, 2]
                        + zld[j-1, i-1, 3, m-1] * v[k-2, j+1, i+1, 3]
                        + zld[j-1, i-1, 4, m-1] * v[k-2, j+1, i+1, 4]);
                    }
                }
            }
            for(j=jst; j<=jend; j++) {
                for(i = ist; i<= iend; i++) {
                    for(m = 1; m<= 5; m++) {
                        v[k-1, j+1, i+1, m-1] =  v[k-1, j+1, i+1, m-1]
                        - omega * ( yld[j-1, i-1, 0, m-1] * v[k-1, j, i+1, 0]
                                  + xld[j-1, i-1, 0, m-1] * v[k-1, j+1, i, 0]
                                  + yld[j-1, i-1, 1, m-1] * v[k-1, j, i+1, 1]
                                  + xld[j-1, i-1, 1, m-1] * v[k-1, j+1, i, 1]
                                  + yld[j-1, i-1, 2, m-1] * v[k-1, j, i+1, 2]
                                  + xld[j-1, i-1, 2, m-1] * v[k-1, j+1, i, 2]
                                  + yld[j-1, i-1, 3, m-1] * v[k-1, j, i+1, 3]
                                  + xld[j-1, i-1, 3, m-1] * v[k-1, j+1, i, 3]
                                  + yld[j-1, i-1, 4, m-1] * v[k-1, j, i+1, 4]
                                  + xld[j-1, i-1, 4, m-1] * v[k-1, j+1, i, 4]);
                    }
                    //---------------------------------------------------------------------
                    //   diagonal block inversion
                    //
                    //   forward elimination
                    //---------------------------------------------------------------------
                    for(m = 0; m< 5; m++) {
                        tmat[0, m] = d[j-1, i-1, 0, m];
                        tmat[1, m] = d[j-1, i-1, 1, m];
                        tmat[2, m] = d[j-1, i-1, 2, m];
                        tmat[3, m] = d[j-1, i-1, 3, m];
                        tmat[4, m] = d[j-1, i-1, 4, m];
                    }
                    tmp1 = 1.0d /tmat[0, 0];
                    tmp = tmp1 * tmat[0, 1];
                    tmat[1, 1] =  tmat[1, 1] - tmp * tmat[1, 0];
                    tmat[2, 1] =  tmat[2, 1] - tmp * tmat[2, 0];
                    tmat[3, 1] =  tmat[3, 1] - tmp * tmat[3, 0];
                    tmat[4, 1] =  tmat[4, 1] - tmp * tmat[4, 0];
                    v[k-1, j+1, i+1, 1] = v[k-1, j+1, i+1, 1] - v[k-1, j+1, i+1, 0] * tmp;

                    tmp = tmp1 * tmat[0, 2];
                    tmat[1, 2] =  tmat[1, 2] - tmp * tmat[1, 0];
                    tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 0];
                    tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 0];
                    tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 0];
                    v[k-1, j+1, i+1, 2] = v[k-1, j+1, i+1, 2] - v[k-1, j+1, i+1, 0] * tmp;

                    tmp = tmp1 * tmat[0, 3];
                    tmat[1, 3] =  tmat[1, 3] - tmp * tmat[1, 0];
                    tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 0];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 0];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 0];
                    v[k-1, j+1, i+1, 3] = v[k-1, j+1, i+1, 3] - v[k-1, j+1, i+1, 0] * tmp;

                    tmp = tmp1 * tmat[0, 4];
                    tmat[1, 4] =  tmat[1, 4] - tmp * tmat[1, 0];
                    tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 0];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 0];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 0];
                    v[k-1, j+1, i+1, 4] = v[k-1, j+1, i+1, 4] - v[k-1, j+1, i+1, 0] * tmp;

                    tmp1 = 1.0d /tmat[1, 1];
                    tmp = tmp1 * tmat[1, 2];
                    tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 1];
                    tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 1];
                    tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 1];
                    v[k-1, j+1, i+1, 2] = v[k-1, j+1, i+1, 2] - v[k-1, j+1, i+1, 1] * tmp;

                    tmp = tmp1 * tmat[1, 3];
                    tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 1];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 1];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 1];
                    v[k-1, j+1, i+1, 3] = v[k-1, j+1, i+1, 3] - v[k-1, j+1, i+1, 1] * tmp;

                    tmp = tmp1 * tmat[1, 4];
                    tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 1];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 1];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 1];
                    v[k-1, j+1, i+1, 4] = v[k-1, j+1, i+1, 4] - v[k-1, j+1, i+1, 1] * tmp;

                    tmp1 = 1.0d /tmat[2, 2];
                    tmp = tmp1 * tmat[2, 3];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 2];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 2];
                    v[k-1, j+1, i+1, 3] = v[k-1, j+1, i+1, 3] - v[k-1, j+1, i+1, 2] * tmp;

                    tmp = tmp1 * tmat[2, 4];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 2];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 2];
                    v[k-1, j+1, i+1, 4] = v[k-1, j+1, i+1, 4] - v[k-1, j+1, i+1, 2] * tmp;

                    tmp1 = 1.0d /tmat[3, 3];
                    tmp = tmp1 * tmat[3, 4];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 3];
                    v[k-1, j+1, i+1, 4] = v[k-1, j+1, i+1, 4] - v[k-1, j+1, i+1, 3] * tmp;

                    //---------------------------------------------------------------------
                    //   back substitution
                    //---------------------------------------------------------------------

                    v[k-1, j+1, i+1, 4] = v[k-1, j+1, i+1, 4]/ tmat[4, 4];
                    v[k-1, j+1, i+1, 3] = v[k-1, j+1, i+1, 3]- tmat[4, 3] * v[k-1, j+1, i+1, 4];
                    v[k-1, j+1, i+1, 3] = v[k-1, j+1, i+1, 3]/ tmat[3, 3];
                    v[k-1, j+1, i+1, 2] = v[k-1, j+1, i+1, 2] -tmat[3, 2] * v[k-1, j+1, i+1, 3] - tmat[4, 2] * 
                                                                       v[k-1, j+1, i+1, 4];
                    v[k-1, j+1, i+1, 2] = v[k-1, j+1, i+1, 2] /tmat[2, 2];
                    v[k-1, j+1, i+1, 1] = v[k-1, j+1, i+1, 1]- tmat[2, 1] * v[k-1, j+1, i+1, 2]-tmat[3, 1]*
                                                                       v[k-1, j+1, i+1, 3]-tmat[4, 1]*
                                                                       v[k-1, j+1, i+1, 4];
                    v[k-1, j+1, i+1, 1] = v[k-1, j+1, i+1, 1] /tmat[1, 1];
                    v[k-1, j+1, i+1, 0] = v[k-1, j+1, i+1, 0] -tmat[1, 0] * v[k-1, j+1, i+1, 1]-
                                                           tmat[2, 0] * v[k-1, j+1, i+1, 2]-
                                                           tmat[3, 0] * v[k-1, j+1, i+1, 3]-
                                                           tmat[4, 0] * v[k-1, j+1, i+1, 4];
                    v[k-1, j+1, i+1, 0] = v[k-1, j+1, i+1, 0] /tmat[0, 0];

                }
            }
            //---------------------------------------------------------------------
            //   send data to east and south
            //---------------------------------------------------------------------
            iex = 2;
            exchange_1(v, k, iex);
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
                    worldcomm.Receive<double>(north, from_n, ref dum1);
                    for(j=jst; j<=jend; j++) {
                        g[k-1, j+1, 1, 0] = dum1[0+idx];//g[1,0,j,k] = dum1[1,j];
                        g[k-1, j+1, 1, 1] = dum1[1+idx];//g[2,0,j,k] = dum1[2,j];
                        g[k-1, j+1, 1, 2] = dum1[2+idx];//g[3,0,j,k] = dum1[3,j];
                        g[k-1, j+1, 1, 3] = dum1[3+idx];//g[4,0,j,k] = dum1[4,j];
                        g[k-1, j+1, 1, 4] = dum1[4+idx];//g[5,0,j,k] = dum1[5,j];
                        idx = idx + 5;
                    }
                }
                if(west != -1) {
                    double[] dum1 = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    //call MPI_RECV[ dum1[1,ist],5*[iend-ist+1],dp_type,west,from_w,MPI_COMM_WORLD,status,IERROR ]
                    worldcomm.Receive<double>(west, from_w, ref dum1);
                    for(i=ist; i<=iend; i++) {
                        g[k-1, 1, i+1, 0] = dum1[0+idx];//g[1,i,0,k] = dum1[1,i];
                        g[k-1, 1, i+1, 1] = dum1[1+idx];//g[2,i,0,k] = dum1[2,i];
                        g[k-1, 1, i+1, 2] = dum1[2+idx];//g[3,i,0,k] = dum1[3,i];
                        g[k-1, 1, i+1, 3] = dum1[3+idx];//g[4,i,0,k] = dum1[4,i];
                        g[k-1, 1, i+1, 4] = dum1[4+idx];//g[5,i,0,k] = dum1[5,i];
                        idx = idx + 5;
                    }
                }
            }
            else if(iex == 1) {
                if(south != -1) {
                    double[] dum1 = new double[(5*(jend-jst+1))];
                    int idx = 0;
                    //call MPI_RECV[ dum1[1,jst],5*[jend-jst+1],dp_type,south,from_s,MPI_COMM_WORLD,status,IERROR ]
                    worldcomm.Receive<double>(south, from_s, ref dum1);
                    for(j=jst; j<=jend; j++) {
                        g[k-1, j+1, nx+2, 0] = dum1[0+idx];//g[1,nx+1,j,k] = dum1[1,j];
                        g[k-1, j+1, nx+2, 1] = dum1[1+idx];//g[2,nx+1,j,k] = dum1[2,j];
                        g[k-1, j+1, nx+2, 2] = dum1[2+idx];//g[3,nx+1,j,k] = dum1[3,j];
                        g[k-1, j+1, nx+2, 3] = dum1[3+idx];//g[4,nx+1,j,k] = dum1[4,j];
                        g[k-1, j+1, nx+2, 4] = dum1[4+idx];//g[5,nx+1,j,k] = dum1[5,j];
                        idx = idx + 5;
                    }
                }
                if(east != -1) {
                    double[] dum1 = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    //call MPI_RECV[ dum1[1,ist],5*[iend-ist+1],dp_type,east,from_e,MPI_COMM_WORLD,status,IERROR ]
                    worldcomm.Receive<double>(east, from_e, ref dum1);
                    for(i=ist; i<=iend; i++) {
                        g[k-1, ny+2, i+1, 0] = dum1[0+idx];//g[1,i,ny+1,k] = dum1[1,i];
                        g[k-1, ny+2, i+1, 1] = dum1[1+idx];//g[2,i,ny+1,k] = dum1[2,i];
                        g[k-1, ny+2, i+1, 2] = dum1[2+idx];//g[3,i,ny+1,k] = dum1[3,i];
                        g[k-1, ny+2, i+1, 3] = dum1[3+idx];//g[4,i,ny+1,k] = dum1[4,i];
                        g[k-1, ny+2, i+1, 4] = dum1[4+idx];//g[5,i,ny+1,k] = dum1[5,i];
                        idx = idx + 5;
                    }
                }
            }
            else if(iex == 2) {
                if(south != -1) {
                    double[] dum = new double[5*(jend-jst+1)];
                    int idx = 0;
                    for(j=jst; j<=jend; j++) {
                        dum[0+idx] = g[k-1, j+1, nx+1, 0];//dum[1,j] = g[1,nx,j,k];
                        dum[1+idx] = g[k-1, j+1, nx+1, 1];//dum[2,j] = g[2,nx,j,k];
                        dum[2+idx] = g[k-1, j+1, nx+1, 2];//dum[3,j] = g[3,nx,j,k];
                        dum[3+idx] = g[k-1, j+1, nx+1, 3];//dum[4,j] = g[4,nx,j,k];
                        dum[4+idx] = g[k-1, j+1, nx+1, 4];//dum[5,j] = g[5,nx,j,k];
                        idx = idx + 5;
                    }
                    //call MPI_SEND[ dum[1,jst], 5*[jend-jst+1], dp_type, south, from_n, MPI_COMM_WORLD, IERROR ]
                    worldcomm.Send<double>(dum, south, from_n);
                }
                if(east != -1) {
                    double[] dum = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    for(i=ist; i<=iend; i++) {
                        dum[0+idx] = g[k-1, ny+1, i+1, 0];//dum[1,i] = g[1,i,ny,k];
                        dum[1+idx] = g[k-1, ny+1, i+1, 1];//dum[2,i] = g[2,i,ny,k];
                        dum[2+idx] = g[k-1, ny+1, i+1, 2];//dum[3,i] = g[3,i,ny,k];
                        dum[3+idx] = g[k-1, ny+1, i+1, 3];//dum[4,i] = g[4,i,ny,k];
                        dum[4+idx] = g[k-1, ny+1, i+1, 4];//dum[5,i] = g[5,i,ny,k];
                        idx = idx + 5;
                    }
                    //call MPI_SEND[ dum[1,ist], 5*[iend-ist+1], dp_type, east, from_w, MPI_COMM_WORLD, IERROR ]
                    worldcomm.Send<double>(dum, east, from_w);
                }
            }
            else {
                if(north != -1) {
                    double[] dum = new double[(5*(jend-jst+1))];
                    int idx = 0;
                    for(j=jst; j<=jend; j++) {
                        dum[0+idx] = g[k-1, j+1, 2, 0];//dum[1,j] = g[1,1,j,k];
                        dum[1+idx] = g[k-1, j+1, 2, 1];//dum[2,j] = g[2,1,j,k];
                        dum[2+idx] = g[k-1, j+1, 2, 2];//dum[3,j] = g[3,1,j,k];
                        dum[3+idx] = g[k-1, j+1, 2, 3];//dum[4,j] = g[4,1,j,k];
                        dum[4+idx] = g[k-1, j+1, 2, 4];//dum[5,j] = g[5,1,j,k];
                        idx = idx + 5;
                    }
                    //call MPI_SEND[ dum[1,jst], 5*[jend-jst+1], dp_type, north, from_s, MPI_COMM_WORLD, IERROR ]
                    worldcomm.Send<double>(dum, north, from_s);
                }
                if(west != -1) {
                    double[] dum = new double[(5*(iend-ist+1))];
                    int idx = 0;
                    for(i=ist; i<=iend; i++) {
                        dum[0+idx] = g[k-1, 2, i+1, 0];//dum[1,i] = g[1,i,1,k];
                        dum[1+idx] = g[k-1, 2, i+1, 1];//dum[2,i] = g[2,i,1,k];
                        dum[2+idx] = g[k-1, 2, i+1, 2];//dum[3,i] = g[3,i,1,k];
                        dum[3+idx] = g[k-1, 2, i+1, 3];//dum[4,i] = g[4,i,1,k];
                        dum[4+idx] = g[k-1, 2, i+1, 4];//dum[5,i] = g[5,i,1,k];
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

            r43 = (4.0d / 3.0d);
            c1345 = c1 * c3 * c4 * c5;
            c34 = c3 * c4;

            for(j = jst; j<= jend; j++) {
                for(i = ist; i<= iend; i++) {
                    //---------------------------------------------------------------------
                    //   form the block daigonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k-1, j+1, i+1, 0];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    d[j-1, i-1, 0, 0] =  1.0d+ dt * 2.0d * (tx1 * dx1+ ty1 * dy1+ tz1 * dz1);
                    d[j-1, i-1, 1, 0] =  0.0d;
                    d[j-1, i-1, 2, 0] =  0.0d;
                    d[j-1, i-1, 3, 0] =  0.0d;
                    d[j-1, i-1, 4, 0] =  0.0d;

                    d[j-1, i-1, 0, 1] =  dt*2.0d*(tx1*(-r43*c34*tmp2*u[k-1, j+1, i+1, 1])+ty1*(-c34*tmp2*u[k-1, j+1, i+1, 1])+tz1*(-c34*tmp2*u[k-1, j+1, i+1, 1]));
                    d[j-1, i-1, 1, 1] =  1.0d+dt*2.0d*(tx1*r43*c34*tmp1+ty1*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx2+ty1*dy2+tz1*dz2);
                    d[j-1, i-1, 2, 1] = 0.0d;
                    d[j-1, i-1, 3, 1] = 0.0d;
                    d[j-1, i-1, 4, 1] = 0.0d;

                    d[j-1, i-1, 0, 2] = dt * 2.0d*(tx1*(-c34*tmp2*u[k-1, j+1, i+1, 2])+ty1*(-r43*c34*tmp2*u[k-1, j+1, i+1, 2])+tz1*(-c34*tmp2 * u[k-1, j+1, i+1, 2]));
                    d[j-1, i-1, 1, 2] = 0.0d;
                    d[j-1, i-1, 2, 2] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*r43*c34*tmp1+tz1*c34*tmp1)+dt*2.0d*(tx1*dx3+ty1*dy3+tz1*dz3);
                    d[j-1, i-1, 3, 2] = 0.0d;
                    d[j-1, i-1, 4, 2] = 0.0d;

                    d[j-1, i-1, 0, 3] = dt * 2.0d*(tx1*(-c34*tmp2*u[k-1, j+1, i+1, 3])+ty1*(-c34*tmp2*u[k-1, j+1, i+1, 3])+tz1*(-r43*c34*tmp2*u[k-1, j+1, i+1, 3]));
                    d[j-1, i-1, 1, 3] = 0.0d;
                    d[j-1, i-1, 2, 3] = 0.0d;
                    d[j-1, i-1, 3, 3] = 1.0d+dt*2.0d*(tx1*c34*tmp1+ty1*c34*tmp1+tz1*r43*c34*tmp1)+dt*2.0d*(tx1*dx4+ty1*dy4+tz1*dz4);
                    d[j-1, i-1, 4, 3] = 0.0d;

                    d[j-1, i-1, 0, 4] = dt * 2.0d*(tx1*  (-(r43*c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 1]))- (c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 2]))- (c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 3]))- (c1345) * tmp2 * u[k-1, j+1, i+1, 4])+ ty1 * (-(c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 1]))- (r43*c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 2]))- (c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 3]))- (c1345) * tmp2 * u[k-1, j+1, i+1, 4])+ tz1 * (-(c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 1]))- (c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 2]))- (r43*c34 - c1345)*tmp3*(pow2(u[k-1, j+1, i+1, 3]))- (c1345) * tmp2 * u[k-1, j+1, i+1, 4]));
                    d[j-1, i-1, 1, 4] = dt * 2.0d  * (tx1 * (r43*c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 1]+ ty1 * (c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 1]+ tz1 * (c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 1]);
                    d[j-1, i-1, 2, 4] = dt * 2.0d * (tx1 * (c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 2]+ ty1 * (r43*c34 -c1345) * tmp2 * u[k-1, j+1, i+1, 2]+ tz1 * (c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 2]);
                    d[j-1, i-1, 3, 4] = dt * 2.0d* (tx1 * (c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 3]+ ty1 * (c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 3]+ tz1 * (r43*c34 - c1345) * tmp2 * u[k-1, j+1, i+1, 3]);
                    d[j-1, i-1, 4, 4] = 1.0d+dt*2.0d*(tx1*c1345*tmp1+ty1*c1345*tmp1+tz1*c1345*tmp1)+dt*2.0d*(tx1*dx5+ty1*dy5+tz1*dz5);
                    //---------------------------------------------------------------------
                    //   form the first block sub-diagonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k-1, j+1, i+2, 0]; //tmp1 = 1.0d / u[1,i+1,j,k];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    a[j-1, i-1, 0, 0] = -dt * tx1 * dx1;
                    a[j-1, i-1, 1, 0] =   dt * tx2;
                    a[j-1, i-1, 2, 0] =   0.0d;
                    a[j-1, i-1, 3, 0] =   0.0d;
                    a[j-1, i-1, 4, 0] =   0.0d;

                    a[j-1, i-1, 0, 1] =  dt * tx2*(-pow2((u[k-1, j+1, i+2, 1]*tmp1))+ c2 * 0.50d * (u[k-1, j+1, i+2, 1] * u[k-1, j+1, i+2, 1]+ u[k-1, j+1, i+2, 2] * u[k-1, j+1, i+2, 2]+ u[k-1, j+1, i+2, 3] * u[k-1, j+1, i+2, 3]) * tmp2)-dt*tx1*(-r43*c34*tmp2*u[k-1, j+1, i+2, 1]);
                    a[j-1, i-1, 1, 1] =  dt*tx2*((2.0d-c2)*(u[k-1, j+1, i+2, 1]*tmp1))-dt*tx1*(r43*c34*tmp1)-dt*tx1*dx2;
                    a[j-1, i-1, 2, 1] =  dt* tx2*(-c2*(u[k-1, j+1, i+2, 2] * tmp1));
                    a[j-1, i-1, 3, 1] =  dt* tx2*(-c2*(u[k-1, j+1, i+2, 3] * tmp1));
                    a[j-1, i-1, 4, 1] =  dt * tx2 * c2;

                    a[j-1, i-1, 0, 2] =  dt * tx2*(-(u[k-1, j+1, i+2, 1] * u[k-1, j+1, i+2, 2]) * tmp2)- dt * tx1 * (-c34 * tmp2 * u[k-1, j+1, i+2, 2]);
                    a[j-1, i-1, 1, 2] =  dt * tx2 * (u[k-1, j+1, i+2, 2] * tmp1);
                    a[j-1, i-1, 2, 2] =  dt * tx2 * (u[k-1, j+1, i+2, 1] * tmp1)- dt * tx1 * (c34 * tmp1)- dt * tx1 * dx3;
                    a[j-1, i-1, 3, 2] = 0.0d;
                    a[j-1, i-1, 4, 2] = 0.0d;

                    a[j-1, i-1, 0, 3] = dt * tx2*(-(u[k-1, j+1, i+2, 1]*u[k-1, j+1, i+2, 3]) * tmp2)- dt * tx1 * (-c34 * tmp2 * u[k-1, j+1, i+2, 3]);
                    a[j-1, i-1, 1, 3] = dt * tx2 * (u[k-1, j+1, i+2, 3] * tmp1);
                    a[j-1, i-1, 2, 3] = 0.0d;
                    a[j-1, i-1, 3, 3] = dt * tx2 * (u[k-1, j+1, i+2, 1] * tmp1)- dt * tx1 * (c34 * tmp1)- dt * tx1 * dx4;
                    a[j-1, i-1, 4, 3] = 0.0d;

                    a[j-1, i-1, 0, 4] = dt*tx2*((c2 * (u[k-1, j+1, i+2, 1] * u[k-1, j+1, i+2, 1]+ u[k-1, j+1, i+2, 2] * u[k-1, j+1, i+2, 2]+ u[k-1, j+1, i+2, 3] * u[k-1, j+1, i+2, 3]) * tmp2- c1 * (u[k-1, j+1, i+2, 4] * tmp1))* (u[k-1, j+1, i+2, 1] * tmp1))- dt * tx1* (-(r43*c34 - c1345) * tmp3 * (pow2(u[k-1, j+1, i+2, 1]))- (c34 - c1345) * tmp3 * (pow2(u[k-1, j+1, i+2, 2]))- (c34 - c1345) * tmp3 * (pow2(u[k-1, j+1, i+2, 3]))- c1345 * tmp2 * u[k-1, j+1, i+2, 4]);
                    a[j-1, i-1, 1, 4] = dt * tx2*(c1 * (u[k-1, j+1, i+2, 4]*tmp1)-0.50d*c2*((3.0d*u[k-1, j+1, i+2, 1]*u[k-1, j+1, i+2, 1]+ u[k-1, j+1, i+2, 2]*u[k-1, j+1, i+2, 2]+ u[k-1, j+1, i+2, 3]*u[k-1, j+1, i+2, 3]) * tmp2))- dt * tx1*(r43*c34-c1345)*tmp2*u[k-1, j+1, i+2, 1];
                    a[j-1, i-1, 2, 4] = dt*tx2*(-c2*(u[k-1, j+1, i+2, 2]*u[k-1, j+1, i+2, 1])*tmp2)-dt*tx1*(c34-c1345)*tmp2* u[k-1, j+1, i+2, 2];
                    a[j-1, i-1, 3, 4] = dt*tx2*(-c2*(u[k-1, j+1, i+2, 3]*u[k-1, j+1, i+2, 1]) * tmp2)- dt * tx1* (c34 - c1345) * tmp2 * u[k-1, j+1, i+2, 3];
                    a[j-1, i-1, 4, 4] = dt*tx2*(c1*(u[k-1, j+1, i+2, 1] * tmp1))- dt * tx1 * c1345 * tmp1- dt * tx1 * dx5;
                    //---------------------------------------------------------------------
                    //   form the second block sub-diagonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k-1, j+2, i+1, 0]; //tmp1 = 1.0d / u[1,i,j+1,k];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    b[j-1, i-1, 0, 0] = -dt * ty1 * dy1;
                    b[j-1, i-1, 1, 0] =   0.0d;
                    b[j-1, i-1, 2, 0] =  dt * ty2;
                    b[j-1, i-1, 3, 0] =   0.0d;
                    b[j-1, i-1, 4, 0] =   0.0d;

                    b[j-1, i-1, 0, 1] = dt*ty2*(-(u[k-1, j+2, i+1, 1]*u[k-1, j+2, i+1, 2]) * tmp2)- dt * ty1 * (-c34 * tmp2 * u[k-1, j+2, i+1, 1]);
                    b[j-1, i-1, 1, 1] =  dt*ty2*(u[k-1, j+2, i+1, 2] * tmp1)- dt * ty1 * (c34 * tmp1)- dt * ty1 * dy2;
                    b[j-1, i-1, 2, 1] =  dt*ty2*(u[k-1, j+2, i+1, 1] * tmp1);
                    b[j-1, i-1, 3, 1] = 0.0d;
                    b[j-1, i-1, 4, 1] = 0.0d;

                    b[j-1, i-1, 0, 2] =  dt * ty2* (-pow2((u[k-1, j+2, i+1, 2]*tmp1))+0.50d*c2*((u[k-1, j+2, i+1, 1] * u[k-1, j+2, i+1, 1]+ u[k-1, j+2, i+1, 2] * u[k-1, j+2, i+1, 2]+ u[k-1, j+2, i+1, 3] * u[k-1, j+2, i+1, 3])* tmp2))-dt*ty1*(-r43*c34*tmp2* u[k-1, j+2, i+1, 2]);
                    b[j-1, i-1, 1, 2] =  dt * ty2*(-c2*(u[k-1, j+2, i+1, 1] * tmp1));
                    b[j-1, i-1, 2, 2] =  dt*ty2*((2.0d-c2)*(u[k-1, j+2, i+1, 2] * tmp1))- dt * ty1 * (r43 * c34 * tmp1)- dt * ty1 * dy3;
                    b[j-1, i-1, 3, 2] =  dt*ty2*(-c2*(u[k-1, j+2, i+1, 3] * tmp1));
                    b[j-1, i-1, 4, 2] =  dt * ty2 * c2;

                    b[j-1, i-1, 0, 3] =  dt * ty2* (-(u[k-1, j+2, i+1, 2]*u[k-1, j+2, i+1, 3]) * tmp2)- dt * ty1 * (-c34 * tmp2 * u[k-1, j+2, i+1, 3]);
                    b[j-1, i-1, 1, 3] = 0.0d;
                    b[j-1, i-1, 2, 3] =  dt * ty2 *    (u[k-1, j+2, i+1, 3] * tmp1);
                    b[j-1, i-1, 3, 3] =  dt * ty2 *    (u[k-1, j+2, i+1, 2] * tmp1)- dt * ty1 * (c34 * tmp1)- dt * ty1 * dy4;
                    b[j-1, i-1, 4, 3] = 0.0d;

                    b[j-1, i-1, 0, 4] =  dt * ty2* ((c2*(u[k-1, j+2, i+1, 1] * u[k-1, j+2, i+1, 1]+ u[k-1, j+2, i+1, 2] * u[k-1, j+2, i+1, 2]+ u[k-1, j+2, i+1, 3] *u[k-1, j+2, i+1, 3]) * tmp2- c1 * (u[k-1, j+2, i+1, 4]*tmp1))*(u[k-1, j+2, i+1, 2] * tmp1))- dt * ty1* (-(c34 - c1345)*tmp3*   (pow2(u[k-1, j+2, i+1, 1]))- (r43*c34 - c1345)*tmp3*   (pow2(u[k-1, j+2, i+1, 2]))- (c34 - c1345)*tmp3*   (pow2(u[k-1, j+2, i+1, 3]))- c1345*tmp2*u[k-1, j+2, i+1, 4]);
                    b[j-1, i-1, 1, 4] =  dt * ty2* (-c2*(u[k-1, j+2, i+1, 1]*u[k-1, j+2, i+1, 2]) * tmp2)- dt * ty1* (c34 - c1345) * tmp2 * u[k-1, j+2, i+1, 1];
                    b[j-1, i-1, 2, 4] =  dt * ty2* (c1 * (u[k-1, j+2, i+1, 4] * tmp1)- 0.50d * c2 * ((u[k-1, j+2, i+1, 1]*u[k-1, j+2, i+1, 1]+ 3.0d * u[k-1, j+2, i+1, 2]*u[k-1, j+2, i+1, 2]+ u[k-1, j+2, i+1, 3]*u[k-1, j+2, i+1, 3]) * tmp2))- dt * ty1* (r43*c34 - c1345) * tmp2 * u[k-1, j+2, i+1, 2];
                    b[j-1, i-1, 3, 4] =  dt * ty2* (-c2*(u[k-1, j+2, i+1, 2]*u[k-1, j+2, i+1, 3]) * tmp2)- dt * ty1 * (c34 - c1345) * tmp2 * u[k-1, j+2, i+1, 3];
                    b[j-1, i-1, 4, 4] =  dt * ty2* (c1 * (u[k-1, j+2, i+1, 2] * tmp1))- dt * ty1 * c1345 * tmp1- dt * ty1 * dy5;
                    //---------------------------------------------------------------------
                    //   form the third block sub-diagonal
                    //---------------------------------------------------------------------
                    tmp1 = 1.0d / u[k, j+1, i+1, 0];   //tmp1 = 1.0d / u[1,i,j,k+1];
                    tmp2 = tmp1 * tmp1;
                    tmp3 = tmp1 * tmp2;

                    c[j-1, i-1, 0, 0] = -dt * tz1 * dz1;
                    c[j-1, i-1, 1, 0] =   0.0d;
                    c[j-1, i-1, 2, 0] =   0.0d;
                    c[j-1, i-1, 3, 0] = dt * tz2;
                    c[j-1, i-1, 4, 0] =   0.0d;

                    c[j-1, i-1, 0, 1] = dt * tz2* (-(u[k, j+1, i+1, 1]*u[k, j+1, i+1, 3]) * tmp2)- dt * tz1 * (-c34 * tmp2 * u[k, j+1, i+1, 1]);
                    c[j-1, i-1, 1, 1] = dt * tz2 *    (u[k, j+1, i+1, 3] * tmp1)- dt * tz1 * c34 * tmp1- dt * tz1 * dz2;
                    c[j-1, i-1, 2, 1] = 0.0d;
                    c[j-1, i-1, 3, 1] = dt * tz2 *    (u[k, j+1, i+1, 1] * tmp1);
                    c[j-1, i-1, 4, 1] = 0.0d;

                    c[j-1, i-1, 0, 2] = dt * tz2* (-(u[k, j+1, i+1, 2]*u[k, j+1, i+1, 3]) * tmp2)- dt * tz1 * (-c34 * tmp2 * u[k, j+1, i+1, 2]);
                    c[j-1, i-1, 1, 2] = 0.0d;
                    c[j-1, i-1, 2, 2] = dt * tz2 *    (u[k, j+1, i+1, 3] * tmp1)- dt * tz1 * (c34 * tmp1)- dt * tz1 * dz3;
                    c[j-1, i-1, 3, 2] = dt * tz2 *    (u[k, j+1, i+1, 2] * tmp1);
                    c[j-1, i-1, 4, 2] = 0.0d;

                    c[j-1, i-1, 0, 3] = dt * tz2     * (-pow2((u[k, j+1, i+1, 3]*tmp1)) + 0.50d * c2 * ((u[k, j+1, i+1, 1] * u[k, j+1, i+1, 1]+ u[k, j+1, i+1, 2] * u[k, j+1, i+1, 2]+ u[k, j+1, i+1, 3] * u[k, j+1, i+1, 3]) * tmp2))- dt * tz1 * (-r43 * c34 * tmp2 * u[k, j+1, i+1, 3]);
                    c[j-1, i-1, 1, 3] = dt * tz2* (-c2 * (u[k, j+1, i+1, 1] * tmp1));
                    c[j-1, i-1, 2, 3] = dt * tz2* (-c2 * (u[k, j+1, i+1, 2] * tmp1));
                    c[j-1, i-1, 3, 3] = dt * tz2*(2.0d-c2)*(u[k, j+1, i+1, 3] * tmp1)- dt * tz1 * (r43 * c34 * tmp1)- dt * tz1 * dz4;
                    c[j-1, i-1, 4, 3] = dt * tz2 * c2;

                    c[j-1, i-1, 0, 4] = dt * tz2*((c2 * (u[k, j+1, i+1, 1] * u[k, j+1, i+1, 1]+ u[k, j+1, i+1, 2] * u[k, j+1, i+1, 2]+ u[k, j+1, i+1, 3] * u[k, j+1, i+1, 3]) * tmp2- c1 * (u[k, j+1, i+1, 4] * tmp1))* (u[k, j+1, i+1, 3] * tmp1)) - dt * tz1*(-(c34-c1345)*tmp3*(pow2(u[k, j+1, i+1, 1]))-(c34-c1345)*tmp3 * (pow2(u[k, j+1, i+1, 2]))-(r43*c34-c1345)*tmp3*(pow2(u[k, j+1, i+1, 3]))- c1345 * tmp2 * u[k, j+1, i+1, 4]);
                    c[j-1, i-1, 1, 4] = dt * tz2* (-c2*(u[k, j+1, i+1, 1]*u[k, j+1, i+1, 3]) * tmp2)- dt * tz1 * (c34 - c1345) * tmp2 * u[k, j+1, i+1, 1];
                    c[j-1, i-1, 2, 4] = dt * tz2*(-c2* (u[k, j+1, i+1, 2]*u[k, j+1, i+1, 3]) * tmp2)- dt * tz1 * (c34 - c1345) * tmp2 * u[k, j+1, i+1, 2];
                    c[j-1, i-1, 3, 4] = dt * tz2* (c1* (u[k, j+1, i+1, 4] * tmp1) - 0.50d * c2* ((u[k, j+1, i+1, 1]*u[k, j+1, i+1, 1]+ u[k, j+1, i+1, 2]*u[k, j+1, i+1, 2]+ 3.0d*u[k, j+1, i+1, 3]*u[k, j+1, i+1, 3]) * tmp2))- dt * tz1 * (r43*c34 - c1345) * tmp2 * u[k, j+1, i+1, 3];
                    c[j-1, i-1, 4, 4] = dt * tz2* (c1* (u[k, j+1, i+1, 3] * tmp1))- dt * tz1 * c1345 * tmp1- dt * tz1 * dz5;
                }
            }
        }
        // end jacu.f
        // buts.f
        public void buts(int ldmx, int ldmy, int ldmz, int nx, int ny, int nz, int k, double omega, double[, , ,] v, double[, ,] tv, 
                         double[, , ,] d, double[, , ,] xud, double[, , ,] yud, double[, , ,] zud, 
                         int ist, int iend, int jst, int jend, int nx0, int ny0, int ipt, int jpt) {
            //---------------------------------------------------------------------
            //   compute the regular-sparse, block upper triangular solution:
            //                     v <-- [ U-inv ] * v
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //int ldmx, ldmy, ldmz, nx, ny, nz, k; double  omega;
            //double  v[ 5, -1:ldmx+2, -1:ldmy+2, *], tv[5, ldmx, ldmy], 
            //        d[ 5, 5, ldmx, ldmy], udx[ 5, 5, ldmx, ldmy], udy[ 5, 5, ldmx, ldmy], udz[ 5, 5, ldmx, ldmy ];
            //int ist, iend,jst, jend,nx0, ny0,ipt, jpt;
            //---------------------------------------------------------------------
            //  local variables
            //---------------------------------------------------------------------
            int i, j, m, iex;
            double tmp, tmp1;
            double[,] tmat = new double[5, 5];//tmat[5,5] 
            //---------------------------------------------------------------------
            //   receive data from south and east
            //---------------------------------------------------------------------
            iex = 1;
            exchange_1(v, k, iex);
            //Debug 
            //if ((isiz1 + 2) != (ldmx + 2) || (isiz2 + 2) != (ldmy + 2)) {
            //    throw new ArgumentException("Look this code: vetor v");
            //} end debug 

            for(j = jend; j>= jst; j--) { //for(j = jend, jst, -1;
                for(i = iend; i>= ist; i--) { //for(i = iend, ist, -1;
                    for(m = 1; m<= 5; m++) {//tv[ m, i, j ] = 
                        tv[j, i, m] = 
                                         omega * (  zud[j-1, i-1, 0, m-1] * v[k, j+1, i+1, 0]
                                                  + zud[j-1, i-1, 1, m-1] * v[k, j+1, i+1, 1]
                                                  + zud[j-1, i-1, 2, m-1] * v[k, j+1, i+1, 2]
                                                  + zud[j-1, i-1, 3, m-1] * v[k, j+1, i+1, 3]
                                                  + zud[j-1, i-1, 4, m-1] * v[k, j+1, i+1, 4]);
                    }
                }
            }
            for(j = jend; j>=jst; j--) {   //for(j = jend,jst,-1;
                for(i = iend; i>=ist; i--) { //for(i = iend,ist,-1;
                    for(m = 1; m<= 5; m++) {
                        tv[j, i, m] = tv[j, i, m]
                                            + omega * ( yud[j-1, i-1, 0, m-1] * v[k-1, j+2, i+1, 0]
                                                      + xud[j-1, i-1, 0, m-1] * v[k-1, j+1, i+2, 0]
                                                      + yud[j-1, i-1, 1, m-1] * v[k-1, j+2, i+1, 1]
                                                      + xud[j-1, i-1, 1, m-1] * v[k-1, j+1, i+2, 1]
                                                      + yud[j-1, i-1, 2, m-1] * v[k-1, j+2, i+1, 2]
                                                      + xud[j-1, i-1, 2, m-1] * v[k-1, j+1, i+2, 2]
                                                      + yud[j-1, i-1, 3, m-1] * v[k-1, j+2, i+1, 3]
                                                      + xud[j-1, i-1, 3, m-1] * v[k-1, j+1, i+2, 3]
                                                      + yud[j-1, i-1, 4, m-1] * v[k-1, j+2, i+1, 4]
                                                      + xud[j-1, i-1, 4, m-1] * v[k-1, j+1, i+2, 4]);
                    }
                    //---------------------------------------------------------------------
                    //   diagonal block inversion
                    //---------------------------------------------------------------------
                    for(m = 0; m< 5; m++) {
                        tmat[0, m] = d[j-1, i-1, 0, m];
                        tmat[1, m] = d[j-1, i-1, 1, m];
                        tmat[2, m] = d[j-1, i-1, 2, m];
                        tmat[3, m] = d[j-1, i-1, 3, m];
                        tmat[4, m] = d[j-1, i-1, 4, m];
                    }

                    tmp1 = 1.0d /tmat[0, 0];
                    tmp = tmp1 * tmat[0, 1];
                    tmat[1, 1] =  tmat[1, 1] - tmp * tmat[1, 0];
                    tmat[2, 1] =  tmat[2, 1] - tmp * tmat[2, 0];
                    tmat[3, 1] =  tmat[3, 1] - tmp * tmat[3, 0];
                    tmat[4, 1] =  tmat[4, 1] - tmp * tmat[4, 0];
                    tv[j, i, 2] = tv[j, i, 2]- tv[j, i, 1] * tmp;

                    tmp = tmp1 * tmat[0, 2];
                    tmat[1, 2] =  tmat[1, 2] - tmp * tmat[1, 0];
                    tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 0];
                    tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 0];
                    tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 0];
                    tv[j, i, 3] = tv[j, i, 3] - tv[j, i, 1] * tmp;

                    tmp = tmp1 * tmat[0, 3];
                    tmat[1, 3] =  tmat[1, 3] - tmp * tmat[1, 0];
                    tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 0];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 0];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 0];
                    tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 1] * tmp;

                    tmp = tmp1 * tmat[0, 4];
                    tmat[1, 4] =  tmat[1, 4] - tmp * tmat[1, 0];
                    tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 0];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 0];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 0];
                    tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 1] * tmp;

                    tmp1 = 1.0d /tmat[1, 1];
                    tmp = tmp1 * tmat[1, 2];
                    tmat[2, 2] =  tmat[2, 2] - tmp * tmat[2, 1];
                    tmat[3, 2] =  tmat[3, 2] - tmp * tmat[3, 1];
                    tmat[4, 2] =  tmat[4, 2] - tmp * tmat[4, 1];
                    tv[j, i, 3] = tv[j, i, 3] - tv[j, i, 2] * tmp;

                    tmp = tmp1 * tmat[1, 3];
                    tmat[2, 3] =  tmat[2, 3] - tmp * tmat[2, 1];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 1];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 1];
                    tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 2] * tmp;

                    tmp = tmp1 * tmat[1, 4];
                    tmat[2, 4] =  tmat[2, 4] - tmp * tmat[2, 1];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 1];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 1];
                    tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 2] * tmp;

                    tmp1 = 1.0d /tmat[2, 2];
                    tmp = tmp1 * tmat[2, 3];
                    tmat[3, 3] =  tmat[3, 3] - tmp * tmat[3, 2];
                    tmat[4, 3] =  tmat[4, 3] - tmp * tmat[4, 2];
                    tv[j, i, 4] = tv[j, i, 4] - tv[j, i, 3] * tmp;

                    tmp = tmp1 * tmat[2, 4];
                    tmat[3, 4] =  tmat[3, 4] - tmp * tmat[3, 2];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 2];
                    tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 3] * tmp;

                    tmp1 = 1.0d /tmat[3, 3];
                    tmp = tmp1 * tmat[3, 4];
                    tmat[4, 4] =  tmat[4, 4] - tmp * tmat[4, 3];
                    tv[j, i, 5] = tv[j, i, 5] - tv[j, i, 4] * tmp;
                    //---------------------------------------------------------------------
                    //   back substitution
                    //---------------------------------------------------------------------
                    tv[j, i, 5] = tv[j, i, 5]/ tmat[4, 4];
                    tv[j, i, 4] = tv[j, i, 4]- tmat[4, 3] * tv[j, i, 5];
                    tv[j, i, 4] = tv[j, i, 4]/ tmat[3, 3];
                    tv[j, i, 3] = tv[j, i, 3]- tmat[3, 2] * tv[j, i, 4]
                                           - tmat[4, 2] * tv[j, i, 5];
                    tv[j, i, 3] = tv[j, i, 3]/ tmat[2, 2];
                    tv[j, i, 2] = tv[j, i, 2]- tmat[2, 1] * tv[j, i, 3]
                                           - tmat[3, 1] * tv[j, i, 4]
                                           - tmat[4, 1] * tv[j, i, 5];
                    tv[j, i, 2] = tv[j, i, 2]/ tmat[1, 1];
                    tv[j, i, 1] = tv[j, i, 1]- tmat[1, 0]*tv[j, i, 2]
                                           - tmat[2, 0]*tv[j, i, 3]
                                           - tmat[3, 0]*tv[j, i, 4]
                                           - tmat[4, 0]*tv[j, i, 5];
                    tv[j, i, 1] = tv[j, i, 1] /tmat[0, 0];

                    v[k-1, j+1, i+1, 0] = v[k-1, j+1, i+1, 0] - tv[j, i, 1];
                    v[k-1, j+1, i+1, 1] = v[k-1, j+1, i+1, 1] - tv[j, i, 2];
                    v[k-1, j+1, i+1, 2] = v[k-1, j+1, i+1, 2] - tv[j, i, 3];
                    v[k-1, j+1, i+1, 3] = v[k-1, j+1, i+1, 3] - tv[j, i, 4];
                    v[k-1, j+1, i+1, 4] = v[k-1, j+1, i+1, 4] - tv[j, i, 5];
                }
            }
            //---------------------------------------------------------------------
            //   send data to north and west
            //---------------------------------------------------------------------
            iex = 3;
            exchange_1(v, k, iex);
        }
        // end buts.f
        //end ssor.f
        //error.f
        public void error() {
            //---------------------------------------------------------------------
            //   compute the solution error
            //---------------------------------------------------------------------
            int i, j, k, m, iglob, jglob;
            double tmp;
            double[,,,] u000ijk = new double[1, 1, 1, 5]; //u000ijk[5]
            double[] dummy = new double[5]; //dummy[5]
            //int IERROR
            for(m = 0; m< 5; m++) {
                errnm[m] = 0.0d;
                dummy[m] = 0.0d;
            }
            for(k = 2; k<= nz-1; k++) {
                for(j = jst; j<= jend; j++) {
                    jglob = jpt + j;
                    for(i = ist; i<= iend; i++) {
                        iglob = ipt + i;
                        exact(iglob, jglob, k, u000ijk, 0, 0, 0);
                        for(m = 0; m< 5; m++) {
                            tmp = (u000ijk[0, 0, 0, m] - u[k-1, j+1, i+1, m]); //tmp = (u000ijk[m] - u[m,i,j,k]);
                            dummy[m] = dummy[m] + pow2(tmp);
                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //   compute the global sum of individual contributions to dot product.
            //---------------------------------------------------------------------
            //call MPI_ALLREDUCE[ dummy,errnm,5,dp_type,MPI_SUM,MPI_COMM_WORLD,IERROR ];
            worldcomm.Allreduce<double>(dummy, MPI.Operation<double>.Add, ref errnm);

            for(m = 0; m< 5; m++) {
                errnm[m] = Math.Sqrt(errnm[m]/((nx0-2)*(ny0-2)*(nz0-2)));   //sqrt (errnm[m]/((nx0-2)*(ny0-2)*(nz0-2)));
            }
        }
        //exact.f
        public void exact(int i, int j, int k, double[, , ,] u000ijk, int i1, int i2, int i3) {
            //---------------------------------------------------------------------
            //
            //   compute the exact solution at [i,j,k]
            //
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            int m;
            double xi, eta, zeta;

            xi   = ((double)(i-1))/(nx0-1); //( dble ( i - 1 ) ) / ( nx0 - 1 );
            eta  = ((double)(j-1))/(ny0-1); //( dble ( j - 1 ) ) / ( ny0 - 1 );
            zeta = ((double)(k-1))/(nz -1);  //( dble ( k - 1 ) ) / ( nz - 1 );
            for(m = 0; m< 5; m++) {
                u000ijk[i1, i2, i3, m] = ce[0, m] 
                    + ce[1, m]*xi 
                    + ce[2, m]*eta 
                    + ce[3, m]*zeta 
                    + ce[4, m]*xi*xi 
                    + ce[5, m]*eta*eta
                    + ce[6, m]*zeta*zeta 
                    + ce[7, m]*xi*xi*xi 
                    + ce[8, m]*eta*eta*eta 
                    + ce[9, m]*zeta*zeta*zeta 
                    + ce[10, m]*xi*xi*xi*xi
                    + ce[11, m]*eta*eta*eta*eta 
                    + ce[12, m]*zeta*zeta*zeta*zeta;
            }
        }
        //end exact.f
        //end error.f
        //pintgr.f
        public void pintgr() {
            //---------------------------------------------------------------------
            //  local variables
            //---------------------------------------------------------------------
            int i, j, k, ibeg, ifin, ifin1, jbeg, jfin, jfin1, iglob, iglob1, iglob2, jglob, jglob1, jglob2, ind1, ind2;
            double[,] phi1 = new double[isiz3+2, isiz2+2]; //phi1[0:isiz2+1,0:isiz3+1]
            double[,] phi2 = new double[isiz3+2, isiz2+2]; //phi2[0:isiz2+1,0:isiz3+1]
            double frc1, frc2, frc3, dummy;
            //int IERROR
            //---------------------------------------------------------------------
            //   set up the sub-domains for intation in each processor
            //---------------------------------------------------------------------
            ibeg = nx + 1;
            ifin = 0;
            iglob1 = ipt + 1;
            iglob2 = ipt + nx;
            if(iglob1>=ii1 && iglob2<ii2+nx)
                ibeg = 1;
            if(iglob1>ii1-nx && iglob2<=ii2)
                ifin = nx;
            if(ii1>=iglob1 && ii1<=iglob2)
                ibeg = ii1 - ipt;
            if(ii2>=iglob1 && ii2<=iglob2)
                ifin = ii2 - ipt;
            jbeg = ny + 1;
            jfin = 0;
            jglob1 = jpt + 1;
            jglob2 = jpt + ny;
            if(jglob1>=ji1 && jglob2<ji2+ny)
                jbeg = 1;
            if(jglob1>ji1-ny && jglob2<=ji2)
                jfin = ny;
            if(ji1>=jglob1 && ji1<=jglob2)
                jbeg = ji1 - jpt;
            if(ji2>=jglob1 && ji2<=jglob2)
                jfin = ji2 - jpt;
            ifin1 = ifin;
            jfin1 = jfin;
            if(ipt+ifin1 == ii2)
                ifin1 = ifin - 1;
            if(jpt+jfin1 == ji2)
                jfin1 = jfin - 1;
            //---------------------------------------------------------------------
            //   initialize
            //---------------------------------------------------------------------
            for(i = 0; i<=isiz2+1; i++) {
                for(k = 0; k<=isiz3+1; k++) {
                    phi1[k, i] = 0.0;//phi1[i,k] = 0.0; //Obs: change index
                    phi2[k, i] = 0.0;//phi2[i,k] = 0.0; //Obs: change index
                }
            }
            for(j = jbeg; j<=jfin; j++) {
                jglob = jpt + j;
                for(i = ibeg; i<=ifin; i++) {
                    iglob = ipt + i;
                    k = ki1;
                    phi1[j, i] = c2*(u[k-1, j+1, i+1, 4]   //phi1[i,j] = c2*( u[5,i,j,k]
                      -0.50d*(pow2(u[k-1, j+1, i+1, 1])
                            + pow2(u[k-1, j+1, i+1, 2])
                            + pow2(u[k-1, j+1, i+1, 3]))
                                 / u[k-1, j+1, i+1, 0]);
                    k = ki2;
                    phi2[j, i] = c2*(u[k-1, j+1, i+1, 4]
                      -0.50d*(pow2(u[k-1, j+1, i+1, 1])
                            + pow2(u[k-1, j+1, i+1, 2])
                            + pow2(u[k-1, j+1, i+1, 3]))
                                 / u[k-1, j+1, i+1, 0]);
                }
            }
            //---------------------------------------------------------------------
            //  communicate in i and j directions
            //---------------------------------------------------------------------
            exchange_4(phi1, phi2, ibeg, ifin1, jbeg, jfin1); // Obs: phi1 e phi2 will change order index

            frc1 = 0.0d;

            for(j = jbeg; j<=jfin1; j++) {
                for(i = ibeg; i<= ifin1; i++) {
                    frc1 = frc1 + (phi1[j, i]    //frc1=frc1+(phi1[i   , j  ]
                                + phi1[j, i+1]             //+ phi1[i+1 , j  ]  
                                + phi1[j+1, i]             //+ phi1[i   , j+1]  
                                + phi1[j+1, i+1]             //+ phi1[i+1 , j+1]  
                                + phi2[j, i]             //+ phi2[i   , j  ]  
                                + phi2[j, i+1]             //+ phi2[i+1 , j  ]  
                                + phi2[j+1, i]             //+ phi2[i   , j+1]  
                                + phi2[j+1, i+1]);          //+ phi2[i+1 , j+1] );
                }
            }
            //---------------------------------------------------------------------
            //  compute the global sum of individual contributions to frc1
            //---------------------------------------------------------------------
            dummy = frc1;
            //call MPI_ALLREDUCE[ dummy,frc1,1,dp_type,MPI_SUM,MPI_COMM_WORLD,IERROR ];
            frc1 = worldcomm.Allreduce<double>(dummy, MPI.Operation<double>.Add);
            frc1 = dxi * deta * frc1;
            //---------------------------------------------------------------------
            //   initialize
            //---------------------------------------------------------------------
            for(i = 0; i<=isiz2+1; i++) {
                for(k = 0; k<=isiz3+1; k++) {
                    phi1[k, i] = 0.0;  //phi1[i,k] = 0.0;
                    phi2[k, i] = 0.0;  //phi2[i,k] = 0.0;
                }
            }
            jglob = jpt + jbeg;
            ind1 = 0;
            if(jglob==ji1) {
                ind1 = 1;
                for(k = ki1; k<= ki2; k++) {
                    for(i = ibeg; i<= ifin; i++) {
                        iglob = ipt + i;
                        phi1[k, i] = c2 * (u[k-1, jbeg+1, i+1, 4]  //phi1[i,k] = c2 * (u[5,i,jbeg,k]
                        - 0.50d*(pow2(u[k-1, jbeg+1, i+1, 1])
                               + pow2(u[k-1, jbeg+1, i+1, 2])
                               + pow2(u[k-1, jbeg+1, i+1, 3]))
                                    / u[k-1, jbeg+1, i+1, 0]);
                    }
                }
            }
            jglob = jpt + jfin;
            ind2 = 0;
            if(jglob==ji2) {
                ind2 = 1;
                for(k = ki1; k<= ki2; k++) {
                    for(i = ibeg; i<= ifin; i++) {
                        iglob = ipt + i;
                        phi2[k, i] = c2*(u[k-1, jfin+1, i+1, 4]        //phi2[i,k] = c2*( u[5,i,jfin,k] 
                        -0.50d*(pow2(u[k-1, jfin+1, i+1, 1])
                              + pow2(u[k-1, jfin+1, i+1, 2])
                              + pow2(u[k-1, jfin+1, i+1, 3]))
                                   / u[k-1, jfin+1, i+1, 0]);
                    }
                }
            }
            //---------------------------------------------------------------------
            //  communicate in i direction
            //---------------------------------------------------------------------
            if(ind1==1) {
                computeRightSideSouthToNorth(phi1, ibeg, ifin1);
            }
            if(ind2==1) {
                computeRightSideSouthToNorth(phi2, ibeg, ifin1);
            }
            frc2 = 0.0d;
            for(k = ki1; k<= ki2-1; k++) {
                for(i = ibeg; i<= ifin1; i++) {
                    frc2 = frc2 + (phi1[k, i]    //frc2=frc2+(phi1[i   , k  ]    
                                + phi1[k, i+1]             //+ phi1[i+1 , k  ]    
                                + phi1[k+1, i]             //+ phi1[i   , k+1]    
                                + phi1[k+1, i+1]             //+ phi1[i+1 , k+1]    
                                + phi2[k, i]             //+ phi2[i   , k  ]    
                                + phi2[k, i+1]             //+ phi2[i+1 , k  ]    
                                + phi2[k+1, i]             //+ phi2[i   , k+1]    
                                + phi2[k+1, i+1]);          //+ phi2[i+1 , k+1] ); 
                }
            }
            //---------------------------------------------------------------------
            //  compute the global sum of individual contributions to frc2
            //---------------------------------------------------------------------
            dummy = frc2;
            //call MPI_ALLREDUCE[ dummy,frc2,1,dp_type,MPI_SUM,MPI_COMM_WORLD,IERROR ];
            frc2 = worldcomm.Allreduce<double>(dummy, MPI.Operation<double>.Add);
            frc2 = dxi * dzeta * frc2;
            //---------------------------------------------------------------------
            //   initialize
            //---------------------------------------------------------------------
            for(i = 0; i<=isiz2+1; i++) {
                for(k = 0; k<=isiz3+1; k++) {
                    phi1[k, i] = 0.0;//phi1[i,k] = 0.0;
                    phi2[k, i] = 0.0;//phi2[i,k] = 0.0;
                }
            }
            iglob = ipt + ibeg;
            ind1 = 0;
            if(iglob==ii1) {
                ind1 = 1;
                for(k = ki1; k<= ki2; k++) {
                    for(j = jbeg; j<= jfin; j++) {
                        jglob = jpt + j;
                        phi1[k, j] = c2*(u[k-1, j+1, ibeg+1, 4]    //phi1[j,k] = c2*( u[5,ibeg,j,k]
                        -0.50d*(pow2(u[k-1, j+1, ibeg+1, 1])
                              + pow2(u[k-1, j+1, ibeg+1, 2])
                              + pow2(u[k-1, j+1, ibeg+1, 3]))
                                   / u[k-1, j+1, ibeg+1, 0]);
                    }
                }
            }
            iglob = ipt + ifin;
            ind2 = 0;
            if(iglob==ii2) {
                ind2 = 1;
                for(k = ki1; k<= ki2; k++) {
                    for(j = jbeg; j<= jfin; j++) {
                        jglob = jpt + j;
                        phi2[k, j] = c2*(u[k-1, j+1, ifin+1, 4]            //phi2[j,k] = c2*( u[5,ifin,j,k]
                        -0.50d*(pow2(u[k-1, j+1, ifin+1, 1])
                              + pow2(u[k-1, j+1, ifin+1, 2])
                              + pow2(u[k-1, j+1, ifin+1, 3]))
                                   / u[k-1, j+1, ifin+1, 0]);
                    }
                }
            }
            //---------------------------------------------------------------------
            //  communicate in j direction
            //---------------------------------------------------------------------
            if(ind1==1) {
                computeRightSideEastToWest(phi1, jbeg, jfin1);
            }
            if(ind2==1) {
                computeRightSideEastToWest(phi2, jbeg, jfin1);
            }
            frc3 = 0.0d;
            for(k = ki1; k<= ki2-1; k++) {
                for(j = jbeg; j<= jfin1; j++) {
                    frc3 = frc3 + (phi1[k, j]    //frc3 = frc3 + ( phi1[j   , k  ]    
                                + phi1[k, j+1]                  //+ phi1[j+1 , k  ]    
                                + phi1[k+1, j]                  //+ phi1[j   , k+1]    
                                + phi1[k+1, j+1]                  //+ phi1[j+1 , k+1]    
                                + phi2[k, j]                  //+ phi2[j   , k  ]    
                                + phi2[k, j+1]                  //+ phi2[j+1 , k  ]    
                                + phi2[k+1, j]                  //+ phi2[j   , k+1]    
                                + phi2[k+1, j+1]);               //+ phi2[j+1 , k+1] ); 
                }
            }
            //---------------------------------------------------------------------
            //  compute the global sum of individual contributions to frc3
            //---------------------------------------------------------------------
            dummy = frc3;
            //call MPI_ALLREDUCE[ dummy,frc3,1,dp_type,MPI_SUM,MPI_COMM_WORLD,IERROR ];
            frc3 = worldcomm.Allreduce<double>(dummy, MPI.Operation<double>.Add);
            frc3 = deta * dzeta * frc3;
            frc = 0.25d * (frc1 + frc2 + frc3);
        }
        //exchange_4.f
        public void exchange_4(double[,] g, double[,] h, int ibeg, int ifin1, int jbeg, int jfin1) {
            //---------------------------------------------------------------------
            //   compute the right hand side based on exact solution
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //Fortran: double  g[0:isiz2+1,0:isiz3+1],h[0:isiz2+1,0:isiz3+1]
            //C#     : double  g[0:isiz3+1,0:isiz2+1],h[0:isiz3+1,0:isiz2+1]
            //int ibeg, ifin1
            //int jbeg, jfin1
            //---------------------------------------------------------------------
            //  local variables
            //---------------------------------------------------------------------
            int i, j, ny2;
            //MPI.Request[] msgid1 = new MPI.Request[1];  //int msgid1;
            //MPI.Request[] msgid3 = new MPI.Request[1];  //int msgid3;
            //double[] dum = new double[1025];//dum[1024]
            //int STATUS[MPI_STATUS_SIZE]
            //int IERROR
            ny2 = ny + 2;
            //---------------------------------------------------------------------
            //   communicate in the east and west directions
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //   receive from east
            //---------------------------------------------------------------------
            if(jfin1==ny) {
                MPI.Request[] msgid3 = new MPI.Request[1];  //int msgid3;
                double[] dum = new double[2*nx];
                //call MPI_IRECV[ dum,2*nx,dp_type,MPI_ANY_SOURCE,from_e,MPI_COMM_WORLD,msgid3,IERROR ];
                msgid3[0] = worldcomm.ImmediateReceive<double>(east, from_e, dum);
                msgid3[0].Wait(); //call MPI_WAIT[ msgid3, STATUS, IERROR ];
                for(i = 1; i<=nx; i++) {
                    g[ny+1, i] = dum[i-1];     //g[i,ny+1] = dum[i];
                    h[ny+1, i] = dum[i+nx-1];  //h[i,ny+1] = dum[i+nx];
                }
            }
            //---------------------------------------------------------------------
            //   send west
            //---------------------------------------------------------------------
            if(jbeg==1) {
                double[] dum = new double[2*nx];
                for(i = 1; i<=nx; i++) {
                    dum[i-1]    = g[1, i];    //dum[i] = g[i,1];
                    dum[i+nx-1] = h[1, i]; //dum[i+nx] = h[i,1];
                }
                //call MPI_SEND[ dum,2*nx,dp_type,west,from_e,MPI_COMM_WORLD,IERROR ];
                worldcomm.Send<double>(dum, west, from_e);
            }
            //---------------------------------------------------------------------
            //   communicate in the south and north directions
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //   receive from south
            //---------------------------------------------------------------------
            if(ifin1==nx) {
                MPI.Request[] msgid1 = new MPI.Request[1];  //int msgid1;
                double[] dum = new double[2*ny2];
                //call MPI_IRECV[ dum, 2*ny2, dp_type,MPI_ANY_SOURCE,from_s,MPI_COMM_WORLD,msgid1,IERROR ];
                msgid1[0] = worldcomm.ImmediateReceive<double>(south, from_s, dum);
                msgid1[0].Wait(); //call MPI_WAIT[ msgid1, STATUS, IERROR ];

                for(j = 0; j<=ny+1; j++) {
                    g[j, nx+1] = dum[j];      //g[nx+1,j] = dum[j+1];
                    h[j, nx+1] = dum[j+ny2];  //h[nx+1,j] = dum[j+ny2+1];
                }
            }
            //---------------------------------------------------------------------
            //   send north
            //---------------------------------------------------------------------
            if(ibeg==1) {
                double[] dum = new double[2*ny2];
                for(j = 0; j<=ny+1; j++) {
                    dum[j]     = g[j, 1];     //dum[j+1] = g[1,j];
                    dum[j+ny2] = h[j, 1];     //dum[j+ny2+1] = h[1,j];
                }
                //call MPI_SEND[ dum, 2*ny2, dp_type,north,from_s,MPI_COMM_WORLD,IERROR ];
                worldcomm.Send<double>(dum, north, from_s);
            }
        }
        //end exchange_4.f
        // exchange_5.f
        public void computeRightSideSouthToNorth(double[,] g, int ibeg, int ifin1) {
            //---------------------------------------------------------------------
            //   compute the right hand side based on exact solution
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //double  g[0:isiz2+1,0:isiz3+1];
            //int ibeg, ifin1;
            //---------------------------------------------------------------------
            //  local variables
            //---------------------------------------------------------------------
            int k;
            //double[] dum = new double[1025]; //dum[1024]
            //int msgid1
            //int STATUS[MPI_STATUS_SIZE]
            //int IERROR
            //---------------------------------------------------------------------
            //   communicate in the south and north directions
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //   receive from south
            //---------------------------------------------------------------------
            if(ifin1==nx) {
                MPI.Request[] msgid1 = new MPI.Request[1];
                double[] dum = new double[nz];
                //call MPI_IRECV[ dum, nz, dp_type, MPI_ANY_SOURCE, from_s, MPI_COMM_WORLD, msgid1, IERROR ];
                msgid1[0] = worldcomm.ImmediateReceive<double>(south, from_s, dum);
                msgid1[0].Wait();    //call MPI_WAIT[ msgid1, STATUS, IERROR ]
                for(k = 1; k<=nz; k++) {
                    g[k, nx+1] = dum[k-1];  //g[nx+1,k] = dum[k];
                }
            }
            //---------------------------------------------------------------------
            //   send north
            //---------------------------------------------------------------------
            if(ibeg==1) {
                double[] dum = new double[nz];
                for(k = 1; k<=nz; k++) {
                    dum[k-1] = g[k, 1];  //dum[k] = g[1,k];
                }
                //call MPI_SEND[ dum, nz, dp_type, north, from_s, MPI_COMM_WORLD, IERROR ];
                worldcomm.Send<double>(dum, north, from_s);
            }
        }
        //end exchange_5.f
        // exchange_6.f
        public void computeRightSideEastToWest(double[,] g, int jbeg, int jfin1) {
            //---------------------------------------------------------------------
            //   compute the right hand side based on exact solution
            //---------------------------------------------------------------------
            //  input parameters
            //---------------------------------------------------------------------
            //double  g[0:isiz2+1,0:isiz3+1]
            //int jbeg, jfin1
            //---------------------------------------------------------------------
            //  local parameters
            //---------------------------------------------------------------------
            int k;
            //double[] dum = new double[1025]; //dum[1024]
            //int msgid3
            //int STATUS[MPI_STATUS_SIZE]
            //int IERROR
            //---------------------------------------------------------------------
            //   communicate in the east and west directions
            //---------------------------------------------------------------------
            //   receive from east
            //---------------------------------------------------------------------
            if(jfin1==ny) {
                double[] dum = new double[nz];
                MPI.Request[] msgid3 = new MPI.Request[1];
                //call MPI_IRECV[ dum, nz,dp_type,MPI_ANY_SOURCE,from_e,MPI_COMM_WORLD,msgid3,IERROR ];
                msgid3[0] = worldcomm.ImmediateReceive<double>(east, from_e, dum);
                msgid3[0].Wait(); //call MPI_WAIT[ msgid3, STATUS, IERROR ]

                for(k = 1; k<=nz; k++) {
                    g[k, ny+1] = dum[k-1];  //g[ny+1,k] = dum[k];
                }
            }
            //---------------------------------------------------------------------
            //   send west
            //---------------------------------------------------------------------
            if(jbeg==1) {
                double[] dum = new double[nz];
                for(k = 1; k<=nz; k++) {
                    dum[k-1] = g[k, 1];  //dum[k] = g[1,k];
                }
                //call MPI_SEND[ dum, nz, dp_type, west, from_e, MPI_COMM_WORLD, IERROR ];
                worldcomm.Send<double>(dum, west, from_e);
            }
        }
        // end exchange_6.f
        //end pintgr.f
        // verify.f
        public int verify(double[] xcr, double[] xce, double xci) {
            //---------------------------------------------------------------------
            //  verification routine                         
            //---------------------------------------------------------------------
            ////double xcr[5], xce[5], xci
            double[] xcrref = new double[5];//xcrref[5];
            double[] xceref = new double[5];//xceref[5];
            double[] xcrdif = new double[5];//xcrdif[5];
            double[] xcedif = new double[5];//xcedif[5];
            double epsilon, dtref=0.0, xciref, xcidif;
            int m;
            //---------------------------------------------------------------------
            //   tolerance level
            //---------------------------------------------------------------------
            epsilon = 1.0E-08;

            char clss = 'U';
            int verified = 1;

            for(m = 0; m<5; m++) {
                xcrref[m] = 1.0;
                xceref[m] = 1.0;
            }
            xciref = 1.0;

            if((nx0 == 12) && (ny0 == 12) && (nz0 == 12) && (itmax == 50)) {
                clss = 'S';
                dtref = 5.0E-1;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of residual, for the [12X12X12] grid,
                //   after 50 time steps, with  DT = 5.0d-01
                //---------------------------------------------------------------------
                xcrref[0] = 1.6196343210976702E-02;
                xcrref[1] = 2.1976745164821318E-03;
                xcrref[2] = 1.5179927653399185E-03;
                xcrref[3] = 1.5029584435994323E-03;
                xcrref[4] = 3.4264073155896461E-02;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of solution error, for the [12X12X12] grid,
                //   after 50 time steps, with  DT = 5.0d-01
                //---------------------------------------------------------------------
                xceref[0] = 6.4223319957960924E-04;
                xceref[1] = 8.4144342047347926E-05;
                xceref[2] = 5.8588269616485186E-05;
                xceref[3] = 5.8474222595157350E-05;
                xceref[4] = 1.3103347914111294E-03;
                //---------------------------------------------------------------------
                //   Reference value of surface integral, for the [12X12X12] grid,
                //   after 50 time steps, with DT = 5.0d-01
                //---------------------------------------------------------------------
                xciref = 7.8418928865937083E+00;
            }
            else if((nx0 == 33) && (ny0 == 33) && (nz0 == 33) && (itmax == 300)) {
                clss = 'W';   //!SPEC95fp size;
                dtref = 1.5E-3;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of residual, for the [33x33x33] grid,
                //   after 300 time steps, with  DT = 1.5d-3
                //---------------------------------------------------------------------
                xcrref[0] =   0.1236511638192E+02;
                xcrref[1] =   0.1317228477799E+01;
                xcrref[2] =   0.2550120713095E+01;
                xcrref[3] =   0.2326187750252E+01;
                xcrref[4] =   0.2826799444189E+02;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of solution error, for the [33X33X33] grid,
                //---------------------------------------------------------------------
                xceref[0] =   0.4867877144216E+00;
                xceref[1] =   0.5064652880982E-01;
                xceref[2] =   0.9281818101960E-01;
                xceref[3] =   0.8570126542733E-01;
                xceref[4] =   0.1084277417792E+01;
                //---------------------------------------------------------------------
                //   Reference value of surface integral, for the [33X33X33] grid,
                //   after 300 time steps, with  DT = 1.5d-3
                //---------------------------------------------------------------------
                xciref    =   0.1161399311023E+02;
            }
            else if((nx0 == 64) && (ny0 == 64) && (nz0 == 64) && (itmax == 250)) {
                clss = 'A';
                dtref = 2.0E+0;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of residual, for the [64X64X64] grid,
                //   after 250 time steps, with  DT = 2.0d+00
                //---------------------------------------------------------------------
                xcrref[0] = 7.7902107606689367E+02;
                xcrref[1] = 6.3402765259692870E+01;
                xcrref[2] = 1.9499249727292479E+02;
                xcrref[3] = 1.7845301160418537E+02;
                xcrref[4] = 1.8384760349464247E+03;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of solution error, for the [64X64X64] grid,
                //   after 250 time steps, with  DT = 2.0d+00
                //---------------------------------------------------------------------
                xceref[0] = 2.9964085685471943E+01;
                xceref[1] = 2.8194576365003349E+00;
                xceref[2] = 7.3473412698774742E+00;
                xceref[3] = 6.7139225687777051E+00;
                xceref[4] = 7.0715315688392578E+01;
                //---------------------------------------------------------------------
                //   Reference value of surface integral, for the [64X64X64] grid,
                //   after 250 time steps, with DT = 2.0d+00
                //---------------------------------------------------------------------
                xciref = 2.6030925604886277E+01;
            }
            else if((nx0 == 102) && (ny0 == 102) && (nz0 == 102) && (itmax == 250)) {
                clss = 'B';
                dtref = 2.0E+0;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of residual, for the [102X102X102] grid,
                //   after 250 time steps, with  DT = 2.0d+00
                //---------------------------------------------------------------------
                xcrref[0] = 3.5532672969982736E+03;
                xcrref[1] = 2.6214750795310692E+02;
                xcrref[2] = 8.8333721850952190E+02;
                xcrref[3] = 7.7812774739425265E+02;
                xcrref[4] = 7.3087969592545314E+03;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of solution error, for the [102X102X102] 
                //   grid, after 250 time steps, with  DT = 2.0d+00
                //---------------------------------------------------------------------
                xceref[0] = 1.1401176380212709E+02;
                xceref[1] = 8.1098963655421574E+00;
                xceref[2] = 2.8480597317698308E+01;
                xceref[3] = 2.5905394567832939E+01;
                xceref[4] = 2.6054907504857413E+02;
                //---------------------------------------------------------------------
                //   Reference value of surface integral, for the [102X102X102] grid,
                //   after 250 time steps, with DT = 2.0d+00
                //---------------------------------------------------------------------
                xciref = 4.7887162703308227E+01;
            }
            else if((nx0 == 162) && (ny0 == 162) && (nz0 == 162) && (itmax == 250)) {
                clss = 'C';
                dtref = 2.0E+0;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of residual, for the [162X162X162] grid,
                //   after 250 time steps, with  DT = 2.0d+00
                //---------------------------------------------------------------------
                xcrref[0] = 1.03766980323537846E+04;
                xcrref[1] = 8.92212458801008552E+02;
                xcrref[2] = 2.56238814582660871E+03;
                xcrref[3] = 2.19194343857831427E+03;
                xcrref[4] = 1.78078057261061185E+04;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of solution error, for the [162X162X162] 
                //   grid, after 250 time steps, with  DT = 2.0d+00
                //---------------------------------------------------------------------
                xceref[0] = 2.15986399716949279E+02;
                xceref[1] = 1.55789559239863600E+01;
                xceref[2] = 5.41318863077207766E+01;
                xceref[3] = 4.82262643154045421E+01;
                xceref[4] = 4.55902910043250358E+02;
                //---------------------------------------------------------------------
                //   Reference value of surface integral, for the [162X162X162] grid,
                //   after 250 time steps, with DT = 2.0d+00
                //---------------------------------------------------------------------
                xciref = 6.66404553572181300E+01;
            }
            else if((nx0 == 408) && (ny0 == 408) && (nz0 == 408) && (itmax == 300)) {
                clss = 'D';
                dtref = 1.0E+0;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of residual, for the [408X408X408] grid,
                //   after 300 time steps, with  DT = 1.0d+00
                //---------------------------------------------------------------------
                xcrref[0] = 0.4868417937025E+05;
                xcrref[1] = 0.4696371050071E+04;
                xcrref[2] = 0.1218114549776E+05;
                xcrref[3] = 0.1033801493461E+05;
                xcrref[4] = 0.7142398413817E+05;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of solution error, for the [408X408X408] 
                //   grid, after 300 time steps, with  DT = 1.0d+00
                //---------------------------------------------------------------------
                xceref[0] = 0.3752393004482E+03;
                xceref[1] = 0.3084128893659E+02;
                xceref[2] = 0.9434276905469E+02;
                xceref[3] = 0.8230686681928E+02;
                xceref[4] = 0.7002620636210E+03;
                //---------------------------------------------------------------------
                //   Reference value of surface integral, for the [408X408X408] grid,
                //   after 300 time steps, with DT = 1.0d+00
                //---------------------------------------------------------------------
                xciref =    0.8334101392503E+02;
            }
            else if((nx0 == 1020) && (ny0 == 1020) && (nz0 == 1020) && (itmax == 300)) {
                clss = 'E';
                dtref = 0.5E+0;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of residual, for the [1020X1020X1020] grid,
                //   after 300 time steps, with  DT = 0.5d+00
                //---------------------------------------------------------------------
                xcrref[0] = 0.2099641687874E+06;
                xcrref[1] = 0.2130403143165E+05;
                xcrref[2] = 0.5319228789371E+05;
                xcrref[3] = 0.4509761639833E+05;
                xcrref[4] = 0.2932360006590E+06;
                //---------------------------------------------------------------------
                //   Reference values of RMS-norms of solution error, for the [1020X1020X1020] 
                //   grid, after 300 time steps, with  DT = 0.5d+00
                //---------------------------------------------------------------------
                xceref[0] = 0.4800572578333E+03;
                xceref[1] = 0.4221993400184E+02;
                xceref[2] = 0.1210851906824E+03;
                xceref[3] = 0.1047888986770E+03;
                xceref[4] = 0.8363028257389E+03;
                //---------------------------------------------------------------------
                //   Reference value of surface integral, for the [1020X1020X1020] grid,
                //   after 300 time steps, with DT = 0.5d+00
                //---------------------------------------------------------------------
                xciref =    0.9512163272273E+02;
            }
            else {
                verified = 0;
            }
            //---------------------------------------------------------------------
            //    verification test for residuals if gridsize is one of 
            //    the defined grid sizes above [class != 'U']
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //    Compute the difference of solution values and the known reference values.
            //---------------------------------------------------------------------
            for(m = 0; m< 5; m++) {
                xcrdif[m] = Math.Abs((xcr[m]-xcrref[m])/xcrref[m]); //xcrdif[m] = dabs((xcr[m]-xcrref[m])/xcrref[m]);
                xcedif[m] = Math.Abs((xce[m]-xceref[m])/xceref[m]); //xcedif[m] = dabs((xce[m]-xceref[m])/xceref[m]);
            }
            xcidif = Math.Abs((xci - xciref)/xciref); //xcidif = dabs((xci - xciref)/xciref);
            //---------------------------------------------------------------------
            //    Output the comparison of computed results to known cases.
            //---------------------------------------------------------------------
            if(clss != 'U') {
                Console.WriteLine(" Verification being performed for class " + clss);//   write[*, 1990] class //1990 format[/, ' Verification being performed for class ', a]
                Console.WriteLine(" Accuracy setting for epsilon = " + epsilon); //   write [*,2000] epsilon //        2000      format[' Accuracy setting for epsilon = ', E20.13]
                verified = (Math.Abs(dt-dtref) <= epsilon)?1:0; //verified = (dabs(dt-dtref) <= epsilon);
                if(!(verified==1)) {
                    clss = 'U';
                    Console.WriteLine(" DT does not match the reference value of " + dtref); //write [*,1000] dtref //1000         format[' DT does not match the reference value of ', E15.8]
                }
            }
            else {
                Console.WriteLine(" Unknown class"); //   write[*, 1995]  //    1995      format[' Unknown class']
            }
            if(clss != 'U') {
                Console.WriteLine(" Comparison of RMS-norms of residual");//   write [*,2001] // 2001   format[' Comparison of RMS-norms of residual']
            }
            else {
                Console.WriteLine(" RMS-norms of residual");//   write [*, 2005]// 2005   format[' RMS-norms of residual']
            }
            for(m = 0; m< 5; m++) {
                if(clss == 'U') {
                    Console.WriteLine("          " + m + " " + xcr[m]);//write[*, 2015] m, xcr[m];
                }
                else if(xcrdif[m] <= epsilon) {
                    Console.WriteLine("          " + m + " " + xcr[m] + " " + xcrref[m] + " " + xcrdif[m]);
                    //      write [*,2011] m,xcr[m],xcrref[m],xcrdif[m];
                }
                else {
                    Console.WriteLine(" FAILURE: "+m+" "+xcr[m]+" "+xcrref[m]+" "+xcrdif[m]);//write [*,2010] m,xcr[m],xcrref[m],xcrdif[m];
                    verified = 0;
                }
            }
            if(clss != 'U') {
                Console.WriteLine(" Comparison of RMS-norms of solution error");//write [*,2002];// 2002   format[' Comparison of RMS-norms of solution error'];
            }
            else {
                Console.WriteLine(" RMS-norms of solution error");//write [*,2006];// 2006   format[' RMS-norms of solution error'];
            }
            for(m = 0; m< 5; m++) {
                if(clss == 'U') {
                    Console.WriteLine("          " + m + " " + xce[m]);//write[*, 2015] m, xce[m];
                }
                else if(xcedif[m] <= epsilon) {
                    Console.WriteLine("          "+m+" "+xce[m]+" "+xceref[m]+" "+xcedif[m]);//write [*,2011] m,xce[m],xceref[m],xcedif[m];
                }
                else {
                    verified = 0;
                    Console.WriteLine(" FAILURE: "+m+" "+xce[m]+" "+xceref[m]+" "+xcedif[m]);//write [*,2010] m,xce[m],xceref[m],xcedif[m];
                }
            }
            if(clss != 'U') {
                Console.WriteLine(" Comparison of surface integral");//write [*,2025];// 2025   format[' Comparison of surface integral'];
            }
            else {
                Console.WriteLine(" Surface integral");//   write [*,2026];// 2026   format[' Surface integral'];
            }
            if(clss == 'U') {
                Console.WriteLine("          "+xci);//   write[*, 2030] xci;
            }
            else if(xcidif <= epsilon) {
                Console.WriteLine("          " + xci + " " + xciref + " " + xcidif);//   write[*, 2032] xci, xciref, xcidif;
            }
            else {
                verified = 0;
                Console.WriteLine(" FAILURE: " + xci + " " + xciref + " " + xcidif);//write[*, 2031] xci, xciref, xcidif;
            }
            if(clss == 'U') {
                Console.WriteLine("' No reference values provided");//   write[*, 2022]//    2022      format[' No reference values provided']
                Console.WriteLine(" No verification performed");//   write[*, 2023]//    2023      format[' No verification performed']
            }
            else if(verified==1) {
                Console.WriteLine(" Verification Successful");//   write[*, 2020]//2020      format[' Verification Successful']
            }
            else {
                Console.WriteLine(" Verification failed");//   write[*, 2021]//    2021      format[' Verification failed']
            }
            return verified;
        }
        // end verify.f
    }
}