using System;
using System.IO;
using NPB3_0_JAV;
using NPB3_0_JAV.BMInOut;

namespace NPB {
    public class BT:BTBase {

        public BT(char c):base(c){}

        static void Main(String[] argv) {

            BT bt = null;
            BTBase.debug = false;

            try {
                string param = argv[0];
            }
            catch (Exception) {
                argv = new String[1];
                argv[0] = "CLASS=S"; // CLASS DEFAULT, IF USER NOT TYPE CLASS=S IN COMMAND-LINE ARGUMENT
            }
            char paramClass;
            if (!BTBase.debug) {
                BMArgs.ParseCmdLineArgs(argv, BMName);
                paramClass = BMArgs.CLASS;
            }
            else {
                paramClass = 'K';  //DEBUG: CHANGE TO [K=(S and 4 PROCESSORS)] OR [S=(S and 1 PROCESSOR)]
            }                      //DEBUG: OR [T=(A and 4 PROCESSORS)] OR [I=(B and 4 PROCESSORS)]

            try {
                bt = new BT(paramClass);
            }
            catch (OutOfMemoryException e) {
                Console.WriteLine(e.ToString());
                Environment.Exit(0);
            }
            bt.runBenchMark();
        }

        public void runBenchMark() {
            if (!active) goto GoToEnd;
            /*-----------------------------------------------------------------------------------------
                  Root node reads input file [if it exists] } else { takes defaults from parameters
            -----------------------------------------------------------------------------------------*/
            if (node == root) {
                string[] vetTemp = new string[9];
                Console.WriteLine(" NAS Parallel Benchmarks "+npbversion+" -- BT Benchmark ");
                try {
                    Console.Write("Trying Read from input file inputbt.data: ");
                    int[] conf = { 1,1,3,2,2};
                    vetTemp = MPIIO.readFileData("inputbt.data",conf);
                    niter          = int.Parse(vetTemp[0]); 
                    dt             = double.Parse(vetTemp[1]);
                    grid_points[1] = int.Parse(vetTemp[2]); 
                    grid_points[2] = int.Parse(vetTemp[3]);
                    grid_points[3] = int.Parse(vetTemp[4]);
                } catch (System.IO.FileNotFoundException) {
                    Console.WriteLine("inputbt.data not found");
                    fstatus = 1;
                }
                rd_interval = 0;
                if (fstatus == 0) { //This if, should retranslation from Fortran Code. Here, we used just Default values.
                    Console.WriteLine(" Reading from input file inputbt.data");
                    if (iotype != 0) {//read [2,'[A]'] cbuff
                        wr_interval = int.Parse(vetTemp[5]);
                        rd_interval = int.Parse(vetTemp[6]); //read [cbuff,*,iostat=i] wr_interval, rd_interval
                        if (i != 0) rd_interval = 0;
                        if (wr_interval <= 0) wr_interval = wr_default;
                        throw new System.ArgumentException("Check this: comparing with fortran.");
                    }
                    if (iotype == 1) {//read [2,*] collbuf_nodes, collbuf_size //write[*,*] 'collbuf_nodes ', collbuf_nodes //write[*,*] 'collbuf_size  ', collbuf_size
                        throw new System.ArgumentException("Check this: comparing with fortran.");
                    }
                } else { //Default values.
                    Console.WriteLine(" No input file inputbt.data. Using compiled defaults");
                    niter = niter_default;
                    dt    = dt_default;
                    grid_points[1] = problem_size; 
                    grid_points[2] = problem_size; 
                    grid_points[3] = problem_size; 
                    wr_interval = wr_default;      
                    if (iotype == 1) {             //set number of nodes involved in collective buffering to 4,
                        collbuf_nodes = 0;         //unless total number of nodes is smaller than that.
                        collbuf_size = 1000000;    //set buffer size for collective buffering to 1MB per node
                                                   //collbuf_nodes = min[4,no_nodes]
                    }                              //set default to No-File-Hints with a value of 0
                }                                  
                Console.WriteLine("Size: "+grid_points[1]+" x "+grid_points[2]+" x "+grid_points[3]);
                Console.WriteLine("Iterations: " + niter + " dt: " + dt);
                if (no_nodes != total_nodes) Console.WriteLine("Total number of processes: " + total_nodes);
                if (no_nodes != maxcells*maxcells) Console.WriteLine("WARNING: compiled for "+ maxcells*maxcells + " processes ");
                Console.WriteLine(" Number of active processes: " + no_nodes);
                if (iotype == 1) Console.WriteLine("BTIO -- FULL MPI-IO"   + " write interval: "+wr_interval);
                if (iotype == 2) Console.WriteLine("BTIO -- SIMPLE MPI-IO" + " write interval: "+wr_interval);
                if (iotype == 3) Console.WriteLine("BTIO -- EPIO"          + " write interval: "+wr_interval);
                if (iotype == 4) Console.WriteLine("BTIO -- FORTRAN IO"    + " write interval: "+wr_interval);
            }
            worldcomm.Broadcast<int>(ref niter, root);          //call mpi_bcast[niter, 1, MPI_INTEGER,root, comm_setup, error]
            worldcomm.Broadcast<double>(ref dt, root);          //call mpi_bcast[dt, 1, dp_type, root, comm_setup, error]
            worldcomm.Broadcast<int>(ref grid_points[1], root); 
            worldcomm.Broadcast<int>(ref grid_points[2], root); //call mpi_bcast[grid_points[1], 3, MPI_INTEGER, root, comm_setup, error]
            worldcomm.Broadcast<int>(ref grid_points[3], root);
            worldcomm.Broadcast<int>(ref wr_interval, root);    //call mpi_bcast[wr_interval, 1, MPI_INTEGER,root, comm_setup, error]
            worldcomm.Broadcast<int>(ref rd_interval, root);    //call mpi_bcast[rd_interval, 1, MPI_INTEGER,root, comm_setup, error]
            make_set();
            for( c = 1; c<= maxcells; c++){
                if ((cell_size[1,c] > IMAX) || (cell_size[2,c] > JMAX) || (cell_size[3,c] > KMAX) ) {
                    Console.WriteLine(node+" "+c+" ");
                    Console.Write(" "+(cell_size[i,c])+" "); Console.Write(" "+(cell_size[i,c])+" "); Console.Write(" "+(cell_size[i,c])+" ");
                    i=3+1;
                    Console.WriteLine(" Problem size too big for compiled array sizes");
                    goto GoToEnd;
                }
            }
            set_constants();
            initialize();
            idump = 0;
            lhsinit();
            exact_rhs();
            compute_buffer_size(5);
            //---------------------------------------------------------------------
            //      for(one time step to touch all code, and reinitialize
            //---------------------------------------------------------------------
            adi();
            initialize(); 
            timer.resetTimer(2);
            //---------------------------------------------------------------------
            //      Synchronize before placing time stamp
            //---------------------------------------------------------------------
            comm_setup.Barrier();
            timer.resetTimer(1);
            timer.start(1);
            for( step = 1; step<= niter; step++){
                if (node == root) {
                    if (mod(step, 20) == 0 || step == niter || step == 1) {
                        Console.WriteLine(" Time step "+step);
                    }
                }
                adi();
                if (iotype != 0) {
                    timer.start(2);
                    if (mod(step, wr_interval)==0 || step == niter) {
                        if (node == root) {
                            Console.WriteLine("Writing data set, time step "+step);
                        }
                        if (step == niter && rd_interval > 1) {
                            rd_interval = 1;
                        }
                        idump = idump + 1;
                    }
                    timer.stop(2);
                }
            }
            timer.stop(1);
            t = timer.readTimer(1);
            verify(niter, clss, ref verified);
            tmax = comm_setup.Reduce<double>(t, MPI.Operation<double>.Max, root);//mpi_reduce(t,tmax,1,dp_type,MPI_MAX,root,comm_setup,error)
            if (iotype != 0) {
                t = timer.readTimer(2);
                if (t != 0.0) t = 1.0d / t;
                tiominv = comm_setup.Reduce<double>(t, MPI.Operation<double>.Add, root);//mpi_reduce(t,tiominv,1,dp_type,MPI_SUM,root,comm_setup,error)
            }
            if( node == root ) {
                n3 = 1.0d*grid_points[1]*grid_points[2]*grid_points[3];
                navg = (grid_points[1]+grid_points[2]+grid_points[3])/3.0;
                if( tmax != 0.0 ) {
                    mflops = (1.0/1000000)*niter*(3478.8*n3-17655.7*pow2(navg)+28023.7*navg)/tmax;//mflops = 1.0e-6*float[niter]*(3478.8*n3-17655.7*navg**2+28023.7*navg)/tmax;
                } else {
                    mflops = 0.0;
                }
                if (iotype != 0) {
                    mbytes = n3 * 40.0 * idump * (1.0d/1000000);  //mbytes = n3 * 40.0 * idump * (1.0d-6);
                    tiominv = tiominv / no_nodes;
                    t = 0.0;
                    if (tiominv != 0.0) t = 1.0 / tiominv;
                    tpc = 0.0;
                    if (tmax != 0.0) tpc = t * 100.0 / tmax;
                    Console.WriteLine(" BTIO -- statistics:"+" "+"   I/O timing in seconds   : "+t+"   I/O timing percentage   : "+
                        tpc+"   Total data written [MB] : "+mbytes+"   I/O data rate  [MB/sec] : "+mbytes*tiominv);
                }
                MPIIO.print_results("BT", clss, grid_points[1], grid_points[2], grid_points[3], niter, maxcells*maxcells, total_nodes, tmax, mflops, 
                    "          floating point", verified, npbversion);
            }
            GoToEnd: {
                worldcomm.Barrier(); 
                mpi.Dispose();       
            }
        }

        public static double mod(double a, double b) { return (a % b); }

        public double min(int n1, int n2) { return n1<n2?n1:n2; }

        public double max(double n1, double n2) { return n1>n2?n1:n2; }

        public double pow2(double p) { return p * p; }

        // make_set.f
        public void make_set() {
            //---------------------------------------------------------------------
            //     This function allocates space for a set of cells and fills the set     
            //     such that communication between cells on different nodes is only
            //     nearest neighbor                                                   
            //---------------------------------------------------------------------
            int p, i, j, c, dir, size, excess;//, ierr, ierrcode;
            //---------------------------------------------------------------------
            //     compute square root; add small number to allow for roundoff
            //     [note: this is computed in setup_mpi.f also, but prefer to do
            //     it twice because of some include file problems].
            //---------------------------------------------------------------------
            ncells = Convert.ToInt32(Math.Sqrt(no_nodes) + 0.00001d);  // ncells = dint(dsqrt(dble[no_nodes] + 0.00001d0));
            //---------------------------------------------------------------------
            //     this makes coding easier
            //---------------------------------------------------------------------
            p = ncells;
            //---------------------------------------------------------------------
            //     determine the location of the cell at the bottom of the 3D 
            //     array of cells
            //---------------------------------------------------------------------
            cell_coord[1, 1] = (int) mod(node, p);
            cell_coord[2,1] = node/p;
            cell_coord[3,1] = 0;
            //---------------------------------------------------------------------
            //     set the cell_coords for cells in the rest of the z-layers; 
            //     this comes down to a simple linear numbering in the z-direct-
            //     ion, and to the doubly-cyclic numbering in the other dirs     
            //---------------------------------------------------------------------
            for(c=2; c<=p;c++) {
                cell_coord[1,c] = (int) mod(cell_coord[1,c-1]+1,p);
                cell_coord[2,c] = (int) mod(cell_coord[2,c-1]-1+p,p);
                cell_coord[3,c] = c - 1;
            }
            //---------------------------------------------------------------------
            //     offset all the coordinates by 1 to adjust for Fortran arrays
            //---------------------------------------------------------------------
            for(dir = 1; dir<=3; dir++) {
                for(c = 1; c<=p; c++) {
                    cell_coord[dir, c] = cell_coord[dir, c] + 1;
                }
            }
            //---------------------------------------------------------------------
            //     slice[dir,n] contains the sequence number of the cell that is in
            //     coordinate plane n in the dir direction
            //---------------------------------------------------------------------
            for(dir = 1; dir<=3; dir++){
                for(c = 1; c<= p; c++){
                    slice[dir, cell_coord[dir, c]] = c;
                }
            }
            //---------------------------------------------------------------------
            //     fill the predecessor and successor entries, using the indices 
            //     of the bottom cells [they are the same at each level of k 
            //     anyway] acting as if full periodicity pertains; note that p is
            //     added to those arguments to the mod functions that might
            //     otherwise return wrong values when using the modulo function
            //---------------------------------------------------------------------
            i = cell_coord[1,1]-1;
            j = cell_coord[2,1]-1;

            predecessor[1] = ((int)mod(i - 1 + p, p)) + p*j;
            predecessor[2] = i + p*((int)mod(j - 1 + p, p));
            predecessor[3] = ((int)mod(i + 1, p)) + p*((int)mod(j - 1 + p, p));
            successor[1] = ((int)mod(i + 1, p)) + p*j;
            successor[2] = i + p*((int)mod(j + 1, p));
            successor[3] = ((int)mod(i - 1 + p, p)) + p*((int)mod(j + 1, p));

            //---------------------------------------------------------------------
            //     now compute the sizes of the cells                                    
            //---------------------------------------------------------------------
            for(dir= 1; dir<= 3; dir++){
                //---------------------------------------------------------------------
                //     set cell_coord range for each direction                            
                //---------------------------------------------------------------------
                size   = grid_points[dir]/p;
                excess = (int) mod(grid_points[dir],p);
                for(c=1; c<= ncells; c++){
                    if (cell_coord[dir,c] <= excess) {
                        cell_size[dir,c] = size+1;
                        cell_low[dir,c] = (cell_coord[dir,c]-1)*(size+1);
                        cell_high[dir,c] = cell_low[dir,c]+size;
                    } else {
                        cell_size[dir,c] = size;
                        cell_low[dir,c]  = excess*(size+1)+(cell_coord[dir,c]-excess-1)*size;
                        cell_high[dir,c] = cell_low[dir,c]+size-1;
                    }
                    if (cell_size[dir, c] <= 2) {
                        Console.WriteLine("Error: Cell size too small. Min size is 3");
                        worldcomm.Dispose();
                        mpi.Dispose();
                        //throw new System.ArgumentException("Error: Cell size too small. Min size is 3");
                        System.Environment.Exit(0);
                        //call MPI_Abort[mpi_comm_world,ierrcode,ierr];
                        //stop;
                    }
                }
            }
        }
        // End make_set.f

        // set_constants.f
        public void set_constants() {
            ce[1,1]  = 2.0d;
            ce[1,2]  = 0.0d;
            ce[1,3]  = 0.0d;
            ce[1,4]  = 4.0d;
            ce[1,5]  = 5.0d;
            ce[1,6]  = 3.0d;
            ce[1,7]  = 0.5d;
            ce[1,8]  = 0.02d;
            ce[1,9]  = 0.01d;
            ce[1,10] = 0.03d;
            ce[1,11] = 0.5d;
            ce[1,12] = 0.4d;
            ce[1,13] = 0.3d;

            ce[2,1]  = 1.0d;
            ce[2,2]  = 0.0d;
            ce[2,3]  = 0.0d;
            ce[2,4]  = 0.0d;
            ce[2,5]  = 1.0d;
            ce[2,6]  = 2.0d;
            ce[2,7]  = 3.0d;
            ce[2,8]  = 0.01d;
            ce[2,9]  = 0.03d;
            ce[2,10] = 0.02d;
            ce[2,11] = 0.4d;
            ce[2,12] = 0.3d;
            ce[2,13] = 0.5d;

            ce[3,1]  = 2.0d;
            ce[3,2]  = 2.0d;
            ce[3,3]  = 0.0d;
            ce[3,4]  = 0.0d;
            ce[3,5]  = 0.0d;
            ce[3,6]  = 2.0d;
            ce[3,7]  = 3.0d;
            ce[3,8]  = 0.04d;
            ce[3,9]  = 0.03d;
            ce[3,10] = 0.05d;
            ce[3,11] = 0.3d;
            ce[3,12] = 0.5d;
            ce[3,13] = 0.4d;

            ce[4,1]  = 2.0d;
            ce[4,2]  = 2.0d;
            ce[4,3]  = 0.0d;
            ce[4,4]  = 0.0d;
            ce[4,5]  = 0.0d;
            ce[4,6]  = 2.0d;
            ce[4,7]  = 3.0d;
            ce[4,8]  = 0.03d;
            ce[4,9]  = 0.05d;
            ce[4,10] = 0.04d;
            ce[4,11] = 0.2d;
            ce[4,12] = 0.1d;
            ce[4,13] = 0.3d;

            ce[5,1]  = 5.0d;
            ce[5,2]  = 4.0d;
            ce[5,3]  = 3.0d;
            ce[5,4]  = 2.0d;
            ce[5,5]  = 0.1d;
            ce[5,6]  = 0.4d;
            ce[5,7]  = 0.3d;
            ce[5,8]  = 0.05d;
            ce[5,9]  = 0.04d;
            ce[5,10] = 0.03d;
            ce[5,11] = 0.1d;
            ce[5,12] = 0.3d;
            ce[5,13] = 0.2d;

            c1 = 1.4d;
            c2 = 0.4d;
            c3 = 0.1d;
            c4 = 1.0d;
            c5 = 1.4d;

            bt = Math.Sqrt(0.5d);  //bt = dsqrt[0.5d0];

            dnxm1 = 1.0d/(grid_points[1]-1); //dnxm1 = 1.0d/dble[grid_points[1]-1];
            dnym1 = 1.0d/(grid_points[2]-1); //dnym1 = 1.0d/dble[grid_points[2]-1];
            dnzm1 = 1.0d/(grid_points[3]-1); //dnzm1 = 1.0d/dble[grid_points[3]-1];

            c1c2 = c1 * c2;
            c1c5 = c1 * c5;
            c3c4 = c3 * c4;
            c1345 = c1c5 * c3c4;

            conz1 = (1.0d-c1c5);

            tx1 = 1.0d/(dnxm1*dnxm1);
            tx2 = 1.0d/(2.0d*dnxm1);
            tx3 = 1.0d/dnxm1;

            ty1 = 1.0d / (dnym1*dnym1);
            ty2 = 1.0d / (2.0d*dnym1);
            ty3 = 1.0d / dnym1;

            tz1 = 1.0d / (dnzm1 * dnzm1);
            tz2 = 1.0d / (2.0d * dnzm1);
            tz3 = 1.0d / dnzm1;

            dx1 = 0.75d;
            dx2 = 0.75d;
            dx3 = 0.75d;
            dx4 = 0.75d;
            dx5 = 0.75d;

            dy1 = 0.75d;
            dy2 = 0.75d;
            dy3 = 0.75d;
            dy4 = 0.75d;
            dy5 = 0.75d;

            dz1 = 1.0d;
            dz2 = 1.0d;
            dz3 = 1.0d;
            dz4 = 1.0d;
            dz5 = 1.0d;

            dxmax = max(dx3, dx4); //dxmax = dmax1[dx3, dx4];
            dymax = max(dy2, dy4); //dymax = dmax1[dy2, dy4];
            dzmax = max(dz2, dz3); //dzmax = dmax1[dz2, dz3];

            dssp = 0.25d*(max(dx1, max(dy1, dz1)));

            c4dssp = 4.0d * dssp;
            c5dssp = 5.0d * dssp;

            dttx1 = dt*tx1;
            dttx2 = dt*tx2;
            dtty1 = dt*ty1;
            dtty2 = dt*ty2;
            dttz1 = dt*tz1;
            dttz2 = dt*tz2;

            c2dttx1 = 2.0d*dttx1;
            c2dtty1 = 2.0d*dtty1;
            c2dttz1 = 2.0d*dttz1;

            dtdssp = dt*dssp;

            comz1  = dtdssp;
            comz4  = 4.0d*dtdssp;
            comz5  = 5.0d*dtdssp;
            comz6  = 6.0d*dtdssp;

            c3c4tx3 = c3c4*tx3;
            c3c4ty3 = c3c4*ty3;
            c3c4tz3 = c3c4*tz3;

            dx1tx1 = dx1*tx1;
            dx2tx1 = dx2*tx1;
            dx3tx1 = dx3*tx1;
            dx4tx1 = dx4*tx1;
            dx5tx1 = dx5*tx1;

            dy1ty1 = dy1*ty1;
            dy2ty1 = dy2*ty1;
            dy3ty1 = dy3*ty1;
            dy4ty1 = dy4*ty1;
            dy5ty1 = dy5*ty1;

            dz1tz1 = dz1*tz1;
            dz2tz1 = dz2*tz1;
            dz3tz1 = dz3*tz1;
            dz4tz1 = dz4*tz1;
            dz5tz1 = dz5*tz1;

            c2iv  = 2.5d;
            con43 = 4.0d/3.0d;
            con16 = 1.0d/6.0d;

            xxcon1 = c3c4tx3*con43*tx3;
            xxcon2 = c3c4tx3*tx3;
            xxcon3 = c3c4tx3*conz1*tx3;
            xxcon4 = c3c4tx3*con16*tx3;
            xxcon5 = c3c4tx3*c1c5*tx3;

            yycon1 = c3c4ty3*con43*ty3;
            yycon2 = c3c4ty3*ty3;
            yycon3 = c3c4ty3*conz1*ty3;
            yycon4 = c3c4ty3*con16*ty3;
            yycon5 = c3c4ty3*c1c5*ty3;

            zzcon1 = c3c4tz3*con43*tz3;
            zzcon2 = c3c4tz3*tz3;
            zzcon3 = c3c4tz3*conz1*tz3;
            zzcon4 = c3c4tz3*con16*tz3;
            zzcon5 = c3c4tz3*c1c5*tz3;
        }
        // End set_constants.f

        // initialize.f
        public void initialize() {
            //---------------------------------------------------------------------
            //     This public void initializes the field variable u using 
            //     tri-linear transfinite interpolation of the boundary values     
            //---------------------------------------------------------------------
            //include 'header.h'
            int c, i, j, k, m, ii, jj, kk, ix, iy, iz;
            double  xi, eta, zeta, Pxi, Peta, Pzeta;
            double[] temp = new double[5+le];        //temp[5]
            double[,,] Pface = new double[2+le,3+le,5+le]; //Pface[5,3,2]
            //---------------------------------------------------------------------
            //  Later [in compute_rhs] we compute 1/u for every element. A few of 
            //  the corner elements are not used, but it convenient [and faster] 
            //  to compute the whole thing with a simple loop. Make sure those 
            //  values are nonzero by initializing the whole thing here. 
            //---------------------------------------------------------------------
            for(c = 1; c<= ncells; c++){ //for(c = 1; c<= ncells; c++){
                for(kk = -1; kk<= KMAX; kk++){ //for(kk = -1; kk<= KMAX; kk++)
                    for(jj = -1; jj<= JMAX; jj++){ //for(jj = -1; jj<= JMAX; jj++)
                        for(ii = -1; ii<= IMAX; ii++){ //for(ii = -1; ii<= IMAX; ii++)
                            for(m = 1; m<= 5; m++){  //for(m = 1; m<= 5; m++){
                                u[c, kk+2, jj+2, ii+2, m] = 1.0; //u[m, ii, jj, kk, c] = 1.0;
                            }
                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //     first store the "interpolated" values everywhere on the grid    
            //---------------------------------------------------------------------
            for(c=1; c<= ncells; c++){
                kk = 0;
                for(k = cell_low[3,c]; k<= cell_high[3,c]; k++){
                    zeta = k*dnzm1;
                    jj = 0;
                    for(j = cell_low[2,c]; j<= cell_high[2,c]; j++){
                        eta = j*dnym1;
                        ii = 0;
                        for(i = cell_low[1,c]; i<= cell_high[1,c]; i++){
                            xi = i*dnxm1;
                            for(ix = 1; ix<= 2; ix++){
                                exact_solution3((ix - 1), eta, zeta, ref Pface, ix, 1); //call exact_solution[dble[ix-1], eta, zeta, Pface[ix,1,1]];
                            }
                            for(iy = 1; iy<= 2; iy++){
                                exact_solution3(xi, (iy-1) , zeta, ref Pface, iy, 2); //call exact_solution[xi, dble[iy-1] , zeta, Pface[iy,2,1]];
                            }
                            for(iz = 1; iz<= 2; iz++){
                                exact_solution3(xi, eta, (iz-1), ref Pface, iz,3); //call exact_solution[xi, eta, dble[iz-1], Pface[iz,3,1]];
                            }
                            for(m = 1; m<= 5; m++){
                                Pxi   = xi   * Pface[2,1,m] + (1.0d-xi)   * Pface[1,1,m];
                                Peta  = eta  * Pface[2,2,m] + (1.0d-eta)  * Pface[1,2,m];
                                Pzeta = zeta * Pface[2,3,m] + (1.0d-zeta) * Pface[1,3,m];
                                u[c,kk+2,jj+2,ii+2,m] = Pxi + Peta + Pzeta - Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + Pxi*Peta*Pzeta;
                            }
                            ii = ii + 1;
                        }
                        jj = jj + 1;
                    }
                    kk = kk+1;
                }
            }
            //---------------------------------------------------------------------
            //     now store the exact values on the boundaries        
            //---------------------------------------------------------------------
            //---------------------------------------------------------------------
            //     west face                                                  
            //---------------------------------------------------------------------
            c = slice[1,1];
            ii = 0;
            xi = 0.0d;
            kk = 0;
            for(k = cell_low[3,c]; k<= cell_high[3,c]; k++){
                zeta = k*dnzm1;
                jj = 0;
                for(j = cell_low[2,c]; j<= cell_high[2,c]; j++){
                    eta = j*dnym1;
                    exact_solution1(xi, eta, zeta, ref temp);
                    for(m = 1; m<= 5; m++){
                        u[c,kk+2,jj+2,ii+2,m] = temp[m]; //u[m,ii,jj,kk,c] = temp[m];
                    }
                    jj = jj + 1;
                }
                kk = kk + 1;
            }
            //---------------------------------------------------------------------
            //     east face                                                      
            //---------------------------------------------------------------------
            c  = slice[1,ncells];
            ii = cell_size[1,c]-1;
            xi = 1.0d;
            kk = 0;
            for(k = cell_low[3,c]; k<= cell_high[3,c]; k++){
                zeta = k*dnzm1;
                jj = 0;
                for(j = cell_low[2,c]; j<= cell_high[2,c]; j++){
                    eta = j*dnym1;
                    exact_solution1(xi, eta, zeta, ref temp);
                    for(m = 1; m<= 5; m++){
                        u[c,kk+2,jj+2,ii+2,m] = temp[m]; //u[m,ii,jj,kk,c] = temp[m];
                    }
                    jj = jj + 1;
                }
                kk = kk + 1;
            }
            //---------------------------------------------------------------------
            //     south face                                                 
            //---------------------------------------------------------------------
            c = slice[2,1];
            jj = 0;
            eta = 0.0d;
            kk = 0;
            for(k = cell_low[3,c]; k<= cell_high[3,c]; k++){
                zeta = k*dnzm1;
                ii = 0;
                for(i = cell_low[1,c]; i<= cell_high[1,c]; i++){
                    xi = i*dnxm1;
                    exact_solution1(xi, eta, zeta, ref temp);
                    for(m = 1; m<= 5; m++){
                        u[c,kk+2,jj+2,ii+2,m] = temp[m];  //u[m,ii,jj,kk,c] = temp[m];
                    }
                    ii = ii + 1;
                }
                kk = kk + 1;
            }
            //---------------------------------------------------------------------
            //     north face                                    
            //---------------------------------------------------------------------
            c = slice[2,ncells];
            jj = cell_size[2,c]-1;
            eta = 1.0d;
            kk = 0;
            for(k = cell_low[3,c]; k<= cell_high[3,c]; k++){
                zeta = k*dnzm1;
                ii = 0;
                for(i = cell_low[1,c]; i<= cell_high[1,c]; i++){
                    xi = i*dnxm1;
                    exact_solution1(xi, eta, zeta, ref temp);
                    for(m = 1; m<= 5; m++){
                        u[c,kk+2,jj+2,ii+2,m] = temp[m];  // u[m,ii,jj,kk,c] = temp[m];
                    }
                    ii = ii + 1;
                }
                kk = kk + 1;
            }
            //---------------------------------------------------------------------
            //     bottom face                                       
            //---------------------------------------------------------------------
            c = slice[3,1];
            kk = 0;
            zeta = 0.0d;
            jj = 0;
            for(j = cell_low[2,c]; j<= cell_high[2,c]; j++){
                eta = j*dnym1;
                ii = 0;
                for(i =cell_low[1,c]; i<= cell_high[1,c]; i++){
                    xi = i*dnxm1;
                    exact_solution1(xi, eta, zeta, ref temp);
                    for(m = 1; m<= 5; m++){
                        u[c,kk+2,jj+2,ii+2,m] = temp[m]; //u[m,ii,jj,kk,c] = temp[m];
                    }
                    ii = ii + 1;
                }
                jj = jj + 1;
            }
            //---------------------------------------------------------------------
            //     top face     
            //---------------------------------------------------------------------
            c = slice[3,ncells];
            kk = cell_size[3,c]-1;
            zeta = 1.0d;
            jj = 0;
            for(j = cell_low[2,c]; j<= cell_high[2,c]; j++){
                eta = j*dnym1;
                ii = 0;
                for(i =cell_low[1,c]; i<= cell_high[1,c]; i++){
                    xi = i*dnxm1;
                    exact_solution1(xi, eta, zeta, ref temp);
                    for(m = 1; m<= 5; m++){
                        u[c,kk+2,jj+2,ii+2,m] = temp[m];  //u[m,ii,jj,kk,c] = temp[m];
                    }
                    ii = ii + 1;
                }
                jj = jj + 1;
            }
        }

        public void lhsinit() {
            int i, j, k, d, c, m, n;
            //---------------------------------------------------------------------
            //     loop over all cells                                       
            //---------------------------------------------------------------------
            for(c = 1; c<= ncells; c++){
                //---------------------------------------------------------------------
                //     first, initialize the start and end arrays
                //---------------------------------------------------------------------
                for(d = 1; d<= 3; d++){
                    if (cell_coord[d,c] == 1) {
                        start[d,c] = 1;
                    } else { 
                        start[d,c] = 0;
                    }
                    if (cell_coord[d,c] == ncells) {
                        end[d,c] = 1;
                    } else {
                        end[d,c] = 0;
                    }
                }
                //---------------------------------------------------------------------
                //     zero the whole left hand side for starters
                //---------------------------------------------------------------------
                for(k = 0; k<= cell_size[3,c]-1; k++){
                    for(j = 0; j<= cell_size[2,c]-1; j++){
                        for(i = 0; i<= cell_size[1,c]-1; i++){
                            for(m = 1; m<=5; m++){
                                for(n = 1; n<= 5; n++){
                                    lhsc[c,k+1,j+1,i+1,n,m] = 0.0d; //lhsc[m,n,i,j,k,c] = 0.0d;
                                }
                            }
                        }
                    }
                }
            }
        }

        public void lhsabinit(ref double[,,] lhsa, ref double[,,] lhsb, int size) {
            //lhsa[5, 5, -1:size], 
            //lhsb[5, 5, -1:size];
            int i, m, n;
            if (MAX_CELL_DIM < size) { throw new ArgumentException("Instruction Debug: Verificar size vetor lhsb"); }
            //---------------------------------------------------------------------
            //     next, set all diagonal values to 1. This is overkill, but convenient
            //---------------------------------------------------------------------
            for(i = 0; i<= size; i++){
                for(m = 1; m<= 5; m++){
                    for(n = 1; n<= 5; n++){
                        lhsa[i+1,n,m] = 0.0d; //lhsa[m,n,i] = 0.0d;
                        lhsb[i+1,n,m] = 0.0d; //lhsb[m,n,i] = 0.0d;
                    }
                    lhsb[i+1,m,m] = 1.0d; //lhsb[m,m,i] = 1.0d;
                }
            }
        }

        public void exact_solution3(double xi, double eta, double zeta, ref double[,,] dtemp, int d1, int d2) {
            //---------------------------------------------------------------------
            //     this function returns the exact solution at point xi, eta, zeta  
            //---------------------------------------------------------------------
            //dtemp(5)
            int m;
            for(m = 1; m<= 5;m++){
                dtemp[d1,d2,m] =  ce[m,1] + xi*(ce[m,2] + xi*(ce[m,5] + xi*(ce[m,8] + xi*ce[m,11]))) + 
                    eta*(ce[m,3] + eta*(ce[m,6] + eta*(ce[m,9] + eta*ce[m,12])))+ 
                    zeta*(ce[m,4] + zeta*(ce[m,7] + zeta*(ce[m,10] + zeta*ce[m,13])));
            }
        }

        public void exact_solution1(double xi, double eta, double zeta, ref double[] dtemp) {
            //---------------------------------------------------------------------
            //     this function returns the exact solution at point xi, eta, zeta  
            //---------------------------------------------------------------------
            //dtemp[5]
            int m;
            for (m = 1; m <= 5; m++) {
                dtemp[m] = ce[m, 1] + xi * (ce[m, 2] + xi * (ce[m, 5] + xi * (ce[m, 8] + xi * ce[m, 11]))) + 
                    eta * (ce[m, 3] + eta * (ce[m, 6] + eta * (ce[m, 9] + eta * ce[m, 12]))) + 
                    zeta * (ce[m, 4] + zeta * (ce[m, 7] + zeta * (ce[m, 10] + zeta * ce[m, 13])));
            }
        }
        // End initialize.f

        // exact_rhs.f
        public void exact_rhs() {
            //---------------------------------------------------------------------
            //     compute the right hand side based on exact solution
            //---------------------------------------------------------------------
            double[] dtemp = new double[5+le];
            double xi, eta, zeta, dtpp;
            int c, m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;
            //---------------------------------------------------------------------
            //     loop over all cells owned by this node                   
            //---------------------------------------------------------------------
            for(c = 1; c<= ncells; c++){
                //---------------------------------------------------------------------
                //     initialize                                  
                //---------------------------------------------------------------------
                for(k= 0; k<= cell_size[3,c]-1; k++){
                    for(j = 0; j<= cell_size[2,c]-1; j++){
                        for(i = 0; i<= cell_size[1,c]-1; i++){ //forcing[5,   0:IMAX-1, 0:JMAX-1, 0:KMAX-1, maxcells]
                            for(m = 1; m<= 5; m++){             //forcing new double[maxcells+1,KMAX,JMAX,IMAX,6];
                                forcing[c,k,j,i,m] = 0.0d; //forcing[m,i,j,k,c] = 0.0d;
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     xi-direction flux differences                      
                //---------------------------------------------------------------------
                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    zeta = (k+cell_low[3,c])*dnzm1;  //zeta = dble[k+cell_low[3,c]] * dnzm1;
                    for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                        eta = (j+cell_low[2,c])*dnym1;
                        for(i=-2*(1-start[1,c]); i<= cell_size[1,c]+1-2*end[1,c]; i++){
                            xi = (i+cell_low[1,c])*dnxm1;
                            exact_solution1(xi, eta, zeta, ref dtemp);
                            for(m = 1; m<= 5; m++){
                                ue[2+i,m] = dtemp[m];
                            }
                            dtpp = 1.0d / dtemp[1];
                            for(m = 2; m<= 5; m++){
                                buf[2+i,m] = dtpp * dtemp[m];
                            }
                            cuf[2+i]   = buf[2+i,2] * buf[2+i,2];
                            buf[2+i,1] = cuf[2+i] + buf[2+i,3] * buf[2+i,3] + buf[2+i,4] * buf[2+i,4];
                            q[2+i] = 0.5d*(buf[2+i,2]*ue[2+i,2] + buf[2+i,3]*ue[2+i,3] + buf[2+i,4]*ue[2+i,4]);
                        }
                        for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                            im1 = i-1;
                            ip1 = i+1;
                            forcing[c,k,j,i,1] = forcing[c,k,j,i,1] - tx2*(ue[2+ip1,2]-ue[2+im1,2])+dx1tx1*(ue[2+ip1,1]-2.0d*ue[2+i,1]+ue[2+im1,1]);
                            forcing[c,k,j,i,2] = forcing[c,k,j,i,2] - tx2*((ue[2+ip1,2]*buf[2+ip1,2]+c2*(ue[2+ip1,5]-q[2+ip1]))-(ue[2+im1,2]*buf[2+im1,2]+c2*(ue[2+im1,5]-q[2+im1])))+xxcon1*(buf[2+ip1,2]-2.0d*buf[2+i,2]+buf[2+im1,2])+dx2tx1*(ue[2+ip1,2]-2.0d* ue[2+i,2]+ue[2+im1,2]);
                            forcing[c,k,j,i,3] = forcing[c,k,j,i,3] - tx2*(ue[2+ip1,3]*buf[2+ip1,2]-ue[2+im1,3]*buf[2+im1,2])+xxcon2*(buf[2+ip1,3]-2.0d*buf[2+i,3]+buf[2+im1,3])+dx3tx1*( ue[2+ip1,3]-2.0d*ue[2+i,3] +ue[2+im1,3]);
                            forcing[c,k,j,i,4] = forcing[c,k,j,i,4] - tx2*(ue[2+ip1,4]*buf[2+ip1,2]-ue[2+im1,4]*buf[2+im1,2])+xxcon2*(buf[2+ip1,4]-2.0d*buf[2+i,4]+buf[2+im1,4])+dx4tx1*( ue[2+ip1,4]-2.0d* ue[2+i,4]+ ue[2+im1,4]);
                            forcing[c,k,j,i,5] = forcing[c,k,j,i,5] - tx2*(buf[2+ip1,2]*(c1*ue[2+ip1,5]-c2*q[2+ip1])-buf[2+im1,2]*(c1*ue[2+im1,5]-c2*q[2+im1]))+0.5d*xxcon3*(buf[2+ip1,1]-2.0d*buf[2+i,1]+buf[2+im1,1])+xxcon4*(cuf[2+ip1]-2.0d*cuf[2+i]+cuf[2+im1])+xxcon5*(buf[2+ip1,5]-2.0d*buf[2+i,5]+buf[2+im1,5])+dx5tx1*( ue[2+ip1,5]-2.0d* ue[2+i,5]+ ue[2+im1,5]);
                        }
                        //---------------------------------------------------------------------
                        //      Fourth-order dissipation                         
                        //---------------------------------------------------------------------
                        if (start[1,c] > 0) {
                            for(m = 1; m<= 5; m++){
                                i = 1;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp*(5.0d*ue[2+i,m]-4.0d*ue[2+i+1,m]+ue[2+i+2,m]);
                                i = 2;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp*(-4.0d*ue[2+i-1,m] + 6.0d*ue[2+i,m] - 4.0d*ue[2+i+1,m] + ue[2+i+2,m]);
                            }
                        }
                        for(i = start[1,c]*3; i<= cell_size[1,c]-3*end[1,c]-1; i++){
                            for(m = 1; m<= 5; m++){
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp*(ue[2+i-2,m] - 4.0d*ue[2+i-1,m] + 6.0d*ue[2+i,m] - 4.0d*ue[2+i+1,m] + ue[2+i+2,m]);
                            }
                        }
                        if (end[1,c] > 0) {
                            for(m = 1; m<= 5; m++){
                                i = cell_size[1,c]-3;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (ue[2+i-2,m] - 4.0d*ue[2+i-1,m] + 6.0d*ue[2+i,m] - 4.0d*ue[2+i+1,m]);
                                i = cell_size[1,c]-2;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (ue[2+i-2,m] - 4.0d*ue[2+i-1,m] + 5.0d*ue[2+i,m]);
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     eta-direction flux differences             
                //---------------------------------------------------------------------
                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    zeta = (k+cell_low[3,c]) * dnzm1;
                    for(i=start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        xi = (i+cell_low[1,c]) * dnxm1;
                        for(j=-2*(1-start[2,c]); j<= cell_size[2,c]+1-2*end[2,c]; j++){
                            eta = (j+cell_low[2,c]) * dnym1;
                            exact_solution1(xi, eta, zeta, ref dtemp);
                            for(m = 1; m<= 5; m++){
                                ue[2+j,m] = dtemp[m];
                            }
                            dtpp = 1.0d/dtemp[1];
                            for(m = 2; m<= 5; m++){
                                buf[2+j,m] = dtpp * dtemp[m];
                            }
                            cuf[2+j]   = buf[2+j,3] * buf[2+j,3];
                            buf[2+j,1] = cuf[2+j] + buf[2+j,2] * buf[2+j,2] + buf[2+j,4] * buf[2+j,4];
                            q[2+j] = 0.5d*(buf[2+j,2]*ue[2+j,2] + buf[2+j,3]*ue[2+j,3] + buf[2+j,4]*ue[2+j,4]);
                        }
                        for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                            jm1 = j-1;
                            jp1 = j+1;
                            forcing[c,k,j,i,1] = forcing[c,k,j,i,1] - ty2*( ue[2+jp1,3]-ue[2+jm1,3] )+ dy1ty1*(ue[2+jp1,1]-2.0d*ue[2+j,1]+ue[2+jm1,1]);
                            forcing[c,k,j,i,2] = forcing[c,k,j,i,2] - ty2*(ue[2+jp1,2]*buf[2+jp1,3]-ue[2+jm1,2]*buf[2+jm1,3])+ yycon2*(buf[2+jp1,2]-2.0d*buf[2+j,2]+buf[2+jm1,2])+ dy2ty1*( ue[2+jp1,2]-2.0* ue[2+j,2]+ ue[2+jm1,2]);
                            forcing[c,k,j,i,3] = forcing[c,k,j,i,3] - ty2*((ue[2+jp1,3]*buf[2+jp1,3]+c2*(ue[2+jp1,5]-q[2+jp1]))- (ue[2+jm1,3]*buf[2+jm1,3]+c2*(ue[2+jm1,5]-q[2+jm1])))+ yycon1*(buf[2+jp1,3]-2.0d*buf[2+j,3]+buf[2+jm1,3])+ dy3ty1*( ue[2+jp1,3]-2.0d*ue[2+j,3] +ue[2+jm1,3]);
                            forcing[c,k,j,i,4] = forcing[c,k,j,i,4] - ty2*(ue[2+jp1,4]*buf[2+jp1,3]-ue[2+jm1,4]*buf[2+jm1,3])+yycon2*(buf[2+jp1,4]-2.0d*buf[2+j,4]+buf[2+jm1,4])+dy4ty1*(ue[2+jp1,4]-2.0d*ue[2+j,4]+ ue[2+jm1,4]);
                            forcing[c,k,j,i,5] = forcing[c,k,j,i,5] - ty2*(buf[2+jp1,3]*(c1*ue[2+jp1,5]-c2*q[2+jp1])-buf[2+jm1,3]*(c1*ue[2+jm1,5]-c2*q[2+jm1]))+0.5d*yycon3*(buf[2+jp1,1]-2.0d*buf[2+j,1]+buf[2+jm1,1])+yycon4*(cuf[2+jp1]-2.0d*cuf[2+j]+cuf[2+jm1])+yycon5*(buf[2+jp1,5]-2.0d*buf[2+j,5]+buf[2+jm1,5])+dy5ty1*(ue[2+jp1,5]-2.0d*ue[2+j,5]+ue[2+jm1,5]);
                        }
                        //---------------------------------------------------------------------
                        //     Fourth-order dissipation                      
                        //---------------------------------------------------------------------
                        if (start[2,c] > 0) {
                            for(m = 1; m<= 5; m++){
                                j = 1;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp *(5.0d*ue[2+j,m] - 4.0d*ue[2+j+1,m] +ue[2+j+2,m]);
                                j = 2;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp*(-4.0d*ue[2+j-1,m] + 6.0d*ue[2+j,m] - 4.0d*ue[2+j+1,m] + ue[2+j+2,m]);
                            }
                        }
                        for(j = start[2,c]*3; j<= cell_size[2,c]-3*end[2,c]-1; j++){
                            for(m = 1; m<= 5; m++){
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp*(ue[2+j-2,m] - 4.0d*ue[2+j-1,m] + 6.0d*ue[2+j,m] - 4.0d*ue[2+j+1,m] + ue[2+j+2,m]);
                            }
                        }
                        if (end[2,c] > 0) {
                            for(m = 1; m<= 5; m++){
                                j = cell_size[2,c]-3;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (ue[2+j-2,m] - 4.0d*ue[2+j-1,m] + 6.0d*ue[2+j,m] - 4.0d*ue[2+j+1,m]);
                                j = cell_size[2,c]-2;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (ue[2+j-2,m] - 4.0d*ue[2+j-1,m] + 5.0d*ue[2+j,m]);
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     zeta-direction flux differences                      
                //---------------------------------------------------------------------
                for(j=start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    eta = (j+cell_low[2,c]) * dnym1;
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        xi = (i+cell_low[1,c]) * dnxm1;
                        for(k=-2*(1-start[3,c]); k<= cell_size[3,c]+1-2*end[3,c]; k++){
                            zeta = (k+cell_low[3,c]) * dnzm1;
                            exact_solution1(xi, eta, zeta, ref dtemp);
                            for(m = 1; m<= 5; m++){
                                ue[2+k,m] = dtemp[m];
                            }
                            dtpp = 1.0d/dtemp[1];
                            for(m = 2; m<= 5; m++){
                                buf[2+k,m] = dtpp * dtemp[m];
                            }
                            cuf[2+k]   = buf[2+k,4] * buf[2+k,4];
                            buf[2+k,1] = cuf[2+k] + buf[2+k,2] * buf[2+k,2] + buf[2+k,3] * buf[2+k,3];
                            q[2+k] = 0.5d*(buf[2+k,2]*ue[2+k,2] + buf[2+k,3]*ue[2+k,3] + buf[2+k,4]*ue[2+k,4]);
                        }
                        for(k=start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                            km1 = k-1;
                            kp1 = k+1;
                            forcing[c,k,j,i,1] = forcing[c,k,j,i,1] - tz2*( ue[2+kp1,4]-ue[2+km1,4] )+dz1tz1*(ue[2+kp1,1]-2.0d*ue[2+k,1]+ue[2+km1,1]);
                            forcing[c,k,j,i,2] = forcing[c,k,j,i,2] - tz2 * (ue[2+kp1,2]*buf[2+kp1,4]-ue[2+km1,2]*buf[2+km1,4])+ zzcon2*(buf[2+kp1,2]-2.0d*buf[2+k,2]+buf[2+km1,2])+ dz2tz1*( ue[2+kp1,2]-2.0d* ue[2+k,2]+ ue[2+km1,2]);
                            forcing[c,k,j,i,3] = forcing[c,k,j,i,3] - tz2 * (ue[2+kp1,3]*buf[2+kp1,4]-ue[2+km1,3]*buf[2+km1,4])+ zzcon2*(buf[2+kp1,3]-2.0d*buf[2+k,3]+buf[2+km1,3])+ dz3tz1*(ue[2+kp1,3]-2.0d*ue[2+k,3]+ue[2+km1,3]);
                            forcing[c,k,j,i,4] = forcing[c,k,j,i,4] - tz2 * ((ue[2+kp1,4]*buf[2+kp1,4]+c2*(ue[2+kp1,5]-q[2+kp1]))-(ue[2+km1,4]*buf[2+km1,4]+c2*(ue[2+km1,5]-q[2+km1])))+zzcon1*(buf[2+kp1,4]-2.0d*buf[2+k,4]+buf[2+km1,4])+dz4tz1*( ue[2+kp1,4]-2.0d*ue[2+k,4] +ue[2+km1,4]);
                            forcing[c,k,j,i,5] = forcing[c,k,j,i,5] - tz2*(buf[2+kp1,4]*(c1*ue[2+kp1,5]-c2*q[2+kp1])-buf[2+km1,4]*(c1*ue[2+km1,5]-c2*q[2+km1]))+0.5d*zzcon3*(buf[2+kp1,1]-2.0d*buf[2+k,1]+buf[2+km1,1])+zzcon4*(cuf[2+kp1]-2.0d*cuf[2+k]+cuf[2+km1])+zzcon5*(buf[2+kp1,5]-2.0d*buf[2+k,5]+buf[2+km1,5])+dz5tz1*( ue[2+kp1,5]-2.0d*ue[2+k,5]+ ue[2+km1,5]);
                        }
                        //---------------------------------------------------------------------
                        //     Fourth-order dissipation                        
                        //---------------------------------------------------------------------
                        if (start[3,c] > 0) {
                            for(m = 1; m<= 5; m++){
                                k = 1;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (5.0d*ue[2+k,m] - 4.0d*ue[2+k+1,m] +ue[2+k+2,m]);
                                k = 2;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (-4.0d*ue[2+k-1,m] + 6.0d*ue[2+k,m] - 4.0d*ue[2+k+1,m] + ue[2+k+2,m]);
                            }
                        }
                        for(k = start[3,c]*3; k<= cell_size[3,c]-3*end[3,c]-1; k++){
                            for(m = 1; m<= 5; m++){
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp*(ue[2+k-2,m] - 4.0d*ue[2+k-1,m] + 6.0d*ue[2+k,m] - 4.0d*ue[2+k+1,m] + ue[2+k+2,m]);
                            }
                        }
                        if (end[3,c] > 0) {
                            for(m = 1; m<= 5; m++){
                                k = cell_size[3,c]-3;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (ue[2+k-2,m] - 4.0d*ue[2+k-1,m] + 6.0d*ue[2+k,m] - 4.0d*ue[2+k+1,m]);
                                k = cell_size[3,c]-2;
                                forcing[c,k,j,i,m] = forcing[c,k,j,i,m] - dssp * (ue[2+k-2,m] - 4.0d*ue[2+k-1,m] + 5.0d*ue[2+k,m]);
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     now change the sign of the forcing function, 
                //---------------------------------------------------------------------
                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                        for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                            for(m = 1; m<= 5; m++){
                                forcing[c,k,j,i,m] = -1.0 * forcing[c,k,j,i,m];
                            }
                        }
                    }
                }
            }
        }

        // End exact_rhs.f

        // Define.f
        public void compute_buffer_size(int dim) {
            int  c, face_size;
            if (ncells == 1) return;
            //---------------------------------------------------------------------
            //     compute the actual sizes of the buffers; note that there is 
            //     always one cell face that doesn't need buffer space, because it 
            //     is at the boundary of the grid
            //---------------------------------------------------------------------

            west_size = 0;
            east_size = 0;
            for(  c = 1; c<= ncells; c++){
                face_size = cell_size[2,c] * cell_size[3,c] * dim * 2;
                if (cell_coord[1,c]!=1) west_size = west_size + face_size;
                if (cell_coord[1,c]!=ncells) east_size = east_size + face_size;
            }

            north_size = 0;
            south_size = 0;
            for(  c = 1; c<= ncells; c++){
                face_size = cell_size[1,c]*cell_size[3,c] * dim * 2;
                if (cell_coord[2,c]!=1) south_size = south_size + face_size;
                if (cell_coord[2,c]!=ncells) north_size = north_size + face_size;
            }

            top_size = 0;
            bottom_size = 0;
            for(  c = 1; c<= ncells; c++){
               face_size = cell_size[1,c] * cell_size[2,c] * dim * 2;
               if (cell_coord[3,c]!=1) bottom_size = bottom_size + face_size;
               if (cell_coord[3,c]!=ncells) top_size = top_size + face_size;   
            }

            start_send_west   = 1;
            start_send_east   = start_send_west   + west_size;
            start_send_south  = start_send_east   + east_size;
            start_send_north  = start_send_south  + south_size;
            start_send_bottom = start_send_north  + north_size;
            start_send_top = start_send_bottom + bottom_size;
            start_recv_west = 1;
            start_recv_east = start_recv_west + west_size;
            start_recv_south = start_recv_east + east_size;
            start_recv_north = start_recv_south + south_size;
            start_recv_bottom = start_recv_north + north_size;
            start_recv_top = start_recv_bottom + bottom_size;
        }
        // End Define.f

        // Adi.f
        public void adi() {
            copy_faces();
            x_solve();
            y_solve();
            z_solve();
            add();
        }

        public void copy_faces() {
            //---------------------------------------------------------------------     
            // This function copies the face values of a variable defined on a set 
            // of cells to the overlap locations of the adjacent sets of cells. 
            // Because a set of cells interfaces in each direction with exactly one 
            // other set, we only need to fill six different buffers. We could try to 
            // overlap communication with computation, by computing
            // some internal values while communicating boundary values, but this
            // adds so much overhead that it's not clearly useful. 
            //---------------------------------------------------------------------
            //int requests[0:11], b_size[0:5], ss[0:5], sr[0:5],statuses[MPI_STATUS_SIZE, 0:11];
            int i, j, k, c, m, p0, p1, p2, p3, p4, p5;
            MPI.Request[] requests = new MPI.Request[12];
            double[][] in_buffer = new double[6][]; 
            double[][] out_buffer = new double[6][];
            int[] b_size = new int[6];
            //int[] ss = new int[5+le];
            //int[] sr = new int[5+le];
            //int[,] statuses= new int[MPI_STATUS_SIZE+1, 12];
            //---------------------------------------------------------------------
            //     exit immediately if there are no faces to be copied           
            //---------------------------------------------------------------------
            if (no_nodes == 1) {
                compute_rhs();
                return;
            }/*
            ss[0] = start_send_east;
            ss[1] = start_send_west;
            ss[2] = start_send_north;
            ss[3] = start_send_south;
            ss[4] = start_send_top;
            ss[5] = start_send_bottom;

            sr[0] = start_recv_east;
            sr[1] = start_recv_west;
            sr[2] = start_recv_north;
            sr[3] = start_recv_south;
            sr[4] = start_recv_top;
            sr[5] = start_recv_bottom; */

            b_size[0] = east_size;
            b_size[1] = west_size;
            b_size[2] = north_size;
            b_size[3] = south_size;
            b_size[4] = top_size;
            b_size[5] = bottom_size;

            for (i = 0; i < 6; i++) {
                out_buffer[i] = new double[b_size[i]];
                in_buffer[i] = new double[b_size[i]];
            }
            //---------------------------------------------------------------------
            //     because the difference stencil for the diagonalized scheme is 
            //     orthogonal, we for(not have to perform the staged copying of faces, 
            //     but can send all face information simultaneously to the neighboring 
            //     cells in all directions          
            //---------------------------------------------------------------------
            p0 = 0;
            p1 = 0;
            p2 = 0;
            p3 = 0;
            p4 = 0;
            p5 = 0;
            for (c = 1; c <= ncells; c++) {
                //---------------------------------------------------------------------
                //     fill the buffer to be sent to eastern neighbors [i-dir]
                //---------------------------------------------------------------------
                if (cell_coord[1, c] != ncells) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = cell_size[1, c] - 2; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    out_buffer[0][p0] = u[c, k + 2, j + 2, i + 2, m];//u[m,i,j,k,c]; //out_buffer[ss[0] + p0] 
                                    p0 = p0 + 1;
                                }
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     fill the buffer to be sent to western neighbors 
                //---------------------------------------------------------------------
                if (cell_coord[1, c] != 1) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = 0; i <= 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    out_buffer[1][p1] = u[c, k + 2, j + 2, i + 2, m]; //u[m,i,j,k,c]; //out_buffer[ss[1] + p1] 
                                    p1 = p1 + 1;
                                }
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     fill the buffer to be sent to northern neighbors [j_dir]
                //---------------------------------------------------------------------
                if (cell_coord[2, c] != ncells) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = cell_size[2, c] - 2; j <= cell_size[2, c] - 1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    out_buffer[2][p2] = u[c, k + 2, j + 2, i + 2, m];  //u[m,i,j,k,c]; //out_buffer[ss[2] + p2] 
                                    p2 = p2 + 1;
                                }
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     fill the buffer to be sent to southern neighbors 
                //---------------------------------------------------------------------
                if (cell_coord[2, c] != 1) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = 0; j <= 1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    out_buffer[3][p3] = u[c, k + 2, j + 2, i + 2, m]; //u[m,i,j,k,c]; //out_buffer[ss[3] + p3] 
                                    p3 = p3 + 1;
                                }
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     fill the buffer to be sent to top neighbors [k-dir]
                //---------------------------------------------------------------------
                if (cell_coord[3, c] != ncells) {
                    for (k = cell_size[3, c] - 2; k <= cell_size[3, c] - 1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    out_buffer[4][p4] = u[c, k + 2, j + 2, i + 2, m];  //u[m,i,j,k,c]; //out_buffer[ss[4] + p4] 
                                    p4 = p4 + 1;
                                }
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     fill the buffer to be sent to bottom neighbors
                //---------------------------------------------------------------------
                if (cell_coord[3, c] != 1) {
                    for (k = 0; k <= 1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    out_buffer[5][p5] = u[c, k + 2, j + 2, i + 2, m]; //u[m,i,j,k,c]; //out_buffer[ss[5] + p5] = 
                                    p5 = p5 + 1;
                                }
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     cell loop
                //---------------------------------------------------------------------
            }
            {   MPI.RequestList requestList = new MPI.RequestList();

                requests[0] = comm_rhs.ImmediateReceive<double>(  successor[1], WEST,   in_buffer[0]);//call mpi_irecv[in_buffer[sr[0]], b_size[0], dp_type, successor[1], WEST,  comm_rhs, requests[0], error];
                requests[1] = comm_rhs.ImmediateReceive<double>(predecessor[1], EAST,   in_buffer[1]);//call mpi_irecv[in_buffer[sr[1]], b_size[1], dp_type, predecessor[1], EAST, comm_rhs, requests[1], error];
                requests[2] = comm_rhs.ImmediateReceive<double>(  successor[2], SOUTH,  in_buffer[2]);//call mpi_irecv[in_buffer[sr[2]], b_size[2], dp_type, successor[2], SOUTH, comm_rhs, requests[2], error];
                requests[3] = comm_rhs.ImmediateReceive<double>(predecessor[2], NORTH,  in_buffer[3]);//call mpi_irecv[in_buffer[sr[3]], b_size[3], dp_type, predecessor[2], NORTH, comm_rhs, requests[3], error];
                requests[4] = comm_rhs.ImmediateReceive<double>(  successor[3], BOTTOM, in_buffer[4]);//call mpi_irecv[in_buffer[sr[4]], b_size[4], dp_type, successor[3], BOTTOM,  comm_rhs, requests[4], error];
                requests[5] = comm_rhs.ImmediateReceive<double>(predecessor[3], TOP,    in_buffer[5]);//call mpi_irecv[in_buffer[sr[5]], b_size[5], dp_type, predecessor[3], TOP, comm_rhs, requests[5], error];

                requests[6]  = comm_rhs.ImmediateSend<double>(out_buffer[0],   successor[1], EAST);//call mpi_isend[out_buffer[ss[0]], b_size[0], dp_type, successor[1],   EAST, comm_rhs, requests[6], error];
                requests[7]  = comm_rhs.ImmediateSend<double>(out_buffer[1], predecessor[1], WEST);//call mpi_isend[out_buffer[ss[1]], b_size[1], dp_type, predecessor[1], WEST, comm_rhs, requests[7], error];
                requests[8]  = comm_rhs.ImmediateSend<double>(out_buffer[2],   successor[2], NORTH);//call mpi_isend[out_buffer[ss[2]], b_size[2], dp_type,successor[2],   NORTH, comm_rhs, requests[8], error];
                requests[9]  = comm_rhs.ImmediateSend<double>(out_buffer[3], predecessor[2], SOUTH);//call mpi_isend[out_buffer[ss[3]], b_size[3], dp_type,predecessor[2], SOUTH, comm_rhs, requests[9], error];
                requests[10] = comm_rhs.ImmediateSend<double>(out_buffer[4],   successor[3], TOP);//call mpi_isend[out_buffer[ss[4]], b_size[4], dp_type,successor[3],   TOP,   comm_rhs,   requests[10], error];
                requests[11] = comm_rhs.ImmediateSend<double>(out_buffer[5], predecessor[3], BOTTOM);//call mpi_isend[out_buffer[ss[5]], b_size[5], dp_type,predecessor[3], BOTTOM, comm_rhs,requests[11], error];

                foreach (MPI.Request request in requests) {
                    requestList.Add(request);
                }
                requestList.WaitAll();//call mpi_waitall[12, requests, statuses, error];
            }
            //---------------------------------------------------------------------
            //     unpack the data that has just been received;             
            //---------------------------------------------------------------------
            p0 = 0;
            p1 = 0;
            p2 = 0;
            p3 = 0;
            p4 = 0;
            p5 = 0;

            for (c = 1; c <= ncells; c++) {
                if (cell_coord[1, c] != 1) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = -2; i <= -1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    u[c, k + 2, j + 2, i + 2, m] = in_buffer[1][p0]; //in_buffer[sr[1] + p0];
                                    p0 = p0 + 1;
                                }
                            }
                        }
                    }
                }

                if (cell_coord[1, c] != ncells) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = cell_size[1, c]; i <= cell_size[1, c] + 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    u[c, k + 2, j + 2, i + 2, m] = in_buffer[0][p1];//in_buffer[sr[0] + p1];
                                    p1 = p1 + 1;
                                }
                            }
                        }
                    }
                }

                if (cell_coord[2, c] != 1) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = -2; j <= -1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    u[c, k + 2, j + 2, i + 2, m] = in_buffer[3][p2];//in_buffer[sr[3] + p2];
                                    p2 = p2 + 1;
                                }
                            }
                        }
                    }
                }

                if (cell_coord[2, c] != ncells) {
                    for (k = 0; k <= cell_size[3, c] - 1; k++) {
                        for (j = cell_size[2, c]; j <= cell_size[2, c] + 1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    u[c, k + 2, j + 2, i + 2, m] = in_buffer[2][p3];//in_buffer[sr[2] + p3];
                                    p3 = p3 + 1;
                                }
                            }
                        }
                    }
                }

                if (cell_coord[3, c] != 1) {
                    for (k = -2; k <= -1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    u[c, k + 2, j + 2, i + 2, m] = in_buffer[5][p4];//in_buffer[sr[5] + p4];
                                    p4 = p4 + 1;
                                }
                            }
                        }
                    }
                }

                if (cell_coord[3, c] != ncells) {
                    for (k = cell_size[3, c]; k <= cell_size[3, c] + 1; k++) {
                        for (j = 0; j <= cell_size[2, c] - 1; j++) {
                            for (i = 0; i <= cell_size[1, c] - 1; i++) {
                                for (m = 1; m <= 5; m++) {
                                    u[c, k + 2, j + 2, i + 2, m] = in_buffer[4][p5];//in_buffer[sr[4] + p5];
                                    p5 = p5 + 1;
                                }
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     cells loop
                //---------------------------------------------------------------------
            }
            //---------------------------------------------------------------------
            //     for(the rest of the rhs that uses the copied face values          
            //---------------------------------------------------------------------
            compute_rhs();
        }

              //rhs.f
        public void compute_rhs() {
            int c, i, j, k, m;
            double rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;
            //---------------------------------------------------------------------
            //     loop over all cells owned by this node                           
            //---------------------------------------------------------------------
            for(c = 1; c<= ncells; c++){
                //---------------------------------------------------------------------
                //     compute the reciprocal of density, and the kinetic energy, 
                //     and the speed of sound.
                //---------------------------------------------------------------------
                //Add in index arrays: square(+0, +1, +1 + 1); qs(+0, +1, +1, +1); rho_i(+0,+1,+1,+1); u(+0,+2,+2,+2,+0); 
                for (k = -1; k <= cell_size[3, c]; k++) {
                    for (j = -1; j <= cell_size[2, c]; j++) {
                        for (i = -1; i <= cell_size[1, c]; i++) {
                            rho_inv = 1.0d / u[c, k+2, j+2, i+2, 1]; //u[1,i,j,k,c]
                            rho_i[c, k+1, j+1, i+1] = rho_inv;                        //rho_i[i, j, k, c] = rho_inv;
                            us[c, k+1, j+1, i+1] = u[c, k+2, j+2, i+2, 2] * rho_inv;        //us[i, j, k, c] = u[2, i, j, k, c]
                            vs[c, k+1, j+1, i+1] = u[c, k+2, j+2, i+2, 3] * rho_inv;        //vs[i, j, k, c] = u[3, i, j, k, c] 
                            ws[c, k+1, j+1, i+1] = u[c, k+2, j+2, i+2, 4] * rho_inv;        //ws[i, j, k, c] = u[4, i, j, k, c] 
                            square[c, k+1, j+1, i+1] = 0.5d * (                       /*square[i, j, k, c]*/
                                u[c, k+2, j+2, i+2, 2] * u[c, k+2, j+2, i+2, 2] +           /*u[2, i, j, k, c] * u[2, i, j, k, c]*/
                                u[c, k+2, j+2, i+2, 3] * u[c, k+2, j+2, i+2, 3] +           /*u[3, i, j, k, c] * u[3, i, j, k, c]*/
                                u[c, k+2, j+2, i+2, 4] * u[c, k+2, j+2, i+2, 4]) * rho_inv; /*u[4, i, j, k, c] * u[4, i, j, k, c])*/
                            qs[c, k+1, j+1, i+1] = square[c, k+1, j+1, i+1] * rho_inv;      /*qs[i, j, k, c] = square[i, j, k, c]*/
                        }
                    }
                }
                //---------------------------------------------------------------------
                // copy the exact forcing term to the right hand side;  because 
                // this forcing term is known, we can store it on the whole of every 
                // cell,  including the boundary                   
                //---------------------------------------------------------------------
                for(k = 0; k<= cell_size[3,c]-1; k++){
                    for(j = 0; j<= cell_size[2,c]-1; j++){
                        for(i = 0; i<= cell_size[1,c]-1; i++){
                            for(m = 1; m<= 5; m++){  //rhs(+0, +1, +1, +1, +0);
                                rhs[c,k+1,j+1,i+1,m] = forcing[c,k,j,i,m];  //rhs[m,i,j,k,c] = forcing[m,i,j,k,c];
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     compute xi-direction fluxes 
                //---------------------------------------------------------------------
                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                        for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){  //us(+0, +1, +1, +1);
                            uijk = us[c,k+1,j+1,i+1];   //uijk = us[i,j,k,c]; 
                            up1  = us[c,k+1,j+1,i+2]; //up1  = us[i+1,j,k,c];
                            um1  = us[c,k+1,j+1,i]; //um1  = us[i-1,j,k,c];
                            //u(+0, +2, +2, +2, +0); rhs(+0, +1, +1, +1, +0); square(+0, +1, +1 + 1); /*vs(+0, +1, +1, +1);*/
                            rhs[c,k+1,j+1,i+1,1] = rhs[c,k+1,j+1,i+1,1] + dx1tx1*(u[c,k+2,j+2,i+1+2,1] - /*rhs[1,i,j,k,c]=rhs[1,i,j,k,c]+dx1tx1*(u[1,i+1,j,k,c]*/
                                             2.0d*u[c,k+2,j+2,i+2,1] + u[c,k+2,j+2,i-1+2,1]) - /*2.0d*u[1,i,j,k,c] + u[1,i-1,j,k,c])*/
                                             tx2*(u[c,k+2,j+2,i+1+2,2] - u[c,k+2,j+2,i-1+2,2]); /*tx2*(u[2,i+1,j,k,c] - u[2,i-1,j,k,c])*/
                            rhs[c,k+1,j+1,i+1,2] = rhs[c,k+1,j+1,i+1,2] + dx2tx1 * /*rhs[2, i, j, k, c] = rhs[2, i, j, k, c] + dx2tx1 */
                                (u[c,k+2,j+2,i+1+2,2] - 2.0d*u[c,k+2,j+2,i+2,2] +  /*(u[2, i + 1, j, k, c] - 2.0d * u[2, i, j, k, c] +*/
                                u[c,k+2,j+2,i-1+2,2]) + /* u[2, i - 1, j, k, c]) +*/
                                xxcon2*con43 * (up1 - 2.0d*uijk + um1) - /*xxcon2 * con43 * (up1 - 2.0d * uijk + um1) -*/
                                tx2 * (u[c,k+2,j+2,i+1+2,2]*up1 -  /*tx2 * (u[2, i + 1, j, k, c] * up1 -*/
                                u[c,k+2,j+2,i-1+2,2]*um1 + /*u[2, i - 1, j, k, c] * um1 +*/
                                (u[c,k+2,j+2,i+1+2,5]- square[c,k+1,j+1,i+1+1]- /*(u[5, i + 1, j, k, c] - square[i + 1, j, k, c] -*/
                                u[c,k+2,j+2,i-1+2,5]+ square[c,k+1,j+1,i-1+1])* /*u[5, i - 1, j, k, c] + square[i - 1, j, k, c]) **/
                                c2);
                            rhs[c,k+1,j+1,i+1,3] = rhs[c,k+1,j+1,i+1,3] + dx3tx1 * /*rhs[3,i,j,k,c] = rhs[3,i,j,k,c] + dx3tx1 * */
                                (u[c,k+2,j+2,i+1+2,3] - 2.0d*u[c,k+2,j+2,i+2,3] +  /*[u[3,i+1,j,k,c] - 2.0d0*u[3,i,j,k,c] +*/
                                u[c,k+2,j+2,i-1+2,3]) + /*u[3,i-1,j,k,c]] +*/
                                xxcon2 * (vs[c,k+1,j+1,i+1+1] - 2.0d*vs[c,k+1,j+1,i+1] + /*xxcon2 * [vs[i+1,j,k,c] - 2.0d0*vs[i,j,k,c] +*/
                                vs[c,k+1,j+1,i-1+1]) - /*vs[i-1,j,k,c]] -*/
                                tx2 * (u[c,k+2,j+2,i+1+2,3]*up1 - /* tx2 * [u[3,i+1,j,k,c]*up1 - */
                                u[c,k+2,j+2,i-1+2,3]*um1); /*u[3,i-1,j,k,c]*um1];*/

                            rhs[c,k+1,j+1,i+1,4] = rhs[c,k+1,j+1,i+1,4] + dx4tx1 *  /*rhs[4,i,j,k,c] = rhs[4,i,j,k,c] + dx4tx1 * */
                                (u[c,k+2,j+2,i+1+2,4] - 2.0d*u[c,k+2,j+2,i+2,4] + /*[u[4,i+1,j,k,c] - 2.0d0*u[4,i,j,k,c] +*/
                                u[c,k+2,j+2,i-1+2,4]) +  /*u[4,i-1,j,k,c]] +*/
                                xxcon2 * (ws[c,k+1,j+1,i+1+1] - 2.0d*ws[c,k+1,j+1,i+1] + /*xxcon2 * [ws[i+1,j,k,c] - 2.0d0*ws[i,j,k,c] +*/
                                ws[c,k+1,j+1,i-1+1]) - /*ws[i-1,j,k,c]] -*/
                                tx2 * (u[c,k+2,j+2,i+1+2,4]*up1 - /*tx2 * [u[4,i+1,j,k,c]*up1 - */
                                u[c,k+2,j+2,i-1+2,4]*um1); /*u[4,i-1,j,k,c]*um1];*/

                            rhs[c,k+1,j+1,i+1,5] = rhs[c,k+1,j+1,i+1,5] + dx5tx1 *            /*rhs[5,i,j,k,c] = rhs[5,i,j,k,c] + dx5tx1 * */
                                (u[c,k+2,j+2,i+1+2,5] - 2.0d*u[c,k+2,j+2,i+2,5] +             /*[u[5,i+1,j,k,c] - 2.0d0*u[5,i,j,k,c] +*/
                                u[c,k+2,j+2,i-1+2,5]) +                                       /*u[5,i-1,j,k,c]] +*/
                                xxcon3 * (qs[c,k+1,j+1,i+1+1] - 2.0d*qs[c,k+1,j+1,i+1] +      /*xxcon3 * [qs[i+1,j,k,c] - 2.0d0*qs[i,j,k,c] +*/
                                qs[c,k+1,j+1,i-1+1]) +                                        /*qs[i-1,j,k,c]] +*/ 
                                xxcon4 * (up1*up1 -2.0d*uijk*uijk +                           /*xxcon4 * [up1*up1 -       2.0d0*uijk*uijk + */
                                um1*um1) +                                                    /*um1*um1] +*/
                                xxcon5 * (u[c,k+2,j+2,i+1+2,5]*rho_i[c,k+1,j+1,i+1+1] -       /*xxcon5 * [u[5,i+1,j,k,c]*rho_i[i+1,j,k,c] - */
                                2.0d*u[c,k+2,j+2,i+2,5]*rho_i[c,k+1,j+1,i+1] +                /*2.0d0*u[5,i,j,k,c]*rho_i[i,j,k,c] +*/
                                u[c,k+2,j+2,i-1+2,5]*rho_i[c,k+1,j+1,i-1+1]) -                /*u[5,i-1,j,k,c]*rho_i[i-1,j,k,c]] -*/
                                tx2 * ( (c1*u[c,k+2,j+2,i+1+2,5] -                            /*tx2 * [ [c1*u[5,i+1,j,k,c] - */
                                c2*square[c,k+1,j+1,i+1+1])*up1 -                             /*c2*square[i+1,j,k,c]]*up1 -*/
                                (c1*u[c,k+2,j+2,i-1+2,5] -                                    /*[c1*u[5,i-1,j,k,c] - */
                                c2*square[c,k+1,j+1,i-1+1])*um1 );                            /*c2*square[i-1,j,k,c]]*um1 ];*/
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     add fourth order xi-direction dissipation               
                //---------------------------------------------------------------------
                if (start[1,c] > 0) {
                    for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                        for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                            i = 1;
                            for(m = 1; m<= 5; m++){
                                rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m]- dssp*(5.0d*u[c,k+2,j+2,i+2,m] - 
                                                       4.0d*u[c,k+2,j+2,i+1+2,m] + u[c,k+2,j+2,i+2+2,m]);
                            }
                            i = 2;
                            for(m = 1; m<= 5; m++){
                                rhs[c,k+1,j+1,i+1,m]=rhs[c,k+1,j+1,i+1,m]-dssp*(-4.0d*u[c,k+2,j+2,i-1+2,m]+6.0d*u[c,k+2,j+2,i+2,m]-
                                                     4.0d*u[c,k+2,j+2,i+1+2,m]+u[c,k+2,j+2,i+2+2,m]);
                            }
                        }
                    }
                }

                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                        for(i = 3*start[1,c]; i<= cell_size[1,c]-3*end[1,c]-1; i++){
                            for(m = 1; m<= 5; m++){
                                rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - dssp * 
                                   (u[c,k+2,j+2,i-2+2,m] - 4.0d*u[c,k+2,j+2,i-1+2,m] + 
                                   6.0*u[c,k+2,j+2,i+2,m] - 4.0d*u[c,k+2,j+2,i+1+2,m] + 
                                   u[c,k+2,j+2,i+2+2,m]);
                            }
                        }
                    }
                }
                if (end[1,c] > 0) {
                    for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                        for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                            i = cell_size[1,c]-3;
                            for(m = 1; m<= 5; m++){
                                rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                                dssp * (u[c,k+2,j+2,i-2+2,m] - 4.0d*u[c,k+2,j+2,i-1+2,m] + 6.0d*u[c,k+2,j+2,i+2,m] - 4.0d*u[c,k+2,j+2,i+1+2,m]);
                            }
                            i = cell_size[1,c]-2;
                            for(m = 1; m<= 5; m++){
                                rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - dssp * (u[c,k+2,j+2,i-2+2,m] - 4.0*u[c,k+2,j+2,i-1+2,m] + 5.0*u[c,k+2,j+2,i+2,m]);
                            }
                        }
                    }
                }
                //---------------------------------------------------------------------
                //     compute eta-direction fluxes 
                //---------------------------------------------------------------------
                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                        for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                            vijk = vs[c,k+1,j+1,i+1];   //vijk = vs[i,j,k,c];
                            vp1  = vs[c,k+1,j+1+1,i+1]; //vp1  = vs[i,j+1,k,c];
                            vm1  = vs[c,k+1,j-1+1,i+1]; //vm1  = vs[i,j-1,k,c];

                            rhs[c,k+1,j+1,i+1,1] = rhs[c,k+1,j+1,i+1,1] + dy1ty1 *       /*rhs[1,i,j,k,c] = rhs[1,i,j,k,c] + dy1ty1 */
                                (u[c,k+2,j+1+2,i+2,1] - 2.0d*u[c,k+2,j+2,i+2,1] +        /*(u[1,i,j+1,k,c] - 2.0d*u[1,i,j,k,c] + */
                                u[c,k+2,j-1+2,i+2,1]) -                                  /*u[1,i,j-1,k,c]) -*/
                                ty2 * (u[c,k+2,j+1+2,i+2,3] - u[c,k+2,j-1+2,i+2,3]);     /*ty2 * (u[3,i,j+1,k,c] - u[3,i,j-1,k,c]);*/

                            rhs[c,k+1,j+1,i+1,2] = rhs[c,k+1,j+1,i+1,2] + dy2ty1 *       /*rhs[2,i,j,k,c] = rhs[2,i,j,k,c] + dy2ty1 * */       
                                (u[c,k+2,j+1+2,i+2,2] - 2.0d*u[c,k+2,j+2,i+2,2] +        /*(u[2,i,j+1,k,c] - 2.0d*u[2,i,j,k,c] + */      
                                u[c,k+2,j-1+2,i+2,2]) +                                  /*u[2,i,j-1,k,c]) + */                         
                                yycon2 * (us[c,k+1,j+1+1,i+1] - 2.0d*us[c,k+1,j+1,i+1] + /*yycon2 * (us[i,j+1,k,c] - 2.0d*us[i,j,k,c] + */                      
                                us[c,k+1,j-1+1,i+1]) -                                   /*us[i,j-1,k,c]) -   */    
                                ty2 * (u[c,k+2,j+1+2,i+2,2]*vp1 -                        /*ty2 * (u[2,i,j+1,k,c]*vp1 - */     
                                u[c,k+2,j-1+2,i+2,2]*vm1);                               /*u[2,i,j-1,k,c]*vm1);*/  

                            rhs[c,k+1,j+1,i+1,3] = rhs[c,k+1,j+1,i+1,3] + dy3ty1 *       /*rhs[3,i,j,k,c] = rhs[3,i,j,k,c] + dy3ty1 * */                          
                                (u[c,k+2,j+1+2,i+2,3] - 2.0d*u[c,k+2,j+2,i+2,3] +        /*(u[3,i,j+1,k,c] - 2.0d*u[3,i,j,k,c] +*/                           
                                u[c,k+2,j-1+2,i+2,3]) +                                  /*u[3,i,j-1,k,c]) +*/                 
                                yycon2*con43 * (vp1 - 2.0d*vijk + vm1) -                 /*yycon2*con43 * (vp1 - 2.0d*vijk + vm1) -*/                   
                                ty2 * (u[c,k+2,j+1+2,i+2,3]*vp1 -                        /*ty2 * (u[3,i,j+1,k,c]*vp1 -*/                   
                                u[c,k+2,j-1+2,i+2,3]*vm1 +                               /*u[3,i,j-1,k,c]*vm1 +*/                        
                                (u[c,k+2,j+1+2,i+2,5] - square[c,k+1,j+1+1,i+1] -        /*(u[5,i,j+1,k,c] - square[i,j+1,k,c] - */                    
                                u[c,k+2,j-1+2,i+2,5] + square[c,k+1,j-1+1,i+1])          /*u[5,i,j-1,k,c] + square[i,j-1,k,c])*/                  
                                *c2);                                                   

                            rhs[c,k+1,j+1,i+1,4] = rhs[c,k+1,j+1,i+1,4] + dy4ty1 *       /*rhs[4,i,j,k,c] = rhs[4,i,j,k,c] + dy4ty1 * */                          
                                (u[c,k+2,j+1+2,i+2,4] - 2.0d*u[c,k+2,j+2,i+2,4] +        /*(u[4,i,j+1,k,c] - 2.0d*u[4,i,j,k,c] + */                        
                                u[c,k+2,j-1+2,i+2,4]) +                                  /*u[4,i,j-1,k,c]) + */                   
                                yycon2 * (ws[c,k+1,j+1+1,i+1] - 2.0d*ws[c,k+1,j+1,i+1] + /*yycon2 * (ws[i,j+1,k,c] - 2.0d*ws[i,j,k,c] + */                                  
                                ws[c,k+1,j-1+1,i+1]) -                                   /*ws[i,j-1,k,c]) */                         
                                ty2 * (u[c,k+2,j+1+2,i+2,4]*vp1 -                        /*ty2 * (u[4,i,j+1,k,c]*vp1 */                              
                                u[c,k+2,j-1+2,i+2,4]*vm1);                               /*u[4,i,j-1,k,c]*vm1)*/    
          
                            rhs[c,k+1,j+1,i+1,5] = rhs[c,k+1,j+1,i+1,5] + dy5ty1 *       /*rhs[5,i,j,k,c] = rhs[5,i,j,k,c] + dy5ty1 */                                  
                                (u[c,k+2,j+1+2,i+2,5] - 2.0d*u[c,k+2,j+2,i+2,5] +        /*(u[5,i,j+1,k,c] - 2.0d*u[5,i,j,k,c] */                                    
                                u[c,k+2,j-1+2,i+2,5]) +                                  /*u[5,i,j-1,k,c]) */                                
                                yycon3 * (qs[c,k+1,j+1+1,i+1] - 2.0d*qs[c,k+1,j+1,i+1] + /*yycon3 * (qs[i,j+1,k,c] - 2.0d*qs[i,j,k,c] */                     
                                qs[c,k+1,j-1+1,i+1]) +                                   /*qs[i,j-1,k,c]) */                  
                                yycon4 * (vp1*vp1 - 2.0d*vijk*vijk +                                         
                                vm1*vm1) +                                                             
                                yycon5 * (u[c,k+2,j+1+2,i+2,5]*rho_i[c,k+1,j+1+1,i+1] -  /*yycon5 * (u[5,i,j+1,k,c]*rho_i[i,j+1,k,c] */                      
                                2.0d*u[c,k+2,j+2,i+2,5]*rho_i[c,k+1,j+1,i+1] +           /*2.0d*u[5,i,j,k,c]*rho_i[i,j,k,c] */                         
                                u[c,k+2,j-1+2,i+2,5]*rho_i[c,k+1,j-1+1,i+1]) -           /*u[5,i,j-1,k,c]*rho_i[i,j-1,k,c]) */                         
                                ty2 * ((c1*u[c,k+2,j+1+2,i+2,5] -                        /*ty2 * ((c1*u[5,i,j+1,k,c] */                    
                                c2*square[c,k+1,j+1+1,i+1]) * vp1 -                      /*c2*square[i,j+1,k,c]) * vp1 */                            
                                (c1*u[c,k+2,j-1+2,i+2,5] -                               /*(c1*u[5,i,j-1,k,c] */                             
                                c2*square[c,k+1,j-1+1,i+1]) * vm1);                      /*c2*square[i,j-1,k,c]) * vm1)*/                        
                     }
                }
            }
            //---------------------------------------------------------------------
            //     add fourth order eta-direction dissipation         
            //---------------------------------------------------------------------
            if (start[2,c] > 0) {
                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    j = 1;
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){ //rhs[m,i,j,k,c]=rhs[m,i,j,k,c]-dssp*(5.0d*u[m,i,j,k,c]-4.0d*u[m,i,j+1,k,c]+u[m,i,j+2,k,c]);
                            rhs[c,k+1,j+1,i+1,m]=rhs[c,k+1,j+1,i+1,m]-dssp*(5.0d*u[c,k+2,j+2,i+2,m]-4.0d*u[c,k+2,j+1+2,i+2,m]+u[c,k+2,j+2+2,i+2,m]);
                        }
                    }
                    j = 2;
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                            dssp * (-4.0d*u[c,k+2,j-1+2,i+2,m] + 6.0d*u[c,k+2,j+2,i+2,m] - 4.0d*u[c,k+2,j+1+2,i+2,m] + u[c,k+2,j+2+2,i+2,m]);
                        }
                    }
                }
            }
            for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                for(j = 3*start[2,c]; j<= cell_size[2,c]-3*end[2,c]-1; j++){
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - dssp*(u[c,k+2,j-2+2,i+2,m] - 4.0d*u[c,k+2,j-1+2,i+2,m] + 
                                6.0*u[c,k+2,j+2,i+2,m] - 4.0d*u[c,k+2,j+1+2,i+2,m] + u[c,k+2,j+2+2,i+2,m]);
                        }
                    }
                }
            }
            if (end[2,c] > 0) {
                for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                    j = cell_size[2,c]-3;
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){//rhs[m,i,j,k,c]=rhs[m,i,j,k,c]-dssp*(u[m,i,j-2,k,c]-4.0d*u[m,i,j-1,k,c]+6.0d*u[m,i,j,k,c]-4.0d*u[m,i,j+1,k,c]);
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                            dssp*(u[c,k+2,j-2+2,i+2,m]-4.0d*u[c,k+2,j-1+2,i+2,m]+6.0d*u[c,k+2,j+2,i+2,m]-4.0d*u[c,k+2,j+1+2,i+2,m]);
                        }
                    }
                    j = cell_size[2,c]-2;
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){//rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - dssp*(u[m,i,j-2,k,c]-4.0*u[m,i,j-1,k,c]+5.0*u[m,i,j,k,c] );
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                            dssp*(u[c,k+2,j-2+2,i+2,m]-4.0*u[c,k+2,j-1+2,i+2,m]+5.0*u[c,k+2,j+2,i+2,m]);
                        }
                    }
                }
            }
            //---------------------------------------------------------------------
            //     compute zeta-direction fluxes 
            //---------------------------------------------------------------------
            for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        wijk = ws[c,k+1  ,j+1,i+1]; //wijk = ws[i,j,k,c];
                        wp1  = ws[c,k+1+1,j+1,i+1]; //wp1  = ws[i,j,k+1,c];
                        wm1  = ws[c,k-1+1,j+1,i+1]; //wm1  = ws[i,j,k-1,c];
                        rhs[c,k+1,j+1,i+1,1] = rhs[c,k+1,j+1,i+1,1] + dz1tz1 *        /*rhs[1,i,j,k,c] = rhs[1,i,j,k,c] + dz1tz1 * */
                            (u[c,k+1+2,j+2,i+2,1] - 2.0d*u[c,k+2,j+2,i+2,1] +         /*(u[1,i,j,k+1,c] - 2.0d*u[1,i,j,k,c] + */
                            u[c,k-1+2,j+2,i+2,1]) -                                   /*u[1,i,j,k-1,c]) - */
                            tz2 * (u[c,k+1+2,j+2,i+2,4] - u[c,k-1+2,j+2,i+2,4]);      /*tz2 * (u[4,i,j,k+1,c] - u[4,i,j,k-1,c]);*/

                        rhs[c,k+1,j+1,i+1,2] = rhs[c,k+1,j+1,i+1,2] + dz2tz1 *        /*rhs[2,i,j,k,c] = rhs[2,i,j,k,c] + dz2tz1 **/
                            (u[c,k+1+2,j+2,i+2,2] - 2.0d*u[c,k+2,j+2,i+2,2] +         /*(u[2,i,j,k+1,c] - 2.0d*u[2,i,j,k,c] */
                            u[c,k-1+2,j+2,i+2,2]) +                                   /*u[2,i,j,k-1,c]) */
                            zzcon2 * (us[c,k+1+1,j+1,i+1] - 2.0d*us[c,k+1,j+1,i+1] +  /*zzcon2 * (us[i,j,k+1,c] - 2.0d*us[i,j,k,c] */
                            us[c,k-1+1,j+1,i+1]) -                                    /*us[i,j,k-1,c]) */
                            tz2 * (u[c,k+1+2,j+2,i+2,2]*wp1 -                         /*tz2 * (u[2,i,j,k+1,c]*wp1 */
                            u[c,k-1+2,j+2,i+2,2]*wm1);                                /*u[2,i,j,k-1,c]*wm1)*/

                        rhs[c,k+1,j+1,i+1,3] = rhs[c,k+1,j+1,i+1,3] + dz3tz1 *        /*rhs[3,i,j,k,c] = rhs[3,i,j,k,c] + dz3tz1 */
                            (u[c,k+1+2,j+2,i+2,3] - 2.0d*u[c,k+2,j+2,i+2,3] +         /*(u[3,i,j,k+1,c] - 2.0d*u[3,i,j,k,c] */
                            u[c,k-1+2,j+2,i+2,3]) +                                   /*u[3,i,j,k-1,c]) */
                            zzcon2 * (vs[c,k+1+1,j+1,i+1] - 2.0d*vs[c,k+1,j+1,i+1] +  /*zzcon2 * (vs[i,j,k+1,c] - 2.0d*vs[i,j,k,c] */
                            vs[c,k-1+1,j+1,i+1]) -                                    /*vs[i,j,k-1,c]) */
                            tz2 * (u[c,k+1+2,j+2,i+2,3]*wp1 -                         /*tz2 * (u[3,i,j,k+1,c]*wp1 */
                            u[c,k-1+2,j+2,i+2,3]*wm1);                                /*u[3,i,j,k-1,c]*wm1)*/

                        rhs[c,k+1,j+1,i+1,4] = rhs[c,k+1,j+1,i+1,4] + dz4tz1 *        /*rhs[4,i,j,k,c] = rhs[4,i,j,k,c] */
                            (u[c,k+1+2,j+2,i+2,4] - 2.0d*u[c,k+2,j+2,i+2,4] +         /*(u[4,i,j,k+1,c] - 2.0d*u[4,i,j,k,c] */
                            u[c,k-1+2,j+2,i+2,4]) +                                   /*u[4,i,j,k-1,c]) */
                            zzcon2*con43 * (wp1 - 2.0d*wijk + wm1) -            
                            tz2 * (u[c,k+1+2,j+2,i+2,4]*wp1 -                         /*tz2 * (u[4,i,j,k+1,c]*wp1 */
                            u[c,k-1+2,j+2,i+2,4]*wm1 +                                /*u[4,i,j,k-1,c]*wm1 */
                            (u[c,k+1+2,j+2,i+2,5] - square[c,k+1+1,j+1,i+1] -         /*(u[5,i,j,k+1,c] - square[i,j,k+1,c] */
                            u[c,k-1+2,j+2,i+2,5] + square[c,k-1+1,j+1,i+1])           /*u[5,i,j,k-1,c] + square[i,j,k-1,c]*/
                            *c2);

                        rhs[c,k+1,j+1,i+1,5] = rhs[c,k+1,j+1,i+1,5] + dz5tz1 *              /*rhs[5,i,j,k,c] = rhs[5,i,j,k,c] + */
                            (u[c,k+1+2,j+2,i+2,5] - 2.0d*u[c,k+2,j+2,i+2,5] +               /*(u[5,i,j,k+1,c] - 2.0d*u[5,i,j,k,c] */
                            u[c,k-1+2,j+2,i+2,5]) +                                   /*u[5,i,j,k-1,c]) */
                            zzcon3 * (qs[c,k+1+1,j+1,i+1] - 2.0d*qs[c,k+1,j+1,i+1] +        /*zzcon3 * (qs[i,j,k+1,c] - 2.0d*qs[i,j,k,c] */
                            qs[c,k-1+1,j+1,i+1]) +                                    /*qs[i,j,k-1,c]) */
                            zzcon4 * (wp1*wp1 - 2.0d*wijk*wijk +                
                            wm1*wm1) +                                          
                            zzcon5 * (u[c,k+1+2,j+2,i+2,5]*rho_i[c,k+1+1,j+1,i+1] -         /*zzcon5 * (u[5,i,j,k+1,c]*rho_i[i,j,k+1,c] */
                            2.0d*u[c,k+2,j+2,i+2,5]*rho_i[c,k+1,j+1,i+1] +                  /*2.0d*u[5,i,j,k,c]*rho_i[i,j,k,c] */
                            u[c,k-1+2,j+2,i+2,5]*rho_i[c,k-1+1,j+1,i+1]) -                  /*u[5,i,j,k-1,c]*rho_i[i,j,k-1,c]) */
                            tz2*((c1*u[c,k+1+2,j+2,i+2,5] -                        /*tz2 * ( (c1*u[5,i,j,k+1,c] */
                            c2*square[c,k+1+1,j+1,i+1])*wp1 -                         /*c2*square[i,j,k+1,c])*wp1 */
                            (c1*u[c,k-1+2,j+2,i+2,5] -                                /*(c1*u[5,i,j,k-1,c] */
                            c2*square[c,k-1+1,j+1,i+1])*wm1);                         /*c2*square[i,j,k-1,c])*wm1)*/
                    }
                }
            }
            //---------------------------------------------------------------------
            //     add fourth order zeta-direction dissipation                
            //---------------------------------------------------------------------
            if (start[3,c] > 0) {
                k = 1;
                for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){ //rhs[m,i,j,k,c] = rhs[m,i,j,k,c]-dssp*(5.0d*u[m,i,j,k,c]-4.0d*u[m,i,j,k+1,c] + u[m,i,j,k+2,c]);
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m]- 
                            dssp * (5.0d*u[c,k+2,j+2,i+2,m] - 4.0d*u[c,k+1+2,j+2,i+2,m] + u[c,k+2+2,j+2,i+2,m]);
                        }
                    }
                }
                k = 2;
                for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){//rhs[m,i,j,k,c]=rhs[m,i,j,k,c]-dssp*(-4.0d*u[m,i,j,k-1,c]+6.0d*u[m,i,j,k,c]-4.0d*u[m,i,j,k+1,c]+u[m,i,j,k+2,c]);
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                            dssp * (-4.0d*u[c,k-1+2,j+2,i+2,m] + 6.0d*u[c,k+2,j+2,i+2,m] - 4.0d*u[c,k+1+2,j+2,i+2,m] + u[c,k+2+2,j+2,i+2,m]);
                        }
                    }
                }
            }
            for(k = 3*start[3,c]; k<= cell_size[3,c]-3*end[3,c]-1; k++){
                for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    for(i = start[1,c]; i<=cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                            dssp*(u[c,k-2+2,j+2,i+2,m]-4.0d*u[c,k-1+2,j+2,i+2,m]+
                            6.0*u[c,k+2,j+2,i+2,m]-4.0d*u[c,k+1+2,j+2,i+2,m]+u[c,k+2+2,j+2,i+2,m]);
                        }
                    }
                }
            }

            if (end[3,c] > 0) {
                k = cell_size[3,c]-3;
                for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                            dssp*(u[c,k-2+2,j+2,i+2,m] - 4.0d*u[c,k-1+2,j+2,i+2,m] + 6.0d*u[c,k+2,j+2,i+2,m] - 4.0d*u[c,k+1+2,j+2,i+2,m]);
                        }
                    }
                }
                k = cell_size[3,c]-2;
                for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - 
                            dssp*(u[c,k-2+2,j+2,i+2,m] - 4.0*u[c,k-1+2,j+2,i+2,m] + 5.0*u[c,k+2,j+2,i+2,m]);
                        }
                    }
                }
            }
            for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                    for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){
                            rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] * dt;
                        }
                    }
                }
            }
        }
    }
              //End rhs.f
        // End adi.f

        // x_solve.f
        public void x_solve() {
            //---------------------------------------------------------------------
            //     
            //     Performs line solves in X direction by first factoring
            //     the block-tridiagonal matrix into an upper triangular matrix, 
            //     and { performing back substitution to solve for the unknow
            //     vectors of each line.  
            //     
            //     Make sure we treat elements zero to cell_size in the direction
            //     of the sweep.
            //     
            //---------------------------------------------------------------------
            int  c, stage, first, last, isize,jsize,ksize,buffer_size; //r_status[MPI_STATUS_SIZE];
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE);
            MPI.Request[] recv_id = new MPI.Request[1];
            MPI.Request[] send_id = new MPI.Request[1];
            double[] out_buffer_x = new double[buffer_size];
            //istart = 0;
            //---------------------------------------------------------------------
            //     in our terminology stage is the number of the cell in the x-direction
            //     i.e. stage = 1 means the start of the line stage=ncells means end
            //---------------------------------------------------------------------
            for(stage = 1; stage<=ncells; stage++){
               c = slice[1,stage];
               isize = cell_size[1,c] - 1;
               jsize = cell_size[2,c] - 1;
               ksize = cell_size[3,c] - 1;
               //---------------------------------------------------------------------
               //     set last-cell flag
               //---------------------------------------------------------------------
               if (stage == ncells) {
                  last = 1;
               } else {
                  last = 0;
               }
               if (stage == 1) {
                  //---------------------------------------------------------------------
                  //     This is the first cell, so solve without receiving data
                  //---------------------------------------------------------------------
                  first = 1;
                  //          //c            call lhsx[c];
                  x_solve_cell(first,last,c);
               } else {
                   //---------------------------------------------------------------------
                   //     Not the first cell of this line, so receive info from
                   //     processor working on preceeding cell
                   //---------------------------------------------------------------------
                   first = 0;
                   x_receive_solve_info(out_buffer_x, recv_id, c);  //x_receive_solve_info(recv_id,c);
                   //---------------------------------------------------------------------
                   //     overlap computations and communications
                   //---------------------------------------------------------------------
                   //            call lhsx[c]
                   //---------------------------------------------------------------------
                   //     wait for completion
                   //---------------------------------------------------------------------
                   send_id[0].Wait();
                   recv_id[0].Wait();
                   //Fortran: call mpi_wait[send_id,r_status,error]
                   //Fortran: call mpi_wait[recv_id,r_status,error]
                   //---------------------------------------------------------------------
                   //     install C'[istart] and rhs'[istart] to be used in this cell
                   //---------------------------------------------------------------------
                   x_unpack_solve_info(out_buffer_x, c);
                   x_solve_cell(first,last,c);
               }
               if (last == 0) x_send_solve_info(send_id, c); //x_send_solve_info(send_id,c);
            }
            out_buffer_x = null;
            buffer_size = MAX_CELL_DIM * MAX_CELL_DIM * BLOCK_SIZE;
            out_buffer_x = new double[buffer_size];
            //---------------------------------------------------------------------
            //     now perform backsubstitution in reverse direction
            //---------------------------------------------------------------------
            for (stage = ncells; stage >= 1; stage--) {  //for(stage = ncells, 1, -1;
                c = slice[1,stage];
                first = 0;
                last = 0;
                if (stage == 1) first = 1;
                if (stage == ncells) {
                   last = 1; //---------------------------------------------------------------------
                             //     last cell, so perform back substitute without waiting
                             //---------------------------------------------------------------------
                   x_backsubstitute(first, last,c); //call x_backsubstitute[first, last,c];
                } else {
                   x_receive_backsub_info(out_buffer_x,recv_id, c);  //      call x_receive_backsub_info[recv_id,c];
                   send_id[0].Wait();
                   recv_id[0].Wait();
                   //      call mpi_wait[send_id,r_status,error];
                   //      call mpi_wait[recv_id,r_status,error];
                   x_unpack_backsub_info(out_buffer_x, c); 
                   x_backsubstitute(first,last,c);   //      call x_backsubstitute[first,last,c];
                }
                if (first == 0) x_send_backsub_info(send_id, c);  //call x_send_backsub_info[send_id,c];
            }
        }

        public void x_unpack_solve_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack C'[-1] and rhs'[-1] for
            //     all j and k
            //---------------------------------------------------------------------
            int j,k,m,n,ptr,istart;
            istart = 0;
            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
                for(j=0; j<=JMAX-1; j++){
                    for(m=1; m<=BLOCK_SIZE; m++){
                        for(n=0; n<BLOCK_SIZE; n++){
                            lhsc[c,k+1,j+1,istart,n+1,m] = out_buffer_x[ptr+n]; //lhsc[m,n,istart-1,j,k,c] = out_buffer[ptr+n];
                        }
                        ptr = ptr+BLOCK_SIZE;
                    }
                    for(n=0; n<BLOCK_SIZE; n++){
                        rhs[c,k+1,j+1,istart,n+1] = out_buffer_x[ptr+n]; //rhs[n,istart-1,j,k,c] = out_buffer[ptr+n];
                    }
                    ptr = ptr+BLOCK_SIZE;
                }
            }
        }

        public void x_send_solve_info(MPI.Request[] send_id, int c) {
            //---------------------------------------------------------------------
            //     pack up and send C'[iend] and rhs'[iend] for
            //     all j and k
            //---------------------------------------------------------------------
            int j,k,m,n,isize,ptr,jp,kp;
            int buffer_size;

            isize = cell_size[1,c]-1;
            jp = cell_coord[2,c] - 1;
            kp = cell_coord[3,c] - 1;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE);
            double[] in_buffer_x = new double[buffer_size];
            //---------------------------------------------------------------------
            //     pack up buffer
            //---------------------------------------------------------------------
            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
               for(j=0; j<=JMAX-1; j++){
                  for(m=1; m<=BLOCK_SIZE; m++){
                     for(n=0; n<BLOCK_SIZE; n++){
                        in_buffer_x[ptr+n] = lhsc[c,k+1,j+1,isize+1,n+1,m];  //in_buffer[ptr+n] = lhsc[m,n,isize,j,k,c];
                     }
                     ptr = ptr+BLOCK_SIZE;
                  }
                  for(n=0; n<BLOCK_SIZE; n++){
                     in_buffer_x[ptr+n] = rhs[c,k+1,j+1,isize+1,n+1];  //in_buffer[ptr+n] = rhs[n,isize,j,k,c];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
            //---------------------------------------------------------------------
            //     send buffer 
            //---------------------------------------------------------------------
            //call mpi_isend[in_buffer, buffer_size,dp_type, successor[1],WEST+jp+kp*NCELLS, comm_solve, send_id,error];
            send_id[0] = comm_solve.ImmediateSend<double>(in_buffer_x, successor[1], WEST + jp + kp * ncells);
        }

        public void x_send_backsub_info(MPI.Request[] send_id, int c) {
            //---------------------------------------------------------------------
            //     pack up and send U[istart] for all j and k
            //---------------------------------------------------------------------
            int j,k,n,ptr,istart,jp,kp;
            int buffer_size;
            //---------------------------------------------------------------------
            //     Send element 0 to previous processor
            //---------------------------------------------------------------------
            istart = 0;
            jp = cell_coord[2,c]-1;
            kp = cell_coord[3,c]-1;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE;

            double[] in_buffer_x = new double[buffer_size];

            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
               for(j=0; j<=JMAX-1; j++){
                  for(n=0; n<BLOCK_SIZE; n++){
                     in_buffer_x[ptr+n] = rhs[c,k+1,j+1,istart+1,n+1]; //in_buffer[ptr+n] = rhs[n,istart,j,k,c];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
            //call mpi_isend[in_buffer, buffer_size,dp_type, predecessor[1], EAST+jp+kp*NCELLS, comm_solve, send_id,error];
            send_id[0] = comm_solve.ImmediateSend<double>(in_buffer_x, predecessor[1], EAST+jp+kp*ncells);
        }

        public void x_unpack_backsub_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack U[isize] for all j and k
            //---------------------------------------------------------------------
            int j,k,n,ptr;
            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
               for(j=0; j<=JMAX-1; j++){
                  for(n=0; n<BLOCK_SIZE; n++){
                     backsub_info[c,k,j,n+1] = out_buffer_x[ptr+n];  //backsub_info[n,j,k,c] = out_buffer[ptr+n];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
        }

        public void x_receive_backsub_info(double[] out_buffer_x, MPI.Request[] recv_id, int c) {
            //---------------------------------------------------------------------
            //     post mpi receives
            //---------------------------------------------------------------------
            int jp,kp;
            jp = cell_coord[2,c] - 1;
            kp = cell_coord[3,c] - 1;
            //call mpi_irecv[out_buffer, buffer_size, dp_type, successor[1], EAST+jp+kp*NCELLS, comm_solve, recv_id, error];

            recv_id[0] = comm_solve.ImmediateReceive<double>(successor[1], EAST+jp+kp*ncells, out_buffer_x);
        }

        public void x_receive_solve_info(double[] out_buffer_x, MPI.Request[] recv_id, int c) {
            //---------------------------------------------------------------------
            //     post mpi receives 
            //---------------------------------------------------------------------
            int jp,kp;
            jp = cell_coord[2,c] - 1;
            kp = cell_coord[3,c] - 1;
            //call mpi_irecv[out_buffer, buffer_size, dp_type, predecessor[1], WEST+jp+kp*NCELLS,  comm_solve, recv_id, error];

            recv_id[0] = comm_solve.ImmediateReceive<double>(predecessor[1], WEST + jp + kp * ncells, out_buffer_x);
        }

        public void x_backsubstitute(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     back solve: if last cell, { generate U[isize]=rhs[isize]
            //     } else { assume U[isize] is loaded in un pack backsub_info
            //     so just use it
            //     after call u[istart] will be sent to next cell
            //---------------------------------------------------------------------
            int i, j, k;
            int m,n,isize,jsize,ksize,istart;

            istart = 0;
            isize = cell_size[1,c]-1;
            jsize = cell_size[2,c]-end[2,c]-1;
            ksize = cell_size[3,c]-end[3,c]-1;
            if (last == 0) {
               for(k=start[3,c]; k<=ksize; k++){
                  for(j=start[2,c]; j<=jsize; j++){
                      //---------------------------------------------------------------------
                      //     U[isize] uses info from previous cell if not last cell
                      //---------------------------------------------------------------------
                     for(m=1; m<=BLOCK_SIZE; m++){
                        for(n=1; n<=BLOCK_SIZE; n++){//rhs[m,isize,j,k,c] = rhs[m,isize,j,k,c] - lhsc[m,n,isize,j,k,c]*backsub_info[n,j,k,c]
                           rhs[c,k+1,j+1,isize+1,m] = rhs[c,k+1,j+1,isize+1,m] - lhsc[c,k+1,j+1,isize+1,n,m]*backsub_info[c,k,j,n];
                        }
                     }
                  }
               }
            }
            for(k=start[3,c]; k<=ksize; k++){
               for(j=start[2,c]; j<=jsize; j++){
                  for(i=isize-1; i>=istart; i--){  //for(i=isize-1,istart,-1;
                     for(m=1; m<=BLOCK_SIZE; m++){
                        for(n=1; n<=BLOCK_SIZE; n++){ //rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - lhsc[m,n,i,j,k,c]*rhs[n,i+1,j,k,c];
                           rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - lhsc[c,k+1,j+1,i+1,n,m]*rhs[c,k+1,j+1,i+2,n];
                        }
                     }
                  }
               }
            }
        }

        public void x_solve_cell(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     performs guaussian elimination on this cell.
            //     
            //     assumes that unpacking routines for non-first cells 
            //     preload C' and rhs' from previous cell.
            //     
            //     assumed send happens outside this routine, but that
            //     c'[IMAX] and rhs'[IMAX] will be sent to next cell
            //---------------------------------------------------------------------
            int i,j,k,isize,ksize,jsize,istart;

            istart = 0;
            isize = cell_size[1,c]-1;
            jsize = cell_size[2,c]-end[2,c]-1;
            ksize = cell_size[3,c]-end[3,c]-1;

            lhsabinit(ref lhsa, ref lhsb, isize);
            for(k=start[3,c]; k<=ksize; k++){
               for(j=start[2,c]; j<=jsize; j++){
                  //---------------------------------------------------------------------
                  //     This function computes the left hand side in the xi-direction
                  //---------------------------------------------------------------------
                  //---------------------------------------------------------------------
                  //     determine a [labeled f] and n jacobians for cell c
                  //---------------------------------------------------------------------
                  for(i = start[1,c]-1; i<= cell_size[1,c] - end[1,c]; i++){
                     tmp1 = rho_i[c,k+1,j+1,i+1];  //rho_i[i,j,k,c];
                     tmp2 = tmp1 * tmp1;
                     tmp3 = tmp1 * tmp2;

                     fjac[2+i, 1, 1] = 0.0d; // fjac[1, 1, i] = 0.0d;
                     fjac[2+i, 2, 1] = 1.0d; // fjac[1, 2, i] = 1.0d;
                     fjac[2+i, 3, 1] = 0.0d; // fjac[1, 3, i] = 0.0d;
                     fjac[2+i, 4, 1] = 0.0d; // fjac[1, 4, i] = 0.0d;
                     fjac[2+i, 5, 1] = 0.0d; // fjac[1, 5, i] = 0.0d;

                     fjac[2+i,1,2] = -(u[c,k+2,j+2,i+2,2] * tmp2 * u[c,k+2,j+2,i+2,2]) + c2 * qs[c,k+1,j+1,i+1];
                     fjac[2+i,2,2] = ( 2.0d - c2 ) * ( u[c,k+2,j+2,i+2,2] * tmp1 );
                     fjac[2+i,3,2] = - c2 * ( u[c,k+2,j+2,i+2,3] * tmp1 );
                     fjac[2+i,4,2] = - c2 * ( u[c,k+2,j+2,i+2,4] * tmp1 );
                     fjac[2+i,5,2] = c2;

                     fjac[2+i,1,3] = - ( u[c,k+2,j+2,i+2,2]*u[c,k+2,j+2,i+2,3] ) * tmp2;
                     fjac[2+i,2,3] = u[c,k+2,j+2,i+2,3] * tmp1;
                     fjac[2+i,3,3] = u[c,k+2,j+2,i+2,2] * tmp1;
                     fjac[2+i,4,3] = 0.0d;
                     fjac[2+i,5,3] = 0.0d;

                     fjac[2+i,1,4] = - ( u[c,k+2,j+2,i+2,2]*u[c,k+2,j+2,i+2,4] ) * tmp2;
                     fjac[2+i,2,4] = u[c,k+2,j+2,i+2,4] * tmp1;
                     fjac[2+i,3,4] = 0.0d;
                     fjac[2+i,4,4] = u[c,k+2,j+2,i+2,2] * tmp1;
                     fjac[2+i,5,4] = 0.0d;

                     fjac[2+i,1,5] = ( c2 * 2.0d * qs[c,k+1,j+1,i+1] - c1 * ( u[c,k+2,j+2,i+2,5] * tmp1 ) ) * ( u[c,k+2,j+2,i+2,2] * tmp1 );
                     fjac[2+i,2,5] = c1 *  u[c,k+2,j+2,i+2,5] * tmp1 - c2 * ( u[c,k+2,j+2,i+2,2]*u[c,k+2,j+2,i+2,2] * tmp2 + qs[c,k+1,j+1,i+1]);
                     fjac[2+i,3,5] = - c2 * ( u[c,k+2,j+2,i+2,3]*u[c,k+2,j+2,i+2,2] ) * tmp2;
                     fjac[2+i,4,5] = - c2 * ( u[c,k+2,j+2,i+2,4]*u[c,k+2,j+2,i+2,2] ) * tmp2;
                     fjac[2+i,5,5] = c1 * ( u[c,k+2,j+2,i+2,2] * tmp1 );

                     njac[2+i,1,1] = 0.0d;
                     njac[2+i,2,1] = 0.0d;
                     njac[2+i,3,1] = 0.0d;
                     njac[2+i,4,1] = 0.0d;
                     njac[2+i,5,1] = 0.0d;

                     njac[2+i,1,2] = - con43 * c3c4 * tmp2 * u[c,k+2,j+2,i+2,2];  //njac[2 + i, 1, 2] = -con43 * c3c4 * tmp2 * u[2, i, j, k, c];
                     njac[2+i,2,2] =   con43 * c3c4 * tmp1;
                     njac[2+i,3,2] =   0.0d;
                     njac[2+i,4,2] =   0.0d;
                     njac[2+i,5,2] =   0.0d;

                     njac[2+i,1,3] = - c3c4 * tmp2 * u[c,k+2,j+2,i+2,3];  //c3c4 * tmp2 * u[3,i,j,k,c];
                     njac[2+i,2,3] =   0.0d;
                     njac[2+i,3,3] =   c3c4 * tmp1;
                     njac[2+i,4,3] =   0.0d;
                     njac[2+i,5,3] =   0.0d;

                     njac[2+i,1,4] = - c3c4 * tmp2 * u[c,k+2,j+2,i+2,4];  //c3c4 * tmp2 * u[4,i,j,k,c];
                     njac[2+i,2,4] =   0.0d;
                     njac[2+i,3,4] =   0.0d;
                     njac[2+i,4,4] =   c3c4 * tmp1;
                     njac[2+i,5,4] =   0.0d;

                     njac[2+i,1,5] = - ( con43 * c3c4
                         - c1345 ) * tmp3 * (pow2(u[c,k+2,j+2,i+2,2]))           /*- c1345 ) * tmp3 * (u[2,i,j,k,c]**2) */
                         - ( c3c4 - c1345 ) * tmp3 * (pow2(u[c,k+2,j+2,i+2,3]))  /*- ( c3c4 - c1345 ) * tmp3 * (u[3,i,j,k,c]**2) */
                         - ( c3c4 - c1345 ) * tmp3 * (pow2(u[c,k+2,j+2,i+2,4]))  /*- ( c3c4 - c1345 ) * tmp3 * (u[4,i,j,k,c]**2) */
                         - c1345 * tmp2 * u[c,k+2,j+2,i+2,5];  //- c1345 * tmp2 * u[5,i,j,k,c];

                     njac[2+i,2,5] = ( con43 * c3c4 - c1345 ) * tmp2 * u[c,k+2,j+2,i+2,2];  // u[2,i,j,k,c];
                     njac[2+i,3,5] = ( c3c4 - c1345 ) * tmp2 * u[c,k+2,j+2,i+2,3];          // u[3,i,j,k,c];
                     njac[2+i,4,5] = ( c3c4 - c1345 ) * tmp2 * u[c,k+2,j+2,i+2,4];          // u[4,i,j,k,c];
                     njac[2+i,5,5] = ( c1345 ) * tmp1;
                  }
                  //---------------------------------------------------------------------
                  //     now jacobians set, so form left hand side in x direction
                  //---------------------------------------------------------------------
                  for(i = start[1,c]; i<= isize - end[1,c]; i++){

                      tmp1 = dt * tx1;
                      tmp2 = dt * tx2;

                      lhsa[1+i, 1, 1] = -tmp2 * fjac[1+i, 1, 1] - tmp1 * njac[1+i, 1, 1] - tmp1 * dx1;//2+i-1
                      lhsa[1+i, 2, 1] = -tmp2 * fjac[1+i, 2, 1] - tmp1 * njac[1+i, 2, 1];
                      lhsa[1+i, 3, 1] = -tmp2 * fjac[1+i, 3, 1] - tmp1 * njac[1+i, 3, 1];
                      lhsa[1+i, 4, 1] = -tmp2 * fjac[1+i, 4, 1] - tmp1 * njac[1+i, 4, 1];
                      lhsa[1+i, 5, 1] = -tmp2 * fjac[1+i, 5, 1] - tmp1 * njac[1+i, 5, 1];

                      lhsa[1+i, 1, 2] = -tmp2 * fjac[1+i, 1, 2] - tmp1 * njac[1+i, 1, 2];
                      lhsa[1+i, 2, 2] = -tmp2 * fjac[1+i, 2, 2] - tmp1 * njac[1+i, 2, 2] - tmp1 * dx2;
                      lhsa[1+i, 3, 2] = -tmp2 * fjac[1+i, 3, 2] - tmp1 * njac[1+i, 3, 2];
                      lhsa[1+i, 4, 2] = -tmp2 * fjac[1+i, 4, 2] - tmp1 * njac[1+i, 4, 2];
                      lhsa[1+i, 5, 2] = -tmp2 * fjac[1+i, 5, 2] - tmp1 * njac[1+i, 5, 2];

                      lhsa[1+i, 1, 3] = -tmp2 * fjac[1+i, 1, 3] - tmp1 * njac[1+i, 1, 3];
                      lhsa[1+i, 2, 3] = -tmp2 * fjac[1+i, 2, 3] - tmp1 * njac[1+i, 2, 3];
                      lhsa[1+i, 3, 3] = -tmp2 * fjac[1+i, 3, 3] - tmp1 * njac[1+i, 3, 3] - tmp1 * dx3;
                      lhsa[1+i, 4, 3] = -tmp2 * fjac[1+i, 4, 3] - tmp1 * njac[1+i, 4, 3];
                      lhsa[1+i, 5, 3] = -tmp2 * fjac[1+i, 5, 3] - tmp1 * njac[1+i, 5, 3];

                      lhsa[1+i, 1, 4] = -tmp2 * fjac[1+i, 1, 4] - tmp1 * njac[1+i, 1, 4];
                      lhsa[1+i, 2, 4] = -tmp2 * fjac[1+i, 2, 4] - tmp1 * njac[1+i, 2, 4];
                      lhsa[1+i, 3, 4] = -tmp2 * fjac[1+i, 3, 4] - tmp1 * njac[1+i, 3, 4];
                      lhsa[1+i, 4, 4] = -tmp2 * fjac[1+i, 4, 4] - tmp1 * njac[1+i, 4, 4] - tmp1 * dx4;
                      lhsa[1+i, 5, 4] = -tmp2 * fjac[1+i, 5, 4] - tmp1 * njac[1+i, 5, 4];

                      lhsa[1+i, 1, 5] = -tmp2 * fjac[1+i, 1, 5] - tmp1 * njac[1+i, 1, 5];
                      lhsa[1+i, 2, 5] = -tmp2 * fjac[1+i, 2, 5] - tmp1 * njac[1+i, 2, 5];
                      lhsa[1+i, 3, 5] = -tmp2 * fjac[1+i, 3, 5] - tmp1 * njac[1+i, 3, 5];
                      lhsa[1+i, 4, 5] = -tmp2 * fjac[1+i, 4, 5] - tmp1 * njac[1+i, 4, 5];
                      lhsa[1+i, 5, 5] = -tmp2 * fjac[1+i, 5, 5] - tmp1 * njac[1+i, 5, 5] - tmp1 * dx5;

                      lhsb[1 + i, 1, 1] = 1.0d + tmp1 * 2.0d * njac[2 + i, 1, 1] + tmp1 * 2.0d * dx1;
                      lhsb[1 + i, 2, 1] = tmp1 * 2.0d * njac[2 + i, 2, 1];
                      lhsb[1 + i, 3, 1] = tmp1 * 2.0d * njac[2 + i, 3, 1];
                      lhsb[1 + i, 4, 1] = tmp1 * 2.0d * njac[2 + i, 4, 1];
                      lhsb[1 + i, 5, 1] = tmp1 * 2.0d * njac[2 + i, 5, 1];

                      lhsb[1 + i, 1, 2] = tmp1 * 2.0d * njac[2 + i, 1, 2];
                      lhsb[1 + i, 2, 2] = 1.0d + tmp1 * 2.0d * njac[2 + i, 2, 2] + tmp1 * 2.0d * dx2;
                      lhsb[1 + i, 3, 2] = tmp1 * 2.0d * njac[2 + i, 3, 2];
                      lhsb[1 + i, 4, 2] = tmp1 * 2.0d * njac[2 + i, 4, 2];
                      lhsb[1 + i, 5, 2] = tmp1 * 2.0d * njac[2 + i, 5, 2];

                      lhsb[1 + i, 1, 3] = tmp1 * 2.0d * njac[2 + i, 1, 3];
                      lhsb[1 + i, 2, 3] = tmp1 * 2.0d * njac[2 + i, 2, 3];
                      lhsb[1 + i, 3, 3] = 1.0d + tmp1 * 2.0d * njac[2 + i, 3, 3] + tmp1 * 2.0d * dx3;
                      lhsb[1 + i, 4, 3] = tmp1 * 2.0d * njac[2 + i, 4, 3];
                      lhsb[1 + i, 5, 3] = tmp1 * 2.0d * njac[2 + i, 5, 3];

                      lhsb[1 + i, 1, 4] = tmp1 * 2.0d * njac[2 + i, 1, 4];
                      lhsb[1 + i, 2, 4] = tmp1 * 2.0d * njac[2 + i, 2, 4];
                      lhsb[1 + i, 3, 4] = tmp1 * 2.0d * njac[2 + i, 3, 4];
                      lhsb[1 + i, 4, 4] = 1.0d + tmp1 * 2.0d * njac[2 + i, 4, 4] + tmp1 * 2.0d * dx4;
                      lhsb[1 + i, 5, 4] = tmp1 * 2.0d * njac[2 + i, 5, 4];

                      lhsb[1 + i, 1, 5] = tmp1 * 2.0d * njac[2 + i, 1, 5];
                      lhsb[1 + i, 2, 5] = tmp1 * 2.0d * njac[2 + i, 2, 5];
                      lhsb[1 + i, 3, 5] = tmp1 * 2.0d * njac[2 + i, 3, 5];
                      lhsb[1 + i, 4, 5] = tmp1 * 2.0d * njac[2 + i, 4, 5];
                      lhsb[1 + i, 5, 5] = 1.0d + tmp1 * 2.0d * njac[2 + i, 5, 5] + tmp1 * 2.0d * dx5;

                      lhsc[c, k+1, j+1, i+1, 1, 1] = tmp2 * fjac[3+i, 1, 1] - tmp1 * njac[3+i, 1, 1] - tmp1 * dx1; //2+i+1
                      lhsc[c, k+1, j+1, i+1, 2, 1] = tmp2 * fjac[3+i, 2, 1] - tmp1 * njac[3+i, 2, 1];
                      lhsc[c, k+1, j+1, i+1, 3, 1] = tmp2 * fjac[3+i, 3, 1] - tmp1 * njac[3+i, 3, 1];
                      lhsc[c, k+1, j+1, i+1, 4, 1] = tmp2 * fjac[3+i, 4, 1] - tmp1 * njac[3+i, 4, 1];
                      lhsc[c, k+1, j+1, i+1, 5, 1] = tmp2 * fjac[3+i, 5, 1] - tmp1 * njac[3+i, 5, 1];

                      lhsc[c, k+1, j+1, i+1, 1, 2] = tmp2 * fjac[3+i, 1, 2] - tmp1 * njac[3+i, 1, 2];
                      lhsc[c, k+1, j+1, i+1, 2, 2] = tmp2 * fjac[3+i, 2, 2] - tmp1 * njac[3+i, 2, 2] - tmp1 * dx2;
                      lhsc[c, k+1, j+1, i+1, 3, 2] = tmp2 * fjac[3+i, 3, 2] - tmp1 * njac[3+i, 3, 2];
                      lhsc[c, k+1, j+1, i+1, 4, 2] = tmp2 * fjac[3+i, 4, 2] - tmp1 * njac[3+i, 4, 2];
                      lhsc[c, k+1, j+1, i+1, 5, 2] = tmp2 * fjac[3+i, 5, 2] - tmp1 * njac[3+i, 5, 2];

                      lhsc[c, k+1, j+1, i+1, 1, 3] = tmp2 * fjac[3+i, 1, 3] - tmp1 * njac[3+i, 1, 3];
                      lhsc[c, k+1, j+1, i+1, 2, 3] = tmp2 * fjac[3+i, 2, 3] - tmp1 * njac[3+i, 2, 3];
                      lhsc[c, k+1, j+1, i+1, 3, 3] = tmp2 * fjac[3+i, 3, 3] - tmp1 * njac[3+i, 3, 3] - tmp1 * dx3;
                      lhsc[c, k+1, j+1, i+1, 4, 3] = tmp2 * fjac[3+i, 4, 3] - tmp1 * njac[3+i, 4, 3];
                      lhsc[c, k+1, j+1, i+1, 5, 3] = tmp2 * fjac[3+i, 5, 3] - tmp1 * njac[3+i, 5, 3];

                      lhsc[c, k+1, j+1, i+1, 1, 4] = tmp2 * fjac[3+i, 1, 4] - tmp1 * njac[3+i, 1, 4];
                      lhsc[c, k+1, j+1, i+1, 2, 4] = tmp2 * fjac[3+i, 2, 4] - tmp1 * njac[3+i, 2, 4];
                      lhsc[c, k+1, j+1, i+1, 3, 4] = tmp2 * fjac[3+i, 3, 4] - tmp1 * njac[3+i, 3, 4];
                      lhsc[c, k+1, j+1, i+1, 4, 4] = tmp2 * fjac[3+i, 4, 4] - tmp1 * njac[3+i, 4, 4] - tmp1 * dx4;
                      lhsc[c, k+1, j+1, i+1, 5, 4] = tmp2 * fjac[3+i, 5, 4] - tmp1 * njac[3+i, 5, 4];

                      lhsc[c, k+1, j+1, i+1, 1, 5] = tmp2 * fjac[3+i, 1, 5] - tmp1 * njac[3+i, 1, 5];
                      lhsc[c, k+1, j+1, i+1, 2, 5] = tmp2 * fjac[3+i, 2, 5] - tmp1 * njac[3+i, 2, 5];
                      lhsc[c, k+1, j+1, i+1, 3, 5] = tmp2 * fjac[3+i, 3, 5] - tmp1 * njac[3+i, 3, 5];
                      lhsc[c, k+1, j+1, i+1, 4, 5] = tmp2 * fjac[3+i, 4, 5] - tmp1 * njac[3+i, 4, 5];
                      lhsc[c, k+1, j+1, i+1, 5, 5] = tmp2 * fjac[3+i, 5, 5] - tmp1 * njac[3+i, 5, 5] - tmp1 * dx5;
                  }
                  //---------------------------------------------------------------------
                  //     outer most for(loops - sweeping in i direction
                  //---------------------------------------------------------------------
                  if (first == 1) { 
                      //---------------------------------------------------------------------
                      //     multiply c[istart,j,k] by b_inverse and copy back to c
                      //     multiply rhs[istart] by b_inverse[istart] and copy to rhs
                      //---------------------------------------------------------------------
                               /* call binvcrhs(lhsb[1,1,istart], lhsc[1,1,istart,j,k,c], rhs[1,istart,j,k,c]); */
                      binvcrhs(ref lhsb, ref lhsc, ref rhs,   istart+1,      c,k+1,j+1,istart+1,      c,k+1,j+1,istart+1);
                      //binvcrhs(lhsb[istart+1, 1, 1], lhsc[c, k+1, j+1, istart+1, 1, 1], rhs[c, k+1, j+1, istart+1, 1]);  
                  }
                  //---------------------------------------------------------------------
                  //     begin inner most for(loop
                  //     for(all the elements of the cell unless last 
                  //---------------------------------------------------------------------
                  for(i=istart+first; i<=isize-last; i++){
                        //---------------------------------------------------------------------
                        //     rhs[i] = rhs[i] - A*rhs[i-1]
                        //---------------------------------------------------------------------
                            /*Fortran:  matvec_sub(lhsa[1,1,i], rhs[1,i-1,j,k,c],rhs[1,i,j,k,c]);*/
                            /*Padrao c# matvec_sub(lhsa[i+1,1,1], rhs[c,k+1,j+1,i-1+1,1],rhs[c,k+1,j+1,i+1,1]);*/
                        matvec_sub(ref lhsa,ref rhs,ref rhs,    (i+1),       c,(k+1),(j+1),(i-1+1),           c,(k+1),(j+1),(i+1)); 
                        //---------------------------------------------------------------------
                        //     B[i] = B[i] - C[i-1]*A[i]
                        //---------------------------------------------------------------------
                        /*Fortran:   matmul_sub(lhsa[1,1,i], lhsc[1,1,i-1,j,k,c], lhsb[1,1,i]);*/
                        /*Padrao c#: matmul_sub(lhsa[i+1,1,1], lhsc[c,k+1,j+1,i-1+1,1,1], lhsb[i+1,1,1]);*/
                        matmul_sub(ref lhsa,ref lhsc,ref lhsb,    (i+1),    c,(k+1),(j+1),(i-1+1),    (i+1));
                        //---------------------------------------------------------------------
                        //     multiply c[i,j,k] by b_inverse and copy back to c
                        //     multiply rhs[1,j,k] by b_inverse[1,j,k] and copy to rhs
                        //---------------------------------------------------------------------
                          /*Fortran binvcrhs( lhsb[1,1,i], lhsc[1,1,i,j,k,c], rhs[1,i,j,k,c] );*/
                        binvcrhs(ref lhsb,ref lhsc, ref rhs,   i+1,    c, k+1, j+1, i+1,    c,k+1,j+1,i+1);
                  }
                  //---------------------------------------------------------------------
                  //     Now finish up special cases for last cell
                  //---------------------------------------------------------------------
                  if (last == 1) {
                      //---------------------------------------------------------------------
                      //     rhs[isize] = rhs[isize] - A*rhs[isize-1]
                      //---------------------------------------------------------------------
                      /*Fortran: matvec_sub[lhsa[1,1,isize], rhs[1,isize-1,j,k,c],rhs[1,isize,j,k,c]];*/
                      /*C#:      matvec_sub[lhsa[isize+1,1,1], rhs[c,k+1,j+1,isize-1+1,1],rhs[c,k+1,j+1,isize+1,1]];*/
                      matvec_sub(ref lhsa, ref rhs, ref rhs, (isize + 1), c, (k + 1), (j + 1), (isize), c, (k + 1), (j + 1), (isize + 1));
                      //---------------------------------------------------------------------
                      //     B[isize] = B[isize] - C[isize-1]*A[isize]
                      //---------------------------------------------------------------------
                      /*Fortran: matmul_sub[lhsa[1,1,isize], lhsc[1,1,isize-1,j,k,c], lhsb[1,1,isize]];*/
                      /*C#:      matmul_sub(lhsa[isize+1,1,1], lhsc[c,k+1,j+1,isize-1+1,1,1], lhsb[isize+1,1,1]);*/
                      matmul_sub(ref lhsa, ref lhsc, ref lhsb, (isize + 1), c, (k + 1), (j + 1), isize, (isize + 1));
                      //---------------------------------------------------------------------
                      //     multiply rhs[] by b_inverse[] and copy to rhs
                      //---------------------------------------------------------------------
                      /*Fortran: binvrhs[ lhsb[1,1,isize], rhs[1,isize,j,k,c] ];*/
                      /*C#:      binvrhs( lhsb[isize+1,1,1], rhs[c,k+1,j+1,isize+1,1] ); */
                      binvrhs(ref lhsb, ref rhs, (isize + 1), c, (k + 1), (j + 1), (isize + 1));
                  }
               }
            }
        }
                   // solve_subs.f
        public void binvcrhs(ref double[,,] lhs, ref double[,,,,,] c, ref double[,,,,] r, int l1, int c1,int c2,int c3,int c4, int r1,int r2,int r3,int r4) {
            double pivot, coeff; //dimension lhs[5,5]; //double c[5,5], r[5];
            pivot = 1.00d/lhs[l1,1,1];
            lhs[l1,2,1] = lhs[l1,2,1]*pivot;
            lhs[l1,3,1] = lhs[l1,3,1]*pivot;
            lhs[l1,4,1] = lhs[l1,4,1]*pivot;
            lhs[l1,5,1] = lhs[l1,5,1]*pivot;

            c[c1,c2,c3,c4,1,1] = c[c1,c2,c3,c4,1,1] * pivot;
            c[c1,c2,c3,c4,2,1] = c[c1,c2,c3,c4,2,1] * pivot;
            c[c1,c2,c3,c4,3,1] = c[c1,c2,c3,c4,3,1] * pivot;
            c[c1,c2,c3,c4,4,1] = c[c1,c2,c3,c4,4,1] * pivot;
            c[c1,c2,c3,c4,5,1] = c[c1,c2,c3,c4,5,1] * pivot;

            r[r1,r2,r3,r4,1]   = r[r1,r2,r3,r4,1]  *pivot;
            coeff = lhs[l1,1,2];
            lhs[l1,2, 2] = lhs[l1,2, 2] - coeff * lhs[l1,2, 1];
            lhs[l1,3, 2] = lhs[l1,3, 2] - coeff * lhs[l1,3, 1];
            lhs[l1,4, 2] = lhs[l1,4, 2] - coeff * lhs[l1,4, 1];
            lhs[l1,5, 2] = lhs[l1,5, 2] - coeff * lhs[l1,5, 1];

            c[c1,c2,c3,c4,1,2] = c[c1,c2,c3,c4,1,2] - coeff * c[c1,c2,c3,c4,1,1];
            c[c1,c2,c3,c4,2,2] = c[c1,c2,c3,c4,2,2] - coeff * c[c1,c2,c3,c4,2,1];
            c[c1,c2,c3,c4,3,2] = c[c1,c2,c3,c4,3,2] - coeff * c[c1,c2,c3,c4,3,1];
            c[c1,c2,c3,c4,4,2] = c[c1,c2,c3,c4,4,2] - coeff * c[c1,c2,c3,c4,4,1];
            c[c1,c2,c3,c4,5,2] = c[c1,c2,c3,c4,5,2] - coeff * c[c1,c2,c3,c4,5,1];

            r[r1,r2,r3,r4,2]   = r[r1,r2,r3,r4,2]   - coeff*r[r1,r2,r3,r4,1];
            coeff = lhs[l1,1, 3];
            lhs[l1,2, 3] = lhs[l1,2, 3] - coeff * lhs[l1,2, 1];
            lhs[l1,3, 3] = lhs[l1,3, 3] - coeff * lhs[l1,3, 1];
            lhs[l1,4, 3] = lhs[l1,4, 3] - coeff * lhs[l1,4, 1];
            lhs[l1,5, 3] = lhs[l1,5, 3] - coeff * lhs[l1,5, 1];

            c[c1,c2,c3,c4,1,3] = c[c1,c2,c3,c4,1, 3] - coeff * c[c1,c2,c3,c4,1, 1];
            c[c1,c2,c3,c4,2,3] = c[c1,c2,c3,c4,2, 3] - coeff * c[c1,c2,c3,c4,2, 1];
            c[c1,c2,c3,c4,3,3] = c[c1,c2,c3,c4,3, 3] - coeff * c[c1,c2,c3,c4,3, 1];
            c[c1,c2,c3,c4,4,3] = c[c1,c2,c3,c4,4, 3] - coeff * c[c1,c2,c3,c4,4, 1];
            c[c1,c2,c3,c4,5,3] = c[c1,c2,c3,c4,5, 3] - coeff * c[c1,c2,c3,c4,5, 1];

            r[r1,r2,r3,r4,3]   = r[r1,r2,r3,r4,3]   - coeff*r[r1,r2,r3,r4,1];
            coeff = lhs[l1,1, 4];
            lhs[l1,2, 4] = lhs[l1,2, 4] - coeff * lhs[l1,2, 1];
            lhs[l1,3, 4] = lhs[l1,3, 4] - coeff * lhs[l1,3, 1];
            lhs[l1,4, 4] = lhs[l1,4, 4] - coeff * lhs[l1,4, 1];
            lhs[l1,5, 4] = lhs[l1,5, 4] - coeff * lhs[l1,5, 1];

            c[c1,c2,c3,c4,1, 4] = c[c1,c2,c3,c4,1, 4] - coeff * c[c1,c2,c3,c4,1, 1];
            c[c1,c2,c3,c4,2, 4] = c[c1,c2,c3,c4,2, 4] - coeff * c[c1,c2,c3,c4,2, 1];
            c[c1,c2,c3,c4,3, 4] = c[c1,c2,c3,c4,3, 4] - coeff * c[c1,c2,c3,c4,3, 1];
            c[c1,c2,c3,c4,4, 4] = c[c1,c2,c3,c4,4, 4] - coeff * c[c1,c2,c3,c4,4, 1];
            c[c1,c2,c3,c4,5, 4] = c[c1,c2,c3,c4,5, 4] - coeff * c[c1,c2,c3,c4,5, 1];

            r[r1,r2,r3,r4,4]   = r[r1,r2,r3,r4,4]   - coeff*r[r1,r2,r3,r4,1];
            coeff = lhs[l1,1, 5];
            lhs[l1,2, 5] = lhs[l1,2, 5] - coeff * lhs[l1,2, 1];
            lhs[l1,3, 5] = lhs[l1,3, 5] - coeff * lhs[l1,3, 1];
            lhs[l1,4, 5] = lhs[l1,4, 5] - coeff * lhs[l1,4, 1];
            lhs[l1,5, 5] = lhs[l1,5, 5] - coeff * lhs[l1,5, 1];

            c[c1,c2,c3,c4,1, 5] = c[c1,c2,c3,c4,1, 5] - coeff * c[c1,c2,c3,c4,1, 1];
            c[c1,c2,c3,c4,2, 5] = c[c1,c2,c3,c4,2, 5] - coeff * c[c1,c2,c3,c4,2, 1];
            c[c1,c2,c3,c4,3, 5] = c[c1,c2,c3,c4,3, 5] - coeff * c[c1,c2,c3,c4,3, 1];
            c[c1,c2,c3,c4,4, 5] = c[c1,c2,c3,c4,4, 5] - coeff * c[c1,c2,c3,c4,4, 1];
            c[c1,c2,c3,c4,5, 5] = c[c1,c2,c3,c4,5, 5] - coeff * c[c1,c2,c3,c4,5, 1];

            r[r1,r2,r3,r4,5]   = r[r1,r2,r3,r4,5]   - coeff*r[r1,r2,r3,r4,1];
            pivot = 1.00d/lhs[l1,2,2];
            lhs[l1,3,2] = lhs[l1,3,2]*pivot;
            lhs[l1,4,2] = lhs[l1,4,2]*pivot;
            lhs[l1,5,2] = lhs[l1,5,2]*pivot;

            c[c1,c2,c3,c4,1, 2] = c[c1,c2,c3,c4,1, 2] * pivot;
            c[c1,c2,c3,c4,2, 2] = c[c1,c2,c3,c4,2, 2] * pivot;
            c[c1,c2,c3,c4,3, 2] = c[c1,c2,c3,c4,3, 2] * pivot;
            c[c1,c2,c3,c4,4, 2] = c[c1,c2,c3,c4,4, 2] * pivot;
            c[c1,c2,c3,c4,5, 2] = c[c1,c2,c3,c4,5, 2] * pivot;

            r[r1,r2,r3,r4,2]   = r[r1,r2,r3,r4,2]  *pivot;
            coeff = lhs[l1,2,1];
            lhs[l1,3,1]= lhs[l1,3,1] - coeff*lhs[l1,3,2];
            lhs[l1,4,1]= lhs[l1,4,1] - coeff*lhs[l1,4,2];
            lhs[l1,5,1]= lhs[l1,5,1] - coeff*lhs[l1,5,2];

            c[c1,c2,c3,c4,1, 1] = c[c1,c2,c3,c4,1, 1] - coeff * c[c1,c2,c3,c4,1, 2];
            c[c1,c2,c3,c4,2, 1] = c[c1,c2,c3,c4,2, 1] - coeff * c[c1,c2,c3,c4,2, 2];
            c[c1,c2,c3,c4,3, 1] = c[c1,c2,c3,c4,3, 1] - coeff * c[c1,c2,c3,c4,3, 2];
            c[c1,c2,c3,c4,4, 1] = c[c1,c2,c3,c4,4, 1] - coeff * c[c1,c2,c3,c4,4, 2];
            c[c1,c2,c3,c4,5, 1] = c[c1,c2,c3,c4,5, 1] - coeff * c[c1,c2,c3,c4,5, 2];

            r[r1,r2,r3,r4,1]   = r[r1,r2,r3,r4,1]   - coeff*r[r1,r2,r3,r4,2];
            coeff = lhs[l1,2,3];
            lhs[l1,3,3]= lhs[l1,3,3] - coeff*lhs[l1,3,2];
            lhs[l1,4,3]= lhs[l1,4,3] - coeff*lhs[l1,4,2];
            lhs[l1,5,3]= lhs[l1,5,3] - coeff*lhs[l1,5,2];

            c[c1,c2,c3,c4,1, 3] = c[c1,c2,c3,c4,1, 3] - coeff * c[c1,c2,c3,c4,1, 2];
            c[c1,c2,c3,c4,2, 3] = c[c1,c2,c3,c4,2, 3] - coeff * c[c1,c2,c3,c4,2, 2];
            c[c1,c2,c3,c4,3, 3] = c[c1,c2,c3,c4,3, 3] - coeff * c[c1,c2,c3,c4,3, 2];
            c[c1,c2,c3,c4,4, 3] = c[c1,c2,c3,c4,4, 3] - coeff * c[c1,c2,c3,c4,4, 2];
            c[c1,c2,c3,c4,5, 3] = c[c1,c2,c3,c4,5, 3] - coeff * c[c1,c2,c3,c4,5, 2];

            r[r1,r2,r3,r4,3]   = r[r1,r2,r3,r4,3]   - coeff*r[r1,r2,r3,r4,2];
            coeff = lhs[l1,2,4];
            lhs[l1,3,4]= lhs[l1,3,4] - coeff*lhs[l1,3,2];
            lhs[l1,4,4]= lhs[l1,4,4] - coeff*lhs[l1,4,2];
            lhs[l1,5,4]= lhs[l1,5,4] - coeff*lhs[l1,5,2];

            c[c1,c2,c3,c4,1, 4] = c[c1,c2,c3,c4,1, 4] - coeff * c[c1,c2,c3,c4,1, 2];
            c[c1,c2,c3,c4,2, 4] = c[c1,c2,c3,c4,2, 4] - coeff * c[c1,c2,c3,c4,2, 2];
            c[c1,c2,c3,c4,3, 4] = c[c1,c2,c3,c4,3, 4] - coeff * c[c1,c2,c3,c4,3, 2];
            c[c1,c2,c3,c4,4, 4] = c[c1,c2,c3,c4,4, 4] - coeff * c[c1,c2,c3,c4,4, 2];
            c[c1,c2,c3,c4,5, 4] = c[c1,c2,c3,c4,5, 4] - coeff * c[c1,c2,c3,c4,5, 2];

            r[r1,r2,r3,r4,4]   = r[r1,r2,r3,r4,4]   - coeff*r[r1,r2,r3,r4,2];
            coeff = lhs[l1,2,5];
            lhs[l1,3,5]= lhs[l1,3,5] - coeff*lhs[l1,3,2];
            lhs[l1,4,5]= lhs[l1,4,5] - coeff*lhs[l1,4,2];
            lhs[l1,5,5]= lhs[l1,5,5] - coeff*lhs[l1,5,2];

            c[c1,c2,c3,c4,1, 5] = c[c1,c2,c3,c4,1, 5] - coeff * c[c1,c2,c3,c4,1, 2];
            c[c1,c2,c3,c4,2, 5] = c[c1,c2,c3,c4,2, 5] - coeff * c[c1,c2,c3,c4,2, 2];
            c[c1,c2,c3,c4,3, 5] = c[c1,c2,c3,c4,3, 5] - coeff * c[c1,c2,c3,c4,3, 2];
            c[c1,c2,c3,c4,4, 5] = c[c1,c2,c3,c4,4, 5] - coeff * c[c1,c2,c3,c4,4, 2];
            c[c1,c2,c3,c4,5, 5] = c[c1,c2,c3,c4,5, 5] - coeff * c[c1,c2,c3,c4,5, 2];

            r[r1,r2,r3,r4,5]   = r[r1,r2,r3,r4,5]   - coeff*r[r1,r2,r3,r4,2];
            pivot = 1.00d/lhs[l1,3,3];
            lhs[l1,4,3] = lhs[l1,4,3]*pivot;
            lhs[l1,5,3] = lhs[l1,5,3]*pivot;

            c[c1,c2,c3,c4,1, 3] = c[c1,c2,c3,c4,1, 3] * pivot;
            c[c1,c2,c3,c4,2, 3] = c[c1,c2,c3,c4,2, 3] * pivot;
            c[c1,c2,c3,c4,3, 3] = c[c1,c2,c3,c4,3, 3] * pivot;
            c[c1,c2,c3,c4,4, 3] = c[c1,c2,c3,c4,4, 3] * pivot;
            c[c1,c2,c3,c4,5, 3] = c[c1,c2,c3,c4,5, 3] * pivot;

            r[r1,r2,r3,r4,3]   = r[r1,r2,r3,r4,3]  *pivot;
            coeff = lhs[l1,3,1];
            lhs[l1,4,1]= lhs[l1,4,1] - coeff*lhs[l1,4,3];
            lhs[l1,5,1]= lhs[l1,5,1] - coeff*lhs[l1,5,3];

            c[c1,c2,c3,c4,1, 1] = c[c1,c2,c3,c4,1, 1] - coeff * c[c1,c2,c3,c4,1, 3];
            c[c1,c2,c3,c4,2, 1] = c[c1,c2,c3,c4,2, 1] - coeff * c[c1,c2,c3,c4,2, 3];
            c[c1,c2,c3,c4,3, 1] = c[c1,c2,c3,c4,3, 1] - coeff * c[c1,c2,c3,c4,3, 3];
            c[c1,c2,c3,c4,4, 1] = c[c1,c2,c3,c4,4, 1] - coeff * c[c1,c2,c3,c4,4, 3];
            c[c1,c2,c3,c4,5, 1] = c[c1,c2,c3,c4,5, 1] - coeff * c[c1,c2,c3,c4,5, 3];

            r[r1,r2,r3,r4,1]   = r[r1,r2,r3,r4,1]   - coeff*r[r1,r2,r3,r4,3];
            coeff = lhs[l1,3,2];
            lhs[l1,4,2]= lhs[l1,4,2] - coeff*lhs[l1,4,3];
            lhs[l1,5,2]= lhs[l1,5,2] - coeff*lhs[l1,5,3];

            c[c1,c2,c3,c4,1, 2] = c[c1,c2,c3,c4,1, 2] - coeff * c[c1,c2,c3,c4,1, 3];
            c[c1,c2,c3,c4,2, 2] = c[c1,c2,c3,c4,2, 2] - coeff * c[c1,c2,c3,c4,2, 3];
            c[c1,c2,c3,c4,3, 2] = c[c1,c2,c3,c4,3, 2] - coeff * c[c1,c2,c3,c4,3, 3];
            c[c1,c2,c3,c4,4, 2] = c[c1,c2,c3,c4,4, 2] - coeff * c[c1,c2,c3,c4,4, 3];
            c[c1,c2,c3,c4,5, 2] = c[c1,c2,c3,c4,5, 2] - coeff * c[c1,c2,c3,c4,5, 3];

            r[r1,r2,r3,r4,2]   = r[r1,r2,r3,r4,2]   - coeff*r[r1,r2,r3,r4,3];
            coeff = lhs[l1,3,4];
            lhs[l1,4,4]= lhs[l1,4,4] - coeff*lhs[l1,4,3];
            lhs[l1,5,4]= lhs[l1,5,4] - coeff*lhs[l1,5,3];

            c[c1,c2,c3,c4,1, 4] = c[c1,c2,c3,c4,1, 4] - coeff * c[c1,c2,c3,c4,1, 3];
            c[c1,c2,c3,c4,2, 4] = c[c1,c2,c3,c4,2, 4] - coeff * c[c1,c2,c3,c4,2, 3];
            c[c1,c2,c3,c4,3, 4] = c[c1,c2,c3,c4,3, 4] - coeff * c[c1,c2,c3,c4,3, 3];
            c[c1,c2,c3,c4,4, 4] = c[c1,c2,c3,c4,4, 4] - coeff * c[c1,c2,c3,c4,4, 3];
            c[c1,c2,c3,c4,5, 4] = c[c1,c2,c3,c4,5, 4] - coeff * c[c1,c2,c3,c4,5, 3];

            r[r1,r2,r3,r4,4]   = r[r1,r2,r3,r4,4]   - coeff*r[r1,r2,r3,r4,3];
            coeff = lhs[l1,3,5];
            lhs[l1,4,5]= lhs[l1,4,5] - coeff*lhs[l1,4,3];
            lhs[l1,5,5]= lhs[l1,5,5] - coeff*lhs[l1,5,3];

            c[c1,c2,c3,c4,1, 5] = c[c1,c2,c3,c4,1, 5] - coeff * c[c1,c2,c3,c4,1, 3];
            c[c1,c2,c3,c4,2, 5] = c[c1,c2,c3,c4,2, 5] - coeff * c[c1,c2,c3,c4,2, 3];
            c[c1,c2,c3,c4,3, 5] = c[c1,c2,c3,c4,3, 5] - coeff * c[c1,c2,c3,c4,3, 3];
            c[c1,c2,c3,c4,4, 5] = c[c1,c2,c3,c4,4, 5] - coeff * c[c1,c2,c3,c4,4, 3];
            c[c1,c2,c3,c4,5, 5] = c[c1,c2,c3,c4,5, 5] - coeff * c[c1,c2,c3,c4,5, 3];

            r[r1,r2,r3,r4,5]   = r[r1,r2,r3,r4,5]   - coeff*r[r1,r2,r3,r4,3];
            pivot = 1.00d/lhs[l1,4,4];
            lhs[l1,5,4] = lhs[l1,5,4]*pivot;

            c[c1,c2,c3,c4,1, 4] = c[c1,c2,c3,c4,1, 4] * pivot;
            c[c1,c2,c3,c4,2, 4] = c[c1,c2,c3,c4,2, 4] * pivot;
            c[c1,c2,c3,c4,3, 4] = c[c1,c2,c3,c4,3, 4] * pivot;
            c[c1,c2,c3,c4,4, 4] = c[c1,c2,c3,c4,4, 4] * pivot;
            c[c1,c2,c3,c4,5, 4] = c[c1,c2,c3,c4,5, 4] * pivot;

            r[r1,r2,r3,r4,4]   = r[r1,r2,r3,r4,4]  *pivot;
            coeff = lhs[l1,4,1];
            lhs[l1,5,1]= lhs[l1,5,1] - coeff*lhs[l1,5,4];

            c[c1,c2,c3,c4,1, 1] = c[c1,c2,c3,c4,1, 1] - coeff * c[c1,c2,c3,c4,1, 4];
            c[c1,c2,c3,c4,2, 1] = c[c1,c2,c3,c4,2, 1] - coeff * c[c1,c2,c3,c4,2, 4];
            c[c1,c2,c3,c4,3, 1] = c[c1,c2,c3,c4,3, 1] - coeff * c[c1,c2,c3,c4,3, 4];
            c[c1,c2,c3,c4,4, 1] = c[c1,c2,c3,c4,4, 1] - coeff * c[c1,c2,c3,c4,4, 4];
            c[c1,c2,c3,c4,5, 1] = c[c1,c2,c3,c4,5, 1] - coeff * c[c1,c2,c3,c4,5, 4];

            r[r1,r2,r3,r4,1]   = r[r1,r2,r3,r4,1]   - coeff*r[r1,r2,r3,r4,4];
            coeff = lhs[l1,4,2];
            lhs[l1,5,2]= lhs[l1,5,2] - coeff*lhs[l1,5,4];

            c[c1,c2,c3,c4,1, 2] = c[c1,c2,c3,c4,1, 2] - coeff * c[c1,c2,c3,c4,1, 4];
            c[c1,c2,c3,c4,2, 2] = c[c1,c2,c3,c4,2, 2] - coeff * c[c1,c2,c3,c4,2, 4];
            c[c1,c2,c3,c4,3, 2] = c[c1,c2,c3,c4,3, 2] - coeff * c[c1,c2,c3,c4,3, 4];
            c[c1,c2,c3,c4,4, 2] = c[c1,c2,c3,c4,4, 2] - coeff * c[c1,c2,c3,c4,4, 4];
            c[c1,c2,c3,c4,5, 2] = c[c1,c2,c3,c4,5, 2] - coeff * c[c1,c2,c3,c4,5, 4];

            r[r1,r2,r3,r4,2]   = r[r1,r2,r3,r4,2]   - coeff*r[r1,r2,r3,r4,4];
            coeff = lhs[l1,4,3];
            lhs[l1,5,3]= lhs[l1,5,3] - coeff*lhs[l1,5,4];

            c[c1,c2,c3,c4,1, 3] = c[c1,c2,c3,c4,1, 3] - coeff * c[c1,c2,c3,c4,1, 4];
            c[c1,c2,c3,c4,2, 3] = c[c1,c2,c3,c4,2, 3] - coeff * c[c1,c2,c3,c4,2, 4];
            c[c1,c2,c3,c4,3, 3] = c[c1,c2,c3,c4,3, 3] - coeff * c[c1,c2,c3,c4,3, 4];
            c[c1,c2,c3,c4,4, 3] = c[c1,c2,c3,c4,4, 3] - coeff * c[c1,c2,c3,c4,4, 4];
            c[c1,c2,c3,c4,5, 3] = c[c1,c2,c3,c4,5, 3] - coeff * c[c1,c2,c3,c4,5, 4];

            r[r1,r2,r3,r4,3]   = r[r1,r2,r3,r4,3]   - coeff*r[r1,r2,r3,r4,4];
            coeff = lhs[l1,4,5];
            lhs[l1,5,5]= lhs[l1,5,5] - coeff*lhs[l1,5,4];

            c[c1,c2,c3,c4,1, 5] = c[c1,c2,c3,c4,1, 5] - coeff * c[c1,c2,c3,c4,1, 4];
            c[c1,c2,c3,c4,2, 5] = c[c1,c2,c3,c4,2, 5] - coeff * c[c1,c2,c3,c4,2, 4];
            c[c1,c2,c3,c4,3, 5] = c[c1,c2,c3,c4,3, 5] - coeff * c[c1,c2,c3,c4,3, 4];
            c[c1,c2,c3,c4,4, 5] = c[c1,c2,c3,c4,4, 5] - coeff * c[c1,c2,c3,c4,4, 4];
            c[c1,c2,c3,c4,5, 5] = c[c1,c2,c3,c4,5, 5] - coeff * c[c1,c2,c3,c4,5, 4];

            r[r1,r2,r3,r4,5]   = r[r1,r2,r3,r4,5]   - coeff*r[r1,r2,r3,r4,4];
            pivot = 1.00d/lhs[l1,5,5];
            c[c1,c2,c3,c4,1, 5] = c[c1,c2,c3,c4,1, 5] * pivot;
            c[c1,c2,c3,c4,2, 5] = c[c1,c2,c3,c4,2, 5] * pivot;
            c[c1,c2,c3,c4,3, 5] = c[c1,c2,c3,c4,3, 5] * pivot;
            c[c1,c2,c3,c4,4, 5] = c[c1,c2,c3,c4,4, 5] * pivot;
            c[c1,c2,c3,c4,5, 5] = c[c1,c2,c3,c4,5, 5] * pivot;

            r[r1,r2,r3,r4,5]   = r[r1,r2,r3,r4,5]  *pivot;
            coeff = lhs[l1,5,1];
            c[c1,c2,c3,c4,1, 1] = c[c1,c2,c3,c4,1, 1] - coeff * c[c1,c2,c3,c4,1, 5];
            c[c1,c2,c3,c4,2, 1] = c[c1,c2,c3,c4,2, 1] - coeff * c[c1,c2,c3,c4,2, 5];
            c[c1,c2,c3,c4,3, 1] = c[c1,c2,c3,c4,3, 1] - coeff * c[c1,c2,c3,c4,3, 5];
            c[c1,c2,c3,c4,4, 1] = c[c1,c2,c3,c4,4, 1] - coeff * c[c1,c2,c3,c4,4, 5];
            c[c1,c2,c3,c4,5, 1] = c[c1,c2,c3,c4,5, 1] - coeff * c[c1,c2,c3,c4,5, 5];

            r[r1,r2,r3,r4,1]   = r[r1,r2,r3,r4,1]   - coeff*r[r1,r2,r3,r4,5];
            coeff = lhs[l1,5,2];
            c[c1,c2,c3,c4,1, 2] = c[c1,c2,c3,c4,1, 2] - coeff * c[c1,c2,c3,c4,1, 5];
            c[c1,c2,c3,c4,2, 2] = c[c1,c2,c3,c4,2, 2] - coeff * c[c1,c2,c3,c4,2, 5];
            c[c1,c2,c3,c4,3, 2] = c[c1,c2,c3,c4,3, 2] - coeff * c[c1,c2,c3,c4,3, 5];
            c[c1,c2,c3,c4,4, 2] = c[c1,c2,c3,c4,4, 2] - coeff * c[c1,c2,c3,c4,4, 5];
            c[c1,c2,c3,c4,5, 2] = c[c1,c2,c3,c4,5, 2] - coeff * c[c1,c2,c3,c4,5, 5];

            r[r1,r2,r3,r4,2]   = r[r1,r2,r3,r4,2]   - coeff*r[r1,r2,r3,r4,5];
            coeff = lhs[l1,5,3];
            c[c1,c2,c3,c4,1, 3] = c[c1,c2,c3,c4,1, 3] - coeff * c[c1,c2,c3,c4,1, 5];
            c[c1,c2,c3,c4,2, 3] = c[c1,c2,c3,c4,2, 3] - coeff * c[c1,c2,c3,c4,2, 5];
            c[c1,c2,c3,c4,3, 3] = c[c1,c2,c3,c4,3, 3] - coeff * c[c1,c2,c3,c4,3, 5];
            c[c1,c2,c3,c4,4, 3] = c[c1,c2,c3,c4,4, 3] - coeff * c[c1,c2,c3,c4,4, 5];
            c[c1,c2,c3,c4,5, 3] = c[c1,c2,c3,c4,5, 3] - coeff * c[c1,c2,c3,c4,5, 5];

            r[r1,r2,r3,r4,3]   = r[r1,r2,r3,r4,3]   - coeff*r[r1,r2,r3,r4,5];
            coeff = lhs[l1,5,4];
            c[c1,c2,c3,c4,1, 4] = c[c1,c2,c3,c4,1, 4] - coeff * c[c1,c2,c3,c4,1, 5];
            c[c1,c2,c3,c4,2, 4] = c[c1,c2,c3,c4,2, 4] - coeff * c[c1,c2,c3,c4,2, 5];
            c[c1,c2,c3,c4,3, 4] = c[c1,c2,c3,c4,3, 4] - coeff * c[c1,c2,c3,c4,3, 5];
            c[c1,c2,c3,c4,4, 4] = c[c1,c2,c3,c4,4, 4] - coeff * c[c1,c2,c3,c4,4, 5];
            c[c1,c2,c3,c4,5, 4] = c[c1,c2,c3,c4,5, 4] - coeff * c[c1,c2,c3,c4,5, 5];

            r[r1,r2,r3,r4,4]   = r[r1,r2,r3,r4,4]   - coeff*r[r1,r2,r3,r4,5];
        }

        public void matvec_sub(ref double[, ,] ablock, ref double[, , , ,] avec, ref double[, , , ,] bvec, int ab1, int av1, int av2, int av3, int av4, int bv1, int bv2, int bv3, int bv4) {
            //---------------------------------------------------------------------
            //     subtracts bvec=bvec - ablock*avec
            //---------------------------------------------------------------------
            //double ablock,avec,bvec
            //dimension ablock[5,5],avec[5],bvec[5]
            //---------------------------------------------------------------------
            //            rhs[i,ic,jc,kc,ccell] = rhs[i,ic,jc,kc,ccell] 
            //     $           - lhs[i,1,ablock,ia,ja,ka,acell]*
            //---------------------------------------------------------------------
            bvec[bv1, bv2, bv3, bv4, 1] = bvec[bv1, bv2, bv3, bv4, 1] - ablock[ab1, 1, 1] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 1] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 1] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 1] * avec[av1, av2, av3, av4, 4] - ablock[ab1, 5, 1] * avec[av1, av2, av3, av4, 5];
            bvec[bv1, bv2, bv3, bv4, 2] = bvec[bv1, bv2, bv3, bv4, 2] - ablock[ab1, 1, 2] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 2] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 2] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 2] * avec[av1, av2, av3, av4, 4] - ablock[ab1, 5, 2] * avec[av1, av2, av3, av4, 5];
            bvec[bv1, bv2, bv3, bv4, 3] = bvec[bv1, bv2, bv3, bv4, 3] - ablock[ab1, 1, 3] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 3] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 3] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 3] * avec[av1, av2, av3, av4, 4] - ablock[ab1, 5, 3] * avec[av1, av2, av3, av4, 5];
            bvec[bv1, bv2, bv3, bv4, 4] = bvec[bv1, bv2, bv3, bv4, 4] - ablock[ab1, 1, 4] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 4] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 4] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 4] * avec[av1, av2, av3, av4, 4] - ablock[ab1, 5, 4] * avec[av1, av2, av3, av4, 5];
            bvec[bv1, bv2, bv3, bv4, 5] = bvec[bv1, bv2, bv3, bv4, 5] - ablock[ab1, 1, 5] * avec[av1, av2, av3, av4, 1] - ablock[ab1, 2, 5] * avec[av1, av2, av3, av4, 2] - ablock[ab1, 3, 5] * avec[av1, av2, av3, av4, 3] - ablock[ab1, 4, 5] * avec[av1, av2, av3, av4, 4] - ablock[ab1, 5, 5] * avec[av1, av2, av3, av4, 5];
        }

        public void matmul_sub(ref double[, ,] ablock, ref double[, , , , ,] bblock, ref double[, ,] cblock, int a1, int b1, int b2, int b3, int b4, int c1) {
            //---------------------------------------------------------------------
            //     subtracts a[i,j,k] X b[i,j,k] from c[i,j,k]
            //---------------------------------------------------------------------
            //      double ablock, bblock, cblock
            //      dimension ablock[5,5], bblock[5,5], cblock[5,5]
            cblock[c1, 1, 1] = cblock[c1, 1, 1] - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 1, 4]
                                     - ablock[a1, 5, 1] * bblock[b1, b2, b3, b4, 1, 5];

            cblock[c1, 1, 2] = cblock[c1, 1, 2] - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 1, 4]
                                     - ablock[a1, 5, 2] * bblock[b1, b2, b3, b4, 1, 5];

            cblock[c1, 1, 3] = cblock[c1, 1, 3] - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 1, 4]
                                     - ablock[a1, 5, 3] * bblock[b1, b2, b3, b4, 1, 5];

            cblock[c1, 1, 4] = cblock[c1, 1, 4] - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 1, 4]
                                     - ablock[a1, 5, 4] * bblock[b1, b2, b3, b4, 1, 5];

            cblock[c1, 1, 5] = cblock[c1, 1, 5] - ablock[a1, 1, 5] * bblock[b1, b2, b3, b4, 1, 1]
                                     - ablock[a1, 2, 5] * bblock[b1, b2, b3, b4, 1, 2]
                                     - ablock[a1, 3, 5] * bblock[b1, b2, b3, b4, 1, 3]
                                     - ablock[a1, 4, 5] * bblock[b1, b2, b3, b4, 1, 4]
                                     - ablock[a1, 5, 5] * bblock[b1, b2, b3, b4, 1, 5];

            cblock[c1, 2, 1] = cblock[c1, 2, 1] - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 2, 4]
                                     - ablock[a1, 5, 1] * bblock[b1, b2, b3, b4, 2, 5];

            cblock[c1, 2, 2] = cblock[c1, 2, 2] - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 2, 4]
                                     - ablock[a1, 5, 2] * bblock[b1, b2, b3, b4, 2, 5];

            cblock[c1, 2, 3] = cblock[c1, 2, 3] - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 2, 4]
                                     - ablock[a1, 5, 3] * bblock[b1, b2, b3, b4, 2, 5];

            cblock[c1, 2, 4] = cblock[c1, 2, 4] - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 2, 4]
                                     - ablock[a1, 5, 4] * bblock[b1, b2, b3, b4, 2, 5];

            cblock[c1, 2, 5] = cblock[c1, 2, 5] - ablock[a1, 1, 5] * bblock[b1, b2, b3, b4, 2, 1]
                                     - ablock[a1, 2, 5] * bblock[b1, b2, b3, b4, 2, 2]
                                     - ablock[a1, 3, 5] * bblock[b1, b2, b3, b4, 2, 3]
                                     - ablock[a1, 4, 5] * bblock[b1, b2, b3, b4, 2, 4]
                                     - ablock[a1, 5, 5] * bblock[b1, b2, b3, b4, 2, 5];

            cblock[c1, 3, 1] = cblock[c1, 3, 1] - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 3, 4]
                                     - ablock[a1, 5, 1] * bblock[b1, b2, b3, b4, 3, 5];

            cblock[c1, 3, 2] = cblock[c1, 3, 2] - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 3, 4]
                                     - ablock[a1, 5, 2] * bblock[b1, b2, b3, b4, 3, 5];

            cblock[c1, 3, 3] = cblock[c1, 3, 3] - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 3, 4]
                                     - ablock[a1, 5, 3] * bblock[b1, b2, b3, b4, 3, 5];

            cblock[c1, 3, 4] = cblock[c1, 3, 4] - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 3, 4]
                                     - ablock[a1, 5, 4] * bblock[b1, b2, b3, b4, 3, 5];

            cblock[c1, 3, 5] = cblock[c1, 3, 5] - ablock[a1, 1, 5] * bblock[b1, b2, b3, b4, 3, 1]
                                     - ablock[a1, 2, 5] * bblock[b1, b2, b3, b4, 3, 2]
                                     - ablock[a1, 3, 5] * bblock[b1, b2, b3, b4, 3, 3]
                                     - ablock[a1, 4, 5] * bblock[b1, b2, b3, b4, 3, 4]
                                     - ablock[a1, 5, 5] * bblock[b1, b2, b3, b4, 3, 5];

            cblock[c1, 4, 1] = cblock[c1, 4, 1] - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 4, 4]
                                     - ablock[a1, 5, 1] * bblock[b1, b2, b3, b4, 4, 5];

            cblock[c1, 4, 2] = cblock[c1, 4, 2] - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 4, 4]
                                     - ablock[a1, 5, 2] * bblock[b1, b2, b3, b4, 4, 5];

            cblock[c1, 4, 3] = cblock[c1, 4, 3] - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 4, 4]
                                     - ablock[a1, 5, 3] * bblock[b1, b2, b3, b4, 4, 5];

            cblock[c1, 4, 4] = cblock[c1, 4, 4] - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 4, 4]
                                     - ablock[a1, 5, 4] * bblock[b1, b2, b3, b4, 4, 5];

            cblock[c1, 4, 5] = cblock[c1, 4, 5] - ablock[a1, 1, 5] * bblock[b1, b2, b3, b4, 4, 1]
                                     - ablock[a1, 2, 5] * bblock[b1, b2, b3, b4, 4, 2]
                                     - ablock[a1, 3, 5] * bblock[b1, b2, b3, b4, 4, 3]
                                     - ablock[a1, 4, 5] * bblock[b1, b2, b3, b4, 4, 4]
                                     - ablock[a1, 5, 5] * bblock[b1, b2, b3, b4, 4, 5];

            cblock[c1, 5, 1] = cblock[c1, 5, 1] - ablock[a1, 1, 1] * bblock[b1, b2, b3, b4, 5, 1]
                                     - ablock[a1, 2, 1] * bblock[b1, b2, b3, b4, 5, 2]
                                     - ablock[a1, 3, 1] * bblock[b1, b2, b3, b4, 5, 3]
                                     - ablock[a1, 4, 1] * bblock[b1, b2, b3, b4, 5, 4]
                                     - ablock[a1, 5, 1] * bblock[b1, b2, b3, b4, 5, 5];

            cblock[c1, 5, 2] = cblock[c1, 5, 2] - ablock[a1, 1, 2] * bblock[b1, b2, b3, b4, 5, 1]
                                     - ablock[a1, 2, 2] * bblock[b1, b2, b3, b4, 5, 2]
                                     - ablock[a1, 3, 2] * bblock[b1, b2, b3, b4, 5, 3]
                                     - ablock[a1, 4, 2] * bblock[b1, b2, b3, b4, 5, 4]
                                     - ablock[a1, 5, 2] * bblock[b1, b2, b3, b4, 5, 5];

            cblock[c1, 5, 3] = cblock[c1, 5, 3] - ablock[a1, 1, 3] * bblock[b1, b2, b3, b4, 5, 1]
                                     - ablock[a1, 2, 3] * bblock[b1, b2, b3, b4, 5, 2]
                                     - ablock[a1, 3, 3] * bblock[b1, b2, b3, b4, 5, 3]
                                     - ablock[a1, 4, 3] * bblock[b1, b2, b3, b4, 5, 4]
                                     - ablock[a1, 5, 3] * bblock[b1, b2, b3, b4, 5, 5];

            cblock[c1, 5, 4] = cblock[c1, 5, 4] - ablock[a1, 1, 4] * bblock[b1, b2, b3, b4, 5, 1]
                                     - ablock[a1, 2, 4] * bblock[b1, b2, b3, b4, 5, 2]
                                     - ablock[a1, 3, 4] * bblock[b1, b2, b3, b4, 5, 3]
                                     - ablock[a1, 4, 4] * bblock[b1, b2, b3, b4, 5, 4]
                                     - ablock[a1, 5, 4] * bblock[b1, b2, b3, b4, 5, 5];

            cblock[c1, 5, 5] = cblock[c1, 5, 5] - ablock[a1, 1, 5] * bblock[b1, b2, b3, b4, 5, 1]
                                     - ablock[a1, 2, 5] * bblock[b1, b2, b3, b4, 5, 2]
                                     - ablock[a1, 3, 5] * bblock[b1, b2, b3, b4, 5, 3]
                                     - ablock[a1, 4, 5] * bblock[b1, b2, b3, b4, 5, 4]
                                     - ablock[a1, 5, 5] * bblock[b1, b2, b3, b4, 5, 5];
        }

        public void binvrhs(ref double[, ,] lhs, ref double[, , , ,] r, int l1, int r1, int r2, int r3, int r4) {
            double pivot, coeff; // dimension lhs[5,5]; r[5];

            pivot = 1.00d / lhs[l1, 1, 1];
            lhs[l1, 2, 1] = lhs[l1, 2, 1] * pivot;
            lhs[l1, 3, 1] = lhs[l1, 3, 1] * pivot;
            lhs[l1, 4, 1] = lhs[l1, 4, 1] * pivot;
            lhs[l1, 5, 1] = lhs[l1, 5, 1] * pivot;
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] * pivot;

            coeff = lhs[l1, 1, 2];
            lhs[l1, 2, 2] = lhs[l1, 2, 2] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 1];
            lhs[l1, 5, 2] = lhs[l1, 5, 2] - coeff * lhs[l1, 5, 1];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 1];

            coeff = lhs[l1, 1, 3];
            lhs[l1, 2, 3] = lhs[l1, 2, 3] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 1];
            lhs[l1, 5, 3] = lhs[l1, 5, 3] - coeff * lhs[l1, 5, 1];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 1];

            coeff = lhs[l1, 1, 4];
            lhs[l1, 2, 4] = lhs[l1, 2, 4] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 1];
            lhs[l1, 5, 4] = lhs[l1, 5, 4] - coeff * lhs[l1, 5, 1];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 1];

            coeff = lhs[l1, 1, 5];
            lhs[l1, 2, 5] = lhs[l1, 2, 5] - coeff * lhs[l1, 2, 1];
            lhs[l1, 3, 5] = lhs[l1, 3, 5] - coeff * lhs[l1, 3, 1];
            lhs[l1, 4, 5] = lhs[l1, 4, 5] - coeff * lhs[l1, 4, 1];
            lhs[l1, 5, 5] = lhs[l1, 5, 5] - coeff * lhs[l1, 5, 1];
            r[r1, r2, r3, r4, 5] = r[r1, r2, r3, r4, 5] - coeff * r[r1, r2, r3, r4, 1];

            pivot = 1.00d / lhs[l1, 2, 2];
            lhs[l1, 3, 2] = lhs[l1, 3, 2] * pivot;
            lhs[l1, 4, 2] = lhs[l1, 4, 2] * pivot;
            lhs[l1, 5, 2] = lhs[l1, 5, 2] * pivot;
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] * pivot;

            coeff = lhs[l1, 2, 1];
            lhs[l1, 3, 1] = lhs[l1, 3, 1] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 2];
            lhs[l1, 5, 1] = lhs[l1, 5, 1] - coeff * lhs[l1, 5, 2];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 2];

            coeff = lhs[l1, 2, 3];
            lhs[l1, 3, 3] = lhs[l1, 3, 3] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] - coeff * lhs[l1, 4, 2];
            lhs[l1, 5, 3] = lhs[l1, 5, 3] - coeff * lhs[l1, 5, 2];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 2];

            coeff = lhs[l1, 2, 4];
            lhs[l1, 3, 4] = lhs[l1, 3, 4] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 2];
            lhs[l1, 5, 4] = lhs[l1, 5, 4] - coeff * lhs[l1, 5, 2];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 2];

            coeff = lhs[l1, 2, 5];
            lhs[l1, 3, 5] = lhs[l1, 3, 5] - coeff * lhs[l1, 3, 2];
            lhs[l1, 4, 5] = lhs[l1, 4, 5] - coeff * lhs[l1, 4, 2];
            lhs[l1, 5, 5] = lhs[l1, 5, 5] - coeff * lhs[l1, 5, 2];
            r[r1, r2, r3, r4, 5] = r[r1, r2, r3, r4, 5] - coeff * r[r1, r2, r3, r4, 2];

            pivot = 1.00d / lhs[l1, 3, 3];
            lhs[l1, 4, 3] = lhs[l1, 4, 3] * pivot;
            lhs[l1, 5, 3] = lhs[l1, 5, 3] * pivot;
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] * pivot;

            coeff = lhs[l1, 3, 1];
            lhs[l1, 4, 1] = lhs[l1, 4, 1] - coeff * lhs[l1, 4, 3];
            lhs[l1, 5, 1] = lhs[l1, 5, 1] - coeff * lhs[l1, 5, 3];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 3];

            coeff = lhs[l1, 3, 2];
            lhs[l1, 4, 2] = lhs[l1, 4, 2] - coeff * lhs[l1, 4, 3];
            lhs[l1, 5, 2] = lhs[l1, 5, 2] - coeff * lhs[l1, 5, 3];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 3];

            coeff = lhs[l1, 3, 4];
            lhs[l1, 4, 4] = lhs[l1, 4, 4] - coeff * lhs[l1, 4, 3];
            lhs[l1, 5, 4] = lhs[l1, 5, 4] - coeff * lhs[l1, 5, 3];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 3];

            coeff = lhs[l1, 3, 5];
            lhs[l1, 4, 5] = lhs[l1, 4, 5] - coeff * lhs[l1, 4, 3];
            lhs[l1, 5, 5] = lhs[l1, 5, 5] - coeff * lhs[l1, 5, 3];
            r[r1, r2, r3, r4, 5] = r[r1, r2, r3, r4, 5] - coeff * r[r1, r2, r3, r4, 3];

            pivot = 1.00d / lhs[l1, 4, 4];
            lhs[l1, 5, 4] = lhs[l1, 5, 4] * pivot;
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] * pivot;

            coeff = lhs[l1, 4, 1];
            lhs[l1, 5, 1] = lhs[l1, 5, 1] - coeff * lhs[l1, 5, 4];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 4];

            coeff = lhs[l1, 4, 2];
            lhs[l1, 5, 2] = lhs[l1, 5, 2] - coeff * lhs[l1, 5, 4];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 4];

            coeff = lhs[l1, 4, 3];
            lhs[l1, 5, 3] = lhs[l1, 5, 3] - coeff * lhs[l1, 5, 4];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 4];

            coeff = lhs[l1, 4, 5];
            lhs[l1, 5, 5] = lhs[l1, 5, 5] - coeff * lhs[l1, 5, 4];
            r[r1, r2, r3, r4, 5] = r[r1, r2, r3, r4, 5] - coeff * r[r1, r2, r3, r4, 4];

            pivot = 1.00d / lhs[l1, 5, 5];
            r[r1, r2, r3, r4, 5] = r[r1, r2, r3, r4, 5] * pivot;

            coeff = lhs[l1, 5, 1];
            r[r1, r2, r3, r4, 1] = r[r1, r2, r3, r4, 1] - coeff * r[r1, r2, r3, r4, 5];

            coeff = lhs[l1, 5, 2];
            r[r1, r2, r3, r4, 2] = r[r1, r2, r3, r4, 2] - coeff * r[r1, r2, r3, r4, 5];

            coeff = lhs[l1, 5, 3];
            r[r1, r2, r3, r4, 3] = r[r1, r2, r3, r4, 3] - coeff * r[r1, r2, r3, r4, 5];

            coeff = lhs[l1, 5, 4];
            r[r1, r2, r3, r4, 4] = r[r1, r2, r3, r4, 4] - coeff * r[r1, r2, r3, r4, 5];
        }
                   // end solve_subs.f
        // end x_solve.f

        // y_solve.f
        public void y_solve() {
            //---------------------------------------------------------------------
            //     Performs line solves in Y direction by first factoring
            //     the block-tridiagonal matrix into an upper triangular matrix, 
            //     and { performing back substitution to solve for the unknow
            //     vectors of each line.  
            //     
            //     Make sure we treat elements zero to cell_size in the direction
            //     of the sweep.
            //---------------------------------------------------------------------
            int c, stage, first, last, isize,jsize,ksize,buffer_size;//int r_status[MPI_STATUS_SIZE];
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE);
            MPI.Request[] recv_id = new MPI.Request[1];
            MPI.Request[] send_id = new MPI.Request[1];
            double[] out_buffer_x = new double[buffer_size];
            //---------------------------------------------------------------------
            //     in our terminology stage is the number of the cell in the y-direction
            //     i.e. stage = 1 means the start of the line stage=ncells means end
            //---------------------------------------------------------------------
            for(stage = 1; stage<=ncells; stage++){
               c = slice[2,stage];
               isize = cell_size[1,c] - 1;
               jsize = cell_size[2,c] - 1;
               ksize = cell_size[3,c] - 1;
               //---------------------------------------------------------------------
               //     set last-cell flag
               //---------------------------------------------------------------------
               if (stage == ncells) {
                  last = 1;
               } else {
                  last = 0;
               }
               if (stage == 1) {
                  //---------------------------------------------------------------------
                  //     This is the first cell, so solve without receiving data
                  //---------------------------------------------------------------------
                  first = 1;
                        //c            call lhsy[c]
                  y_solve_cell(first,last,c); //call y_solve_cell[first,last,c];
               } else {
                  //c---------------------------------------------------------------------
                  //c     Not the first cell of this line, so receive info from
                  //c     processor working on preceeding cell
                  //c---------------------------------------------------------------------
                  first = 0;
                  y_receive_solve_info(out_buffer_x, recv_id, c); //call y_receive_solve_info[recv_id,c];
                  //      c---------------------------------------------------------------------
                  //      c     overlap computations and communications
                  //      c---------------------------------------------------------------------
                  //      c            call lhsy[c]
                  //      c---------------------------------------------------------------------
                  //      c     wait for completion
                  //      c---------------------------------------------------------------------
                  send_id[0].Wait();
                  recv_id[0].Wait();
                  //call mpi_wait[send_id,r_status,error];
                  //call mpi_wait[recv_id,r_status,error];
                  //      c---------------------------------------------------------------------
                  //      c     install C'[jstart+1] and rhs'[jstart+1] to be used in this cell
                  //      c---------------------------------------------------------------------
                  y_unpack_solve_info(out_buffer_x, c);
                  y_solve_cell(first,last,c); //call y_solve_cell[first,last,c];
               }
               if (last == 0) y_send_solve_info(send_id, c);  //call y_send_solve_info[send_id,c];
            }
            out_buffer_x = null;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE;
            out_buffer_x = new double[buffer_size];
            //---------------------------------------------------------------------
            //     now perform backsubstitution in reverse direction
            //---------------------------------------------------------------------
            for(stage = ncells; stage>= 1; stage--){  //for(stage = ncells, 1, -1
               c = slice[2,stage];
               first = 0;
               last = 0;
               if (stage == 1) first = 1;
               if (stage == ncells) {
                  last = 1;
                  //---------------------------------------------------------------------
                  //     last cell, so perform back substitute without waiting
                  //---------------------------------------------------------------------
                  y_backsubstitute(first, last,c);     //call y_backsubstitute[first, last,c];
               } else {
                  y_receive_backsub_info(out_buffer_x, recv_id, c);   //call y_receive_backsub_info[recv_id,c];
                  send_id[0].Wait();
                  recv_id[0].Wait();
                          //call mpi_wait[send_id,r_status,error];
                          //call mpi_wait[recv_id,r_status,error];
                  y_unpack_backsub_info(out_buffer_x, c);
                  y_backsubstitute(first,last,c);      //call y_backsubstitute[first,last,c];
               }
               if (first == 0) y_send_backsub_info(send_id, c);  //call y_send_backsub_info[send_id,c];
            }
        }

        public void y_unpack_solve_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack C'[-1] and rhs'[-1] for
            //     all i and k
            //---------------------------------------------------------------------
            int i,k,m,n,ptr,jstart;
            jstart = 0;
            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
               for(i=0; i<=IMAX-1; i++){
                  for(m=1; m<=BLOCK_SIZE; m++){
                     for(n=0; n<BLOCK_SIZE; n++){
                        lhsc[c,k+1,jstart,i+1,n+1,m] = out_buffer_x[ptr+n]; //lhsc[m,n,i,jstart-1,k,c] = out_buffer_x[ptr+n];
                     }
                     ptr = ptr+BLOCK_SIZE;
                  }
                  for(n=0; n<BLOCK_SIZE; n++){
                     rhs[c,k+1,jstart,i+1,n+1] = out_buffer_x[ptr+n];  // rhs[n,i,jstart-1,k,c] = out_buffer_x[ptr+n];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
        }

        public void y_send_solve_info(MPI.Request[] send_id, int c) {
            //---------------------------------------------------------------------
            //     pack up and send C'[jend] and rhs'[jend] for
            //     all i and k
            //---------------------------------------------------------------------
            int i,k,m,n,jsize,ptr,ip,kp;
            int buffer_size;

            jsize = cell_size[2,c]-1;
            ip = cell_coord[1,c] - 1;
            kp = cell_coord[3,c] - 1;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE);
            double[] in_buffer_x = new double[buffer_size];

            //---------------------------------------------------------------------
            //     pack up buffer
            //---------------------------------------------------------------------
            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
               for(i=0; i<=IMAX-1; i++){
                  for(m=1; m<=BLOCK_SIZE; m++){
                     for(n=0; n<BLOCK_SIZE; n++){
                        in_buffer_x[ptr+n] = lhsc[c,k+1,jsize+1,i+1,n+1,m];  //in_buffer[ptr+n] = lhsc[m,n,i,jsize,k,c]
                     }
                     ptr = ptr+BLOCK_SIZE;
                  }
                  for(n=0; n<BLOCK_SIZE; n++){
                     in_buffer_x[ptr+n] = rhs[c,k+1,jsize+1,i+1,n+1];       //in_buffer[ptr+n] = rhs[n,i,jsize,k,c] 
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
            //---------------------------------------------------------------------
            //     send buffer 
            //---------------------------------------------------------------------
            //call mpi_isend[in_buffer, buffer_size, dp_type, successor[2], SOUTH+ip+kp*NCELLS, comm_solve, send_id,error];

            send_id[0] = comm_solve.ImmediateSend<double>(in_buffer_x, successor[2], SOUTH+ip+kp*ncells);
        }

        public void y_send_backsub_info(MPI.Request[] send_id, int c) {
            //---------------------------------------------------------------------
            //     pack up and send U[jstart] for all i and k
            //---------------------------------------------------------------------
            int i,k,n,ptr,jstart,ip,kp;
            int buffer_size;
            //---------------------------------------------------------------------
            //     Send element 0 to previous processor
            //---------------------------------------------------------------------
            jstart = 0;
            ip = cell_coord[1,c]-1;
            kp = cell_coord[3,c]-1;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE;

            double[] in_buffer_x = new double[buffer_size];

            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
               for(i=0; i<=IMAX-1; i++){
                  for(n=0; n<BLOCK_SIZE; n++){
                     in_buffer_x[ptr+n] = rhs[c,k+1,jstart+1,i+1,n+1];  //in_buffer[ptr+n] = rhs[n,i,jstart,k,c];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
            //call mpi_isend[in_buffer, buffer_size, dp_type, predecessor[2], NORTH+ip+kp*NCELLS, comm_solve, send_id,error];

            send_id[0] = comm_solve.ImmediateSend<double>(in_buffer_x, predecessor[2], NORTH+ip+kp*ncells);
        }

        public void y_unpack_backsub_info(double[] out_buffer_x, int c) {
            //  c---------------------------------------------------------------------
            //  c     unpack U[jsize] for all i and k
            //  c---------------------------------------------------------------------
            int i,k,n,ptr;

            ptr = 0;
            for(k=0; k<=KMAX-1; k++){
               for(i=0; i<=IMAX-1; i++){
                  for(n=0; n<BLOCK_SIZE; n++){
                     backsub_info[c,k,i,n+1] = out_buffer_x[ptr+n];  //backsub_info[n,i,k,c] = out_buffer[ptr+n];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
        }

        public void y_receive_backsub_info(double[] out_buffer_x,MPI.Request[] recv_id,int c) {
            //---------------------------------------------------------------------
            //     post mpi receives
            //---------------------------------------------------------------------
            int ip,kp;
            ip = cell_coord[1,c] - 1;
            kp = cell_coord[3,c] - 1;
            //call mpi_irecv[out_buffer, buffer_size, dp_type, successor[2], NORTH+ip+kp*NCELLS, comm_solve, recv_id, error];

            recv_id[0] = comm_solve.ImmediateReceive<double>(successor[2], NORTH+ip+kp*ncells, out_buffer_x);
        }

        public void y_receive_solve_info(double[] out_buffer_x, MPI.Request[] recv_id, int c) {
            //---------------------------------------------------------------------
            //     post mpi receives 
            //---------------------------------------------------------------------
            int ip,kp;
            ip = cell_coord[1,c] - 1;
            kp = cell_coord[3,c] - 1;
            //call mpi_irecv[out_buffer, buffer_size, dp_type, predecessor[2], SOUTH+ip+kp*NCELLS,  comm_solve, recv_id, error];

            recv_id[0] = comm_solve.ImmediateReceive<double>(predecessor[2], SOUTH+ip+kp*ncells, out_buffer_x);
        }

        public void y_backsubstitute(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     back solve: if last cell, { generate U[jsize]=rhs[jsize]
            //     } else { assume U[jsize] is loaded in un pack backsub_info
            //     so just use it
            //     after call u[jstart] will be sent to next cell
            //---------------------------------------------------------------------
            int i, k;
            int m,n,j,jsize,isize,ksize,jstart;
            jstart = 0;
            isize = cell_size[1,c]-end[1,c]-1;
            jsize = cell_size[2,c]-1;
            ksize = cell_size[3,c]-end[3,c]-1;
            if (last == 0) {
               for(k=start[3,c]; k<=ksize; k++){
                  for(i=start[1,c]; i<=isize; i++){
                     //---------------------------------------------------------------------
                     //     U[jsize] uses info from previous cell if not last cell
                     //---------------------------------------------------------------------
                     for(m=1; m<=BLOCK_SIZE; m++){
                        for(n=1; n<=BLOCK_SIZE; n++){ //rhs[m,i,jsize,k,c]=rhs[m,i,jsize,k,c]-lhsc[m,n,i,jsize,k,c]*backsub_info[n,i,k,c];
                           rhs[c,k+1,jsize+1,i+1,m] = rhs[c,k+1,jsize+1,i+1,m] - lhsc[c,k+1,jsize+1,i+1,n,m]*backsub_info[c,k,i,n];
                        }
                     }
                  }
               }
            }
            for(k=start[3,c]; k<=ksize; k++){
               for(j=jsize-1; j>=jstart; j--){//for(j=jsize-1,jstart,-1;
                  for(i=start[1,c]; i<=isize; i++){
                     for(m=1; m<=BLOCK_SIZE; m++){
                        for(n=1; n<=BLOCK_SIZE; n++){  //rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - lhsc[m,n,i,j,k,c]*rhs[n,i,j+1,k,c];
                           rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - lhsc[c,k+1,j+1,i+1,n,m]*rhs[c,k+1,j+2,i+1,n];
                        }
                     }
                  }
               }
            }
        }

        public void y_solve_cell(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     performs guaussian elimination on this cell.
            //     
            //     assumes that unpacking routines for non-first cells 
            //     preload C' and rhs' from previous cell.
            //     
            //     assumed send happens outside this routine, but that
            //     c'[JMAX] and rhs'[JMAX] will be sent to next cell
            //---------------------------------------------------------------------
            int i,j,k,isize,ksize,jsize,jstart;
            double[,] utmp = new double[JMAX+4,7];   //double utmp[6,-2:JMAX+1];

            jstart = 0;
            isize = cell_size[1,c]-end[1,c]-1;
            jsize = cell_size[2,c]-1;
            ksize = cell_size[3,c]-end[3,c]-1;

            lhsabinit(ref lhsa, ref lhsb, jsize);

            for(k=start[3,c]; k<=ksize; k++){
               for(i=start[1,c]; i<=isize; i++){
                    //---------------------------------------------------------------------
                    //     This function computes the left hand side for the three y-factors   
                    //---------------------------------------------------------------------
                    //---------------------------------------------------------------------
                    //     Compute the indices for storing the tri-diagonal matrix;
                    //     determine a [labeled f] and n jacobians for cell c
                    //---------------------------------------------------------------------
                    for(j = start[2,c]-1; j<= cell_size[2,c]-end[2,c]; j++){
                        utmp[2+j,1] = 1.0d / u[c,k+2,j+2,i+2,1];  //u[1,i,j,k,c]
                        utmp[2+j,2] = u[c,k+2,j+2,i+2,2];         //u[2,i,j,k,c]; 
                        utmp[2+j,3] = u[c,k+2,j+2,i+2,3];         //u[3,i,j,k,c]; 
                        utmp[2+j,4] = u[c,k+2,j+2,i+2,4];         //u[4,i,j,k,c]; 
                        utmp[2+j,5] = u[c,k+2,j+2,i+2,5];         //u[5,i,j,k,c]; 
                        utmp[2+j,6] =qs[c,k+1,j+1,i+1];           //qs[i,j,k,c]; 
                    }

                    for(j = start[2,c]-1; j<= cell_size[2,c]-end[2,c]; j++){
                        tmp1 = utmp[2+j,1];
                        tmp2 = tmp1 * tmp1;
                        tmp3 = tmp1 * tmp2;

                        fjac[2+j,1,1] = 0.0d;
                        fjac[2+j,2,1] = 0.0d;
                        fjac[2+j,3,1] = 1.0d;
                        fjac[2+j,4,1] = 0.0d;
                        fjac[2+j,5,1] = 0.0d;

                        fjac[2+j,1,2] = - ( utmp[2+j,2]*utmp[2+j,3] ) * tmp2;
                        fjac[2+j,2,2] = utmp[2+j,3] * tmp1;
                        fjac[2+j,3,2] = utmp[2+j,2] * tmp1;
                        fjac[2+j,4,2] = 0.0d;
                        fjac[2+j,5,2] = 0.0d;

                        fjac[2+j,1,3] = - ( utmp[2+j,3]*utmp[2+j,3]*tmp2) + c2 * utmp[2+j,6];
                        fjac[2+j,2,3] = - c2 *  utmp[2+j,2] * tmp1;
                        fjac[2+j,3,3] = ( 2.0d - c2 ) *  utmp[2+j,3] * tmp1;
                        fjac[2+j,4,3] = - c2 * utmp[2+j,4] * tmp1;
                        fjac[2+j,5,3] = c2;

                        fjac[2+j,1,4] = - ( utmp[2+j,3]*utmp[2+j,4] ) * tmp2;
                        fjac[2+j,2,4] = 0.0d;
                        fjac[2+j,3,4] = utmp[2+j,4] * tmp1;
                        fjac[2+j,4,4] = utmp[2+j,3] * tmp1;
                        fjac[2+j,5,4] = 0.0d;

                        fjac[2+j,1,5] = ( c2 * 2.0d * utmp[2+j,6] - c1 * utmp[2+j,5] * tmp1 ) * utmp[2+j,3] * tmp1;
                        fjac[2+j,2,5] = - c2 * utmp[2+j,2]*utmp[2+j,3] * tmp2;
                        fjac[2+j,3,5] = c1 * utmp[2+j,5] * tmp1 - c2 * ( utmp[2+j,6] + utmp[2+j,3]*utmp[2+j,3] * tmp2 );
                        fjac[2+j,4,5] = - c2 * ( utmp[2+j,3]*utmp[2+j,4] ) * tmp2;
                        fjac[2+j,5,5] = c1 * utmp[2+j,3] * tmp1;

                        njac[2+j,1,1] = 0.0d;
                        njac[2+j,2,1] = 0.0d;
                        njac[2+j,3,1] = 0.0d;
                        njac[2+j,4,1] = 0.0d;
                        njac[2+j,5,1] = 0.0d;

                        njac[2+j,1,2] = - c3c4 * tmp2 * utmp[2+j,2];
                        njac[2+j,2,2] =   c3c4 * tmp1;
                        njac[2+j,3,2] =   0.0d;
                        njac[2+j,4,2] =   0.0d;
                        njac[2+j,5,2] =   0.0d;

                        njac[2+j,1,3] = - con43 * c3c4 * tmp2 * utmp[2+j,3];
                        njac[2+j,2,3] =   0.0d;
                        njac[2+j,3,3] =   con43 * c3c4 * tmp1;
                        njac[2+j,4,3] =   0.0d;
                        njac[2+j,5,3] =   0.0d;

                        njac[2+j,1,4] = - c3c4 * tmp2 * utmp[2+j,4];
                        njac[2+j,2,4] =   0.0d;
                        njac[2+j,3,4] =   0.0d;
                        njac[2+j,4,4] =   c3c4 * tmp1;
                        njac[2+j,5,4] =   0.0d; 

                        njac[2+j,1,5] = - (c3c4-c1345)*tmp3*(pow2(utmp[2+j,2]))-
                                        (con43*c3c4-c1345)*tmp3*(pow2(utmp[2+j,3]))-(c3c4-c1345)*tmp3*(pow2(utmp[2+j,4]))-c1345*tmp2*utmp[2+j,5];

                        njac[2+j,2,5] = (  c3c4 - c1345 ) * tmp2 * utmp[2+j,2];
                        njac[2+j,3,5] = ( con43 * c3c4 - c1345 ) * tmp2 * utmp[2+j,3];
                        njac[2+j,4,5] = ( c3c4 - c1345 ) * tmp2 * utmp[2+j,4];
                        njac[2+j,5,5] = ( c1345 ) * tmp1;
                    }
                    //---------------------------------------------------------------------
                    //     now joacobians set, so form left hand side in y direction
                    //---------------------------------------------------------------------
                    for(j = start[2,c]; j<= jsize-end[2,c]; j++){
                        tmp1 = dt * ty1;
                        tmp2 = dt * ty2;//lhsa[1,1,j]=-tmp2*fjac[1,1,j-1]-tmp1*njac[1,1,j-1]-tmp1*dy1;
                        lhsa[1+j,1,1] = - tmp2 * fjac[j+1,1,1] - tmp1 * njac[j+1,1,1] - tmp1 * dy1;                        
                        lhsa[1+j,2,1] = - tmp2 * fjac[j+1,2,1] - tmp1 * njac[j+1,2,1];//lhsa[1,2,j]=-tmp2*fjac[1,2,j-1]-tmp1*njac[1,2,j-1];
                        lhsa[1+j,3,1] = - tmp2 * fjac[j+1,3,1] - tmp1 * njac[j+1,3,1];
                        lhsa[1+j,4,1] = - tmp2 * fjac[j+1,4,1] - tmp1 * njac[j+1,4,1];
                        lhsa[1+j,5,1] = - tmp2 * fjac[j+1,5,1] - tmp1 * njac[j+1,5,1]; //Obs: fjac[2+j-1,5,1] - tmp1 * njac[2+j-1,5,1];

                        lhsa[1+j,1,2] = - tmp2 * fjac[j+1,1,2] - tmp1 * njac[j+1,1,2];
                        lhsa[1+j,2,2] = - tmp2 * fjac[j+1,2,2] - tmp1 * njac[j+1,2,2] - tmp1 * dy2;
                        lhsa[1+j,3,2] = - tmp2 * fjac[j+1,3,2] - tmp1 * njac[j+1,3,2];
                        lhsa[1+j,4,2] = - tmp2 * fjac[j+1,4,2] - tmp1 * njac[j+1,4,2];
                        lhsa[1+j,5,2] = - tmp2 * fjac[j+1,5,2] - tmp1 * njac[j+1,5,2];

                        lhsa[1+j,1,3] = - tmp2 * fjac[j+1,1,3] - tmp1 * njac[j+1,1,3];
                        lhsa[1+j,2,3] = - tmp2 * fjac[j+1,2,3] - tmp1 * njac[j+1,2,3];
                        lhsa[1+j,3,3] = - tmp2 * fjac[j+1,3,3] - tmp1 * njac[j+1,3,3] - tmp1 * dy3;
                        lhsa[1+j,4,3] = - tmp2 * fjac[j+1,4,3] - tmp1 * njac[j+1,4,3];
                        lhsa[1+j,5,3] = - tmp2 * fjac[j+1,5,3] - tmp1 * njac[j+1,5,3];

                        lhsa[1+j,1,4] = - tmp2 * fjac[j+1,1,4] - tmp1 * njac[j+1,1,4];
                        lhsa[1+j,2,4] = - tmp2 * fjac[j+1,2,4] - tmp1 * njac[j+1,2,4];
                        lhsa[1+j,3,4] = - tmp2 * fjac[j+1,3,4] - tmp1 * njac[j+1,3,4];
                        lhsa[1+j,4,4] = - tmp2 * fjac[j+1,4,4] - tmp1 * njac[j+1,4,4] - tmp1 * dy4;
                        lhsa[1+j,5,4] = - tmp2 * fjac[j+1,5,4] - tmp1 * njac[j+1,5,4];

                        lhsa[1+j,1,5] = - tmp2 * fjac[j+1,1,5] - tmp1 * njac[j+1,1,5];
                        lhsa[1+j,2,5] = - tmp2 * fjac[j+1,2,5] - tmp1 * njac[j+1,2,5];
                        lhsa[1+j,3,5] = - tmp2 * fjac[j+1,3,5] - tmp1 * njac[j+1,3,5];
                        lhsa[1+j,4,5] = - tmp2 * fjac[j+1,4,5] - tmp1 * njac[j+1,4,5];
                        lhsa[1+j,5,5] = - tmp2 * fjac[j+1,5,5] - tmp1 * njac[j+1,5,5] - tmp1 * dy5;

                        lhsb[1+j,1,1] = 1.0d + tmp1 * 2.0d * njac[2+j,1,1] + tmp1 * 2.0d * dy1;
                        lhsb[1+j,2,1] = tmp1 * 2.0d * njac[2+j,2,1];  //lhsb[1,2,j] = tmp1 * 2.0d * njac[1,2,j];
                        lhsb[1+j,3,1] = tmp1 * 2.0d * njac[2+j,3,1];  //lhsb[1,3,j] = tmp1 * 2.0d * njac[1,3,j];
                        lhsb[1+j,4,1] = tmp1 * 2.0d * njac[2+j,4,1];  //lhsb[1,4,j] = tmp1 * 2.0d * njac[1,4,j];
                        lhsb[1+j,5,1] = tmp1 * 2.0d * njac[2+j,5,1];  //lhsb[1,5,j] = tmp1 * 2.0d * njac[1,5,j];

                        lhsb[1+j,1,2] = tmp1 * 2.0d * njac[2+j,1,2];
                        lhsb[1+j,2,2] = 1.0d + tmp1 * 2.0d * njac[2+j,2,2] + tmp1 * 2.0d * dy2;
                        lhsb[1+j,3,2] = tmp1 * 2.0d * njac[2+j,3,2];
                        lhsb[1+j,4,2] = tmp1 * 2.0d * njac[2+j,4,2];
                        lhsb[1+j,5,2] = tmp1 * 2.0d * njac[2+j,5,2];

                        lhsb[1+j,1,3] = tmp1 * 2.0d * njac[2+j,1,3];
                        lhsb[1+j,2,3] = tmp1 * 2.0d * njac[2+j,2,3];
                        lhsb[1+j,3,3] = 1.0d + tmp1 * 2.0d * njac[2+j,3,3] + tmp1 * 2.0d * dy3;
                        lhsb[1+j,4,3] = tmp1 * 2.0d * njac[2+j,4,3];
                        lhsb[1+j,5,3] = tmp1 * 2.0d * njac[2+j,5,3];

                        lhsb[1+j,1,4] = tmp1 * 2.0d * njac[2+j,1,4];
                        lhsb[1+j,2,4] = tmp1 * 2.0d * njac[2+j,2,4];
                        lhsb[1+j,3,4] = tmp1 * 2.0d * njac[2+j,3,4];
                        lhsb[1+j,4,4] = 1.0d + tmp1 * 2.0d * njac[2+j,4,4] + tmp1 * 2.0d * dy4;
                        lhsb[1+j,5,4] = tmp1 * 2.0d * njac[2+j,5,4];

                        lhsb[1+j,1,5] = tmp1 * 2.0d * njac[2+j,1,5];
                        lhsb[1+j,2,5] = tmp1 * 2.0d * njac[2+j,2,5];
                        lhsb[1+j,3,5] = tmp1 * 2.0d * njac[2+j,3,5];
                        lhsb[1+j,4,5] = tmp1 * 2.0d * njac[2+j,4,5];
                        lhsb[1+j,5,5] = 1.0d + tmp1 * 2.0d * njac[2+j,5,5] + tmp1 * 2.0d * dy5;

                        lhsc[c,k+1,j+1,i+1,1,1] =  tmp2 * fjac[j+3,1,1] - tmp1 * njac[j+3,1,1] - tmp1 * dy1;
                        lhsc[c,k+1,j+1,i+1,2,1] =  tmp2 * fjac[j+3,2,1] - tmp1 * njac[j+3,2,1];//lhsc[1,2,i,j,k,c]=tmp2*fjac[1,2,j+1]-tmp1*njac[1,2,j+1];
                        lhsc[c,k+1,j+1,i+1,3,1] =  tmp2 * fjac[j+3,3,1] - tmp1 * njac[j+3,3,1];//lhsc[1,3,i,j,k,c]=tmp2*fjac[1,3,j+1]-tmp1*njac[1,3,j+1];
                        lhsc[c,k+1,j+1,i+1,4,1] =  tmp2 * fjac[j+3,4,1] - tmp1 * njac[j+3,4,1];//lhsc[1,4,i,j,k,c]=tmp2*fjac[1,4,j+1]-tmp1*njac[1,4,j+1];
                        lhsc[c,k+1,j+1,i+1,5,1] =  tmp2 * fjac[j+3,5,1] - tmp1 * njac[j+3,5,1];//lhsc[1,5,i,j,k,c]=tmp2*fjac[1,5,j+1]-tmp1*njac[1,5,j+1];

                        lhsc[c,k+1,j+1,i+1,1,2] =  tmp2 * fjac[j+3,1,2] - tmp1 * njac[j+3,1,2];
                        lhsc[c,k+1,j+1,i+1,2,2] =  tmp2 * fjac[j+3,2,2] - tmp1 * njac[j+3,2,2] - tmp1 * dy2;
                        lhsc[c,k+1,j+1,i+1,3,2] =  tmp2 * fjac[j+3,3,2] - tmp1 * njac[j+3,3,2];
                        lhsc[c,k+1,j+1,i+1,4,2] =  tmp2 * fjac[j+3,4,2] - tmp1 * njac[j+3,4,2];
                        lhsc[c,k+1,j+1,i+1,5,2] =  tmp2 * fjac[j+3,5,2] - tmp1 * njac[j+3,5,2];

                        lhsc[c,k+1,j+1,i+1,1,3] =  tmp2 * fjac[j+3,1,3] - tmp1 * njac[j+3,1,3];
                        lhsc[c,k+1,j+1,i+1,2,3] =  tmp2 * fjac[j+3,2,3] - tmp1 * njac[j+3,2,3];
                        lhsc[c,k+1,j+1,i+1,3,3] =  tmp2 * fjac[j+3,3,3] - tmp1 * njac[j+3,3,3] - tmp1 * dy3;
                        lhsc[c,k+1,j+1,i+1,4,3] =  tmp2 * fjac[j+3,4,3] - tmp1 * njac[j+3,4,3];
                        lhsc[c,k+1,j+1,i+1,5,3] =  tmp2 * fjac[j+3,5,3] - tmp1 * njac[j+3,5,3];

                        lhsc[c,k+1,j+1,i+1,1,4] =  tmp2 * fjac[j+3,1,4] - tmp1 * njac[j+3,1,4];
                        lhsc[c,k+1,j+1,i+1,2,4] =  tmp2 * fjac[j+3,2,4] - tmp1 * njac[j+3,2,4];
                        lhsc[c,k+1,j+1,i+1,3,4] =  tmp2 * fjac[j+3,3,4] - tmp1 * njac[j+3,3,4];
                        lhsc[c,k+1,j+1,i+1,4,4] =  tmp2 * fjac[j+3,4,4] - tmp1 * njac[j+3,4,4] - tmp1 * dy4;
                        lhsc[c,k+1,j+1,i+1,5,4] =  tmp2 * fjac[j+3,5,4] - tmp1 * njac[j+3,5,4];

                        lhsc[c,k+1,j+1,i+1,1,5] =  tmp2 * fjac[j+3,1,5] - tmp1 * njac[j+3,1,5];
                        lhsc[c,k+1,j+1,i+1,2,5] =  tmp2 * fjac[j+3,2,5] - tmp1 * njac[j+3,2,5];
                        lhsc[c,k+1,j+1,i+1,3,5] =  tmp2 * fjac[j+3,3,5] - tmp1 * njac[j+3,3,5];
                        lhsc[c,k+1,j+1,i+1,4,5] =  tmp2 * fjac[j+3,4,5] - tmp1 * njac[j+3,4,5];
                        lhsc[c,k+1,j+1,i+1,5,5] =  tmp2 * fjac[j+3,5,5] - tmp1 * njac[j+3,5,5] - tmp1 * dy5;
                  }
                  //---------------------------------------------------------------------
                  //     outer most for(loops - sweeping in i direction
                  //---------------------------------------------------------------------
                  if (first == 1) { 
                        //c---------------------------------------------------------------------
                        //c     multiply c[i,jstart,k] by b_inverse and copy back to c
                        //c     multiply rhs[jstart] by b_inverse[jstart] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvcrhs[ lhsb[1,1,jstart],lhsc[1,1,i,jstart,k,c], rhs[1,i,jstart,k,c] ];
                        //C#:           binvcrhs(lhsb[jstart+1, 1, 1], lhsc[c, k+1, jstart+1, i+1, 1, 1], rhs[c, k+1, jstart+1, i+1, 1]);
                        binvcrhs(ref lhsb,ref lhsc,ref rhs,      (jstart+1),     c,(k+1),(jstart+1),(i+1),     c,(k+1),(jstart+1),(i+1));
                  }
                  //---------------------------------------------------------------------
                  //     begin inner most for(loop
                  //     for(all the elements of the cell unless last 
                  //---------------------------------------------------------------------
                  for(j=jstart+first; j<=jsize-last; j++){
                        //c---------------------------------------------------------------------
                        //c     subtract A*lhs_vector[j-1] from lhs_vector[j]
                        //c     
                        //c     rhs[j] = rhs[j] - A*rhs[j-1]
                        //c---------------------------------------------------------------------
                        //Fortran: call matvec_sub[lhsa[1,1,j], rhs[1,i,j-1,k,c],rhs[1,i,j,k,c]];
                        //C#:           matvec_sub(lhsa[j+1,1,1], rhs[c,k+1,j,i+1,1],rhs[c,k+1,j+1,i+1,1]);
                        matvec_sub(ref lhsa,ref rhs,ref rhs,         (j+1),     c,(k+1),j,(i+1),      c,(k+1),(j+1),(i+1));
                        //c---------------------------------------------------------------------
                        //c     B[j] = B[j] - C[j-1]*A[j]
                        //c---------------------------------------------------------------------
                        //call matmul_sub[lhsa[1,1,j],lhsc[1,1,i,j-1,k,c],lhsb[1,1,j]];
                        //C#:  matmul_sub(lhsa[j+1,1,1],lhsc[c,k+1,j-1+1,i+1,1,1],lhsb[j+1,1,1]);
                        matmul_sub(ref lhsa,ref lhsc,ref lhsb,      (j+1),      c,(k+1),j,(i+1),      (j+1));
                        //c---------------------------------------------------------------------
                        //c     multiply c[i,j,k] by b_inverse and copy back to c
                        //c     multiply rhs[i,1,k] by b_inverse[i,1,k] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvcrhs[ lhsb[1,1,j],lhsc[1,1,i,j,k,c], rhs[1,i,j,k,c] ];
                        //C#:           binvcrhs( lhsb[j+1,1,1],lhsc[c,k+1,j+1,i+1,1,1], rhs[c,k+1,j+1,i+1,1] );
                        binvcrhs(ref lhsb,ref lhsc,ref rhs,          (j+1),       c,(k+1),(j+1),(i+1),        c,(k+1),(j+1),(i+1));
                  }
                  //c---------------------------------------------------------------------
                  //c     Now finish up special cases for last cell
                  //c---------------------------------------------------------------------
                  if (last == 1) {
                        //c---------------------------------------------------------------------
                        //c     rhs[jsize] = rhs[jsize] - A*rhs[jsize-1]
                        //c---------------------------------------------------------------------
                        //Fortran: call matvec_sub[lhsa[1,1,jsize], rhs[1,i,jsize-1,k,c],rhs[1,i,jsize,k,c]];
                        //C#:           matvec_sub(lhsa[jsize+1,1,1], rhs[c,k+1,jsize-1+1,i+1,1],rhs[c,k+1,jsize+1,i+1,1]);
                        matvec_sub(ref lhsa,ref rhs,ref rhs,         (jsize+1),        c,(k+1),jsize,(i+1),      c,(k+1),(jsize+1),(i+1));
                        //c---------------------------------------------------------------------
                        //c     B[jsize] = B[jsize] - C[jsize-1]*A[jsize]
                        //c     call matmul_sub[aa,i,jsize,k,c,
                        //c     $              cc,i,jsize-1,k,c,bb,i,jsize,k,c]
                        //c---------------------------------------------------------------------
                        //Fortran: call matmul_sub[lhsa[1,1,jsize], lhsc[1,1,i,jsize-1,k,c], lhsb[1,1,jsize]];
                        //C#:           matmul_sub(lhsa[jsize+1,1,1], lhsc[c,k+1,jsize-1+1,i+1,1,1], lhsb[jsize+1,1,1]);
                        matmul_sub(ref lhsa,ref lhsc,ref lhsb,     (jsize+1),       c,(k+1),jsize,(i+1),   (jsize+1));
                        //c---------------------------------------------------------------------
                        //c     multiply rhs[jsize] by b_inverse[jsize] and copy to rhs
                        //c---------------------------------------------------------------------
                        //Fortran: call binvrhs[ lhsb[1,1,jsize],rhs[1,i,jsize,k,c] ];
                        //C#:           binvrhs(lhsb[jsize+1,1,1],rhs[c,k+1,jsize+1,i+1,1] );
                        binvrhs(ref lhsb, ref rhs,     (jsize+1),   c,(k+1),(jsize+1),(i+1));
                  }
               }
            }
        }
        // end y_solve.f

        // z_solve.f
        public void z_solve() {
            //---------------------------------------------------------------------
            //     Performs line solves in Z direction by first factoring
            //     the block-tridiagonal matrix into an upper triangular matrix, 
            //     and { performing back substitution to solve for the unknow
            //     vectors of each line.  
            //     
            //     Make sure we treat elements zero to cell_size in the direction
            //     of the sweep.
            //---------------------------------------------------------------------
            int c, stage, first, last, isize,jsize,ksize,buffer_size; //int[] r_status[MPI_STATUS_SIZE];
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE);
            MPI.Request[] recv_id = new MPI.Request[1];
            MPI.Request[] send_id = new MPI.Request[1];
            double[] out_buffer_x = new double[buffer_size];
            //---------------------------------------------------------------------
            //     in our terminology stage is the number of the cell in the y-direction
            //     i.e. stage = 1 means the start of the line stage=ncells means end
            //---------------------------------------------------------------------
            for(stage = 1; stage<=ncells; stage++){
               c = slice[3,stage];
               isize = cell_size[1,c] - 1;
               jsize = cell_size[2,c] - 1;
               ksize = cell_size[3,c] - 1;
               //---------------------------------------------------------------------
               //     set last-cell flag
               //---------------------------------------------------------------------
               if (stage == ncells) {
                  last = 1;
               } else {
                  last = 0;
               }

               if (stage == 1) {
                  //---------------------------------------------------------------------
                  //     This is the first cell, so solve without receiving data
                  //---------------------------------------------------------------------
                  first = 1;
                  z_solve_cell(first,last,c);  //call z_solve_cell[first,last,c];
               } else {
                  //---------------------------------------------------------------------
                  //     Not the first cell of this line, so receive info from
                  //     processor working on preceeding cell
                  //---------------------------------------------------------------------
                  first = 0;
                  z_receive_solve_info(out_buffer_x,recv_id, c);//call z_receive_solve_info[recv_id,c];
                  //  c---------------------------------------------------------------------
                  //  c     overlap computations and communications
                  //  c---------------------------------------------------------------------
                  //      c      call lhsz[c]
                  //  c---------------------------------------------------------------------
                  //  c     wait for completion
                  //  c---------------------------------------------------------------------
                  send_id[0].Wait();
                  recv_id[0].Wait();
                  //Fortran: call mpi_wait[send_id,r_status,error];
                  //Fortran: call mpi_wait[recv_id,r_status,error];
                  //  c---------------------------------------------------------------------
                  //  c     install C'[kstart+1] and rhs'[kstart+1] to be used in this cell
                  //  c---------------------------------------------------------------------
                  z_unpack_solve_info(out_buffer_x, c);
                  z_solve_cell(first,last,c);  //call z_solve_cell[first,last,c];
               }
               if (last == 0) z_send_solve_info(send_id, c);  //call z_send_solve_info[send_id,c];
            }
            out_buffer_x = null;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE;
            out_buffer_x = new double[buffer_size];
            //---------------------------------------------------------------------
            //     now perform backsubstitution in reverse direction
            //---------------------------------------------------------------------
            for(stage = ncells; stage>= 1; stage--){  //for(stage = ncells, 1, -1;
               c = slice[3,stage];
               first = 0;
               last = 0;
               if (stage == 1) first = 1;
               if (stage == ncells) {
                  last = 1;
                  //---------------------------------------------------------------------
                  //     last cell, so perform back substitute without waiting
                  //---------------------------------------------------------------------
                  z_backsubstitute(first, last,c); //call z_backsubstitute[first, last,c];
               } else {
                  z_receive_backsub_info(out_buffer_x, recv_id, c);            //call z_receive_backsub_info[recv_id,c];
                  send_id[0].Wait();
                  recv_id[0].Wait();
                  // Fortran: call mpi_wait[send_id,r_status,error]; 
                  // Fortran: call mpi_wait[recv_id,r_status,error];
                  z_unpack_backsub_info(out_buffer_x, c);                  
                  z_backsubstitute(first,last,c); //call z_backsubstitute[first,last,c];
               }
               if (first == 0) z_send_backsub_info(send_id, c);  //call z_send_backsub_info[send_id,c];  
            }
        }

        public void z_unpack_solve_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack C'[-1] and rhs'[-1] for
            //     all i and j
            //---------------------------------------------------------------------
            int i,j,m,n,ptr,kstart;
            kstart = 0;
            ptr = 0;
            for(j=0; j<=JMAX-1; j++){
               for(i=0; i<=IMAX-1; i++){
                  for(m=1; m<=BLOCK_SIZE; m++){
                     for(n=0; n<BLOCK_SIZE; n++){  //lhsc[m,n,i,j,kstart-1,c] = out_buffer[ptr+n];
                        lhsc[c,kstart,j+1,i+1,n+1,m] = out_buffer_x[ptr+n];
                     }
                     ptr = ptr+BLOCK_SIZE;
                  }
                  for(n=0; n<BLOCK_SIZE; n++){ //rhs[n,i,j,kstart-1,c] = out_buffer[ptr+n];
                     rhs[c,kstart,j+1,i+1,n+1] = out_buffer_x[ptr+n];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
        }

        public void z_send_solve_info(MPI.Request[] send_id,int c) {
            //---------------------------------------------------------------------
            //     pack up and send C'[kend] and rhs'[kend] for
            //     all i and j
            //---------------------------------------------------------------------
            int i,j,m,n,ksize,ptr,ip,jp;
            int buffer_size;

            ksize = cell_size[3,c]-1;
            ip = cell_coord[1,c] - 1;
            jp = cell_coord[2,c] - 1;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE);
            double[] in_buffer_x = new double[buffer_size];

            //c---------------------------------------------------------------------
            //c     pack up buffer
            //c---------------------------------------------------------------------
            ptr = 0;
            for(j=0; j<=JMAX-1; j++){
               for(i=0; i<=IMAX-1; i++){
                  for(m=1; m<=BLOCK_SIZE; m++){
                     for(n=0; n<BLOCK_SIZE; n++){  //in_buffer[ptr+n] = lhsc[m,n,i,j,ksize,c]
                        in_buffer_x[ptr+n] = lhsc[c,ksize+1,j+1,i+1,n+1,m];
                     }
                     ptr = ptr+BLOCK_SIZE;
                  }
                  for(n=0; n<BLOCK_SIZE; n++){  //in_buffer[ptr+n] = rhs[n,i,j,ksize,c];
                     in_buffer_x[ptr+n] = rhs[c,ksize+1,j+1,i+1,n+1];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }

            //---------------------------------------------------------------------
            //     send buffer 
            //---------------------------------------------------------------------
            //call mpi_isend[in_buffer, buffer_size, dp_type, successor[3], BOTTOM+ip+jp*NCELLS, comm_solve, send_id,error];
            send_id[0] = comm_solve.ImmediateSend<double>(in_buffer_x, successor[3], BOTTOM+ip+jp*ncells);
        }

        public void z_send_backsub_info(MPI.Request[] send_id,int c) {
            //---------------------------------------------------------------------
            //     pack up and send U[jstart] for all i and j
            //---------------------------------------------------------------------
            int i,j,n,ptr,kstart,ip,jp;
            int buffer_size;
            //---------------------------------------------------------------------
            //     Send element 0 to previous processor
            //---------------------------------------------------------------------
            kstart = 0;
            ip = cell_coord[1,c]-1;
            jp = cell_coord[2,c]-1;
            buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE;

            double[] in_buffer_x = new double[buffer_size];

            ptr = 0;
            for(j=0; j<=JMAX-1; j++){
               for(i=0; i<=IMAX-1; i++){
                  for(n=0; n<BLOCK_SIZE; n++){  //in_buffer[ptr+n] = rhs[n,i,j,kstart,c];
                     in_buffer_x[ptr+n] = rhs[c,kstart+1,j+1,i+1,n+1];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
            //call mpi_isend[in_buffer, buffer_size, dp_type, predecessor[3], TOP+ip+jp*NCELLS, comm_solve, send_id,error];

            send_id[0] = comm_solve.ImmediateSend<double>(in_buffer_x, predecessor[3], TOP+ip+jp*ncells);
        }

        public void z_unpack_backsub_info(double[] out_buffer_x, int c) {
            //---------------------------------------------------------------------
            //     unpack U[ksize] for all i and j
            //---------------------------------------------------------------------
            int i,j,n,ptr;
            ptr = 0;
            for(j=0; j<=JMAX-1; j++){
               for(i=0; i<=IMAX-1; i++){
                  for(n=0; n<BLOCK_SIZE; n++){  //backsub_info[n,i,j,c] = out_buffer[ptr+n];
                     backsub_info[c,j,i,n+1] = out_buffer_x[ptr+n];
                  }
                  ptr = ptr+BLOCK_SIZE;
               }
            }
        }

        public void z_receive_backsub_info(double[] out_buffer_x,MPI.Request[] recv_id,int c) {
            //---------------------------------------------------------------------
            //     post mpi receives
            //---------------------------------------------------------------------
            int ip,jp;
            ip = cell_coord[1,c] - 1;
            jp = cell_coord[2,c] - 1;
            //call mpi_irecv[out_buffer, buffer_size, dp_type, successor[3], TOP+ip+jp*NCELLS, comm_solve, recv_id, error];

            recv_id[0] = comm_solve.ImmediateReceive<double>(successor[3], TOP+ip+jp*ncells, out_buffer_x);
        }

        public void z_receive_solve_info(double[] out_buffer_x,MPI.Request[] recv_id,int c) {
            //---------------------------------------------------------------------
            //     post mpi receives 
            //---------------------------------------------------------------------
            int ip,jp;
            ip = cell_coord[1,c] - 1;
            jp = cell_coord[2,c] - 1;
            //call mpi_irecv[out_buffer, buffer_size, dp_type, predecessor[3], BOTTOM+ip+jp*NCELLS, comm_solve, recv_id, error];

            recv_id[0] = comm_solve.ImmediateReceive<double>(predecessor[3], BOTTOM+ip+jp*ncells, out_buffer_x);
        }

        public void z_backsubstitute(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     back solve: if last cell, { generate U[ksize]=rhs[ksize]
            //     } else { assume U[ksize] is loaded in un pack backsub_info
            //     so just use it
            //     after call u[kstart] will be sent to next cell
            //---------------------------------------------------------------------
            int i, k;
            int m,n,j,jsize,isize,ksize,kstart;

            kstart = 0;
            isize = cell_size[1,c]-end[1,c]-1;
            jsize = cell_size[2,c]-end[2,c]-1;
            ksize = cell_size[3,c]-1;
            if (last == 0) {
               for(j=start[2,c]; j<=jsize; j++){
                  for(i=start[1,c]; i<=isize; i++){
                     //---------------------------------------------------------------------
                     //     U[jsize] uses info from previous cell if not last cell
                     //---------------------------------------------------------------------
                     for(m=1; m<=BLOCK_SIZE; m++){
                        for(n=1; n<=BLOCK_SIZE; n++){ //rhs[m,i,j,ksize,c] = rhs[m,i,j,ksize,c] - lhsc[m,n,i,j,ksize,c]*backsub_info[n,i,j,c]
                           rhs[c,ksize+1,j+1,i+1,m] = rhs[c,ksize+1,j+1,i+1,m] - lhsc[c,ksize+1,j+1,i+1,n,m]* backsub_info[c,j,i,n];
                        }
                     }
                  }
               }
            }
            for(k=ksize-1; k>=kstart; k--){  //for(k=ksize-1,kstart,-1;
               for(j=start[2,c]; j<=jsize; j++){
                  for(i=start[1,c]; i<=isize; i++){
                     for(m=1; m<=BLOCK_SIZE; m++){
                        for(n=1; n<=BLOCK_SIZE; n++){  //rhs[m,i,j,k,c] = rhs[m,i,j,k,c] - lhsc[m,n,i,j,k,c]*rhs[n,i,j,k+1,c];
                           rhs[c,k+1,j+1,i+1,m] = rhs[c,k+1,j+1,i+1,m] - lhsc[c,k+1,j+1,i+1,n,m]*rhs[c,k+2,j+1,i+1,n];
                        }
                     }
                  }
               }
            }
        }

        public void z_solve_cell(int first, int last, int c) {
            //---------------------------------------------------------------------
            //     performs guaussian elimination on this cell.
            //     
            //     assumes that unpacking routines for non-first cells 
            //     preload C' and rhs' from previous cell.
            //     
            //     assumed send happens outside this routine, but that
            //     c'[KMAX] and rhs'[KMAX] will be sent to next cell.
            //---------------------------------------------------------------------
            int i,j,k,isize,ksize,jsize,kstart;
            double[,] utmp = new double[KMAX+4,7];   //double utmp[6,-2:KMAX+1];
            kstart = 0;
            isize = cell_size[1,c]-end[1,c]-1;
            jsize = cell_size[2,c]-end[2,c]-1;
            ksize = cell_size[3,c]-1;

            lhsabinit(ref lhsa, ref lhsb, ksize); //call lhsabinit[lhsa, lhsb, ksize];

            for(j=start[2,c]; j<=jsize; j++){
               for(i=start[1,c]; i<=isize; i++){
                  //---------------------------------------------------------------------
                  //     This function computes the left hand side for the three z-factors   
                  //---------------------------------------------------------------------
                  //---------------------------------------------------------------------
                  //     Compute the indices for storing the block-diagonal matrix;
                  //     determine c [labeled f] and s jacobians for cell c
                  //---------------------------------------------------------------------
                  for(k = start[3,c]-1; k<= cell_size[3,c]-end[3,c]; k++){
                     utmp[2+k,1] = 1.0d/u[c,k+2,j+2,i+2,1];  //utmp[1,k] = 1.0d0 / u[1,i,j,k,c]
                     utmp[2+k,2] =      u[c,k+2,j+2,i+2,2];  //utmp[2,k] =         u[2,i,j,k,c]
                     utmp[2+k,3] =      u[c,k+2,j+2,i+2,3];  //utmp[3,k] =         u[3,i,j,k,c]
                     utmp[2+k,4] =      u[c,k+2,j+2,i+2,4];  //utmp[4,k] =         u[4,i,j,k,c]
                     utmp[2+k,5] =      u[c,k+2,j+2,i+2,5];  //utmp[5,k] =         u[5,i,j,k,c]
                     utmp[2+k,6] =       qs[c,k+1,j+1,i+1];  //utmp[6,k] =          qs[i,j,k,c]
                  }
                  for(k = start[3,c]-1; k<= cell_size[3,c]-end[3,c]; k++){

                     tmp1 = utmp[2+k,1];
                     tmp2 = tmp1 * tmp1;
                     tmp3 = tmp1 * tmp2;

                     fjac[2+k,1,1] = 0.0d;   //fjac[1,1,k] = 0.0d+00
                     fjac[2+k,2,1] = 0.0d;
                     fjac[2+k,3,1] = 0.0d;
                     fjac[2+k,4,1] = 1.0d;
                     fjac[2+k,5,1] = 0.0d;

                     fjac[2+k,1,2] = - ( utmp[2+k,2]*utmp[2+k,4] ) * tmp2;
                     fjac[2+k,2,2] = utmp[2+k,4] * tmp1;
                     fjac[2+k,3,2] = 0.0d;
                     fjac[2+k,4,2] = utmp[2+k,2] * tmp1;
                     fjac[2+k,5,2] = 0.0d;

                     fjac[2+k,1,3] = - ( utmp[2+k,3]*utmp[2+k,4] ) * tmp2;
                     fjac[2+k,2,3] = 0.0d;
                     fjac[2+k,3,3] = utmp[2+k,4] * tmp1;
                     fjac[2+k,4,3] = utmp[2+k,3] * tmp1;
                     fjac[2+k,5,3] = 0.0d;

                     fjac[2+k,1,4] = - (utmp[2+k,4]*utmp[2+k,4] * tmp2 ) + c2 * utmp[2+k,6];
                     fjac[2+k,2,4] = - c2 *  utmp[2+k,2] * tmp1;
                     fjac[2+k,3,4] = - c2 *  utmp[2+k,3] * tmp1;
                     fjac[2+k,4,4] = ( 2.0d - c2 ) *  utmp[2+k,4] * tmp1;
                     fjac[2+k,5,4] = c2;

                     fjac[2+k,1,5] = ( c2 * 2.0d * utmp[2+k,6] - c1 * ( utmp[2+k,5] * tmp1 ) ) * ( utmp[2+k,4] * tmp1 );
                     fjac[2+k,2,5] = - c2 * ( utmp[2+k,2]*utmp[2+k,4] ) * tmp2;
                     fjac[2+k,3,5] = - c2 * ( utmp[2+k,3]*utmp[2+k,4] ) * tmp2;
                     fjac[2+k,4,5] = c1 * ( utmp[2+k,5] * tmp1 ) - c2 * ( utmp[2+k,6] + utmp[2+k,4]*utmp[2+k,4] * tmp2 );
                     fjac[2+k,5,5] = c1 * utmp[2+k,4] * tmp1;

                     njac[2+k,1,1] = 0.0d;
                     njac[2+k,2,1] = 0.0d;
                     njac[2+k,3,1] = 0.0d;
                     njac[2+k,4,1] = 0.0d;
                     njac[2+k,5,1] = 0.0d;

                     njac[2+k,1,2] = - c3c4 * tmp2 * utmp[2+k,2];
                     njac[2+k,2,2] =   c3c4 * tmp1;
                     njac[2+k,3,2] =   0.0d;
                     njac[2+k,4,2] =   0.0d;
                     njac[2+k,5,2] =   0.0d;

                     njac[2+k,1,3] = - c3c4 * tmp2 * utmp[2+k,3];
                     njac[2+k,2,3] =   0.0d;
                     njac[2+k,3,3] =   c3c4 * tmp1;
                     njac[2+k,4,3] =   0.0d;
                     njac[2+k,5,3] =   0.0d;

                     njac[2+k,1,4] = - con43 * c3c4 * tmp2 * utmp[2+k,4];
                     njac[2+k,2,4] =   0.0d;
                     njac[2+k,3,4] =   0.0d;
                     njac[2+k,4,4] =   con43 * c3 * c4 * tmp1;
                     njac[2+k,5,4] =   0.0d;

                     njac[2+k,1,5] = - (c3c4 - c1345)*tmp3*(pow2(utmp[2+k,2]))-
                                     (c3c4-c1345)*tmp3*(pow2(utmp[2+k,3]))-(con43*c3c4-c1345)*tmp3*(pow2(utmp[2+k,4]))-c1345*tmp2*utmp[2+k,5];

                     njac[2+k,2,5] = (  c3c4 - c1345 ) * tmp2 * utmp[2+k,2];
                     njac[2+k,3,5] = (  c3c4 - c1345 ) * tmp2 * utmp[2+k,3];
                     njac[2+k,4,5] = ( con43 * c3c4 - c1345 ) * tmp2 * utmp[2+k,4];
                     njac[2+k,5,5] = ( c1345 )* tmp1;
                  }
                  //---------------------------------------------------------------------
                  //     now joacobians set, so form left hand side in z direction
                  //---------------------------------------------------------------------
                  for(k = start[3,c]; k<= ksize-end[3,c]; k++){
                     tmp1 = dt * tz1;
                     tmp2 = dt * tz2;  
                     lhsa[k+1,1,1] = - tmp2 * fjac[k+1,1,1] - tmp1 * njac[k+1,1,1] - tmp1*dz1;//lhsa[1,1,k]=-tmp2*fjac[1,1,k-1]-tmp1*njac[1,1,k-1]
                     lhsa[k+1,2,1] = - tmp2 * fjac[k+1,2,1] - tmp1 * njac[k+1,2,1];    //fjac[k-1+2,2,1] njac[k-1+2,2,1]
                     lhsa[k+1,3,1] = - tmp2 * fjac[k+1,3,1] - tmp1 * njac[k+1,3,1];
                     lhsa[k+1,4,1] = - tmp2 * fjac[k+1,4,1] - tmp1 * njac[k+1,4,1];
                     lhsa[k+1,5,1] = - tmp2 * fjac[k+1,5,1] - tmp1 * njac[k+1,5,1];

                     lhsa[k+1,1,2] = - tmp2 * fjac[k+1,1,2] - tmp1 * njac[k+1,1,2];
                     lhsa[k+1,2,2] = - tmp2 * fjac[k+1,2,2] - tmp1 * njac[k+1,2,2] - tmp1 * dz2;
                     lhsa[k+1,3,2] = - tmp2 * fjac[k+1,3,2] - tmp1 * njac[k+1,3,2];
                     lhsa[k+1,4,2] = - tmp2 * fjac[k+1,4,2] - tmp1 * njac[k+1,4,2];
                     lhsa[k+1,5,2] = - tmp2 * fjac[k+1,5,2] - tmp1 * njac[k+1,5,2];

                     lhsa[k+1,1,3] = - tmp2 * fjac[k+1,1,3] - tmp1 * njac[k+1,1,3];
                     lhsa[k+1,2,3] = - tmp2 * fjac[k+1,2,3] - tmp1 * njac[k+1,2,3];
                     lhsa[k+1,3,3] = - tmp2 * fjac[k+1,3,3] - tmp1 * njac[k+1,3,3] - tmp1 * dz3;
                     lhsa[k+1,4,3] = - tmp2 * fjac[k+1,4,3] - tmp1 * njac[k+1,4,3];
                     lhsa[k+1,5,3] = - tmp2 * fjac[k+1,5,3] - tmp1 * njac[k+1,5,3];

                     lhsa[k+1,1,4] = - tmp2 * fjac[k+1,1,4] - tmp1 * njac[k+1,1,4];
                     lhsa[k+1,2,4] = - tmp2 * fjac[k+1,2,4] - tmp1 * njac[k+1,2,4];
                     lhsa[k+1,3,4] = - tmp2 * fjac[k+1,3,4] - tmp1 * njac[k+1,3,4];
                     lhsa[k+1,4,4] = - tmp2 * fjac[k+1,4,4] - tmp1 * njac[k+1,4,4] - tmp1 * dz4;
                     lhsa[k+1,5,4] = - tmp2 * fjac[k+1,5,4] - tmp1 * njac[k+1,5,4];

                     lhsa[k+1,1,5] = - tmp2 * fjac[k+1,1,5] - tmp1 * njac[k+1,1,5];
                     lhsa[k+1,2,5] = - tmp2 * fjac[k+1,2,5] - tmp1 * njac[k+1,2,5];
                     lhsa[k+1,3,5] = - tmp2 * fjac[k+1,3,5] - tmp1 * njac[k+1,3,5];
                     lhsa[k+1,4,5] = - tmp2 * fjac[k+1,4,5] - tmp1 * njac[k+1,4,5];
                     lhsa[k+1,5,5] = - tmp2 * fjac[k+1,5,5] - tmp1 * njac[k+1,5,5] - tmp1 * dz5;

                     lhsb[k+1,1,1] = 1.0d + tmp1 * 2.0d * njac[2+k,1,1] + tmp1 * 2.0d * dz1; //lhsb[1,1,k]=1.0d+tmp1*2.0d+00*njac[1,1,k]
                     lhsb[k+1,2,1] = tmp1 * 2.0d * njac[2+k,2,1];
                     lhsb[k+1,3,1] = tmp1 * 2.0d * njac[2+k,3,1];
                     lhsb[k+1,4,1] = tmp1 * 2.0d * njac[2+k,4,1];
                     lhsb[k+1,5,1] = tmp1 * 2.0d * njac[2+k,5,1];

                     lhsb[k+1,1,2] = tmp1 * 2.0d * njac[2+k,1,2];
                     lhsb[k+1,2,2] = 1.0d + tmp1 * 2.0d * njac[2+k,2,2] + tmp1 * 2.0d * dz2;
                     lhsb[k+1,3,2] = tmp1 * 2.0d * njac[2+k,3,2];
                     lhsb[k+1,4,2] = tmp1 * 2.0d * njac[2+k,4,2];
                     lhsb[k+1,5,2] = tmp1 * 2.0d * njac[2+k,5,2];

                     lhsb[k+1,1,3] = tmp1 * 2.0d * njac[2+k,1,3];
                     lhsb[k+1,2,3] = tmp1 * 2.0d * njac[2+k,2,3];
                     lhsb[k+1,3,3] = 1.0d + tmp1 * 2.0d * njac[2+k,3,3] + tmp1 * 2.0d * dz3;
                     lhsb[k+1,4,3] = tmp1 * 2.0d * njac[2+k,4,3];
                     lhsb[k+1,5,3] = tmp1 * 2.0d * njac[2+k,5,3];

                     lhsb[k+1,1,4] = tmp1 * 2.0d * njac[2+k,1,4];
                     lhsb[k+1,2,4] = tmp1 * 2.0d * njac[2+k,2,4];
                     lhsb[k+1,3,4] = tmp1 * 2.0d * njac[2+k,3,4];
                     lhsb[k+1,4,4] = 1.0d + tmp1 * 2.0d * njac[2+k,4,4] + tmp1 * 2.0d * dz4;
                     lhsb[k+1,5,4] = tmp1 * 2.0d * njac[2+k,5,4];

                     lhsb[k+1,1,5] = tmp1 * 2.0d * njac[2+k,1,5];
                     lhsb[k+1,2,5] = tmp1 * 2.0d * njac[2+k,2,5];
                     lhsb[k+1,3,5] = tmp1 * 2.0d * njac[2+k,3,5];
                     lhsb[k+1,4,5] = tmp1 * 2.0d * njac[2+k,4,5];
                     lhsb[k+1,5,5] = 1.0d + tmp1 * 2.0d * njac[2+k,5,5] + tmp1 * 2.0d * dz5;
                     
                     lhsc[c,k+1,j+1,i+1,1,1] =  tmp2 * fjac[k+3,1,1] - tmp1 * njac[k+3,1,1]-tmp1*dz1;//lhsc[1,1,i,j,k,c] fjac[1,1,k+1] njac[1,1,k+1]
                     lhsc[c,k+1,j+1,i+1,2,1] =  tmp2 * fjac[k+3,2,1] - tmp1 * njac[k+3,2,1];//fjac[k+1+2,2,1] njac[k+1+2,2,1]
                     lhsc[c,k+1,j+1,i+1,3,1] =  tmp2 * fjac[k+3,3,1] - tmp1 * njac[k+3,3,1];
                     lhsc[c,k+1,j+1,i+1,4,1] =  tmp2 * fjac[k+3,4,1] - tmp1 * njac[k+3,4,1];
                     lhsc[c,k+1,j+1,i+1,5,1] =  tmp2 * fjac[k+3,5,1] - tmp1 * njac[k+3,5,1];

                     lhsc[c,k+1,j+1,i+1,1,2] =  tmp2 * fjac[k+3,1,2] - tmp1 * njac[k+3,1,2];
                     lhsc[c,k+1,j+1,i+1,2,2] =  tmp2 * fjac[k+3,2,2] - tmp1 * njac[k+3,2,2] - tmp1 * dz2;
                     lhsc[c,k+1,j+1,i+1,3,2] =  tmp2 * fjac[k+3,3,2] - tmp1 * njac[k+3,3,2];
                     lhsc[c,k+1,j+1,i+1,4,2] =  tmp2 * fjac[k+3,4,2] - tmp1 * njac[k+3,4,2];
                     lhsc[c,k+1,j+1,i+1,5,2] =  tmp2 * fjac[k+3,5,2] - tmp1 * njac[k+3,5,2];

                     lhsc[c,k+1,j+1,i+1,1,3] =  tmp2 * fjac[k+3,1,3] - tmp1 * njac[k+3,1,3];
                     lhsc[c,k+1,j+1,i+1,2,3] =  tmp2 * fjac[k+3,2,3] - tmp1 * njac[k+3,2,3];
                     lhsc[c,k+1,j+1,i+1,3,3] =  tmp2 * fjac[k+3,3,3] - tmp1 * njac[k+3,3,3] - tmp1 * dz3;
                     lhsc[c,k+1,j+1,i+1,4,3] =  tmp2 * fjac[k+3,4,3] - tmp1 * njac[k+3,4,3];
                     lhsc[c,k+1,j+1,i+1,5,3] =  tmp2 * fjac[k+3,5,3] - tmp1 * njac[k+3,5,3];

                     lhsc[c,k+1,j+1,i+1,1,4] =  tmp2 * fjac[k+3,1,4] - tmp1 * njac[k+3,1,4];
                     lhsc[c,k+1,j+1,i+1,2,4] =  tmp2 * fjac[k+3,2,4] - tmp1 * njac[k+3,2,4];
                     lhsc[c,k+1,j+1,i+1,3,4] =  tmp2 * fjac[k+3,3,4] - tmp1 * njac[k+3,3,4];
                     lhsc[c,k+1,j+1,i+1,4,4] =  tmp2 * fjac[k+3,4,4] - tmp1 * njac[k+3,4,4] - tmp1 * dz4;
                     lhsc[c,k+1,j+1,i+1,5,4] =  tmp2 * fjac[k+3,5,4] - tmp1 * njac[k+3,5,4];

                     lhsc[c,k+1,j+1,i+1,1,5] =  tmp2 * fjac[k+3,1,5] - tmp1 * njac[k+3,1,5];
                     lhsc[c,k+1,j+1,i+1,2,5] =  tmp2 * fjac[k+3,2,5] - tmp1 * njac[k+3,2,5];
                     lhsc[c,k+1,j+1,i+1,3,5] =  tmp2 * fjac[k+3,3,5] - tmp1 * njac[k+3,3,5];
                     lhsc[c,k+1,j+1,i+1,4,5] =  tmp2 * fjac[k+3,4,5] - tmp1 * njac[k+3,4,5];
                     lhsc[c,k+1,j+1,i+1,5,5] =  tmp2 * fjac[k+3,5,5] - tmp1 * njac[k+3,5,5] - tmp1 * dz5;
                  }
                  //---------------------------------------------------------------------
                  //     outer most for(loops - sweeping in i direction
                  //---------------------------------------------------------------------
                  if (first == 1) { 
                     //---------------------------------------------------------------------
                     //     multiply c[i,j,kstart] by b_inverse and copy back to c
                     //     multiply rhs[kstart] by b_inverse[kstart] and copy to rhs
                     //---------------------------------------------------------------------
                     //Fortran: call binvcrhs[ lhsb[1,1,kstart], lhsc[1,1,i,j,kstart,c], rhs[1,i,j,kstart,c] ];
                     //C#:           binvcrhs( lhsb[kstart+1,1,1], lhsc[c,kstart+1,j+1,i+1,1,1], rhs[c,kstart+1,j+1,i+1,1] );
                     binvcrhs(ref lhsb,ref lhsc,ref rhs,           (kstart+1),     c,(kstart+1),(j+1),(i+1),     c,(kstart+1),(j+1),(i+1));
                  }
                  //c---------------------------------------------------------------------
                  //c     begin inner most for(loop
                  //c     for(all the elements of the cell unless last 
                  //c---------------------------------------------------------------------
                  for(k=kstart+first; k<=ksize-last; k++){
                     //c---------------------------------------------------------------------
                     //c     subtract A*lhs_vector[k-1] from lhs_vector[k]
                     //c     
                     //c     rhs[k] = rhs[k] - A*rhs[k-1]
                     //c---------------------------------------------------------------------
                     //Fortran: call matvec_sub[lhsa[1,1,k], rhs[1,i,j,k-1,c],rhs[1,i,j,k,c]];
                     //C#: matvec_sub(lhsa[k+1,1,1], rhs[c,k-1+1,j+1,i+1,1],rhs[c,k+1,j+1,i+1,1]);
                     matvec_sub(ref lhsa,ref rhs,ref rhs,      (k+1),     c,(k),(j+1),(i+1),      c,(k+1),(j+1),(i+1));
                     //c---------------------------------------------------------------------
                     //c     B[k] = B[k] - C[k-1]*A[k]
                     //c     call matmul_sub[aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k,c]
                     //c---------------------------------------------------------------------
                     //Fortran: call matmul_sub[lhsa[1,1,k], lhsc[1,1,i,j,k-1,c], lhsb[1,1,k]];
                     //C#:           matmul_sub(lhsa[k+1,1,1], lhsc[c,k-1+1,j+1,i+1,1,1], lhsb[k+1,1,1]);
                     matmul_sub(ref lhsa,ref lhsc,ref lhsb,             (k+1),     c,k,(j+1),(i+1),    (k+1));
                     //c---------------------------------------------------------------------
                     //c     multiply c[i,j,k] by b_inverse and copy back to c
                     //c     multiply rhs[i,j,1] by b_inverse[i,j,1] and copy to rhs
                     //c---------------------------------------------------------------------
                     //Fortran: call binvcrhs[ lhsb[1,1,k], lhsc[1,1,i,j,k,c], rhs[1,i,j,k,c] ];
                     //C#:           binvcrhs( lhsb[k+1,1,1],    lhsc[c,k+1,j+1,i+1,1,1],    rhs[c,k+1,j+1,i+1,1] );
                     binvcrhs(ref lhsb,ref lhsc,ref rhs,      (k+1),      c,(k+1),(j+1),(i+1),     c,(k+1),(j+1),(i+1));
                  }
                  //c---------------------------------------------------------------------
                  //c     Now finish up special cases for last cell
                  //c---------------------------------------------------------------------
                  if (last == 1) {
                     //c---------------------------------------------------------------------
                     //c     rhs[ksize] = rhs[ksize] - A*rhs[ksize-1]
                     //c---------------------------------------------------------------------
                     //Fortran: call matvec_sub[lhsa[1,1,ksize], rhs[1,i,j,ksize-1,c],rhs[1,i,j,ksize,c]];
                     //C#:           matvec_sub(lhsa[ksize+1,1,1], rhs[c,ksize-1+1,j+1,i+1,1],rhs[c,ksize+1,j+1,i+1,1]);
                     matvec_sub(ref lhsa,ref rhs,ref rhs,          (ksize+1),      c,(ksize),(j+1),(i+1),     c,(ksize+1),(j+1),(i+1));
                     //c---------------------------------------------------------------------
                     //c     B[ksize] = B[ksize] - C[ksize-1]*A[ksize]
                     //c     call matmul_sub[aa,i,j,ksize,c,
                     //c     $              cc,i,j,ksize-1,c,bb,i,j,ksize,c]
                     //c---------------------------------------------------------------------
                     //Fortran: call matmul_sub[lhsa[1,1,ksize],lhsc[1,1,i,j,ksize-1,c],lhsb[1,1,ksize]];
                     //C#:           matmul_sub[lhsa[ksize+1,1,1],lhsc[c,ksize-1+1,j+1,i+1,1,1],lhsb[ksize+1,1,1]];
                     matmul_sub(ref lhsa,ref lhsc,ref lhsb,         (ksize+1),     c,(ksize),(j+1),(i+1),      (ksize+1));
                     //c---------------------------------------------------------------------
                     //c     multiply rhs[ksize] by b_inverse[ksize] and copy to rhs
                     //c---------------------------------------------------------------------
                     //Fortran: call binvrhs[ lhsb[1,1,ksize], rhs[1,i,j,ksize,c] ];
                     //C#:           binvrhs(lhsb[ksize+1,1,1], rhs[c,ksize+1,j+1,i+1,1] );
                     binvrhs(ref lhsb,ref rhs,        (ksize+1),     c,(ksize+1),(j+1),(i+1));
                  }
               }
            }
        }
        // end z_solve.f

        // add.f
        public void add() {
            //---------------------------------------------------------------------
            //     addition of update to the vector u
            //---------------------------------------------------------------------
            int  c, i, j, k, m;
            for(c = 1; c<= ncells; c++){
               for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                  for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                     for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){  //u[m,i,j,k,c] = u[m,i,j,k,c] + rhs[m,i,j,k,c];
                           u[c,k+2,j+2,i+2,m] = u[c,k+2,j+2,i+2,m] + rhs[c,k+1,j+1,i+1,m];
                        }
                     }
                  }
               }
            }
        }
        // end add.f

        // verify.f
        public void verify(int no_time_steps, char clss, ref bool verified) {
            //---------------------------------------------------------------------
            //  verification routine                         
            //---------------------------------------------------------------------
            double[] xcrref = new double[6];
            double[] xceref = new double[6];
            double[] xcrdif = new double[6];
            double[] xcedif = new double[6];
            double[] xce    = new double[6];
            double[] xcr    = new double[6];
            double dtref=0, epsilon;
            int m;
            //---------------------------------------------------------------------
            //   tolerance level
            //---------------------------------------------------------------------
            epsilon = 1.0 * Math.Pow(.1, 8);//1.0d/100000000; //1.0d-08;
            verified = true;
            //---------------------------------------------------------------------
            //   compute the error norm and the residual norm, and exit if not printing
            //---------------------------------------------------------------------
            if (iotype != 0) {
                //call accumulate_norms[xce]; Function null in Fortran code
            } else {
                error_norm(xce);//call error_norm[xce];
            }
            copy_faces();
            rhs_norm(xcr);
            for(m = 1; m<= 5; m++){
                xcr[m] = xcr[m] / dt;
            }
            if (node != 0) return;
            clss = 'U';
            for(m = 1; m<=5; m++){
                xcrref[m] = 1.0;
                xceref[m] = 1.0;
            }
            //---------------------------------------------------------------------
            //    reference data for 12X12X12 grids after 60 time steps, with DT = 1.0d-02
            //---------------------------------------------------------------------
            if ((grid_points[1]  == 12) && (grid_points[2]  == 12) && (grid_points[3]  == 12) && (no_time_steps   == 60))  {
                clss = 'S';
                dtref = 0.01; // 1.0d / 100;//1.0d-2;
                //c---------------------------------------------------------------------
                //c  Reference values of RMS-norms of residual.
                //c---------------------------------------------------------------------
                xcrref[1] = 1.7034283709541311E-01d;//xcrref[1] = 1.7034283709541311d-01;
                xcrref[2] = 1.2975252070034097E-02d;
                xcrref[3] = 3.2527926989486055E-02d;
                xcrref[4] = 2.6436421275166801E-02d;
                xcrref[5] = 1.9211784131744430E-01d;
                //c---------------------------------------------------------------------
                //c  Reference values of RMS-norms of solution error.
                //c---------------------------------------------------------------------
                if (iotype == 0) {
                    xceref[1] = 4.9976913345811579E-04d;
                    xceref[2] = 4.5195666782961927E-05d;
                    xceref[3] = 7.3973765172921357E-05d;
                    xceref[4] = 7.3821238632439731E-05d;
                    xceref[5] = 8.9269630987491446E-04d;
                } else {
                    xceref[1] = 0.1149036328945E+02d;
                    xceref[2] = 0.9156788904727E+00d;
                    xceref[3] = 0.2857899428614E+01d;
                    xceref[4] = 0.2598273346734E+01d;
                    xceref[5] = 0.2652795397547E+02d;
                }
                //   c---------------------------------------------------------------------
                //   c    reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
                //   c---------------------------------------------------------------------
            } else if ((grid_points[1] == 24) && (grid_points[2] == 24) && (grid_points[3] == 24) && (no_time_steps == 200)) {
                clss = 'W';
                dtref = 0.0008;//0.8d/1000; //0.8d-3
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[1] = 0.1125590409344E+03;
                xcrref[2] = 0.1180007595731E+02;
                xcrref[3] = 0.2710329767846E+02;
                xcrref[4] = 0.2469174937669E+02;
                xcrref[5] = 0.2638427874317E+03;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                if (iotype == 0) {
                    xceref[1] = 0.4419655736008E+01;
                    xceref[2] = 0.4638531260002;
                    xceref[3] = 0.1011551749967E+01;
                    xceref[4] = 0.9235878729944;
                    xceref[5] = 0.1018045837718E+02;
                } else {
                    xceref[1] = 0.6729594398612E+02;
                    xceref[2] = 0.5264523081690E+01;
                    xceref[3] = 0.1677107142637E+02;
                    xceref[4] = 0.1508721463436E+02;
                    xceref[5] = 0.1477018363393E+03;
                }
                //   c---------------------------------------------------------------------
                //   c    reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
                //   c---------------------------------------------------------------------
            } else if ((grid_points[1] == 64) && (grid_points[2] == 64) && (grid_points[3] == 64) && (no_time_steps == 200)) {
                clss = 'A';
                dtref = 0.0008;//0.8d/1000; //0.8d-3
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[1] = 1.0806346714637264E+02d;
                xcrref[2] = 1.1319730901220813E+01d;
                xcrref[3] = 2.5974354511582465E+01d;
                xcrref[4] = 2.3665622544678910E+01d;
                xcrref[5] = 2.5278963211748344E+02d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                if (iotype == 0) {
                    xceref[1] = 4.2348416040525025E+00d;
                    xceref[2] = 4.4390282496995698E-01d;
                    xceref[3] = 9.6692480136345650E-01d;
                    xceref[4] = 8.8302063039765474E-01d;
                    xceref[5] = 9.7379901770829278E+00d;
                } else {
                    xceref[1] = 0.6482218724961E+02d;
                    xceref[2] = 0.5066461714527E+01d;
                    xceref[3] = 0.1613931961359E+02d;
                    xceref[4] = 0.1452010201481E+02d;
                    xceref[5] = 0.1420099377681E+03d;
                }
                //   c---------------------------------------------------------------------
                //   c    reference data for 102X102X102 grids after 200 time steps,
                //   c    with DT = 3.0d-04
                //   c---------------------------------------------------------------------
            } else if ( (grid_points[1] == 102) && (grid_points[2] == 102) && (grid_points[3] == 102) && (no_time_steps == 200) ) {
                clss = 'B';
                dtref = .0003;//3.0d/10000;  //3.0d-4
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[1] = 1.4233597229287254E+03d;
                xcrref[2] = 9.9330522590150238E+01d;
                xcrref[3] = 3.5646025644535285E+02d;
                xcrref[4] = 3.2485447959084092E+02d;
                xcrref[5] = 3.2707541254659363E+03d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                if (iotype == 0) {
                    xceref[1] = 5.2969847140936856E+01d;
                    xceref[2] = 4.4632896115670668E+00d;
                    xceref[3] = 1.3122573342210174E+01d;
                    xceref[4] = 1.2006925323559144E+01d;
                    xceref[5] = 1.2459576151035986E+02d;
                } else {
                    xceref[1] = 0.1477545106464E+03d;
                    xceref[2] = 0.1108895555053E+02d;
                    xceref[3] = 0.3698065590331E+02d;
                    xceref[4] = 0.3310505581440E+02d;
                    xceref[5] = 0.3157928282563E+03d;
                }
                //   c---------------------------------------------------------------------
                //   c    reference data for 162X162X162 grids after 200 time steps,
                //   c    with DT = 1.0d-04
                //   c---------------------------------------------------------------------
            } else if ( (grid_points[1] == 162) && (grid_points[2] == 162) && (grid_points[3] == 162) && (no_time_steps == 200) ) {
                clss = 'C';
                dtref = .0001;//1.0d/10000; //1.0d-4
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[1] = 0.62398116551764615E+04d;
                xcrref[2] = 0.50793239190423964E+03d;
                xcrref[3] = 0.15423530093013596E+04d;
                xcrref[4] = 0.13302387929291190E+04d;
                xcrref[5] = 0.11604087428436455E+05d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                if (iotype == 0) {
                    xceref[1] = 0.16462008369091265E+03d;
                    xceref[2] = 0.11497107903824313E+02d;
                    xceref[3] = 0.41207446207461508E+02d;
                    xceref[4] = 0.37087651059694167E+02d;
                    xceref[5] = 0.36211053051841265E+03d;
                } else {
                    xceref[1] = 0.2597156483475E+03d;
                    xceref[2] = 0.1985384289495E+02d;
                    xceref[3] = 0.6517950485788E+02d;
                    xceref[4] = 0.5757235541520E+02d;
                    xceref[5] = 0.5215668188726E+03d;
                }
                //   c---------------------------------------------------------------------
                //   c    reference data for 408x408x408 grids after 250 time steps,
                //   c    with DT = 0.2d-04
                //   c---------------------------------------------------------------------
            } else if ( (grid_points[1] == 408) && (grid_points[2] == 408) && (grid_points[3] == 408) && (no_time_steps == 250) ) {
                clss = 'D';
                dtref = 0.2/10000; //0.2d-4
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[1] = 0.2533188551738E+05d;
                xcrref[2] = 0.2346393716980E+04d;
                xcrref[3] = 0.6294554366904E+04d;
                xcrref[4] = 0.5352565376030E+04d;
                xcrref[5] = 0.3905864038618E+05d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                if (iotype == 0) {
                    xceref[1] = 0.3100009377557E+03d;
                    xceref[2] = 0.2424086324913E+02d;
                    xceref[3] = 0.7782212022645E+02d;
                    xceref[4] = 0.6835623860116E+02d;
                    xceref[5] = 0.6065737200368E+03d;
                } else {
                    xceref[1] = 0.3813781566713E+03d;
                    xceref[2] = 0.3160872966198E+02d;
                    xceref[3] = 0.9593576357290E+02d;
                    xceref[4] = 0.8363391989815E+02d;
                    xceref[5] = 0.7063466087423E+03d;
                }
                //   c---------------------------------------------------------------------
                //   c    reference data for 1020x1020x1020 grids after 250 time steps,
                //   c    with DT = 0.4d-05
                //   c---------------------------------------------------------------------
            } else if ( (grid_points[1] == 1020) && (grid_points[2] == 1020) && (grid_points[3] == 1020) && (no_time_steps == 250) ) {
                clss = 'E';
                dtref = 0.4/100000; //0.4d-5;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of residual.
                //       c---------------------------------------------------------------------
                xcrref[1] = 0.9795372484517E+05d;
                xcrref[2] = 0.9739814511521E+04d;
                xcrref[3] = 0.2467606342965E+05d;
                xcrref[4] = 0.2092419572860E+05d;
                xcrref[5] = 0.1392138856939E+06d;
                //       c---------------------------------------------------------------------
                //       c  Reference values of RMS-norms of solution error.
                //       c---------------------------------------------------------------------
                if (iotype == 0) {
                    xceref[1] = 0.4327562208414E+03d;
                    xceref[2] = 0.3699051964887E+02d;
                    xceref[3] = 0.1089845040954E+03d;
                    xceref[4] = 0.9462517622043E+02d;
                    xceref[5] = 0.7765512765309E+03d;
                } else {
                    xceref[1] = 0.4729898413058E+03d;
                    xceref[2] = 0.4145899331704E+02d;
                    xceref[3] = 0.1192850917138E+03d;
                    xceref[4] = 0.1032746026932E+03d;
                    xceref[5] = 0.8270322177634E+03d;
                }
            } else {
                verified = false;
            }
            //---------------------------------------------------------------------
            //    verification test for residuals if gridsize is one of 
            //    the defined grid sizes above [class != 'U']
            //---------------------------------------------------------------------
            //    Compute the difference of solution values and the known reference 
            //    values.
            //---------------------------------------------------------------------
            for(m = 1; m<= 5; m++){
                xcrdif[m] = Math.Abs((xcr[m]-xcrref[m])/xcrref[m]); //dabs((xcr[m]-xcrref[m])/xcrref[m]);
                xcedif[m] = Math.Abs((xce[m]-xceref[m])/xceref[m]); //dabs((xce[m]-xceref[m])/xceref[m]);
            }
            //---------------------------------------------------------------------
            //    Output the comparison of computed results to known cases.
            //---------------------------------------------------------------------
            if (clss != 'U') {
                Console.WriteLine(" Verification being performed for class "+clss);
                Console.WriteLine(" accuracy setting for epsilon = "+epsilon);
                verified = (Math.Abs(dt-dtref) <= epsilon); 
                if (!verified) {
                    verified = false;
                    clss = 'U';
                    Console.WriteLine(" DT does not match the reference value of "+dtref);
                }
            } else {
                Console.WriteLine(" Unknown class");
            }
            if (clss != 'U') {
                Console.WriteLine(" Comparison of RMS-norms of residual"); 
            } else {
                Console.WriteLine(" RMS-norms of residual");
            }
            for(m = 1; m<= 5; m++){
                if (clss == 'U') {
                    Console.WriteLine("           "+m+" "+xcr[m]);
                } else if (xcrdif[m] <= epsilon) {
                    Console.WriteLine("     "+m+" "+xcr[m]+" "+xcrref[m]+" "+xcrdif[m]+" ");
                } else { 
                    verified = false;
                    Console.WriteLine(" FAILURE: "+m+" "+xcr[m]+" "+xcrref[m]+" "+xcrdif[m]+" ");
                }
            }
            if (clss != 'U') {
                Console.WriteLine(" Comparison of RMS-norms of solution error");
            } else {
                Console.WriteLine(" RMS-norms of solution error");
            }
            for(m = 1; m<= 5; m++){
                if (clss == 'U') {
                    Console.WriteLine("         " + m + " " + xce[m] + " ");
                } else if (xcedif[m] <= epsilon) {
                    Console.WriteLine("      "+m+" "+xce[m]+" "+xceref[m]+" "+xcedif[m]+" "); 
                } else {
                    verified = false;
                    Console.WriteLine(" FAILURE: "+m+" "+xce[m]+" "+xceref[m]+" "+xcedif[m]+" "); 
                }
            }
            if (clss == 'U') {
                Console.WriteLine(" No reference values provided");
                Console.WriteLine(" No verification performed");
            } else if (verified) {
                Console.WriteLine(" Verification Successful");
            } else {
                Console.WriteLine(" Verification failed");
            }
        }
              //error.f
        public void error_norm(double[] rms) {
            //---------------------------------------------------------------------
            //     this function computes the norm of the difference between the
            //     computed solution and the exact solution
            //---------------------------------------------------------------------
            int c, i, j, k, m, ii, jj, kk, d;//double xi, eta, zeta, u_exact[5], rms[5], rms_work[5], add;
            double xi, eta, zeta, add;
            double[] u_exact   = new double[6];
            double[] rms_work  = new double[6];
            for(m = 1; m<= 5; m++){
               rms_work[m] = 0.0; //0.0d;
            }
            for(c = 1; c<= ncells; c++){
               kk = 0;
               for(k = cell_low[3,c]; k<= cell_high[3,c]; k++){
                  zeta = k*dnzm1; 
                  jj = 0;
                  for(j = cell_low[2,c]; j<= cell_high[2,c]; j++){
                     eta = j*dnym1;
                     ii = 0;
                     for(i = cell_low[1,c]; i<= cell_high[1,c]; i++){
                        xi = i*dnxm1;
                        exact_solution1(xi, eta, zeta, ref u_exact);
                        for(m = 1; m<= 5; m++){
                           add = u[c,kk+2,jj+2,ii+2,m]-u_exact[m];
                           rms_work[m] += add*add;
                        }
                        ii = ii + 1;
                     }
                     jj = jj + 1;
                  }
                  kk = kk + 1;
               }
            }
            comm_setup.Allreduce<double>(rms_work, MPI.Operation<double>.Add, ref rms);//call mpi_allreduce[rms_work, rms, 5, dp_type, MPI_SUM, comm_setup, error];
            for(m = 1; m<= 5; m++){
               for(d = 1; d<= 3; d++){
                  rms[m] /= (grid_points[d]-2); 
               }
               rms[m] = Math.Sqrt(rms[m]);               
            }
        }

        public void rhs_norm(double[] rms) {
            int c, i, j, k, d, m;
            double[] rms_work = new double[6];
            double add;

            for(m = 1; m<= 5; m++){
                rms_work[m] = 0.0; //0.0d;
            }
            for(c = 1; c<= ncells; c++){
               for(k = start[3,c]; k<= cell_size[3,c]-end[3,c]-1; k++){
                  for(j = start[2,c]; j<= cell_size[2,c]-end[2,c]-1; j++){
                     for(i = start[1,c]; i<= cell_size[1,c]-end[1,c]-1; i++){
                        for(m = 1; m<= 5; m++){
                           add = rhs[c,k+1,j+1,i+1,m]; 
                           rms_work[m] += add*add; 
                        }
                     }
                  } 
               } 
            }
            comm_setup.Allreduce<double>(rms_work, MPI.Operation<double>.Add, ref rms); //call mpi_allreduce[rms_work, rms, 5, dp_type, MPI_SUM, comm_setup, error];
            for(m = 1; m<= 5; m++){
               for(d = 1; d<= 3; d++){
                  rms[m] = rms[m] / (grid_points[d]-2);
               } 
               rms[m] = Math.Sqrt(rms[m]);             
            } 

        }
              // end error.f
        // end verify.f
    }
}


