/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			C# 	V E R S I O N			  !
!									  !
!                                  F T                                    !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is a serial/multithreaded version of the              !
!    NPB3_0_JAV FT code.                                                  !
!									  !
!    Permission to use, copy, distribute and modify this software	  !
!    for any purpose with or without fee is hereby granted.  We 	  !
!    request, however, that all derived work reference the NAS  	  !
!    Parallel Benchmarks 3.0. This software is provided "as is" 	  !
!    without express or implied warranty.				  !
!									  !
!    Information on NPB 3.0, including the Technical Report NAS-02-008	  !
!    "Implementation of the NAS Parallel Benchmarks in Java",		  !
!    original specifications, source code, results and information	  !
!    on how to submit new results, is available at:			  !
!									  !
!	    http://www.nas.nasa.gov/Software/NPB/			  !
!									  !
!    Send comments or suggestions to  npb@nas.nasa.gov  		  !
!									  !
!	   NAS Parallel Benchmarks Group				  !
!	   NASA Ames Research Center					  !
!	   Mail Stop: T27A-1						  !
!	   Moffett Field, CA   94035-1000				  !
!									  !
!	   E-mail:  npb@nas.nasa.gov					  !
!	   Fax:     (650) 604-3957					  !
!									  !
!-------------------------------------------------------------------------!
! Authors: D. Bailey 					                  !
!	   W. Saphir							  !
! Translation to Java and MultiThreaded Code	 			  !
!	   M. Frumkin							  !
!	   M. Schultz							  !
! Translation to C#	 			  !
!      Cenez Araújo de Rezende, MDCC/UFC
!-------------------------------------------------------------------------!
*/

using System;
using System.IO;
using NPB3_0_JAV.FTThreads;
using NPB3_0_JAV.BMInOut;

namespace NPB3_0_JAV {

    public class FT : FTBase
    {
        public int bid = -1;
        public BMResults results;
        public bool serial = true;
        bool done = false;
        public FT(char clss, int np, bool ser) : base(clss, np, ser)
        {
            //super(clss, np, ser);
            serial = ser;
        }
        static void Main(String[] argv)
        {
            FT ft = null;

            BMArgs.ParseCmdLineArgs(argv, BMName);
            char CLSS = BMArgs.CLASS;
            int np = BMArgs.num_threads;
            bool serial = BMArgs.serial;

            try
            {
                ft = new FT(CLSS, np, serial);
            }
            catch (OutOfMemoryException e)
            {
                BMArgs.outOfMemoryMessage();
                Environment.Exit(0);
            }
            ft.runBenchMark();
        }
        public void runBenchMark()
        {
            BMArgs.Banner(BMName, CLASS, serial, num_threads);
            Console.WriteLine(" Size = " + nx + " X " + ny + " X " + nz
                               + " niter = " + niter_default);
            setTimers();
            timer.resetAllTimers();

            if (serial) appft_serial();

            if (timeron) timer.start(14);
            int verified = verify(4, nx, ny, nz, niter_default, checksum);
            if (timeron) timer.stop(14);
            timer.stop(1);

            double time = timer.readTimer(1);
            results = new BMResults(BMName,
                          CLASS,
                          nx,
                          ny,
                          nz,
                          niter_default,
                          time,
                          getMFLOPS(time, nx, ny, nz),
                          "floating point",
                          verified,
                          serial,
                          num_threads,
                          bid);
            results.print();
            if (timeron) printTimers();
            done = true;
        }

        public void appft_serial() {
            if (timeron) timer.start(2);
            initial_conditions(xtr, ny, nx, nz);
            CompExp(nx, exp1);
            CompExp(ny, exp2);
            CompExp(nz, exp3);
            fftXYZ(1, xtr, exp2, exp1, exp3, ny, nx, nz);
            if (timeron) timer.stop(2);

            timer.start(1);
            if (timeron) timer.start(12);
            initial_conditions(xtr, ny, nx, nz);
            if (timeron) timer.stop(12);
            if (timeron) timer.start(15);
            fftXYZ(1, xtr, exp2, exp1, exp3, ny, nx, nz);
            if (timeron) timer.stop(15);

            double ap = (-4.0 * alpha * pow2(pi));
            int n12 = nx / 2;
            int n22 = ny / 2;
            int n32 = nz / 2;

            for (int it = 0; it < niter_default; it++)
            {
                if (timeron) timer.start(11);

                for (int i = 0; i < nx; i++)
                {
                    int ii = i - ((i) / n12) * nx;
                    int ii2 = ii * ii;
                    for (int k = 0; k < nz; k++)
                    {
                        int kk = k - ((k) / n32) * nz;
                        int ik2 = ii2 + kk * kk;
                        for (int j = 0; j < ny; j++){
                            int jj = j - ((j) / n22) * ny;
                            //xnt[REAL+j*isize4+k*jsize4+i*ksize4] = xtr[REAL+j*isize3+i*jsize3+k*ksize3] * Math.Exp((ap*(jj*jj+ik2))*(it+1));
                            //xnt[IMAG+j*isize4+k*jsize4+i*ksize4] = xtr[IMAG+j*isize3+i*jsize3+k*ksize3] * Math.Exp((ap*(jj*jj+ik2))*(it + 1));
                            xnt[REAL,j,k,i] = xtr[REAL,j,i,k] * Math.Exp((ap*(jj*jj+ik2))*(it+1));
                            xnt[IMAG,j,k,i] = xtr[IMAG,j,i,k] * Math.Exp((ap*(jj*jj+ik2))*(it+1));
                        }
                    }
                }
                if (timeron) timer.stop(11);

                if (timeron) timer.start(15);
                fftXYZ(-1, xnt, exp2, exp3, exp1, ny, nz, nx);
                if (timeron) timer.stop(15);

                if (timeron) timer.start(10);
                CalculateChecksum(checksum,REAL+it,it,xnt,ny,nz,nx);  // CalculateChecksum(checksum,REAL+it*isize2,it,xnt,ny,nz,nx);
                if (timeron) timer.stop(10);
            }
        }
        public void setTimers()
        {
            //File f1 = new File("timer.flag");
            timeron = false;
            if (File.Exists("timer.flag")) timeron = true;
        }
        public void printTimers()
        {
            //DecimalFormat fmt = new DecimalFormat("0.000");
            Console.WriteLine("  SECTION   Time (secs)");
            Console.WriteLine("FT time =		      " + timer.readTimer(1).ToString("N3"));
            Console.WriteLine("WarmUp time =		      " + timer.readTimer(2).ToString("N3"));
            Console.WriteLine("ffXYZ body time =	      " + timer.readTimer(3).ToString("N3"));
            Console.WriteLine("Swarztrauber body time =      " + timer.readTimer(4).ToString("N3"));
            Console.WriteLine("Redistribution time =	      " + timer.readTimer(5).ToString("N3"));
            Console.WriteLine("Transposition time =	      " + timer.readTimer(6).ToString("N3"));
            Console.WriteLine("X time =		      " + timer.readTimer(7).ToString("N3"));
            Console.WriteLine("Y time =		      " + timer.readTimer(8).ToString("N3"));
            Console.WriteLine("Z time =		      " + timer.readTimer(9).ToString("N3"));
            Console.WriteLine("CalculateChecksum =	      " + timer.readTimer(10).ToString("N3"));
            Console.WriteLine("evolve =		      " + timer.readTimer(11).ToString("N3"));
            Console.WriteLine("compute_initial_conditions =  " + timer.readTimer(12).ToString("N3"));
            Console.WriteLine("twiddle =		      " + timer.readTimer(13).ToString("N3"));
            Console.WriteLine("verify =		      " + timer.readTimer(14).ToString("N3"));
            Console.WriteLine("fftXYZ =		      " + timer.readTimer(15).ToString("N3"));
        }

        public double getMFLOPS(double total_time, int nx, int ny, int nz)
        {
            double mflops = 0.0;
            int ntotal = nx * ny * nz;
            if (total_time > 0)
            {
                mflops = 14.8157 + 7.19641 * Math.Log(ntotal) + (5.23518 + 7.21113 * Math.Log(ntotal)) * niter_default;
                mflops *= ntotal / (total_time * 1000000.0);
            }
            return mflops;
        }

        public void CalculateChecksum(double[,] csum, int csmffst, int iterN, double[,,,] u, int d1, int d2, int d3)
        {
            int i, ii, ji, ki;
            int isize3 = 2,
                jsize3 = isize3 * (d1 + 1),
                ksize3 = jsize3 * d2;
            csum[REAL,csmffst] = 0.0;
            csum[IMAG,csmffst] = 0.0;

            double csumr = 0.0, csumi = 0.0;
            for (i = 1; i <= 1024; i++)
            {
                ii = (1 * i) % d3;
                ji = (3 * i) % d1;
                ki = (5 * i) % d2;
                csumr += u[ REAL, ji, ki,ii]; // csumr += u[ REAL, ji, ki,ii];
                csumi += u[ IMAG, ji, ki,ii]; // csumi += u[ IMAG, ji, ki,ii];
            }
            csum[ REAL,csmffst] = csumr / (d1 * d2 * d3);
            csum[ IMAG,csmffst] = csumi / (d1 * d2 * d3);
        }

        public void fftXYZ(int sign, double[,,,] x, double[,] exp1, double[,] exp2, double[,] exp3, int n1, int n2, int n3)
        {
            int i = 0, j = 0, k, log;
            int isize3 = 2, jsize3, ksize3;
            jsize3 = isize3 * (n1 + 1);
            ksize3 = jsize3 * n2;


            if (timeron) timer.start(3);

            log = ilog2(n2);
            if (timeron) timer.start(7);
            for (k = 0; k < n3; k++) Swarztrauber(sign, log, n1, n2, x, k, n1, exp2, scr);
            if (timeron) timer.stop(7);

            log = ilog2(n1);
            if (timeron) timer.start(8);
            for (k = 0; k < n3; k++)
            {
                for (j = 0; j < n2; j++)
                {
                    for (i = 0; i < n1; i++)
                    {
                        plane[ REAL, j, i,0] = x[ REAL, i, j,k];
                        plane[ IMAG, j, i,0] = x[ IMAG, i, j,k];
                    }
                }
                Swarztrauber(sign, log, n2, n1, plane, 0, n2, exp1, scr);
                for (j = 0; j < n2; j++)
                {
                    for (i = 0; i < n1; i++)
                    {
                        x[REAL,i,j,k] = plane[REAL,j,i,0];
                        x[IMAG,i,j,k] = plane[IMAG,j,i,0];
                    }
                }
            }
            if (timeron) timer.stop(8);

            log = ilog2(n3);
            if (timeron) timer.start(9);
            for (k = 0; k < n2; k++)
            {
                for (i = 0; i < n3; i++)
                {
                    for (j = 0; j < n1; j++)
                    {
                        plane[REAL,j,i,0] = x[REAL,j,k,i];
                        plane[IMAG,j,i,0] = x[IMAG,j,k,i];
                    }
                }
                Swarztrauber(sign, log, n1, n3, plane, 0, n1, exp3, scr);
                for (i = 0; i < n3; i++)
                {
                    for (j = 0; j < n1; j++)
                    {
                        x[REAL,j,k,i] = plane[REAL,j,i,0];
                        x[IMAG,j,k,i] = plane[IMAG,j,i,0];
                    }
                }
            }

            if (timeron) timer.stop(9);
            if (timeron) timer.stop(3);
        }

        public int verify(int ires, int n1, int n2, int n3, int nt, double[,] cksum)
        {
            int verified = -1;
            bool[] temp = new bool[niter_default];
            double[,] cexpd = new double[2,21];
            if ((n1 == 64) && (n2 == 64) && (n3 == 64) && (nt == 6)) {
                //
                // Class S reference values.
                //
                cexpd[REAL,0] = 554.6087004964;
                cexpd[REAL,1] = 554.6385409189;
                cexpd[REAL,2] = 554.6148406171;
                cexpd[REAL,3] = 554.5423607415;
                cexpd[REAL,4] = 554.4255039624;
                cexpd[REAL,5] = 554.2683411902;

                cexpd[IMAG,0] = 484.5363331978;
                cexpd[IMAG,1] = 486.5304269511;
                cexpd[IMAG,2] = 488.3910722336;
                cexpd[IMAG,3] = 490.1273169046;
                cexpd[IMAG,4] = 491.7475857993;
                cexpd[IMAG,5] = 493.2597244941;

            }
            else if ((n1 == 128) && (n2 == 128) && (n3 == 32) && (nt == 6)) {
                //
                // Class W reference values.
                //
                cexpd[REAL,0] = 567.3612178944;
                cexpd[REAL,1] = 563.1436885271;
                cexpd[REAL,2] = 559.4024089970;
                cexpd[REAL,3] = 556.0698047020;
                cexpd[REAL,4] = 553.0898991250;
                cexpd[REAL,5] = 550.4159734538;

                cexpd[IMAG,0] = 529.3246849175;
                cexpd[IMAG,1] = 528.2149986629;
                cexpd[IMAG,2] = 527.0996558037;
                cexpd[IMAG,3] = 526.0027904925;
                cexpd[IMAG,4] = 524.9400845633;
                cexpd[IMAG,5] = 523.9212247086;
                //
            }
            else if ((n1 == 256) && (n2 == 256) && (n3 == 128) && (nt == 6)){
                //
                // Class A reference values.
                //
                cexpd[REAL,0] = 504.6735008193;
                cexpd[REAL,1] = 505.9412319734;
                cexpd[REAL,2] = 506.9376896287;
                cexpd[REAL,3] = 507.7892868474;
                cexpd[REAL,4] = 508.5233095391;
                cexpd[REAL,5] = 509.1487099959;

                cexpd[IMAG,0] = 511.4047905510;
                cexpd[IMAG,1] = 509.8809666433;
                cexpd[IMAG,2] = 509.8144042213;
                cexpd[IMAG,3] = 510.1336130759;
                cexpd[IMAG,4] = 510.4914655194;
                cexpd[IMAG,5] = 510.7917842803;
                //
            }
            else if ((n1 == 512) && (n2 == 256) && (n3 == 256) && (nt == 20)){
                //
                // Class B reference values.
                //
                cexpd[REAL,0] = 517.7643571579;
                cexpd[REAL,1] = 515.4521291263;
                cexpd[REAL,2] = 514.6409228649;
                cexpd[REAL,3] = 514.2378756213;
                cexpd[REAL,4] = 513.9626667737;
                cexpd[REAL,5] = 513.7423460082;
                cexpd[REAL,6] = 513.5547056878;
                cexpd[REAL,7] = 513.3910925466;
                cexpd[REAL,8] = 513.2470705390;
                cexpd[REAL,9] = 513.1197729984;
                cexpd[REAL,10] = 513.0070319283;
                cexpd[REAL,11] = 512.9070537032;
                cexpd[REAL,12] = 512.8182883502;
                cexpd[REAL,13] = 512.7393733383;
                cexpd[REAL,14] = 512.6691062020;
                cexpd[REAL,15] = 512.6064276004;
                cexpd[REAL,16] = 512.5504076570;
                cexpd[REAL,17] = 512.5002331720;
                cexpd[REAL,18] = 512.4551951846;
                cexpd[REAL,19] = 512.4146770029;

                cexpd[IMAG,0] = 507.7803458597;
                cexpd[IMAG,1] = 508.8249431599;
                cexpd[IMAG,2] = 509.6208912659;
                cexpd[IMAG,3] = 510.1023387619;
                cexpd[IMAG,4] = 510.3976610617;
                cexpd[IMAG,5] = 510.5948019802;
                cexpd[IMAG,6] = 510.7404165783;
                cexpd[IMAG,7] = 510.8576573661;
                cexpd[IMAG,8] = 510.9577278523;
                cexpd[IMAG,9] = 511.0460304483;
                cexpd[IMAG,10] = 511.1252433800;
                cexpd[IMAG,11] = 511.1968077718;
                cexpd[IMAG,12] = 511.2616233064;
                cexpd[IMAG,13] = 511.3203605551;
                cexpd[IMAG,14] = 511.3735928093;
                cexpd[IMAG,15] = 511.4218460548;
                cexpd[IMAG,16] = 511.4656139760;
                cexpd[IMAG,17] = 511.5053595966;
                cexpd[IMAG,18] = 511.5415130407;
                cexpd[IMAG,19] = 511.5744692211;
                //
            }
            else if ((n1 == 512) && (n2 == 512) &&
                     (n3 == 512) && (nt == 20))
            {
                //
                // Class C reference values.
                //
                cexpd[REAL,0] = 519.5078707457;
                cexpd[REAL,1] = 515.5422171134;
                cexpd[REAL,2] = 514.4678022222;
                cexpd[REAL,3] = 514.0150594328;
                cexpd[REAL,4] = 513.7550426810;
                cexpd[REAL,5] = 513.5811056728;
                cexpd[REAL,6] = 513.4569343165;
                cexpd[REAL,7] = 513.3651975661;
                cexpd[REAL,8] = 513.2955192805;
                cexpd[REAL,9] = 513.2410471738;
                cexpd[REAL,10] = 513.1971141679;
                cexpd[REAL,11] = 513.1605205716;
                cexpd[REAL,12] = 513.1290734194;
                cexpd[REAL,13] = 513.1012720314;
                cexpd[REAL,14] = 513.0760908195;
                cexpd[REAL,15] = 513.0528295923;
                cexpd[REAL,16] = 513.0310107773;
                cexpd[REAL,17] = 513.0103090133;
                cexpd[REAL,18] = 512.9905029333;
                cexpd[REAL,19] = 512.9714421109;

                cexpd[IMAG,0] = 514.9019699238;
                cexpd[IMAG,1] = 512.7578201997;
                cexpd[IMAG,2] = 512.2251847514;
                cexpd[IMAG,3] = 512.1090289018;
                cexpd[IMAG,4] = 512.1143685824;
                cexpd[IMAG,5] = 512.1496764568;
                cexpd[IMAG,6]= 512.1870921893;
                cexpd[IMAG,7] = 512.2193250322;
                cexpd[IMAG,8] = 512.2454735794;
                cexpd[IMAG,9] = 512.2663649603;
                cexpd[IMAG,10] = 512.2830879827;
                cexpd[IMAG,11] = 512.2965869718;
                cexpd[IMAG,12] = 512.3075927445;
                cexpd[IMAG,13] = 512.3166486553;
                cexpd[IMAG,14] = 512.3241541685;
                cexpd[IMAG,15] = 512.3304037599;
                cexpd[IMAG,16] = 512.3356167976;
                cexpd[IMAG,17] = 512.3399592211;
                cexpd[IMAG,18] = 512.3435588985;
                cexpd[IMAG,19] = 512.3465164008;
            }
            double epsilon = 1.0E-12;
            //
            // Verification test for results.
            //
            if (nt <= 0)
            {

            }
            else
            {
                for (int it = 0; it < nt; it++)
                {
                    double csumr = (cksum[REAL,it] - cexpd[REAL,it]) / cexpd[REAL,it];
                    double csumi = (cksum[IMAG,it] - cexpd[IMAG,it]) / cexpd[IMAG,it];
                    if (Math.Abs(csumr) <= epsilon
                     || Math.Abs(csumi) <= epsilon
                   )
                    {
                        if (verified == -1) verified = 1;
                    }
                    else
                    {
                        verified = 0;
                    }
                }
            }
            BMResults.printVerificationStatus(CLASS, verified, BMName);
            return verified;
        }

        public double getTime() { return timer.readTimer(1); }
        public bool isDone() { return done; }
        public void finalize() { // throws Throwable{
            Console.WriteLine("FT: is about to be garbage collected");
            //super.finalize();
        }
    }
}