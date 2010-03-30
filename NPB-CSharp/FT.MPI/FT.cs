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
            //char CLSS = BMArgs.CLASS;
            char CLSS = 'S';
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

            double ap = (-4.0 * alpha * Math.Pow(pi, 2));
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
                            xnt[i,k,j,REAL] = xtr[k,i,j,REAL] * Math.Exp((ap*(jj*jj+ik2))*(it+1));
                            xnt[i,k,j,IMAG] = xtr[k,i,j,IMAG] * Math.Exp((ap*(jj*jj+ik2))*(it+1));
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
            csum[csmffst,REAL] = 0.0;
            csum[csmffst,IMAG] = 0.0;

            double csumr = 0.0, csumi = 0.0;
            for (i = 1; i <= 1024; i++)
            {
                ii = (1 * i) % d3;
                ji = (3 * i) % d1;
                ki = (5 * i) % d2;
                csumr += u[ii, ki, ji, REAL]; // csumr += u[ii, ki, ji, REAL];
                csumi += u[ii, ki, ji, IMAG]; // csumi += u[ii, ki, ji, IMAG];
            }
            csum[csmffst, REAL] = csumr / (d1 * d2 * d3);
            csum[csmffst, IMAG] = csumi / (d1 * d2 * d3);
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
                        plane[0, i, j, REAL] = x[k, j, i, REAL];
                        plane[0, i, j, IMAG] = x[k, j, i, IMAG];
                    }
                }
                Swarztrauber(sign, log, n2, n1, plane, 0, n2, exp1, scr);
                for (j = 0; j < n2; j++)
                {
                    for (i = 0; i < n1; i++)
                    {
                        x[k,j,i,REAL] = plane[0,i,j,REAL];
                        x[k,j,i,IMAG] = plane[0,i,j,IMAG];
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
                        plane[0,i,j,REAL] = x[i,k,j,REAL];
                        plane[0,i,j,IMAG] = x[i,k,j,IMAG];
                    }
                }
                Swarztrauber(sign, log, n1, n3, plane, 0, n1, exp3, scr);
                for (i = 0; i < n3; i++)
                {
                    for (j = 0; j < n1; j++)
                    {
                        x[i,k,j,REAL] = plane[0,i,j,REAL];
                        x[i,k,j,IMAG] = plane[0,i,j,IMAG];
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
            double[,] cexpd = new double[21,2];
            if ((n1 == 64) && (n2 == 64) && (n3 == 64) && (nt == 6)) {
                //
                // Class S reference values.
                //
                cexpd[0,REAL] = 554.6087004964;
                cexpd[1,REAL] = 554.6385409189;
                cexpd[2,REAL] = 554.6148406171;
                cexpd[3,REAL] = 554.5423607415;
                cexpd[4,REAL] = 554.4255039624;
                cexpd[5,REAL] = 554.2683411902;

                cexpd[0,IMAG] = 484.5363331978;
                cexpd[1,IMAG] = 486.5304269511;
                cexpd[2,IMAG] = 488.3910722336;
                cexpd[3,IMAG] = 490.1273169046;
                cexpd[4,IMAG] = 491.7475857993;
                cexpd[5,IMAG] = 493.2597244941;

            }
            else if ((n1 == 128) && (n2 == 128) && (n3 == 32) && (nt == 6)) {
                //
                // Class W reference values.
                //
                cexpd[0,REAL] = 567.3612178944;
                cexpd[1,REAL] = 563.1436885271;
                cexpd[2,REAL] = 559.4024089970;
                cexpd[3,REAL] = 556.0698047020;
                cexpd[4,REAL] = 553.0898991250;
                cexpd[5,REAL] = 550.4159734538;

                cexpd[0,IMAG] = 529.3246849175;
                cexpd[1,IMAG] = 528.2149986629;
                cexpd[2,IMAG] = 527.0996558037;
                cexpd[3,IMAG] = 526.0027904925;
                cexpd[4,IMAG] = 524.9400845633;
                cexpd[5,IMAG] = 523.9212247086;
                //
            }
            else if ((n1 == 256) && (n2 == 256) && (n3 == 128) && (nt == 6)){
                //
                // Class A reference values.
                //
                cexpd[0,REAL] = 504.6735008193;
                cexpd[1,REAL] = 505.9412319734;
                cexpd[2,REAL] = 506.9376896287;
                cexpd[3,REAL] = 507.7892868474;
                cexpd[4,REAL] = 508.5233095391;
                cexpd[5,REAL] = 509.1487099959;

                cexpd[0,IMAG] = 511.4047905510;
                cexpd[1,IMAG] = 509.8809666433;
                cexpd[2,IMAG] = 509.8144042213;
                cexpd[3,IMAG] = 510.1336130759;
                cexpd[4,IMAG] = 510.4914655194;
                cexpd[5,IMAG] = 510.7917842803;
                //
            }
            else if ((n1 == 512) && (n2 == 256) && (n3 == 256) && (nt == 20)){
                //
                // Class B reference values.
                //
                cexpd[0,REAL] = 517.7643571579;
                cexpd[1,REAL] = 515.4521291263;
                cexpd[2,REAL] = 514.6409228649;
                cexpd[3,REAL] = 514.2378756213;
                cexpd[4,REAL] = 513.9626667737;
                cexpd[5,REAL] = 513.7423460082;
                cexpd[6,REAL] = 513.5547056878;
                cexpd[7,REAL] = 513.3910925466;
                cexpd[8,REAL] = 513.2470705390;
                cexpd[9,REAL] = 513.1197729984;
                cexpd[10,REAL] = 513.0070319283;
                cexpd[11,REAL] = 512.9070537032;
                cexpd[12,REAL] = 512.8182883502;
                cexpd[13,REAL] = 512.7393733383;
                cexpd[14,REAL] = 512.6691062020;
                cexpd[15,REAL] = 512.6064276004;
                cexpd[16,REAL] = 512.5504076570;
                cexpd[17,REAL] = 512.5002331720;
                cexpd[18,REAL] = 512.4551951846;
                cexpd[19,REAL] = 512.4146770029;

                cexpd[0,IMAG] = 507.7803458597;
                cexpd[1,IMAG] = 508.8249431599;
                cexpd[2,IMAG] = 509.6208912659;
                cexpd[3,IMAG] = 510.1023387619;
                cexpd[4,IMAG] = 510.3976610617;
                cexpd[5,IMAG] = 510.5948019802;
                cexpd[6,IMAG] = 510.7404165783;
                cexpd[7,IMAG] = 510.8576573661;
                cexpd[8,IMAG] = 510.9577278523;
                cexpd[9,IMAG] = 511.0460304483;
                cexpd[10,IMAG] = 511.1252433800;
                cexpd[11,IMAG] = 511.1968077718;
                cexpd[12,IMAG] = 511.2616233064;
                cexpd[13,IMAG] = 511.3203605551;
                cexpd[14,IMAG] = 511.3735928093;
                cexpd[15,IMAG] = 511.4218460548;
                cexpd[16,IMAG] = 511.4656139760;
                cexpd[17,IMAG] = 511.5053595966;
                cexpd[18,IMAG] = 511.5415130407;
                cexpd[19,IMAG] = 511.5744692211;
                //
            }
            else if ((n1 == 512) && (n2 == 512) &&
                     (n3 == 512) && (nt == 20))
            {
                //
                // Class C reference values.
                //
                cexpd[0,REAL] = 519.5078707457;
                cexpd[1,REAL] = 515.5422171134;
                cexpd[2,REAL] = 514.4678022222;
                cexpd[3,REAL] = 514.0150594328;
                cexpd[4,REAL] = 513.7550426810;
                cexpd[5,REAL] = 513.5811056728;
                cexpd[6,REAL] = 513.4569343165;
                cexpd[7,REAL] = 513.3651975661;
                cexpd[8,REAL] = 513.2955192805;
                cexpd[9,REAL] = 513.2410471738;
                cexpd[10,REAL] = 513.1971141679;
                cexpd[11,REAL] = 513.1605205716;
                cexpd[12,REAL] = 513.1290734194;
                cexpd[13,REAL] = 513.1012720314;
                cexpd[14,REAL] = 513.0760908195;
                cexpd[15,REAL] = 513.0528295923;
                cexpd[16,REAL] = 513.0310107773;
                cexpd[17,REAL] = 513.0103090133;
                cexpd[18,REAL] = 512.9905029333;
                cexpd[19,REAL] = 512.9714421109;

                cexpd[0,IMAG] = 514.9019699238;
                cexpd[1,IMAG] = 512.7578201997;
                cexpd[2,IMAG] = 512.2251847514;
                cexpd[3,IMAG] = 512.1090289018;
                cexpd[4,IMAG] = 512.1143685824;
                cexpd[5,IMAG] = 512.1496764568;
                cexpd[6,IMAG]= 512.1870921893;
                cexpd[7,IMAG] = 512.2193250322;
                cexpd[8,IMAG] = 512.2454735794;
                cexpd[9,IMAG] = 512.2663649603;
                cexpd[10,IMAG] = 512.2830879827;
                cexpd[11,IMAG] = 512.2965869718;
                cexpd[12,IMAG] = 512.3075927445;
                cexpd[13,IMAG] = 512.3166486553;
                cexpd[14,IMAG] = 512.3241541685;
                cexpd[15,IMAG] = 512.3304037599;
                cexpd[16,IMAG] = 512.3356167976;
                cexpd[17,IMAG] = 512.3399592211;
                cexpd[18,IMAG] = 512.3435588985;
                cexpd[19,IMAG] = 512.3465164008;
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
                    double csumr = (cksum[it,REAL] - cexpd[it,REAL]) / cexpd[it,REAL];
                    double csumi = (cksum[it,IMAG] - cexpd[it,IMAG]) / cexpd[it,IMAG];
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