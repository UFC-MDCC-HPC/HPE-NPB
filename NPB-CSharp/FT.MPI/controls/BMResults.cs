/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                           B M R E S A L T S                             !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    BMResults implements Benchmark Result class                          !
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
!     Translation to Java and to MultiThreaded Code:			  !
!     Michael A. Frumkin					          !
!-------------------------------------------------------------------------!
*/

using System;
using System.Text;

namespace NPB3_0_JAV.BMInOut {

    [Serializable()]
    public class BMResults
    {
        public String name;
        public String MachineName;
        public String PrLang;
        public char clss;
        public int n1, n2, n3, niter;
        public double time, acctime, wctime, mops;
        public double tmSent = 0.0, tmReceived = 0.0;
        public int RecArrSize = 0;
        public String optype;
        public int numthreads;
        public bool serial;
        public int pid;
        public int verified;
        public System.IO.StreamWriter out1 = null; 

        public BMResults() { }

        public BMResults(int bid)
        {
            pid = bid;
            clss = 'S';
            optype = "floating point";
        }

        public BMResults(String bname,
                 char CLASS,
                 int bn1,
                 int bn2,
                 int bn3,
                 int bniter,
                 double btime,
                 double bmops,
                 String boptype,
                 int passed_verification,
                 bool bserial,
                     int num_threads,
                     int bid)
        {
            pid = bid;
            name = bname;
            clss = CLASS;
            n1 = bn1;
            n2 = bn2;
            n3 = bn3;
            niter = bniter;
            time = btime;
            mops = bmops;
            optype = boptype;
            verified = passed_verification;
            serial = bserial;
            numthreads = num_threads;
        }

        public void print()
        {
            StringBuilder outbuf = new StringBuilder("                                "
                                                 + "                               *"); 
            int len = outbuf.Length;
            string outline = "***** NAS Parallel Benchmarks" +
                             " Java version (NPB3_0_JAV) " + name + " ****";
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* Class             = " + clss;
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            if (n2 == 0 && n3 == 0)
            {
                outline = "* Size              = " + n1;
            }
            else
            {
                outline = "* Size              = " + n1 + " X " + n2 + " X " + n3;
            }
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* Iterations        = " + niter;
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* Time in seconds   = " + time.ToString("N3");
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* ACCTime           = " + acctime.ToString("N3");
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* Mops total        = " + mops.ToString("N3");
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* Operation type    = " + optype;
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            if (verified == 1) outline ="* Verification      = Successful";
            else if (verified == 0) outline = "* Verification      = Failed";
            else outline ="* Verification      = Not Performed";
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            if (!serial)
            {
                outbuf = new StringBuilder("                          "
                                        + "                         *");
                outline = "* Threads requested = " + numthreads;
                outbuf.Insert(0, outline);
                outbuf.Length = len;
                outbuf.Insert(len - 1, "*");
                Console.WriteLine(outbuf.ToString());
            }

            outbuf = new StringBuilder("*                               "
                                    + "                               *");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* Please send all errors/feedbacks to:";
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* NPB Working Team";
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "* npb@nas.nasa.gov";
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            outbuf = new StringBuilder("                          "
                                    + "                         *");
            outline = "********************************"
                               + "*******************************";
            outbuf.Insert(0, outline);
            outbuf.Length = len;
            outbuf.Insert(len - 1, "*");
            Console.WriteLine(outbuf.ToString());

            if (out1 != null)
            {
                try
                {
                    outline = "***** NAS Parallel Benchmarks Java version (NPB3_0_JAV) "
                                  + name + " Report *****";
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = "Class           = " + clss.ToString();
                    out1.WriteLine(outline, 0, outline.Length);
                    if (n2 == 0 && n3 == 0)
                    {
                        outline = "Size            = " + n1.ToString();
                    }
                    else
                    {
                        outline = "Size            = " + n1.ToString() + " X " +
                                                         n2.ToString() + " X " +
                                                         n3.ToString();
                    }
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = "Iterations      = " + niter.ToString();
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = "Time in seconds = " + time.ToString("N3");
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = "ACCTime         = " + acctime.ToString("N3");
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = "Mops total      = " + mops.ToString("N3");
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = "Operation type  = " + optype.ToString();
                    out1.WriteLine(outline, 0, outline.Length);
                    if (verified == 1) outline = "Verification    = Successful";
                    else if (verified == 0) outline = "Verification Failed";
                    else outline = "Verification Not Performed";
                    out1.WriteLine(outline, 0, outline.Length);

                    outline = "\n Please send all errors/feedbacks to:";
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = " NPB Working Team";
                    out1.WriteLine(outline, 0, outline.Length);
                    outline = " npb@nas.nasa.gov\n";
                    out1.WriteLine(outline, 0, outline.Length);
                    out1.Flush();
                }
                catch (Exception e)
                {
                    Console.Error.WriteLine("Res.print: write file: " + e.ToString());
                }
            }
        }

        public int getFromFile(String filename)
        {
            System.IO.StreamReader in1 = null;
            verified = -1;
            try
            {
                // System.IO.File file = new System.IO.File(filename);
                System.IO.FileStream fs = new System.IO.FileStream(filename, System.IO.FileMode.Open);

                in1 = new System.IO.StreamReader(fs);
            }
            catch (Exception e)
            {
                Console.Error.WriteLine("BMResults.getFromFile: filename " + e.ToString());
                return 0;
            }
            string line;
            string keyword;
            int idx1;
            try
            {
                while ((line = in1.ReadLine()) != null)
                {
                    if (line.IndexOf("Time in seconds =") >= 0)
                    {
                        keyword = "Time in seconds =";
                        idx1 = line.IndexOf(keyword);
                        idx1 += keyword.Length;
                        double dbl;
                        Double.TryParse(line.Substring(idx1), out dbl);
                        acctime = wctime = dbl;
                    }
                    else if (line.IndexOf("Verification    =") >= 0)
                    {
                        verified = 0;
                        if (line.IndexOf("successful") >= 0 && line.IndexOf("successful") < 0
                       || line.IndexOf("SUCCESSFUL") >= 0 && line.IndexOf("UNSUCCESSFUL") < 0)
                            verified = 1;
                    }
                    else if (line.IndexOf("Mop/s total     =") >= 0)
                    {
                        keyword = "Mop/s total     =";
                        idx1 = line.IndexOf(keyword);
                        idx1 += keyword.Length;
                        Double.TryParse(line.Substring(idx1),out mops);
                    }
                }
            }
            catch (Exception e)
            {
                Console.Error.WriteLine("BMResults.getFromFile: " + e.ToString());
                return 0;
            }
            //    print();
            return 1;
        }
        public static void printVerificationStatus(char clss, int verified, String BMName)
        {
            if (clss == 'U' || verified == -1)
            {
                verified = -1;
                Console.WriteLine(" Problem size unknown");
                Console.WriteLine(BMName + "." + clss + ": Verification Not Performed");
            }
            else if (verified == 1)
            {
                Console.WriteLine(BMName + "." + clss + ": Verification Successful");
            }
            else
            {
                Console.WriteLine(BMName + "." + clss + ": Verification Failed");
            }
        }
        public static int printComparisonStatus(char clss, int verified, double epsilon,
                             double[] xcr, double[] xcrref, double[] xcrdif)
        {
            for (int m = 0; m < xcr.Length; m++)
            {
                if (clss == 'U')
                {
                    Console.WriteLine(m + ". " + xcr[m]);
                }
                else
                {
                    if (xcrdif[m] <= epsilon)
                    {
                        if (verified == -1) verified = 1;
                    }
                    else
                    {
                        verified = 0;
                        Console.Write("FAILURE: ");
                    }
                    Console.WriteLine(m + ". " + xcr[m] + " " + xcrref[m] + " " + xcrdif[m]);
                }
            }
            return verified;
        }
        public static int printComparisonStatus(char clss, int verified, double epsilon,
                             double xcr, double xcrref, double xcrdif)
        {
            if (clss == 'U')
            {
                Console.WriteLine(" " + xcr);
            }
            else
            {
                if (xcrdif <= epsilon)
                {
                    if (verified == -1) verified = 1;
                }
                else
                {
                    verified = 0;
                    Console.Write("FAILURE: ");
                }
                Console.WriteLine(xcr + " " + xcrref + " " + xcrdif);
            }
            return verified;
        }
    }
}
