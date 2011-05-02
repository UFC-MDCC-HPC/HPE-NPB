/*
!-------------------------------------------------------------------------!
!									  !
!	 N  A  S     P A R A L L E L	 B E N C H M A R K S  3.0	  !
!									  !
!			J A V A 	V E R S I O N			  !
!									  !
!                               B M A R G S                               !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    BMArgs implements Command Line Benchmark Arguments class             !
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
!     Translation to C# :          			  !
!     Francisco Heron de Carvalho Junior (MDCC/UFC)					          !
!-------------------------------------------------------------------------!
*/

using System;

namespace NPB3_0_JAV.BMInOut {

    [Serializable()]
    public class BMArgs
    {
        public static char CLASS = 'U';
        public static int num_threads = 4;
        public static bool serial = true;
        public BMArgs()
        {
            CLASS = 'U';
            num_threads = 4;
            serial = true;
        }
        static public void ParseCmdLineArgs(string[] argv, string BMName)
        {
            for (int i = 0; i < argv.Length; i++)
            {
                if (argv[i].Equals("SERIAL")
                   || argv[i].Equals("serial")
                   || argv[i].Equals("-serial")
                   || argv[i].Equals("-SERIAL"))
                {
                    serial = true;
                }
                else
                    if (argv[i].StartsWith("class=")
                       || argv[i].StartsWith("CLASS=")
                       || argv[i].StartsWith("-class")
                       || argv[i].StartsWith("-CLASS"))
                    {

                        if (argv[i].Length > 6)
                            CLASS = Char.ToUpper(argv[i][6]);
                        if (CLASS != 'A' && CLASS != 'B' && CLASS != 'C' && CLASS != 'S' && CLASS != 'W')
                        {
                            Console.WriteLine("classes allowed are A,B,C,W and S.");
                            commandLineError(BMName);
                        }
                    }
                    else if (argv[i].StartsWith("np=")
                             || argv[i].StartsWith("NP=")
                             || argv[i].StartsWith("-NP")
                             || argv[i].StartsWith("-np"))
                    {
                        try
                        {
                            if (argv[i].Length > 3)
                                num_threads = Int32.Parse(argv[i].Substring(3));
                            serial = false;
                        }
                        catch (Exception e)
                        {
                            Console.WriteLine("argument to " + argv[i].Substring(0, 3)
                                       + " must be an integer.");
                            commandLineError(BMName);
                        }
                    }
            }
        }
        public static void commandLineError(String BMName)
        {
            Console.WriteLine("synopsis: java " + BMName
                          + " CLASS=[ABCWS] -serial [-NPnnn]");
            Console.WriteLine("[ABCWS] is the size class \n"
                           + "-serial specifies the serial version and\n"
                           + "-NP specifies number of threads where nnn "
                           + "is an integer");
           Environment.Exit(0);
        }
        public static void outOfMemoryMessage()
        {
            Console.WriteLine("The java maximum heap size is "
                           + "to small to run this benchmark class");
            Console.WriteLine("To allocate more memory, use the -mxn option"
                           + " where n is the number of bytes to be allocated");
        }
        public static void Banner(String BMName,
                                  char clss, bool serial, int np)
        {
            Console.WriteLine(" NAS Parallel Benchmarks C# version (NPB3_0_CS)");
            if (serial) Console.WriteLine(" Serial Version " + BMName + "." + clss);
            else Console.WriteLine(" Multithreaded Version " + BMName + "." + clss +
                                    " np=" + np);
        }
    }

}
