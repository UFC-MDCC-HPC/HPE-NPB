using System;

namespace NPB.Ftc {

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
                        catch (Exception e) {
                            Console.WriteLine(e.ToString());
                            Console.WriteLine("argument to " + argv[i].Substring(0, 3) + " must be an integer.");
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
        public static int[] readInputFtData(string path) {
            String s, temp = "";
            int[] vet = new int[4];
            System.IO.StreamReader file = new System.IO.StreamReader(path);

            s = file.ReadLine(); s = s.Trim(' ');
            for (int i = 0; i < s.Length; i++) {
                if (s[i] == ' ') {
                    break;
                }
                temp = temp + s[i];
            }
            vet[0] = int.Parse(temp); temp = "";

            s = file.ReadLine(); s = s.Trim(' ');
            for (int i = 0; i < s.Length; i++) {
                if (s[i] == ' ') {
                    break;
                }
                temp = temp + s[i];
            }
            vet[1] = int.Parse(temp); temp = "";

            s = file.ReadLine(); s = s.Trim(' ');
            for (int i = 0; i < s.Length; i++) {
                if (s[0] == ' ') {
                    break;
                }
                temp = temp + s[0];
                s = s.Substring(1);
            }
            vet[2] = int.Parse(temp); temp = "";

            s = s.Trim(' ');
            for (int i = 0; i < s.Length; i++) {
                if (s[i] == ' ') {
                    break;
                }
                temp = temp + s[i];
            }
            vet[3] = int.Parse(temp);
            return vet;
        }
    }
}
