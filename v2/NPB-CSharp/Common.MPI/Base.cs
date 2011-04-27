using System;
using NPB;

namespace NPB {
    public class Base {

        public static double mod(double a, double b) { return (a % b); }

        public double min(int n1, int n2) { return n1 < n2 ? n1 : n2; }

        public double max(double n1, double n2) { return n1 > n2 ? n1 : n2; }

        public double pow2(double p) { return p * p; }
    }
}
