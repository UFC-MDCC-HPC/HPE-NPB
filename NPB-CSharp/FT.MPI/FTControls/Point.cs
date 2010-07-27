using System;
using System.Collections.Generic;
using System.Text;

namespace NPB.Ftc {
    class Point {

        public static unsafe void setAddress(double[,,,] s, int i, double[, , ,] d,  int j) {
            if (arrayBoundCheck(s, i) && arrayBoundCheck(d, j)) {
                fixed (double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    p2 += j;
                    p1 += i;
                    *((double*)p2) = *((double*)p1);
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe void setAddress(double[] s, int i, double[] d, int j) {
            if (arrayBoundCheck(s, i) && arrayBoundCheck(d, j)) {
                fixed (double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    p2 += j;
                    p1 += i;
                    *((double*)p2) = *((double*)p1);
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe void setAddress(double[,,,] s, int i, double[,,] d, int j) {
            if (arrayBoundCheck(s, i) && arrayBoundCheck(d, j)) {
                fixed (double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    p2 += j;
                    p1 += i;
                    *((double*)p2) = *((double*)p1);
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe void setAddress(double[,,] s, int i, double[,,,] d, int j) {
            if (arrayBoundCheck(s, i) && arrayBoundCheck(d, j)) {
                fixed (double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    p2 += j;
                    p1 += i;
                    *((double*)p2) = *((double*)p1);
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe void setValue(double[, , ,] v, int i, double value) {
            if (arrayBoundCheck(v,i)) {
                fixed (double* pv = v) {
                    double* p = pv;
                    p += i;
                    *((double*)p) = value;
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe void setRealAndImagToVetsD1(double[, , ,] s, double[] r, double[] i) {
            int size1 = s.Length; 
            int size2 = r.Length; 
            int size3 = i.Length; 
            int size = size1/2;
            if (((size) == size2) && ((size) == size3)) {
                fixed (double* _S = s, _R = r, _I = i) {
                    double* p1 = _S; 
                    double* p2 = _R;
                    double* p3 = _S; 
                    double* p4 = _I;
                    for (int n = -1; n < size; n++) {
                        *((double*)p2) = *((double*)p1);
                        *((double*)p4) = *((double*)p3);
                        if (n != -1) {
                            p2 += 1; p1 += 2;
                            p4 += 1; p3 += 2;
                        }
                        else p3 += 1;
                    }
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe double getValue(double[, , ,] s, int i) {
            if (arrayBoundCheck(s, i)) {
                fixed (double* ps = s) {
                    double* p1 = ps;
                    p1 += i;
                    return *((double*)p1);
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe double[] getImagD1(double[, , ,] s) {
            int size = s.Length;
            double[] d = new double[size/2];
            fixed (double* ps = s, pd = d) {
                double* p1 = ps;
                double* p2 = pd;
                for (int n = -1; n < size / 2; n++) {
                    *((double*)p2) = *((double*)p1);
                    if (n != -1) {
                        p2 += 1;
                        p1 += 2;
                    }
                    else p1 += 1;
                }
            }
            return d;
        }

        public static unsafe double[] getRealD1(double[, , ,] s) {
            int size = s.Length;
            double[] d = new double[size / 2];
            fixed (double* ps = s, pd = d) {
                double* p1 = ps;
                double* p2 = pd;
                for (int n = 0; n < size / 2; n++) {
                    *((double*)p2) = *((double*)p1);
                    p2 += 1;
                    p1 += 2;
                }
            }
            return d;
        }
        
        public static unsafe double[] getD1(double[,,,] s) {
            int size = s.Length;
            double[] d = new double[size];
            fixed (double* ps = s, pd = d) {
                double* p1 = ps;
                double* p2 = pd;
                for (int n = 0; n < size / 2; n++) {
                    *((decimal*)p2) = *((decimal*)p1);
                    p2 += 2;
                    p1 += 2;
                }
            }
            return d;
        }

        public static unsafe void setVetor(double[] s, double[,,,] d) {
            int size = s.Length;
            if (size == d.Length) {
                fixed (double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    for (int n = 0; n < size / 2; n++) {
                        *((decimal*)p2) = *((decimal*)p1);
                        p2 += 2;
                        p1 += 2;
                    }
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static unsafe void setVetor(double[,,,] s, double[] d) {
            int size = s.Length;
            if (size == d.Length) {
                fixed (double* ps = s, pd = d) {
                    double* p1 = ps;
                    double* p2 = pd;
                    for (int n = 0; n < size/2; n++) {
                        *((decimal*)p2) = *((decimal*)p1);
                        p2 += 2;
                        p1 += 2;
                    }
                }
            }
            else { throw new IndexOutOfRangeException(); }
        }

        public static bool arrayBoundCheck(double[,,,] v, int i) {
            return ((i < v.Length) && (i >= 0)) ? true : false;
        }
        public static bool arrayBoundCheck(double[,,] v, int i) {
            return ((i < v.Length) && (i >= 0)) ? true : false;
        }
        public static bool arrayBoundCheck(double[,] v, int i) {
            return ((i < v.Length) && (i >= 0)) ? true : false;
        }
        public static bool arrayBoundCheck(double[] v, int i) {
            return ((i < v.Length) && (i >= 0)) ? true : false;
        }

    }
}
