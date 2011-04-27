using System;
using MPI;

namespace NPB{
    public class Timer{
        public static int max_counters = 64;
        double[] start_time = new double[max_counters];
        double[] elapsed_time = new double[max_counters];
        double[] total_time = new double[max_counters];

        public Timer(){
            for (int i = 0; i < max_counters; i++){
                start_time[i] = 0;
                elapsed_time[i] = 0;
                total_time[i] = 0;
            }
        }

        public void start(int n){
            start_time[n] = MPI.Unsafe.MPI_Wtime();

        }

        public void stop(int n){
            double now;
            now = MPI.Unsafe.MPI_Wtime();
            elapsed_time[n] = now - start_time[n];
            total_time[n] += elapsed_time[n];
        }

        public double readTimer(int n){
            return total_time[n];
        }

        public void resetTimer(int n){
            elapsed_time[n] = start_time[n] = total_time[n] = 0;
        }

        public void resetAllTimers(){
            for (int i = 0; i < max_counters; i++) resetTimer(i);
        }
    }
}