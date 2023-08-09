/*
The C# example to call the if97 shared library on Windows x64 

>csc -out:demo.exe demo.cs /platform:"x64"

Author:   Cheng Maohua
Email:    cmh@seu.edu.cn
*/

using System;
using System.Runtime.InteropServices;

public class if97
    {
        [DllImport("../target/release/if97", CallingConvention = CallingConvention.StdCall)]
        public static extern double pt(double p, double t, int o_id);
    
}

namespace demo_if97
{
    class demo_if97
    {

        static void Main(string[] args)
        {
            double p = 16.13;
            double t = 535.0;
            double h, s, v;
            h = if97.pt(p, t, 4);
            s = if97.pt(p, t, 5);
            v = if97.pt(p, t, 3);
            Console.WriteLine("(p,t) h,s,v {0 :.00} {1:.0} {2:.000} {3:.000} {4:.000}", p, t, h, s, v);

          }
    }
}
