/*
The C# example to call the shared library on Windows x64 
  
  winget install Microsoft.DotNet.SDK.8

  Add the path of C# compiler,for example
  
  *  C:\Windows\Microsoft.NET\Framework64\v4.0.30319\ to the system environment variables #Path# 

  csc -out:demo.exe demo.cs /platform:"x64"

Author:   Cheng Maohua
Email:    cmh@seu.edu.cn
*/

using System;
using System.Runtime.InteropServices;

public static class seuif97
    {
        [DllImport("../target/release/seuif97", CallingConvention = CallingConvention.StdCall)]
        public static extern double pt(double p, double t, int o_id);
        [DllImport("../target/release/seuif97", CallingConvention = CallingConvention.StdCall)]
        public static extern double pt2s(double p, double t);
    
}

namespace demo_seuif97
{
    static class demo_seuif97
    {

        static void Main(string[] args)
        {
            double p = 16.13;
            double t = 535.0;
            double h, s;
            // Universal Functions (with o_id parameter))
            h = seuif97.pt(p, t, 4);
            // Direct Property Functions
            s = seuif97.pt2s(p, t);
            Console.WriteLine("(p,t) h,s {0 :.00} {1:.0} {2:.000} {3:.000} ", p, t, h, s);

          }
    }
}
