
/*
   The Java example with Java Native Access(JNA) to call the shared library

     https://github.com/java-native-access/jna
     
   Download JNA:
      https://repo1.maven.org/maven2/net/java/dev/jna/jna/5.13.0/jna-5.13.0.jar

    javac -cp jna.jar demo_jni_cdecl.java
    java -cp .;jna.jar demo_jni_cdecl
 */
import com.sun.jna.Native;
import com.sun.jna.Library;

interface seuif97 extends Library {

    seuif97 lib = (seuif97) Native.load("../target/release/seuif97", seuif97.class);
    // universal functions with o_id parameter
    public double pt(double p, double t, int o_id);
    // direct property functions 
    public double pt2s(double p, double t); 
}

public class demo_jni_cdecl {

    public static void main(String[] args) {
        double p = 16.0;
        double t = 540.0;
        double h,s;
        //  the Universal Function
        h = seuif97.lib.pt(p, t, 4);
        // the Direct Property Function  
        s = seuif97.lib.pt2s(p, t);
        System.out.printf("(p,t)->h,s: (%.1f %.1f) h: %.2f, s: %.2f\n", p, t, h, s);
    }
}