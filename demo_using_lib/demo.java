/*
  The Java example with Java Native Access(JNA)
     https://github.com/java-native-access/jna
  
  javac -cp jna.jar demo.java
  java -cp .;jna.jar demo
 */
import com.sun.jna.Native;
import com.sun.jna.Library;

interface if97 extends Library { 
  
    if97 lib = (if97)Native.load("../target/release/if97",  if97.class);
    public double pt(double p, double t, int o_id);  
   
}

public class demo {
    
    public static void main(String[] args){
        double p=16.0;
        double t=540.0;
        double h;
        h=if97.lib.pt(p,t,4);
        System.out.printf("(p,t)->h: (%.1f %.1f) h: %.2f\n",p,t,h);
    }
}