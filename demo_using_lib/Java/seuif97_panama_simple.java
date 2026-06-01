/**
 * SEUIF97 Panama FFI Demo (No jextract required)
 * 
 * Call Rust library directly using Java 22+ Project Panama FFI
 * No binding code generation required
 * 
 * Usage:
 * 
 * 1. Build Rust library:
 *    cargo build -r --features cdecl --target x86_64-pc-windows-msvc
 * 
 * 2. Compile:
 *    javac seuif97_panama_simple.java
 * 
 * 3. Run:
 *    java seuif97_panama_simple
 */

import java.lang.foreign.*;
import java.lang.invoke.MethodHandle;

public class seuif97_panama_simple {
    
    // Load library in static initializer
    static {
        // Try to load seuif97 library
        try {
            System.loadLibrary("seuif97");
            System.out.println("[OK] Successfully loaded seuif97 library");
        } catch (UnsatisfiedLinkError e) {
            System.err.println("[ERROR] Failed to load seuif97 library: " + e.getMessage());
            System.err.println("        Please build with: cargo build -r --features cdecl --target x86_64-pc-windows-msvc");
            System.exit(1);
        }
    }
    
    // Define function signatures
    private static final Linker LINKER = Linker.nativeLinker();
    private static final SymbolLookup LOOKUP = SymbolLookup.loaderLookup();
    
    // Function handles - manually defined
    private static final MethodHandle PT;
    private static final MethodHandle PH;
    private static final MethodHandle PS;
    
    static {
        try {
            // Find function addresses - Java 22+ uses MemorySegment
            MemorySegment ptSeg = LOOKUP.find("pt").orElseThrow(() -> 
                new RuntimeException("Function 'pt' not found in library"));
            MemorySegment phSeg = LOOKUP.find("ph").orElseThrow(() -> 
                new RuntimeException("Function 'ph' not found in library"));
            MemorySegment psSeg = LOOKUP.find("ps").orElseThrow(() -> 
                new RuntimeException("Function 'ps' not found in library"));
            
            // Define function signature: double pt(double, double, int)
            FunctionDescriptor ptDesc = FunctionDescriptor.of(
                ValueLayout.JAVA_DOUBLE,    // return: double
                ValueLayout.JAVA_DOUBLE,    // param1: p (MPa)
                ValueLayout.JAVA_DOUBLE,    // param2: t (°C)
                ValueLayout.JAVA_INT        // param3: o_id
            );
            
            // Create MethodHandle
            PT = LINKER.downcallHandle(ptSeg, ptDesc);
            PH = LINKER.downcallHandle(phSeg, ptDesc);
            PS = LINKER.downcallHandle(psSeg, ptDesc);
            
        } catch (Exception e) {
            throw new RuntimeException("Failed to initialize native bindings", e);
        }
    }
    
    // Wrapper methods
    public static double pt(double p, double t, int o_id) throws Throwable {
        return (double) PT.invokeExact(p, t, o_id);
    }
    
    public static double ph(double p, double h, int o_id) throws Throwable {
        return (double) PH.invokeExact(p, h, o_id);
    }
    
    public static double ps(double p, double s, int o_id) throws Throwable {
        return (double) PS.invokeExact(p, s, o_id);
    }
    
    public static void main(String[] args) throws Throwable {
        System.out.println("\n=== SEUIF97 Panama FFI Demo (Simple) ===");
        System.out.println();
        
        double p = 16.0;   // MPa
        double t = 540.0;  // °C
        
        System.out.println("Input: p = " + p + " MPa, t = " + t + " °C");
        System.out.println();
        
        // o_id: 4=enthalpy, 5=entropy, 3=volume, 1=temperature
        double h = pt(p, t, 4);
        double s = pt(p, t, 5);
        double v = pt(p, t, 3);
        
        System.out.println("Results from (p, t):");
        System.out.printf("  h = %.3f kJ/kg%n", h);
        System.out.printf("  s = %.5f kJ/(kg·K)%n", s);
        System.out.printf("  v = %.6f m³/kg%n", v);
        
        // Test with (p, h) input
        double p2 = 3.0;
        double h2 = 3000.0;
        double t2 = ph(p2, h2, 1);  // o_id=1 returns temperature
        
        System.out.println();
        System.out.println("Results from (p, h):");
        System.out.printf("  p = %.1f MPa, h = %.0f kJ/kg -> t = %.2f °C%n", p2, h2, t2);
        
        System.out.println();
        System.out.println("=== Demo Complete ===");
    }
}

