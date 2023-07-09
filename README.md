# IF97

[![DOI](https://zenodo.org/badge/110833324.svg)](https://zenodo.org/badge/latestdoi/110833324)

IF97 is the high-speed shared library of IAPWS-IF97. It is suitable for computation-intensive calculations，such as heat cycle calculations, simulations of non-stationary processes, real-time process monitoring and optimizations.   
 
Through the high-speed library, the results of the IAPWS-IF97 are accurately produced at about 3x speed-up compared to the repeated squaring method for fast computation of integer powers, 10x speed-up compared to  the `math.pow()` of the C standard library.   
The following input pairs are implemented: 

```
(p,t) (p,h) (p,s) (p,x) 

(t,x) 

(h,s)  
```

### Usage

The functions are provided in the if97 pacakge:

* ??(in1,in2,propertyID), , e.g:`h=pt(p,t,4)`, the propertyID h is 4 
   * first,second input parameters: the input properties(double)
   * third input parameter: the propertyID of the calculated property(int, 0-29), see Properties in libseuif97
   * the return: the calculated property value(double)

## Properties in if97

| Properties                            |    Unit     | symbol | propertyID |
| ------------------------------------- | :---------: | -----: | ---------: |
| Pressure                              |     MPa     |      p |          0 |
| Temperature                           |     °C      |      t |          1 |
| Density                               |   kg/m^3    |      d |          2 |
| Specific Volume                       |   m^3/kg    |      v |          3 |
| Specific enthalpy                     |    kJ/kg    |      h |          4 |
| Specific entropy                      |  kJ/(kg·K)  |      s |          5 |
| Specific exergy                       |    kJ/kg    |      e |          6 |
| Specific internal energy              |    kJ/kg    |      u |          7 |
| Specific isobaric heat capacity       |  kJ/(kg·K)  |     cp |          8 |
| Specific isochoric heat capacity      |  kJ/(kg·K)  |     cv |          9 |
| Speed of sound                        |     m/s     |      w |         10 |
| Isentropic exponent                   |             |     ks |         11 |
| Specific Helmholtz free energy        |    kJ/kg    |      f |         12 |
| Specific Gibbs free energy            |    kJ/kg    |      g |         13 |
| Compressibility factor                |             |      z |         14 |
| Steam quality                         |             |      x |         15 |
| Region                                |             |      r |         16 |
| Isobaric volume expansion coefficient |     1/K     |     ec |         17 |
| Isothermal compressibility            |    1/MPa    |     kt |         18 |
| Partial derivative (dV/dT)p           |  m3/(kg·K)  |   dvdt |         19 |
| Partial derivative (dV/dP)T           | m3/(kg·MPa) |   dvdp |         20 |
| Partial derivative (dP/dT)v           |    MPa/K    |   dpdt |         21 |
| Isothermal Joule-Thomson coefficient  | kJ/(kg·MPa) |   iJTC |         22 |
| Joule-Thomson coefficient             |    K/MPa    |    JTC |         23 |
| Dynamic viscosity                     |  kg/(m·s)   |     dv |         24 |
| Kinematic viscosity                   |    m^2/s    |     kv |         25 |
| Thermal conductivity                  |   W/(m.K)   |     tc |         26 |
| Thermal diffusivity                   |   um^2/s    |     td |         27 |
| Prandtl number                        |             |     pr |         28 |
| Surface tension                       |    mN/m     |     st |         29 |

## Cite as

* Wang Pei-hong, Jia Jun-ying, Cheng Mao-hua. General Calculating Models of Water and Steam Properties(IAPWS-IF97), Power Engineering, 21 (6), 2001, pp. 1564-1567. (in Chinese)

* Cheng Maohua. SEUIF97: the high-speed shared library of IAPWS-IF97(1.0.1). Zenodo. https://doi.org/10.5281/zenodo.4586961