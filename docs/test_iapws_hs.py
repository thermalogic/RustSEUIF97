from iapws import IAPWS97
steam = IAPWS97(h=1500, s= 4.0)           
print(steam.P,steam.x)
steam = IAPWS97(h=1800, s=5.3)          
print(steam.P,steam.x)
s13 = IAPWS97(T=623.15, P= 100.0)           
print("s13",s13.s)
s13s = IAPWS97(T=623.15, P=16.5291642526045)           
print("s13",s13s.s)

sl = IAPWS97(T=273.15, P=0.000611212677444)           
print("sl",sl.h,sl.s)
