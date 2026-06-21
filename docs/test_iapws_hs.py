from iapws import IAPWS97
steam = IAPWS97(h=1500, s= 4.0)           
print(steam.P)
steam = IAPWS97(h=1800, s=5.3)          
print(steam.P)