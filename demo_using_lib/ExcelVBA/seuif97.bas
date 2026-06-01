Attribute VB_Name = "seuif97"
'  IF97
'   The VBA Module to access seuif97.dll
'
'   License: this code is in the public domain
'
'   Author:   Cheng Maohua
'   Email:    cmh@ .edu.cn
'
' Last modified: 2016.4.20
'
Option Explicit

' API declarations
#If Win64 Then

    Declare PtrSafe Function pt Lib "seuif97" (ByVal p As Double, ByVal t As Double, ByVal OutWhat As Integer) As Double
    Declare PtrSafe Function ph Lib "seuif97" (ByVal p As Double, ByVal h As Double, ByVal OutWhat As Integer) As Double
    Declare PtrSafe Function ps Lib "seuif97" (ByVal p As Double, ByVal s As Double, ByVal OutWhat As Integer) As Double
    Declare PtrSafe Function pv Lib "seuif97" (ByVal p As Double, ByVal v As Double, ByVal OutWhat As Integer) As Double
    
    Declare PtrSafe Function th Lib "seuif97" (ByVal t As Double, ByVal h As Double, ByVal OutWhat As Integer) As Double
    Declare PtrSafe Function ts Lib "seuif97" (ByVal t As Double, ByVal s As Double, ByVal OutWhat As Integer) As Double
    Declare PtrSafe Function tv Lib "seuif97" (ByVal t As Double, ByVal v As Double, ByVal OutWhat As Integer) As Double
    
    Declare PtrSafe Function px Lib "seuif97" (ByVal p As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    Declare PtrSafe Function tx Lib "seuif97" (ByVal t As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    
    Declare PtrSafe Function hx Lib "seuif97" (ByVal h As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    Declare PtrSafe Function sx Lib "seuif97" (ByVal s As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    
    Declare PtrSafe Function hs Lib "seuif97" (ByVal h As Double, ByVal s As Double, ByVal OutWhat As Integer) As Double
     
    

#Else

    Declare Function pt Lib "seuif97" (ByVal p As Double, ByVal t As Double, ByVal OutWhat As Integer) As Double
    Declare Function ph Lib "seuif97" (ByVal p As Double, ByVal h As Double, ByVal OutWhat As Integer) As Double
    Declare Function ps Lib "seuif97" (ByVal p As Double, ByVal s As Double, ByVal OutWhat As Integer) As Double
    Declare Function pv Lib "seuif97" (ByVal p As Double, ByVal v As Double, ByVal OutWhat As Integer) As Double
    
    Declare Function th Lib "seuif97" (ByVal t As Double, ByVal h As Double, ByVal OutWhat As Integer) As Double
    Declare Function ts Lib "seuif97" (ByVal t As Double, ByVal s As Double, ByVal OutWhat As Integer) As Double
    Declare Function tv Lib "seuif97" (ByVal t As Double, ByVal v As Double, ByVal OutWhat As Integer) As Double
    
    Declare Function px Lib "seuif97" (ByVal p As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    Declare Function tx Lib "seuif97" (ByVal t As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    
    Declare Function hx Lib "seuif97" (ByVal h As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    Declare Function sx Lib "seuif97" (ByVal s As Double, ByVal x As Double, ByVal OutWhat As Integer) As Double
    
    Declare Function hs Lib "seuif97" (ByVal h As Double, ByVal s As Double, ByVal OutWhat As Integer) As Double
    
    
#End If

Public Function if97pt(ByVal p As Double, ByVal t As Double, ByVal wp As Integer) As Double
   if97pt = pt(p, t, wp)
End Function

Public Function if97ph(ByVal p As Double, ByVal h As Double, ByVal wp As Integer) As Double
   if97ph = ph(p, h, wp)
End Function

Public Function if97ps(ByVal p As Double, ByVal s As Double, ByVal wp As Integer) As Double
   if97ps = ps(p, s, wp)
End Function

Public Function if97pv(ByVal p As Double, ByVal v As Double, ByVal wp As Integer) As Double
   if97pv = pv(p, v, wp)
End Function

Public Function if97th(ByVal t As Double, ByVal h As Double, ByVal wp As Integer) As Double
   if97th = th(t, h, wp)
End Function

Public Function if97ts(ByVal t As Double, ByVal s As Double, ByVal wp As Integer) As Double
   if97ts = ts(t, s, wp)
End Function

Public Function if97tv(ByVal t As Double, ByVal v As Double, ByVal wp As Integer) As Double
   if97tv = tv(t, v, wp)
End Function

Public Function if97hs(ByVal h As Double, ByVal s As Double, ByVal wp As Integer) As Double
   if97hs = hs(h, s, wp)
End Function

Public Function if97px(ByVal p As Double, ByVal x As Double, ByVal wp As Integer) As Double
   if97px = px(p, x, wp)
End Function

Public Function if97tx(ByVal t As Double, ByVal x As Double, ByVal wp As Integer) As Double
   if97tx = tx(t, x, wp)
End Function

Public Function if97hx(ByVal h As Double, ByVal x As Double, ByVal wp As Integer) As Double
   if97hx = hx(h, x, wp)
End Function

Public Function if97sx(ByVal s As Double, ByVal x As Double, ByVal wp As Integer) As Double
   if97sx = sx(t, x, wp)
End Function

