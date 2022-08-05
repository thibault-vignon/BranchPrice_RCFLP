\ENCODING=ISO-8859-1
\Problem name: RCFLP

Minimize
 obj: V_0 + 2 V_1 + 2.5 V_2 + 1.5 V_3 + 1.5 V_4 + 0.5 V_5
Subject To
 Convexity(0)#0:   V_0 + V_2 + V_4  = 1
 Convexity(1)#1:   V_1 + V_3 + V_5  = 1
 Capacity(0)#2:    - 40 V_0 - 30 V_1 + 100 V_2 + 100 V_3 + 100 V_4 + 100 V_5
                   >= 0
 Capacity(1)#3:    50 V_0 + 50 V_1 + 10 V_2 + 20 V_3 - 40 V_4 - 30 V_5 >= 0
 _Equality(0,0)#4: V_2 - V_3 + V_4 - V_5  = 0
 _Equality(1,0)#5: V_0 - V_1 + V_2 - V_3  = 0
End
