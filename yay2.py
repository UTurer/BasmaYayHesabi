# -*- coding: utf-8 -*-
"""
Basma yayı hesabı


"""

import numpy as np
import math

print("\n"*80)

dspan = np.linspace(0.1,1,10)*1e-3
Dspan = np.linspace(1,10,10)*1e-3
G = 69e9 #Pa, AISI302 stainless steel
# F_max = (200/1000*9.81)*(100/80) # %80 sıkışınca 200 gr kuvvet verecek.
# defl = 100/3*1e-3; # %20-%80 arası sıkışınca 20mm hareket edecek.
nsspan = np.linspace(1,1.5,5)
F1 = 300/1000*9.81
F2 = 400/1000*9.81
k = (5/1000)*9.81/1e-3

NN = len(dspan)*len(Dspan)*len(nsspan)

d = np.zeros(NN)
D = np.zeros(NN)
ns = np.zeros(NN)

for i in np.arange(0,NN):
    i1=i%len(dspan)
    i2=math.floor(i/len(dspan))%len(Dspan)
    i3=math.floor(i/len(dspan)/len(Dspan))%len(nsspan)
    d[i]=dspan[i1]
    D[i]=Dspan[i2]
    ns[i]=nsspan[i3]
    
del i1
del i2
del i3

# C, spring factor
C = D/d

# W, Wahl factor, Shigley 7th Ed, sayfa 511, Eqn 10-5
W = (4*C-1)/(4*C-4)+0.615/C

# Ssy, Pa, yield strength, Shigley 7th Ed. sayfa 517 Table 10-4 ve altındaki paragraf
Ssy1 = 1867/(d*1000)**0.146*0.45
Ssy2 = 2065/(d*1000)**0.263*0.45
Ssy3 = 2911/(d*1000)**0.478*0.45
Ssy = Ssy1
Ssy[(d>=2.5e-3) & (d<5e-3)] = Ssy2[(d>=2.5e-3) & (d<5e-3)]
Ssy[d>5e-3] = Ssy3[d>5e-3] # Mpa
Ssy = Ssy*1e6 # Pa

# ns, safety factor
Tau_max = Ssy/ns

# Tau_max [Pa], Max shear stress, Shigley 7th Ed. sayfa 511 Eqn 10-3 ve sayfa 512 Eqn 10-7 ve altındaki paragraf
F_max = Tau_max*np.pi*(d**3)/8/W/D

# k [N/m], stiffness
defl = F_max/k

# Na, number of active coils, Shigley 7th Ed., sayfa 512, Eq10-9
Na = G*(d**4/D**3)/8/k

# Nt, number of total coils, Shigley 7th Ed., sayfa 513, Table 10-1 (Squared and Ground)
Nt = Na+2

# L_solid [m], tamamen sıkışmış halde yay uzunluğu, Shigley 7th Ed., sayfa 513, Table 10-1 (Squared and Ground)
L_solid = Nt*d

# L_free [m], sıkışmış halde yay uzunluğu, Shigley 7th Ed., sayfa 513, Table 10-1 (Squared and Ground)
L_free = defl+L_solid

# p [m], pitch, Shigley 7th Ed., sayfa 513, Table 10-1 (Squared and Ground)
p = (L_free-2*d)/Na

# CONDITIONLARI Belirleme
cond1 = ns>=0.8
cond2 = (C>=4) & (C<=12)
cond3 = (Na>=3) & (Na<=15);
cond4 = defl>=(20*1e-3*100/60)
cond5 = F_max>(400/1000*9.81)
cond6 = (Na<=30)
#cond = (cond1) & (cond2) & (cond3)
cond = (cond2) & (cond4) & (cond5)

# Uygun çözüm parametrelerini belirleme
d_ok = d[cond]
D_ok = D[cond]
Nt_ok = Nt[cond]
Na_ok = Na[cond]
L_solid_ok = L_solid[cond]
L_free_ok = L_free[cond]
ns_ok = ns[cond]
p_ok = p[cond]
C_ok = C[cond]
defl_ok = defl[cond]
F_max_ok = F_max[cond]

# COST belirleme
W1 = 1/5
W2 = 1
W3 = 2
COST_W1 = (C_ok-8)**2*W1
COST_W2 = (L_solid_ok*1000)*W2
COST_W3 = 1/(ns_ok-0.8)*W3
COST = COST_W1+COST_W2+COST_W3

from prettytable import PrettyTable
tbl = PrettyTable()
tbl.field_names = ["d (mm)","D (mm)","Na","L_solid (mm)","L_free (mm)","Defl (mm)","F (gr)","pitch (mm)","ns","C","COST_W1","COST_W2","COST_W3","COST"]
for i in np.arange(0,len(d_ok)):
    tbl.add_row([d_ok[i]*1000,D_ok[i]*1000,Na_ok[i],L_solid_ok[i]*1000,L_free_ok[i]*1000,defl_ok[i]*1000,F_max_ok[i]/9.81*1000,p_ok[i]*1000,ns_ok[i],C_ok[i],COST_W1[i],COST_W2[i],COST_W3[i],COST[i]])

tbl.sortby = "Na"
tbl.float_format=".2"
print(tbl.get_string(start=0,end=9))


Na=50.88
p=2.31
d=0.5
D=6

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
t=np.linspace(0,Na*2*np.pi,1000)
X1=p/2/np.pi*t
X2=X1+d
to=d*2*np.pi/p
Y1=D/2*np.sin(t)
Y2=D/2*np.sin(t)
ax.plot(X1,Y1,"k")
ax.plot(X2,Y2,"k")
# plt.show()
plt.grid()
# plt.ylim(0,100)
X3=d/2/np.pi*t
X4=X3+d
ax.plot(X3,Y1,"r")
ax.plot(X4,Y2,"r")


# Sonuçları sıralama
# [SORTED,IND_SORT] = sort(COST,2,'ascend');

# Sonuçları yazdırma
#T = table(d_ok(IND_SORT)'*1000,D_ok(IND_SORT)'*1000,Na_ok(IND_SORT)',L_solid_ok(IND_SORT)'*1000,p_ok(IND_SORT)'*1000,ns_ok(IND_SORT)',C_ok(IND_SORT)',COST_W1(IND_SORT)',COST_W2(IND_SORT)',COST(IND_SORT)','VariableNames',{'d_mm','D_mm','Na','L_solid_mm','pitch_mm','ns','C','COST_W1','COST_W2','COST'})