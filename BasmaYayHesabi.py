import numpy as np
import math

print("\n"*80)

d_range = np.array(np.arange(1,6.5,0.1))*1e-3
D_mean_range = np.array(np.arange(18,28.5,0.1))*1e-3
G = 69e9 #Pa, AISI302 stainless steel
# F_max = (200/1000*9.81)*(100/80) # %80 sıkışınca 200 gr kuvvet verecek.
# defl = 100/3*1e-3; # %20-%80 arası sıkışınca 20mm hareket edecek.
ns_range = np.array(np.arange(1.2,1.51,0.1))
F1 = 400/1000*9.81
F2 = 1000/1000*9.81
k = (F2-F1)/4e-3
defl1 = F1/k
defl2 = F2/k

N_range = len(d_range)*len(D_mean_range)*len(ns_range)

d = np.zeros(N_range)
D = np.zeros(N_range)
ns = np.zeros(N_range)

for i in np.arange(0,N_range):
    i1=i%len(d_range)
    i2=math.floor(i/len(d_range))%len(D_mean_range)
    i3=math.floor(i/len(d_range)/len(D_mean_range))%len(ns_range)
    d[i]=d_range[i1]
    D[i]=D_mean_range[i2]
    ns[i]=ns_range[i3]
    
del i1
del i2
del i3
del N_range

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
k = F_max/defl

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
cond3 = (Na>=3) & (Na<=15)
cond4 = defl>=0
cond5 = F_max>0
cond = (cond1) & (cond2) & (cond3) & (cond4) & (cond5)
#cond = (cond3) & (cond4) & (cond5)
#cond = (cond2)

# Uygun çözüm parametrelerini belirleme
d_ok = d[cond]
D_mean_ok = D[cond]
D_in_ok = D_mean_ok - d_ok
D_out_ok = D_mean_ok + d_ok
Nt_ok = Nt[cond]
Na_ok = Na[cond]
L_solid_ok = L_solid[cond]
L_free_ok = L_free[cond]
k_ok = k[cond]
ns_ok = ns[cond]
p_ok = p[cond]
C_ok = C[cond]
defl_ok = defl[cond]
F_max_ok = F_max[cond]

# COST belirleme
W1 = 1
W2 = 1
W3 = 2
COST_W1 = (C_ok-8)**2*W1
COST_W2 = (L_solid_ok*1000)*W2
COST_W3 = 1/(ns_ok-0.8)*W3
COST = COST_W1+COST_W2+COST_W3

from prettytable import PrettyTable
tbl = PrettyTable()
tbl.field_names = ["index",
                   "d (mm)",
                   "D_mean (mm)",
                   "D_in (mm)", 
                   "D_out (mm)",
                   "Na",
                   "L_solid (mm)",
                   "L_free (mm)",
                   "Defl (mm)",
                   "F (gr)",
                   "pitch (mm)",
                   "k (gr/mm)",
                   "ns",
                   "C",
                   "COST_W1",
                   "COST_W2",
                   "COST_W3",
                   "COST"]
for i in np.arange(0,len(d_ok)):
    tbl.add_row([i,
                 d_ok[i]*1000,
                 D_mean_ok[i]*1000,
                 D_in_ok[i]*1000,
                 D_out_ok[i]*1000,
                 Na_ok[i],
                 L_solid_ok[i]*1000,
                 L_free_ok[i]*1000,
                 defl_ok[i]*1000,
                 F_max_ok[i]/9.81*1000,
                 p_ok[i]*1000,
                 k_ok[i]/9.81,
                 ns_ok[i],
                 C_ok[i],
                 COST_W1[i],
                 COST_W2[i],
                 COST_W3[i],
                 COST[i]])

tbl.sortby = "L_free (mm)"
tbl.float_format=".2"
print(tbl)

# # Change some column alignments; default was 'c'
# x.align['column_one'] = 'r'
# x.align['col_two'] = 'r'
# x.align['column_3'] = 'l'


# import matplotlib.pyplot as plt
# # plot
# fig, ax = plt.subplots()
# X=np.arange(0.6,2.1,0.01)
# Y=np.exp(-1*X)
# Y2=1/(X-0.8)
# ax.plot(X,Y2)
# # plt.show()
# plt.grid()
# plt.ylim(0,100)

# Sonuçları sıralama
# [SORTED,IND_SORT] = sort(COST,2,'ascend');

# Sonuçları yazdırma
#T = table(d_ok(IND_SORT)'*1000,D_ok(IND_SORT)'*1000,Na_ok(IND_SORT)',L_solid_ok(IND_SORT)'*1000,p_ok(IND_SORT)'*1000,ns_ok(IND_SORT)',C_ok(IND_SORT)',COST_W1(IND_SORT)',COST_W2(IND_SORT)',COST(IND_SORT)','VariableNames',{'d_mm','D_mm','Na','L_solid_mm','pitch_mm','ns','C','COST_W1','COST_W2','COST'})