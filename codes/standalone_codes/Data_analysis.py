'''
This code was used to generate Regina data. Not used actively.
'''

import numpy as np
import matplotlib.pyplot as plt
import xlrd
import csv

loc=('E:/Python Projects/Sewage Network/Regent Data (selected).xlsx')

# To open Workbook
wb = xlrd.open_workbook(loc)
sheet = wb.sheet_by_index(2)


#line_id=
line_id=[]
line_type={}
for i in range(1,sheet.nrows):
    line_id.append(sheet.cell_value(i,0))

xs={}
ys={}
xe={}
ye={}
for i in line_id:
    j=sheet.row_values(int(i))
    xs[i]=j[1]/10000
    ys[i]=j[2]/10000
    xe[i]=j[3]/10000
    ye[i]=j[4]/10000
    line_type[i]=int(j[5])

# for i in line_id[0:20]:
#     plt.plot([xs[i],xe[i]],[ys[i],ye[i]],"r")
#     print(i)

ax=52.6
bx=52.7
ay=558.6
by=558.71

selected=[]
for i in line_id:
    if xs[i]>=ax and xs[i]<=bx and xe[i]>=ax and xe[i]<=bx:
        if ys[i] >= ay and ys[i] <= by and ye[i] >= ay and ye[i] <= by:
            selected.append(i)

for i in selected:
    if line_type[i]==1:
        plt.plot([xs[i],xe[i]],[ys[i],ye[i]],"r",linewidth=0.5)
    if line_type[i]==2:
        plt.plot([xs[i], xe[i]], [ys[i], ye[i]], "b", linewidth=0.5)
    print(i)

#plt.show()

sheet2=wb.sheet_by_index(3)
#line_id=
pid=[]
for i in range(1,sheet2.nrows):
    pid.append(sheet2.cell_value(i,0))

px={}
py={}
for i in pid:
    j=sheet2.row_values(int(i))
    px[i]=j[1]/10000
    py[i]=j[2]/10000

# for i in line_id[0:20]:
#     plt.plot([xs[i],xe[i]],[ys[i],ye[i]],"r")
#     print(i)

selected=[]
for i in pid:
    if px[i]>=ax and px[i]<=bx:
        if py[i] >= ay and py[i] <= by:
            selected.append(i)

for i in selected:
    plt.plot(px[i],py[i],"go",markersize=1)
    print(i)

plt.show()
