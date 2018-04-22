# check clustered barcodes against actual for each sample

import math

def jaccard_similarity(x,y):
    intersect = len(set.intersection(*[set(x), set(y)]))
    union = len(set.union(*[set(x), set(y)]))
    return intersect/float(union)

direct = '/Users/Caitlin/Documents/Bioinformatics/Summer_Project_2017/analysis/data'

filename_A = '{}/pbmc_donor_A/bc_100k_sorted.txt'.format(direct)
filename_B = '{}/pbmc_donor_B/bc_100k_sorted.txt'.format(direct)
filename_C = '{}/pbmc_donor_C/bc_100k_sorted.txt'.format(direct)

filenames_m_A = ['{}/pbmc_A-C_mix/model_A{}_ABC_100k.txt'.format(direct,_) for _ in range(1,4)]
filenames_m_B = ['{}/pbmc_A-C_mix/model_B{}_ABC_100k.txt'.format(direct,_) for _ in range(1,4)]
filenames_m_C = ['{}/pbmc_A-C_mix/model_C{}_ABC_100k.txt'.format(direct,_) for _ in range(1,4)]

donor_A = []
donor_B = []
donor_C = []
model_A = []
model_B = []
model_C = []

for barcode in open(filename_A, 'r'):
    donor_A.append(barcode.strip())
print('donor_A len: ', len(donor_A))

for barcode in open(filename_B, 'r'):
    donor_B.append(barcode.strip())
print('donor_C len: ', len(donor_B))

for barcode in open(filename_C, 'r'):
    donor_C.append(barcode.strip())
print('donor_C len: ', len(donor_C))


for filename in filenames_m_A:
    for barcode in open(filename, 'r'):
        if barcode.strip() not in model_A:
            model_A.append(barcode.strip())
print('model_A len: ', len(model_A))


for filename in filenames_m_B:
    for barcode in open(filename, 'r'):
        if barcode.strip() not in model_B:
            model_B.append(barcode.strip())
print('model_B len: ', len(model_B))

for filename in filenames_m_C:
    for barcode in open(filename, 'r'):
        if barcode.strip() not in model_C:
            model_C.append(barcode.strip())
print('model_C len: ', len(model_C))


donA_modA = jaccard_similarity(donor_A, model_A)
donA_modB = jaccard_similarity(donor_A, model_B)
donA_modC = jaccard_similarity(donor_A, model_C)
donB_modA = jaccard_similarity(donor_B, model_A)
donB_modB = jaccard_similarity(donor_B, model_B)
donB_modC = jaccard_similarity(donor_B, model_C)
donC_modA = jaccard_similarity(donor_C, model_A)
donC_modB = jaccard_similarity(donor_C, model_B)
donC_modC = jaccard_similarity(donor_C, model_C)

donA_donB = jaccard_similarity(donor_A, donor_C)
donA_donC = jaccard_similarity(donor_A, donor_C)
donB_donC = jaccard_similarity(donor_B, donor_C)

modA_modB = jaccard_similarity(model_A, model_B)
modA_modC = jaccard_similarity(model_A, model_C)
modB_modC = jaccard_similarity(model_B, model_C)


print("similarity donor A to model A: ", donA_modA)
print("similarity donor A to model B: ", donA_modB)
print("similarity donor A to model C: ", donA_modC)

print("similarity donor B to model A: ", donB_modA)
print("similarity donor B to model B: ", donB_modB)
print("similarity donor B to model C: ", donB_modC)

print("similarity donor C to model A: ", donC_modA)
print("similarity donor C to model B: ", donC_modB)
print("similarity donor C to model C: ", donC_modC)

print("similarity donor A to donor B:", donA_donB)
print("similarity donor A to donor C: ", donA_donC)
print("similarity donor B to donor C: ", donB_donC)
print("similarity model A to model B: ", modA_modB)
print("similarity model A to model C: ", modA_modC)
print("similarity model B to model C: ", modB_modC)


a = 0
b = 0
c = 0
x = 0
denom = len(model_A)
for barcode in model_A:
    if barcode in donor_A:
        a += 1
    if barcode in donor_B:
        b += 1
    if barcode in donor_C:
        c += 1
    if (barcode not in donor_A) and (barcode not in donor_B) and (barcode not in donor_C):
        x += 1
print("Barcodes in model A found in donor A = {:.2f}%, donor B = {:.2f%}, "
      "donor C = {:.2f}%, none = {:.2f}%".format((a/denom)*100, (b/denom)*100,
                                                 (c/denom)*100, (x/denom)*100))

a = 0
b = 0
c = 0
x = 0
denom = len(model_B)
for barcode in model_B:
    if barcode in donor_A:
        a += 1
    if barcode in donor_B:
        b += 1
    if barcode in donor_C:
        c += 1
    if (barcode not in donor_A) and (barcode not in donor_B) and (barcode not in donor_C):
        x += 1
print("Barcodes in model B found in donor A = {:.2f}%, donor B = {:.2f%}, "
      "donor C = {:.2f}%, none = {:.2f}%".format((a/denom)*100, (b/denom)*100,
                                                 (c/denom)*100, (x/denom)*100))


a = 0
b = 0
c = 0
x = 0
denom = len(model_C)
for barcode in model_C:
    if barcode in donor_A:
        a += 1
    if barcode in donor_B:
        b += 1
    if barcode in donor_C:
        c += 1
    if (barcode not in donor_A) and (barcode not in donor_B) and (barcode not in donor_C):
        x += 1
print("Barcodes in model C found in donor A = {:.2f}%, donor B = {:.2f%}, "
      "donor C = {:.2f}%, none = {:.2f}%".format((a/denom)*100, (b/denom)*100,
                                                 (c/denom)*100, (x/denom)*100))

y = 0
n = 0
denom = len(model_A)
for barcode in model_A:
    if barcode in model_B:
        y += 1
    else:
        n += 1
print("Barcode in model A also in model B: {:.2f}%".format((y/denom)*100))

y = 0
n = 0
denom = len(model_A)
for barcode in model_A:
    if barcode in model_C:
        y += 1
    else:
        n += 1
print("Barcode in model A also in model C: {:.2f}%".format((y/denom)*100))

y = 0
n = 0
denom = len(model_B)
for barcode in model_B:
    if barcode in model_C:
        y += 1
    else:
        n += 1
print("Barcode in model B also in model C: {:.2f}%".format((y/denom)*100))