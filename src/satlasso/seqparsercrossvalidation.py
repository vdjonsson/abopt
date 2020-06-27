import csv 
import numpy 
import random 

fpath =  './neutdata/'
outp = 'neut/'
fitp ='fits/'
pffile =  'NeutSeqData_'
seqfile = 'MatrixSeq'  
neutfile = 'MatrixNeut_' 
indfile = 'Index_'
virusfile = 'Virus_'
aafile = 'AA'  
ss20 = 'SS20'
ss80 = 'SS80'

bnabs = [ 'NIH45-46.csv', 'NIH45-46G54W.csv', '45-46m2.csv', '45-46m7.csv', '45-46m25.csv', '45-46m28.csv' ,'8ANC195.csv', '10-996.csv','10-1074.csv', '2G12.csv', '2F5.csv']
# bnabs = [ '12A12.csv', 'VRC-PG04b.csv', 'PG9.csv','PGT121.csv', 'VRC-PG04.csv', 'NIH45-46G54W.csv']

bnabstr = bnabs[6]

# Matrix sequence file, same for all bnabs
fw = open(fpath + outp + seqfile + bnabstr, 'rU')
fwn = open(fpath + outp+ neutfile + bnabstr , 'rU')
#fwv = open(fpath + outp+ virusfile + bnabstr , 'rU')

fw80 = open(fpath + outp + ss80 + seqfile + bnabstr, 'w+')
fwn80 = open(fpath + outp + ss80 + neutfile + bnabstr  , 'w+')
fwi80 = open(fpath + outp + ss80 + indfile + bnabstr, 'w+')


fw20 = open(fpath + outp + ss20 + seqfile + bnabstr, 'w+')
fwn20 = open(fpath + outp + ss20 + neutfile + bnabstr, 'w+')
fwi20 = open(fpath + outp + ss20 + indfile + bnabstr, 'w+')

fwaa = open(fpath + fitp + aafile + bnabstr, 'w+')
fwaa80 = open(fpath + fitp + ss80 + aafile + bnabstr, 'w+')
fwaa20 = open(fpath + fitp + ss20 + aafile + bnabstr, 'w+')

fwreader = csv.reader(fw)
fwnreader = csv.reader(fwn)
#fwvreader = csv.reader(fwv)

i = 1

aalist = ['A', 'R', 'N', 'D','C','Q','E','G','H', 'I','L','K','M', 'F','P','S', 'T', 'W', 'Y' ,'V']

hivseqmatrix = []
hivseqneut = [] 

m=1

# Get all viruses and neutralization values into a list 
for row in fwreader:	
	hivseqmatrix.append(row)
	
for row in fwnreader: 
	hivseqneut.append(row)
	
lenseqmatrix = len(hivseqmatrix)

size20 = int(numpy.floor(lenseqmatrix*0.2))
randind = random.sample(range(lenseqmatrix), size20)
	
for i in range(lenseqmatrix):
	if i in randind:
		fw20.write(str(hivseqmatrix[i])) 
		fw20.write('\n')
		fwn20.write(str(hivseqneut[i]))
		fwn20.write('\n')
		fwi20.write(str(i) + '\n')
	else:
		fw80.write(str(hivseqmatrix[i])) 
		fw80.write('\n')
		fwn80.write(str(hivseqneut[i]))
		fwn80.write('\n')
		fwi80.write(str(i) + '\n')
	

uniqueaa = []
uniqueaa20 = [] 
uniqueaa80 = [] 

rowlen = 0



# ALL SEQUENCES
# Find list of unique amino acids per position 
for i in range(len(hivseq[0])):
	aas = 	[]
	j=0
	rowstr =''
	for j in  range(len(hivseq)):
			
		aa = hivseq[j][i]
		if aa not in aas and aa !='-' and aa !='#':
			aas.append(aa)	
			rowstr = rowstr + str(aa) +','
			rowlen= rowlen +1	 
	
	uniqueaa.append(list(aas))
	rowstr = rowstr.strip(',')
	fwaa.write(rowstr + '\n')
	
print(uniqueaa)


# ALL SEQUENCES 
aamatrix = [] 

# Now take all and convert to matrix form 
for i in range(len(hivseq)):
	aas = hivseq[i]
	aarow =[]	
	for j in  range(len(aas)):
		aa = hivseq[i][j]
		aachoices = uniqueaa[j]
		aavector = numpy.zeros(len(aachoices))
			
		if aa !='-' and aa != '#':
			indexaa = aachoices.index(aa)
			aavector[indexaa] = 1
		aarow.append(aavector)
				
	aamatrix.append(aarow)
	
# Convert this matrix to a csv format 
for k in range(len(aamatrix)):
	rowstr =''
	for m in range(len(aamatrix[k])):
		for l in range(len(aamatrix[k][m])):
			rowstr = rowstr + str(int(aamatrix[k][m][l])) + str(',') 
	
	# take the last comma off 
	rowstr = rowstr.strip(',')
	fw.write(rowstr + '\n')
	
	



		
fwn.close()
fw20.close()
fw80.close()
fw.close()
fwaa.close()