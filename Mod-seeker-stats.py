#!/usr/bin/env python

#Mod-seq data analysis pipeline
#Copyright (C) 2014  Yizhu Lin
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import csv, copy, datetime
import math

Q = 0.05 #FDR for BHcontrol
OR_TH = 1.5
TAIL = 45
filelist = []

class geneModCount(object):
    
    def __init__(self, name, chrms, strain, featureType, posStart, posEnd):
        self.name = name
        self.chrms = chrms
        self.strain = strain
        self.featureType = featureType
        self.Start = int(posStart)
        self.End = int(posEnd)
        self.length = self.End - self.Start+1
        self.countArray = [0 for i in xrange(self.length)]
        self.count = "NA"
        self.countPerNt = "NA"
        self.testArray = ['NA' for i in xrange(self.length)]
        self.stat()
        self.getDescription()
        
    def stat(self):
        self.count = sum(self.countArray)
        self.countPerNt = float(self.count)/self.length
    
    def getDescription(self):
        self.stat()
        self.description = [self.name, self.chrms, self.strain, self.featureType, self.Start, self.End, self.length, self.count, self.countPerNt]
        return(self.description)    
    
    def copyGene(self, newName):
        newGene = geneModCount(newName, self.chrms, self.strain,self.featureType, self.Start, self.End)
        newGene.length = self.length
        newGene.countArray = copy.deepcopy(self.countArray)
        newGene.count = self.count
        newGene.countPerNt = self.countPerNt
        newGene.testArray = copy.deepcopy(self.testArray)
        return newGene
        
    def mergeGenes(self, other, newName):
        newGene = self.copyGene(newName)
        for i in xrange(min(newGene.length, other.length)):
            newGene.countArray[i] = int(self.countArray[i])+int(other.countArray[i])
            newGene.count = int(self.count) + int(other.count)
        return newGene  



def resetSetting():
    #TODO: this func is not being used
    try:
        with open("stats_settings.txt", 'w') as fh:
            fh.write("""FDR = 0.05
Odds_Ratio_Threshold = 1.5
Number_Of_Replicates = 2

File_List:
Treated_1 = filename1
Control_1 = filename2
Treated_2 = filename3
Control_2 = filename4
""")
            print "stats_settings.txt file has been reset."
    except IOError:
        print "Failed to reset 'stats_settings.txt'."

def parseSetting(fh):
    global Q  #FDR for BHcontrol
    global OR_TH
    global TAIL
    
    reps = 2
    #import from setting.txt
    lines = fh.readlines()
    if len(lines) < 3:
        return False

    for line in lines:
        if line.startswith("FDR"):
            Q = float(line.split('=')[1].strip())
        if line.startswith("Odds_Ratio_Threshold"):
            OR_TH = float(line.split('=')[1].strip())
        if line.startswith("Number_Of_Replicates"):
            reps = int(line.split('=')[1].strip())
    for i in xrange(2*reps):
        filelist.append('')
    for i in xrange(reps):
        for line in lines:
            if line.startswith("Treated_%d" % (i+1)):
                filelist[2*i] = line.split('=')[1].strip()
            if line.startswith("Control_%d" % (i+1)):
                filelist[2*i+1] = line.split('=')[1].strip()           
    
    #check all settings:
    #TODO: check ref files exist, check adapters have correct format
    if Q and type(Q)== float and OR_TH and type(OR_TH) == float:
        print "FDR:", Q
        print "Odds Ratio Threshold:", OR_TH
    else:
        print "wrong format!"
        return False
    count = 0
    for filename in filelist:
        if filename:
            try:
                with open(filename, 'r') as fh:
                    fh.close()
                    count += 1
            except:
                print "Failed to read file ", "'"+filename+"'"
                return False
    if count != 2*reps:
        print "%d files is expected with %d replicates. Only %d valid files were found." %(2*reps, reps, count)
    else:
        for i in xrange(reps):
            print "Treated", i+1, ':', filelist[2*i]
            print "Control", i+1, ':', filelist[2*i+1]
        return True


def setting():
    print "Importing from 'stats_settings.txt ...'"
    try:
        with open("stats_settings.txt", 'r') as fh:
            success = parseSetting(fh)
            if not success:
                print "ERROR: setting.txt has wrong format."
                resetSetting()
                exit()
    except IOError:
        print "ERROR: Failed to open 'setting.txt'."
        resetSetting()
        exit()
    
def readData(filename):
    fh = open(filename, "r")
    print "From", filename,"reading..."
    data = dict()
    line = fh.readline()
    while(line):
        if line.startswith(">"): #gene discription
            line = line.split()
            gene = geneModCount(*line[1:7])
            gene.length = int(line[7])
            gene.count = int(line[8])
            gene.countPerNt = float(line[9])
            gene.countArray = fh.readline().split()#read new line
            data[gene.name]=gene #add gene to dict
        line = fh.readline()#read another line, should starts with '>'
    fh.close()
    print "finished."
    return data

def mergeWrap(data, names, newGeneName):
    l = len(names)
    if names[0] in data:
        newGene = data[names[0]]
    for i in xrange(l-1):
        if names[i+1] in data:
            #print 'merge',names[i+1], 'length', data[names[i+1]].length
            newGene = newGene.mergeGenes(data[names[i+1]],newGeneName)
    data[newGeneName] = newGene
    return data
        
def pchisq(x):
    #return p-value of chi-sq test with df=1
    return 1 - math.erf(math.sqrt(0.5*float(x)))

def testChisq(array):
    #calculate chi-value of 2X2 table
    # a c
    # b d
    #Pearson chi-sq with Yates' continuity
    (a, b, c, d)= tuple(array)
    n = float(sum([a, b, c, d]))
    Ea = (a+b)*(a+c)/n
    Eb = (a+b)*(b+d)/n
    Ec = (a+c)*(c+d)/n
    Ed = (b+d)*(c+d)/n
    #chi = (a-Ea)**2/Ea + (b-Eb)**2/Eb +(c-Ec)**2/Ec +(d-Ed)**2/Ed
    chi = (abs(a-Ea)-0.5)**2/Ea + (abs(b-Eb)-0.5)**2/Eb +(abs(c-Ec)-0.5)**2/Ec +(abs(d-Ed)-0.5)**2/Ed
    OR = a*d/float(b*c)
    SE = math.sqrt(1/float(a)+1/float(b)+1/float(c)+1/float(d))
    z = 1.959964
    ORL = OR*math.exp(SE*(-z))
    ORU = OR*math.exp(SE*z)    
    return (chi, pchisq(chi), OR, ORL, ORU)

def testCMH(array):
    
    #type(array) == list
    #len(array) must be multiple of 4
    assert len(array)%4 == 0
    reps = len(array)/4
    chit1Sum = 0
    chit2Sum = 0
    ORt1Sum = 0
    ORt2Sum = 0
    SEt1_nSum = 0
    SEt1_dSum = 0
    SEt2_nSum = 0
    SEt2_dSum = 0
    SEt3_nSum = 0
    for i in xrange(reps):
        a = array[4*i]
        b = array[4*i+1]
        c = array[4*i+2]
        d = array[4*i+3]
        n = float(a+b+c+d)
         
        chit1 = a-(a+b)*(a+c)/n
        chit2 = (a+b)*(a+c)*(b+d)*(c+d)/(n**3-n**2)
        chit1Sum+=chit1
        chit2Sum+=chit2

        ORt1 = a*d/n
        ORt2 = b*c/n
        ORt1Sum+=ORt1
        ORt2Sum+=ORt2
        
        SEt1_n = (a+d)*a*d/n**2
        SEt1_d = a*d/n
        SEt2_n = (b+c)*b*c/n**2
        SEt2_d = b*c/n                
        SEt3_n = ((a+d)*b*c+(b+c)*a*d)/n**2
        SEt1_nSum+=SEt1_n
        SEt1_dSum+=SEt1_d
        SEt2_nSum+=SEt2_n
        SEt2_dSum+=SEt2_d
        SEt3_nSum+=SEt3_n
        
    chi = (abs(chit1Sum)-0.5)**2/chit2Sum
    OR = ORt1Sum/ORt2Sum
    var = SEt1_nSum/(2*SEt1_dSum**2)+SEt2_nSum/(2*SEt2_dSum**2)+SEt3_nSum/(2*SEt1_dSum*SEt2_dSum)
    SE = math.sqrt(var)
    z = 1.959964
    ORL = OR*math.exp(SE*(-z))
    ORU = OR*math.exp(SE*z)
    return (chi, pchisq(chi), OR, ORL, ORU)

def checkGeneExist(gene, datalist):
    for data in datalist:
        if gene not in data:
            return False
    return True

def test_Wrap(datalist):
    assert len(datalist)%2 == 0
    reps = len(datalist)/2    
    output = dict()

    for gene in datalist[0]:
        
        if checkGeneExist(gene, datalist):
            output[gene]=datalist[0][gene]
            output[gene].testArray = ["NA" for i in xrange(output[gene].length)]
            pvalues = []
            for i in xrange(datalist[0][gene].length):
                counts = []
                array = []
                for j in xrange(len(datalist)):
                    array.append(max(int(datalist[j][gene].countArray[i]), 1))
                    array.append(max(int(datalist[j][gene].count), 1))
                    counts.append(int(datalist[j][gene].countArray[i]))
                    counts.append(int(datalist[j][gene].count))
                if len(array) == 4:
                    (chi, p, OR, ORL, ORU) = testChisq(array)
                elif (len(array)>=8 and len(array)%4 == 0):
                    (chi, p, OR, ORL, ORU) = testCMH(array)
                else:
                    print "Number of input files should be mulitple of 4."
                    exit()
                output[gene].testArray[i] = list((chi, p, 'NA', OR, ORL, ORU))+counts #'NA' is for p_adjusted
                pvalues.append(p)

            p_sig = BHcontrol(pvalues) #list pvalues has been sorted in BHcontrol
            #if gene == 'RDN18' or gene == 'RDN25':
                #print pvalues
            if p_sig != 'NA':
                for i in xrange(output[gene].length):
                    if output[gene].testArray[i][1] <= p_sig:
                        j = pvalues.index(output[gene].testArray[i][1])+1 #rank
                        output[gene].testArray[i][2] = output[gene].testArray[i][1]*output[gene].length/float(j)
    return output

def BHcontrol(pvalues):
    pvalues.sort()
    m = float(len(pvalues))
    p = 0
    i = 0
    while p<=((i+1)/m)*Q and i < m :
        p = pvalues[i]
        i += 1
    p_sig =  pvalues[max(i-1, 0)]
    if p_sig <= Q:
        #print "BHcontrol:", p_sig, m, i-1
        return p_sig
    else: return 'NA'

def writeData(outputfile, output):
    fh = open(outputfile, 'w')
    mywriter = csv.writer(fh, delimiter='\t', lineterminator = '\n')
    fileheader = []
    for i in xrange(len(filelist)):
        fileheader.append(filelist[i])
        fileheader.append(filelist[i]+'_SUM')
    header = ['geneName', 'chrom', 'strand', 'type', 'start', 'end', 'length', 'pos_on_gene',
              'chisq', 'p_origin', 'p_adjusted', 'odds_ratio', 'OR_lower', 'OR_upper']+fileheader
    mywriter.writerow(header)
    for gene in output:
        for i in xrange(max(output[gene].length-TAIL, 0)):
            p_adjusted = output[gene].testArray[i][2]
            OR = output[gene].testArray[i][3]
            if p_adjusted != "NA" and p_adjusted <= Q and OR > OR_TH:
                #refomat chi, p_origin, p_adjusted, OR, OR_L, OR_U
                chi = '{:.3f}'.format(output[gene].testArray[i][0])
                p_origin = '{:.4g}'.format(output[gene].testArray[i][1])
                p_adjusted = '{:.4g}'.format(output[gene].testArray[i][2])
                OR = '{:.3f}'.format(output[gene].testArray[i][3])
                OR_l = '{:.3f}'.format(output[gene].testArray[i][4])
                OR_u = '{:.3f}'.format(output[gene].testArray[i][5])
                row = output[gene].description[:7]+[i+1]+[chi, p_origin, p_adjusted, OR, OR_l, OR_u]+output[gene].testArray[i][6:]
                mywriter.writerow(row)
    fh.close()   

def main(filelist, outputfile):
    setting()
    dataDictsList = []
    for filename in filelist:
        data = readData(filename)
        dataDictsList.append(data)
    output = test_Wrap(dataDictsList)
    writeData(outputfile, output)



start = datetime.datetime.now()
outputfile = 'Stats_output_'+start.strftime("%Y-%m-%d_%H-%M-%S.txt")
main(filelist,  outputfile)
end = datetime.datetime.now()
print "Overall runtime:", end - start