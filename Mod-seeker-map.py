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

import os
import sys
import subprocess
import time
import datetime
import csv
import copy
import threading
import Queue

#initiate args
inputpath = "InputFiles/"
refpath = ''
ref = ''
intersectRef = ''
intersectRefType = ''
adapter5 = ''
adapter3 = ''
now = datetime.datetime.now()
outputpath = "Results_"+now.strftime("%Y-%m-%d_%H-%M-%S/")
bowtie_p = 16
Max_threads = 6

class geneModCount(object):
    
    def __init__(self, name, chrms, strain, featureType, posStart, posEnd):
        self.name = name
        self.chrms = chrms
        self.strain = strain
        self.featureType = featureType
        self.Start = int(posStart)
        self.End = int(posEnd)
        self.length = self.End - self.Start+1
        self.count = "NA"
        self.countPerNt = "NA"
        self.countArray = [0 for i in xrange(self.length)]
        self.description = [self.name, self.chrms, self.strain, self.featureType, self.Start, self.End, self.length, self.count, self.countPerNt]   
    
    def getDescription(self):
        self.description = [self.name, self.chrms, self.strain, self.featureType, self.Start, self.End, self.length, self.count, self.countPerNt]
        
        
def welcome():
    print "This is a mod-seq data analysis software. Run on linux."
    return

#import setting args from 'map_settings.txt'
def resetSetting():
    #TODO: this func is not being used
    try:
        with open("map_settings.txt", 'w') as fh:
            fh.write("""ref_path = ~/bowtie_indexes/
ref_genome = S_cer
GeneAnnotationFile = ~/SGD_SacCer3_filtered.gff
adapter5 = ^ATCGTAGGCACCTGAAA
adapter3 = CTGTAGGCACCATCAAT
output_path = default
bowtie_alignment_threads = 16
max_python_threads = 6
""")
            print "map_settings.txt file has been reset."
    except IOError:
        print "Failed to reset 'map_settings.txt'."

def parseSetting(fh):
        
    global refpath
    global ref
    global intersectRef
    global intersectRefType
    global adapter5
    global adapter3
    global now
    global outputpath
    global bowtie_p
    global Max_threads
    
    #import from map_settings.txt
    lines = fh.readlines()
    if len(lines) < 8:
        return False
    if lines[0].startswith("ref_path"):
        refpath = lines[0].split('=')[1].strip()
    if lines[1].startswith("ref_genome"):
        ref = lines[1].split('=')[1].strip()
    if lines[2].startswith("GeneAnnotationFile"):
        intersectRef = lines[2].split('=')[1].strip()
        intersectRefType = intersectRef.split('.')[-1].strip()
    if lines[3].startswith("adapter5"):
        adapter5 = lines[3].split('=')[1].strip().upper()
    if lines[4].startswith("adapter3"):
        adapter3 = lines[4].split('=')[1].strip().upper()
    if lines[5].startswith("output_path"):
        outputpath_new = lines[5].split('=')[1].strip()
        if outputpath_new != 'default':
            outputpath = outputpath_new
        subprocess.check_output("mkdir "+outputpath, shell=True)
    if lines[6].startswith("bowtie_alignment_threads"):
        bowtie_p = int(lines[6].split('=')[1].strip())
    if lines[7].startswith("max_python_threads"):
        Max_threads = int(lines[7].split('=')[1].strip())
    
    #check all settings:
    #TODO: check ref files exist, check adapters have correct format
    if refpath and ref and intersectRef and outputpath:
        if not (intersectRefType == 'gff' or intersectRefType=='bed'):
            print "intersectRef file %s is invalid." % intersectRef
            print intersectRefType
            print "Please choose a bed file or a gff file."
            return False
        print "bowtie index location: "+refpath
        print "Aligned to: "+ref
        print "intersectBed reference: "+intersectRef
        print "intersectBed reference type: "+intersectRefType
    else:
        return False    
    if adapter5:
        print "5' adapter: "+adapter5
    else:
        print "Warning: no 5' adapter."
    if adapter3:
        print "3' adapter: "+adapter3
    else:
        print "Warning: no 3' adapter."
    print "All settings imported."
    return True

def setting():
    print "Importing from 'map_settings.txt ...'"
    try:
        with open("map_settings.txt", 'r') as fh:
            success = parseSetting(fh)
            if not success:
                print "ERROR: map_settings.txt has wrong format."
                resetSetting()
                exit()
    except IOError:
        print "ERROR: Failed to open 'map_settings.txt'."
        resetSetting()
        exit()

def envCheck():
    try:
        version = subprocess.check_output("bowtie --version", shell=True)
        print "bowtie is installed."
    except:
        print "ERROR: bowtie is not properly installed."
        exit()
    try:
        version = subprocess.check_output("cutadapt --version", shell=True)
        print "cutadapt is installed."
    except:
        print "ERROR: cutadapt is not properly installed."
        exit()
    try:
        version = subprocess.check_output("bedtools --version", shell=True)
        print "bedtools is installed"
    except:
        print "ERROR: bedtools is not properly installed."
        exit()
    #TODO: test samtools
    return

def test():
    #TODO
    return

#get filenames in path
def listFiles(path):
    if (os.path.isdir(path) == False):
        print 'not folder'
        return [path]
    else:
        files = []
        samples = []
        for filename in os.listdir(path):
            sampleName = filename.split('.')[0]
            if filename.endswith(".fastq.gz") or filename.endswith(".fastq"):
                if sampleName not in samples:
                    if filename.endswith(".fastq"):
                        subprocess.check_output("gzip -c "+inputpath+filename, shell=True)
                    files.append(filename)
                    samples.append(sampleName)
                else:
                    print ("Warning: Multiple files named "+sampleName+
                           ".fastq(.gz) were detected. '"+filename+"' is ignored.")
        return samples

def printList(samples):
    print "Files in path:"
    for sampleName in samples:
        print sampleName
    return

def bamToSam(filename):
    assert(filename.endswith(".bam"))
    outfilename = filename[:-3]+".sam"
    cmd = "samtools view -bS "
    mycmd = cmd + filename + " > "+ outfilename
    subprocess.check_output(mycmd, shell=True)
    return(outfilename)

#if input is a sam file, convert to bam
def fixWrap(inputfile, outputfile, mismatchfile):
    if inputfile.endswith(".bam"):
        inputfile = bamToSam(inputfile)
    fhIn = open(inputfile, "r")
    fhOut = open(outputfile, "w")
    fhOut.close()
    fhOut = open(outputfile, "a")
    fhMis = open(mismatchfile, "w")
    fhMis.close()
    fhMis = open(mismatchfile, "a")

    i = 0
    mis = 0
    writerOut = csv.writer(fhOut, delimiter = "\t")
    writerMis = csv.writer(fhMis, delimiter = "\t")
    #read lines, loop through file
    line = fhIn.readline()   
    while line:
        lineOrigin = line
        lineStr,misCount = fix(line)
        if line :           
            #writerOut.writerow(line)
            fhOut.write(lineStr)
        if misCount:
            #writerMis.writerow(lineOrigin)
            fhMis.write(lineOrigin)
            mis += 1
        line = fhIn.readline()
        i += 1
    """    
    print 'Fix 5p mismatch: '+inputfile,
    print "Number of reads:",i,';',
    print "Number of reads with 5p mismatch:", mis
    """
    fhIn.close()
    fhOut.close()
    fhMis.close()


def fix(line):
    #sam file format example:
    #HISEQ:108:H7N5WADXX:1:1107:18207:18800  0       gi|207113128|ref|NR_002819.2|   553     255     51M     *       0       0       AAAATTTCCGTGCGGGCCGTGGGGGGCTGGCGGCAACTGGGGGGCCGCAGA     BBBFFFFFFFFFFIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFF   XA:i:0  MD:Z:51 NM:i:0
    #QNAME FLAG RNAME POS(leftmost) MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    misCount = 0
    lineStr = line
    line = line.split()
    if len(line) >= 13:
        m = line[12].split(':')
        #if(line[1]=='0' or line[1]=='16') and len(m[2]) > 2: #mismatch detected
        #if (line[1]=='0' or line[1]=='16'):
        #if m[1] != 'Z':

        if line[1]=='0': #postive strand, remove 5' end

            while m[2].startswith('0'):#have a mismatch
                m[2] = m[2][2:] #strip first two character
                misCount += 1 #mismatch count +1
            if misCount != 0 :
                line[3] = str(int(line[3])+misCount) #move 5'start to 3' direction
                line[5] = str(int(line[5].rstrip('M'))-misCount)+'M' #change seq length
                line[9] = line[9][misCount:]
                line[10] = line[10][misCount:]
                line[12] = 'MD:Z:'+str(len(line[9]))
            lineStr = '\t'.join(line)
            lineStr = lineStr+'\n'              
        elif line[1]=='16': #negtive strand, remove 3' end
            while (m[2][-1]=='0' and m[2][-2] in "ATCG"):#have a mismatch
                m[2] = m[2][:-2] #strip last two character
                misCount += 1 #mismatch count +1
            if misCount != 0 :
                line[5] = str(int(line[5].rstrip('M'))-misCount)+'M' #change seq length
                line[9] = line[9][:-misCount]
                line[10] = line[10][:-misCount]
                line[12] = 'MD:Z:'+str(len(line[9]))
            lineStr = '\t'.join(line)
            lineStr = lineStr+'\n'
    return (lineStr, misCount)

def cutadapters(sampleName):
    #cut 3' adapter
    cmd1 = "cutadapt -m 25 -a "+adapter3
    inputfile1 = inputpath+sampleName+".fastq.gz"
    outputfile1 = outputpath+"TRIMMED_"+sampleName+".fastq"
    logfile1 = outputpath+"TRIMMED_"+sampleName+".log"          
    mycmd1 = cmd1+" "+inputfile1+" > "+outputfile1+" 2> "+logfile1
    #print "Call command:\n"+mycmd1
    subprocess.check_output(mycmd1, shell=True)
    
    #cut 5'adapter and seperate modstop from adapter stop
    cmd2 = "cutadapt -g "+adapter5+" --untrimmed-output "
    inputfile2 = outputfile1
    modStopOutput = outputpath+"ModStop_"+sampleName+".fastq"
    adaptStopOutput = outputpath+"AdaptStop_"+sampleName+".fastq"
    logfile2 = outputpath+"Mod-Adapt_"+sampleName+".log"
    mycmd2 = cmd2+modStopOutput+" "+inputfile2+" > "+adaptStopOutput+" 2> "+logfile2
    #print "Call command:\n"+mycmd2
    subprocess.check_output(mycmd2, shell=True)        
    subprocess.check_output("gzip "+adaptStopOutput, shell=True)
    subprocess.check_output("gzip "+outputfile1, shell=True)
    
def bowtieAlign(sampleName):

    #use bowtie for alignment
    cmd1 = "bowtie --best --chunkmbs 500 -p %d -t -S " % bowtie_p
    inputfile1 = outputpath+"ModStop_"+sampleName+".fastq"
    outputfile1 = outputpath+"Aligned_"+ref+'_'+sampleName+".sam"
    mycmd1 = cmd1+refpath+ref+" "+inputfile1+" "+outputfile1
    #print "Call command:\n"+mycmd1
    subprocess.check_output(mycmd1, shell=True)
    subprocess.check_output("gzip "+inputfile1, shell=True)

def fivePrimeFix(sampleName):
    #5p fix
    inputfile2 = outputpath+"Aligned_"+ref+'_'+sampleName+".sam"
    outputfile2 = outputpath+"5pFixed_"+ref+'_'+sampleName+".sam"
    mismatchfile = outputpath+"5pMisMatch_"+ref+'_'+sampleName+".sam"
    fixWrap(inputfile2, outputfile2, mismatchfile)
    subprocess.check_output("rm -f "+inputfile2, shell=True)

def samToBam(sampleName):    
    #convert sam file to bam file
    cmd3 = "samtools view -bS "
    inputfile3 = outputpath+"5pFixed_"+ref+'_'+sampleName+".sam"
    outputfile3 = outputpath+"5pFixed_"+ref+'_'+sampleName+".bam"
    mycmd3 = cmd3 + inputfile3 + " > "+ outputfile3
    subprocess.check_output(mycmd3, shell=True)
    subprocess.check_output("rm -f "+inputfile3, shell=True)

def sortBam(sampleName):    
    #sort
    cmd4 = "samtools sort -m 5000000000 "
    inputfile4 = outputpath+"5pFixed_"+ref+'_'+sampleName+".bam"
    outputfile4 = outputpath+"SORTED_"+ref+'_'+sampleName
    mycmd4 = cmd4+inputfile4+" "+outputfile4
    subprocess.check_output(mycmd4, shell=True)
    subprocess.check_output("rm -f "+inputfile4, shell=True)
    #output file name: "SORTED_foo.bam"

def intersect(sampleName):

    inputfile1 = outputpath+"SORTED_"+ref+'_'+sampleName+'.bam'
    inputfile2 = intersectRef         
    cmd1 = "intersectBed -s -wo -split -bed -abam "
    if intersectRefType == 'gff':
        cmd2 = """| awk 'BEGIN { OFS = "\t";} {split($21,a,";"); sub(/ID=/,"",a[1]); print $1, $2, $3, $4, $6, $15, $16, $17, a[1]}' """
    elif intersectRefType == 'bed':
        cmd2 = """| awk 'BEGIN { OFS = "\t";} {print $1, $2, $3, $4, $6, "NA", $14, $15, $16}' """        
    outputfile1 = outputpath+"Intersect_"+ref+'_'+sampleName+'.tab'
    mycmd1 = cmd1+inputfile1+' -b '+inputfile2+cmd2+" > "+outputfile1
    #print "Call Command: "+mycmd1
    subprocess.check_output(mycmd1, shell=True)    

def modCount(sampleName):
        
    #count ModStops on each position
    inputfile = outputpath+"Intersect_"+ref+'_'+sampleName+'.tab'
    outputfile = outputpath+"CountMod_"+ref+'_'+sampleName+".tab"
    fh = open(inputfile, 'r')
    genes = dict()
    
    """.tab file format:
    0       1       2       3       4       5       6       7       8
    chr     rStart  rEnd    rName   +/-     type    geneS   geneEnd geneName
    """
            
    line = fh.readline()
    while (line): 
        line = line.split()
        gene = line[8]
        if (int(line[1]) > int(line[6]) and int(line[2]) < int(line[7])):#if mod-stop is inside gene
            if gene not in genes:
                genes[gene]=geneModCount(line[8], line[0], *line[4:8])
            if line[4] == '+':
                if intersectRefType == 'gff':
                    genes[gene].countArray[int(line[1])-int(line[6])]+=1
                elif intersectRefType == 'bed':
                    genes[gene].countArray[int(line[1])-int(line[6])-1]+=1
            elif line[4] == '-':
                genes[gene].countArray[int(line[7])-int(line[2])-1]+=1
                
        line=fh.readline()
    #print genes[gene]
    fh.close()
    
    for gene in genes:           
        genes[gene].count = sum(genes[gene].countArray)
        genes[gene].countPerNt = float(genes[gene].count)/genes[gene].length
        genes[gene].getDescription()
        
    fh = open(outputfile, 'w')
    mywriter = csv.writer(fh, delimiter=' ')
    for gene in genes:
        #if genes[gene].countPerNt > 1:
        mywriter.writerow([">"]+genes[gene].description)
        mywriter.writerow(genes[gene].countArray)
    fh.close()

def run(sampleName):
    try:    
        cutadapters(sampleName)
    except:
        print "ERROR: Failed to run cutadapt."
        exit()
    try:
        bowtieAlign(sampleName)
    except:
        print "ERROR: Failed to run bowtie."
        exit()
    fivePrimeFix(sampleName)
    try:
        samToBam(sampleName)
        sortBam(sampleName)
    except:
        print "ERROR: Failed to run samtools."
        exit()
    try: 
        intersect(sampleName)
    except:
        print "ERROR: Failed to run intersectBed."
        exit()
    modCount(sampleName)
    #print sampleName + ' completed'


def runThreads():
    while True:
        sampleName = q.get()
        run(sampleName)
        q.task_done()

def main():
    welcome()
    setting()
    envCheck()
    test()
    
    print 'running analysis...'
    samples = listFiles(inputpath)
    if samples:
        print str(len(samples))+'files to process:'
        #printList(samples)
    else:
        print "Can not find any input fastq(.gz) files"
        exit() 
    try:
        for i in xrange(Max_threads):
            t = threading.Thread(target=runThreads)
            t.daemon = True
            t.start()
        
        for sampleName in samples:
            q.put(sampleName)
            
        q.join()
        print "done"
    except:
        exit()
    end = datetime.datetime.now()
    print "Total time: "
    print end - now
    return
    
if __name__ == '__main__':
    q = Queue.Queue()
    main() 