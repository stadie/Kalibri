#! /usr/bin/python

import os
import sys

algo_types_all = ['ak5Calo_','ak5FastCalo_','ak5FastPF_','ak5JPT_','ak7PF_','ak7Calo_','ak5PFCHS_']

length=len(sys.argv)
samplelist=[]

for i in range(1,length):
#	print str(i) + ": " + sys.argv[i]
	samplelist.append(sys.argv[i])

targetdir=samplelist[length-2]


print "Use this script to set links to n-tuple files from multiple samples in a common targetdir."
print "Usage: ./setlinks.py sample1 sample2 ... targetdir"
print "------------------------------------------------------------------------------------------"
print "number of arguments: " + str(length-1) + "; will set links to " + str(length-2) + " samples in targetdir"
print "Warning: targetdir " + targetdir + " will be deleted before setting new links"


if len(sys.argv)<3:
	print "error occured..."
	print "usage: ./setlinks.py sample1 sample2 ... targetdir"
	sys.exit()

rep = ''
while not rep in [ 'c', 'C' ]:
	for j, sample in enumerate(samplelist):
		if j<(length-2):
			print "copying " + sample + " to " + targetdir
	rep = raw_input( 'enter "c" to continue ("Ctrl+c" to abort): ' )
	if 1 < len(rep):
		rep = rep[0]


print "targetdir: " +  targetdir

if os.path.exists(targetdir):
    os.system("rm -r " + targetdir)
    os.system("mkdir " + targetdir)
else:
    os.system("mkdir " + targetdir)

for algo_types in algo_types_all:
	for j, sample in enumerate(samplelist):
		if j<(length-2):
			print "setting links to " + algo_types+ " " + sample + " in " + targetdir
			for i in range(0,10):
#				print i
				os.system("ln -s " +sample +"/"+algo_types+str(i)+".root " + targetdir +"/"+algo_types+str(i)+"_sam"+str(j)+".root")

