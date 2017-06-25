#!/usr/bin/env python
# -*- coding: utf8 -*-
from __future__ import print_function
import subprocess as sub
import re, os

if __name__ == "__main__":
	EX1_PRED = ["nottaken", "bimod", "2lev", "comb", "perfect"]
	BENCHMARK = ["fibonacci", "matmul", "pi", "whetstone", "memcopy"]
	FILE_LOC = "/home/yuting/Desktop/SimpleScalar/benchmarks"
	EX2_OPTIONS = ["1", "2", "4", "8"]
	EX3_REG = [["16", "8"], ["32", "8"], ["32", "16"], ["64", "16"], ["64", "32"]]
	
	exercise1_file = open("%s/Exercise.txt" %FILE_LOC, 'w+')

	exercise1_file.write("Exercise1:\n")
	exercise1_file.write('["nottaken", "bimod", "2lev", "comb", "perfect"]\n')
	print ("Exercise1:")
	print ('["nottaken", "bimod", "2lev", "comb", "perfect"]')
	for ben_name in BENCHMARK:
		temp = []
		exercise1_file.write(" "+ben_name+":\n [")
		print (" "+ben_name+":\n [", end='')
		for pre in EX1_PRED:
			cmd = ("%s/sim-outorder -fetch:ifqsize 1 -decode:width 1 -issue:width 1 -commit:width 1 -issue:inorder true -res:ialu 1 -res:imult 1 -res:fpalu 1 -res:fpmult 1 -bpred %s %s/%s 2>&1 | grep 'sim_IPC'" %(FILE_LOC, pre, FILE_LOC, ben_name))
			temp.append(sub.check_output(cmd, shell=True).split()[1])
			#print (sub.check_output(cmd, shell=True).split()[1], end='')
		exercise1_file.write("%s]\n" %(', '.join(temp)))
		print ("%s]" %(', '.join(temp)))
	#exercise1_file.close()

	print ("\n========================================================\n")
	exercise1_file.write("\n========================================================\n")
	exercise1_file.write("Exercise2.1 <in-order>:\n")
	exercise1_file.write('["1", "2", "4", "8" "investigate"]\n')
	print ("Exercise2.1 <in-order>:")
	print ('["1", "2", "4", "8" "investigate"]')
	for ben_name in BENCHMARK:
		temp = []
		exercise1_file.write(" "+ben_name+":\n [")
		print (" "+ben_name+":\n [", end='')
		for _option in EX2_OPTIONS:
			cmd = ("%s/sim-outorder -fetch:ifqsize %s -decode:width %s -issue:width %s -commit:width %s -issue:inorder true -res:ialu %s -res:imult %s -res:fpalu %s -res:fpmult %s -bpred 2lev %s/%s 2>&1 | grep 'sim_IPC'" %(FILE_LOC, _option, _option, _option, _option, _option, _option, _option, _option, FILE_LOC, ben_name))
			temp.append(sub.check_output(cmd, shell=True).split()[1])
		cmd = ("%s/sim-outorder -fetch:ifqsize 8 -decode:width 8 -issue:width 8 -commit:width 8 -issue:inorder true -res:ialu 4 -res:imult 1 -res:fpalu 4 -res:fpmult 1 -bpred 2lev %s/%s 2>&1 | grep 'sim_IPC'" %(FILE_LOC, FILE_LOC, ben_name))
		temp.append(sub.check_output(cmd, shell=True).split()[1])
		exercise1_file.write("%s]\n" %(', '.join(temp)))
		print ("%s]" %(', '.join(temp)))
	
	
	print ("\n========================================================\n")
	exercise1_file.write("\n========================================================\n")
	exercise1_file.write("Exercise2.2 <out-of-order>:\n")
	exercise1_file.write('["1", "2", "4", "8" "investigate"]\n')
	print ("Exercise2.2 <out-of-order>:")
	print ('["1", "2", "4", "8" "investigate"]')
	for ben_name in BENCHMARK:
		temp = []
		exercise1_file.write(" "+ben_name+":\n [")
		print (" "+ben_name+":\n [", end='')
		for _option in EX2_OPTIONS:
			cmd = ("%s/sim-outorder -fetch:ifqsize %s -decode:width %s -issue:width %s -commit:width %s -issue:inorder false -res:ialu %s -res:imult %s -res:fpalu %s -res:fpmult %s -bpred 2lev %s/%s 2>&1 | grep 'sim_IPC'" %(FILE_LOC, _option, _option, _option, _option, _option, _option, _option, _option, FILE_LOC, ben_name))
			temp.append(sub.check_output(cmd, shell=True).split()[1])
		cmd = ("%s/sim-outorder -fetch:ifqsize 8 -decode:width 8 -issue:width 8 -commit:width 8 -issue:inorder false -res:ialu 4 -res:imult 1 -res:fpalu 4 -res:fpmult 1 -bpred 2lev %s/%s 2>&1 | grep 'sim_IPC'" %(FILE_LOC, FILE_LOC, ben_name))
		temp.append(sub.check_output(cmd, shell=True).split()[1])
		exercise1_file.write("%s]\n" %(', '.join(temp)))
		print ("%s]" %(', '.join(temp)))


	print ("\n========================================================\n")
	exercise1_file.write("\n========================================================\n")
	exercise1_file.write("Exercise3.1 <default amount>:\n")
	exercise1_file.write('["ruu=16, lsq=8", "ruu=32, lsq=8", "ruu=32, lsq=16", "ruu=64, lsq=16", "ruu=64, lsq=32"]\n')
	print ("Exercise3.1 <default amount>:")
	print ('["ruu=16, lsq=8", "ruu=32, lsq=8", "ruu=32, lsq=16", "ruu=64, lsq=16", "ruu=64, lsq=32"]')
	
	for ben_name in BENCHMARK:
		temp = []
		exercise1_file.write(" "+ben_name+":\n [")
		print (" "+ben_name+":\n [", end='')
		for _reg in EX3_REG:
			cmd = ("%s/sim-outorder -fetch:ifqsize 8 -decode:width 8 -issue:width 8 -commit:width 8 -issue:inorder false -res:ialu 4 -res:imult 1 -res:fpalu 4 -res:fpmult 1 -ruu:size %s -lsq:size %s -bpred 2lev %s/%s 2>&1 | grep 'sim_IPC'" %(FILE_LOC, _reg[0], _reg[1], FILE_LOC, ben_name))
			temp.append(sub.check_output(cmd, shell=True).split()[1])
		exercise1_file.write("%s]\n" %(', '.join(temp)))
		print ("%s]" %(', '.join(temp)))


	print ("\n========================================================\n")
	exercise1_file.write("\n========================================================\n")
	exercise1_file.write("Exercise3.2 <double amount>:\n")
	exercise1_file.write('["ruu=16, lsq=8", "ruu=32, lsq=8", "ruu=32, lsq=16", "ruu=64, lsq=16", "ruu=64, lsq=32"]\n')
	print ("Exercise3.2 <double amount>:")
	print ('["ruu=16, lsq=8", "ruu=32, lsq=8", "ruu=32, lsq=16", "ruu=64, lsq=16", "ruu=64, lsq=32"]')
	for ben_name in BENCHMARK:
		temp = []
		exercise1_file.write(" "+ben_name+":\n [")
		print (" "+ben_name+":\n [", end='')
		for _reg in EX3_REG:
			
			cmd = ("%s/sim-outorder -fetch:ifqsize 8 -decode:width 8 -issue:width 8 -commit:width 8 -issue:inorder false -res:ialu 8 -res:imult 2 -res:fpalu 8 -res:fpmult 2 -ruu:size %s -lsq:size %s -bpred 2lev %s/%s 2>&1 | grep 'sim_IPC'" %(FILE_LOC, _reg[0], _reg[1], FILE_LOC, ben_name))
			temp.append(sub.check_output(cmd, shell=True).split()[1])
		exercise1_file.write("%s]\n" %(', '.join(temp)))
		print ("%s]" %(', '.join(temp)))
	exercise1_file.close()








			
	
