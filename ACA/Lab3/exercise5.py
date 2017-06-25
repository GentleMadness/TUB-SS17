#!/usr/bin/env python
# -*- coding: utf8 -*-
from __future__ import print_function
import subprocess as sub
import re, os

if __name__ == "__main__":
	FILE_LOC = "/home/shankuan/ACA/SimpleScalar/benchmarks"
	DEFAULT_MEMORY_LATENCY = 18
	INCREASE_MEMORY_LATENCY = 200
	SETS_CACHE = ["1024", "2048", "4096", "8192"]
	ASSOCIATIVITY_CACHE = ["2", "4", "8", "16", "32"]
	L2_LATENCY = ["6", "7", "8"]
	exercise_file = open("%s/Exercise5.txt" %FILE_LOC, 'w+')
	print ("AMAT (Hit Rate + Miss Rate * Miss Penalty)")
	exercise_file.write("AMAT (Hit Rate + Miss Rate * Miss Penalty)\n")
	exercise_file.write("Default Memory latency <18 cycles>\n")
	print("Default Memory latency <18 cycles>")
	
	cmd = ("%s/sim-outorder -issue:width 4 -issue:inorder false -res:ialu 4 -res:imult 1 -res:fpalu 4 -res:fpmult 1 -ruu:size 32 -lsq:size 16 -bpred comb %s/memcopy 2>&1 | egrep 'sim_IPC|ul2.accesses|ul2.hits|ul2.miss_rate'" %(FILE_LOC, FILE_LOC))
	result = sub.check_output(cmd, shell=True).splitlines()
	dict_amat = {}
	for _list_line in result:
		_list = _list_line.split()
		if _list[0] == "sim_IPC":
			print (" IPC:", end='')
			print (" "+_list[1])
			exercise_file.write(" IPC:")
			exercise_file.write(" "+_list[1]+"\n")
		else:
			dict_amat[_list[0]] = float(_list[1])
	amat = dict_amat["ul2.hits"]/dict_amat["ul2.accesses"]+dict_amat["ul2.miss_rate"]*DEFAULT_MEMORY_LATENCY
	print(" AMAT: %lf\n" %amat)
	exercise_file.write(" AMAT: %lf\n\n" %amat)

	print("========================================================\n")
	exercise_file.write("========================================================\n")

	exercise_file.write("Increase Memory latency <200 cycles>\n")
	print("Increase Memory latency <200 cycles>")
	cmd = ("%s/sim-outorder -issue:width 4 -issue:inorder false -res:ialu 4 -res:imult 1 -res:fpalu 4 -res:fpmult 1 -ruu:size 32 -lsq:size 16 -bpred comb -mem:lat 200 2 %s/memcopy 2>&1 | egrep 'sim_IPC|ul2.accesses|ul2.hits|ul2.miss_rate'" %(FILE_LOC, FILE_LOC))
	result = sub.check_output(cmd, shell=True).splitlines()
	dict_amat = {}
	for _list_line in result:
		_list = _list_line.split()
		if _list[0] == "sim_IPC":
			print (" IPC:", end='')
			print (" "+_list[1])
			exercise_file.write(" IPC:")
			exercise_file.write(" "+_list[1]+"\n")
		else:
			dict_amat[_list[0]] = float(_list[1])
	amat = dict_amat["ul2.hits"]/dict_amat["ul2.accesses"]+dict_amat["ul2.miss_rate"]*INCREASE_MEMORY_LATENCY
	print(" AMAT: %lf\n" %amat)
	exercise_file.write(" AMAT: %lf\n\n" %amat)


	for _l2_latency in L2_LATENCY:
		print("========================================================\n")
		exercise_file.write("========================================================\n")
		exercise_file.write("L2 Latency <%s cycles>\n" %_l2_latency)
		print("L2 Latency <%s cycles>\n" %_l2_latency)
		for _setcache in SETS_CACHE:
			for _assoccache in ASSOCIATIVITY_CACHE:
				print("<L2 Cache Latency: %s, Sets in the cache: %s, Associativity: %s>" %(_l2_latency, _setcache, _assoccache))
				exercise_file.write("<L2 Cache Latency: %s, Sets in the cache: %s, Associativity: %s>\n" %(_l2_latency, _setcache, _assoccache))
				
				cmd = ("%s/sim-outorder -issue:width 4 -issue:inorder false -res:ialu 4 -res:imult 1 -res:fpalu 4 -res:fpmult 1 -ruu:size 32 -lsq:size 16 -bpred comb -mem:lat 200 2 -cache:il1 il1:128:64:1:l -cache:il2 dl2 -cache:dl1 dl1:256:32:1:l -cache:dl2 ul2:%s:64:%s:l -cache:dl2lat %s -cache:il2lat %s %s/memcopy 2>&1 | egrep 'sim_IPC|ul2.accesses|ul2.hits|ul2.miss_rate'" %(FILE_LOC, _setcache, _assoccache, _l2_latency, _l2_latency, FILE_LOC))
				result = sub.check_output(cmd, shell=True).splitlines()
				dict_amat = {}
				for _list_line in result:
					_list = _list_line.split()
					if _list[0] == "sim_IPC":
						print (" IPC:", end='')
						print (" "+_list[1])
						exercise_file.write(" IPC:")
						exercise_file.write(" "+_list[1]+"\n")
					else:
						dict_amat[_list[0]] = float(_list[1])
				amat = dict_amat["ul2.hits"]/dict_amat["ul2.accesses"]+dict_amat["ul2.miss_rate"]*INCREASE_MEMORY_LATENCY
				print(" AMAT: %lf\n" %amat)
				exercise_file.write(" AMAT: %lf\n\n" %amat)
	exercise_file.close()
