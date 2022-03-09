#!/usr/bin/env python

import sys
import numpy as np
import argparse

def main():

 parser = argparse.ArgumentParser()
 parser.add_argument('-i',type=str)
 args=parser.parse_args()

 in_name = args.i+".dat"
 out_name = args.i
 in_data = open(in_name,"r")
 dist_data = np.loadtxt(in_data, comments=["#","@"], dtype=float)
 in_data.close()
 
 length_data = dist_data.shape[0]
 bin_data = np.zeros(length_data)
 
# for i in range(length_data):
#     if dist_data[i,1] < 0.375:
#         bin_data[i] = 0
#     else:
#         bin_data[i] = 1
 
 for i in range(length_data):
     if dist_data[i,1] <= 0.3:
         bin_data[i] = 1
     elif dist_data[i,1] >= 0.9:
         bin_data[i] = 2
     else:
         bin_data[i] = 0

# remove unsuccesfull jumps

 for i in range(length_data):
     if i == 0:
         bin_data[i] = bin_data[i]
     else:
         if bin_data[i] == 0 and bin_data[i-1] != 0:
             bin_data[i] = bin_data[i-1]
         else:
             bin_data[i] = bin_data[i]
 
# bin_result = open("./100fs_statebin_1_0.34to0.4_2.dat","w")
 
# calculate the waiting times

 left_to_right = []
 right_to_left = []

 count_left = 0
 count_right = 0

 for i in range(length_data):
     if i == 0:
        if bin_data[i] == 1:
            count_left += 1
        else:
            count_right +=1
     else:
        if bin_data[i] == 1 and bin_data[i-1] == 1:
            count_left += 1
        elif bin_data[i] == 2 and bin_data[i-1] == 2:
            count_right += 1
        elif bin_data[i] == 1 and bin_data[i-1] == 2:
            left_to_right = np.append(left_to_right,count_left)
            count_left = 0
        elif bin_data[i] == 2 and bin_data[i-1] == 1:
            right_to_left = np.append(right_to_left,count_right)
            count_right = 0

 #bin_result.write("#t  state (0: max; 1: bound; 2: unbound)  \n")
 
 #for i in range(length_data):
 #    bin_result.write("{:15.8f} {:1.0f}\n".format(dist_data[i,0],bin_data[i]))
 #    bin_result.write("{:1.0f}\n".format(bin_data[i]))
 
 #bin_result.close()
 
 outfile_lefttoright = out_name+"_left_to_right.dat"
 outfile_righttoleft = out_name+"_right_to_left.dat"

 np.savetxt(outfile_lefttoright,left_to_right,delimiter=',')
 np.savetxt(outfile_righttoleft,right_to_left,delimiter=',')

if __name__ == "__main__":
        main()

