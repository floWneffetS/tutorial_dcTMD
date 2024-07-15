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

# set cores
 
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
 
# calculate the waiting times
# left: bound state
# right: bound state

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
            right_to_left = np.append(right_to_left,count_right)
            count_right = 0
        elif bin_data[i] == 2 and bin_data[i-1] == 1:
            left_to_right = np.append(left_to_right,count_left)
            count_left = 0
 
 outfile_lefttoright = out_name+"_unbindingtimes.dat"
 outfile_righttoleft = out_name+"_bindingtimes.dat"

 np.savetxt(outfile_lefttoright,left_to_right,delimiter=',',fmt='%10.3f')
 np.savetxt(outfile_righttoleft,right_to_left,delimiter=',',fmt='%10.3f')

if __name__ == "__main__":
        main()

