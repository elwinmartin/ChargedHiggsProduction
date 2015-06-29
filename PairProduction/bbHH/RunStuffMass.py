'''
Created on Apr 2, 2013
Mass loop, woo!
@author: elwin
'''
open('gutg.dat', 'w').close()
import os
for c in range(0, 100):
    var = 165 + .15 * c 
    var = str(var)
    string = "14000 \n 16 \n 0 \n 5.0d0 \n 8000.d0 \n" + var + "d0" + "\n 0.254d0 \n 1"
 #   print string
    with open("gutg.dat", "a") as myfile:
        myfile.write(string)
    os.system("./gutg.x")
    with open("outputmass.dat", "r") as myfile2:
        output = myfile2.read()
    with open("datums.data","a") as myfile3:
        myfile3.write(output)
    os.system("rm gutg.dat")
