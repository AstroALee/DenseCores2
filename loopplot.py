import sys
import os



for i in range(int(sys.argv[1])):
	fname = "Output_" + str(i)
	os.system("python oneDplot.py " + fname + " " + str(sys.argv[2]) + " n z 0 17 32")
        os.system("mv OneDPlot.pdf OneDPlot_" + str(sys.argv[2]) + "_" + str(i) + ".pdf")
	os.system("python contour.py " + fname + " " + str(sys.argv[2]) + " n n 6 n 6")
	os.system("mv Contour.pdf Contour_" + str(sys.argv[2]) + "_" + str(i) + ".pdf")


