import scipy.stats.mstats as mstats
import numpy
import os 
import glob

MCLAB = os.getenv('MCLAB')
os.chdir(MCLAB + "/arbeitman_fru_network/data/for_wiggles")

OUT = open('./background.txt', 'w')

OUT.write("rep,sum,mean,q0.5,q0.75,q0.9\n")

for myfile in glob.glob("*.csv"):

    count = []
    currep = open(myfile, 'r')
    next(currep)

    for line in currep:
        count.append(float(line.split(',')[2]))

    mysum = sum(count)
    mymean = round(numpy.mean(count) ,2)
    myquant = ','.join(map(str,mstats.mquantiles(count, prob=(.50, .75, .90))))

    myout = ','.join(map(str,(myfile,mysum,mymean,myquant)))

    OUT.write(myout + '\n')

OUT.close()
