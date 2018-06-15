# script that takes the fourth row of a series 'csv' files with spaces instead of commas and then stores those rows as a csv file
import sys

if len(sys.argv) < 5:
    print("Error! you must give me an input of the following form: [splitter2.py] [fileNForm] [noFiles] [fileExt] [noTemp]")
    sys.exit(0)

fileNForm = sys.argv[1]                                                    # finds file name
noFiles =   int(sys.argv[2])
fileExt =   sys.argv[3]
noTemp =    int(sys.argv[4])                                               # finds number of discrete temperature values
outputFile = ''.join([fileNForm, '.csv'])                                  # creates name of output
outputFileOpen = open(outputFile, 'w')


jList = [[] for i in range(noTemp)]
n = 0                                                                      # n = number of the file
while n < noFiles:
    inputFile = ''.join([fileNForm, str(n), fileExt])                      # creates name of input file
    lines = [line.rstrip('\n') for line in open(inputFile)]                # splits lines into individual variables

    x=2                                                                    # x = row of file
    while x< noTemp + 2:                                                   # replaces all spaces with commas
        lineSplit = lines[x].split(' ')
        heatCap = lineSplit[3]
        jList[x-2].append(heatCap)
        x=x+1
    n=n+1

noTempCount = 0
while noTempCount < noTemp:
        activeRowList = jList[noTempCount]
        newline = ', '.join(activeRowList)
        outputFileOpen.write(newline+'\n')
        noTempCount = noTempCount+1

outputFileOpen.close()
