# script that takes the fourth row of a series 'csv' files with spaces instead of commas and then stores those rows as a csv file

fileNForm = input('\nEnter generic form of file name: ')                   # finds file name
noFiles = input('\nEnter number of files: ')
fileExt = input('\nEnter file extension: ')
noTemp = input('\nEnter number of discrete temperature values: ')          # finds number of discrete temperature values
outputFile = ''.join([fileNForm, '.csv'])                                  # creates name of output
outputFileOpen = open(outputFile, 'w')


jList = [[] for i in range(int(noTemp))]
n = 0                                                                      # n = number of the file
while n < int(noFiles):
    inputFile = ''.join([fileNForm, str(n), fileExt])                      # creates name of input file
    lines = [line.rstrip('\n') for line in open(inputFile)]                # splits lines into individual variables

    x=2                                                                    # x = row of file
    while x<int(noTemp)+2:                                                 # replaces all spaces with commas
        lineSplit = lines[x].split(' ')
        heatCap = lineSplit[3]
        jList[x-2].append(heatCap)
        x=x+1
    n=n+1

noTempCount = 0
while noTempCount < int(noTemp):
        activeRowList = jList[noTempCount]
        newline = ', '.join(activeRowList)
        outputFileOpen.write(newline+'\n')
        noTempCount = noTempCount+1

outputFileOpen.close()
