# pulls heat capacity and temperature data out of file containing data seperated by spaces

userInput = input('\n[File Name] [No. Files] [File Ext] [ No. Temps]')  # accepts input for specifications
inputList = userInput.split(' ')
fileNForm = inputList[0]                                                # file name
noFiles = inputList[1]                                                  # number of files
fileExt = inputList[2]                                                  # file extension
noTemp = inputList[3]                                                   # number of discrete temperature values
outputFile = ''.join([fileNForm, '.csv'])                               # creates name of output
outputFileOpen = open(outputFile, 'w')


jList = [[] for i in range(int(noTemp))]                                # creates a list of rows (blank)

# generates first row: Temperatures
inputFile = ''.join([fileNForm, str(0), fileExt])                       # creates name of input file
lines = [line.rstrip('\n') for line in open(inputFile)]                 # splits lines into individual variables
x = 2                                                                   # x = row of file
while x < int(noTemp) + 2:                                              # replaces all spaces with commas
    lineSplit = lines[x].split(' ')
    temp = lineSplit[0]
    jList[x - 2].append(temp)
    x = x + 1

# iterates through files
n = 0                                                                   # n = number of the file
while n < int(noFiles):
    inputFile = ''.join([fileNForm, str(n), fileExt])                   # creates name of input file
    lines = [line.rstrip('\n') for line in open(inputFile)]             # splits lines into individual variables

    # iterates through row of file
    x=2                                                                 # x = row of file
    while x<int(noTemp)+2:
        lineSplit = lines[x].split(' ')                                 # separates values in row by spaces
        heatCap = lineSplit[3]
        jList[x-2].append(heatCap)
        x=x+1
    n=n+1

# writes rows to .csv
noTempCount = 0
while noTempCount < int(noTemp):
        activeRowList = jList[noTempCount]
        newline = ', '.join(activeRowList)
        outputFileOpen.write(newline+'\n')
        noTempCount = noTempCount+1

outputFileOpen.close()
