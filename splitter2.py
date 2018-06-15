fileNForm = input('\nEnter generic form of file name: ')                    # finds file name
noFiles = input('\nEnter number of files: ')
fileExt = input('\nEnter file extension: ')
outputFile = ''.join([fileNForm, '.csv'])                                   # creates name of output file
outputFileOpen = open(outputFile, 'w')

jList = []
n = 0
while n < int(noFiles):
    inputFile = ''.join([fileNForm, n, fileExt])                       # creates name of input file
    lines = [line.rstrip('\n') for line in open(inputFile)]                 # splits lines into individual variables
    jList.append([lines[3]])


    x=0
    while x<203:                                                            # replaces all spaces with commas
        lineSplit = lines[x].split(' ')
        heatCap = lineSplit[3]
        jList.append([heatCap])
        x=x+1
    n=n+1


outputFileOpen.close()
