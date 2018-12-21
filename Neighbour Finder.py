# script that finds the nearest neighbours in a cluster to determine the average number of nearest neighbours, the average bond angle, the average ring size, the largest ring size, and the number of islands

import sys
from math import sqrt, acos, degrees

#inputFile = sys.argv[1]                                                     # finds name of file it acts upon
inputFile = 'tip3p_n7_g.snapshot.1490000.0.extxyz'

linesTraj = [line.rstrip('\n') for line in open(inputFile)]
noLines = len(linesTraj)
noAt = int(linesTraj[0])                                                    # finds the number of atoms in the cluster
noMol = int(noAt/3)
cellSize = float(((linesTraj[1].split(' '))[0].split('"'))[1])              # finds lattice dimensions
centre = cellSize/2

def coordMover(Prime, mol2, jList, centreDist):                             # moves cluster to centre
    x1old = float(jList[Prime][0])                                          # finds x coord of prime molecule
    y1old = float(jList[Prime][1])
    z1old = float(jList[Prime][2])
    x2old = float(jList[mol2][0])                                           # finds x coord of distance checking molecule
    y2old = float(jList[mol2][1])
    z2old = float(jList[mol2][2])

    x1 = x1old - centreDist[0]                                              # moves coords
    if x1 > cellSize:
        x1 = x1 - cellSize
    if x1 < 0:
        x1 = cellSize + x1
    y1 = y1old - centreDist[1]
    if y1 > cellSize:
        y1 = y1 - cellSize
    if y1 < 0:
        y1 = cellSize + y1
    z1 = z1old - centreDist[2]
    if z1 > cellSize:
        z1 = z1 - cellSize
    if z1 < 0:
        z1 = cellSize + z1
    x2 = x2old - centreDist[0]
    if x2 > cellSize:
        x2 = x2 - cellSize
    if x2 < 0:
        x2 = cellSize + x2
    y2 = y2old - centreDist[1]
    if y2 > cellSize:
        y2 = y2 - cellSize
    if y2 < 0:
        y2 = cellSize + y2
    z2 = z2old - centreDist[2]
    if z2 > cellSize:
        z2 = z2 - cellSize
    if z2 < 0:
        z2 = cellSize + z2
    return x1, y1, z1, x2, y2, z2

def centreDrift(lines) :                                                    # finds distance cluster has drifted from centre
    lineSplit = lines[2].split()
    xcoord = float(lineSplit[1])                                            # finds x coordinate
    ycoord = float(lineSplit[2])                                            # finds y coordinate
    zcoord = float(lineSplit[3])                                            # finds z coordinate
    xDist = xcoord - centre
    yDist = ycoord - centre
    zDist = zcoord - centre

    return xDist, yDist, zDist


def distanceCalc(Prime, mol2, jList, centreDist):                           # finds the distance between molecules
    newCoords = coordMover(Prime, mol2, jList, centreDist)

    x1 = newCoords[0]                                                       # finds x coord of prime molecule
    y1 = newCoords[1]
    z1 = newCoords[2]
    x2 = newCoords[3]                                                       # finds x coord of distance checking molecule
    y2 = newCoords[4]
    z2 = newCoords[5]
    xcomp = (x2 - x1) ** 2
    ycomp = (y2 - y1) ** 2
    zcomp = (z2 - z1) ** 2
    distance = sqrt(xcomp + ycomp + zcomp)                                  # finds distance between molecules
    return distance


def angleFinder1(lines, centreDist):                                        # calculates angle for single config
    jList = [[] for i in range(noMol)]                                      # creates a list of molecules (blank)
    angleList = []
    neighbourListPerm = []

    n = 2                                                                   # iterates through molecules to make list of coords
    n2 = 0                                                                  # increments through where coords are printed
    while n < noAt + 2:
        lineSplit = lines[n].split()
        xcoord = lineSplit[1]                                               # finds x coordinate
        ycoord = lineSplit[2]                                               # finds y coordinate
        zcoord = lineSplit[3]                                               # finds z coordinate
        jList[n2].append(xcoord)                                            # adds coordinates as list to in list of molecules
        jList[n2].append(ycoord)
        jList[n2].append(zcoord)
        n = n + 3
        n2 = n2 + 1

    Prime = 0                                                               # iterates through molecules to find angles
    graph = {}
    while Prime < noMol:
        neighbourList = []
        mol2 = 0
        while mol2 < noMol:
            if mol2 - Prime != 0:
                if distanceCalc(Prime, mol2, jList, centreDist) <= 3.1:
                    neighbourList.append(mol2)
            mol2 = mol2 + 1
        noNeighbour = len(neighbourList)                                    # finds the number of neighbours
        neighbourListPerm.append(noNeighbour)
        graph[Prime] = neighbourList

        if noNeighbour > 2:                                                 # checks if there's enough neighbours to find angle
            secondIt = 0                                                    # loops through first set of neighbours
            while secondIt < noNeighbour:
                second = neighbourList[secondIt]                            # finds molecule for angle
                thirdIt = secondIt + 1
                while thirdIt < noNeighbour:
                    third = neighbourList[thirdIt]                          # finds molecule for angle
                    dist1 = distanceCalc(Prime, second, jList, centreDist)  # finds first triangle length
                    dist2 = distanceCalc(Prime, third, jList, centreDist)   # finds second triangle length
                    dist3 = distanceCalc(second, third, jList, centreDist)  # finds third triangle length
                    cosAngle = (dist1 ** 2 + dist2 ** 2 - dist3 ** 2) / (2 * dist1 * dist2)
                    angle = degrees(acos(cosAngle))
                    thirdIt = thirdIt + 1
                    angleList.append(angle)
                secondIt = secondIt + 1
        Prime = Prime + 1

    noAngles = len(angleList)                                                       # finds number of angles
    angleSum = 0
    n = 0
    while n < noAngles:                                                             # summs all the angles
        angleSum = angleSum + angleList[n]
        n = n + 1
    if noAngles > 0:                                                                # tests if there are any angles
        angleAvg = angleSum / noAngles                                              # finds average angle
    else:
        angleAvg = 0

    noNeighbours = len(neighbourListPerm)                                           # finds number of angles
    neighbourSum = 0
    n = 0
    while n < noNeighbours:
        neighbourSum = neighbourSum + neighbourListPerm[n]
        n = n +1
    if noNeighbours > 0:                                                            # tests if there are any angles
        aveNN = neighbourSum / noMol                                                # finds average angle
    else:
        aveNN = 0

    return angleAvg, aveNN, graph


def find_path(graph, start, end, path=[]):                                          # finds paths within a graph
    path = path + [start]
    if start == end:
        return [path]
    if not start in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_path(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


def sort_paths2(foundPaths, noPaths):
    shortestPathLength = 999998                                                     # tracks shortest path
    secondShortestPathLength = 999999                                               # tracks second shortest path
    pathIteration = 0                                                               # iterates through all paths to find shortest
    while pathIteration < noPaths:                                                  # tracks whether path is forbidden
        currentPath = foundPaths[pathIteration]
        currentPathLength = len(currentPath)
        if currentPathLength < secondShortestPathLength:
            if currentPathLength < shortestPathLength:
                shortestPathList = []
                secondShortestPathList = []
                shortestPathList.append(currentPath)
                secondShortestPathList.append(currentPath)
                secondShortestPathLength = shortestPathLength
                shortestPathLength = currentPathLength
            else:
                if currentPathLength == shortestPathLength:
                    shortestPathList.append(currentPath)
                else:
                    secondShortestPathList = []
                    secondShortestPathList.append(currentPath)
                    secondShortestPathLength = currentPathLength
        pathIteration = pathIteration + 1
    noShortestPaths = len(shortestPathList)
    if noShortestPaths == 1 and len(foundPaths) > 1:
        shortestPathList.append(secondShortestPathList[0])
    return shortestPathList


def sort_paths(foundPaths, noPaths, forbiddenPaths):
    shortestPathLength = 999999                                                     # tracks shortest path
    pathIteration = 0                                                               # iterates through all paths to find shortest
    noForbiddenPaths = len(forbiddenPaths)
    shortestPathList =[]
    while pathIteration < noPaths:
        forbidden = 0                                                               # tracks whether path is forbidden
        currentPath = foundPaths[pathIteration]
        forbiddenPathsIt = 0
        while forbiddenPathsIt < noForbiddenPaths:                                  # iterates through forbidden paths to check if path matches any in list
            if forbiddenPaths[forbiddenPathsIt] == currentPath:
                forbidden = 1
                forbiddenPathsIt = noForbiddenPaths
            else:
                forbiddenPathsIt = forbiddenPathsIt + 1
        if forbidden == 1:                                                          # if path is forbidden: skip this path
            pathIteration = pathIteration + 1
        else:
            currentPathLength = len(currentPath)
            if currentPathLength < shortestPathLength:
                shortestPathList = []
                shortestPathList.append(currentPath)
                shortestPathLength = currentPathLength
            else:
                if currentPathLength == shortestPathLength:
                    shortestPathList.append(currentPath)
            pathIteration = pathIteration + 1
    return shortestPathList


def ring_processor(allRingsTotal):
    if len(allRingsTotal) > 0:
        largestRingSize = 0
        allProcessed = 0                                                                # keeps track of whether all rings have been processed
        allRingsFinal = []
        while allProcessed == 0:
            noRings = len(allRingsTotal)
            comparisonRingIt = 1
            primeRing = allRingsTotal[0]
            primeRingLen = len(primeRing)
            if primeRingLen > largestRingSize:
                largestRingSize = primeRingLen
            allRingsFinal.append(primeRing)
            while comparisonRingIt < noRings:                                            # iterates through all other rings to check for duplicates
                comparisonRing = allRingsTotal[comparisonRingIt]
                if primeRingLen == len(comparisonRing):                                 # if prime ring and comparison ring are same length, check if they contain the same atom
                    commonAtoms = 0                                                     # keeps track of number of common atoms in both lists
                    primeRingElementIt = 0
                    while primeRingElementIt < primeRingLen:
                        primeRingElement = primeRing[primeRingElementIt]
                        if primeRingElement in comparisonRing:
                            commonAtoms = commonAtoms + 1
                        primeRingElementIt = primeRingElementIt + 1
                    if commonAtoms == primeRingLen:
                        del allRingsTotal[comparisonRingIt]
                    else:
                        comparisonRingIt = comparisonRingIt + 1
                else:
                    comparisonRingIt = comparisonRingIt + 1
                noRings = len(allRingsTotal)
            del allRingsTotal[0]
            noRings = len(allRingsTotal)
            if noRings == 0:
                allProcessed = 1

        noFinalRings = len(allRingsFinal)
        totalRingSize = 0
        allRingsFinalIt = 0
        while allRingsFinalIt < noFinalRings:                                           # iterates through the final list of rings and sums the ring size
            totalRingSize = totalRingSize + len(allRingsFinal[allRingsFinalIt])
            allRingsFinalIt = allRingsFinalIt + 1
        averageRingSize = totalRingSize / float(noFinalRings)
    else:
        allRingsFinal = []
        averageRingSize = 0
        largestRingSize = 0

    return averageRingSize, largestRingSize


config = 0                                                                          # iterates through configurations
counter = 0
while counter < noLines:
    allRingsTotal = []                                                              # creates a list of rings which are append as they're found
    lineFirst = int(linesTraj[counter])                                             # finds value for first row of config (should be == noMol)
    lines = linesTraj[counter:counter+((noAt)+2)]                                   # creates list for given configuration
    centreDist = centreDrift(lines)                                                 # finds distance of cluster from centre
    angleAvg2 = angleFinder1(lines, centreDist)                                     # finds average angle of config
    graph = angleAvg2[2]


    totalRingSize = 0                                                               # keeps track of the total loop length
    atomNumber = 0                                                                  # finds loops for atom in position[atomNumber]
    while atomNumber < noMol:
        graphs = {}                                                                 # creates new graph without subject atom (blank at first)
        if atomNumber == 0:                                                         # ensures the following process includes all atoms except the subject atom
            otherAtoms = 1
        else:
            otherAtoms = 0
        while otherAtoms < noMol:                                                   # removes atom in position[atomNumber] so it is not included in loops
            atomsNeighList = graph[otherAtoms]
            atomsNeighlistLen = len(atomsNeighList)
            atomsNeighListMod = []                                                  # list of neighbours with subject atom removed
            for neighbour in atomsNeighList:
                if neighbour != str(atomNumber):
                    atomsNeighListMod.append(neighbour)
            graphs[otherAtoms] = atomsNeighListMod

            otherAtoms = otherAtoms + 1
            if otherAtoms == atomNumber:
                graphs[otherAtoms] = []
                otherAtoms = otherAtoms + 1

        atomsNeighList = graph[atomNumber]
        noNeighbours = len(atomsNeighList)
        allPathsList = []
        allShortestPathsList = []
        if noNeighbours < 2:                                                        # tests if atom has enough neighbours to form loop
            actualpath = []

        if noNeighbours == 2:                                                       # atoms with 2 neighbours must be handled
            neighPathComMax = 1
            startPoint = atomsNeighList[0]
            endPoint = atomsNeighList[1]
            foundPaths = find_path(graphs, startPoint, endPoint)
            allPathsList.append(foundPaths)
            noPaths = len(foundPaths)
            if noPaths > 0:
                sort_paths_output = sort_paths2(foundPaths, noPaths)
                shortestPathList = sort_paths_output
                allShortestPathsList.append(shortestPathList)
            else:
                shortestPathList = []

        else:
            neighSPos = 0                                                           # position in list of starting neighbour
            neighPathCom = 0                                                        # keeps track of the combinations of neighbours that a path has been found between
            while neighSPos < noNeighbours:
                startPoint = atomsNeighList[neighSPos]
                neighEPos = neighSPos + 1                                           # position in list of ending neighbour
                while neighEPos < noNeighbours:
                    endPoint = atomsNeighList[neighEPos]
                    forbiddenPaths = []
                    foundPaths = find_path(graphs, startPoint, endPoint)            # finds shortest path between two points
                    allPathsList.append(foundPaths)
                    noPaths = len(foundPaths)                                       # finds the number of paths found
                    if noPaths > 0:
                        sort_paths_output = sort_paths(foundPaths, noPaths, forbiddenPaths)
                        shortestPathList = sort_paths_output
                        allShortestPathsList.append(shortestPathList)
                    else:
                        allShortestPathsList.append([])

                    neighPathCom = neighPathCom + 1
                    neighEPos = neighEPos + 1

                neighSPos = neighSPos + 1

            neighPathComMax = neighPathCom                                         # finds the total number of neighbour combinations
            neighPathCom = 0
            shortestPathLength = 999998
            secondShortestPathLength = 999999
            allPathsWithAltRoutes = []
            loopbackCheckThreshhold = 999999
            while neighPathCom < neighPathComMax:                                  # iterates through all ring sizes to find shortest to help determine if there are any loopbacks
                currentPathSet = allShortestPathsList[neighPathCom]
                currentPathSetSize = len(currentPathSet)
                currentPathSetIt = 0
                while currentPathSetIt < currentPathSetSize:                        # iterates through all paths for one combination in order to add them to list for loopback testing
                    currentPath = currentPathSet[currentPathSetIt]
                    currentPathReverse = currentPath[::-1]                          # finds reverse of current path
                    allPathsWithAltRoutes.append(currentPath)
                    allPathsWithAltRoutes.append(currentPathReverse)
                    currentPathLength = len(currentPath)
                    if currentPathLength <= secondShortestPathLength:
                        if currentPathLength < shortestPathLength:
                            secondShortestPathLength = shortestPathLength
                            shortestPathLength = currentPathLength
                        else:
                            secondShortestPathLength = currentPathLength

                    currentPathSetIt = currentPathSetIt + 1

                neighPathCom = neighPathCom + 1

            loopbackCheckThreshhold = shortestPathLength + secondShortestPathLength - 1
            if loopbackCheckThreshhold < 4:                                         # 4 is smallest possible loopback
                loopbackCheckThreshhold = 999999

            numberOfPaths = len(allPathsWithAltRoutes)
            neighPathCom = 0
            forbiddenPaths = []
            while neighPathCom < neighPathComMax:
                currentPathSet = allShortestPathsList[neighPathCom]
                if len(currentPathSet) == 0:
                    neighPathCom = neighPathCom + 1
                else:
                    currentPathLength = len(currentPathSet[0])
                    if currentPathLength >= loopbackCheckThreshhold and len(forbiddenPaths) == 0:
                        expandedPathIt = 0
                        while expandedPathIt < numberOfPaths:                                           # iterates through all alternate routes to see if a loopback can be made
                            forbiddenPathAddition = allPathsWithAltRoutes[expandedPathIt]
                            forbiddenPaths.append(forbiddenPathAddition)
                            noForbiddenPaths = len(forbiddenPaths)
                            forbiddenPathsIt = 0
                            while forbiddenPathsIt < noForbiddenPaths:                                  # iterates through existing forbidden paths and appends new path to it
                                currentForbiddenPath = forbiddenPaths[forbiddenPathsIt]
                                currentForbiddenPathLength = len(currentForbiddenPath)
                                if currentForbiddenPathLength <= noMol:                                 # checks path isn't already larger than number of molecules in system
                                    currentForbiddenPathEnd = currentForbiddenPath[currentForbiddenPathLength-1]
                                    if forbiddenPathAddition[0] == currentForbiddenPathEnd:             # checks to see if new path starts with the same node as the end of the last path
                                        forbiddenPathAdditionNew = forbiddenPathAddition.copy()
                                        del forbiddenPathAdditionNew[0]
                                        newForbiddenPath = currentForbiddenPath + forbiddenPathAdditionNew
                                        forbiddenPaths.append(newForbiddenPath)
                                forbiddenPathsIt = forbiddenPathsIt + 1
                            expandedPathIt = expandedPathIt + 1

                    if currentPathLength >= loopbackCheckThreshhold:                                    # only checks for loopbacks if path is longer than threshold
                        shortestPaths = allShortestPathsList[neighPathCom]
                        degeneracy = len(shortestPaths)
                        noForbiddenPaths = len(forbiddenPaths)
                        shortestPathsIt = 0
                        while shortestPathsIt < degeneracy:                                             # iterates through all paths to remove any loopbacks
                            forbiddenPathsIt = 0
                            while forbiddenPathsIt < noForbiddenPaths:
                                if shortestPaths[0] == forbiddenPaths[forbiddenPathsIt]:
                                    del shortestPaths[0]
                                    forbiddenPathsIt = noForbiddenPaths
                                else:
                                    forbiddenPathsIt = forbiddenPathsIt + 1
                            shortestPathsIt = shortestPathsIt + 1
                        allShortestPathsList[neighPathCom] = shortestPaths                              # updates list of shortest paths to include new paths
                        currentPathSet = allShortestPathsList[neighPathCom]

                    degeneracy = len(currentPathSet)
                    currentPathIt = 0
                    while currentPathIt < degeneracy:
                        currentPathSet[currentPathIt].append(atomNumber)
                        allRingsTotal.append(currentPathSet[currentPathIt])
                        currentPathIt = currentPathIt + 1

                    neighPathCom = neighPathCom + 1

        atomNumber = atomNumber + 1

    ring_processor_out = ring_processor(allRingsTotal)
    averageRing = ring_processor_out[0]
    largestRingSize = ring_processor_out[1]

    atomNumber = 0                                                                                  # finds loops for atom in position[atomNumber]
    noIslands = 1                                                                                   # keeps track of the number of islands
    while atomNumber < noMol:
        isolatedAt = []
        otherAtom = atomNumber + 1
        while otherAtom < noMol:
            actualpath = find_path(graph, atomNumber, otherAtom)                                    # finds shortest path between two points
            if len(actualpath) == 0:
                isolatedAt.append(otherAtom)
            otherAtom = otherAtom + 1
        noIsolatedAt = len(isolatedAt)
        if noIsolatedAt > 0:
            atomNumber = int(isolatedAt[0])
            noIslands = noIslands + 1
        else:
            atomNumber = noAt


    newLine = str(angleAvg2[1]) + ' ' + str(angleAvg2[0]) + ' ' + str(noIslands) + ' ' + str(averageRing) + ' ' + str(largestRingSize) # generates row of data

    print(newLine)
    config = config + 1
    counter = config * (noAt+2)
