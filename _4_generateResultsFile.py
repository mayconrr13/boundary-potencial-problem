def generateResultsFile(PotentialVector, FlowVector, internalPoints, IntPotentialResults, IntFlowResults):
    boundaryResults = [PotentialVector, FlowVector]
    internalResults = []
    for i in range(len(internalPoints)):
        internalResults.append([IntPotentialResults[i], []])

        for j in range(2):
            internalResults[i][1].append(IntFlowResults[2 * i + j])

    file = open('potentialProblem-' + str(1) + '.txt', "a")

    file.write('- - - - - - - - - - - - - - - - - - - - - -\n')
    file.write('Boundary Element Method - Potential Problem\n')
    file.write('- - - - - - - - - - - - - - - - - - - - - -\n\n')

    file.write('- - - - - - - - - - - - - -\n')
    file.write('Results on colocation nodes\n')
    file.write('- - - - - - - - - - - - - -\n')
    file.write('point, pontential, flow\n\n')

    for i in range(len(boundaryResults[0])):
        file.write('%.0f' % i + ', ' + '%.3f' % boundaryResults[0][i] + ', ' + '%.3f' % boundaryResults[1][i] + '\n')

    file.write('\n')
    file.write('- - - - - - - - - - - - - -\n')
    file.write('Results on internal points\n')
    file.write('- - - - - - - - - - - - - -\n')
    file.write('point, pontential, flowX, flowY\n\n')

    for i in range(len(internalResults)):
        file.write('I-%.0f' % i + ', ' + '%.3f' % internalResults[i][0] + ', ' + '%.3f' % internalResults[i][1][0] + ', ' + '%.3f' % internalResults[i][1][1] + '\n')
        
    file.close()

    print("|||||||||||||")
    print("||| E N D |||")
    print("|||||||||||||")

# print("Boundary Potential: ", boundaryResults[0])
# print("Boundary Flow: ", boundaryResults[1])

# print("Internal Points Potential: ", IntPotentialResults)
# print("Internal Points Flow: ", internalResults[1])