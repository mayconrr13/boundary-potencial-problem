def generateResultsFile(boundaryResults, internalResults):
    file = open('potentialProblem-' + str(1) + '.txt', "a")

    file.write('- - - - - - - - - - - - - - - - - - - - - -\n')
    file.write('Boundary Element Method - Potential Problem\n')
    file.write('- - - - - - - - - - - - - - - - - - - - - -\n\n')

    file.write('- - - - - - - - - - - - - -\n')
    file.write('Results on colocation nodes\n')
    file.write('- - - - - - - - - - - - - -\n')
    file.write('point, pontential, flow\n\n')

    for i in range(len(boundaryResults[0])):
        file.write('%.0f' % i + ', ' + '%.6f' % boundaryResults[0][i] + ', ' + '%.6f' % boundaryResults[1][i] + '\n')

    file.write('\n')
    file.write('- - - - - - - - - - - - - -\n')
    file.write('Results on internal points\n')
    file.write('- - - - - - - - - - - - - -\n')
    file.write('point, pontential, flowX, flowY\n\n')

    for i in range(len(internalResults)):
        file.write('I-%.0f' % i + ', ' + '%.6f' % internalResults[i][0] + ', ' + '%.6f' % internalResults[i][1][0] + ', ' + '%.6f' % internalResults[i][1][1] + '\n')
        
    file.close()

# print("Boundary Potential: ", boundaryResults[0])
# print("Boundary Flow: ", boundaryResults[1])

# print("Internal Points Potential: ", IntPotentialResults)
# print("Internal Points Flow: ", internalResults[1])