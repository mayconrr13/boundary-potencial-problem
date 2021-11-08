import numpy as np
from shapeFunctions import getShapeFunctionValueOnNode

class Element:
    def __init__(self, nodeList: list):
        self.nodeList = nodeList
        self.order = len(nodeList) - 1

    def getElementNodesRealCoordinates(self, geometricNodes):
        coordinates = np.zeros(len(self.nodeList), dtype=list)

        for i in range(len(self.nodeList)):
            coordinates[i] = geometricNodes[self.nodeList[i]]

        return coordinates

    def getDimensionlessPointsBasedOnGeometricCoordinates(self):
        dimensionlessPoints = np.zeros(len(self.nodeList), dtype=float)

        for i in range(len(self.nodeList)):
            if i == 0:
                dimensionlessPoints[i] = -1

            else:
                ksiValue = (2 / (len(self.nodeList) - 1)) * i - 1
                dimensionlessPoints[i] = ksiValue

        return dimensionlessPoints

    def handleElementType(self, duplicatedNodes):
        amountOfDuplicatedNodes = 0

        for i in range(len(self.nodeList)):
            if self.nodeList[i] in duplicatedNodes:
                amountOfDuplicatedNodes += 1

        if amountOfDuplicatedNodes == 2:
            return "discontinuous"
        elif amountOfDuplicatedNodes == 1:
            return "semicontinuous"
        else:
            return "continuous"

    def getDimensionlessPointsBasedOnElementContinuity(self, duplicatedNodes):
        dimensionlessPoints = self.getDimensionlessPointsBasedOnGeometricCoordinates()
        continuityBasedDimensionlessPoints = np.zeros(len(dimensionlessPoints))

        for j in range(len(self.nodeList)):
            if self.nodeList[j] in duplicatedNodes:
                discontinuousPoint = dimensionlessPoints[j]

                if j == 0:
                    nextPoint = dimensionlessPoints[j + 1]

                    ksiValue = discontinuousPoint + (nextPoint - discontinuousPoint) * 0.25

                    continuityBasedDimensionlessPoints[j] = ksiValue

                else:
                    previousPoint = dimensionlessPoints[j - 1]

                    ksiValue = discontinuousPoint + (previousPoint - discontinuousPoint) * 0.25

                    continuityBasedDimensionlessPoints[j] = ksiValue

            else:
                continuityBasedDimensionlessPoints[j] = dimensionlessPoints[j]

        return continuityBasedDimensionlessPoints

    def getAuxiliaryNodesCoordinates(self, duplicatedNodes, geometricNodes):
        dimensionlessPoints = self.getDimensionlessPointsBasedOnGeometricCoordinates()
        continuityBaseddimensionlessPoints = self.getDimensionlessPointsBasedOnElementContinuity(duplicatedNodes)
        auxiliaryNodes = []

        for i in range(len(continuityBaseddimensionlessPoints)):
            xCoordinate = 0
            yCoordinate = 0

            for j in range(len(self.nodeList)):
                xCoordinate += getShapeFunctionValueOnNode(
                    continuityBaseddimensionlessPoints[i], j, dimensionlessPoints) * geometricNodes[self.nodeList[j]][0]
                yCoordinate += getShapeFunctionValueOnNode(
                    continuityBaseddimensionlessPoints[i], j, dimensionlessPoints) * geometricNodes[self.nodeList[j]][1]

            auxiliaryNodes.append([xCoordinate, yCoordinate])

        return auxiliaryNodes

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)
