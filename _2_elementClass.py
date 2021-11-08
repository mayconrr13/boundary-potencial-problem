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

    def getAdimensionalPointsBasedOnGeometricCoordinates(self):
        adimensionalPoints = np.zeros(len(self.nodeList), dtype=float)

        for i in range(len(self.nodeList)):
            if i == 0:
                adimensionalPoints[i] = -1

            else:
                ksiValue = (2 / (len(self.nodeList) - 1)) * i - 1
                adimensionalPoints[i] = ksiValue

        return adimensionalPoints

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

    def getAdimensionalPointsBasedOnElementContinuity(self, duplicatedNodes):
        adimensionalPoints = self.getAdimensionalPointsBasedOnGeometricCoordinates()
        continuityBasedAdimensionalPoints = np.zeros(len(adimensionalPoints))

        for j in range(len(self.nodeList)):
            if self.nodeList[j] in duplicatedNodes:
                discontinuousPoint = adimensionalPoints[j]

                if j == 0:
                    nextPoint = adimensionalPoints[j + 1]

                    ksiValue = discontinuousPoint + \
                        (nextPoint - discontinuousPoint) * 0.25

                    continuityBasedAdimensionalPoints[j] = ksiValue

                else:
                    previousPoint = adimensionalPoints[j - 1]

                    ksiValue = discontinuousPoint + \
                        (previousPoint - discontinuousPoint) * 0.25

                    continuityBasedAdimensionalPoints[j] = ksiValue

            else:
                continuityBasedAdimensionalPoints[j] = adimensionalPoints[j]

        return continuityBasedAdimensionalPoints

    def getColocationNodesCoordinates(self, duplicatedNodes, geometricNodes):
        adimensionalPoints = self.getAdimensionalPointsBasedOnGeometricCoordinates()
        continuityBasedAdimensionalPoints = self.getAdimensionalPointsBasedOnElementContinuity(
            duplicatedNodes)
        # colocationNodes = np.zeros(len(self.nodeList), dtype=list)
        colocationNodes = []
        # colocationNodes = np.array([], dtype=int)

        for i in range(len(continuityBasedAdimensionalPoints)):
            xCoordinate = 0
            yCoordinate = 0

            for j in range(len(self.nodeList)):
                xCoordinate += getShapeFunctionValueOnNode(
                    continuityBasedAdimensionalPoints[i], j, adimensionalPoints) * geometricNodes[self.nodeList[j]][0]
                yCoordinate += getShapeFunctionValueOnNode(
                    continuityBasedAdimensionalPoints[i], j, adimensionalPoints) * geometricNodes[self.nodeList[j]][1]

            # colocationNodes = np.append(colocationNodes, ([xCoordinate, yCoordinate]))
            colocationNodes.append([xCoordinate, yCoordinate])
            # colocationNodes[i] = [xCoordinate, yCoordinate]

        # duplicatedNodes = np.array(colocationNodes, dtype=list)

        return colocationNodes

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)
