import math
import numpy as np
from shapeFunctions import getShapeFunctionValueOnNode
from tangentNormalAndJacobian import getPointsPropertiesOnElement

integrationPoints = [-0.9815606342467192506906, -0.9041172563704748566785, -0.769902674194304687037, -0.5873179542866174472967, -0.3678314989981801937527, -0.1252334085114689154724, 0.1252334085114689154724, 0.3678314989981801937527, 0.5873179542866174472967, 0.7699026741943046870369, 0.9041172563704748566785, 0.9815606342467192506906]
weights = [0.0471753363865118271946, 0.1069393259953184309603, 0.1600783285433462263347, 0.2031674267230659217491, 0.233492536538354808761, 0.2491470458134027850006, 0.2491470458134027850006, 0.233492536538354808761, 0.203167426723065921749, 0.160078328543346226335, 0.1069393259953184309603, 0.0471753363865118271946]

# element, value
# u = [[1, 0], [3, 10]]
# q = [[0, 0], [2, 0]]
# node, value
u = [[5, 0], [6, 0], [9, 10], [10, 10]]
q = [[0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [7, 0], [8, 0]]

geometricNodes = [
    [2, 0],
    [3.5, 0],
    [5, 0],
    [6.5, 0],
    [8, 0],
    [8, 0],
    [8, 4],
    [8, 4],
    [2, 4],
    [2, 4],
    [2, 0]
]

elements = [
    [0,1],
    [1,2],
    [2,3],
    [3,4],
    [5,6],
    [7,8],
    [9,10]
]

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

        # return amountOfDuplicatedNodes
    
    def getAdimensionalPointsBasedOnElementContinuity(self, duplicatedNodes):
        adimensionalPoints = self.getAdimensionalPointsBasedOnGeometricCoordinates()
        continuityBasedAdimensionalPoints = np.zeros(len(adimensionalPoints))      

        for j in range(len(self.nodeList)):
            if self.nodeList[j] in duplicatedNodes:
                discontinuousPoint = adimensionalPoints[j]

                if j == 0:
                    nextPoint = adimensionalPoints[j + 1]

                    ksiValue = discontinuousPoint + (nextPoint - discontinuousPoint) * 0.25
                    
                    continuityBasedAdimensionalPoints[j] = ksiValue
                
                else:
                    previousPoint = adimensionalPoints[j - 1]

                    ksiValue = discontinuousPoint + (previousPoint - discontinuousPoint) * 0.25
                    
                    continuityBasedAdimensionalPoints[j] = ksiValue

            else:
                continuityBasedAdimensionalPoints[j] = adimensionalPoints[j]
        
        return continuityBasedAdimensionalPoints

    def getColocationNodesCoordinates(self, duplicatedNodes, geometricNodes):
        adimensionalPoints = self.getAdimensionalPointsBasedOnGeometricCoordinates()
        continuityBasedAdimensionalPoints = self.getAdimensionalPointsBasedOnElementContinuity(duplicatedNodes)
        # colocationNodes = np.zeros(len(self.nodeList), dtype=list)
        colocationNodes = []
        # colocationNodes = np.array([], dtype=int)

        for i in range(len(continuityBasedAdimensionalPoints)):
            xCoordinate = 0
            yCoordinate = 0

            for j in range(len(self.nodeList)):
                xCoordinate += getShapeFunctionValueOnNode(continuityBasedAdimensionalPoints[i], j, adimensionalPoints) * geometricNodes[self.nodeList[j]][0]
                yCoordinate += getShapeFunctionValueOnNode(continuityBasedAdimensionalPoints[i], j, adimensionalPoints) * geometricNodes[self.nodeList[j]][1]

            # colocationNodes = np.append(colocationNodes, ([xCoordinate, yCoordinate]))
            colocationNodes.append([xCoordinate, yCoordinate])
            # colocationNodes[i] = [xCoordinate, yCoordinate]

        # duplicatedNodes = np.array(colocationNodes, dtype=list)

        return colocationNodes

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

elementsList = np.zeros(len(elements), dtype=list)
for el in range(len(elements)):
    elementsList[el] = Element(elements[el])


# process
def getDuplicatedNodes(nodeList: list):
    # newList = []
    duplicatedNodes = np.array([], dtype=int)

    for i in range (len(nodeList)):
        nodeToBeChecked = nodeList[i]
        nodeList.count(nodeToBeChecked)

        if nodeList.count(nodeToBeChecked) > 1:
            # duplicatedNodes.append(i)
            duplicatedNodes = np.append(duplicatedNodes, i)
    
    return duplicatedNodes

def generateColocationMesh(elementsList, duplicatedNodes, geometricNodes):
    colocationMesh = []

    for i in range(len(elementsList)):
        elementColocationCoordinates = elementsList[i].getColocationNodesCoordinates(duplicatedNodes, geometricNodes)
        
        for j in range(len(elementColocationCoordinates)):            
            if elementColocationCoordinates[j] not in colocationMesh:
                colocationMesh.append(elementColocationCoordinates[j])

    # colocationMesh = np.array(colocationMesh, dtype=list)   
    return colocationMesh

def getSourcePoints(duplicatedNodes, geometricNodes):
    sourcePointsList = []

    for i in range(len(elementsList)):
        elementColocationNodes = elementsList[i].getColocationNodesCoordinates(duplicatedNodes, geometricNodes)  
        adimentionalPoints = elementsList[i].getAdimensionalPointsBasedOnGeometricCoordinates()
        nodesProperties = getPointsPropertiesOnElement(elementsList[i].nodeList, elementColocationNodes, adimentionalPoints)
        
        for k in range(len(elementsList[i].nodeList)):
            xCoordinate = nodesProperties[2][k][0] * 1 + elementColocationNodes[k][0]
            yCoordinate = nodesProperties[2][k][1] * 1 + elementColocationNodes[k][1]

            sourcePointCoordinates = [xCoordinate, yCoordinate]

            if sourcePointsList.count(sourcePointCoordinates) < 1:
                sourcePointsList.append(sourcePointCoordinates)

    sourcePointsList = np.array(sourcePointsList)

    return sourcePointsList     

duplicatedNodes = getDuplicatedNodes(geometricNodes)
colocationMesh = generateColocationMesh(elementsList, duplicatedNodes, geometricNodes)
# sourcePoints = getSourcePoints(duplicatedNodes, geometricNodes)
sourcePoints = colocationMesh

def getIntegrationPointCoordinates(points: float, elementNodes: list, adimentionalPoints: list):
    elementPointsCoordinates = []
    xComponent = 0
    yComponent = 0

    for j in range(len(elementNodes)):
        xComponent += getShapeFunctionValueOnNode(points, j, adimentionalPoints) * elementNodes[j][0]
        yComponent += getShapeFunctionValueOnNode(points, j, adimentionalPoints) * elementNodes[j][1]
    
    elementPointsCoordinates = [xComponent, yComponent]

    elementPointsCoordinates = np.array(elementPointsCoordinates)

    return elementPointsCoordinates

def getIntegrationPointsCoordinatesPerElement(element: list, geometricNodes: list):
    nodesCoordinatesPerElement = []
    elementNodes = element.getElementNodesRealCoordinates(geometricNodes)
    adimentionalNodes = element.getAdimensionalPointsBasedOnGeometricCoordinates()

    for ip in range(len(integrationPoints)):
        coordinate = getIntegrationPointCoordinates(integrationPoints[ip], elementNodes, adimentionalNodes)
        nodesCoordinatesPerElement.append(coordinate)
    
    return nodesCoordinatesPerElement

def getRadius(sourcePoint: list, integrationPointCoordinates: list):
    
    xComponent = integrationPointCoordinates[0] - sourcePoint[0]
    yComponent = integrationPointCoordinates[1] - sourcePoint[1]
    radius = (xComponent ** 2 + yComponent ** 2) ** (1/2)

    return [[xComponent, yComponent], radius]

def handleColocationNodeOnElement(element, sourcePointCoordinate):
    colocationNodesCoordinates = element.getColocationNodesCoordinates(duplicatedNodes, geometricNodes)
    
    # if colocationNodesCoordinates.count(sorcePointCoordinate) == 1:
    if sourcePointCoordinate in colocationNodesCoordinates:
        return True
    else: 
        return False    

def UContribuitionWithSingularitySubtraction(elementType, jacobian, sourcePoint, fieldPoint, radius, integrationPointIndex):
    UInitialContribuition = (-1 / (2 * math.pi)) * math.log(radius, math.e) * jacobian * weights[integrationPointIndex]
    
    radiusCheck = jacobian * (fieldPoint - sourcePoint)
    USecondTerm = (-1 / (2 * math.pi)) * math.log(abs(radiusCheck), math.e) * jacobian * weights[integrationPointIndex]

    cauchyPrincipalValue = 0

    if elementType == "semicontinuous" or elementType == "discontinuous":
        secondTermCPV = (1 + fieldPoint) * math.log(jacobian * (1 + fieldPoint), math.e) + (1 - fieldPoint) * math.log10(jacobian * (1 - fieldPoint)) - (1 + fieldPoint) - (1 - fieldPoint)

        cauchyPrincipalValue += (-1 / (2 * math.pi)) * secondTermCPV
    else: 
        secondTermCPV = (1 + fieldPoint) * math.log(jacobian * (1 + fieldPoint), math.e) + (1 - fieldPoint) * math.log10(jacobian * (1 - fieldPoint)) - (1 + fieldPoint) - (1 - fieldPoint)
        cauchyPrincipalValue += (-1 / (2 * math.pi)) * secondTermCPV

    U = UInitialContribuition - USecondTerm + cauchyPrincipalValue

    return U

def getIndex(newList, parameter):
    for i in range(len(newList)):
        if parameter == newList[i]:
            return i

def getHandGMatrices(sourcePoints: list, elementsList: list, geometricNodes):
    HMatrix = np.zeros((len(sourcePoints), len(sourcePoints)))
    GMatrix = np.zeros((len(sourcePoints), len(sourcePoints)))
    correctionMatrix = np.identity(len(sourcePoints)) * 1/2

    for sp in range(len(sourcePoints)):
        
        for el in range(len(elementsList)):
            integrationPointsRealCoordinates = getIntegrationPointsCoordinatesPerElement(elementsList[el], geometricNodes)
            elementNodes = elementsList[el].getElementNodesRealCoordinates(geometricNodes)
            elementType = elementsList[el].handleElementType(duplicatedNodes)
            colocationCoordinates = elementsList[el].getColocationNodesCoordinates(duplicatedNodes, geometricNodes)
            adimentionalPoints = elementsList[el].getAdimensionalPointsBasedOnGeometricCoordinates()
            colocationAdimentionalPoints = elementsList[el].getAdimensionalPointsBasedOnElementContinuity(duplicatedNodes)
            [_, jacobian, normalVector] = getPointsPropertiesOnElement(integrationPoints, elementNodes, adimentionalPoints)
            sourcePointIsOnElement = handleColocationNodeOnElement(elementsList[el], sourcePoints[sp])
            nodeIndex = getIndex(list(colocationCoordinates), sourcePoints[sp])

            if type(nodeIndex) is int:
                sourcePointAdimentionalCoordinate = colocationAdimentionalPoints[nodeIndex]
                                   
            for ip in range(len(integrationPoints)):
                integrationPointRadius = getRadius(sourcePoints[sp], integrationPointsRealCoordinates[ip])
                radiusComponent1 = integrationPointRadius[0][0] / integrationPointRadius[1]
                radiusComponent2 = integrationPointRadius[0][1] / integrationPointRadius[1]

                DRDN = radiusComponent1 * normalVector[ip][0] + radiusComponent2 * normalVector[ip][1]
                
                Q = (-1 / (2 * math.pi * integrationPointRadius[1])) * DRDN * jacobian[ip] * weights[ip]

                # checagem se o elemento contem o ponto fonte e se é continuo/semi/descontinuo e
                # avaliar a contribuição no fluxo.
                if sourcePointIsOnElement:
                    U = UContribuitionWithSingularitySubtraction(elementType, jacobian[ip], sourcePointAdimentionalCoordinate, integrationPoints[ip], integrationPointRadius[1], ip)
                                        
                else:
                    U = (-1 / (2 * math.pi)) * math.log(integrationPointRadius[1], math.e) * jacobian[ip] * weights[ip]

                for en in range(len(elementNodes)):
                    shapeFunctionValueOnIP = getShapeFunctionValueOnNode(integrationPoints[ip], en, adimentionalPoints)

                    Qcontribution = Q * shapeFunctionValueOnIP
                    Ucontribution = U * shapeFunctionValueOnIP

                    HMatrix[sp][elements[el][en]] += Qcontribution
                    GMatrix[sp][elements[el][en]] += Ucontribution

    HMatrix = np.array(HMatrix) + correctionMatrix
    GMatrix = np.array(GMatrix)
    
    return HMatrix, GMatrix

HMatrix, GMatrix = getHandGMatrices(sourcePoints, elementsList, geometricNodes)

def getFinalComponents(HMatrix, GMatrix, u, q, colocationMesh):
    FVector = np.zeros([len(colocationMesh)], dtype=float)

    for j in range(len(u)):
        HMatrix[:, u[j][0]] = - GMatrix[:, u[j][0]]
        GMatrix[:, u[j][0]] = - HMatrix[:, u[j][0]]

        FVector[u[j][0]] = u[j][1]

    for k in range(len(q)):
        FVector[q[k][0]] = q[k][1]

    return HMatrix, GMatrix, FVector

HMatrix, GMatrix, FVector = getFinalComponents(HMatrix, GMatrix, u, q, colocationMesh)

results = np.linalg.solve(HMatrix, np.dot(GMatrix, FVector))

print(results)

def checkHMatrixRowsValues(HMatrix):
    sum = np.zeros(len(sourcePoints))
    for j in range(len(sum)):
        for i in range(len(HMatrix)):
            sum[j] += HMatrix[j][i]

    print(sum)
# checkHMatrixRowsValues(HMatrix)