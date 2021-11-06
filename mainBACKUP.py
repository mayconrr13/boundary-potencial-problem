import math
import numpy as np
from numpy.core.arrayprint import printoptions
from shapeFunctions import getShapeFunctionValueOnNode
from tangentNormalAndJacobian import getPointsPropertiesOnElement
# data input
## nodes, elements - GEOMETRIC MESH

geometricNodes = [
    [2, 0],
    [5, 0],
    [8, 0],
    [8, 0],
    [8, 4],
    [8, 4],
    [2, 0]
]

elements = [
    [0, 1],
    [1, 2],
    [3, 4],
    [5, 6]
]
4

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
        colocationNodes = np.zeros(len(self.nodeList), dtype=list)

        for i in range(len(continuityBasedAdimensionalPoints)):
            xCoordinate = 0
            yCoordinate = 0

            for j in range(len(self.nodeList)):
                xCoordinate += getShapeFunctionValueOnNode(continuityBasedAdimensionalPoints[i], j, adimensionalPoints) * geometricNodes[self.nodeList[j]][0]
                yCoordinate += getShapeFunctionValueOnNode(continuityBasedAdimensionalPoints[i], j, adimensionalPoints) * geometricNodes[self.nodeList[j]][1]

            colocationNodes[i] = [xCoordinate, yCoordinate]

        return colocationNodes

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

elementsList = np.zeros(len(elements), dtype=list)
for el in range(len(elements)):
    elementsList[el] = Element(elements[el])

# process
## determinar malha de colocação
# checagem de nós duplicados
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

duplicatedNodes = getDuplicatedNodes(geometricNodes)

# criação da malha
# def generateColocationNodes(elementsList: list, duplicatedNodes):
#     colocationNodes = []

#     for i in range(len(elementsList)):
#         element = elementsList[i]
        
#         for j in range(len(element)):
#             if element[j] in duplicatedNodes:
#                 discontinuousNode = geometricNodes[element[j]]

#                 if j == 0:
#                     nextNode = geometricNodes[element[j + 1]]

#                     xColocationCoordinate = discontinuousNode[0] + (nextNode[0] - discontinuousNode[0]) * 0.25
#                     yColocationCoordinate = discontinuousNode[1] + (nextNode[1] - discontinuousNode[1]) * 0.25
                    
#                     colocationNodes.append([xColocationCoordinate, yColocationCoordinate])
                
#                 else:
#                     previousNode = geometricNodes[element[j - 1]]

#                     xColocationCoordinate = discontinuousNode[0] + (previousNode[0] - discontinuousNode[0]) * 0.25
#                     yColocationCoordinate = discontinuousNode[1] + (previousNode[1] - discontinuousNode[1]) * 0.25
                    
#                     colocationNodes.append([xColocationCoordinate, yColocationCoordinate])

#             else:
#                 continuousNode = geometricNodes[element[j]]

#                 if continuousNode not in colocationNodes:
#                     colocationNodes.append(continuousNode)
                         
#     return colocationNodes

# colocationNodes = generateColocationNodes(elements, duplicatedNodes)

# def getElementColocationNodes(element: list, colocationNodes: list): 
#     elementColocationNodes = []

#     for j in range(len(element)):
#         elementColocationNodes.append(colocationNodes[element[j]])

#     return elementColocationNodes

# ## criação dos pontos fontes
# def getSourcePoints(colocationNodes: list, elementsList: list):
#     sourcePointsList = []

#     for i in range(len(elementsList)):
#         elementColocationNodes = getElementColocationNodes(elementsList[i], colocationNodes)  
#         nodesProperties = getPointsPropertiesOnElement(elementsList[i], elementColocationNodes, [-1, 1])
        
#         for k in range(len(elementsList[i])):
#             xCoordinate = nodesProperties[2][k][0] * 1 + elementColocationNodes[k][0]
#             yCoordinate = nodesProperties[2][k][1] * 1 + elementColocationNodes[k][1]

#             sourcePointCoordinates = [xCoordinate, yCoordinate]

#             if sourcePointsList.count(sourcePointCoordinates) < 1:
#                 sourcePointsList.append(sourcePointCoordinates)

#     return sourcePointsList

# sourcePoints = getSourcePoints(colocationNodes, elements)

# integrationPoints = [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116]
# weights = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]

# ## para cara ponto fonte - para cada elemento - para cada ponto de integração

#     ## coordenada real do ponto de integração
# def getIntegrationPointCoordinates(points: float, elementNodes: list, adimentionalPoints: list):
#     elementPointsCoordinates = []
#     xComponent = 0
#     yComponent = 0

#     for j in range(len(elementNodes)):
#         xComponent += getShapeFunctionValueOnNode(points, j, adimentionalPoints) * elementNodes[j][0]
#         yComponent += getShapeFunctionValueOnNode(points, j, adimentionalPoints) * elementNodes[j][1]
    
#     elementPointsCoordinates = [xComponent, yComponent]

#     return elementPointsCoordinates
  
# def getElementNodes(element: list, geometricNodes: list):
#     elementNodesCoordinates = []

#     for i in range(len(element)):
#         elementNodesCoordinates.append(geometricNodes[element[i]])

#     return elementNodesCoordinates

# def getIntegrationPointsCoordinatesPerElement(element: list):
#     nodesCoordinatesPerElement = []
#     elementNodes = getElementNodes(element, geometricNodes)

#     for ip in range(len(integrationPoints)):
#         result = getIntegrationPointCoordinates(integrationPoints[ip], elementNodes, [-1, 1])
#         nodesCoordinatesPerElement.append(result)
    
#     return nodesCoordinatesPerElement

#     ## determinar raio

# def getRadius(sourcePoint: list, integrationPointCoordinates: list):
#     xComponent = integrationPointCoordinates[0] - sourcePoint[0]
#     yComponent = integrationPointCoordinates[1] - sourcePoint[1]
#     radius = (xComponent ** 2 + yComponent ** 2) ** (1/2)

#     return [[xComponent, yComponent], radius]

#     ## preencher matrizes H e G


# HMatrix = np.zeros((len(sourcePoints), len(sourcePoints)))
# GMatrix = np.zeros((len(sourcePoints), len(sourcePoints)))

# for sp in range(len(sourcePoints)):
#     for el in range(len(elements)):
#         integrationPointsRealCoordinates = getIntegrationPointsCoordinatesPerElement(elements[el])
#         elementNodes = getElementNodes(elements[el], geometricNodes)
#         [_, jacobian, normalVector] = getPointsPropertiesOnElement(integrationPoints, elementNodes, [-1, 1])
        
#         for ip in range(len(integrationPoints)):
#             integrationPointRadius = getRadius(sourcePoints[sp], integrationPointsRealCoordinates[ip])

#             radiusComponent1 = integrationPointRadius[0][0] / integrationPointRadius[1]
#             radiusComponent2 = integrationPointRadius[0][1] / integrationPointRadius[1]

#             DRDN = radiusComponent1 * normalVector[ip][0] + radiusComponent2 * normalVector[ip][1]
            
#             Q = (-1 / (2 * math.pi * integrationPointRadius[1])) * DRDN * jacobian[ip] * weights[ip]
#             U = (-1 / (2 * math.pi)) * math.log10(integrationPointRadius[1]) * jacobian[ip] * weights[ip]

#             for en in range(len(elementNodes)):
#                 shapeFunctionValueOnIP = getShapeFunctionValueOnNode(integrationPoints[ip], en, [-1, 1])

#                 Qcontribution = Q * shapeFunctionValueOnIP
#                 Ucontribution = U * shapeFunctionValueOnIP

#                 HMatrix[sp][elements[el][en]] += Qcontribution
#                 GMatrix[sp][elements[el][en]] += Ucontribution
    

# # results output
# ## retornar incognitas