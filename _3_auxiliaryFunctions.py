import math
import numpy as np

from _2_elementClass import Element
from integrationsPointsAndWeights import integrationPoints, weights
from shapeFunctions import getShapeFunctionValueOnNode
from tangentNormalAndJacobian import getPointProperty, getPointsPropertiesOnElement
from _1_readInputFile import readInputFile

def getElementsList(elements: list): 
    elementsList = np.zeros(len(elements), dtype=list)
    
    for el in range(len(elements)):
        elementsList[el] = Element(elements[el])

    return elementsList

def getDuplicatedNodes(nodeList: list):
    duplicatedNodes = np.array([], dtype=int)

    for i in range(len(nodeList)):
        nodeToBeChecked = nodeList[i]
        nodeList.count(nodeToBeChecked)

        if nodeList.count(nodeToBeChecked) > 1:
            duplicatedNodes = np.append(duplicatedNodes, i)

    return duplicatedNodes

def generateColocationMesh(elementsList, duplicatedNodes, geometricNodes):
    colocationMesh = []

    for i in range(len(elementsList)):
        elementColocationCoordinates = elementsList[i].getColocationNodesCoordinates(duplicatedNodes, geometricNodes)

        for j in range(len(elementColocationCoordinates)):
            if elementColocationCoordinates[j] not in colocationMesh:
                colocationMesh.append(elementColocationCoordinates[j])

    return colocationMesh

def getSourcePoints(duplicatedNodes, geometricNodes, elementsList):
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

    return sourcePointsList

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


def handleColocationNodeOnElement(element, sourcePointCoordinate, duplicatedNodes, geometricNodes):
    colocationNodesCoordinates = element.getColocationNodesCoordinates(duplicatedNodes, geometricNodes)

    if sourcePointCoordinate in colocationNodesCoordinates:
        return True
    else:
        return False


def getCauchyPricipalValue(jacobian, sourcePoint):
    cauchyPrincipalValue = 0

    if abs(sourcePoint - 1) <= (10 ** -12) or abs(sourcePoint + 1) <= (10 ** -12):
        secondTermCPV = 2 * math.log((jacobian * 2), math.e) - 2

        cauchyPrincipalValue = (-1 / (2 * math.pi)) * jacobian * secondTermCPV

    else:
        p1 = (1 + sourcePoint) * \
            math.log((jacobian * (1 + sourcePoint)), math.e)
        p2 = (1 - sourcePoint) * \
            math.log((jacobian * (1 - sourcePoint)), math.e)
        p3 = - (1 + sourcePoint) - (1 - sourcePoint)

        cauchyPrincipalValue = (-1 / (2 * math.pi)) * jacobian * (p1+p2+p3)

    return cauchyPrincipalValue


def getIndex(newList, parameter):
    for i in range(len(newList)):
        if parameter == newList[i]:
            return i

def getHandGMatrices(sourcePoints: list, colocationMesh: list, elementsList: list, geometricNodes, duplicatedNodes, elements):
        HMatrix = np.zeros((len(sourcePoints), len(colocationMesh)))
        GMatrix = np.zeros((len(sourcePoints), len(colocationMesh)))

        for sp in range(len(sourcePoints)):
            for el in range(len(elementsList)):
                integrationPointsRealCoordinates = getIntegrationPointsCoordinatesPerElement(
                    elementsList[el], geometricNodes)
                elementNodes = elementsList[el].getElementNodesRealCoordinates(
                    geometricNodes)
                elementType = elementsList[el].handleElementType(duplicatedNodes)
                colocationCoordinates = elementsList[el].getColocationNodesCoordinates(
                    duplicatedNodes, geometricNodes)
                adimentionalPoints = elementsList[el].getAdimensionalPointsBasedOnGeometricCoordinates(
                )
                colocationAdimentionalPoints = elementsList[el].getAdimensionalPointsBasedOnElementContinuity(
                    duplicatedNodes)
                [_, jacobian, normalVector] = getPointsPropertiesOnElement(
                    integrationPoints, elementNodes, adimentionalPoints)
                sourcePointIsOnElement = handleColocationNodeOnElement(
                    elementsList[el], sourcePoints[sp], duplicatedNodes, geometricNodes)
                nodeIndex = getIndex(list(colocationCoordinates), sourcePoints[sp])

                for en in range(len(elementNodes)):
                    Q = 0
                    U = 0
                    cauchyPrincipalValue = 0

                    if sourcePointIsOnElement:
                        sourcePointAdimentionalCoordinate = colocationAdimentionalPoints[nodeIndex]
                        sourcePointJacobian = getPointProperty(
                            sourcePointAdimentionalCoordinate, elementNodes, adimentionalPoints)

                        shapeFunctionValueOnSource = getShapeFunctionValueOnNode(
                            sourcePointAdimentionalCoordinate, en, colocationAdimentionalPoints)
                        cauchyPrincipalValue = getCauchyPricipalValue(
                            sourcePointJacobian, sourcePointAdimentionalCoordinate)

                        GMatrix[sp][elements[el][en]] += cauchyPrincipalValue * \
                            shapeFunctionValueOnSource

                    for ip in range(len(integrationPoints)):
                        integrationPointRadius = getRadius(
                            sourcePoints[sp], integrationPointsRealCoordinates[ip])
                        radiusComponent1 = integrationPointRadius[0][0] / \
                            integrationPointRadius[1]
                        radiusComponent2 = integrationPointRadius[0][1] / \
                            integrationPointRadius[1]

                        DRDN = radiusComponent1 * \
                            normalVector[ip][0] + \
                            radiusComponent2 * normalVector[ip][1]
                        shapeFunctionValueOnIP = getShapeFunctionValueOnNode(
                            integrationPoints[ip], en, colocationAdimentionalPoints)

                        Q += (-1 / (2 * math.pi * integrationPointRadius[1])) * \
                            DRDN * jacobian[ip] * weights[ip] * \
                            shapeFunctionValueOnIP
                        U += (-1 / (2 * math.pi)) * math.log(
                            integrationPointRadius[1], math.e) * jacobian[ip] * weights[ip] * shapeFunctionValueOnIP

                        if sourcePointIsOnElement:
                            radiusCheck = jacobian[ip] * (
                                integrationPoints[ip] - sourcePointAdimentionalCoordinate)

                            U -= (-1 / (2 * math.pi)) * math.log(abs(radiusCheck), math.e) * \
                                sourcePointJacobian * \
                                weights[ip] * shapeFunctionValueOnSource

                    HMatrix[sp][elements[el][en]] += Q
                    GMatrix[sp][elements[el][en]] += U

        return HMatrix, GMatrix

def getFinalComponents(HMatrix, GMatrix, u, q, sourcePoints, colocationMesh):
        FHMatrix = np.zeros((len(sourcePoints), len(sourcePoints)), dtype=float)
        FGMatrix = np.zeros((len(sourcePoints), len(sourcePoints)), dtype=float)
        FVector = np.zeros(len(sourcePoints), dtype=float)

        if sourcePoints == colocationMesh:
            HMatrix += np.identity(len(sourcePoints)) * 1/2

        for j in range(len(u)):
            FHMatrix[:, u[j][0]] += - GMatrix[:, u[j][0]]
            FGMatrix[:, u[j][0]] += - HMatrix[:, u[j][0]]

            FVector[u[j][0]] += u[j][1]

        for k in range(len(q)):
            FHMatrix[:, q[k][0]] += HMatrix[:, q[k][0]]
            FGMatrix[:, q[k][0]] += GMatrix[:, q[k][0]]

            FVector[q[k][0]] += q[k][1]

        return FHMatrix, FGMatrix, FVector

def getPotentialAndFlowVector(results, u, q, sourcePoints):
        PotentialVector = np.zeros(len(sourcePoints), dtype=float)
        FlowVector = np.zeros(len(sourcePoints), dtype=float)

        for j in range(len(q)):
            FlowVector[q[j][0]] += q[j][1]
            PotentialVector[q[j][0]] += results[q[j][0]]

        for k in range(len(u)):
            FlowVector[u[k][0]] += results[u[k][0]]
            PotentialVector[u[k][0]] += u[k][1]

        return PotentialVector, FlowVector


def getInternalFluxMatrices(sourcePoints: list, colocationMesh: list, elementsList: list, geometricNodes, duplicatedNodes, elements):
    DMatrix = np.zeros((2 * len(sourcePoints), len(colocationMesh)))
    SMatrix = np.zeros((2 * len(sourcePoints), len(colocationMesh)))

    for sp in range(len(sourcePoints)):

        for el in range(len(elementsList)):
            integrationPointsRealCoordinates = getIntegrationPointsCoordinatesPerElement(elementsList[el], geometricNodes)
            elementNodes = elementsList[el].getElementNodesRealCoordinates(geometricNodes)
            adimentionalPoints = elementsList[el].getAdimensionalPointsBasedOnGeometricCoordinates()
            colocationAdimentionalPoints = elementsList[el].getAdimensionalPointsBasedOnElementContinuity(duplicatedNodes)
            [_, jacobian, normalVector] = getPointsPropertiesOnElement(integrationPoints, elementNodes, adimentionalPoints)

            for ip in range(len(integrationPoints)):
                integrationPointRadius = getRadius(sourcePoints[sp], integrationPointsRealCoordinates[ip])
                radiusComponent1 = integrationPointRadius[0][0] / integrationPointRadius[1]
                radiusComponent2 = integrationPointRadius[0][1] / integrationPointRadius[1]

                DRDN = radiusComponent1 * normalVector[ip][0] + radiusComponent2 * normalVector[ip][1]

                partialD = (1 / (2 * math.pi * integrationPointRadius[1])) * jacobian[ip] * weights[ip]
                D = partialD * np.array(integrationPointRadius[0]) / integrationPointRadius[1]

                partialS = (1 / (2 * math.pi * (integrationPointRadius[1]) ** 2)) * jacobian[ip] * weights[ip]
                S = partialS * (np.array(normalVector[ip] - 2 * (np.array(integrationPointRadius[0]) / integrationPointRadius[1]) * DRDN))

                for en in range(len(elementNodes)):
                    shapeFunctionValueOnIP = getShapeFunctionValueOnNode(integrationPoints[ip], en, colocationAdimentionalPoints)

                    Dcontribution = D * shapeFunctionValueOnIP
                    Scontribution = S * shapeFunctionValueOnIP

                    for i in range(2):

                        DMatrix[sp * 2 + i][elements[el][en]] += Dcontribution[i]
                        SMatrix[sp * 2 + i][elements[el][en]] += Scontribution[i]

    DMatrix = np.array(DMatrix)
    SMatrix = np.array(SMatrix)

    return DMatrix, SMatrix
