from shapeFunctions import getDiffShapeFunction

# property POINTS could be an integration points list or an element nodes list
def getPointTangentVector(points: list, elementNodes: list, adimentionalPoints: list):
    tangentVectors = []

    for i in range(len(points)):
        xTangentComponent = 0
        yTangentComponent = 0

        for j in range(len(elementNodes)):
            xTangentComponent += getDiffShapeFunction(points[i], j, adimentionalPoints) * elementNodes[j][0]
            yTangentComponent += getDiffShapeFunction(points[i], j, adimentionalPoints) * elementNodes[j][1]
        
        tangentVectors.append([xTangentComponent, yTangentComponent])

    return tangentVectors

def getPointJacobian(tangentVectors: list):
    jacobians = []

    for i in range(len(tangentVectors)):
        jacobian = (tangentVectors[i][0] ** 2 + tangentVectors[i][1] ** 2) ** (1 / 2)
        jacobians.append(jacobian)

    return jacobians

def getPointNormalVetor(tangentVectors: list, jacobians: list):
    normalVectors = []

    for i in range(len(tangentVectors)):
        xNormalComponent = tangentVectors[i][1] / jacobians[i]
        yNormalComponent = - tangentVectors[i][0] / jacobians[i]

        normalVector = [xNormalComponent, yNormalComponent]

        normalVectors.append(normalVector)
    
    return normalVectors

def getPointsPropertiesOnElement(points: list, elementNodes: list, adimentionalPoints: list):
    tangentVectors = getPointTangentVector(points, elementNodes, adimentionalPoints)
    jacobians = getPointJacobian(tangentVectors)
    normalVectors = getPointNormalVetor(tangentVectors, jacobians)

    elementPointsProperties = [tangentVectors, jacobians, normalVectors]

    return elementPointsProperties

# results = getPointsPropertiesOnElement([7, 15], [[2, 0], [5, 0]], [-1, 1])
# print(results)

def getPointProperty(ksi, elementNodes, adimentionalPoints):
    xTangentComponent = 0
    yTangentComponent = 0

    for j in range(len(elementNodes)):
        xTangentComponent += getDiffShapeFunction(ksi, j, adimentionalPoints) * elementNodes[j][0]
        yTangentComponent += getDiffShapeFunction(ksi, j, adimentionalPoints) * elementNodes[j][1]
    
    jacobian = (xTangentComponent ** 2 + yTangentComponent ** 2) ** (1 / 2)

    xNormalComponent = yTangentComponent / jacobian
    yNormalComponent = - xTangentComponent / jacobian

    return jacobian