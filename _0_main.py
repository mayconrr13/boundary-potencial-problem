import numpy as np

from _1_readInputFile import u, q, geometricNodes, internalPoints, elements
from _3_auxiliaryFunctions import generateColocationMesh, getDuplicatedNodes, getElementsList, getFinalComponents, getHandGMatrices, getInternalFluxMatrices, getPotentialAndFlowVector
from _4_generateResultsFile import generateResultsFile

def potentialBEMProcess():
    # Leitura do arquivo de entrada
    # desenvolver função

    # Criação da lista de lementos com base na classe Elementos
    elementsList = getElementsList(elements)

    # Geração da malha com pontos fontes no contorno
    duplicatedNodes = getDuplicatedNodes(geometricNodes)
    colocationMesh = generateColocationMesh(elementsList, duplicatedNodes, geometricNodes)
    sourcePoints = colocationMesh

    # Obtenção das matrizes globais H e G
    HMatrix, GMatrix = getHandGMatrices(sourcePoints, colocationMesh, elementsList, geometricNodes, duplicatedNodes)
   
    # Troca de colunas com base nas grandezas prescritas
    FHMatrix, FGMatrix, FVector = getFinalComponents(HMatrix, GMatrix, u, q, sourcePoints, colocationMesh)
    
    # Resolução do problema e geração dos vetores de potencial e fluxo no contorno
    resultsVector = np.linalg.solve(FHMatrix, np.dot(FGMatrix, FVector))
    PotentialVector, FlowVector = getPotentialAndFlowVector(resultsVector, u, q, sourcePoints)

    # Determinação do potencial nos pontos internos
    IntHMatrix, IntGMatrix = getHandGMatrices(internalPoints, colocationMesh, elementsList, geometricNodes, duplicatedNodes)
    IntPotentialResults = np.dot(IntGMatrix, FlowVector) - np.dot(IntHMatrix, PotentialVector)

    # Determinação do fluxo nos pontos internos
    DMatrix, SMatrix = getInternalFluxMatrices(internalPoints, colocationMesh, elementsList, geometricNodes, duplicatedNodes)
    IntFlowResults = np.dot(DMatrix, FlowVector) - np.dot(SMatrix, PotentialVector)

    # Organização e geração do arquivos de resultados
    boundaryResults = [PotentialVector, FlowVector]
    internalResults = []
    for i in range(len(internalPoints)):
        internalResults.append([IntPotentialResults[i], []])

        for j in range(2):
            internalResults[i][1].append(IntFlowResults[2 * i + j])

    generateResultsFile(boundaryResults, internalResults)


potentialBEMProcess()