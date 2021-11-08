import numpy as np

from _1_readInputFile import readInputFile
from _3_auxiliaryFunctions import generateAuxiliaryMesh, getDuplicatedNodes, getElementsList, getFinalComponents, getHandGMatrices, getInternalFluxMatrices, getPotentialAndFlowVector
from _4_generateResultsFile import generateResultsFile

def potentialBEMProcess():
    # Leitura do arquivo de entrada
        # Artigo - 3 elementos lineares por lado
    u, q, geometricNodes, internalPoints, elements = readInputFile("ex1_imputFile.txt")

        # Artigo - 1 elemento cúbico por lado
    # u, q, geometricNodes, internalPoints, elements = readInputFile("ex2_imputFile.txt")

    # Criação da lista de lementos com base na classe Elementos
    elementsList = getElementsList(elements)

    # Geração da malha com pontos fontes no contorno
    duplicatedNodes = getDuplicatedNodes(geometricNodes)
    auxiliaryMesh = generateAuxiliaryMesh(elementsList, duplicatedNodes, geometricNodes)
    sourcePoints = auxiliaryMesh

    # Obtenção das matrizes globais H e G
    HMatrix, GMatrix = getHandGMatrices(sourcePoints, auxiliaryMesh, elementsList, geometricNodes, duplicatedNodes, elements)
   
    # Troca de colunas com base nas grandezas prescritas
    FHMatrix, FGMatrix, FVector = getFinalComponents(HMatrix, GMatrix, u, q, sourcePoints, auxiliaryMesh)
    
    # Resolução do problema e geração dos vetores de potencial e fluxo no contorno
    resultsVector = np.linalg.solve(FHMatrix, np.dot(FGMatrix, FVector))
    PotentialVector, FlowVector = getPotentialAndFlowVector(resultsVector, u, q, sourcePoints)

    # Determinação do potencial nos pontos internos
    IntHMatrix, IntGMatrix = getHandGMatrices(internalPoints, auxiliaryMesh, elementsList, geometricNodes, duplicatedNodes, elements)
    IntPotentialResults = np.dot(IntGMatrix, FlowVector) - np.dot(IntHMatrix, PotentialVector)

    # Determinação do fluxo nos pontos internos
    DMatrix, SMatrix = getInternalFluxMatrices(internalPoints, auxiliaryMesh, elementsList, geometricNodes, duplicatedNodes, elements)
    IntFlowResults = np.dot(DMatrix, FlowVector) - np.dot(SMatrix, PotentialVector)

    # Geração do arquivos de resultados
    generateResultsFile(PotentialVector, FlowVector, internalPoints, IntPotentialResults, IntFlowResults)

potentialBEMProcess()