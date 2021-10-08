import math
from matplotlib import style    
import matplotlib.pyplot as plt 
import numpy as np
import time
import copy
##############################################################
#graph sets of xy coordinates
def graph_coords(x2, y2, min_dist, generation):
    #Define graph style
    style.use('dark_background')
    plt.clf()
    # plotting the points
    plt.plot(x2, y2,'yo-', label="Optimum Path")
    for i in range(len(x2) - 1):
        plt.annotate(i + 1, (x2[i], y2[i]), textcoords="offset points", xytext=(0,5), ha = 'center')
        
    # naming the axes 
    plt.xlabel('x - axis') 
    plt.ylabel('y - axis') 
    plt.legend()
    # giving a title to my graph 
    plt.title(("Optimum Distance : " + str(round(min_dist, 2)) + " Gen: " + str(generation)))
      
    # function to show the plot 
    plt.pause(.05)
    plt.show()    

    return

########################################
class node():
    def __init__(self, name, x, y):
        self.name = name
        self.x = x
        self.y = y
########################################
class chromosome():
    def __init__(self, name, nodelist):
        self.name = name
        self.genes = nodelist
        self.distance = math.inf
        self.currentAge = 0
        self.maxAge = ageMax
        
    def listNames(self):
        names = []
        for node in self.genes:
            names.append(node.name)
        return names
            
########################################
#read and parse data file
def read_datafile(path):
    #Import data file
    i = 0
    x = []
    y = []
    with open((path), "r") as file:
      for line in file:
        split_line = line.strip().split(" ")
        
        #Track line number to remove header info
        if i > 6:
            #Populate x,y coordinate pairs into arrays
            x.append(float(split_line[1]))
            y.append(float(split_line[2]))
        #increment line counter
        i += 1    
    return x, y
#########################################
#Calculate distance for the trip		
def calculate_trip_dist(trip, nodes):		
    dist = 0		
    for j in range(len(trip)-1):		
        node1 = nodes[trip[j]   -1]		
        node2 = nodes[trip[j+1] -1]		
        		
        dist = dist + (math.hypot(node1.x - node2.x, node1.y - node2.y))			
    		
    return dist  
##############################################################
def getXY(trip):
    x = []
    y = []
    for gene in trip:
        x.append(nodes[gene -1].x)
        y.append(nodes[gene -1].y)
                      
    return x, y
##############################################################
def generateSplitPoint(nodes):		
    a = np.random.randint(1,len(nodes)+1)		
    b = a		
    while(a == b and abs(a-b) < 1):		
        b = np.random.randint(1,len(nodes)+ 1)   		
    		
    splitPoint1 = min(a, b)		
    splitPoint2 = max(a, b) 		
    return splitPoint1, splitPoint2	
##############################################################		
def Crossover(parent1genes, parent2genes, splitPoint1, splitPoint2, nodes):		
    tempChromosome = []		
    #First section		
        		
    for x in range(splitPoint1):		
        tempChromosome.append(nodes[parent1genes[x] - 1])		
    #Last Section    		
    for x in range(splitPoint2, len(parent1genes)):		
        if(x <= len(nodes)):		
            tempChromosome.append(nodes[parent1genes[x] - 1])		
        else:		
            tempChromosome.append(nodes[parent1genes[0] - 1])		
    #Mid Section		
    for gene in parent2genes:		
        found = False		
        for gene2 in tempChromosome:		
            if(gene2.name == nodes[gene -1].name):		
                found = True		
                break		
        if(not found):		
            tempChromosome.insert(splitPoint1,(nodes[gene -1]))		
    if(len(tempChromosome) <101):
        tempChromosome.append(tempChromosome[0])
    
    return tempChromosome
##############################################################
def deviation(populationArr):
    best = populationArr[0].distance
    devArr = []
    for chromo in populationArr:
        devArr.append(chromo.distance - best)
    
    dev = np.sum(devArr)/ len(devArr)
    
    return dev   

'''#############################################################
def heatmap(matrix):
    plt.clf()
    plt.imshow(matrix, cmap='gray', interpolation='None')
    plt.colorbar()
    plt.clim(800, 4000)
    plt.pause(.05)
    plt.show()     
  
##############################################################'''
    
x     = []
y     = []
nodes = []
#######
#INPUT
######
    
#data file path
file_path = str(r'C:\Users\burkh\OneDrive\Desktop\CECS\AI\Project4\Random100.tsp')
#used to read and parse the tsp file
x, y = read_datafile(file_path)

#create array of nodes
for i in range(len(x)):
    n = node((i + 1), x[i], y[i])
    nodes.append(n)


#Number of times to run experiment
for q in range(65):
        
    ####################
    #BEGIN GA PROCESSES
    ###################
    #Population Size
    matrixSize = 10  #NxN
    gaPopSize = matrixSize * matrixSize#round(matrixSize * matrixSize /3)

    ageMax = 25 #99999
    i = 0
    deathCount = 0
    tempArr = []
    
    populationTorus = np.empty((matrixSize,matrixSize), dtype=chromosome) #Originally populated with "None"
    prevdist = math.inf
    '''graphTorus = np.zeros((matrixSize,matrixSize))'''
    #Create Chromosomes and add to population
    
    ###########################
    #CREATE INITIAL POPULATION
    ##########################
    i = 0
    y = 0
    for row in populationTorus:
        x = 0
        for cell in row:
            arr = np.random.permutation(nodes)	
            populationTorus[x][y] = copy.deepcopy(chromosome(i, np.append(arr, arr[0])))
            i +=1
            x +=1
        y +=1
    
    
    """for i in range(0,gaPopSize-1):
        #Randomly order the nodes to create a chromosome
        arr = np.random.permutation(nodes)		
        inserted = False
        
        #Randomly place in torus
        while(inserted == False):
            x = np.random.randint(0,len(populationTorus) - 1)
            y = np.random.randint(0,len(populationTorus) - 1)
            if(populationTorus[x][y] == None):
                populationTorus[x][y] = chromosome(i, np.append(arr, arr[0]))
                inserted = True"""
        
    mutationChance = 10000 #10000
    numGenertations = 4500
    generationNumber = 0
    noChange = 0
    #numNoChangeGen = 2000
    worstDist = math.inf
    chromoCount = 0
    chromoCount += gaPopSize
    
    start = time.time()
    
    #print("End Time" + "\t" + "Best Dist" + "\t" + "gen number" + "\t" + "When Found")
    
    #Loop through generations
    stopping = math.inf
    bestDist = math.inf
    #################################
    # Initial TEST FITNESS OF EACH CHROMOSOME
    ################################
    #Measure trip distance
    i = 0
    for row in populationTorus:
        j = 0
        for cell in row:
            if(cell != None):
                trip = cell.listNames()
                distance = calculate_trip_dist(trip, nodes)
                cell.distance = copy.deepcopy(distance)
                best = copy.deepcopy(cell)
                

                ########################
                
                #STORE BEST FITNESS 
                
                ########################
                
                if(distance < bestDist):
                    bestDist = distance
                    change = True
                    noChange = 0
                    whenFound = time.time() - start
                    genFound = copy.deepcopy(generationNumber)
                    #cell.currentAge - 1
                    #print(bestDist)    
        j +=1
    i +=1   
    
    
    while((generationNumber < numGenertations)):
        #print("gen: " + str(generationNumber))
        
        '''heatmap(graphTorus)'''
        change = False
        #loop through matrix
        i = 0
        for row in populationTorus:
            j = 0
            for cell in row:
                '''if(cell != None):
                    graphTorus[i][j] = populationTorus[i][j].distance
                else:
                    graphTorus[i][j] = 0'''
                if(cell != None):
                    #age the chromosome by 1 iteration
                    cell.currentAge += 1
                    
                    if((cell.maxAge - cell.currentAge) <= 0):
                        #Cell dies
                        deathCount +=1
                        populationTorus[i][j] = None
                        j +=1
                        continue
                    #Cell searches to reproduce locally
                    #i = column
                    #j = row
                    #(i-1,j-1)   (i,j-1)   (i+1,j-1)
                    #(i-1,j)     (i,j)     (i+1,j)
                    #(i-1,j+1)   (i,j+1)   (i+1,j+1)
                    
                    #select randomly which of the 8 sorrounding spaces to reproduce with
                    x = np.random.randint(-1,2) + j
                    y = np.random.randint(-1,2) + i
                    
                    #If indice goes off matrix - flip to other edge to form torus
                    if(x > len(populationTorus) - 1):
                        x = 0
                    if(x < 0):
                        x = len(populationTorus) - 1
                    if(y > len(populationTorus) - 1):
                        y = 0 
                    if(y < 0):
                        y = len(populationTorus) - 1
                        
                    #See if cell is None - asexual reproduction (with mutation) if so, crossover if not
                    if(populationTorus[x][y] == None):
                        
                        ##########
                        #ASEXUAL
                        #########
                        
                        #cell is empty - asexually fill with mutated version of self
                        populationTorus[x][y] = copy.deepcopy(populationTorus[i][j])
                        mutationIndex1 = np.random.randint(0,len(nodes)-1)
                        mutationIndex2 = np.random.randint(0,len(nodes)-1)
                        
                        if(mutationIndex1 != mutationIndex2):
                            #swap genes
                            tempChromo = copy.deepcopy(populationTorus[x][y].genes[mutationIndex1])
                            populationTorus[x][y].genes[mutationIndex1] = copy.deepcopy(populationTorus[x][y].genes[mutationIndex2])
                            populationTorus[x][y].genes[mutationIndex2] = copy.deepcopy(tempChromo)
                            
                            #If index was start of trip. Make end of trip the same as start
                            if(mutationIndex1 == 0 or mutationIndex2 == 0):
                                populationTorus[x][y].genes[len(nodes)] = populationTorus[x][y].genes[0]

                        #re-calculate fitness
                        trip = populationTorus[x][y].listNames()
                        distance = calculate_trip_dist(trip, nodes)
                        populationTorus[x][y].distance = copy.deepcopy(distance)                            
                        populationTorus[x][y].currentAge = 0
                        populationTorus[x][y].maxAge = ageMax   
                        #Increment pop size count
                        chromoCount +=1
                        

                    else:
                        #cell has a chromosome in it - Crossover to replace the cell with lower fitness
                        
                        ##########
                        #SEXUAL
                        #########
                            #######################
                            #CROSSOVER - TWO POINT
                            ######################

                        #index of chromosome to be split at for crossover
                        splitPoint1, splitPoint2 = generateSplitPoint(nodes)
                        
                        #Get parent genes
                        parent1genes = populationTorus[i][j].listNames()
                        parent2genes = populationTorus[x][y].listNames()
                        
                        #pull first n nodes from parent1 til min split point
                        tempChromosome1 = Crossover(parent1genes, parent2genes, splitPoint1, splitPoint2, nodes)

                        #Replace less fit parent with offspring
                            #pick parent to replace
                            
                        if(populationTorus[i][j].distance < populationTorus[x][y].distance):
                            populationTorus[x][y].genes = tempChromosome1
                            #re-calculate fitness
                            trip = populationTorus[x][y].listNames()
                            distance = calculate_trip_dist(trip, nodes)
                            populationTorus[x][y].distance = copy.deepcopy(distance)                            
                            populationTorus[x][y].currentAge = 0
                            populationTorus[x][y].maxAge = ageMax
                          
                            
                            #############
                            #MUTATION
                            #############
                            
                            mutationIndicator = np.random.randint(1,mutationChance)
                            
                            if(mutationIndicator == 1):
                                #Mutation occurred - select gene index randomly to swap with
                                mutationIndex1 = np.random.randint(0,len(nodes)-1)
                                mutationIndex2 = np.random.randint(0,len(nodes)-1)
                                
                                if(mutationIndex1 != mutationIndex2):
                                    #swap genes
                                    tempChromo = copy.deepcopy(populationTorus[x][y].genes[mutationIndex1])
                                    populationTorus[x][y].genes[mutationIndex1] = copy.deepcopy(populationTorus[x][y].genes[mutationIndex2])
                                    populationTorus[x][y].genes[mutationIndex2] = copy.deepcopy(tempChromo)
                                    
                                    #If index was start of trip. Make end of trip the same as start
                                    if(mutationIndex1 == 0 or mutationIndex2 == 0):
                                        populationTorus[x][y].genes[len(nodes)] = populationTorus[x][y].genes[0]
        
                                #re-calculate fitness
                                trip = populationTorus[x][y].listNames()
                                distance = calculate_trip_dist(trip, nodes)
                                populationTorus[x][y].distance = copy.deepcopy(distance)                            
                                populationTorus[x][y].currentAge = 0
                                populationTorus[x][y].maxAge = ageMax  
                                #increment pop size counter
                                chromoCount +=1
                                if(distance < bestDist):
                                    bestDist = distance
                                    currTime = time.time() - start
                                    whenFound = time.time() - start
                                    genFound = copy.deepcopy(generationNumber)                                    
                                    #print(bestDist)                                         
                            
                            
                            
                        else:
                            populationTorus[i][j].genes = tempChromosome1
                            #re-calculate fitness
                            trip = populationTorus[i][j].listNames()
                            distance = calculate_trip_dist(trip, nodes)
                            populationTorus[i][j].distance = copy.deepcopy(distance)                                
                            populationTorus[i][j].currentAge = 0
                            populationTorus[i][j].maxAge = ageMax  
                            #increment pop size counter
                            chromoCount +=1                            
                            if(distance < bestDist):
                                bestDist = distance
                                whenFound = time.time() - start
                                genFound = copy.deepcopy(generationNumber)                                
                                #print(bestDist)                        
                            #############
                            #MUTATION
                            #############
                            
                            mutationIndicator = np.random.randint(1,mutationChance)
                            
                            if(mutationIndicator == 1):
                                #Mutation occurred - select gene index randomly to swap with
                                mutationIndex1 = np.random.randint(0,len(nodes)-1)
                                mutationIndex2 = np.random.randint(0,len(nodes)-1)
                                
                                if(mutationIndex1 != mutationIndex2):
                                    #swap genes
                                    tempChromo = copy.deepcopy(populationTorus[i][j].genes[mutationIndex1])
                                    populationTorus[i][j].genes[mutationIndex1] = copy.deepcopy(populationTorus[i][j].genes[mutationIndex2])
                                    populationTorus[i][j].genes[mutationIndex2] = copy.deepcopy(tempChromo)
                                    
                                    #If index was start of trip. Make end of trip the same as start
                                    if(mutationIndex1 == 0 or mutationIndex2 == 0):
                                        populationTorus[i][j].genes[len(nodes)] = populationTorus[i][j].genes[0]
        
                                #re-calculate fitness
                                trip = populationTorus[i][j].listNames()
                                distance = calculate_trip_dist(trip, nodes)
                                populationTorus[i][j].distance = copy.deepcopy(distance)                            
                                populationTorus[i][j].currentAge = 0
                                populationTorus[i][j].maxAge = ageMax                           
                                if(distance < bestDist):
                                    bestDist = distance
                                    whenFound = time.time() - start
                                    genFound = copy.deepcopy(generationNumber)                                    
                                    #print(bestDist)                    
                j +=1
            i +=1        
        if(change == False):
            noChange +=1
        generationNumber += 1 #Iterate generation number an repeat process until stopping criteria
    #Time at which best solution was found    
    #print(str(round(time.time() - start,2)) + "\t" + str(round(bestDist,2)) + "\t" + str(generationNumber) + "\t" + str(whenFound))     
    #print("deaths: " + str(deathCount) + " Total Chromos: ")
    print(str(round(time.time() - start ,2)) + "\t" + str(round(bestDist,2)) + "\t" + str(round(genFound,2)) + "\t" + str(round(whenFound,2)) + "\t" + str(deathCount) + "\t" + str(chromoCount))  