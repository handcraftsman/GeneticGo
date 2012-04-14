package genetic

import (
	"bytes"
	"strings"
)

func (solver *Solver) getMutationStrategyResultChannel() chan sequenceInfo {
	mutateStrategyResults := solver.strategies[0].results
	for _, strategy := range solver.strategies {
		if strings.Contains(strategy.name, "mutate") {
			mutateStrategyResults = strategy.results
			break
		}
	}
	return mutateStrategyResults
}

func (solver *Solver) crossover(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)
	mutateStrategyResults := solver.getMutationStrategyResultChannel()

	for {
		select {
		case <-solver.quit:
			solver.quit <- true
			return
		default:
			break
		}

		parentA := solver.pool[random.Intn(len(solver.pool))].genes
		parentB := solver.pool[0].genes

		if len(parentA) == 1 || len(parentB) == 1 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		sourceStart := random.Intn((len(parentB)-1)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		destinationStart := random.Intn((len(parentA)-1)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		maxLength := min(len(parentA)-destinationStart, len(parentB)-sourceStart)
		length := (1 + random.Intn(maxLength/numberOfGenesPerChromosome-1)) * numberOfGenesPerChromosome

		childGenes := bytes.NewBuffer(make([]byte, 0, max(len(parentA), len(parentB))))

		if destinationStart > 0 {
			childGenes.WriteString(parentA[0:destinationStart])
		}

		childGenes.WriteString(parentB[sourceStart : sourceStart+length])

		if childGenes.Len() < len(parentA) {
			childGenes.WriteString(parentA[childGenes.Len():len(parentA)])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy}
		child.fitness = getFitness(child.genes)

		select {
		case strategy.results <- child:
			continue
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) mutate(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)

	for {
		select {
		case <-solver.quit:
			solver.quit <- true
			return
		default:
			break
		}

		useBestParent := random.Intn(100) == 0
		parent := solver.pool[0].genes
		if !useBestParent {
			parent = solver.pool[random.Intn(len(solver.pool))].genes
		}

		parentIndex := random.Intn(len(parent))

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parent)))
		if parentIndex > 0 {
			childGenes.WriteString(parent[:parentIndex])
		}

		currentGene := parent[parentIndex : parentIndex+1]

		for {
			gene := <-solver.nextGene
			if len(gene) != 1 {
				return
			}
			if gene != currentGene {
				childGenes.WriteString(gene)
				break
			}
		}

		if parentIndex+1 < len(parent) {
			childGenes.WriteString(parent[parentIndex+1:])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy}
		child.fitness = getFitness(child.genes)

		select {
		case strategy.results <- child:
			continue
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) reverse(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)
	mutateStrategyResults := solver.getMutationStrategyResultChannel()

	for {
		select {
		case <-solver.quit:
			solver.quit <- true
			return
		default:
			break
		}

		useBestParent := random.Intn(100) == 0
		parent := solver.pool[0].genes
		if !useBestParent {
			parent = solver.pool[random.Intn(len(solver.pool))].genes
		}

		if len(parent) == 1 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		reversePointA := random.Intn(len(parent)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		reversePointB := random.Intn(len(parent)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		for ; reversePointA == reversePointB; reversePointB = random.Intn(len(parent)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome {
		}

		min, max := sort(reversePointA, reversePointB)

		fragments := make([]string, max-min)
		for i := min; i < max; i += numberOfGenesPerChromosome {
			fragments[i-min] = parent[i : i+numberOfGenesPerChromosome]
		}

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parent)))
		if min > 0 {
			childGenes.WriteString(parent[0:min])
		}

		reverseArray(fragments)
		for _, fragment := range fragments {
			childGenes.WriteString(fragment)
		}

		if childGenes.Len() < len(parent) {
			childGenes.WriteString(parent[childGenes.Len():len(parent)])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy}
		child.fitness = getFitness(child.genes)

		select {
		case strategy.results <- child:
			continue
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) shift(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)
	mutateStrategyResults := solver.getMutationStrategyResultChannel()

	for {
		select {
		case <-solver.quit:
			solver.quit <- true
			return
		default:
			break
		}

		useBestParent := random.Intn(100) == 0
		parent := solver.pool[0].genes
		if !useBestParent {
			parent = solver.pool[random.Intn(len(solver.pool))].genes
		}

		if len(parent) < numberOfGenesPerChromosome+1 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}
		shiftRight := random.Intn(2) == 1

		segmentStart := random.Intn((len(parent)-numberOfGenesPerChromosome)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		segmentCount := 1 + random.Intn((len(parent)-segmentStart+numberOfGenesPerChromosome)/numberOfGenesPerChromosome-1)

		// +2 because first and last will be empty to leave room
		fragments := make([]string, segmentCount+2)
		for i := 0; i < segmentCount; i++ {
			start := segmentStart + i*numberOfGenesPerChromosome
			fragments[i+1] = parent[start : start+numberOfGenesPerChromosome]
		}
		start := 1
		end := segmentCount

		if shiftRight {
			start = 0
			fragments[0] = fragments[end]
			end--
		} else {
			end++
			fragments[end] = fragments[0]
			start++
		}

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parent)))
		if segmentStart > 0 {
			childGenes.WriteString(parent[0:segmentStart])
		}

		for i := start; i <= end; i++ {
			childGenes.WriteString(fragments[i])
		}

		if childGenes.Len() < len(parent) {
			childGenes.WriteString(parent[childGenes.Len():len(parent)])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy}
		child.fitness = getFitness(child.genes)

		select {
		case strategy.results <- child:
			continue
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) swap(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)
	mutateStrategyResults := solver.getMutationStrategyResultChannel()

	for {
		select {
		case <-solver.quit:
			solver.quit <- true
			return
		default:
			break
		}

		useBestParent := random.Intn(100) == 0
		parent := solver.pool[0].genes
		if !useBestParent {
			parent = solver.pool[random.Intn(len(solver.pool))].genes
		}

		if len(parent) == 1 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		parentIndexA := random.Intn(len(parent))
		parentIndexB := random.Intn(len(parent))

		for tries := 0; parentIndexA == parentIndexB && tries < 10; parentIndexB = random.Intn(len(parent)) {
			tries++
		}

		parentIndexA, parentIndexB = sort(parentIndexA, parentIndexB)

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parent)))
		if parentIndexA > 0 {
			childGenes.WriteString(parent[:parentIndexA])
		}

		childGenes.WriteString(parent[parentIndexB : parentIndexB+1])

		if parentIndexB-parentIndexA > 1 {
			childGenes.WriteString(parent[parentIndexA+1 : parentIndexB])
		}

		childGenes.WriteString(parent[parentIndexA : parentIndexA+1])

		if parentIndexB+1 < len(parent) {
			childGenes.WriteString(parent[parentIndexB+1:])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy}
		child.fitness = getFitness(child.genes)

		select {
		case strategy.results <- child:
			continue
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) initializeStrategies(numberOfGenesPerChromosome int, getFitness func(string) int) {
	solver.strategies = []strategyInfo{
		{name: "crossover ", start: func(strategyIndex int) {
			solver.crossover(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan sequenceInfo, 1)},
		{name: "mutate    ", start: func(strategyIndex int) {
			solver.mutate(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan sequenceInfo, 1)},
		{name: "reverse   ", start: func(strategyIndex int) {
			solver.reverse(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan sequenceInfo, 1)},
		{name: "shift     ", start: func(strategyIndex int) {
			solver.shift(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan sequenceInfo, 1)},
		{name: "swap      ", start: func(strategyIndex int) {
			solver.swap(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan sequenceInfo, 1)},
	}
	solver.maxStrategySuccess = 1

	for i, _ := range solver.strategies {
		solver.strategies[i].index = i
		go func(index int) { solver.strategies[index].start(index) }(i)
	}
}

func (solver *Solver) inPool(childGenes string) bool {
	solver.poolLock.Lock()
	if solver.distinctPool[childGenes] {
		solver.poolLock.Unlock()
		return true
	}
	solver.distinctPool[childGenes] = true
	solver.poolLock.Unlock()
	return false
}
