package genetic

import (
	"bytes"
	"strings"
)

func (solver *Solver) getMutationStrategyResultChannel() chan *sequenceInfo {
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

		parentA := <-solver.randomParent
		parentAgenes := (*parentA).genes
		parentB := <-solver.randomParent
		parentBgenes := (*parentB).genes

		if len(parentAgenes) == 1 || len(parentBgenes) == 1 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		sourceStart := random.Intn((len(parentBgenes)-1)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		destinationStart := random.Intn((len(parentAgenes)-1)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		maxLength := min(len(parentAgenes)-destinationStart, len(parentBgenes)-sourceStart)
		length := (1 + random.Intn(maxLength/numberOfGenesPerChromosome-1)) * numberOfGenesPerChromosome

		childGenes := bytes.NewBuffer(make([]byte, 0, max(len(parentAgenes), len(parentBgenes))))

		if destinationStart > 0 {
			childGenes.WriteString(parentAgenes[0:destinationStart])
		}

		childGenes.WriteString(parentBgenes[sourceStart : sourceStart+length])

		if childGenes.Len() < len(parentAgenes) {
			childGenes.WriteString(parentAgenes[childGenes.Len():len(parentAgenes)])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy, parent: parentA}

		select {
		case strategy.results <- &child:
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

		parent := <-solver.randomParent
		parentGenes := (*parent).genes

		parentIndex := random.Intn(len(parentGenes))

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if parentIndex > 0 {
			childGenes.WriteString(parentGenes[:parentIndex])
		}

		currentGene := parentGenes[parentIndex : parentIndex+1]

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

		if parentIndex+1 < len(parentGenes) {
			childGenes.WriteString(parentGenes[parentIndex+1:])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy, parent: parent}

		select {
		case strategy.results <- &child:
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

		parent := <-solver.randomParent
		parentGenes := (*parent).genes

		if len(parent.genes) == 1 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		reversePointA := random.Intn(len(parentGenes)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		reversePointB := random.Intn(len(parentGenes)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		for ; reversePointA == reversePointB; reversePointB = random.Intn(len(parentGenes)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome {
		}

		min, max := sort(reversePointA, reversePointB)

		fragments := make([]string, max-min)
		for i := min; i < max; i += numberOfGenesPerChromosome {
			fragments[i-min] = parentGenes[i : i+numberOfGenesPerChromosome]
		}

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if min > 0 {
			childGenes.WriteString(parentGenes[0:min])
		}

		reverseArray(fragments)
		for _, fragment := range fragments {
			childGenes.WriteString(fragment)
		}

		if childGenes.Len() < len(parentGenes) {
			childGenes.WriteString(parentGenes[childGenes.Len():len(parentGenes)])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy, parent: parent}

		select {
		case strategy.results <- &child:
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

		parent := <-solver.randomParent
		parentGenes := (*parent).genes

		if len(parent.genes) < numberOfGenesPerChromosome+1 {
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

		segmentStart := random.Intn((len(parentGenes)-numberOfGenesPerChromosome)/numberOfGenesPerChromosome) * numberOfGenesPerChromosome
		segmentCount := 1 + random.Intn((len(parentGenes)-segmentStart+numberOfGenesPerChromosome)/numberOfGenesPerChromosome-1)

		// +2 because first and last will be empty to leave room
		fragments := make([]string, segmentCount+2)
		for i := 0; i < segmentCount; i++ {
			start := segmentStart + i*numberOfGenesPerChromosome
			fragments[i+1] = parentGenes[start : start+numberOfGenesPerChromosome]
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

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if segmentStart > 0 {
			childGenes.WriteString(parentGenes[0:segmentStart])
		}

		for i := start; i <= end; i++ {
			childGenes.WriteString(fragments[i])
		}

		if childGenes.Len() < len(parentGenes) {
			childGenes.WriteString(parentGenes[childGenes.Len():len(parentGenes)])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy, parent: parent}

		select {
		case strategy.results <- &child:
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

		parent := <-solver.randomParent
		parentGenes := (*parent).genes

		if len(parentGenes) == 1 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		parentIndexA := random.Intn(len(parentGenes))
		parentIndexB := random.Intn(len(parentGenes))

		for tries := 0; parentIndexA == parentIndexB && tries < 10; parentIndexB = random.Intn(len(parentGenes)) {
			tries++
		}

		parentIndexA, parentIndexB = sort(parentIndexA, parentIndexB)

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if parentIndexA > 0 {
			childGenes.WriteString(parentGenes[:parentIndexA])
		}

		childGenes.WriteString(parentGenes[parentIndexB : parentIndexB+1])

		if parentIndexB-parentIndexA > 1 {
			childGenes.WriteString(parentGenes[parentIndexA+1 : parentIndexB])
		}

		childGenes.WriteString(parentGenes[parentIndexA : parentIndexA+1])

		if parentIndexB+1 < len(parentGenes) {
			childGenes.WriteString(parentGenes[parentIndexB+1:])
		}

		child := sequenceInfo{genes: childGenes.String(), strategy: strategy, parent: parent}

		select {
		case strategy.results <- &child:
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
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
		{name: "mutate    ", start: func(strategyIndex int) {
			solver.mutate(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
		{name: "reverse   ", start: func(strategyIndex int) {
			solver.reverse(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
		{name: "shift     ", start: func(strategyIndex int) {
			solver.shift(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
		{name: "swap      ", start: func(strategyIndex int) {
			solver.swap(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
	}
	solver.maxStrategySuccess = 1

	for i, _ := range solver.strategies {
		solver.strategies[i].index = i
		go func(index int) { solver.strategies[index].start(index) }(i)
	}
}

func (solver *Solver) inPool(childGenes string) bool {
	solver.distinctPoolLock.Lock()
	if solver.distinctPool[childGenes] {
		solver.distinctPoolLock.Unlock()
		return true
	}
	solver.distinctPool[childGenes] = true
	solver.distinctPoolLock.Unlock()
	return false
}
