package genetic

import (
	"bytes"
	"strings"
)

func (solver *Solver) getStrategyResultChannel(name string) chan *sequenceInfo {
	strategyResults := solver.strategies[0].results
	for _, strategy := range solver.strategies {
		if strings.Contains(strategy.name, name) {
			strategyResults = strategy.results
			break
		}
	}
	return strategyResults
}

func (solver *Solver) add(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)
	crossoverStrategyResults := solver.getStrategyResultChannel("crossover")

	for {
		if !solver.isHillClimbing ||
			random.Intn(100) != 0 {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-crossoverStrategyResults:
				strategy.results <- child
				continue
			}
		}

		parentA := <-solver.randomParent
		parentAgenes := (*parentA).genes
		parentB := <-solver.randomParent
		parentBgenes := (*parentB).genes
		for parentBgenes == parentAgenes {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case parentB = <-solver.randomParent:
				parentBgenes = (*parentB).genes
			}
		}

		childGenes := parentAgenes + parentBgenes[len(parentBgenes)-numberOfGenesPerChromosome:]

		child := sequenceInfo{genes: childGenes, strategy: strategy, parent: parentA}

		select {
		case strategy.results <- &child:
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) crossover(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)
	mutateStrategyResults := solver.getStrategyResultChannel("mutate")

	for {
		parentA := <-solver.randomParent
		parentAgenes := (*parentA).genes
		parentB := <-solver.randomParent
		parentBgenes := (*parentB).genes

		if len(parentAgenes) == numberOfGenesPerChromosome || len(parentBgenes) == numberOfGenesPerChromosome {
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
		maxLength := min(len(parentAgenes)-destinationStart, len(parentBgenes)-sourceStart) / numberOfGenesPerChromosome * numberOfGenesPerChromosome
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
	mutateOneGene := func(parentGenes string) (string, bool) {
		parentIndex := random.Intn(len(parentGenes))

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if parentIndex > 0 {
			childGenes.WriteString(parentGenes[:parentIndex])
		}

		currentGene := parentGenes[parentIndex : parentIndex+1]

		gene := currentGene
		for gene == currentGene {
			select {
			case <-solver.quit:
				solver.quit <- true
				return "", true
			case gene = <-solver.nextGene:
				if len(gene) != 1 {
					return "", true
				}
			}
		}
		childGenes.WriteString(gene)

		if parentIndex+1 < len(parentGenes) {
			childGenes.WriteString(parentGenes[parentIndex+1:])
		}
		return childGenes.String(), false
	}

	for {
		parent := <-solver.randomParent
		childGenes, quit := mutateOneGene((*parent).genes)
		if quit {
			return
		}
		if random.Intn(2) == 1 {
			childGenes, quit = mutateOneGene(childGenes)
			if quit {
				return
			}
		}

		child := sequenceInfo{genes: childGenes, strategy: strategy, parent: parent}

		select {
		case strategy.results <- &child:
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) remove(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator(solver.RandSeed)
	mutateStrategyResults := solver.getStrategyResultChannel("mutate")
	swapStrategyResults := solver.getStrategyResultChannel("swap")

	for {
		if !solver.isHillClimbing {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-swapStrategyResults:
				strategy.results <- child
				continue
			}
		}

		parent := <-solver.randomParent
		if len(parent.genes) <= numberOfGenesPerChromosome {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		// prefer removing from the end
		parentGenes := (*parent).genes
		index := len(parentGenes) - numberOfGenesPerChromosome
		for ; index > 0; index -= numberOfGenesPerChromosome {
			if random.Intn(len(parentGenes)/numberOfGenesPerChromosome) == 0 {
				break
			}
		}

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if index > 0 {
			childGenes.WriteString(parentGenes[0:index])
		}
		if index < len(parentGenes)-numberOfGenesPerChromosome {
			childGenes.WriteString(parentGenes[index+numberOfGenesPerChromosome:])
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
	mutateStrategyResults := solver.getStrategyResultChannel("mutate")

	for {
		parent := <-solver.randomParent
		parentGenes := (*parent).genes

		if len(parent.genes) == numberOfGenesPerChromosome {
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
	mutateStrategyResults := solver.getStrategyResultChannel("mutate")

	for {
		parent := <-solver.randomParent
		parentGenes := (*parent).genes

		numberOfChromosomesInParent := len(parent.genes) / numberOfGenesPerChromosome
		if numberOfChromosomesInParent < 2 {
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

		segmentStart := random.Intn(numberOfChromosomesInParent - 1)
		if segmentStart > 0 {
			segmentStart--
		}
		segmentCount := 2
		if numberOfChromosomesInParent > 2+segmentStart {
			segmentCount = 2 + random.Intn(numberOfChromosomesInParent-(1+segmentStart))
		}

		segmentOffset := numberOfGenesPerChromosome * segmentStart
		segmentLength := numberOfGenesPerChromosome * segmentCount

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if segmentStart > 0 {
			childGenes.WriteString(parentGenes[0:segmentOffset])
		}
		if shiftRight {
			childGenes.WriteString(parentGenes[segmentOffset+segmentLength-numberOfGenesPerChromosome : segmentOffset+segmentLength])
			childGenes.WriteString(parentGenes[segmentOffset : segmentOffset+segmentLength-numberOfGenesPerChromosome])
		} else {
			childGenes.WriteString(parentGenes[segmentOffset+numberOfGenesPerChromosome : segmentOffset+segmentLength])
			childGenes.WriteString(parentGenes[segmentOffset : segmentOffset+numberOfGenesPerChromosome])
		}
		if segmentOffset+segmentLength < len(parentGenes) {
			childGenes.WriteString(parentGenes[segmentOffset+segmentLength : len(parentGenes)])
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
	mutateStrategyResults := solver.getStrategyResultChannel("mutate")

	for {
		parent := <-solver.randomParent
		parentGenes := (*parent).genes

		swapLength := numberOfGenesPerChromosome
		if random.Intn(2) == 0 {
			swapLength = 1
		}

		if len(parentGenes) == swapLength {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case child := <-mutateStrategyResults:
				strategy.results <- child
				continue
			}
		}

		parentIndexA := random.Intn(len(parentGenes)/swapLength) * swapLength
		parentIndexB := random.Intn(len(parentGenes)/swapLength) * swapLength
		if parentIndexA == parentIndexB {
			parentIndexB += swapLength
			parentIndexB %= len(parentGenes)
		}

		parentIndexA, parentIndexB = sort(parentIndexA, parentIndexB)

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))
		if parentIndexA > 0 {
			childGenes.WriteString(parentGenes[:parentIndexA])
		}

		childGenes.WriteString(parentGenes[parentIndexB : parentIndexB+swapLength])

		if parentIndexB-parentIndexA > swapLength {
			childGenes.WriteString(parentGenes[parentIndexA+swapLength : parentIndexB])
		}

		childGenes.WriteString(parentGenes[parentIndexA : parentIndexA+swapLength])

		if parentIndexB+swapLength < len(parentGenes) {
			childGenes.WriteString(parentGenes[parentIndexB+swapLength:])
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
		{name: "add       ", start: func(strategyIndex int) {
			solver.add(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
		{name: "crossover ", start: func(strategyIndex int) {
			solver.crossover(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
		{name: "mutate    ", start: func(strategyIndex int) {
			solver.mutate(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: 0, results: make(chan *sequenceInfo, 1)},
		{name: "remove    ", start: func(strategyIndex int) {
			solver.remove(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
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
