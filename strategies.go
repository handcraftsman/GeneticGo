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
	random := createRandomNumberGenerator()
	crossoverStrategyResults := solver.getStrategyResultChannel("crossover")

	for {
		if !solver.isHillClimbing ||
			numberOfGenesPerChromosome > 1 && random.Intn(100) != 0 {
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
	random := createRandomNumberGenerator()
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

func (solver *Solver) flutter(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator()
	for {
		parent := <-solver.randomParent
		parentGenes := (*parent).genes
		chromosomeIndex := chooseWeightedChromosome(len(parentGenes), numberOfGenesPerChromosome, random)

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))

		numberOfGenesToFlutter := 1 + random.Intn(numberOfGenesPerChromosome)
		start := chromosomeIndex
		if numberOfGenesToFlutter < numberOfGenesPerChromosome {
			start += random.Intn(numberOfGenesPerChromosome - numberOfGenesToFlutter + 1)
		}

		if start > 0 {
			childGenes.WriteString(parentGenes[:start])
		}
		anyChanged := false
		for i := 0; i < numberOfGenesToFlutter; i++ {
			modifier := random.Intn(5) - 2
			if modifier == 0 {
				if anyChanged {
					childGenes.WriteString(parentGenes[start+i : start+i+1])
					continue
				}
				modifier++
				anyChanged = true
			}
			geneSetIndex := strings.Index(solver.geneSet, parentGenes[start+i:start+i+1])
			geneSetIndex += modifier
			if geneSetIndex < 0 {
				geneSetIndex += len(solver.geneSet)
			} else if geneSetIndex >= len(solver.geneSet) {
				geneSetIndex -= len(solver.geneSet)
			}
			childGenes.WriteString(solver.geneSet[geneSetIndex : geneSetIndex+1])
		}

		if start+numberOfGenesToFlutter < len(parentGenes) {
			childGenes.WriteString(parentGenes[start+numberOfGenesToFlutter:])
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

func (solver *Solver) mutate(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator()
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

func (solver *Solver) rand(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	for {
		parent := <-solver.randomParent
		parentLen := len((*parent).genes)

		childGenes := bytes.NewBuffer(make([]byte, 0, parentLen))
		length := 0
		for length < parentLen {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case chromosome := <-solver.nextChromosome:
				if len(chromosome) != numberOfGenesPerChromosome {
					return
				}
				childGenes.WriteString(chromosome)
				length += len(chromosome)
			}
		}
		child := sequenceInfo{genes: childGenes.String(), strategy: strategy}
		child.parent = &child

		select {
		case strategy.results <- &child:
		case <-solver.quit:
			solver.quit <- true
			return
		}
	}
}

func (solver *Solver) remove(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator()
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
		index := chooseWeightedChromosome(len(parentGenes), numberOfGenesPerChromosome, random)

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

func (solver *Solver) replace(strategy strategyInfo, numberOfGenesPerChromosome int, getFitness func(string) int) {
	random := createRandomNumberGenerator()
	mutateStrategyResults := solver.getStrategyResultChannel("mutate")

	for {
		parent := <-solver.randomParent
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

		parentGenes := (*parent).genes
		chromosomeIndex := chooseWeightedChromosome(len(parentGenes), numberOfGenesPerChromosome, random)

		childGenes := bytes.NewBuffer(make([]byte, 0, len(parentGenes)))

		numberOfGenesToMutate := 1 + random.Intn(numberOfGenesPerChromosome)
		start := 0
		if numberOfGenesToMutate < numberOfGenesPerChromosome {
			start = random.Intn(numberOfGenesPerChromosome - numberOfGenesToMutate + 1)
		}

		childGenes.WriteString(parentGenes[:chromosomeIndex+start])

		for i := 0; i < numberOfGenesToMutate; i++ {
			select {
			case <-solver.quit:
				solver.quit <- true
				return
			case gene := <-solver.nextGene:
				if len(gene) != 1 {
					return
				}
				childGenes.WriteString(gene)
			}
		}

		if chromosomeIndex+start+numberOfGenesToMutate < len(parentGenes) {
			childGenes.WriteString(parentGenes[chromosomeIndex+start+numberOfGenesToMutate:])
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
	random := createRandomNumberGenerator()
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
	random := createRandomNumberGenerator()
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
	random := createRandomNumberGenerator()
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

const initialStrategySuccess = 3

func (solver *Solver) initializeStrategies(numberOfGenesPerChromosome int, getFitness func(string) int) {
	solver.strategies = []strategyInfo{
		{name: "add       ", start: func(strategyIndex int) {
			solver.add(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "crossover ", start: func(strategyIndex int) {
			solver.crossover(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "flutter   ", start: func(strategyIndex int) {
			solver.flutter(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "mutate    ", start: func(strategyIndex int) {
			solver.mutate(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "random    ", start: func(strategyIndex int) {
			solver.rand(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "remove    ", start: func(strategyIndex int) {
			solver.remove(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "replace   ", start: func(strategyIndex int) {
			solver.replace(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "reverse   ", start: func(strategyIndex int) {
			solver.reverse(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "shift     ", start: func(strategyIndex int) {
			solver.shift(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
		{name: "swap      ", start: func(strategyIndex int) {
			solver.swap(solver.strategies[strategyIndex], numberOfGenesPerChromosome, getFitness)
		}, successCount: initialStrategySuccess, results: make(chan *sequenceInfo, 1)},
	}
	solver.maxStrategySuccess = initialStrategySuccess

	for i, _ := range solver.strategies {
		solver.strategies[i].index = i
		go func(index int) { solver.strategies[index].start(index) }(i)
	}
}

func chooseWeightedChromosome(lenParentGenes, numberOfGenesPerChromosome int, random randomSource) int {
	// prefer chromosomes near the end
	numberOfChromosomes := lenParentGenes / numberOfGenesPerChromosome
	index := lenParentGenes - numberOfGenesPerChromosome
	for ; index > 0 && random.Intn(numberOfChromosomes) != 0; numberOfChromosomes, index = numberOfChromosomes-1, numberOfGenesPerChromosome-1 {
	}
	return index
}
