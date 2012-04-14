package genetic

import (
	"math/rand"
	"sync"
)

type Solver struct {
	RandSeed                          int64
	MaxSecondsToRunWithoutImprovement float64
	MaxRoundsWithoutImprovement       int
	LowerFitnessesAreBetter           bool
	PrintStrategyUsage                bool
	PrintDiagnosticInfo               bool

	childFitnessIsBetter, childFitnessIsSameOrBetter func(child, other sequenceInfo) bool

	quit                     chan bool
	nextGene, nextChromosome chan string

	strategies         []strategyInfo
	maxStrategySuccess int

	needNewlineBeforeDisplay bool

	maxPoolSize  int
	pool         []sequenceInfo
	distinctPool map[string]bool
	poolLock     sync.Mutex

	rand rand.Rand
}

type sequenceInfo struct {
	genes    string
	fitness  int
	strategy strategyInfo
}

type strategyInfo struct {
	name         string
	start        func(strategyIndex int)
	successCount int
	results      chan sequenceInfo
	index        int
}
