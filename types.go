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
	randomParent             chan *sequenceInfo

	strategies                     []strategyInfo
	maxStrategySuccess             int
	numberOfImprovements           int
	successParentIsBestParentCount int

	needNewlineBeforeDisplay bool

	maxPoolSize      int
	pool             []sequenceInfo
	poolLock         sync.Mutex
	distinctPool     map[string]bool
	distinctPoolLock sync.Mutex

	random rand.Rand
}

type sequenceInfo struct {
	genes    string
	fitness  int
	strategy strategyInfo
	parent   *sequenceInfo
}

type strategyInfo struct {
	name         string
	start        func(strategyIndex int)
	successCount int
	results      chan *sequenceInfo
	index        int
}
