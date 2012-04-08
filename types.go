package genetic

type Solver struct {
	RandSeed                          int64
	MaxSecondsToRunWithoutImprovement float64
	LowerFitnessesAreBetter           bool
	PrintStrategyUsage                bool
	PrintDiagnosticInfo               bool
}

type sequenceInfo struct {
	genes   string
	fitness int
}

type strategyInfo struct {
	name                         string
	function                     func(parentA, parentB, geneSet string, numberOfGenesPerChromosome int, nextGene chan string, useBestParent bool) string
	count                        int
	incrementParentSuccessCounts bool
}
