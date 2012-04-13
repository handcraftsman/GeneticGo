package main

import (
	"fmt"
	genetic "github.com/handcraftsman/GeneticGo"
	"time"
)

func main() {
	const genes = " abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!."
	target := "Not all those who wander are lost."
	calc := func(current string) int {
		return calculate(target, current)
	}

	start := time.Now()

	disp := func(current string) {
		fmt.Print(current)
		fmt.Print("\t")
		fmt.Print(calc(current))
		fmt.Print("\t")
		fmt.Println(time.Since(start))
	}

	var solver = new(genetic.Solver)
	solver.MaxSecondsToRunWithoutImprovement = 1

	var best = solver.GetBest(calc, disp, genes, len(target), 1)
	fmt.Println(best)

	fmt.Print("Total time: ")
	fmt.Println(time.Since(start))
}

func calculate(target, current string) int {
	differenceCount := 0
	for i := 0; i < len(target); i++ {
		if target[i] != current[i] {
			differenceCount++
		}
	}
	return len(target) - differenceCount
}
