package main

import (
	"fmt";
	"math/rand";
	"time"
)

func main() {
	rand.Seed(time.Now().UnixNano())
	const genes = " abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!."
	target := "Not all those who wander are lost."
	calc := func (current string) int {
		return calculate(target, current)
	}
	
	start := time.Now()
	
	disp := func (current string) {
		fmt.Print(current)
		fmt.Print("\t")
		fmt.Print(calc(current))
		fmt.Print("\t")
		fmt.Println(time.Since(start))
	}
	
	var best = getBest(calc, disp, genes, len(target))
	println(best)
	
	fmt.Print("Total time: ")
	fmt.Println(time.Since(start))
}

func getBest(calculate func(string) int, disp func(string), genes string, length int) string {

	var bestGenes = generateParent(genes, length)
	value := calculate(bestGenes)
	var bestValue = value
	
	for bestValue < length {
		current := mutateParent(bestGenes, genes)
		value := calculate(current)
		if value > bestValue {
			disp(current)
			bestValue = value
			bestGenes = current
		}
	}

	return bestGenes
}

func mutateParent(parent, genes string) string {
	geneIndex := rand.Intn(len(genes))
	parentIndex := rand.Intn(len(parent))
	current := ""
	if parentIndex > 0 {
		current += parent[:parentIndex]
	}
	current += genes[geneIndex:1+geneIndex]
	if parentIndex+1 < len(parent) {
		current += parent[parentIndex+1:]
	}
	return current
}

func generateParent(genes string, length int) string {
	s := ""
	for i := 0; i < length; i++ {
		index := rand.Intn(len(genes))
		s += genes[index:1+index]
	}
	return s
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