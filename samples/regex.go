package main

import (
	".."
	"fmt"
	"math"
	"regexp"
	"strings"
	"time"
)

const regexSpecials = "[]()|?*+"

func main() {
	wanted := []string{"AL","AK","AS","AZ","AR"}
	unwanted := []string{"AA"}

	geneSet := getUniqueCharacters(wanted) + regexSpecials

	calc := func(genes string) int {
		return calculate(wanted, unwanted, geneSet, genes)
	}
	start := time.Now()

	disp := func(genes string) {
		fmt.Println(genes,
		"\t",
		calc(genes),
		"\t",
		time.Since(start))
	}

	var solver = new(genetic.Solver)
	solver.MaxSecondsToRunWithoutImprovement = 3
	solver.LowerFitnessesAreBetter = true

	var best = solver.GetBestUsingHillClimbing(calc, disp, geneSet, 10, 1, 0)
	fitness := calc(best)

	if fitness == 0 {
		println("\nsolved with: " + best)
	} else {
		println("\nfailed to find a solution")
	}

	fmt.Print("Total time: ")
	fmt.Println(time.Since(start))
}

func getUniqueCharacters(wanted []string) string {
	uniqueCharacters := make(map[string]bool)

	characters := ""
	for _, item := range wanted {
		for i := 0; i < len(item); i++ {
			token := item[i : i+1]
			if !uniqueCharacters[token] {
				characters += token
				uniqueCharacters[token] = true
			}
		}
	}
	return characters
}

func calculate(wanted, unwanted []string, geneSet string, genes string) int {
	if !isValidRegex(genes) {
		return math.MaxInt32
	}

	regex := regexp.MustCompile("^(" + genes + ")$")
	fitness := 0
	for _, item := range wanted {
		if !regex.MatchString(item) {
			fitness++
		}
	}

	for _, item := range unwanted {
		if regex.MatchString(item) {
			fitness += 10
		}
	}

	return fitness
}

func isValidRegex(genes string) bool {
	if strings.Contains(genes, "()") || strings.Contains(genes, "??") {
		return false
	}

	_, err := regexp.Compile(genes)
	return err == nil
}
