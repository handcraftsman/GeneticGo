package main

import (
	"fmt"
	genetic "github.com/handcraftsman/GeneticGo"
	"math"
	"regexp"
	"strings"
	"time"
)

const regexSpecials = "[]()|?*+"

func main() {
	wanted := []string{"AL", "AK", "AS", "AZ", "AR"}
	unwanted := []string{"AA"}

	geneSet := getUniqueCharacters(wanted) + regexSpecials

	calc := func(candidate string) int {
		return calculate(wanted, unwanted, geneSet, candidate)
	}
	start := time.Now()

	disp := func(candidate string) {
		fmt.Println(candidate,
			"\t",
			calc(candidate),
			"\t",
			time.Since(start))
	}

	var solver = new(genetic.Solver)
	solver.MaxSecondsToRunWithoutImprovement = .3
	solver.LowerFitnessesAreBetter = true
	solver.MaxRoundsWithoutImprovement = 10

	var best = solver.GetBestUsingHillClimbing(calc, disp, geneSet, 10, 1, 0)
	fitness := calc(best)

	if fitness == 0 {
		fmt.Println("\nsolved with: " + best)
	} else {
		fmt.Println("\nfailed to find a solution")
		fmt.Println("consider increasing the following:")
		fmt.Println("\tsolver.MaxSecondsToRunWithoutImprovement")
		fmt.Println("\tsolver.MaxRoundsWithoutImprovement")
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

func calculate(wanted, unwanted []string, geneSet, candidate string) int {
	if !isValidRegex(candidate) {
		return math.MaxInt32
	}

	regex := regexp.MustCompile("^(" + candidate + ")$")
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

func isValidRegex(candidate string) bool {
	if strings.Contains(candidate, "()") || strings.Contains(candidate, "??") {
		return false
	}

	_, err := regexp.Compile(candidate)
	return err == nil
}
