package main

import (
	"fmt"
	genetic "github.com/handcraftsman/GeneticGo"
	"strconv"
	"strings"
	"time"
)

type Direction struct {
	xdiff int
	ydiff int
}

var North = Direction{0, -1}
var NorthEast = Direction{1, -1}
var East = Direction{1, 0}
var SouthEast = Direction{1, 1}
var South = Direction{0, 1}
var SouthWest = Direction{-1, 1}
var West = Direction{-1, 0}
var NorthWest = Direction{-1, -1}

const boardWidthHeight = 8

func main() {
	genes := ""
	for i := 0; i < boardWidthHeight; i++ {
		genes += strconv.Itoa(i)
	}

	start := time.Now()

	calc := func(current string) int {
		return getFitness(current, boardWidthHeight)
	}

	disp := func(current string) {
		display(current, boardWidthHeight)
		fmt.Print(current)
		fmt.Print("\t")
		fmt.Print(getFitness(current, boardWidthHeight))
		fmt.Print("\t")
		fmt.Println(time.Since(start))
	}

	var solver = new(genetic.Solver)
	solver.MaxSecondsToRunWithoutImprovement = 1

	var best = solver.GetBest(calc, disp, genes, boardWidthHeight, 2)
	disp(best)
	fmt.Print("Total time: ")
	fmt.Println(time.Since(start))
}

func display(current string, boardWidthHeight int) {
	board := convertGenesToBoard(current)
	fmt.Println()
	for y := 0; y < boardWidthHeight; y++ {
		for x := 0; x < boardWidthHeight; x++ {
			key := strconv.Itoa(x) + "," + strconv.Itoa(y)
			if board[key] {
				fmt.Print("Q ")
			} else {
				fmt.Print(". ")
			}
		}
		fmt.Println()
	}
}

func getFitness(current string, boardWidthHeight int) int {
	distinctX := make(map[int]bool)
	distinctY := make(map[int]bool)
	board := convertGenesToBoard(current)

	safeQueens := 0
	for coordinate, _ := range board {
		parts := strings.Split(coordinate, ",")

		x, err := strconv.Atoi(parts[0])
		if err != nil {
			panic(err)
		}
		distinctX[x] = true

		y, err := strconv.Atoi(parts[1])
		if err != nil {
			panic(err)
		}
		distinctY[y] = true

		nextPosition := make(chan string)
		defer func() { nextPosition = nil }()

		quit := false
		go getAttackablePositions(x, y, boardWidthHeight, nextPosition, &quit)

		isValid := true
		for n := range nextPosition {
			if board[n] {
				quit = true
				<-nextPosition
				isValid = false
				break
			}
		}
		if isValid {
			safeQueens++
		}
	}
	fitness := 1000*len(board) + safeQueens*100 + len(distinctX)*len(distinctY)

	return fitness
}

func getAttackablePositions(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generators := []func(x, y, boardWidthHeight int, nextPosition chan string, quit *bool){
		generatePositions(North), generatePositions(NorthEast),
		generatePositions(East), generatePositions(SouthEast),
		generatePositions(South), generatePositions(SouthWest),
		generatePositions(West), generatePositions(NorthWest)}

	for _, generator := range generators {
		if *quit {
			break
		}
		generator(x, y, boardWidthHeight, nextPosition, quit)
	}

	close(nextPosition)
}

func generatePositions(direction Direction) func(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	return func(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
		x += direction.xdiff
		y += direction.ydiff
		for y >= 0 && y < boardWidthHeight && x >= 0 && x < boardWidthHeight {
			if *quit {
				return
			}
			nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
			x += direction.xdiff
			y += direction.ydiff
		}
	}
}

func convertGenesToBoard(genes string) map[string]bool {
	board := make(map[string]bool)
	for i := 0; i < len(genes); i += 2 {
		coordinate := genes[i:i+1] + "," + genes[i+1:i+2]
		board[coordinate] = true
	}
	return board
}
