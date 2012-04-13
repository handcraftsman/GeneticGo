package main

import (
	"fmt"
	genetic "github.com/handcraftsman/GeneticGo"
	"strconv"
	"strings"
	"time"
)

const North, West, South, East, Same int = -1, -1, 1, 1, 0
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
		print(current)
		print("\t")
		print(getFitness(current, boardWidthHeight))
		print("\t")
		fmt.Println(time.Since(start))
	}

	var solver = new(genetic.Solver)
	solver.MaxSecondsToRunWithoutImprovement = 1

	var best = solver.GetBest(calc, disp, genes, boardWidthHeight, 2)
	disp(best)
	print("Total time: ")
	fmt.Println(time.Since(start))
}

func display(current string, boardWidthHeight int) {
	board := convertGenesToBoard(current)
	println()
	for y := 0; y < boardWidthHeight; y++ {
		for x := 0; x < boardWidthHeight; x++ {
			key := strconv.Itoa(x) + "," + strconv.Itoa(y)
			if board[key] {
				print("Q ")
			} else {
				print(". ")
			}
		}
		println()
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
		generatePositionsNorth, generatePositionsNorthEast,
		generatePositionsEast, generatePositionsSouthEast,
		generatePositionsSouth, generatePositionsSouthWest,
		generatePositionsWest, generatePositionsNorthWest}

	for _, generator := range generators {
		if *quit {
			break
		}
		generator(x, y, boardWidthHeight, nextPosition, quit)
	}

	close(nextPosition)
}

/*
func generatePositions(direction, x, y, boardWidthHeight int, nextPosition chan string, quit *bool)
	return func(x, y, boardWidthHeight, nextPosition, quit) {
		if direction ==
		generatePositions(x, y, )
	}
}
*/

func generatePositionsNorth(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, North, Same, boardWidthHeight, nextPosition, quit)
}

func generatePositionsSouth(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, South, Same, boardWidthHeight, nextPosition, quit)
}

func generatePositionsWest(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, Same, West, boardWidthHeight, nextPosition, quit)
}

func generatePositionsEast(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, Same, East, boardWidthHeight, nextPosition, quit)
}

func generatePositionsNorthEast(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, North, East, boardWidthHeight, nextPosition, quit)
}

func generatePositionsNorthWest(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, North, West, boardWidthHeight, nextPosition, quit)
}

func generatePositionsSouthEast(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, South, East, boardWidthHeight, nextPosition, quit)
}

func generatePositionsSouthWest(x, y, boardWidthHeight int, nextPosition chan string, quit *bool) {
	generatePositions(x, y, South, West, boardWidthHeight, nextPosition, quit)
}

func generatePositions(x, y, yDifference, xDifference, boardWidthHeight int, nextPosition chan string, quit *bool) {
	x += xDifference
	y += yDifference
	for y >= 0 && y < boardWidthHeight && x >= 0 && x < boardWidthHeight {
		if *quit {
			return
		}
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
		x += xDifference
		y += yDifference
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
