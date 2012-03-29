package main

import (
	"fmt";
	"strconv";
	"strings";
	"time";
	"./genetic"
)

func main() {
	const boardWidthHeight = 8
	genes := ""
	for i := 0; i < boardWidthHeight; i++ {
		genes += strconv.Itoa(i)
	}
	
	start := time.Now()
	
	calc := func (current string) int {
		return getFitness(current, boardWidthHeight)
	}
	
	disp := func (current string) {
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
	for y := 0; y < boardWidthHeight; y++ {
		for x:= 0; x < boardWidthHeight; x++ {
			key := strconv.Itoa(x) + "," + strconv.Itoa(y)
			if board[key] {
				fmt.Print("Q ")
			} else {
				fmt.Print(". ")
			}
		}
		fmt.Println("")
	}
}

func getFitness(current string, boardWidthHeight int) int {
	fitness := 0
	board := convertGenesToBoard(current)
	for coordinate, _ := range ( board ) {
	
		nextPosition := make(chan string)
		go getAttackablePositions(coordinate, boardWidthHeight, nextPosition)
		n := <- nextPosition
        for {
			if len(n) == 0 {
				fitness++
				break;
			}
			if board[n] {
				break
			}
			n = <- nextPosition
		}
		nextPosition = nil
    }
	fitness += 10 * len(board)
	return fitness
}

func getAttackablePositions(coordinate string, boardWidthHeight int, nextPosition chan string) {
	parts := strings.Split(coordinate, ",")
	x, err := strconv.Atoi(parts[0])
	if err != nil { panic(err) }
	y, err := strconv.Atoi(parts[1])
	if err != nil { panic(err) }
	getPositionsNorth(x, y, boardWidthHeight, nextPosition)
	getPositionsNorthEast(x, y, boardWidthHeight, nextPosition)
	getPositionsEast(x, y, boardWidthHeight, nextPosition)
	getPositionsSouthEast(x, y, boardWidthHeight, nextPosition)
	getPositionsSouth(x, y, boardWidthHeight, nextPosition)
	getPositionsSouthWest(x, y, boardWidthHeight, nextPosition)
	getPositionsWest(x, y, boardWidthHeight, nextPosition)
	getPositionsNorthWest(x, y, boardWidthHeight, nextPosition)
	nextPosition <- ""
}

func getPositionsNorth(x, y, boardWidthHeight int, nextPosition chan string) {
	for y-=1 ; y >= 0; y-- {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
	}
}

func getPositionsSouth(x, y, boardWidthHeight int, nextPosition chan string) {
	for y+=1 ; y < boardWidthHeight; y++ {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
	}
}

func getPositionsWest(x, y, boardWidthHeight int, nextPosition chan string) {
	for x-=1 ; x >= 0; x-- {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
	}
}

func getPositionsEast(x, y, boardWidthHeight int, nextPosition chan string) {
	for x+=1 ; x < boardWidthHeight; x++ {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
	}
}

func getPositionsNorthEast(x, y, boardWidthHeight int, nextPosition chan string) {
	y--
	x++
	for ; y >= 0 && x < boardWidthHeight; y-- {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
		x++
	}
}

func getPositionsNorthWest(x, y, boardWidthHeight int, nextPosition chan string) {
	y--
	x--
	for ; y >= 0 && x >= 0; y-- {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
		x--
	}
}

func getPositionsSouthEast(x, y, boardWidthHeight int, nextPosition chan string) {
	y++
	x++
	for ; y < boardWidthHeight && x < boardWidthHeight; y++ {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
		x++ 
	}
}

func getPositionsSouthWest(x, y, boardWidthHeight int, nextPosition chan string) {
	y++ 
	x--
	for ; y < boardWidthHeight && x >= 0; y++ {
		nextPosition <- strconv.Itoa(x) + "," + strconv.Itoa(y)
		x--
	}
}

func convertGenesToBoard(genes string) map[string] bool {
	board := make(map[string]bool)
	for i := 0; i < len(genes); i+=2 {
		coordinate := genes[i:i+1] + "," + genes[i+1:i+2]
		board[coordinate] = true
	}
	return board
}


