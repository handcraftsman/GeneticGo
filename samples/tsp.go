package main

import (
	"../genetic"
	"flag"
	"fmt"
	"github.com/handcraftsman/File"
	"math"
	"strconv"
	"strings"
	"time"
)

const genericGeneSet string = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"

func main() {
	flag.Parse()
	if flag.NArg() != 1 {
		println("Usage: go run samples/tsp.go ROUTEFILEPATH")
		return
	}
	var routeFileName = flag.Arg(0)
	if !File.Exists(routeFileName) {
		println("file " + routeFileName + " does not exist.")
		return
	}
	println("using route file: " + routeFileName)

	idToPointLookup := readPoints(routeFileName)
	println("read " + strconv.Itoa(len(idToPointLookup)) + " points...")

	calc := func(genes string) int {
		points := genesToPoints(genes, idToPointLookup)
		return getFitness(genes, points)
	}

	if File.Exists(routeFileName + ".opt.tour") {
		println("found optimal solution file: " + routeFileName + ".opt")
		optimalRoute := readOptimalRoute(routeFileName+".opt.tour", len(idToPointLookup))
		println("read " + strconv.Itoa(len(optimalRoute)) + " segments in the optimal route")
		points := getPointsInOptimalOrder(idToPointLookup, optimalRoute)
		genes := genericGeneSet[0:len(idToPointLookup)]
		idToPointLookup = make(map[string]Point, len(idToPointLookup))
		for i, v := range points {
			idToPointLookup[genericGeneSet[i:i+1]] = v
		}
		print("optimal route: " + genes)
		print("\t")
		println(getFitness(genes, points))
	}

	genes := genericGeneSet[0:len(idToPointLookup)]

	start := time.Now()

	disp := func(genes string) {
		points := genesToPoints(genes, idToPointLookup)
		print(genes)
		print("\t")
		print(getFitness(genes, points))
		print("\t")
		fmt.Println(time.Since(start))
	}

	var solver = new(genetic.Solver)
	solver.MaxSecondsToRunWithoutImprovement = 20
	solver.LowerFitnessesAreBetter = true

	var best = solver.GetBest(calc, disp, genes, len(idToPointLookup), 1)
	disp(best)
	print("Total time: ")
	fmt.Println(time.Since(start))
}

func genesToPoints(genes string, idToPointLookup map[string]Point) []Point {
	points := make([]Point, len(idToPointLookup))
	if len(genes) != len(idToPointLookup) {
		panic(genes + " - incorrect gene sequence length")
	}
	for i := 0; i < len(genes); i++ {
		geneId := genes[i : i+1]
		point := idToPointLookup[geneId]
		points[i] = point
	}
	return points
}

func getFitness(genes string, points []Point) int {
	distinctPoints := make(map[string]bool)

	for i := 0; i < len(genes); i++ {
		distinctPoints[genes[i:i+1]] = true
	}

	fitness := getDistance(points[0], points[len(points)-1])
	for i := 0; i < len(points)-1; i++ {
		fitness += getDistance(points[i], points[i+1])
	}
	if len(distinctPoints) != len(genes) {
		fitness += 10000 * (len(genes) - len(distinctPoints))
	}
	return fitness
}

func getDistance(pointA, pointB Point) int {
	sideA := float64(pointA.row - pointB.row)
	sideB := float64(pointA.col - pointB.col)
	sideC := math.Sqrt(sideA*sideA + sideB*sideB)
	return int(.5 + sideC)
}

func convertGenesToBoard(genes string) map[string]bool {
	board := make(map[string]bool)
	for i := 0; i < len(genes); i += 2 {
		coordinate := genes[i:i+1] + "," + genes[i+1:i+2]
		board[coordinate] = true
	}
	return board
}

type Point struct {
	row int
	col int
}

func getPointsInOptimalOrder(idToPointLookup map[string]Point, optimalRoute []string) []Point {
	points := make([]Point, len(optimalRoute))
	i := 0
	for _, v := range optimalRoute {
		points[i] = idToPointLookup[v]
		i++
	}
	return points
}

func readPoints(routeFileName string) map[string]Point {
	pointLines := make(chan string)

	lineHandler := tspRouteFileHeader
	go func() {
		for line := range File.EachLine(routeFileName) {
			lineHandler = lineHandler(line, pointLines)
		}
		close(pointLines)
	}()

	points := make(map[string]Point)

	for pointLine := range pointLines {
		parts := strings.Split(pointLine, " ")

		geneIndex, err := strconv.Atoi(parts[0])
		if err != nil {
			panic(err)
		}

		x, err := strconv.Atoi(parts[1])
		if err != nil {
			panic(err)
		}

		y, err := strconv.Atoi(parts[2])
		if err != nil {
			panic(err)
		}

		point := Point{col: x, row: y}

		geneId := genericGeneSet[geneIndex : geneIndex+1]

		points[geneId] = point
	}

	return points
}
func readOptimalRoute(optimalRouteFileName string, numberExpected int) []string {
	pointLines := make(chan string)

	lineHandler := tspOptimalRouteFileHeader
	go func() {
		for line := range File.EachLine(optimalRouteFileName) {
			lineHandler = lineHandler(line, pointLines)
		}
		close(pointLines)
	}()

	pointIds := make([]string, numberExpected)

	i := 0
	for pointLine := range pointLines {
		x, err := strconv.Atoi(pointLine)
		if err != nil {
			panic(err)
		}

		pointIds[i] = genericGeneSet[x : x+1]

		i++
	}

	return pointIds
}

func values(m map[string]Point) []Point {
	list := make([]Point, len(m))
	i := 0
	for _, v := range m {
		list[i] = v
		i++
	}
	return list
}

type tspLineFn func(line string, pointLines chan string) tspLineFn

func tspRouteFileHeader(line string, pointLines chan string) tspLineFn {
	if line != "NODE_COORD_SECTION" {
		return tspRouteFileHeader
	}
	return tspRoutePoint
}

func tspRoutePoint(line string, pointLines chan string) tspLineFn {
	if line != "EOF" {
		pointLines <- line
		return tspRoutePoint
	}
	return tspRouteFileDone
}

func tspRouteFileDone(line string, pointLines chan string) tspLineFn {
	return tspRouteFileDone
}

func tspOptimalRouteFileHeader(line string, pointLines chan string) tspLineFn {
	if line != "TOUR_SECTION" {
		return tspOptimalRouteFileHeader
	}
	return tspOptimalRoutePoint
}

func tspOptimalRoutePoint(line string, pointLines chan string) tspLineFn {
	if line != "-1" {
		pointLines <- line
		return tspOptimalRoutePoint
	}
	return tspRouteFileDone
}
