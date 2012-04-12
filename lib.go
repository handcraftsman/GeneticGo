package genetic

import (
	"math/rand"
	"time"
)

func createRandomNumberGenerator(seed int64) rand.Rand {
	if seed == 0 {
		seed = time.Now().UnixNano()
	}
	return *rand.New(rand.NewSource(seed))
}

func insertionSort(items []sequenceInfo, compare func(sequenceInfo, sequenceInfo) bool, index int) {
	if index < 1 || index > len(items) {
		return
	}
	for i := index; i > 0; i-- {
		if compare(items[i], items[i-1]) {
			items[i], items[i-1] = items[i-1], items[i]
			continue
		}
		break
	}
}

func max(a, b int) int {
	if a >= b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a <= b {
		return a
	}
	return b
}

func reverseArray(a []string) {
	for i, j := 0, len(a)-1; i < j; i, j = i+1, j-1 {
		a[i], a[j] = a[j], a[i]
	}
}

func sort(a, b int) (int, int) {
	if a < b {
		return a, b
	}
	return b, a
}
