package genetic

import (
	"math/rand"
	"time"
)

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

func min(a, b int) int {
	if a <= b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a >= b {
		return a
	}
	return b
}

func reverseArray(a []string) {
	for i, j := 0, len(a)-1; i < j; i, j = i+1, j-1 {
		a[i], a[j] = a[j], a[i]
	}
}

func seedRandomNumberGenerator(seed int64) {
	if seed > 0 {
		rand.Seed(seed)
	} else {
		rand.Seed(time.Now().UnixNano())
	}
}

func sort(a, b int) (int, int) {
	if a < b {
		return a, b
	}
	return b, a
}
