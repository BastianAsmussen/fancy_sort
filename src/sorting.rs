use rand::prelude::SliceRandom;

/// Different sorting algorithms implementations.
///
/// # Implementations
/// - [Bubble sort](enum.SortType.html#variant.Bubble)
/// - [Insertion sort](enum.SortType.html#variant.Insertion)
/// - [Selection sort](enum.SortType.html#variant.Selection)
/// - [Merge sort](enum.SortType.html#variant.Merge)
/// - [Quicksort](enum.SortType.html#variant.Quick)
/// - [Heapsort](enum.SortType.html#variant.Heap)
/// - [Shell sort](enum.SortType.html#variant.Shell)
/// - [Cocktail sort](enum.SortType.html#variant.Cocktail)
/// - [Comb sort](enum.SortType.html#variant.Comb)
/// - [Counting sort](enum.SortType.html#variant.Counting)
/// - [Radix sort](enum.SortType.html#variant.Radix)
/// - [Bogo sort](enum.SortType.html#variant.Bogo)
/// - [Stalin sort](enum.SortType.html#variant.Stalin)
pub enum SortType {
    /// Bubble sort implementation.
    /// * Runtime: O(n^2)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Bubble sort, sometimes referred to as sinking sort, is a simple sorting algorithm that repeatedly steps through the list,
    /// compares adjacent elements and swaps them if they are in the wrong order.
    /// The pass through the list is repeated until the list is sorted.
    /// The algorithm, which is a comparison sort, is named for the way smaller or larger elements "bubble" to the top of the list.
    /// Although the algorithm is simple, it is too slow and impractical for most problems even when compared to insertion sort.
    /// Bubble sort can be practical if the input is in mostly sorted order with some out-of-order elements nearly in position.
    ///
    /// # Logic
    /// 1. Iterate through the list of elements.
    /// 2. Compare each element with the next one.
    /// 3. If the current element is greater than the next one, swap them.
    /// 4. Repeat until the list is sorted.
    /// 5. Return the sorted list.
    Bubble,
    /// Insertion sort implementation.
    /// * Runtime: O(n^2)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Insertion sort is a simple sorting algorithm that builds the final sorted array (or list) one item at a time.
    /// It is much less efficient on large lists than more advanced algorithms such as quicksort, heapsort, or merge sort.
    ///
    /// # Logic
    /// 1. Iterate through the list of elements.
    /// 2. Compare each element with the previous one.
    /// 3. If the current element is smaller than the previous one, swap them.
    /// 4. Repeat until the list is sorted.
    /// 5. Return the sorted list.
    Insertion,
    /// Selection sort implementation.
    /// * Runtime: O(n^2)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Selection sort is an in-place comparison sorting algorithm.
    /// It has an O(n2) time complexity, which makes it inefficient on large lists,
    /// and generally performs worse than the similar insertion sort.
    /// Selection sort is noted for its simplicity, and it has performance advantages over more complicated algorithms in certain situations,
    /// particularly where auxiliary memory is limited.
    /// The algorithm divides the input list into two parts: a sorted sublist of items which is built up from left to right at the front (left) of the list
    /// and a sublist of the remaining unsorted items that occupy the rest of the list.
    /// Initially, the sorted sublist is empty and the unsorted sublist is the entire input list.
    /// The algorithm proceeds by finding the smallest (or largest, depending on sorting order) element in the unsorted sublist,
    /// exchanging (swapping) it with the leftmost unsorted element (putting it in sorted order),
    /// and moving the sublist boundaries one element to the right.
    ///
    /// # Logic
    /// 1. Iterate through the list of elements.
    /// 2. Find the smallest element in the list.
    /// 3. Swap it with the first element.
    /// 4. Repeat until the list is sorted.
    /// 5. Return the sorted list.
    Selection,
    /// Merge sort implementation.
    /// * Runtime: O(n log n)
    /// * Memory: O(n)
    ///
    /// # Info
    /// Merge sort is an efficient, general-purpose, comparison-based sorting algorithm.
    /// Most implementations produce a stable sort, which means that the order of equal elements is the same in the input and output.
    /// Merge sort is a divide and conquer algorithm that was invented by John von Neumann in 1945.
    /// A detailed description and analysis of bottom-up mergesort appeared in a report by Goldstine and von Neumann as early as 1948.
    ///
    /// # Logic
    /// 1. Divide the unsorted list into n sublists, each containing one element (a list of one element is considered sorted).
    /// 2. Repeatedly merge sublists to produce new sorted sublists until there is only one sublist remaining. This will be the sorted list.
    /// 3. Return the sorted list.
    Merge,
    /// Quick sort implementation.
    /// * Runtime: O(n log n)
    /// * Memory: O(log n)
    ///
    /// # Info
    /// Quicksort (sometimes called partition-exchange sort) is an efficient sorting algorithm,
    /// serving as a systematic method for placing the elements of an array in order.
    /// Developed by Tony Hoare in 1959 and published in 1961, it is still a commonly used algorithm for sorting.
    /// When implemented well, it can be about two or three times faster than its main competitors, merge sort and heapsort.
    ///
    /// # Logic
    /// 1. Pick an element, called a pivot, from the list.
    /// 2. Reorder the list so that all elements with values less than the pivot come before the pivot,
    /// while all elements with values greater than the pivot come after it (equal values can go either way).
    /// After this partitioning, the pivot is in its final position.
    /// This is called the partition operation.
    /// 3. Recursively apply the above steps to the sub-list of elements with smaller values and separately
    /// to the sub-list of elements with greater values.
    /// 4. Return the sorted list.
    Quick,
    /// Heap sort implementation.
    /// * Runtime: O(n log n)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Heapsort is a comparison-based sorting algorithm.
    /// Heapsort can be thought of as an improved selection sort: like that algorithm,
    /// it divides its input into a sorted and an unsorted region, and it iteratively shrinks the unsorted region by extracting the largest element
    /// and moving that to the sorted region.
    /// The improvement consists of the use of a heap data structure rather than a linear-time search to find the maximum.
    ///
    /// # Logic
    /// 1. Call the build_max_heap() function on the list.
    /// Also referred to as heapify(), this builds a heap from a list in O(n) operations.
    /// 2. Swap the first element of the list with the final element.
    /// Decrease the considered range of the list by one.
    /// 3. Call the sift_down() function on the list to sift the new first element to its appropriate index in the heap.
    /// 4. Go to step 2 unless the considered range of the list is one element.
    /// 5. Return the sorted list.
    Heap,
    /// Shell sort implementation.
    /// * Runtime: O(n log n)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Shellsort, also known as Shell sort or Shell's method, is an in-place comparison sort.
    /// It can be seen as either a generalization of sorting by exchange (bubble sort) or sorting by insertion (insertion sort).
    /// The method starts by sorting pairs of elements far apart from each other, then progressively reducing the gap between elements to be compared.
    /// Starting with far apart elements, it can move some out-of-place elements into position faster than a simple nearest neighbor exchange.
    /// Donald Shell published the first version of this sort in 1959.
    /// The running time of Shellsort is heavily dependent on the gap sequence it uses.
    /// For many practical variants, determining their time complexity remains an open problem.
    ///
    /// # Logic
    /// 1. Rearrange elements at each n/2, n/4, n/8, ... intervals.
    /// 2. Repeat until the interval is 1.
    /// 3. Return the sorted list.
    Shell,
    /// Cocktail sort implementation.
    /// * Runtime: O(n^2)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Cocktail shaker sort is a variation of bubble sort that is both a stable sorting algorithm and a comparison sort.
    /// The algorithm differs from a bubble sort in that it sorts in both directions on each pass through the list.
    /// This sorting algorithm is only marginally more difficult to implement than a bubble sort,
    /// and solves the problem of turtles in bubble sorts.
    /// It provides only marginal performance improvements, and does not improve asymptotic performance;
    /// like the bubble sort, it is not of practical interest (insertion sort is preferred for simple sorts),
    /// though it finds some use in education.
    /// It performs slightly better than bubble sort, and is comparable to insertion sort.
    ///
    /// # Logic
    /// 1. Iterate through the list of elements.
    /// 2. If the current element is greater than the next one, swap them.
    /// 3. If the current element is less than the previous one, swap them.
    /// 4. Repeat until the list is sorted.
    /// 5. Return the sorted list.
    Cocktail,
    /// Comb sort implementation.
    /// * Runtime: O(n^2)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Comb sort is a relatively simple sorting algorithm originally designed by Włodzimierz Dobosiewicz in 1980.
    /// Later it was rediscovered by Stephen Lacey and Richard Box in 1991.
    /// Comb sort improves on bubble sort.
    /// The basic idea is to eliminate turtles, or small values near the end of the list,
    /// since in a bubble sort these slow the sorting down tremendously.
    /// Rabbits, large values around the beginning of the list, do not pose a problem in bubble sort.
    /// In bubble sort, when any two elements are compared, they always have a gap (distance from each other) of 1.
    /// The basic idea of comb sort is that the gap can be much more than 1.
    /// The inner loop of bubble sort, which does the actual swap, is modified such that gap between swapped elements goes down (for each iteration of outer loop) in steps of a "shrink factor" k: [ n/k, n/k2, n/k3, ..., 1 ].
    /// Unlike in bubble sort, where the gap is constant, here the gap is a decreasing function of i.
    ///
    /// # Logic
    /// 1. Initialize gap size for swapping.
    /// 2. While gap size is greater than 1, iterate through the list of elements.
    /// 3. If the current element is greater than the next one, swap them.
    /// 4. Repeat until the list is sorted.
    /// 5. Return the sorted list.
    Comb,
    /// Counting sort implementation.
    /// * Runtime: O(n + k)
    /// * Memory: O(k) where k is the range of the non-negative key values.
    ///
    /// # Info
    /// It is not a comparison sort.
    /// In computer science, counting sort is an algorithm for sorting a collection of objects according to keys that are small integers;
    /// that is, it is an integer sorting algorithm.
    /// It operates by counting the number of objects that have each distinct key value,
    /// and using arithmetic on those counts to determine the positions of each key value in the output sequence.
    /// Its running time is linear in the number of items and the difference between the maximum and minimum key values,
    /// so it is only suitable for direct use in situations where the variation in keys is not significantly greater than the number of items.
    /// However, it is often used as a subroutine in another sorting algorithm, radix sort,
    /// that can handle larger keys more efficiently.
    /// Because counting sort uses key values as indexes into an array, it is not a comparison sort,
    /// and the Ω(n log n) lower bound for comparison sorting does not apply to it.
    /// Bucket sort may be used for many of the same tasks as counting sort,
    /// with a similar time analysis; however, compared to counting sort, bucket sort requires linked lists,
    /// dynamic arrays or a large amount of preallocated memory to hold the sets of items within each bucket,
    /// whereas counting sort instead stores a single number (the count of items) per bucket.
    ///
    /// # Logic
    /// 1. Find the maximum value in the list.
    /// 2. Initialize an array of the size of the maximum value plus one.
    /// 3. Iterate through the list of elements.
    /// 4. Increment the value at the index of the current element in the array.
    /// 5. Iterate through the array.
    /// 6. If the current element is greater than the next one, swap them.
    /// 7. Repeat until the list is sorted.
    /// 8. Return the sorted list.
    Counting,
    /// Radix sort implementation.
    /// * Runtime: O(nk)
    /// * Memory: O(n + k) where k is the length of the longest key.
    ///
    /// # Info
    /// It is not a comparison sort.
    /// In computer science, radix sort is a non-comparative sorting algorithm.
    /// It avoids comparison by creating and distributing elements into buckets according to their radix.
    /// For elements with more than one significant digit, this bucketing process is repeated for each digit, while preserving the ordering of the prior step,
    /// until all digits have been considered.
    /// For this reason, radix sort has also been called bucket sort and digital sort.
    /// Radix sort can be applied to data that can be sorted lexicographically, be they integers, words, punch cards, playing cards, or the mail.
    ///
    /// # Logic
    /// 1. Find the maximum value in the list.
    /// 2. Initialize an array of the size of the maximum value plus one.
    /// 3. Iterate through the list of elements.
    /// 4. Increment the value at the index of the current element in the array.
    /// 5. Iterate through the array.
    /// 6. If the current element is greater than the next one, swap them.
    /// 7. Repeat until the list is sorted.
    /// 8. Return the sorted list.
    Radix,
    /// Bogo sort implementation.
    /// * Runtime: O((n+1)!)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Bogo sort is a highly ineffective sorting function based on the generate and test paradigm.
    //  The function successively generates permutations of its input until it finds one that is sorted.
    //  It is not useful for sorting, but may be used for educational purposes, to contrast it with more efficient algorithms.
    //  Its average run-time is unbounded, but its worst-case time complexity is O((n+1)!).
    //  The probability of the algorithm terminating on a randomly ordered list approaches zero as the list grows.
    //  Bogo sort is equivalent to generating random permutations of the input list until it finds one that is sorted.
    //
    /// # Logic
    /// 7. Randomize the array.
    /// 8. Check if the array is sorted.
    /// 9. If not, return to step 1.
    Bogo,
    /// Stalin sort implementation.
    /// * Runtime: O(n)
    /// * Memory: O(1)
    ///
    /// # Info
    /// Stalin Sort is an efficient sorting algorithm,
    /// serving as a systematic method
    /// for placing the elements of a random access file or an array in order.
    /// Stalin Sort is also know as the best sorting algorithm of all times because of its AMAZING capacity
    /// of always ordering an array with an O(n)
    /// performance.
    ///
    /// # Logic
    /// 1. Iterate through the list of elements.
    /// 2. If the current element is greater than the previous one: keep the element.
    /// 3. If not: remove the element.
    /// 4. Repeat until the list is sorted.
    /// 5. Return the sorted list.
    Stalin,
}

/// Sorts the given list of numbers using the given algorithm.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
/// * `algorithm` - The algorithm to use for sorting.
///
/// # Returns
/// The sorted list of numbers.
pub fn sort(numbers: &[f64], algorithm: SortType) -> Vec<f64> {
    match algorithm {
        SortType::Bubble => bubble_sort(numbers),
        SortType::Insertion => insertion_sort(numbers),
        SortType::Selection => selection_sort(numbers),
        SortType::Shell => shell_sort(numbers),
        SortType::Merge => merge_sort(numbers),
        SortType::Quick => quick_sort(numbers),
        SortType::Heap => heap_sort(numbers),
        SortType::Cocktail => cocktail_sort(numbers),
        SortType::Comb => comb_sort(numbers),
        SortType::Counting => counting_sort(numbers),
        SortType::Radix => radix_sort(numbers),
        SortType::Bogo => bogo_sort(numbers),
        SortType::Stalin => stalin_sort(numbers),
    }
}

/// Bubble sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn bubble_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();
    let mut swapped = true;
    let mut i = 0;

    while swapped {
        swapped = false;
        i += 1;

        for j in 0..(sorted.len() - i) {
            if sorted[j] > sorted[j + 1] {
                sorted.swap(j, j + 1);
                swapped = true;
            }
        }
    }

    sorted
}

/// Insertion sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn insertion_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();

    for i in 1..sorted.len() {
        let mut j = i;

        while j > 0 && sorted[j - 1] > sorted[j] {
            sorted.swap(j - 1, j);
            j -= 1;
        }
    }

    sorted
}

/// Selection sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn selection_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();

    for i in 0..(sorted.len() - 1) {
        let mut min = i;

        for j in (i + 1)..sorted.len() {
            if sorted[j] < sorted[min] {
                min = j;
            }
        }

        if min != i {
            sorted.swap(i, min);
        }
    }

    sorted
}

/// Shell sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn shell_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();
    let mut gap = sorted.len() / 2;

    while gap > 0 {
        for i in gap..sorted.len() {
            let mut j = i;

            while j >= gap && sorted[j - gap] > sorted[j] {
                sorted.swap(j - gap, j);
                j -= gap;
            }
        }

        gap /= 2;
    }

    sorted
}

/// Merge sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn merge_sort(numbers: &[f64]) -> Vec<f64> {
    let sorted = numbers.to_vec();

    if sorted.len() <= 1 {
        return sorted;
    }

    let mid = sorted.len() / 2;
    let left = merge_sort(&sorted[..mid]);
    let right = merge_sort(&sorted[mid..]);

    merge(&left, &right)
}

/// Merge helper function for merge sort.
///
/// # Arguments
/// * `left` - The left list of numbers to merge.
/// * `right` - The right list of numbers to merge.
///
/// # Returns
/// The merged list of numbers.
fn merge(left: &[f64], right: &[f64]) -> Vec<f64> {
    let mut merged = Vec::with_capacity(left.len() + right.len());
    let mut i = 0;
    let mut j = 0;

    while i < left.len() && j < right.len() {
        if left[i] < right[j] {
            merged.push(left[i]);
            i += 1;
        } else {
            merged.push(right[j]);
            j += 1;
        }
    }

    if i < left.len() {
        merged.extend_from_slice(&left[i..]);
    }

    if j < right.len() {
        merged.extend_from_slice(&right[j..]);
    }

    merged
}

/// Quick sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn quick_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();

    if sorted.len() <= 1 {
        return sorted;
    }

    let pivot = sorted.pop().unwrap();
    let mut left = Vec::new();
    let mut right = Vec::new();

    for number in sorted {
        if number < pivot {
            left.push(number);
        } else {
            right.push(number);
        }
    }

    let mut sorted_left = quick_sort(&left);
    let mut sorted_right = quick_sort(&right);

    sorted_left.push(pivot);
    sorted_left.append(&mut sorted_right);

    sorted_left
}

/// Heap sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn heap_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();
    let len = sorted.len();

    for i in (0..len / 2).rev() {
        heapify(&mut sorted, len, i);
    }

    for i in (1..len).rev() {
        sorted.swap(0, i);
        heapify(&mut sorted, i, 0);
    }

    sorted
}

/// Heapify helper function for heap sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to heapify.
/// * `len` - The length of the list of numbers.
/// * `i` - The index of the list of numbers.
fn heapify(numbers: &mut [f64], len: usize, i: usize) {
    let mut largest = i;
    let left = 2 * i + 1;
    let right = 2 * i + 2;

    if left < len && numbers[left] > numbers[largest] {
        largest = left;
    }

    if right < len && numbers[right] > numbers[largest] {
        largest = right;
    }

    if largest != i {
        numbers.swap(i, largest);
        heapify(numbers, len, largest);
    }
}

/// Cocktail shaker sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn cocktail_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();
    let mut swapped = true;
    let mut start = 0;
    let mut end = sorted.len() - 1;

    while swapped {
        swapped = false;

        for i in start..end {
            if sorted[i] > sorted[i + 1] {
                sorted.swap(i, i + 1);
                swapped = true;
            }
        }

        if !swapped {
            break;
        }

        swapped = false;
        end -= 1;

        for i in (start..end).rev() {
            if sorted[i] > sorted[i + 1] {
                sorted.swap(i, i + 1);
                swapped = true;
            }
        }

        start += 1;
    }

    sorted
}

/// Comb sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn comb_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();
    let mut gap = sorted.len();
    let shrink = 1.3;
    let mut swapped = true;

    while swapped || gap > 1 {
        gap = (gap as f64 / shrink) as usize;

        if gap < 1 {
            gap = 1;
        }

        swapped = false;

        for i in 0..sorted.len() - gap {
            if sorted[i] > sorted[i + gap] {
                sorted.swap(i, i + gap);
                swapped = true;
            }
        }
    }

    sorted
}

/// Counting sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn counting_sort(numbers: &[f64]) -> Vec<f64> {
    let sorted = numbers.to_vec();

    let mut count = vec![0; sorted.len()];
    let mut output = vec![0.0; sorted.len()];

    for i in 0..sorted.len() {
        let mut smaller = 0;

        for j in 0..sorted.len() {
            if sorted[j] < sorted[i] {
                smaller += 1;
            } else if sorted[j] == sorted[i] {
                count[i] += 1;
            }
        }

        for j in 0..count[i] {
            output[smaller + j] = sorted[i];
        }
    }

    output
}

/// Radix sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
pub fn radix_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();
    let mut max = sorted[0];

    for (i, ..) in sorted.iter().enumerate() {
        if sorted[i] > max {
            max = sorted[i];
        }
    }

    let mut exp = 1.0;

    while max / exp > 0.0 {
        counting_sort_by_digit(&mut sorted, exp);
        exp *= 10.0;
    }

    sorted
}

/// Counting sort by digit helper function for radix sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
/// * `exp` - The exponent.
fn counting_sort_by_digit(numbers: &mut [f64], exp: f64) {
    let mut output = vec![0.0; numbers.len()];
    let mut count = vec![0; 10];

    for i in 0..numbers.len() {
        count[((numbers[i] / exp) % 10.0) as usize] += 1;
    }

    for i in 1..10 {
        count[i] += count[i - 1];
    }

    for i in (0..numbers.len()).rev() {
        output[count[((numbers[i] / exp) % 10.0) as usize] - 1] = numbers[i];
        count[((numbers[i] / exp) % 10.0) as usize] -= 1;
    }

    numbers.copy_from_slice(&output[..numbers.len()]);
}

/// Bogo sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
///
/// # Remarks
/// This is a joke algorithm. It is not meant to be used in production.
pub fn bogo_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = numbers.to_vec();

    let mut rng = rand::thread_rng();
    while !is_sorted(&sorted) {
        sorted.shuffle(&mut rng);
    }

    sorted
}

/// Is sorted helper function for bogo sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to check.
///
/// # Returns
/// True if the list of numbers is sorted; otherwise, false.
fn is_sorted(numbers: &[f64]) -> bool {
    for i in 1..numbers.len() {
        if numbers[i] < numbers[i - 1] {
            return false;
        }
    }

    true
}

/// Stalin sort.
///
/// # Arguments
/// * `numbers` - The list of numbers to sort.
///
/// # Returns
/// The sorted list of numbers.
///
/// # Remarks
/// This is a joke algorithm. It is not meant to be used in production.
pub fn stalin_sort(numbers: &[f64]) -> Vec<f64> {
    let mut sorted = vec![numbers[0]];

    // Loop through the numbers in the list, but skip the first one since it is already sorted.
    for i in numbers.iter().skip(1) {
        if i >= sorted.last().unwrap() {
            sorted.push(*i);
        }
    }

    sorted
}
