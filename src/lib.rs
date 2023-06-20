mod sorting;

#[cfg(test)]
mod tests {
    use crate::sorting;
    use crate::sorting::SortType;

    pub const UNSORTED_NUMBERS: [f64; 10] = [9.0, 3.0, 1.0, 4.0, 2.0, 8.0, 7.0, 6.0, 5.0, 0.0];
    pub const SORTED_NUMBERS: [f64; 10] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

    #[test]
    fn test_bubble_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Bubble);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_insertion_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Insertion);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_selection_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Selection);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_merge_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Merge);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_quick_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Quick);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_heap_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Heap);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_shell_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Shell);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_cocktail_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Cocktail);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_comb_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Comb);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_counting_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Counting);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_radix_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Radix);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_bogo_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Bogo);

        assert_eq!(sorted, SORTED_NUMBERS);
    }

    #[test]
    fn test_stalin_sort() {
        let sorted = sorting::sort(&UNSORTED_NUMBERS, SortType::Stalin);

        assert_eq!(sorted, [9.0]);
    }
}
