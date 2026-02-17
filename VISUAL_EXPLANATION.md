# Visual Explanation of the LIANA SCE Fix

## The Problem Visualized

### Scenario: Original data has 10 cells, we filter to 5 cells

```
ORIGINAL DATA (10 cells):
Cells:   [1] [2] [3] [4] [5] [6] [7] [8] [9] [10]
Labels:   A   B   A   C   A   B   C   A   B   C
Counts:   10  0   5   8   0   12  0   9   0   6
                ↓
Filter: Keep only cells with counts > 0
                ↓
FILTERED DATA (5 cells):
Cells:   [1]     [3] [4]     [6]     [8]     [10]
Labels:   A       A   C       B       A       C
Counts:   10      5   8       12      9       6
```

### BUGGY APPROACH ❌

```r
# Step 1: Create boolean filter from FILTERED data
nonzero_cells <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
#                  [1]   [2]   [3]   [4]   [5]   [6]   [7]   [8]   [9]   [10]

# Step 2: Filter the matrix
filtered_cnt <- original_cnt[, nonzero_cells]
# Result: 5 cells (indices 1, 3, 4, 6, 8, 10 from original)

# Step 3: BUG! Apply the SAME boolean filter to ORIGINAL labels
buggy_labels <- original_labels[nonzero_cells]
# This takes ORIGINAL positions 1, 3, 4, 6, 8, 10
# Which gives: A, A, C, B, A, C ✓ (happens to work in this example)

# BUT if we filter again...
filtered_cnt2 <- filtered_cnt[rowSums(filtered_cnt) > 0, ]  # More filtering
nonzero_cells2 <- colSums(filtered_cnt2) > 0  # Boolean from FILTERED data

# Step 4: BUG APPEARS! Apply boolean from FILTERED data to ORIGINAL labels
buggy_labels2 <- original_labels[nonzero_cells2]
# nonzero_cells2 is length 5 (from filtered data)
# But it's used to index original_labels (length 10)
# DIMENSION MISMATCH! ❌
```

### FIXED APPROACH ✓

```r
# Step 1: Track original indices
cell_indices <- seq_len(10)  # [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Step 2: Create boolean filter
nonzero_cells <- c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)

# Step 3: Filter the matrix
filtered_cnt <- original_cnt[, nonzero_cells]
# Result: 5 cells

# Step 4: Update tracked indices (KEY FIX!)
cell_indices <- cell_indices[nonzero_cells]
# cell_indices = [1, 3, 4, 6, 8, 10]  # Tracks which original cells remain

# Step 5: Use tracked indices for labels
fixed_labels <- original_labels[cell_indices]
# Takes original positions 1, 3, 4, 6, 8, 10
# Gives: A, A, C, B, A, C ✓

# Step 6: Filter again (second filtering)
nonzero_genes <- rowSums(filtered_cnt) > 0
filtered_cnt2 <- filtered_cnt[nonzero_genes, ]
# cell_indices stays [1, 3, 4, 6, 8, 10] because we didn't filter cells

# If we filter cells again:
nonzero_cells2 <- c(TRUE, TRUE, FALSE, TRUE, TRUE)  # From 5 cells, keep 4
filtered_cnt3 <- filtered_cnt2[, nonzero_cells2]

# Step 7: Update indices again
cell_indices <- cell_indices[nonzero_cells2]
# cell_indices = [1, 3, 6, 8]  # Still tracks original positions

# Step 8: Labels still align perfectly
fixed_labels2 <- original_labels[cell_indices]
# Takes original positions 1, 3, 6, 8
# Gives: A, A, B, A ✓
# Length = 4 = ncol(filtered_cnt3) ✓ PERFECT MATCH!
```

## Step-by-Step Comparison

### BUGGY: Using Boolean Filter Directly

```
Original:     [1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9]  [10]
Labels:        A    B    A    C    A    B    C    A    B    C
                    ↓ FILTER (keep 1,3,4,6,8,10)
Boolean:      [T]  [F]  [T]  [T]  [F]  [T]  [F]  [T]  [F]  [T]
Filtered:     [A]       [A]  [C]       [B]       [A]       [C]
                    ↓ FILTER AGAIN (keep 1,2,4,5)
Boolean:      [T]  [T]  [F]  [T]  [T]
Apply to:      A    B    A    C    A  ← WRONG! These are original labels
Result:        A    B    A    C    A  ← Length 5 but matrix has 4 cells!
                                        DIMENSION MISMATCH! ❌
```

### FIXED: Tracking Indices

```
Original:     [1]  [2]  [3]  [4]  [5]  [6]  [7]  [8]  [9]  [10]
Labels:        A    B    A    C    A    B    C    A    B    C
Indices:       1    2    3    4    5    6    7    8    9    10
                    ↓ FILTER (keep 1,3,4,6,8,10)
Boolean:      [T]  [F]  [T]  [T]  [F]  [T]  [F]  [T]  [F]  [T]
Filtered:     [A]       [A]  [C]       [B]       [A]       [C]
Indices:       1         3    4         6         8         10  ← TRACKED!
                    ↓ FILTER AGAIN (keep 1,2,4,5 of filtered)
Boolean:      [T]  [T]  [F]  [T]  [T]
Indices:       1    3         6    8   ← UPDATED!
Apply to:      A    A         B    A   ← Use tracked indices
Result:        A    A    B    A         ← Length 4 = matrix cells!
                                         PERFECT ALIGNMENT! ✓
```

## The Key Insight

**BUGGY APPROACH:**
- Creates boolean filter from FILTERED data
- Applies boolean filter to ORIGINAL labels
- ❌ Dimensions don't match after multiple filtering steps

**FIXED APPROACH:**
- Tracks which ORIGINAL indices remain
- Updates tracked indices after each filter
- Uses tracked indices to get correct labels
- ✓ Dimensions always match perfectly

## Code Pattern

### Buggy Pattern ❌
```r
# Filter 1
bool1 <- some_condition(data)
data <- data[, bool1]

# Filter 2  
bool2 <- another_condition(data)  # Boolean from FILTERED data
data <- data[, bool2]

# Labels - BUG!
labels <- original_labels[bool2]  # Boolean from FILTERED data on ORIGINAL labels
# Length mismatch! ❌
```

### Fixed Pattern ✓
```r
# Initialize
indices <- seq_len(ncol(data))

# Filter 1
bool1 <- some_condition(data)
data <- data[, bool1]
indices <- indices[bool1]  # UPDATE INDICES

# Filter 2
bool2 <- another_condition(data)
data <- data[, bool2]
indices <- indices[bool2]  # UPDATE INDICES AGAIN

# Labels - CORRECT!
labels <- original_labels[indices]  # Tracked indices always align
# Perfect match! ✓
```

## Summary

The fix is simple but crucial:
1. **Track** original cell indices from the start
2. **Update** tracked indices after EVERY filtering operation
3. **Use** tracked indices (not boolean filters) to extract labels
4. **Verify** dimensions match with `stopifnot()`

This ensures expression matrices and metadata labels remain aligned through all filtering operations, preventing the SingleCellExperiment dimension mismatch error.
